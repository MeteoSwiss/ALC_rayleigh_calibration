#!/usr/bin/env python3
"""
Download CAMS global atmospheric composition FORECAST aerosol optical profiles
from the Atmosphere Data Store (ADS), compute model-level altitude, and write a
netCDF whose variable names match the existing CAMS_Beta_*.nc files.

Variables retrieved (TOA backscatter intentionally skipped):
    aerext355/532/1064            aerosol extinction coefficient        [m-1]
    aerbackscatgnd355/532/1064    attenuated backscatter from ground    [m-1 sr-1]
    t, q, lnsp, z                 temperature, specific humidity,
                                  ln(surface pressure), surface geopotential

Computed and added:
    z         geopotential on each model level   [m2 s-2]   (overwrites the
              surface-only z that ADS delivers broadcast across levels)
    altitude       geometric height above sea level   [m]
    altitude_agl   geometric height above model orography [m]

Usage
-----
    python download_cams_beta.py 20260601                    # one day -> CAMS_Beta_20260601.nc
    python download_cams_beta.py 20260601 20260605           # inclusive range -> one combined file
    python download_cams_beta.py 202606                      # whole month -> CAMS_Beta_202606.nc
    python download_cams_beta.py --daily 20260601 20260605   # one CAMS_Beta_YYYYMMDD.nc per day
    python download_cams_beta.py --out D:/CAMS 202606         # write into a chosen folder

All forms include the full variable set (aerosol extinction + ground backscatter at
355/532/1064 nm, plus t/q/z/lnsp). The calibration auto-download is the only path that
trims to t/q/z/lnsp; for aerosol work use this CLI (or download_daily_files()).

Notes
-----
* These optical fields exist only as type=forecast (no analysis version), so the
  request always uses type=forecast, the 00 UTC run, steps 3..24 (8 valid times
  per day -> matches the 248 steps in a 31-day file).
* Requests are chunked by day (one ADS request per date) to stay under the ADS
  per-request cost limit, with an automatic variable-split fallback for big days.
* ADS pre-interpolates to a regular 0.4 deg grid and ignores the 'grid' keyword.
  Set REGRID_TO_1DEG = True to bilinearly resample onto the historical 1 deg grid
  (-27..45 E, 27.5..73.5 N) so the series matches older MARS-derived files.
* GRIB -> netCDF uses cfgrib (pure Python; the eccodes binary ships in the pip
  wheel, no conda), or the eccodes `grib_to_netcdf` CLI if it is on PATH.
"""

import os
import sys
import shutil
import subprocess
import tempfile
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
# NOTE: ``cdsapi`` and ``cfgrib`` are imported lazily (inside download() /
# grib_to_netcdf()) so this module stays importable — for build_output, the
# calibration's auto-download wiring, and the tests — on machines where they are
# not installed.

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
DATASET = "cams-global-atmospheric-composition-forecasts"
# CAMS lives on the Atmosphere Data Store (ADS), NOT the Climate Data Store (CDS) that
# ~/.cdsapirc usually points at. We force the ADS endpoint here and reuse the unified
# ECMWF token (override with the ADS_API_URL / ADS_API_KEY environment variables). One-
# time setup: log in at https://ads.atmosphere.copernicus.eu and accept the site licences.
ADS_URL = "https://ads.atmosphere.copernicus.eu/api"
OUTPUT_DIR = "."
# North, West, South, East. Europe + Arctic box covering 421/427 E-PROFILE census
# stations at 0.4 deg (the north edge reaches 80 N to include Hopen 76.5 N and
# Bjornoya 74.5 N). 6 far-flung affiliates (Canada, Bonaire, New Zealand) fall
# outside this box and are not served by the European CAMS download. Override via
# the ALC_CAMS_AREA env var ("N,W,S,E") for a different domain.
AREA = [float(x) for x in os.environ.get("ALC_CAMS_AREA", "80,-30,27,45").split(",")]
RUN_TIME = ["00:00"]
LEADTIME = ["3", "6", "9", "12", "15", "18", "21", "24"]
MODEL_LEVELS = ["1"] + [str(k) for k in range(38, 138)]   # 1 (for z) + 38..137

REGRID_TO_1DEG = False                 # resample 0.4 deg -> 1 deg to match old files
KEEP_INTERMEDIATE = False             # keep the raw .grib and grib_to_netcdf .nc

VARIABLES = [
    "aerosol_extinction_coefficient_355nm",
    "aerosol_extinction_coefficient_532nm",
    "aerosol_extinction_coefficient_1064nm",
    "attenuated_backscatter_due_to_aerosol_355nm_from_ground",
    "attenuated_backscatter_due_to_aerosol_532nm_from_ground",
    "attenuated_backscatter_due_to_aerosol_1064nm_from_ground",
    "temperature",
    "specific_humidity",
    "logarithm_of_surface_pressure",
    "geopotential",
]

# The 4 model-level fields the calibration actually reads (molecular T/p, water-vapour q,
# surface z/lnsp). The aerosol optical fields above are NOT used by either calibration, so
# the auto-download omits them — that cuts the ADS request cost ~60% (the per-request limit
# scales with variables x levels x steps x dates) and keeps the daily files small.
CALIBRATION_VARIABLES = [
    "temperature",
    "specific_humidity",
    "logarithm_of_surface_pressure",
    "geopotential",
]

# Lean set for OmB + water-vapour work: aerosol backscatter at the ceilometer
# wavelengths (532 nm = Mini-MPL native, 1064 nm = CHM15k native, 910 nm derived
# via the Angstrom exponent between 532 and 1064) PLUS the molecular/WV fields.
# Skips 355 nm (no E-PROFILE instrument) and the extinction fields (OmB uses
# backscatter), so it is ~40% cheaper than the full VARIABLES set while still
# carrying everything OmB and the WV correction need.
OMB_VARIABLES = [
    "attenuated_backscatter_due_to_aerosol_532nm_from_ground",
    "attenuated_backscatter_due_to_aerosol_1064nm_from_ground",
    "temperature",
    "specific_humidity",
    "logarithm_of_surface_pressure",
    "geopotential",
]

# Physical constants (ECMWF IFS values)
G0 = 9.80665        # m s-2
RD = 287.0597       # J kg-1 K-1, dry-air gas constant
RE = 6371229.0      # m, Earth radius used by the IFS

# L137 half-level coefficients, index k = 0 (model top) .. 137 (surface).
# p_half(k) = A_HALF[k] + B_HALF[k] * surface_pressure
A_HALF = [
    0.000000, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966,
    13.605424, 18.608931, 24.985718, 32.985710, 42.879242, 54.955463,
    69.520576, 86.895882, 107.415741, 131.425507, 159.279404, 191.338562,
    227.968948, 269.539581, 316.420746, 368.982361, 427.592499, 492.616028,
    564.413452, 643.339905, 729.744141, 823.967834, 926.344910, 1037.201172,
    1156.853638, 1285.610352, 1423.770142, 1571.622925, 1729.448975, 1897.519287,
    2076.095947, 2265.431641, 2465.770508, 2677.348145, 2900.391357, 3135.119385,
    3381.743652, 3640.468262, 3911.490479, 4194.930664, 4490.817383, 4799.149414,
    5119.895020, 5452.990723, 5798.344727, 6156.074219, 6526.946777, 6911.870605,
    7311.869141, 7727.412109, 8159.354004, 8608.525391, 9076.400391, 9562.682617,
    10065.978516, 10584.631836, 11116.662109, 11660.067383, 12211.547852, 12766.873047,
    13324.668945, 13881.331055, 14432.139648, 14975.615234, 15508.256836, 16026.115234,
    16527.322266, 17008.789063, 17467.613281, 17901.621094, 18308.433594, 18685.718750,
    19031.289063, 19343.511719, 19620.042969, 19859.390625, 20059.931641, 20219.664063,
    20337.863281, 20412.308594, 20442.078125, 20425.718750, 20361.816406, 20249.511719,
    20087.085938, 19874.025391, 19608.572266, 19290.226563, 18917.460938, 18489.707031,
    18006.925781, 17471.839844, 16888.687500, 16262.046875, 15596.695313, 14898.453125,
    14173.324219, 13427.769531, 12668.257813, 11901.339844, 11133.304688, 10370.175781,
    9617.515625, 8880.453125, 8163.375000, 7470.343750, 6804.421875, 6168.531250,
    5564.382813, 4993.796875, 4457.375000, 3955.960938, 3489.234375, 3057.265625,
    2659.140625, 2294.242188, 1961.500000, 1659.476563, 1387.546875, 1143.250000,
    926.507813, 734.992188, 568.062500, 424.414063, 302.476563, 202.484375,
    122.101563, 62.781250, 22.835938, 3.757813, 0.000000, 0.000000,
]
B_HALF = [
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000700, 0.00002400, 0.00005900, 0.00011200, 0.00019900,
    0.00034000, 0.00056200, 0.00089000, 0.00135300, 0.00199200, 0.00285700,
    0.00397100, 0.00537800, 0.00713300, 0.00926100, 0.01180600, 0.01481600,
    0.01831800, 0.02235500, 0.02696400, 0.03217600, 0.03802600, 0.04454800,
    0.05177300, 0.05972800, 0.06844800, 0.07795800, 0.08828600, 0.09946200,
    0.11150500, 0.12444800, 0.13831300, 0.15312500, 0.16891000, 0.18568900,
    0.20349100, 0.22233300, 0.24224400, 0.26324200, 0.28535400, 0.30859800,
    0.33293900, 0.35825400, 0.38436300, 0.41112500, 0.43839100, 0.46600300,
    0.49380000, 0.52161900, 0.54930100, 0.57669200, 0.60364800, 0.63003600,
    0.65573600, 0.68064300, 0.70466900, 0.72773900, 0.74979700, 0.77079800,
    0.79071700, 0.80953600, 0.82725600, 0.84388100, 0.85943200, 0.87392900,
    0.88740800, 0.89990000, 0.91144800, 0.92209600, 0.93188100, 0.94086000,
    0.94906400, 0.95655000, 0.96335200, 0.96951300, 0.97507800, 0.98007200,
    0.98454200, 0.98850000, 0.99198400, 0.99500300, 0.99763000, 1.00000000,
]


# ----------------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------------
def parse_args(argv):
    """Return (date_list, label) from the command-line tokens."""
    if len(argv) == 1 and len(argv[0]) == 6:
        y, m = int(argv[0][:4]), int(argv[0][4:6])
        start = date(y, m, 1)
        end = date(y + (m == 12), (m % 12) + 1, 1) - timedelta(days=1)
        label = f"{y:04d}{m:02d}"
    elif len(argv) == 1 and len(argv[0]) == 8:
        start = end = _d(argv[0])
        label = argv[0]
    elif len(argv) == 2 and len(argv[0]) == 8 and len(argv[1]) == 8:
        start, end = _d(argv[0]), _d(argv[1])
        if end < start:
            start, end = end, start
        label = f"{argv[0]}_{argv[1]}"
    else:
        sys.exit(
            "Usage:\n"
            "  download_cams_beta.py YYYYMMDD            (single day)\n"
            "  download_cams_beta.py YYYYMMDD YYYYMMDD   (inclusive range)\n"
            "  download_cams_beta.py YYYYMM              (whole month)"
        )
    n = (end - start).days + 1
    dates = [(start + timedelta(days=i)).strftime("%Y-%m-%d") for i in range(n)]
    return dates, label


def _d(s):
    return date(int(s[:4]), int(s[4:6]), int(s[6:8]))


# ----------------------------------------------------------------------------
# Download + GRIB->netCDF
# ----------------------------------------------------------------------------
def _ads_client():
    """A cdsapi client pointed at the ADS (where CAMS lives), reusing the unified ECMWF
    token. The URL is forced to the ADS because ~/.cdsapirc usually carries the CDS URL;
    the key comes from ADS_API_KEY, else ~/.cdsapirc, else cdsapi's own resolution."""
    import cdsapi

    url = os.environ.get("ADS_API_URL") or ADS_URL
    key = os.environ.get("ADS_API_KEY")
    if not key:
        rcfile = Path.home() / ".cdsapirc"
        if rcfile.exists():
            for line in rcfile.read_text().splitlines():
                if ":" in line:
                    k, v = line.split(":", 1)
                    if k.strip() == "key":
                        key = v.strip()
                        break
    return cdsapi.Client(url=url, key=key) if key else cdsapi.Client(url=url)


def _request_too_large(exc) -> bool:
    """True if an ADS error means the request exceeded the per-request cost limit."""
    s = str(exc).lower()
    return ("cost limit" in s) or ("too large" in s) or ("reduce your selection" in s)


def _retrieve_day(client, day, variables, out_handle, grib_path):
    """Retrieve one day's model-level GRIB and append it to *out_handle*. The ADS caps the
    cost of a single request (~ variables x levels x steps); if a day is rejected as too
    large, split the variable list in half and recurse — adapting to the limit without
    hard-coding it. GRIB messages are self-describing, so appended parts form one valid
    multi-message file that cfgrib / grib_to_netcdf read as a whole."""
    part = f"{grib_path}.{day}.{len(variables)}.part"
    try:
        client.retrieve(
            DATASET,
            {
                "date": [day],
                "time": RUN_TIME,
                "leadtime_hour": LEADTIME,
                "type": ["forecast"],
                "variable": variables,
                "model_level": MODEL_LEVELS,
                "area": AREA,
                "data_format": "grib",
                "download_format": "unarchived",
            },
            part,
        )
    except Exception as exc:  # noqa: BLE001
        if _request_too_large(exc) and len(variables) > 1:
            mid = len(variables) // 2
            print(f"[download]   {day}: request too large; splitting {len(variables)} vars")
            _retrieve_day(client, day, variables[:mid], out_handle, grib_path)
            _retrieve_day(client, day, variables[mid:], out_handle, grib_path)
            return
        raise
    with open(part, "rb") as f:
        shutil.copyfileobj(f, out_handle)
    os.remove(part)


def download(dates, grib_path, variables=None):
    """Retrieve the model-level GRIB for *dates* into *grib_path*.

    Chunked by day (one ADS request per date) to stay under the ADS per-request cost
    limit, with an automatic variable-split fallback for any day that is still too large.
    *variables* defaults to the full archive set; the calibration auto-download passes the
    smaller ``CALIBRATION_VARIABLES`` (t/q/z/lnsp) it actually needs.
    """
    if variables is None:
        variables = VARIABLES
    print(f"[download] {len(dates)} day(s) x {len(variables)} var(s): {dates[0]} .. {dates[-1]}")
    c = _ads_client()
    with open(grib_path, "wb") as out:
        for d in dates:
            print(f"[download]   {d}")
            _retrieve_day(c, d, list(variables), out, grib_path)


def grib_to_netcdf(grib_path, nc_path):
    """Convert the downloaded GRIB to a netCDF shaped like the existing CAMS_Beta_*.nc
    (``time/level/latitude/longitude`` + ``t/q/z/lnsp``).

    Prefers the eccodes ``grib_to_netcdf`` CLI when it is on PATH (e.g. a conda env),
    matching the historical pipeline bit-for-bit. Otherwise falls back to a pure-Python
    conversion via :func:`_cfgrib_to_dataset` — ``cfgrib`` bundles the eccodes binary
    wheel (incl. on Windows), so no conda / no system eccodes / no CLI is required.
    """
    if shutil.which("grib_to_netcdf") is not None:
        print("[convert] grib_to_netcdf (eccodes CLI)")
        subprocess.run(["grib_to_netcdf", "-o", nc_path, grib_path], check=True)
        return
    print("[convert] cfgrib (pure-Python, no eccodes CLI)")
    ds = _cfgrib_to_dataset(grib_path)
    try:
        ds.to_netcdf(nc_path)
    finally:
        ds.close()


def _cfgrib_to_dataset(grib_path):
    """Read a CAMS model-level GRIB into a Dataset shaped like grib_to_netcdf output:
    dims ``(time, level, latitude, longitude)``, variables ``t``/``q``/``z``/``lnsp``
    (plus whatever aerosol fields cfgrib exposes), with ``z``/``lnsp`` broadcast across
    levels as surface fields (the layout :func:`build_output` and the calibration
    readers expect).

    ``cfgrib.open_datasets`` splits the heterogeneous GRIB into consistent hypercubes;
    ``time_dims=['valid_time']`` collapses the forecast (reference_time x step) grid to a
    single valid-time axis. NOTE: validated structurally against the archive layout; run
    one real download to confirm end-to-end once the ADS licences are accepted.
    """
    import cfgrib

    dss = cfgrib.open_datasets(
        grib_path, backend_kwargs={"time_dims": ["valid_time"], "indexpath": ""}
    )
    if not dss:
        raise RuntimeError(f"cfgrib found no readable GRIB messages in {grib_path}")

    norm = []
    for d in dss:
        ren = {}
        if "valid_time" in d.variables:
            ren["valid_time"] = "time"
        for lv in ("hybrid", "hybrid_level", "model_level", "modelLevel"):
            if lv in d.dims:
                ren[lv] = "level"
                break
        norm.append(d.rename(ren))

    ds = xr.merge(norm, compat="override", join="outer")
    if "level" not in ds.dims:
        raise RuntimeError(
            "cfgrib output has no model-level axis; unexpected GRIB layout "
            f"(dims={dict(ds.sizes)})"
        )
    for name in ("t", "q", "z", "lnsp"):
        if name not in ds:
            raise RuntimeError(f"cfgrib output is missing required variable '{name}'")

    # z and lnsp are surface fields; broadcast them across the model-level axis so
    # build_output/compute_geopotential read the surface value at any level slot.
    for name in ("z", "lnsp"):
        da = ds[name]
        surf = da.max("level", skipna=True) if "level" in da.dims else da
        ds[name] = surf.broadcast_like(ds["t"])

    ds = ds.sortby("level").transpose("time", "level", "latitude", "longitude")
    return ds


# ----------------------------------------------------------------------------
# Geopotential / altitude on model levels
# ----------------------------------------------------------------------------
def compute_geopotential(ds):
    """Hydrostatic integration from the surface upward over the contiguous
    model-level block. Returns z [m2 s-2], altitude ASL [m], altitude AGL [m].
    Levels that are not contiguous with the surface (e.g. the lone level 1) stay
    NaN, since the integration cannot bridge the missing 2..37."""
    levels = ds["level"].values.astype(int)
    nt, nl, ny, nx = ds["t"].shape
    A = np.asarray(A_HALF)
    B = np.asarray(B_HALF)

    ibot = int(np.where(levels == levels.max())[0][0])    # model level 137
    sp = np.exp(ds["lnsp"].isel(level=ibot).values).astype(np.float64)   # Pa
    zs = ds["z"].isel(level=ibot).values.astype(np.float64)              # m2 s-2 (orography)

    t = ds["t"].values.astype(np.float64)
    q = ds["q"].values.astype(np.float64)

    z_full = np.full((nt, nl, ny, nx), np.nan, np.float64)
    order = np.argsort(-levels)         # bottom (137) -> top
    z_h = zs.copy()                     # half-level geopotential at the surface
    prev = None
    broken = False
    for oi in order:
        lev = int(levels[oi])
        if prev is not None and lev != prev - 1:
            broken = True               # gap: cannot integrate further up
        if broken:
            prev = lev
            continue
        tv = t[:, oi] * (1.0 + 0.609133 * q[:, oi])
        if lev == 1:
            ph_lp = A[1] + B[1] * sp
            dlogp = np.log(ph_lp / 0.1)
            alpha = np.log(2.0)
        else:
            ph_l = A[lev - 1] + B[lev - 1] * sp
            ph_lp = A[lev] + B[lev] * sp
            dlogp = np.log(ph_lp / ph_l)
            alpha = 1.0 - (ph_l / (ph_lp - ph_l)) * dlogp
        z_full[:, oi] = z_h + RD * tv * alpha
        z_h = z_h + RD * tv * dlogp
        prev = lev

    height = z_full / G0                      # geopotential height
    altitude = RE * height / (RE - height)    # geometric altitude ASL
    altitude_agl = altitude - (zs / G0)[:, None, :, :]
    return z_full, altitude, altitude_agl


# ----------------------------------------------------------------------------
# Assemble + write
# ----------------------------------------------------------------------------
def build_output(ds, out_path):
    z_full, altitude, altitude_agl = compute_geopotential(ds)
    dims = ds["t"].dims
    nt, nl, ny, nx = ds["t"].shape

    keep = [
        "aerext355", "aerext532", "aerext1064",
        "aerbackscatgnd355", "aerbackscatgnd532", "aerbackscatgnd1064",
        "t", "q",
    ]
    out = ds[[v for v in keep if v in ds]].copy()

    # z and lnsp are SURFACE fields. The calibration readers parse them exactly as the
    # existing CAMS_Beta_*.nc files store them: the surface value at model level 1
    # (level index 0), NaN on every other level. _cams_levels (rayleigh/WV) reads
    # z_raw[isfinite][0]; _cams_levels_all_times (cloud) reads level index 0. Writing
    # the full model-level geopotential into z instead corrupts the surface reference
    # (the reader then picks an upper-level value ~24 km too high). The per-level
    # geopotential height is provided separately as `altitude` (geometric ASL).
    levels = ds["level"].values.astype(int)
    ibot = int(np.where(levels == levels.max())[0][0])    # surface = model level 137

    def _surface_only(name):
        full = np.full((nt, nl, ny, nx), np.nan, np.float64)
        full[:, 0] = np.asarray(ds[name].isel(level=ibot).values, np.float64)
        return full

    out["z"] = (dims, _surface_only("z"))
    out["z"].attrs = {"units": "m**2 s**-2",
                      "long_name": "Surface geopotential (orography)",
                      "standard_name": "surface_geopotential"}
    out["lnsp"] = (dims, _surface_only("lnsp"))
    out["lnsp"].attrs = {"units": "1",
                         "long_name": "Logarithm of surface pressure",
                         "standard_name": "lnsp"}

    out["altitude"] = (dims, altitude.astype(np.float32))
    out["altitude"].attrs = {"units": "m",
                             "long_name": "Geometric altitude above sea level"}
    out["altitude_agl"] = (dims, altitude_agl.astype(np.float32))
    out["altitude_agl"].attrs = {"units": "m",
                                 "long_name": "Geometric altitude above model surface"}

    if REGRID_TO_1DEG:
        lon = np.arange(AREA[1], AREA[3] + 0.001, 1.0)
        lat = np.arange(AREA[0], AREA[2] - 0.001, -1.0)
        print(f"[regrid] -> 1 deg ({lat.size} x {lon.size})")
        out = out.interp(latitude=lat, longitude=lon, method="linear")

    out.attrs["processing"] = (
        "Downloaded from CAMS ADS (forecast, 00 UTC, steps 3-24). "
        "z and altitude computed by hydrostatic integration on L137 model levels."
    )

    enc = {}
    for v in out.data_vars:
        enc[v] = {"zlib": True, "complevel": 4, "dtype": "float32"}
    enc["time"] = {"units": "hours since 1900-01-01 00:00:00", "calendar": "gregorian"}

    print(f"[write] {out_path}")
    out.to_netcdf(out_path, encoding=enc)


# ----------------------------------------------------------------------------
def download_to_netcdf(dates, out_path, variables=None, keep_intermediate=KEEP_INTERMEDIATE):
    """Download *dates* from the ADS and write the calibration-ready CAMS_Beta netCDF.

    This is the importable core of the CLI (and what the calibration's CAMS auto-download
    calls): ADS retrieve -> GRIB->netCDF (cfgrib or the eccodes CLI) -> build_output.
    *dates* is a list of "YYYY-MM-DD" strings (as produced by parse_args). *variables*
    defaults to the full archive set; pass ``CALIBRATION_VARIABLES`` for the lean t/q/z/lnsp
    file the calibration needs. Returns *out_path*.
    """
    out_path = str(out_path)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    workdir = tempfile.mkdtemp(prefix="cams_")
    grib_path = os.path.join(workdir, "cams.grib")
    raw_nc = os.path.join(workdir, "cams_raw.nc")

    download(dates, grib_path, variables=variables)
    grib_to_netcdf(grib_path, raw_nc)

    with xr.open_dataset(raw_nc) as ds:
        build_output(ds, out_path)

    if keep_intermediate:
        print(f"[info] intermediate files kept in {workdir}")
    else:
        for p in (grib_path, raw_nc):
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(workdir)
        except OSError:
            pass
    print("[done]", out_path)
    return out_path


def download_daily_files(dates, out_dir=OUTPUT_DIR, variables=None):
    """Download each day in *dates* to its OWN ``CAMS_Beta_YYYYMMDD.nc`` (one file per
    day), with the full variable set by default (aerosol optical fields + t/q/z/lnsp).

    *dates* is a list of "YYYY-MM-DD" strings (as produced by parse_args). Days whose
    file already exists are skipped, so an interrupted run resumes cleanly. Returns the
    list of written paths.
    """
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for d in dates:
        ymd = d.replace("-", "")
        out_path = os.path.join(out_dir, f"CAMS_Beta_{ymd}.nc")
        if os.path.exists(out_path):
            print(f"[skip] {out_path} already exists")
        else:
            download_to_netcdf([d], out_path, variables=variables)
        written.append(out_path)
    return written


def main():
    argv = list(sys.argv[1:])
    out_dir = OUTPUT_DIR
    if "--out" in argv:
        i = argv.index("--out")
        out_dir = argv[i + 1]
        del argv[i:i + 2]
    per_day = "--daily" in argv
    argv = [a for a in argv if a != "--daily"]

    dates, label = parse_args(argv)
    if per_day:
        download_daily_files(dates, out_dir)
    else:
        download_to_netcdf(dates, os.path.join(out_dir, f"CAMS_Beta_{label}.nc"))


if __name__ == "__main__":
    main()
