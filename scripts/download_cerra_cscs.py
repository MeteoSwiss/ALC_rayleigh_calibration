#!/usr/bin/env python3
"""Download CERRA pressure-level profiles at the E-PROFILE ALC stations (for the WV correction).

CERRA is the Copernicus regional reanalysis for Europe at **5.5 km** (Lambert conformal,
~1069x1069 grid) — far finer than CAMS (0.4 deg ~= 44 km), so its T/q profiles resolve
Alpine and coastal orography much better where the CAMS water-vapour correction is weakest.

Why pressure levels (not model levels): CERRA ``reanalysis-cerra-model-levels`` only carries
t/q/u/v — there is **no geopotential** there, so altitude would have to be reconstructed from
the 106-level HARMONIE hybrid coefficients + surface pressure + orography. The pressure-level
product natively carries **geopotential(z) + temperature(t) + relative_humidity(r)**, giving
altitude, T, and (via r,T,p) specific humidity q in one clean download. 1000..1 hPa (29 levels).

Why full-domain download + local subset: CERRA lives ONLY on the public CDS
(``reanalysis-cerra-pressure-levels``). Polytope feature-extraction does NOT serve CERRA
(class=rr is rejected by every feature/MARS datasource on the ECMWF endpoint), and the CDS
cannot subset the projected Lambert grid by lat/lon. So we fetch the whole grid per day, then
extract the nearest grid column at each census station and DISCARD the full grid — ~540 GB of
transfer collapses to ~5 GB of kept station profiles.

Output: one netCDF per month ``CERRA_stations_<YYYYMM>.nc`` with dims (time, level, station)
and variables z/t/r/q plus geopotential height, and per-station coords lat/lon/alt/site/wmo/type.

Resumable: a month whose ``CERRA_stations_<YYYYMM>.nc`` already exists is skipped; within a
month, per-day station-profile parts are cached so an interrupted month resumes day-by-day.

WHERE TO RUN (CSCS): the CDS download needs outbound HTTPS, which balfrin COMPUTE nodes lack.
Run on a login node inside tmux (network-bound, not compute-bound), e.g.:

    tmux new-session -d -s cerra 'bash run_cerra.sh; exec bash'

where run_cerra.sh calls this script with the base-env python (which has cdsapi/cfgrib):

    /users/mhrvo/miniforge3/bin/python -u scripts/download_cerra_cscs.py \
        --start 202501 --end 202603 \
        --census validation/scope_l1_2026_census.json \
        --out /scratch/mch/mhrvo/CERRA_stations

The CDS queue is rate-limited; the full 15-month run can take days. The per-month/per-day
resume makes that safe across interruptions.
"""
from __future__ import annotations

import argparse
import calendar
import json
import os
import sys
import traceback
from datetime import date, timedelta
from pathlib import Path

import numpy as np

DATASET = "reanalysis-cerra-pressure-levels"
CDS_URL = "https://cds.climate.copernicus.eu/api"

# All 29 CERRA pressure levels (hPa), surface -> top. The water-vapour correction only needs
# the troposphere (>= ~200 hPa); pass --levels tropo to drop the 10 stratospheric levels and
# cut the transfer by ~1/3. Default keeps the full set for generality.
LEVELS_ALL = ["1000", "975", "950", "925", "900", "875", "850", "825", "800", "750",
              "700", "600", "500", "400", "300", "250", "200", "150", "100", "70",
              "50", "30", "20", "10", "7", "5", "3", "2", "1"]
LEVELS_TROPO = ["1000", "975", "950", "925", "900", "875", "850", "825", "800", "750",
                "700", "600", "500", "400", "300", "250", "200"]

# CERRA HRES analyses are 3-hourly (8 per day).
TIMES = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]

VARIABLES = ["geopotential", "temperature", "relative_humidity"]
G0 = 9.80665  # m s-2, for geopotential -> geopotential height


def month_list(start: str, end: str) -> list[str]:
    """Inclusive list of 'YYYYMM' from start to end."""
    y, m = int(start[:4]), int(start[4:6])
    ey, em = int(end[:4]), int(end[4:6])
    out = []
    while (y, m) <= (ey, em):
        out.append(f"{y:04d}{m:02d}")
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return out


def month_days(ym: str, last_allowed: date) -> list[date]:
    """Every day of the month as date objects, capped at last_allowed."""
    y, m = int(ym[:4]), int(ym[4:6])
    n = calendar.monthrange(y, m)[1]
    out = []
    for d in range(1, n + 1):
        day = date(y, m, d)
        if day <= last_allowed:
            out.append(day)
    return out


def load_stations(census_path: Path) -> dict:
    """Read the census JSON into parallel arrays of station lat/lon/metadata.

    The census can list the same site under several idents; we keep unique (wmo, ident)
    rows so each physical stream maps to one column (they may share a grid cell, which is
    fine — the nearest-neighbour lookup is cheap)."""
    rows = json.loads(census_path.read_text())
    seen = set()
    lats, lons, alts, sites, wmos, idents, types = [], [], [], [], [], [], []
    for r in rows:
        key = (r.get("wmo"), r.get("ident"))
        if key in seen:
            continue
        seen.add(key)
        lats.append(float(r["lat"]))
        lons.append(float(r["lon"]))
        alts.append(float(r.get("alt", np.nan)))
        sites.append(str(r.get("site", "")))
        wmos.append(str(r.get("wmo", "")))
        idents.append(str(r.get("ident", "")))
        types.append(str(r.get("type", "")))
    return {
        "lat": np.array(lats), "lon": np.array(lons), "alt": np.array(alts),
        "site": sites, "wmo": wmos, "ident": idents, "type": types,
        "n": len(lats),
    }


def _cds_client():
    import cdsapi
    return cdsapi.Client(url=os.environ.get("CDS_API_URL", CDS_URL),
                         quiet=False, wait_until_complete=True)


def download_day_grib(client, day: date, levels, grib_path: Path) -> None:
    """Retrieve one day's full-domain CERRA pressure-level GRIB (all 8 analysis times)."""
    client.retrieve(
        DATASET,
        {
            "variable": VARIABLES,
            "pressure_level": levels,
            "data_type": ["reanalysis"],
            "product_type": ["analysis"],
            "year": [f"{day.year:04d}"],
            "month": [f"{day.month:02d}"],
            "day": [f"{day.day:02d}"],
            "time": TIMES,
            "data_format": "grib",
        },
        str(grib_path),
    )


def _nearest_indices(grid_lat, grid_lon, st_lat, st_lon):
    """Nearest flat grid index for each station via a KD-tree on the unit sphere.

    Works directly on the CERRA 2D lat/lon arrays (no projection maths needed); great-circle
    nearest neighbour is exact enough at 5.5 km for a point column. The tree is built once per
    run because the CERRA grid is static."""
    from scipy.spatial import cKDTree

    def to_xyz(lat, lon):
        la, lo = np.deg2rad(lat), np.deg2rad(lon)
        return np.stack([np.cos(la) * np.cos(lo),
                         np.cos(la) * np.sin(lo),
                         np.sin(la)], axis=-1)

    gxyz = to_xyz(np.asarray(grid_lat).ravel(), np.asarray(grid_lon).ravel())
    sxyz = to_xyz(np.asarray(st_lat), np.asarray(st_lon))
    tree = cKDTree(gxyz)
    _, idx = tree.query(sxyz, k=1)
    return idx


def subset_day(grib_path: Path, stations: dict, levels, nn_cache: dict):
    """Extract per-station column profiles from a day's GRIB -> Dataset(time, level, station).

    Reads t/r/z (short names), finds the nearest grid column for every station (cached across
    days), derives geopotential height and specific humidity, and returns an xarray Dataset."""
    import xarray as xr

    # cfgrib returns t/r/z on isobaricInhPa with 2D latitude/longitude coords. Read the three
    # short names together; fall back to per-variable reads if cfgrib refuses the mix.
    def _open(**kw):
        return xr.open_dataset(grib_path, engine="cfgrib",
                               backend_kwargs={"indexpath": "", **kw})

    try:
        ds = _open()
    except Exception:
        parts = []
        for sn in ("t", "r", "z"):
            parts.append(_open(filter_by_keys={"shortName": sn}))
        ds = xr.merge(parts)

    lat2d = ds["latitude"].values
    lon2d = ds["longitude"].values
    # CERRA longitudes are 0..360; census lon is -180..180 -> wrap for the KD-tree.
    lon2d = np.where(lon2d > 180.0, lon2d - 360.0, lon2d)

    if "idx" not in nn_cache:
        nn_cache["idx"] = _nearest_indices(lat2d, lon2d, stations["lat"], stations["lon"])
        nn_cache["shape"] = lat2d.shape
        nn_cache["glat"] = lat2d.ravel()[nn_cache["idx"]]
        nn_cache["glon"] = lon2d.ravel()[nn_cache["idx"]]
    idx = nn_cache["idx"]

    lev = ds["isobaricInhPa"].values.astype("float64")
    time = ds["time"].values if ds["time"].ndim else np.array([ds["time"].values])
    nt = time.size

    out = {}
    for sn in ("t", "r", "z"):
        if sn not in ds:
            continue
        arr = ds[sn].values  # (time, level, y, x) or (level, y, x) if single time
        if arr.ndim == 3:
            arr = arr[None, ...]
        nt_, nl_, ny, nx = arr.shape
        flat = arr.reshape(nt_, nl_, ny * nx)
        out[sn] = flat[:, :, idx]  # (time, level, station)

    coords = {
        "time": time,
        "level": lev,
        "station": np.arange(stations["n"]),
    }
    data_vars = {}
    if "t" in out:
        data_vars["t"] = (("time", "level", "station"), out["t"])
    if "r" in out:
        data_vars["r"] = (("time", "level", "station"), out["r"])
    if "z" in out:
        data_vars["z"] = (("time", "level", "station"), out["z"])
        data_vars["gph"] = (("time", "level", "station"), out["z"] / G0)  # geopotential height [m]
    # specific humidity from RH, T, p (Magnus over water); p from the pressure-level coordinate.
    if "t" in out and "r" in out:
        T = out["t"]
        RH = out["r"]
        p = (lev[None, :, None] * 100.0)  # hPa -> Pa, broadcast (1, level, 1)
        esat = 611.2 * np.exp(17.62 * (T - 273.15) / (T - 30.03))  # Pa
        e = np.clip(RH, 0.0, None) / 100.0 * esat
        q = 0.622 * e / (p - 0.378 * e)
        data_vars["q"] = (("time", "level", "station"), q.astype("float32"))

    result = xr.Dataset(data_vars, coords=coords)
    result = result.assign_coords(
        st_lat=("station", stations["lat"]),
        st_lon=("station", stations["lon"]),
        st_alt=("station", stations["alt"]),
        grid_lat=("station", nn_cache["glat"]),
        grid_lon=("station", nn_cache["glon"]),
        site=("station", np.array(stations["site"], dtype=object)),
        wmo=("station", np.array(stations["wmo"], dtype=object)),
        ident=("station", np.array(stations["ident"], dtype=object)),
        type=("station", np.array(stations["type"], dtype=object)),
    )
    ds.close()
    return result


def process_month(client, ym: str, days, stations, levels, out_dir: Path, tmp_dir: Path):
    """Download+subset every day of a month, concatenate, write CERRA_stations_<YYYYMM>.nc."""
    import xarray as xr

    out_path = out_dir / f"CERRA_stations_{ym}.nc"
    if out_path.exists():
        print(f"[skip] {out_path.name} (exists)")
        return "skip"

    nn_cache: dict = {}
    day_parts = []
    for day in days:
        ds_part = tmp_dir / f"CERRA_stations_{day:%Y%m%d}.nc"
        if ds_part.exists():
            print(f"[cache] {ds_part.name}")
            day_parts.append(ds_part)
            continue
        grib_path = tmp_dir / f"cerra_{day:%Y%m%d}.grib"
        print(f"[get ] {day:%Y-%m-%d}  (8 times x {len(levels)} lvl x {len(VARIABLES)} var, full domain)",
              flush=True)
        try:
            download_day_grib(client, day, levels, grib_path)
            prof = subset_day(grib_path, stations, levels, nn_cache)
            prof.to_netcdf(ds_part)
            prof.close()
            day_parts.append(ds_part)
            print(f"[sub ] {ds_part.name}  ({stations['n']} stations)", flush=True)
        except Exception as exc:  # noqa: BLE001 — keep going; the run is resumable
            print(f"[FAIL] {day:%Y-%m-%d}: {exc}")
            traceback.print_exc()
        finally:
            if grib_path.exists():
                grib_path.unlink()  # discard the full-domain grid; keep only station columns

    if not day_parts:
        print(f"[FAIL] {ym}: no days retrieved")
        return "fail"
    month_ds = xr.open_mfdataset([str(p) for p in day_parts], combine="by_coords")
    month_ds.to_netcdf(out_path)
    month_ds.close()
    for p in day_parts:
        p.unlink()
    print(f"[ok  ] {out_path.name}  ({len(day_parts)} days)", flush=True)
    return "ok"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="202501", help="first month YYYYMM")
    ap.add_argument("--end", default="202603", help="last month YYYYMM")
    ap.add_argument("--census", default=os.environ.get("ALC_CENSUS",
                    "validation/scope_l1_2026_census.json"),
                    help="census JSON with station lat/lon")
    ap.add_argument("--out", default=os.environ.get("ALC_CERRA_DIR", "."),
                    help="output folder (default $ALC_CERRA_DIR or .)")
    ap.add_argument("--levels", choices=["all", "tropo"], default="all",
                    help="'all' 29 levels (1000..1 hPa) or 'tropo' 17 levels (1000..200 hPa)")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_parts"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    stations = load_stations(Path(args.census))
    levels = LEVELS_ALL if args.levels == "all" else LEVELS_TROPO
    yesterday = date.today() - timedelta(days=1)

    print(f"[cerra] dataset={DATASET}  {len(levels)} levels ({args.levels})  "
          f"{stations['n']} stations -> {out_dir}", flush=True)
    client = _cds_client()
    months = month_list(args.start, args.end)
    n_ok = n_skip = n_fail = 0
    for ym in months:
        days = month_days(ym, yesterday)
        if not days:
            print(f"[skip] {ym}: no days on/before {yesterday}")
            continue
        status = process_month(client, ym, days, stations, levels, out_dir, tmp_dir)
        n_ok += status == "ok"
        n_skip += status == "skip"
        n_fail += status == "fail"
    print(f"[cerra] done: {n_ok} months written, {n_skip} skipped, {n_fail} failed "
          f"of {len(months)}", flush=True)


if __name__ == "__main__":
    main()
