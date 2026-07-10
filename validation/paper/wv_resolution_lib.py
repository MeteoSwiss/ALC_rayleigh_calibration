"""Shared library for the CAMS-resolution / water-vapour (WV) sensitivity study.

The operational WV correction of the 910 nm cloud calibration reads humidity from CAMS at a
coarse resolution: 0.4 deg horizontal (~44 km, smoothed orography), 3-hourly, ~101 model levels.
This module builds the WV number-density profile n_wv(z) from three independent sources and turns
each into the two-way transmission T2_wv(z) and its impact on the cloud calibration:

  - CAMS 0.4 deg 3-hourly  (operational baseline; D:/CAMS_Monthly_04)
  - ECMWF/IFS at the point  (Cloudnet ..._ecmwf.nc; native ~9 km, HOURLY, 137 levels)
  - Payerne radiosonde      (in-situ vertical truth, 00/12 UT)
  - ERA5 via earthkit       (optional cross-check; 0.25 deg, hourly, 37 pressure levels)

The cloud calibration integrates B = int(beta'/T2_wv) dz up to cloud base, and C = 1/(2 B S_liq).
The integrand is dominated by the cloud-base peak, so to first order C is proportional to
T2_wv(cbh): the relative calibration error from a WV source is (T2_source(cbh)/T2_truth(cbh) - 1).
We report that, and also run the real pipeline (Payerne) to confirm the magnitude and sign.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from calibration.water_vapor_correction.water_vapor import (  # noqa: E402
    _cams_levels, two_way_wv_transmission, DEFAULT_ABS_CROSS_SECTION,
    LASER_SPECTRUM, cams_nearest_offset_deg, KB, EPS, G0)
from calibration.cloud.calibration import _nw_from_T_RH  # noqa: E402
from calibration.water_vapor_correction.water_vapor import cams_levels_all_times  # noqa: E402

# ---------------------------------------------------------------- paths / data
CAMS_04 = Path("D:/CAMS_Monthly_04")          # operational 0.4 deg, 3-hourly, monthly files
CAMS_1D = Path("D:/CAMS")                      # legacy 1 deg archive (L137 for recent years)
IFS_DIR = Path("C:/Users/hervo/Downloads/cloudnet-collection-00529733605c43b6")  # Payerne IFS
SND_DIR = Path("D:/Soundings")                 # Payerne radiosonde sounding_pay_YYYY.csv
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/wv_resolution")
REPO_FIGS = Path(__file__).resolve().parents[2] / "doc" / "reports" / "figs_wv_resolution"

# 10 stations: 5 flat + 5 mountain/complex terrain. key = <wmo>_<ident> (dashboard style).
STATIONS = [
    # --- flat / plateau / coastal ---
    dict(key="0-20000-0-06610_C", name="Payerne",   ctry="CH", inst="CL61", lat=46.813, lon=6.943, alt=490.0,  terrain="flat"),
    dict(key="0-20000-0-06447_B", name="Uccle",      ctry="BE", inst="CL61", lat=50.798, lon=4.359, alt=100.0,  terrain="flat"),
    dict(key="0-276-13-19987_A",  name="Koln",       ctry="DE", inst="CL51", lat=50.928, lon=6.928, alt=50.0,   terrain="flat"),
    dict(key="0-20000-0-07643_A", name="Montpellier",ctry="FR", inst="CL31", lat=43.577, lon=3.963, alt=1.0,    terrain="flat"),
    dict(key="0-20000-0-02615_A", name="Falsterbo",  ctry="SE", inst="CL31", lat=55.383, lon=12.816,alt=2.0,    terrain="flat"),
    # --- mountain / complex terrain ---
    dict(key="0-380-5-1_B",       name="Aosta",      ctry="IT", inst="CL61", lat=45.742, lon=7.357, alt=560.0,  terrain="mountain"),
    dict(key="0-20000-0-11208_A", name="Radstadt",   ctry="AT", inst="CL51", lat=47.383, lon=13.439,alt=835.0,  terrain="mountain"),
    dict(key="0-20000-0-06735_A", name="Adelboden",  ctry="CH", inst="CL31", lat=46.492, lon=7.561, alt=1327.0, terrain="mountain"),
    dict(key="0-20000-0-06792_A", name="Samedan",    ctry="CH", inst="CL31", lat=46.526, lon=9.879, alt=1709.0, terrain="mountain"),
    dict(key="0-20000-0-06751_A", name="Robiei",     ctry="CH", inst="CL31", lat=46.442, lon=8.515, alt=1896.0, terrain="mountain"),
]

PAYERNE = STATIONS[0]


def laser(inst: str):
    return LASER_SPECTRUM.get(inst, (910.0, 3.4))


# ---------------------------------------------------------------- CAMS (0.4 deg, 3-hourly)
def cams_month_file(day8: str) -> Path:
    return CAMS_04 / f"CAMS_Beta_{day8[:6]}.nc"


def cams1_month_file(day8: str) -> Path:
    """Legacy 1 deg CAMS monthly file (filename case varies: CAMS_Beta_ / CAMS_beta_)."""
    for name in (f"CAMS_Beta_{day8[:6]}.nc", f"CAMS_beta_{day8[:6]}.nc"):
        p = CAMS_1D / name
        if p.is_file():
            return p
    return CAMS_1D / f"CAMS_Beta_{day8[:6]}.nc"


def nwv_cams(cams_file: Path, lat: float, lon: float, when: np.datetime64, win_min: int = 90):
    """CAMS n_wv(z ASL) near time `when` (mean of the steps within +/-win_min; nearest if none)."""
    t0 = when - np.timedelta64(win_min, "m")
    t1 = when + np.timedelta64(win_min, "m")
    lev = _cams_levels(cams_file, lat, lon, t0, t1)
    if lev is None:
        return None
    H, _T, _P, n = lev
    return np.asarray(H, float), np.asarray(n, float)


_CAMS_DT_CACHE: dict = {}


def _cams_dt(cams_file: Path, lat: float, lon: float):
    """Cached (datetime64 steps, z[lev,t] ASL, nw[lev,t]) at the nearest CAMS grid point."""
    k = (str(cams_file), round(lat, 4), round(lon, 4))
    if k not in _CAMS_DT_CACHE:
        tnum, z, _T, nw = cams_levels_all_times(str(cams_file), lat, lon)
        dt = np.datetime64("1970-01-01") + ((tnum - 719529) * 86400).round().astype("timedelta64[s]")
        _CAMS_DT_CACHE[k] = (dt, np.asarray(z, float), np.asarray(nw, float))
    return _CAMS_DT_CACHE[k]


def nwv_cams_fast(cams_file: Path, lat: float, lon: float, when: np.datetime64):
    """CAMS n_wv(z ASL) at the nearest 3-hourly step to `when` (cached whole-file read)."""
    dt, z, nw = _cams_dt(cams_file, lat, lon)
    ti = int(np.abs(dt - when).argmin())
    H, n = z[:, ti], nw[:, ti]
    o = np.argsort(H)
    return H[o], n[o]


def cams_surface_alt(cams_file: Path, lat: float, lon: float) -> float:
    """Model-orography altitude [m] of the CAMS grid cell nearest the station (surface geopotential)."""
    import xarray as xr
    with xr.open_dataset(cams_file) as ds:
        sub = ds.sel(latitude=lat, longitude=lon, method="nearest")
        z = np.asarray(sub["z"].values, float)
        z = z[np.isfinite(z)]
    return float(z.flat[0]) / G0 if z.size else np.nan


def cams_grid(cams_file: Path, lat: float, lon: float):
    """(offset_deg to nearest grid point, grid spacing_deg)."""
    return cams_nearest_offset_deg(cams_file, lat, lon)


_ORO_CACHE: dict = {}


def cams_orography_field(cams_file: Path):
    """(lats, lons, surface_altitude[lat,lon]) of the CAMS 0.4 deg grid (model orography)."""
    import xarray as xr
    k = str(cams_file)
    if k in _ORO_CACHE:
        return _ORO_CACHE[k]
    with xr.open_dataset(cams_file) as ds:
        lats = np.asarray(ds["latitude"].values, float)
        lons = np.asarray(ds["longitude"].values, float)
        alt = np.asarray(ds["altitude"].isel(time=0, level=-1).values, float)   # level 137 = surface
    _ORO_CACHE[k] = (lats, lons, alt)
    return lats, lons, alt


def altitude_matched_cell(cams_file: Path, station: dict, radius_deg: float = 1.0):
    """Nearby CAMS cell (within radius) whose model orography best matches the true station
    altitude -> (lat, lon, cell_alt). A sensitivity probe for the orography error: for a valley
    station the nearest cell is a smoothed high average, and a lower cell nearby is a better proxy
    for the sub-cloud air the beam actually traverses."""
    lats, lons, alt = cams_orography_field(cams_file)
    la, lo = station["lat"], station["lon"]
    ila = np.where(np.abs(lats - la) <= radius_deg)[0]
    ilo = np.where(np.abs(((lons - lo + 180) % 360) - 180) <= radius_deg)[0]
    best, bd = None, 1e9
    for i in ila:
        for j in ilo:
            a = alt[i, j]
            if np.isfinite(a) and abs(a - station["alt"]) < bd:
                bd, best = abs(a - station["alt"]), (float(lats[i]), float(lons[j]), float(a))
    return best


# ---------------------------------------------------------------- ERA5 (0.25 deg, independent)
def nwv_era5(era5_file: Path, lat: float, lon: float, when: np.datetime64):
    """ERA5 n_wv(z ASL) at the point/time nearest (lat,lon,when) from pressure-level q/t/z."""
    import xarray as xr
    with xr.open_dataset(era5_file) as ds:
        sub = ds.sel(latitude=lat, longitude=lon, method="nearest")
        tname = "valid_time" if "valid_time" in sub.coords or "valid_time" in sub.dims else "time"
        t = np.asarray(sub[tname].values)
        ti = int(np.abs(t - when).argmin())
        sub = sub.isel({tname: ti})
        P = np.asarray(sub["pressure_level"].values, float) * 100.0    # hPa -> Pa
        q = np.asarray(sub["q"].values, float)
        T = np.asarray(sub["t"].values, float)
        z = np.asarray(sub["z"].values, float)                        # geopotential m2/s2
    h = z / G0
    Pw = q * P / (EPS + (1.0 - EPS) * q)
    n = Pw / (KB * T)
    good = np.isfinite(h) & np.isfinite(n)
    h, n = h[good], n[good]
    o = np.argsort(h)
    return h[o], n[o]


# ---------------------------------------------------------------- IFS (Cloudnet, hourly, at point)
def ifs_file_payerne(day8: str) -> Path:
    return IFS_DIR / f"{day8}_payerne_ecmwf.nc"


def _ifs_surface_alt(ds) -> float:
    return float(np.asarray(ds["sfc_geopotential"].values).flat[0]) / G0


def nwv_ifs(ifs_file: Path, when: np.datetime64):
    """IFS n_wv(z ASL) at the hour nearest `when` from Cloudnet q/T/pressure (height is AGL)."""
    import xarray as xr
    with xr.open_dataset(ifs_file) as ds:
        t = np.asarray(ds["time"].values)
        ti = int(np.abs(t - when).argmin())
        q = np.asarray(ds["q"].values[ti], float)              # kg/kg
        T = np.asarray(ds["temperature"].values[ti], float)    # K
        P = np.asarray(ds["pressure"].values[ti], float)       # Pa
        h_agl = np.asarray(ds["height"].values[ti], float)     # m AGL
        sfc = _ifs_surface_alt(ds)
    Pw = q * P / (EPS + (1.0 - EPS) * q)
    n = Pw / (KB * T)
    h_asl = h_agl + sfc
    order = np.argsort(h_asl)
    return h_asl[order], n[order]


def nwv_ifs_all_hours(ifs_file: Path):
    """All hourly IFS n_wv(z ASL) profiles for a day: (times, list[(h_asl,n_wv)])."""
    import xarray as xr
    out_t, out_p = [], []
    with xr.open_dataset(ifs_file) as ds:
        t = np.asarray(ds["time"].values)
        sfc = _ifs_surface_alt(ds)
        for ti in range(t.size):
            q = np.asarray(ds["q"].values[ti], float)
            T = np.asarray(ds["temperature"].values[ti], float)
            P = np.asarray(ds["pressure"].values[ti], float)
            h = np.asarray(ds["height"].values[ti], float) + sfc
            Pw = q * P / (EPS + (1.0 - EPS) * q)
            n = Pw / (KB * T)
            o = np.argsort(h)
            out_t.append(t[ti]); out_p.append((h[o], n[o]))
    return np.asarray(out_t), out_p


# ---------------------------------------------------------------- Radiosonde (Payerne)
_SND_CACHE: dict = {}


def _load_sounding_year(year: int):
    import pandas as pd
    if year in _SND_CACHE:
        return _SND_CACHE[year]
    fp = SND_DIR / f"sounding_pay_{year}.csv"
    if not fp.is_file():
        _SND_CACHE[year] = None
        return None
    s = pd.read_csv(fp)
    if "RH" in s.columns and "rh" not in s.columns:
        s = s.rename(columns={"RH": "rh"})
    s["dt"] = pd.to_datetime(s["t"] - 719529, unit="D")   # MATLAB datenum -> UTC
    _SND_CACHE[year] = s
    return s


def nwv_sounding(day8: str, hour: int = 0, alt: float = 490.0):
    """Payerne radiosonde n_wv(z ASL) for the launch nearest `day8`T`hour`Z (T,rh -> n_wv)."""
    import pandas as pd
    year = int(day8[:4])
    s = _load_sounding_year(year)
    if s is None:
        return None
    target = pd.Timestamp(f"{day8[:4]}-{day8[4:6]}-{day8[6:8]}T{hour:02d}:00:00")
    win = s[(s["dt"] >= target - pd.Timedelta(hours=1.5)) & (s["dt"] <= target + pd.Timedelta(hours=3.0))]
    win = win.dropna(subset=["z", "T", "rh"])
    win = win[(win["z"] >= alt - 25) & (win["z"] <= 12000)].sort_values("z")
    if len(win) < 20:
        return None
    # keep a single launch (the first within the window)
    t0 = win["dt"].min()
    win = win[win["dt"] <= t0 + pd.Timedelta(hours=2.5)]
    z = win["z"].to_numpy(float)
    n = _nw_from_T_RH(win["T"].to_numpy(float), win["rh"].to_numpy(float))
    # de-duplicate identical altitudes (interp needs strictly increasing)
    z, uidx = np.unique(z, return_index=True)
    return z, n[uidx]


# ---------------------------------------------------------------- T2_wv and impact
def zgrid_asl(alt: float, top_agl: float = 6000.0, dz: float = 15.0):
    return alt + np.arange(0.0, top_agl + dz, dz)


def t2_of(h_asl, n_wv, alt: float, inst: str, zg=None):
    """Two-way WV transmission T2_wv on an ASL grid for the station instrument's laser spectrum."""
    if zg is None:
        zg = zgrid_asl(alt)
    lam, fwhm = laser(inst)
    t2 = two_way_wv_transmission(zg, alt, np.asarray(h_asl, float), np.asarray(n_wv, float),
                                 DEFAULT_ABS_CROSS_SECTION, lam, fwhm)
    return zg, np.asarray(t2, float)


def t2_at(zg, t2, alt: float, cbh_agl: float) -> float:
    """T2_wv at a cloud base height (AGL) by interpolation on the ASL grid."""
    return float(np.interp(alt + cbh_agl, zg, t2))


def iwv_mm(h_asl, n_wv) -> float:
    """Integrated water-vapour column [kg/m2 = mm] from a number-density profile (rho = n*m_H2O)."""
    m_h2o = 2.9915e-26  # kg per H2O molecule
    h = np.asarray(h_asl, float); n = np.asarray(n_wv, float)
    o = np.argsort(h); h, n = h[o], n[o]
    good = np.isfinite(h) & np.isfinite(n)
    return float(np.trapz(n[good] * m_h2o, h[good]))


def dC_pct(t2_source_cbh: float, t2_truth_cbh: float) -> float:
    """Relative cloud-calibration error (%) from a WV source vs truth: C ~ T2_wv(cbh)."""
    if not (np.isfinite(t2_source_cbh) and np.isfinite(t2_truth_cbh) and t2_truth_cbh > 0):
        return np.nan
    return 100.0 * (t2_source_cbh / t2_truth_cbh - 1.0)


if __name__ == "__main__":
    # smoke test: all three sources for one winter day at Payerne
    day = "20250115"
    p = PAYERNE
    when = np.datetime64(f"2025-01-15T00:00:00")
    c = nwv_cams(cams_month_file(day), p["lat"], p["lon"], when)
    i = nwv_ifs(ifs_file_payerne(day), when)
    s = nwv_sounding(day, 0, p["alt"])
    print("CAMS surface alt:", round(cams_surface_alt(cams_month_file(day), p["lat"], p["lon"]), 0),
          "  station alt:", p["alt"], "  grid offset/spacing:", cams_grid(cams_month_file(day), p["lat"], p["lon"]))
    for lab, prof in (("CAMS", c), ("IFS", i), ("SONDE", s)):
        if prof is None:
            print(f"  {lab}: None"); continue
        zg, t2 = t2_of(prof[0], prof[1], p["alt"], p["inst"])
        print(f"  {lab:6s}: IWV={iwv_mm(*prof):5.1f} mm  T2@1km={t2_at(zg,t2,p['alt'],1000):.3f}  "
              f"T2@2km={t2_at(zg,t2,p['alt'],2000):.3f}  T2@3km={t2_at(zg,t2,p['alt'],3000):.3f}")
