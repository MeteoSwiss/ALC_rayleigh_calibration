# -*- coding: utf-8 -*-
"""q17a - CAMS aerosol tracers for the advection test (dust / smoke) and the free-tropospheric climatology.

Downloads from the ADS (cams-global-atmospheric-composition-forecasts, 00 UTC run, lead 0 h = the night
centre of European calibration nights), one request per month and product, over the CHM15k box
(N61 W-10 S36 E28, native 0.4 deg):
  * 'sl': single-level AOD at 550 nm by species (dust, organic matter, black carbon, sulphate, sea salt,
    total) + total AOD at 1240 nm;
  * 'ml': aerosol extinction at 1064 nm on 9 model levels spanning ~1.8-7.8 km ASL (113..89, every 3rd).
Then extracts the nearest grid point of every CHM15k station into one compact npz.

GRIB parts are kept under q17/q17a_cams_grib/ (total ~0.4 GB). Re-running skips months already present.
Never kill this while a request is running server-side: it keeps running and the queue slot is lost.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from calibration.io.download_cams_beta import _ads_client  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
CENSUS = REPO / "validation/scope_l1_2026_census.json"
OUT = Path("C:/DATA/Projects/202610_CHM15k_stability/q17")
GRIB = OUT / "q17a_cams_grib"
DATASET = "cams-global-atmospheric-composition-forecasts"
AREA = [61, -10, 36, 28]
SL_VARS = ["dust_aerosol_optical_depth_550nm", "organic_matter_aerosol_optical_depth_550nm",
           "black_carbon_aerosol_optical_depth_550nm", "sulphate_aerosol_optical_depth_550nm",
           "sea_salt_aerosol_optical_depth_550nm", "total_aerosol_optical_depth_550nm",
           "total_aerosol_optical_depth_1240nm"]
ML_LEVELS = [str(k) for k in range(89, 114, 3)]          # 89 (~7.8 km) ... 113 (~1.8 km ASL)
# nominal geometric altitude ASL of those levels (from a CAMS profile at 49 N, ps ~ 1000 hPa)
ML_ALT_KM = {89: 7.77, 92: 6.89, 95: 6.00, 98: 5.13, 101: 4.30, 104: 3.54, 107: 2.87, 110: 2.30, 113: 1.81}
MONTHS = [f"{y}-{m:02d}" for y in (2025, 2026) for m in range(1, 13) if (y, m) <= (2026, 7)]


def month_dates(ym):
    import calendar
    y, m = int(ym[:4]), int(ym[5:])
    n = calendar.monthrange(y, m)[1]
    return f"{ym}-01/{ym}-{n:02d}"


def download(products, months=None):
    GRIB.mkdir(parents=True, exist_ok=True)
    c = _ads_client()
    for prod in products:
        for ym in (months or MONTHS):
            f = GRIB / f"cams_{prod}_{ym}.grib"
            if (f.exists() and f.stat().st_size > 0) or Path(str(f) + ".part").exists():
                continue
            req = {"date": [month_dates(ym)], "time": ["00:00"], "leadtime_hour": ["0"], "type": ["forecast"],
                   "area": AREA, "data_format": "grib", "download_format": "unarchived"}
            if prod == "sl":
                req["variable"] = SL_VARS
            else:
                req["variable"] = ["aerosol_extinction_coefficient_1064nm"]
                req["model_level"] = ML_LEVELS
            t0 = time.time()
            print(f"request {prod} {ym}", flush=True)
            c.retrieve(DATASET, req, str(f) + ".part")
            Path(str(f) + ".part").rename(f)
            print(f"  got {f.name}: {f.stat().st_size / 1e6:.1f} MB in {time.time() - t0:.0f} s", flush=True)


def stations():
    rows = []
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["type"].upper().startswith("CHM"):
            rows.append((f"{s['wmo']}_{s['ident']}", float(s["lat"]), float(s["lon"])))
    return rows


def extract():
    """Nearest-grid-point series per CHM15k station -> q17a_cams_points.npz."""
    import xarray as xr
    st = stations()
    keys = np.array([s[0] for s in st])
    lat = np.array([s[1] for s in st])
    lon = np.array([s[2] for s in st])
    inside = (lat <= AREA[0]) & (lat >= AREA[2]) & (lon >= AREA[1]) & (lon <= AREA[3])
    out = {}
    for prod in ("sl", "ml"):
        times, data = [], []
        for ym in MONTHS:
            f = GRIB / f"cams_{prod}_{ym}.grib"
            if not f.exists():
                continue
            dss = __import__("cfgrib").open_datasets(str(f), backend_kwargs={"indexpath": ""})
            ds = xr.merge(dss, compat="override")
            ilat = xr.DataArray(np.where(inside, lat, np.nan), dims="station")
            ilon = xr.DataArray(np.where(inside, lon, np.nan), dims="station")
            sub = ds.sel(latitude=ilat.fillna(45.0), longitude=ilon.fillna(5.0), method="nearest")
            tt = (sub["time"].values + (sub["step"].values if "step" in sub.coords else np.timedelta64(0, "h")))
            tt = (np.atleast_1d(tt).astype("datetime64[s]").astype("f8") / 86400.0)
            names = sorted(v for v in sub.data_vars)
            arr = np.stack([np.asarray(sub[v].transpose("time", "station", *[d for d in sub[v].dims
                                                                         if d not in ("time", "station")]).values)
                            for v in names], axis=-1)
            times.append(tt)
            data.append(arr)
            out[f"{prod}_names"] = np.array(names)
            if prod == "ml" and "hybrid" in sub.coords:
                out["ml_levels"] = np.asarray(sub["hybrid"].values)
            print(f"extract {prod} {ym}: {arr.shape} vars {names}", flush=True)
        if times:
            out[f"{prod}_t"] = np.concatenate(times)
            d = np.concatenate(data, axis=0).astype("f4")
            d[:, ~inside] = np.nan
            out[f"{prod}_data"] = d
    out["keys"], out["lat"], out["lon"] = keys, lat, lon
    np.savez_compressed(OUT / "q17a_cams_points.npz", **out)
    print("wrote q17a_cams_points.npz", {k: v.shape for k, v in out.items()}, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("download", "extract"))
    ap.add_argument("--products", default="sl,ml")
    ap.add_argument("--months", default="", help="comma list YYYY-MM (parallel workers); default all")
    a = ap.parse_args()
    if a.cmd == "download":
        download(a.products.split(","), [m for m in a.months.split(",") if m] or None)
    else:
        extract()


if __name__ == "__main__":
    sys.exit(main())
