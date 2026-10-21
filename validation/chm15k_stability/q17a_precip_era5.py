# -*- coding: utf-8 -*-
"""q17a - daily ERA5 total precipitation at every CHM15k station (washout test of the aerosol hypothesis).

Source: DestinE Earth Data Hub, ERA5 single levels (ARCO Zarr, hourly, 0.25 deg), lazy nearest-grid-point
access - no domain download (each point pulls only the few compressed chunks it needs).
Needs zarr + aiohttp (run with the system Python 3.12, the repo .venv has no zarr) and the EDH token in
~/.destinationearth.

ERA5 'tp' at valid_time h is the accumulation over the previous hour, so the UTC-day total of day D is the
sum over valid_time in (D 00:00, D+1 00:00]; we shift by one hour before grouping.

Output: C:/DATA/Projects/202610_CHM15k_stability/q17/q17a_era5_tp_daily.csv  (date, stream, tp_mm, cp_mm)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

REPO = Path(__file__).resolve().parents[2]
CENSUS = REPO / "validation/scope_l1_2026_census.json"
OUT = Path("C:/DATA/Projects/202610_CHM15k_stability/q17")
EDH_URL = "https://edh:{tok}@api.earthdatahub.destine.eu/era5/reanalysis-era5-single-levels-v0.zarr"
START, END = "2024-12-01T01:00", "2026-08-01T00:00"


def chm15k_points():
    """CHM15k streams with their 0.25 deg grid point (several streams share one point)."""
    rows = []
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["type"].upper().startswith("CHM"):
            lat, lon = float(s["lat"]), float(s["lon"])
            glat = round(lat * 4) / 4.0
            glon = round((lon % 360.0) * 4) / 4.0 % 360.0
            rows.append(dict(stream=f"{s['wmo']}_{s['ident']}", site=s["site"], lat=lat, lon=lon, glat=glat, glon=glon))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tok = Path.home().joinpath(".destinationearth").read_text().strip()
    ds = xr.open_dataset(EDH_URL.format(tok=tok), chunks={}, engine="zarr")
    pts = chm15k_points()
    grid = pts[["glat", "glon"]].drop_duplicates().reset_index(drop=True)
    print(f"{len(pts)} CHM15k streams on {len(grid)} ERA5 grid points; ERA5 ends {ds.valid_time.values[-1]}", flush=True)
    out_f = OUT / "q17a_era5_tp_daily.csv"
    done = set()
    if out_f.exists():
        prev = pd.read_csv(out_f)
        done = set(zip(prev.glat.round(2), prev.glon.round(2)))
    t0 = time.time()
    for k, g in grid.iterrows():
        if (round(g.glat, 2), round(g.glon, 2)) in done:
            continue
        for attempt in range(4):
            try:
                sub = ds[["tp", "cp"]].sel(latitude=g.glat, longitude=g.glon, method="nearest")
                sub = sub.sel(valid_time=slice(START, END)).load()
                break
            except Exception as exc:  # noqa: BLE001 - network hiccups: retry
                print(f"  retry {attempt} {g.glat},{g.glon}: {exc!r}"[:200], flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            continue
        vt = pd.to_datetime(sub.valid_time.values) - pd.Timedelta(hours=1)
        df = pd.DataFrame({"tp": sub.tp.values * 1000.0, "cp": sub.cp.values * 1000.0}, index=vt)
        daily = df.groupby(df.index.floor("D")).agg(["sum", "count"])
        daily = daily[daily[("tp", "count")] == 24]
        res = pd.DataFrame({"date": daily.index.strftime("%Y-%m-%d"), "glat": g.glat, "glon": g.glon,
                            "tp_mm": daily[("tp", "sum")].round(3).values, "cp_mm": daily[("cp", "sum")].round(3).values})
        res.to_csv(out_f, mode="a", header=not out_f.exists(), index=False)
        print(f"[{k + 1}/{len(grid)}] {g.glat:.2f},{g.glon:.2f}: {len(res)} days, total {res.tp_mm.sum():.0f} mm "
              f"({time.time() - t0:.0f} s)", flush=True)
    pts.to_csv(OUT / "q17a_era5_points.csv", index=False)
    print("done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
