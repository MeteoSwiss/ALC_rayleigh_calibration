# -*- coding: utf-8 -*-
"""q16e check 1: do masked (_FillValue = -999.9) rcs_0 cells reach the Rayleigh fit?

calibration.io.data_loader.load_l1_data concatenates the per-file masked arrays with np.concatenate,
which DROPS the mask; average_ceilometer_data then takes np.asarray(...) of the raw data, so a masked
cell enters the mean as -999.9. This script counts how often that happens in CHM15k L1 files
(night hours, 1.5-8 km range) for a geographically spread sample of units and months.
"""
from __future__ import annotations

import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import L1_ROOT, OUT, chm15k_streams  # noqa: E402

Q16 = OUT / "q16"
DAYS = ["20250115", "20250415", "20250615", "20250915", "20251215", "20260315", "20260615"]


def scan(args):
    import netCDF4
    wmo, ident, day = args
    f = L1_ROOT / wmo / day[:4] / day[4:6] / f"L1_{wmo}_{ident}{day}.nc"
    if not f.exists():
        return None
    try:
        with netCDF4.Dataset(f) as d:
            rng = np.asarray(d["range"][:], float)
            t = np.asarray(d["time"][:], float)
            r = d["rcs_0"][:]
            m = np.ma.getmaskarray(r)
            raw = np.asarray(r.data, float)
            lon = float(d["station_longitude"][:])
    except Exception as exc:  # noqa: BLE001
        return dict(stream=f"{wmo}_{ident}", day=day, err=str(exc)[:60])
    sh = np.mod((t % 1) * 24 + lon / 15, 24)
    night = (sh >= 22) | (sh < 2)
    band = (rng >= 1500) & (rng <= 8000)
    sub = m[np.ix_(night, band)]
    return dict(stream=f"{wmo}_{ident}", day=day, n_prof=int(night.sum()),
                frac_masked_all=float(m.mean()), frac_masked_night_band=float(sub.mean()) if sub.size else np.nan,
                prof_with_mask_night=float(np.any(sub, axis=1).mean()) if sub.size else np.nan,
                fill_value_seen=float(np.unique(raw[m])[0]) if m.any() else np.nan,
                n_nan_unmasked=int(np.isnan(raw[~m]).sum()))


def main():
    streams = chm15k_streams()
    rng = np.random.default_rng(16)
    pick = [streams[i] for i in rng.choice(len(streams), size=min(30, len(streams)), replace=False)]
    jobs = [(w, i, d) for (w, i, *_rest) in pick for d in DAYS]
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = [r for r in ex.map(scan, jobs) if r]
    keys = sorted({k for r in rows for k in r})
    with open(Q16 / "q16e_masked_fill.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    ok = [r for r in rows if "err" not in r]
    fm = np.array([r["frac_masked_night_band"] for r in ok], float)
    pm = np.array([r["prof_with_mask_night"] for r in ok], float)
    print(f"files {len(ok)}; files with any masked night/1.5-8 km cell: {(fm > 0).sum()}; "
          f"max masked fraction {np.nanmax(fm):.2e}; max share of night profiles touched {np.nanmax(pm):.3f}")
    for r in ok:
        if r["frac_masked_night_band"] > 0:
            print(r)


if __name__ == "__main__":
    main()
