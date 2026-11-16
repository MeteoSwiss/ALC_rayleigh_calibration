# -*- coding: utf-8 -*-
"""Seasonal behaviour of the raw Lufft housekeeping nn1 / nn2 / nn3 / scaling at night (report 19, q17).

Wiegner & Geiss (2012, AMT 5, 1953) document that the CHM15k(x) firmware changes the APD high voltage (nn1,
0.1 V units) automatically with the background light and resets it at 01:00 UTC; each 0.5 V step changed the
lidar constant by ~24 %. This script reads the raw daily Lufft NetCDF archive (A:/CHM15k/<site>), every 3rd
day, and keeps night medians (00-03 UTC) of nn1, nn2, nn3, scaling, p_calc and temp_int.
Output: C:/DATA/Projects/202610_CHM15k_stability/q17/q17_nn1_<site>.csv
"""
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

OUT = Path("C:/DATA/Projects/202610_CHM15k_stability/q17")
VARS = ("nn1", "nn2", "nn3", "scaling", "p_calc", "temp_int", "temp_lom", "temp_det")


def night(f):
    import netCDF4
    try:
        with netCDF4.Dataset(f) as d:
            t = netCDF4.num2date(d.variables["time"][:], d.variables["time"].units)
            hr = np.array([x.hour for x in t])
            m = hr < 3
            if m.sum() < 10:
                return None
            row = {"file": Path(f).name, "date": f"{t[0]:%Y-%m-%d}",
                   "module": getattr(d, "serlom", getattr(d, "serial_lom", "")),
                   "firmware": getattr(d, "software_version", "")}
            for v in VARS:
                if v in d.variables:
                    x = np.asarray(d.variables[v][:], float)
                    row[v] = float(np.nanmedian(x[m])) if x.ndim == 1 else np.nan
                    if v == "nn1":
                        row["nn1_day_min"] = float(np.nanmin(x)) if x.ndim == 1 else np.nan
                        row["nn1_changes"] = int(np.sum(np.diff(x) != 0)) if x.ndim == 1 else -1
                else:
                    row[v] = np.nan
            return row
    except Exception:
        return None


def main(site):
    files = sorted(Path(f"A:/CHM15k/{site}").rglob("*.nc"))[::3]
    with ProcessPoolExecutor(8) as ex:
        rows = [r for r in ex.map(night, files, chunksize=8) if r]
    OUT.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with open(OUT / f"q17_nn1_{site}.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(site, len(rows))


if __name__ == "__main__":
    for s in sys.argv[1:]:
        main(s)
