# -*- coding: utf-8 -*-
"""q16g (task 1a, side table): firmware (Lufft software_version) and calibration-software version per
station-year of the 11-year archive, read from the MeteoSwiss MATLAB output A:/CHM15k/<STN>/CeiloCalib_*.nc.
Lets q16g_archive_years.csv be split by firmware generation. Output: q16g_archive_firmware.csv
"""
from __future__ import annotations

import csv
import glob
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

Q16 = Path("C:/DATA/Projects/202610_CHM15k_stability/q16")


def main():
    rows = []
    for f in sorted(glob.glob("A:/CHM15k/*/CeiloCalib_*.nc")):
        stn = Path(f).parent.name.lower()
        with netCDF4.Dataset(f) as d:
            t = np.ma.filled(d["start_time"][:].astype("f8"), np.nan)
            sm = d["software_version_measurements"][:]
            sw = netCDF4.chartostring(sm) if sm.dtype.kind == "S" else np.asarray(sm)
            sw = np.array([str(x).replace("\x00", "").strip() for x in np.atleast_1d(sw)])
            if sw.size != t.size:          # stored (soft_mes, time) in some files
                sw = netCDF4.chartostring(np.asarray(sm).T)
                sw = np.array([str(x).replace("\x00", "").strip() for x in sw])
            vc = np.ma.filled(d["software_version_calibration"][:].astype("f8"), np.nan)
        years = np.array([(datetime(1904, 1, 1) + timedelta(seconds=float(x))).year if np.isfinite(x) else -1
                          for x in t])
        for y in np.unique(years[years > 0]):
            s = years == y
            fw, n = np.unique(sw[s], return_counts=True)
            rows.append(dict(station=stn, file=Path(f).name, year=int(y), n=int(s.sum()),
                             firmware="; ".join(f"{a or 'n/a'} ({b})" for a, b in zip(fw, n)),
                             calib_software=";".join(f"{v:.3f}" for v in np.unique(vc[s][np.isfinite(vc[s])]))))
    with open(Q16 / "q16g_archive_firmware.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
