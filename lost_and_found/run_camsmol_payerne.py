"""run_camsmol_payerne.py

Re-run the Payerne CL61 (0-20000-0-06610_C) corrected-Rayleigh calibration with
molecular_source = 'standard' and 'cams', over the analysis period, to measure the
impact of building the molecular reference from CAMS T/p instead of US Std 1976.

Both runs use the PRODUCTION water-vapor correction (apply_wv_correction=1 from
options.json) and the CORRECT Payerne coordinates. NOTE: the L2 file and the
station manifest carry station_latitude/longitude = 0 for this CL61 (same class of
bug as the station_altitude=0 we fixed in RAW2L2); left as-is the WV/CAMS read
would sample the Gulf of Guinea. We override the coordinates here so both runs are
physically correct and differ ONLY in the molecular source.

Writes, per source:
  C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_<tag>/0-20000-0-06610_C/0-20000-0-06610_C_cl.csv
(same schema as run_all_l2monthly.py) so load_rayleigh_python_kalman.m can read it.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import calendar as cal
import csv
import logging
import warnings
from pathlib import Path

import numpy as np

from rayleigh_calibration import (
    calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel,
)
from rayleigh_calibration.config import InstrumentType

BASE = Path("A:/E-PROFILE_L2_monthly")
WMO, IDENT = "0-20000-0-06610", "C"
KEY = f"{WMO}_{IDENT}"
# Correct Payerne coordinates (file/manifest carry 0,0; see module docstring).
LAT, LON, ALT = 46.8137, 6.9425, 491.0
FIELDS = ["date", "flag", "lidar_constant", "uncertainty",
          "bottom_height", "top_height", "message"]


def dates_2026(months=(2, 3, 4, 5, 6)):
    out = []
    for mo in months:
        for day in range(1, cal.monthrange(2026, mo)[1] + 1):
            out.append(f"2026{mo:02d}{day:02d}")
    return out


def run(source):
    tag = "camsmol" if source == "cams" else "stdmol_check"
    outdir = Path(f"C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_{tag}/{KEY}")
    outdir.mkdir(parents=True, exist_ok=True)

    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = BASE
    o.data_level = DataLevel.L2_MONTHLY
    o.molecular_source = source
    o.plot_all = o.plot_main = False
    o.folder_output = outdir

    info = InstrumentInfo(
        site_name=WMO, wmo_id=WMO, identifier=IDENT,
        instrument_type=InstrumentType.CL61,
        latitude=LAT, longitude=LON, altitude=ALT,
    )

    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)

    rows = []
    for ds in dates_2026():
        try:
            r = calibrate_rayleigh(ds, info, o)
            rows.append(dict(date=ds, flag=r.flag, lidar_constant=r.lidar_constant,
                             uncertainty=r.uncertainty,
                             bottom_height=r.calibration_bottom_height,
                             top_height=r.calibration_top_height, message=r.message))
        except Exception as exc:
            rows.append(dict(date=ds, flag=-99, lidar_constant=-1, uncertainty=0,
                             bottom_height=None, top_height=None,
                             message=f"{type(exc).__name__}: {exc}"))
    with open(outdir / f"{KEY}_cl.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    return rows


def ok_by_date(rows):
    return {r["date"]: float(r["lidar_constant"])
            for r in rows if str(r["flag"]) in ("1", "1.0", "0.5")}


print("Payerne CL61 Rayleigh re-run: standard vs CAMS molecular "
      "(WV on, correct coords)", flush=True)
res = {}
for src in ("standard", "cams"):
    rows = run(src)
    cl = list(ok_by_date(rows).values())
    res[src] = ok_by_date(rows)
    if cl:
        print(f"  {src:9s}: {len(cl):2d} ok nights, median CL = {np.median(cl):.5f}",
              flush=True)
    else:
        print(f"  {src:9s}: no ok nights", flush=True)

common = sorted(set(res["standard"]) & set(res["cams"]))
if common:
    ratios = np.array([res["cams"][k] / res["standard"][k] for k in common])
    med_std = np.median(list(res["standard"].values()))
    med_cam = np.median(list(res["cams"].values()))
    print(f"\nPaired nights N={len(common)}: CL(CAMS)/CL(std) "
          f"mean={ratios.mean():.4f}  median={np.median(ratios):.4f}  "
          f"std={ratios.std():.4f}")
    print(f"median-CL ratio (CAMS/std) = {med_cam / med_std:.4f}  "
          f"({(med_cam / med_std - 1) * 100:+.2f} %)")
    print("Production fullcal_all median CL = 0.641 (correct-coords baseline).")
print("CAMSMOL_RUN_DONE", flush=True)
