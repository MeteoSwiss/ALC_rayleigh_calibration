"""run_rayleigh_diag_dates.py

Regenerate the Rayleigh-calibration diagnostics for SPECIFIC Payerne CL61 nights
(default: the 2026-03-11/12/13 window that exposed the low-R² selection), so the
new molecular-window gates can be compared night-by-night against the old figure.

Same configuration as run_rayleigh_diag_payerne.py (reads options.json -> inherits
the new min_window_start_m / min_window_r2 / max_window_rel_error gates).
"""
from __future__ import annotations
import os
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import logging
import warnings
from pathlib import Path

from calibration import (
    calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel,
)
from calibration.config import InstrumentType

BASE = Path("A:/E-PROFILE_L2_monthly")
WMO, IDENT = "0-20000-0-06610", "C"
LAT, LON, ALT = 46.8137, 6.9425, 491.0
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")

DATES = sys.argv[1:] or ["20260311", "20260312", "20260313"]

o = CalibrationOptions.from_json(Path("options.json"))
o.folder_root = BASE
o.data_level = DataLevel.L2_MONTHLY
o.molecular_source = "standard"
o.plot_main = True
o.plot_all = True
o.folder_output = OUT

print(f"window gates: start>={o.min_window_start_m:.0f} m  R2>={o.min_window_r2:.2f}  "
      f"rel_err<={o.max_window_rel_error:.0f}%")

info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61,
                      latitude=LAT, longitude=LON, altitude=ALT)

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

for ds in DATES:
    try:
        r = calibrate_rayleigh(ds, info, o)
    except Exception as exc:
        print(f"  {ds}: EXCEPTION {exc}")
        continue
    bot = r.calibration_bottom_height
    top = r.calibration_top_height
    win = f"{bot:.0f}-{top:.0f} m" if bot is not None and top is not None else "n/a"
    print(f"  {ds}: flag={r.flag} CL={r.lidar_constant:.4g}  window {win}  {r.message}")

print(f"\nplots -> {OUT / 'plots' / WMO / '2026'}")
print("RAYLEIGH_DIAG_DATES_DONE")
