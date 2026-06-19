"""run_rayleigh_diag_payerne.py

Save the per-night Rayleigh-calibration DIAGNOSTICS for the Payerne CL61
(0-20000-0-06610_C): the molecular-fit diagnostic (RCS vs molecular reference + fit
window), the window search, the lidar-constant panel and the RCS time series. These
let us inspect whether the molecular fit window is aerosol-contaminated and whether
the fit is clean.

Runs with the correct Payerne coordinates and the production WV correction; plots are
written to figs_paper_validation/rayleigh_diag/plots/<wmo>/<year>/.
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
import calendar as cal
import logging
import warnings
from pathlib import Path

from rayleigh_calibration import (
    calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel,
)
from rayleigh_calibration.config import InstrumentType

BASE = Path("A:/E-PROFILE_L2_monthly")
WMO, IDENT = "0-20000-0-06610", "C"
LAT, LON, ALT = 46.8137, 6.9425, 491.0
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")

o = CalibrationOptions.from_json(Path("options.json"))
o.folder_root = BASE
o.data_level = DataLevel.L2_MONTHLY
o.molecular_source = "standard"
o.plot_main = True
o.plot_all = True                      # adds the window-search panel
o.folder_output = OUT

info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61,
                      latitude=LAT, longitude=LON, altitude=ALT)

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

# Try March nights; save diagnostics for the first few clear (flag 1 / 0.5) nights.
saved = 0
for day in range(1, cal.monthrange(2026, 3)[1] + 1):
    ds = f"202603{day:02d}"
    try:
        r = calibrate_rayleigh(ds, info, o)
    except Exception as exc:
        continue
    if r.flag in (1, 1.0, 0.5):
        print(f"  {ds}: flag={r.flag} CL={r.lidar_constant:.4g} "
              f"window {r.calibration_bottom_height:.0f}-{r.calibration_top_height:.0f} m  {r.message}")
        saved += 1
        if saved >= 3:
            break

print(f"\nsaved diagnostics for {saved} clear night(s) to: {OUT / 'plots' / WMO / '2026'}")
print("RAYLEIGH_DIAG_DONE")
