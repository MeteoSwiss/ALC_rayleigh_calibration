"""Sensitivity of the Payerne CL61 calibrated-night count to the molecular-window
gates. Re-runs March 2026 for several min_window_start_m / min_window_r2 settings
(no plots, writable output) so we can judge whether the 2000 m / R2>=0.5 defaults
are too strict. ONE read per (setting, day)."""
from __future__ import annotations
import os, sys, calendar as cal, logging, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
from pathlib import Path
from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from rayleigh_calibration.config import InstrumentType

WMO, IDENT = "0-20000-0-06610", "C"
LAT, LON, ALT = 46.8137, 6.9425, 491.0
YEAR, MONTH = 2026, 3
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")

info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

SETTINGS = [
    (1000, 0.5), (1500, 0.5), (2000, 0.5), (2500, 0.5),
    (1500, 0.4), (2000, 0.4),
]
ndays = cal.monthrange(YEAR, MONTH)[1]
print(f"=== Payerne CL61 {YEAR}-{MONTH:02d}: calibrated nights vs gate (rel<=50%, {ndays} days) ===")
for start_m, r2 in SETTINGS:
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = Path("A:/E-PROFILE_L2_monthly")
    o.data_level = DataLevel.L2_MONTHLY
    o.molecular_source = "standard"
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    o.min_window_start_m = float(start_m)
    o.min_window_r2 = float(r2)
    n_ok = 0
    for day in range(1, ndays + 1):
        ds = f"{YEAR}{MONTH:02d}{day:02d}"
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception:
            continue
        if r.flag in (1, 1.0, 0.5):
            n_ok += 1
    print(f"  start>={start_m:5d} m  R2>={r2:.2f} :  {n_ok:2d} / {ndays} nights")
print("SWEEP_DONE")
