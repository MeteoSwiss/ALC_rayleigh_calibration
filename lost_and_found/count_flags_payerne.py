"""Quick flag tally for Payerne CL61 over a month range, NO plots.
Reports how many nights calibrate vs. are rejected under the new window gates,
and (for comparison) the best-window R2/start of each night.
"""
from __future__ import annotations
import os, sys, calendar as cal, logging, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
from pathlib import Path
from collections import Counter
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType

WMO, IDENT = "0-20000-0-06610", "C"
LAT, LON, ALT = 46.8137, 6.9425, 491.0
YEAR, MONTH = 2026, int(sys.argv[1]) if len(sys.argv) > 1 else 3

o = CalibrationOptions.from_json(Path("options.json"))
o.folder_root = Path("A:/E-PROFILE_L2_monthly")
o.data_level = DataLevel.L2_MONTHLY
o.molecular_source = "standard"
o.plot_main = False
o.plot_all = False
o.folder_output = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")
info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

flags = Counter()
ok = []
for day in range(1, cal.monthrange(YEAR, MONTH)[1] + 1):
    ds = f"{YEAR}{MONTH:02d}{day:02d}"
    try:
        r = calibrate_rayleigh(ds, info, o)
    except Exception:
        flags["exception"] += 1
        continue
    flags[r.flag] += 1
    if r.flag in (1, 1.0, 0.5):
        ok.append((ds, r.lidar_constant, r.calibration_bottom_height, r.calibration_top_height))

print(f"=== Payerne CL61 {YEAR}-{MONTH:02d}  gates: start>={o.min_window_start_m:.0f} R2>={o.min_window_r2:.2f} rel<={o.max_window_rel_error:.0f}% ===")
for k in sorted(flags, key=lambda x: str(x)):
    print(f"  flag {k}: {flags[k]}")
n_ok = len(ok)
print(f"  -> calibrated nights: {n_ok}")
for ds, cl, b, t in ok:
    print(f"     {ds}  CL={cl:.4g}  window {b:.0f}-{t:.0f} m")
print("COUNT_DONE")
