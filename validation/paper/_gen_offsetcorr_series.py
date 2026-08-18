# -*- coding: utf-8 -*-
"""Generate the offset-corrected Rayleigh calibration series (daily -> Kalman) for the Payerne
CHM15k (NEW) and CL61 (refresh), by recalibrating the full record with the measured hood offset
subtracted from rcs_0 in memory (patch load_data), then reusing the operational Kalman bridge and
the calib CSV format. Writes <key>_offsetcorr_{L1,L2,}.csv into the CALIB dir so
run_paper_validation can plot them as extra channels. The offset is applied in rcs_0 units:
CL61 b_phys is Mm^-1 sr^-1 (x1e-6 -> SI); CHM15k b_phys is already in rcs_0 (counts/s.m^2)."""
import sys, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import calibration.rayleigh.calibration as RC
from validation.paper.calib_benchmark import raw_rayleigh, kalman, OUT, BENCHMARK

D = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/"
PAY = BENCHMARK["payerne"]; START, END = PAY["start"], PAY["end"]
CHM = next(c for c in PAY["channels"] if c["itype"] == "CHM15k" and c["calib"] == "rayleigh")
CL61 = next(c for c in PAY["channels"] if c["itype"] == "CL61" and c["calib"] == "rayleigh")
print(f"period {START}..{END}")

# offset correction applied in memory (rcs_0 units) via load_data wrapper
STATE = {"rng": None, "corr": None}
_orig = RC.load_data
def _patched(file_list, instrument_type, data_level):
    d = _orig(file_list, instrument_type, data_level)
    if d is not None and STATE["corr"] is not None:
        d.rcs = d.rcs - np.interp(d.range_alc, STATE["rng"], STATE["corr"])[None, :]
    return d
RC.load_data = _patched

def write_series(key, dates, C, Cstd):
    res = kalman(dates, C, Cstd, normalise=True)
    if res is None:
        print(f"  {key}: no Kalman series"); return
    for suffix in ("_L1", "_L2", ""):
        with open(OUT / f"{key}_offsetcorr{suffix}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
            w.writerows(res)
    nk = sum(1 for r in res if np.isfinite(r[3]))
    ck = [r[3] for r in res if np.isfinite(r[3])]
    print(f"  {key}_offsetcorr: {len(dates)} nights -> {len(res)} rows ({nk} Kalman), "
          f"C_kalman {min(ck):.3f}..{max(ck):.3f}")

# ---- CHM15k (NEW): b_phys already in rcs_0 units ----
chm = np.load(D + "chm15k_b_dark.npz")
STATE["rng"], STATE["corr"] = chm["rng"].astype(float), chm["b_phys"].astype(float)
dc, Cc, Sc = raw_rayleigh(CHM, START, END, level="L1")
print(f"CHM15k corrected: {len(Cc)} nights, median C={np.median(Cc):.3e}" if len(Cc) else "CHM15k: 0 nights")
STATE["corr"] = None
dn, Cn, Sn = raw_rayleigh(CHM, START, END, level="L1")     # native, for the shift report
if len(Cn) and len(Cc):
    print(f"  CHM15k native median C={np.median(Cn):.3e}  -> offset shift {100*(np.median(Cc)/np.median(Cn)-1):+.1f}%")
STATE["rng"], STATE["corr"] = chm["rng"].astype(float), chm["b_phys"].astype(float)
write_series("0-20000-0-06610_A_rayleigh", dc, Cc, Sc)

# ---- CL61 (refresh with the current physical model): b_phys Mm^-1 sr^-1 -> x1e-6 ----
cl = np.load(D + "cl61_b_dark.npz")
STATE["rng"], STATE["corr"] = cl["rng"].astype(float), cl["b_phys"].astype(float) * 1e-6
dc2, Cc2, Sc2 = raw_rayleigh(CL61, START, END, level="L1")
print(f"CL61 corrected: {len(Cc2)} nights, median C={np.median(Cc2):.3f}" if len(Cc2) else "CL61: 0 nights")
write_series("0-20000-0-06610_C_rayleigh", dc2, Cc2, Sc2)
print("GEN_OFFSETCORR_DONE")
