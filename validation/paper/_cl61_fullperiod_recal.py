# -*- coding: utf-8 -*-
"""Full-period Payerne CL61 Rayleigh recalibration, WITH and WITHOUT the measured hood-offset
correction, over the whole CL61 record (2026-02-24 .. 2026-06-30). The offset (physical AC-coupling
model b_phys, cl61_b_dark.npz) is subtracted from rcs_0 IN MEMORY by patching load_l1_data (avoids
copying 24 GB of L1). Subtraction happens before the WV division (correct order: the offset is an
additive term in the raw rcs_0). Writes an incremental CSV so progress is visible."""
import sys, warnings, io, contextlib
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import calibration.rayleigh.calibration as RC
from calibration.config import CalibrationOptions, InstrumentInfo, InstrumentType, DataLevel
from validation.paper.calib_benchmark import RAYLEIGH_METHOD, L1_ROOT, CAMS, REPO

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "C", 46.813, 6.943, 491.0
CL_CLOUD = 1.4252
OUT_CSV = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_fullperiod_recal.csv")

# offset correction (physical AC-coupling model), Mm-1 sr-1 -> SI for rcs_0
dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
rng_d = dk["rng"].astype(float); b_SI = dk["b_phys"].astype(float) * 1e-6

STATE = {"correct": False}
# calibrate_rayleigh loads via load_data(file_list, instrument_type, data_level) (line ~312),
# so we must wrap that reference in the rayleigh module namespace, not load_l1_data.
_orig_load_data = RC.load_data
def _patched_load_data(file_list, instrument_type, data_level):
    d = _orig_load_data(file_list, instrument_type, data_level)
    if d is not None and STATE["correct"]:
        corr = np.interp(d.range_alc, rng_d, b_SI)
        d.rcs = d.rcs - corr[None, :]
    return d
RC.load_data = _patched_load_data

def opts():
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = L1_ROOT; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = True
    import tempfile
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    return o

INFO = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)

days = sorted(p.name[len(f"L1_{WMO}_{IDENT}"):len(f"L1_{WMO}_{IDENT}") + 8]
              for p in (L1_ROOT / WMO).glob(f"2026/*/L1_{WMO}_{IDENT}2026*.nc"))
print(f"CL61 record: {len(days)} days  {days[0]} .. {days[-1]}")

def cal(ds, correct):
    STATE["correct"] = correct
    o = opts()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            r = RC.calibrate_rayleigh(ds, INFO, o)
    except Exception:
        return np.nan, -99
    cl = float(r.lidar_constant) if (np.isfinite(r.lidar_constant) and r.lidar_constant > 0) else np.nan
    return cl, int(r.flag) if np.isfinite(r.flag) else -99

with open(OUT_CSV, "w") as f:
    f.write("date,CL_orig,flag_orig,CL_corr,flag_corr\n")
    for k, ds in enumerate(days):
        c0, f0 = cal(ds, False)
        c1, f1 = cal(ds, True)
        f.write(f"{ds},{c0:.6g},{f0},{c1:.6g},{f1}\n"); f.flush()
        if (k + 1) % 10 == 0 or k == len(days) - 1:
            print(f"  {k+1}/{len(days)} done ({ds})", flush=True)

# summary
import pandas as pd
df = pd.read_csv(OUT_CSV)
ok = df[(df.flag_orig.isin([1, 0.5])) & (df.flag_corr.isin([1, 0.5])) & df.CL_orig.notna() & df.CL_corr.notna()]
o3 = ok[(ok.CL_orig < 3 * ok.CL_orig.median()) & (ok.CL_corr < 3 * ok.CL_corr.median())]
print(f"\ncalibrated days: orig {df.CL_orig.notna().sum()}, corr {df.CL_corr.notna().sum()}, "
      f"both {len(ok)} ({len(o3)} after outlier drop)")
print(f"MEDIAN C_L   orig {o3.CL_orig.median():.3f}   corr {o3.CL_corr.median():.3f}   "
      f"shift {100*(o3.CL_corr.median()/o3.CL_orig.median()-1):+.1f}%")
print(f"gap to cloud {CL_CLOUD}:  orig {100*(o3.CL_orig.median()/CL_CLOUD-1):+.1f}%  ->  "
      f"corr {100*(o3.CL_corr.median()/CL_CLOUD-1):+.1f}%")
print(f"night-to-night scatter (robust):  orig {100*1.4826*(o3.CL_orig-o3.CL_orig.median()).abs().median()/o3.CL_orig.median():.1f}%"
      f"   corr {100*1.4826*(o3.CL_corr-o3.CL_corr.median()).abs().median()/o3.CL_corr.median():.1f}%")
print("FULLPERIOD_RECAL_DONE saved", OUT_CSV)
