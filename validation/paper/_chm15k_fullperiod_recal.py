# -*- coding: utf-8 -*-
"""Full-period Payerne CHM15k Rayleigh recalibration WITH and WITHOUT the measured hood-offset
correction (counterpart of _cl61_fullperiod_recal). The CHM15k offset b_phys is already in rcs_0
units (counts/s.m^2), subtracted in memory via the load_data patch. Photon-counting CHM15k has no
cloud-method reference, so the figure compares the corrected constant to the NATIVE median: the
correction is expected to be small (within the ~13 % night scatter) -> confirms it is not material."""
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

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "A", 46.813, 6.943, 491.0
OUT_CSV = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/chm15k_fullperiod_recal.csv")

dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/chm15k_b_dark.npz")
rng_d = dk["rng"].astype(float); b_off = dk["b_phys"].astype(float)   # already rcs_0 units

STATE = {"correct": False}
_orig = RC.load_data
def _patched(file_list, instrument_type, data_level):
    d = _orig(file_list, instrument_type, data_level)
    if d is not None and STATE["correct"]:
        d.rcs = d.rcs - np.interp(d.range_alc, rng_d, b_off)[None, :]
    return d
RC.load_data = _patched

def opts():
    import tempfile
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = L1_ROOT; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = False                    # CHM15k 1064 nm -> no WV
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    return o

INFO = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CHM15k, latitude=LAT, longitude=LON, altitude=ALT)
days = sorted(p.name[len(f"L1_{WMO}_{IDENT}"):len(f"L1_{WMO}_{IDENT}") + 8]
              for p in (L1_ROOT / WMO).glob(f"2026/*/L1_{WMO}_{IDENT}2026*.nc"))
print(f"CHM15k record: {len(days)} days {days[0]}..{days[-1]}")

def cal(ds, correct):
    STATE["correct"] = correct
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            r = RC.calibrate_rayleigh(ds, INFO, opts())
    except Exception:
        return np.nan, -99
    cl = float(r.lidar_constant) if (np.isfinite(r.lidar_constant) and r.lidar_constant > 0) else np.nan
    return cl, int(r.flag) if np.isfinite(r.flag) else -99

with open(OUT_CSV, "w") as f:
    f.write("date,CL_orig,flag_orig,CL_corr,flag_corr\n")
    for k, ds in enumerate(days):
        c0, f0 = cal(ds, False); c1, f1 = cal(ds, True)
        f.write(f"{ds},{c0:.6g},{f0},{c1:.6g},{f1}\n"); f.flush()
        if (k + 1) % 20 == 0 or k == len(days) - 1:
            print(f"  {k+1}/{len(days)} done", flush=True)

import pandas as pd
df = pd.read_csv(OUT_CSV)
for fl, tag in ([1], "flag=1"), ([1, 0.5, 0], "flag>=0"):
    o = df[df.flag_orig.isin(fl) & df.CL_orig.notna()]; c = df[df.flag_corr.isin(fl) & df.CL_corr.notna()]
    if len(o) and len(c):
        print(f"{tag}: orig n={len(o)} med={o.CL_orig.median():.3e} | corr n={len(c)} med={c.CL_corr.median():.3e} "
              f"| shift {100*(c.CL_corr.median()/o.CL_orig.median()-1):+.1f}%")
print("CHM_FULLPERIOD_DONE")
