"""CL61 offset-correction experiment (Payerne).

Procedure (equivalent to: un-range-correct, remove the measured dark offset, re-range-correct —
the two are identical because range correction is the fixed z^2 map: (P - P_dark) z^2 =
beta - beta_dark):
  1. b_dark(z) = median beta_att profile of the three covered-telescope windows (May/June 2026),
     smoothed 300 m (measured directly, cl61_b_dark.npz);
  2. corrected copies of the nightly L1 files with rcs_0' = rcs_0 - b_dark(z);
  3. calibrate_rayleigh (eprof_v2, WV on) on original vs corrected archives;
  4. compare the per-night lidar constants with the cloud-method constant C_L = 1.425
     (the CHM15k-anchored reference).
"""
import sys, shutil, tempfile, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

from calibration.config import CalibrationOptions, InstrumentInfo, InstrumentType, DataLevel
from calibration.rayleigh.calibration import calibrate_rayleigh
from validation.paper.calib_benchmark import RAYLEIGH_METHOD, L1_ROOT, CAMS, REPO

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "C", 46.813, 6.943, 491.0
NIGHTS = ["20260225", "20260316", "20260328", "20260402", "20260407", "20260410", "20260411",
          "20260417", "20260420", "20260422", "20260423", "20260424", "20260428", "20260521", "20260602"]
CL_CLOUD = 1.4252

dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
rng_d, b_dark = dk["rng"], dk["b_dark"]
# smooth 300 m (boxcar) and zero the first 300 m (near-range dark is contaminated by the hood)
b_s = dk["b_smooth"]   # robust per-gate median (330 m running median), see _cl61_dark_windows
print(f"b_dark (smoothed): 3 km {b_s[np.argmin(np.abs(rng_d-3000))]:+.4f}  "
      f"5 km {b_s[np.argmin(np.abs(rng_d-5000))]:+.4f}  8 km {b_s[np.argmin(np.abs(rng_d-8000))]:+.4f} Mm-1 sr-1")

# corrected archive (mirror layout)
CORR = Path(tempfile.mkdtemp(prefix="cl61_corr_"))
for ds in NIGHTS:
    src = L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
    if not src.exists():
        continue
    dst = CORR / WMO / ds[:4] / ds[4:6] / src.name
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    with Dataset(dst, "r+") as nc:
        r = np.asarray(nc.variables["range"][:], "f8")
        corr = np.interp(r, rng_d, b_s) * 1e-6            # Mm^-1 sr^-1 -> SI
        v = nc.variables["rcs_0"]
        x = np.asarray(v[:], "f8")
        if x.shape[-1] == r.size:
            v[:] = x - corr[None, :]
        else:
            v[:] = (x.T - corr[None, :]).T

def run(root, tag):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = root; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = True
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                          instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
    out = {}
    for ds in NIGHTS:
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception:
            continue
        if r.flag in (1, 0.5) and np.isfinite(r.lidar_constant) and r.lidar_constant > 0:
            out[ds] = float(r.lidar_constant)
    return out

orig = run(L1_ROOT, "original")
corr = run(CORR, "corrected")
print("CORRALL", {k: round(v, 4) for k, v in sorted(corr.items())})
common = sorted(set(orig) & set(corr))
co = np.array([orig[d] for d in common]); cc = np.array([corr[d] for d in common])
ok = (co < 3 * np.median(co)) & (cc < 3 * np.median(cc))
print(f"\nnights calibrated: original {len(orig)}, corrected {len(corr)}, common {len(common)} ({ok.sum()} after outlier drop)")
print("night     C_L orig   C_L corrected   change")
for d, a, b in zip(common, co, cc):
    print(f"{d}   {a:7.3f}    {b:7.3f}       {100*(b/a-1):+6.1f}%")
print(f"\nMEDIAN    {np.median(co[ok]):7.3f}    {np.median(cc[ok]):7.3f}       {100*(np.median(cc[ok])/np.median(co[ok])-1):+6.1f}%")
print(f"gap to C_L(cloud)={CL_CLOUD}:  original {100*(np.median(co[ok])/CL_CLOUD-1):+.1f}%  ->  "
      f"corrected {100*(np.median(cc[ok])/CL_CLOUD-1):+.1f}%")
print("OFFSET_CORRECTION_DONE")
