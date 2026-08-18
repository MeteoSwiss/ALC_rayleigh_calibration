# -*- coding: utf-8 -*-
"""Verify the CHM15k (anchor) hood offset is a REAL calibration bias, not a hood artifact:
characterise b_dark(z) from the four covered-telescope windows, subtract it from the clean-night
L1 rcs_0, and re-run the CHM15k Rayleigh calibration. If C_L shifts by ~the fractional bias
(+~25%), the offset biases the CHM15k molecular calibration exactly as it does the CL61.
Same procedure proven on CL61 (_cl61_offset_correction_experiment)."""
import sys, shutil, tempfile, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.signal import medfilt

from calibration.config import CalibrationOptions, InstrumentInfo, InstrumentType, DataLevel
from calibration.rayleigh.calibration import calibrate_rayleigh
from validation.paper.calib_benchmark import RAYLEIGH_METHOD, L1_ROOT, CAMS, REPO

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "A", 46.813, 6.943, 491.0
HOOD = [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
        ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")]
NIGHTS = ["20260225", "20260316", "20260328", "20260402", "20260407", "20260410", "20260411",
          "20260417", "20260420", "20260422", "20260423", "20260424", "20260428", "20260521", "20260602"]

def load(t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, rng = [], None
    for ds in sorted(days):
        f = L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            tv = np.asarray(nc.variables["time"][:], "f8")
            tt = np.array([datetime(1970, 1, 1) + timedelta(days=x) for x in tv])
            r = np.asarray(nc.variables["range"][:], "f8")
            x = np.asarray(nc.variables["rcs_0"][:], "f8")
            if x.shape[0] != tt.size:
                x = x.T
        s = (tt >= t1) & (tt <= t2)
        if s.any():
            X.append(x[s]); rng = r
    return (np.vstack(X), rng) if X else (None, None)

# ---- b_dark(z): pooled hood median in non-range-corrected space, smoothed, x z^2 ----
pool, rng = [], None
for s1, s2 in HOOD:
    X, rng = load(datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
    if X is not None:
        pool.append(X)
Xh = np.vstack(pool); zkm = rng / 1000.0
P = Xh / zkm[None, :] ** 2
Pm = np.nanmedian(P, axis=0)
dr = np.median(np.diff(rng)); kk = int(round(300 / dr)); kk += 1 - kk % 2
b_dark = medfilt(np.nan_to_num(Pm), kk) * zkm ** 2          # back to rcs_0 (range-corrected) units
b_dark[rng < 1000] = 0.0                                    # near-field/overlap: leave alone
print(f"CHM15k b_dark (rcs_0 units): 3 km {b_dark[np.argmin(abs(rng-3000))]:+.4g}  "
      f"4 km {b_dark[np.argmin(abs(rng-4000))]:+.4g}  5 km {b_dark[np.argmin(abs(rng-5000))]:+.4g}  "
      f"(N hood profiles {Xh.shape[0]})")

# ---- corrected archive: rcs_0' = rcs_0 - b_dark(z) ----
CORR = Path(tempfile.mkdtemp(prefix="chm_corr_"))
for ds in NIGHTS:
    src = L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
    if not src.exists():
        continue
    dst = CORR / WMO / ds[:4] / ds[4:6] / src.name
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    with Dataset(dst, "r+") as nc:
        r = np.asarray(nc.variables["range"][:], "f8")
        corr = np.interp(r, rng, b_dark)
        v = nc.variables["rcs_0"]; x = np.asarray(v[:], "f8")
        v[:] = (x - corr[None, :]) if x.shape[-1] == r.size else (x.T - corr[None, :]).T

def run(root):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = root; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = True
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                          instrument_type=InstrumentType.CHM15k, latitude=LAT, longitude=LON, altitude=ALT)
    out = {}
    for ds in NIGHTS:
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception:
            continue
        if r.flag in (1, 0.5) and np.isfinite(r.lidar_constant) and r.lidar_constant > 0:
            out[ds] = float(r.lidar_constant)
    return out

orig = run(L1_ROOT); corr = run(CORR)
common = sorted(set(orig) & set(corr))
print(f"\nnights calibrated: original {len(orig)}, corrected {len(corr)}, common {len(common)}")
print("night      C_L orig     C_L corrected   change")
for d in common:
    print(f"{d}   {orig[d]:.4g}    {corr[d]:.4g}     {100*(corr[d]/orig[d]-1):+6.1f}%")
co = np.array([orig[d] for d in common]); cc = np.array([corr[d] for d in common])
ok = (co < 3 * np.median(co)) & (cc < 3 * np.median(cc))
print(f"\nMEDIAN   {np.median(co[ok]):.4g}    {np.median(cc[ok]):.4g}     "
      f"{100*(np.median(cc[ok])/np.median(co[ok])-1):+.1f}%   (fractional-bias prediction ~ +15..+25%)")
print("CHM_VERIFY_DONE")
