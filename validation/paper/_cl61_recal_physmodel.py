# -*- coding: utf-8 -*-
"""Payerne CL61 Rayleigh recalibration with the PHYSICAL offset model, vs original and
vs the empirical linear ramp. Proves the AC-high-pass model (fig_cl61_offset_physical_model)
is not just a better curve fit but a better calibration.

Correction applied to L1 rcs_0 (gate-by-gate, beta-space):
  (0) original           : no correction
  (1) linear ramp        : min(P_lin*z^2, 0), valid 0.3-12 km  (the current report value)
  (2) physical model     : (Ap e^-r/Lp - Au e^-r/Lu) * z^2, DC floor removed, valid 0.35-14 km
Compared against the CHM15k-anchored cloud constant C_L = 1.4252.
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
NPZ = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz"

d = np.load(NPZ)
rng_d = d["rng"].astype(float); zkm = rng_d / 1000.0
# (1) linear ramp (as in the report): clamp <=0, valid 0.3-12 km
b_lin = np.minimum(d["b_linfit"].astype(float), 0.0)
b_lin[(rng_d < 300) | (rng_d > 12000)] = 0.0
# (2) physical model, DC floor removed so the far field -> 0 (no z^2 blow-up of b_inf)
Ap, Lp, Au, Lu, binf = [float(x) for x in d["phys_params"]]
P_phys = Ap * np.exp(-rng_d / Lp) - Au * np.exp(-rng_d / Lu)      # b_inf removed
b_phys = P_phys * zkm ** 2
b_phys[(rng_d < 350) | (rng_d > 14000)] = 0.0

def build(corr_beta, prefix):
    root = Path(tempfile.mkdtemp(prefix=prefix))
    for ds in NIGHTS:
        src = L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not src.exists():
            continue
        dst = root / WMO / ds[:4] / ds[4:6] / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, dst)
        with Dataset(dst, "r+") as nc:
            r = np.asarray(nc.variables["range"][:], "f8")
            corr = np.interp(r, rng_d, corr_beta) * 1e-6
            v = nc.variables["rcs_0"]; x = np.asarray(v[:], "f8")
            v[:] = (x - corr[None, :]) if x.shape[-1] == r.size else (x.T - corr[None, :]).T
    return root

ROOT_LIN = build(b_lin, "cl61_lin_")
ROOT_PHY = build(b_phys, "cl61_phy_")

def run(root):
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

orig = run(L1_ROOT); lin = run(ROOT_LIN); phy = run(ROOT_PHY)
common = sorted(set(orig) & set(lin) & set(phy))
print("night      original   linear-ramp   physical")
for ds in common:
    print(f"{ds}   {orig[ds]:8.3f}   {lin[ds]:8.3f}    {phy[ds]:8.3f}")

def summ(name, dct):
    v = np.array([dct[k] for k in common])
    v = v[v < 3 * np.median(v)]
    med = np.median(v); sc = 1.4826 * np.median(np.abs(v - med))
    print(f"  {name:14s} n={v.size}  median C_L={med:.3f}  scatter={100*sc/med:4.1f}%  "
          f"gap to cloud({CL_CLOUD}) = {100*(med/CL_CLOUD-1):+.1f}%")
print("\nsummary (common nights):")
summ("original", orig); summ("linear ramp", lin); summ("physical", phy)
print("RECAL_PHYS_DONE")
