# -*- coding: utf-8 -*-
"""Step 1 (Rayleigh v2) — two intercomparison figures: OPERATIONAL v1.0 vs v2 Rayleigh calibration.

Legacy = the operational calibration NetCDFs the user provided (ALC_calibration_..._A/C 2026.nc),
Kalman-smoothed with the SAME operational bridge the v2 series uses (calib_benchmark.kalman), written
to a calib_legacy/ dir. v2 = the existing calib/ series. Only the Rayleigh channels switch source;
the cloud channels are identical in both figures. Reuses the ladder pipeline (2x2 layout, fixed axes).
"""
import sys, warnings, csv
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset

import presentation.make_ladder_figures as M
from validation.paper import calib_benchmark as CB

NC = {"A": r"C:/Users/hervo/OneDrive/Documents/ALC_calibration_0-20000-0-06610_A2026.nc",
      "C": r"C:/Users/hervo/OneDrive/Documents/ALC_calibration_0-20000-0-06610_C2026.nc"}
LEG = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib_legacy")
LEG.mkdir(parents=True, exist_ok=True)
OUT = Path(r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/presentation/figs_intercomp/ladder")


def nc_to_kalman_csv(ident, key):
    """Read the operational NetCDF lidar constants, Kalman-smooth (operational bridge), write CSV."""
    with Dataset(NC[ident]) as nc:
        t = np.ma.filled(nc.variables["time"][:].astype("f8"), np.nan)
        lc = np.ma.filled(nc.variables["lidar_constant"][:].astype("f8"), np.nan)
        u = np.ma.filled(nc.variables["lidar_constant_uncertainty"][:].astype("f8"), np.nan)
    dates = [datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t]
    ok = np.isfinite(lc) & (lc > 0)
    dates = [d for d, o in zip(dates, ok) if o]
    C, Cstd = lc[ok], u[ok]
    print(f"  {key}: {len(C)} operational nights  median C_L={np.median(C):.3e}  "
          f"span {dates[0].date()}..{dates[-1].date()}" if len(C) else f"  {key}: no nights")
    res = CB.kalman(dates, C, Cstd, normalise=True)         # operational Kalman bridge (>=3 nights)
    if res is None:
        print(f"    -> Kalman returned None ({len(C)} nights < 3) — no legacy series for {key}")
        return False
    with open(LEG / f"{key}_L1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
        w.writerows(res)
    nk = sum(1 for r in res if np.isfinite(r[3]))
    print(f"    -> {LEG.name}/{key}_L1.csv  ({len(res)} daily rows, {nk} finite Kalman)")
    return True


print("== Kalman-smoothing the operational calibration NetCDFs ==")
have_chm = nc_to_kalman_csv("A", "0-20000-0-06610_A_rayleigh")
have_cl61 = nc_to_kalman_csv("C", "0-20000-0-06610_C_rayleigh")

# ---- render the two figures (reuse the ladder pipeline) ----
cache = M.load_l1_cache(M.L1_NPZ)
cfg = M._cfg(wl="none", wv=False)              # WV OFF, wavelength conversion OFF (pure Rayleigh, no harmonisation)

CHM = M._ch("A", "CHM15k", "rayleigh", "CHM15k (Rayleigh)")
CL31 = M._ch("B", "CL31", "cloud", "CL31 (cloud)")
CL61C = M._ch("C", "CL61", "cloud", "CL61 (cloud)")
CL61R = M._ch("C", "CL61", "rayleigh", "CL61 (Rayleigh)")


def render(tag, title, chm_dir, cl61r_dir):
    chans = [M.cl(CHM, chm_dir)]          # cloud calibration OFF -> only the Rayleigh channels remain
    if cl61r_dir is not None and (Path(cl61r_dir) / "0-20000-0-06610_C_rayleigh_L1.csv").exists():
        chans.append(M.cl(CL61R, cl61r_dir))
    else:
        print(f"  [{tag}] CL61 (Rayleigh) absent for this source -> CHM15k only")
    R = M.build_R(cache, chans, cfg)
    if R is None:
        print(f"  [{tag}] build failed"); return
    M.fig_ladder(R, OUT / f"{tag}.png",
                 f"Payerne (0-20000-0-06610)  —  {title}  —  2026-03-01 to 2026-06-30")
    for k, ch in enumerate(chans):
        st = R["stats"][k]
        print("     %-22s med %+6.1f%%  logr %.2f  N=%d" % (ch["label"], st["medrelbias_pct"], st["r_log"], st["n"]))
    print(f"  -> {tag}.png")


print("\n== rendering legacy (operational v1.0) — WV off, wavelength off, cloud off ==")
render("step1c_operational_raw", "Rayleigh v1.0 (operational)  —  no WV, no wavelength conv, no cloud",
       LEG, LEG if have_cl61 else None)
print("\n== rendering v2 — WV off, wavelength off, cloud off ==")
render("step1d_v2_raw", "Rayleigh v2  —  no WV, no wavelength conv, no cloud",
       M.CALIB_NEW, M.CALIB_NEW)
print("STEP1_DONE")
