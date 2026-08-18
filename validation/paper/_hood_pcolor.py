# -*- coding: utf-8 -*-
"""Time-height pcolor of the covered-telescope sessions (CL61 + CHM15k) to check for aerosol/cloud
contamination or hood on/off transitions inside the 'dark' windows. beta_att in Mm^-1 sr^-1
(CL61: rcs_0*1e6; CHM15k: rcs_0/5e11*1e6), diverging scale +-0.5 as in the operator MATLAB fig1:
a covered telescope should be uniform noise; any bright layer/period is contamination."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
SESS = [("12 May", ("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-12 09:24", "2026-05-12 14:53")),
        ("26-27 May (24 h)", ("2026-05-26 11:45", "2026-05-27 13:15"), ("2026-05-26 12:00", "2026-05-27 13:15")),
        ("9 Jun", ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-09 09:16", "2026-06-09 11:55")),
        ("23 Jun", ("2026-06-23 10:00", "2026-06-23 12:25"), ("2026-06-23 12:32", "2026-06-23 15:09"))]

def load(ident, t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, T, rng = [], [], None
    for ds in sorted(days):
        f = L1 / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_{ident}{ds}.nc"
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
            X.append(x[s]); T.append(tt[s]); rng = r
    return (np.vstack(X), np.concatenate(T), rng) if X else (None, None, None)

fig, ax = plt.subplots(2, 4, figsize=(21, 10))
print("contamination check: %% of profiles with median beta_att(1-4 km) > 0.15 Mm-1 sr-1")
for col, (lab, wcl61, wchm) in enumerate(SESS):
    for row, (inst, ident, CL, win) in enumerate([("CL61", "C", 1.0, wcl61), ("CHM15k", "A", 5e11, wchm)]):
        t1 = datetime.strptime(win[0], "%Y-%m-%d %H:%M"); t2 = datetime.strptime(win[1], "%Y-%m-%d %H:%M")
        X, T, rng = load(ident, t1, t2)
        a = ax[row, col]
        if X is None:
            a.set_visible(False); continue
        beta = (X / CL) * 1e6                                  # Mm^-1 sr^-1
        hrs = np.array([(t - T[0]).total_seconds() / 3600.0 for t in T])
        pm = a.pcolormesh(hrs, rng, beta.T, cmap="RdBu_r", vmin=-0.3, vmax=0.3, shading="auto")
        # contamination flag: a real return at 1-4 km (above the near-field, below the noise floor)
        band = (rng >= 1000) & (rng <= 4000)
        prof15 = np.nanmedian(beta[:, band], axis=1)
        contam = prof15 > 0.15
        a.plot(hrs[contam], np.full(contam.sum(), 700), "s", ms=3, color="lime", mec="k", mew=0.2)
        a.set_ylim(0, 6000)
        a.set_title(f"{inst} — {lab}  ({100*contam.mean():.0f}% contaminated)", fontsize=10, fontweight="bold")
        if col == 0:
            a.set_ylabel("Range [m]")
        if row == 1:
            a.set_xlabel(f"hours since {T[0]:%m-%d %H:%M} UTC")
        if col == 3:
            fig.colorbar(pm, ax=a, label=r"$\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]")
        print(f"  {inst:7s} {lab:18s}: {100*contam.mean():5.1f}%  (n={X.shape[0]}, {contam.sum()} flagged)")
fig.suptitle("Terminal-hood sessions — time-height β_att, zoom 0-6 km (green = profiles with a real "
             "1-4 km return > 0.15 Mm⁻¹sr⁻¹ = aerosol/cloud or hood on/off contamination)",
             fontweight="bold", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_pcolor.png"
fig.savefig(OUT, dpi=140); print("saved", OUT)
