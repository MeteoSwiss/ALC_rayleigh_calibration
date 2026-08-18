# -*- coding: utf-8 -*-
"""First look at the covered-telescope (hood) offset for ALL THREE Payerne instruments
(CHM15k _A, CL31 _B, CL61 _C), from the operator log of four hood sessions. Different rcs_0
units per instrument (counts/s.m2 / V.m2 / m-1sr-1) -> compare SHAPE and SIGN, not magnitude.
Non-range-corrected space P = rcs_0 / z^2 (homoscedastic). 23 Jun _A/_B not yet local."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
# per-instrument hood windows (UTC) from the operator log; 26 May = ~24 h (26->27 May)
SESS = {
    "CHM15k": ("A", "#d62728", [("2026-05-12 09:24", "2026-05-12 14:53"),
                                 ("2026-05-26 12:00", "2026-05-27 13:15"),
                                 ("2026-06-09 09:16", "2026-06-09 11:55")]),
    "CL31":   ("B", "#9467bd", [("2026-05-12 09:33", "2026-05-12 14:55"),
                                 ("2026-05-26 12:00", "2026-05-27 13:10"),
                                 ("2026-06-09 09:16", "2026-06-09 11:57")]),
    "CL61":   ("C", "#1f77b4", [("2026-05-12 09:17", "2026-05-12 14:58"),
                                 ("2026-05-26 11:45", "2026-05-27 13:15"),
                                 ("2026-06-09 09:10", "2026-06-09 11:51")]),
}
SLAB = ["12 May", "26-27 May (24 h)", "9 Jun"]

def load(ident, t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, rng = [], None
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
            X.append(x[s]); rng = r
    return (np.vstack(X), rng) if X else (None, None)

fig, ax = plt.subplots(1, 3, figsize=(17, 6.4))
print("instrument  session            N     P@0.5km    P@1km     P@3km   (non-range-corr, native units)")
for k, (inst, (ident, col, wins)) in enumerate(SESS.items()):
    a = ax[k]; a.axhline(0, color="0.6", lw=0.8)
    for i, (s1, s2) in enumerate(wins):
        t1 = datetime.strptime(s1, "%Y-%m-%d %H:%M"); t2 = datetime.strptime(s2, "%Y-%m-%d %H:%M")
        X, rng = load(ident, t1, t2)
        if X is None:
            continue
        zkm = rng / 1000.0
        P = X / zkm[None, :] ** 2                    # non-range-corrected
        Pm = np.nanmedian(P, axis=0)
        dr = np.median(np.diff(rng)); kk = int(round(300 / dr)); kk += 1 - kk % 2
        Ps = medfilt(np.nan_to_num(Pm), kk)
        # normalise each session's profile to its |value| at 0.5 km so shapes overlay
        ref = abs(Ps[np.argmin(np.abs(rng - 500))]) or 1.0
        a.plot(Ps / ref, rng, lw=1.6, alpha=0.85,
               color=plt.cm.viridis(i / 2.0), label=SLAB[i])
        at = lambda z: Ps[np.argmin(np.abs(rng - z))]
        print(f"{inst:9s}  {SLAB[i]:17s} {X.shape[0]:5d}  {at(500):+9.3g} {at(1000):+9.3g} {at(3000):+9.3g}")
    a.set_title(f"{inst}  (ident {ident})"); a.set_xlabel("P / |P(0.5 km)|  (shape, sign)")
    a.set_ylabel("range [m]"); a.set_ylim(0, 15000 if inst != "CL31" else 8000)
    a.set_xlim(-3, 3); a.grid(alpha=0.3); a.legend(fontsize=9, title="hood session")
fig.suptitle("Payerne covered-telescope offset SHAPE per instrument (non-range-corrected, "
             "normalised to |P(0.5 km)|) — does the anchor (CHM15k) and the CL31 show the CL61 undershoot?",
             fontweight="bold", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_multi_instrument.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
