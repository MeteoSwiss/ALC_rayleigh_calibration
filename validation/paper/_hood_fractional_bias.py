# -*- coding: utf-8 -*-
"""Does the hood offset MATTER? Per instrument, express the pooled hood offset as a FRACTION of
that instrument's own clear-night molecular signal at the same heights (self-consistent rcs_0
units -> the ratio is dimensionless and comparable across CHM15k / CL31 / CL61).
frac(z) = P_hood(z) / P_clear(z) is the fractional signal bias = the C_L bias for a fit at z.
All four hood sessions pooled (per-instrument windows); clear reference = a deep clear night."""
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
HOOD = {
    "CHM15k": ("A", "#d62728", [("2026-05-12 09:24", "2026-05-12 14:53"),
                                 ("2026-05-26 12:00", "2026-05-27 13:15"),
                                 ("2026-06-09 09:16", "2026-06-09 11:55"),
                                 ("2026-06-23 12:32", "2026-06-23 15:09")]),
    "CL31":   ("B", "#9467bd", [("2026-05-12 09:33", "2026-05-12 14:55"),
                                 ("2026-05-26 12:00", "2026-05-27 13:10"),
                                 ("2026-06-09 09:16", "2026-06-09 11:57"),
                                 ("2026-06-23 10:12", "2026-06-23 12:13")]),
    "CL61":   ("C", "#1f77b4", [("2026-05-12 09:17", "2026-05-12 14:58"),
                                 ("2026-05-26 11:45", "2026-05-27 13:15"),
                                 ("2026-06-09 09:10", "2026-06-09 11:51"),
                                 ("2026-06-23 10:00", "2026-06-23 12:25")]),
}
# deep clear night (low background, molecular return at 3-5 km); fallbacks if a file is absent
CLEAR = [("2026-04-22 20:00", "2026-04-23 03:30"), ("2026-04-24 20:00", "2026-04-25 03:30"),
         ("2026-04-02 20:00", "2026-04-03 03:30")]

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

def med_P(ident, wins):
    pool, rng = [], None
    for s1, s2 in wins:
        X, rng = load(ident, datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
        if X is not None:
            pool.append(X)
    if not pool:
        return None, None
    X = np.vstack(pool); zkm = rng / 1000.0
    return np.nanmedian(X / zkm[None, :] ** 2, axis=0), rng

def clear_ref(ident):
    """Median molecular P over ALL available clear nights + per-gate robust noise (MAD)."""
    prof, rng = [], None
    for c1, c2 in CLEAR:
        Xc, rc = load(ident, datetime.strptime(c1, "%Y-%m-%d %H:%M"), datetime.strptime(c2, "%Y-%m-%d %H:%M"))
        if Xc is not None:
            zkm = rc / 1000.0; prof.append(Xc / zkm[None, :] ** 2); rng = rc
    if not prof:
        return None, None, None
    A = np.vstack(prof)
    Pc = np.nanmedian(A, axis=0)
    noise = 1.4826 * np.nanmedian(np.abs(A - Pc[None, :]), axis=0) / np.sqrt(A.shape[0])
    return Pc, noise, rng

fig, ax = plt.subplots(1, 3, figsize=(17, 6.2))
print("instrument   hood N-sessions   frac bias @3km  @4km  @5km   (offset / clear molecular)")
for k, (inst, (ident, col, wins)) in enumerate(HOOD.items()):
    Ph, rng = med_P(ident, wins)
    Pc, noise, _ = clear_ref(ident)
    dr = np.median(np.diff(rng)); kk = int(round(300 / dr)); kk += 1 - kk % 2
    Phs = medfilt(np.nan_to_num(Ph), kk); Pcs = medfilt(np.nan_to_num(Pc), kk)
    a = ax[k]; a.axhline(0, color="0.6", lw=0.8); a.axhspan(3000, 5000, color="orange", alpha=0.10)
    # where the molecular signal falls below ~3x its own noise, a fractional bias is meaningless
    usable = np.abs(Pcs) > 3 * medfilt(np.nan_to_num(noise), kk)
    frac = 100.0 * Phs / np.where(usable, Pcs, np.nan)
    a.plot(frac, rng, "-", color=col, lw=1.8)
    # mark the altitude above which the molecular signal is lost in noise (short-range instruments)
    lost = rng[usable].max() if usable.any() else 0
    if lost < 11000:
        a.axhspan(lost, 12000, color="0.6", alpha=0.12)
        a.text(0.5, 0.90, f"molecular < 3σ above {lost/1000:.1f} km\n→ no Rayleigh fit possible",
               transform=a.transAxes, ha="center", va="top", fontsize=9, color="0.25")
    a.set_xlim(-60, 60); a.set_ylim(0, 8000 if inst == "CL31" else 12000)
    a.set_title(f"{inst}  (ident {ident})"); a.set_xlabel("hood offset / clear molecular  [%]")
    a.set_ylabel("range [m]"); a.grid(alpha=0.3)
    band = (rng >= 3000) & (rng <= 5000) & usable
    b35 = np.nanmean(frac[band]) if band.any() else float("nan")
    txt35 = f"3-5 km bias: {b35:+.0f}%" if band.any() else "3-5 km: no molecular signal"
    a.text(0.5, 0.02, txt35, transform=a.transAxes, ha="center", fontsize=11, fontweight="bold",
           bbox=dict(boxstyle="round", fc="white", ec=col))
    print(f"{inst:9s}   {len(wins)} sessions   " + (f"{b35:+.0f}% (3-5 km)" if band.any() else "no molecular signal at 3-5 km"))
fig.suptitle("Does the hood offset matter? Offset as a fraction of each instrument's own clear-night "
             "molecular signal (median of 3 nights). Orange = Rayleigh window 3-5 km; grey = molecular lost in noise.",
             fontweight="bold", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_fractional_bias.png"
fig.savefig(OUT, dpi=150); print("saved", OUT)
