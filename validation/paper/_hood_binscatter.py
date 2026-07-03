# -*- coding: utf-8 -*-
"""Reproduce (and improve) the MATLAB terminal-hood dark binscatter: for CHM15k and CL61, the
2-D histogram of beta/r^2 (non-range-corrected) vs range with the per-gate median overlaid, plus
the histogram of the per-gate medians. Pooled over the four 2026 hood sessions (more robust than a
single session). beta_att = rcs_0 / C_L (CL61: rcs_0 already m^-1 sr^-1 -> C_L=1; CHM15k: C_L=5e11
per the Payerne unit, as in the MATLAB). beta/r^2 in Mm^-1 sr^-1 m^-2 to match the MATLAB axes.
Improvement: per-gate median precision band (MAD/sqrt N) and the 3-5 km Rayleigh window marked."""
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
from matplotlib.colors import Normalize

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
INST = {"chm15k": ("A", 5e11, [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
                                ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")]),
        "CL61":   ("C", 1.0,  [("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-26 11:45", "2026-05-27 13:15"),
                                ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-23 10:00", "2026-06-23 12:25")])}

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

fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.28, wspace=0.22)
XLIM = 2.5e-7
for k, (inst, (ident, CL, wins)) in enumerate(INST.items()):
    pool, rng = [], None
    for s1, s2 in wins:
        X, rng = load(ident, datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
        if X is not None:
            pool.append(X)
    X = np.vstack(pool)
    betar2 = (X / CL) * 1e6 / (rng[None, :] ** 2)          # Mm^-1 sr^-1 m^-2 (MATLAB convention)
    med = np.nanmedian(betar2, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(betar2 - med[None, :]), axis=0)
    se = mad / np.sqrt(np.isfinite(betar2).sum(axis=0))     # per-gate median precision
    dr = np.median(np.diff(rng)); kk = int(round(300 / dr)); kk += 1 - kk % 2
    meds = medfilt(np.nan_to_num(med), kk)

    ax = fig.add_subplot(gs[0, k])
    xx = betar2.ravel(); yy = np.repeat(rng[None, :], X.shape[0], axis=0).ravel()
    m = np.isfinite(xx) & (np.abs(xx) <= XLIM) & (yy <= 15000)
    h = ax.hist2d(xx[m], yy[m], bins=[np.linspace(-XLIM, XLIM, 90), np.linspace(0, 15000, 90)],
                  cmap="Blues")
    ax.axvline(0, color="0.5", lw=0.8)
    ax.axhspan(3000, 5000, color="orange", alpha=0.12)                 # Rayleigh fit window
    ax.plot(med, rng, ".", ms=1.4, color="k")                          # per-gate median (MATLAB)
    ax.fill_betweenx(rng, meds - 3 * se, meds + 3 * se, color="#d62728", alpha=0.35, lw=0)  # +-3 SE
    ax.plot(meds, rng, "-", color="#d62728", lw=1.4, label="median +- 3x precision")
    ax.set_xlim(-XLIM, XLIM); ax.set_ylim(0, 15000)
    ax.set_title(f"{inst}   (pooled {X.shape[0]} dark profiles)", fontweight="bold")
    ax.set_xlabel(r"$\beta/r^2$  [Mm$^{-1}$sr$^{-1}$ m$^{-2}$]"); ax.set_ylabel("Range [m]")
    fig.colorbar(h[3], ax=ax, label="Bin counts")
    ax.legend(loc="upper right", fontsize=8)

    axh = fig.add_subplot(gs[1, k])
    HX = 2.5e-9                       # zoom ~100x: the offset lives here, far below the per-sample noise
    # STACKED histogram of the per-gate medians, split by 1 km altitude layer (which altitudes carry
    # the negative offset?). Colour = altitude; near-field gates (|median|>HX) go off-scale.
    edges = np.arange(0, 15001, 1000); nb = len(edges) - 1
    cmap = matplotlib.cm.turbo
    stacks, cols = [], []
    for i in range(nb):
        ma = (rng >= edges[i]) & (rng < edges[i + 1])
        v = med[ma]; v = v[np.isfinite(v) & (np.abs(v) <= HX)]
        stacks.append(v); cols.append(cmap(i / (nb - 1)))
    axh.hist(stacks, bins=np.linspace(-HX, HX, 70), stacked=True, color=cols, edgecolor="none")
    axh.axvline(0, color="0.5", lw=0.8)
    m35 = np.nanmedian(meds[(rng >= 3000) & (rng <= 5000)])
    axh.axvline(m35, color="k", lw=1.8, ls="--", label=f"3-5 km median = {m35:+.2e}")
    noff = int((np.abs(med[np.isfinite(med)]) > HX).sum())
    axh.set_xlim(-HX, HX); axh.set_xlabel(r"median($\beta/r^2$)  [Mm$^{-1}$sr$^{-1}$ m$^{-2}$]  (100$\times$ zoom)")
    axh.set_ylabel("Count of range gates")
    axh.set_title(f"per-gate medians, stacked by 1 km layer ({noff} near-field gates off-scale)", fontsize=9)
    axh.legend(fontsize=8, loc="upper left")
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 15))
    fig.colorbar(sm, ax=axh, label="altitude [km]", pad=0.01)

fig.suptitle("Terminal hood dark measurements — Payerne CHM15k & CL61 (4 sessions pooled, 2026)\n"
             "per-sample beta/r^2 (2-D histogram) with the per-gate median (black) and its +-3x precision band (red)",
             fontweight="bold", fontsize=12)
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_binscatter.png"
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("saved", OUT, "; per-gate median at 3-5 km resolves at >3x precision where the red band clears 0")
