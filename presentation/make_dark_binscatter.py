# -*- coding: utf-8 -*-
"""Deck figure: terminal-hood dark binscatter for the THREE Payerne co-located instruments
(CHM15k A, CL31 B, CL61 C).  For each instrument the 2-D histogram of beta/r^2 (non-range-
corrected, calibrated) vs altitude is drawn with the per-gate median overlaid; a zoom row shows
the per-gate median offset shape (the electronic distortion invisible in a single profile).

Extends validation/paper/_hood_binscatter.py (which covered only CHM15k + CL61) with CL31.
Altitude on the Y axis; beta/r^2 on X.  Pooled over the 4 covered-telescope sessions (2026).
"""
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

# ident, calibration constant C_L (rcs_0 -> beta_att Mm^-1 sr^-1), per-panel x half-window, colour,
# and the 4 covered-telescope windows for that instrument.
INST = {
    "CHM15k (1064 nm)": ("A", 5e11, 1.5e-7, "#c02a2a", [
        ("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
        ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")]),
    "CL31 (910 nm)": ("B", 3.37e7, 6.0e-7, "#e07b1a", [
        ("2026-05-12 09:33", "2026-05-12 14:55"), ("2026-05-26 12:00", "2026-05-27 13:10"),
        ("2026-06-09 09:16", "2026-06-09 11:57"), ("2026-06-23 10:12", "2026-06-23 12:13")]),
    "CL61 (910 nm)": ("C", 1.0, 1.5e-7, "#2166ac", [
        ("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-26 11:45", "2026-05-27 13:15"),
        ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-23 10:00", "2026-06-23 12:25")]),
}
ZMAX = 12000.0


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


fig = plt.figure(figsize=(16, 8.6))
gs = fig.add_gridspec(2, 3, height_ratios=[2.3, 1.0], hspace=0.30, wspace=0.24)

for k, (inst, (ident, CL, xlim, col, wins)) in enumerate(INST.items()):
    pool, rng = [], None
    for s1, s2 in wins:
        X, rng = load(ident, datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
        if X is not None:
            pool.append(X)
    if not pool:
        print(f"  {inst}: NO hood profiles found -> skipped")
        continue
    X = np.vstack(pool)
    betar2 = (X / CL) * 1e6 / (rng[None, :] ** 2)                 # Mm^-1 sr^-1 m^-2
    med = np.nanmedian(betar2, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(betar2 - med[None, :]), axis=0)
    se = mad / np.sqrt(np.maximum(np.isfinite(betar2).sum(axis=0), 1))
    dr = np.median(np.diff(rng)); kk = int(round(300 / dr)); kk += 1 - kk % 2
    meds = medfilt(np.nan_to_num(med), max(kk, 3))

    # --- top: per-sample binscatter (2-D histogram) with the per-gate median ---
    ax = fig.add_subplot(gs[0, k])
    xx = betar2.ravel(); yy = np.repeat(rng[None, :], X.shape[0], axis=0).ravel()
    m = np.isfinite(xx) & (np.abs(xx) <= xlim) & (yy <= ZMAX)
    h = ax.hist2d(xx[m], yy[m], bins=[np.linspace(-xlim, xlim, 90), np.linspace(0, ZMAX, 90)],
                  cmap="Blues", cmin=1)
    ax.axvline(0, color="0.45", lw=0.8)
    ax.axhspan(3000, 5000, color="orange", alpha=0.13, lw=0)                     # Rayleigh fit window
    ax.plot(med, rng, ".", ms=1.5, color="0.25")                                # per-gate median
    ax.plot(meds, rng, "-", color=col, lw=1.7)
    ax.set_xlim(-xlim, xlim); ax.set_ylim(0, ZMAX)
    ax.set_title(f"{inst}\n{X.shape[0]} pooled dark profiles", fontweight="bold", fontsize=11)
    ax.set_xlabel(r"$\beta/r^2$  [Mm$^{-1}$sr$^{-1}$ m$^{-2}$]")
    if k == 0:
        ax.set_ylabel("Altitude a.g.l. [m]")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))

    # --- bottom: ZOOM on the per-gate median (the electronic-distortion shape) +/- 3x precision ---
    axz = fig.add_subplot(gs[1, k])
    hx = np.nanpercentile(np.abs(meds[(rng >= 400) & (rng <= ZMAX)]), 99) * 1.4 + 1e-12
    axz.axvline(0, color="0.45", lw=0.8)
    axz.axhspan(3000, 5000, color="orange", alpha=0.13, lw=0)
    axz.fill_betweenx(rng, meds - 3 * se, meds + 3 * se, color=col, alpha=0.30, lw=0)
    axz.plot(meds, rng, "-", color=col, lw=1.7, label="per-gate median")
    m35 = np.nanmedian(meds[(rng >= 3000) & (rng <= 5000)])
    axz.set_xlim(-hx, hx); axz.set_ylim(0, ZMAX)
    axz.set_xlabel(r"median $\beta/r^2$ (zoom)")
    if k == 0:
        axz.set_ylabel("Altitude a.g.l. [m]")
    axz.set_title(f"3–5 km median = {m35:+.1e}", fontsize=9)
    axz.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    print(f"  {inst}: {X.shape[0]} profiles, 3-5 km median beta/r2 = {m35:+.2e}")

fig.suptitle("Payerne terminal-hood dark measurements — per-sample $\\beta/r^2$ binscatter with per-gate median "
             "(3 co-located instruments, 4 pooled sessions 2026)\n"
             "top: symmetric per-sample noise dominates (median ~0);  bottom: the per-gate median exposes the "
             "instrument-specific electronic distortion (invisible in a single profile)",
             fontweight="bold", fontsize=12)
OUT = "presentation/figs_intercomp/fig_dark_binscatter_payerne.png"
fig.savefig(OUT, dpi=160, bbox_inches="tight")
print("saved", OUT)
