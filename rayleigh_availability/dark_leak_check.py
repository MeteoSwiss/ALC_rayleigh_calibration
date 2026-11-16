# -*- coding: utf-8 -*-
"""Is the measured Payerne dark b(z) contaminated by atmospheric signal leaking through the hood?

For each session (data-verified windows, dark_windows_detected.json) and instrument (CHM15k A, CL61 C):
  d(z) = median rcs_0 over the covered window                      (the "dark")
  a(z) = median rcs_0 over the 60 min before + 60 min after it      (the atmosphere around it)
A leak makes d(z) contain k*a(z), so (i) d(z) follows a(z) session by session -- in particular its
structure moves with each session's boundary-layer top -- and (ii) the regression of d on a over
150-3000 m gives a consistent positive k. A pure electronic baseline is the same in every session
whatever the atmosphere. The 24-h session (26-27 May) is also split into day and night halves.

Run:  python rayleigh_availability/dark_leak_check.py   -> figure + table on stdout
"""
from __future__ import annotations
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
OUT_FIG = REPO / "doc" / "reports" / "figs_remote_dark" / "dark_leak_check.png"
AROUND = timedelta(minutes=60)
GAP = timedelta(minutes=5)
NAMES = {"A": "CHM15k", "C": "CL61"}


def read(ident, t0, t1):
    T, X, rng = [], [], None
    day = t0.date()
    while day <= t1.date():
        f = L1 / f"{day:%Y}" / f"{day:%m}" / f"L1_0-20000-0-06610_{ident}{day:%Y%m%d}.nc"
        if f.exists():
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                x = np.ma.filled(ds["rcs_0"][:].astype("f8"), np.nan)
            if x.shape[0] != t.size:
                x = x.T
            tt = np.array([datetime(1970, 1, 1) + timedelta(days=float(v)) for v in t])
            m = (tt >= t0) & (tt <= t1)
            T.append(tt[m])
            X.append(x[m])
        day += timedelta(days=1)
    return np.concatenate(T), np.vstack(X), rng


def smooth(z, y, w_m=150.0):
    dz = float(np.median(np.diff(z)))
    k = max(3, int(round(w_m / dz)) | 1)
    yy = np.where(np.isfinite(y), y, 0.0)
    return np.convolve(yy, np.ones(k) / k, mode="same")


def analyse(ident, w0, w1):
    t0 = datetime.strptime(w0, "%Y-%m-%d %H:%M:%S")
    t1 = datetime.strptime(w1, "%Y-%m-%d %H:%M:%S")
    t, x, z = read(ident, t0 - AROUND - GAP, t1 + AROUND + GAP)
    dark = (t >= t0) & (t <= t1)
    atm = ((t >= t0 - AROUND - GAP) & (t <= t0 - GAP)) | ((t >= t1 + GAP) & (t <= t1 + AROUND + GAP))
    d = smooth(z, np.nanmedian(x[dark], axis=0))
    a = smooth(z, np.nanmedian(x[atm], axis=0))
    band = (z >= 150) & (z <= 3000)
    A = np.vstack([a[band], np.ones(band.sum())]).T
    k, c = np.linalg.lstsq(A, d[band], rcond=None)[0]                 # d = k*a + c over 150-3000 m
    r = np.corrcoef(a[band], d[band])[0, 1]
    # boundary-layer top proxy of the surrounding atmosphere: height of the strongest negative
    # gradient of a(z) between 300 and 3000 m
    g = np.gradient(a, z)
    sel = (z >= 300) & (z <= 3000)
    zbl = float(z[sel][np.argmin(g[sel])])
    gd = np.gradient(d, z)
    zd = float(z[sel][np.argmin(gd[sel])]) if ident == "C" else np.nan    # dark's own sharpest drop
    return dict(z=z, d=d, a=a, k=k, c=c, r=r, zbl=zbl, zd=zd, n=int(dark.sum()), start=w0[:16])


def main():
    win = json.loads((DATA / "dark_windows_detected.json").read_text())
    res = {}
    for ident in ("A", "C"):
        res[ident] = []
        print(f"== {NAMES[ident]}")
        print("   session            n    k(leak)     r(d,a)  BL-top(atm)  dark drop")
        for w0, w1 in win[ident]:
            s = analyse(ident, w0, w1)
            res[ident].append(s)
            print(f"   {s['start']}  {s['n']:5d}  {s['k']:+.2e}  {s['r']:+.2f}   {s['zbl']:7.0f} m   "
                  f"{s['zd']:7.0f} m")
    # 24-h session split day / night (CL61 and CHM15k)
    split = {}
    for ident in ("A", "C"):
        w0, w1 = win[ident][1]
        t0 = datetime.strptime(w0, "%Y-%m-%d %H:%M:%S")
        t1 = datetime.strptime(w1, "%Y-%m-%d %H:%M:%S")
        t, x, z = read(ident, t0, t1)
        h = np.array([v.hour for v in t])
        night = (h >= 21) | (h < 3)
        day = (h >= 9) & (h < 15)
        split[ident] = (z, smooth(z, np.nanmedian(x[night], axis=0)), smooth(z, np.nanmedian(x[day], axis=0)))
    figure(res, split)


def figure(res, split):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(22, 10))
    for row, ident in enumerate(("A", "C")):
        cmap = plt.get_cmap("viridis")
        ss = res[ident]
        ax0, ax1, ax2, ax3 = axs[row]
        for i, s in enumerate(ss):
            col = cmap(i / max(len(ss) - 1, 1))
            m = s["z"] <= 6000
            ax0.plot(s["d"][m], s["z"][m], color=col, lw=1.2, label=s["start"][:10])
            ax1.plot(s["a"][m] / np.nanmax(s["a"][m]), s["z"][m], color=col, lw=1.2)
            ax1.axhline(s["zbl"], color=col, lw=0.6, ls=":")
            ax0.axhline(s["zbl"], color=col, lw=0.6, ls=":")
        ax0.axvline(0, color="k", lw=0.8)
        ax0.set_title(f"{NAMES[ident]}: dark d(z) per session (dotted: BL top of the surrounding atmosphere)")
        ax0.set_xlabel("rcs_0 (instrument units)")
        ax0.set_ylabel("Altitude a.g.l. [m]")
        ax0.legend(fontsize=7, ncol=2)
        ax1.set_title("atmosphere 1 h before/after (normalised)")
        ax1.set_xlabel("a(z) / max")
        # scatter of the CL61/CHM 'dark drop' vs BL top
        ax2.plot([s["zbl"] for s in ss], [s["k"] for s in ss], "o")
        ax2.set_xlabel("BL top of surrounding atmosphere [m]")
        ax2.set_ylabel("leak coefficient k (d = k a + c, 150-3000 m)")
        ax2.axhline(0, color="k", lw=0.8)
        ax2.set_title("leak coefficient per session")
        z, dn, dd = split[ident]
        m = z <= 6000
        ax3.plot(dn[m], z[m], label="night 21-03 UTC")
        ax3.plot(dd[m], z[m], label="day 09-15 UTC")
        ax3.axvline(0, color="k", lw=0.8)
        ax3.set_title("24-h session 26-27 May: night vs day dark")
        ax3.legend(fontsize=8)
        for ax in (ax0, ax1, ax3):
            ax.set_ylim(0, 6000)
            ax.grid(alpha=0.3)
        ax2.grid(alpha=0.3)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=110)
    print(f"-> {OUT_FIG}")


if __name__ == "__main__":
    main()
