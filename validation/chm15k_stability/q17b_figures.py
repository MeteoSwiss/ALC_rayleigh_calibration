# -*- coding: utf-8 -*-
"""q17b (figures): time series, height dependence and scaling scatter of the CHM15k / CL61 test.

Inputs: q17/q17b_nights.csv, q17b_bands.csv (q17b_nights.py, q17b_analyse.py). Landscape PNGs, altitude on Y.
"""
from __future__ import annotations

import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from q17b_extract import Q17, SITES  # noqa: E402
from q17b_nights import BAND_NAMES, BANDS  # noqa: E402
from q17b_analyse import R_MID, cl61_valid, flag_iqr, load_kappa_in, monthly, running_median  # noqa: E402

warnings.filterwarnings("ignore")
EPOCH = datetime(1970, 1, 1)


def series(df, site, bn, kap):
    d = df[df.site == site]
    t = d["t"].to_numpy(float)
    x = d[f"chm_{bn}"].to_numpy(float)
    y = d[f"cl_wv_{bn}"].to_numpy(float)
    y0 = d[f"cl_nowv_{bn}"].to_numpy(float)
    T = d["T_int"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(y0) & cl61_valid(site, t)
    t, x, y, y0, T = t[ok], x[ok], y[ok], y0[ok], T[ok]
    bad = flag_iqr(x) | flag_iqr(y)
    t, x, y, y0, T = t[~bad], x[~bad], y[~bad], y0[~bad], T[~bad]
    xk = x - kap / 100 * (T - np.median(T)) if np.isfinite(kap) else x
    return t, x, xk, y, y0


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = pd.read_csv(Q17 / "q17b_nights.csv")
    bands = pd.read_csv(Q17 / "q17b_bands.csv")
    k_in = load_kappa_in()
    sites = list(SITES)

    # 1. time series ------------------------------------------------------------------------
    fig, axs = plt.subplots(2, 4, figsize=(22, 8), sharex=False)
    for j, site in enumerate(sites):
        wmo, ich = SITES[site][0], SITES[site][1]
        kap = k_in.get(f"{wmo}_{ich}", np.nan)
        for i, bn in enumerate(("3-4", "6-8")):
            ax = axs[i, j]
            t, x, xk, y, y0 = series(df, site, bn, kap)
            if len(t) < 20:
                continue
            dt = [EPOCH + timedelta(days=float(v)) for v in t]
            for v, lab, c, ls in ((x, "CHM15k raw", "0.6", "-"), (xk, "CHM15k kappa-norm.", "k", "-"),
                                  (y, "CL61 WV-corrected", "C0", "-"), (y0, "CL61 no WV corr.", "C9", "--")):
                yy = 100 * (v - np.median(v))
                ax.plot(dt, yy, ".", ms=2, color=c, alpha=0.25)
                ax.plot(dt, 100 * running_median(t, v - np.median(v)), ls, color=c, lw=1.8, label=lab)
            ax.plot(dt, 100 * R_MID * running_median(t, xk - np.median(xk)), ":", color="C1", lw=2,
                    label=f"aerosol prediction {R_MID:.2f} x CHM15k")
            ax.axhline(0, color="k", lw=0.5)
            ax.set_ylim(-45, 45)
            ax.set_title(f"{site}: clear-sky ln X, {bn} km above the instrument", fontsize=9)
            if j == 0:
                ax.set_ylabel("ln X - median [%] (31-d running median)")
            if i == 0 and j == 0:
                ax.legend(fontsize=7, loc="lower left")
    fig.suptitle("q17b - CHM15k (1064 nm) vs CL61 (910.55 nm) on the same clear-sky night blocks, "
                 "CAMS molecular, CL61 WV-corrected", fontsize=10)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(Q17 / "q17b_series.png", dpi=110)
    plt.close(fig)

    # 2. height dependence --------------------------------------------------------------------
    zmid = np.array([0.5 * (a + b) / 1000 for a, b in BANDS])
    fig, axs = plt.subplots(1, 4, figsize=(22, 6.5), sharey=True)
    for j, site in enumerate(sites):
        ax = axs[j]
        b = bands[(bands.site == site) & (bands.ycol == "cl_wv")].set_index("band")
        b0 = bands[(bands.site == site) & (bands.ycol == "cl_nowv")].set_index("band")
        if b.empty:
            continue
        def col(frame, c):
            return np.array([frame[c].get(bn, np.nan) if c in frame else np.nan for bn in BAND_NAMES])
        ax_ = col(b, "amp_x")
        axk = col(b, "k_amp_x")
        # signed in-phase CL61 response on the CHM15k phase
        yin = col(b, "amp_y_inphase")
        yin0 = col(b0, "amp_y_inphase")
        bg = bands[(bands.site == site) & (bands.ycol == "cl_wv+cloudgain")].set_index("band")
        ba = bands[(bands.site == site) & (bands.ycol == "cl_wv+aeronet")].set_index("band")
        yg = col(bg, "amp_y_inphase")
        ya = col(ba, "amp_y_inphase")
        ax.plot(ax_, zmid, "o-", color="0.6", label="CHM15k raw")
        ax.plot(axk, zmid, "o-", color="k", label="CHM15k kappa-normalised")
        ax.plot(yin, zmid, "s-", color="C0", label="CL61 WV-corr., in phase with CHM15k")
        ax.plot(yin0, zmid, "s--", color="C9", label="CL61 no WV corr., in phase")
        ax.plot(yg, zmid, "^-", color="C2", label="CL61 WV-corr. / own cloud constant, in phase")
        ax.plot(ya, zmid, "v-", color="C4", label="CL61 WV + AERONET transmission corr., in phase")
        ax.plot(ax_ - yg / R_MID, zmid, "d-", color="C3", alpha=0.6,
                label="CHM15k gain part = CHM15k - CL61(cloud-gain)/r")
        ax.plot(R_MID * ax_, zmid, ":", color="C1", lw=2, label=f"aerosol prediction ({R_MID:.2f} x CHM15k raw)")
        ax.axvline(0, color="C3", lw=1, label="gain prediction (CL61 = 0)")
        ax.set_title(site, fontsize=10)
        ax.set_xlabel("annual amplitude of ln X [%]")
        ax.set_xlim(-15, 40)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("height above instrument [km]")
            ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("q17b - height dependence of the annual cycle (common clear nights)", fontsize=10)
    fig.tight_layout()
    fig.savefig(Q17 / "q17b_height.png", dpi=110)
    plt.close(fig)

    # 3. scaling scatter ------------------------------------------------------------------------
    fig, axs = plt.subplots(2, 4, figsize=(22, 9))
    for j, site in enumerate(sites):
        wmo, ich = SITES[site][0], SITES[site][1]
        kap = k_in.get(f"{wmo}_{ich}", np.nan)
        ax = axs[0, j]
        for bn, c in (("2-3", "C2"), ("4-5", "C4"), ("6-8", "C3")):
            t, x, xk, y, y0 = series(df, site, bn, kap)
            if len(t) < 20:
                continue
            mo = monthly(t, [x, y])
            if len(mo) < 3:
                continue
            xm = 100 * (mo[:, 1] - mo[:, 1].mean())
            ym = 100 * (mo[:, 2] - mo[:, 2].mean())
            ax.plot(xm, ym, "o", color=c, label=f"{bn} km monthly medians")
        lim = np.array([-35, 35])
        ax.plot(lim, R_MID * lim, ":", color="C1", lw=2, label=f"aerosol slope {R_MID:.2f}")
        ax.plot(lim, 0 * lim, "-", color="C3", lw=1, label="gain slope 0")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_title(f"{site}: seasonal (monthly) CL61 vs CHM15k", fontsize=9)
        ax.set_xlabel("CHM15k ln X anomaly [%]")
        if j == 0:
            ax.set_ylabel("CL61 ln X anomaly [%]")
            ax.legend(fontsize=7)
        ax = axs[1, j]
        for bn, c in (("2-3", "C2"), ("6-8", "C3")):
            t, x, xk, y, y0 = series(df, site, bn, kap)
            if len(t) < 20:
                continue
            a = 100 * (x - running_median(t, x))
            bb = 100 * (y - running_median(t, y))
            ax.plot(a, bb, ".", ms=3, color=c, alpha=0.5, label=f"{bn} km nightly anomalies")
        ax.plot(lim, R_MID * lim, ":", color="C1", lw=2)
        ax.plot(lim, 0 * lim, "-", color="C3", lw=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_title(f"{site}: synoptic anomalies (vs 31-d running median)", fontsize=9)
        ax.set_xlabel("CHM15k ln X anomaly [%]")
        if j == 0:
            ax.set_ylabel("CL61 ln X anomaly [%]")
            ax.legend(fontsize=7)
    fig.suptitle("q17b - wavelength scaling: aerosol predicts slope ~0.6, a CHM15k gain predicts 0", fontsize=10)
    fig.tight_layout()
    fig.savefig(Q17 / "q17b_scatter.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
