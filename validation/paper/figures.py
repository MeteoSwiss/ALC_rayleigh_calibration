"""
figures.py — Python reproductions of the MATLAB paper-validation figures
(make_validation_figures.m / paper_val_figure_payerne.m / paper_val_earlinet_figure.m):

  fig_calib_timeseries : one grid of ALL calibrated channels, raw daily (x) + Kalman (line +/-1 sigma)
  fig_multi_alc        : per multi-ALC site, the 3x3 layout
                         (a) median+/-IQR profiles | (b) scatter vs ref | (c) histogram of differences
                         + the four channel pcolors (d-g) in the lower-right 2x2 block
  fig_earlinet         : per EARLINET site, the 2x2 layout
                         (a) median matched profile +/-IQR | (b) density scatter | (c) EARLINET | (d) CHM curtains

Same subplot structure, colours, comparison band and noise-floor truncation as the MATLAB.
"""
from __future__ import annotations
import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

# MATLAB channel_colors base (blue, red, green, purple, cyan, dark-red, gold, grey)
COLORS = np.array([[0.00, 0.45, 0.74], [0.85, 0.33, 0.10], [0.47, 0.67, 0.19], [0.49, 0.18, 0.56],
                   [0.30, 0.75, 0.93], [0.64, 0.08, 0.18], [0.93, 0.69, 0.13], [0.25, 0.25, 0.25]])
CLIM = (-2.0, 1.0)
BETALIM = (1e-2, 1e2)


def _col(k):
    return COLORS[k % len(COLORS)]


def _median_iqr(B):
    """median, 25th, 75th percentile over axis 0 (profiles), NaN-aware."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(B, axis=0)
        q1 = np.nanpercentile(B, 25, axis=0)
        q3 = np.nanpercentile(B, 75, axis=0)
        nz = np.sum(np.isfinite(B), axis=0)
    return med, q1, q3, nz


def _truncate_noise_floor(med, nz, z, zmin, nprof, minfrac=0.05):
    """MATLAB: keep up to the first altitude above zMin where the median is no longer positive
    (or coverage too thin). Above that the screened signal is noise."""
    good = np.isfinite(med) & (med > 0) & (nz >= max(10, minfrac * nprof))
    keep = np.ones(med.size, bool)
    bad = np.where(~good & (z > zmin))[0]
    if bad.size:
        keep[bad[0]:] = False
    return keep & good


# ---------------------------------------------------------------------------
#  Multi-ALC 3x3 figure (profile | scatter | hist + 4 pcolors)
# ---------------------------------------------------------------------------
def fig_multi_alc(R, cfg, out_png, title, zmax_plot=6000):
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]      # m AGL
    zmask = (z >= 0) & (z <= zmax_plot)
    zc = z[zmask]
    iref = cfg["referenceChannel"]
    zmin, zmax = cfg["zMin"], cfg["zMax"]
    band = (z >= zmin) & (z <= zmax)
    nch = len(R["channels"])
    tx = mdates.date2num(np.asarray(R["time_sync"]).astype("datetime64[s]").astype(datetime))

    fig = plt.figure(figsize=(19, 10))
    gs = GridSpec(3, 3, figure=fig, hspace=0.28, wspace=0.26,
                  left=0.05, right=0.93, top=0.92, bottom=0.07)

    # (a) profile median +/- IQR — left column
    axp = fig.add_subplot(gs[:, 0])
    for k in range(nch):
        B = R["beta"][k][:, zmask]
        med, q1, q3, nz = _median_iqr(B)
        keep = _truncate_noise_floor(med, nz, zc, zmin, B.shape[0])
        c = _col(k)
        m1 = keep & np.isfinite(q1) & (q1 > 0)
        axp.plot(q1[m1], zc[m1], ":", color=c, lw=1.0)
        m3 = keep & np.isfinite(q3) & (q3 > 0)
        axp.plot(q3[m3], zc[m3], ":", color=c, lw=1.0)
        axp.plot(med[keep], zc[keep], "-", color=c, lw=1.8, label=R["channels"][k]["label"])
    axp.axhline(zmin, ls="--", color="k", lw=0.8); axp.axhline(zmax, ls="--", color="k", lw=0.8)
    axp.set_xscale("log"); axp.set_xlim(*BETALIM); axp.set_ylim(0, zmax_plot)
    axp.grid(alpha=0.3); axp.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]"); axp.set_ylabel("Altitude a.g.l. [m]")
    axp.set_title("(a) Median (solid) $\\pm$ IQR (dotted)  %s-%s"
                  % (_fmt_my(cfg["start"]), _fmt_my(cfg["end"])), fontsize=10)
    axp.legend(loc="upper right", fontsize=8)

    # (b) scatter vs reference
    axs = fig.add_subplot(gs[0, 1])
    ref = R["beta"][iref][:, band].ravel()
    allv = []
    for k in range(nch):
        if k == iref:
            continue
        cur = R["beta"][k][:, band].ravel()
        m = np.isfinite(cur) & np.isfinite(ref) & (cur > 0) & (ref > 0)
        a, b = ref[m], cur[m]
        if a.size > 6000:
            sel = np.random.default_rng(k).permutation(a.size)[:6000]; a, b = a[sel], b[sel]
        st = R["stats"][k]
        axs.scatter(a, b, 4, color=_col(k), alpha=0.25, edgecolors="none",
                    label="%s (r=%.2f, %+.0f%%)" % (R["channels"][k]["label"], st["r"], st["relbias_pct"]))
        allv.append(a); allv.append(b)
    if allv:
        allv = np.concatenate(allv)
        lim = (max(np.nanmin(allv), 1e-2), np.nanpercentile(allv, 99.8))
        axs.plot(lim, lim, "k--", lw=1.0); axs.set_xlim(*lim); axs.set_ylim(*lim)
    axs.set_xscale("log"); axs.set_yscale("log"); axs.grid(alpha=0.3)
    axs.set_xlabel("%s [Mm$^{-1}$ sr$^{-1}$]" % R["channels"][iref]["label"]); axs.set_ylabel("channel [Mm$^{-1}$ sr$^{-1}$]")
    axs.set_title("(b) Scatter vs %s (%.0f-%.0f m)" % (R["channels"][iref]["label"], zmin, zmax), fontsize=10)
    axs.legend(loc="lower right", fontsize=7)

    # (c) histogram of (channel - reference) differences
    axh = fig.add_subplot(gs[0, 2])
    refb = R["beta"][iref][:, band]
    dall = []
    for k in range(nch):
        if k == iref:
            continue
        d = (R["beta"][k][:, band] - refb).ravel(); dall.append(d[np.isfinite(d)])
    if dall:
        xmax = np.nanpercentile(np.abs(np.concatenate(dall)), 99)
        edges = np.linspace(-xmax, xmax, 61)
        for k in range(nch):
            if k == iref:
                continue
            d = (R["beta"][k][:, band] - refb).ravel(); d = d[np.isfinite(d)]
            axh.hist(d, edges, density=True, histtype="step", color=_col(k), lw=1.6,
                     label="%s (med %+.2f)" % (R["channels"][k]["label"], np.median(d)))
        axh.axvline(0, ls="--", color="k", lw=0.8); axh.set_xlim(-xmax, xmax)
    axh.grid(alpha=0.3); axh.set_xlabel(r"$\beta_{att}$ difference [Mm$^{-1}$ sr$^{-1}$]"); axh.set_ylabel("pdf")
    axh.set_title("(c) Difference vs %s" % R["channels"][iref]["label"], fontsize=10)
    axh.legend(loc="upper right", fontsize=7)

    # (d-g) four channel pcolors (lower-right 2x2)
    pc_pos = [(1, 1), (1, 2), (2, 1), (2, 2)]
    letters = "defg"
    last = None
    for k in range(min(nch, 4)):
        r_, c_ = pc_pos[k]
        ax = fig.add_subplot(gs[r_, c_])
        B = R["beta_disp"][k][:, zmask].T.copy()
        B[B < 1e-3] = 1e-3
        pcm = ax.pcolormesh(tx, zc, np.log10(np.abs(B)), shading="auto", vmin=CLIM[0], vmax=CLIM[1], cmap="viridis")
        cbh = np.asarray(R["cbh"][k]) - R["station"]["altitude"] if R.get("cbh") else None
        if cbh is not None and np.isfinite(cbh).any():
            ax.plot(tx, np.where((cbh > 0) & (cbh < zmax_plot), cbh, np.nan), ".", color="k", ms=2)
        ax.set_ylim(0, zmax_plot); ax.set_title("(%s) %s" % (letters[k], R["channels"][k]["label"]), fontsize=10)
        ax.set_ylabel("Alt a.g.l. [m]")
        if r_ < 2:
            ax.set_xticklabels([])
        else:
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); ax.set_xlabel("Date")
        last = pcm
    if last is not None:
        cax = fig.add_axes([0.945, 0.07, 0.012, 0.55])
        cb = fig.colorbar(last, cax=cax); cb.set_label(r"log$_{10}\beta_{att}$")

    fig.suptitle(title, fontweight="bold", fontsize=12)
    fig.savefig(out_png, dpi=200); plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  EARLINET 2x2 figure
# ---------------------------------------------------------------------------
def fig_earlinet(code, label, betaE, betaC, grid, times, stats, matlab, out_png, zmin=500, zmax=5000, zmax_plot=6000):
    zmask = (grid >= 0) & (grid <= zmax_plot)
    z = grid[zmask]
    order = np.argsort(times)
    bE = betaE[order][:, zmask]; bC = betaC[order][:, zmask]; ts = np.asarray(times)[order]
    colE, colC = COLORS[0], COLORS[1]

    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 2, figure=fig, hspace=0.22, wspace=0.2, left=0.07, right=0.95, top=0.91, bottom=0.08)

    # (a) median matched profile +/- IQR
    ax1 = fig.add_subplot(gs[0, 0])
    for B, c, nm in ((bE, colE, "EARLINET (%s)" % label), (bC, colC, "CHM15k (Rayleigh)")):
        med, q1, q3, _ = _median_iqr(B)
        good = np.isfinite(med)
        q1c = np.where(q1 > 0, q1, np.nan); q3c = np.where(q3 > 0, q3, np.nan)
        v = good & np.isfinite(q1c) & np.isfinite(q3c)
        if v.any():
            ax1.fill_betweenx(z[v], q1c[v], q3c[v], color=c, alpha=0.15)
        ax1.plot(med[good], z[good], "-", color=c, lw=2.0, label=nm)
    ax1.set_xscale("log"); ax1.set_xlim(*BETALIM); ax1.set_ylim(0, zmax_plot); ax1.grid(alpha=0.3)
    ax1.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]"); ax1.set_ylabel("Altitude a.g.l. [m]")
    ax1.set_title(r"(a) Median matched profile ($\pm$ IQR)", fontsize=10); ax1.legend(loc="upper right", fontsize=8)

    # (b) density scatter (CHM y vs EARLINET x), log-log
    ax2 = fig.add_subplot(gs[0, 1])
    band = (z >= zmin) & (z <= zmax)
    a = bC[:, band].ravel(); b = bE[:, band].ravel()
    m = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    a, b = a[m], b[m]
    if a.size:
        lim = (np.log10(max(np.nanmin(np.r_[a, b]), 1e-2)), np.log10(np.nanpercentile(np.r_[a, b], 99.8)))
        hb = ax2.hexbin(np.log10(b), np.log10(a), gridsize=55, bins="log", cmap="Blues", mincnt=1,
                        extent=(lim[0], lim[1], lim[0], lim[1]))
        ax2.plot(lim, lim, "k--", lw=1.2); ax2.set_xlim(*lim); ax2.set_ylim(*lim)
        ticks = np.arange(np.ceil(lim[0]), np.floor(lim[1]) + 1)
        ax2.set_xticks(ticks); ax2.set_yticks(ticks)
        ax2.set_xticklabels([r"10$^{%d}$" % t for t in ticks]); ax2.set_yticklabels([r"10$^{%d}$" % t for t in ticks])
        cb = fig.colorbar(hb, ax=ax2); cb.set_label("counts")
    ax2.grid(alpha=0.3); ax2.set_xlabel("EARLINET (%s) [Mm$^{-1}$ sr$^{-1}$]" % label); ax2.set_ylabel("CHM15k (Rayleigh) [Mm$^{-1}$ sr$^{-1}$]")
    ax2.set_title("(b) Density (%.0f-%.0f m): r=%.2f, bias=%+.0f%%, N=%d"
                  % (zmin, zmax, stats["r"], stats["relbias_pct"], stats["n"]), fontsize=10)

    # (c) EARLINET curtain ; (d) CHM curtain (profile index x-axis, date ticks)
    npr = bE.shape[0]
    tlbl = [np.datetime64(t, "D").astype(datetime).strftime("%y-%m-%d") for t in ts]
    ti = np.round(np.linspace(0, npr - 1, min(5, npr))).astype(int)
    for tile, (B, nm) in ((gs[1, 0], (bE, "EARLINET (%s)" % label)), (gs[1, 1], (bC, "CHM15k (Rayleigh)"))):
        ax = fig.add_subplot(tile)
        Bp = B.T.copy(); Bp[Bp < 1e-3] = 1e-3
        pcm = ax.pcolormesh(np.arange(npr), z, np.log10(np.abs(Bp)), shading="auto", vmin=CLIM[0], vmax=CLIM[1], cmap="viridis")
        ax.set_xlim(0, max(npr - 1, 1)); ax.set_xticks(ti); ax.set_xticklabels([tlbl[i] for i in ti])
        ax.set_ylabel("Alt. a.g.l. [m]"); ax.set_xlabel("Matched profile (by date)")
        ax.set_title("(%s) %s" % ("c" if tile == gs[1, 0] else "d", nm), fontsize=10)
        cb = fig.colorbar(pcm, ax=ax); cb.set_label(r"log$_{10}\beta_{att}$")

    fig.suptitle("%s — EARLINET vs CHM15k (Rayleigh)  (%d matched)" % (label, npr),
                 fontweight="bold", fontsize=12)
    fig.savefig(out_png, dpi=200); plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  Combined calibration time-series grid (all channels)
# ---------------------------------------------------------------------------
def fig_calib_timeseries(channels, calib_dir, out_png, ncol=3):
    """channels: list of dict(key, title, unit). Reads <calib_dir>/<key>.csv."""
    have = [c for c in channels if (Path(calib_dir) / f"{c['key']}.csv").is_file()]
    n = len(have); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 2.5 * nrow), squeeze=False)
    for i, c in enumerate(have):
        ax = axes[i // ncol][i % ncol]
        t, cd, ck, cks = [], [], [], []
        with open(Path(calib_dir) / f"{c['key']}.csv", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    t.append(datetime.strptime(row["time"][:10], "%Y-%m-%d"))
                except Exception:
                    continue
                cd.append(_f(row.get("C_daily"))); ck.append(_f(row.get("C_kalman"))); cks.append(_f(row.get("C_kalman_std")))
        t = np.array(t); cd = np.array(cd); ck = np.array(ck); cks = np.array(cks)
        ax.plot(t, cd, "x", color="0.45", ms=4, mew=0.8, label="raw daily")
        good = np.isfinite(ck)
        if good.any():
            ax.plot(t[good], ck[good], "-", color=COLORS[1], lw=1.6, label="Kalman")
            sg = good & np.isfinite(cks)
            if sg.any():
                ax.fill_between(t[sg], (ck - cks)[sg], (ck + cks)[sg], color=COLORS[1], alpha=0.2)
        # focus the y-axis on the actual coefficient values (the Kalman ±1σ band on sparse
        # channels can be much larger than the spread and would otherwise flatten the panel).
        vals = np.concatenate([cd[np.isfinite(cd)], ck[good]])
        if vals.size:
            lo, hi = np.nanpercentile(vals, 2), np.nanpercentile(vals, 98)
            pad = 0.15 * (hi - lo) if hi > lo else 0.1 * abs(hi) + 1e-12
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(c["title"], fontsize=9); ax.grid(alpha=0.3)
        ax.set_ylabel(c.get("unit", "C$_L$ [a.u.]"), fontsize=8)
        ax.tick_params(labelsize=7)
        for lab in ax.get_xticklabels():
            lab.set_rotation(0); lab.set_fontsize(7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper left")
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Calibration coefficient time series — raw daily (x) and Kalman estimate (line, $\\pm1\\sigma$)",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out_png, dpi=150); plt.close(fig)
    return out_png


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _fmt_my(yyyymmdd):
    return datetime.strptime(yyyymmdd, "%Y%m%d").strftime("%b %Y")
