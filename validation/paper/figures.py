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

# General colour rules (2026-07-02): CHM15k / Rayleigh-calibrated reference = red; CL61 = blue
# (Rayleigh) / dark grey (cloud); Mini-MPL = green; CL31 = orange; CL51 = purple; EARLINET = black.
# Exception: several units of the SAME type at one station (Amsterdam 4x CHM15k) -> default palette.
TYPE_COLORS = {"CHM15k": "#d62728", "CHM8k": "#d62728", "CL31": "#ff7f0e", "CL51": "#9467bd",
               "Mini-MPL": "#2ca02c", "MPL": "#2ca02c", "EARLINET": "#000000"}
CL61_COLORS = {"rayleigh": "#1f77b4", "cloud": "#404040"}


def _col(k):
    return COLORS[k % len(COLORS)]


def channel_colors(channels):
    """Per-channel colours by instrument type + calibration method (see rules above).
    channels: list of dicts with itype + calib (falls back to the palette when a type appears
    more than once at the station, or the type is unknown)."""
    from collections import Counter
    cnt = Counter(c.get("itype", "") for c in channels)
    cols = []
    seen_cl61 = {}
    for k, c in enumerate(channels):
        it = c.get("itype", "")
        if it == "CL61":
            cal = c.get("calib", "")
            if seen_cl61.get(cal):          # 2nd CL61 entry of the same method (e.g. offset-corr)
                cols.append("#17BECF")       # -> cyan to stay distinguishable
            else:
                cols.append(CL61_COLORS.get(cal, "#404040"))
            seen_cl61[cal] = True
        elif it in TYPE_COLORS and cnt[it] == 1:
            cols.append(TYPE_COLORS[it])
        else:
            cols.append(_col(k))
    return cols


def _median_iqr(B):
    """median, 25th, 75th percentile over axis 0 (profiles), NaN-aware."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(B, axis=0)
        q1 = np.nanpercentile(B, 25, axis=0)
        q3 = np.nanpercentile(B, 75, axis=0)
        nz = np.sum(np.isfinite(B), axis=0)
    return med, q1, q3, nz


def _truncate_noise_floor(med, nz, z, zmin, nprof, minfrac=0.20):
    """Keep up to the first altitude above zMin where the median stops being defined or the
    coverage gets too thin. On the LINEAR profile axis the median is shown even where it goes
    NEGATIVE (a channel whose median dips below zero has hit its noise/offset floor — that is the
    information, not something to hide); only the >20 %-coverage safeguard truncates the curve, so a
    handful of profiles can never define it."""
    good = np.isfinite(med) & (nz >= max(10, minfrac * nprof))
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
    cols = channel_colors(R["channels"])
    tx = mdates.date2num(np.asarray(R["time_sync"]).astype("datetime64[s]").astype(datetime))

    fig = plt.figure(figsize=(19, 10))
    gs = GridSpec(3, 3, figure=fig, hspace=0.28, wspace=0.26,
                  left=0.05, right=0.93, top=0.92, bottom=0.07)

    # (a) profile median +/- IQR — left column.
    # Common-time restriction: only hours where EVERY channel has a valid screened profile enter
    # the medians, so the per-channel medians describe the same atmospheric sample (a sparse
    # channel would otherwise be compared against medians of different weather).
    axp = fig.add_subplot(gs[:, 0])
    have = [np.any(np.isfinite(R["beta"][k][:, zmask]), axis=1) for k in range(nch)]
    common = np.logical_and.reduce(have)
    ncom = int(common.sum())
    synced = ncom > 0
    xmax_seen = 0.0
    xmin_seen = 0.0
    for k in range(nch):
        B = R["beta"][k][:, zmask]
        if synced:
            B = B[common]
        med, q1, q3, nz = _median_iqr(B)
        keep = _truncate_noise_floor(med, nz, zc, zmin, B.shape[0])
        c = cols[k]
        m1 = keep & np.isfinite(q1) & (q1 > 0)
        axp.plot(q1[m1], zc[m1], ":", color=c, lw=1.0)
        m3 = keep & np.isfinite(q3) & (q3 > 0)
        axp.plot(q3[m3], zc[m3], ":", color=c, lw=1.0)
        axp.plot(med[keep], zc[keep], "-", color=c, lw=1.8, label=R["channels"][k]["label"])
        if keep.any() and np.isfinite(med[keep]).any():
            xmin_seen = min(xmin_seen, float(np.nanmin(med[keep])))   # let the axis show negative medians
        # overlap-UNcorrected twin (CHM15k): dashed median in the same colour, same common-hour
        # sample and truncation, so the effect of the temperature-dependent overlap correction
        # (below ~700 m) is directly visible against the corrected solid line.
        Bn = R.get("beta_noovl", [None] * nch)[k]
        if Bn is not None:
            Bn = Bn[:, zmask]
            if synced:
                Bn = Bn[common]
            medn, _, _, nzn = _median_iqr(Bn)
            keepn = _truncate_noise_floor(medn, nzn, zc, zmin, Bn.shape[0])
            axp.plot(medn[keepn], zc[keepn], "--", color=c, lw=1.4,
                     label=R["channels"][k]["label"] + " (no overlap corr)")
        if m3.any():
            xmax_seen = max(xmax_seen, float(np.nanpercentile(q3[m3], 98)))
    axp.axhline(zmin, ls="--", color="k", lw=0.8); axp.axhline(zmax, ls="--", color="k", lw=0.8)
    # clamp the profile axis to +/-2 Mm^-1 sr^-1: keep the median's zero-crossing / negative dip visible
    # without letting an extreme excursion (e.g. the CL31 diving to ~-4) stretch the whole panel.
    axp.set_xlim(max(-2.0, xmin_seen * 1.05), min(2.0, xmax_seen * 1.05) if xmax_seen > 0 else 1.0)
    if xmin_seen < 0:
        axp.axvline(0, color="0.6", lw=0.8, ls="-")     # mark beta_att = 0 when negative medians are shown
    axp.set_ylim(0, zmax_plot)
    axp.grid(alpha=0.3); axp.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]"); axp.set_ylabel("Altitude a.g.l. [m]")
    axp.set_title("(a) Median (solid) $\\pm$ IQR (dotted)  %s-%s\n%s"
                  % (_fmt_my(cfg["start"]), _fmt_my(cfg["end"]),
                     "N=%d common hours (all channels)" % ncom if synced
                     else "no common hours - per-channel sampling"), fontsize=10)
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
        axs.scatter(a, b, 4, color=cols[k], alpha=0.25, edgecolors="none",
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
            axh.hist(d, edges, density=True, histtype="step", color=cols[k], lw=1.6,
                     label="%s (med %+.2f)" % (R["channels"][k]["label"], np.median(d)))
        axh.axvline(0, ls="--", color="k", lw=0.8); axh.set_xlim(-xmax, xmax)
    axh.grid(alpha=0.3); axh.set_xlabel(r"$\beta_{att}$ difference [Mm$^{-1}$ sr$^{-1}$]"); axh.set_ylabel("pdf")
    axh.set_title("(c) Difference vs %s" % R["channels"][iref]["label"], fontsize=10)
    axh.legend(loc="upper right", fontsize=7)

    # (d-g) four channel pcolors (lower-right 2x2) — OmB style: ALL data in greyscale, only the
    # KEPT gates (screening + SNR + coverage; what the medians/statistics use) in colour on top,
    # so everything flagged stays grey.
    pc_pos = [(1, 1), (1, 2), (2, 1), (2, 2)]
    letters = "defg"
    last = None
    for k in range(min(nch, 4)):
        r_, c_ = pc_pos[k]
        ax = fig.add_subplot(gs[r_, c_])
        B = R["beta_disp"][k][:, zmask].T.copy()
        B[B < 1e-3] = 1e-3
        ax.pcolormesh(tx, zc, np.log10(np.abs(B)), shading="auto", vmin=CLIM[0], vmax=CLIM[1],
                      cmap="gray_r")
        K = R["beta"][k][:, zmask].T.copy()
        K[K < 1e-3] = 1e-3
        pcm = ax.pcolormesh(tx, zc, np.log10(np.abs(K)), shading="auto", vmin=CLIM[0], vmax=CLIM[1],
                            cmap="viridis")
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
        fig.text(0.951, 0.655, "colour = kept\ngrey = flagged", fontsize=7, ha="left")

    fig.suptitle(title, fontweight="bold", fontsize=12)
    fig.savefig(out_png, dpi=200); plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  EARLINET 2x2 figure
# ---------------------------------------------------------------------------
def fig_earlinet(code, label, betaE, betaC, grid, times, stats, matlab, out_png, zmin=500, zmax=5000,
                 zmax_plot=6000, betaC_raw=None, instr="CHM15k (Rayleigh)", itype="CHM15k"):
    zmask = (grid >= 0) & (grid <= zmax_plot)
    z = grid[zmask]
    order = np.argsort(times)
    bE = betaE[order][:, zmask]; bC = betaC[order][:, zmask]; ts = np.asarray(times)[order]
    bCraw = betaC_raw[order][:, zmask] if betaC_raw is not None else None
    colE = TYPE_COLORS["EARLINET"]                                 # black
    colC = TYPE_COLORS.get(itype, TYPE_COLORS["CHM15k"])           # per-type (CHM red, MPL green)

    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 2, figure=fig, hspace=0.22, wspace=0.2, left=0.07, right=0.95, top=0.91, bottom=0.08)

    # (a) median matched profile +/- IQR. A gate's median is only shown when MORE THAN 20 % of
    # the matched profiles contribute — otherwise a single profile could define the curve.
    ax1 = fig.add_subplot(gs[0, 0])
    nprof = bE.shape[0]
    minn = max(10, int(round(0.20 * nprof)))
    xmax_seen = 0.0
    for B, c, nm in ((bE, colE, "EARLINET (%s)" % label), (bC, colC, instr)):
        med, q1, q3, nz = _median_iqr(B)
        good = np.isfinite(med) & (nz >= minn)
        q1c = np.where(q1 > 0, q1, np.nan); q3c = np.where(q3 > 0, q3, np.nan)
        v = good & np.isfinite(q1c) & np.isfinite(q3c)
        if v.any():
            ax1.fill_betweenx(z[v], q1c[v], q3c[v], color=c, alpha=0.15)
            xmax_seen = max(xmax_seen, float(np.nanpercentile(q3c[v], 98)))
        ax1.plot(med[good], z[good], "-", color=c, lw=2.0, label=nm)
    ax1.set_xlim(0, xmax_seen * 1.05 if xmax_seen > 0 else 1.0); ax1.set_ylim(0, zmax_plot); ax1.grid(alpha=0.3)
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
    ax2.grid(alpha=0.3); ax2.set_xlabel("EARLINET (%s) [Mm$^{-1}$ sr$^{-1}$]" % label); ax2.set_ylabel("%s [Mm$^{-1}$ sr$^{-1}$]" % instr)
    ax2.set_title("(b) Density (%.0f-%.0f m): r=%.2f (log r=%.2f), bias=%+.0f%% (med %+.0f%%), N=%d"
                  % (zmin, zmax, stats["r"], stats.get("r_log", np.nan),
                     stats["relbias_pct"], stats.get("medrelbias_pct", np.nan), stats["n"]), fontsize=10)

    # (c) EARLINET curtain ; (d) CHM curtain (profile index x-axis, date ticks). The CHM curtain
    # is OmB-style: the unscreened data in greyscale, only the KEPT gates in colour on top —
    # everything flagged (clouds/fog/qf/SNR) stays grey.
    npr = bE.shape[0]
    tlbl = [np.datetime64(t, "D").astype(datetime).strftime("%y-%m-%d") for t in ts]
    ti = np.round(np.linspace(0, npr - 1, min(5, npr))).astype(int)
    for tile, (B, Braw, nm) in ((gs[1, 0], (bE, None, "EARLINET (%s)" % label)),
                                (gs[1, 1], (bC, bCraw, instr))):
        ax = fig.add_subplot(tile)
        if Braw is not None:
            Bg = Braw.T.copy(); Bg[Bg < 1e-3] = 1e-3
            ax.pcolormesh(np.arange(npr), z, np.log10(np.abs(Bg)), shading="auto",
                          vmin=CLIM[0], vmax=CLIM[1], cmap="gray_r")
        Bp = B.T.copy(); Bp[Bp < 1e-3] = 1e-3
        pcm = ax.pcolormesh(np.arange(npr), z, np.log10(np.abs(Bp)), shading="auto", vmin=CLIM[0], vmax=CLIM[1], cmap="viridis")
        ax.set_xlim(0, max(npr - 1, 1)); ax.set_xticks(ti); ax.set_xticklabels([tlbl[i] for i in ti])
        ax.set_ylabel("Alt. a.g.l. [m]"); ax.set_xlabel("Matched profile (by date)")
        ax.set_title("(%s) %s%s" % ("c" if tile == gs[1, 0] else "d", nm,
                                    "" if Braw is None else "  (grey = flagged)"), fontsize=10)
        cb = fig.colorbar(pcm, ax=ax); cb.set_label(r"log$_{10}\beta_{att}$")

    fig.suptitle("%s — EARLINET vs %s  (%d matched)" % (label, instr, npr),
                 fontweight="bold", fontsize=12)
    fig.savefig(out_png, dpi=200); plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
#  Combined calibration time-series grid (all channels)
# ---------------------------------------------------------------------------
def fig_calib_timeseries(channels, calib_dir, out_png, ncol=5):
    """One panel per INSTRUMENT showing the absolute lidar constant C_L: raw daily (x) + Kalman
    (line +/- 1 sigma) for every calibration method available (CL61: Rayleigh AND cloud in the
    same panel). Colours follow the general per-type rules (CHM15k red, CL31 orange, CL51 purple,
    Mini-MPL green; CL61 blue=Rayleigh / dark grey=cloud). channels: list of
    dict(title, unit, series=[dict(key, calib, itype)]). Landscape."""
    have = []
    for c in channels:
        ser = [s for s in c["series"] if (Path(calib_dir) / f"{s['key']}.csv").is_file()]
        if ser:
            have.append(dict(c, series=ser))
    n = len(have); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 2.6 * nrow), squeeze=False)
    for i, c in enumerate(have):
        ax = axes[i // ncol][i % ncol]
        vals = []
        for s in c["series"]:
            col = (CL61_COLORS.get(s.get("calib", ""), "#404040") if s.get("itype") == "CL61"
                   else TYPE_COLORS.get(s.get("itype", ""), "0.45"))
            t, cd, ck, cks = _read_calib_csv(Path(calib_dir) / f"{s['key']}.csv")
            ax.plot(t, cd, "x", color=col, ms=4, mew=0.8, alpha=0.65)
            good = np.isfinite(ck)
            if good.any():
                ax.plot(t[good], ck[good], "-", color=col, lw=1.7, label=s["calib"])
                sg = good & np.isfinite(cks)
                if sg.any():
                    ax.fill_between(t[sg], (ck - cks)[sg], (ck + cks)[sg], color=col, alpha=0.18)
            vals.append(cd[np.isfinite(cd)]); vals.append(ck[good])
        # focus the y-axis on the actual constant values (the Kalman ±1σ band on sparse channels
        # can be much larger than the spread and would otherwise flatten the panel).
        vals = np.concatenate(vals) if vals else np.array([])
        if vals.size:
            lo, hi = np.nanpercentile(vals, 2), np.nanpercentile(vals, 98)
            pad = 0.15 * (hi - lo) if hi > lo else 0.1 * abs(hi) + 1e-12
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(c["title"], fontsize=9); ax.grid(alpha=0.3)
        ax.set_ylabel(c.get("unit", "C$_L$"), fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc="upper left")
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Lidar constant C$_L$ time series — raw daily (x) and Kalman estimate (line, $\\pm1\\sigma$); "
                 "CHM15k red, CL31 orange, CL51 purple, Mini-MPL green, CL61 blue=Rayleigh / dark grey=cloud",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out_png, dpi=150); plt.close(fig)
    return out_png


def _read_calib_csv(path):
    t, cd, ck, cks = [], [], [], []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                t.append(datetime.strptime(row["time"][:10], "%Y-%m-%d"))
            except Exception:
                continue
            cd.append(_f(row.get("C_daily"))); ck.append(_f(row.get("C_kalman"))); cks.append(_f(row.get("C_kalman_std")))
    return np.array(t), np.array(cd), np.array(ck), np.array(cks)


def fig_calib_l1l2(channels, calib_dir, out_png, ncol=3):
    """One grid of all channels overlaying the L1 (blue) and L2 (orange) calibration series:
    raw daily (markers) + Kalman (line). Shows that L1 (binned) and L2 give the same coefficient."""
    cd_dir = Path(calib_dir)
    have = [c for c in channels if (cd_dir / f"{c['key']}_L1.csv").is_file() or (cd_dir / f"{c['key']}_L2.csv").is_file()]
    n = len(have); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 2.5 * nrow), squeeze=False)
    # Each series is normalised by its own Kalman median, so L1 and L2 share a y-scale and the
    # comparison is of the *temporal pattern* (the absolute L1/L2 ratio is reported in the table:
    # Rayleigh lidar constants match to ~1 %, the cloud coefficient is on a different input scale).
    styles = [("L2", COLORS[1], "x"), ("L1", COLORS[0], "+")]
    for i, c in enumerate(have):
        ax = axes[i // ncol][i % ncol]
        for level, color, marker in styles:
            f = cd_dir / f"{c['key']}_{level}.csv"
            if not f.is_file():
                continue
            t, cd, ck, cks = _read_calib_csv(f)
            if t.size == 0:
                continue
            med = np.nanmedian(ck[np.isfinite(ck)]) if np.isfinite(ck).any() else np.nan
            if not np.isfinite(med) or med == 0:
                continue
            ax.plot(t, cd / med, marker, color=color, ms=4, mew=0.9, alpha=0.55, label=f"{level} daily")
            g = np.isfinite(ck)
            if g.any():
                ax.plot(t[g], ck[g] / med, "-", color=color, lw=1.5, label=f"{level} Kalman")
        ax.axhline(1.0, color="0.6", lw=0.6, ls=":")
        ax.set_ylim(0.4, 1.6)
        ax.set_title(c["title"], fontsize=9); ax.grid(alpha=0.3)
        ax.set_ylabel("C / median", fontsize=8); ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=6.5, loc="best", ncol=2)
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Calibration coefficient L1 (binned to L2 grid) vs L2 — normalised to each median (eprof_v2)",
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
