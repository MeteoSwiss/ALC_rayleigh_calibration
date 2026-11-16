# -*- coding: utf-8 -*-
"""q20 stage 3 - figures of the per-unit correction budget (report 19, q20). Landscape PNGs, English.

Run after q20_corrections.py:
    .venv/Scripts/python.exe validation/chm15k_stability/q20_figures.py
Outputs in C:/DATA/Projects/202610_CHM15k_stability/q20/:
    q20_unit_<label>.png        one panel per unit: raw vs fully corrected nightly C_L + paper-Kalman lines,
                                amplitude / scatter after every cumulative step
    q20_units_grid_1.png, _2.png  all 20 units, normalised Kalman lines (overview)
    q20_waterfall.png           network amplitude and scatter waterfall + monthly composite of what remains
    q20_unit_steps_heatmap.png  amplitude and scatter per unit and step
    q20_residual_map.png        map and bar chart of the residual amplitude
    q20_twins.png               co-located twin ratios and the seasonal-kappa fraction
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.io as sio  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q20_lib import ALL_UNITS, PAYERNE, Q20, STEPS, UNITS, to_datetime  # noqa: E402

warnings.filterwarnings("ignore")
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
COAST = Path("C:/Program Files/MATLAB/R2025b/toolbox/map/mapdata/coastlines.mat")
C_RAW, C_SYN, C_FULL, C_TWIN = "0.35", "#1f77b4", "#d95f02", "#7570b3"
STEP_SHORT = {"raw": "operational", "kappa": "1 κ synoptic", "aerosol": "2 aerosol", "window": "3 window",
              "module": "4 module", "pcal": "5 p_calc", "molecular": "6 molecular", "tier2": "7 tier-2",
              "other": "8 sampling"}


def load():
    s = pd.read_csv(Q20 / "q20_series.csv", dtype={"date": str, "module": str})
    st = pd.read_csv(Q20 / "q20_unit_steps.csv")
    terms = pd.read_csv(Q20 / "q20_terms.csv")
    units = pd.read_csv(Q20 / "q20_units.csv")
    net = pd.read_csv(Q20 / "q20_network.csv")
    tw = pd.read_csv(Q20 / "q20_twins.csv")
    kf = pd.read_csv(Q20 / "q20_kappa_fraction.csv")
    kal = dict(np.load(Q20 / "q20_kalman.npz"))
    return s, st, terms, units, net, tw, kf, kal


def metric(st, key, version, col):
    x = st[(st.stream == key) & (st.version == version)][col]
    return float(x.iloc[0]) if len(x) else np.nan


def kal_line(kal, key, version):
    k = f"{key}__{version}"
    if f"{k}__grid" not in kal:
        return None, None
    return to_datetime(kal[f"{k}__grid"]), kal[f"{k}__state"]


def fname(label):
    return label.replace(" ", "_")


# ---------------------------------------------------------------------------------------------
# Per-unit panel
# ---------------------------------------------------------------------------------------------
def unit_figure(s, st, units, kal, key):
    g = s[s.stream == key].sort_values("tmid")
    u = units[units.stream == key].iloc[0]
    label = u.label
    fig = plt.figure(figsize=(16, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.7, 1.0], wspace=0.25)
    ax = fig.add_subplot(gs[0])
    tt = to_datetime(g.tmid.values)
    t2 = g.is_tier2.values
    ax.scatter(tt, g.c / 1e11, s=9, color=C_RAW, alpha=0.45, lw=0, label="operational nightly C_L")
    ax.scatter(np.array(tt)[t2], g.c[t2] / 1e11, s=22, facecolor="none", edgecolor="tab:red", lw=0.8,
               label="tier-2 (fallback) nights")
    corr = np.exp(g.lnc_other.values) / 1e11
    ax.scatter(np.array(tt)[~t2], corr[~t2], s=9, color=C_SYN, alpha=0.55, lw=0,
               label="fully corrected nightly C_L (κ synoptic)")
    for v, col, ls, lab in (("raw", "k", "-", "Kalman, operational"),
                            ("cum_other", C_SYN, "-", "Kalman, fully corrected (κ synoptic)"),
                            ("final_kfull", C_FULL, "--", "Kalman, fully corrected (κ on full T_int)")):
        gx, gy = kal_line(kal, key, v)
        if gx is not None:
            ax.plot(gx, np.asarray(gy) / 1e11, color=col, ls=ls, lw=1.6, label=lab)
    lo, hi = np.nanpercentile(np.r_[g.c / 1e11, corr], [1, 99])
    ax.set_ylim(lo * 0.85, hi * 1.25)
    mods = g.module.values
    for i in np.where(mods[1:] != mods[:-1])[0]:
        ax.axvline(tt[i + 1], color="tab:purple", ls=":", lw=1.2)
        ax.text(tt[i + 1], lo * 0.88, f" module TUB{mods[i]} → TUB{mods[i + 1]}", color="tab:purple",
                fontsize=8, va="bottom")
    ax.set_ylabel("C_L [10$^{11}$]")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.5, ncol=2, framealpha=0.85)
    a0, p0, n0 = (metric(st, key, "raw", c) for c in ("amp", "peak", "n2n"))
    a1, p1, n1 = (metric(st, key, "cum_other", c) for c in ("amp", "peak", "n2n"))
    a2, p2 = (metric(st, key, "final_kfull", c) for c in ("amp", "peak"))
    note = ""
    if key == PAYERNE[0]:
        a3, p3 = (metric(st, key, "final_t2scale", c) for c in ("amp", "peak"))
        note = (f"\nABNORMAL UNIT: κ unreliable (noisy receiver); {int(t2.sum())} of {len(g)} nights are tier-2; "
                f"with tier-2 nights scaled (-{9.4:.1f} %) instead of dropped: ±{a3:.1f} % (DOY {p3:.0f})")
    ax.set_title(f"{label} ({u.country}, {key}) - module(s) TUB{u.modules.replace('/', '/TUB')}, "
                 f"κ = {u.kappa:+.2f} ± {u.kappa_se:.2f} %/K (intra-night, q2a)\n"
                 f"annual amplitude ±{a0:.1f} % (peak DOY {p0:.0f}) → ±{a1:.1f} % (DOY {p1:.0f}) "
                 f"[full κ: ±{a2:.1f} %, DOY {p2:.0f}];  night-to-night {n0:.1f} → {n1:.1f} %" + note, fontsize=9.5)
    # right: cumulative steps
    axb = fig.add_subplot(gs[1])
    vers = ["raw"] + [f"cum_{x}" for x in STEPS[1:]] + ["final_kfull"]
    amps = [metric(st, key, v, "amp") for v in vers]
    n2ns = [metric(st, key, v, "n2n") for v in vers]
    x = np.arange(len(vers))
    cols = [C_RAW] + [C_SYN] * (len(vers) - 2) + [C_FULL]
    axb.bar(x, amps, color=cols, alpha=0.85)
    axb.set_xticks(x)
    axb.set_xticklabels([STEP_SHORT.get(v.replace("cum_", ""), v) if v != "final_kfull" else "all, κ full T"
                         for v in vers], rotation=60, ha="right", fontsize=8)
    axb.set_ylabel("annual amplitude of ln C_L [%]")
    axb.grid(axis="y", alpha=0.3)
    ax2 = axb.twinx()
    ax2.plot(x, n2ns, "kD", ms=5, label="night-to-night scatter")
    ax2.set_ylabel("night-to-night scatter [%]")
    ax2.set_ylim(0, max(np.nanmax(n2ns) * 1.3, 1))
    axb.set_ylim(0, max(np.nanmax(amps) * 1.25, 1))
    ax2.legend(loc="upper right", fontsize=7.5)
    axb.set_title("after each cumulative correction", fontsize=9.5)
    fig.savefig(Q20 / f"q20_unit_{fname(label)}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Overview grids
# ---------------------------------------------------------------------------------------------
def grid_figures(s, st, units, kal):
    keys = [k for k, *_ in UNITS]
    for part, chunk in enumerate((keys[:10], keys[10:]), start=1):
        fig, axs = plt.subplots(2, 5, figsize=(22, 8.4), sharey=True)
        for ax, key in zip(axs.ravel(), chunk):
            g = s[s.stream == key].sort_values("tmid")
            u = units[units.stream == key].iloc[0]
            ref = np.median(g.c)
            t2 = g.is_tier2.values
            ax.scatter(np.array(to_datetime(g.tmid.values))[~t2], np.exp(g.lnc_other.values[~t2]) / ref, s=4,
                       color=C_SYN, alpha=0.35, lw=0)
            for v, col, ls in (("raw", "k", "-"), ("cum_other", C_SYN, "-"), ("final_kfull", C_FULL, "--")):
                gx, gy = kal_line(kal, key, v)
                if gx is not None:
                    ax.plot(gx, np.asarray(gy) / ref, color=col, ls=ls, lw=1.5)
            a0, p0 = metric(st, key, "raw", "amp"), metric(st, key, "raw", "peak")
            a1, p1 = metric(st, key, "cum_other", "amp"), metric(st, key, "cum_other", "peak")
            a2 = metric(st, key, "final_kfull", "amp")
            ax.set_title(f"{u.label} ({u.country}), κ {u.kappa:+.1f} %/K", fontsize=9.5)
            ax.text(0.02, 0.03, f"±{a0:.0f} % (DOY {p0:.0f}) → ±{a1:.0f} % (DOY {p1:.0f})\nfull κ: ±{a2:.0f} %",
                    transform=ax.transAxes, fontsize=8, va="bottom",
                    bbox=dict(facecolor="white", alpha=0.75, lw=0))
            ax.set_ylim(0.5, 1.7)
            ax.grid(alpha=0.3)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%y-%m"))
            ax.tick_params(axis="x", labelsize=7.5)
        for ax in axs[:, 0]:
            ax.set_ylabel("C_L / operational median")
        h = [plt.Line2D([], [], color="k", lw=1.5), plt.Line2D([], [], color=C_SYN, lw=1.5),
             plt.Line2D([], [], color=C_FULL, lw=1.5, ls="--"),
             plt.Line2D([], [], color=C_SYN, marker="o", lw=0, alpha=0.5, ms=3)]
        fig.legend(h, ["Kalman, operational", "Kalman, all corrections (κ on synoptic T_int deviations)",
                       "Kalman, all corrections (κ on the full T_int)", "corrected nightly C_L"],
                   loc="upper center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, 1.01))
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(Q20 / f"q20_units_grid_{part}.png", dpi=110)
        plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Network waterfall + monthly composite of the residual
# ---------------------------------------------------------------------------------------------
def monthly_composite(s, col, drop_t2):
    out = []
    for key, *_ in UNITS:
        g = s[s.stream == key]
        if drop_t2:
            g = g[~g.is_tier2]
        y = g[col].values
        m = pd.to_datetime(g.tmid.values, unit="D").month
        mm = pd.Series(y).groupby(m).median()
        out.append(100 * (mm - mm.mean()).reindex(range(1, 13)))
    return pd.concat(out, axis=1).median(axis=1)


def waterfall(s, st, terms, net):
    fig, axs = plt.subplots(1, 3, figsize=(20, 6.2), gridspec_kw=dict(width_ratios=[1.35, 1.0, 1.0]))
    ax = axs[0]
    vers = ["raw"] + [f"cum_{x}" for x in STEPS[1:]]
    sub = st[st.stream.isin([k for k, *_ in UNITS])]
    med = [np.nanmedian(sub[sub.version == v].amp) for v in vers]
    q1 = [np.nanpercentile(sub[sub.version == v].amp, 25) for v in vers]
    q3 = [np.nanpercentile(sub[sub.version == v].amp, 75) for v in vers]
    x = np.arange(len(vers))
    # waterfall: change of the network median at each step
    for i in range(len(vers)):
        if i == 0:
            ax.bar(i, med[0], color=C_RAW)
        else:
            d = med[i] - med[i - 1]
            ax.bar(i, d, bottom=med[i - 1], color="tab:green" if d < 0 else "tab:red", alpha=0.8)
            ax.text(i, max(med[i], med[i - 1]) + 0.4, f"{d:+.1f}", ha="center", fontsize=8)
    ax.bar(len(vers), med[-1], color=C_SYN)
    for extra, col, lab in (("final_ktwin", C_TWIN, "residual, κ twin-f"), ("final_kfull", C_FULL, "residual, κ full T")):
        v = np.nanmedian(sub[sub.version == extra].amp)
        ax.bar(len(vers) + (1 if extra == "final_ktwin" else 2), v, color=col)
        ax.text(len(vers) + (1 if extra == "final_ktwin" else 2), v + 0.4, f"{v:.1f}", ha="center", fontsize=8)
    ax.errorbar(x, med, yerr=[np.array(med) - q1, np.array(q3) - med], fmt="none", ecolor="k", lw=0.8, capsize=3)
    ax.text(len(vers), med[-1] + 0.4, f"{med[-1]:.1f}", ha="center", fontsize=8)
    ax.text(0, med[0] + 0.4, f"{med[0]:.1f}", ha="center", fontsize=8)
    ax.set_xticks(np.arange(len(vers) + 3))
    ax.set_xticklabels([STEP_SHORT[v.replace("cum_", "")] for v in vers] + ["residual\n(κ synoptic)", "residual\n(κ twin f)",
                                                                              "residual\n(κ full T)"],
                       rotation=50, ha="right", fontsize=8.5)
    ax.set_ylabel("annual amplitude of ln C_L, network median [%]")
    ax.set_title("(a) annual amplitude: operational → each cumulative correction → residual\n"
                 "(bars = change of the median of 20 units; whiskers = IQR after the step)", fontsize=9.5)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, max(q3) * 1.15)
    # (b) one-at-a-time in-phase content and n2n
    ax = axs[1]
    tsub = terms[terms.stream.isin([k for k, *_ in UNITS])]
    order = STEPS[1:]
    ip = [np.nanmedian(tsub[tsub.term == t].inphase) for t in order]
    ip_q1 = [np.nanpercentile(tsub[tsub.term == t].inphase, 25) for t in order]
    ip_q3 = [np.nanpercentile(tsub[tsub.term == t].inphase, 75) for t in order]
    y = np.arange(len(order))
    ax.barh(y, ip, color=["tab:red" if v < 0 else "tab:green" for v in ip], alpha=0.8)
    ax.errorbar(ip, y, xerr=[np.array(ip) - ip_q1, np.array(ip_q3) - ip], fmt="none", ecolor="k", lw=0.8, capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([STEP_SHORT[t] for t in order])
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("in-phase annual content of the correction term [%]\n(> 0: removing it lowers the cycle)")
    ax.set_title("(b) what each correction removes, one at a time\n(median and IQR over the 20 units)", fontsize=9.5)
    ax.grid(axis="x", alpha=0.3)
    ax3 = ax.twiny()
    n2 = [np.nanmedian(tsub[tsub.term == t].n2n) for t in order]
    ax3.plot(n2, y, "kD", ms=5)
    ax3.set_xlabel("night-to-night content of the term [%] (◆)")
    ax3.set_xlim(0, max(n2) * 1.4)
    # (c) monthly composite of what remains
    ax = axs[2]
    months = np.arange(1, 13)
    for col, drop, lab, c, ls in (("lnc_raw", False, "operational", "k", "-"),
                                  ("lnc_other", True, "all corrections, κ synoptic", C_SYN, "-"),
                                  ("lnc_final_ktwin", True, "all corrections, κ twin f", C_TWIN, "-."),
                                  ("lnc_final_kfull", True, "all corrections, κ on full T_int", C_FULL, "--")):
        comp = monthly_composite(s, col, drop)
        ax.plot(months, comp.values, color=c, ls=ls, marker="o", ms=4, lw=1.6, label=lab)
    ax.axvline(6.7, color="goldenrod", lw=1, ls=":")
    ax.text(6.8, ax.get_ylim()[1] * 0.92, "solstice", color="goldenrod", fontsize=8, rotation=90, va="top")
    ax.set_xticks(months)
    ax.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
    ax.set_ylabel("monthly median of ln C_L minus its annual mean [%]")
    ax.set_title("(c) what remains: monthly composite, median of 20 units", fontsize=9.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower center")
    # night-to-night scatter after each cumulative step (right axis of panel a)
    n2 = [np.nanmedian(sub[sub.version == v].n2n) for v in vers]
    n2x = [np.nanmedian(sub[sub.version == v].n2n) for v in ("final_ktwin", "final_kfull")]
    axn = axs[0].twinx()
    axn.plot(np.arange(len(vers) + 3), n2 + [n2[-1]] + n2x, "kD", ms=5, label="night-to-night scatter (median)")
    for i, v in enumerate(n2 + [n2[-1]] + n2x):
        axn.text(i + 0.15, v, f"{v:.1f}", fontsize=7, va="center")
    axn.set_ylim(0, max(n2) * 2.6)
    axn.set_ylabel("night-to-night scatter, network median [%] (◆)")
    axn.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(Q20 / "q20_waterfall.png", dpi=125)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Heatmap per unit and step
# ---------------------------------------------------------------------------------------------
def heatmap(st, units):
    keys = [k for k, *_ in UNITS] + [PAYERNE[0]]
    labs = [units.set_index("stream").loc[k, "label"] for k in keys]
    vers = ["raw"] + [f"cum_{x}" for x in STEPS[1:]] + ["final_ktwin", "final_kfull"]
    vlab = [STEP_SHORT[v.replace("cum_", "")] if v.startswith("cum_") or v == "raw" else
            {"final_ktwin": "all, κ twin f", "final_kfull": "all, κ full T"}[v] for v in vers]
    fig, axs = plt.subplots(2, 1, figsize=(19, 10.5))
    for ax, col, title, cmap in ((axs[0], "amp", "annual amplitude of ln C_L [%] after each cumulative step", "viridis"),
                                 (axs[1], "n2n", "night-to-night scatter [%] after each cumulative step", "magma")):
        M = np.array([[metric(st, k, v, col) for k in keys] for v in vers])
        im = ax.imshow(M, aspect="auto", cmap=cmap)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=7,
                            color="w" if M[i, j] < np.nanpercentile(M, 60) else "k")
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels(labs, rotation=40, ha="right", fontsize=8)
        ax.set_yticks(range(len(vers)))
        ax.set_yticklabels(vlab, fontsize=8)
        ax.axvline(len(keys) - 1.5, color="w", lw=2)
        ax.set_title(title + " (Payerne, abnormal unit, at the right)", fontsize=9.5)
        plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    fig.tight_layout()
    fig.savefig(Q20 / "q20_unit_steps_heatmap.png", dpi=115)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Map
# ---------------------------------------------------------------------------------------------
def residual_map(st, units):
    fig = plt.figure(figsize=(17, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0])
    if COAST.exists():
        z = sio.loadmat(str(COAST))
        ax.plot(z["coastlon"].ravel(), z["coastlat"].ravel(), color="0.5", lw=0.6)
    keys = [k for k, *_ in ALL_UNITS]
    U = units.set_index("stream")
    amps = np.array([metric(st, k, "cum_other", "amp") for k in keys])
    peaks = np.array([metric(st, k, "cum_other", "peak") for k in keys])
    lat = np.array([U.loc[k, "lat"] for k in keys])
    lon = np.array([U.loc[k, "lon"] for k in keys])
    # spread the co-located twins so that they are visible
    jit = {"0-20000-0-06240_A": (-0.9, 0.5), "0-20000-0-06240_B": (0.9, 0.5), "0-20000-0-06240_C": (-0.9, -0.5),
           "0-20000-0-06240_D": (0.9, -0.5), "0-20000-0-12843_A": (-0.6, 0), "0-20000-0-12843_B": (0.6, 0),
           "0-20000-0-06610_A": (-1.6, 0.6)}
    lon_p = lon + np.array([jit.get(k, (0, 0))[0] for k in keys])
    lat_p = lat + np.array([jit.get(k, (0, 0))[1] for k in keys])
    sc = ax.scatter(lon_p, lat_p, s=12 * amps, c=peaks, cmap="twilight_shifted", vmin=100, vmax=250,
                    edgecolor="k", lw=0.7, zorder=3)
    twin_sites = {"0-20000-0-06240": "Schiphol A/B/C/D", "0-20000-0-12843": "Budapest A/B"}
    done = set()
    for k, x, y, a in zip(keys, lon_p, lat_p, amps):
        site = k.rsplit("_", 1)[0]
        if site in twin_sites:
            if site in done:
                continue
            done.add(site)
            members = [kk for kk in keys if kk.startswith(site)]
            txt = f"{twin_sites[site]} ±" + "/".join(f"{amps[keys.index(kk)]:.0f}" for kk in members)
            ax.text(max(lon_p[keys.index(kk)] for kk in members) + 0.7, max(lat_p[keys.index(kk)] for kk in members) + 0.4,
                    txt, fontsize=7.5, zorder=4)
            continue
        if k == PAYERNE[0]:
            ax.text(x - 0.6, y + 0.35, f"{U.loc[k, 'label']} ±{a:.0f}", fontsize=7.5, zorder=4, color="tab:red",
                    ha="right")
            continue
        ax.text(x + 0.5, y + 0.35, f"{U.loc[k, 'label']} ±{a:.0f}", fontsize=7.5, zorder=4)
    cb = plt.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.05, pad=0.09, aspect=40)
    cb.set_label("peak day of year of the residual cycle (solstice = DOY 172)")
    ax.set_xlim(-12, 31)
    ax.set_ylim(35, 62)
    ax.set_aspect(1 / np.cos(np.deg2rad(48)))
    ax.set_xlabel("longitude [°E]")
    ax.set_ylabel("latitude [°N]")
    ax.set_title("Residual annual amplitude after all corrections (κ synoptic); marker area ∝ amplitude\n"
                 "(Schiphol A-D and Budapest A/B spread for visibility; Payerne in red = abnormal unit)", fontsize=9.5)
    ax.grid(alpha=0.3)
    # bar chart sorted by latitude
    ax = fig.add_subplot(gs[1])
    order = np.argsort(lat)
    y = np.arange(len(keys))
    a_raw = np.array([metric(st, k, "raw", "amp") for k in keys])[order]
    a_syn = amps[order]
    a_tw = np.array([metric(st, k, "final_ktwin", "amp") for k in keys])[order]
    a_full = np.array([metric(st, k, "final_kfull", "amp") for k in keys])[order]
    ax.barh(y - 0.3, a_raw, 0.2, color=C_RAW, label="operational")
    ax.barh(y - 0.1, a_syn, 0.2, color=C_SYN, label="all corrections, κ synoptic")
    ax.barh(y + 0.1, a_tw, 0.2, color=C_TWIN, label="all corrections, κ twin f")
    ax.barh(y + 0.3, a_full, 0.2, color=C_FULL, label="all corrections, κ full T_int")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{U.loc[keys[i], 'label']} ({lat[i]:.0f}°N)" for i in order], fontsize=8)
    ax.set_xlabel("annual amplitude of ln C_L [%]")
    ax.set_title("Per unit, sorted by latitude", fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    fig.savefig(Q20 / "q20_residual_map.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Twins
# ---------------------------------------------------------------------------------------------
def twins_figure(tw, kf, units):
    fig, axs = plt.subplots(1, 3, figsize=(22, 5.8), gridspec_kw=dict(width_ratios=[1.5, 1.0, 1.1]))
    ax = axs[0]
    pairs = list(dict.fromkeys(zip(tw.site, tw.pair)))
    vers = [("raw", C_RAW, "operational"), ("final", C_SYN, "all corrections, κ synoptic"),
            ("final_ktwin", C_TWIN, "all corrections, κ twin f"), ("final_kfull", C_FULL, "all corrections, κ full T_int")]
    x = np.arange(len(pairs))
    for j, (v, col, lab) in enumerate(vers):
        a = [tw[(tw.site == s_) & (tw.pair == p_) & (tw.version == v)].amp.iloc[0] for s_, p_ in pairs]
        ax.bar(x + (j - 1.5) * 0.2, a, 0.2, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([("Schiphol " if "06240" in s_ else "Budapest ") + p_ for s_, p_ in pairs])
    ax.set_ylabel("annual amplitude of the same-night ratio ln(C_1/C_2) [%]")
    ax.set_title("Co-located twins: atmosphere and any common gain cancel in the ratio", fontsize=9.5)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    ax = axs[1]
    y = np.arange(len(kf))
    ax.errorbar(kf.f, y, xerr=[kf.f - kf.ci_lo, kf.ci_hi - kf.f], fmt="o", color="k", capsize=4)
    ax.axvline(0, color=C_SYN, lw=1.2, ls="-", label="f = 0: κ on synoptic T only (main chain)")
    ax.axvline(1, color=C_FULL, lw=1.2, ls="--", label="f = 1: κ on the full T_int")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.sites}\n({r.nights} nights)" for r in kf.itertuples()])
    ax.set_xlabel("fraction f of the seasonal κ·T_int term expressed between twins")
    ax.set_title("Seasonal κ fraction from the twins (night + unit fixed effects,\nmonth-block bootstrap 95 % CI)",
                 fontsize=9.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.6, len(kf) - 0.4)
    # (c) between-unit test
    ax = axs[2]
    bu = pd.read_csv(Q20 / "q20_between_units.csv")
    ph = pd.read_csv(Q20 / "q20_between_units_phasors.csv")
    free = bu[bu.model == "free f"].iloc[0]
    from q20_lib import W
    e = np.exp(1j * W * free.common_peak)
    P = ph.P_amp.values * np.exp(1j * W * ph.P_peak.values)
    S = ph.S_amp.values * np.exp(1j * W * ph.S_peak.values)
    x = np.real(S * np.conj(e))
    y = np.real(P * np.conj(e))
    U = units.set_index("stream")
    ax.scatter(x, y, color=C_SYN, s=30, zorder=3)
    for k, xx, yy in zip(ph.stream, x, y):
        ax.text(xx + 0.2, yy + 0.2, U.loc[k, "label"], fontsize=7)
    xl = np.array([min(0, x.min()) - 1, x.max() + 2])
    ax.plot(xl, free.common_amp + free.f * xl, color="k", lw=1.2,
            label=f"fit: G = ±{free.common_amp:.1f} % (DOY {free.common_peak:.0f}) + f·S, f = {free.f:.2f} "
                  f"[{free.f_lo:.2f}, {free.f_hi:.2f}]")
    ax.plot(xl, free.common_amp + xl, color=C_FULL, lw=1, ls="--", label="slope 1 (κ on the full T_int)")
    ax.axhline(free.common_amp, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("seasonal κ·T_int term S of the unit, in-phase amplitude [%]")
    ax.set_ylabel("residual cycle after all corrections (κ synoptic), in-phase [%]")
    ax.set_title("Between units: the residual grows with the unit's κ × seasonal T_int swing\n"
                 "(projection on the phase of the common term G; 20 units, unit bootstrap)", fontsize=9.5)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q20 / "q20_twins.png", dpi=125)
    plt.close(fig)


def main():
    s, st, terms, units, net, tw, kf, kal = load()
    for key, *_ in ALL_UNITS:
        if key in set(units.stream):
            unit_figure(s, st, units, kal, key)
    grid_figures(s, st, units, kal)
    waterfall(s, st, terms, net)
    heatmap(st, units)
    residual_map(st, units)
    twins_figure(tw, kf, units)
    print("figures ->", Q20)


if __name__ == "__main__":
    main()
