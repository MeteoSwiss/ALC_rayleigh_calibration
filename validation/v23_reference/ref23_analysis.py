# -*- coding: utf-8 -*-
"""ref23 step 5: the reference-site agreement analysis -- tables, statistics, figures.

  python validation/v23_reference/ref23_run.py run --workers 12     # harness (reads L1 once)
  python validation/v23_reference/ref23_tables.py                    # per-night tables
  python validation/v23_reference/ref23_analysis.py                  # this script

Outputs (C:/DATA/Projects/202610_v23/reference_agreement/):
  pairs.csv              one row per site x pair x night (ln rho per state / band / alpha / window,
                         every sigma term)
  summary.csv            statistics per site x pair x version (+ sensitivities)
  summary_tier.csv, summary_season.csv, summary_pooled.csv
  ref23_tables.md        markdown tables quoted in the report
  ref23_*.png            figures (landscape; height always on the Y axis)
"""
from __future__ import annotations

import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402
import ref23_pairs as P  # noqa: E402
import ref23_stats as S  # noqa: E402
import ref23_cloud as CL  # noqa: E402

OUT = RC.OUT
SITE_LABEL = {k: v["name"] for k, v in RC.SITES.items()}
TYPE = {(sk, i): v[0] for sk, s in RC.SITES.items() for i, v in s["streams"].items()}
# Cloud-uncertainty model of the reference channels (see ref23_cloud and the report): the published
# per-day number (std / median of the per-profile coefficients) -- checked against the observed
# day-to-day reproducibility; "empirical" (the channel's lag-1 reproducibility) and "sem" (standard
# error of the median over independent scenes) are the sensitivities.
CLOUD_MODEL = "reported"
CLOUD_ALT = ("empirical", "sem")


def pair_label(site, pair):
    (x, mx), (r, mr) = [p.split(":") for p in pair.split("/")]
    nm = {"ray": "Rayleigh", "cloud": "cloud"}
    tx, tr = TYPE[(site, x)], TYPE[(site, r)]
    return f"{tx} {x} ({nm[mx]}) / {tr} {r} ({nm[mr]})"


def short_label(site, pair):
    (x, mx), (r, mr) = [p.split(":") for p in pair.split("/")]
    nm = {"ray": "R", "cloud": "c"}
    return f"{SITE_LABEL[site]} {TYPE[(site, x)]}-{x}{nm[mx]} / {TYPE[(site, r)]}-{r}{nm[mr]}"


# ------------------------------------------------------------------------------- cloud sigma
def cloud_sig_models(models, clouds, kind=CLOUD_MODEL):
    """Per cloud channel: day(v) -> 1-sigma of one day value (v = the day rows used);
    drift(L) -> the extra sigma of a value carried L days (ref23_cloud.channel_model)."""
    out = {}
    for key, m in models.items():
        sw = m["sig_w"]
        if kind == "empirical":
            day = (lambda v, sw=sw: sw)
        elif kind == "reported":
            day = (lambda v: float(np.sqrt(np.nanmean(np.square(v.sig_rep)))))
        else:
            day = (lambda v: float(np.sqrt(np.nanmean(np.square(v.sig_sem)))))
        out[key] = dict(day=day, drift=m["drift"])
    return out


# ------------------------------------------------------------------------------- statistics
def site_pair_stats(pt):
    rows = []
    for (site, pair), g in pt.groupby(["site", "pair"], sort=False):
        kind = g.kind.iloc[0]
        vers = ["cloud"] if kind == "CC" else ["v22pub", "v22", "v22b", "v23", "v23pre"]
        for v in vers:
            for tag, kw in (("main", {}), ("with_mol", dict(with_mol=True)),
                            ("no_conv", dict(with_conv=False)), ("hi_band", dict(band="_hi")),
                            ("w60", dict(w="_w60")),
                            ("clean", dict(w="_clean")),
                            ("clean_hi", dict(band="_hi", w="_clean")),
                            ("dark", dict(w="_dk")),
                            ("dark_hi", dict(band="_hi", w="_dk")),
                            ("cloud_empirical", dict(cloud_model="empirical")),
                            ("cloud_sem", dict(cloud_model="sem")),
                            ("cloud_local", dict(ref="_cl"))):
                if tag == "with_mol" and kind not in ("RR_twin", "RR_mixed"):
                    continue
                if tag in ("no_conv", "clean", "clean_hi") and kind not in ("RC", "RR_mixed"):
                    continue
                if kind == "self" and tag in ("hi_band", "w60"):
                    continue
                if tag.startswith("cloud_") and kind not in ("RC", "self", "CC"):
                    continue
                s = S.summarize(g, v, **kw)
                if s is None:
                    continue
                rows.append(dict(site=site, pair=pair, kind=kind, version=v, variant=tag, **s))
            if kind != "CC" and "lnr_v22pub" in g and "lnr_v23" in g:
                # like-for-like: the nights both the published v2.2 and the v2.3 series calibrate
                gc = g[np.isfinite(g["lnr_v22pub"]) & np.isfinite(g["lnr_v23"])]
                s = S.summarize(gc, v)
                if s is not None:
                    rows.append(dict(site=site, pair=pair, kind=kind, version=v, variant="common",
                                     **s))
    return pd.DataFrame(rows)


def cl31_era_stats(pt):
    """Payerne CL31 pairs per optical-block era (swap 2026-07-07), dark-free and with the measured
    dark (the covered-telescope campaign predates the swap, so the dark is only valid before)."""
    rows = []
    sel = pt[(pt.site == "payerne") & ((pt.x == "B") | (pt.r == "B"))].copy()
    sel["era"] = np.where(sel.date < P.CL31_STEP_DAY, "before 2026-07-07", "after 2026-07-07")
    for (pair, era), g in sel.groupby(["pair", "era"]):
        v = "cloud" if g.kind.iloc[0] == "CC" else "v23"
        for tag, kw in (("main", {}), ("dark", dict(w="_dk"))):
            s = S.summarize(g, v, min_n=6, **kw)
            if s is not None and s.get("n", 0) >= 6:
                rows.append(dict(pair=pair, era=era, variant=tag, **s))
    return pd.DataFrame(rows)


def split_stats(pt, by):
    rows = []
    for (site, pair), g in pt.groupby(["site", "pair"], sort=False):
        kind = g.kind.iloc[0]
        if kind == "CC":
            continue
        for key, gg in g.groupby(by):
            for v in ("v22pub", "v23", "v23pre"):
                s = S.summarize(gg, v, min_n=6)
                if s is None or s.get("n", 0) < 6:
                    continue
                rows.append(dict(site=site, pair=pair, kind=kind, split=str(key), version=v, **s))
    return pd.DataFrame(rows)


def wv_transmission_stats(ing):
    """Two-way water-vapour transmission of the 910 nm units at 1.5 and 3 km (per-night CAMS, the
    unit's own laser line): how much a humidity error could move the comparison."""
    rows = []
    z = RC.Z_AGL
    i15, i30 = int(np.argmin(abs(z - 1500))), int(np.argmin(abs(z - 3000)))
    for site, ident in (("payerne", "C"), ("lindenberg", "C"), ("aosta", "B"), ("payerne", "B")):
        T = [np.asarray(v["corr"][ident]["t2"], float) for (s_, d_), v in ing.items()
             if s_ == site and ident in v["corr"] and v["corr"][ident].get("t2") is not None]
        if not T:
            continue
        T = np.vstack(T)
        rows.append(dict(site=site, ident=ident, type=TYPE[(site, ident)], nights=len(T),
                         t2_15_med=float(np.nanmedian(T[:, i15])),
                         t2_30_med=float(np.nanmedian(T[:, i30])),
                         t2_30_p05=float(np.nanpercentile(T[:, i30], 5)),
                         # a +-20 % humidity error scales ln T2 by +-20 %
                         dlnr_20pct_30_p05=100 * 0.2 * abs(np.log(np.nanpercentile(T[:, i30], 5)))))
    return pd.DataFrame(rows)


def fate_table(pt):
    """What v2.3 did to the nights of the v2.2 series, judged by the independent reference.

    Each night of the Rayleigh unit X (v2.2 rerun = the published constant) is classed by its fate
    in v2.3: strict tier kept; strict tier rejected by the new relative QC (u/C > 60 %); fallback
    tier kept on the same window; fallback tier moved to another window; fallback tier rejected.
    For each class: the median ln rho with the v2.2 constant and with the v2.3 constant."""
    rows = []
    sel = pt[pt.kind.isin(["RC", "RR_mixed", "RR_twin", "self"]) & np.isfinite(pt.get("lnr_v22",
                                                                                      np.nan))]
    for (site, pair), g in sel.groupby(["site", "pair"], sort=False):
        t22 = g.tier_x_v22.values
        f23 = g.v23_flag_x.values
        same = np.isclose(g.get("lncx_v23", pd.Series(np.nan, index=g.index)).values,
                          g.get("lncx_v22", pd.Series(np.nan, index=g.index)).values,
                          rtol=0, atol=1e-12)
        ok23 = np.isin(f23, RC.SUCCESS)
        cls = np.full(len(g), "other", object)
        cls[(t22 == 1) & ok23] = "strict, kept"
        cls[(t22 == 1) & (f23 == -6)] = "strict, rejected by u/C > 60 %"
        cls[(t22 == 2) & ok23 & same] = "fallback, same window"
        cls[(t22 == 2) & ok23 & ~same] = "fallback, window moved"
        cls[(t22 == 2) & ~ok23] = "fallback, rejected"
        for c in np.unique(cls):
            m = cls == c
            y22 = 100 * g.lnr_v22.values[m]
            y23 = 100 * g.get("lnr_v23", pd.Series(np.nan, index=g.index)).values[m]
            rows.append(dict(site=site, pair=pair, fate=c, n=int(m.sum()),
                             med_v22=float(np.nanmedian(y22)), scat_v22=S.robust_std(y22),
                             med_v23=float(np.nanmedian(y23)) if np.isfinite(y23).any() else np.nan,
                             n_v23=int(np.isfinite(y23).sum())))
    return pd.DataFrame(rows)


def pooled_groups(pt):
    """Groups pooled across sites for the pull / reliability figures."""
    g = {}
    pay_a = (pt.site == "payerne") & ((pt.x == "A") | (pt.r == "A"))
    g["CHM15k Rayleigh vs CL61 cloud"] = pt[(pt.kind == "RC") & (pt.r != "A") & (pt.x != "B")
                                            & ~pay_a]
    g["CHM15k Rayleigh vs CL61 Rayleigh"] = pt[(pt.kind == "RR_mixed") & ~pay_a]
    g["CL61 Rayleigh vs CL61 cloud"] = pt[pt.kind == "self"]
    g["CHM15k twins"] = pt[pt.kind == "RR_twin"]
    g["Payerne CHM15k A (abnormal) vs CL61"] = pt[pay_a & pt.kind.isin(["RC", "RR_mixed"])
                                                  & (pt.x == "A") & (pt.r == "C")]
    return g


# ------------------------------------------------------------------------------- figures
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_timeseries(pt, pairs, fname, ncols=4, title="", ylim=80):
    plt = _mpl()
    n = len(pairs)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 3.6 * nrows + 0.6), squeeze=False)
    for k, (site, pair) in enumerate(pairs):
        a = axes[k // ncols][k % ncols]
        g = pt[(pt.site == site) & (pt.pair == pair)].copy()
        if g.empty:
            a.set_title(short_label(site, pair) + " (no data)", fontsize=9)
            a.axis("off")
            continue
        g["t"] = pd.to_datetime(g.date)
        cc = g.kind.iloc[0] == "CC"
        vlist = ((("cloud", "#8e44ad", 0.0, 4, "cloud / cloud"),) if cc else
                 (("v22pub", "0.55", 0.35, 3, "v2.2 published"), ("v23", "#c0392b", 0.0, 4, "v2.3")))
        for v, col, dx, ms, lab in vlist:
            ln = S.lnr_col(g, v)
            if ln not in g:
                continue
            y = 100 * g[ln].values
            sp = 100 * S.sigma_pred(g, v)
            ok = np.isfinite(y) & np.isfinite(sp)
            t = g.t.values[ok] + np.timedelta64(int(dx * 24), "h")
            a.errorbar(t, y[ok], yerr=2 * sp[ok], fmt="none", ecolor=col, alpha=0.45, lw=0.8)
            if v == "cloud":
                a.plot(t, y[ok], "D", color=col, ms=ms - 1, label=lab)
            elif v == "v23":
                tier = g.tier_x.values[ok]
                t2 = (tier == 2)
                if g.kind.iloc[0] in ("RR_twin", "RR_mixed"):
                    t2 = t2 | (g.tier_r.values[ok] == 2)
                a.plot(t[~t2], y[ok][~t2], "o", color=col, ms=ms, label=f"{lab} (strict tier)")
                a.plot(t[t2], y[ok][t2], "o", mfc="none", mec=col, ms=ms,
                       label=f"{lab} (fallback tier)")
            else:
                a.plot(t, y[ok], "s", color=col, ms=ms - 1, label=lab)
        a.axhline(0, color="k", lw=0.8)
        med = np.nanmedian(100 * g.get(S.lnr_col(g, "cloud" if cc else "v23"),
                                       pd.Series(dtype=float)))
        if np.isfinite(med):
            a.axhline(med, color="#c0392b", lw=0.8, ls="--")
        a.set_ylim(-ylim, ylim)
        a.set_title(short_label(site, pair), fontsize=9)
        a.set_ylabel("ln rho [%]  (+-2 sigma)")
        a.grid(alpha=0.3)
        a.tick_params(axis="x", labelsize=7)
        if k == 0:
            a.legend(fontsize=7, loc="upper left")
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=120)
    plt.close(fig)


def fig_pulls(pt):
    plt = _mpl()
    groups = pooled_groups(pt)
    names = list(groups)
    fig, axes = plt.subplots(2, len(names), figsize=(4.6 * len(names), 7.6), squeeze=False)
    bins = np.linspace(-6, 6, 49)
    xx = np.linspace(-6, 6, 300)
    gauss = np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi)
    for j, nm in enumerate(names):
        g = groups[nm]
        for i, centred in enumerate((False, True)):
            a = axes[i][j]
            for v, col in (("v22pub", "0.45"), ("v23", "#c0392b")):
                y = g[S.lnr_col(g, v)].values if S.lnr_col(g, v) in g else np.array([])
                sp = S.sigma_pred(g, v)
                ok = np.isfinite(y) & np.isfinite(sp)
                if ok.sum() < 5:
                    continue
                yy = y[ok]
                if centred:
                    # centred per site x pair (the bias of each pair removed)
                    yy = yy - g[ok].groupby(["site", "pair"])[S.lnr_col(g, v)].transform("median").values
                pu = np.clip(yy / sp[ok], -5.99, 5.99)
                st = 1.4826 * np.median(np.abs(yy / sp[ok]))
                frac = 100 * np.mean(np.abs(yy / sp[ok]) < 2)
                a.hist(pu, bins=bins, density=True, histtype="step", lw=1.6, color=col,
                       label=f"{S.VERSION_LABEL[v]}: robust std {st:.2f}, {frac:.0f} % in 2 sigma "
                             f"(n={ok.sum()})")
            a.plot(xx, gauss, "k:", lw=1, label="N(0, 1)")
            a.set_title(f"{nm}\n{'centred on each pair median' if centred else 'uncentred'}",
                        fontsize=9)
            a.set_xlabel("pull = ln rho / sigma_pred" if not centred
                         else "(ln rho - median) / sigma_pred")
            a.legend(fontsize=6.5, loc="upper left")
            a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "ref23_pulls.png", dpi=120)
    plt.close(fig)


def fig_reliability(pt):
    plt = _mpl()
    groups = pooled_groups(pt)
    names = [n for n in groups if len(groups[n]) >= 20]
    fig, axes = plt.subplots(1, len(names), figsize=(4.8 * len(names), 4.8), squeeze=False)
    rel_rows = []
    for j, nm in enumerate(names):
        g = groups[nm]
        a = axes[0][j]
        top = 0
        for v, col, mk in (("v22pub", "0.45", "s"), ("v23pre", "#e59866", "v"),
                           ("v22b", "#2980b9", "^"), ("v23", "#c0392b", "o")):
            col_l = S.lnr_col(g, v)
            if col_l not in g:
                continue
            y = g[col_l].values
            sp = S.sigma_pred(g, v)
            ok = np.isfinite(y) & np.isfinite(sp)
            if ok.sum() < 20:
                continue
            yy = y[ok] - g[ok].groupby(["site", "pair"])[col_l].transform("median").values
            ss = sp[ok]
            nb = int(np.clip(ok.sum() // 25, 3, 8))
            edges = np.quantile(ss, np.linspace(0, 1, nb + 1))
            xs, ys, es = [], [], []
            for b in range(nb):
                m = (ss >= edges[b]) & (ss <= edges[b + 1])
                if m.sum() < 8:
                    continue
                rms = 1.4826 * np.median(np.abs(yy[m]))
                xs.append(100 * np.median(ss[m]))
                ys.append(100 * rms)
                es.append(100 * rms / np.sqrt(2 * (m.sum() - 1)))
                rel_rows.append(dict(group=nm, version=v, bin=b, n=int(m.sum()),
                                     sig_pred=100 * np.median(ss[m]), obs_rms=100 * rms))
            a.errorbar(xs, ys, yerr=es, fmt=mk + "-", color=col, ms=6, capsize=3,
                       label=S.VERSION_LABEL[v])
            top = max(top, max(xs + ys) if xs else 0)
        lim = 1.1 * top if top > 0 else 30
        a.plot([0, lim], [0, lim], "k:", lw=1, label="1:1 (calibrated)")
        a.set_xlim(0, lim)
        a.set_ylim(0, lim)
        a.set_xlabel("predicted 1 sigma of ln rho (bin median) [%]")
        a.set_ylabel("observed robust RMS of ln rho - pair median [%]")
        a.set_title(nm, fontsize=9)
        a.grid(alpha=0.3)
        a.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "ref23_reliability.png", dpi=120)
    plt.close(fig)
    return pd.DataFrame(rel_rows)


def fig_profiles(profs, pt):
    plt = _mpl()
    panels = [("payerne", ["A:ray/C:cloud", "A:ray/C:ray", "B:cloud/C:cloud"]),
              ("lindenberg", ["0:ray/C:cloud", "0:ray/C:ray"]),
              ("aosta", ["0:ray/B:cloud", "0:ray/B:ray"]),
              ("schiphol", [p for p in pt[pt.site == "schiphol"].pair.unique()]),
              ("budapest", ["B:ray/A:ray"]), ("magurele", ["B:ray/A:ray"])]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels), 6.0), sharey=True)
    z = RC.Z_AGL / 1000.0
    rows = []
    for a, (site, plist) in zip(axes, panels):
        for k, pair in enumerate(plist):
            state = "cloud" if pair.startswith("B:cloud") else "v23"
            arr = [v[0] for (s_, p_, d_, st_), v in profs.items()
                   if s_ == site and p_ == pair and st_ == state]
            if len(arr) < 5:
                continue
            A = 100 * np.vstack(arr).astype(float)
            n = np.sum(np.isfinite(A), axis=0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                med = np.nanmedian(A, axis=0)
                q1, q3 = np.nanpercentile(A, 25, axis=0), np.nanpercentile(A, 75, axis=0)
            keep = n >= max(5, 0.2 * len(arr))
            med[~keep] = np.nan
            q1[~keep] = np.nan
            q3[~keep] = np.nan
            col = plt.cm.tab10(k)
            a.plot(med, z, "-", color=col, lw=1.6, label=f"{pair} (n={len(arr)})")
            a.fill_betweenx(z, q1, q3, color=col, alpha=0.15)
            for zz, mm in zip(RC.Z_AGL, med):
                rows.append(dict(site=site, pair=pair, z=zz, med=mm))
        a.axvline(0, color="k", lw=0.8)
        for z0 in (0.5, 1.5, 3.0):
            a.axhline(z0, color="0.6", lw=0.6, ls=":")
        a.set_xlim(-60, 60)
        a.set_ylim(0, 5)
        a.set_title(SITE_LABEL[site], fontsize=10)
        a.set_xlabel("ln rho(z) [%]  (median, IQR over nights)")
        a.grid(alpha=0.3)
        a.legend(fontsize=6.5, loc="upper right")
    axes[0].set_ylabel("height above ground [km]")
    fig.tight_layout()
    fig.savefig(OUT / "ref23_profiles.png", dpi=120)
    plt.close(fig)
    return pd.DataFrame(rows)


def fig_lambda(pt):
    plt = _mpl()
    sel = pt[pt.kind.isin(["RC", "RR_mixed"])]
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2))
    a = axes[0]
    rows = []
    for k, ((site, pair), g) in enumerate(sel.groupby(["site", "pair"])):
        ys, lo, hi = [], [], []
        for al in RC.ALPHAS:
            col = "lnr_v23" if al == RC.ALPHA_NOMINAL else f"lnr_v23_a{al:g}"
            y = 100 * g[col].values if col in g else np.array([np.nan])
            ys.append(np.nanmedian(y))
            lo.append(np.nanpercentile(y, 25) if np.isfinite(y).any() else np.nan)
            hi.append(np.nanpercentile(y, 75) if np.isfinite(y).any() else np.nan)
            rows.append(dict(site=site, pair=pair, alpha=al, med=ys[-1]))
        a.plot(RC.ALPHAS, ys, "o-", color=plt.cm.tab10(k % 10), label=short_label(site, pair))
    a.axvline(RC.ALPHA_NOMINAL, color="k", lw=0.8, ls=":")
    a.set_xlabel("Angstrom exponent of the 910 -> 1064 nm conversion")
    a.set_ylabel("median ln rho [%]")
    a.set_title("lambda-conversion sensitivity (v2.3 constants, 500-3000 m)", fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=6.5)
    a = axes[1]
    data, labels = [], []
    for (site, pair), g in sel.groupby(["site", "pair"]):
        data.append(100 * S.sigma_conv(g, "v23"))
        (x, mx), (r, mr) = [p.split(":") for p in pair.split("/")]
        nm = {"ray": "R", "cloud": "c"}
        labels.append(f"{SITE_LABEL[site]} {x}{nm[mx]}/{r}{nm[mr]}")
    a.boxplot([d[np.isfinite(d)] for d in data], vert=True, showfliers=False)
    a.set_xticks(range(1, len(labels) + 1))
    a.set_xticklabels(labels, fontsize=7, rotation=25, ha="right")
    a.set_ylabel("sigma_conv per night (alpha = 1 +- 0.5) [%]")
    a.set_title("per-night lambda-conversion term", fontsize=10)
    a.grid(alpha=0.3)
    a = axes[2]
    for k, ((site, pair), g) in enumerate(sel.groupby(["site", "pair"])):
        sx = 100 * g.get("sx_v23", pd.Series(np.nan, index=g.index)).values
        a.plot(sx, 100 * S.sigma_conv(g, "v23"), ".", ms=3, color=plt.cm.tab10(k % 10),
               label=short_label(site, pair))
    a.plot([0, 40], [0, 40], "k:", lw=1)
    a.set_xlim(0, 40)
    a.set_ylim(0, 20)
    a.set_xlabel("v2.3 1 sigma of the CHM15k constant [%]")
    a.set_ylabel("sigma_conv [%]")
    a.set_title("conversion term vs calibration uncertainty", fontsize=10)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "ref23_lambda.png", dpi=120)
    plt.close(fig)
    return pd.DataFrame(rows)


def fig_budget(rays):
    """Median v2.3 budget components per Rayleigh unit (calibrated nights), v2.2 ensemble for
    reference -- what the predicted sigma of each pair is made of."""
    plt = _mpl()
    r = rays[(rays.variant == "v23") & rays.flag.isin(RC.SUCCESS)].copy()
    comps = list(P.COMPONENTS)
    rows = []
    for (site, ident), g in r.groupby(["site", "ident"]):
        tot = np.sqrt(np.sum(np.square(g[comps].values), axis=1))
        d = dict(site=site, ident=ident, label=f"{SITE_LABEL[site]} {TYPE[(site, ident)]} {ident}",
                 n=len(g), total=100 * np.nanmedian(tot),
                 frac_tier2=100 * float(np.mean(g.tier == 2)))
        for c in comps:
            d[c] = 100 * float(np.nanmedian(g[c]))
        rows.append(d)
    b = pd.DataFrame(rows).sort_values("total")
    fig, ax = plt.subplots(figsize=(16, 5.6))
    x = np.arange(len(b))
    w = 0.14
    cols = {"sig_ens": "0.5", "sig_vert": "#c0392b", "sig_mol": "#2980b9", "sig_aer": "#8e44ad",
            "sig_noise": "#27ae60"}
    names = {"sig_ens": "ensemble (= v2.2)", "sig_vert": "vertical consistency",
             "sig_mol": "molecular model", "sig_aer": "residual aerosol", "sig_noise": "photon noise"}
    for k, c in enumerate(comps):
        ax.bar(x + (k - 2.5) * w, b[c], w, color=cols[c], label=names[c])
    ax.bar(x + 2.5 * w, b.total, w, color="k", label="v2.3 total (quadrature)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{l.replace('Amsterdam ', '')} (n={n}, {f:.0f} % fallback)"
                        for l, n, f in zip(b.label, b.n, b.frac_tier2)],
                       fontsize=7.5, rotation=25, ha="right")
    ax.set_ylabel("median 1 sigma of the nightly constant [%]")
    ax.set_title("Reviewed v2.3 uncertainty budget of the Rayleigh units of this study "
                 "(calibrated v2.3 nights)", fontsize=10)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "ref23_budget.png", dpi=120)
    plt.close(fig)
    return b


# ------------------------------------------------------------------------------- tables
def fmt_row(s22, s23, site, pair, kind, s23_extra=None):
    def f(x, nd=1, sgn=False):
        return "–" if x is None or not np.isfinite(x) else (f"{x:+.{nd}f}" if sgn else f"{x:.{nd}f}")
    return (f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {s23['n']} | "
            f"{f(s23['bias'], 1, True)} [{f(s23['ci_lo'], 1, True)}, {f(s23['ci_hi'], 1, True)}] | "
            f"{f(s23['scatter'])} | {f(s22['sig_pred']) if s22 else '–'} → {f(s23['sig_pred'])} | "
            f"{f(s22['pull_c'], 2) if s22 else '–'} → {f(s23['pull_c'], 2)} | "
            f"{f(s22['pull_u'], 2) if s22 else '–'} → {f(s23['pull_u'], 2)} | "
            f"{f(s22['in2_u'], 0) if s22 else '–'} → {f(s23['in2_u'], 0)} % | "
            f"{f(s23['spearman'], 2)} | {S.verdict(s23)} |")


def reproduction(rays, pub):
    """The v2.2 rerun against the published constants (calibration_version 220), per stream."""
    rows = []
    rs = rays[rays.variant == "v22"].copy()
    for (site, ident), g in rs.groupby(["site", "ident"]):
        ok = g.flag.isin(RC.SUCCESS) & np.isfinite(g.c) & (g.c > 0) & ~(g.u_ens > g.c)
        mine = g[ok].set_index("date")
        pb = pub[(pub.site == site) & (pub.ident == ident) & (pub.method == "ray")].set_index("date")
        both = mine.index.intersection(pb.index)
        r = mine.loc[both, "c"].values / pb.loc[both, "c"].values
        ur = (mine.loc[both, "u_ens"].values / pb.loc[both, "u"].values) if len(both) else []
        rows.append(dict(site=site, ident=ident, type=TYPE[(site, ident)],
                         nights_rerun=len(g), ok_rerun=int(ok.sum()), published=len(pb),
                         both=len(both), only_rerun=int(len(mine.index.difference(pb.index))),
                         only_published=int(len(pb.index.difference(mine.index))),
                         identical_pct=(100 * float(np.mean(np.abs(r - 1) < 1e-9)) if len(both)
                                        else np.nan),
                         within_1pct=(100 * float(np.mean(np.abs(r - 1) < 0.01)) if len(both)
                                      else np.nan),
                         u_identical_pct=(100 * float(np.mean(np.abs(np.asarray(ur) - 1) < 1e-6))
                                          if len(both) else np.nan)))
    return pd.DataFrame(rows)


def nonidentical_nights(rays, pub):
    """The common nights whose v2.2 rerun constant differs from the published one (same selection as
    `reproduction`), one row per unit-night with the ratio rerun / published."""
    rows = []
    rs = rays[rays.variant == "v22"]
    for (site, ident), g in rs.groupby(["site", "ident"]):
        ok = g.flag.isin(RC.SUCCESS) & np.isfinite(g.c) & (g.c > 0) & ~(g.u_ens > g.c)
        mine = g[ok].set_index("date")
        pb = pub[(pub.site == site) & (pub.ident == ident) & (pub.method == "ray")].set_index("date")
        both = mine.index.intersection(pb.index)
        ratio = mine.loc[both, "c"].values / pb.loc[both, "c"].values
        for day, x in zip(both, ratio):
            if abs(x - 1) >= 1e-9:
                rows.append(dict(site=site, ident=ident, date=day, ratio=float(x)))
    return pd.DataFrame(rows, columns=["site", "ident", "date", "ratio"])


def main():
    warnings.filterwarnings("ignore")
    rays, clouds, pub, ing = P.load_tables()
    rep = reproduction(rays, pub)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    nonidentical_nights(rays, pub).to_csv(OUT / "reproduction_nights.csv", index=False)
    print(rep.to_string())
    ctab, sfa, cmp, models = CL.run(make_figure=True)
    sigm = cloud_sig_models(models, clouds, CLOUD_MODEL)
    alts = {k: cloud_sig_models(models, clouds, k) for k in CLOUD_ALT}
    pt, profs = P.build_pair_table(rays, clouds, pub, ing, sigm, alt_models=alts)
    pt.to_csv(OUT / "pairs.csv", index=False)
    with open(OUT / "tables" / "profiles.pkl", "wb") as f:
        pickle.dump(profs, f, protocol=pickle.HIGHEST_PROTOCOL)
    st = site_pair_stats(pt)
    st.to_csv(OUT / "summary.csv", index=False)
    tier = split_stats(pt, pt.tier_x.map({1: "strict", 2: "fallback"}).fillna("other"))
    tier.to_csv(OUT / "summary_tier.csv", index=False)
    seas = split_stats(pt, "season")
    seas.to_csv(OUT / "summary_season.csv", index=False)
    cl31_era_stats(pt).to_csv(OUT / "summary_cl31_era.csv", index=False)
    pooled = []
    for nm, g in pooled_groups(pt).items():
        g = g.assign(sp_key=g.site + "|" + g.pair)
        for v in ("v22pub", "v22b", "v23", "v23pre"):
            for tag, kw in (("main", {}), ("hi_band", dict(band="_hi")),
                            ("with_mol", dict(with_mol=True)), ("no_conv", dict(with_conv=False))):
                s = S.summarize(g, v, center_by="sp_key", **kw)
                if s:
                    pooled.append(dict(group=nm, version=v, variant=tag, **s))
        # like-for-like: nights where both v2.2 published and v2.3 exist
        gc = g[np.isfinite(g.get("lnr_v22pub", np.nan)) & np.isfinite(g.get("lnr_v23", np.nan))]
        for v in ("v22pub", "v23"):
            s = S.summarize(gc, v, center_by="sp_key")
            if s:
                pooled.append(dict(group=nm, version=v, variant="common", **s))
    pd.DataFrame(pooled).to_csv(OUT / "summary_pooled.csv", index=False)
    # figures
    cl61_pairs = [(s_, p_) for s_ in ("payerne", "lindenberg", "aosta")
                  for p_ in [f"{a}:{ma}/{b}:{mb}" for (a, ma), (b, mb) in RC.PAIRS[s_]]
                  if not (s_ == "payerne" and p_ == "A:ray/B:cloud")]
    twin_pairs = [(s_, f"{a}:{ma}/{b}:{mb}") for s_ in ("schiphol", "budapest", "magurele")
                  for (a, ma), (b, mb) in RC.PAIRS[s_]]
    fig_timeseries(pt, cl61_pairs, "ref23_timeseries_cl61.png", ncols=4,
                   title="ln rho per night, v2.3 (red, +-2 sigma) vs v2.2 published (grey)")
    fig_timeseries(pt, twin_pairs, "ref23_timeseries_twins.png", ncols=4,
                   title="CHM15k twins: ln rho per night, v2.3 (red, +-2 sigma) vs v2.2 published (grey)",
                   ylim=50)
    fig_pulls(pt)
    rel = fig_reliability(pt)
    rel.to_csv(OUT / "reliability.csv", index=False)
    prof = fig_profiles(profs, pt)
    prof.to_csv(OUT / "profiles_median.csv", index=False)
    lam = fig_lambda(pt)
    lam.to_csv(OUT / "lambda_alpha.csv", index=False)
    bud = fig_budget(rays)
    bud.to_csv(OUT / "budget_units.csv", index=False)
    fate = fate_table(pt)
    fate.to_csv(OUT / "fate_v22_to_v23.csv", index=False)
    wv_transmission_stats(ing).to_csv(OUT / "wv_transmission.csv", index=False)
    # markdown summary table
    lines = ["| site | pair (X / R) | n | bias ln rho [%] [95 % CI] | observed sigma [%] | "
             "predicted sigma [%] v2.2 → v2.3 | centred pull v2.2 → v2.3 | uncentred pull v2.2 → v2.3 | "
             "within 2 sigma v2.2 → v2.3 | Spearman v2.3 | verdict (v2.3) |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    m = st[st.variant == "main"]
    for (site, pair), g in m.groupby(["site", "pair"], sort=False):
        kind = g.kind.iloc[0]
        if kind == "CC":
            r = g[g.version == "cloud"]
            if len(r) and r.iloc[0].n >= 8:
                s = r.iloc[0].to_dict()
                lines.append(fmt_row(s, s, site, pair, kind))
            continue
        r22 = g[g.version == "v22pub"]
        r23 = g[g.version == "v23"]
        if not len(r23) or r23.iloc[0].n < 8:
            continue
        lines.append(fmt_row(r22.iloc[0].to_dict() if len(r22) and r22.iloc[0].n >= 8 else None,
                             r23.iloc[0].to_dict(), site, pair, kind))
    (OUT / "ref23_tables.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
