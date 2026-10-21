# -*- coding: utf-8 -*-
"""Figures of the independent review (u23r_review.md). Landscape PNGs in the review output folder.

  u23r_fig_truth.png        pulls under the truth variants (leakage, truth error)
  u23r_fig_reliability.png  conditional calibration: observed vs reported by decile (time series),
                            by quintile (twins), and the vertical term vs |g| (shipped vs revised)
  u23r_fig_noise.png        split-night test of the photon-noise term per instrument type
  u23r_fig_bracket.png      variance shares and the pulls a perfect budget would give
  u23r_fig_molecular.png    molecular climatology out of sample (latitude clip)
  u23r_fig_kalman.png       Kalman weighting: permutation null and the QC / cap variants
Usage: python u23r_figures.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
from u23r_truth import prepare  # noqa: E402
from u23r_saturation import saturate  # noqa: E402

warnings.filterwarnings("ignore")
OUT = C.OUT
CAP = 0.15


def _revised_sigma(d, cap=CAP):
    sv = saturate(d.rev_sig_vert.values, cap)
    return np.sqrt(d.sig_ens.values ** 2 + sv ** 2 + np.nan_to_num(d.rev_sig_noise.values) ** 2
                   + np.nan_to_num(d.rev_sig_mol.values) ** 2 + np.nan_to_num(d.def_sig_aer.values) ** 2)


def fig_truth(plt):
    import pandas as pd
    t = pd.read_csv(OUT / "u23r_truth.csv")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.2))
    order = ["u23 truth (published nights)", "exclude +-3 d", "exclude +-7 d", "causal (past 31 d)",
             "no kappa removal"]
    a = ax[0]
    m = t[t.population == "CHM15k (Payerne excl.)"].set_index("variant").reindex(order)
    x = np.arange(len(order))
    a.bar(x - 0.2, m.pull23, 0.4, color="tab:blue", label="raw anomaly / reported sigma")
    a.bar(x + 0.2, m.pull23_secorr, 0.4, color="tab:cyan",
          label="truth error removed: anomaly / sqrt(sigma^2 + se_median^2)")
    a.axhline(1.0, color="k", lw=1)
    a.set_xticks(x)
    a.set_xticklabels(["u23 truth\n(night excluded)", "+-3 d\nexcluded", "+-7 d\nexcluded",
                       "causal\n(past 31 d)", "no kappa\nremoval"], fontsize=8)
    a.set_ylabel("robust pull std (CHM15k, v2.3 shipped)")
    a.set_title("Time-series truth: leakage of neighbouring nights vs the truth's own error", fontsize=10)
    a.legend(fontsize=8)
    a.grid(alpha=0.3, axis="y")
    a = ax[1]
    tw = t[t.population == "twins"].set_index("variant")
    labs = ["u23 truth", "exclude +-3 d", "exclude +-7 d", "no kappa removal",
            "u23 truth, median error removed"]
    vals = [tw.loc[v, "pull23"] if v in tw.index else np.nan for v in labs]
    a.bar(np.arange(len(labs)), vals, 0.6, color="tab:orange")
    a.axhline(1.0, color="k", lw=1)
    a.set_xticks(np.arange(len(labs)))
    a.set_xticklabels(["u23 truth", "+-3 d\nexcluded", "+-7 d\nexcluded", "no kappa\nremoval",
                       "median error\nremoved"], fontsize=8)
    a.set_ylabel("robust pull std (twins, v2.3 shipped)")
    a.set_title("Same-night twin truth", fontsize=10)
    a.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_truth.png", dpi=130)
    plt.close(fig)


def _deciles(a, s, nq=10):
    import pandas as pd
    q = pd.qcut(s, nq, labels=False)
    rep = np.array([np.median(s[q == k]) for k in range(nq)])
    obs = np.array([C.rstd(a[q == k]) for k in range(nq)])
    return rep, obs


def fig_reliability(plt):
    import pandas as pd
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6))
    a = ax[0]
    styles = [("orig", "eprof_v2.2", "o", "0.5"), ("orig", "eprof_v2.3", "s", "tab:blue"),
              ("fresh", "eprof_v2.3", "^", "tab:red")]
    for set_name, method, mk, col in styles:
        try:
            d = prepare(set_name, method)["df"]
        except FileNotFoundError:
            continue
        d = d[d.ok & np.isfinite(d.anom) & np.isfinite(d.sig23) & (d.type == "CHM15k")
              & (d.key != C.PAYERNE)]
        rep, obs = _deciles(d.anom.values, d.sig23.values)
        a.plot(100 * rep, 100 * obs, mk + "-", color=col, label=f"{set_name} units, {method}: shipped")
        if "rev_sig_vert" in d and np.isfinite(d.rev_sig_vert).mean() > 0.9:
            dd = d[np.isfinite(d.rev_sig_vert)]
            rep, obs = _deciles(dd.anom.values, _revised_sigma(dd))
            a.plot(100 * rep, 100 * obs, mk + "--", color=col, mfc="none",
                   label=f"{set_name} units, {method}: revised")
    lim = [2, 45]
    a.plot(lim, lim, "k:", lw=1, label="1:1")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(lim)
    a.set_ylim(4, 30)
    a.set_xlabel("reported 1-sigma, decile median [%]")
    a.set_ylabel("observed anomaly scatter (robust) [%]")
    a.set_title("Time-series truth (CHM15k): observed vs reported per decile", fontsize=10)
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    a = ax[1]
    p = prepare("origx", "eprof_v2.3")          # twin units re-run with every field (shipped + revised)
    d, tw = p["df"], p["tw"]
    s23 = C.twin_sigma(tw, d, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
    ok = np.isfinite(s23)
    rep, obs = _deciles(tw.anom.values[ok], s23[ok], 5)
    a.plot(100 * rep, 100 * obs, "s-", color="tab:blue", label="shipped (eprof_v2.3)")
    if "rev_sig_vert" in d and np.isfinite(d.rev_sig_vert[d.ok]).mean() > 0.9:
        dd = d[d.ok].copy()
        dd["sig_rev_nomol"] = np.sqrt(np.clip(_revised_sigma(dd) ** 2 - dd.rev_sig_mol ** 2
                                              - dd.def_sig_aer ** 2, 0, None))
        sr = C.twin_sigma(tw, dd, "sig_rev_nomol")
        okr = np.isfinite(sr)
        rep, obs = _deciles(tw.anom.values[okr], sr[okr], 5)
        a.plot(100 * rep, 100 * obs, "s--", color="tab:blue", mfc="none", label="revised")
    lim = [3, 40]
    a.plot(lim, lim, "k:", lw=1, label="1:1")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(lim)
    a.set_ylim(3, 20)
    a.set_xlabel("reported pair sigma (molecular + aerosol removed), quintile median [%]")
    a.set_ylabel("observed twin-ratio anomaly scatter [%]")
    a.set_title("Twin truth (real gain changes cancel), design units", fontsize=10)
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    a = ax[2]
    for set_name, method, mk, col in styles[1:]:
        try:
            d = prepare(set_name, method)["df"]
        except FileNotFoundError:
            continue
        d = d[d.ok & np.isfinite(d.anom) & np.isfinite(d.sig23) & (d.type == "CHM15k")
              & (d.key != C.PAYERNE) & np.isfinite(d.def_g)]
        q = pd.qcut(np.abs(d.def_g), 8, labels=False)
        gq = np.array([np.median(np.abs(d.def_g[q == k])) for k in range(8)])
        obs = np.array([C.rstd(d.anom[q == k]) for k in range(8)])
        vert = np.array([np.median(d.def_sig_vert[q == k]) for k in range(8)])
        a.plot(100 * gq, 100 * obs, mk + "-", color=col, label=f"{set_name}: observed scatter")
        a.plot(100 * gq, 100 * vert, mk + ":", color=col, mfc="none",
               label=f"{set_name}: vertical term (shipped, D = 1 km)")
        a.plot(100 * gq, 100 * saturate(vert, CAP), mk + "--", color=col, mfc="none", alpha=0.7,
               label=f"{set_name}: vertical term, saturated at {100 * CAP:.0f} %")
    a.set_xscale("log")
    a.set_xlabel("|d ln C_L / dz| over the probed range, octile median [%/km]")
    a.set_ylabel("[%] (1 sigma)")
    a.set_title("The vertical term grows like |g|, the error does not (eprof_v2.3)", fontsize=10)
    a.legend(fontsize=7)
    a.grid(alpha=0.3, which="both")
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter
    ticks = [0.5, 1, 2, 3, 5, 10, 20, 30, 40]
    for j, a in enumerate(ax):
        axes = (a.xaxis, a.yaxis) if j < 2 else (a.xaxis,)
        for axis in axes:
            axis.set_major_locator(FixedLocator(ticks))
            axis.set_major_formatter(FormatStrFormatter("%g"))
            axis.set_minor_formatter(NullFormatter())
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_reliability.png", dpi=130)
    plt.close(fig)


def fig_noise(plt):
    import pandas as pd
    rows = []
    for set_name in ("fresh", "origx"):
        for method in ("eprof_v2.2",):
            try:
                d = prepare(set_name, method)["df"]
            except Exception:  # noqa: BLE001
                continue
            d = d[d.ok & np.isfinite(d.c_odd) & np.isfinite(d.c_even)]
            for typ, g in d.groupby("type"):
                n = 0.5 * (g.c_odd - g.c_even)
                h = 0.5 * (g.c_h1 - g.c_h2)
                rows.append(dict(set=set_name, type=typ, n=len(g),
                                 per_gate=C.rstd(n / g.def_sig_noise),
                                 per_gate_measured=C.rstd(n / g.sig_noise_meas),
                                 window_mean=C.rstd(n / g.sig_noise_wm),
                                 halves_window_mean=C.rstd(h / g.sig_noise_wm),
                                 ratio_prop_meas=float(np.nanmedian(g.noise_ratio_36))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "u23r_noise_summary.csv", index=False)
    print(t.round(3).to_string(index=False))
    fig, ax = plt.subplots(figsize=(15, 5))
    x = np.arange(len(t))
    for j, (col, lab, c) in enumerate((("per_gate", "per-gate, propagated / sqrt(k) (shipped)", "0.5"),
                                        ("per_gate_measured", "per-gate, measured on the binned stack", "tab:orange"),
                                        ("window_mean", "window mean, measured on the binned stack (revised)", "tab:blue"))):
        ax.bar(x + (j - 1) * 0.27, t[col], 0.27, color=c, label=lab)
    ax.axhline(1.0, color="k", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.type}\n{r.set} units (n {r.n})" for r in t.itertuples()], fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("robust std of (C_odd - C_even)/2 / sigma_noise")
    ax.set_title("Split-night test of the photon-noise term (odd vs even 300-s profiles; 1 = the "
                 "noise term is right)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_noise.png", dpi=130)
    plt.close(fig)


def fig_bracket(plt):
    import pandas as pd
    sh = pd.read_csv(OUT / "u23r_bracket_shares.csv")
    ex = pd.read_csv(OUT / "u23r_bracket_expected.csv")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.2))
    a = ax[0]
    y = np.arange(len(sh))
    cols = ["tab:blue" if k == "twin" else ("tab:red" if f else "tab:green")
            for k, f in zip(sh.kind, sh.fresh)]
    a.errorbar(sh.share_a, y, xerr=[sh.share_a - sh.share_a_lo, sh.share_a_hi - sh.share_a], fmt="none",
               ecolor="0.6")
    a.scatter(sh.share_a, y, c=cols, zorder=3)
    a.set_yticks(y)
    a.set_yticklabels([f"{r.kind}: {r.a.split('_')[0][-5:]}_{r.a[-1]} / {r.b[-1]}" + (" (fresh)" if r.fresh else "")
                       for r in sh.itertuples()], fontsize=7)
    a.set_xlabel("share of the CHM15k anomaly variance shared with the co-located unit (95 % bootstrap)")
    a.set_title("Twins share atmosphere + gain, a CL61 only the atmosphere", fontsize=10)
    a.grid(alpha=0.3)
    a = ax[1]
    a.fill_between([ex.ts_pull.min(), ex.ts_pull.max()], ex.twin_pull.min(), ex.twin_pull.max(),
                   color="0.85", label="what a perfect budget would give (share CIs, w = 0.63-1)")
    t = pd.read_csv(OUT / "u23r_truth.csv")
    def pick(variant, pop, col="pull23"):
        r = t[(t.variant == variant) & (t.population == pop)]
        return float(r[col].iloc[0]) if len(r) else np.nan
    pts = [("v2.2 ensemble, u23 truth", pick("u23 truth (published nights)", "CHM15k (Payerne excl.)", "pull22"),
            pick("u23 truth", "twins", "pull22"), "0.4"),
           ("v2.3 shipped, u23 truth", pick("u23 truth (published nights)", "CHM15k (Payerne excl.)"),
            pick("u23 truth", "twins"), "tab:blue"),
           ("v2.3 shipped, +-3 d truth", pick("exclude +-3 d", "CHM15k (Payerne excl.)"),
            pick("exclude +-3 d", "twins"), "tab:red")]
    for lab, xx, yy, c in pts:
        a.plot(xx, yy, "o", ms=9, color=c, label=lab)
    a.set_xscale("log")
    a.set_xlim(0.8, 5)
    a.set_ylim(0.3, 2.0)
    a.set_yscale("log")
    a.set_xlabel("time-series pull (CHM15k)")
    a.set_ylabel("twin pull")
    a.set_title("The bracket is wide: +-30 % on the budget level", fontsize=10)
    a.legend(fontsize=7, loc="upper left")
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_bracket.png", dpi=130)
    plt.close(fig)


def fig_molecular(plt):
    import pandas as pd
    d = pd.read_csv(OUT / "u23r_molecular_oos.csv")
    g = d.groupby(["key", "lat", "in_fit"]).apply(
        lambda x: pd.Series(dict(err=np.mean(x.b - x.clim_asl), sd=np.std(x.b - x.clim_asl),
                                 cover=np.mean(np.abs(x.b) <= 2 * np.sqrt(x.clim_agl ** 2 + 1.35 ** 2)),
                                 n=len(x)))).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    for fit, c, lab in ((True, "tab:blue", "in the fit"), (False, "tab:red", "not in the fit")):
        x = g[g.in_fit == fit]
        ax[0].errorbar(x.lat, x.err, yerr=x.sd, fmt="o", color=c, label=lab, alpha=0.8)
        ax[1].plot(x.lat, 100 * x.cover, "o", color=c, label=lab)
    for a in ax:
        a.axvspan(37.2, 55.0, color="0.9", label="fitted latitudes")
        a.axvline(62.0, color="k", ls="--", lw=1, label="clip of the module (62 N)")
        a.grid(alpha=0.3)
        a.set_xlabel("station latitude [deg N]")
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_ylabel("realised b - climatology [%] (mean +- sd per unit)")
    ax[0].set_title("US-1976 bias of C_L: CAMS per night vs the climatology", fontsize=10)
    ax[1].axhline(95.4, color="k", lw=0.8)
    ax[1].set_ylabel("nights with |b| <= 2 sigma_mol [%]")
    ax[1].set_title("Coverage of the molecular term", fontsize=10)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_molecular.png", dpi=130)
    plt.close(fig)


def fig_kalman(plt):
    import pandas as pd
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    for j, method in enumerate(("eprof_v2.2", "eprof_v2.3")):
        f = OUT / f"u23r_kalman_qc_orig_{method}.csv"
        if not f.exists():
            continue
        q = pd.read_csv(f)
        q = q[q.variant != "flat, all nights"]
        x = np.arange(len(q))
        ax[j].bar(x - 0.2, q.twin_rms, 0.4, color="tab:blue", label="twin best-estimate rms (calibration error)")
        ax[j].bar(x + 0.2, q.jitter, 0.4, color="tab:orange", label="day-to-day jitter")
        ax[j].axhline(0, color="k", lw=1)
        ax[j].set_xticks(x)
        ax[j].set_xticklabels([v.replace(", ", ",\n").replace(" + ", "\n+ ") for v in q.variant], fontsize=7)
        ax[j].set_ylabel("median change vs equal weights, all nights [%] (negative = better)")
        ax[j].set_title(f"Kalman best estimate, design units, {method}", fontsize=10)
        ax[j].legend(fontsize=7)
        ax[j].grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "u23r_fig_kalman.png", dpi=130)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    which = sys.argv[1:] or ["truth", "reliability", "noise", "bracket", "molecular", "kalman"]
    for w in which:
        globals()[f"fig_{w}"](plt)
        print("done", w)


if __name__ == "__main__":
    main()
