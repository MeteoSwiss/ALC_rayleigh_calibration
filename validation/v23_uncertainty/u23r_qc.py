# -*- coding: utf-8 -*-
"""Review attack: what does the relative QC (u/C_L > 0.6, 2 sigma) cost, and is 0.6 well chosen?

* "bad" night = |time-series anomaly| > 25 % (the doc's definition), "good" = |anomaly| < 10 %.
* ROC of the gate (fraction of bad nights rejected vs fraction of good nights rejected) for the full
  budget and for single ingredients (vertical term, |g|, noise, the v2.2 ensemble) and a random gate.
* Cost per unit (rejected fraction, good nights lost), Payerne A.
* eprof_v2.3 interplay: is the gate redundant with the repaired tier 2 (rejects the same nights) or
  complementary? Needs both review runs (same nights).

Usage: python u23r_qc.py [orig|fresh]
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
from u23r_truth import prepare  # noqa: E402

warnings.filterwarnings("ignore")
BAD = np.log(1.25)
GOOD = np.log(1.10)


def roc(score, bad, good, thresholds):
    tpr = np.array([np.mean(score[bad] > t) for t in thresholds])
    fpr = np.array([np.mean(score[good] > t) for t in thresholds])
    # AUC = P(score_bad > score_good) (Mann-Whitney)
    from scipy.stats import mannwhitneyu
    sb, sg = score[bad], score[good]
    sb, sg = sb[np.isfinite(sb)], sg[np.isfinite(sg)]
    auc = mannwhitneyu(sb, sg).statistic / (sb.size * sg.size) if sb.size and sg.size else np.nan
    return tpr, fpr, float(auc)


def main():
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    d22 = prepare(set_name, "eprof_v2.2")["df"]
    d23 = prepare(set_name, "eprof_v2.3")["df"]
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    d = d22[d22.ok & np.isfinite(d22.anom) & np.isfinite(d22.sig23)].copy()
    d["u_rel"] = 2 * d.sig23
    a = np.abs(d.anom.values)
    bad, good = a > BAD, a < GOOD
    say(f"[{set_name}, eprof_v2.2] {len(d)} nights with a truth: bad (>25 % off) {bad.mean() * 100:.1f} %,"
        f" good (<10 % off) {good.mean() * 100:.1f} %")
    thr = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0])
    rows = []
    for t in thr:
        rej = d.u_rel.values > t
        rows.append(dict(limit=t, rejected=rej.mean(), bad_rejected=rej[bad].mean(),
                         good_rejected=rej[good].mean(),
                         bad_rate_rejected=bad[rej].mean() if rej.any() else np.nan,
                         bad_rate_kept=bad[~rej].mean(),
                         enrichment=(bad[rej].mean() / bad[~rej].mean()) if rej.any() else np.nan,
                         good_among_rejected=good[rej].mean() if rej.any() else np.nan,
                         payerne_rejected=rej[d.key.values == C.PAYERNE].mean()))
    qt = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    say(qt.round(3).to_string(index=False))
    qt.to_csv(C.OUT / f"u23r_qc_{set_name}.csv", index=False)

    # ROC of alternative scores
    scores = {"v2.3 budget u/C": d.u_rel.values,
              "vertical term only": d.def_sig_vert.values,
              "|gradient|": np.abs(d.def_g.values),
              "photon noise only": d.def_sig_noise.values,
              "v2.2 ensemble u/C": 2 * d.sig_ens.values,
              "random": np.random.default_rng(0).random(len(d))}
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6))
    ts = np.unique(np.concatenate([np.quantile(v[np.isfinite(v)], np.linspace(0, 1, 201))
                                   for v in scores.values()]))
    for name, s in scores.items():
        tq = np.quantile(s[np.isfinite(s)], np.linspace(0, 1, 201))
        tpr, fpr, auc = roc(s, bad, good, tq)
        ax[0].plot(100 * fpr, 100 * tpr, "-", label=f"{name} (AUC {auc:.2f})")
        say(f"ROC {name:22s}: AUC(bad vs good) {auc:.3f}")
    rej = d.u_rel.values > 0.6
    ax[0].plot(100 * rej[good].mean(), 100 * rej[bad].mean(), "k*", ms=14, label="u/C > 0.6 (shipped)")
    ax[0].plot([0, 100], [0, 100], "k:", lw=1)
    ax[0].set_xlabel("good nights rejected (|anomaly| < 10 %) [%]")
    ax[0].set_ylabel("bad nights rejected (|anomaly| > 25 %) [%]")
    ax[0].set_title(f"ROC of the QC score ({set_name} units, eprof_v2.2)", fontsize=10)
    ax[0].legend(fontsize=7, loc="lower right")
    ax[0].grid(alpha=0.3)

    # per unit cost
    rows = []
    for key, g in d.groupby("key"):
        r = g.u_rel.values > 0.6
        aa = np.abs(g.anom.values)
        rows.append(dict(key=key, site=g.site.iloc[0] if "site" in g else key, type=g.type.iloc[0],
                         n=len(g), rejected=r.mean(), good_lost=np.mean(r & (aa < GOOD)),
                         bad_caught=(np.mean(r[aa > BAD]) if np.any(aa > BAD) else np.nan),
                         bad_rate=np.mean(aa > BAD)))
    ut = pd.DataFrame(rows).sort_values("rejected", ascending=False)
    ut.to_csv(C.OUT / f"u23r_qc_units_{set_name}.csv", index=False)
    say(ut.round(3).head(15).to_string(index=False))
    say(f"units with > 20 % of nights rejected: {int((ut.rejected > 0.2).sum())} of {len(ut)}; median "
        f"unit rejection {100 * ut.rejected.median():.1f} %")
    a1 = ax[1]
    a1.bar(np.arange(len(ut)), 100 * ut.rejected, color="tab:blue", alpha=0.6, label="rejected")
    a1.bar(np.arange(len(ut)), 100 * ut.good_lost, color="tab:green", alpha=0.8,
           label="of which good (|anomaly| < 10 %)")
    a1.plot(np.arange(len(ut)), 100 * ut.bad_rate, "k_", ms=10, mew=2, label="share of bad nights")
    a1.set_xticks(np.arange(len(ut)))
    a1.set_xticklabels([k.split("_")[0][-5:] + "_" + k.rsplit("_", 1)[1] for k in ut.key], rotation=80,
                       fontsize=6)
    a1.set_ylabel("[% of the unit's nights]")
    a1.set_title("Cost per unit of u/C > 0.6", fontsize=10)
    a1.legend(fontsize=7)
    a1.grid(alpha=0.3, axis="y")

    # eprof_v2.3 interplay (same nights)
    m22 = d22.set_index(["key", "date"])
    m23 = d23.set_index(["key", "date"])
    common = m22.index.intersection(m23.index)
    x22 = m22.loc[common]
    x23 = m23.loc[common]
    ok22 = x22.ok.values
    ok23 = x23.ok.values
    t2 = x22.tier.values == 2.0
    same_c = ok22 & ok23 & np.isclose(x22.c.values, x23.c.values, rtol=1e-12)
    say(f"eprof_v2.3 vs v2.2 on {ok22.sum()} v2.2 successes: tier 2 {100 * np.mean(t2[ok22]):.1f} %; "
        f"v2.3 keeps {100 * np.mean(ok23[ok22]):.1f} % (tier 1 identical: "
        f"{100 * np.mean(same_c[ok22 & ~t2]):.2f} %); tier-2 nights: v2.3 rejects "
        f"{100 * np.mean(~ok23[ok22 & t2]):.1f} %, moves the window of "
        f"{100 * np.mean((ok23 & ~same_c)[ok22 & t2]):.1f} %")
    qc22 = ok22 & (x22.def_exceeds.values == 1)
    qc23 = ok23 & (x23.def_exceeds.values == 1)
    say(f"u/C gate on v2.2 successes: rejects {qc22.sum()} ({100 * qc22.sum() / ok22.sum():.1f} %), of "
        f"which tier 2: {100 * np.mean(t2[qc22]):.0f} %; of those tier-2 gate rejections, v2.3's "
        f"vertical-consistency test already rejects {100 * np.mean(~ok23[qc22 & t2]):.0f} %")
    say(f"u/C gate on v2.3 successes: rejects {qc23.sum()} ({100 * qc23.sum() / max(ok23.sum(), 1):.1f} %)"
        f"; of them tier 1 {100 * np.mean(x23.tier.values[qc23] == 1):.0f} %")
    # enrichment of bad nights among the gate's rejections, v2.3
    e = d23[d23.ok & np.isfinite(d23.anom) & np.isfinite(d23.sig23)]
    r = (2 * e.sig23.values) > 0.6
    aa = np.abs(e.anom.values)
    say(f"eprof_v2.3: gate rejects {100 * r.mean():.1f} % of {len(e)} nights; bad rate rejected "
        f"{100 * np.mean(aa[r] > BAD):.1f} % vs kept {100 * np.mean(aa[~r] > BAD):.1f} % (enrichment "
        f"{np.mean(aa[r] > BAD) / max(np.mean(aa[~r] > BAD), 1e-9):.1f} x); good among rejected "
        f"{100 * np.mean(aa[r] < GOOD):.0f} %")
    # what v2.3's test removed vs what the gate removes, in terms of the truth
    a2 = ax[2]
    vt = x22[ok22 & t2]
    rej_vc = ~x23.loc[vt.index].ok.values
    aa2 = np.abs(vt.anom.values)
    fin = np.isfinite(aa2)
    bins = np.linspace(0, 0.6, 25)
    a2.hist(100 * aa2[fin & rej_vc], bins=100 * bins, alpha=0.5, density=True,
            label=f"tier-2 nights v2.3 rejects (n {int((fin & rej_vc).sum())})")
    a2.hist(100 * aa2[fin & ~rej_vc], bins=100 * bins, alpha=0.5, density=True,
            label=f"tier-2 nights v2.3 keeps (n {int((fin & ~rej_vc).sum())})")
    g_rej = qc22 & ~t2
    a2.hist(100 * np.abs(x22.anom.values[g_rej & np.isfinite(x22.anom.values)]), bins=100 * bins,
            histtype="step", lw=2, color="k", density=True,
            label=f"tier-1 nights the u/C gate rejects (n {int((g_rej & np.isfinite(x22.anom.values)).sum())})")
    a2.set_xlabel("|time-series anomaly| [%]")
    a2.set_ylabel("density")
    a2.set_title("What each gate removes (truth = 31-day median)", fontsize=10)
    a2.legend(fontsize=7)
    a2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(C.OUT / f"u23r_qc_{set_name}.png", dpi=130)
    plt.close(fig)
    (C.OUT / f"u23r_qc_{set_name}.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
