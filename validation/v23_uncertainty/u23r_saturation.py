# -*- coding: utf-8 -*-
"""Review: a saturating vertical term, sigma_v' = sigma_v / sqrt(1 + (sigma_v / c)^2).

Physical rationale: the window search keeps the cleanest part of the column, so the error of the
window constant does not grow without bound with the gradient measured ACROSS the probed range
(the steepest nights scatter by 12-19 % on the time-series truth, an upper bound of their
calibration error, while |g| x D reports 25-35 %). One parameter, the cap c (1 sigma).

For each cap: pooled time-series pull, pulls in the bottom / top decile of the reported sigma, the
spread of log(pull) over the deciles, the twin pull, the share of nights over the 0.6 (2 sigma) QC
limit, the AUC of the full budget as a QC score (bad > 25 % vs good < 10 %), and the unit Spearman.
Noise: "shipped" = per-gate (def_*), "revised" = window-mean estimator + molecular ASL (rev_*, review
harness >= fresh run). Splits: all nights, 2025, 2026.

Usage: python u23r_saturation.py <set> <method> [shipped|revised]
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
CAPS = (np.inf, 0.30, 0.20, 0.15, 0.12, 0.10)


def saturate(sv, cap):
    sv = np.asarray(sv, float)
    return sv if not np.isfinite(cap) else sv / np.sqrt(1.0 + (sv / cap) ** 2)


def components(d, noise):
    if noise == "revised":
        return d.rev_sig_vert.values, d.rev_sig_noise.values, d.rev_sig_mol.values
    return d.def_sig_vert.values, d.def_sig_noise.values, d.def_sig_mol.values


def metrics(d, tw, noise, cap):
    import pandas as pd
    from scipy.stats import mannwhitneyu, spearmanr
    sv, sn, sm = components(d, noise)
    sv = saturate(sv, cap)
    base = d.sig_ens.values ** 2 + np.nan_to_num(sn) ** 2 + sv ** 2
    sig = np.sqrt(base + np.nan_to_num(sm) ** 2 + np.nan_to_num(d.def_sig_aer.values) ** 2)
    out = {}
    chm = ((d.type == "CHM15k") & (d.key != C.PAYERNE) & np.isfinite(d.anom)).values
    a, s = d.anom.values[chm], sig[chm]
    st = C.pull_stats(a / s)
    q = pd.qcut(s, 10, labels=False)
    dec = np.array([C.pull_stats(a[q == k] / s[q == k])["rstd"] for k in range(10)])
    out.update(n=st["n"], ts_pull=st["rstd"], ts_out2=st["out2"], bottom=dec[0], top=dec[-1],
               decile_log_sd=float(np.std(np.log(dec))), med_sig=100 * np.median(s))
    allm = np.isfinite(d.anom.values)
    u = 2 * sig[allm]
    aa = np.abs(d.anom.values[allm])
    bad, good = aa > np.log(1.25), aa < np.log(1.10)
    out["qc06_rejected"] = float(np.mean(u > 0.6))
    out["qc06_good_among_rej"] = float(np.mean(good[u > 0.6])) if np.any(u > 0.6) else np.nan
    out["auc"] = float(mannwhitneyu(u[bad], u[good]).statistic / (bad.sum() * good.sum()))
    for typ in ("CL61", "Mini-MPL"):
        m = ((d.type == typ) & np.isfinite(d.anom)).values
        if m.sum() >= 10:
            out[f"pull_{typ}"] = C.pull_stats(d.anom.values[m] / sig[m])["rstd"]
    # unit ranking
    g = pd.DataFrame(dict(key=d.key.values[allm], a=d.anom.values[allm], s=sig[allm]))
    ut = g.groupby("key").agg(obs=("a", C.rstd), rep=("s", "median"), n=("a", "size"))
    ut = ut[ut.n >= 8]
    out["spearman_units"] = float(spearmanr(ut.rep, ut.obs)[0])
    if tw is not None and len(tw):
        nomol = dict(zip(zip(d.key, d.date), np.sqrt(base)))
        sa = np.array([nomol.get((x, y), np.nan) for x, y in zip(tw.ka, tw.date)])
        sb = np.array([nomol.get((x, y), np.nan) for x, y in zip(tw.kb, tw.date)])
        tp = C.pull_stats(tw.anom.values / np.hypot(sa, sb))
        out.update(twin_pull=tp["rstd"], twin_out2=tp["out2"])
    return out


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    noise = sys.argv[3] if len(sys.argv) > 3 else "shipped"
    p = prepare(set_name, method)
    d = p["df"]
    d = d[d.ok & np.isfinite(d.sig23)].reset_index(drop=True)
    if noise == "revised":
        d = d[np.isfinite(d.rev_sig_vert)].reset_index(drop=True)
    tw = p.get("tw")
    rows = []
    for split in ("all", "2025", "2026"):
        dd = d if split == "all" else d[[x[:4] == split for x in d.date]].reset_index(drop=True)
        tt = tw if (tw is None or split == "all") else tw[[x[:4] == split for x in tw.date]]
        for cap in CAPS:
            r = metrics(dd, tt, noise, cap)
            rows.append(dict(set=set_name, method=method, noise=noise, split=split,
                             cap=100 * cap if np.isfinite(cap) else np.inf, **r))
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / f"u23r_saturation_{set_name}_{method}_{noise}.csv", index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
