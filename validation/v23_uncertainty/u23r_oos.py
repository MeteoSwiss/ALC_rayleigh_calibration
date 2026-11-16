# -*- coding: utf-8 -*-
"""Review: before / after on the SAME nights, out of sample.

Budgets compared on every night of a review run (exact code values of the review harness):
  v2.2       the ensemble alone (published uncertainty)
  shipped    v2.3 as delivered (def_*: per-gate noise, molecular climatology at the window centre AGL)
  revised    the review's revision: window-mean noise + molecular climatology ASL (rev_*) and the
             vertical term saturated at 15 % (1 sigma)
Truths: the u23 time-series truth (night excluded) and the leakage-controlled one (+-3 d excluded);
twins where the set has them. Per population: pooled pull and share outside 2 sigma, pulls in the
bottom / top decile of the reported sigma, the decile spread of log(pull), QC rejection at 0.6.

Usage: python u23r_oos.py <set> [eprof_v2.2|eprof_v2.3]
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
CAP = 0.15


def sigmas(d):
    rev_v = saturate(d.rev_sig_vert.values, CAP)
    rev = np.sqrt(d.sig_ens.values ** 2 + rev_v ** 2 + np.nan_to_num(d.rev_sig_noise.values) ** 2
                  + np.nan_to_num(d.rev_sig_mol.values) ** 2 + np.nan_to_num(d.def_sig_aer.values) ** 2)
    rev_nomol = np.sqrt(d.sig_ens.values ** 2 + rev_v ** 2 + np.nan_to_num(d.rev_sig_noise.values) ** 2)
    shp_nomol = np.sqrt(np.clip(d.sig23.values ** 2 - d.def_sig_mol.values ** 2
                                - d.def_sig_aer.values ** 2, 0, None))
    return {"v2.2": (d.sig_ens.values, d.sig_ens.values), "shipped": (d.sig23.values, shp_nomol),
            "revised": (rev, rev_nomol)}


def stats(a, s):
    import pandas as pd
    ok = np.isfinite(a) & np.isfinite(s) & (s > 0)
    a, s = a[ok], s[ok]
    st = C.pull_stats(a / s)
    out = dict(n=st["n"], pull=st["rstd"], out2=st["out2"], med_sig=100 * np.median(s),
               obs=100 * C.rstd(a))
    if a.size >= 100:
        q = pd.qcut(s, 10, labels=False, duplicates="drop")
        dec = np.array([C.pull_stats(a[q == k] / s[q == k])["rstd"] for k in np.unique(q)])
        out.update(bottom=dec[0], top=dec[-1], decile_log_sd=float(np.std(np.log(dec))))
    return out


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "fresh"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.3"
    p = prepare(set_name, method)
    d = p["df"]
    d = d[d.ok & np.isfinite(d.sig23) & np.isfinite(d.rev_sig_vert)].reset_index(drop=True)
    T = C.Truth()
    d = C.add_truth(d, T, exclude_days=3.0, col="anom3")
    S = sigmas(d)
    pops = {"CHM15k (Payerne excl.)": (d.type == "CHM15k") & (d.key != C.PAYERNE),
            "CL61": d.type == "CL61", "Mini-MPL": d.type == "Mini-MPL"}
    rows = []
    for truth, col in (("u23 truth", "anom"), ("+-3 d excluded", "anom3")):
        for pop, m in pops.items():
            m = m.values & np.isfinite(d[col].values)
            if m.sum() < 10:
                continue
            for b, (s, _) in S.items():
                r = stats(d[col].values[m], s[m])
                qc = float(np.mean(2 * s[m] > 0.6))
                rows.append(dict(set=set_name, method=method, truth=truth, population=pop, budget=b,
                                 qc06=qc, **r))
    tw = p.get("tw")
    if tw is not None and len(tw):
        for b, (_, snm) in S.items():
            smap = dict(zip(zip(d.key, d.date), snm))
            sa = np.array([smap.get((x, y), np.nan) for x, y in zip(tw.ka, tw.date)])
            sb = np.array([smap.get((x, y), np.nan) for x, y in zip(tw.kb, tw.date)])
            r = stats(tw.anom.values, np.hypot(sa, sb))
            rows.append(dict(set=set_name, method=method, truth="twins", population="twins", budget=b, **r))
    co = p.get("coloc")
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / f"u23r_oos_{set_name}_{method}.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    print(out.round(3).to_string(index=False))
    # per unit (u23 truth): pulls shipped vs revised
    rows = []
    for key, g in d[np.isfinite(d.anom)].groupby("key"):
        if len(g) < 10:
            continue
        idx = g.index.values
        r = dict(key=key, site=g.site.iloc[0], type=g.type.iloc[0], n=len(g), obs=100 * C.rstd(g.anom))
        for b, (s, _) in S.items():
            r[f"rep_{b}"] = 100 * np.median(s[idx])
            r[f"pull_{b}"] = C.pull_stats(g.anom.values / s[idx])["rstd"]
        rows.append(r)
    ut = pd.DataFrame(rows).sort_values("obs", ascending=False)
    ut.to_csv(C.OUT / f"u23r_oos_units_{set_name}_{method}.csv", index=False)
    from scipy.stats import spearmanr
    print(ut.round(2).to_string(index=False))
    for b in S:
        print(f"Spearman(unit reported, observed) {b}: {spearmanr(ut[f'rep_{b}'], ut.obs)[0]:+.2f}; "
              f"sd(log unit pull) {np.std(np.log(ut[f'pull_{b}'])):.3f}")
    if co is not None and len(co):
        print(f"fresh co-located CHM15k/CL61 (Camborne): {len(co)} same-night pairs, ratio anomaly "
              f"robust sigma {100 * C.rstd(co.anom):.1f} %")


if __name__ == "__main__":
    main()
