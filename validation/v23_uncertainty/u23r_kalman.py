# -*- coding: utf-8 -*-
"""Review attack (claim 6): is the Kalman improvement from the v2.3 weighting significant?

Metrics as in u23_analysis.cmd_kalman (twin best-estimate consistency rms / p95, day-to-day jitter,
one-step-ahead innovation of the more certain half). Two null tests:

* permutation: the SAME weights shuffled across the nights of each unit (same distribution, no
  information on which night is uncertain). If shuffled weights give the same gain, the gain comes
  from the weighting mechanics (fewer effective observations -> smoother), not from the budget;
* bootstrap over units / pairs of the median relative change.

Usage: python u23r_kalman.py [orig] [eprof_v2.2|eprof_v2.3] [n_perm]
"""
from __future__ import annotations

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
import u23_analysis as A  # noqa: E402
from u23r_truth import prepare  # noqa: E402

warnings.filterwarnings("ignore")


def series_metrics(times, vals, unc, sig):
    days, st = A._kalman_series(times, vals, unc)
    tday = np.array([(t - C.EPOCH).days for t in times])
    e = A._innovations(days, st, tday, np.log(vals))
    good = sig <= np.nanmedian(sig)
    return (days, st), dict(innov_good=C.pull_stats(e[good])["rstd"],
                            innov_all=C.pull_stats(e)["rstd"],
                            jitter=float(np.nanstd(np.diff(np.log(st)))))


def twin_metric(sa, sb):
    (da, xa), (db, xb) = sa, sb
    common, ia, ib = np.intersect1d(da, db, return_indices=True)
    if common.size < 60:
        return np.nan, np.nan
    r = np.log(xa[ia] / xb[ib])
    dev = r - A.running_median(common.astype(float), r)
    return float(np.sqrt(np.nanmean(dev ** 2))), float(np.nanpercentile(np.abs(dev), 95))


def run_all(units, unc_mode, seed=None):
    """Metrics for every unit and twin pair. unc_mode: 'flat' | 'w23' | 'perm'."""
    rs = np.random.default_rng(seed)
    met, ser = {}, {}
    for key, (times, vals, sig) in units.items():
        if unc_mode == "flat":
            unc = None
        else:
            s = sig if unc_mode == "w23" else rs.permutation(sig)
            unc = 2 * s * vals
        ser[key], met[key] = series_metrics(times, vals, unc, sig)
    tw = {}
    for wmo, idents in C.TWINS.items():
        for i, a in enumerate(idents):
            for b in idents[i + 1:]:
                ka, kb = f"{wmo}_{a}", f"{wmo}_{b}"
                if ka in ser and kb in ser:
                    tw[f"{ka}/{kb}"] = twin_metric(ser[ka], ser[kb])
    return met, tw


def summary(base_m, base_t, m, t):
    """Median relative change vs the flat filter (%), per metric."""
    out = {}
    for k in ("jitter", "innov_good", "innov_all"):
        v = [100 * (m[u][k] / base_m[u][k] - 1) for u in base_m if np.isfinite(base_m[u][k])]
        out[k] = float(np.median(v))
    for j, k in enumerate(("twin_rms", "twin_p95")):
        v = [100 * (t[p][j] / base_t[p][j] - 1) for p in base_t if np.isfinite(base_t[p][j])]
        out[k] = float(np.median(v))
    return out


def _perm_worker(args):
    units, seed = args
    return run_all(units, "perm", seed)


def qc_variants(ok):
    """Downstream value of the relative QC and of a saturated vertical term: Kalman best estimates
    (weighted) with / without the u/C > 0.6 rejection, against the equal-weight filter on ALL nights.
    Twin consistency is the metric that isolates the calibration error (common gain cancels)."""
    import pandas as pd
    from u23r_saturation import saturate
    out = []
    variants = {}
    specs = [("flat, all nights", None, False, False), ("weighted", np.inf, False, False),
             ("weighted + QC 0.6", np.inf, True, False),
             ("weighted, vertical capped 15 %", 0.15, False, False),
             ("weighted, capped 15 % + QC 0.6 on the uncapped budget", 0.15, True, False)]
    if "rev_sig_vert" in ok and np.isfinite(ok.rev_sig_vert).mean() > 0.95:
        specs.append(("weighted, revised budget (noise + ASL + cap 15 %)", 0.15, False, True))
    for lab, cap, qc, revised in specs:
        units = {}
        for key, g in ok.groupby("key"):
            g = g.sort_values("tday")
            if qc:
                g = g[~(2 * g.sig23 > 0.6)]
            if len(g) < 30:
                continue
            if revised:
                sv = saturate(g.rev_sig_vert.values, cap)
                sig = np.sqrt(g.sig_ens.values ** 2 + sv ** 2 + np.nan_to_num(g.rev_sig_noise.values) ** 2
                              + np.nan_to_num(g.rev_sig_mol.values) ** 2
                              + np.nan_to_num(g.def_sig_aer.values) ** 2)
            else:
                sv = g.def_sig_vert.values if cap is None else saturate(g.def_sig_vert.values, cap)
                sig = np.sqrt(g.sig23.values ** 2 - g.def_sig_vert.values ** 2 + sv ** 2)
            times = [C.EPOCH + timedelta(days=int(x)) for x in g.tday]
            units[key] = (times, g.c.values, sig)
        variants[lab] = run_all(units, "flat" if cap is None else "w23")
    bm, bt = variants["flat, all nights"]
    for lab, (m, t) in variants.items():
        s = summary(bm, bt, m, t) if lab != "flat, all nights" else {k: 0.0 for k in
                                                                    ("jitter", "innov_good", "innov_all",
                                                                     "twin_rms", "twin_p95")}
        out.append(dict(variant=lab, **s))
    return pd.DataFrame(out)


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    n_perm = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    df = prepare(set_name, method)["df"]
    ok = df[df.ok & np.isfinite(df.sig23)]
    units = {}
    for key, g in ok.groupby("key"):
        g = g.sort_values("tday")
        if len(g) < 30:
            continue
        times = [C.EPOCH + timedelta(days=int(x)) for x in g.tday]
        units[key] = (times, g.c.values, g.sig23.values)
    bm, bt = run_all(units, "flat")
    wm, wt = run_all(units, "w23")
    obs = summary(bm, bt, wm, wt)
    with ProcessPoolExecutor(max_workers=24) as ex:
        perms = list(ex.map(_perm_worker, [(units, s) for s in range(n_perm)]))
    null = pd.DataFrame([summary(bm, bt, m, t) for m, t in perms])
    rows = []
    for k, v in obs.items():
        nv = null[k].values
        # one-sided: an improvement is a NEGATIVE change; p = share of shuffled weights doing as well
        rows.append(dict(metric=k, observed_change_pct=v, null_median=np.median(nv),
                         null_p05=np.percentile(nv, 5), null_p95=np.percentile(nv, 95),
                         p_one_sided=float(np.mean(nv <= v))))
    # bootstrap over units / pairs of the observed median change
    rs = np.random.default_rng(7)
    ukeys = list(bm)
    pkeys = list(bt)
    for k in ("jitter", "innov_good"):
        v = np.array([100 * (wm[u][k] / bm[u][k] - 1) for u in ukeys])
        bs = [np.median(rs.choice(v, v.size)) for _ in range(2000)]
        rows.append(dict(metric=k + " (bootstrap over units)", observed_change_pct=float(np.median(v)),
                         null_p05=np.percentile(bs, 2.5), null_p95=np.percentile(bs, 97.5),
                         n=len(v), better=int(np.sum(v < 0))))
    for j, k in enumerate(("twin_rms", "twin_p95")):
        v = np.array([100 * (wt[p][j] / bt[p][j] - 1) for p in pkeys if np.isfinite(bt[p][j])])
        bs = [np.median(rs.choice(v, v.size)) for _ in range(2000)]
        rows.append(dict(metric=k + " (bootstrap over pairs)", observed_change_pct=float(np.median(v)),
                         null_p05=np.percentile(bs, 2.5), null_p95=np.percentile(bs, 97.5),
                         n=len(v), better=int(np.sum(v < 0))))
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / f"u23r_kalman_{set_name}_{method}.csv", index=False)
    pd.set_option("display.width", 220)
    print(f"{len(units)} units, {len(bt)} twin pairs, {n_perm} permutations")
    print(out.round(3).to_string(index=False))
    qv = qc_variants(ok)
    qv.to_csv(C.OUT / f"u23r_kalman_qc_{set_name}_{method}.csv", index=False)
    print("median change vs the equal-weight filter on all nights [%] (negative = better):")
    print(qv.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
