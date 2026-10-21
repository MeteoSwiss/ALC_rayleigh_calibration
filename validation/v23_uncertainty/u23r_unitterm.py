# -*- coding: utf-8 -*-
"""Review attack: is something structurally missing for units like the Mini-MPL or Schleiz?

Candidate: a unit-level repeatability term estimated CAUSALLY from the unit's own recent history,
  sigma_rep(t)^2 = max(0, rstd(anomalies of the unit's nights in [t-120 d, t-1 d])^2
                        - median(budget sigma^2 of those nights)),
added in quadrature to the nightly budget. It is evaluated only on nights with >= 10 past nights,
so every number is out of sample in time. Reported: pooled and per-unit pulls before/after (time-series
truth), the per-unit spread, and the twin pull (where the term must NOT be needed: it would then
over-cover, because it also absorbs real gain changes common to co-located units).

Usage: python u23r_unitterm.py [orig|fresh] [eprof_v2.2|eprof_v2.3]
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
PAST_DAYS = 120
MIN_PAST = 10


def add_rep(d, sig_col="sig23"):
    d = d.sort_values(["key", "tday"]).copy()
    d["sig_rep"] = np.nan
    d["n_past"] = 0
    for key, g in d.groupby("key"):
        t = g.tday.values.astype(float)
        a = g.anom.values
        s = g[sig_col].values
        rep = np.full(len(g), np.nan)
        npast = np.zeros(len(g), int)
        for i in range(len(g)):
            m = (t >= t[i] - PAST_DAYS) & (t <= t[i] - 1) & np.isfinite(a) & np.isfinite(s)
            npast[i] = int(m.sum())
            if m.sum() >= MIN_PAST:
                so = C.rstd(a[m])
                sb2 = np.median(s[m] ** 2)
                rep[i] = np.sqrt(max(so ** 2 - sb2, 0.0))
        d.loc[g.index, "sig_rep"] = rep
        d.loc[g.index, "n_past"] = npast
    return d


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    p = prepare(set_name, method)
    df = p["df"]
    d = df[df.ok & np.isfinite(df.anom) & np.isfinite(df.sig23)].copy()
    d = add_rep(d)
    e = d[np.isfinite(d.sig_rep)].copy()
    e["sig_new"] = np.hypot(e.sig23, e.sig_rep)
    rows = []
    for lab, m in (("CHM15k (Payerne excl.)", (e.type == "CHM15k") & (e.key != C.PAYERNE)),
                   ("Payerne A", e.key == C.PAYERNE), ("CL61", e.type == "CL61"),
                   ("Mini-MPL", e.type == "Mini-MPL"), ("Schleiz", e.key == C.SCHLEIZ)):
        g = e[m]
        if len(g) < 10:
            continue
        a, b = C.pull_stats(g.anom / g.sig23), C.pull_stats(g.anom / g.sig_new)
        per_a = g.groupby("key").apply(lambda x: C.pull_stats(x.anom / x.sig23)["rstd"])
        per_b = g.groupby("key").apply(lambda x: C.pull_stats(x.anom / x.sig_new)["rstd"])
        rows.append(dict(population=lab, n=a["n"], units=g.key.nunique(),
                         pull_budget=a["rstd"], out2_budget=a["out2"], pull_with_rep=b["rstd"],
                         out2_with_rep=b["out2"], unit_pull_sd_budget=float(np.std(np.log(per_a))),
                         unit_pull_sd_with_rep=float(np.std(np.log(per_b))),
                         med_sig=100 * g.sig23.median(), med_sig_rep=100 * g.sig_rep.median()))
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False))
    out.to_csv(C.OUT / f"u23r_unitterm_{set_name}_{method}.csv", index=False)
    tw = p.get("tw")
    if tw is not None and len(tw):
        s_new = e.set_index(["key", "date"]).sig_new
        s_old = d.set_index(["key", "date"])
        nomol = np.sqrt(np.clip(s_old.sig23 ** 2 - s_old.def_sig_mol ** 2 - s_old.def_sig_aer ** 2, 0, None))
        rep = e.set_index(["key", "date"]).sig_rep
        za, zb = [], []
        for r in tw.itertuples():
            ka, kb = (r.ka, r.date), (r.kb, r.date)
            if ka in rep.index and kb in rep.index:
                sa, sb = nomol.get(ka, np.nan), nomol.get(kb, np.nan)
                za.append(r.anom / np.hypot(sa, sb))
                zb.append(r.anom / np.hypot(np.hypot(sa, rep[ka]), np.hypot(sb, rep[kb])))
        a, b = C.pull_stats(np.array(za)), C.pull_stats(np.array(zb))
        print(f"twins ({a['n']} pair-nights with a causal term for both units): pull {a['rstd']:.3f} "
              f"-> {b['rstd']:.3f}, outside 2 sigma {100 * a['out2']:.1f} -> {100 * b['out2']:.1f} %")


if __name__ == "__main__":
    main()
