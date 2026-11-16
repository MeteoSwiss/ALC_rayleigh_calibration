# -*- coding: utf-8 -*-
"""Review attack (claim 3): does the budget rank the units, and is 0.58 -> 0.67 a real gain?

* Spearman(reported unit median sigma, observed unit scatter) with bootstrap CIs (nights resampled
  within units, then units resampled) and the CI of the v2.3 - v2.2 difference;
* a benchmark predictor that needs no budget at all: the unit's own observed scatter in the OTHER
  year (2025 predicts 2026 and vice versa) -- i.e. a unit-level repeatability term;
* out of sample: the same on the fresh units.

Usage: python u23r_ranking.py [orig|fresh] [eprof_v2.2|eprof_v2.3]
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


def unit_stats(d):
    import pandas as pd
    rows = []
    for key, g in d.groupby("key"):
        if len(g) < 8:
            continue
        rows.append(dict(key=key, site=g.site.iloc[0], type=g.type.iloc[0], n=len(g),
                         obs=C.rstd(g.anom), rep22=g.sig_ens.median(), rep23=g.sig23.median(),
                         gabs=np.abs(g.def_g).median()))
    return pd.DataFrame(rows)


def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y)[0])


def boot_rho(d, n_boot=1000, seed=3):
    """Two-stage bootstrap: units resampled, nights resampled within each unit."""
    rs = np.random.default_rng(seed)
    groups = {k: g for k, g in d.groupby("key") if len(g) >= 8}
    keys = list(groups)
    out = []
    for _ in range(n_boot):
        ks = rs.choice(keys, len(keys))
        obs, r22, r23 = [], [], []
        for k in ks:
            g = groups[k]
            ii = rs.integers(0, len(g), len(g))
            gg = g.iloc[ii]
            obs.append(C.rstd(gg.anom))
            r22.append(gg.sig_ens.median())
            r23.append(gg.sig23.median())
        a, b = spearman(r22, obs), spearman(r23, obs)
        out.append((a, b, b - a))
    return np.array(out)


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    df = prepare(set_name, method)["df"]
    d = df[df.ok & np.isfinite(df.anom) & np.isfinite(df.sig23)].copy()
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    for label, sub in (("all types", d), ("CHM15k", d[d.type == "CHM15k"])):
        ut = unit_stats(sub)
        r22, r23 = spearman(ut.rep22, ut.obs), spearman(ut.rep23, ut.obs)
        bs = boot_rho(sub)
        say(f"[{set_name}/{method}, {label}] {len(ut)} units: Spearman v2.2 {r22:+.2f} "
            f"[{np.percentile(bs[:, 0], 2.5):+.2f}, {np.percentile(bs[:, 0], 97.5):+.2f}], v2.3 "
            f"{r23:+.2f} [{np.percentile(bs[:, 1], 2.5):+.2f}, {np.percentile(bs[:, 1], 97.5):+.2f}], "
            f"difference {r23 - r22:+.2f} [{np.percentile(bs[:, 2], 2.5):+.2f}, "
            f"{np.percentile(bs[:, 2], 97.5):+.2f}], P(diff <= 0) = {np.mean(bs[:, 2] <= 0):.3f}; "
            f"median |g| alone {spearman(ut.gabs, ut.obs):+.2f}")
        # own history: observed scatter in the other year
        sub = sub.copy()
        sub["year"] = [int(x[:4]) for x in sub.date]
        rows = []
        for key, g in sub.groupby("key"):
            a, b = g[g.year == 2025], g[g.year == 2026]
            if len(a) >= 8 and len(b) >= 8:
                rows.append(dict(key=key, obs25=C.rstd(a.anom), obs26=C.rstd(b.anom),
                                 rep23_25=a.sig23.median(), rep23_26=b.sig23.median(),
                                 rep22_25=a.sig_ens.median(), rep22_26=b.sig_ens.median()))
        h = pd.DataFrame(rows)
        if len(h) >= 6:
            say(f"   year split ({len(h)} units with >= 8 nights in both years): predict 2026 scatter "
                f"from 2025 scatter {spearman(h.obs25, h.obs26):+.2f}, from v2.3 2026 budget "
                f"{spearman(h.rep23_26, h.obs26):+.2f}, from v2.2 2026 {spearman(h.rep22_26, h.obs26):+.2f}; "
                f"predict 2025 from 2026 scatter {spearman(h.obs26, h.obs25):+.2f}, from v2.3 2025 "
                f"budget {spearman(h.rep23_25, h.obs25):+.2f}")
        if label == "all types":
            ut["rank23"] = ut.rep23.rank(ascending=False).astype(int)
            ut["rank_obs"] = ut.obs.rank(ascending=False).astype(int)
            ut.sort_values("rank23").to_csv(C.OUT / f"u23r_ranking_{set_name}_{method}.csv",
                                            index=False)
            top = ut.sort_values("rank23").head(8)
            say("   top-8 by v2.3 reported sigma (rank v2.3 / observed): " + "; ".join(
                f"{r.site[:14]} {r.key.rsplit('_', 1)[1]} {r.rank23}/{r.rank_obs}" for r in top.itertuples()))
            topo = ut.sort_values("rank_obs").head(8)
            say("   top-8 by observed scatter (rank observed / v2.3): " + "; ".join(
                f"{r.site[:14]} {r.key.rsplit('_', 1)[1]} {r.rank_obs}/{r.rank23}" for r in topo.itertuples()))
    (C.OUT / f"u23r_ranking_{set_name}_{method}.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
