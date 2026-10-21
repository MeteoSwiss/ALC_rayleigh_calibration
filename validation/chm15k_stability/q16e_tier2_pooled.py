# -*- coding: utf-8 -*-
"""q16e check 5b: pooled effect of the eprof_v2.2 fallback tier on the annual cycle (from q16e_tier2_selection.csv).

Per unit, ln C is centred on the unit's tier-1 mean; the network-pooled monthly mean is computed with all
published nights and with tier-1 nights only, and a first harmonic is fitted to each. The difference is
the in-phase contribution of the fallback tier (share x offset) to the annual amplitude.
"""
import numpy as np
import pandas as pd

from common import OUT

d = pd.read_csv(OUT / "q16" / "q16e_tier2_selection.csv")
d["tier"] = np.where(d["c_eprof_v2"] > 0, 1, np.where(d.get("c_eprof_v2.2", np.nan) > 0, 2, 0))
d = d[d.tier > 0].copy()
d["t"] = pd.to_datetime(d.date.astype(str))
d["m"] = d.t.dt.month
d["y"] = np.log(d.c_pub)
d["y"] -= d.groupby("stream").y.transform(lambda s: s[d.loc[s.index, "tier"] == 1].mean())


def harm(m, y):
    w = 2 * np.pi / 12
    A = np.c_[np.ones_like(m), np.cos(w * (m - 0.5)), np.sin(w * (m - 0.5))]
    p, *_ = np.linalg.lstsq(A, y, rcond=None)
    return 100 * np.hypot(p[1], p[2]), (np.arctan2(p[2], p[1]) / w) * 30.44 % 365


m = np.arange(1, 13, dtype=float)
all_ = d.groupby("m").y.mean().reindex(range(1, 13)).values
t1 = d[d.tier == 1].groupby("m").y.mean().reindex(range(1, 13)).values
share = d.groupby("m").tier.apply(lambda s: np.mean(s == 2)).reindex(range(1, 13)).values
for name, y in (("all published", all_), ("tier-1 only", t1), ("difference", all_ - t1)):
    a, p = harm(m, y)
    print(f"{name:14s} amp {a:5.2f} %  peak DOY {p:4.0f}   monthly: " + " ".join(f"{100 * v:+.1f}" for v in y))
print("tier-2 share   " + " ".join(f"{100 * v:.0f}" for v in share))
# leave-one-unit-out spread of the difference amplitude
amps = []
for u in d.stream.unique():
    e = d[d.stream != u]
    ya = e.groupby("m").y.mean().reindex(range(1, 13)).values
    yb = e[e.tier == 1].groupby("m").y.mean().reindex(range(1, 13)).values
    amps.append(harm(m, ya)[0] - harm(m, yb)[0])
print(f"amplitude difference (all - tier1), leave-one-unit-out range {min(amps):+.2f} .. {max(amps):+.2f} %-pt")
