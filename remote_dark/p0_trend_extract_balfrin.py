# -*- coding: utf-8 -*-
"""Runs ON balfrin — per-stream noise-floor level + robust trend from the P0 caches
(the second triage axis: rising floor = ageing laser = 1/laser-inflated bound)."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np

OUT = Path("/scratch/mch/mhrvo/remote_dark_net")

rows = []
for f in sorted(OUT.glob("0-*.npz")):
    try:
        z = np.load(f, allow_pickle=True)
        sig, n_prof, rng, dates = z["sigma"], z["n_prof"], z["rng"], z["dates"]
    except Exception:                                                       # noqa: BLE001
        continue
    z2 = np.where(rng > 0, rng, np.nan) ** 2
    top = rng >= np.nanpercentile(rng, 92)
    fl = np.nanmedian((sig * np.sqrt(np.maximum(n_prof.astype(float), 1.0))[:, None]
                       / z2[None, :])[:, top], axis=1)
    ok = np.isfinite(fl) & (fl > 0)
    if ok.sum() < 20:
        continue
    t = np.array([(datetime.strptime(str(int(d)), "%Y%m%d")
                   - datetime(2025, 1, 1)).days for d in dates[ok]], float) / 365.25
    y = np.log(fl[ok])
    A = np.vstack([t, np.ones_like(t)]).T
    beta = np.linalg.lstsq(A, y, rcond=None)[0]
    rows.append(f"{f.stem},{int(ok.sum())},{np.median(fl[ok]):.4e},"
                f"{100 * (np.exp(beta[0]) - 1):.1f}")

with open(OUT / "p0_floor_trends.csv", "w") as fh:
    fh.write("key,n_nights,floor_med,trend_pct_per_yr\n")
    fh.write("\n".join(rows) + "\n")
print(f"TRENDS_DONE {len(rows)}")
