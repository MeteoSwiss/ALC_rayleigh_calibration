# -*- coding: utf-8 -*-
"""Is the night-to-night variance shared by co-located CHM15k twins a RETRIEVAL error (atmosphere:
aerosol in/below the window, molecular model) or a real CHM15k GAIN change?

Test: the 31-day-running-median anomalies of ln C_L (kappa-corrected) of co-located units, from the
published (all-night) series. Twins (same type, same wavelength) share both the atmosphere and any
CHM15k-specific environmental gain forcing; a co-located CL61 (Rayleigh, 910 nm) shares the
atmosphere but not the CHM15k gain. If the CHM15k-CL61 anomalies correlate as much as the twins, the
common part is atmospheric (a retrieval error that belongs in the uncertainty); if they do not, it is
a CHM15k gain variation (real, not a calibration error).
Output: C:/DATA/Projects/202610_v23/uncertainty/u23_common_mode.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from u23_analysis import (OUT, Published, census, kappa_table, running_median,  # noqa: E402
                          unit_kappa, pull_stats)

PAIRS = [
    # twins (CHM15k / CHM15k)
    ("0-20000-0-06240_A", "0-20000-0-06240_B"), ("0-20000-0-06240_A", "0-20000-0-06240_C"),
    ("0-20000-0-06240_A", "0-20000-0-06240_D"), ("0-20000-0-06240_B", "0-20000-0-06240_C"),
    ("0-20000-0-06240_B", "0-20000-0-06240_D"), ("0-20000-0-06240_C", "0-20000-0-06240_D"),
    ("0-20000-0-12843_A", "0-20000-0-12843_B"), ("0-20008-0-INO_A", "0-20008-0-INO_B"),
    # CHM15k / CL61, same site
    ("0-20000-0-10393_0", "0-20000-0-10393_C"), ("0-380-5-1_0", "0-380-5-1_B"),
    ("0-20000-0-06610_A", "0-20000-0-06610_C"),
]


def anomalies(P, kappa):
    lncc = P.corrected(kappa)
    base = running_median(P.t, lncc, exclude_self=True)
    return dict(zip(P.dates, lncc - base))


def main():
    import pandas as pd
    cen = census()
    ktab = kappa_table()
    cache = {}
    rows = []
    for a, b in PAIRS:
        an = {}
        for k in (a, b):
            if k not in cache:
                wmo, ident = k.rsplit("_", 1)
                s = cen[(wmo, ident)]
                P = Published(wmo, ident, s["lon"])
                cache[k] = (anomalies(P, unit_kappa(k, s["type"], ktab)) if P.ok else {}, s)
            an[k] = cache[k][0]
        common = sorted(set(an[a]) & set(an[b]))
        x = np.array([an[a][d] for d in common])
        y = np.array([an[b][d] for d in common])
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        if x.size < 20:
            continue
        # robust correlation: clip at 4 robust sd, then Pearson
        def clip(v):
            s = pull_stats(v)["rstd"]
            return np.abs(v - np.median(v)) <= 4 * s
        m = clip(x) & clip(y)
        r = float(np.corrcoef(x[m], y[m])[0, 1])
        sa, sb = pull_stats(x)["rstd"], pull_stats(y)["rstd"]
        sd = pull_stats(x - y)["rstd"]
        # shared variance estimate: cov / (sa sb) with robust scales
        rows.append(dict(a=a, b=b, type_a=cache[a][1]["type"], type_b=cache[b][1]["type"],
                         n=int(x.size), sigma_a=100 * sa, sigma_b=100 * sb,
                         sigma_ratio=100 * sd, corr=r,
                         unit_specific_a=100 * sa * np.sqrt(max(1 - r, 0.0)),
                         shared_sigma=100 * np.sqrt(max(r, 0.0) * sa * sb)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "u23_common_mode.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
