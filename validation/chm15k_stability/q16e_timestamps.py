# -*- coding: utf-8 -*-
"""q16e check 2: how the published nightly Rayleigh constant is time-stamped.

output.write_calibration_result stamps a night with floor(max(profile time)) and overwrites any
same-method row with the same floor. A night whose usable profiles all precede 00 UTC is therefore
stamped D-1 (and can overwrite night D-1). Count, per month, the share of CHM15k nights whose
calibration period ends before 00 UTC, plus the night-centre local solar hour, to see whether
the stamping or the averaged period is seasonal.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import OUT, chm15k_streams, days_to_datetime, load_constants  # noqa: E402


def main():
    rows = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        c = load_constants(wmo, ident)
        if c is None:
            continue
        ok = np.isfinite(c["t"]) & (c["c"] > 0)
        t, t0, t1 = c["t"][ok], c["t0"][ok], c["t1"][ok]
        mon = np.array([d.month for d in days_to_datetime(t)])
        evening_only = np.floor(t1) == np.floor(t0)          # whole period before 00 UTC (or after)
        stamp_prev = np.floor(t1) < np.floor(t0 + 0.5)        # stamped on the evening day
        mid_solar = np.mod(((t0 + t1) / 2 % 1) * 24 + lon / 15, 24)
        dur_h = (t1 - t0) * 24
        dup = len(t) - len(np.unique(np.floor(t)))
        for m in range(1, 13):
            s = mon == m
            if s.any():
                rows.append((m, s.sum(), evening_only[s].sum(), stamp_prev[s].sum(),
                             np.median(np.where(mid_solar[s] > 12, mid_solar[s] - 24, mid_solar[s])),
                             np.median(dur_h[s]), dup))
    R = np.array(rows, float)
    print("month  n_nights  evening_only[%]  stamped_on_evening_day[%]  median_mid_solar_h  median_period_h")
    for m in range(1, 13):
        s = R[:, 0] == m
        n = R[s, 1].sum()
        print(f"{m:5d} {n:9.0f} {100 * R[s, 2].sum() / n:15.1f} {100 * R[s, 3].sum() / n:25.1f} "
              f"{np.median(R[s, 4]):19.2f} {np.median(R[s, 5]):16.2f}")


if __name__ == "__main__":
    main()
