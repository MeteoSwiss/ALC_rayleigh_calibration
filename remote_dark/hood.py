# -*- coding: utf-8 -*-
"""Payerne hood sessions: the ground truth, in two forms.

1. `truth(ident)` — the pooled median dark profile b(z) from `dark_profiles_payerne.npz`
   (what M1's retrieved profile is judged against).
2. `session_frames(ident)` — the TIME-RESOLVED dark of each logbook session, with the
   housekeeping modulators alongside. This is what M3 is judged against: the two multi-hour
   sessions (CL31 26–27 May ~25 h and 09–10 Jul ~21 h) span full diurnal cycles UNDER the hood,
   so the dark's own response to internal temperature and daylight background is in the data.

The logbook (window list per instrument) is imported from rayleigh_availability.dark_profiles —
one source of truth for what was measured when, including the CL31 optic-block swap of 07 Jul
13:00 that splits the CL31 record into two different instruments.
"""
from __future__ import annotations

from datetime import timedelta

import numpy as np

from remote_dark.common import HOOD_NPZ, PAYERNE, read_day
from rayleigh_availability.dark_profiles import CL31_SWAP, WINDOWS, _parse

#: Kotthaus 2016 (applied by Looschelders 2025): the first ~20 min of a hood session show a
#: near-range transient while the sensor's compensation settles. Discarded everywhere.
SETTLE_MIN = 20.0


def truth(ident: str):
    """(rng, b, sem, b_raw) — the pooled hood dark for one Payerne instrument, rcs_0 units.

    `b` is the campaign-smoothed profile (median 150 m + mean 450 m, the subtractable product);
    `b_raw` is the unsmoothed pooled median — the one to fit parametric FAMILIES against, since
    smoothing at 450 m would erase exactly the resonant structure the CL31 family models.
    """
    z = np.load(HOOD_NPZ, allow_pickle=True)
    return (np.asarray(z[f"{ident}_range"], float), np.asarray(z[f"{ident}_b"], float),
            np.asarray(z[f"{ident}_sem"], float), np.asarray(z[f"{ident}_b_rcs"], float))


def session_frames(ident: str, min_hours: float = 0.0):
    """Yield one dict per logbook session: time-resolved dark + modulators.

    dark(t, z) stays in rcs_0 units. Sessions shorter than *min_hours* are skipped (M3 wants the
    multi-hour ones; M1's pooled truth already integrates the short ones).
    """
    wmo = PAYERNE["wmo"]
    for t0s, t1s in WINDOWS[ident]:
        t0, t1 = _parse(t0s), _parse(t1s)
        if (t1 - t0).total_seconds() / 3600.0 < min_hours:
            continue
        t0e = t0 + timedelta(minutes=SETTLE_MIN)
        frames, hks, times = [], {k: [] for k in ("bckgrd", "t_int", "laser", "window")}, []
        rng = None
        day = t0.replace(hour=0, minute=0, second=0, microsecond=0)
        while day <= t1:
            d = read_day(wmo, ident, day)
            day += timedelta(days=1)
            if d is None:
                continue
            m = (d["times"] >= t0e) & (d["times"] <= t1)
            if not m.any():
                continue
            rng = d["rng"]
            frames.append(d["rcs"][m])
            times.append(d["times"][m])
            for k in hks:
                hks[k].append(d["hk"][k][m])
        if not frames:
            continue
        yield {
            "ident": ident, "t0": t0, "t1": t1,
            "era": ("pre_swap" if (ident == "B" and t1 <= CL31_SWAP) else
                    "post_swap" if ident == "B" else "single"),
            "times": np.concatenate(times),
            "rng": rng,
            "dark": np.vstack(frames),
            "hk": {k: np.concatenate(v) for k, v in hks.items()},
        }


__all__ = ["truth", "session_frames", "WINDOWS", "CL31_SWAP", "SETTLE_MIN"]
