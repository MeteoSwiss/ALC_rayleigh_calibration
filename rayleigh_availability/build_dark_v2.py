# -*- coding: utf-8 -*-
"""Payerne dark profiles v2: the full campaign (10 sessions, May-September 2026) for the CHM15k (A)
and the CL61 (C), using the covered intervals DETECTED in the L1 data (dark_windows_detect.py, already
trimmed inside the hood-on/off edges) instead of the approximate logbook times.

Same processing as dark_profiles.py (time median per session in rcs_0 units, range smoothing, median
over sessions). Written to a NEW file so that runs made with the v1 profile (6 sessions, logbook
windows) stay reproducible. The CL31 (B) is included with its data-confirmed windows (29 Sep: hood
removed 09:39-09:40 UTC, not logged); dark_profiles.pool keeps only its post-optical-module-swap
sessions (after 2026-07-07 13:00), since the swap changed the CL31 gain.

Run:  python rayleigh_availability/build_dark_v2.py
"""
from __future__ import annotations
import json
from datetime import datetime

import numpy as np

import dark_profiles as DP   # run from rayleigh_availability/ (or with it on sys.path)

WINDOWS_JSON = DP.DATA / "dark_windows_detected.json"
OUT = DP.DATA / "dark_profiles_payerne_v2.npz"


def main():
    win = json.loads(WINDOWS_JSON.read_text())
    DP.MARGIN_S = 0.0                                # windows are already trimmed from the data
    DP._parse = lambda s: datetime.strptime(s, "%Y-%m-%d %H:%M:%S")
    out = {}
    for ident, name in (("A", "CHM15k"), ("B", "CL31"), ("C", "CL61")):
        DP.WINDOWS[ident] = [tuple(w) for w in win[ident]]
        print(f"\n== {name} ({ident}): {len(DP.WINDOWS[ident])} sessions")
        rng, per_window = DP.collect(ident)
        pooled = DP.pool(ident, rng, per_window)
        b = pooled["b"]
        band = (rng >= 2500) & (rng <= 6500)
        print(f"   pooled over {pooled['n_windows']} sessions, {pooled['n_profiles']} profiles; "
              f"<b> 2.5-6.5 km = {np.nanmean(b[band]):+.4e}")
        out[f"{ident}_range"] = rng
        out[f"{ident}_b_rcs"] = b
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"{ident}_b"] = np.where(rng > 0, b / rng ** 2, np.nan)
        out[f"{ident}_sem"] = pooled["sem"]
        out[f"{ident}_nwin"] = pooled["n_windows"]
        out[f"{ident}_nprof"] = pooled["n_profiles"]
    np.savez(OUT, **out)
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
