# -*- coding: utf-8 -*-
"""Build the Payerne CL31 (B) night cache with the co-located CHM15k's clear-night list.

Why borrowed nights. The v4 night selection takes clear nights from the stream's own Rayleigh
successes — and the CL31 has none in the archive (it is not a Rayleigh-calibrated type), which is
why v4 recorded Payerne B as untestable. But "the night was clear at Payerne" is a property of the
SKY, not of the instrument that certified it: the CHM15k ten metres away is the certifier, and its
218 clear nights apply to the CL31 verbatim.

Era discipline: the optic-block swap of 2026-07-07 ~13:00 changes the instrument (gain ×2.255);
this builds the PRE-swap cache (20250101..20260707), the era with six hood sessions to verify
against. The output uses the exact v4 cache schema so m1_fit consumes it unchanged.
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from rayleigh_availability import dark_from_clearsky as v4              # noqa: E402
from remote_dark.common import V4_DIR                                   # noqa: E402

WMO = "0-20000-0-06610"
ERA = ("20250101", "20260707")           # pre-swap; the hood pool of this era has 6 sessions


def main():
    dates = v4.clear_nights_from_csv(v4.DEFAULT_CSV, WMO, "A", *ERA)
    print(f"CHM15k clear-night list: {len(dates)} nights {ERA[0]}..{ERA[1]}")
    jobs = [dict(date=d, wmo=WMO, ident="B", l1_root=str(v4.DEFAULT_L1),
                 cams_folder=v4.DEFAULT_CAMS) for d in dates]
    nights = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for i, res in enumerate(ex.map(v4.process_night, jobs, chunksize=2)):
            if res is not None:
                nights.append(res)
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(jobs)} processed, {len(nights)} usable", flush=True)
    if not nights:
        print("no usable night — nothing written")
        return
    rng = nights[0]["rng"]
    keep = [n for n in nights if n["rng"].size == rng.size]
    out = V4_DIR / "cache" / f"{WMO}_B_borrowA_{ERA[0]}_{ERA[1]}_220.npz"
    np.savez(out,
             S=np.vstack([n["S"] for n in keep]),
             sigma=np.vstack([n["sigma"] for n in keep]),
             M=np.vstack([n["m"] for n in keep]),
             Aer=np.vstack([n["a"] for n in keep]),
             dates=np.array([int(n["date"]) for n in keep]),
             n_prof=np.array([n["n_prof"] for n in keep]),
             rng=rng, itype=np.array("CL31"), wl_nm=np.array(910.0))
    print(f"{len(keep)} nights -> {out}")


if __name__ == "__main__":
    main()
