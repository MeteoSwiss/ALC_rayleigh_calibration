# -*- coding: utf-8 -*-
"""Runs ON balfrin: extend the Schiphol night caches over Jun–Aug 2026.

Why: the local L1 mirror ends 2026-07-12 and unit D swapped its optical module
(TUB150037 -> TUB160055) on 2026-07-11 — one local night of the new module. The post-swap
nights live only in the balfrin archive. This produces small per-unit extension npz files
(same schema as the v4 cache) to be pulled back and merged locally.
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "alc_v22_code"))

from rayleigh_availability import dark_from_clearsky as v4          # noqa: E402

WMO = "0-20000-0-06240"
L1 = "/scratch/mch/mhrvo/E-PROFILE_L1"
CAMS = "/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS"
CSV = Path("/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel")
OUT = Path("/scratch/mch/mhrvo/remote_dark_ext")
ERA = ("20260601", "20260813")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for ident in "ABCD":
        dates = v4.clear_nights_from_csv(CSV, WMO, ident, *ERA)
        print(f"{ident}: {len(dates)} clear nights {ERA[0]}..{ERA[1]}", flush=True)
        jobs = [dict(date=d, wmo=WMO, ident=ident, l1_root=L1, cams_folder=CAMS)
                for d in dates]
        nights = []
        with ProcessPoolExecutor(max_workers=8) as ex:
            for res in ex.map(v4.process_night, jobs, chunksize=1):
                if res is not None:
                    nights.append(res)
        if not nights:
            print(f"  {ident}: nothing usable", flush=True)
            continue
        rng = nights[0]["rng"]
        keep = [n for n in nights if n["rng"].size == rng.size]
        np.savez(OUT / f"{WMO}_{ident}_ext_{ERA[0]}_{ERA[1]}.npz",
                 S=np.vstack([n["S"] for n in keep]),
                 sigma=np.vstack([n["sigma"] for n in keep]),
                 M=np.vstack([n["m"] for n in keep]),
                 Aer=np.vstack([n["a"] for n in keep]),
                 dates=np.array([int(n["date"]) for n in keep]),
                 n_prof=np.array([n["n_prof"] for n in keep]),
                 rng=rng, itype=np.array("CHM15k"), wl_nm=np.array(1064.0))
        post = sum(1 for n in keep if int(n["date"]) >= 20260711)
        print(f"  {ident}: {len(keep)} nights cached ({post} post-swap-date)", flush=True)
    print("EXT_DONE", flush=True)


if __name__ == "__main__":
    main()
