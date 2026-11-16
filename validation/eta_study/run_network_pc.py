"""Home-PC share of the full-network eta-study run: the same extraction + step-2 fit as the CSCS
array (validation/eta_study/cscs/extract_fit_array.sbatch), for a list of units taken off the
Linard queue (see C:/DATA/Projects/202611_eta_study/network/split.csv).

Pool over unit-months (largest unit first); when all months of a unit are done the unit is merged
and fitted (working gates), exactly like `extract_incloud_profiles.py --unit` + `fit_alpha_shape.py --unit`.
Outputs go to a separate tree so nothing mixes with the 43-unit study run:
  C:/DATA/Projects/202611_eta_study/network/pc/incloud/   (<key>.npz, units/<key>.json)
  C:/DATA/Projects/202611_eta_study/network/pc/step2/profiles/<key>_working.parquet

Run:  python validation/eta_study/run_network_pc.py --keys <file> --workers 8
"""
from __future__ import annotations

import os
# set BEFORE importing the study modules (they read the paths at import; spawned workers re-import
# this module first, so they see the same values)
os.environ.setdefault("ETA_OUT", "C:/DATA/Projects/202611_eta_study/network/pc/incloud")
os.environ.setdefault("ETA_BASE", "C:/DATA/Projects/202611_eta_study/network/pc")
os.environ.setdefault("ETA_LUT", "C:/DATA/Projects/202611_eta_study/forward_lut_ext.npz")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import sys
import time as time_mod
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_incloud_profiles as E  # noqa: E402
import fit_alpha_shape as F  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", required=True, help="file with one unit key per line")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    keys = [l.strip() for l in open(a.keys, encoding="utf-8") if l.strip()]
    units = E.select_units(network=True)
    sub = {k: units[k] for k in keys}
    (E.OUT / "units").mkdir(parents=True, exist_ok=True)
    for k, u in sub.items():
        with open(E.OUT / "units" / f"{k}.json", "w", encoding="utf-8") as f:
            json.dump({kk: vv for kk, vv in u.items() if kk != "days"}, f, indent=1)
    js = E.jobs_for(sub)
    js.sort(key=lambda j: (-j[0]["n_cloud_days"], -len(j[2])))           # largest unit first
    remaining = defaultdict(int)
    for u, _, _ in js:
        remaining[u["key"]] += 1
    print(f"{len(sub)} units, {len(js)} unit-months, {sum(len(j[2]) for j in js)} days -> {a.workers} processes",
          flush=True)
    t0 = time_mod.time()
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(E.process_unit_month, j) for j in js]
        for f in as_completed(futs):
            key, ym, msgs, n, dt = f.result()
            done += 1
            remaining[key] -= 1
            if remaining[key] == 0:
                n_m = E.merge_unit(sub[key])
                fit = F.fit_unit((key, sub[key], "working")) if n_m else None
                print(f"[{done}/{len(js)} {time_mod.time() - t0:6.0f} s] {key}: merged {n_m} profiles, "
                      f"fit {fit[2] if fit else 0} valid", flush=True)
    print(f"total {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
