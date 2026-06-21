"""
run_network_v2_vs_v11.py — network-wide comparison of the optimized E-PROF v2 (config C8, now the
eprof_v2 default) against E-PROF v1.1, on EVERY CHM15k / CL61 / Mini-MPL stream in the 2026 archive,
at BOTH L1 and L2. Per clear (fit-reaching) night we record, for each method, whether it yields a
valid calibration (window passes the pipeline QC rel_error<=15%) and the resulting lidar constant.

One load per night (calibrate_rayleigh -> fit_inputs); both methods are then evaluated on that one
prepared profile. Output: figs_paper_validation/network_v2_v11/net_<level>_<label>.json
    = {ds: {"v2": [ok,cl,rel], "v1.1": [ok,cl,rel]}}.
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
from calibration.rayleigh.molecular_methods import select_molecular_window

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/network_v2_v11")
OUT.mkdir(parents=True, exist_ok=True)
HALF = tuple(range(250, 2000, 240))
QC = 15.0
METHODS = {"v2": "eprof_v2", "v1.1": "eprof_v1.1"}   # v2 uses the optimized C8 defaults

ROOTS = {"L1": (Path("D:/E-PROFILE_L1_2026"), DataLevel.L1),
         "L2": (Path("D:/E-PROFILE_L2_2026"), DataLevel.L2_DAILY)}


def base_options(level):
    root, dl = ROOTS[level]
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = root
    o.data_level = dl
    o.apply_wv_correction = True
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def date_strs(first, last):
    d0, d1 = datetime.strptime(first, "%Y%m%d"), datetime.strptime(last, "%Y%m%d")
    step = int(os.environ.get("RC_DAYSTEP", "1"))   # >1 = sample every Nth day (re-run speedup)
    d = d0
    while d <= d1:
        yield d.strftime("%Y%m%d")
        d += timedelta(days=step)


def eval_night(signal, p_mol, rng, stack):
    out = {}
    for tag, method in METHODS.items():
        try:
            w = select_molecular_window(method, signal, p_mol, rng, HALF, signal_stack=stack)
            ok = bool(w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC
                      and np.isfinite(w.cl) and w.cl > 0)
            out[tag] = [ok, float(w.cl), float(w.rel_error)]
        except Exception:
            out[tag] = [False, float("nan"), float("nan")]
    return out


def run_one(args):
    level, inst = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    o = base_options(level)
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    per_night = {}
    for ds in date_strs(inst["first"], inst["last"]):
        fin = {}
        try:
            calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
        except Exception:
            continue
        if not fin:
            continue
        try:
            per_night[ds] = eval_night(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
        except Exception:
            continue
    (OUT / f"net_{level}_{inst['label']}.json").write_text(json.dumps(per_night), encoding="utf-8")
    n2 = sum(v["v2"][0] for v in per_night.values())
    n1 = sum(v["v1.1"][0] for v in per_night.values())
    return level, inst["label"], inst["group"], len(per_night), n2, n1


def run_slice(level, slice_i, n_slices):
    """Single-process: handle MANIFEST[slice_i::n_slices] for one level (Windows-safe — no
    ProcessPoolExecutor orphans; parallelism comes from launching N slice processes)."""
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.ERROR)
    insts = MANIFEST[slice_i::n_slices]
    print(f"slice {slice_i}/{n_slices} {level}: {len(insts)} streams", flush=True)
    for k, inst in enumerate(insts):
        try:
            _, label, group, n_fit, n2, n1 = run_one((level, inst))
            print(f"  [{slice_i}:{k+1}/{len(insts)}] {level} {label:22s} ({group:8s}) "
                  f"fit={n_fit:3d} v2={n2:3d} v11={n1:3d}", flush=True)
        except Exception as e:
            print(f"  [{slice_i}] FAILED {inst['label']}: {e}", flush=True)
    print(f"NETSLICE_{slice_i}_DONE", flush=True)


def main():
    # slice mode: <level> <slice_i> <n_slices>  (single-process, shell-parallel)
    if len(sys.argv) >= 4:
        run_slice(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
        return
    only_level = sys.argv[1] if len(sys.argv) > 1 else None
    jobs = [(lvl, inst) for lvl in ROOTS for inst in MANIFEST
            if only_level is None or lvl == only_level]
    workers = int(os.environ.get("RC_WORKERS", "14"))
    print(f"network v2-vs-v1.1: {len(jobs)} (level x stream) jobs, {workers} workers")
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            try:
                level, label, group, n_fit, n2, n1 = fut.result()
                done += 1
                if done % 20 == 0 or n_fit:
                    print(f"  [{done}/{len(jobs)}] {level} {label:22s} ({group:8s}) "
                          f"{n_fit:3d} fit-nights | v2 {n2:3d} valid, v1.1 {n1:3d} valid")
            except Exception as e:
                print(f"  FAILED: {e}")
    print("NETWORK_DONE")


if __name__ == "__main__":
    main()
