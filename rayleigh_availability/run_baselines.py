# -*- coding: utf-8 -*-
"""Phase 0 — re-run the v1.1 and v2(C8) baselines on the availability corpus.

The existing per-night archive CANNOT be used as a baseline: it mixes two code vintages at
20260601 (the -2 message text and the -9 frequency both change there), so any "v2 today" number
read from it is a blend of two builds. Everything is therefore recomputed here with one build.

This runs the FULL pipeline per night per method (not just the window selection) because the
baseline of interest is the real flag histogram -- v1.1 and v2 differ in the post-fit gates
(-3 method disagreement, -6 uncertainty, -9 lower-signal layer) as well as in window eligibility.

No diagnostic plots are produced (plot_main/plot_all off): the run is for numbers only, and the
PNGs dominate the wall clock otherwise.

Out: <DATA>/rayleigh_availability/baselines/base_<method>_<label>.json
     = {date: [flag, lidar_constant, uncertainty, bottom_height, top_height, message]}
       (the fields of calibration.CalibrationResult)
Resumable: an existing output file is skipped unless --force.
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

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel  # noqa: E402
from calibration.config import InstrumentType  # noqa: E402

MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/baselines")
OUT.mkdir(parents=True, exist_ok=True)
L1_ROOT = Path("D:/E-PROFILE_L1_2026")

METHODS = ("eprof_v1.1", "eprof_v2")


def base_options(method):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = method
    o.molecular_source = "standard"     # CHM15k is 1064 nm; keeps CAMS out of the baseline
    o.apply_wv_correction = True
    o.plot_main = False                 # numbers only -- no diagnostics
    o.plot_all = False
    o.folder_output = OUT
    return o


def dates_of(inst):
    for first, last in inst["periods"]:
        d = datetime.strptime(first, "%Y%m%d")
        d1 = datetime.strptime(last, "%Y%m%d")
        while d <= d1:
            yield d.strftime("%Y%m%d")
            d += timedelta(days=1)


def run_one(args):
    method, inst = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    out_f = OUT / f"base_{method}_{inst['label']}.json"
    # Date-level resume: keep the nights already computed and process only the ones the corpus has
    # gained. Re-running a whole stream to add a few months would otherwise dominate the cost.
    rec = {}
    if out_f.exists() and "--force" not in sys.argv:
        rec = json.loads(out_f.read_text())
        todo = [d for d in dates_of(inst) if d not in rec]
        if not todo:
            return method, inst["label"], len(rec), _n_ok(rec), True

    o = base_options(method)
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    for ds in dates_of(inst):
        if ds in rec:
            continue
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception as exc:                       # keep the night, record the failure
            rec[ds] = [-99, None, None, None, None, f"EXC {type(exc).__name__}: {exc}"[:120]]
            continue
        if r is None:
            continue
        rec[ds] = [_f(r.flag), _f(r.lidar_constant), _f(r.uncertainty),
                   _f(r.calibration_bottom_height), _f(r.calibration_top_height),
                   str(r.message or "")[:120]]
    out_f.write_text(json.dumps(rec), encoding="utf-8")
    return method, inst["label"], len(rec), _n_ok(rec), False


def _f(v):
    try:
        if v is None:
            return None
        f = float(v)
        return f if f == f else None                   # NaN -> None
    except (TypeError, ValueError):
        return None


def _n_ok(rec):
    return sum(1 for v in rec.values() if v[0] in (1.0, 0.5))


def main():
    only = [a for a in sys.argv[1:] if not a.startswith("--")]
    insts = [i for i in MANIFEST if not only or i["label"] in only or i["group"] in only]
    jobs = [(m, i) for m in METHODS for i in insts]
    workers = int(os.environ.get("RA_WORKERS", str(max(1, (os.cpu_count() or 8) - 2))))
    print(f"baselines: {len(jobs)} (method x stream) jobs over {workers} workers; "
          f"{len(insts)} streams, methods {METHODS}", flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            j = futs[fut]
            done += 1
            try:
                method, label, n, n_ok, cached = fut.result()
                print(f"  [{done}/{len(jobs)}] {method:11s} {label:28s} {n:4d} nights, "
                      f"{n_ok:3d} valid{' (cached)' if cached else ''}", flush=True)
            except Exception as e:
                print(f"  [{done}/{len(jobs)}] FAILED {j[0]} {j[1]['label']}: {e}", flush=True)
    print("BASELINES_DONE", flush=True)


if __name__ == "__main__":
    main()
