# -*- coding: utf-8 -*-
"""Phase 2 — run candidate eprof_v2.2 configurations over the availability corpus.

Same harness shape as run_baselines.py (full pipeline per night, so the recorded flag is the real
QC outcome, not just a window-selection verdict), but sweeping the noise-aware gate parameters.

eprof_v2.2 is a strict SUPERSET of eprof_v2 by construction: the noise-relative tier is only
consulted when the strict v2 gates leave no eligible window at all. A candidate can therefore only
ADD nights -- never move or lose one -- which is what keeps the continuity indicator at zero and
confines the risk to "is each ADDED night trustworthy?".

No plots.

Out: <DATA>/rayleigh_availability/candidates/cand_<config>_<label>.json
     = {date: [flag, lidar_constant, uncertainty, bottom, top, message]}
Resumable; --force to recompute.
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
from calibration.io.classification import (  # noqa: E402
    classification_files_for_night, read_classification_curtain)

MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
CLASSIF = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/classification")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/candidates")
OUT.mkdir(parents=True, exist_ok=True)
L1_ROOT = Path("D:/E-PROFILE_L1_2026")

# Candidate gate settings. N* sweeps the chi-square tolerance: how far a window may depart from a
# pure Rayleigh shape, in units of its OWN measured noise. 1.5 is strict (a window must be almost
# perfectly noise-consistent), 4.0 is permissive. R* varies the de-biased scattering reference.
CONFIGS = {
    "N1.5":     dict(max_chi2red=1.5),
    "N2.0":     dict(max_chi2red=2.0),
    "N2.5":     dict(max_chi2red=2.5),      # the registered eprof_v2.2 default
    "N3.0":     dict(max_chi2red=3.0),
    "N4.0":     dict(max_chi2red=4.0),
    # reference-percentile variants at the chosen tolerance
    "N2.5_p05": dict(max_chi2red=2.5, ref_pct=5.0),
    "N2.5_p20": dict(max_chi2red=2.5, ref_pct=20.0),
    # keep the raw (min-based) scattering reference: isolates how much the de-biasing contributes
    "N2.5_rawref": dict(max_chi2red=2.5, ref_chi2max=-1.0),
}

# Phase 3: the same gates PLUS the target-classification pre-fit cell mask (contaminated cells,
# aerosol included, removed before the window search). Only the streams in run_classification.SUBSET
# have a classification, so these configs are only meaningful there. "Mv2" applies the mask to the
# CURRENT operational gates -- the mask's effect on its own, with no noise-relative gate involved.
MASK_CONFIGS = {
    "Mv2":   ("eprof_v2",   {}),
    "M2.5":  ("eprof_v2.2", dict(max_chi2red=2.5)),
}
CONFIGS.update({k: v[1] for k, v in MASK_CONFIGS.items()})


def base_options(cfg, params):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = MASK_CONFIGS[cfg][0] if cfg in MASK_CONFIGS else "eprof_v2.2"
    o.molecular_params = dict(params)
    o.plot_main = False
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


def _f(v):
    try:
        if v is None:
            return None
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def run_one(args):
    cfg, inst = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    out_f = OUT / f"cand_{cfg}_{inst['label']}.json"
    # Date-level resume (see run_baselines.py): only the nights the corpus has gained are computed.
    rec = {}
    if out_f.exists() and "--force" not in sys.argv:
        rec = json.loads(out_f.read_text())
        if not [d for d in dates_of(inst) if d not in rec]:
            return (cfg, inst["label"], len(rec),
                    sum(1 for v in rec.values() if v[0] in (1.0, 0.5)), True)
    o = base_options(cfg, CONFIGS[cfg])
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    masked = cfg in MASK_CONFIGS
    if masked and not (CLASSIF / inst["label"]).is_dir():
        return cfg, inst["label"], 0, 0, True          # not in the classified subset -> nothing to do
    for ds in dates_of(inst):
        if ds in rec:
            continue
        try:
            cls = None
            if masked:
                # BOTH days or nothing: a night covered by d alone would be masked over its morning
                # half only, and the record would look complete while the evening went unscreened.
                files = classification_files_for_night(CLASSIF, inst["label"], inst["wmo"], ds)
                if len(files) < 2:
                    continue
                cls = read_classification_curtain(files)
                if cls is None:
                    continue                            # no classification -> not a comparable night
            r = calibrate_rayleigh(ds, info, o, classification=cls)
        except Exception as exc:
            rec[ds] = [-99, None, None, None, None, f"EXC {type(exc).__name__}: {exc}"[:120]]
            continue
        if r is None:
            continue
        rec[ds] = [_f(r.flag), _f(r.lidar_constant), _f(r.uncertainty),
                   _f(r.calibration_bottom_height), _f(r.calibration_top_height),
                   str(r.message or "")[:120]]
    out_f.write_text(json.dumps(rec), encoding="utf-8")
    return cfg, inst["label"], len(rec), sum(1 for v in rec.values() if v[0] in (1.0, 0.5)), False


def main():
    only_cfg = [a for a in sys.argv[1:] if a in CONFIGS]
    pats = [a.upper() for a in sys.argv[1:] if a not in CONFIGS and not a.startswith("--")]
    cfgs = only_cfg or list(CONFIGS)
    insts = [i for i in MANIFEST
             if not pats or any(p in i["label"].upper() for p in pats)]
    jobs = [(c, i) for c in cfgs for i in insts]
    workers = int(os.environ.get("RA_WORKERS", str(max(1, (os.cpu_count() or 8) - 2))))
    print(f"candidates: {len(jobs)} (config x stream) jobs over {workers} workers; "
          f"configs {cfgs}; {len(insts)} streams", flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            j = futs[fut]
            done += 1
            try:
                cfg, label, n, n_ok, cached = fut.result()
                print(f"  [{done}/{len(jobs)}] {cfg:12s} {label:28s} {n:4d} nights, {n_ok:3d} valid"
                      f"{' (cached)' if cached else ''}", flush=True)
            except Exception as e:
                print(f"  [{done}/{len(jobs)}] FAILED {j[0]} {j[1]['label']}: {e}", flush=True)
    print("CANDIDATES_DONE", flush=True)


if __name__ == "__main__":
    main()
