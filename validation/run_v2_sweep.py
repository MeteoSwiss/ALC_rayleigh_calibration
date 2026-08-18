"""
run_v2_sweep.py — optimize the E-PROF v2 (optimal) molecular-window method: sweep candidate gate
configurations over clear-sky nights for CHM15k + Mini-MPL + CL61, on BOTH L1 and L2, and measure
the valid-calibration fraction and the short-term variability (sigma_SD) of each.

Also runs a leave-one-gate-out diagnostic on the baseline (relax each gate to infinity in turn) to
reveal WHY clear nights fail per instrument type — i.e. which gate is the binding constraint.

The per-night cost is dominated by the single load+prepare (calibrate_rayleigh -> fit_inputs); all
configs are then evaluated on that one prepared profile (cheap), exactly like the method comparison.
A "clear-sky night" = a fit-night (the profile reached the molecular fit; cloudy/fog nights return
before it). "valid" = the v2 window passes the pipeline QC (rel_error <= 15 %).

Output: figs_paper_validation/v2_sweep/sweep_<level>_<label>.json
        = {ds: {config: [ok, cl, rel_error, temporal_cv, start_m]}}.
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
from calibration.rayleigh.molecular_methods import (
    select_molecular_window, compute_window_grid, flag_contaminated_cells)

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/v2_sweep")
OUT.mkdir(parents=True, exist_ok=True)
HALF = tuple(range(250, 2000, 240))
QC = 15.0

ROOTS = {"L1": (Path("D:/E-PROFILE_L1_2026"), DataLevel.L1),
         "L2": (Path("D:/E-PROFILE_L2_2026"), DataLevel.L2_DAILY)}

INF = 1e9
# Candidate v2 configurations (gate overrides; max_rel_error stays 15 = the pipeline QC, since a
# window above it can never be a valid calibration). C0 is the current production v2.
CONFIGS = {
    "C0_baseline":    dict(),
    "C1_tcv0.8":      dict(max_temporal_cv=0.8),
    "C2_tcv1.2":      dict(max_temporal_cv=1.2),
    "C3_shape":       dict(max_residual_pct=20.0, max_scattering_ratio=1.15, max_ratio_std=0.40),
    "C4_r2_0.35":     dict(min_r2=0.35),
    "C5_start1200":   dict(min_window_start_m=1200.0),
    "C6_balanced":    dict(max_temporal_cv=0.8, max_residual_pct=16.0, min_r2=0.40,
                           min_window_start_m=1500.0, max_ratio_std=0.40, max_scattering_ratio=1.12),
    "C7_aggressive":  dict(max_temporal_cv=1.5, max_residual_pct=25.0, min_r2=0.30,
                           min_window_start_m=1000.0, max_scattering_ratio=1.25, max_ratio_std=0.50),
    # C8 = investigation-informed recommendation: relax the binding scattering-ratio gate to 1.15
    # (the dominant cause of clear-night failures, esp. CL61) + the moderate C6 relaxations.
    "C8_recommended": dict(max_temporal_cv=0.8, max_residual_pct=16.0, min_r2=0.40,
                           min_window_start_m=1500.0, max_ratio_std=0.40, max_scattering_ratio=1.15),
    # leave-one-gate-out from the baseline (relax ONE gate) -> failure diagnostic
    "LOO_temporal":   dict(max_temporal_cv=INF),
    "LOO_r2":         dict(min_r2=0.0),
    "LOO_residual":   dict(max_residual_pct=INF),
    "LOO_scattering": dict(max_scattering_ratio=INF),
    "LOO_ratiostd":   dict(max_ratio_std=INF),
    "LOO_start":      dict(min_window_start_m=0.0),
}


def base_options(level):
    root, dl = ROOTS[level]
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = root
    o.data_level = dl
    o.molecular_source = "standard"
    o.apply_wv_correction = True
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def date_strs(first, last):
    d0, d1 = datetime.strptime(first, "%Y%m%d"), datetime.strptime(last, "%Y%m%d")
    d = d0
    while d <= d1:
        yield d.strftime("%Y%m%d")
        d += timedelta(days=1)


def eval_night(signal, p_mol, rng, stack):
    """Run every config's v2 selection on one prepared profile -> {config: [ok,cl,rel,tcv,start]}.

    The gate sweep does not touch v2's time-resolved flagging (flag_nmad/min_excess), so the
    cleaned-mean profile and its window grid are identical across configs: build them ONCE and
    re-run only the cheap eligibility+selection per config (grid passed -> flagging skipped).
    """
    stack = np.asarray(stack, float)
    if stack.ndim == 2 and stack.shape[0] >= 5:
        flag = flag_contaminated_cells(stack, p_mol, rng, nmad=4.0, min_excess=0.25)
        masked = np.where(flag, np.nan, stack)
        with np.errstate(invalid="ignore"):
            sig = np.nanmean(masked, axis=0)
        grid = compute_window_grid(sig, p_mol, rng, HALF, range_start_m=2000, range_end_m=6000,
                                   increment_bins=8, signal_stack=masked)
    else:
        sig = signal
        grid = compute_window_grid(sig, p_mol, rng, HALF, range_start_m=2000, range_end_m=6000,
                                   increment_bins=8, signal_stack=stack)
    out = {}
    for name, params in CONFIGS.items():
        try:
            w = select_molecular_window("eprof_v2", sig, p_mol, rng, HALF, grid=grid, **params)
            ok = bool(w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC and np.isfinite(w.cl) and w.cl > 0)
            out[name] = [ok, float(w.cl), float(w.rel_error), float(w.temporal_cv), float(w.start_m)]
        except Exception:
            out[name] = [False, float("nan"), float("nan"), float("nan"), float("nan")]
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
    (OUT / f"sweep_{level}_{inst['label']}.json").write_text(json.dumps(per_night), encoding="utf-8")
    n_base = sum(v["C0_baseline"][0] for v in per_night.values())
    return level, inst["label"], inst["group"], len(per_night), n_base


def main():
    only_level = sys.argv[1] if len(sys.argv) > 1 else None
    jobs = [(lvl, inst) for lvl in ROOTS for inst in MANIFEST
            if only_level is None or lvl == only_level]
    workers = int(os.environ.get("RC_WORKERS", "14"))
    print(f"v2 sweep: {len(jobs)} (level x instrument) jobs, {len(CONFIGS)} configs each, {workers} workers")
    print("  configs:", ", ".join(CONFIGS))
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            j = futs[fut]
            try:
                level, label, group, n_fit, n_base = fut.result()
                done += 1
                print(f"  [{done}/{len(jobs)}] {level} {label:26s} ({group:8s}) "
                      f"{n_fit:3d} fit-nights, baseline {n_base:3d} valid")
            except Exception as e:
                print(f"  FAILED {j[0]} {j[1]['label']}: {e}")
    print("V2_SWEEP_DONE")


if __name__ == "__main__":
    main()
