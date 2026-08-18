"""
run_l1_2026_variability.py — Rayleigh-calibrate the selected CHM15k + Mini-MPL + all CL61
over the E-PROFILE **L1 2026** archive (D:/E-PROFILE_L1_2026), every available night, with
the molecular-window methods (calipso dropped). Saves per-instrument per-night JSON so the
variability metrics can be computed without re-running.

For each night the L1 file is loaded once (cloud-screened, molecular profile built, WV
correction applied on the 910 nm units), then ALL methods are evaluated on the one shared
prepared profile (`run_methods`) — the per-night cost is dominated by the single load, not by
the number of methods. Instruments run in parallel (one process each).

Instruments come from validation/scope_l1_2026.json (built by scope_l1_2026.py from the L1
`instrument_type` metadata). Output JSON format matches longrun_methods.mw_to_list so the
existing precision metrics apply unchanged.
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
from calibration.rayleigh.molecular_methods import METHODS, compute_window_grid, select_molecular_window

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "validation" / "scope_l1_2026.json"
ROOT = Path("D:/E-PROFILE_L1_2026")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/l1_2026_variability")
OUT.mkdir(parents=True, exist_ok=True)

HALF = tuple(range(250, 2000, 240))
QC_THR = 15.0                       # pipeline threshold_quality (rel_error %) -> flag -2


def base_options():
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = ROOT
    o.data_level = DataLevel.L1
    o.molecular_source = "standard"
    o.apply_wv_correction = True     # 910 nm (CL61) gets the bundled WV LUT + CAMS
    # The pipeline result r reproduces E-PROF v1.0 (sign error): the legacy main window
    # (eprof_v1.1) run through the FULL pipeline with the historical Klett sign error. The
    # prepared profile (fit_inputs_out) is built before window selection and is method- and
    # sign-independent, so the 6 live methods are evaluated on it unaffected -> one load/night
    # gives both v1.0 (from r) and the 6 method proxies (from the profile).
    o.molecular_method = "eprof_v1.1"
    o.sign_error_v10 = True
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def run_methods(signal, p_mol, rng, stack):
    """All methods on one prepared profile (optimal additionally gets the time stack)."""
    grid = compute_window_grid(signal, p_mol, rng, HALF, range_start_m=2000,
                               range_end_m=6000, increment_bins=8, signal_stack=stack)
    out = {}
    for m in METHODS:
        if m == "eprof_v2":
            out[m] = select_molecular_window("eprof_v2", signal, p_mol, rng, HALF, signal_stack=stack)
        else:
            out[m] = select_molecular_window(m, signal, p_mol, rng, HALF, grid=grid)
    return out


def calibrates(w):
    return bool(w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC_THR)


def mw_to_list(w):
    return [bool(w.ok), float(w.cl), float(w.cl_err), float(w.rel_error), float(w.r2),
            float(w.temporal_cv), float(w.scattering_ratio), float(w.start_m), float(w.end_m),
            float(w.center_m)]


def result_to_list(r):
    """E-PROF v1.0 from the full-pipeline CalibrationResult (same 10-field layout)."""
    ok = r.flag in (1, 1.0, 0.5)
    cl = float(r.lidar_constant)
    err = float(r.uncertainty)
    rel = abs(100.0 * err / cl) if cl > 0 else 999.0
    nan = float("nan")
    bot = float(getattr(r, "bottom_height", 0.0) or 0.0)
    top = float(getattr(r, "top_height", 0.0) or 0.0)
    return [bool(ok), cl, err, rel, nan, nan, nan, bot, top, (bot + top) / 2.0]


def date_strs(first, last):
    d0 = datetime.strptime(first, "%Y%m%d")
    d1 = datetime.strptime(last, "%Y%m%d")
    d = d0
    while d <= d1:
        yield d.strftime("%Y%m%d")
        d += timedelta(days=1)


def run_instrument(inst):
    """Process one instrument over its full 2026 date span. Saves results_<label>.json."""
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    o = base_options()
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=InstrumentType(inst["type"]),
                          latitude=inst["lat"], longitude=inst["lon"], altitude=inst["alt"])
    per_night = {}
    n_cal = 0
    for ds in date_strs(inst["first"], inst["last"]):
        fin = {}
        r = None
        try:
            r = calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)   # r = E-PROF v1.0 (sign error)
        except Exception:
            continue
        if not fin:
            continue
        try:
            res = run_methods(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
        except Exception:
            continue
        row = {m: mw_to_list(res[m]) for m in METHODS}
        if r is not None:
            try:
                row["eprof_v1.0"] = result_to_list(r)
            except Exception:
                pass
        per_night[ds] = row
        if any(calibrates(res[m]) for m in METHODS):
            n_cal += 1
    (OUT / f"results_{inst['label']}.json").write_text(json.dumps(per_night), encoding="utf-8")
    return inst["label"], inst["group"], len(per_night), n_cal


def main():
    insts = json.loads(MANIFEST.read_text())
    only = sys.argv[1] if len(sys.argv) > 1 else None
    if only:
        insts = [i for i in insts if only.lower() in i["label"].lower() or only == i["group"]]
    workers = int(os.environ.get("RC_WORKERS", "12"))
    print(f"Running {len(insts)} instruments over L1 2026 with {workers} workers...")
    print("  methods:", ", ".join(METHODS))
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_instrument, inst): inst for inst in insts}
        for fut in as_completed(futs):
            inst = futs[fut]
            try:
                label, group, n_fit, n_cal = fut.result()
                done += 1
                print(f"  [{done}/{len(insts)}] {label:24s} ({group:8s}) "
                      f"{n_fit:3d} fit-nights, {n_cal:3d} calibrated")
            except Exception as e:
                print(f"  FAILED {inst['label']}: {e}")
    print("L1_2026_RUN_DONE")


if __name__ == "__main__":
    main()
