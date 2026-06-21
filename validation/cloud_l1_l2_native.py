"""
cloud_l1_l2_native.py — definitive L1-vs-L2 / averaging / config comparison for the liquid-cloud
calibration. For streams that have BOTH L1 and L2, compute valid% and sigma_SD over the cross product:

  variants : L1_native (L1, native time, 30 m range)   <- finest
             L1_300s   (L1 averaged to 300 s x 30 m)    <- L2 grid, from L1
             L2        (L2 as delivered, 300 s x 30 m)   <- cleaned attenuated backscatter

  configs  : K0 (MATLAB baseline)   K6 (balanced)   K7 (aggressive)

This answers (a) does changing the gate config K0 -> K6/K7 add value, and (b) is native-resolution
L1 actually better than L2 (the averaging finding implies it should be). The expensive WV correction
is recomputed per variant (the grid changes). Range fixed at 30 m so only the time axis / data level
varies. sigma_SD is convention-independent, reported as % of median C_L.

Usage:  python cloud_l1_l2_native.py <slice_i> <n_slices>
Output: figs_paper_validation/cloud_l1l2/ll_<label>.json = {ds: {variant: {config: [cal_median, n]}}}
"""
from __future__ import annotations
import json, logging, os, sys, warnings
from dataclasses import replace
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

from calibration.cloud.calibration import (
    set_defaults, read_ceilometer_data, average_ceilo_data,
    compute_wv_transmission, liquid_cloud_calibration_from_data)
from validation.run_cloud_sweep import base_config, CONFIGS

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_l1l2")
OUT.mkdir(parents=True, exist_ok=True)
ROOTS = {"L1": Path("D:/E-PROFILE_L1_2026"), "L2": Path("D:/E-PROFILE_L2_2026")}
CFG = {"K0": CONFIGS["K0_baseline"], "K6": CONFIGS["K6_balanced"], "K7": CONFIGS["K7_aggressive"]}
# variant -> (level, average_time_s, average_range_m)
VARIANTS = {
    "L1_native": ("L1", None, 30.0),
    "L1_300s":   ("L1", 300.0, 30.0),
    "L2":        ("L2", None, None),
}
NTYPE = int(os.environ.get("RC_NTYPE", "4"))
DAYSTEP = int(os.environ.get("RC_DAYSTEP", "10"))


def select():
    sub = []
    for t in ("CL31", "CL51", "CL61"):
        sub += sorted([m for m in MANIFEST if m["group"] == t], key=lambda m: -m["n_days"])[:NTYPE]
    return sub


def path_for(level, inst, ds):
    prefix = level
    return ROOTS[level] / inst["wmo"] / "2026" / ds[4:6] / f"{prefix}_{inst['wmo']}_{inst['ident']}{ds}.nc"


def eval_variant(fp, inst, avg_s, avg_r):
    """Read fp, (optionally) average, WV-correct once, then run K0/K6/K7. Returns {config:[cal_median,n]}."""
    base = set_defaults(base_config(fp, inst["type"], inst["lat"], inst["lon"]))
    data, status = read_ceilometer_data(base.nc_file, base)
    if status != 0:
        return None
    cfg = replace(base, average_time_s=avg_s, average_range_m=avg_r)
    if avg_s or avg_r:
        data = average_ceilo_data(data, cfg)
    trans2 = np.asarray(compute_wv_transmission(data, cfg), float)
    if trans2.size == 0 or not np.any(np.isfinite(trans2)) or np.all(trans2 == 1):
        return {name: [float("nan"), 0] for name in CFG}
    data.beta = data.beta / trans2
    data.trans2_wv = trans2
    out = {}
    for name, gates in CFG.items():
        cfg2 = replace(cfg, apply_wv_correction=False, average_time_s=None, average_range_m=None, **gates)
        try:
            res = liquid_cloud_calibration_from_data(data, cfg2)
            out[name] = [float(res.cal_median), int(res.n_profiles)]
        except Exception:
            out[name] = [float("nan"), 0]
    return out


def eval_day(inst, ds):
    day = {}
    for variant, (level, avg_s, avg_r) in VARIANTS.items():
        fp = path_for(level, inst, ds)
        if not fp.exists():
            continue
        try:
            r = eval_variant(fp, inst, avg_s, avg_r)
        except Exception:
            r = None
        if r is not None:
            day[variant] = r
    return day


def run_slice(slice_i, n_slices):
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.CRITICAL)
    insts = select()[slice_i::n_slices]
    print(f"l1l2 slice {slice_i}/{n_slices}: {len(insts)} streams, variants={list(VARIANTS)}", flush=True)
    for k, inst in enumerate(insts):
        d0 = datetime.strptime(inst["first"], "%Y%m%d"); d1 = datetime.strptime(inst["last"], "%Y%m%d")
        per_day = {}; d = d0
        while d <= d1:
            ds = d.strftime("%Y%m%d"); d += timedelta(days=DAYSTEP)
            day = eval_day(inst, ds)
            if day:
                per_day[ds] = day
        (OUT / f"ll_{inst['label']}.json").write_text(json.dumps(per_day), encoding="utf-8")
        print(f"  [{slice_i}:{k+1}/{len(insts)}] {inst['label']} ({inst['group']}) days={len(per_day)}", flush=True)
    print(f"LL_SLICE_{slice_i}_DONE", flush=True)


if __name__ == "__main__":
    run_slice(int(sys.argv[1]), int(sys.argv[2]))
