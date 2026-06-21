"""
run_cloud_avg_sweep.py — find the optimal pre-averaging for the liquid-cloud calibration. On NATIVE L1
data (which can be averaged to any resolution), sweep average_time_s and measure the number of valid
cloud calibrations and the coefficient's short-term variability (sigma_SD), with the recommended K6
gates fixed. Range averaging is fixed at 30 m (matching L2). The expensive CAMS read is cached, but the
WV transmission is re-interpolated per averaging level (the grid changes), so each level pays a WV cost.

Usage:  python run_cloud_avg_sweep.py <slice_i> <n_slices>   (single-process slice; shell-parallel)
Output: figs_paper_validation/cloud_avg/avg_<label>.json = {ds: {avg_s: [cal_median, n_profiles]}}
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
    set_defaults, read_ceilometer_data, average_ceilo_data, compute_wv_transmission,
    liquid_cloud_calibration_from_data)
from validation.run_cloud_sweep import base_config, CONFIGS

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_avg")
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("D:/E-PROFILE_L1_2026")
AVG_LEVELS = [30, 60, 120, 300, 600, 900, 1200]     # average_time_s to test
K6 = CONFIGS["K6_balanced"]                          # recommended gates, fixed
NTYPE = int(os.environ.get("RC_NTYPE", "3"))
DAYSTEP = int(os.environ.get("RC_DAYSTEP", "5"))


def select():
    sub = []
    for t in ("CL31", "CL51", "CL61"):
        sub += sorted([m for m in MANIFEST if m["group"] == t], key=lambda m: -m["n_days"])[:NTYPE]
    return sub


def eval_day_avg(nc_file, inst):
    """Read native ONCE, then for each averaging level: average + WV + K6 gates."""
    base = set_defaults(base_config(nc_file, inst["type"], inst["lat"], inst["lon"]))
    data0, status = read_ceilometer_data(base.nc_file, base)
    if status != 0:
        return None
    out = {}
    for avg_s in AVG_LEVELS:
        cfg = replace(base, average_time_s=float(avg_s), average_range_m=30.0)
        try:
            data = average_ceilo_data(data0, cfg)
            trans2 = compute_wv_transmission(data, cfg)
            trans2 = np.asarray(trans2, float)
            if trans2.size == 0 or not np.any(np.isfinite(trans2)) or np.all(trans2 == 1):
                out[avg_s] = [float("nan"), 0]; continue
            data.beta = data.beta / trans2; data.trans2_wv = trans2
            cfg2 = replace(cfg, apply_wv_correction=False, average_time_s=None, average_range_m=None, **K6)
            res = liquid_cloud_calibration_from_data(data, cfg2)
            out[avg_s] = [float(res.cal_median), int(res.n_profiles)]
        except Exception:
            out[avg_s] = [float("nan"), 0]
    return out


def run_slice(slice_i, n_slices):
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.CRITICAL)
    insts = select()[slice_i::n_slices]
    print(f"avg-sweep slice {slice_i}/{n_slices}: {len(insts)} streams, levels={AVG_LEVELS}", flush=True)
    for k, inst in enumerate(insts):
        d0 = datetime.strptime(inst["first"], "%Y%m%d"); d1 = datetime.strptime(inst["last"], "%Y%m%d")
        per_day = {}; d = d0
        while d <= d1:
            ds = d.strftime("%Y%m%d"); d += timedelta(days=DAYSTEP)
            fp = ROOT / inst["wmo"] / "2026" / ds[4:6] / f"L1_{inst['wmo']}_{inst['ident']}{ds}.nc"
            if not fp.exists():
                continue
            try:
                r = eval_day_avg(fp, inst)
            except Exception:
                continue
            if r is not None:
                per_day[ds] = r
        (OUT / f"avg_{inst['label']}.json").write_text(json.dumps(per_day), encoding="utf-8")
        print(f"  [{slice_i}:{k+1}/{len(insts)}] {inst['label']} ({inst['group']}) days={len(per_day)}", flush=True)
    print(f"AVGSLICE_{slice_i}_DONE", flush=True)


if __name__ == "__main__":
    run_slice(int(sys.argv[1]), int(sys.argv[2]))
