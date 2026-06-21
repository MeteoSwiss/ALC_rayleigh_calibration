"""
run_cloud_sweep.py — optimize the liquid-cloud (O'Connor/Hopkin) calibration for CL31/CL51/CL61:
sweep candidate gate configurations to maximize the number of valid cloud calibrations while keeping
the night-to-night variability (sigma_SD) of the calibration coefficient low. Runs on BOTH L1 and L2.

Efficiency: the expensive per-day work (read + pre-average + water-vapour correction) is done ONCE;
each config is then evaluated by liquid_cloud_calibration_from_data on the already-WV-corrected data
(apply_wv_correction=False, averaging disabled) — only the cheap instrument/cloud/consistency gates
re-run per config.

Usage:  python run_cloud_sweep.py <phase> [level]
        phase 1 = top-10 streams per type (CL31/CL51/CL61); phase 2 = full network.
        level omitted = both L1 and L2.

Output: figs_paper_validation/cloud_sweep/cloud_<phase>_<level>_<label>.json
        = {ds: {config: [cal_median, n_profiles]}}.
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    CloudCalConfig, set_defaults, read_ceilometer_data, average_ceilo_data,
    compute_wv_transmission, liquid_cloud_calibration_from_data)

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_sweep")
OUT.mkdir(parents=True, exist_ok=True)
CAMS = "D:/CAMS"
WV_LUT = str(REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc")
ROOTS = {"L1": Path("D:/E-PROFILE_L1_2026"), "L2": Path("D:/E-PROFILE_L2_2026")}

# Candidate cloud-calibration configs (gate overrides; defaults: n_consec=5, consistency_range=10,
# ratio_filter=0.05, cbh_max=2400, cal_max=2400, temp=-20, attenuation_factor=20).
CONFIGS = {
    "K0_baseline":    dict(),
    "K1_consec3":     dict(n_consecutive=3),
    "K2_consist20":   dict(consistency_range=20.0),
    "K3_ratio0.1":    dict(ratio_filter=0.10),
    "K4_cbh3000":     dict(cbh_maxheight=3000.0, cal_maxheight=3000.0),
    "K5_tempcold":    dict(temp_threshold=-25.0),
    "K6_balanced":    dict(n_consecutive=3, consistency_range=15.0, ratio_filter=0.08,
                           cbh_maxheight=3000.0, cal_maxheight=3000.0),
    "K7_aggressive":  dict(n_consecutive=3, consistency_range=25.0, ratio_filter=0.15,
                           cbh_maxheight=3500.0, cal_maxheight=3500.0, attenuation_factor=10.0),
}


def base_config(nc_file, inst, lat, lon):
    return CloudCalConfig(
        nc_file=str(nc_file), instrument=inst, apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=CAMS, abs_cs_lookup_table=WV_LUT,
        station_latitude=lat, station_longitude=lon,
        average_time_s=300.0, average_range_m=10.0)


def eval_day(nc_file, inst, lat, lon):
    """Read + average + WV ONCE, then evaluate every config on the WV-corrected data."""
    base = set_defaults(base_config(nc_file, inst, lat, lon))
    data, status = read_ceilometer_data(base.nc_file, base)
    if status != 0:
        return None
    data = average_ceilo_data(data, base)
    in_band = inst.upper() in ("CL31", "CL51", "CL61")
    if in_band:
        trans2 = compute_wv_transmission(data, base)          # raises if CAMS missing
        trans2 = np.asarray(trans2, float)
        if trans2.size == 0 or not np.any(np.isfinite(trans2)) or np.all(trans2 == 1):
            return None
        data.beta = data.beta / trans2
        data.trans2_wv = trans2
    out = {}
    for name, gates in CONFIGS.items():
        cfg = replace(base, apply_wv_correction=False, average_time_s=None, average_range_m=None, **gates)
        try:
            res = liquid_cloud_calibration_from_data(data, cfg)
            out[name] = [float(res.cal_median), int(res.n_profiles)]
        except Exception:
            out[name] = [float("nan"), 0]
    return out


def run_one(args):
    level, inst = args
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    root = ROOTS[level]
    prefix = "L1" if level == "L1" else "L2"
    d0 = datetime.strptime(inst["first"], "%Y%m%d")
    d1 = datetime.strptime(inst["last"], "%Y%m%d")
    per_day = {}
    d = d0
    while d <= d1:
        ds = d.strftime("%Y%m%d"); d += timedelta(days=1)
        fp = root / inst["wmo"] / "2026" / ds[4:6] / f"{prefix}_{inst['wmo']}_{inst['ident']}{ds}.nc"
        if not fp.exists():
            continue
        try:
            r = eval_day(fp, inst["type"], inst["lat"], inst["lon"])
        except Exception:
            continue
        if r is not None:
            per_day[ds] = r
    tag = f"{sys.argv[1]}"
    (OUT / f"cloud_{tag}_{level}_{inst['label']}.json").write_text(json.dumps(per_day), encoding="utf-8")
    nval = sum(1 for v in per_day.values()
               if np.isfinite(v["K0_baseline"][0]) and v["K0_baseline"][0] > 0 and v["K0_baseline"][1] >= 1)
    return level, inst["label"], inst["group"], len(per_day), nval


def select(phase):
    if phase == "2":
        return list(MANIFEST)
    sub = []
    for t in ("CL31", "CL51", "CL61"):
        sub += sorted([m for m in MANIFEST if m["group"] == t], key=lambda m: -m["n_days"])[:10]
    return sub


def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "1"
    only_level = sys.argv[2] if len(sys.argv) > 2 else None
    insts = select(phase)
    jobs = [(lvl, m) for lvl in ROOTS for m in insts if only_level is None or lvl == only_level]
    workers = int(os.environ.get("RC_WORKERS", "12"))
    print(f"cloud sweep phase {phase}: {len(jobs)} (level x stream) jobs, {len(CONFIGS)} configs each, {workers} workers")
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            try:
                level, label, group, ndays, nval = fut.result()
                done += 1
                print(f"  [{done}/{len(jobs)}] {level} {label:22s} ({group}) days={ndays:3d} baseline-valid={nval:3d}")
            except Exception as e:
                print(f"  FAILED: {e}")
    print("CLOUD_SWEEP_DONE")


if __name__ == "__main__":
    main()
