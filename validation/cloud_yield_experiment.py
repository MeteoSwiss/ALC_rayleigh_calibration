"""cloud_yield_experiment.py -- test 5 propositions to raise the number of DASHBOARD-VALID liquid-cloud
calibrations from L1 data (valid = a positive lidar constant C_L, the dashboard's success criterion,
NOT just a finite O'Connor coefficient).

Why the distinction matters: CL61 L1 carries no applied calibration constant, so C_L = const / coef
is NaN even when a good coefficient is found -> 0 valid CL61 on the dashboard. CL31/CL51 default to
const=1e8 and are fine. So yield levers help CL31/CL51 directly, but CL61 ALSO needs a C_L fallback.

Per station-day we read L1 once per AVERAGING variant, water-vapour-correct once, then run the cheap
gate configs. Stored: {ds: {variant: {config: [coef, n_profiles]}}}. The analysis step turns these
into per-(type, proposition) valid% and sigma_SD using the type C_L constant (with/without the CL61
fallback). Slice-parallel:  python cloud_yield_experiment.py <slice_i> <n_slices>
"""
from __future__ import annotations
import json
import os
import sys
import warnings
from dataclasses import replace
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root for calibration/validation
from calibration.cloud.calibration import (
    set_defaults, read_ceilometer_data, average_ceilo_data,
    compute_wv_transmission, liquid_cloud_calibration_from_data)
from validation.run_cloud_sweep import base_config, CONFIGS, MANIFEST, ROOTS

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_yield")
OUT.mkdir(parents=True, exist_ok=True)

# Averaging variants (time_s, range_m). The chosen production change is 30 s / 10 m vs the old 300 s.
VARIANTS = {"base_300s": (300.0, 10.0), "fine_30s": (30.0, 10.0)}
# Keep the LITERATURE-default gates (K0); K1/K7 retained only for context.
CFG = {"K0": CONFIGS["K0_baseline"], "K1": CONFIGS["K1_consec3"], "K7": CONFIGS["K7_aggressive"]}

START, END = datetime(2026, 1, 1), datetime(2026, 5, 31)
DAYSTEP = int(os.environ.get("CY_DAYSTEP", "3"))     # sample every Nth day
MAXDAYS = int(os.environ.get("CY_MAXDAYS", "32"))    # cap days-with-data per station


def select():
    """All CL61 + top-10 CL31 + top-10 CL51 by coverage."""
    sub = []
    for t, n in (("CL61", 999), ("CL31", 10), ("CL51", 10)):
        sub += sorted([m for m in MANIFEST if m["group"] == t], key=lambda m: -m["n_days"])[:n]
    return sub


def l1_path(inst, ds):
    return ROOTS["L1"] / inst["wmo"] / "2026" / ds[4:6] / f"L1_{inst['wmo']}_{inst['ident']}{ds}.nc"


def eval_day(fp, inst):
    """Per averaging variant: read + average + WV once, then run each config. {variant:{config:[coef,n]}}."""
    out = {}
    for vname, (avg_s, avg_r) in VARIANTS.items():
        base = set_defaults(base_config(fp, inst["type"], inst["lat"], inst["lon"]))
        data, status = read_ceilometer_data(base.nc_file, base)
        if status != 0:
            return None
        cfg = replace(base, average_time_s=avg_s, average_range_m=avg_r)
        if avg_s or avg_r:
            data = average_ceilo_data(data, cfg)
        try:
            trans2 = np.asarray(compute_wv_transmission(data, cfg), float)
        except Exception:
            return None
        if trans2.size == 0 or not np.any(np.isfinite(trans2)) or np.all(trans2 == 1):
            out[vname] = {c: [float("nan"), 0] for c in CFG}
            continue
        data.beta = data.beta / trans2
        data.trans2_wv = trans2
        res_v = {}
        for cname, gates in CFG.items():
            c2 = replace(cfg, apply_wv_correction=False, average_time_s=None, average_range_m=None, **gates)
            try:
                r = liquid_cloud_calibration_from_data(data, c2)
                res_v[cname] = [float(r.cal_median), int(r.n_profiles)]
            except Exception:
                res_v[cname] = [float("nan"), 0]
        out[vname] = res_v
    return out


def main():
    si, ns = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
    streams = select()[si::ns]
    for inst in streams:
        outfp = OUT / f"cy_{inst['label']}.json"
        if outfp.exists():
            continue
        days = {}
        d, ndata = START, 0
        while d <= END and ndata < MAXDAYS:
            ds = d.strftime("%Y%m%d")
            d += timedelta(days=DAYSTEP)
            fp = l1_path(inst, ds)
            if not fp.exists():
                continue
            try:
                r = eval_day(fp, inst)
            except Exception:
                r = None
            if r is not None:
                days[ds] = r
                ndata += 1
        outfp.write_text(json.dumps(days))
        print(f"{inst['label']} ({inst['group']}): {len(days)} days", flush=True)
    print("SLICE DONE", flush=True)


if __name__ == "__main__":
    main()
