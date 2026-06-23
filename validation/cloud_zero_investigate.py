"""Why do some streams return 0 valid cloud calibrations under the new 30 s / K0 / fallback version?
Identify the zero-yield streams from the experiment output, then re-run them at 30 s + K0 and tabulate
WHY each day failed: genuinely no liquid cloud (-1) vs the dominant filter rejection (window / energy /
peak-above / peak-below / aerosol-ratio / cbh / consistency) vs no usable data.
Run: python validation/cloud_zero_investigate.py
"""
import glob
import json
import sys
from collections import Counter
from datetime import datetime, timedelta
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from validation.run_cloud_sweep import base_config, MANIFEST, ROOTS, CONFIGS
from calibration.cloud.calibration import (
    set_defaults, read_ceilometer_data, average_ceilo_data, compute_wv_transmission,
    liquid_cloud_calibration_from_data)
from calibration.flags import dominant_cloud_reject_flag

CYD = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_yield")
META = {m["label"]: m for m in MANIFEST}
START = datetime(2026, 1, 1)


def zero_streams():
    """Streams whose fine_30s/K0 yielded 0 valid days (n>=1, coef>0)."""
    zero, good = [], []
    for fp in glob.glob(str(CYD / "cy_*.json")):
        lab = Path(fp).stem[3:]
        if lab not in META:
            continue
        days = json.loads(Path(fp).read_text())
        nv = sum(1 for ds in days
                 if (lambda c: c and c[1] >= 1 and np.isfinite(c[0]) and c[0] > 0)(days[ds].get("fine_30s", {}).get("K0")))
        (zero if nv == 0 else good).append((lab, META[lab]["group"], len(days), nv))
    return zero, good


def l1(inst, ds):
    return ROOTS["L1"] / inst["wmo"] / "2026" / ds[4:6] / f"L1_{inst['wmo']}_{inst['ident']}{ds}.nc"


def investigate(lab, ndays=24):
    m = META[lab]
    reasons = Counter()
    d, seen = START, 0
    while seen < ndays and d <= datetime(2026, 5, 31):
        ds = d.strftime("%Y%m%d"); d += timedelta(days=3)
        fp = l1(m, ds)
        if not fp.exists():
            continue
        base = set_defaults(replace(base_config(fp, m["type"], m["lat"], m["lon"]),
                                    average_time_s=30.0, average_range_m=10.0))
        try:
            data, st = read_ceilometer_data(base.nc_file, base)
            if st != 0 or data is None or not np.any(np.isfinite(np.asarray(getattr(data, "beta", np.nan), float))):
                reasons["no_data"] += 1; seen += 1; continue
            data = average_ceilo_data(data, base)
            tr = np.asarray(compute_wv_transmission(data, base), float)
            if tr.size == 0 or np.all(tr == 1):
                reasons["no_wv/cams"] += 1; seen += 1; continue
            data.beta = data.beta / tr
            res = liquid_cloud_calibration_from_data(data, base)
        except Exception as e:
            reasons[f"exc:{type(e).__name__}"] += 1; seen += 1; continue
        seen += 1
        if int(getattr(res, "n_profiles", 0)) >= 1:
            reasons["VALID"] += 1
        else:
            _, reason, counts = dominant_cloud_reject_flag(
                getattr(res, "filter_stats", None), getattr(res, "cloud_stats", None),
                getattr(res, "consistency_stats", None))
            reasons[reason] += 1
    print(f"  {lab:24s} {m['group']:5s} days={seen:2d}  {dict(reasons)}")


def main():
    zero, good = zero_streams()
    print(f"ZERO-yield streams (fine_30s/K0): {len(zero)}  |  good: {len(good)}")
    print("=== WHY (zero-yield streams) ===")
    for lab, t, nd, nv in sorted(zero, key=lambda x: x[1]):
        investigate(lab)
    print("=== contrast: 2 good streams per type ===")
    bytype = {}
    for lab, t, nd, nv in good:
        bytype.setdefault(t, []).append(lab)
    for t, labs in bytype.items():
        for lab in labs[:2]:
            investigate(lab)
    print("DONE")


if __name__ == "__main__":
    main()
