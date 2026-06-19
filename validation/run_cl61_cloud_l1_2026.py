"""
run_cl61_cloud_l1_2026.py — liquid-water-cloud calibration (O'Connor/Hopkin) of the CL61 over
the L1 2026 archive, as an INDEPENDENT cross-check of the Rayleigh calibration.

For every CL61 (from validation/scope_l1_2026.json, group==CL61) and every day in its span we
run the cloud calibration on the daily file (WV correction mandatory at 910 nm) and store the
per-day calibration coefficient. The cloud method needs liquid clouds in view, not a clean
molecular night, so it samples a different population of days than the Rayleigh fit — which is
exactly what makes it a useful cross-check.

Data level: the cloud reader is bit-for-bit validated on E-PROFILE **monthly L2** (the O'Connor
method needs a full month of profiles to accumulate enough liquid-cloud returns — daily files are
too sparse and yield no calibration). We therefore run the cross-check on the **monthly L2 archive**
(A:/E-PROFILE_L2_monthly), one calibration coefficient per month per CL61, using the same full
spectral WV absorption LUT the module was validated against (abs_cross_647_full_levels_1000.nc).

Output: figs_paper_validation/l1_2026_variability/cloud_<label>.json  = {YYYYMM: [ok, C, std, n_prof]}.
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

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
OPTS = json.loads((REPO / "options.json").read_text())
ROOT = Path("A:/E-PROFILE_L2_monthly")    # monthly L2: enough cloud profiles for the O'Connor method
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/l1_2026_variability")
OUT.mkdir(parents=True, exist_ok=True)
LR_AEROSOL = float(OPTS.get("LRaer", 52))
# The cloud WV correction needs the full spectral absorption LUT the module was validated against
# (the bundled 910 nm LUT is NOT equivalent for the cloud path). options.json may override.
_VALIDATED_LUT = "C:/Users/hervo/OneDrive/Documents/MATLAB/MDA/monitoring_alc_monthly/abs_cross_647_full_levels_1000.nc"
WV_LUT = OPTS.get("abs_cs_lookup_table", "") or _VALIDATED_LUT


def file_path(wmo, ident, ym):
    return ROOT / wmo / ym[:4] / f"L2_{wmo}_{ident}{ym}.nc"


def months_in(first, last):
    """YYYYMM strings spanning the instrument's first..last daily dates."""
    y0, m0 = int(first[:4]), int(first[4:6])
    y1, m1 = int(last[:4]), int(last[4:6])
    out, y, m = [], y0, m0
    while (y, m) <= (y1, m1):
        out.append(f"{y:04d}{m:02d}")
        m += 1
        if m > 12:
            y, m = y + 1, 1
    return out


def run_instrument(inst):
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    from calibration.cloud import liquid_cloud_calibration, CloudCalConfig
    out = {}
    n_ok = 0
    for ym in months_in(inst["first"], inst["last"]):
        f = file_path(inst["wmo"], inst["ident"], ym)
        if not f.is_file():
            continue
        cfg = CloudCalConfig(
            nc_file=str(f),
            instrument="CL61",
            apply_wv_correction=True,
            # The optional above-cloud aerosol-transmission refinement currently collapses the
            # median coefficient to 0 on the 2026 monthly L2 (a regression in that step); the base
            # O'Connor coefficient (without it) is valid and is what we cross-check against.
            apply_transmission_correction=False,
            cams_folder=OPTS.get("cams_folder", "D:/CAMS/"),
            abs_cs_lookup_table=WV_LUT,
            station_latitude=inst["lat"],
            station_longitude=inst["lon"],
            aerosol_lidar_ratio=LR_AEROSOL,
        )
        try:
            res = liquid_cloud_calibration(cfg)
            n_prof = int(getattr(res, "n_profiles", 0))
            coef = float(res.cal_mean) if n_prof > 0 else float("nan")
            std = float(res.cal_std) if n_prof > 0 else float("nan")
            ok = bool(n_prof > 0 and np.isfinite(coef) and coef > 0)
        except Exception:
            ok, coef, std, n_prof = False, float("nan"), float("nan"), 0
        out[ym] = [ok, coef, std, n_prof]
        n_ok += int(ok)
    (OUT / f"cloud_{inst['label']}.json").write_text(json.dumps(out), encoding="utf-8")
    return inst["label"], len(out), n_ok


def main():
    cl61 = [m for m in MANIFEST if m["group"] == "CL61"]
    only = sys.argv[1] if len(sys.argv) > 1 else None
    if only:
        cl61 = [m for m in cl61 if only.lower() in m["label"].lower()]
    workers = int(os.environ.get("RC_WORKERS", "10"))
    print(f"CL61 cloud calibration cross-check: {len(cl61)} instruments, {workers} workers")
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_instrument, m): m for m in cl61}
        for fut in as_completed(futs):
            m = futs[fut]
            try:
                label, n_days, n_ok = fut.result()
                done += 1
                print(f"  [{done}/{len(cl61)}] {label:28s} {n_days:3d} days, {n_ok:3d} cloud-calibrated")
            except Exception as e:
                print(f"  FAILED {m['label']}: {e}")
    print("CL61_CLOUD_DONE")


if __name__ == "__main__":
    main()
