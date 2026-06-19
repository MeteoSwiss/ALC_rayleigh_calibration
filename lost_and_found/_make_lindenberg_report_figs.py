"""
_make_lindenberg_report_figs.py — regenerate the Lindenberg CL61 diagnostic figures and
Kalman CSVs from the CURRENT calibration CSVs, without recomputing any calibration.

Reuses every plotting / Kalman helper from run_lindenberg_cl61_cal.py.  Snapshots the CSVs
first so a still-running calibration process cannot corrupt the read.  Safe to run while the
main script is still filling the cloud CSV: the main run regenerates these same outputs (with
complete data) when it finishes.
"""
from __future__ import annotations
import os, sys, shutil, warnings, logging
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import types
from pathlib import Path
import numpy as np

# The heavy rayleigh_calibration package is only needed by run_lindenberg_cl61_cal for
# calibrate_rayleigh (which we never call here — we only reuse its plotting/Kalman helpers).
# OneDrive is mid-sync on that package, so stub it out to keep the import clean.
for _name, _attrs in {
    "rayleigh_calibration": ["calibrate_rayleigh", "CalibrationOptions",
                             "InstrumentInfo", "DataLevel"],
    "rayleigh_calibration.config": ["InstrumentType"],
}.items():
    _m = types.ModuleType(_name)
    for _a in _attrs:
        setattr(_m, _a, object)
    sys.modules[_name] = _m

import run_lindenberg_cl61_cal as L


def _snapshot(src: Path) -> Path:
    dst = src.with_name(src.stem + "_snapshot.csv")
    shutil.copy2(src, dst)
    return dst


def main():
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    OUT = L.OUT

    ray_csv = _snapshot(OUT / "rayleigh_lindenberg_cl61.csv")
    cloud_csv = _snapshot(OUT / "cloud_lindenberg_cl61.csv")

    ray_rows = L.load_csv(ray_csv)
    cloud_rows = L._load_cloud_csv_migrated(cloud_csv)

    ray_series = L.extract_rayleigh(ray_rows)
    cloud_series = L.extract_cloud(cloud_rows)

    for tag in ("wv", "nowv"):
        t, v = ray_series.get(tag, (np.array([]), np.array([])))
        tc, vc = cloud_series.get(tag, (np.array([]), np.array([])))
        if v.size:
            print(f"  Rayleigh {tag}: {v.size} nights, median CL = {np.median(v):.3f}")
        if vc.size:
            print(f"  Cloud    {tag}: {vc.size} days,   median CL = {np.median(vc):.3f}")

    print("Building diagnostic plots...")
    ray_k = L.plot_series(
        ray_series, "rayleigh",
        ylabel=r"Rayleigh lidar constant  $C_L$  (Wiegner 2014)  [—]",
        title="Lindenberg CL61 — Rayleigh calibration (improved method) + E-PROFILE Kalman",
        outfile=OUT / "lindenberg_cl61_rayleigh_diag.png",
    )
    cloud_k = L.plot_series(
        cloud_series, "cloud",
        ylabel=r"Cloud lidar constant  $C_L = 1/C$  (Wiegner 2014)  [—]",
        title="Lindenberg CL61 — Liquid-cloud calibration ($C_L = 1/C$) + E-PROFILE Kalman",
        outfile=OUT / "lindenberg_cl61_cloud_diag.png",
    )
    L.plot_wv_impact(ray_k, cloud_k, OUT / "lindenberg_cl61_wv_impact.png")
    L._plot_rayleigh_vs_cloud(ray_k, cloud_k, OUT / "lindenberg_cl61_ray_vs_cloud.png")

    L.save_kalman_csv(ray_k, OUT / "rayleigh_lindenberg_cl61_kalman.csv")
    L.save_kalman_csv(cloud_k, OUT / "cloud_lindenberg_cl61_kalman.csv")

    # cross-check: median product of the two Kalman series on common dates
    if "wv" in ray_k and "wv" in cloud_k:
        gt_r, gs_r, _ = ray_k["wv"]
        gt_c, gs_c, _ = cloud_k["wv"]
        common = np.intersect1d(gt_r, gt_c)
        if common.size:
            sr = gs_r[np.isin(gt_r, common)]
            sc = gs_c[np.isin(gt_c, common)]
            print(f"  Cross-check: median CL_Ray = {np.median(sr):.3f}, "
                  f"median CL_Cloud = {np.median(sc):.3f}, "
                  f"median ratio = {np.median(sr / sc):.3f}")

    print("REPORT_FIGS_DONE.")


if __name__ == "__main__":
    main()
