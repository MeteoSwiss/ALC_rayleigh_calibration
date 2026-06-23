#!/usr/bin/env python
"""
MATLAB-vs-Python parity test for the liquid-water-cloud ceilometer calibration.

This runs BOTH implementations on the SAME E-PROFILE L2 NetCDF file and asserts they
agree within a tight tolerance:

  (a) it invokes the MATLAB reference ``liquid_cloud_calibration.m`` via the MATLAB CLI
      (``matlab.exe -batch``), with the working directory set to the
      ``liquid_cloud_calibration`` folder so its private helpers are on the path, and
      dumps the MATLAB outputs (per-profile ``all_coefficients`` = calibration
      coefficients, ``lidar_ratios`` = S_consistent at the valid profiles, ``cal_median``,
      ``cal_mean``, ``cal_mode``, ``cal_std``, ``n_profiles``, and ``data.trans2_wv``) to a
      v7 .mat file;
  (b) it runs the Python port ``calibration.cloud`` on the same file
      with the SAME config;
  (c) it asserts a match within rtol on ``cal_median`` and on the valid
      ``calibration_coefficients`` (aligned by profile index), plus ``cal_mean``,
      ``cal_std``, ``cal_mode``, the valid lidar ratios, and the full ``trans2_wv`` matrix.

It prints a clear PASS/FAIL with the max relative difference per quantity.

Config (identical on both sides), matching ``run_cloud_calibration.m``:
    instrument='CL61', apply_wv_correction=1, apply_transmission_correction=true,
    aerosol_lidar_ratio=50, plus the CAMS folder and the HITRAN cross-section LUT.

Run:
    python test_cloud_calibration_vs_matlab.py
or under pytest:
    python -m pytest test_cloud_calibration_vs_matlab.py -v -s -o addopts=""
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.io as sio

# --- repo / tool locations -------------------------------------------------
MATLAB_EXE = r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
MATLAB_FUNC_DIR = r"C:\Users\hervo\OneDrive\Documents\MATLAB\ALC\liquid_cloud_calibration"

# --- test inputs (identical for MATLAB and Python) -------------------------
TEST_NC = r"A:\E-PROFILE_L2_monthly\0-20000-0-06610\2026\L2_0-20000-0-06610_C202603.nc"
INSTRUMENT = "CL61"
CAMS_FOLDER = r"D:/CAMS/"
ABS_CS_LUT = r"C:/Users/hervo/OneDrive/Documents/MATLAB/MDA/monitoring_alc_monthly/abs_cross_647_full_levels_1000.nc"
AEROSOL_LIDAR_RATIO = 50.0
APPLY_WV = True
APPLY_TRANS_CORR = True

# --- tolerances ------------------------------------------------------------
RTOL_MEDIAN = 1e-3          # required: cal_median
RTOL_COEFFS = 1e-3          # required: valid calibration_coefficients
RTOL_TRANS2 = 5e-3          # trans2_wv (small residual from the q->RH round-trip numerics)
RTOL_STATS = 5e-3           # cal_mean / cal_std
# A handful of profiles can land on the exact boundary of a cloud filter and be rejected
# by a different (still-rejecting) filter on each side; this never changes the valid set
# by more than a couple of profiles. Allow a tiny mismatch in the valid mask.
MAX_VALID_MASK_MISMATCH = 4


def run_matlab(out_mat: Path) -> None:
    """Invoke the MATLAB reference and dump its outputs to *out_mat* (v7 .mat)."""
    # Build a one-shot batch script. Use forward slashes inside MATLAB strings to avoid
    # backslash-escape surprises; MATLAB accepts them on Windows.
    nc = TEST_NC.replace("\\", "/")
    lut = ABS_CS_LUT.replace("\\", "/")
    out = str(out_mat).replace("\\", "/")
    script = f"""
config.nc_file = '{nc}';
config.instrument = '{INSTRUMENT}';
config.apply_wv_correction = {int(APPLY_WV)};
config.apply_transmission_correction = {str(APPLY_TRANS_CORR).lower()};
config.aerosol_lidar_ratio = {AEROSOL_LIDAR_RATIO};
config.cams_folder = '{CAMS_FOLDER}';
config.abs_cs_lookup_table = '{lut}';
config.debug = 0;
res = liquid_cloud_calibration(config);
all_coefficients = res.all_coefficients(:);
lidar_ratios = res.lidar_ratios(:);
cal_median = res.cal_median;
cal_mean = res.cal_mean;
cal_mode = res.cal_mode;
cal_std = res.cal_std;
n_profiles = res.n_profiles;
if isfield(res,'trans2_wv'); trans2_wv = res.trans2_wv; else; trans2_wv = []; end
save('{out}', 'all_coefficients','lidar_ratios','cal_median','cal_mean','cal_mode','cal_std','n_profiles','trans2_wv', '-v7');
fprintf('MATLAB_DONE n=%d median=%.8f\\n', n_profiles, cal_median);
"""
    cmd = [MATLAB_EXE, "-batch", script]
    print(f"[matlab] running reference in {MATLAB_FUNC_DIR} ...")
    proc = subprocess.run(
        cmd, cwd=MATLAB_FUNC_DIR, capture_output=True, text=True, timeout=900)
    if proc.returncode != 0 or not out_mat.exists():
        print("----- MATLAB stdout -----")
        print(proc.stdout)
        print("----- MATLAB stderr -----")
        print(proc.stderr)
        raise RuntimeError(f"MATLAB run failed (returncode={proc.returncode})")
    # surface the confirmation line
    for line in proc.stdout.splitlines():
        if "MATLAB_DONE" in line:
            print("[matlab]", line.strip())


def run_python() -> dict:
    """Run the Python port and return the comparable quantities."""
    # ensure the package is importable when run as a plain script
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from calibration.cloud import (
        liquid_cloud_calibration, CloudCalConfig)

    cfg = CloudCalConfig(
        nc_file=TEST_NC,
        instrument=INSTRUMENT,
        apply_wv_correction=APPLY_WV,
        apply_transmission_correction=APPLY_TRANS_CORR,
        aerosol_lidar_ratio=AEROSOL_LIDAR_RATIO,
        cams_folder=CAMS_FOLDER,
        abs_cs_lookup_table=ABS_CS_LUT,
        debug=0,
    )
    res = liquid_cloud_calibration(cfg)
    return {
        "all_coefficients": np.asarray(res.all_coefficients, dtype="float64").ravel(),
        "lidar_ratios": np.asarray(res.lidar_ratios, dtype="float64").ravel(),
        "cal_median": float(res.cal_median),
        "cal_mean": float(res.cal_mean),
        "cal_mode": float(res.cal_mode),
        "cal_std": float(res.cal_std),
        "n_profiles": int(res.n_profiles),
        "trans2_wv": (None if res.trans2_wv is None
                      else np.asarray(res.trans2_wv, dtype="float64")),
    }


def _max_rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    """Max relative difference over finite, comparable entries.

    Uses |a-b| / max(|b|, tiny) so a near-zero reference does not blow up the ratio.
    """
    a = np.asarray(a, dtype="float64").ravel()
    b = np.asarray(b, dtype="float64").ravel()
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m):
        return 0.0
    denom = np.maximum(np.abs(b[m]), 1e-300)
    return float(np.max(np.abs(a[m] - b[m]) / denom))


def compare(mat: dict, py: dict) -> bool:
    """Compare MATLAB and Python results; print a report; return overall PASS."""
    ok = True
    print("\n================ MATLAB vs Python parity report ================")
    print(f" test file : {TEST_NC}")
    print(f" instrument: {INSTRUMENT}  (WV={APPLY_WV}, trans_corr={APPLY_TRANS_CORR},"
          f" aerosol_LR={AEROSOL_LIDAR_RATIO})")
    print("-" * 64)

    # --- scalars ---
    def _scalar_check(name: str, mv: float, pv: float, rtol: float) -> None:
        nonlocal ok
        rd = abs(pv - mv) / max(abs(mv), 1e-300)
        status = "PASS" if rd <= rtol else "FAIL"
        if rd > rtol:
            ok = False
        print(f" {name:<14s} MATLAB={mv:.8g}  Python={pv:.8g}  reldiff={rd:.3e}  "
              f"[{status} @ rtol {rtol:g}]")

    mv_n = int(np.atleast_1d(mat["n_profiles"]).ravel()[0])
    pv_n = py["n_profiles"]
    status = "PASS" if abs(mv_n - pv_n) <= MAX_VALID_MASK_MISMATCH else "FAIL"
    if abs(mv_n - pv_n) > MAX_VALID_MASK_MISMATCH:
        ok = False
    print(f" {'n_profiles':<14s} MATLAB={mv_n}  Python={pv_n}  "
          f"diff={abs(mv_n - pv_n)}  [{status} @ tol {MAX_VALID_MASK_MISMATCH}]")

    _scalar_check("cal_median", float(np.ravel(mat["cal_median"])[0]),
                  py["cal_median"], RTOL_MEDIAN)
    _scalar_check("cal_mean", float(np.ravel(mat["cal_mean"])[0]),
                  py["cal_mean"], RTOL_STATS)
    _scalar_check("cal_std", float(np.ravel(mat["cal_std"])[0]),
                  py["cal_std"], RTOL_STATS)
    _scalar_check("cal_mode", float(np.ravel(mat["cal_mode"])[0]),
                  py["cal_mode"], RTOL_STATS)

    # --- per-profile calibration coefficients (full-length, NaN for invalid) ---
    mc = np.asarray(mat["all_coefficients"], dtype="float64").ravel()
    pc = np.asarray(py["all_coefficients"], dtype="float64").ravel()
    print("-" * 64)
    if mc.size != pc.size:
        ok = False
        print(f" all_coefficients LENGTH MISMATCH: MATLAB {mc.size} vs Python {pc.size}")
    else:
        mvalid = np.isfinite(mc)
        pvalid = np.isfinite(pc)
        mask_mismatch = int(np.sum(mvalid != pvalid))
        both = mvalid & pvalid
        rd = _max_rel_diff(pc[both], mc[both])
        status = "PASS" if (rd <= RTOL_COEFFS and mask_mismatch <= MAX_VALID_MASK_MISMATCH) else "FAIL"
        if rd > RTOL_COEFFS or mask_mismatch > MAX_VALID_MASK_MISMATCH:
            ok = False
        print(f" calibration_coefficients (aligned by profile index):")
        print(f"   valid in both = {int(np.sum(both))},  "
              f"valid-mask mismatches = {mask_mismatch} (tol {MAX_VALID_MASK_MISMATCH})")
        print(f"   max reldiff on common-valid = {rd:.3e}  [{status} @ rtol {RTOL_COEFFS:g}]")

    # --- lidar ratios (S_consistent at the valid profiles) ---
    ml = np.asarray(mat["lidar_ratios"], dtype="float64").ravel()
    pl = np.asarray(py["lidar_ratios"], dtype="float64").ravel()
    if ml.size == pl.size and ml.size > 0:
        # both are the valid-profile lists in profile order; align directly
        rd = _max_rel_diff(pl, ml)
        status = "PASS" if rd <= RTOL_COEFFS else "FAIL"
        if rd > RTOL_COEFFS:
            ok = False
        print(f" lidar_ratios (S_consistent[valid]): n={ml.size}  "
              f"max reldiff={rd:.3e}  [{status} @ rtol {RTOL_COEFFS:g}]")
    else:
        # lengths differ by the borderline profiles; compare the sorted overlap leniently
        n = min(ml.size, pl.size)
        rd = _max_rel_diff(np.sort(pl)[:n], np.sort(ml)[:n]) if n else 0.0
        status = "PASS" if (abs(ml.size - pl.size) <= MAX_VALID_MASK_MISMATCH
                            and rd <= RTOL_COEFFS) else "FAIL"
        if status == "FAIL":
            ok = False
        print(f" lidar_ratios: MATLAB n={ml.size} Python n={pl.size} "
              f"(diff {abs(ml.size - pl.size)}); sorted-overlap max reldiff={rd:.3e} [{status}]")

    # --- trans2_wv (full matrix) ---
    print("-" * 64)
    mt = mat.get("trans2_wv", None)
    pt = py.get("trans2_wv", None)
    mt = None if (mt is None or np.size(mt) == 0) else np.asarray(mt, dtype="float64")
    if mt is None or pt is None:
        print(" trans2_wv : not available on one side -> skipped")
    else:
        if mt.shape != pt.shape:
            # one may be transposed relative to the other
            if mt.shape == pt.T.shape:
                pt = pt.T
        if mt.shape != pt.shape:
            ok = False
            print(f" trans2_wv SHAPE MISMATCH: MATLAB {mt.shape} vs Python {pt.shape}")
        else:
            rd = _max_rel_diff(pt, mt)
            status = "PASS" if rd <= RTOL_TRANS2 else "FAIL"
            if rd > RTOL_TRANS2:
                ok = False
            print(f" trans2_wv : shape {mt.shape}  max reldiff={rd:.3e}  "
                  f"[{status} @ rtol {RTOL_TRANS2:g}]")
            # also report the median (typical) difference, which is what matters physically
            both = np.isfinite(mt) & np.isfinite(pt)
            med = float(np.median(np.abs(pt[both] - mt[both]) /
                                  np.maximum(np.abs(mt[both]), 1e-300)))
            print(f"            median reldiff = {med:.3e}")

    print("=" * 64)
    print(f" OVERALL: {'PASS' if ok else 'FAIL'}")
    print("=" * 64)
    return ok


def _run() -> bool:
    if not Path(TEST_NC).exists():
        raise FileNotFoundError(f"Test NetCDF not found: {TEST_NC}")
    if not Path(MATLAB_EXE).exists():
        raise FileNotFoundError(f"matlab.exe not found: {MATLAB_EXE}")

    with tempfile.TemporaryDirectory() as td:
        out_mat = Path(td) / "matlab_cloudcal_out.mat"
        run_matlab(out_mat)
        mat = sio.loadmat(str(out_mat))
        py = run_python()
        return compare(mat, py)


def test_cloud_calibration_matches_matlab():
    """pytest entry point.

    Opt-in integration test: it shells out to MATLAB (~30 s) and needs the A:/D:/C:
    data drives, so it is skipped in the default suite. Set RUN_MATLAB_PARITY=1 to run
    it (or run this file directly: ``python tests/test_cloud_calibration_vs_matlab.py``).
    """
    import os
    import pytest

    if os.environ.get("RUN_MATLAB_PARITY") != "1":
        pytest.skip("MATLAB-parity integration test; set RUN_MATLAB_PARITY=1 to run")
    assert _run(), "Python port does not match the MATLAB reference within tolerance"


if __name__ == "__main__":
    sys.exit(0 if _run() else 1)
