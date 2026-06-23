"""
Tests for the cloud-calibration aerosol two-way transmission correction
(``apply_transmission_correction``) and its scale-bug fix.

The bug: the below-cloud aerosol optical depth was computed from the *uncalibrated* stored
attenuated backscatter, ``B_aerosol = sum(beta)*range_resol``. For L2 input (stored in
``1E-6*1/(m*sr)`` units, ~O(1)) this overshoots by ~1e6, so ``T2 = exp(-2*LR*B_aerosol)``
underflowed to 0 and the median calibration coefficient collapsed to 0. The fix converts the
integral to a *physical* optical depth with the per-profile coefficient ``C = S/S_THEORETICAL``:
``AOD = LR * C * sum(beta)*range_resol`` — which is correct for ANY beta scale (L1 physical or
L2 stored), because ``C`` absorbs the scale.

Two layers of tests:
  * a fast, deterministic UNIT test of ``apply_transmission_correction`` at both the L1 (physical)
    and L2 (stored, ×1e6) beta scales — always runs;
  * end-to-end INTEGRATION tests of ``liquid_cloud_calibration`` on real L1 and L2 CL61 files
    (skipped when the CAMS files needed for the mandatory 910 nm water-vapour correction are
    absent, e.g. in CI).
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from calibration.cloud.calibration import (
    apply_transmission_correction, S_THEORETICAL, CloudCalConfig, set_defaults)

REPO = Path(__file__).resolve().parents[1]
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"   # bundled full-spectral LUT
CAMS_DIR = Path("D:/CAMS")
# Sample CL61 files bundled in the repo, plus a real L1 (cloudy day) from the 2026 archive.
BUNDLED_L2 = REPO / "examples/data/L2/0-756-4-EERLCL61/2026/04/L2_0-756-4-EERLCL61_A20260414.nc"
BUNDLED_L1 = REPO / "examples/data/L1/0-756-4-EERLCL61/2026/03/L1_0-756-4-EERLCL61_A20260304.nc"
REAL_L1 = Path("D:/E-PROFILE_L1_2026/0-756-4-EERLCL61/2026/04/L1_0-756-4-EERLCL61_A20260414.nc")
SION = dict(station_latitude=46.2246, station_longitude=7.3607)


def _cams_for(yyyymm: str) -> bool:
    return (CAMS_DIR / f"CAMS_Beta_{yyyymm}.nc").is_file()


def _synthetic_profile(scale: float, n_range: int = 100, dz: float = 30.0):
    """A single-layer below-cloud aerosol + an opaque cloud peak, range x 1 profile.

    ``scale`` multiplies the physical backscatter (scale=1 -> L1 physical ~1/(m*sr);
    scale=1e6 -> L2 stored in 1E-6 units). The matching apparent lidar ratio S is divided by
    the same scale, so the calibration coefficient C=S/S_THEORETICAL maps the scaled beta back
    to physical — i.e. the *physical* atmosphere is identical at both scales.
    """
    rng = np.arange(n_range, dtype="float64") * dz
    beta_phys = np.zeros(n_range)
    beta_phys[4:36] = 2.0e-6        # aerosol below the cloud (physical 1/(m*sr))
    beta_phys[40] = 1.0e-3          # opaque liquid-cloud peak (argmax -> max_idx = 40)
    beta = (beta_phys * scale)[:, None]                 # (n_range, 1)
    S = np.array([S_THEORETICAL / scale])               # so C_base = S/S_THEORETICAL = 1/scale
    data = SimpleNamespace(range=rng)
    cfg = set_defaults(CloudCalConfig(instrument="CL61", aerosol_lidar_ratio=50.0))
    return beta, data, S, cfg


# ---------------------------------------------------------------------------
# Unit test of the fixed function (always runs)
# ---------------------------------------------------------------------------
def _apply_tc(scale: float):
    """Run apply_transmission_correction on the synthetic profile -> (C_corrected, C_base, T2)."""
    beta, data, S, cfg = _synthetic_profile(scale)
    C_corrected, C_low, C_high = apply_transmission_correction(beta, data, S, cfg)
    C_base = float(S[0] / S_THEORETICAL)
    c = float(C_corrected[0])
    T2 = c / C_base
    return c, C_base, T2, float(C_low[0]), float(C_high[0])


@pytest.mark.parametrize("scale,label", [(1.0, "L1-physical"), (1.0e6, "L2-stored")])
def test_transmission_correction_meaningful(scale, label):
    c, C_base, T2, c_low, c_high = _apply_tc(scale)
    # 1) meaningful: finite, positive, NOT collapsed to zero (the bug)
    assert np.isfinite(c) and c > 0, f"{label}: C_corrected not positive-finite ({c})"
    # 2) the correction is a real two-way transmission in (0, 1]
    assert 0.0 < T2 <= 1.0, f"{label}: implied T2={T2} outside (0,1]"
    # 3) the implied below-cloud AOD is physical (~0.1 here), NOT the ~1e4 of the bug
    aod = -0.5 * np.log(T2)
    assert 1e-3 < aod < 1.0, f"{label}: implied AOD={aod} not physical"
    # low LR -> less attenuation -> larger C; high LR -> smaller C
    assert c_high <= c <= c_low


def test_transmission_correction_is_scale_invariant():
    """The relative correction (T2) must be identical at the L1 and L2 beta scales — i.e. the
    fix gives the same physical answer whether the input beta is L1-physical or L2-stored."""
    _, _, t2_l1, _, _ = _apply_tc(1.0)
    _, _, t2_l2, _, _ = _apply_tc(1.0e6)
    assert t2_l1 == pytest.approx(t2_l2, rel=1e-9), (t2_l1, t2_l2)
    assert 0.5 < t2_l1 < 0.95     # a clearly non-trivial but moderate correction


def test_regression_old_formula_underflows_on_l2_scale():
    """Guard the specific bug: the OLD integral (no C factor) underflows T2->0 at the L2 scale,
    while the fixed code keeps it finite."""
    beta, data, S, cfg = _synthetic_profile(1.0e6)        # L2 stored scale
    range_resol = float(data.range[1] - data.range[0])
    seg = beta[4:(40 - 5) + 1, 0]
    B_old = np.nansum(seg) * range_resol                  # OLD: uncalibrated integral (~960)
    T2_old = np.exp(-2.0 * cfg.aerosol_lidar_ratio * B_old)
    assert T2_old == 0.0, f"expected the old formula to underflow, got T2={T2_old}"
    # the fixed function does NOT underflow on the same input
    C_corrected, _, _ = apply_transmission_correction(beta, data, S, cfg)
    assert C_corrected[0] > 0


# ---------------------------------------------------------------------------
# End-to-end integration tests on real L1 / L2 CL61 (need CAMS for the WV correction)
# ---------------------------------------------------------------------------
def _run_cloud(nc_file: str, tc: bool):
    cfg = CloudCalConfig(
        nc_file=str(nc_file), instrument="CL61", apply_wv_correction=True,
        apply_transmission_correction=tc, aerosol_lidar_ratio=50.0,
        cams_folder=str(CAMS_DIR), abs_cs_lookup_table=str(WV_LUT), **SION)
    from calibration.cloud import liquid_cloud_calibration
    return liquid_cloud_calibration(cfg)


@pytest.mark.skipif(not (BUNDLED_L2.is_file() and WV_LUT.is_file() and _cams_for("202604")),
                    reason="needs the bundled L2 sample, the WV LUT and CAMS_Beta_202604.nc")
def test_l2_cloud_calibration_with_transmission_correction():
    """L2 input: the transmission correction must produce a meaningful (non-zero) coefficient,
    slightly below the uncorrected one."""
    on = _run_cloud(BUNDLED_L2, tc=True)
    off = _run_cloud(BUNDLED_L2, tc=False)
    assert on.n_profiles > 0
    assert np.isfinite(on.cal_median) and on.cal_median > 0
    ratio = on.cal_median / off.cal_median
    assert 0.80 <= ratio <= 1.0, f"transmission correction effect off-range: {ratio:.3f}"


@pytest.mark.skipif(not (REAL_L1.is_file() and WV_LUT.is_file() and _cams_for("202604")),
                    reason="needs the L1 2026 archive, the WV LUT and CAMS_Beta_202604.nc")
def test_l1_cloud_calibration_with_transmission_correction():
    """L1 input (physical 1/(m*sr) beta, ~O(1) coefficient): same meaningful, non-zero result."""
    on = _run_cloud(REAL_L1, tc=True)
    off = _run_cloud(REAL_L1, tc=False)
    assert on.n_profiles > 0
    assert np.isfinite(on.cal_median) and on.cal_median > 0
    ratio = on.cal_median / off.cal_median
    assert 0.80 <= ratio <= 1.0, f"transmission correction effect off-range: {ratio:.3f}"


@pytest.mark.skipif(not (BUNDLED_L1.is_file() and WV_LUT.is_file() and _cams_for("202603")),
                    reason="needs the bundled L1 sample, the WV LUT and CAMS_Beta_202603.nc")
def test_l1_input_is_read_without_error():
    """The reader must accept the L1 physical units ('m^-1.sr^-1') and run end-to-end even on a
    cloudless day (n_profiles may be 0), without raising or crashing."""
    res = _run_cloud(BUNDLED_L1, tc=True)
    assert res.n_profiles >= 0          # ran to completion; cloudless day -> 0 is acceptable
