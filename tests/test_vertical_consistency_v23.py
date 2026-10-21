"""eprof_v2.3: the repaired fallback tier of eprof_v2.2 (q19; doc/v2.3/vertical_consistency.md).

* eprof_v2.2 is unchanged: its parameters, and every night its strict tier calibrates (tier 1 never
  reads the noise), which come out bit-identical under eprof_v2.3;
* the night-mean noise measured on the binned stack recovers sigma / sqrt(n) on white noise;
* the vertical-consistency test accepts a clean column and a clean window under a brighter elevated
  layer, and rejects aerosol in the window and a window sitting in a decaying tail; when every tier-2
  candidate fails, the night fails instead of keeping an aerosol-loaded window;
* instrument scope: measured noise for CHM15k and CL61, the test for CHM15k only, the v2.2 fallback
  tier for every other type (an explicit molecular_params entry wins).
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType, DataLevel
from calibration.io.output import version_code
from calibration.rayleigh import calibration as C
from calibration.rayleigh import molecular_methods as M

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "examples" / "data"
OPTIONS = REPO / "options.json"

# Fixture nights the strict tier calibrates, one per instrument family.
TIER1_NIGHTS = [
    ("0-20000-0-06610", "20260225", InstrumentType.CHM15k),
    ("0-20000-0-07014", "20260424", InstrumentType.MINI_MPL),
    ("0-756-4-EERLCL61", "20260304", InstrumentType.CL61),
]
FORCE_TIER2 = {"min_r2": 1.01}          # no window reaches R2 > 1: every night goes to the fallback tier


def _info(wmo, itype, date):
    f = DATA / "L2" / wmo / date[:4] / date[4:6] / f"L2_{wmo}_A{date}.nc"
    with Dataset(f) as ds:
        lat = float(ds.variables["station_latitude"][...])
        lon = float(ds.variables["station_longitude"][...])
        alt = float(ds.variables["station_altitude"][...])
    return InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier="A", instrument_type=itype,
                          latitude=lat, longitude=lon, altitude=alt)


def _options(method, **params):
    o = CalibrationOptions.from_json(OPTIONS)
    o.data_level = DataLevel.L2_DAILY
    o.folder_root = DATA / "L2"
    o.folder_output = Path(tempfile.mkdtemp())
    o.plot_main = False
    o.molecular_method = method
    o.molecular_params = dict(params)
    return o


# --------------------------------------------------------------------------------------------
# Registry: v2.2 untouched, v2.3 = v2.2 + the q19 switches
# --------------------------------------------------------------------------------------------
def test_v22_preset_untouched_and_v23_extends_it():
    v22, v23 = M.DEFAULT_PARAMS["eprof_v2.2"], M.DEFAULT_PARAMS["eprof_v2.3"]
    assert "use_vc" not in v22 and "measured_noise" not in v22
    assert {k: v23[k] for k in v22} == v22
    assert v23["measured_noise"] is True and v23["use_vc"] is True
    assert (v23["vc_k"], v23["vc_eps_tol"], v23["vc_tol_lo"], v23["vc_block_m"]) == (2.0, 0.1, 0.1, 750.0)
    assert "eprof_v2.3" in M.METHODS
    assert version_code("eprof_v2.3") == 230


# --------------------------------------------------------------------------------------------
# The two building blocks on synthetic data
# --------------------------------------------------------------------------------------------
def test_binned_mean_noise_recovers_sigma_over_sqrt_n():
    rng = np.random.default_rng(0)
    stack = 5.0 + rng.normal(0.0, 1.0, size=(64, 400))      # 64 binned profiles, unit noise
    sigma = M._binned_mean_noise(stack)
    assert np.nanmedian(sigma) == pytest.approx(1.0 / np.sqrt(64), rel=0.10)
    assert M._binned_mean_noise(stack[:5]) is None           # too short to measure


Z = np.arange(30.0, 12000.0, 30.0)
P_MOL = np.exp(-Z / 8000.0) / (Z / 1000.0) ** 2              # any smooth positive molecular profile
SIGMA = 0.002 * P_MOL                                        # 0.2 % per gate: measured up to the 9 km cap


def _grid(ratio, sigma=SIGMA):
    """Window grid of the profile signal = ratio(z) * p_mol (ratio = scattering ratio up to C_L)."""
    return M.compute_window_grid(ratio * P_MOL, P_MOL, Z, (300.0, 600.0), range_start_m=2000.0,
                                 range_end_m=6000.0, sigma_signal=sigma)


def _cell(g, center_m, half_m):
    return int(np.argmin(np.abs(g.center_m - center_m))), int(np.argmin(np.abs(g.half_m - half_m)))


def test_vertical_consistency_verdicts():
    clean = np.ones_like(Z)
    g = _grid(clean)
    assert M._vertical_consistency(g, *_cell(g, 3800, 600)) == (True, "accept")

    # aerosol in the window (R 1.3 over 3.0-4.7 km, clean above): the window is brighter than the column
    layer = np.where((Z >= 3000) & (Z <= 4700), 1.3, 1.0)
    g = _grid(layer)
    assert M._vertical_consistency(g, *_cell(g, 3800, 600)) == (False, "aerosol")

    # a window just above a decaying tail: the 1 km below it is much brighter
    tail = np.where(Z >= 3000, 1.0 + 4.0 * np.exp(-(Z - 3000.0) / 400.0), 5.0)
    g = _grid(tail)
    assert M._vertical_consistency(g, *_cell(g, 4800, 300)) == (False, "tail")

    # a BRIGHTER layer above a clean window (cirrus, an elevated plume) never rejects it
    cirrus = np.where((Z >= 7000) & (Z <= 7600), 3.0, 1.0)
    g = _grid(cirrus)
    assert M._vertical_consistency(g, *_cell(g, 3800, 600)) == (True, "accept")

    # nothing measured -> behave exactly as before
    g = M.compute_window_grid(clean * P_MOL, P_MOL, Z, (300.0, 600.0), range_start_m=2000.0,
                              range_end_m=6000.0)
    assert g.sigma is None
    assert M._vertical_consistency(g, *_cell(g, 3800, 600)) == (True, "no-noise")


def test_fallback_tier_fails_when_every_candidate_is_aerosol_loaded():
    # a deep layer (R 1.3 up to 6.8 km) covers the whole search range; the clean air is above it
    g = _grid(np.where(Z < 6800, 1.3, 1.0))
    v22 = M._select_optimal(g, **dict(M.DEFAULT_PARAMS["eprof_v2.2"], **FORCE_TIER2))
    v23 = M._select_optimal(g, **dict(M.DEFAULT_PARAMS["eprof_v2.3"], **FORCE_TIER2))
    assert v22.ok                                             # v2.2 rescues a window inside the layer
    assert not v23.ok
    assert "vertical-consistency" in v23.message and "first candidate: aerosol" in v23.message


def test_fallback_tier_keeps_a_consistent_window_above_a_layer():
    # aerosol below 3.5 km, clean above: the kept tier-2 window passes the test, above the layer
    g = _grid(np.where(Z < 3500, 1.3, 1.0))
    mw = M._select_optimal(g, **dict(M.DEFAULT_PARAMS["eprof_v2.3"], **FORCE_TIER2))
    assert mw.ok and "vertical consistency" in mw.message
    assert mw.start_m >= 3500
    assert M._vertical_consistency(g, *_cell(g, mw.center_m, mw.half_m)) == (True, "accept")


# --------------------------------------------------------------------------------------------
# Instrument scope (calibrate_rayleigh._window_search_setup)
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("itype,noise,vc", [
    (InstrumentType.CHM15k, True, True),
    (InstrumentType.CL61, True, False),
    (InstrumentType.CL31, False, False),
    (InstrumentType.CL51, False, False),
    (InstrumentType.MINI_MPL, False, False),
], ids=lambda x: getattr(x, "value", str(x)))
def test_scope_by_instrument_type(itype, noise, vc):
    o = _options("eprof_v2.3")
    sigma = np.full(10, 4.0)
    params, sigma_search = C._window_search_setup(o, itype, sigma, native_per_bin=16.0)
    merged = dict(M.DEFAULT_PARAMS["eprof_v2.3"], **(params or {}))
    assert merged["measured_noise"] is noise
    assert merged["use_vc"] is vc
    # the propagated sigma is corrected for the native profiles in each bin only where it is used
    assert np.allclose(sigma_search, 1.0 if noise else 4.0)
    assert np.allclose(sigma, 4.0)                           # the caller's sigma_signal is not modified


def test_scope_leaves_other_methods_and_explicit_params_alone():
    sigma = np.full(10, 4.0)
    o = _options("eprof_v2.2", max_chi2red=3.0)
    params, sigma_search = C._window_search_setup(o, InstrumentType.CHM15k, sigma, 16.0)
    assert params == {"max_chi2red": 3.0} and sigma_search is sigma
    o = _options("eprof_v2.3", use_vc=True)                  # an explicit entry wins over the scope
    params, _ = C._window_search_setup(o, InstrumentType.CL61, sigma, 16.0)
    assert params["use_vc"] is True


# --------------------------------------------------------------------------------------------
# End to end on the fixtures
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("wmo,date,itype", TIER1_NIGHTS, ids=["CHM15k", "Mini-MPL", "CL61"])
def test_tier1_nights_bit_identical(wmo, date, itype):
    results = {}
    for method in ("eprof_v2.2", "eprof_v2.3"):
        o = _options(method)
        results[method] = calibrate_rayleigh(date, _info(wmo, itype, date), o)
        files = list(Path(o.folder_output).rglob("*.nc"))
        with Dataset(files[0]) as ds:
            assert int(ds["calibration_version"][-1]) == version_code(method)
    r22, r23 = results["eprof_v2.2"], results["eprof_v2.3"]
    assert r22.flag == 1 and r23.flag == 1
    assert r23.lidar_constant == r22.lidar_constant
    assert r23.uncertainty == r22.uncertainty


@pytest.mark.parametrize("wmo,date,itype,noise,vc", [
    ("0-20000-0-06610", "20260225", InstrumentType.CHM15k, True, True),
    ("0-756-4-EERLCL61", "20260304", InstrumentType.CL61, True, False),
    ("0-20000-0-07014", "20260424", InstrumentType.MINI_MPL, False, False),
], ids=["CHM15k", "CL61", "Mini-MPL"])
def test_forced_tier2_runs_only_the_validated_parts(wmo, date, itype, noise, vc, monkeypatch):
    calls = {"noise": 0, "vc": 0}
    measure_noise = M._binned_mean_noise
    check_column = M._vertical_consistency

    def counted_noise(*args, **kwargs):
        calls["noise"] += 1
        return measure_noise(*args, **kwargs)

    def counted_vc(*args, **kwargs):
        calls["vc"] += 1
        return check_column(*args, **kwargs)

    monkeypatch.setattr(M, "_binned_mean_noise", counted_noise)
    monkeypatch.setattr(M, "_vertical_consistency", counted_vc)
    r23 = calibrate_rayleigh(date, _info(wmo, itype, date), _options("eprof_v2.3", **FORCE_TIER2))
    assert (calls["noise"] > 0) is noise
    assert (calls["vc"] > 0) is vc
    assert r23.flag == 1
    if not noise:                         # out of scope: exactly the eprof_v2.2 fallback tier
        r22 = calibrate_rayleigh(date, _info(wmo, itype, date), _options("eprof_v2.2", **FORCE_TIER2))
        assert r23.lidar_constant == r22.lidar_constant
