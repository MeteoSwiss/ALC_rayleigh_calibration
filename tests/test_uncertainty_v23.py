"""v2.3 uncertainty budget (calibration/rayleigh/uncertainty.py).

* option OFF (default 'ensemble'): the v2.2 result and NetCDF are untouched -- the v2.3 code is not
  even entered, the constant/uncertainty equal the pre-v2.3 snapshot and no new variable is written;
* option ON: same constant, a larger uncertainty, finite non-negative components, the components in
  the NetCDF, and the relative QC gate;
* a synthetic molecular profile with an injected range-dependent distortion (additive B r^2 on the
  range-corrected signal, i.e. a constant electronic offset) gives a much larger vertical term than
  the clean profile.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType, DataLevel
from calibration.rayleigh import uncertainty as U
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.rayleigh.rayleigh_fit import RayleighFitResult

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "examples" / "data"
OPTIONS = REPO / "options.json"

# Self-contained Rayleigh fixtures (no CAMS needed) and the pre-v2.3 (git 35e7c23) results on them.
FIXTURES = [
    ("0-20000-0-06610", "20260225", InstrumentType.CHM15k, 492307124267.0247, 31669922985.400066),
    ("0-20000-0-07014", "20260424", InstrumentType.MINI_MPL, 262415.5271364225, 16687.475435226414),
]
V23_VARS = ("lidar_constant_uncertainty_ensemble", "lidar_constant_uncertainty_vertical",
            "lidar_constant_uncertainty_molecular", "lidar_constant_uncertainty_aerosol",
            "lidar_constant_uncertainty_noise", "lidar_constant_height_gradient")


def _info(wmo, itype, date):
    f = DATA / "L2" / wmo / "2026" / date[4:6] / f"L2_{wmo}_A{date}.nc"
    with Dataset(f) as ds:
        lat = float(ds.variables["station_latitude"][...])
        lon = float(ds.variables["station_longitude"][...])
        alt = float(ds.variables["station_altitude"][...])
    return InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier="A", instrument_type=itype,
                          latitude=lat, longitude=lon, altitude=alt)


def _options(model=None, **params):
    o = CalibrationOptions.from_json(OPTIONS)
    o.data_level = DataLevel.L2_DAILY
    o.folder_root = DATA / "L2"
    o.folder_output = Path(tempfile.mkdtemp())
    o.plot_main = False
    if model is not None:
        o.uncertainty_model = model
    if params:
        o.uncertainty_params = dict(params)
    return o


def _written(o):
    files = list(Path(o.folder_output).rglob("*.nc"))
    assert len(files) == 1
    return files[0]


# --------------------------------------------------------------------------------------------
# Option OFF: bit-identical to v2.2
# --------------------------------------------------------------------------------------------
def test_default_option_is_ensemble():
    assert CalibrationOptions().uncertainty_model == "ensemble"
    assert CalibrationOptions.from_json(OPTIONS).uncertainty_model == "ensemble"
    assert CalibrationOptions.from_json(OPTIONS).kalman_uncertainty_weighting is False


@pytest.mark.parametrize("wmo,date,itype,c_ref,u_ref", FIXTURES, ids=["CHM15k", "Mini-MPL"])
def test_option_off_is_bit_identical(wmo, date, itype, c_ref, u_ref, monkeypatch):
    # the v2.3 hook must not even be entered with the default option
    def _boom(**_):
        raise AssertionError("v2.3 budget computed with uncertainty_model='ensemble'")
    monkeypatch.setattr(U, "budget_for_night", _boom)
    o = _options()
    r = calibrate_rayleigh(date, _info(wmo, itype, date), o)
    assert r.flag in (1, 0.5)
    # pre-v2.3 snapshot (exact on the reference machine; rel 1e-12 tolerates BLAS/platform rounding)
    assert r.lidar_constant == pytest.approx(c_ref, rel=1e-12)
    assert r.uncertainty == pytest.approx(u_ref, rel=1e-12)
    assert r.uncertainty_components is None
    assert "v2.3" not in r.message
    with Dataset(_written(o)) as ds:
        assert not any(v in ds.variables for v in V23_VARS)
        unc = ds["lidar_constant_uncertainty"]
        assert unc.comments == "(np.std(CL_all) + error_fit * np.median(CL_all))*2"
        assert "coverage_factor" not in unc.ncattrs()


# --------------------------------------------------------------------------------------------
# Option ON
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("wmo,date,itype,c_ref,u_ref", FIXTURES, ids=["CHM15k", "Mini-MPL"])
def test_option_on_same_constant_full_budget(wmo, date, itype, c_ref, u_ref):
    o = _options("v2.3")
    fi = {}
    r = calibrate_rayleigh(date, _info(wmo, itype, date), o, fit_inputs_out=fi)
    assert r.flag in (1, 0.5)
    assert r.lidar_constant == pytest.approx(c_ref, rel=1e-12)      # C_L unchanged
    b = fi["uncertainty_budget"]
    assert b.u_ensemble == pytest.approx(u_ref, rel=1e-12)          # (a) = the v2.2 value
    for name in ("sigma_ensemble_rel", "sigma_vertical_rel", "sigma_molecular_rel",
                 "sigma_aerosol_rel", "sigma_noise_rel"):
        v = getattr(b, name)
        assert np.isfinite(v) and v >= 0, name
    # the residual-aerosol allowance scales with the wavelength (532 nm: 1/8 of the 1064 nm value)
    assert b.sigma_aerosol_rel == pytest.approx(0.02 * (itype.wavelength_nm / 1064.0) ** 3)
    assert r.uncertainty == pytest.approx(b.u_total)
    assert r.uncertainty >= b.u_ensemble
    total = 2 * r.lidar_constant * np.sqrt(b.sigma_ensemble_rel ** 2 + b.sigma_vertical_rel ** 2
                                           + b.sigma_molecular_rel ** 2 + b.sigma_aerosol_rel ** 2
                                           + b.sigma_noise_rel ** 2)
    assert r.uncertainty == pytest.approx(total, rel=1e-9)
    assert b.n_subwindows >= 3 and np.isfinite(b.gradient_per_km)
    with Dataset(_written(o)) as ds:
        for v in V23_VARS:
            assert v in ds.variables, v
        assert ds["lidar_constant_uncertainty"].coverage_factor == 2
        assert ds["lidar_constant_uncertainty_ensemble"][0] == pytest.approx(b.u_ensemble)
        assert ds["lidar_constant_uncertainty_vertical"][0] == pytest.approx(b.u_vertical)
        assert ds["lidar_constant_uncertainty"][0] == pytest.approx(b.u_total)


def test_relative_qc_gate():
    wmo, date, itype = FIXTURES[0][:3]
    r = calibrate_rayleigh(date, _info(wmo, itype, date), _options("v2.3", max_rel_uncertainty=0.01))
    assert r.flag == -6
    assert "Relative uncertainty exceeds limit" in r.message
    r = calibrate_rayleigh(date, _info(wmo, itype, date), _options("v2.3", max_rel_uncertainty=10.0))
    assert r.flag in (1, 0.5)


def test_netcdf_rows_mixing_v22_and_v23(tmp_path):
    """A file can hold v2.2 and v2.3 rows: the components exist only once a v2.3 row is written,
    a v2.2 row (re)written later carries fill values (never stale components), and
    strip_calibration_method keeps the extra variables when it rewrites the file."""
    from calibration.config import CalibrationResult
    from calibration.io.output import write_calibration_result, strip_calibration_method
    info = InstrumentInfo(site_name="X", wmo_id="0-0-0-TEST", identifier="A",
                          instrument_type=InstrumentType.CHM15k, latitude=47.0, longitude=7.0,
                          altitude=500.0)
    hk = {}
    comps = {v: float(i + 1) for i, v in enumerate(V23_VARS)}

    def write(day, components):
        r = CalibrationResult(lidar_constant=1e12, flag=1, uncertainty=5e10,
                              uncertainty_components=components)
        return write_calibration_result(tmp_path, info, r, date_epoch=20000.0 + day,
                                        time_start=20000.0 + day - 0.3,
                                        time_end=20000.0 + day + 0.2, wavelength_nm=1064.0,
                                        housekeeping=hk)

    f = write(0, None)                                        # v2.2 row: plain file
    with Dataset(f) as ds:
        assert not any(v in ds.variables for v in V23_VARS)
    write(1, comps)                                           # v2.3 row: variables appear
    write(0, None)                                            # rewrite the v2.2 row
    with Dataset(f) as ds:
        ds.set_auto_mask(False)
        for v in V23_VARS:
            vals = ds[v][:]
            assert vals[0] == pytest.approx(-999.9) and vals[1] == pytest.approx(comps[v])
        assert ds["lidar_constant_uncertainty"].coverage_factor == 2
    removed = strip_calibration_method(tmp_path, info, method=0, t_start=19999.5, t_end=20000.5)
    assert removed == 1
    with Dataset(f) as ds:
        ds.set_auto_mask(False)
        assert len(ds["time"]) == 1
        for v in V23_VARS:
            assert ds[v][0] == pytest.approx(comps[v])


def test_budget_failure_never_loses_the_night(monkeypatch):
    def _boom(**_):
        raise RuntimeError("synthetic failure")
    monkeypatch.setattr(U, "uncertainty_budget", _boom)
    wmo, date, itype, c_ref, u_ref = FIXTURES[0]
    fi = {}
    r = calibrate_rayleigh(date, _info(wmo, itype, date), _options("v2.3"), fit_inputs_out=fi)
    assert r.flag in (1, 0.5)
    assert fi["uncertainty_budget"].model == "v2.3-fallback"
    assert r.uncertainty == pytest.approx(u_ref, rel=1e-12)       # the ensemble value


def test_unknown_param_is_ignored():
    o = CalibrationOptions()
    o.uncertainty_params = {"not_a_param": 1, "window_displacement_m": 750.0}
    p = U.resolve_params(o)
    assert p["window_displacement_m"] == 750.0 and "not_a_param" not in p


# --------------------------------------------------------------------------------------------
# Synthetic profiles: an electronic distortion must inflate the vertical term
# --------------------------------------------------------------------------------------------
def _synthetic(distortion=0.0, noise_rel_5km=0.0, seed=0):
    """Pure molecular CHM15k-like profile on a 30 m grid (C_L = 1e12), optional additive B r^2."""
    rng = np.arange(15.0, 12000.0, 30.0)
    atm = load_standard_atmosphere(None, rng + 500.0)
    mp = calculate_molecular_properties(atm.temperature, atm.pressure, rng, 1064e-9)
    c_true = 1e12
    rcs = c_true * mp.beta_mol * mp.transmission                    # range-corrected signal
    signal = rcs / rng ** 2
    sigma = None
    if noise_rel_5km > 0:
        i5 = int(np.argmin(np.abs(rng - 5000.0)))
        sigma = np.full(rng.size, noise_rel_5km * signal[i5])        # background-limited noise
        rcs = rcs + np.random.default_rng(seed).normal(0.0, 1.0, rng.size) * sigma * rng ** 2
    if distortion:
        # constant offset in the background-subtracted power: rcs += B r^2, B = distortion x the
        # molecular power at 5 km (a 5 % offset there grows to ~20 % at 7 km)
        i5 = int(np.argmin(np.abs(rng - 5000.0)))
        rcs = rcs + distortion * signal[i5] * rng ** 2
    fit = RayleighFitResult(slope=c_true, intercept=0.0, r_squared=1.0, std_error=0.0,
                            p_value=0.0, center_range_m=3500.0, half_length_m=500.0,
                            range_start_m=3000.0, range_end_m=4000.0, altitude_start=3500.0,
                            altitude_end=4500.0, relative_error=0.0)
    return rng, rcs, mp, fit, sigma, c_true


def _budget(rng, rcs, mp, fit, sigma, cl):
    o = CalibrationOptions()
    o.uncertainty_model = "v2.3"
    return U.uncertainty_budget(
        cl=cl, u_ensemble=0.02 * cl, rcs_mean=rcs, range_alc=rng, p_mol=mp.p_mol,
        fit_result=fit, sigma_signal=sigma,
        cl_in_window=U.operational_cl_in_window(rcs, rng, mp.beta_mol, fit, o),
        options=o, latitude=47.0, station_altitude=500.0, date_str="20250415")


def test_synthetic_clean_profile_is_vertically_consistent():
    rng, rcs, mp, fit, sigma, c = _synthetic()
    b = _budget(rng, rcs, mp, fit, sigma, c)
    assert b.n_subwindows >= 3
    # every sub-window recovers the true constant
    assert np.allclose(b.sub_cl / c, 1.0, rtol=2e-3)
    assert abs(b.gradient_per_km) < 2e-3
    assert b.sigma_vertical_rel < 2e-3


def test_synthetic_distortion_inflates_vertical_term():
    rng, rcs, mp, fit, sigma, c = _synthetic(noise_rel_5km=0.05, seed=1)
    clean = _budget(rng, rcs, mp, fit, sigma, c)
    rng, rcs_d, mp, fit, sigma, c = _synthetic(distortion=0.05, noise_rel_5km=0.05, seed=1)
    dist = _budget(rng, rcs_d, mp, fit, sigma, c)
    assert np.isfinite(clean.sigma_vertical_rel) and np.isfinite(dist.sigma_vertical_rel)
    # positive offset -> C_L grows with height -> positive gradient, much larger vertical term
    # (the offset is 5 % of the molecular signal at 5 km: +1 % at 2 km, +4 % at 4.5 km)
    assert dist.gradient_per_km > 0.008
    assert dist.sigma_vertical_rel > 5.0 * max(clean.sigma_vertical_rel, 1e-3)
    # a negative offset (over-subtracted background) gives the mirror image
    rng, rcs_n, mp, fit, sigma, c = _synthetic(distortion=-0.05, noise_rel_5km=0.05, seed=1)
    neg = _budget(rng, rcs_n, mp, fit, sigma, c)
    assert neg.gradient_per_km < -0.008
    assert neg.sigma_vertical_rel > 5.0 * max(clean.sigma_vertical_rel, 1e-3)
    assert dist.u_total > clean.u_total
    for b in (clean, dist):
        for v in (b.sigma_vertical_rel, b.sigma_noise_rel, b.sigma_molecular_rel, b.u_total):
            assert np.isfinite(v) and v >= 0


def test_fit_log_gradient_recovers_slope():
    z = np.arange(2000.0, 6001.0, 500.0)
    cl = 1e12 * np.exp(-0.10 * (z - 4000.0) / 1000.0)                # -10 %/km
    f = U.fit_log_gradient(z, cl, np.full(z.size, 0.01), 4000.0, noise_floor_rel=0.0)
    assert f["g"] == pytest.approx(-0.10, rel=1e-9)
    assert f["resid_rel"] == pytest.approx(0.0, abs=1e-9)


def test_runner_kalman_weighting_option(monkeypatch):
    """_kalman_rows: unweighted by default (operational); weighted on request, and weight-neutral
    when every night carries the same relative uncertainty."""
    import importlib
    import sys
    sys.path.insert(0, str(REPO / "scripts"))
    monkeypatch.delenv("ALC_KALMAN_UNCERTAINTY_WEIGHTING", raising=False)
    rnc = importlib.import_module("run_network_calibration")
    rs = np.random.default_rng(3)
    rows = []
    for k in range(60):
        c = 1e12 * (1.0 + 0.05 * rs.standard_normal())
        rel = 0.04 if k % 3 else 0.40                    # every third night much less certain
        if k % 3 == 0:
            c *= 1.3                                     # ... and biased
        rows.append(dict(date=f"2025{1 + k // 28:02d}{1 + k % 28:02d}", method="rayleigh", flag=1,
                         cal_value=c, uncertainty=rel * c))
    plain = rnc._kalman_rows(rows)
    assert plain == rnc._kalman_rows(rows, weighted=False)
    weighted = rnc._kalman_rows(rows, weighted=True)
    assert len(weighted) == len(plain) and weighted != plain
    # the biased, uncertain nights pull the weighted estimate less
    last_plain = float(plain[-1]["kalman"])
    last_w = float(weighted[-1]["kalman"])
    assert abs(last_w / 1e12 - 1) < abs(last_plain / 1e12 - 1)
    same = [dict(r, uncertainty=0.1 * float(r["cal_value"])) for r in rows]
    assert rnc._kalman_rows(same, weighted=True) == rnc._kalman_rows(same, weighted=False)
    monkeypatch.setenv("ALC_KALMAN_UNCERTAINTY_WEIGHTING", "1")
    assert rnc._kalman_rows(rows) == weighted


def test_molecular_climatology_bounds():
    vals = [U.molecular_bias_climatology(lat, doy, zc)
            for lat in (10, 40, 47, 55, 75, -45) for doy in (15, 105, 196, 288)
            for zc in (2.5, 4.5, 6.5)]
    assert np.all(np.isfinite(vals)) and np.max(np.abs(vals)) < 0.06
    # summer: the US-1976 model is too dense -> C_L biased low (negative ln-bias)
    assert U.molecular_bias_climatology(47.0, 196, 4.5) < 0 < U.molecular_bias_climatology(47.0, 30, 4.5)
    # southern hemisphere = mirrored season
    assert U.molecular_bias_climatology(-47.0, 196, 4.5) == pytest.approx(
        U.molecular_bias_climatology(47.0, 196 + 365.25 / 2, 4.5))
