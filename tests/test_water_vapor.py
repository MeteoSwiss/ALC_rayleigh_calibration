"""
Validation tests for the water-vapor correction (calibration.water_vapor_correction.water_vapor).

Cross-validates the Python implementation against:
  1. the MATLAB reference code (wv_t2eff.m, compute_wv_transmission.m,
     get_water_vapor_number_concentration_from_RH.m) via hardcoded reference values
     produced by running the MATLAB on identical inputs;
  2. the independent ACTRIS-Cloudnet ``atmoslib`` thermodynamics library
     (vapor_pressure, saturation_vapor_pressure, absolute_humidity).

Run:  pytest tests/test_water_vapor.py -v
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from calibration.water_vapor_correction.water_vapor import (
    wv_t2eff_core,
    two_way_wv_transmission,
    cams_water_vapor_profile,
    load_abs_cross_section,
    in_water_vapor_band,
    laser_spectrum_for,
    KB,
    EPS,
)

# Optional external reference (ACTRIS-Cloudnet)
try:
    import atmoslib
    HAVE_ATMOSLIB = True
except Exception:
    HAVE_ATMOSLIB = False

try:
    import xarray as xr  # noqa: F401
    HAVE_XARRAY = True
except Exception:
    HAVE_XARRAY = False

from calibration.water_vapor_correction.water_vapor import G0

N_A = 6.02214076e23      # Avogadro [1/mol]
M_H2O = 0.01801528       # molar mass of water [kg/mol]

# Local data files (integration tests skip if absent)
LUT = Path(r"C:\Users\hervo\OneDrive\Documents\MATLAB\MDA\monitoring_alc_monthly\abs_cross_647_full_levels_1000.nc")
CAMS_FEB = Path(r"D:\CAMS\CAMS_Beta_202602.nc")


# --------------------------------------------------------------------------- #
# 1. Transmission core vs MATLAB wv_t2eff.m (synthetic inputs, hardcoded refs)
# --------------------------------------------------------------------------- #
def _synthetic_grids():
    wl = np.arange(905.0, 915.0 + 1e-9, 0.05)   # MATLAB 905:0.05:915
    rng = np.arange(0.0, 3000.0 + 1e-9, 100.0)  # MATLAB 0:100:3000
    idx = [0, 10, 20, 30]                        # 0, 1000, 2000, 3000 m
    return wl, rng, idx


def test_core_vs_matlab_constant():
    """Constant abscs (1e-24 cm^2) + constant n_wv (1e23 m^-3), lambda0=910, FWHM=2."""
    wl, rng, idx = _synthetic_grids()
    abscs = np.full((wl.size, rng.size), 1e-24)
    nw = np.full(rng.size, 1e23)
    t2 = wv_t2eff_core(910.0, 2.0, wl, abscs, rng, nw)
    matlab = np.array([1.000000, 0.980257, 0.960905, 0.941936])   # from wv_t2eff.m
    assert np.allclose(t2[idx], matlab, atol=2e-4), f"{t2[idx]} vs {matlab}"


def test_core_vs_matlab_varying():
    """Wavelength-varying abscs and decaying n_wv profile."""
    wl, rng, idx = _synthetic_grids()
    abscs = (1e-24 * (1 + 0.5 * np.sin(wl - 910.0)))[:, None] * np.ones((1, rng.size))
    nw = 1e23 * np.exp(-rng / 2000.0)
    t2 = wv_t2eff_core(910.0, 2.0, wl, abscs, rng, nw)
    matlab = np.array([1.000000, 0.984439, 0.975130, 0.969532])   # from wv_t2eff.m
    assert np.allclose(t2[idx], matlab, atol=2e-4), f"{t2[idx]} vs {matlab}"


def test_core_bounds_and_monotonic():
    """T^2 must stay in (0, 1] and decrease with range for a positive WV profile."""
    wl, rng, _ = _synthetic_grids()
    abscs = np.full((wl.size, rng.size), 1e-24)
    nw = np.full(rng.size, 1e23)
    t2 = wv_t2eff_core(910.0, 2.0, wl, abscs, rng, nw)
    assert t2[0] == pytest.approx(1.0, abs=1e-6)
    assert np.all(t2 > 0) and np.all(t2 <= 1.0 + 1e-9)
    assert np.all(np.diff(t2) <= 1e-9)        # non-increasing


def test_no_absorption_gives_unity():
    wl, rng, _ = _synthetic_grids()
    abscs = np.zeros((wl.size, rng.size))
    nw = np.full(rng.size, 1e23)
    t2 = wv_t2eff_core(910.0, 2.0, wl, abscs, rng, nw)
    assert np.allclose(t2, 1.0, atol=1e-9)


# --------------------------------------------------------------------------- #
# 2. Water-vapor number density: MATLAB convention + atmoslib (Cloudnet)
# --------------------------------------------------------------------------- #
def _n_wv_from_q(p_pa, q, t_k):
    """Replicates water_vapor.cams_water_vapor_profile's n_wv computation."""
    pw = q * p_pa / (EPS + (1.0 - EPS) * q)
    return pw / (KB * t_k)


def test_n_wv_constant_matches_matlab():
    """Python n_wv = Pw/(kB*T); MATLAB get_water_vapor...RH uses nw = 7.25e22*Pw[Pa]/T.
    The two number-density conventions must agree (1/kB = 7.243e22)."""
    assert (1.0 / KB) == pytest.approx(7.25e22, rel=2e-3)


@pytest.mark.skipif(not HAVE_ATMOSLIB, reason="atmoslib not installed")
def test_vapor_pressure_matches_atmoslib():
    """My Pw(q,P) must match ACTRIS-Cloudnet atmoslib.vapor_pressure.
    rel=1e-3: same formula, but atmoslib uses a slightly different Mw/Md constant
    than our EPS=0.621981 (agreement is ~2 ppm — well below any physical effect)."""
    for p_pa, q in [(85000.0, 0.005), (95000.0, 0.001), (70000.0, 0.012)]:
        pw_mine = q * p_pa / (EPS + (1.0 - EPS) * q)
        pw_ref = float(atmoslib.vapor_pressure(p_pa, q))
        assert pw_mine == pytest.approx(pw_ref, rel=1e-3)


@pytest.mark.skipif(not HAVE_ATMOSLIB, reason="atmoslib not installed")
def test_n_wv_matches_atmoslib():
    """My n_wv must match atmoslib absolute_humidity converted to number density."""
    for p_pa, q, t_k in [(85000.0, 0.005, 280.0), (70000.0, 0.012, 290.0)]:
        n_mine = _n_wv_from_q(p_pa, q, t_k)
        pw = float(atmoslib.vapor_pressure(p_pa, q))
        rho_wv = float(atmoslib.absolute_humidity(t_k, pw))   # kg/m^3
        n_ref = rho_wv * N_A / M_H2O
        assert n_mine == pytest.approx(n_ref, rel=2e-3)


@pytest.mark.skipif(not HAVE_ATMOSLIB, reason="atmoslib not installed")
def test_n_wv_from_RH_matches_matlab():
    """RH -> n_wv: atmoslib saturation pressure + ideal gas vs MATLAB
    get_water_vapor_number_concentration_from_RH (IAPWS). Different saturation
    formulas, so allow ~1%."""
    cases = [(300.0, 0.50, 4.273534e23), (273.15, 0.80, 1.297834e23)]   # MATLAB refs
    for t_k, rh, n_matlab in cases:
        pws = float(atmoslib.saturation_vapor_pressure(t_k))
        pw = pws * rh
        n_mine = pw / (KB * t_k)
        assert n_mine == pytest.approx(n_matlab, rel=1e-2), f"T={t_k} RH={rh}: {n_mine:.4e} vs {n_matlab:.4e}"


# --------------------------------------------------------------------------- #
# 3. Helpers
# --------------------------------------------------------------------------- #
def test_band_and_laser_spectrum():
    assert in_water_vapor_band(910.74) and in_water_vapor_band(905.0)
    assert not in_water_vapor_band(1064.47)
    assert laser_spectrum_for("CL61", 910.0) == (910.74, 1.0)
    assert laser_spectrum_for("CHM15k", 1064.0)[0] == pytest.approx(1064.47)


# --------------------------------------------------------------------------- #
# 4. Integration vs MATLAB compute_wv_transmission.m (needs CAMS + LUT)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not (LUT.exists() and CAMS_FEB.exists()),
                    reason="CAMS/LUT data files not available")
def test_full_profile_vs_matlab_payerne():
    """Full Python chain (CAMS q -> n_wv -> T^2_wv) vs the MATLAB
    compute_wv_transmission reference for Payerne CL61, Feb 2026."""
    rng = np.arange(0.0, 6001.0, 30.0)
    salt = 491.0
    alt = rng + salt
    prof = cams_water_vapor_profile(CAMS_FEB, 46.81, 6.94,
                                    np.datetime64("2026-02-26T00:00"),
                                    np.datetime64("2026-02-26T04:00"))
    assert prof is not None
    h_wv, n_wv = prof
    # n_wv in the boundary layer must be physically sane (~1e22-1e24 m^-3)
    bl = (h_wv >= salt) & (h_wv <= salt + 2000)
    assert 1e22 < np.nanmedian(n_wv[bl]) < 1e24

    t2 = two_way_wv_transmission(alt, salt, h_wv, n_wv, LUT, 910.74, 1.0)
    def at(z):
        return t2[np.argmin(np.abs(rng - z))]
    # MATLAB compute_wv_transmission reference (CL61 910.74/1.0, Payerne):
    assert at(1000) == pytest.approx(0.9025, abs=0.012)
    assert at(3000) == pytest.approx(0.8201, abs=0.012)
    assert at(6000) == pytest.approx(0.8051, abs=0.012)
    assert np.all((t2 > 0) & (t2 <= 1.0 + 1e-9))


@pytest.mark.skipif(not LUT.exists(), reason="HITRAN LUT not available")
def test_lut_sanity():
    wl, h_m, abscs = load_abs_cross_section(LUT)
    assert wl.min() < 910 < wl.max()
    assert abscs.shape[0] == wl.size
    band = (wl > 909) & (wl < 912)
    assert np.nanmax(abscs[band]) > 0


# --------------------------------------------------------------------------- #
# 5. Regression test for the top-of-atmosphere / level-order bug (audit F3)
# --------------------------------------------------------------------------- #
def _synthetic_cams(tmp_path, level_ascending=True):
    """Write a minimal 137-level CAMS-like NetCDF (1x1 grid, 2 time steps).

    z and lnsp are surface fields broadcast over levels (as in real CAMS, only a
    surface value is meaningful). Profiles are smooth: cold/dry at the top
    (level 1), warm/moist at the surface (level 137)."""
    nlev = 137
    lev = np.arange(1, nlev + 1)
    frac = (lev - 1.0) / (nlev - 1.0)            # 0 at top -> 1 at surface
    t_prof = 220.0 + 70.0 * frac                 # 220 K (top) -> 290 K (surface)
    q_prof = 1e-6 + 8e-3 * frac ** 3             # ~0 (top) -> 8 g/kg (surface)
    sp = 95000.0
    z_surf = 491.0 * G0                          # surface geopotential [m^2/s^2]

    if not level_ascending:
        order = np.arange(nlev - 1, -1, -1)
        lev = lev[order]
        t_prof = t_prof[order]
        q_prof = q_prof[order]

    times = np.array(["2026-02-26T00:00", "2026-02-26T03:00"], dtype="datetime64[ns]")
    ntime = times.size
    t4 = np.broadcast_to(t_prof[None, :, None, None], (ntime, nlev, 1, 1)).copy()
    q4 = np.broadcast_to(q_prof[None, :, None, None], (ntime, nlev, 1, 1)).copy()
    z4 = np.full((ntime, nlev, 1, 1), z_surf)
    lnsp4 = np.full((ntime, nlev, 1, 1), np.log(sp))

    dims = ("time", "level", "latitude", "longitude")
    ds = xr.Dataset(
        {"t": (dims, t4), "q": (dims, q4), "z": (dims, z4), "lnsp": (dims, lnsp4)},
        coords={"time": times, "level": lev,
                "latitude": [46.81], "longitude": [6.94]},
    )
    suffix = "asc" if level_ascending else "desc"
    path = tmp_path / f"cams_{suffix}.nc"
    ds.to_netcdf(path)
    return path


@pytest.mark.skipif(not HAVE_XARRAY, reason="xarray not installed")
def test_geopotential_invariant_to_level_order(tmp_path):
    """The CAMS profile must be identical whether the level axis is stored
    ascending or descending. Guards against the shared latent bug where the
    top-of-atmosphere special case keyed off the loop index (i==0) instead of the
    physical top half-level (idx==0) -- a reordered axis would then divide by a
    nonzero Ph_lev at the wrong level (or by zero), corrupting the integration."""
    t0 = np.datetime64("2026-02-26T00:00")
    t1 = np.datetime64("2026-02-26T04:00")
    h_asc, n_asc = cams_water_vapor_profile(_synthetic_cams(tmp_path, True),
                                            46.81, 6.94, t0, t1)
    h_desc, n_desc = cams_water_vapor_profile(_synthetic_cams(tmp_path, False),
                                              46.81, 6.94, t0, t1)
    assert np.allclose(h_asc, h_desc, rtol=0, atol=1e-6)
    assert np.allclose(n_asc, n_desc, rtol=1e-10, atol=0)
    # sanity: finite, monotonic ascending height, physically plausible surface n_wv
    assert np.all(np.isfinite(h_asc)) and np.all(np.isfinite(n_asc))
    assert np.all(np.diff(h_asc) > 0)
    assert 1e23 < n_asc[0] < 5e23          # ~290 K, 8 g/kg surface layer
