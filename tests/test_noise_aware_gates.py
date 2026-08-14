# -*- coding: utf-8 -*-
"""Regression tests for the noise-aware molecular-window gates (eprof_v2.2).

The properties asserted here are the ones the availability study depends on, and each would be
silently lost by a plausible refactor:

  * eprof_v2.2 is a strict SUPERSET of eprof_v2 -- if the strict gates admit any window, v2.2 must
    return exactly v2's answer. This is what keeps existing calibrations bit-identical.
  * the noise-relative statistics stay NaN when no noise profile is supplied, so every other method
    is unaffected by their existence.
  * the de-biased scattering reference is only used by the fallback tier.
  * method params actually reach the selector (they used to be dropped).
  * a non-finite rel_error is rejected by gated methods instead of passing QC as if it were 0.
"""
from __future__ import annotations

import numpy as np
import pytest

from calibration.rayleigh.molecular_methods import (
    DEFAULT_PARAMS, METHODS, compute_window_grid, select_molecular_window)
from calibration.rayleigh.rayleigh_fit import find_optimal_molecular_window

HALF = tuple(range(250, 2000, 240))


def _profile(n=220, dz=30.0, seed=0, noise=0.02, aerosol=False):
    """A synthetic molecular profile: signal proportional to p_mol, plus optional aerosol below."""
    rng = np.random.default_rng(seed)
    rng_m = (np.arange(n) + 1) * dz
    p_mol = np.exp(-rng_m / 8000.0) / rng_m ** 0
    signal = 2.0e11 * p_mol
    if aerosol:                                   # a layer that breaks proportionality low down
        signal = signal * (1.0 + 3.0 * np.exp(-((rng_m - 2000.0) / 700.0) ** 2))
    sigma = noise * np.median(signal)
    signal = signal + rng.normal(0, sigma, n)
    stack = signal[None, :] + rng.normal(0, sigma * np.sqrt(8), (8, n))
    return signal, p_mol, rng_m, stack, np.full(n, sigma)


def test_v22_registered():
    assert "eprof_v2.2" in METHODS
    p = DEFAULT_PARAMS["eprof_v2.2"]
    assert p["use_noise"] is True
    # the aerosol threshold must NOT be relaxed relative to v2 -- only the statistic is de-biased
    assert p["max_scattering_ratio"] == DEFAULT_PARAMS["eprof_v2"]["max_scattering_ratio"]


def test_noise_stats_are_nan_without_sigma():
    """Without a noise profile the new arrays stay NaN, so other methods cannot be affected."""
    signal, p_mol, rng_m, stack, _ = _profile()
    g = compute_window_grid(signal, p_mol, rng_m, HALF, signal_stack=stack)
    assert np.all(np.isnan(g.chi2red))
    assert np.all(np.isnan(g.scattering_ratio_dbz))
    assert np.any(np.isfinite(g.scattering_ratio))          # the raw one is still computed


def test_v22_matches_v2_when_strict_gates_pass():
    """On a clean profile v2 finds a window, so v2.2 must return THE SAME one."""
    signal, p_mol, rng_m, stack, sigma = _profile(noise=0.01)
    a = select_molecular_window("eprof_v2", signal, p_mol, rng_m, HALF, signal_stack=stack)
    b = select_molecular_window("eprof_v2.2", signal, p_mol, rng_m, HALF, signal_stack=stack,
                                sigma_signal=sigma)
    if not a.ok:
        pytest.skip("synthetic profile does not pass the strict v2 gates")
    assert b.ok
    assert b.start_m == a.start_m and b.end_m == a.end_m
    assert b.cl == pytest.approx(a.cl, rel=1e-12)


def test_v22_never_loses_a_v2_window():
    """Across several noise levels: whenever v2 succeeds, v2.2 succeeds identically."""
    for seed in range(6):
        for noise in (0.005, 0.02, 0.05):
            signal, p_mol, rng_m, stack, sigma = _profile(seed=seed, noise=noise)
            a = select_molecular_window("eprof_v2", signal, p_mol, rng_m, HALF, signal_stack=stack)
            if not a.ok:
                continue
            b = select_molecular_window("eprof_v2.2", signal, p_mol, rng_m, HALF,
                                        signal_stack=stack, sigma_signal=sigma)
            assert b.ok, f"v2.2 lost a window v2 accepted (seed={seed}, noise={noise})"
            assert b.cl == pytest.approx(a.cl, rel=1e-12)


def test_chi2red_is_about_one_for_pure_noise():
    """A window whose only departure from Rayleigh is its own noise should have chi2red ~ 1."""
    signal, p_mol, rng_m, stack, sigma = _profile(noise=0.03, seed=3)
    g = compute_window_grid(signal, p_mol, rng_m, HALF, signal_stack=stack, sigma_signal=sigma)
    ok = np.isfinite(g.chi2red)
    assert ok.any()
    assert 0.2 < float(np.nanmedian(g.chi2red[ok])) < 5.0


def test_aerosol_still_rejected_by_chi2():
    """chi2red must stay large where a real layer breaks proportionality (not just noise)."""
    _, p_mol, rng_m, _, _ = _profile()
    sig_a, p_mol, rng_m, stack_a, sigma = _profile(noise=0.01, aerosol=True)
    g = compute_window_grid(sig_a, p_mol, rng_m, HALF, signal_stack=stack_a, sigma_signal=sigma)
    low = g.start_m < 2500                       # windows inside the aerosol layer
    vals = g.chi2red[low & np.isfinite(g.chi2red)]
    if vals.size:
        assert np.nanmax(vals) > 5.0, "aerosol curvature should be far above the noise level"


def test_method_params_reach_the_selector():
    """find_optimal_molecular_window used to DROP extra params for the pluggable methods."""
    signal, p_mol, rng_m, stack, _ = _profile(noise=0.01)
    wide = find_optimal_molecular_window(signal, p_mol, rng_m, HALF, method="eprof_v2",
                                         signal_stack=stack)
    blocked = find_optimal_molecular_window(signal, p_mol, rng_m, HALF, method="eprof_v2",
                                            signal_stack=stack,
                                            method_params=dict(min_window_start_m=1e9))
    assert not np.isfinite(blocked.slope) or np.isinf(blocked.relative_error), \
        "an impossible min_window_start_m must reject every window"
    assert wide is not None


def test_nonfinite_rel_error_rejected_for_gated_methods():
    """A window whose rel_error cannot be evaluated must not pass QC as if it were perfect."""
    class _MW:
        ok = True
        rel_error = np.nan
        slope = intercept = r2 = std_err = p_value = 1.0
        center_m = half_m = start_m = end_m = 1000.0
        grid = None
    from calibration.rayleigh.rayleigh_fit import _result_from_method_window
    assert _result_from_method_window(_MW(), gated=True).relative_error == np.inf
    assert _result_from_method_window(_MW(), gated=False).relative_error == 0.0
