# -*- coding: utf-8 -*-
"""Tests for the target-classification pre-fit cell mask and the repaired -11 night denominator.

The two screens read the SAME curtain with different code sets and for different purposes:

  * the pre-fit mask (aerosol included) removes contaminated cells before the window search, and
    must survive a classification grid finer than the fit grid -- one contaminated native profile
    contaminates the binned profile it lands in, so the reduction is ANY, not nearest;
  * the -11 veto (cloud/ice only) rejects the night, and its fraction must be computed over the
    NIGHT. Its denominator used to be the whole classification file: a cloud filling the fit window
    for a whole 7 h night scored ~15 % of 48 h and never reached the 30 % threshold, which is why
    the flag fired zero times across the network.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from calibration.rayleigh.calibration import (
    PREFIT_CONTAM_CODES, VETO_CONTAM_CODES, _classification_on_fit_grid)
from calibration.rayleigh.molecular_methods import select_molecular_window

HALF = tuple(range(250, 2000, 240))


def _curtain(day="2026-03-04", dt_s=30, n_range=1024, dz=15.0):
    """A 48 h all-clear curtain at 30 s x 15 m, ending at the end of ``day``."""
    t1 = np.datetime64(day, "s") + np.timedelta64(1, "D")
    n_t = int(2 * 86400 / dt_s)
    times = t1 - np.arange(n_t, 0, -1) * np.timedelta64(dt_s, "s")
    rng = (np.arange(n_range) + 1) * dz
    return times, rng, np.zeros((n_t, n_range), np.int8)


def _night_times(day="2026-03-04", hours=(22, 5), step_s=300):
    """Fit-profile times: 22:00 on d-1 to 05:00 on d, on the 300 s L2 grid."""
    t0 = np.datetime64(day, "s") - np.timedelta64(24 - hours[0], "h")
    n = int((24 - hours[0] + hours[1]) * 3600 / step_s)
    return [(t0 + np.timedelta64(i * step_s, "s")).astype("datetime64[s]").item()
            for i in range(n)]


def test_mask_maps_a_layer_to_the_right_cells():
    times, rng, codes = _curtain()
    layer = (rng >= 3000) & (rng < 4000)
    codes[:, layer] = 5                                    # aerosol, all 48 h
    night = _night_times()
    fit_rng = np.arange(1, 201) * 30.0
    m = _classification_on_fit_grid((times, rng, codes), night, fit_rng, PREFIT_CONTAM_CODES)
    assert m is not None and m.shape == (len(night), fit_rng.size)
    in_layer = (fit_rng >= 3000) & (fit_rng < 4000)
    assert m[:, in_layer].all()
    assert not m[:, ~in_layer].any()


def test_mask_is_any_not_nearest_when_classification_is_finer():
    """One contaminated 30 s profile inside a 300 s bin contaminates that bin."""
    times, rng, codes = _curtain()
    night = _night_times()
    t0 = np.datetime64(night[3], "s")
    hit = np.abs((times - t0).astype("int64")) < 20        # a single 30 s profile in bin 3
    assert hit.sum() == 1
    codes[np.ix_(hit, (rng >= 3000) & (rng < 3200))] = 1   # droplet
    fit_rng = np.arange(1, 201) * 30.0
    m = _classification_on_fit_grid((times, rng, codes), night, fit_rng, PREFIT_CONTAM_CODES)
    band = (fit_rng >= 3000) & (fit_rng < 3200)
    assert m[3, band].all(), "a contaminated sub-profile must contaminate its bin"
    assert not m[np.arange(len(night)) != 3][:, band].any()


def test_aerosol_is_in_the_prefit_mask_but_never_in_the_veto():
    times, rng, codes = _curtain()
    codes[:, (rng >= 3000) & (rng < 4000)] = 5             # aerosol only
    night = _night_times()
    fit_rng = np.arange(1, 201) * 30.0
    pre = _classification_on_fit_grid((times, rng, codes), night, fit_rng, PREFIT_CONTAM_CODES)
    veto = _classification_on_fit_grid((times, rng, codes), night, fit_rng, VETO_CONTAM_CODES)
    assert pre.any() and not veto.any()


def test_veto_fraction_is_over_the_night_not_the_file():
    """A cloud filling the fit window all night: ~100 % of the night, ~15 % of the 48 h file."""
    times, rng, codes = _curtain()
    night = _night_times()
    t_lo, t_hi = np.datetime64(night[0], "s"), np.datetime64(night[-1], "s")
    in_night = (times >= t_lo) & (times <= t_hi)
    band = (rng >= 3000) & (rng < 4000)
    codes[np.ix_(in_night, band)] = 3                      # ice, the whole night
    fit_rng = np.arange(1, 201) * 30.0
    veto = _classification_on_fit_grid((times, rng, codes), night, fit_rng, VETO_CONTAM_CODES)
    frac_night = veto[:, (fit_rng >= 3000) & (fit_rng < 4000)].mean()
    frac_file = np.isin(codes[:, band], VETO_CONTAM_CODES).mean()
    assert frac_night > 0.95, "over the night the window is fully contaminated"
    assert frac_file < 0.30, "over the file it never reaches the veto threshold"


def test_unclassified_cells_stay_unmasked():
    """A curtain that stops short in range must not mask the gates above it."""
    times, rng, codes = _curtain(n_range=200, dz=15.0)     # tops out at 3000 m
    codes[:] = 1
    night = _night_times()
    fit_rng = np.arange(1, 201) * 30.0                     # up to 6000 m
    m = _classification_on_fit_grid((times, rng, codes), night, fit_rng, PREFIT_CONTAM_CODES)
    assert m[:, fit_rng <= 3000].all()
    assert not m[:, fit_rng > 3100].any(), "unknown is not contaminated"


def test_extra_cell_mask_reaches_the_selector():
    """The mask must be unioned with the temporal screen, not ignored or substituted."""
    rng_m = (np.arange(220) + 1) * 30.0
    p_mol = np.exp(-rng_m / 8000.0)
    rs = np.random.default_rng(0)
    stack = 2.0e11 * p_mol[None, :] * (1.0 + 0.01 * rs.standard_normal((40, rng_m.size)))
    signal = stack.mean(axis=0)
    extra = np.zeros(stack.shape, bool)
    extra[:, (rng_m >= 3000) & (rng_m < 4000)] = True
    a = select_molecular_window("eprof_v2", signal, p_mol, rng_m, HALF, signal_stack=stack)
    b = select_molecular_window("eprof_v2", signal, p_mol, rng_m, HALF, signal_stack=stack,
                                extra_cell_mask=extra)
    assert b.n_clean_frac < a.n_clean_frac
    assert b.cell_flag[:, (rng_m >= 3000) & (rng_m < 4000)].all()


def test_mask_of_wrong_shape_is_ignored():
    rng_m = (np.arange(220) + 1) * 30.0
    p_mol = np.exp(-rng_m / 8000.0)
    stack = np.tile(2.0e11 * p_mol, (40, 1))
    a = select_molecular_window("eprof_v2", stack.mean(axis=0), p_mol, rng_m, HALF,
                                signal_stack=stack)
    b = select_molecular_window("eprof_v2", stack.mean(axis=0), p_mol, rng_m, HALF,
                                signal_stack=stack, extra_cell_mask=np.zeros((3, 3), bool))
    assert a.ok == b.ok and a.start_m == b.start_m
