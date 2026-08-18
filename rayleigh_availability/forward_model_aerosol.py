# -*- coding: utf-8 -*-
"""MECHANISM SCAN 2 -- aerosol INSIDE the molecular fit window, and what the gates do about it.

WHY
---
The observational study shows the retrieved Rayleigh constant is not independent of where the
molecular window sits (Payerne CHM15k -22.5 % between v2.0-kept and v2.2-recovered nights,
Aosta +32 %). Two families of mechanisms can do that:

  * CONTAMINATION -- aerosol INSIDE the window inflates signal/p_mol locally, so C_L jumps when
    the window overlaps the layer and returns to normal when it does not: a LOCALISED bump;
  * TRANSMISSION -- aerosol BELOW (or spread through) the window attenuates everything above it,
    which the Klett correction is supposed to undo: a MONOTONIC slope.

This module quantifies the first one on synthetic profiles built from a KNOWN ``C_true`` and
pushed through the SHIPPED retrieval (imported from :mod:`rayleigh_availability.forward_model`,
which itself calls the production functions -- nothing is reimplemented here), and it runs the
REAL v2.0 / v2.2 gates on the same profiles.

THE FORWARD EQUATION (identical to forward_model.py, repeated here for the record)
----------------------------------------------------------------------------------
    beta_tot(z) = beta_mol(z) + beta_aer(z)                [m^-1 sr^-1]
    ext_tot(z)  = beta_mol(z) * 8*pi/3 + beta_aer(z)*S_aer [m^-1]
    T2(z)       = exp(-2 * INT_0^z ext_tot dz')            [-]
    rcs(z)      = C_true * beta_tot(z) * T2(z)             (overlap = 1, no background)

with S_aer = 52 sr (the operational ``lidar_ratio_aerosol``), z in m AGL. The aerosol layer is a
smooth top-hat (Hann-tapered edges, ``forward_model.aerosol_layer``) whose PEAK backscatter is set
from a target backscatter ratio at the layer centre:

    R_peak = beta_tot / beta_mol  at  z = layer centre   ->   beta_peak = (R_peak-1)*beta_mol(centre)

NOISE. The night is built as a STACK of ``n_profiles`` binned profiles (300 s x 30 m), each with
1-sigma ``sigma_night * sqrt(n_profiles)`` on the range-normalised signal, so the stack mean has
exactly the measured night-mean noise of a clean Payerne CHM15k night (1.2e-4 signal units, see
``forward_model.PAYERNE_CHM15K``). The window search receives that stack -- as production does --
so the temporal-variability gate and the contaminated-cell screen are exercised for real.

THE SIGMA BUG. ``calibration.rayleigh.calibration._sigma_on_fit_grid`` divides the native per-gate
noise by sqrt(n_BINNED profiles) instead of sqrt(n_NATIVE profiles), i.e. it omits
sqrt(avg_time/dt_native) = sqrt(300/15) = 4.472. The v2.2 chi2red gate therefore sees a sigma
4.472x too large and a chi2red 20.0x too small. Every v2.2 case below is run BOTH ways
(``sigma_shipped = 4.472 * sigma_true`` and ``sigma_true``) so the question "does fixing the bug
change the science?" gets a number.

OUTPUTS
-------
    doc/reports/figs_altitude_audit/fwd_aer_ladders.png     (a) the C_L(window) ladder
    doc/reports/figs_altitude_audit/fwd_aer_gates.png       (b)+(c) what the gates accept
    doc/reports/figs_altitude_audit/fwd_aer_signature.png   (d) bump vs slope discriminator
    doc/reports/figs_altitude_audit/fwd_aer_results.json    every number in the report

Run with::

    python rayleigh_availability/forward_model_aerosol.py
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from calibration.rayleigh.molecular_methods import (          # noqa: E402
    DEFAULT_PARAMS,
    compute_window_grid,
    select_molecular_window,
)
from rayleigh_availability.forward_model import (              # noqa: E402
    PAYERNE_CHM15K,
    aerosol_layer,
    default_grid,
    default_options,
    exponential_aerosol,
    forward_signal,
    make_atmosphere,
    retrieve,
    scan_windows,
)

FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"

# --- experiment constants (all explicit, all SI unless stated) ------------------------------
C_TRUE = float(PAYERNE_CHM15K["C_true"])          # imposed lidar constant, instrument units
STATION_ALT_M = float(PAYERNE_CHM15K["station_alt_m"])
WAVELENGTH_NM = float(PAYERNE_CHM15K["wavelength_nm"])
S_AER = 52.0                                      # aerosol lidar ratio, sr (options.json LRaer)
N_PROFILES = int(PAYERNE_CHM15K["n_profiles_binned"])          # 135 binned profiles = one night
SIGMA_NIGHT = float(PAYERNE_CHM15K["sigma_night_clean"])       # 1.2e-4, night-MEAN signal noise
SIGMA_BUG_FACTOR = float(np.sqrt(300.0 / 15.0))                # 4.4721: the omitted binning factor

LADDER_HALF_M = 490.0                             # a real half_length_options_m entry
LADDER_CENTRES_M = np.arange(2250.0, 6501.0, 250.0)
LAYER_THICKNESS_M = 750.0                         # "thin lofted layer", 0.5-1 km as specified


# ---------------------------------------------------------------------------
# 1. Case construction
# ---------------------------------------------------------------------------
def build_atmosphere(z: Optional[np.ndarray] = None) -> tuple[np.ndarray, dict]:
    """Payerne CHM15k grid + US-Std molecular atmosphere (the production molecular routine)."""
    z = default_grid() if z is None else np.asarray(z, float)
    atm = make_atmosphere(z, station_alt_m=STATION_ALT_M, wavelength_nm=WAVELENGTH_NM)
    return z, atm


def beta_peak_for_ratio(z: np.ndarray, beta_mol: np.ndarray, centre_m: float,
                        r_peak: float) -> float:
    """Peak aerosol backscatter [m^-1 sr^-1] giving ``beta_tot/beta_mol = r_peak`` at ``centre_m``."""
    i = int(np.argmin(np.abs(z - float(centre_m))))
    return float((float(r_peak) - 1.0) * beta_mol[i])


def make_case(
    z: np.ndarray,
    atm: dict,
    *,
    kind: str = "layer",
    layer_centre_m: float = 3500.0,
    thickness_m: float = LAYER_THICKNESS_M,
    r_peak: float = 1.0,
    scale_height_m: float = 2000.0,
    haze_top_m: Optional[float] = None,
    seed: int = 0,
    noise: bool = True,
) -> dict:
    """Build one synthetic night: aerosol field, night-mean rcs and the per-profile stack.

    ``kind`` is ``"layer"`` (Hann-tapered top-hat of ``thickness_m`` centred on
    ``layer_centre_m``, peak set by ``r_peak``) or ``"haze"`` (exponential
    ``beta_surface*exp(-z/scale_height_m)`` optionally truncated at ``haze_top_m``; ``r_peak``
    then sets the ratio at ``layer_centre_m``, so the two families are directly comparable).

    Returns a dict with ``rcs_mean`` (night mean, the array the constant is computed from),
    ``signal_stack`` (n_profiles x n_range, what the window search consumes), ``beta_aer``,
    ``ext_aer``, ``aod`` (aerosol optical depth, one-way, whole column) and the case metadata.
    """
    beta_mol = atm["beta_mol"]
    if kind == "layer":
        base = float(layer_centre_m) - 0.5 * float(thickness_m)
        top = float(layer_centre_m) + 0.5 * float(thickness_m)
        peak = beta_peak_for_ratio(z, beta_mol, layer_centre_m, r_peak)
        beta_aer, ext_aer = aerosol_layer(z, base, top, peak, S_AER)
        geom = dict(base_m=base, top_m=top)
    elif kind == "haze":
        i = int(np.argmin(np.abs(z - float(layer_centre_m))))
        target = (float(r_peak) - 1.0) * beta_mol[i]
        beta_surface = float(target * np.exp(z[i] / float(scale_height_m)))
        beta_aer, ext_aer = exponential_aerosol(z, beta_surface, scale_height_m, S_AER,
                                                top_m=haze_top_m)
        geom = dict(beta_surface=beta_surface, scale_height_m=float(scale_height_m),
                    haze_top_m=haze_top_m)
    else:
        raise ValueError(f"unknown case kind {kind!r}")

    rcs_true = forward_signal(z, C_TRUE, beta_mol, beta_aer=beta_aer, ext_aer=ext_aer)
    signal_true = rcs_true / z ** 2

    if noise:
        rng = np.random.default_rng(int(seed))
        sigma_profile = SIGMA_NIGHT * np.sqrt(N_PROFILES)
        stack = signal_true[None, :] + rng.normal(0.0, sigma_profile, (N_PROFILES, z.size))
    else:
        stack = np.repeat(signal_true[None, :], N_PROFILES, axis=0)
    signal_mean = stack.mean(axis=0)
    rcs_mean = signal_mean * z ** 2

    aod = float(np.trapezoid(ext_aer, z)) if hasattr(np, "trapezoid") else float(np.trapz(ext_aer, z))
    return dict(kind=kind, layer_centre_m=float(layer_centre_m), thickness_m=float(thickness_m),
                r_peak=float(r_peak), seed=int(seed), noise=bool(noise),
                beta_aer=beta_aer, ext_aer=ext_aer, aod=aod,
                rcs_true=rcs_true, rcs_mean=rcs_mean, signal_stack=stack, **geom)


# ---------------------------------------------------------------------------
# 2. (a) the C_L(window) ladder
# ---------------------------------------------------------------------------
def ladder(z, case, atm, centres_m=LADDER_CENTRES_M, half_m=LADDER_HALF_M) -> np.ndarray:
    """Retrieved C_L for FORCED windows [centre-half, centre+half], via the shipped chain.

    Returns the relative deviation (C_L/C_true - 1) * 100 in %, one value per centre.
    """
    values = np.array([v for _, v in scan_windows(z, case["rcs_mean"], atm, centres_m, half_m)],
                      float)
    return (values / C_TRUE - 1.0) * 100.0


def ladder_metrics(centres_m, dev_pct) -> dict:
    """Split a ladder into its MONOTONIC part and its LOCALISED part -- the (d) discriminator.

    * ``slope_pct_per_km``  least-squares slope of dev vs centre over the whole ladder;
    * ``linear_span_pct``   |slope| * (top - bottom), i.e. how much of the excursion a straight
                            line explains;
    * ``bump_pct``          max |dev - linear fit|, the part a straight line canNOT explain;
    * ``localisation``      bump / (bump + linear_span) -- 1 = pure localised bump,
                            0 = pure monotonic slope;
    * ``bump_centre_km``    the window centre where |residual| peaks (where the layer is);
    * ``fwhm_km``           full width at half maximum of |residual| around that peak.
    """
    c = np.asarray(centres_m, float) / 1000.0
    d = np.asarray(dev_pct, float)
    ok = np.isfinite(d)
    if ok.sum() < 4:
        return dict(slope_pct_per_km=np.nan, linear_span_pct=np.nan, bump_pct=np.nan,
                    localisation=np.nan, bump_centre_km=np.nan, fwhm_km=np.nan)
    coeff = np.polyfit(c[ok], d[ok], 1)
    resid = d - np.polyval(coeff, c)
    resid[~ok] = np.nan
    span = float(c[ok].max() - c[ok].min())
    linear_span = abs(float(coeff[0])) * span
    k = int(np.nanargmax(np.abs(resid)))
    bump = float(abs(resid[k]))
    half = 0.5 * bump
    above = np.abs(resid) >= half
    fwhm = float(np.nansum(above) * np.median(np.diff(c))) if bump > 0 else np.nan
    denom = bump + linear_span
    # Below 0.05 % total excursion the ladder is flat to within the estimator's own closure
    # (0.005 %) and the noise floor: the split into "bump" and "slope" is then meaningless.
    if denom < 0.05:
        return dict(slope_pct_per_km=float(coeff[0]), linear_span_pct=float(linear_span),
                    bump_pct=bump, localisation=np.nan, bump_centre_km=np.nan, fwhm_km=np.nan)
    return dict(slope_pct_per_km=float(coeff[0]), linear_span_pct=float(linear_span),
                bump_pct=bump, localisation=float(bump / denom) if denom > 0 else np.nan,
                bump_centre_km=float(c[k]), fwhm_km=fwhm)


# ---------------------------------------------------------------------------
# 3. (b)+(c) the REAL gates
# ---------------------------------------------------------------------------
def run_gates(z, case, atm, method: str, sigma_signal=None) -> dict:
    """Run the SHIPPED window search (``select_molecular_window``) and price the window it picks.

    ``method`` is ``"eprof_v2"`` or ``"eprof_v2.2"``; both dispatch to ``_select_optimal`` with
    their registered ``DEFAULT_PARAMS`` -- exactly what production does (options.json's
    min_window_* are ignored by every method except v1.2, see rayleigh_fit.py).

    The C_L of the accepted window is then computed by forcing that window through
    ``forward_model.retrieve`` (Klett + ``calculate_lidar_constant``), on the FULL night-mean rcs
    -- production computes the constant from ``rcs_mean``, not from the cell-cleaned mean.
    """
    opts = default_options(method)
    stack = case["signal_stack"]
    mw = select_molecular_window(
        method, stack.mean(axis=0), atm["p_mol"], z, opts.half_length_options_m,
        range_start_m=opts.range_start_m, range_end_m=opts.range_end_m,
        increment_bins=opts.fit_range_increment_bins,
        signal_stack=stack, sigma_signal=sigma_signal,
    )
    out = dict(method=method, ok=bool(mw.ok), message=str(mw.message),
               n_eligible=int(np.nansum(mw.eligible)) if mw.eligible is not None else 0,
               start_m=float(mw.start_m), end_m=float(mw.end_m), centre_m=float(mw.center_m),
               r2=float(mw.r2), residual_pct=float(mw.residual_pct),
               scattering_ratio=float(mw.scattering_ratio), rel_error=float(mw.rel_error),
               temporal_cv=float(mw.temporal_cv), C_L=np.nan, bias_pct=np.nan)
    if not mw.ok:
        return out
    res = retrieve(z, case["rcs_mean"], atm, options=opts,
                   window=(float(mw.start_m), float(mw.end_m)))
    if res["ok"]:
        out["C_L"] = float(res["C_L"])
        out["bias_pct"] = float(res["C_L"] / C_TRUE - 1.0) * 100.0
    return out


def run_tier2_probe(loads=(1.02, 1.05, 1.10, 1.20, 1.30, 1.50, 1.80, 2.20),
                    kind="haze", layer_centre_m=3500.0, seed=3, verbose=True) -> dict:
    """Isolate the sigma bug: FORCE the v2.2 noise tier and compare shipped vs corrected sigma.

    In these synthetic cases the strict tier almost never empties, so the noise tier -- the only
    part of v2.2 that reads ``sigma_signal`` -- is rarely exercised, and the sigma bug looks
    harmless. To price it we make the strict tier empty ON PURPOSE by passing an impossible
    ``min_r2`` (2.0) to the SHIPPED selector: everything else (the chi2red gate at 2.5, the
    de-biased scattering ratio, the excess-CV gates, the composite score) is untouched, and the
    fallback branch runs exactly as it would on a real flag -2 night.

    Reported per load: how many windows the noise tier accepts and which one it picks, with the
    sigma the pipeline currently computes (4.472x too large -> chi2red 20x too small) and with the
    correct one.
    """
    z, atm = build_atmosphere()
    opts = default_options("eprof_v2.2")
    sigma_true = np.full(z.size, SIGMA_NIGHT)
    rows = []
    for r in loads:
        kwargs = (dict(kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0, r_peak=r)
                  if kind == "haze" else
                  dict(kind="layer", layer_centre_m=layer_centre_m,
                       thickness_m=LAYER_THICKNESS_M, r_peak=r))
        case = make_case(z, atm, seed=seed, **kwargs)
        entry = dict(r_peak=float(r), kind=kind, aod=case["aod"])
        for label, sigma in (("shipped", sigma_true * SIGMA_BUG_FACTOR), ("fixed", sigma_true)):
            mw = select_molecular_window(
                "eprof_v2.2", case["signal_stack"].mean(axis=0), atm["p_mol"], z,
                opts.half_length_options_m, range_start_m=opts.range_start_m,
                range_end_m=opts.range_end_m, increment_bins=opts.fit_range_increment_bins,
                signal_stack=case["signal_stack"], sigma_signal=sigma,
                min_r2=2.0,                       # impossible -> the strict tier is empty
            )
            d = dict(ok=bool(mw.ok), n_eligible=int(np.nansum(mw.eligible)) if mw.eligible
                     is not None else 0, start_m=float(mw.start_m), end_m=float(mw.end_m),
                     bias_pct=np.nan)
            if mw.ok:
                res = retrieve(z, case["rcs_mean"], atm, options=opts,
                               window=(float(mw.start_m), float(mw.end_m)))
                if res["ok"]:
                    d["bias_pct"] = float(res["C_L"] / C_TRUE - 1.0) * 100.0
            entry[label] = d
        rows.append(entry)
        if verbose:
            s, f = entry["shipped"], entry["fixed"]
            print(f"  {kind} R={r:4.2f} (AOD {entry['aod']:.4f}): shipped sigma -> "
                  f"{s['n_eligible']:3d} windows, {s['start_m']/1000:4.2f}-{s['end_m']/1000:4.2f} km,"
                  f" bias {s['bias_pct']:+7.2f} %   |   fixed sigma -> {f['n_eligible']:3d} windows,"
                  f" {f['start_m']/1000:4.2f}-{f['end_m']/1000:4.2f} km, bias {f['bias_pct']:+7.2f} %")
    return dict(kind=kind, seed=int(seed), rows=rows)


def chi2red_stats(z, case, atm, sigma_signal) -> dict:
    """Min / median reduced chi-square over the window grid, for the sigma-bug discussion."""
    opts = default_options("eprof_v2.2")
    g = compute_window_grid(case["signal_stack"].mean(axis=0), atm["p_mol"], z,
                            opts.half_length_options_m, range_start_m=opts.range_start_m,
                            range_end_m=opts.range_end_m,
                            increment_bins=opts.fit_range_increment_bins,
                            signal_stack=case["signal_stack"], sigma_signal=sigma_signal)
    c = g.chi2red[np.isfinite(g.chi2red)]
    if c.size == 0:
        return dict(min=np.nan, median=np.nan, frac_below_2p5=np.nan)
    return dict(min=float(np.min(c)), median=float(np.median(c)),
                frac_below_2p5=float(np.mean(c <= 2.5)))


def _optical_depth(ext: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Cumulative one-way optical depth from the FIRST GATE, the convention of
    ``rayleigh_fit.calculate_lidar_constant`` (trapezoid, no rectangle below z[0])."""
    out = np.zeros_like(z)
    out[1:] = np.cumsum(0.5 * np.diff(z) * (ext[:-1] + ext[1:]))
    return out


def decompose_bias(z, case, atm, window) -> dict:
    """Split the C_L error into its BACKSCATTER and TRANSMISSION halves, and test the law.

    ``calculate_lidar_constant`` forms ``cl(z) = rcs(z) / beta_tot_ret(z) * exp(2*OD_ret(z))``.
    With ``rcs = C_true * beta_tot_true * exp(-2*OD_true)`` this is exactly

        cl(z) / C_true = [beta_tot_true(z) / beta_tot_ret(z)] * exp(2*(OD_ret(z) - OD_true(z)))
                         \\________ backscatter term ________/  \\____ transmission term ____/

    so the two error sources are separable with no modelling assumption. The Klett reference is
    ``nanmean(beta_att/beta_mol)`` over the window and ``beta_att`` is normalised by the fitted
    slope, so the aerosol INSIDE the window is absorbed into the constant by construction: the
    prediction is the window-MEAN backscatter ratio,

        C_L / C_true - 1  ~  <beta_tot/beta_mol>_window - 1                                [law]

    Returned: the measured bias, the two terms, and the law's prediction (all in %).
    """
    from calibration.rayleigh.atmosphere import MOLECULAR_LIDAR_RATIO
    res = retrieve(z, case["rcs_mean"], atm, window=window)
    if not res["ok"]:
        return dict(window_m=[float(window[0]), float(window[1])], ok=False)
    beta_mol = atm["beta_mol"]
    beta_tot_true = beta_mol + case["beta_aer"]
    ext_true = beta_mol * MOLECULAR_LIDAR_RATIO + case["ext_aer"]
    m = (z >= window[0]) & (z <= window[1])
    backscatter = float(np.median(beta_tot_true[m] / res["beta_tot"][m]))
    transmission = float(np.median(np.exp(2.0 * (_optical_depth(res["ext_tot"], z)[m]
                                                 - _optical_depth(ext_true, z)[m]))))
    ratio_mean = float(np.mean(beta_tot_true[m] / beta_mol[m]))
    aer_recovered = float(np.median(res["beta_tot"][m] - beta_mol[m])
                          / np.median(case["beta_aer"][m])) if np.median(case["beta_aer"][m]) > 0 \
        else np.nan
    return dict(window_m=[float(window[0]), float(window[1])], ok=True,
                bias_pct=float(res["C_L"] / C_TRUE - 1.0) * 100.0,
                backscatter_term_pct=(backscatter - 1.0) * 100.0,
                transmission_term_pct=(transmission - 1.0) * 100.0,
                law_prediction_pct=(ratio_mean - 1.0) * 100.0,
                aerosol_recovered_fraction=aer_recovered)


def run_decomposition(verbose=True) -> dict:
    """Where the bias comes from, on representative (case, window) pairs."""
    z, atm = build_atmosphere()
    cases = [
        ("thin layer 3.5 km, R=1.50", dict(kind="layer", layer_centre_m=3500.0, r_peak=1.50,
                                           noise=False), (3010.0, 3990.0)),
        ("thin layer 3.5 km, R=1.50", dict(kind="layer", layer_centre_m=3500.0, r_peak=1.50,
                                           noise=False), (1740.0, 5570.0)),
        ("haze R(3 km)=1.10", dict(kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                                   r_peak=1.10, noise=False), (3010.0, 3990.0)),
        ("haze R(3 km)=1.50", dict(kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                                   r_peak=1.50, noise=False), (3010.0, 3990.0)),
        ("haze R(3 km)=1.50", dict(kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                                   r_peak=1.50, noise=False), (3420.0, 7250.0)),
        ("residual layer 0.1-1.7 km, R=3.0", dict(kind="layer", layer_centre_m=900.0,
                                                  thickness_m=1600.0, r_peak=3.0, noise=False),
         (3010.0, 3990.0)),
    ]
    out = []
    for label, kwargs, win in cases:
        case = make_case(z, atm, **kwargs)
        d = decompose_bias(z, case, atm, win)
        d["label"] = label
        out.append(d)
        if verbose:
            print(f"  {label:34s} window {win[0]/1000:4.2f}-{win[1]/1000:4.2f} km: bias "
                  f"{d['bias_pct']:+7.2f} % = backscatter {d['backscatter_term_pct']:+7.2f} % x "
                  f"transmission {d['transmission_term_pct']:+6.2f} %  |  law predicts "
                  f"{d['law_prediction_pct']:+7.2f} %  |  aerosol recovered in window: "
                  f"{d['aerosol_recovered_fraction']:+.3f}")
    return out


def _gate_case(args) -> dict:
    """One (case, seed) through v2, v2.2-with-shipped-sigma and v2.2-with-corrected-sigma."""
    family, case_kwargs, seed = args
    z, atm = build_atmosphere()
    case = make_case(z, atm, seed=seed, **case_kwargs)
    sigma_true = np.full(z.size, SIGMA_NIGHT)
    sigma_ship = sigma_true * SIGMA_BUG_FACTOR
    sig = case["rcs_true"] / z ** 2
    band = (z >= 2000.0) & (z <= 6000.0)
    ratio_slope = float(np.polyfit(z[band] / 1000.0,
                                   np.log(sig[band] / atm["p_mol"][band]), 1)[0] * 100.0)
    row = dict(family=str(family), seed=int(seed), aod=case["aod"],
               layer_centre_m=float(case_kwargs.get("layer_centre_m", np.nan)),
               r_peak=float(case_kwargs.get("r_peak", np.nan)),
               thickness_m=float(case_kwargs.get("thickness_m", LAYER_THICKNESS_M)),
               ratio_slope_pct_per_km=ratio_slope)
    row["v2"] = run_gates(z, case, atm, "eprof_v2")
    row["v22_shipped"] = run_gates(z, case, atm, "eprof_v2.2", sigma_signal=sigma_ship)
    row["v22_fixed"] = run_gates(z, case, atm, "eprof_v2.2", sigma_signal=sigma_true)
    return row


# ---------------------------------------------------------------------------
# 4. Studies
# ---------------------------------------------------------------------------
def run_closure(n_seeds: int = 20, verbose=True) -> dict:
    """Closure of THIS scan's chain: R_peak = 1 (no aerosol) through the identical code path.

    Two numbers: the forced-window ladder on a NOISELESS pure-molecular profile (must return
    C_true; this is the estimator closure) and the searched-window constant with the real gates on
    NOISY 135-profile stacks (the closure of the whole gate+retrieval chain used in (b)/(c)),
    averaged over ``n_seeds`` noise realisations so the bias is separable from the scatter.
    """
    z, atm = build_atmosphere()
    clean = make_case(z, atm, kind="layer", r_peak=1.0, noise=False)
    dev = ladder(z, clean, atm)
    m = ladder_metrics(LADDER_CENTRES_M, dev)

    sigma_true = np.full(z.size, SIGMA_NIGHT)
    variants = dict(v2=("eprof_v2", None),
                    v22_shipped=("eprof_v2.2", sigma_true * SIGMA_BUG_FACTOR),
                    v22_fixed=("eprof_v2.2", sigma_true))
    gated = {}
    for key, (mth, sg) in variants.items():
        bias, starts, ends, ok = [], [], [], []
        for seed in range(int(n_seeds)):
            noisy = make_case(z, atm, kind="layer", r_peak=1.0, seed=100 + seed, noise=True)
            g = run_gates(z, noisy, atm, mth, sigma_signal=sg)
            ok.append(g["ok"])
            if g["ok"]:
                bias.append(g["bias_pct"])
                starts.append(g["start_m"])
                ends.append(g["end_m"])
        arr = np.asarray(bias, float)
        gated[key] = dict(accept_rate=float(np.mean(ok)), n_seeds=int(n_seeds),
                          bias_mean_pct=float(np.mean(arr)) if arr.size else np.nan,
                          bias_sem_pct=float(np.std(arr) / np.sqrt(arr.size)) if arr.size else np.nan,
                          bias_std_pct=float(np.std(arr)) if arr.size else np.nan,
                          median_start_m=float(np.median(starts)) if starts else np.nan,
                          median_end_m=float(np.median(ends)) if ends else np.nan)

    out = dict(max_abs_dev_pct=float(np.nanmax(np.abs(dev))),
               ladder_slope_pct_per_km=m["slope_pct_per_km"],
               dev_pct=dev.tolist(), centres_m=LADDER_CENTRES_M.tolist(), gated=gated)
    if verbose:
        print("CLOSURE (R_peak = 1, no aerosol)")
        print(f"  forced-window ladder, noiseless : max |C_L/C_true - 1| = "
              f"{out['max_abs_dev_pct']:.4f} %, slope = {m['slope_pct_per_km']:+.5f} %/km")
        for k, v in gated.items():
            print(f"  searched window + gates, {n_seeds} noisy nights, {k:12s}: accept "
                  f"{v['accept_rate']*100:.0f} %, window {v['median_start_m']:.0f}-"
                  f"{v['median_end_m']:.0f} m, bias {v['bias_mean_pct']:+.3f} "
                  f"+- {v['bias_sem_pct']:.3f} % (night-to-night sd {v['bias_std_pct']:.3f} %)")
    return out


def run_ladder_study(verbose=True) -> dict:
    """(a) ladder shape vs backscatter ratio (layer fixed) and vs layer altitude (R fixed)."""
    z, atm = build_atmosphere()
    r_list = [1.05, 1.10, 1.20, 1.50, 2.00, 3.00]
    alt_list = [2000.0, 3000.0, 4000.0, 5000.0]

    by_ratio = []
    for r in r_list:
        case = make_case(z, atm, kind="layer", layer_centre_m=3500.0, r_peak=r, noise=False)
        dev = ladder(z, case, atm)
        by_ratio.append(dict(r_peak=r, aod=case["aod"], dev_pct=dev.tolist(),
                             **ladder_metrics(LADDER_CENTRES_M, dev)))
        if verbose:
            m = by_ratio[-1]
            print(f"  layer 3.5 km, R={r:4.2f} (AOD {case['aod']:.4f}): peak excursion "
                  f"{np.nanmax(dev):+7.2f} %, bump {m['bump_pct']:.2f} % at "
                  f"{m['bump_centre_km']:.2f} km, FWHM {m['fwhm_km']:.2f} km, "
                  f"slope {m['slope_pct_per_km']:+.2f} %/km")

    by_altitude = []
    for a in alt_list:
        case = make_case(z, atm, kind="layer", layer_centre_m=a, r_peak=1.50, noise=False)
        dev = ladder(z, case, atm)
        by_altitude.append(dict(layer_centre_m=a, aod=case["aod"], dev_pct=dev.tolist(),
                                **ladder_metrics(LADDER_CENTRES_M, dev)))

    # How the bump dilutes with window WIDTH (all are real half_length_options_m entries): a wide
    # window averages the layer with clean air, so the same layer moves C_L less -- which is why
    # the amplitude of the excursion is a property of the (layer, window) pair, not of the layer.
    by_width = []
    case = make_case(z, atm, kind="layer", layer_centre_m=3500.0, r_peak=1.50, noise=False)
    for half in (250.0, 490.0, 970.0, 1450.0, 1930.0):
        dev = ladder(z, case, atm, half_m=half)
        mm = ladder_metrics(LADDER_CENTRES_M, dev)
        by_width.append(dict(half_m=half, peak_pct=float(np.nanmax(dev)), **mm))
        if verbose:
            print(f"  window half-width {half:6.0f} m ({2*half:5.0f} m wide): peak excursion "
                  f"{np.nanmax(dev):+7.2f} %, FWHM {mm['fwhm_km']:.2f} km")

    # thin (500 m) vs thick (1500 m) layer at fixed R -- the thick layer carries more aerosol, so
    # it both fills the window better and attenuates more.
    by_thickness = []
    for thick in (500.0, 750.0, 1500.0):
        case = make_case(z, atm, kind="layer", layer_centre_m=3500.0, thickness_m=thick,
                         r_peak=1.50, noise=False)
        dev = ladder(z, case, atm)
        by_thickness.append(dict(thickness_m=thick, aod=case["aod"], peak_pct=float(np.nanmax(dev)),
                                 **ladder_metrics(LADDER_CENTRES_M, dev)))

    return dict(centres_m=LADDER_CENTRES_M.tolist(), half_m=LADDER_HALF_M,
                thickness_m=LAYER_THICKNESS_M, by_ratio=by_ratio, by_altitude=by_altitude,
                by_width=by_width, by_thickness=by_thickness)


def run_gate_study(layer_centres=(2000., 2500., 3000., 3500., 4000., 4500., 5000.),
                   ratios=(1.05, 1.10, 1.20, 1.35, 1.50, 2.00, 2.50, 3.00),
                   haze_ratios=(1.02, 1.05, 1.10, 1.20, 1.30, 1.50, 1.80, 2.20, 2.60, 3.00),
                   thick_ratios=(1.05, 1.10, 1.20, 1.35, 1.50, 2.00, 2.50, 3.00),
                   seeds=(0, 1, 2, 3, 4, 5, 6, 7), workers=6, verbose=True) -> dict:
    """(b)+(c) the real gates, with noise, over three contamination families.

    * ``thin``  -- a 750 m layer at one altitude, the rest of the band clean. A clean window
                   ALWAYS exists somewhere in 1.5-6 km, so this family tests whether the gates
                   FIND it.
    * ``thick`` -- a 2500 m layer centred at 3.5 km: it covers most of the search band, so an
                   uncontaminated window is much harder to find.
    * ``haze``  -- exponential haze (scale height 2.5 km) reaching THROUGH the whole band: there
                   is NO clean window at all. This is the regime the v2.2-recovered nights are in
                   (their signal/p_mol falls ~19 %/km between 2 and 6 km).
    """
    jobs = []
    for a in layer_centres:
        for r in ratios:
            for s in seeds:
                jobs.append(("thin", dict(kind="layer", layer_centre_m=a,
                                          thickness_m=LAYER_THICKNESS_M, r_peak=r), s))
    for r in thick_ratios:
        for s in seeds:
            jobs.append(("thick", dict(kind="layer", layer_centre_m=3500.0,
                                       thickness_m=2500.0, r_peak=r), s))
    for r in haze_ratios:
        for s in seeds:
            jobs.append(("haze", dict(kind="haze", layer_centre_m=3000.0,
                                      scale_height_m=2500.0, r_peak=r), s))
    if verbose:
        print(f"gate study: {len(jobs)} cases on {workers} workers")
    with ProcessPoolExecutor(max_workers=int(workers)) as pool:
        rows = list(pool.map(_gate_case, jobs, chunksize=2))
    return dict(layer_centres_m=list(layer_centres), ratios=list(ratios),
                haze_ratios=list(haze_ratios), thick_ratios=list(thick_ratios),
                seeds=list(seeds), thickness_m=LAYER_THICKNESS_M, rows=rows)


def run_signature_study(verbose=True) -> dict:
    """(d) localised in-window contamination vs monotonic below-window / through-window haze."""
    z, atm = build_atmosphere()
    scenarios = []

    # 1. thin lofted layer INSIDE the search band
    for r in (1.20, 1.50):
        scenarios.append((f"thin layer 3.1-3.9 km, R={r:.2f}",
                          dict(kind="layer", layer_centre_m=3500.0, thickness_m=750.0,
                               r_peak=r, noise=False)))
    # 2. thick residual layer entirely BELOW the search band (pure transmission)
    for r in (1.50, 3.00):
        scenarios.append((f"residual layer 0.1-1.7 km, R={r:.2f}",
                          dict(kind="layer", layer_centre_m=900.0, thickness_m=1600.0,
                               r_peak=r, noise=False)))
    # 3. exponential haze reaching THROUGH the band (transmission + contamination)
    for r in (1.10, 1.30):
        scenarios.append((f"haze H=2.5 km, R(3 km)={r:.2f}",
                          dict(kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                               r_peak=r, noise=False)))

    out = []
    for label, kwargs in scenarios:
        case = make_case(z, atm, **kwargs)
        dev = ladder(z, case, atm)
        m = ladder_metrics(LADDER_CENTRES_M, dev)
        # signal/p_mol slope between 2 and 6 km -- the observational discriminator (%/km)
        sig = case["rcs_true"] / z ** 2
        ratio = sig / atm["p_mol"]
        band = (z >= 2000.0) & (z <= 6000.0)
        slope_ratio = float(np.polyfit(z[band] / 1000.0, np.log(ratio[band]), 1)[0] * 100.0)
        out.append(dict(label=label, kind=case["kind"], aod=case["aod"], case_kwargs=kwargs,
                        dev_pct=dev.tolist(), ratio_slope_pct_per_km=slope_ratio, **m))
        if verbose:
            loc = m["localisation"]
            loc_s = "  n/a" if not np.isfinite(loc) else f"{loc:5.2f}"
            print(f"  {label:34s} AOD {case['aod']:.4f}  ladder slope {m['slope_pct_per_km']:+8.4f} "
                  f"%/km  bump {m['bump_pct']:8.4f} %  localisation {loc_s}  "
                  f"signal/p_mol slope {slope_ratio:+6.2f} %/km")
    return dict(centres_m=LADDER_CENTRES_M.tolist(), scenarios=out)


# ---------------------------------------------------------------------------
# 4b. The observational comparison, reproduced
# ---------------------------------------------------------------------------
# Measured medians of (recovered night) / (kept night) - 1, per stream, from the availability
# study (rayleigh_availability/): these are the targets this scan must be compared against.
OBSERVED_RECOVERED_VS_KEPT = {
    "Payerne CHM15k": -22.5, "Payerne CL61": -17.5, "Lindenberg CHM15k": +23.2,
    "Aosta CHM15k": +31.7, "Aosta CL61": +32.1, "SIRTA CHM15k": +7.4,
}
OBSERVED_RATIO_SLOPE_PCT_PER_KM = -19.0   # signal/p_mol, 2-6 km, on recovered nights
OBSERVED_RATIO_SLOPE_KEPT_PCT_PER_KM = -4.0


def run_recovery_proxy(loads=(1.05, 1.10, 1.20, 1.30, 1.50, 1.80, 2.20, 2.60, 3.00),
                       low_centre_m=3000.0, high_centre_m=4500.0, verbose=True) -> dict:
    """Reproduce the OBSERVED comparison: a night fitted 1.5 km higher than the usual window.

    The availability study compares v2.2-recovered nights (which fit 1-2 km higher, because no
    low window passes the strict gates) with v2.0-kept nights. Two different quantities are mixed
    in that comparison and this function separates them on the same synthetic haze:

      * ``delta_height_pct``  = C_L(window at ``high_centre_m``) / C_L(window at ``low_centre_m``)
        - 1, i.e. the SAME night fitted higher -> NEGATIVE, because the window-mean scattering
        ratio decreases upward in a haze;
      * ``dev_low_pct`` / ``dev_high_pct`` = the bias of each window against C_true, i.e. a hazy
        night compared with a CLEAN night at the same window height -> POSITIVE.

    The sign a station actually shows is the sum of the two, so a site whose recovered nights fit
    much higher lands negative (Payerne) and a site whose recovered nights are mostly just hazier
    lands positive (Aosta, Lindenberg).
    """
    z, atm = build_atmosphere()
    centres = list(LADDER_CENTRES_M)
    i_low, i_high = centres.index(float(low_centre_m)), centres.index(float(high_centre_m))
    rows = []
    for r in loads:
        case = make_case(z, atm, kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                         r_peak=r, noise=False)
        dev = ladder(z, case, atm)
        m = ladder_metrics(LADDER_CENTRES_M, dev)
        sig = (case["rcs_true"] / z ** 2) / atm["p_mol"]
        band = (z >= 2000.0) & (z <= 6000.0)
        ratio_slope = float(np.polyfit(z[band] / 1000.0, np.log(sig[band]), 1)[0] * 100.0)
        delta = float(((1 + dev[i_high] / 100.0) / (1 + dev[i_low] / 100.0) - 1.0) * 100.0)
        rows.append(dict(r_peak=float(r), aod=case["aod"], ratio_slope_pct_per_km=ratio_slope,
                         ladder_slope_pct_per_km=m["slope_pct_per_km"],
                         dev_low_pct=float(dev[i_low]), dev_high_pct=float(dev[i_high]),
                         delta_height_pct=delta, dev_pct=dev.tolist()))
        if verbose:
            print(f"  R(3 km)={r:4.2f}  AOD {case['aod']:.4f}  signal/p_mol slope "
                  f"{ratio_slope:+6.2f} %/km  |  same night fitted "
                  f"{low_centre_m/1000:.1f}->{high_centre_m/1000:.1f} km: {delta:+7.2f} %  |  "
                  f"vs a CLEAN night: {dev[i_low]:+7.2f} % (low) / {dev[i_high]:+7.2f} % (high)")
    return dict(low_centre_m=float(low_centre_m), high_centre_m=float(high_centre_m),
                centres_m=LADDER_CENTRES_M.tolist(), rows=rows,
                observed=OBSERVED_RECOVERED_VS_KEPT,
                observed_ratio_slope_pct_per_km=OBSERVED_RATIO_SLOPE_PCT_PER_KM)


# ---------------------------------------------------------------------------
# 5. Figures (altitude ALWAYS on Y, landscape, dpi >= 130)
# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_ladders(res_ladder, path: Path) -> None:
    plt = _mpl()
    z, atm = build_atmosphere()
    centres = np.asarray(res_ladder["centres_m"], float) / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6))

    ax = axes[0]
    for entry in res_ladder["by_ratio"]:
        case = make_case(z, atm, kind="layer", layer_centre_m=3500.0,
                         r_peak=entry["r_peak"], noise=False)
        ratio = (case["rcs_true"] / z ** 2) / atm["p_mol"] / C_TRUE
        ax.plot(ratio, z / 1000.0, lw=1.4, label=f"R={entry['r_peak']:.2f}")
    ax.axhspan(3.125, 3.875, color="0.85", zorder=0)
    ax.set_xlim(0.9, 3.2)
    ax.set_ylim(0, 8)
    ax.set_xlabel("(signal / $p_{mol}$) / $C_{true}$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("Synthetic layer, 750 m thick, centred 3.5 km\n(grey band = the layer)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    ax.axvline(0.0, color="C3", ls="--", lw=1.4, label="$C_{true}$")
    for entry in res_ladder["by_ratio"]:
        ax.plot(entry["dev_pct"], centres, "o-", ms=3.5, lw=1.4,
                label=f"R={entry['r_peak']:.2f} (AOD {entry['aod']:.3f})")
    ax.axhspan(3.125, 3.875, color="0.85", zorder=0)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("(a) $C_L$ ladder, layer at 3.5 km\nforced windows $\\pm$490 m", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[2]
    ax.axvline(0.0, color="C3", ls="--", lw=1.4)
    for entry in res_ladder["by_altitude"]:
        zc = entry["layer_centre_m"] / 1000.0
        ax.plot(entry["dev_pct"], centres, "o-", ms=3.5, lw=1.4,
                label=f"layer at {zc:.1f} km")
        ax.axhline(zc, color="0.7", lw=0.7, ls=":")
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("(a) the bump follows the layer\nR = 1.50, 750 m thick", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Aerosol INSIDE the molecular window: the $C_L$(window) ladder "
                 f"(Payerne CHM15k grid, $C_{{true}}$={C_TRUE:.2e}, $S_{{aer}}$=52 sr)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _grid_from_rows(res_gate, key, field, family="thin"):
    """(n_altitudes, n_ratios) median over seeds of ``rows[key][field]`` (NaN where none)."""
    alts = res_gate["layer_centres_m"]
    ratios = res_gate["ratios"]
    out = np.full((len(alts), len(ratios)), np.nan)
    for i, a in enumerate(alts):
        for j, r in enumerate(ratios):
            vals = [row[key][field] for row in res_gate["rows"]
                    if row["family"] == family and row["layer_centre_m"] == a
                    and row["r_peak"] == r]
            vals = ([float(v) for v in vals] if field == "ok"
                    else [v for v in vals if v is not None and np.isfinite(v)])
            if vals:
                out[i, j] = float(np.median(vals))
    return np.asarray(alts, float), np.asarray(ratios, float), out


def family_summary(res_gate, family: str) -> list:
    """Per-ratio summary of one family: acceptance, chosen window, bias, for the three variants."""
    ratios = sorted({row["r_peak"] for row in res_gate["rows"] if row["family"] == family})
    out = []
    for r in ratios:
        rows = [row for row in res_gate["rows"]
                if row["family"] == family and row["r_peak"] == r]
        entry = dict(r_peak=float(r), n=len(rows),
                     aod=float(np.median([row["aod"] for row in rows])),
                     ratio_slope_pct_per_km=float(np.median(
                         [row["ratio_slope_pct_per_km"] for row in rows])))
        for key in ("v2", "v22_shipped", "v22_fixed"):
            ok = [bool(row[key]["ok"]) for row in rows]
            acc = [row[key] for row in rows if row[key]["ok"]]
            # The window is reported as the MODE over seeds, never as (median start, median end):
            # when the seeds split between two very different windows the two medians describe a
            # window no seed ever chose, and its bias is not the bias of any real retrieval.
            windows = [(round(a["start_m"]), round(a["end_m"])) for a in acc]
            mode_win, n_mode = (max(set(windows), key=windows.count),
                                windows.count(max(set(windows), key=windows.count))) \
                if windows else ((np.nan, np.nan), 0)
            biases = [a["bias_pct"] for a in acc]
            entry[key] = dict(
                accept=float(np.mean(ok)),
                n_eligible=float(np.median([row[key]["n_eligible"] for row in rows])),
                start_m=float(mode_win[0]), end_m=float(mode_win[1]),
                n_mode=int(n_mode), n_distinct_windows=int(len(set(windows))),
                bias_pct=float(np.median(biases)) if biases else np.nan,
                bias_min_pct=float(np.min(biases)) if biases else np.nan,
                bias_max_pct=float(np.max(biases)) if biases else np.nan,
                bias_spread_pct=float(np.std(biases)) if len(biases) > 1 else np.nan)
        # Did the v2.2 noise tier actually fire? It can only fire on a night v2 rejected.
        entry["tier2_fired_shipped"] = float(np.mean(
            [(not row["v2"]["ok"]) and row["v22_shipped"]["ok"] for row in rows]))
        entry["tier2_fired_fixed"] = float(np.mean(
            [(not row["v2"]["ok"]) and row["v22_fixed"]["ok"] for row in rows]))
        out.append(entry)
    return out


def fig_gates(res_gate, path: Path) -> None:
    plt = _mpl()
    keys = [("v2", "(b) E-PROF v2.0 gates"),
            ("v22_shipped", "(c) v2.2, sigma AS SHIPPED (4.47x too large)"),
            ("v22_fixed", "(c) v2.2, sigma CORRECTED")]
    fig, axes = plt.subplots(2, 3, figsize=(17.0, 9.0), sharey=True)
    alts, ratios, _ = _grid_from_rows(res_gate, "v2", "ok")
    y = alts / 1000.0
    xt = np.arange(len(ratios))

    for col, (key, title) in enumerate(keys):
        _, _, acc = _grid_from_rows(res_gate, key, "ok")
        _, _, bias = _grid_from_rows(res_gate, key, "bias_pct")
        _, _, wcen = _grid_from_rows(res_gate, key, "centre_m")

        _, _, nel = _grid_from_rows(res_gate, key, "n_eligible")
        ax = axes[0, col]
        im = ax.pcolormesh(xt, y, nel, cmap="RdYlGn", vmin=0, vmax=110, shading="nearest")
        for i in range(len(y)):
            for j in range(len(ratios)):
                if np.isfinite(nel[i, j]):
                    ax.text(xt[j], y[i], f"{nel[i, j]:.0f}", ha="center", va="center",
                            fontsize=7.5, color="0.15")
                    ax.text(xt[j], y[i] - 0.19, f"{acc[i, j]*100:.0f} %", ha="center", va="center",
                            fontsize=6, color="0.35")
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{r:.2f}" for r in ratios], fontsize=8)
        ax.set_xlabel("backscatter ratio $R$ at layer centre")
        if col == 0:
            ax.set_ylabel("layer centre altitude  [km AGL]")
        ax.set_title(title + "\neligible windows out of 136 (small text = acceptance)",
                     fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)

        ax = axes[1, col]
        vmax = 40.0
        im = ax.pcolormesh(xt, y, bias, cmap="coolwarm", vmin=-vmax, vmax=vmax, shading="nearest")
        for i in range(len(y)):
            for j in range(len(ratios)):
                if np.isfinite(bias[i, j]):
                    ax.text(xt[j], y[i], f"{bias[i, j]:+.0f}", ha="center", va="center",
                            fontsize=7.5, color="0.1")
                    if np.isfinite(wcen[i, j]):
                        ax.text(xt[j], y[i] - 0.19, f"{wcen[i, j]/1000:.1f}km", ha="center",
                                va="center", fontsize=6, color="0.35")
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{r:.2f}" for r in ratios], fontsize=8)
        ax.set_xlabel("backscatter ratio $R$ at layer centre")
        if col == 0:
            ax.set_ylabel("layer centre altitude  [km AGL]")
        ax.set_title("bias of the ACCEPTED $C_L$ [%]\n(small text = centre of the chosen window)",
                     fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)

    n_seeds = len(res_gate.get("seeds", []))
    fig.suptitle("What the shipped gates actually do with a THIN aerosol layer in the fit band "
                 f"(750 m thick, 135-profile nights, $\\sigma_{{night}}$=1.2e-4, {n_seeds} noise "
                 "seeds per cell)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_saturation(res_gate, path: Path) -> None:
    """Where the strict tier finally breaks: contamination the gates cannot dodge.

    Panels 1-2 keep altitude on Y (the atmospheres, and the window the gates end up choosing);
    panel 3 is the acceptance/bias curve versus aerosol load, which has no altitude axis.
    """
    plt = _mpl()
    z, atm = build_atmosphere()
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6))

    haze = family_summary(res_gate, "haze")
    thick = family_summary(res_gate, "thick")

    ax = axes[0]
    for e in haze:
        case = make_case(z, atm, kind="haze", layer_centre_m=3000.0, scale_height_m=2500.0,
                         r_peak=e["r_peak"], noise=False)
        ratio = (case["rcs_true"] / z ** 2) / atm["p_mol"] / C_TRUE
        ax.plot(ratio, z / 1000.0, lw=1.2,
                label=f"R(3 km)={e['r_peak']:.2f}, {e['ratio_slope_pct_per_km']:+.0f} %/km")
    ax.axhspan(1.5, 6.0, color="C2", alpha=0.06)
    ax.set_xlim(0.8, 3.0)
    ax.set_ylim(0, 8)
    ax.set_xlabel("(signal / $p_{mol}$) / $C_{true}$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("haze family: no clean window exists\n(green = the 1.5-6 km search band)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")

    ax = axes[1]
    for label, fam, colour in (("haze", haze, "C0"), ("thick layer 2.25-4.75 km", thick, "C1")):
        for e in fam:
            for key, mark, off in (("v2", "o", 0.0), ("v22_fixed", "x", 0.03)):
                if np.isfinite(e[key]["start_m"]) and e[key]["accept"] > 0:
                    ax.plot([e["r_peak"] + off] * 2,
                            [e[key]["start_m"] / 1000.0, e[key]["end_m"] / 1000.0],
                            colour, lw=2.0 if key == "v2" else 1.0, alpha=0.8)
                    ax.plot(e["r_peak"] + off, 0.5 * (e[key]["start_m"] + e[key]["end_m"]) / 1000.0,
                            mark, color=colour, ms=5)
                elif key == "v2":
                    ax.plot(e["r_peak"], 3.5, "v", color="C3", ms=8)
    ax.set_ylim(0, 7)
    ax.set_xlabel("backscatter ratio $R$ (at 3 km for haze, at layer centre for the thick layer)")
    ax.set_ylabel("accepted molecular window  [km AGL]")
    ax.set_title("the window the gates end up choosing\n(bars = window span; red triangle = "
                 "night rejected)", fontsize=11)
    ax.grid(alpha=0.3)

    ax = axes[2]
    for label, fam, colour in (("haze", haze, "C0"), ("thick layer", thick, "C1")):
        r = [e["r_peak"] for e in fam]
        ax.plot(r, [e["v2"]["accept"] * 100 for e in fam], "o-", color=colour,
                label=f"{label}: v2 acceptance [%]")
        ax.plot(r, [e["v2"]["bias_pct"] for e in fam], "s--", color=colour, alpha=0.6,
                label=f"{label}: bias of accepted $C_L$ [%]")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel("backscatter ratio $R$")
    ax.set_ylabel("acceptance [%]   /   $C_L$ bias [%]")
    ax.set_title("(b) where the strict tier breaks", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")

    fig.suptitle("Contamination the gates cannot dodge: thick layer and column-filling haze",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_recovery(res_rec, path: Path) -> None:
    """The observed recovered-vs-kept offsets, reproduced by in-window haze."""
    plt = _mpl()
    centres = np.asarray(res_rec["centres_m"], float) / 1000.0
    lo, hi = res_rec["low_centre_m"] / 1000.0, res_rec["high_centre_m"] / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6))

    ax = axes[0]
    for row in res_rec["rows"]:
        ax.plot(row["dev_pct"], centres, "o-", ms=3.0, lw=1.3,
                label=f"{row['ratio_slope_pct_per_km']:+5.1f} %/km (AOD {row['aod']:.3f})")
    ax.axhline(lo, color="0.4", ls=":", lw=1.2)
    ax.axhline(hi, color="0.4", ls="--", lw=1.2)
    ax.axvline(0.0, color="C3", ls="--", lw=1.2)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]  (symlog above 10 %)")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title(f"$C_L$ ladders for column-filling haze\n(dotted {lo:.1f} km = a kept night's "
                 f"window, dashed {hi:.1f} km = a recovered one)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", title="signal/$p_{mol}$ slope, 2-6 km",
              title_fontsize=7)

    ax = axes[1]
    slope = [row["ratio_slope_pct_per_km"] for row in res_rec["rows"]]
    ax.plot(slope, [row["delta_height_pct"] for row in res_rec["rows"]], "o-", color="C0", lw=1.8,
            label=f"same night fitted {lo:.1f} -> {hi:.1f} km")
    ax.plot(slope, [row["dev_low_pct"] for row in res_rec["rows"]], "s--", color="C1", lw=1.4,
            label=f"hazy vs CLEAN night, both at {lo:.1f} km")
    colours = plt.cm.tab10(np.linspace(0, 1, len(res_rec["observed"])))
    # Stagger the labels in x: several stations sit within a few tenths of a percent of each
    # other (Aosta CHM15k +31.7 vs CL61 +32.1) and would overprint on a shared left margin.
    for k, ((name, value), colour) in enumerate(zip(res_rec["observed"].items(), colours)):
        ax.axhline(value, color=colour, lw=1.0, ls=":")
        ax.text(-19.5 + 5.0 * (k % 2), value + 0.6, f"{name} {value:+.1f} %", fontsize=7,
                color=colour, va="bottom")
    ax.axvline(res_rec["observed_ratio_slope_pct_per_km"], color="0.3", lw=1.4)
    ax.text(res_rec["observed_ratio_slope_pct_per_km"], -35,
            " measured on recovered\n nights: $-$19 %/km", fontsize=7.5, color="0.3")
    ax.set_ylim(-40, 60)
    ax.set_xlabel("signal / $p_{mol}$ slope over 2-6 km   [%/km]  (the measurable aerosol load)")
    ax.set_ylabel("$C_L$ offset   [%]")
    ax.set_title("(a) the two halves of the observed offset\nhave OPPOSITE signs", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[2]
    ax.plot(slope, [row["ladder_slope_pct_per_km"] for row in res_rec["rows"]], "o-", color="C2",
            lw=1.8)
    ax.axvline(res_rec["observed_ratio_slope_pct_per_km"], color="0.3", lw=1.4)
    ax.axhline(-15.0, color="C3", ls="--", lw=1.2)
    ax.text(-18.5, -15.5, "Payerne CHM15k measured $dC_L/dz$ = $-$15 %/km", fontsize=7.5,
            color="C3")
    ax.set_xlabel("signal / $p_{mol}$ slope over 2-6 km   [%/km]")
    ax.set_ylabel("$dC_L/dz$ of the ladder   [%/km]")
    ax.set_title("the simulated $C_L$ altitude gradient\nvs the measured aerosol load",
                 fontsize=11)
    ax.grid(alpha=0.3)

    fig.suptitle("Reproducing the observed recovered-vs-kept offsets with aerosol in the window",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_signature(res_sig, path: Path) -> None:
    plt = _mpl()
    z, atm = build_atmosphere()
    centres = np.asarray(res_sig["centres_m"], float) / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.6))

    styles = {"layer": "-", "haze": "--"}
    ax = axes[0]
    for s in res_sig["scenarios"]:
        kind = s["kind"]
        case = make_case(z, atm, **s["case_kwargs"])
        ratio = (case["rcs_true"] / z ** 2) / atm["p_mol"] / C_TRUE
        ax.plot(ratio, z / 1000.0, styles[kind], lw=1.3, label=s["label"])
    ax.set_xlim(0.8, 3.2)
    ax.set_ylim(0, 8)
    ax.set_xlabel("(signal / $p_{mol}$) / $C_{true}$   [-]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("the six synthetic atmospheres", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    ax = axes[1]
    ax.axvline(0.0, color="C3", ls="--", lw=1.4)
    for s in res_sig["scenarios"]:
        ax.plot(s["dev_pct"], centres, styles[s["kind"]], marker="o", ms=3.2, lw=1.3,
                label=s["label"])
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("(d) ladders: bump vs slope", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right")

    ax = axes[2]
    for s in res_sig["scenarios"]:
        ax.scatter(s["localisation"], s["bump_pct"], s=70,
                   marker="o" if s["kind"] == "layer" else "^", label=s["label"])
    ax.axvline(0.5, color="0.5", ls=":")
    ax.set_xlim(0, 1)
    ax.set_xlabel("localisation index  bump / (bump + linear span)   [-]")
    ax.set_ylabel("bump amplitude  [%]")
    ax.set_yscale("log")
    ax.set_title("(d) the discriminator\n(right = localised contamination,\nleft = monotonic "
                 "transmission)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")
    flat = [s for s in res_sig["scenarios"] if not np.isfinite(s["localisation"])]
    if flat:
        ax.text(0.03, 0.97, "not plotted: the two residual layers BELOW the window\n"
                            "(ladder flat to 0.0013 %, both terms zero ->\n"
                            "no localisation index is defined)",
                transform=ax.transAxes, fontsize=7.5, va="top",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    fig.suptitle("Signature separating in-window contamination from below-window transmission",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main(workers: int = 6) -> dict:
    print("=" * 96)
    results = dict(constants=dict(C_true=C_TRUE, S_aer=S_AER, sigma_night=SIGMA_NIGHT,
                                  n_profiles=N_PROFILES, sigma_bug_factor=SIGMA_BUG_FACTOR,
                                  ladder_half_m=LADDER_HALF_M,
                                  layer_thickness_m=LAYER_THICKNESS_M,
                                  v2_params=DEFAULT_PARAMS["eprof_v2"],
                                  v22_params=DEFAULT_PARAMS["eprof_v2.2"]))
    results["closure"] = run_closure()
    print("\n(a) LADDER STUDY")
    results["ladder"] = run_ladder_study()
    print("\nBIAS DECOMPOSITION (backscatter vs transmission, and the window-mean-ratio law)")
    results["decomposition"] = run_decomposition()
    print("\n(d) SIGNATURE STUDY")
    results["signature"] = run_signature_study()
    print("\nRECOVERED-vs-KEPT PROXY (the observed comparison, reproduced)")
    results["recovery"] = run_recovery_proxy()
    print("\n(b)+(c) GATE STUDY")
    results["gates"] = run_gate_study(workers=workers)
    for fam in ("thin", "thick", "haze"):
        results["gates"][f"summary_{fam}"] = family_summary(results["gates"], fam)
        print(f"\n  family '{fam}'  (R | AOD | signal/p_mol slope | v2: accept, n_eligible, modal "
              f"window, bias [min..max] | v2.2 shipped-sigma: accept, bias | v2.2 fixed-sigma: "
              f"accept, bias | tier-2 fired)")
        for e in results["gates"][f"summary_{fam}"]:
            v2, vs, vf = e["v2"], e["v22_shipped"], e["v22_fixed"]
            print(f"    {e['r_peak']:4.2f} {e['aod']:.4f} {e['ratio_slope_pct_per_km']:+6.1f} %/km"
                  f" | v2 {v2['accept']*100:3.0f} % {v2['n_eligible']:4.0f} "
                  f"{v2['start_m']/1000:4.2f}-{v2['end_m']/1000:4.2f}km({v2['n_mode']}/{e['n']}) "
                  f"{v2['bias_pct']:+7.2f} [{v2['bias_min_pct']:+6.1f}..{v2['bias_max_pct']:+6.1f}]"
                  f" | v2.2s {vs['accept']*100:3.0f} % {vs['bias_pct']:+7.2f}"
                  f" | v2.2f {vf['accept']*100:3.0f} % {vf['bias_pct']:+7.2f}"
                  f" | tier2 {e['tier2_fired_shipped']*100:3.0f}/{e['tier2_fired_fixed']*100:3.0f} %")

    print("\n(c) FORCED NOISE-TIER PROBE -- haze family")
    results["tier2_haze"] = run_tier2_probe(kind="haze")
    print("\n(c) FORCED NOISE-TIER PROBE -- thin layer at 3.5 km")
    results["tier2_layer"] = run_tier2_probe(kind="layer",
                                             loads=(1.05, 1.10, 1.20, 1.35, 1.50, 2.00, 3.00))

    z, atm = build_atmosphere()
    sigma_true = np.full(z.size, SIGMA_NIGHT)
    chi2 = {}
    for label, r in (("clean", 1.0), ("R1.20", 1.20), ("R1.50", 1.50), ("R2.00", 2.0)):
        case = make_case(z, atm, kind="layer", layer_centre_m=3500.0, r_peak=r, seed=7)
        chi2[label] = dict(shipped=chi2red_stats(z, case, atm, sigma_true * SIGMA_BUG_FACTOR),
                           fixed=chi2red_stats(z, case, atm, sigma_true))
    results["chi2red"] = chi2
    print("\nchi2red over the window grid (median | fraction <= 2.5):")
    for label, d in chi2.items():
        print(f"   {label:6s} shipped sigma: {d['shipped']['median']:9.3f} | "
              f"{d['shipped']['frac_below_2p5']*100:5.1f} %   "
              f"fixed sigma: {d['fixed']['median']:9.3f} | {d['fixed']['frac_below_2p5']*100:5.1f} %")

    fig_ladders(results["ladder"], FIGDIR / "fwd_aer_ladders.png")
    fig_gates(results["gates"], FIGDIR / "fwd_aer_gates.png")
    fig_saturation(results["gates"], FIGDIR / "fwd_aer_saturation.png")
    fig_signature(results["signature"], FIGDIR / "fwd_aer_signature.png")
    fig_recovery(results["recovery"], FIGDIR / "fwd_aer_recovery.png")

    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, float) and not np.isfinite(o):
            return None
        return o

    out = FIGDIR / "fwd_aer_results.json"
    out.write_text(json.dumps(_clean(results), indent=1), encoding="utf-8")
    print(f"\nJSON  -> {out}")
    for name in ("fwd_aer_ladders.png", "fwd_aer_gates.png", "fwd_aer_saturation.png",
                 "fwd_aer_signature.png", "fwd_aer_recovery.png"):
        print(f"figure -> {FIGDIR / name}")
    return results


if __name__ == "__main__":
    main()
