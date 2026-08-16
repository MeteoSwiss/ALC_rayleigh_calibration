# -*- coding: utf-8 -*-
"""MECHANISM SCAN 4 -- additive signal residual (b) and molecular-model error.

Both mechanisms tested here are "instrument or model" rather than "atmosphere": neither needs a
single aerosol particle to make the retrieved lidar constant depend on the height of the molecular
window. They are run through the SHIPPED retrieval via ``forward_model.retrieve`` -- nothing is
re-implemented -- on a synthetic profile built from a KNOWN ``C_true``.

MECHANISM A -- ADDITIVE RESIDUAL b
----------------------------------
The window fit is ``signal = a * p_mol + b`` with a FREE intercept. ``b`` is whatever additive
residual survives the instrument's own background/afterpulse/dark subtraction, in the units of
``signal = rcs / z^2``. With ``options.subtract_background = False`` (the operational default,
``options.json``) ``rayleigh_fit.calculate_lidar_constant`` divides the RAW ``rcs_mean`` by
``beta_tot``, so the fitted ``b`` is never removed and the pointwise constant carries

    C_L(z) / C_true - 1  =  b / (a * p_mol(z))                                          [A1]

Because ``p_mol(z) = beta_mol(z) * T2(z) / z^2`` falls by a factor ~5.9 between 3 and 6 km
(0.687 from the exp(-z/H) density decay, 4.0 from the 1/z^2, the rest transmission), a residual
worth -3.8 % of the signal at 3 km is worth -20 % at 6 km. That is an ALTITUDE-DEPENDENT
calibration built into the estimator, with the sign of ``b``.

MECHANISM B -- MOLECULAR MODEL ERROR
------------------------------------
The retrieval's molecular reference ``p_mol`` comes from a model. Operationally that model is the
US Standard Atmosphere 1976 (``options.json``: ``molecular_source = 'standard'``), NOT the CAMS
T/p of the night. Any error in ``p_mol`` that is not a pure multiplicative constant maps into an
altitude-dependent constant, since ``C_L(z) ~ rcs(z) / beta_tot(z)``. A pure multiplicative error
(e.g. a surface-pressure bias propagated hydrostatically) is degenerate with ``C_L`` itself and
produces NO slope -- only a SHAPE error does. The four probes are the ones the operational chain
can actually suffer: US-Std instead of CAMS, a +-2 K temperature bias, a +-5 hPa surface-pressure
error, and the AGL-vs-ASL height-reference confusion.

WHAT IS MEASURED
----------------
For every configuration the molecular window is FORCED at 13 centres from 3000 to 6000 m AGL
(half-width 490 m = 33 gates on the Payerne CHM15k L2 grid) and the ladder

    L(z_c) = C_L(z_c) / C_true - 1        [%]

is fitted in log space against the window centre, giving ``dC/dz`` in %/km over 3-6 km. That is
directly comparable with the measured network gradients (Payerne CHM15k -15 %/km).

Run::

    python rayleigh_availability/fwd_additive_molecular.py            # full scan (~2 min)
    python rayleigh_availability/fwd_additive_molecular.py --quick    # fewer noise seeds

Outputs: four PNGs in ``doc/reports/figs_altitude_audit/`` and a JSON of every number in
``<DATA>/rayleigh_availability/fwd_additive_molecular.json``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from calibration.rayleigh.atmosphere import (                       # noqa: E402
    DEFAULT_STANDARD_ATMOSPHERE,
    MOLECULAR_LIDAR_RATIO,
    calculate_molecular_properties,
    load_standard_atmosphere,
)
from rayleigh_availability.forward_model import (                   # noqa: E402
    PAYERNE_CHM15K,
    default_grid,
    default_options,
    exponential_aerosol,
    forward_signal,
    make_atmosphere,
    retrieve,
)

FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
CAMS_DIR = Path("D:/CAMS_run_v20")
OUT_JSON = DATA / "fwd_additive_molecular.json"

# Payerne 0-20000-0-06610 (the site the observational study anchors on).
PAY_LAT, PAY_LON, PAY_ALT = 46.8130, 6.9440, 490.0

# The ladder: window centres (m AGL) and half-width. 3-6 km is the band the task asks for and the
# band the operational search (range_start_m=2000, range_end_m=6000) can actually reach.
CENTRES_M = np.arange(3000.0, 6001.0, 250.0)
HALF_M = 490.0

# Measured Payerne CHM15k intercepts (median over nights, intercept_test_payerne.json,
# 36 v2.0-kept + 80 v2.2-recovered nights). Units of signal = rcs / z^2.
B_KEPT_MEASURED = -1.835e-4
B_RECOVERED_MEASURED = -6.375e-5


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def atmosphere_from_tp(z, temperature, pressure, *, wavelength_nm=1064.0, station_alt_m=PAY_ALT):
    """Wrap an arbitrary (T, p) pair into the dict :func:`forward_model.retrieve` consumes.

    Calls the repo's own ``calculate_molecular_properties`` (Bucholtz 1995), so ``p_mol`` and
    ``beta_mol`` are bit-for-bit what the pipeline would build from the same T/p.

    Parameters
    ----------
    z : ndarray            range grid, m AGL
    temperature : ndarray  K
    pressure : ndarray     Pa
    """
    mol = calculate_molecular_properties(np.asarray(temperature, float),
                                         np.asarray(pressure, float),
                                         np.asarray(z, float), wavelength_nm * 1e-9)
    return dict(z=z, altitude_asl=station_alt_m + z, temperature=np.asarray(temperature, float),
                pressure=np.asarray(pressure, float), beta_mol=mol.beta_mol,
                alpha_mol=mol.alpha_mol, p_mol=mol.p_mol, transmission_mol=mol.transmission,
                beta_att_mol=mol.beta_att_mol, wavelength_nm=float(wavelength_nm),
                station_alt_m=float(station_alt_m),
                lidar_ratio_mol=float(MOLECULAR_LIDAR_RATIO))


def ladder(z, rcs, atm_retrieval, *, subtract_background=False, centres=CENTRES_M, half=HALF_M):
    """Force the molecular window at each centre and return the retrieved C_L (shipped chain)."""
    options = default_options()
    options.subtract_background = bool(subtract_background)
    values = np.full(len(centres), np.nan)
    for k, centre in enumerate(centres):
        result = retrieve(z, rcs, atm_retrieval, options=options,
                          window=(float(centre) - half, float(centre) + half))
        if result["ok"]:
            values[k] = result["C_L"]
    return values


def ladder_slope(values, c_true, centres=CENTRES_M):
    """Log-space slope of the ladder in %/km (the metric the observational study reports)."""
    ratio = np.asarray(values, float) / float(c_true)
    good = np.isfinite(ratio) & (ratio > 0)
    if good.sum() < 3:
        return np.nan
    return float(np.polyfit(np.asarray(centres)[good] / 1000.0,
                            np.log(ratio[good]), 1)[0] * 100.0)


def dev_pct(values, c_true):
    return (np.asarray(values, float) / float(c_true) - 1.0) * 100.0


def at(centres, target):
    return int(np.argmin(np.abs(np.asarray(centres) - target)))


# ---------------------------------------------------------------------------
# 0. closure on THIS ladder grid
# ---------------------------------------------------------------------------
def run_closure(z, atm, c_true):
    """Pure molecular, noiseless, no residual, correct model: the ladder must return C_true."""
    rcs = forward_signal(z, c_true, atm["beta_mol"])
    values = ladder(z, rcs, atm)
    deviation = dev_pct(values, c_true)
    return dict(dev_pct=deviation.tolist(),
                max_abs_dev_pct=float(np.nanmax(np.abs(deviation))),
                slope_pct_per_km=ladder_slope(values, c_true))


# ---------------------------------------------------------------------------
# A. additive residual b
# ---------------------------------------------------------------------------
def scan_additive(z, atm, c_true, fractions=(-0.05, -0.02, -0.005, -0.001,
                                             0.001, 0.005, 0.02, 0.05)):
    """Ladder versus additive residual b, expressed as a fraction of the signal at 3 km.

    ``b`` is passed to ``forward_signal(..., background=b)``: a constant offset of the
    range-normalised signal ``rcs / z^2``, i.e. exactly the quantity the window fit returns as its
    intercept. Both the shipped default (``subtract_background=False``) and the existing fix
    (``True``) are run on the SAME synthetic profile.
    """
    i3 = at(z, 3000.0)
    signal_3km = float(c_true * atm["p_mol"][i3])
    out = []
    for fraction in fractions:
        b = fraction * signal_3km
        rcs = forward_signal(z, c_true, atm["beta_mol"], background=b)
        raw = ladder(z, rcs, atm, subtract_background=False)
        fix = ladder(z, rcs, atm, subtract_background=True)
        # closed-form prediction [A1] evaluated at each window centre
        pred = np.array([b / (c_true * atm["p_mol"][at(z, c)]) * 100.0 for c in CENTRES_M])
        out.append(dict(
            fraction_at_3km=float(fraction), b=float(b),
            raw_dev_pct=dev_pct(raw, c_true).tolist(),
            fix_dev_pct=dev_pct(fix, c_true).tolist(),
            pred_dev_pct=pred.tolist(),
            raw_slope=ladder_slope(raw, c_true), fix_slope=ladder_slope(fix, c_true),
            raw_at_3km=float(dev_pct(raw, c_true)[at(CENTRES_M, 3000)]),
            raw_at_6km=float(dev_pct(raw, c_true)[at(CENTRES_M, 6000)]),
            fix_at_3km=float(dev_pct(fix, c_true)[at(CENTRES_M, 3000)]),
            fix_at_6km=float(dev_pct(fix, c_true)[at(CENTRES_M, 6000)]),
            pred_slope=ladder_slope(c_true * (1 + pred / 100.0), c_true),
        ))
    return dict(signal_3km=signal_3km, cases=out)


def anchor_measured_b(z, atm, c_true):
    """The same ladder at the b actually measured on Payerne CHM15k nights."""
    i3 = at(z, 3000.0)
    signal_3km = float(c_true * atm["p_mol"][i3])
    out = {}
    for name, b in (("kept", B_KEPT_MEASURED), ("recovered", B_RECOVERED_MEASURED)):
        rcs = forward_signal(z, c_true, atm["beta_mol"], background=b)
        raw = ladder(z, rcs, atm, subtract_background=False)
        fix = ladder(z, rcs, atm, subtract_background=True)
        out[name] = dict(b=float(b), fraction_at_3km=float(b / signal_3km),
                         raw_dev_pct=dev_pct(raw, c_true).tolist(),
                         fix_dev_pct=dev_pct(fix, c_true).tolist(),
                         raw_slope=ladder_slope(raw, c_true),
                         fix_slope=ladder_slope(fix, c_true))
    return out


# ---------------------------------------------------------------------------
# A-bis. failure modes of subtract_background = True
# ---------------------------------------------------------------------------
def noisy_intercept(z, atm, c_true, n_seeds=300, sigma=None, b=0.0, seed=20260815):
    """Failure mode 1: the fitted intercept is itself a noisy estimate.

    With ``subtract_background=True`` the retrieval subtracts ``b_hat``, an estimate whose error
    is amplified by 1/p_mol exactly like the true b would be. Truth here has ``b = 0`` (nothing to
    correct), so any extra scatter is pure cost. Bias and scatter are reported SEPARATELY per
    window height, with the s.e.m. of the bias so significance is answerable.
    """
    sigma = PAYERNE_CHM15K["sigma_night_clean"] if sigma is None else sigma
    rng = np.random.default_rng(seed)
    raw = np.full((n_seeds, len(CENTRES_M)), np.nan)
    fix = np.full_like(raw, np.nan)
    for s in range(n_seeds):
        rcs = forward_signal(z, c_true, atm["beta_mol"], background=b,
                             noise_sigma=sigma, rng=rng)
        raw[s] = ladder(z, rcs, atm, subtract_background=False)
        fix[s] = ladder(z, rcs, atm, subtract_background=True)
    def summarise(arr):
        bias = dev_pct(np.nanmean(arr, axis=0), c_true)
        scatter = np.nanstd(arr, axis=0) / c_true * 100.0
        n = np.sum(np.isfinite(arr), axis=0)
        return dict(bias_pct=bias.tolist(), scatter_pct=scatter.tolist(),
                    sem_pct=(scatter / np.sqrt(np.maximum(n, 1))).tolist(),
                    slope=ladder_slope(np.nanmean(arr, axis=0), c_true))
    return dict(sigma=float(sigma), b=float(b), n_seeds=int(n_seeds),
                raw=summarise(raw), fix=summarise(fix))


def aerosol_curvature(z, atm, c_true, beta_surface=2.0e-6, scale_height_m=1500.0,
                      lidar_ratio_sr=52.0, b=0.0):
    """Failure mode 2: a real aerosol gradient pushes curvature into the free intercept.

    With an exponential haze the true signal is NOT a straight line in p_mol, so ``linregress``
    returns a non-zero ``b_hat`` that is pure aerosol, not instrument. Subtracting it removes REAL
    signal. Truth has no additive residual (``b = 0`` by default), so with a perfect estimator both
    modes would agree; whatever separates them is the cost of the fix.
    """
    beta_aer, ext_aer = exponential_aerosol(z, beta_surface, scale_height_m, lidar_ratio_sr)
    rcs = forward_signal(z, c_true, atm["beta_mol"], beta_aer=beta_aer, ext_aer=ext_aer,
                         background=b)
    raw = ladder(z, rcs, atm, subtract_background=False)
    fix = ladder(z, rcs, atm, subtract_background=True)
    signal = rcs / z ** 2
    band = (z >= 2000) & (z <= 6000)
    ratio_slope = float(np.polyfit(z[band] / 1000.0,
                                   np.log(signal[band] / atm["p_mol"][band]), 1)[0] * 100.0)

    # What does the FREE window fit report as its intercept when the truth has b = 0? This is the
    # quantity intercept_test.py measures on real nights as `rel_pred = median(b/(a*p_mol))`.
    options = default_options()
    fitted = []
    for centre in CENTRES_M:
        r = retrieve(z, rcs, atm, options=options,
                     window=(float(centre) - HALF_M, float(centre) + HALF_M))
        if not r["ok"]:
            fitted.append((np.nan, np.nan))
            continue
        win = (z >= centre - HALF_M) & (z <= centre + HALF_M)
        rel_pred = float(np.median(r["intercept"] / (r["slope"] * atm["p_mol"][win])))
        fitted.append((float(r["intercept"]), rel_pred))
    return dict(beta_surface=float(beta_surface), scale_height_m=float(scale_height_m),
                b=float(b),
                fitted_intercept=[f[0] for f in fitted],
                rel_pred_pct=[f[1] * 100.0 for f in fitted],
                scattering_ratio_3km=float(1.0 + beta_aer[at(z, 3000)] / atm["beta_mol"][at(z, 3000)]),
                signal_over_pmol_slope_2_6km=ratio_slope,
                raw_dev_pct=dev_pct(raw, c_true).tolist(),
                fix_dev_pct=dev_pct(fix, c_true).tolist(),
                raw_slope=ladder_slope(raw, c_true), fix_slope=ladder_slope(fix, c_true))


def intercept_vs_aerosol(z, atm, c_true, beta_list=(5.0e-9, 5.0e-8, 1.0e-7, 1.5e-7, 2.5e-7,
                                                    5.0e-7)):
    """Is the intercept measured on real nights instrumental at all?

    ``intercept_test.py`` measures ``rel_pred = median(b_hat/(a*p_mol))`` = -6.32 % (kept, window
    mid 3656 m) and -6.94 % (recovered, 4855 m) on Payerne CHM15k. Here the truth contains NO
    additive residual whatsoever (``b_true = 0``) and only an exponential haze, so every non-zero
    ``b_hat`` the free fit returns is ATMOSPHERE that the two-parameter model cannot represent.
    Matching the measured ``rel_pred`` with a plausible haze would mean the "instrumental residual"
    is mostly a misread aerosol.
    """
    out = []
    options = default_options()
    for beta_surface in beta_list:
        beta_aer, ext_aer = exponential_aerosol(z, beta_surface, 1500.0, 52.0)
        rcs = forward_signal(z, c_true, atm["beta_mol"], beta_aer=beta_aer, ext_aer=ext_aer)
        signal = rcs / z ** 2
        band = (z >= 2000) & (z <= 6000)
        rel_pred, values = [], []
        for centre in CENTRES_M:
            r = retrieve(z, rcs, atm, options=options,
                         window=(float(centre) - HALF_M, float(centre) + HALF_M))
            win = (z >= centre - HALF_M) & (z <= centre + HALF_M)
            rel_pred.append(float(np.median(r["intercept"] / (r["slope"] * atm["p_mol"][win])))
                            * 100.0 if r["ok"] else np.nan)
            values.append(r["C_L"] if r["ok"] else np.nan)
        out.append(dict(
            beta_surface=float(beta_surface),
            scattering_ratio_3km=float(1.0 + beta_aer[at(z, 3000)] / atm["beta_mol"][at(z, 3000)]),
            signal_over_pmol_slope_2_6km=float(
                np.polyfit(z[band] / 1000.0, np.log(signal[band] / atm["p_mol"][band]), 1)[0] * 100.0),
            rel_pred_pct=rel_pred,
            ladder_slope=ladder_slope(np.array(values), c_true)))
    return out


def wide_range_intercept(z, rcs, atm, *, z_lo=3000.0, z_hi=15000.0):
    """PROPOSED estimator (not shipped): estimate b on a WIDE range instead of inside the window.

    ``linregress`` inside a 980 m window extrapolates to ``p_mol = 0`` from a lever arm of only a
    factor ~2 in ``p_mol``, which is why ``b_hat`` is so noisy. Regressing over 3-15 km spans three
    decades of ``p_mol`` and pins the intercept on the far range, where the molecular signal is
    negligible and whatever is left IS the additive residual. The offset is then removed from
    ``rcs`` BEFORE the shipped retrieval runs (which therefore stays at
    ``subtract_background=False`` and needs no code change inside ``calculate_lidar_constant``).
    """
    signal = rcs / z ** 2
    band = (z >= z_lo) & (z <= z_hi) & np.isfinite(signal)
    x, y = atm["p_mol"][band], signal[band]
    a_hat, b_hat = np.polyfit(x, y, 1)
    return float(b_hat), float(a_hat)


def decision_experiment(z, atm, c_true, *, b_true, sigma, n_seeds=300, seed=99,
                        beta_surface=0.0):
    """Bias / scatter / RMSE of the three options, at one (b, sigma, aerosol) operating point.

    Options compared at every ladder height:
      raw        the shipped default (``subtract_background=False``)
      fix        the existing switch (``subtract_background=True``, in-window intercept)
      wide       proposed: intercept fitted over 3-15 km, removed from ``rcs`` beforehand
      far        proposed: intercept fitted over 8-15 km only (aerosol-free far range)

    RMSE = sqrt(bias^2 + scatter^2) over the seeds is the number a calibration owner should act on:
    the switch is only worth flipping where it lowers RMSE.
    """
    if beta_surface:
        beta_aer, ext_aer = exponential_aerosol(z, beta_surface, 1500.0, 52.0)
    else:
        beta_aer = ext_aer = None
    rng = np.random.default_rng(seed)
    keys = ("raw", "fix", "wide", "far")
    store = {k: np.full((n_seeds, len(CENTRES_M)), np.nan) for k in keys}
    b_hats = {"wide": np.full(n_seeds, np.nan), "far": np.full(n_seeds, np.nan)}
    for s in range(n_seeds):
        rcs = forward_signal(z, c_true, atm["beta_mol"], beta_aer=beta_aer, ext_aer=ext_aer,
                             background=b_true, noise_sigma=sigma, rng=rng)
        store["raw"][s] = ladder(z, rcs, atm, subtract_background=False)
        store["fix"][s] = ladder(z, rcs, atm, subtract_background=True)
        for key, z_lo in (("wide", 3000.0), ("far", 8000.0)):
            b_hat, _ = wide_range_intercept(z, rcs, atm, z_lo=z_lo)
            b_hats[key][s] = b_hat
            store[key][s] = ladder(z, rcs - b_hat * z ** 2, atm, subtract_background=False)
    out = dict(b_true=float(b_true), sigma=float(sigma), n_seeds=int(n_seeds),
               beta_surface=float(beta_surface),
               b_hat_wide_mean=float(np.nanmean(b_hats["wide"])),
               b_hat_wide_std=float(np.nanstd(b_hats["wide"])),
               b_hat_far_mean=float(np.nanmean(b_hats["far"])),
               b_hat_far_std=float(np.nanstd(b_hats["far"])))
    for key, arr in store.items():
        bias = dev_pct(np.nanmean(arr, axis=0), c_true)
        scatter = np.nanstd(arr, axis=0) / c_true * 100.0
        out[key] = dict(bias_pct=bias.tolist(), scatter_pct=scatter.tolist(),
                        rmse_pct=np.sqrt(bias ** 2 + scatter ** 2).tolist(),
                        slope=ladder_slope(np.nanmean(arr, axis=0), c_true))
    return out


# ---------------------------------------------------------------------------
# B. molecular model error
# ---------------------------------------------------------------------------
def molecular_model_scan(z, c_true, night, cams_file, *, wavelength_nm=1064.0):
    """Truth = CAMS T/p of a real Payerne night; retrieve with a WRONG molecular model.

    The truth signal is pure molecular (no aerosol, no additive residual, no noise) so the ladder
    isolates the model error alone. Each retrieval atmosphere is a plausible operational error:

      cams        the same CAMS profile          -> closure control, must return ~0
      usstd       US Standard 1976 at ASL        -> the OPERATIONAL default (molecular_source='standard')
      T+2K/T-2K   truth temperature +- 2 K       -> analysis temperature bias
      p+5/-5hPa   truth pressure * (1 +- 5/p_sfc) -> surface-pressure error, hydrostatic (uniform)
      agl         US Standard sampled at z AGL   -> station altitude forgotten (490 m at Payerne)
      agl_cams    CAMS interpolated at z AGL     -> same confusion on the CAMS path
    """
    t_start = np.datetime64(f"{night}T20:00")
    t_end = np.datetime64(night) + np.timedelta64(1, "D") + np.timedelta64(4, "h")
    truth = make_atmosphere(z, station_alt_m=PAY_ALT, wavelength_nm=wavelength_nm,
                            cams_file=cams_file, latitude=PAY_LAT, longitude=PAY_LON,
                            t_start=t_start, t_end=t_end)
    rcs = forward_signal(z, c_true, truth["beta_mol"])

    p_sfc_hpa = float(truth["pressure"][0]) / 100.0
    std_asl = load_standard_atmosphere(DEFAULT_STANDARD_ATMOSPHERE, PAY_ALT + z)
    std_agl = load_standard_atmosphere(DEFAULT_STANDARD_ATMOSPHERE, z)
    cams_agl = make_atmosphere(z, station_alt_m=0.0, wavelength_nm=wavelength_nm,
                               cams_file=cams_file, latitude=PAY_LAT, longitude=PAY_LON,
                               t_start=t_start, t_end=t_end)

    variants = {
        "cams (closure)": truth,
        "US-Std 1976 (operational)": atmosphere_from_tp(z, std_asl.temperature, std_asl.pressure,
                                                        wavelength_nm=wavelength_nm),
        "T +2 K": atmosphere_from_tp(z, truth["temperature"] + 2.0, truth["pressure"],
                                     wavelength_nm=wavelength_nm),
        "T -2 K": atmosphere_from_tp(z, truth["temperature"] - 2.0, truth["pressure"],
                                     wavelength_nm=wavelength_nm),
        "p_sfc +5 hPa": atmosphere_from_tp(z, truth["temperature"],
                                           truth["pressure"] * (1.0 + 5.0 / p_sfc_hpa),
                                           wavelength_nm=wavelength_nm),
        "p_sfc -5 hPa": atmosphere_from_tp(z, truth["temperature"],
                                           truth["pressure"] * (1.0 - 5.0 / p_sfc_hpa),
                                           wavelength_nm=wavelength_nm),
        "US-Std at AGL (no station alt)": atmosphere_from_tp(z, std_agl.temperature,
                                                            std_agl.pressure,
                                                            wavelength_nm=wavelength_nm),
        "CAMS at AGL (no station alt)": cams_agl,
    }

    out = {}
    for name, atm_retr in variants.items():
        values = ladder(z, rcs, atm_retr)
        deviation = dev_pct(values, c_true)
        ratio = atm_retr["p_mol"] / truth["p_mol"]
        band = (z >= 3000) & (z <= 6000)
        out[name] = dict(
            dev_pct=deviation.tolist(),
            slope_pct_per_km=ladder_slope(values, c_true),
            dev_at_3km=float(deviation[at(CENTRES_M, 3000)]),
            dev_at_6km=float(deviation[at(CENTRES_M, 6000)]),
            pmol_ratio_3km=float(ratio[at(z, 3000)]), pmol_ratio_6km=float(ratio[at(z, 6000)]),
            pmol_ratio_slope_pct_per_km=float(
                np.polyfit(z[band] / 1000.0, np.log(ratio[band]), 1)[0] * 100.0),
            pmol_ratio_profile=ratio.tolist(),
        )
    return dict(night=night, cams_file=str(cams_file), p_sfc_hpa=p_sfc_hpa,
                truth_T_3km=float(truth["temperature"][at(z, 3000)]),
                variants=out)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_additive(res, z, atm, c_true, path):
    plt = _mpl()
    y = CENTRES_M / 1000.0
    fig, axes = plt.subplots(1, 3, figsize=(16.8, 6.2))

    ax = axes[0]
    i3 = at(z, 3000.0)
    signal = c_true * atm["p_mol"]
    ax.plot(signal, z / 1000.0, color="C0", lw=1.8, label="signal moléculaire $C_{true}\\,p_{mol}$")
    for frac, colour in ((0.001, "C2"), (0.005, "C1"), (0.02, "C4"), (0.05, "C3")):
        ax.axvline(frac * signal[i3], color=colour, ls=":", lw=1.2,
                   label=f"|b| = {frac*100:g} % du signal à 3 km")
    ax.axvline(abs(B_KEPT_MEASURED), color="k", ls="--", lw=1.6,
               label=f"|b| mesuré Payerne (méd. nuits gardées)\n= {abs(B_KEPT_MEASURED):.2e} "
                     f"= {abs(B_KEPT_MEASURED)/signal[i3]*100:.1f} % à 3 km")
    ax.set_xscale("log")
    ax.set_xlim(1e-6, 1e-1)
    ax.set_ylim(0, 8)
    ax.axhspan(3.0, 6.0, color="C2", alpha=0.07)
    ax.set_xlabel("signal $rcs/z^2$  [unités de signal]")
    ax.set_ylabel("portée au-dessus de l'instrument  [km AGL]")
    ax.set_title("Où se situe le résidu additif b\npar rapport au signal moléculaire", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.4, loc="upper right")

    ax = axes[1]
    ax.axvline(0.0, color="k", lw=1.0)
    palette = {0.001: "C2", 0.005: "C1", 0.02: "C4", 0.05: "C3"}
    for case in res["additive"]["cases"]:
        frac = case["fraction_at_3km"]
        colour = palette[abs(frac)]
        ls = "-" if frac > 0 else "--"
        ax.plot(case["raw_dev_pct"], y, ls, color=colour, lw=1.7, marker="o", ms=3.5,
                label=f"b = {frac*100:+g} % @3 km : {case['raw_slope']:+.2f} %/km")
    anchor = res["anchor"]["kept"]
    ax.plot(anchor["raw_dev_pct"], y, "-", color="k", lw=2.4, marker="s", ms=4,
            label=f"b mesuré Payerne ({anchor['fraction_at_3km']*100:+.1f} % @3 km) : "
                  f"{anchor['raw_slope']:+.2f} %/km")
    ax.set_xlim(-62, 45)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
    ax.set_title("A. Chaîne livrée (subtract_background = False)\n"
                 "le résidu b n'est PAS retiré du signal", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.0, loc="lower left")

    ax = axes[2]
    ax.axvline(0.0, color="k", lw=1.0)
    for case in res["additive"]["cases"]:
        frac = case["fraction_at_3km"]
        ax.plot(case["fix_dev_pct"], y, "-" if frac > 0 else "--",
                color=palette[abs(frac)], lw=1.7, marker="o", ms=3.5,
                label=f"b = {frac*100:+g} % @3 km : {case['fix_slope']:+.3f} %/km")
    ax.plot(anchor["fix_dev_pct"], y, "-", color="k", lw=2.4, marker="s", ms=4,
            label=f"b mesuré Payerne : {anchor['fix_slope']:+.3f} %/km")
    ax.set_xlim(-1.2, 1.2)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
    ax.set_title("A-bis. Avec subtract_background = True\n(noter l'échelle $\\pm$1,2 % au lieu de "
                 "$\\pm$45 %)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.0, loc="lower left")

    fig.suptitle("Mécanisme A — un résidu additif b biaise $C_L$ de $b/(a\\,p_{mol}(z))$ : "
                 "l'erreur croît avec l'altitude de la fenêtre", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_failure_modes(res, path):
    plt = _mpl()
    y = CENTRES_M / 1000.0
    fig, axes = plt.subplots(1, 4, figsize=(21.5, 5.9))

    for ax, key, title in ((axes[0], "clean", "nuit propre  $\\sigma$ = 1,2e-4"),
                           (axes[1], "noisy", "nuit bruitée  $\\sigma$ = 5,0e-4")):
        block = res["noisy_intercept"][key]
        for mode, colour, label in (("raw", "C0", "subtract_background = False (livré)"),
                                    ("fix", "C3", "subtract_background = True")):
            bias = np.array(block[mode]["bias_pct"])
            scatter = np.array(block[mode]["scatter_pct"])
            sem = np.array(block[mode]["sem_pct"])
            offset = 0.0 if mode == "raw" else 0.03
            ax.errorbar(bias, y + offset, xerr=scatter, fmt="none", ecolor=colour, alpha=0.35,
                        capsize=3, lw=1.0)
            ax.errorbar(bias, y + offset, xerr=sem, fmt="o", color=colour, ms=4, capsize=2,
                        lw=1.4, label=f"{label}\n  biais {block[mode]['slope']:+.2f} %/km, "
                                      f"dispersion {scatter[0]:.2f}→{scatter[-1]:.2f} %")
        ax.axvline(0.0, color="k", lw=1.0)
        ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
        ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
        ax.set_title(f"Mode d'échec 1 — l'ordonnée à l'origine estimée est bruitée\n"
                     f"{title}, {block['n_seeds']} tirages, vérité b = 0", fontsize=10.5)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.2, loc="lower left")

    ax = axes[2]
    ax.axvline(0.0, color="k", lw=1.0)
    for case, colour in zip(res["aerosol_curvature"], ("C0", "C1", "C3")):
        lab = (f"$\\beta_{{aer}}(0)$ = {case['beta_surface']:.0e} m$^{{-1}}$sr$^{{-1}}$, "
               f"R(3 km) = {case['scattering_ratio_3km']:.2f}")
        ax.plot(case["raw_dev_pct"], y, "-", color=colour, lw=1.8, marker="o", ms=3.5,
                label=f"{lab}\n  False : {case['raw_slope']:+.2f} %/km")
        ax.plot(case["fix_dev_pct"], y, "--", color=colour, lw=1.8, marker="s", ms=3.5,
                label=f"  True  : {case['fix_slope']:+.2f} %/km")
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]   (échelle symlog)")
    ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
    ax.set_title("Mode d'échec 2 — un gradient d'aérosol réel\nmet de la courbure dans b "
                 "(vérité : b = 0)", fontsize=10.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.8, loc="lower left")

    ax = axes[3]
    colours = plt.cm.viridis(np.linspace(0.05, 0.9, len(res["intercept_vs_aerosol"])))
    for case, colour in zip(res["intercept_vs_aerosol"], colours):
        ax.plot(case["rel_pred_pct"], y, "-", color=colour, lw=1.7, marker="o", ms=3.0,
                label=f"R(3 km) = {case['scattering_ratio_3km']:.2f}, "
                      f"signal/$p_{{mol}}$ {case['signal_over_pmol_slope_2_6km']:+.1f} %/km")
    ax.plot([-6.32], [3.656], "*", color="#d62728", ms=17, mec="k", mew=0.6,
            label="Payerne mesuré, nuits gardées (n = 36)")
    ax.plot([-6.94], [4.855], "P", color="#d62728", ms=11, mec="k", mew=0.6,
            label="Payerne mesuré, nuits récupérées (n = 80)")
    ax.axvline(0.0, color="k", lw=1.0)
    ax.set_xlim(-32, 4)
    ax.set_xlabel("ordonnée à l'origine ajustée, en relatif\n"
                  "médiane $\\hat b/(a\\,p_{mol})$ dans la fenêtre   [%]")
    ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
    ax.set_title("L'« ordonnée à l'origine » est-elle instrumentale ?\n"
                 "vérité : b = 0, uniquement de la brume", fontsize=10.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.6, loc="lower left")

    fig.suptitle("Modes d'échec de subtract_background = True : bruit sur $\\hat b$, absorption "
                 "d'un vrai gradient d'aérosol — et $\\hat b$ n'est même pas instrumental",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_decision(res, path):
    """RMSE of the three options versus window height -- the figure the recommendation rests on."""
    plt = _mpl()
    y = CENTRES_M / 1000.0
    tags = list(res["decision"])
    fig, axes = plt.subplots(1, len(tags), figsize=(4.3 * len(tags), 6.0), sharey=True)
    style = {"raw": ("C0", "-", "livré (subtract_background = False)"),
             "fix": ("C3", "--", "commutateur existant (True, b ajusté dans la fenêtre)"),
             "wide": ("C2", "-", "proposé : b ajusté sur 3-15 km, retiré de rcs"),
             "far": ("C4", "-.", "proposé : b ajusté sur 8-15 km (loin, sans aérosol)")}
    for ax, tag in zip(np.atleast_1d(axes), tags):
        block = res["decision"][tag]
        for key, (colour, ls, label) in style.items():
            ax.plot(block[key]["rmse_pct"], y, ls, color=colour, lw=1.9, marker="o", ms=3.5,
                    label=label)
            ax.plot(np.abs(block[key]["bias_pct"]), y, ls, color=colour, lw=0.9, alpha=0.45)
        ax.set_xscale("log")
        ax.set_xlim(0.02, 200)
        ax.set_xlabel("erreur sur $C_L$   [%]\n(trait épais RMSE, trait fin |biais|)")
        ax.set_title(f"{tag}\n$\\sigma$ = {block['sigma']:.1e}, "
                     f"b = {block['b_true']:.2e}", fontsize=10)
        ax.grid(alpha=0.3, which="both")
    np.atleast_1d(axes)[0].set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
    np.atleast_1d(axes)[0].legend(fontsize=7.2, loc="lower right")
    fig.suptitle(f"Décision : corriger le résidu additif abaisse le biais mais fait exploser la "
                 f"variance si b est ajusté DANS la fenêtre "
                 f"({res['decision'][tags[0]]['n_seeds']} tirages)", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_molecular(res, z, path):
    plt = _mpl()
    y = CENTRES_M / 1000.0
    nights = res["molecular"]
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 6.2))

    style = {
        "cams (closure)": ("k", "-"),
        "US-Std 1976 (operational)": ("C3", "-"),
        "T +2 K": ("C0", "-"),
        "T -2 K": ("C0", "--"),
        "p_sfc +5 hPa": ("C2", "-"),
        "p_sfc -5 hPa": ("C2", "--"),
        "US-Std at AGL (no station alt)": ("C4", "-"),
        "CAMS at AGL (no station alt)": ("C1", "-"),
    }

    # --- panel 1: the p_mol shape error of every probe, as a profile ------------------------
    ax = axes[0]
    winter = nights["winter"]
    for name, block in winter["variants"].items():
        colour, ls = style[name]
        ax.plot(np.array(block["pmol_ratio_profile"]), z / 1000.0, ls, color=colour, lw=1.6,
                label=f"{name} : $dp_{{mol}}/dz$ = "
                      f"{block['pmol_ratio_slope_pct_per_km']:+.2f} %/km")
    summer = nights["summer"]["variants"]["US-Std 1976 (operational)"]
    ax.plot(np.array(summer["pmol_ratio_profile"]), z / 1000.0, ":", color="C3", lw=2.0,
            label=f"US-Std, été {nights['summer']['night']} : "
                  f"{summer['pmol_ratio_slope_pct_per_km']:+.2f} %/km")
    ax.text(0.02, 0.985, "$C_L \\propto 1/\\beta_{mol}$ : la pente de $C_L$\nest l'opposé de "
                         "celle de ce rapport", transform=ax.transAxes, fontsize=7.6, va="top",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.axvline(1.0, color="k", lw=1.0)
    ax.axhspan(3.0, 6.0, color="C2", alpha=0.08)
    ax.set_xlim(0.90, 1.10)
    ax.set_ylim(0, 8)
    ax.set_xlabel("$p_{mol}^{\\rm modèle\\ faux} / p_{mol}^{\\rm vérité\\ CAMS}$   [-]")
    ax.set_ylabel("altitude  [km AGL]")
    ax.set_title("Erreur de modèle moléculaire, en profil\n(vérité = CAMS Payerne, hiver "
                 "sauf mention)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.6, loc="lower left")

    for ax, (tag, night) in zip(axes[1:], nights.items()):
        ax.axvline(0.0, color="k", lw=1.0)
        for name, block in night["variants"].items():
            colour, ls = style[name]
            ax.plot(block["dev_pct"], y, ls, color=colour, lw=1.8, marker="o", ms=3.0,
                    label=f"{name} : {block['slope_pct_per_km']:+.2f} %/km")
        ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
        ax.set_ylabel("centre de la fenêtre moléculaire  [km AGL]")
        ax.set_title(f"B. Vérité = CAMS {night['night']} ({tag})\nrestitution avec un modèle faux",
                     fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.8, loc="lower left")

    fig.suptitle("Mécanisme B — erreur du modèle moléculaire : seule une erreur de FORME crée une "
                 "pente ; une erreur multiplicative est dégénérée avec $C_L$", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def fig_summary(res, path):
    """Every mechanism on one axis against the measured network gradients."""
    plt = _mpl()
    entries = []
    anchor = res["anchor"]["kept"]
    req = res["required_b"]
    entries.append(("A — résidu additif b, amplitude MESURÉE à Payerne\n"
                    f"(b = {anchor['b']:.2e} = {anchor['fraction_at_3km']*100:.1f} % du signal @3 km)",
                    anchor["raw_slope"], "#c0392b"))
    for case in res["additive"]["cases"]:
        if case["fraction_at_3km"] in (-0.05, 0.05):
            entries.append((f"A — b = {case['fraction_at_3km']*100:+g} % du signal @3 km",
                            case["raw_slope"], "#e67e22"))
    entries.append((f"A — b qu'il FAUDRAIT pour {req['target_slope_pct_per_km']:+.1f} %/km "
                    f"({req['ratio_to_measured']:.1f}x le b mesuré)",
                    req["target_slope_pct_per_km"], "#7f8c8d"))
    entries.append(("A-bis — même b, subtract_background = True", anchor["fix_slope"], "#27ae60"))
    for case in res["intercept_vs_aerosol"]:
        if 1.05 < case["scattering_ratio_3km"] < 1.15 or case["scattering_ratio_3km"] > 1.9:
            entries.append((f"contrôle AÉROSOL (scan 2/3) — brume R(3 km) = "
                            f"{case['scattering_ratio_3km']:.2f}, "
                            f"signal/$p_{{mol}}$ {case['signal_over_pmol_slope_2_6km']:+.0f} %/km",
                            case["ladder_slope"], "#8e44ad"))
    for tag, night in res["molecular"].items():
        entries.append((f"B — US-Std au lieu de CAMS ({tag} {night['night']}) = configuration "
                        f"OPÉRATIONNELLE",
                        night["variants"]["US-Std 1976 (operational)"]["slope_pct_per_km"], "#2980b9"))
    night = res["molecular"]["winter"]
    for name, label in (("T +2 K", "B — biais de température +2 K (hiver)"),
                        ("p_sfc +5 hPa", "B — erreur de pression de surface +5 hPa (hiver)"),
                        ("US-Std at AGL (no station alt)",
                         "B — altitude de station oubliée, AGL au lieu d'ASL (hiver)")):
        entries.append((label, night["variants"][name]["slope_pct_per_km"], "#16a085"))

    labels = [e[0] for e in entries][::-1]
    values = [e[1] for e in entries][::-1]
    colours = [e[2] for e in entries][::-1]

    fig, ax = plt.subplots(figsize=(15.5, 7.4))
    positions = np.arange(len(values))
    ax.barh(positions, values, color=colours, alpha=0.88, height=0.62)
    for p, v in zip(positions, values):
        ax.text(v + (0.35 if v >= 0 else -0.35), p, f"{v:+.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8.6, fontweight="bold")
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8.3)
    ax.set_xlim(-26, 14)
    for target, lab, colour in ((-15.3, "Payerne CHM15k mesuré (−15,3 %/km, n = 116)", "#b30000"),
                                (+7.0, "Aosta / Lindenberg : signe inverse", "#00609b")):
        ax.axvline(target, color=colour, ls="--", lw=1.5)
        ax.text(target - 0.35, len(values) - 0.55, lab, color=colour, fontsize=8.4,
                rotation=90, va="top", ha="right")
    ax.axvline(0.0, color="k", lw=1.0)
    ax.set_xlabel("pente de l'échelle  $d\\ln C_L/dz$  entre 3 et 6 km   [%/km]   "
                  "(fenêtres forcées, demi-largeur 490 m, grille Payerne CHM15k)")
    ax.set_title("Bilan du scan 4 : amplitude simulée de chaque mécanisme, face au gradient mesuré",
                 fontsize=12.5)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    quick = "--quick" in sys.argv
    n_seeds = 60 if quick else 300

    z = default_grid()
    c_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAY_ALT, wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])

    res = {}
    res["grid"] = dict(centres_m=CENTRES_M.tolist(), half_width_m=HALF_M, C_true=c_true,
                       n_gates_window=int(np.sum((z >= 3000 - HALF_M) & (z <= 3000 + HALF_M))))

    print("=== CLOSURE on the 3-6 km ladder grid ===")
    res["closure"] = run_closure(z, atm, c_true)
    print(f"  max |C_L/C_true - 1| = {res['closure']['max_abs_dev_pct']:.5f} %   "
          f"residual slope = {res['closure']['slope_pct_per_km']:+.5f} %/km")

    print("\n=== A. additive residual b ===")
    res["additive"] = scan_additive(z, atm, c_true)
    print(f"  signal(3 km) = {res['additive']['signal_3km']:.4e}")
    print(f"  {'b/signal@3km':>13} {'b':>11} {'C@3km':>9} {'C@6km':>9} {'slope':>9} "
          f"{'pred slope':>11} | {'fix@3km':>8} {'fix@6km':>8} {'fix slope':>10}")
    for case in res["additive"]["cases"]:
        print(f"  {case['fraction_at_3km']*100:+12.1f}% {case['b']:+11.3e} "
              f"{case['raw_at_3km']:+8.2f}% {case['raw_at_6km']:+8.2f}% "
              f"{case['raw_slope']:+8.3f} {case['pred_slope']:+10.3f} | "
              f"{case['fix_at_3km']:+7.3f}% {case['fix_at_6km']:+7.3f}% {case['fix_slope']:+9.4f}")

    res["anchor"] = anchor_measured_b(z, atm, c_true)
    print("\n  measured Payerne CHM15k intercepts:")
    for name, block in res["anchor"].items():
        print(f"   {name:>10}: b = {block['b']:+.3e} = {block['fraction_at_3km']*100:+.2f} % of "
              f"signal@3km -> raw slope {block['raw_slope']:+.2f} %/km, "
              f"fix slope {block['fix_slope']:+.4f} %/km")

    # How large would b have to be to reproduce the MEASURED Payerne gradient of -15.3 %/km?
    i3 = at(z, 3000.0)
    signal_3km = float(c_true * atm["p_mol"][i3])
    target = -15.3
    lo, hi = -0.40, -0.0001          # b as a fraction of the signal at 3 km
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        rcs = forward_signal(z, c_true, atm["beta_mol"], background=mid * signal_3km)
        if ladder_slope(ladder(z, rcs, atm), c_true) < target:
            lo = mid
        else:
            hi = mid
    res["required_b"] = dict(target_slope_pct_per_km=target,
                             fraction_at_3km=float(mid), b=float(mid * signal_3km),
                             ratio_to_measured=float(mid * signal_3km / B_KEPT_MEASURED))
    print(f"  to reach the measured {target:+.1f} %/km you would need b = "
          f"{mid*signal_3km:+.3e} = {mid*100:+.1f} % of the signal at 3 km, i.e. "
          f"{mid*signal_3km/B_KEPT_MEASURED:.2f}x the measured Payerne intercept")

    print("\n=== A-bis. failure modes of subtract_background = True ===")
    res["noisy_intercept"] = {
        "clean": noisy_intercept(z, atm, c_true, n_seeds=n_seeds,
                                 sigma=PAYERNE_CHM15K["sigma_night_clean"]),
        "noisy": noisy_intercept(z, atm, c_true, n_seeds=n_seeds,
                                 sigma=PAYERNE_CHM15K["sigma_night_noisy"], seed=7),
    }
    for key, block in res["noisy_intercept"].items():
        for mode in ("raw", "fix"):
            sc = block[mode]["scatter_pct"]
            print(f"  {key:>5} sigma={block['sigma']:.1e} {mode:>3}: bias slope "
                  f"{block[mode]['slope']:+.3f} %/km, scatter {sc[0]:.3f} % @3 km -> "
                  f"{sc[-1]:.3f} % @6 km")

    res["aerosol_curvature"] = [
        aerosol_curvature(z, atm, c_true, beta_surface=b0)
        for b0 in (5.0e-8, 1.5e-7, 5.0e-7)
    ]
    i36 = (at(CENTRES_M, 3500), at(CENTRES_M, 5000))
    for case in res["aerosol_curvature"]:
        print(f"  aerosol beta0={case['beta_surface']:.1e} (R(3km)={case['scattering_ratio_3km']:.2f}, "
              f"signal/p_mol {case['signal_over_pmol_slope_2_6km']:+.1f} %/km 2-6 km): "
              f"raw {case['raw_slope']:+.2f} %/km, fix {case['fix_slope']:+.2f} %/km")
        print(f"     free-fit intercept with TRUE b = 0: b_hat = "
              f"{case['fitted_intercept'][i36[0]]:+.3e} @3.5 km, "
              f"{case['fitted_intercept'][i36[1]]:+.3e} @5 km  ->  rel_pred = "
              f"{case['rel_pred_pct'][i36[0]]:+.2f} % @3.5 km, "
              f"{case['rel_pred_pct'][i36[1]]:+.2f} % @5 km "
              f"(mesuré Payerne : -6,3 % @3,66 km, -6,9 % @4,86 km)")

    res["intercept_vs_aerosol"] = intercept_vs_aerosol(z, atm, c_true)
    print("  is the fitted intercept instrumental? (truth: b = 0, only haze)")
    i1, i2 = at(CENTRES_M, 3660), at(CENTRES_M, 4860)
    for case in res["intercept_vs_aerosol"]:
        print(f"    beta0={case['beta_surface']:.1e}  R(3km)={case['scattering_ratio_3km']:.2f}  "
              f"signal/p_mol {case['signal_over_pmol_slope_2_6km']:+6.1f} %/km  ->  rel_pred "
              f"{case['rel_pred_pct'][i1]:+6.2f} % @3.66 km, {case['rel_pred_pct'][i2]:+6.2f} % "
              f"@4.86 km   (mesuré : -6,32 % / -6,94 %)")

    print("\n=== A-ter. decision: raw vs in-window fix vs wide-range intercept ===")
    res["decision"] = {}
    for tag, b_true, sigma, beta0 in (
            ("clean, b mesuré", B_KEPT_MEASURED, PAYERNE_CHM15K["sigma_night_clean"], 0.0),
            ("noisy, b mesuré", B_KEPT_MEASURED, PAYERNE_CHM15K["sigma_night_noisy"], 0.0),
            ("clean, b=0", 0.0, PAYERNE_CHM15K["sigma_night_clean"], 0.0),
            ("clean, b mesuré + brume légère R(3km)~1.2", B_KEPT_MEASURED,
             PAYERNE_CHM15K["sigma_night_clean"], 1.5e-7),
            ("clean, b mesuré + brume forte R(3km)~2", B_KEPT_MEASURED,
             PAYERNE_CHM15K["sigma_night_clean"], 5.0e-7)):
        block = decision_experiment(z, atm, c_true, b_true=b_true, sigma=sigma,
                                    n_seeds=n_seeds, beta_surface=beta0)
        res["decision"][tag] = block
        i3, i6 = at(CENTRES_M, 3000), at(CENTRES_M, 6000)
        print(f"  --- {tag}  (b_hat 3-15 km = {block['b_hat_wide_mean']:+.3e} "
              f"+- {block['b_hat_wide_std']:.1e}, b_hat 8-15 km = {block['b_hat_far_mean']:+.3e} "
              f"+- {block['b_hat_far_std']:.1e}, vrai b = {b_true:+.3e})")
        for key in ("raw", "fix", "wide", "far"):
            k = block[key]
            print(f"    {key:>4}: slope {k['slope']:+7.2f} %/km | @3 km biais "
                  f"{k['bias_pct'][i3]:+7.2f} disp {k['scatter_pct'][i3]:6.2f} RMSE "
                  f"{k['rmse_pct'][i3]:6.2f} | @6 km biais {k['bias_pct'][i6]:+7.2f} disp "
                  f"{k['scatter_pct'][i6]:6.2f} RMSE {k['rmse_pct'][i6]:6.2f}   [%]")

    print("\n=== B. molecular model error ===")
    res["molecular"] = {}
    for tag, night, cams in (("winter", "2025-01-15", "CAMS_Beta_202501.nc"),
                             ("summer", "2026-05-22", "CAMS_Beta_202605.nc")):
        block = molecular_model_scan(z, c_true, night, CAMS_DIR / cams)
        res["molecular"][tag] = block
        print(f"  --- {tag} {night} (CAMS p_sfc = {block['p_sfc_hpa']:.1f} hPa, "
              f"T(3 km) = {block['truth_T_3km']:.1f} K)")
        for name, v in block["variants"].items():
            print(f"    {name:<32} slope {v['slope_pct_per_km']:+7.3f} %/km   "
                  f"C@3km {v['dev_at_3km']:+7.3f} %  C@6km {v['dev_at_6km']:+7.3f} %   "
                  f"p_mol ratio {v['pmol_ratio_3km']:.4f}->{v['pmol_ratio_6km']:.4f}")

    print("\n=== figures ===")
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig_additive(res, z, atm, c_true, FIGDIR / "fwd_additive_residual_ladder.png")
    fig_failure_modes(res, FIGDIR / "fwd_additive_subtract_failure_modes.png")
    fig_decision(res, FIGDIR / "fwd_additive_fix_decision.png")
    fig_molecular(res, z, FIGDIR / "fwd_molecular_model_error.png")
    fig_summary(res, FIGDIR / "fwd_scan4_amplitude_summary.png")
    for name in ("fwd_additive_residual_ladder.png", "fwd_additive_subtract_failure_modes.png",
                 "fwd_additive_fix_decision.png",
                 "fwd_molecular_model_error.png", "fwd_scan4_amplitude_summary.png"):
        print(f"  -> {FIGDIR / name}")

    # drop the bulky per-gate ratio profiles before serialising (figures already used them)
    for block in res["molecular"].values():
        for variant in block["variants"].values():
            variant.pop("pmol_ratio_profile", None)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"\nnumbers -> {OUT_JSON}")


if __name__ == "__main__":
    main()
