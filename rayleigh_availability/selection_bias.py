# -*- coding: utf-8 -*-
"""MECHANISM SCAN 5 -- is the molecular-window SELECTOR a biased estimator on noisy data?

THE QUESTION
------------
In E-PROF v2 / v2.2 the molecular window is not fixed: it is CHOSEN, out of ~136 candidate
(centre, half-length) cells, as the one that maximises a composite quality score (R2-led) among
the cells that pass a set of gates. Two of those gates and one score term are built on a
statistic that is itself an EXTREMUM over the same 136 candidates:

    scattering_ratio(i,j) = ratio_med(i,j) / min over "clean" windows of ratio_med
    ratio_med(i,j)        = median over the window of signal / p_mol   ==  the C_L PROXY

Choosing an extremum over N noisy estimates is the textbook biased-estimator setup, and the code
itself documents it (molecular_methods.compute_window_grid, the note above `clean = ...`): the
minimum of N noisy ratio_med is biased LOW by an amount that GROWS with the noise. Every window
is then measured against that too-low reference, the gate `scattering_ratio <= 1.15` keeps only
the windows nearest the (downward-biased) minimum, and the score term `-w_ratio*|SR-1|` pulls in
the same direction. If that survives into the constant, then a noisy night returns a LOW C_L for a
purely instrumental reason -- and since v2.2 recovers exactly the noisy nights, it would explain
part of the recovered-vs-kept offset WITHOUT any aerosol.

This module measures that, end to end, on a PERFECTLY CLEAN simulated atmosphere where the true
answer is known, so anything that comes out is estimator bias by construction.

WHAT IS SIMULATED (forward model)
---------------------------------
Pure molecular air, no aerosol at all: beta_aer = 0, ext_aer = 0 everywhere. On the real Payerne
CHM15k L2 grid (512 gates, dz = 29.97 m, first gate 22.4775 m), 1064 nm, station 490 m ASL,

    ext_tot(z) = beta_mol(z) * 8*pi/3              [m^-1]
    T2(z)      = exp(-2 * INT_0^z ext_tot dz')     [-]
    rcs_i(z)   = C_true * beta_mol(z) * T2(z) + eps_i(z) * z^2

with C_true = 6.5e11 (the recorded Payerne CHM15k constant) and eps_i ~ N(0, sigma_prof) the
photon noise OF THE RANGE-NORMALISED SIGNAL of ONE binned (300 s x 30 m) profile. A night is
``n_prof`` such profiles; the night mean therefore carries

    sigma_mean = sigma_prof / sqrt(n_prof)         [signal units]

and ``sigma_mean`` is the single knob of the SNR ladder. Everything else is bit-identical between
the rungs -- same atmosphere, same C_true, same grid, same number of profiles.

THE SNR LADDER IS MEASURED, NOT ASSUMED
---------------------------------------
``sigma_mean`` was measured on 39 real Payerne CHM15k nights (scratchpad/measure_snr.py) as the
range first-difference scatter of the night-mean signal over 11-15 km (pure-noise altitudes),
after running the production front end (calibrate_rayleigh(..., fit_inputs_out=...)):

    v2.0-KEPT nights       median 1.61e-4   (p10 1.30e-4, p90 3.98e-4), n_prof median 102
    v2.2-RECOVERED nights  median 4.64e-4   (p10 3.58e-4, p90 5.86e-4), n_prof median  89

i.e. the nights v2.2 recovers are 2.9x NOISIER than the nights v2.0 keeps. The same run confirms
the sigma bug quantitatively: the pipeline's own sigma_signal / the measured one is 4.10 (kept)
and 4.33 (recovered), against sqrt(avg_time/dt_native) = sqrt(300/15) = 4.472.

THE THREE CONFIGURATIONS COMPARED
---------------------------------
  v2.0              molecular_method 'eprof_v2'   -- no noise input at all
  v2.2 sigma BUGGY  'eprof_v2.2' fed sigma_mean * sqrt(20)  (what the shipped code passes today)
  v2.2 sigma FIXED  'eprof_v2.2' fed sigma_mean             (the missing time-binning factor)

WHAT IS RUN (the SHIPPED chain, nothing reimplemented)
------------------------------------------------------
  molecular_methods.select_molecular_window   window search + cell flagging + gates + score
  rayleigh_fit._result_from_method_window     the same mapping calibrate_rayleigh uses
  calibration._compute_cl_for_perturbation    Klett + calculate_lidar_constant, once per
                                              (time subset x lidar ratio x altitude shift)
and the reported constant is the production best estimate: the MEDIAN over
6 time samples (full night + 4 random 70 % subsets, rng seed 42 as in production)
x 5 lidar ratios (52 +- 10, +- 20 sr) x 5 window shifts (0, +-100, +-200 m) = 150 values.
The production post-fit gates are applied too (flag -2 when no window passes or when
relative_error > threshold_quality = 15 %).

Run:  python rayleigh_availability/selection_bias.py [--seeds N] [--workers N] [--quick]
Out:  doc/reports/figs_altitude_audit/sel_bias_*.png  +  selection_bias_results.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from calibration.rayleigh.calibration import (            # noqa: E402
    ALT_SHIFTS_M, LR_DELTAS, N_TIME_SAMPLES, TIME_SUBSET_FRACTION,
    _compute_cl_for_perturbation,
)
from calibration.rayleigh.molecular_methods import select_molecular_window  # noqa: E402
from calibration.rayleigh.rayleigh_fit import _result_from_method_window    # noqa: E402
from rayleigh_availability.forward_model import (          # noqa: E402
    PAYERNE_CHM15K, default_grid, default_options, forward_signal, make_atmosphere,
)

FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"
OUTJSON = REPO / "rayleigh_availability" / "selection_bias_results.json"

# ---------------------------------------------------------------------------
# Measured SNR ladder (scratchpad/measure_snr.py, 39 real Payerne CHM15k nights).
# sigma_mean = 1-sigma photon noise of the NIGHT-MEAN range-normalised signal, signal units.
# ---------------------------------------------------------------------------
SIGMA_KEPT_MED = 1.61e-4        # v2.0-kept nights, median
SIGMA_REC_MED = 4.64e-4         # v2.2-recovered nights, median
SIGMA_LADDER = (1.20e-4, 1.61e-4, 2.40e-4, 3.20e-4, 4.64e-4, 6.40e-4,
                9.00e-4, 1.30e-3, 1.80e-3, 2.60e-3)
N_PROF = 100                    # binned 300 s profiles per night (real median 102 / 89)
TIME_BIN_FACTOR = np.sqrt(300.0 / 15.0)   # avg_time / dt_native for a CHM15k = the missing factor

# ``params`` are method_params overrides handed to select_molecular_window.
#   max_residual_pct = -1  -> the STRICT tier can never qualify a window, so the eprof_v2.2
#                             noise tier is the one actually exercised. Nothing else changes:
#                             tier 2 does not read max_residual_pct. This is the only way to
#                             characterise tier 2 on a clean atmosphere, because on clean air
#                             the strict tier almost always succeeds and tier 2 is never entered.
#   flag_nmad = 1e9        -> switches OFF the time-altitude contaminated-cell screen
#                             (flag_contaminated_cells), an UPPER-TAIL clip that can only lower
#                             the time-cleaned mean the fit slope is measured on.
CONFIGS = (
    # (key, label, molecular_method, sigma mode, params)
    ("v20",        "v2.0 (eprof_v2)",             "eprof_v2",   "none",  {}),
    ("v20_noflag", "v2.0 sans ecretage de cellules", "eprof_v2", "none", dict(flag_nmad=1e9)),
    ("v22_bug",    "v2.2, sigma tel que livre",   "eprof_v2.2", "buggy", {}),
    ("v22_fix",    "v2.2, sigma corrige",         "eprof_v2.2", "fixed", {}),
    ("t2_bug",     "tier 2 force, sigma livre",   "eprof_v2.2", "buggy", dict(max_residual_pct=-1.0)),
    ("t2_fix",     "tier 2 force, sigma corrige", "eprof_v2.2", "fixed", dict(max_residual_pct=-1.0)),
    # The recovered nights fit 1-2 km HIGHER than the kept ones. On clean air the selector never
    # goes up (see the results), so the height is imposed by the atmosphere, not chosen. These two
    # rungs answer the follow-up question: ONCE the search is confined to the altitudes the
    # recovered nights actually use (window starting above 3 km, centres 4-6 km), where the SNR is
    # 4-5x worse, how big is the selection bias there? Same clean atmosphere, same C_true.
    ("high_v20",   "v2.0, recherche >3 km",       "eprof_v2",   "none",
     dict(min_window_start_m=3000.0, range_start_m=4000.0)),
    ("high_v22",   "v2.2 sigma corrige, >3 km",   "eprof_v2.2", "fixed",
     dict(min_window_start_m=3000.0, range_start_m=4000.0)),
)


# ---------------------------------------------------------------------------
# One synthetic night through the shipped selector + constant
# ---------------------------------------------------------------------------
def _atmosphere_cache(_cache={}):
    """Grid + molecular atmosphere + noiseless rcs (identical for every rung; built once)."""
    if not _cache:
        z = default_grid()
        atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                              wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
        rcs_true = forward_signal(z, PAYERNE_CHM15K["C_true"], atm["beta_mol"])
        _cache.update(z=z, atm=atm, rcs_true=rcs_true)
    return _cache


def run_night(sigma_mean: float, seed: int, method: str, sigma_mode: str,
              n_prof: int = N_PROF, options=None, params: dict | None = None) -> dict:
    """Simulate one pure-molecular night at noise ``sigma_mean`` and run the shipped retrieval.

    Returns a dict with the retrieved constant, the selected window, the eligibility tier and
    the production flag. ``C_L`` is the median over the full production perturbation grid.
    """
    c = _atmosphere_cache()
    z, atm, rcs_true = c["z"], c["atm"], c["rcs_true"]
    options = options or default_options(method)
    C_true = PAYERNE_CHM15K["C_true"]

    rng = np.random.default_rng(seed)
    sigma_prof = float(sigma_mean) * np.sqrt(n_prof)          # noise of ONE binned profile
    z2 = z ** 2
    rcs_stack = rcs_true[None, :] + rng.normal(0.0, sigma_prof, size=(n_prof, z.size)) * z2[None, :]
    rcs_mean = rcs_stack.mean(axis=0)
    signal = rcs_mean / z2
    signal_stack = rcs_stack / z2[None, :]

    if sigma_mode == "none":
        sigma_signal = None
    elif sigma_mode == "buggy":
        sigma_signal = np.full(z.size, float(sigma_mean) * TIME_BIN_FACTOR)
    elif sigma_mode == "fixed":
        sigma_signal = np.full(z.size, float(sigma_mean))
    else:
        raise ValueError(sigma_mode)

    out = dict(seed=int(seed), sigma_mean=float(sigma_mean), method=method,
               sigma_mode=sigma_mode, ok=False, flag=-2, C_L=np.nan, C_L_nominal=np.nan,
               ratio_med=np.nan, slope=np.nan, r2=np.nan, rel_error=np.nan,
               start_m=np.nan, end_m=np.nan, centre_m=np.nan, half_m=np.nan,
               scattering_ratio=np.nan, tier="none", window_ok=False, n_clean_frac=np.nan)

    p = dict(params or {})
    kw = dict(range_start_m=p.pop("range_start_m", options.range_start_m),
              range_end_m=p.pop("range_end_m", options.range_end_m),
              increment_bins=options.fit_range_increment_bins,
              signal_stack=signal_stack, sigma_signal=sigma_signal, **p)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mw = select_molecular_window(method, signal, atm["p_mol"], z,
                                     options.half_length_options_m, **kw)
        # Which TIER produced this window? The selector's message does not say so in the default
        # (score) objective, and tier 2 is by construction only reached when the strict v2 gates
        # qualify nothing -- so the honest test is to ask the strict selector itself on the same
        # data. Deterministic, same noise, so this is the exact tier-1 outcome.
        strict_ok = mw.ok
        if method == "eprof_v2.2":
            strict_ok = select_molecular_window("eprof_v2", signal, atm["p_mol"], z,
                                                options.half_length_options_m, **kw).ok
    out["n_clean_frac"] = float(mw.n_clean_frac) if np.isfinite(mw.n_clean_frac) else np.nan
    if not mw.ok:
        out["message"] = mw.message
        return out
    out.update(window_ok=True,
               tier=("strict" if strict_ok else "fallback"),
               ratio_med=float(mw.cl), slope=float(mw.slope), r2=float(mw.r2),
               start_m=float(mw.start_m), end_m=float(mw.end_m),
               centre_m=float(mw.center_m), half_m=float(mw.half_m),
               scattering_ratio=float(mw.scattering_ratio), message=mw.message)

    fit_result = _result_from_method_window(mw, gated=True)
    out["rel_error"] = float(fit_result.relative_error)
    # Production post-fit gates (calibrate_rayleigh step 7).
    if not np.isfinite(fit_result.slope) or not np.isfinite(fit_result.relative_error):
        return out
    if fit_result.slope <= 0:
        out["flag"] = -7
        return out
    if abs(fit_result.intercept) >= fit_result.slope:
        out["flag"] = -8
        return out
    if fit_result.relative_error > options.threshold_quality:
        return out

    # Production perturbation grid: time subsets x lidar ratio x altitude shift.
    rng_sub = np.random.default_rng(42)                    # production uses this fixed seed
    n_keep = min(n_prof, max(3, int(n_prof * TIME_SUBSET_FRACTION)))
    means = [rcs_mean]
    for _ in range(N_TIME_SAMPLES):
        idx = np.sort(rng_sub.choice(n_prof, size=n_keep, replace=False))
        means.append(rcs_stack[idx].mean(axis=0))

    values, nominal = [], np.nan
    for t_idx, rcs_t in enumerate(means):
        for d in LR_DELTAS:
            lr = options.lidar_ratio_aerosol + d
            if lr <= 0:
                continue
            for shift in ALT_SHIFTS_M:
                pr = _compute_cl_for_perturbation(
                    rcs_mean=rcs_t, range_alc=z, beta_mol=atm["beta_mol"],
                    fit_result=fit_result, lidar_ratio_aerosol=lr, altitude_shift_m=shift,
                    subtract_background=options.subtract_background,
                    consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
                    sign_error_v10=getattr(options, "sign_error_v10", False))
                if pr is None:
                    continue
                values.append(pr.lidar_constant)
                if t_idx == 0 and lr == options.lidar_ratio_aerosol and shift == 0:
                    nominal = float(pr.lidar_constant)
    if not values:
        out["flag"] = -5
        return out
    out.update(ok=True, flag=1, C_L=float(np.median(values)), C_L_nominal=nominal,
               C_true=C_true)
    return out


# ---------------------------------------------------------------------------
# Batch driver
# ---------------------------------------------------------------------------
def _job(args):
    sigma, seed, key, method, mode, params = args
    r = run_night(sigma, seed, method, mode, params=params)
    r["config"] = key
    r.pop("message", None)
    return r


def run_scan(sigmas=SIGMA_LADDER, n_seeds: int = 400, workers: int = 6,
             configs=CONFIGS) -> list:
    jobs = [(s, seed, key, method, mode, params)
            for (key, _lab, method, mode, params) in configs
            for s in sigmas
            for seed in range(n_seeds)]
    print(f"selection-bias scan: {len(jobs)} synthetic nights "
          f"({len(configs)} configs x {len(sigmas)} SNR rungs x {n_seeds} seeds), "
          f"{workers} workers", flush=True)
    res = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_job, j) for j in jobs]
        for k, f in enumerate(as_completed(futs)):
            res.append(f.result())
            if (k + 1) % 500 == 0:
                print(f"   {k + 1}/{len(jobs)}", flush=True)
    return res


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def summarise(res: list, sigmas=SIGMA_LADDER, configs=CONFIGS) -> dict:
    C_true = PAYERNE_CHM15K["C_true"]
    out = {}
    for key, label, method, mode, _params in configs:
        rows = []
        for s in sigmas:
            sub = [r for r in res if r["config"] == key and abs(r["sigma_mean"] - s) < 1e-12]
            n = len(sub)
            ok = [r for r in sub if r["ok"]]
            cl = np.array([r["C_L"] for r in ok], float)
            dev = (cl / C_true - 1.0) * 100.0 if cl.size else np.array([])
            start = np.array([r["start_m"] for r in ok], float)
            end = np.array([r["end_m"] for r in ok], float)
            cen = np.array([r["centre_m"] for r in ok], float)
            fb = sum(1 for r in sub if r["window_ok"] and r["tier"] == "fallback")
            ncf = np.array([r["n_clean_frac"] for r in ok], float)
            rows.append(dict(
                sigma_mean=float(s), n=n, n_ok=len(ok),
                success_pct=100.0 * len(ok) / n if n else np.nan,
                fallback_pct=100.0 * fb / n if n else np.nan,
                clean_frac_median=float(np.nanmedian(ncf)) if ncf.size else np.nan,
                bias_mean_pct=float(np.mean(dev)) if dev.size else np.nan,
                bias_median_pct=float(np.median(dev)) if dev.size else np.nan,
                bias_sem_pct=float(np.std(dev, ddof=1) / np.sqrt(dev.size)) if dev.size > 1 else np.nan,
                scatter_pct=float(np.std(dev, ddof=1)) if dev.size > 1 else np.nan,
                start_median_m=float(np.median(start)) if start.size else np.nan,
                start_p10_m=float(np.percentile(start, 10)) if start.size else np.nan,
                start_p90_m=float(np.percentile(start, 90)) if start.size else np.nan,
                centre_median_m=float(np.median(cen)) if cen.size else np.nan,
                top_median_m=float(np.median(end)) if end.size else np.nan,
                width_median_m=float(np.median(end - start)) if end.size else np.nan,
            ))
        out[key] = dict(label=label, method=method, sigma_mode=mode, rows=rows)
    return out


def _slope_per_decade(x, y):
    """Least-squares slope of y versus log10(x); returns NaN with < 3 finite points."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    g = np.isfinite(x) & np.isfinite(y) & (x > 0)
    if g.sum() < 3:
        return np.nan
    return float(np.polyfit(np.log10(x[g]), y[g], 1)[0])


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
COL = {"v20": "#444444", "v20_noflag": "#9467bd", "v22_bug": "#d62728", "v22_fix": "#1f77b4",
       "t2_bug": "#ff7f0e", "t2_fix": "#17becf", "high_v20": "#8c564b", "high_v22": "#e377c2"}
MAIN = ("v20", "v22_bug", "v22_fix")     # the three production-relevant configurations


def figure_bias(summary: dict, path: Path) -> None:
    """Panel A: C_L bias vs noise. Panel B: selected window (altitude on Y) vs noise.
    Panel C: success rate. All three share the SNR axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17.0, 6.2))
    ax = axes[0]
    for key in MAIN:
        r = summary[key]["rows"]
        s = np.array([q["sigma_mean"] for q in r])
        b = np.array([q["bias_mean_pct"] for q in r])
        e = np.array([q["bias_sem_pct"] for q in r])
        ax.errorbar(s, b, yerr=e, fmt="o-", color=COL[key], ms=5, lw=1.6, capsize=3,
                    label=summary[key]["label"])
    ax.axhline(0, color="0.4", lw=1.2, ls="--")
    ax.axvspan(1.30e-4, 3.98e-4, color="0.75", alpha=0.25)
    ax.axvspan(3.58e-4, 5.86e-4, color="#1f77b4", alpha=0.12)
    ax.axvline(SIGMA_KEPT_MED, color="0.35", lw=1.0, ls=":")
    ax.axvline(SIGMA_REC_MED, color="#1f77b4", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"bruit de la moyenne de nuit  $\sigma_{mean}$  [unites signal]")
    ax.set_ylabel(r"$C_L / C_{true} - 1$   [%]   (moyenne $\pm$ s.e.m.)")
    ax.set_title("A. Biais du constant sur atmosphere PUREMENT moleculaire\n"
                 "(aucun aerosol : tout ecart = biais d'estimateur)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    # Zoom on the band the real instruments actually occupy -- the whole decision is made there,
    # and on the full axis it is a flat line pinned to zero.
    axin = ax.inset_axes((0.50, 0.55, 0.47, 0.42))
    for key in MAIN:
        r = [q for q in summary[key]["rows"] if q["sigma_mean"] <= 6.5e-4]
        s = np.array([q["sigma_mean"] for q in r])
        axin.errorbar(s, [q["bias_mean_pct"] for q in r], yerr=[q["bias_sem_pct"] for q in r],
                      fmt="o-", color=COL[key], ms=3.5, lw=1.3, capsize=2)
    axin.axhline(0, color="0.4", lw=1.0, ls="--")
    axin.axvspan(1.30e-4, 3.98e-4, color="0.75", alpha=0.25)
    axin.axvspan(3.58e-4, 5.86e-4, color="#1f77b4", alpha=0.12)
    axin.set_xscale("log")
    axin.tick_params(labelsize=7)
    axin.set_title("zoom : bruit reel mesure", fontsize=7.5)
    axin.grid(alpha=0.3)

    ax = axes[1]
    for key in MAIN:
        r = summary[key]["rows"]
        s = np.array([q["sigma_mean"] for q in r])
        lo = np.array([q["start_p10_m"] for q in r]) / 1000.0
        hi = np.array([q["start_p90_m"] for q in r]) / 1000.0
        md = np.array([q["start_median_m"] for q in r]) / 1000.0
        tp = np.array([q["top_median_m"] for q in r]) / 1000.0
        ax.fill_betweenx(md, s * 0 + s.min(), s * 0 + s.max(), alpha=0)   # keep Y = altitude
        ax.plot(s, md, "o-", color=COL[key], ms=5, lw=1.8, label=f"{summary[key]['label']} - bas")
        ax.plot(s, tp, "s--", color=COL[key], ms=4, lw=1.2, alpha=0.65,
                label=f"{summary[key]['label']} - haut")
        ax.fill_between(s, lo, hi, color=COL[key], alpha=0.12)
    ax.set_xscale("log")
    ax.set_xlabel(r"bruit de la moyenne de nuit  $\sigma_{mean}$  [unites signal]")
    ax.set_ylabel("altitude de la fenetre moleculaire  [km AGL]")
    ax.set_title("B. Ou le selecteur place la fenetre\n(trait plein = bas, tirets = haut, "
                 "bande = p10-p90 du bas)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="center left", ncol=1)
    ax.text(0.97, 0.97, "le bas ne bouge JAMAIS (1738 m, la\nfenetre la plus basse autorisee) ;\n"
                        "seul le haut descend = fenetre plus etroite",
            transform=ax.transAxes, fontsize=7.5, va="top", ha="right",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    ax = axes[2]
    for key in MAIN:
        r = summary[key]["rows"]
        s = np.array([q["sigma_mean"] for q in r])
        ax.plot(s, [q["success_pct"] for q in r], "o-", color=COL[key], ms=5, lw=1.8,
                label=summary[key]["label"])
        if key.startswith("v22"):
            ax.plot(s, [q["fallback_pct"] for q in r], "^:", color=COL[key], ms=4, lw=1.0,
                    alpha=0.7, label=f"{summary[key]['label']} - % via tier 2")
    ax.set_xscale("log")
    ax.set_ylim(-3, 103)
    ax.set_xlabel(r"bruit de la moyenne de nuit  $\sigma_{mean}$  [unites signal]")
    ax.set_ylabel("nuits calibrees  [%]")
    ax.set_title("C. Disponibilite sur ciel parfaitement clair\n"
                 "(une nuit rejetee ici est un faux negatif pur)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")

    fig.suptitle("Scan 5 - biais de SELECTION de la fenetre moleculaire sur atmosphere purement "
                 f"moleculaire ($C_{{true}}$ = {PAYERNE_CHM15K['C_true']:.2e}, grille Payerne "
                 f"CHM15k, {N_PROF} profils/nuit, 400 tirages/point)\n"
                 "bande grise = bruit des nuits GARDEES par v2.0, bande bleue = bruit des nuits "
                 "RECUPEREES par v2.2 (p10-p90 mesures sur 39 nuits reelles)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    import matplotlib.pyplot as _plt
    _plt.close(fig)


def figure_mechanism(res: list, summary: dict, path: Path) -> None:
    """Why the bias exists: the scattering-ratio reference is a MIN over 136 noisy windows.

    Left: C_L bias versus the selected window's start altitude (altitude on Y), at the
    recovered-night noise level. Right: the min-of-N reference deficit versus noise.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from calibration.rayleigh.molecular_methods import compute_window_grid

    C_true = PAYERNE_CHM15K["C_true"]
    c = _atmosphere_cache()
    z, atm, rcs_true = c["z"], c["atm"], c["rcs_true"]
    options = default_options("eprof_v2")

    fig, axes = plt.subplots(1, 3, figsize=(17.0, 6.0))

    # --- panel 1: bias vs the selected window TOP (the only edge that moves) -----------
    # The window BOTTOM is pinned at 1738 m on clean air at every noise level, so the height
    # degree of freedom the selector actually exercises is the top (= the window width).
    ax = axes[0]
    for s, col, lab in ((SIGMA_KEPT_MED, "0.35", "bruit nuits gardees"),
                        (SIGMA_REC_MED, "#1f77b4", "bruit nuits recuperees"),
                        (9.0e-4, "#d62728", r"$\sigma$ = 9.0e-4")):
        sub = [r for r in res
               if r["config"] == "v22_fix" and r["ok"] and abs(r["sigma_mean"] - s) < 1e-12]
        if not sub:
            continue
        dev = np.array([r["C_L"] / C_true - 1.0 for r in sub]) * 100.0
        top = np.array([r["end_m"] for r in sub]) / 1000.0
        jit = np.random.default_rng(7).normal(0, 0.035, size=top.size)
        ax.plot(dev, top + jit, "o", color=col, ms=2.6, alpha=0.30)
        tops = np.unique(np.round(top, 2))                     # window tops are 10 m multiples
        mid = [t for t in tops if (np.abs(top - t) < 0.02).sum() >= 8]
        med = [float(np.median(dev[np.abs(top - t) < 0.02])) for t in mid]
        ax.plot(med, mid, "o-", color=col, lw=2.2, ms=6,
                label=f"{lab} ($\\sigma$={s:.2e}, n={len(sub)})")
    ax.axvline(0, color="0.3", lw=1.3, ls="--")
    ax.set_xlim(-6, 6)
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("HAUT de la fenetre selectionnee  [km AGL]")
    ax.set_title("A. Le seul bord que le selecteur deplace, c'est le haut\n"
                 "(le bas reste a 1738 m ; points = tirages, trait = mediane)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")

    # --- panel 2: the min-of-N reference deficit ---------------------------------------
    ax = axes[1]
    sig_scan = np.array([1.2e-4, 1.61e-4, 2.4e-4, 3.2e-4, 4.64e-4, 6.4e-4, 9.0e-4])
    n_rep = 60
    deficits = np.full((sig_scan.size, n_rep), np.nan)
    spread = np.full((sig_scan.size, n_rep), np.nan)
    for i, s in enumerate(sig_scan):
        for k in range(n_rep):
            rng = np.random.default_rng(90000 + 977 * i + k)
            sp = s * np.sqrt(N_PROF)
            stack = rcs_true[None, :] + rng.normal(0, sp, size=(N_PROF, z.size)) * (z ** 2)[None, :]
            sig = stack.mean(axis=0) / z ** 2
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                g = compute_window_grid(sig, atm["p_mol"], z, options.half_length_options_m,
                                        range_start_m=options.range_start_m,
                                        range_end_m=options.range_end_m,
                                        increment_bins=options.fit_range_increment_bins)
            clean = np.isfinite(g.ratio_med) & (g.ratio_med > 0) & (g.slope > 0) & (g.r2 >= 0.5)
            if not np.any(clean):
                continue
            c_min = float(np.min(g.ratio_med[clean]))
            deficits[i, k] = (c_min / C_true - 1.0) * 100.0
            spread[i, k] = (np.nanmedian(g.ratio_med[clean]) / C_true - 1.0) * 100.0
    ax.errorbar(sig_scan, np.nanmean(deficits, axis=1),
                yerr=np.nanstd(deficits, axis=1) / np.sqrt(n_rep),
                fmt="o-", color="#8c564b", ms=5, lw=1.8, capsize=3,
                label=r"reference $c_{min}$ = min sur les fenetres")
    ax.errorbar(sig_scan, np.nanmean(spread, axis=1),
                yerr=np.nanstd(spread, axis=1) / np.sqrt(n_rep),
                fmt="s--", color="#2ca02c", ms=5, lw=1.4, capsize=3,
                label=r"mediane des $ratio\_med$ (non biaisee)")
    ax.axhline(0, color="0.3", lw=1.2, ls="--")
    ax.axvline(SIGMA_KEPT_MED, color="0.35", lw=1.0, ls=":")
    ax.axvline(SIGMA_REC_MED, color="#1f77b4", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma_{mean}$  [unites signal]")
    ax.set_ylabel(r"ecart a $C_{true}$   [%]")
    ax.set_title("B. Le denominateur du rapport de diffusion est\n"
                 "un MINIMUM sur ~136 fenetres bruitees", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    ax.text(0.97, 0.05,
            "le biais min-de-N est REEL et croit avec le bruit,\n"
            "mais il n'entre que dans une PORTE (<=1,15) et\n"
            "une penalite de score, jamais dans le constant :\n"
            "il ne se propage donc pas a $C_L$",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.92))

    # --- panel 3: decomposition -- which piece of the selector carries the bias? --------
    ax = axes[2]
    for key in ("v20", "v20_noflag", "t2_bug", "t2_fix", "high_v20", "high_v22"):
        if key not in summary:
            continue
        r = summary[key]["rows"]
        s = np.array([q["sigma_mean"] for q in r])
        b = np.array([q["bias_mean_pct"] for q in r])
        e = np.array([q["bias_sem_pct"] for q in r])
        ax.errorbar(s, b, yerr=e, fmt="o-", color=COL[key], ms=4, lw=1.5, capsize=2,
                    label=summary[key]["label"])
    ax.axhline(0, color="0.3", lw=1.2, ls="--")
    ax.axvline(SIGMA_KEPT_MED, color="0.35", lw=1.0, ls=":")
    ax.axvline(SIGMA_REC_MED, color="#1f77b4", lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma_{mean}$  [unites signal]")
    ax.set_ylabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_title("C. Decomposition : ecretage de cellules,\ntier 2 force, bug de sigma", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")

    fig.suptitle("Scan 5 - d'ou vient le biais : selection d'un extremum sur ~136 fenetres bruitees",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


OBS_DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
OBS_RAW = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/calib_raw")
OBS_VALID = (1.0, 0.5)
OBS_SITES = (
    ("Payerne CHM15k", "PAYERNE_CHM15k_A", 490.0),
    ("Payerne CL61", "__payerne_cl61__", 490.0),
    ("Lindenberg CHM15k", "LINDENBERG_CHM15k_0", 125.0),
    ("Aoste CHM15k", "SAINT-CHRISTOPHE_A_CHM15k_0", 570.0),
    ("Aoste CL61", "SAINT-CHRISTOPHE_A_CL61_B", 570.0),
    ("SIRTA CHM15k", "PALAISEAU_CHM15k_B", 156.0),
)


def _observed_offsets() -> list:
    """(label, dC %, kept window centre km AGL, recovered centre km AGL, n_kept, n_rec) per site."""
    rows = []
    for label, key, alt in OBS_SITES:
        try:
            if key == "__payerne_cl61__":
                base = json.loads((OBS_RAW / "payerne_C_v2.0.json").read_text())
                cand = json.loads((OBS_RAW / "payerne_C_v2.2.json").read_text())
            else:
                base = json.loads((OBS_DATA / "baselines" / f"base_eprof_v2_{key}.json").read_text())
                cand = json.loads((OBS_DATA / "candidates" / f"cand_N2.5_{key}.json").read_text())
        except FileNotFoundError:
            continue
        bv = {d for d, v in base.items() if v[0] in OBS_VALID and v[1]}
        kept = [v for d, v in cand.items() if v[0] in OBS_VALID and v[1] and d in bv]
        rec = [v for d, v in cand.items() if v[0] in OBS_VALID and v[1] and d not in bv]
        if len(kept) < 5 or len(rec) < 5:
            continue
        med = lambda g, i: float(np.median([r[i] for r in g if r[i] is not None]))  # noqa: E731
        ck, cr = med(kept, 1), med(rec, 1)
        zk = 0.5 * (med(kept, 3) + med(kept, 4)) - alt
        zr = 0.5 * (med(rec, 3) + med(rec, 4)) - alt
        rows.append((label, 100.0 * (cr / ck - 1.0), zk / 1000.0, zr / 1000.0, len(kept), len(rec)))
    return rows


def figure_vs_observed(summary: dict, path: Path) -> None:
    """The verdict figure: what the SELECTOR can produce vs what the network actually shows.

    Both panels have the molecular-window centre on Y and the constant offset on X, on the SAME
    x-scale, so the shortfall is read directly off the page. Left = the six observed streams
    (arrow from the v2.0-kept population to the v2.2-recovered one). Right = the simulation on a
    pure molecular atmosphere, moving along the MEASURED noise ladder from the kept-night noise
    (1.61e-4) to the recovered-night noise (4.64e-4).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    obs = _observed_offsets()
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), sharex=True)

    ax = axes[0]
    for k, (label, dc, zk, zr, nk, nr) in enumerate(obs):
        col = f"C{k}"
        ax.annotate("", xy=(dc, zr), xytext=(0.0, zk),
                    arrowprops=dict(arrowstyle="-|>", lw=2.0, color=col, alpha=0.85))
        ax.plot([0.0], [zk], "o", color=col, ms=7)
        ax.plot([dc], [zr], "s", color=col, ms=7,
                label=f"{label}  ({nk} gardees / {nr} recuperees)")
        ax.text(dc, zr + 0.10, f"{dc:+.1f} %", color=col, fontsize=8, ha="center")
    ax.axvline(0, color="0.3", lw=1.3, ls="--")
    ax.set_xlabel(r"$C_L$ des nuits recuperees / $C_L$ des nuits gardees $-$ 1   [%]")
    ax.set_ylabel("centre de la fenetre moleculaire  [km AGL]")
    ax.set_title("OBSERVE - v2.0 gardees $\\rightarrow$ v2.2 recuperees\n"
                 "(medianes par flux, archive 2025-2026)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")

    ax = axes[1]
    for key in MAIN:
        rows = {r["sigma_mean"]: r for r in summary[key]["rows"]}
        a, b = rows.get(SIGMA_KEPT_MED), rows.get(SIGMA_REC_MED)
        if a is None or b is None:
            continue
        dc = b["bias_mean_pct"] - a["bias_mean_pct"]
        sem = float(np.hypot(a["bias_sem_pct"], b["bias_sem_pct"]))
        za, zb = a["centre_median_m"] / 1000.0, b["centre_median_m"] / 1000.0
        ax.annotate("", xy=(dc, zb), xytext=(0.0, za),
                    arrowprops=dict(arrowstyle="-|>", lw=2.0, color=COL[key]))
        ax.errorbar([dc], [zb], xerr=[sem], fmt="s", color=COL[key], ms=7, capsize=3,
                    label=f"{summary[key]['label']}: {dc:+.3f} $\\pm$ {sem:.3f} %")
        ax.plot([0.0], [za], "o", color=COL[key], ms=7)
    ax.axvline(0, color="0.3", lw=1.3, ls="--")
    ax.set_xlabel(r"$C_L(\sigma_{recup}) / C_L(\sigma_{garde}) - 1$   [%]")
    ax.set_ylabel("centre de la fenetre moleculaire  [km AGL]")
    ax.set_title("SIMULE - meme atmosphere PUREMENT moleculaire,\n"
                 f"bruit {SIGMA_KEPT_MED:.2e} $\\rightarrow$ {SIGMA_REC_MED:.2e} "
                 "(l'echelle x est celle de gauche)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    # Same experiment, but with the search confined to the altitudes the recovered nights use.
    if "high_v22" in summary:
        rows = {r["sigma_mean"]: r for r in summary["high_v22"]["rows"]}
        ref = {r["sigma_mean"]: r for r in summary["v22_fix"]["rows"]}[SIGMA_KEPT_MED]
        b = rows[SIGMA_REC_MED]
        dc = b["bias_mean_pct"] - ref["bias_mean_pct"]
        sem = float(np.hypot(ref["bias_sem_pct"], b["bias_sem_pct"]))
        ax.annotate("", xy=(dc, b["centre_median_m"] / 1000.0),
                    xytext=(0.0, ref["centre_median_m"] / 1000.0),
                    arrowprops=dict(arrowstyle="-|>", lw=2.0, color=COL["high_v22"], ls="--"))
        ax.errorbar([dc], [b["centre_median_m"] / 1000.0], xerr=[sem], fmt="D",
                    color=COL["high_v22"], ms=7, capsize=3,
                    label=(f"v2.2 forcee a chercher >3 km (tier 2 a 80 %) : "
                           f"{dc:+.3f} $\\pm$ {sem:.3f} %"))
    if obs:
        worst = max(abs(o[1]) for o in obs)
        keys = list(MAIN) + (["high_v22"] if "high_v22" in summary else [])
        best_sim = max(
            abs({r["sigma_mean"]: r for r in summary[k]["rows"]}[SIGMA_REC_MED]["bias_mean_pct"]
                - {r["sigma_mean"]: r for r in summary["v20"]["rows"]}[SIGMA_KEPT_MED]["bias_mean_pct"])
            for k in keys)
        ax.text(0.97, 0.97,
                f"amplitude observee max : {worst:.1f} %\n"
                f"amplitude simulee max  : {best_sim:.3f} %\n"
                f"facteur manquant       : {worst / max(best_sim, 1e-9):.0f}x\n"
                "et le signe simule est TOUJOURS negatif :\n"
                "il ne peut pas produire les +23 a +32 %",
                transform=ax.transAxes, fontsize=8.5, va="top", ha="right",
                bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.95))
    ax.legend(fontsize=8, loc="lower left")
    ax.set_ylim(axes[0].get_ylim())

    fig.suptitle("Scan 5 - le biais de SELECTION peut-il produire le decalage observe ? "
                 "Non : il en explique moins de 1 %.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def print_table(summary: dict) -> None:
    for key in summary:
        d = summary[key]
        print(f"\n=== {d['label']}  (method={d['method']}, sigma={d['sigma_mode']}) ===")
        print(f"{'sigma_mean':>11} {'ok%':>6} {'tier2%':>7} {'biais moy %':>12} {'s.e.m.':>7} "
              f"{'median %':>9} {'sigma_nuit%':>11} {'bas m':>7} {'haut m':>7} {'cell%':>6}")
        for r in d["rows"]:
            print(f"{r['sigma_mean']:11.2e} {r['success_pct']:6.1f} {r['fallback_pct']:7.1f} "
                  f"{r['bias_mean_pct']:12.3f} {r['bias_sem_pct']:7.3f} "
                  f"{r['bias_median_pct']:9.3f} {r['scatter_pct']:11.2f} "
                  f"{r['start_median_m']:7.0f} {r['top_median_m']:7.0f} "
                  f"{100 * (1 - r['clean_frac_median']):6.2f}")
        s = [r["sigma_mean"] for r in d["rows"]]
        print(f"  pente du biais : {_slope_per_decade(s, [r['bias_mean_pct'] for r in d['rows']]):+.2f} "
              f"%/decade de bruit ; pente du bas de fenetre : "
              f"{_slope_per_decade(s, [r['start_median_m'] for r in d['rows']]):+.0f} m/decade")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=400)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--verdict-only", action="store_true",
                    help="redraw sel_bias_vs_observed.png from the stored results, no simulation")
    a = ap.parse_args()
    if a.verdict_only:
        summary = json.loads(OUTJSON.read_text())["summary"]
        print_table(summary)
        figure_vs_observed(summary, FIGDIR / "sel_bias_vs_observed.png")
        print(f"\nfigure -> {FIGDIR / 'sel_bias_vs_observed.png'}")
        return
    n_seeds = 30 if a.quick else a.seeds
    res = run_scan(n_seeds=n_seeds, workers=a.workers)
    summary = summarise(res)
    print_table(summary)
    OUTJSON.write_text(json.dumps(dict(summary=summary, n_seeds=n_seeds,
                                       sigma_ladder=list(SIGMA_LADDER), n_prof=N_PROF,
                                       C_true=PAYERNE_CHM15K["C_true"]), indent=1))
    # Per-night rows, so any figure can be redrawn without re-simulating.
    OUTJSON.with_name("selection_bias_raw.json").write_text(json.dumps(
        [{k: r[k] for k in ("config", "sigma_mean", "seed", "ok", "C_L", "start_m", "end_m",
                            "ratio_med", "tier", "n_clean_frac")} for r in res]))
    figure_bias(summary, FIGDIR / "sel_bias_ladder.png")
    figure_mechanism(res, summary, FIGDIR / "sel_bias_mechanism.png")
    figure_vs_observed(summary, FIGDIR / "sel_bias_vs_observed.png")
    print(f"\nfigures -> {FIGDIR}")


if __name__ == "__main__":
    main()
