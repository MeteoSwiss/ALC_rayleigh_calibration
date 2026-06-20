"""
Pluggable molecular-window detection for Rayleigh calibration.

Six interchangeable strategies choose the aerosol-free molecular reference window
from one time-collapsed, range-normalized profile (signal = RCS / r^2). They share
a single grid search (:func:`compute_window_grid`) and differ only in the eligibility
mask and the selection rule.

(A seventh "calipso"/CALIOP-style strategy — normalize to the *highest* clean layer —
was removed: it is physically appropriate only for a *down-looking satellite*, which
always has a pure-Rayleigh stratosphere beneath the aerosol, not for a ground-up
ceilometer whose high-altitude signal is photon noise, so it chased noise and was the
least precise strategy. See doc/reports for the drop rationale.)

Method keys are the E-PROF calibration versions (display labels in parentheses):

  "eprof_v1.1"   E-PROF v1.1 (sign cor) — legacy main-branch window rule: center =
              min Σ|intercept| (free intercept), half = max R²; NO R² floor. Kept for
              comparison; known degenerate (selects high-altitude noise where signal→0 ⇒
              intercept→0 and R²→0). E-PROF v1.0 (sign error) is THIS window run through the
              full pipeline with the historical Klett sign error (config.sign_error_v10).
  "eprof_v1.2"   E-PROF v1.2 (improved) — production fix: eligible windows (start above
              aerosol, R²≥floor, slope>0, |b|<a, slope≈median ratio) then MAX R².
  "eprof_v0.25"  E-PROF v0.25 (MATLAB) — Auto_Calib_25/Rayleigh/rayleigh_fit.m (Hervo &
              Poltera 2014): intercept FORCED to 0, center = min Σ RMSE, half = max R²;
              reject R²<min_r2; centers capped at 5000 m.
  "earlinet"  EARLINET / Single Calculus Chain (Mattis, D'Amico, Baars et al. 2016;
              Freudenthaler et al. 2018): Rayleigh-shape residual gate + SNR gate +
              scattering ratio ≤ 1.1, pick the LOWEST qualifying window (best SNR).
  "eprof_v2"  E-PROF v2 (optimal) — best-of: all physical gates + TEMPORAL-variability
              aerosol rejection (molecular scattering is steady in time — it only tracks
              T/p — whereas aerosol advects and fluctuates; the ONLY method that uses the
              full profile time series, not just the night mean) + a composite quality
              score (R², shape residual, scattering-ratio→1, SNR). It additionally FLAGS
              and EXCLUDES contaminated time-altitude cells (aerosol/cloud present only
              part of the night) and fits on the time-cleaned mean profile, so it uses the
              clean part of an otherwise-contaminated night (see flag_contaminated_cells).
  "bellini"   ALICENET / Bellini et al. (2024, AMT 17, 6119). Windows in 3-7 km, widths
              600-3000 m; rejects windows whose fit residuals are autocorrelated
              (Breusch-Godfrey test = undetected aerosol), selects max
              M_Ray = (adjR² + (1-|b|))/std(b), with border-residual-sign and an E_CL
              relative-uncertainty gate (see _select_bellini).

Legacy aliases (main/improved/matlab/optimal/eprof_v10) are still accepted on input.

References: Wiegner & Geiß (2012, AMT 5, 1953); Mattis, D'Amico, Baars et al. (2016,
AMT 9, 3009); Freudenthaler et al. (2018, amt-2017-395); Baars et al. (2016, ACP 16,
5111); Bellini et al. (2024, AMT 17, 6119, ALICENET);
EarthCARE ATLID is HSRL (direct molecular channel — no window search needed; see report).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import linregress


# Live selectable methods, keyed by E-PROF calibration version.
METHODS = ("eprof_v1.1", "eprof_v1.2", "eprof_v0.25", "earlinet", "eprof_v2", "bellini")

# Human-readable display labels (figures, reports, tables).
METHOD_LABELS: Dict[str, str] = {
    "eprof_v1.0": "E-PROF v1.0 (sign error)",
    "eprof_v1.1": "E-PROF v1.1 (sign cor)",
    "eprof_v1.2": "E-PROF v1.2 (improved)",
    "eprof_v0.25": "E-PROF v0.25 (MATLAB)",
    "eprof_v2": "E-PROF v2",
    "earlinet": "EARLINET/SCC",
    "bellini": "Bellini/ALICENET",
}

# Back-compat: old method names accepted on input and mapped to the E-PROF version keys.
ALIASES: Dict[str, str] = {
    "main": "eprof_v1.1", "improved": "eprof_v1.2", "matlab": "eprof_v0.25",
    "optimal": "eprof_v2", "eprof_v10": "eprof_v1.0",
}


def resolve_method(method: str) -> str:
    """Map a method name (E-PROF key or legacy alias) to its canonical E-PROF key."""
    m = method.lower()
    return ALIASES.get(m, m)


# Per-method default parameters. Overridable via select_molecular_window(**params).
DEFAULT_PARAMS: Dict[str, Dict[str, Any]] = {
    "eprof_v1.1":  dict(),
    "eprof_v1.2":  dict(min_window_start_m=2000.0, min_r2=0.5, max_rel_error=50.0),
    "eprof_v0.25": dict(min_r2=0.5, range_end_cap_m=5000.0),
    "earlinet": dict(min_window_start_m=2000.0, max_residual_pct=10.0,
                     max_ratio_std=0.30, max_scattering_ratio=1.1, max_rel_error=15.0),
    # Gates optimized 2026-06 over 24 instruments x both levels x all clear nights (see
    # validation/run_v2_sweep.py + doc/reports/v2_optimization_report.md, config "C8"). The
    # scattering-ratio gate was the dominant cause of clear-night rejections (esp. CL61), so it
    # was relaxed 1.10 -> 1.15; start/r2/residual/ratio_std/temporal were eased moderately. This
    # raises valid-night yield on every instrument type and level while holding sigma_SD ~flat.
    # Previous baseline: start=2000, r2=0.5, residual=12, scattering=1.1, ratio_std=0.30, tcv=0.5.
    "eprof_v2": dict(min_window_start_m=1500.0, min_r2=0.40, max_residual_pct=16.0,
                     max_scattering_ratio=1.15, max_ratio_std=0.40, max_temporal_cv=0.8,
                     max_rel_error=15.0, w_ratio=0.25, w_resid=0.20, w_snr=0.10,
                     w_npts=0.10, w_tvar=0.35, w_rel=0.20,
                     flag_nmad=4.0, flag_min_excess=0.25),
    "bellini":  dict(min_window_start_m=3000.0, max_window_end_m=7000.0, min_width_m=600.0,
                     max_width_m=3000.0, max_rel_error=15.0, max_ecl=0.40, border_m=200.0),
}


@dataclass
class WindowGrid:
    """All per-window statistics of the molecular grid search (n_centers × n_lengths)."""
    center_m: NDArray[np.float64]        # (n_centers,) center range, m AGL
    half_m: NDArray[np.float64]          # (n_lengths,) half-length, m
    start_m: NDArray[np.float64]         # (n_centers, n_lengths) window start, m AGL
    slope: NDArray[np.float64]           # free-fit slope a (signal = a*p_mol + b)
    intercept: NDArray[np.float64]       # free-fit intercept b
    r2: NDArray[np.float64]              # free-fit R²
    std_err: NDArray[np.float64]         # free-fit slope std error
    p_value: NDArray[np.float64]         # free-fit p-value
    slope0: NDArray[np.float64]          # forced-intercept-0 slope (MATLAB)
    rmse0: NDArray[np.float64]           # forced-0 RMSE
    r2_0: NDArray[np.float64]            # forced-0 R²
    ratio_med: NDArray[np.float64]       # median(signal/p_mol) in window = CL proxy
    ratio_std: NDArray[np.float64]       # relative std of signal/p_mol (SNR proxy)
    residual_pct: NDArray[np.float64]    # relative RMSE of forced-0 fit (Rayleigh-shape)
    rel_error: NDArray[np.float64]       # |slope - ratio_med| / ratio_med * 100
    scattering_ratio: NDArray[np.float64]  # ratio_med / clean-min  (1 = aerosol-free)
    n_pts: NDArray[np.float64]           # points per window
    temporal_cv: NDArray[np.float64]     # temporal CV of the window-mean signal/molecular
                                         # ratio across the night (aerosol proxy: molecular
                                         # is steady, aerosol fluctuates). NaN if no stack.
    signal: Optional[NDArray[np.float64]] = None   # the collapsed signal (for re-fits, e.g. bellini)
    p_mol: Optional[NDArray[np.float64]] = None    # the molecular power profile
    range_alc: Optional[NDArray[np.float64]] = None  # range (m)


@dataclass
class MethodWindow:
    """Window chosen by one detection method (+ its eligibility mask for plotting)."""
    method: str
    ok: bool
    message: str = ""
    center_m: float = np.nan
    half_m: float = np.nan
    start_m: float = np.nan
    end_m: float = np.nan
    slope: float = np.nan
    intercept: float = np.nan
    r2: float = np.nan
    std_err: float = np.nan
    p_value: float = np.nan
    cl: float = np.nan               # lidar-constant proxy = median(signal/p_mol) in window
    cl_err: float = np.nan           # absolute uncertainty = in-window std of signal/p_mol
    scattering_ratio: float = np.nan
    residual_pct: float = np.nan
    rel_error: float = np.nan
    n_pts: int = 0
    temporal_cv: float = np.nan
    eligible: Optional[NDArray[np.bool_]] = None
    grid: Optional[WindowGrid] = None
    cell_flag: Optional[NDArray[np.bool_]] = None  # (n_profiles, n_range) aerosol/cloud cells
                                                   # excluded by the "optimal" time-resolved screen
    n_clean_frac: float = np.nan                   # fraction of cells kept after flagging


# ---------------------------------------------------------------------------
# Shared grid search
# ---------------------------------------------------------------------------
def compute_window_grid(
    signal: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    half_length_options_m: tuple,
    range_start_m: float = 2000.0,
    range_end_m: float = 6000.0,
    increment_bins: int = 8,
    signal_stack: Optional[NDArray[np.float64]] = None,
) -> WindowGrid:
    """Compute every candidate window's fit and quality statistics (method-agnostic).

    If ``signal_stack`` (n_profiles × n_range, the per-profile range-normalized signal)
    is given, also compute each window's temporal variability (CV over the night of the
    window-mean signal/molecular ratio) — used by the "optimal" method to reject aerosol
    (steady molecular vs. fluctuating aerosol). Averaging over the window's range bins
    first suppresses photon noise, leaving mostly atmospheric (aerosol) variability.
    """
    dz = float(np.abs(range_alc[1] - range_alc[0])) if len(range_alc) > 1 else 1.0
    half_bins = np.unique(np.floor(np.array(half_length_options_m) / dz)).astype(int)
    half_bins = half_bins[half_bins > 0]
    c0 = int(np.floor(range_start_m / dz))
    c1 = int(np.floor(range_end_m / dz))
    center_bins = np.arange(c0, c1, increment_bins)
    nC, nL = len(center_bins), len(half_bins)

    shape = (nC, nL)
    slope = np.full(shape, np.nan)
    intercept = np.full(shape, np.nan)
    r2 = np.full(shape, np.nan)
    std_err = np.full(shape, np.nan)
    p_value = np.full(shape, np.nan)
    slope0 = np.full(shape, np.nan)
    rmse0 = np.full(shape, np.nan)
    r2_0 = np.full(shape, np.nan)
    ratio_med = np.full(shape, np.nan)
    ratio_std = np.full(shape, np.nan)
    residual_pct = np.full(shape, np.nan)
    n_pts = np.zeros(shape)
    temporal_cv = np.full(shape, np.nan)

    # Per-(time, range) ratio for the temporal-variability metric (optimal method).
    ratio_tz = None
    if signal_stack is not None and np.ndim(signal_stack) == 2 and signal_stack.shape[0] >= 3:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_tz = np.asarray(signal_stack, float) / p_mol[None, :]

    for i, cb in enumerate(center_bins):
        for j, hb in enumerate(half_bins):
            s = cb - hb
            e = min(cb + hb, len(range_alc))
            if s < 0:
                continue
            x = p_mol[s:e]
            y = signal[s:e]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 5:
                continue
            x, y = x[ok], y[ok]
            n_pts[i, j] = x.size
            # Free fit (signal = a*p_mol + b)
            try:
                a, b, r, p, se = linregress(x, y)
            except (ValueError, RuntimeWarning):
                continue
            slope[i, j] = a
            intercept[i, j] = b
            r2[i, j] = r ** 2
            std_err[i, j] = se
            p_value[i, j] = p
            # Forced-intercept-0 fit (MATLAB Auto_Calib_25)
            sxx = float(np.sum(x * x))
            if sxx > 0:
                a0 = float(np.sum(x * y) / sxx)
                slope0[i, j] = a0
                res = y - a0 * x
                sse = float(np.sum(res * res))
                rmse0[i, j] = np.sqrt(sse / x.size)
                sst = float(np.sum((y - np.mean(y)) ** 2))
                r2_0[i, j] = 1.0 - sse / sst if sst > 0 else np.nan
                denom = np.mean(np.abs(a0 * x))
                residual_pct[i, j] = (rmse0[i, j] / denom * 100.0) if denom > 0 else np.nan
            # Ratio statistics (signal / molecular)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratios = y / x
            ratios = ratios[np.isfinite(ratios)]
            if ratios.size:
                med = float(np.median(ratios))
                ratio_med[i, j] = med
                ratio_std[i, j] = float(np.std(ratios) / abs(med)) if med != 0 else np.inf
            # Temporal variability of the window-mean ratio across the night (aerosol proxy)
            if ratio_tz is not None:
                with np.errstate(invalid="ignore"):
                    r_t = np.nanmean(ratio_tz[:, s:e], axis=1)   # one value per profile
                r_t = r_t[np.isfinite(r_t)]
                if r_t.size >= 3:
                    mt = float(np.mean(r_t))
                    if mt != 0:
                        temporal_cv[i, j] = float(np.std(r_t) / abs(mt))

    # Free-fit slope-vs-median consistency (aerosol curvature) — production rel_error.
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_error = np.abs((slope - ratio_med) / ratio_med) * 100.0

    # Scattering ratio: ratio_med normalized by the cleanest (smallest) molecular-only
    # estimate. Reference taken over well-behaved windows so noise can't set a bogus min.
    clean = np.isfinite(ratio_med) & (ratio_med > 0) & (slope > 0) & (r2 >= 0.5)
    if np.any(clean):
        c_min = float(np.min(ratio_med[clean]))
    else:
        pos = ratio_med[np.isfinite(ratio_med) & (ratio_med > 0)]
        c_min = float(np.min(pos)) if pos.size else np.nan
    scattering_ratio = ratio_med / c_min if (c_min and np.isfinite(c_min)) else np.full(shape, np.nan)

    start_m = (center_bins[:, None] - half_bins[None, :]).astype(float) * dz
    return WindowGrid(
        center_m=center_bins.astype(float) * dz,
        half_m=half_bins.astype(float) * dz,
        start_m=start_m,
        slope=slope, intercept=intercept, r2=r2, std_err=std_err, p_value=p_value,
        slope0=slope0, rmse0=rmse0, r2_0=r2_0,
        ratio_med=ratio_med, ratio_std=ratio_std, residual_pct=residual_pct,
        rel_error=rel_error, scattering_ratio=scattering_ratio, n_pts=n_pts,
        temporal_cv=temporal_cv,
        signal=np.asarray(signal, float), p_mol=np.asarray(p_mol, float),
        range_alc=np.asarray(range_alc, float),
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _pack(g: WindowGrid, i: int, j: int, method: str, eligible, msg: str = "") -> MethodWindow:
    """Build a MethodWindow for grid cell (i, j)."""
    half = float(g.half_m[j])
    center = float(g.center_m[i])
    return MethodWindow(
        method=method, ok=True, message=msg,
        center_m=center, half_m=half,
        start_m=float(g.start_m[i, j]), end_m=center + half,
        slope=float(g.slope[i, j]), intercept=float(g.intercept[i, j]),
        r2=float(g.r2[i, j]), std_err=float(g.std_err[i, j]), p_value=float(g.p_value[i, j]),
        cl=float(g.ratio_med[i, j]),
        cl_err=float(g.ratio_std[i, j] * abs(g.ratio_med[i, j])),
        scattering_ratio=float(g.scattering_ratio[i, j]),
        residual_pct=float(g.residual_pct[i, j]), rel_error=float(g.rel_error[i, j]),
        n_pts=int(g.n_pts[i, j]), temporal_cv=float(g.temporal_cv[i, j]),
        eligible=eligible, grid=g,
    )


def _fail(g: WindowGrid, method: str, eligible, msg: str) -> MethodWindow:
    return MethodWindow(method=method, ok=False, message=msg, eligible=eligible, grid=g)


def _argmax2d(score: NDArray[np.float64]):
    flat = int(np.argmax(score))
    i, j = np.unravel_index(flat, score.shape)
    return int(i), int(j)


def flag_contaminated_cells(
    signal_stack: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    nmad: float = 4.0,
    min_excess: float = 0.25,
    min_clean: int = 5,
) -> NDArray[np.bool_]:
    """Flag time-altitude cells contaminated by aerosol or cloud (the "optimal" screen).

    Molecular-only scattering gives signal/p_mol ≈ CL — steady in time at each altitude.
    Aerosol and cloud add backscatter, so their cells sit ABOVE the clean level as a
    *coherent enhancement* (not random noise). For each altitude we take the temporal
    median (the typical/clean value when contamination is intermittent) and a robust spread
    (MAD), and flag cells exceeding ``median + nmad*MAD`` — but only if they also exceed the
    median by at least ``min_excess`` (so altitudes with tiny MAD are not over-flagged).
    This is an *upper-tail outlier* test, so clean molecular noise is left intact (no bias to
    the cleaned mean), while a passing cloud or a layer present only at the start/end of the
    night is removed per cell. (A purely empirical percentile threshold would instead chop the
    upper noise tail at every altitude and bias the cleaned mean low — avoided here.)

    Returns a boolean mask (n_profiles, n_range); True = contaminated (exclude).
    """
    stack = np.asarray(signal_stack, float)
    if stack.ndim != 2:
        return np.zeros_like(stack, dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = stack / p_mol[None, :]          # ≈ CL where molecular
    finite = np.isfinite(ratio)
    masked = np.where(finite, ratio, np.nan)
    n_finite = finite.sum(axis=0)
    with np.errstate(invalid="ignore"):
        med = np.nanmedian(masked, axis=0)
        mad = np.nanmedian(np.abs(masked - med[None, :]), axis=0) * 1.4826
    thr = np.maximum(med + nmad * mad, med * (1.0 + min_excess))
    bad_col = (n_finite < min_clean) | ~np.isfinite(thr) | ~(med > 0)
    thr = np.where(bad_col, np.inf, thr)
    return finite & (ratio > thr[None, :])


# ---------------------------------------------------------------------------
# The selectors
# ---------------------------------------------------------------------------
def _select_main(g: WindowGrid, **_) -> MethodWindow:
    """Legacy: center = min Σ|intercept|, half = max R². No floor (degenerate)."""
    valid_center = np.any(np.isfinite(g.r2), axis=1)
    if not np.any(valid_center):
        return _fail(g, "main", np.zeros_like(g.r2, bool), "no finite fit")
    sum_abs_b = np.where(valid_center, np.nansum(np.abs(g.intercept), axis=1), np.inf)
    i = int(np.argmin(sum_abs_b))
    j = int(np.nanargmax(g.r2[i, :]))
    elig = np.isfinite(g.r2)  # "main" has no real eligibility; show all fitted
    return _pack(g, i, j, "main", elig, "min sum|intercept| then max R2 (legacy)")


def _select_improved(g: WindowGrid, min_window_start_m=2000.0, min_r2=0.5,
                     max_rel_error=50.0, **_) -> MethodWindow:
    """Eligible (above aerosol, R²≥floor, slope>0, |b|<a, slope≈median) then max R²."""
    elig = (
        np.isfinite(g.r2)
        & (g.start_m >= min_window_start_m)
        & (g.r2 >= min_r2)
        & (g.slope > 0)
        & (np.abs(g.intercept) < g.slope)
        & (~np.isfinite(g.rel_error) | (g.rel_error <= max_rel_error))
    )
    if not np.any(elig):
        return _fail(g, "improved", elig, "no eligible window")
    i, j = _argmax2d(np.where(elig, g.r2, -np.inf))
    return _pack(g, i, j, "improved", elig, "max R2 among eligible")


def _select_matlab(g: WindowGrid, min_r2=0.5, range_end_cap_m=5000.0, **_) -> MethodWindow:
    """Auto_Calib_25: forced b=0, center = min Σ RMSE, half = max R², reject R²<floor."""
    in_cap = g.center_m <= range_end_cap_m
    rmse = np.where(in_cap[:, None] & np.isfinite(g.rmse0), g.rmse0, np.nan)
    valid_center = np.any(np.isfinite(rmse), axis=1)
    if not np.any(valid_center):
        return _fail(g, "matlab", np.isfinite(rmse), "no finite forced-0 fit")
    sum_rmse = np.where(valid_center, np.nansum(rmse, axis=1), np.inf)
    i = int(np.argmin(sum_rmse))
    j = int(np.nanargmax(np.where(np.isfinite(rmse[i, :]), g.r2_0[i, :], np.nan)))
    elig = in_cap[:, None] & np.isfinite(g.r2_0) & (g.r2_0 >= min_r2) & (g.slope0 > 0)
    w = _pack(g, i, j, "matlab", elig, "min sum-RMSE then max R2 (forced b=0)")
    # MATLAB rejects the night if the chosen fit is too poor.
    if not (np.isfinite(g.r2_0[i, j]) and g.r2_0[i, j] >= min_r2 and g.slope0[i, j] > 0):
        return _fail(g, "matlab", elig, f"R²₀={g.r2_0[i, j]:.2f} < {min_r2} (rejected)")
    w.r2 = float(g.r2_0[i, j])      # report the forced-0 R² for this method
    w.slope = float(g.slope0[i, j])
    w.intercept = 0.0
    return w


def _select_earlinet(g: WindowGrid, min_window_start_m=2000.0, max_residual_pct=10.0,
                     max_ratio_std=0.30, max_scattering_ratio=1.1, max_rel_error=15.0,
                     **_) -> MethodWindow:
    """EARLINET/SCC: Rayleigh-shape + SNR + ratio≤1.1, pick LOWEST qualifying window."""
    elig = (
        np.isfinite(g.residual_pct)
        & (g.start_m >= min_window_start_m)
        & (g.slope0 > 0)
        & (g.residual_pct <= max_residual_pct)
        & (g.ratio_std <= max_ratio_std)
        & (g.scattering_ratio <= max_scattering_ratio)
        & (np.isfinite(g.rel_error) & (g.rel_error <= max_rel_error))
    )
    if not np.any(elig):
        return _fail(g, "earlinet", elig, "no shape-matching window")
    # Lowest center (Mattis 2016: take the lowest qualifying window = best SNR).
    ci = np.where(np.any(elig, axis=1), g.center_m, np.inf)
    i = int(np.argmin(ci))
    # at that center, smallest shape residual
    resid = np.where(elig[i, :], g.residual_pct[i, :], np.inf)
    j = int(np.argmin(resid))
    return _pack(g, i, j, "earlinet", elig, "lowest qualifying (shape+SNR)")


def _select_optimal(g: WindowGrid, min_window_start_m=2000.0, min_r2=0.5,
                    max_residual_pct=12.0, max_scattering_ratio=1.1, max_ratio_std=0.30,
                    max_temporal_cv=0.5, max_rel_error=15.0, use_bg=False, w_ratio=0.25,
                    w_resid=0.20, w_snr=0.10, w_npts=0.10, w_tvar=0.35, w_rel=0.20, **_) -> MethodWindow:
    """Best-of: physical gates + TEMPORAL-variability aerosol rejection + composite score.

    Molecular scattering is steady in time; aerosol advects and fluctuates. Windows whose
    night-time temporal CV exceeds ``max_temporal_cv`` are rejected as aerosol-contaminated,
    and lower temporal variability is rewarded in the score. When no profile stack was
    supplied (temporal_cv all-NaN) the temporal gate/penalty are skipped (it then falls
    back to the single-profile composite).
    """
    have_tvar = bool(np.any(np.isfinite(g.temporal_cv)))
    elig = (
        np.isfinite(g.r2)
        & (g.start_m >= min_window_start_m)
        & (g.slope > 0)
        & (np.abs(g.intercept) < g.slope)
        & (g.r2 >= min_r2)
        & (np.isfinite(g.residual_pct) & (g.residual_pct <= max_residual_pct))
        & (np.isfinite(g.scattering_ratio) & (g.scattering_ratio <= max_scattering_ratio))
        & (g.ratio_std <= max_ratio_std)
        & (np.isfinite(g.rel_error) & (g.rel_error <= max_rel_error))
    )
    if have_tvar:
        # Reject temporally variable (aerosol) windows; keep NaN-CV windows (can't judge).
        elig = elig & (~np.isfinite(g.temporal_cv) | (g.temporal_cv <= max_temporal_cv))
    if use_bg and g.signal is not None and g.p_mol is not None and g.range_alc is not None:
        # Bellini-style vertical residual-autocorrelation gate (complements the temporal CV:
        # catches a smooth, STEADY aerosol layer the linear fit didn't flag). Applied only to
        # the already-eligible windows for speed.
        rng = g.range_alc
        dz = float(rng[1] - rng[0]) if len(rng) > 1 else 1.0
        for i, j in zip(*np.where(elig)):
            cb = int(round(g.center_m[i] / dz)); hb = int(round(g.half_m[j] / dz))
            s, e = cb - hb, min(cb + hb, len(rng))
            if s < 0:
                elig[i, j] = False
                continue
            x = g.p_mol[s:e]; y = g.signal[s:e]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 8:
                continue
            resid = y[ok] - (g.slope[i, j] * x[ok] + g.intercept[i, j])
            bg_ok, _ = _bg_no_autocorr(resid)
            if not bg_ok:
                elig[i, j] = False
    if not np.any(elig):
        return _fail(g, "optimal", elig, "no eligible window")
    n_norm = g.n_pts / (np.nanmax(g.n_pts) if np.nanmax(g.n_pts) > 0 else 1.0)
    tvar_pen = np.where(np.isfinite(g.temporal_cv), g.temporal_cv, 0.0) if have_tvar else 0.0
    rel_pen = np.where(np.isfinite(g.rel_error), g.rel_error, 0.0) / 100.0
    quality = (
        g.r2
        - w_ratio * np.abs(g.scattering_ratio - 1.0)
        - w_resid * (g.residual_pct / 100.0)
        - w_snr * g.ratio_std
        - w_tvar * tvar_pen
        - w_rel * rel_pen
        + w_npts * n_norm
    )
    i, j = _argmax2d(np.where(elig, quality, -np.inf))
    msg = "max composite quality (R2+shape+purity+SNR" + ("+temporal)" if have_tvar else ")")
    return _pack(g, i, j, "optimal", elig, msg)


def _bg_no_autocorr(resid: NDArray[np.float64], alpha: float = 1.96):
    """Breusch-Godfrey proxy: True if the residuals show no significant (positive) lag-1
    autocorrelation. Random residuals = molecular; coherent (smooth) structure = an
    undetected aerosol layer, which makes residuals positively autocorrelated."""
    n = resid.size
    if n < 8 or np.std(resid) <= 0:
        return False, np.nan
    r1 = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])
    if not np.isfinite(r1):
        return False, np.nan
    return (r1 <= alpha / np.sqrt(n)), r1   # one-sided ~ BG p>0.05 (no autocorrelation)


def _select_bellini(g: WindowGrid, min_window_start_m=3000.0, max_window_end_m=7000.0,
                    min_width_m=600.0, max_width_m=3000.0, max_rel_error=15.0,
                    max_ecl=0.40, border_m=200.0, **_) -> MethodWindow:
    """ALICENET / Bellini et al. (2024, AMT 17, 6119) Rayleigh calibration (Supplement S3).

    Windows lie within 3-7 km with widths 600-3000 m. Per window the signal is linearly fit to
    the molecular profile and rejected unless: the residuals show no coherent structure
    (Breusch-Godfrey no-autocorrelation proxy, QC.CAL1); slope > 0 and intercept ~ 0 (QC.CAL2);
    and the border residuals (+/-200 m of each edge) are predominantly negative (QC.CAL3). Among
    retained windows the one maximising M_Ray = (adjR^2 + (1-|b|))/std(b) is selected, then
    rejected if the relative calibration uncertainty E_CL = err(slope)/slope + std(CL)/median(CL)
    exceeds 0.40 (QC.CAL4). QC.CAL5 ('negative AOD') is enforced by the pipeline's Klett step.
    """
    sig, pmol, rng = g.signal, g.p_mol, g.range_alc
    if sig is None or pmol is None or rng is None:
        return _select_improved(g)   # inputs not stored -> safe fallback
    dz = float(rng[1] - rng[0]) if len(rng) > 1 else 1.0
    width = 2.0 * g.half_m[None, :] * np.ones_like(g.center_m[:, None])
    end_grid = g.center_m[:, None] + g.half_m[None, :]
    cand = (
        np.isfinite(g.r2) & (g.slope > 0)
        & (g.start_m >= min_window_start_m) & (end_grid <= max_window_end_m)
        & (width >= min_width_m) & (width <= max_width_m)
    )
    if not np.any(cand):
        return _fail(g, "bellini", cand, "no candidate window in 3-7 km")
    adjr2 = np.full(g.r2.shape, np.nan)
    b_rel = np.full(g.r2.shape, np.inf)
    keep = np.zeros(g.r2.shape, bool)
    bw = max(int(round(border_m / dz)), 1)
    for i, j in zip(*np.where(cand)):
        cb = int(round(g.center_m[i] / dz))
        hb = int(round(g.half_m[j] / dz))
        s, e = cb - hb, min(cb + hb, len(rng))
        if s < 0:
            continue
        x, y = pmol[s:e], sig[s:e]
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 8:
            continue
        x, y = x[ok], y[ok]
        n = x.size
        a, b = g.slope[i, j], g.intercept[i, j]
        resid = y - (a * x + b)
        adjr2[i, j] = 1.0 - (1.0 - g.r2[i, j]) * (n - 1) / (n - 2) if n > 2 else np.nan
        denom = abs(a) * np.mean(x)
        b_rel[i, j] = abs(b) / denom if denom > 0 else np.inf
        bg_ok, _r1 = _bg_no_autocorr(resid)
        lo = slice(max(s - bw, 0), min(s + bw, len(rng)))   # +/-200 m of each border
        hi = slice(max(e - bw, 0), min(e + bw, len(rng)))
        bx = np.concatenate([pmol[lo], pmol[hi]])
        by = np.concatenate([sig[lo], sig[hi]])
        m2 = np.isfinite(bx) & np.isfinite(by)
        border_ok = (np.sum(np.sign(by[m2] - (a * bx[m2] + b))) < 0) if m2.sum() >= 3 else True
        keep[i, j] = bool(bg_ok and border_ok and np.isfinite(adjr2[i, j]))
    keep = keep & (g.rel_error <= max_rel_error)   # QC.CAL2 intercept ~ 0 (proportionality)
    if not np.any(keep):
        return _fail(g, "bellini", keep, "no window passed QC.CAL1-3")
    std_brel = np.nanstd(b_rel[keep])
    std_brel = std_brel if (np.isfinite(std_brel) and std_brel > 0) else 1.0
    mray = np.where(keep, (adjr2 + (1.0 - b_rel)) / std_brel, -np.inf)
    i, j = _argmax2d(mray)
    cb = int(round(g.center_m[i] / dz)); hb = int(round(g.half_m[j] / dz))
    s, e = cb - hb, min(cb + hb, len(rng))
    x, y = pmol[s:e], sig[s:e]
    ok = np.isfinite(x) & np.isfinite(y)
    ratio = y[ok] / x[ok]
    slope_term = (g.std_err[i, j] / g.slope[i, j]) if g.slope[i, j] > 0 else np.inf
    spread_term = (np.std(ratio) / abs(np.median(ratio))) if (ratio.size and np.median(ratio) != 0) else np.inf
    ecl = slope_term + spread_term
    if ecl > max_ecl:
        return _fail(g, "bellini", keep, f"E_CL={ecl:.2f} > {max_ecl} (QC.CAL4)")
    return _pack(g, i, j, "bellini", keep, "ALICENET: max M_Ray (adjR2 + BG no-autocorr + border-sign)")


# Keyed by E-PROF version; the private selectors keep their algorithm names.
_SELECTORS = {
    "eprof_v1.1": _select_main,
    "eprof_v1.2": _select_improved,
    "eprof_v0.25": _select_matlab,
    "earlinet": _select_earlinet,
    "eprof_v2": _select_optimal,
    "bellini": _select_bellini,
}


def select_molecular_window(
    method: str,
    signal: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    half_length_options_m: tuple,
    range_start_m: float = 2000.0,
    range_end_m: float = 6000.0,
    increment_bins: int = 8,
    grid: Optional[WindowGrid] = None,
    signal_stack: Optional[NDArray[np.float64]] = None,
    **params,
) -> MethodWindow:
    """Detect the molecular window with the chosen ``method``.

    Pass a precomputed ``grid`` to evaluate several methods on one profile cheaply.
    Pass ``signal_stack`` (n_profiles × n_range) to enable the "optimal" method's
    temporal-variability aerosol rejection. Extra keyword args override DEFAULT_PARAMS.
    """
    method = resolve_method(method)
    if method not in _SELECTORS:
        raise ValueError(f"Unknown molecular method '{method}'; choose from {METHODS}")
    merged = dict(DEFAULT_PARAMS.get(method, {}))
    merged.update(params)

    # E-PROF v2 ("optimal") does TIME-RESOLVED aerosol/cloud flagging: it excludes contaminated
    # cells (e.g. aerosol only at the start of the night, a cloud only at the end), then fits the
    # molecular window on the time-cleaned mean profile -> it uses the clean part of an
    # otherwise-contaminated night. The other methods use the full night mean.
    if (method == "eprof_v2" and signal_stack is not None
            and np.ndim(signal_stack) == 2 and np.shape(signal_stack)[0] >= 5):
        flag = flag_contaminated_cells(
            signal_stack, p_mol, range_alc,
            nmad=merged.get("flag_nmad", 4.0),
            min_excess=merged.get("flag_min_excess", 0.25),
        )
        masked = np.where(flag, np.nan, np.asarray(signal_stack, float))
        with np.errstate(invalid="ignore"):
            clean_signal = np.nanmean(masked, axis=0)
        g = compute_window_grid(
            clean_signal, p_mol, range_alc, half_length_options_m,
            range_start_m=range_start_m, range_end_m=range_end_m,
            increment_bins=increment_bins, signal_stack=masked,
        )
        mw = _select_optimal(g, **merged)
        mw.cell_flag = flag
        mw.n_clean_frac = float(1.0 - flag.mean()) if flag.size else np.nan
        return mw

    if grid is None:
        grid = compute_window_grid(
            signal, p_mol, range_alc, half_length_options_m,
            range_start_m=range_start_m, range_end_m=range_end_m,
            increment_bins=increment_bins, signal_stack=signal_stack,
        )
    return _SELECTORS[method](grid, **merged)
