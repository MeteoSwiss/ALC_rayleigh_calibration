"""
Rayleigh fitting algorithms for lidar calibration.

This module contains the core fitting routines for finding the optimal
molecular scattering region and calculating the lidar constant.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.stats import linregress


@dataclass
class WindowSearchDiagnostics:
    """Grid-search arrays produced by the molecular-window optimisation."""
    range_bin_m: NDArray[np.float64]        # Center positions searched (m)
    half_length_m: NDArray[np.float64]      # Half-lengths searched (m)
    slopes: NDArray[np.float64]             # (n_centers, n_lengths)
    intercepts: NDArray[np.float64]         # (n_centers, n_lengths)
    r_squared: NDArray[np.float64]          # (n_centers, n_lengths)
    sum_abs_intercept: NDArray[np.float64]  # (n_centers,)
    valid_window: Optional[NDArray[np.bool_]] = None  # (n_centers, n_lengths) eligible-window mask


@dataclass
class RayleighFitResult:
    """Result of Rayleigh fit optimization."""
    # Fit parameters
    slope: float                  # Linear fit slope (a)
    intercept: float              # Linear fit intercept (b)
    r_squared: float             # Coefficient of determination
    std_error: float             # Standard error of slope
    p_value: float               # p-value of fit

    # Optimal window parameters
    center_range_m: float        # Center of molecular window (m)
    half_length_m: float         # Half-length of window (m)
    range_start_m: float         # Start of molecular window (m)
    range_end_m: float           # End of molecular window (m)

    # Altitude (above sea level)
    altitude_start: float        # Bottom of window (m ASL)
    altitude_end: float          # Top of window (m ASL)

    # Quality metrics
    relative_error: float        # Relative error between fit and median

    # Optional diagnostics for plotting
    search_diagnostics: Optional[WindowSearchDiagnostics] = None

    @property
    def is_valid(self) -> bool:
        """Check if fit result is physically valid."""
        return self.slope > 0 and abs(self.intercept) < self.slope


@dataclass
class CalibrationConstantResult:
    """Result of lidar constant calculation."""
    lidar_constant: float
    uncertainty: float
    cl_profile: NDArray[np.float64]  # CL at each altitude
    molecular_start_idx: int
    molecular_end_idx: int


def _result_from_method_window(mw) -> RayleighFitResult:
    """Map a molecular_methods.MethodWindow onto a RayleighFitResult (+ diagnostics)."""
    g = mw.grid
    diagnostics = None
    if g is not None:
        sum_abs_intercept = np.nansum(np.abs(g.intercept), axis=1)
        diagnostics = WindowSearchDiagnostics(
            range_bin_m=np.asarray(g.center_m, float),
            half_length_m=np.asarray(g.half_m, float),
            slopes=g.slope, intercepts=g.intercept, r_squared=g.r2,
            sum_abs_intercept=sum_abs_intercept, valid_window=mw.eligible,
        )
    if not mw.ok:
        return RayleighFitResult(
            slope=np.nan, intercept=np.nan, r_squared=np.nan, std_error=np.nan,
            p_value=np.nan, center_range_m=np.nan, half_length_m=np.nan,
            range_start_m=np.nan, range_end_m=np.nan, altitude_start=0, altitude_end=0,
            relative_error=np.inf, search_diagnostics=diagnostics,
        )
    rel = mw.rel_error if np.isfinite(mw.rel_error) else 0.0
    return RayleighFitResult(
        slope=mw.slope, intercept=mw.intercept, r_squared=mw.r2, std_error=mw.std_err,
        p_value=mw.p_value, center_range_m=mw.center_m, half_length_m=mw.half_m,
        range_start_m=mw.start_m, range_end_m=mw.end_m, altitude_start=0, altitude_end=0,
        relative_error=rel, search_diagnostics=diagnostics,
    )


def find_optimal_molecular_window(
    signal: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    half_length_options_m: tuple,
    range_start_m: float = 2000.0,
    range_end_m: float = 6000.0,
    increment_bins: int = 8,
    min_window_start_m: float = 2000.0,
    min_r2: float = 0.5,
    max_rel_error: float = 50.0,
    method: str = "improved",
    signal_stack: Optional[NDArray[np.float64]] = None,
) -> RayleighFitResult:
    """
    Find the optimal molecular scattering window using grid search.

    The algorithm searches over different center positions and window sizes,
    then selects, among windows that are *physically* molecular, the one with
    the best linear signal-vs-molecular fit (largest R^2).

    Selection criterion (and why it changed)
    ----------------------------------------
    For each (center, half-length) window we regress the range-normalized signal
    ``y`` on the theoretical molecular power ``x`` (``y = a*x + b``). A genuine
    molecular window has a *high R^2* (signal proportional to molecular), a
    positive slope ``a``, a small background offset ``b``, and lies above the
    boundary-layer aerosol.

    Earlier versions chose the center with the smallest Σ|b| (intercept), then the
    half-length with the largest R^2. With a *free* intercept this is degenerate:
    in the high-altitude noise region the signal → 0, so any line fits it with
    ``b ≈ 0`` (trivially small) *and* ``R^2 ≈ 0`` (no real correlation) — so
    min-Σ|b| systematically selected a low-R^2, noise-dominated window. This is the
    documented EARLINET-SCC minimum-signal failure mode (Mattis, D'Amico, Baars
    et al., AMT 9, 3009, 2016); the original MATLAB ``Auto_Calib_25/rayleigh_fit.m``
    avoided it by forcing ``b = 0``, choosing by RMSE, and rejecting fits with
    ``R^2 < min_r2_rfit (=0.5)``.

    We restore that robustness: keep only windows that pass molecular-validity
    gates, then pick the **largest R^2** among them. R^2 is the right metric here
    precisely because it *collapses* in the noise region (so noise windows are
    rejected, not selected).

    Parameters
    ----------
    signal : ndarray
        Range-normalized signal (RCS / r^2).
    p_mol : ndarray
        Theoretical molecular power profile.
    range_alc : ndarray
        Range bins in meters.
    half_length_options_m : tuple
        Possible half-lengths for the molecular window (m).
    range_start_m : float
        Minimum center range to consider (m).
    range_end_m : float
        Maximum center range to consider (m).
    increment_bins : int
        Step size in bins for center position search.
    min_window_start_m : float
        Minimum altitude (m AGL) at which a window may START — keeps the window
        above the boundary-layer aerosol (rec #1). Note this constrains the window
        start, unlike ``range_start_m`` which constrains the window center.
    min_r2 : float
        Minimum R^2 for a window to be eligible (cf. MATLAB ``min_r2_rfit``).
        Rejects noise-dominated and strongly aerosol-curved windows.
    max_rel_error : float
        Maximum |slope − pointwise-median-ratio| / median-ratio (%). Rejects
        windows where the signal is not truly proportional to molecular, i.e.
        aerosol curvature (rec #2).

    Returns
    -------
    RayleighFitResult
        Optimal fit parameters and window location. If no window passes the
        validity gates, returns a failed fit (relative_error = inf) so the caller
        flags a non-calibration night instead of emitting a spurious constant.

    method : str
        Detection strategy: "improved" (default; the in-line implementation below),
        or one of "main"/"matlab"/"calipso"/"earlinet"/"optimal", dispatched to
        :mod:`rayleigh_calibration.rayleigh.molecular_methods`.
    """
    # Non-"improved" methods are dispatched to the pluggable selectors. "improved"
    # keeps the in-line implementation below (the production path, unchanged).
    if method != "improved":
        from .molecular_methods import select_molecular_window
        mw = select_molecular_window(
            method, signal, p_mol, range_alc, half_length_options_m,
            range_start_m=range_start_m, range_end_m=range_end_m,
            increment_bins=increment_bins, signal_stack=signal_stack,
        )
        return _result_from_method_window(mw)

    dz = np.abs(range_alc[1] - range_alc[0]) if len(range_alc) > 1 else 1.0

    # Convert parameters to bin indices
    half_length_bins = np.unique(np.floor(np.array(half_length_options_m) / dz)).astype(int)
    range_start_bin = int(np.floor(range_start_m / dz))
    range_end_bin = int(np.floor(range_end_m / dz))

    center_bins = np.arange(range_start_bin, range_end_bin, increment_bins)
    n_centers = len(center_bins)
    n_lengths = len(half_length_bins)

    # Pre-allocate result arrays
    slopes = np.full((n_centers, n_lengths), np.nan)
    intercepts = np.full((n_centers, n_lengths), np.nan)
    r_squared = np.full((n_centers, n_lengths), np.nan)
    std_errors = np.full((n_centers, n_lengths), np.nan)
    p_values = np.full((n_centers, n_lengths), np.nan)
    rel_errors = np.full((n_centers, n_lengths), np.nan)  # |slope - median ratio| / median (%)

    # Grid search over center positions and window sizes
    for i, center_bin in enumerate(center_bins):
        for j, half_len in enumerate(half_length_bins):
            start_bin = center_bin - half_len
            end_bin = min(center_bin + half_len, len(range_alc))

            if start_bin < 0:
                continue

            x = p_mol[start_bin:end_bin]
            y = signal[start_bin:end_bin]

            # Skip if any NaN values
            if np.any(np.isnan(x)) or np.any(np.isnan(y)):
                continue

            # Perform linear regression
            try:
                a, b, r, p, se = linregress(x, y)
                slopes[i, j] = a
                intercepts[i, j] = b
                r_squared[i, j] = r ** 2
                std_errors[i, j] = se
                p_values[i, j] = p
                # Slope vs. pointwise-median ratio: large where the signal is not
                # truly proportional to molecular (aerosol curvature) -> rec #2 gate.
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratios = y / x
                ratios = ratios[np.isfinite(ratios)]
                if ratios.size:
                    med_ratio = np.median(ratios)
                    if med_ratio != 0:
                        rel_errors[i, j] = abs((a - med_ratio) / med_ratio * 100.0)
            except (ValueError, RuntimeWarning):
                continue

    # ── Select the optimal molecular window ──────────────────────────────────
    # Keep only physically valid molecular windows, then choose the best linear
    # fit (largest R²) among them. See the function docstring for why min-Σ|b|
    # (the previous criterion) is degenerate at high altitude.
    #   (rec #1) starts above the boundary-layer aerosol : start >= min_window_start_m
    #   (R² fix) genuine signal~molecular correlation    : R² >= min_r2  (cf. MATLAB min_r2_rfit)
    #            signal increasing with molecular         : slope > 0
    #            small background offset vs. slope        : |intercept| < slope
    #   (rec #2) proportional to molecular (no curvature) : relative_error <= max_rel_error
    start_range_grid = (center_bins[:, None] - half_length_bins[None, :]).astype(float) * dz
    valid_window = (
        np.isfinite(r_squared)
        & (start_range_grid >= min_window_start_m)
        & (r_squared >= min_r2)
        & (slopes > 0)
        & (np.abs(intercepts) < slopes)
        & (~np.isfinite(rel_errors) | (rel_errors <= max_rel_error))
    )

    # Per-center sum of |intercept| kept only for the window-search diagnostic plot.
    valid_center = np.any(np.isfinite(r_squared), axis=1)
    sum_abs_intercept = np.where(valid_center, np.nansum(np.abs(intercepts), axis=1), np.inf)

    if not np.any(valid_window):
        # No molecular window passes the validity gates (clouds/aerosol fill the band,
        # or the molecular SNR is too low everywhere). Return a failed fit
        # (relative_error=inf) so the caller flags a non-calibration night, matching the
        # MATLAB R²/RMSE rejection rather than emitting a spurious low-R² constant.
        diagnostics = WindowSearchDiagnostics(
            range_bin_m=center_bins.astype(float) * dz,
            half_length_m=np.array(half_length_options_m, dtype=float),
            slopes=slopes, intercepts=intercepts, r_squared=r_squared,
            sum_abs_intercept=sum_abs_intercept, valid_window=valid_window,
        )
        return RayleighFitResult(
            slope=np.nan, intercept=np.nan, r_squared=np.nan, std_error=np.nan,
            p_value=np.nan, center_range_m=np.nan, half_length_m=np.nan,
            range_start_m=np.nan, range_end_m=np.nan, altitude_start=0, altitude_end=0,
            relative_error=np.inf, search_diagnostics=diagnostics,
        )

    # Best window = largest R² among the valid ones. R² collapses in the noise
    # region, so this can never select the high-altitude degenerate window.
    masked_r2 = np.where(valid_window, r_squared, -np.inf)
    flat_idx = int(np.argmax(masked_r2))
    best_center_idx, best_length_idx = np.unravel_index(flat_idx, masked_r2.shape)
    best_center_idx = int(best_center_idx)
    best_length_idx = int(best_length_idx)

    # Extract best fit parameters
    best_slope = slopes[best_center_idx, best_length_idx]
    best_intercept = intercepts[best_center_idx, best_length_idx]
    best_r2 = r_squared[best_center_idx, best_length_idx]
    best_stderr = std_errors[best_center_idx, best_length_idx]
    best_pvalue = p_values[best_center_idx, best_length_idx]

    best_center_m = center_bins[best_center_idx] * dz
    best_half_m = half_length_options_m[best_length_idx]

    # Calculate relative error between fit slope and median ratio
    center_bin = center_bins[best_center_idx]
    half_bin = half_length_bins[best_length_idx]
    start_bin = center_bin - half_bin
    end_bin = min(center_bin + half_bin, len(range_alc))

    x = p_mol[start_bin:end_bin]
    y = signal[start_bin:end_bin]
    valid = ~(np.isnan(x) | np.isnan(y))

    if np.any(valid):
        median_ratio = np.median(y[valid] / x[valid])
        relative_error = abs((best_slope - median_ratio) / median_ratio * 100)
    else:
        relative_error = np.inf

    # Store grid-search diagnostics for plotting
    diagnostics = WindowSearchDiagnostics(
        range_bin_m=center_bins.astype(float) * dz,
        half_length_m=np.array(half_length_options_m, dtype=float),
        slopes=slopes,
        intercepts=intercepts,
        r_squared=r_squared,
        sum_abs_intercept=sum_abs_intercept,
        valid_window=valid_window,
    )

    return RayleighFitResult(
        slope=best_slope,
        intercept=best_intercept,
        r_squared=best_r2,
        std_error=best_stderr,
        p_value=best_pvalue,
        center_range_m=best_center_m,
        half_length_m=best_half_m,
        range_start_m=best_center_m - best_half_m,
        range_end_m=best_center_m + best_half_m,
        altitude_start=0,  # Will be set by caller
        altitude_end=0,    # Will be set by caller
        relative_error=relative_error,
        search_diagnostics=diagnostics,
    )


def calculate_lidar_constant(
    rcs_mean: NDArray[np.float64],
    beta_tot: NDArray[np.float64],
    ext_tot: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    molecular_mask: NDArray[np.bool_],
    fit_result: RayleighFitResult,
    subtract_background: bool = False,
) -> CalibrationConstantResult:
    """
    Calculate the lidar constant using the Wiegner and Geiss (2012) method.

    Parameters
    ----------
    rcs_mean : ndarray
        Time-averaged range-corrected signal.
    beta_tot : ndarray
        Total (molecular + aerosol) backscatter profile.
    ext_tot : ndarray
        Total extinction profile.
    range_alc : ndarray
        Range bins in meters.
    molecular_mask : ndarray
        Boolean mask indicating molecular region.
    fit_result : RayleighFitResult
        Result from Rayleigh fit optimization.
    subtract_background : bool
        Whether to subtract background from signal.

    Returns
    -------
    CalibrationConstantResult
        Lidar constant and uncertainty.
    """
    mol_indices = np.where(molecular_mask)[0]
    if len(mol_indices) == 0:
        raise ValueError("No molecular region found")

    i_start = mol_indices[0]
    i_end = mol_indices[-1]

    # Optionally subtract background
    if subtract_background:
        signal_norm = rcs_mean / (range_alc ** 2)
        signal_norm = signal_norm - fit_result.intercept
        rcs_corrected = signal_norm * (range_alc ** 2)
    else:
        rcs_corrected = rcs_mean

    # Calculate CL at each altitude using transmission correction.
    # Use cumulative trapezoidal integration instead of per-bin np.trapz.
    cl_profile = np.full_like(range_alc, np.nan)

    # Cumulative optical depth from surface to each bin
    if i_end > 1:
        dz_half = np.diff(range_alc[:i_end]) / 2.0
        mid_sum = ext_tot[:i_end - 1] + ext_tot[1:i_end]
        optical_depth = np.zeros(i_end)
        np.cumsum(dz_half * mid_sum, out=optical_depth[1:])
    else:
        optical_depth = np.zeros(max(i_end, 1))

    inv_transmission = np.exp(2 * optical_depth)

    valid_mask = beta_tot[:i_end] > 0
    cl_profile[:i_end] = np.where(
        valid_mask,
        rcs_corrected[:i_end] / np.where(valid_mask, beta_tot[:i_end], 1.0) * inv_transmission,
        np.nan,
    )

    # Calculate median CL in molecular region
    cl_molecular = cl_profile[i_start:i_end]
    valid_cl = cl_molecular[~np.isnan(cl_molecular)]

    if len(valid_cl) == 0:
        raise ValueError("No valid CL values in molecular region")

    lidar_constant = np.median(valid_cl)

    # Calculate uncertainty
    # Includes both scatter in molecular region and fit error
    fit_error = fit_result.std_error / fit_result.slope if fit_result.slope > 0 else 0
    uncertainty = (np.std(valid_cl) + fit_error * lidar_constant) * 2

    return CalibrationConstantResult(
        lidar_constant=lidar_constant,
        uncertainty=uncertainty,
        cl_profile=cl_profile,
        molecular_start_idx=i_start,
        molecular_end_idx=i_end,
    )


def validate_calibration(
    cl_result: CalibrationConstantResult,
    fit_result: RayleighFitResult,
    inv_transmission_ref: float,
    threshold: float = 15.0,
) -> Tuple[bool, float, str]:
    """
    Validate calibration result using slope method comparison.

    Parameters
    ----------
    cl_result : CalibrationConstantResult
        Result from lidar constant calculation.
    fit_result : RayleighFitResult
        Result from Rayleigh fit.
    inv_transmission_ref : float
        Inverse transmission at reference altitude.
    threshold : float
        Maximum allowed relative error (%).

    Returns
    -------
    tuple
        (is_valid, error_percent, message)
    """
    # Calculate CL using slope method
    cl_slope = fit_result.slope * inv_transmission_ref

    # Relative error between methods
    error_percent = abs((cl_slope - cl_result.lidar_constant) / cl_result.lidar_constant * 100)

    if error_percent > threshold:
        return False, error_percent, f"Method disagreement: {error_percent:.1f}% > {threshold}%"

    if cl_result.uncertainty > cl_result.lidar_constant:
        return False, error_percent, f"Uncertainty exceeds value: {cl_result.uncertainty:.2e} > {cl_result.lidar_constant:.2e}"

    return True, error_percent, "Calibration validated"
