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


def find_optimal_molecular_window(
    signal: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    half_length_options_m: tuple,
    range_start_m: float = 2000.0,
    range_end_m: float = 6000.0,
    increment_bins: int = 8,
) -> RayleighFitResult:
    """
    Find the optimal molecular scattering window using grid search.

    The algorithm searches over different center positions and window sizes
    to find the region where the signal best matches theoretical molecular
    scattering (linear relationship with minimal offset).

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

    Returns
    -------
    RayleighFitResult
        Optimal fit parameters and window location.
    """
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
            except (ValueError, RuntimeWarning):
                continue

    # Find optimal window:
    # 1. First, find center with minimum sum of |intercept| across all window sizes
    sum_abs_intercept = np.nansum(np.abs(intercepts), axis=1)
    best_center_idx = np.nanargmin(sum_abs_intercept)

    # 2. Then, find window size with maximum R² at that center
    best_length_idx = np.nanargmax(r_squared[best_center_idx, :])

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
