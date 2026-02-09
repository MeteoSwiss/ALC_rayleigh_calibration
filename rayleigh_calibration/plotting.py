"""
Plotting utilities for Rayleigh calibration diagnostics.

This module provides optional visualization functions for diagnosing
calibration results. All plotting is optional and can be disabled
via the CalibrationOptions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

# Lazy import matplotlib to avoid dependency if not needed
if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure


def _get_plt():
    """Lazy import of matplotlib."""
    import matplotlib.pyplot as plt
    return plt


def plot_rcs_timeseries(
    hours_since_start: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    rcs: NDArray[np.float64],
    cbh: Optional[NDArray[np.float64]] = None,
    no_cloud_value: float = -9,
    title: str = "",
    save_path: Optional[Path] = None,
) -> "Figure":
    """
    Plot range-corrected signal time series with cloud base heights.
    
    Parameters
    ----------
    hours_since_start : ndarray
        Time axis in hours since start.
    range_alc : ndarray
        Range bins in meters.
    rcs : ndarray
        Range-corrected signal (time x range).
    cbh : ndarray, optional
        Cloud base heights (time x layers).
    no_cloud_value : float
        Value indicating no cloud.
    title : str
        Plot title.
    save_path : Path, optional
        If provided, save figure to this path.
        
    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    plt = _get_plt()
    
    fig, ax = plt.subplots(figsize=(12, 5))
    
    # Plot log10 of positive RCS values
    rcs_pos = rcs.T.copy()
    rcs_pos[rcs_pos <= 0] = np.nan
    log_rcs = np.log10(rcs_pos)
    
    im = ax.pcolormesh(
        hours_since_start, range_alc, log_rcs,
        vmin=0, vmax=6,
        shading='auto',
        cmap='viridis',
    )
    
    cbar = plt.colorbar(im, ax=ax, label='log10(RCS)')
    
    # Overlay cloud base heights
    if cbh is not None:
        cbh_plot = cbh[:, 0].copy()
        cbh_plot[cbh_plot == no_cloud_value] = np.nan
        ax.plot(hours_since_start, cbh_plot, '.', color='white', alpha=0.5, markersize=2)
    
    ax.set_xlabel('Hours since start')
    ax.set_ylabel('Range (m)')
    ax.set_title(title)
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_molecular_fit(
    range_alc: NDArray[np.float64],
    altitude: float,
    rcs_mean: NDArray[np.float64],
    signal_normalized: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    beta_att_mol: NDArray[np.float64],
    fit_altitude_start: float,
    fit_altitude_end: float,
    title: str = "",
    save_path: Optional[Path] = None,
) -> "Figure":
    """
    Plot Rayleigh fit results showing signal vs molecular theory.
    
    Parameters
    ----------
    range_alc : ndarray
        Range bins in meters.
    altitude : float
        Station altitude in meters.
    rcs_mean : ndarray
        Mean range-corrected signal.
    signal_normalized : ndarray
        Normalized signal (after Rayleigh fit).
    p_mol : ndarray
        Theoretical molecular power.
    beta_att_mol : ndarray
        Attenuated molecular backscatter.
    fit_altitude_start : float
        Bottom of molecular fit window (m ASL).
    fit_altitude_end : float
        Top of molecular fit window (m ASL).
    title : str
        Plot title.
    save_path : Path, optional
        If provided, save figure to this path.
        
    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    plt = _get_plt()
    
    z = range_alc + altitude
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 8))
    
    # Panel 1: Raw RCS
    ax1 = axes[0]
    ax1.plot(rcs_mean, z, 'b-', linewidth=0.5)
    ax1.axhline(fit_altitude_start, color='k', linestyle='--', alpha=0.5)
    ax1.axhline(fit_altitude_end, color='k', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Raw RCS')
    ax1.set_ylabel('Altitude (m ASL)')
    ax1.set_title('Range-Corrected Signal')
    
    # Panel 2: Normalized signal vs molecular
    ax2 = axes[1]
    ax2.plot(signal_normalized * 1e6, z, 'b-', label='Normalized signal', linewidth=0.5)
    ax2.plot(p_mol * 1e6, z, 'r-', label='Molecular theory', linewidth=0.5)
    ax2.axhline(fit_altitude_start, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(fit_altitude_end, color='k', linestyle='--', alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_xlabel('Signal (Mm⁻¹)')
    ax2.set_title('Signal Comparison')
    ax2.legend(loc='upper right')
    
    # Panel 3: Attenuated backscatter
    ax3 = axes[2]
    beta_att = signal_normalized * (range_alc ** 2)
    ax3.plot(beta_att * 1e6, z, 'b-', label='Retrieved', linewidth=0.5)
    ax3.plot(beta_att_mol * 1e6, z, 'r-', label='Molecular', linewidth=0.5)
    ax3.axhline(fit_altitude_start, color='k', linestyle='--', alpha=0.5)
    ax3.axhline(fit_altitude_end, color='k', linestyle='--', alpha=0.5)
    ax3.set_xlabel('β_att (Mm⁻¹ sr⁻¹)')
    ax3.set_title('Attenuated Backscatter')
    ax3.legend(loc='upper right')
    
    fig.suptitle(title)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_lidar_constant(
    range_alc: NDArray[np.float64],
    cl_profile: NDArray[np.float64],
    cl_median: float,
    cl_slope: float,
    cl_uncertainty: float,
    fit_range_start: float,
    fit_range_end: float,
    title: str = "",
    save_path: Optional[Path] = None,
) -> "Figure":
    """
    Plot lidar constant profile with comparison methods.
    
    Parameters
    ----------
    range_alc : ndarray
        Range bins in meters.
    cl_profile : ndarray
        Lidar constant at each altitude.
    cl_median : float
        Median lidar constant.
    cl_slope : float
        Lidar constant from slope method.
    cl_uncertainty : float
        Uncertainty on lidar constant.
    fit_range_start : float
        Bottom of molecular fit window (m).
    fit_range_end : float
        Top of molecular fit window (m).
    title : str
        Plot title.
    save_path : Path, optional
        If provided, save figure to this path.
        
    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    plt = _get_plt()
    
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Plot CL profile
    ax.plot(cl_profile, range_alc, 'b-', linewidth=0.5, label='CL(z)')
    
    # Plot median and slope values
    ax.axvline(cl_median, color='r', linestyle='-', label=f'Median: {cl_median:.3e}')
    ax.axvline(cl_slope, color='g', linestyle='--', label=f'Slope: {cl_slope:.3e}')
    
    # Uncertainty bounds
    ax.axvline(cl_median - cl_uncertainty, color='r', linestyle=':', alpha=0.5)
    ax.axvline(cl_median + cl_uncertainty, color='r', linestyle=':', alpha=0.5)
    
    # Molecular window
    ax.axhline(fit_range_start, color='k', linestyle='--', alpha=0.5)
    ax.axhline(fit_range_end, color='k', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Lidar Constant')
    ax.set_ylabel('Range (m)')
    ax.set_title(title)
    ax.legend(loc='upper right')
    
    # Set reasonable x-limits
    valid_cl = cl_profile[~np.isnan(cl_profile)]
    if len(valid_cl) > 0:
        x_margin = (np.max(valid_cl) - np.min(valid_cl)) * 0.2
        ax.set_xlim(np.min(valid_cl) - x_margin, np.max(valid_cl) + x_margin)
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_rayleigh_window_search(
    range_bin_m: NDArray[np.float64],
    half_length_m: NDArray[np.float64],
    slopes: NDArray[np.float64],
    intercepts: NDArray[np.float64],
    r_squared: NDArray[np.float64],
    sum_abs_intercept: NDArray[np.float64],
    best_range_m: float,
    best_half_m: float,
    title: str = "",
    save_path: Optional[Path] = None,
) -> "Figure":
    """
    Plot Rayleigh fit window search results.
    
    Parameters
    ----------
    range_bin_m : ndarray
        Center positions searched (m).
    half_length_m : ndarray
        Half-lengths searched (m).
    slopes : ndarray
        Fit slopes for each window.
    intercepts : ndarray
        Fit intercepts for each window.
    r_squared : ndarray
        R² values for each window.
    sum_abs_intercept : ndarray
        Sum of |intercept| for each center position.
    best_range_m : float
        Optimal center position (m).
    best_half_m : float
        Optimal half-length (m).
    title : str
        Plot title.
    save_path : Path, optional
        If provided, save figure to this path.
        
    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    plt = _get_plt()
    
    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    
    # Panel 1: Slopes
    ax1 = axes[0]
    im1 = ax1.pcolormesh(range_bin_m, half_length_m, slopes.T, shading='auto')
    plt.colorbar(im1, ax=ax1, label='Slope (a)')
    ax1.set_ylabel('Half-length (m)')
    ax1.set_title('Fit Slope')
    
    # Panel 2: Absolute intercepts
    ax2 = axes[1]
    im2 = ax2.pcolormesh(range_bin_m, half_length_m, np.abs(intercepts.T), shading='auto')
    plt.colorbar(im2, ax=ax2, label='|Intercept| (b)')
    ax2.set_ylabel('Half-length (m)')
    ax2.set_title('Fit Intercept (absolute)')
    
    # Panel 3: Sum of |b| vs center position
    ax3 = axes[2]
    ax3.plot(range_bin_m, sum_abs_intercept, 'k-')
    ax3.axvline(best_range_m, color='r', linestyle='--', label=f'Best: {best_range_m:.0f}m')
    ax3.set_ylabel('Sum |b|')
    ax3.set_title('Sum of Intercepts')
    ax3.legend()
    
    # Panel 4: R² with best point marked
    ax4 = axes[3]
    im4 = ax4.pcolormesh(range_bin_m, half_length_m, r_squared.T, shading='auto', vmin=0, vmax=1)
    ax4.scatter([best_range_m], [best_half_m], color='r', s=100, marker='x', linewidths=2)
    plt.colorbar(im4, ax=ax4, label='R²')
    ax4.set_xlabel('Center range (m)')
    ax4.set_ylabel('Half-length (m)')
    ax4.set_title('Coefficient of Determination')
    
    fig.suptitle(title)
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig
