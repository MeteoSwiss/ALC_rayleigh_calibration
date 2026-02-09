"""
Plotting utilities for Rayleigh calibration diagnostics.

This module provides optional visualization functions for diagnosing
calibration results.  Every function:
  * accepts an optional *save_path* -- when given the figure is saved as PNG
    **and closed** automatically so memory is freed during batch runs.
  * returns the ``Figure`` object (useful for interactive inspection).
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, TYPE_CHECKING, List

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib import dates as mdates

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_plt():
    """Lazy import of matplotlib with a non-interactive backend."""
    import matplotlib
    matplotlib.use("Agg")  # headless -- no DISPLAY required
    import matplotlib.pyplot as plt
    return plt


def _save_and_close(fig, save_path: Optional[Path], dpi: int = 150):
    """Save to *save_path* (if given) then close the figure."""
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.info("Plot saved -> %s", save_path)
    plt = _get_plt()
    plt.close(fig)


# =========================================================================
# 1.  RCS time-series
# =========================================================================


def plot_rcs_timeseries(
    hours_since_start: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    rcs: NDArray[np.float64],
    cbh: Optional[NDArray[np.float64]] = None,
    no_cloud_value: float = -9,
    title: str = "",
    save_path: Optional[Path] = None,
    time_datetime: Optional[List[datetime]] = None,
    molecular_window_range: Optional[tuple] = None,
    sensitivity_range: Optional[tuple] = None,
    time_subset_info: Optional[str] = None,
    time_subset_indices: Optional[List[NDArray]] = None,
) -> "Figure":
    """
    Pseudo-colour plot of the range-corrected signal with cloud-base overlay.

    Parameters
    ----------
    hours_since_start : (N,)
        Time axis in hours since the first profile.
    range_alc : (M,)
        Range bins in metres.
    rcs : (N, M)
        Range-corrected signal.
    cbh : (N, L), optional
        Cloud base heights (up to L layers).
    no_cloud_value : float
        Sentinel used for "no cloud".
    time_datetime : (N,), optional
        Datetime objects for better x-axis formatting.
    molecular_window_range : tuple of (start_m, end_m), optional
        Optimal molecular window range to overlay.
    sensitivity_range : tuple of (min_m, max_m), optional
        Full sensitivity test range (with altitude shifts) to overlay.
    time_subset_info : str, optional
        Information about time subsets used in sensitivity analysis.
    time_subset_indices : list of ndarray, optional
        Profile indices for each time subset (for visualization).
    title, save_path : see module docstring.
    """
    plt = _get_plt()

    fig, ax = plt.subplots(figsize=(13, 5))

    # Log-scale colour mesh (mask non-positive values)
    rcs_t = rcs.T.copy().astype(float)
    rcs_t[rcs_t <= 0] = np.nan
    log_rcs = np.log10(rcs_t)

    # Auto colour-limits from 5th / 95th percentile of valid data
    valid = log_rcs[np.isfinite(log_rcs)]
    if len(valid) > 0:
        vmin, vmax = float(np.percentile(valid, 5)), float(np.percentile(valid, 95))
    else:
        vmin, vmax = 0.0, 6.0

    # Use datetime for x-axis if available, otherwise use hours
    use_datetime = time_datetime is not None
    if use_datetime:
        from matplotlib import dates as mdates
        # Convert datetime objects to matplotlib numeric dates
        x_data = np.array(mdates.date2num(time_datetime))
    else:
        x_data = hours_since_start

    im = ax.pcolormesh(
        x_data,
        range_alc * 1e-3,  # km
        log_rcs,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
        cmap="viridis",
    )
    plt.colorbar(im, ax=ax, label=r"log$_{10}$(RCS)", pad=0.02)

    # Visualize time subset coverage at top
    if time_subset_indices is not None:
        y_top = float(range_alc.max()) * 1e-3  # Top of plot in km
        y_spacing = (y_top * 0.04)  # 4% of range for each subset indicator
        colours_subset = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        
        if use_datetime:
            from matplotlib import dates as mdates
            x_min_num = x_data.min()
            x_max_num = x_data.max()
        else:
            x_min_num = x_data.min()
            x_max_num = x_data.max()
        
        for subset_idx, indices in enumerate(time_subset_indices):
            if len(indices) == 0:
                continue
            y_pos = y_top + y_spacing * (subset_idx + 1)
            # Mark profile indices in this subset
            for idx in indices:
                x_val = x_data[idx] if idx < len(x_data) else x_min_num
                ax.plot([x_val, x_val], [y_pos - 0.05, y_pos + 0.05], 
                       color=colours_subset[subset_idx % len(colours_subset)], 
                       linewidth=2, alpha=0.7, zorder=10)

    # Cloud base height overlay (all available layers)
    if cbh is not None:
        n_layers = cbh.shape[1] if cbh.ndim > 1 else 1
        colours = ["white", "silver", "grey"]
        for layer in range(min(n_layers, 3)):
            c = (cbh[:, layer] if cbh.ndim > 1 else cbh).copy().astype(float)
            c[(c == no_cloud_value) | (c <= 0)] = np.nan
            ax.scatter(
                x_data,
                c * 1e-3,
                s=1,
                c=colours[layer % len(colours)],
                alpha=0.6,
                label=f"CBH layer {layer + 1}" if layer == 0 else None,
            )

    # Overlay molecular window and sensitivity ranges
    if sensitivity_range is not None:
        y_sens = np.array(sensitivity_range) * 1e-3  # Convert to km
        ax.axhline(y_sens[0], color='orange', linestyle='--', linewidth=1.5, alpha=0.8,
                   label=f'Sensitivity range ({y_sens[0]:.2f}–{y_sens[1]:.2f} km)')
        ax.axhline(y_sens[1], color='orange', linestyle='--', linewidth=1.5, alpha=0.8)
    
    if molecular_window_range is not None:
        y_mol = np.array(molecular_window_range) * 1e-3  # Convert to km
        ax.axhline(y_mol[0], color='cyan', linestyle='-', linewidth=2.5, 
                   label=f'Molecular window ({y_mol[0]:.2f}–{y_mol[1]:.2f} km)')
        ax.axhline(y_mol[1], color='cyan', linestyle='-', linewidth=2.5)

    # Format x-axis
    if use_datetime:
        from matplotlib import dates as mdates
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        fig.autofmt_xdate(rotation=45, ha="right")
        ax.set_xlabel("Time (UTC)")
    else:
        ax.set_xlabel("Hours since start of window")

    ax.set_ylabel("Range (km)")
    # Extend y-axis if showing subset indicators
    y_max = float(range_alc.max()) * 1e-3
    if time_subset_indices is not None:
        y_max += y_max * 0.25  # Add 25% for subset indicators
    ax.set_ylim(0, y_max)
    ax.set_title(title or "Range-Corrected Signal")
    
    # Add legend if we have overlays (positioned above plot to avoid masking subset indicators)
    if molecular_window_range is not None or sensitivity_range is not None:
        ax.legend(loc='upper center', fontsize=9, framealpha=0.9, 
                 bbox_to_anchor=(0.5, 1.02), ncol=2)
    
    # Add time subset information box
    if time_subset_info is not None:
        ax.text(0.02, 0.8, time_subset_info, transform=ax.transAxes,
                fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    _save_and_close(fig, save_path)
    return fig


# =========================================================================
# 2.  Molecular fit diagnostic (3 panels)
# =========================================================================


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
    Three-panel diagnostic: raw RCS, normalised signal vs molecular theory,
    and attenuated backscatter comparison.
    """
    plt = _get_plt()

    z_km = (range_alc + altitude) * 1e-3  # km ASL
    z_start = fit_altitude_start * 1e-3
    z_end = fit_altitude_end * 1e-3
    z_max = min(z_end + 3.0, float(z_km.max()))

    fig, axes = plt.subplots(1, 3, figsize=(15, 8), sharey=True)

    # -- Panel 1: Raw RCS --
    ax1 = axes[0]
    ax1.plot(rcs_mean, z_km, color="#1f77b4", lw=0.6)
    ax1.axhspan(z_start, z_end, color="gold", alpha=0.15, label="Fit window")
    ax1.axhline(z_start, color="k", ls="--", lw=0.7, alpha=0.5)
    ax1.axhline(z_end, color="k", ls="--", lw=0.7, alpha=0.5)
    ax1.set_xlabel("Raw RCS (a.u.)")
    ax1.set_ylabel("Altitude (km ASL)")
    ax1.set_title("Range-Corrected Signal")
    ax1.set_ylim(0, z_max)
    ax1.grid(True, alpha=0.25)
    ax1.legend(fontsize=8)

    # -- Panel 2: Normalised signal vs molecular --
    ax2 = axes[1]
    scale = 1e6
    ax2.plot(signal_normalized * scale, z_km, color="#1f77b4", lw=0.6,
             label="Normalised signal")
    ax2.plot(p_mol * scale, z_km, color="#d62728", lw=0.8, ls="--",
             label="Molecular theory")
    ax2.axhspan(z_start, z_end, color="gold", alpha=0.15)
    ax2.axhline(z_start, color="k", ls="--", lw=0.7, alpha=0.5)
    ax2.axhline(z_end, color="k", ls="--", lw=0.7, alpha=0.5)
    ax2.set_xscale("log")
    ax2.set_xlabel(r"Signal (Mm$^{-1}$)")
    ax2.set_title("Normalised Signal vs Molecular")
    ax2.set_ylim(0, z_max)
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.25)

    # -- Panel 3: Attenuated backscatter --
    ax3 = axes[2]
    beta_att = signal_normalized * (range_alc ** 2)
    ax3.plot(beta_att * scale, z_km, color="#1f77b4", lw=0.6, label=r"Retrieved $\beta_{att}$")
    ax3.plot(beta_att_mol * scale, z_km, color="#d62728", lw=0.8, ls="--",
             label=r"Molecular $\beta_{att}$")
    ax3.axhspan(z_start, z_end, color="gold", alpha=0.15)
    ax3.axhline(z_start, color="k", ls="--", lw=0.7, alpha=0.5)
    ax3.axhline(z_end, color="k", ls="--", lw=0.7, alpha=0.5)
    ax3.set_xlabel(r"$\beta_{att}$ (Mm$^{-1}$ sr$^{-1}$)")
    ax3.set_title("Attenuated Backscatter")
    ax3.set_ylim(0, z_max)
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.25)

    fig.suptitle(title or "Molecular Fit Diagnostic", fontsize=13)
    fig.tight_layout()

    _save_and_close(fig, save_path)
    return fig


# =========================================================================
# 3.  Lidar-constant profile
# =========================================================================


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
    Altitude profile of the lidar constant with median, slope estimate,
    and uncertainty band.
    """
    plt = _get_plt()

    fig, ax = plt.subplots(figsize=(8, 10))

    r_km = range_alc * 1e-3

    # CL(z)
    ax.plot(cl_profile, r_km, color="#1f77b4", lw=0.6, label="CL(z)")

    # Median + uncertainty band
    ax.axvline(cl_median, color="#d62728", ls="-", lw=1.2,
               label=f"Median: {cl_median:.3e}")
    ax.axvspan(
        cl_median - cl_uncertainty,
        cl_median + cl_uncertainty,
        color="#d62728", alpha=0.10,
        label=f"\u00b12\u03c3: {cl_uncertainty:.2e}",
    )

    # Slope method
    ax.axvline(cl_slope, color="#2ca02c", ls="--", lw=1.0,
               label=f"Slope: {cl_slope:.3e}")

    # Molecular window
    ax.axhline(fit_range_start * 1e-3, color="k", ls="--", lw=0.7, alpha=0.5)
    ax.axhline(fit_range_end * 1e-3, color="k", ls="--", lw=0.7, alpha=0.5)
    ax.axhspan(fit_range_start * 1e-3, fit_range_end * 1e-3,
               color="gold", alpha=0.15, label="Molecular window")

    ax.set_xlabel("Lidar Constant (a.u.)")
    ax.set_ylabel("Range (km)")
    ax.set_title(title or "Lidar Constant Profile")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25)

    # Sensible x-limits
    margin = max(cl_uncertainty * 4, cl_median * 0.3)
    ax.set_xlim(cl_median - margin, cl_median + margin)
    ax.set_ylim(0, fit_range_end * 1.5e-3)

    fig.tight_layout()
    _save_and_close(fig, save_path)
    return fig


# =========================================================================
# 4.  Window search diagnostics
# =========================================================================


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
    """Four-panel diagnostic of the Rayleigh window grid search."""
    plt = _get_plt()

    fig, axes = plt.subplots(4, 1, figsize=(13, 15), sharex=True)

    r_km = range_bin_m * 1e-3
    h_km = half_length_m * 1e-3

    # Panel 1 -- Slopes
    ax1 = axes[0]
    im1 = ax1.pcolormesh(r_km, h_km, slopes.T, shading="auto", cmap="RdBu_r")
    plt.colorbar(im1, ax=ax1, label="Slope (a)", pad=0.02)
    ax1.set_ylabel("Half-length (km)")
    ax1.set_title("Fit Slope")
    ax1.grid(True, alpha=0.15)

    # Panel 2 -- |Intercept|
    ax2 = axes[1]
    im2 = ax2.pcolormesh(r_km, h_km, np.abs(intercepts.T), shading="auto", cmap="magma")
    plt.colorbar(im2, ax=ax2, label="|Intercept| (b)", pad=0.02)
    ax2.set_ylabel("Half-length (km)")
    ax2.set_title("Fit Intercept (absolute)")
    ax2.grid(True, alpha=0.15)

    # Panel 3 -- sum |b| vs centre
    ax3 = axes[2]
    ax3.plot(r_km, sum_abs_intercept, "k-", lw=0.8)
    ax3.axvline(best_range_m * 1e-3, color="#d62728", ls="--", lw=1.0,
                label=f"Best centre: {best_range_m / 1e3:.2f} km")
    ax3.set_ylabel(r"$\Sigma$|b|")
    ax3.set_title("Sum of |Intercepts| (lower is better)")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.25)

    # Panel 4 -- R²
    ax4 = axes[3]
    im4 = ax4.pcolormesh(r_km, h_km, r_squared.T, shading="auto",
                         vmin=0, vmax=1, cmap="YlGnBu")
    ax4.scatter(
        [best_range_m * 1e-3], [best_half_m * 1e-3],
        color="#d62728", s=120, marker="x", linewidths=2.5, zorder=5,
        label=f"Optimum",
    )
    plt.colorbar(im4, ax=ax4, label=r"R$^2$", pad=0.02)
    ax4.set_xlabel("Centre range (km)")
    ax4.set_ylabel("Half-length (km)")
    ax4.set_title("Coefficient of Determination")
    ax4.legend(fontsize=9, loc="upper left")
    ax4.grid(True, alpha=0.15)

    fig.suptitle(title or "Rayleigh Window Search", fontsize=13)
    fig.tight_layout()

    _save_and_close(fig, save_path)
    return fig


# =========================================================================
# 5.  Sensitivity analysis heatmap (NEW)
# =========================================================================


def plot_sensitivity_analysis(
    lr_values: NDArray[np.float64],
    alt_shifts: NDArray[np.float64],
    cl_matrix: NDArray[np.float64],
    cl_median: float,
    cl_uncertainty: float,
    cl_cube: Optional[NDArray[np.float64]] = None,
    time_sample_labels: Optional[List[str]] = None,
    title: str = "",
    save_path: Optional[Path] = None,
) -> "Figure":
    """
    Heatmap of lidar-constant values across the LR x altitude-shift grid,
    with a strip-chart / box-plot summary on the right.

    Parameters
    ----------
    lr_values : (N,)   Lidar-ratio values tested.
    alt_shifts : (M,)  Altitude-window shifts in metres.
    cl_matrix : (N, M) Lidar constant for each (LR, alt) combination
                        (median over time samples).  NaN = failed.
    cl_median : float   Best-estimate (median over all).
    cl_uncertainty : float  Robust 2-sigma uncertainty.
    cl_cube : (N, M, T), optional
        Full 3-D cube of lidar constants (LR × alt-shift × time-sample).
        When provided the right panel shows per-sample spread.
    time_sample_labels : list of str, optional
        Labels for each time sample (length T).
    title, save_path : see module docstring.
    """
    plt = _get_plt()

    n_panels = 3 if cl_cube is not None else 2
    fig, axes = plt.subplots(
        1, n_panels, figsize=(6 * n_panels, 5),
        gridspec_kw={"width_ratios": [2, 1, 1][:n_panels]},
    )
    if n_panels == 2:
        axes = list(axes)

    # -- Left: heatmap (deviation from median %) --
    ax = axes[0]
    rel_dev = (cl_matrix - cl_median) / cl_median * 100

    vabs = max(float(np.nanmax(np.abs(rel_dev))), 1.0)
    im = ax.imshow(
        rel_dev, aspect="auto", origin="lower",
        cmap="RdBu_r", vmin=-vabs, vmax=vabs,
    )
    plt.colorbar(im, ax=ax, label="Deviation from median (%)", pad=0.02)

    # Annotate cells with absolute CL
    for i in range(len(lr_values)):
        for j in range(len(alt_shifts)):
            val = cl_matrix[i, j]
            if np.isfinite(val):
                ax.text(
                    j, i, f"{val:.2e}",
                    ha="center", va="center", fontsize=7,
                    color="black" if abs(rel_dev[i, j]) < vabs * 0.6 else "white",
                )

    ax.set_xticks(range(len(alt_shifts)))
    ax.set_xticklabels([f"{s:+.0f}" for s in alt_shifts])
    ax.set_yticks(range(len(lr_values)))
    ax.set_yticklabels([f"{lr:.0f}" for lr in lr_values])
    ax.set_xlabel("Altitude-window shift (m)")
    ax.set_ylabel("Lidar ratio (sr)")
    ax.set_title("Sensitivity Grid (median over time samples)")

    # -- Middle: overall box-plot / strip-chart --
    ax2 = axes[1]
    if cl_cube is not None:
        all_valid = cl_cube[np.isfinite(cl_cube)].ravel()
    else:
        all_valid = cl_matrix[np.isfinite(cl_matrix)].ravel()

    if len(all_valid) > 0:
        ax2.boxplot(
            all_valid, vert=True, widths=0.4, patch_artist=True,
            boxprops=dict(facecolor="#aec7e8", alpha=0.7),
        )
        rng = np.random.default_rng(42)
        ax2.scatter(
            np.ones_like(all_valid) + rng.uniform(-0.08, 0.08, len(all_valid)),
            all_valid, s=20, alpha=0.4, c="#1f77b4", zorder=3,
        )
        ax2.axhline(cl_median, color="#d62728", ls="-", lw=1.2,
                     label=f"Median: {cl_median:.3e}")
        ax2.axhspan(
            cl_median - cl_uncertainty, cl_median + cl_uncertainty,
            color="#d62728", alpha=0.10,
            label=f"\u00b12\u03c3: {cl_uncertainty:.2e}",
        )
        ax2.legend(fontsize=8, loc="upper right")

    ax2.set_ylabel("Lidar Constant")
    ax2.set_title("All Combinations")
    ax2.set_xticks([])
    ax2.grid(True, axis="y", alpha=0.25)

    # -- Right: per time-sample comparison (only when cube is provided) --
    if cl_cube is not None and n_panels == 3:
        ax3 = axes[2]
        n_t = cl_cube.shape[2]
        labels = time_sample_labels or [str(i) for i in range(n_t)]
        # Collect valid values per time sample
        sample_data = []
        for t in range(n_t):
            vals = cl_cube[:, :, t].ravel()
            sample_data.append(vals[np.isfinite(vals)])

        positions = list(range(1, n_t + 1))
        bp = ax3.boxplot(
            sample_data, positions=positions, vert=True, widths=0.5,
            patch_artist=True,
        )
        colours = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                    "#8c564b", "#e377c2", "#7f7f7f"]
        for i, patch in enumerate(bp["boxes"]):
            patch.set_facecolor(colours[i % len(colours)])
            patch.set_alpha(0.6)

        ax3.axhline(cl_median, color="#d62728", ls="--", lw=1.0,
                     label=f"Median: {cl_median:.3e}")
        ax3.axhspan(
            cl_median - cl_uncertainty, cl_median + cl_uncertainty,
            color="#d62728", alpha=0.08,
        )
        ax3.set_xticks(positions)
        ax3.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax3.set_xlabel("Time sample")
        ax3.set_ylabel("Lidar Constant")
        ax3.set_title("Time-Sample Spread")
        ax3.legend(fontsize=8, loc="upper right")
        ax3.grid(True, axis="y", alpha=0.25)

    fig.suptitle(title or "Sensitivity Analysis", fontsize=13)
    fig.tight_layout()

    _save_and_close(fig, save_path)
    return fig
