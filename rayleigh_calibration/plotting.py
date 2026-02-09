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
    x_data = time_datetime if time_datetime is not None else hours_since_start

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

    # Format x-axis
    if time_datetime is not None:
        from matplotlib import dates as mdates
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        fig.autofmt_xdate(rotation=45, ha="right")
        ax.set_xlabel("Time (UTC)")
    else:
        ax.set_xlabel("Hours since start of window")

    ax.set_ylabel("Range (km)")
    ax.set_ylim(0, min(float(range_alc.max()) * 1e-3, 15.0))
    ax.set_title(title or "Range-Corrected Signal")
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
    cl_matrix : (N, M) Lidar constant for each combination.  NaN = failed.
    cl_median : float   Best-estimate (median).
    cl_uncertainty : float  Robust 2-sigma uncertainty.
    title, save_path : see module docstring.
    """
    plt = _get_plt()

    fig, axes = plt.subplots(
        1, 2, figsize=(14, 5),
        gridspec_kw={"width_ratios": [2, 1]},
    )

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
    ax.set_title("Lidar Constant -- Sensitivity Grid")

    # -- Right: box-plot / strip-chart --
    ax2 = axes[1]
    valid = cl_matrix[np.isfinite(cl_matrix)].ravel()
    if len(valid) > 0:
        ax2.boxplot(
            valid, vert=True, widths=0.4, patch_artist=True,
            boxprops=dict(facecolor="#aec7e8", alpha=0.7),
        )
        rng = np.random.default_rng(42)
        ax2.scatter(
            np.ones_like(valid) + rng.uniform(-0.08, 0.08, len(valid)),
            valid, s=30, alpha=0.6, c="#1f77b4", zorder=3,
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
    ax2.set_title("Distribution")
    ax2.set_xticks([])
    ax2.grid(True, axis="y", alpha=0.25)

    fig.suptitle(title or "Sensitivity Analysis", fontsize=13)
    fig.tight_layout()

    _save_and_close(fig, save_path)
    return fig
