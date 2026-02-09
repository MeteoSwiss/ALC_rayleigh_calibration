"""
Main Rayleigh calibration engine.

This module provides the high-level calibration function that orchestrates
the entire calibration process.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, List
import time as timing

import numpy as np
from numpy.typing import NDArray

from .config import InstrumentInfo, CalibrationOptions, CalibrationResult
from .data_loader import (
    build_file_paths,
    load_l1_data,
    filter_time_range,
    filter_cloudy_profiles,
    CeilometerData,
)
from .atmosphere import (
    load_standard_atmosphere,
    load_ecmwf_profile,
    calculate_molecular_properties,
    klett_inversion,
    MOLECULAR_LIDAR_RATIO,
)
from .rayleigh_fit import (
    find_optimal_molecular_window,
    calculate_lidar_constant,
    validate_calibration,
    RayleighFitResult,
)
from .output import write_calibration_result
from .plotting import (
    plot_rcs_timeseries,
    plot_molecular_fit,
    plot_lidar_constant,
    plot_rayleigh_window_search,
    plot_sensitivity_analysis,
)


# Configure module logger
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Sensitivity perturbations
# ---------------------------------------------------------------------------
LR_DELTAS = (-20, -10, 0, +10, +20)        # Lidar-ratio perturbations (sr)
ALT_SHIFTS_M = (-200, -100, 0, +100, +200)  # Altitude-window shifts (m)


def _plot_dir(options: CalibrationOptions, info: InstrumentInfo, date_str: str) -> Path:
    """Return (and create) the plot output directory for one station/date."""
    d = options.folder_output / "plots" / info.wmo_id / date_str[:4]
    d.mkdir(parents=True, exist_ok=True)
    return d


def _plot_tag(info: InstrumentInfo, date_str: str) -> str:
    """File-name prefix for a plot, e.g. '20260101_0-20000-0-06610'."""
    return f"{date_str}_{info.wmo_id}"


@dataclass
class _PerturbationResult:
    """Lidar constant obtained for one (LR, altitude-shift) combination."""
    lidar_ratio: float
    altitude_shift_m: float
    lidar_constant: float
    uncertainty_single: float  # single-config uncertainty (from CL profile scatter)
    # Optional diagnostic arrays (only stored for nominal run, used for plots)
    cl_profile: Optional[NDArray[np.float64]] = None
    signal_normalized: Optional[NDArray[np.float64]] = None
    beta_att: Optional[NDArray[np.float64]] = None
    beta_tot: Optional[NDArray[np.float64]] = None
    ext_tot: Optional[NDArray[np.float64]] = None


def _compute_cl_for_perturbation(
    rcs_mean: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    beta_mol: NDArray[np.float64],
    fit_result: RayleighFitResult,
    lidar_ratio_aerosol: float,
    altitude_shift_m: float,
    subtract_background: bool,
    consider_points_lower_than_molecular: bool,
    return_diagnostics: bool = False,
) -> Optional[_PerturbationResult]:
    """
    Run Klett inversion + CL calculation for ONE (LR, altitude-shift) pair.

    Returns None if the shifted window is out of range or the inversion fails.
    """
    # --- Shifted molecular window ----------------------------------------
    shifted_start = fit_result.range_start_m + altitude_shift_m
    shifted_end = fit_result.range_end_m + altitude_shift_m

    if shifted_start < range_alc[0] or shifted_end > range_alc[-1]:
        return None

    mol_mask = (
        (range_alc >= shifted_start) &
        (range_alc <= shifted_end) &
        ~np.isnan(rcs_mean)
    )
    mol_indices = np.where(mol_mask)[0]
    if len(mol_indices) < 3:
        return None

    i_start_mol = mol_indices[0]
    i_end_mol = mol_indices[-1]

    # --- Normalize signal ------------------------------------------------
    signal = rcs_mean / (range_alc ** 2)
    if subtract_background:
        signal_normalized = (signal - fit_result.intercept) / fit_result.slope
    else:
        signal_normalized = signal / fit_result.slope

    beta_att = signal_normalized * (range_alc ** 2)

    # --- Reference value in the shifted window ---------------------------
    reference_idx = int((i_start_mol + i_end_mol) / 2)
    ref_vals = beta_att[mol_mask] / beta_mol[mol_mask]
    if len(ref_vals) == 0 or np.all(np.isnan(ref_vals)):
        return None
    reference_value = np.nanmean(ref_vals)

    # --- Klett inversion -------------------------------------------------
    dz = np.abs(range_alc[1] - range_alc[0])
    i_start_ext = 0  # from ground
    i_end_ext = i_end_mol

    try:
        beta_aer, beta_tot, ext_aer = klett_inversion(
            beta_att=beta_att,
            beta_mol=beta_mol,
            range_alc=range_alc,
            reference_index=reference_idx,
            lidar_ratio_aerosol=lidar_ratio_aerosol,
            reference_value=reference_value,
            i_start=i_start_ext,
            i_end=i_end_ext,
        )
    except Exception:
        return None

    ext_tot = ext_aer + beta_mol * MOLECULAR_LIDAR_RATIO

    if not consider_points_lower_than_molecular:
        below_mol = beta_tot < beta_mol
        beta_aer[below_mol] = 0
        beta_tot[below_mol] = beta_mol[below_mol]
        ext_aer[below_mol] = 0

    # --- Lidar constant --------------------------------------------------
    try:
        cl_result = calculate_lidar_constant(
            rcs_mean=rcs_mean,
            beta_tot=beta_tot,
            ext_tot=ext_tot,
            range_alc=range_alc,
            molecular_mask=mol_mask,
            fit_result=fit_result,
            subtract_background=subtract_background,
        )
    except ValueError:
        return None

    if np.isnan(cl_result.lidar_constant) or cl_result.lidar_constant <= 0:
        return None

    return _PerturbationResult(
        lidar_ratio=lidar_ratio_aerosol,
        altitude_shift_m=altitude_shift_m,
        lidar_constant=cl_result.lidar_constant,
        uncertainty_single=cl_result.uncertainty,
        cl_profile=cl_result.cl_profile if return_diagnostics else None,
        signal_normalized=signal_normalized if return_diagnostics else None,
        beta_att=beta_att if return_diagnostics else None,
        beta_tot=beta_tot if return_diagnostics else None,
        ext_tot=ext_tot if return_diagnostics else None,
    )


def calibrate_rayleigh(
    date_str: str,
    info: InstrumentInfo,
    options: CalibrationOptions,
    std_atm_file: Optional[Path] = None,
) -> CalibrationResult:
    """
    Perform Rayleigh calibration for a single instrument on a single date.

    This is the main entry point for calibration. It orchestrates:
    1. Loading L1 data from NetCDF files
    2. Filtering to nighttime and cloud-free profiles
    3. Loading atmospheric model (standard or ECMWF)
    4. Calculating molecular scattering properties
    5. Finding optimal molecular window via Rayleigh fit
    6. Performing Klett inversion for extinction
    7. Calculating lidar constant
    8. Writing results to NetCDF

    Parameters
    ----------
    date_str : str
        Date string in YYYYMMDD format.
    info : InstrumentInfo
        Instrument configuration.
    options : CalibrationOptions
        Calibration options.
    std_atm_file : Path, optional
        Path to standard atmosphere file.

    Returns
    -------
    CalibrationResult
        Calibration result including lidar constant, flag, and uncertainty.
    """
    start_time = timing.time()

    # Default standard atmosphere file
    if std_atm_file is None:
        std_atm_file = Path("standard_atmosphere_US_1976_50km.csv")

    logger.info(f"Starting Rayleigh calibration for {info.site_name} on {date_str}")

    # =========================================================================
    # Step 1: Load L1 data
    # =========================================================================
    file_current, file_previous = build_file_paths(date_str, info, options)

    # Build file list from existing files
    file_list = []
    if file_previous.exists():
        file_list.append(file_previous)
    if file_current.exists():
        file_list.append(file_current)

    if not file_list:
        logger.warning(f"No data files found for {date_str}")
        return CalibrationResult(
            lidar_constant=-1,
            flag=0,
            uncertainty=0,
            message="No data files found",
        )

    logger.info(f"Loading {len(file_list)} files")

    data = load_l1_data(file_list, info.instrument_type)
    if data is None:
        logger.warning("Failed to load data")
        return CalibrationResult(
            lidar_constant=-1,
            flag=0,
            uncertainty=0,
            message="Failed to load data",
        )

    logger.info(f"Loaded {len(data.time)} profiles")

    # =========================================================================
    # Step 2: Filter to nighttime window (solar time)
    # =========================================================================
    solar_offset = data.longitude / 15.0
    logger.info(
        f"Station lon={data.longitude:.2f}° → solar time offset = "
        f"{solar_offset:+.2f} h (solar midnight ≈ "
        f"{(-solar_offset) % 24:05.2f} UTC)"
    )
    data = filter_time_range(data, date_str, options)
    logger.info(f"After time filtering: {len(data.time)} profiles")

    if len(data.time) == 0:
        logger.warning("No profiles in nighttime window")
        return CalibrationResult(
            lidar_constant=-1,
            flag=0,
            uncertainty=0,
            message="No profiles in nighttime window",
        )

    # ── Plot: RCS time-series (before cloud filtering) ──
    if options.plot_main or options.plot_all:
        pdir = _plot_dir(options, info, date_str)
        tag = _plot_tag(info, date_str)
        try:
            plot_rcs_timeseries(
                hours_since_start=data.hours_since_start,
                range_alc=data.range_alc,
                rcs=data.rcs,
                cbh=data.cbh,
                no_cloud_value=info.instrument_type.no_cloud_value,
                title=f"{info.site_name} ({info.wmo_id}) — {date_str} — RCS",
                save_path=pdir / f"{tag}_rcs_timeseries.png",
                time_datetime=data.time_datetime,
            )
        except Exception as exc:
            logger.warning(f"plot_rcs_timeseries failed: {exc}")

    # =========================================================================
    # Step 3: Filter cloudy profiles
    # =========================================================================
    no_cloud_value = info.instrument_type.no_cloud_value
    data, is_clear, is_partial = filter_cloudy_profiles(data, options, no_cloud_value)

    if not is_clear:
        logger.warning("Not a clear night")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-1,
            uncertainty=0,
            message="Not a clear night",
        )

    logger.info(f"After cloud filtering: {len(data.time)} profiles (partial: {is_partial})")

    # =========================================================================
    # Step 4: Load atmospheric model
    # =========================================================================
    logger.info("Loading atmospheric model")

    # Ensure altitude grid is not masked
    altitude_grid = data.altitude_grid
    if np.ma.isMaskedArray(altitude_grid):
        altitude_grid = altitude_grid.data

    if options.use_std_atm:
        atm_profile = load_standard_atmosphere(std_atm_file, altitude_grid)
    else:
        ecmwf_file = options.folder_ecmwf / f"MACC_{date_str}.nc"
        atm_profile = load_ecmwf_profile(
            ecmwf_file, info.latitude, info.longitude, altitude_grid
        )
        if atm_profile is None:
            logger.warning(f"ECMWF file not found: {ecmwf_file}")
            return CalibrationResult(
                lidar_constant=-1,
                flag=-4,
                uncertainty=0,
                message="Missing model data",
            )

    # =========================================================================
    # Step 5: Calculate molecular properties
    # =========================================================================
    logger.info("Calculating molecular scattering properties")

    wavelength_m = info.instrument_type.wavelength_nm * 1e-9

    mol_props = calculate_molecular_properties(
        atm_profile.temperature,
        atm_profile.pressure,
        data.range_alc,
        wavelength_m,
    )

    logger.info(f"Time elapsed: {timing.time() - start_time:.1f}s")

    # =========================================================================
    # Step 6: Find optimal molecular window
    # =========================================================================
    logger.info("Finding optimal molecular window")

    # Calculate mean RCS profile
    rcs_mean = np.nanmean(data.rcs, axis=0)

    if np.all(np.isnan(rcs_mean)):
        logger.warning("RCS contains only NaN")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-5,
            uncertainty=0,
            message="RCS contains only NaN",
        )

    # Range-normalized signal
    signal = rcs_mean / (data.range_alc ** 2)

    fit_result = find_optimal_molecular_window(
        signal=signal,
        p_mol=mol_props.p_mol,
        range_alc=data.range_alc,
        half_length_options_m=options.half_length_options_m,
        range_start_m=options.range_start_m,
        range_end_m=options.range_end_m,
        increment_bins=options.fit_range_increment_bins,
    )

    # Update altitude values
    fit_result.altitude_start = fit_result.range_start_m + data.altitude
    fit_result.altitude_end = fit_result.range_end_m + data.altitude

    logger.info(
        f"Optimal window: {fit_result.range_start_m:.0f}-{fit_result.range_end_m:.0f}m "
        f"(R²={fit_result.r_squared:.4f})"
    )

    # ── Plot: Rayleigh window search diagnostics ──
    if options.plot_all and fit_result.search_diagnostics is not None:
        pdir = _plot_dir(options, info, date_str)
        tag = _plot_tag(info, date_str)
        diag = fit_result.search_diagnostics
        try:
            plot_rayleigh_window_search(
                range_bin_m=diag.range_bin_m,
                half_length_m=diag.half_length_m,
                slopes=diag.slopes,
                intercepts=diag.intercepts,
                r_squared=diag.r_squared,
                sum_abs_intercept=diag.sum_abs_intercept,
                best_range_m=fit_result.center_range_m,
                best_half_m=fit_result.half_length_m,
                title=f"{info.site_name} — {date_str} — Window Search",
                save_path=pdir / f"{tag}_window_search.png",
            )
        except Exception as exc:
            logger.warning(f"plot_rayleigh_window_search failed: {exc}")

    # =========================================================================
    # Step 7: Validate Rayleigh fit
    # =========================================================================
    if not fit_result.is_valid:
        if fit_result.slope <= 0:
            logger.warning("Negative Rayleigh fit slope")
            return CalibrationResult(
                lidar_constant=-1,
                flag=-7,
                uncertainty=0,
                message="Negative Rayleigh fit",
            )
        else:
            logger.warning("Rayleigh fit issue: |b| > a")
            return CalibrationResult(
                lidar_constant=-1,
                flag=-8,
                uncertainty=0,
                message="Rayleigh fit issue: |b| > a",
            )

    if fit_result.relative_error > options.threshold_quality:
        logger.warning(f"Poor Rayleigh fit quality: {fit_result.relative_error:.1f}%")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-2,
            uncertainty=0,
            message="Signal not proportional to molecular scattering",
        )

    # =========================================================================
    # Steps 8-10: Sensitivity analysis over lidar-ratio and altitude window
    # =========================================================================
    # Run Klett inversion + CL calculation for 25 combinations of
    #   LR ∈ {LRaer-20, LRaer-10, LRaer, LRaer+10, LRaer+20}
    #   altitude shift ∈ {-200m, -100m, 0m, +100m, +200m}
    #
    # The best-estimate CL is the median of the 25 values, and the
    # uncertainty is derived from their spread (half IQR × 2 ≈ robust 2σ).
    # =========================================================================
    logger.info(
        "Sensitivity analysis: 5 LR × 5 altitude shifts = 25 combinations"
    )

    perturbation_results: List[_PerturbationResult] = []

    for lr_delta in LR_DELTAS:
        lr = options.lidar_ratio_aerosol + lr_delta
        if lr <= 0:
            continue  # skip non-physical values
        for alt_shift in ALT_SHIFTS_M:
            pr = _compute_cl_for_perturbation(
                rcs_mean=rcs_mean,
                range_alc=data.range_alc,
                beta_mol=mol_props.beta_mol,
                fit_result=fit_result,
                lidar_ratio_aerosol=lr,
                altitude_shift_m=alt_shift,
                subtract_background=options.subtract_background,
                consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
            )
            if pr is not None:
                perturbation_results.append(pr)

    n_ok = len(perturbation_results)
    logger.info(f"Successful perturbation runs: {n_ok}/25")

    if n_ok == 0:
        logger.warning("All perturbation runs failed")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-5,
            uncertainty=0,
            message="All perturbation runs failed",
        )

    cl_values = np.array([pr.lidar_constant for pr in perturbation_results])

    # Robust statistics
    cl_median = float(np.median(cl_values))
    q25, q75 = float(np.percentile(cl_values, 25)), float(np.percentile(cl_values, 75))
    iqr = q75 - q25
    mad = float(np.median(np.abs(cl_values - cl_median)))  # median absolute deviation

    # Robust 2-sigma uncertainty: max of (half-IQR scaled to ~2σ) and (MAD scaled to ~2σ)
    # For a Gaussian: IQR ≈ 1.35σ, MAD ≈ 0.6745σ
    sigma_iqr = iqr / 1.349
    sigma_mad = mad / 0.6745
    robust_sigma = max(sigma_iqr, sigma_mad)
    uncertainty = 2 * robust_sigma  # ~95% confidence interval

    logger.info(
        f"Sensitivity CL: median={cl_median:.4e}, "
        f"IQR=[{q25:.4e}, {q75:.4e}], MAD={mad:.4e}, "
        f"2σ_robust={uncertainty:.4e} ({uncertainty/cl_median*100:.1f}%)"
    )

    # Also compute the nominal (LR=LRaer, shift=0) result for validation
    nominal = [
        pr for pr in perturbation_results
        if pr.lidar_ratio == options.lidar_ratio_aerosol and pr.altitude_shift_m == 0
    ]
    if nominal:
        cl_nominal = nominal[0].lidar_constant
    else:
        cl_nominal = cl_median  # fallback

    # =========================================================================
    # Validation: slope-method comparison (using nominal config)
    # =========================================================================
    # Recompute the nominal ext_tot for transmission at reference
    pr_nominal = _compute_cl_for_perturbation(
        rcs_mean=rcs_mean,
        range_alc=data.range_alc,
        beta_mol=mol_props.beta_mol,
        fit_result=fit_result,
        lidar_ratio_aerosol=options.lidar_ratio_aerosol,
        altitude_shift_m=0,
        subtract_background=options.subtract_background,
        consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
        return_diagnostics=True,   # keep cl_profile etc. for plots
    )

    # Validation using the slope method
    mol_mask_nominal = (
        (data.range_alc >= fit_result.range_start_m) &
        (data.range_alc <= fit_result.range_end_m) &
        ~np.isnan(rcs_mean)
    )
    i_start_mol_nominal = np.where(mol_mask_nominal)[0][0] if np.any(mol_mask_nominal) else 0

    # Slope-based CL estimate for cross-check
    # optical depth from ground to start of molecular window
    signal_check = rcs_mean / (data.range_alc ** 2)
    if options.subtract_background:
        signal_norm_check = (signal_check - fit_result.intercept) / fit_result.slope
    else:
        signal_norm_check = signal_check / fit_result.slope
    beta_att_check = signal_norm_check * (data.range_alc ** 2)

    # Use nominal Klett for transmission estimate
    ref_idx_check = int((np.where(mol_mask_nominal)[0][0] + np.where(mol_mask_nominal)[0][-1]) / 2)
    ref_val_check = np.nanmean(beta_att_check[mol_mask_nominal] / mol_props.beta_mol[mol_mask_nominal])
    _, _, ext_aer_check = klett_inversion(
        beta_att=beta_att_check,
        beta_mol=mol_props.beta_mol,
        range_alc=data.range_alc,
        reference_index=ref_idx_check,
        lidar_ratio_aerosol=options.lidar_ratio_aerosol,
        reference_value=ref_val_check,
        i_start=0,
        i_end=np.where(mol_mask_nominal)[0][-1],
    )
    ext_tot_check = ext_aer_check + mol_props.beta_mol * MOLECULAR_LIDAR_RATIO
    optical_depth_ref = np.trapz(ext_tot_check[:i_start_mol_nominal], data.range_alc[:i_start_mol_nominal])
    inv_trans_ref = np.exp(2 * optical_depth_ref)

    cl_slope = fit_result.slope * inv_trans_ref
    error_pct = abs((cl_slope - cl_median) / cl_median * 100)

    if error_pct > options.threshold_quality:
        logger.warning(f"Method disagreement: {error_pct:.1f}% > {options.threshold_quality}%")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-3,
            uncertainty=0,
            message=f"Method disagreement: {error_pct:.1f}%",
        )

    if uncertainty > cl_median:
        logger.warning(f"Uncertainty exceeds CL value: {uncertainty:.2e} > {cl_median:.2e}")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-6,
            uncertainty=0,
            message=f"Uncertainty exceeds value: {uncertainty:.2e} > {cl_median:.2e}",
        )

    # ── Plot: Molecular fit, CL profile, Sensitivity analysis ──
    if options.plot_main or options.plot_all:
        pdir = _plot_dir(options, info, date_str)
        tag = _plot_tag(info, date_str)
        plot_title_base = f"{info.site_name} ({info.wmo_id}) — {date_str}"

        # (a) Molecular fit diagnostic
        if pr_nominal is not None and pr_nominal.signal_normalized is not None:
            try:
                # Molecular attenuated backscatter = beta_mol × range²
                beta_att_mol = mol_props.beta_mol * (data.range_alc ** 2)
                plot_molecular_fit(
                    range_alc=data.range_alc,
                    altitude=data.altitude,
                    rcs_mean=rcs_mean,
                    signal_normalized=pr_nominal.signal_normalized,
                    p_mol=mol_props.p_mol,
                    beta_att_mol=beta_att_mol,
                    fit_altitude_start=fit_result.altitude_start,
                    fit_altitude_end=fit_result.altitude_end,
                    title=f"{plot_title_base} — Molecular Fit",
                    save_path=pdir / f"{tag}_molecular_fit.png",
                )
            except Exception as exc:
                logger.warning(f"plot_molecular_fit failed: {exc}")

        # (b) Lidar-constant profile
        if pr_nominal is not None and pr_nominal.cl_profile is not None:
            try:
                plot_lidar_constant(
                    range_alc=data.range_alc,
                    cl_profile=pr_nominal.cl_profile,
                    cl_median=cl_median,
                    cl_slope=cl_slope,
                    cl_uncertainty=uncertainty,
                    fit_range_start=fit_result.range_start_m,
                    fit_range_end=fit_result.range_end_m,
                    title=f"{plot_title_base} — Lidar Constant",
                    save_path=pdir / f"{tag}_lidar_constant.png",
                )
            except Exception as exc:
                logger.warning(f"plot_lidar_constant failed: {exc}")

        # (c) Sensitivity analysis heatmap
        try:
            # Build LR × altitude-shift matrix
            lr_used = np.array(sorted({pr.lidar_ratio for pr in perturbation_results}))
            shifts_used = np.array(sorted({pr.altitude_shift_m for pr in perturbation_results}))
            cl_matrix = np.full((len(lr_used), len(shifts_used)), np.nan)
            lr_idx = {v: i for i, v in enumerate(lr_used)}
            sh_idx = {v: i for i, v in enumerate(shifts_used)}
            for pr in perturbation_results:
                cl_matrix[lr_idx[pr.lidar_ratio], sh_idx[pr.altitude_shift_m]] = pr.lidar_constant

            plot_sensitivity_analysis(
                lr_values=lr_used,
                alt_shifts=shifts_used,
                cl_matrix=cl_matrix,
                cl_median=cl_median,
                cl_uncertainty=uncertainty,
                title=f"{plot_title_base} — Sensitivity",
                save_path=pdir / f"{tag}_sensitivity.png",
            )
        except Exception as exc:
            logger.warning(f"plot_sensitivity_analysis failed: {exc}")

    # =========================================================================
    # Step 11: Build result and write to NetCDF
    # =========================================================================
    flag = 0.5 if is_partial else 1

    result = CalibrationResult(
        lidar_constant=cl_median,
        flag=flag,
        uncertainty=uncertainty,
        calibration_bottom_height=fit_result.altitude_start,
        calibration_top_height=fit_result.altitude_end,
        message=f"Successful (method err: {error_pct:.1f}%, sensitivity 2σ: {uncertainty/cl_median*100:.1f}%)",
    )

    # Prepare housekeeping data
    housekeeping = {
        'laser_life_time': np.nanmean(data.laser_life_time),
        'status_detector': np.nanmean(data.status_detector),
        'status_laser': np.nanmean(data.status_laser),
        'temperature_optical_module': np.nanmean(data.temperature_optical_module),
        'window_transmission': np.nanmean(data.window_transmission),
        'optical_module_id': data.optical_module_id,
    }

    # Get time values
    date_epoch = np.floor(np.max(data.time))  # Central date
    time_start = np.min(data.time)
    time_end = np.max(data.time)

    # Write to NetCDF
    output_path = write_calibration_result(
        output_dir=options.folder_output,
        info=info,
        result=result,
        date_epoch=date_epoch,
        time_start=time_start,
        time_end=time_end,
        wavelength_nm=info.instrument_type.wavelength_nm,
        housekeeping=housekeeping,
    )

    logger.info(f"Results written to {output_path}")
    logger.info(f"Lidar constant: {result.lidar_constant:.4e} ± {result.uncertainty:.4e}")
    logger.info(f"Total time: {timing.time() - start_time:.1f}s")

    return result
