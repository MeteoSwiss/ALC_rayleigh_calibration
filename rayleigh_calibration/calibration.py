"""
Main Rayleigh calibration engine.

This module provides the high-level calibration function that orchestrates
the entire calibration process.
"""

from __future__ import annotations

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
import time as timing

import numpy as np

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
)
from .output import write_calibration_result


# Configure module logger
logger = logging.getLogger(__name__)


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
    # Step 8: Normalize signal and calculate beta_att
    # =========================================================================
    if options.subtract_background:
        signal_normalized = (signal - fit_result.intercept) / fit_result.slope
    else:
        signal_normalized = signal / fit_result.slope

    beta_att = signal_normalized * (data.range_alc ** 2)

    # =========================================================================
    # Step 9: Klett inversion for extinction
    # =========================================================================
    logger.info("Performing Klett inversion")

    dz = np.abs(data.range_alc[1] - data.range_alc[0])

    # Define molecular region
    mol_mask = np.logical_and(
        data.range_alc >= fit_result.range_start_m,
        data.range_alc <= fit_result.range_end_m
    ) & ~np.isnan(rcs_mean)

    mol_indices = np.where(mol_mask)[0]
    i_start_mol = mol_indices[0]
    i_end_mol = mol_indices[-1]

    # Reference value at center of molecular region
    reference_idx = int((i_start_mol + i_end_mol) / 2)
    reference_value = np.mean(beta_att[mol_mask] / mol_props.beta_mol[mol_mask])

    # Extinction start
    i_start_ext = np.argmin(np.abs(data.range_alc - options.z_start_ext))

    # Inversion end
    if options.calc_ext_above_molecular:
        i_end_ext = len(data.range_alc)
    else:
        i_end_ext = i_end_mol

    beta_aer, beta_tot, ext_aer = klett_inversion(
        beta_att=beta_att,
        beta_mol=mol_props.beta_mol,
        range_alc=data.range_alc,
        reference_index=reference_idx,
        lidar_ratio_aerosol=options.lidar_ratio_aerosol,
        reference_value=reference_value,
        i_start=i_start_ext,
        i_end=i_end_ext,
    )

    # Total extinction
    ext_tot = ext_aer + mol_props.beta_mol * MOLECULAR_LIDAR_RATIO

    # Consider points lower than molecular
    if not options.consider_points_lower_than_molecular:
        # Set aerosol to zero where total < molecular
        below_mol = beta_tot < mol_props.beta_mol
        beta_aer[below_mol] = 0
        beta_tot[below_mol] = mol_props.beta_mol[below_mol]
        ext_aer[below_mol] = 0

    # =========================================================================
    # Step 10: Calculate lidar constant
    # =========================================================================
    logger.info("Calculating lidar constant")

    cl_result = calculate_lidar_constant(
        rcs_mean=rcs_mean,
        beta_tot=beta_tot,
        ext_tot=ext_tot,
        range_alc=data.range_alc,
        molecular_mask=mol_mask,
        fit_result=fit_result,
        subtract_background=options.subtract_background,
    )

    # Calculate inverse transmission at reference for validation
    optical_depth_ref = np.trapz(ext_tot[:i_start_mol], data.range_alc[:i_start_mol])
    inv_trans_ref = np.exp(2 * optical_depth_ref)

    # Validate
    is_valid, error_pct, message = validate_calibration(
        cl_result=cl_result,
        fit_result=fit_result,
        inv_transmission_ref=inv_trans_ref,
        threshold=options.threshold_quality,
    )

    if not is_valid:
        if "uncertainty" in message.lower():
            flag = -6
        else:
            flag = -3
        logger.warning(f"Calibration validation failed: {message}")
        return CalibrationResult(
            lidar_constant=-1,
            flag=flag,
            uncertainty=0,
            message=message,
        )

    # =========================================================================
    # Step 11: Build result and write to NetCDF
    # =========================================================================
    flag = 0.5 if is_partial else 1

    result = CalibrationResult(
        lidar_constant=cl_result.lidar_constant,
        flag=flag,
        uncertainty=cl_result.uncertainty,
        calibration_bottom_height=fit_result.altitude_start,
        calibration_top_height=fit_result.altitude_end,
        message=f"Successful (error: {error_pct:.1f}%)",
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
