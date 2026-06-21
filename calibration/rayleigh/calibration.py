"""
Main Rayleigh calibration engine.

This module provides the high-level calibration function that orchestrates
the entire calibration process.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, List
import time as timing

import numpy as np
from numpy.typing import NDArray

from ..config import InstrumentInfo, CalibrationOptions, CalibrationResult, DataLevel
from ..io.data_loader import (
    build_file_paths,
    load_l1_data,
    load_data,
    average_ceilometer_data,
    filter_time_range,
    filter_cloudy_profiles,
    CeilometerData,
)
from .atmosphere import (
    DEFAULT_STANDARD_ATMOSPHERE,
    load_standard_atmosphere,
    load_cams_atmosphere,
    calculate_molecular_properties,
    klett_inversion,
    MOLECULAR_LIDAR_RATIO,
)
from ..water_vapor_correction.water_vapor import (
    in_water_vapor_band,
    laser_spectrum_for,
    cams_water_vapor_profile,
    two_way_wv_transmission,
)
from .rayleigh_fit import (
    find_optimal_molecular_window,
    calculate_lidar_constant,
    validate_calibration,
    RayleighFitResult,
)
from ..io.output import write_calibration_result
from ..io.cams import ensure_cams_file
from ..plotting import (
    plot_rcs_timeseries,
    plot_rayleigh_diagnostics_compact,
)


# Configure module logger
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Sensitivity perturbations
# ---------------------------------------------------------------------------
LR_DELTAS = (-20, -10, 0, +10, +20)        # Lidar-ratio perturbations (sr)
ALT_SHIFTS_M = (-200, -100, 0, +100, +200)  # Altitude-window shifts (m)
N_TIME_SAMPLES = 4                          # Number of random time-subsets
TIME_SUBSET_FRACTION = 0.7                  # Fraction of profiles kept per subset


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
    """Lidar constant obtained for one (LR, altitude-shift, time-sample) combination."""
    lidar_ratio: float
    altitude_shift_m: float
    lidar_constant: float
    uncertainty_single: float  # single-config uncertainty (from CL profile scatter)
    time_sample_idx: int = 0   # 0 = all profiles, 1..N = random subset
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
    sign_error_v10: bool = False,
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
            sign_error_v10=sign_error_v10,
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
    fit_inputs_out: Optional[dict] = None,
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
        std_atm_file = DEFAULT_STANDARD_ATMOSPHERE

    logger.info(f"Starting Rayleigh calibration for {info.site_name} on {date_str}")

    # Suitability warning: CL31/CL51 have signal distortion/saturation that makes the molecular
    # Rayleigh fit unreliable (they are normally calibrated by the liquid-cloud method). The fit
    # still RUNS so a whole network can be processed in one pass — warn that the result is only
    # indicative.
    if not info.instrument_type.supports_calibration:
        warnings.warn(
            f"{info.instrument_type.value}: Rayleigh calibration is not well-suited to this "
            "instrument (CL31/CL51 signal distortion/saturation); the result is indicative only.",
            UserWarning)

    # =========================================================================
    # Step 1: Load L1 data
    # =========================================================================
    candidate_files = build_file_paths(date_str, info, options)

    # Keep only files that exist (chronological order preserved)
    file_list = [f for f in candidate_files if f.exists()]

    if not file_list:
        logger.warning(f"No data files found for {date_str}")
        return CalibrationResult(
            lidar_constant=-1,
            flag=0,
            uncertainty=0,
            message="No data files found",
        )

    logger.info(f"Loading {len(file_list)} {options.data_level.value} file(s)")

    data = load_data(file_list, info.instrument_type, options.data_level)
    if data is None:
        logger.warning("Failed to load data")
        return CalibrationResult(
            lidar_constant=-1,
            flag=0,
            uncertainty=0,
            message="Failed to load data",
        )

    logger.info(f"Loaded {len(data.time)} profiles")

    # Optional pre-averaging: reduces the Rayleigh input to coarser blocks before any
    # filtering or fitting, matching the cloud-calibration speedup.
    avg_time_s = getattr(options, "average_time_s", None)
    avg_range_m = getattr(options, "average_range_m", None)
    # Native L1 is on a fine grid (e.g. CHM15k 15 m x 15 s) that the gated molecular methods
    # (v2/earlinet) over-reject even though the signal matches L2's beta_att; bin it to the
    # standard L2 grid (30 m x 300 s) so L1 and L2 calibrate consistently. Only when the level
    # is L1 and no explicit averaging was requested (L2/RAW untouched; a coarser native grid is
    # a no-op). See network_v2_vs_v11_report.md ("L1 vs L2 - the tie on L1 is a native-grid effect").
    if (avg_time_s is None and avg_range_m is None
            and getattr(options, "data_level", None) == DataLevel.L1
            and getattr(options, "l1_bin_to_l2_grid", True)):
        avg_time_s = getattr(options, "l1_grid_time_s", 300.0)
        avg_range_m = getattr(options, "l1_grid_range_m", 30.0)
        logger.info("L1 native grid -> binning to the L2 grid (%.0f s x %.0f m)", avg_time_s, avg_range_m)
    data = average_ceilometer_data(
        data,
        average_time_s=avg_time_s,
        average_range_m=avg_range_m,
    )
    if avg_time_s or avg_range_m:
        logger.info(
            "Averaged Rayleigh input to %s profiles x %s range bins",
            len(data.time), len(data.range_alc),
        )

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

    # ── Plot: RCS time-series (before cloud filtering) — extra detail, gated on plot_all ──
    if options.plot_all:
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

    # Night window (UTC); used for time-windowed model data (the CAMS molecular
    # profile and/or the WV correction below).
    night = np.datetime64(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}")
    t_start = night - np.timedelta64(24 - int(options.hour_min), "h")
    t_end = night + np.timedelta64(int(options.hour_max), "h")

    molecular_source = getattr(options, "molecular_source", "standard")
    if molecular_source == "cams":
        # Build the molecular reference (beta_mol ~ P/T) from the actual CAMS T/p
        # at the site instead of the US Standard 1976 atmosphere. Same CAMS archive
        # as the WV correction; a night without a matching CAMS month cannot use
        # this option and is skipped (consistent with the WV rule).
        cams_mol_file = ensure_cams_file(
            options.cams_folder, date_str,
            auto_download=getattr(options, "auto_download_cams", False),
            scope=getattr(options, "cams_download_scope", "day"),
            log=logger,
        )
        if cams_mol_file is None:
            logger.warning(
                f"No CAMS for {date_str[:6]}; cannot build CAMS molecular profile"
            )
            return CalibrationResult(
                lidar_constant=-1, flag=-4, uncertainty=0,
                message="No CAMS for molecular profile",
            )
        atm_profile = load_cams_atmosphere(
            cams_mol_file, info.latitude, info.longitude, t_start, t_end, altitude_grid
        )
        if atm_profile is None:
            logger.warning(f"CAMS molecular profile unavailable: {cams_mol_file}")
            return CalibrationResult(
                lidar_constant=-1, flag=-4, uncertainty=0,
                message="Missing CAMS molecular data",
            )
        logger.info("Molecular profile from CAMS T/p")
    elif options.use_std_atm:
        atm_profile = load_standard_atmosphere(std_atm_file, altitude_grid)
    else:
        # ECMWF MACC reanalysis path retired (replaced by CAMS, molecular_source='cams',
        # or the US Standard Atmosphere). Fall back to the standard atmosphere.
        atm_profile = load_standard_atmosphere(std_atm_file, altitude_grid)

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

    # --- Water-vapor correction for 910 nm instruments (M. Hervo rule: a 910 nm
    # night is only calibrated when it can be water-vapor corrected). At 905-911 nm
    # H2O absorption attenuates the molecular-reference signal; folding the two-way
    # WV transmission into the molecular model removes the resulting bias in C_L. ---
    nominal_wl_nm = info.instrument_type.wavelength_nm
    if options.apply_wv_correction and in_water_vapor_band(nominal_wl_nm):
        cams_file = ensure_cams_file(
            options.cams_folder, date_str,
            auto_download=getattr(options, "auto_download_cams", False),
            scope=getattr(options, "cams_download_scope", "day"),
            log=logger,
        )
        if cams_file is None:
            logger.warning(f"No CAMS for {date_str[:6]}; skipping WV-required 910 nm night")
            return CalibrationResult(
                lidar_constant=-1, flag=-4, uncertainty=0,
                message="No CAMS for water-vapor correction",
            )
        lam0_nm, fwhm_nm = laser_spectrum_for(info.instrument_type.value, nominal_wl_nm)
        prof = cams_water_vapor_profile(cams_file, info.latitude, info.longitude, t_start, t_end)
        if prof is not None:
            h_wv, n_wv = prof
            wv_alt_grid = info.altitude + data.range_alc        # ASL, aligned to range_alc
            t2_wv = two_way_wv_transmission(
                wv_alt_grid, info.altitude, h_wv, n_wv,
                Path(options.abs_cs_lookup_table), lam0_nm, fwhm_nm,
            )
            # Remove WV absorption from the measured range-corrected signal so the
            # downstream fit / Klett / lidar-constant recover a WV-free CL:
            #   rcs / T2_wv = CL * beta_tot * T2_scattering.
            if t2_wv.shape[0] == data.rcs.shape[1]:
                data.rcs = data.rcs / t2_wv[None, :]
                logger.info(f"WV correction applied to RCS (median T2_wv={np.nanmedian(t2_wv):.3f})")
            else:
                logger.warning("WV transmission length mismatch; correction skipped")

    logger.info(f"Time elapsed: {timing.time() - start_time:.1f}s")

    # =========================================================================
    # Step 6: Find optimal molecular window
    # =========================================================================
    logger.info("Finding optimal molecular window")

    # Screen outlier profiles (residual aerosol/cloud/noise) BEFORE averaging, so the
    # time-mean is robust without the ~1.5x noise penalty of a median profile at the
    # weak-signal molecular altitudes. Metric = median range-normalised signal in the
    # molecular search band; reject profiles > N robust-sigma (MAD) from the median.
    keep_idx = np.arange(data.rcs.shape[0])
    if getattr(options, "screen_profile_outliers", True) and data.rcs.shape[0] >= 10:
        band = (data.range_alc >= options.range_start_m) & (data.range_alc <= options.range_end_m)
        if np.any(band):
            metric = np.nanmedian(data.rcs[:, band] / (data.range_alc[band] ** 2), axis=1)
            finite = np.isfinite(metric)
            if finite.sum() >= 10:
                center = np.median(metric[finite])
                mad = np.median(np.abs(metric[finite] - center))
                thr = getattr(options, "profile_outlier_nmad", 4.0) * 1.4826 * mad
                keep = finite & (np.abs(metric - center) <= max(thr, 1e-30))
                if keep.sum() >= max(3, int(0.5 * data.rcs.shape[0])):
                    keep_idx = np.where(keep)[0]
                    logger.info(f"Outlier screen: kept {keep_idx.size}/{data.rcs.shape[0]} profiles")
    rcs_use = data.rcs[keep_idx]

    # Collapse the kept profiles in time. Mean is efficient for the weak-signal photon
    # noise at 3-6 km; "median" is robust but ~1.5x noisier there (set via options).
    _agg = np.nanmedian if getattr(options, "time_aggregation", "mean") == "median" else np.nanmean
    rcs_mean = _agg(rcs_use, axis=0)

    if np.all(np.isnan(rcs_mean)):
        logger.warning("RCS contains only NaN")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-5,
            uncertainty=0,
            message="RCS contains only NaN",
        )

    # Range-normalized signal (night mean) + per-profile stack. The stack feeds the
    # E-PROF v2 ("optimal") method's temporal-variability aerosol rejection (molecular is
    # steady in time; aerosol fluctuates).
    signal = rcs_mean / (data.range_alc ** 2)
    signal_stack = rcs_use / (data.range_alc[None, :] ** 2)

    # Optional capture hook (used by the method-comparison harness): expose the prepared
    # fit inputs so every method can be evaluated on the identical profile.
    if fit_inputs_out is not None:
        fit_inputs_out.update(
            signal=signal, p_mol=mol_props.p_mol, range_alc=data.range_alc,
            altitude=data.altitude, signal_stack=signal_stack,
            hours=np.asarray(data.hours_since_start)[keep_idx],
        )

    fit_result = find_optimal_molecular_window(
        signal=signal,
        p_mol=mol_props.p_mol,
        range_alc=data.range_alc,
        half_length_options_m=options.half_length_options_m,
        range_start_m=options.range_start_m,
        range_end_m=options.range_end_m,
        increment_bins=options.fit_range_increment_bins,
        min_window_start_m=options.min_window_start_m,
        min_r2=options.min_window_r2,
        max_rel_error=options.max_window_rel_error,
        method=getattr(options, "molecular_method", "eprof_v1.2"),
        signal_stack=signal_stack,
    )

    # Update altitude values
    fit_result.altitude_start = fit_result.range_start_m + data.altitude
    fit_result.altitude_end = fit_result.range_end_m + data.altitude

    logger.info(
        f"Optimal window: {fit_result.range_start_m:.0f}-{fit_result.range_end_m:.0f}m "
        f"(R²={fit_result.r_squared:.4f})"
    )

    # ── Pre-compute time subsets for sensitivity analysis ──
    # Build now so they're available for visualization
    rng = np.random.default_rng(42)
    n_profiles = rcs_use.shape[0]   # screened profile count
    # Cap at n_profiles: on nights with very few profiles, max(3, ...) could exceed the
    # population and rng.choice(replace=False) would raise "larger sample than population".
    n_keep = min(n_profiles, max(3, int(n_profiles * TIME_SUBSET_FRACTION)))

    rcs_mean_samples = [rcs_mean]  # index 0: all (screened) profiles
    time_subset_indices = [keep_idx]  # original indices of the kept profiles (for viz)

    for _ in range(N_TIME_SAMPLES):
        idx = rng.choice(n_profiles, size=n_keep, replace=False)
        idx = np.sort(idx)  # Sort for visualization
        rcs_mean_samples.append(_agg(rcs_use[idx], axis=0))
        time_subset_indices.append(keep_idx[idx])   # map back to original indices

    # ── Plot: RCS with molecular window annotation — extra detail, gated on plot_all ──
    if options.plot_all:
        pdir = _plot_dir(options, info, date_str)
        tag = _plot_tag(info, date_str)
        # Calculate sensitivity range (window ± max altitude shift)
        sens_min = fit_result.range_start_m + min(ALT_SHIFTS_M)
        sens_max = fit_result.range_end_m + max(ALT_SHIFTS_M)
        
        # Build time subset info string
        time_subset_text = (
            f"Time subsets:\n"
            f"  • Sample 0: all {data.rcs.shape[0]} profiles\n"
            f"  • Samples 1–{N_TIME_SAMPLES}: random 70% ({int(data.rcs.shape[0] * TIME_SUBSET_FRACTION)} each)"
        )
        
        try:
            plot_rcs_timeseries(
                hours_since_start=data.hours_since_start,
                range_alc=data.range_alc,
                rcs=data.rcs,
                cbh=data.cbh,
                no_cloud_value=info.instrument_type.no_cloud_value,
                title=f"{info.site_name} ({info.wmo_id}) — {date_str} — RCS + Calibration Regions",
                save_path=pdir / f"{tag}_rcs_annotated.png",
                time_datetime=data.time_datetime,
                molecular_window_range=(fit_result.range_start_m, fit_result.range_end_m),
                sensitivity_range=(sens_min, sens_max),
                time_subset_info=time_subset_text,
                time_subset_indices=time_subset_indices,
            )
        except Exception as exc:
            logger.warning(f"plot_rcs_timeseries (annotated) failed: {exc}")

    # =========================================================================
    # Step 7: Validate Rayleigh fit
    # =========================================================================
    # No eligible molecular window passed the validity gates (start above aerosol,
    # R² >= min_window_r2, slope > 0, |b| < a, slope ~ median ratio). The fit comes
    # back with NaN slope / inf relative_error; flag it as non-proportional rather
    # than emitting a spurious constant.
    if not np.isfinite(fit_result.slope) or not np.isfinite(fit_result.relative_error):
        logger.warning("No eligible molecular window (signal not proportional to molecular)")
        return CalibrationResult(
            lidar_constant=-1,
            flag=-2,
            uncertainty=0,
            message="No molecular window passed the validity gates",
        )

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
    # Steps 8-10: Sensitivity analysis over lidar-ratio, altitude window, and
    #             time-profile subsets
    # =========================================================================
    # Run Klett inversion + CL calculation for combinations of
    #   LR ∈ {LRaer-20, LRaer-10, LRaer, LRaer+10, LRaer+20}
    #   altitude shift ∈ {-200m, -100m, 0m, +100m, +200m}
    #   time sample  ∈ {all profiles, + N_TIME_SAMPLES random subsets}
    #
    # The best-estimate CL is the median of all values, and the
    # uncertainty is derived from their spread (half IQR × 2 ≈ robust 2σ).
    # =========================================================================
    n_time_total = 1 + N_TIME_SAMPLES  # all profiles + random subsets
    n_lr = len([d for d in LR_DELTAS if options.lidar_ratio_aerosol + d > 0])
    n_alt = len(ALT_SHIFTS_M)
    n_combos = n_lr * n_alt * n_time_total
    logger.info(
        f"Sensitivity analysis: {n_lr} LR × {n_alt} altitude shifts "
        f"× {n_time_total} time samples = {n_combos} combinations"
    )

    perturbation_results: List[_PerturbationResult] = []

    for t_idx, rcs_mean_t in enumerate(rcs_mean_samples):
        for lr_delta in LR_DELTAS:
            lr = options.lidar_ratio_aerosol + lr_delta
            if lr <= 0:
                continue  # skip non-physical values
            for alt_shift in ALT_SHIFTS_M:
                pr = _compute_cl_for_perturbation(
                    rcs_mean=rcs_mean_t,
                    range_alc=data.range_alc,
                    beta_mol=mol_props.beta_mol,
                    fit_result=fit_result,
                    lidar_ratio_aerosol=lr,
                    altitude_shift_m=alt_shift,
                    subtract_background=options.subtract_background,
                    consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
                    sign_error_v10=options.sign_error_v10,
                )
                if pr is not None:
                    pr.time_sample_idx = t_idx
                    perturbation_results.append(pr)

    n_ok = len(perturbation_results)
    logger.info(f"Successful perturbation runs: {n_ok}/{n_combos}")

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
        sign_error_v10=options.sign_error_v10,
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
        sign_error_v10=options.sign_error_v10,
    )
    if options.sign_error_v10:
        # E-PROF v1.0 (pre-a4e7140): reference inverse-transmittance from the TOTAL
        # optical depth (aerosol + molecular), which overestimates C_L^slope.
        ext_tot_check = ext_aer_check + mol_props.beta_mol * MOLECULAR_LIDAR_RATIO
        optical_depth_ref = np.trapezoid(
            ext_tot_check[:i_start_mol_nominal], data.range_alc[:i_start_mol_nominal]
        )
    else:
        # The Rayleigh-fit slope corresponds to C_L * T_a^2 (aerosol two-way
        # transmittance), because beta_a ~ 0 in the molecular window. To recover
        # C_L we must divide out only the AEROSOL transmittance T_a^2, not the total
        # transmittance T^2 = T_a^2 * T_m^2 (using the total OD overestimates
        # C_L^slope and can cause false 'method disagreement' rejections).
        optical_depth_ref = np.trapezoid(
            ext_aer_check[:i_start_mol_nominal], data.range_alc[:i_start_mol_nominal]
        )
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

    # ── Plot: compact 4x4 Rayleigh diagnostics dashboard ──
    if options.plot_main:
        pdir = _plot_dir(options, info, date_str)
        tag = _plot_tag(info, date_str)
        plot_title_base = f"{info.site_name} ({info.wmo_id}) — {date_str}"

        # Build LR × altitude-shift × time-sample cube for the dashboard
        lr_used = np.array(sorted({pr.lidar_ratio for pr in perturbation_results}))
        shifts_used = np.array(sorted({pr.altitude_shift_m for pr in perturbation_results}))
        ts_used = np.array(sorted({pr.time_sample_idx for pr in perturbation_results}))
        cl_cube = np.full((len(lr_used), len(shifts_used), len(ts_used)), np.nan)
        lr_idx = {v: i for i, v in enumerate(lr_used)}
        sh_idx = {v: i for i, v in enumerate(shifts_used)}
        ts_idx = {v: i for i, v in enumerate(ts_used)}
        for pr in perturbation_results:
            cl_cube[lr_idx[pr.lidar_ratio], sh_idx[pr.altitude_shift_m], ts_idx[pr.time_sample_idx]] = pr.lidar_constant
        cl_matrix = np.nanmedian(cl_cube, axis=2)

        if (
            pr_nominal is not None
            and pr_nominal.signal_normalized is not None
            and fit_result.search_diagnostics is not None
        ):
            diag = fit_result.search_diagnostics
            try:
                plot_rayleigh_diagnostics_compact(
                    range_alc=data.range_alc,
                    altitude=data.altitude,
                    rcs_mean=rcs_mean,
                    signal_normalized=pr_nominal.signal_normalized,
                    p_mol=mol_props.p_mol,
                    beta_att_mol=mol_props.beta_att_mol,
                    fit_altitude_start=fit_result.altitude_start,
                    fit_altitude_end=fit_result.altitude_end,
                    range_bin_m=diag.range_bin_m,
                    half_length_m=diag.half_length_m,
                    slopes=diag.slopes,
                    intercepts=diag.intercepts,
                    r_squared=diag.r_squared,
                    best_range_m=fit_result.center_range_m,
                    best_half_m=fit_result.half_length_m,
                    lr_values=lr_used,
                    alt_shifts=shifts_used,
                    cl_matrix=cl_matrix,
                    cl_median=cl_median,
                    cl_uncertainty=uncertainty,
                    hours_since_start=data.hours_since_start,
                    rcs=data.rcs,
                    used_profile_indices=np.asarray(keep_idx, dtype=int),
                    cloud_base_height=data.cbh,
                    no_cloud_value=info.instrument_type.no_cloud_value,
                    z_low_cloud=options.z_low_cloud,
                    title=f"{plot_title_base} — Rayleigh diagnostics (compact)",
                    save_path=pdir / f"{tag}_rayleigh_diag_compact.png",
                )
            except Exception as exc:
                logger.warning(f"plot_rayleigh_diagnostics_compact failed: {exc}")

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
