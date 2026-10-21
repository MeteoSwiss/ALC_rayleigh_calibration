"""
NetCDF output handling for calibration results.

This module handles writing calibration results to NetCDF files,
including creating new files and appending to existing ones.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
import time

import numpy as np
from netCDF4 import Dataset

from ..config import InstrumentInfo, CalibrationResult


def create_output_directory(output_path: Path) -> None:
    """Create output directory if it doesn't exist."""
    output_path.mkdir(parents=True, exist_ok=True)


def get_output_filepath(
    output_dir: Path,
    info: InstrumentInfo,
    year: int,
) -> Path:
    """
    Generate the output file path for calibration results.

    Parameters
    ----------
    output_dir : Path
        Base output directory.
    info : InstrumentInfo
        Instrument information.
    year : int
        Calibration year.

    Returns
    -------
    Path
        Full path to output file.
    """
    year_dir = output_dir / str(year)
    year_dir.mkdir(parents=True, exist_ok=True)

    filename = f"ALC_calibration_{info.wmo_id}_{info.identifier}{year}.nc"
    return year_dir / filename


def _create_calibration_file(
    filepath: Path,
    info: InstrumentInfo,
) -> Dataset:
    """
    Create a new calibration NetCDF file with all required variables.

    Parameters
    ----------
    filepath : Path
        Output file path.
    info : InstrumentInfo
        Instrument information.

    Returns
    -------
    Dataset
        Open NetCDF dataset for writing.
    """
    ncid = Dataset(filepath, 'w', format='NETCDF4')

    # Global attributes
    ncid.station_name = info.site_name
    ncid.wigos_station_id = info.wmo_id
    ncid.instrument_id = info.identifier
    ncid.identifier = info.identifier
    ncid.instrument_type = info.instrument_type.value
    ncid.title = f"Calibration data for {info.instrument_type.value} at {info.site_name}"
    ncid.history = f"Created by E-PROFILE ({time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())})"
    ncid.source = "Ground based remote sensing"
    ncid.references = "E-PROFILE calibration description"
    ncid.Conventions = "CF-1.8"

    # Create unlimited time dimension
    ncid.createDimension('time', None)

    # Time variables
    nc_time = ncid.createVariable('time', 'f8', ('time',), zlib=True)
    nc_time.long_name = "Central time (UTC) of the calibration period"
    nc_time.units = "days since 1970-01-01 00:00:00.000"
    nc_time.standard_name = "time"
    nc_time.calendar = "gregorian"

    nc_start = ncid.createVariable('start_time', 'f8', ('time',), zlib=True)
    nc_start.long_name = "Start time (UTC) of the calibration period"
    nc_start.units = "days since 1970-01-01 00:00:00.000"

    nc_end = ncid.createVariable('end_time', 'f8', ('time',), zlib=True)
    nc_end.long_name = "End time (UTC) of the calibration period"
    nc_end.units = "days since 1970-01-01 00:00:00.000"

    # Calibration results
    units = info.instrument_type.lidar_constant_units

    nc_cl = ncid.createVariable('lidar_constant', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_cl.long_name = "Lidar constant C_L (Wiegner & Geiss 2012)"
    nc_cl.definition = ("C_L = RCS / beta_att (range-corrected signal divided by attenuated "
                        "total backscatter); calibrate via beta_att = RCS / C_L")
    nc_cl.units = units

    nc_cl_unc = ncid.createVariable('lidar_constant_uncertainty', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_cl_unc.long_name = "Uncertainties on the lidar constant"
    nc_cl_unc.units = units
    nc_cl_unc.comments = "(np.std(CL_all) + error_fit * np.median(CL_all))*2"

    nc_bottom = ncid.createVariable('calibration_bottom_height', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_bottom.long_name = "Bottom height (above sea level) of calibration vertical range"
    nc_bottom.units = "m"

    nc_top = ncid.createVariable('calibration_top_height', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_top.long_name = "Top height (above sea level) of calibration vertical range"
    nc_top.units = "m"

    nc_method = ncid.createVariable('calibration_method', 'i1', ('time',), zlib=True, fill_value=-127)
    nc_method.long_name = "Method of calibration"
    nc_method.flag_values = np.array([0, 1, 2, 3], dtype=np.int8)
    nc_method.flag_meanings = "Rayleigh Liquid_water_clouds Ground_based_lidar Satellite_lidar"

    # Which ALGORITHM VERSION produced this row. A file accumulates rows over years and can be
    # appended to by a newer release, so the version has to be per calibration, not a global
    # attribute -- otherwise a re-processed night is indistinguishable from an original one and a
    # series that straddles an upgrade cannot be interpreted.
    nc_version = ncid.createVariable('calibration_version', 'i2', ('time',), zlib=True,
                                     fill_value=-999)
    nc_version.long_name = "Algorithm version that produced this calibration"
    nc_version.flag_values = np.array(sorted(VERSION_CODES.values()), dtype=np.int16)
    nc_version.flag_meanings = " ".join(k for k, _ in sorted(VERSION_CODES.items(),
                                                             key=lambda kv: kv[1]))
    nc_version.comment = ("code = major*100 + minor*10 + patch of the calibration method, e.g. "
                          "100 = E-PROF v1.0, 120 = v1.2, 200 = v2 (C8), 220 = v2.2 "
                          "(noise-aware molecular-window gates), 230 = v2.3 (v2.2 with a measured "
                          "noise and a vertical-consistency test in its fallback tier)")

    # Housekeeping variables
    nc_lifetime = ncid.createVariable('laser_life_time', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_lifetime.long_name = "Average laser life time during the calibration (nb of shots since factory)"
    nc_lifetime.units = "1"

    nc_wavelength = ncid.createVariable('laser_wavelength', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_wavelength.long_name = "Wavelength of Laser for channel 0"
    nc_wavelength.units = "nm"

    nc_detector = ncid.createVariable('status_detector', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_detector.long_name = "Average detector status during the calibration"
    nc_detector.units = "%"
    nc_detector.comments = "corresponds to state_detector for Lufft"

    nc_laser = ncid.createVariable('status_laser', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_laser.long_name = "Average laser energy during the calibration"
    nc_laser.units = "%"
    nc_laser.comments = "corresponds to state_laser for Lufft"

    nc_temp = ncid.createVariable('temperature_optical_module', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_temp.long_name = "Average laser temperature during the calibration"
    nc_temp.units = "K"
    nc_temp.comments = "corresponds to temp_lom for Lufft and to temperature_laser for Vaisala"

    nc_om_id = ncid.createVariable('optical_module_id', 'i4', ('time',), zlib=True, fill_value=-999)
    nc_om_id.long_name = "ID of the optical module in place during calibration"
    nc_om_id.comments = "corresponds to the number following 'TUB' for Lufft"

    nc_window = ncid.createVariable('window_transmission', 'f8', ('time',), zlib=True, fill_value=-999.9)
    nc_window.long_name = "Average optical windows transmission during the calibration"
    nc_window.units = "%"
    nc_window.comments = "corresponds to Windows transmission estimate for Vaisala and state_optics for Lufft"

    return ncid


# v2.3 uncertainty budget (calibration/rayleigh/uncertainty.py). These variables are created ONLY when
# a row carrying the budget is written (options.uncertainty_model == 'v2.3'), so a v2.2 file keeps its
# exact v2.2 structure. Rows written without the budget (older versions) hold the fill value.
_V23_UNCERTAINTY_VARS = {
    "lidar_constant_uncertainty_ensemble": dict(
        long_name="Ensemble part of the lidar-constant uncertainty (2 sigma)",
        comments=("2 x robust sigma (max of IQR/1.349, MAD/0.6745) of the 125-member ensemble: aerosol "
                  "lidar ratio +-10/+-20 sr x window shift +-100/+-200 m x 5 time samples (all "
                  "profiles + 4 random 70 % subsets). This is the v2.2 lidar_constant_uncertainty.")),
    "lidar_constant_uncertainty_vertical": dict(
        long_name="Vertical-consistency part of the lidar-constant uncertainty (2 sigma)",
        comments=("2 C_L sqrt((|d ln C_L/dz| D)^2 + s_res^2): d ln C_L/dz from C_L recomputed (same "
                  "Klett/C_L machinery) in sub-windows over the usable molecular range, D = "
                  "representative inter-night window displacement, s_res = scatter of the "
                  "sub-window constants about the line not explained by noise. Large when the "
                  "constant depends on the window height (electronic distortion, faint aerosol, "
                  "non-molecular window).")),
    "lidar_constant_uncertainty_molecular": dict(
        long_name="Molecular-model part of the lidar-constant uncertainty (2 sigma)",
        comments=("US Standard 1976 vs the real temperature/pressure: climatology of the CAMS vs "
                  "US-1976 ln-bias of C_L (latitude, season, window height) plus its synoptic "
                  "residual; a floor only when the molecular reference is CAMS.")),
    "lidar_constant_uncertainty_aerosol": dict(
        long_name="Residual-aerosol part of the lidar-constant uncertainty (2 sigma)",
        comments=("Fixed allowance for a vertically uniform aerosol in the window, invisible to any "
                  "single-profile test: 1 sigma = 2 % at 1064 nm (CAMS in-window aerosol scaled by "
                  "its observed transfer to C_L), scaled as (wavelength / 1064 nm)^3.")),
    "lidar_constant_uncertainty_noise": dict(
        long_name="Photon-noise part of the lidar-constant uncertainty (2 sigma)",
        comments=("Noise of the night-mean signal in the window, from the per-gate noise measured by "
                  "temporal first differences on the native profiles (eprof_v2.2).")),
    "lidar_constant_height_gradient": dict(
        long_name="Height gradient of the lidar constant, d ln C_L / dz",
        units="km-1",
        comments=("Noise-weighted slope of ln C_L against the sub-window height (same machinery as "
                  "lidar_constant_uncertainty_vertical). A true lidar constant does not depend on "
                  "height; a persistent non-zero value flags a profile defect.")),
}


def _ensure_v23_uncertainty_variables(ncid: Dataset, units: str) -> None:
    """Create the v2.3 uncertainty variables if missing and document the 2-sigma convention."""
    for name, meta in _V23_UNCERTAINTY_VARS.items():
        if name in ncid.variables:
            continue
        v = ncid.createVariable(name, 'f8', ('time',), zlib=True, fill_value=-999.9)
        v.long_name = meta["long_name"]
        v.units = meta.get("units", units)
        v.comments = meta["comments"]
    unc = ncid.variables['lidar_constant_uncertainty']
    unc.long_name = "Uncertainty of the lidar constant (2 sigma, coverage factor k = 2)"
    unc.coverage_factor = 2
    unc.comments = ("2-sigma (~95 %) absolute uncertainty. Rows with a valid "
                    "lidar_constant_uncertainty_vertical (uncertainty model v2.3): "
                    "2 C_L sqrt(s_ens^2 + s_vert^2 + s_mol^2 + s_aer^2 + s_noise^2), see the five "
                    "lidar_constant_uncertainty_* components. Other rows (v2.2 and older): "
                    "the ensemble part alone (2 x robust sigma of the 125-member ensemble).")


def _write_v23_uncertainty_row(ncid: Dataset, idx: int, components: Optional[dict]) -> None:
    """Write (or clear, when ``components`` is None) the v2.3 components of row ``idx``."""
    for name in _V23_UNCERTAINTY_VARS:
        if name not in ncid.variables:
            continue
        val = None if components is None else components.get(name)
        ok = val is not None and np.isfinite(val)
        ncid.variables[name][idx] = float(val) if ok else -999.9


# Method name -> numeric version code (major*100 + minor*10 + patch). Numeric so the variable can
# carry flag_values/flag_meanings and be compared arithmetically ("rows older than 200").
VERSION_CODES = {
    "eprof_v1.0": 100, "eprof_v1.1": 110, "eprof_v1.2": 120, "eprof_v0.25": 25,
    "eprof_v2": 200, "eprof_v2p": 205, "eprof_v2.2": 220, "eprof_v2.3": 230,
    "earlinet": 300, "bellini": 310, "main": 90, "improved": 95,
    "cloud_oconnor": 1000,          # the liquid-cloud (O'Connor) retrieval
}


def version_code(method_name: str) -> int:
    """Numeric code for a calibration method name; -999 when it is not a known version."""
    return VERSION_CODES.get(str(method_name), -999)


def write_calibration_result(
    output_dir: Path,
    info: InstrumentInfo,
    result: CalibrationResult,
    date_epoch: float,
    time_start: float,
    time_end: float,
    wavelength_nm: float,
    housekeeping: dict,
    method: int = 0,
    version: int = -999,
) -> Path:
    """
    Write calibration result to NetCDF file.

    Creates a new file if one doesn't exist, or appends to existing file. ``method`` tags
    the row in the ``calibration_method`` variable (0=Rayleigh, 1=Liquid_water_clouds, ...);
    a Rayleigh (night) row and a cloud (day) row for the same calendar day coexist because
    the same-day de-duplication below also keys on the method.

    Parameters
    ----------
    output_dir : Path
        Base output directory.
    info : InstrumentInfo
        Instrument information.
    result : CalibrationResult
        Calibration result to write.
    date_epoch : float
        Central date as days since epoch.
    time_start : float
        Start time as days since epoch.
    time_end : float
        End time as days since epoch.
    wavelength_nm : float
        Laser wavelength in nanometers.
    housekeeping : dict
        Dictionary with housekeeping parameters.

    Returns
    -------
    Path
        Path to output file.
    """
    # Output year from the actual calendar day (NOT a 365.25-day approximation, which
    # buckets some Jan-1 dates into the wrong year and would split a same-day Rayleigh
    # night [integer epoch] and cloud day [+0.5 epoch] into different yearly files).
    year = (datetime(1970, 1, 1) + timedelta(days=float(np.floor(date_epoch)))).year
    filepath = get_output_filepath(output_dir, info, year)

    if not filepath.exists():
        # Create new file
        ncid = _create_calibration_file(filepath, info)
        idx = 0
    else:
        # Open existing file
        ncid = Dataset(filepath, 'a')

        # Overwrite a same-day row ONLY if it is the same method; otherwise append.
        # This lets a Rayleigh (method 0, night) row and a cloud (method 1, day) row for
        # the same calendar day coexist, while a re-run of the same method updates in place.
        existing_times = ncid.variables['time'][:]
        existing_methods = ncid.variables['calibration_method'][:]
        same_day = np.floor(existing_times) == np.floor(date_epoch)
        matching = np.where(same_day & (existing_methods == method))[0]

        if len(matching) > 0:
            idx = int(matching[0])
            print(f"  Overwriting existing calibration (method={method}) for this date")
        else:
            idx = len(existing_times)

    # Write data
    ncid.variables['time'][idx] = date_epoch
    ncid.variables['start_time'][idx] = time_start
    ncid.variables['end_time'][idx] = time_end

    ncid.variables['lidar_constant'][idx] = result.lidar_constant
    ncid.variables['lidar_constant_uncertainty'][idx] = result.uncertainty

    # v2.3 uncertainty components: only when the row carries them (uncertainty_model 'v2.3') or the
    # file already has the variables (then a row without them is cleared to fill, never left stale).
    components = getattr(result, "uncertainty_components", None)
    if components:
        _ensure_v23_uncertainty_variables(ncid, info.instrument_type.lidar_constant_units)
    if components or any(n in ncid.variables for n in _V23_UNCERTAINTY_VARS):
        _write_v23_uncertainty_row(ncid, idx, components)

    if result.calibration_bottom_height is not None:
        ncid.variables['calibration_bottom_height'][idx] = result.calibration_bottom_height
    if result.calibration_top_height is not None:
        ncid.variables['calibration_top_height'][idx] = result.calibration_top_height

    ncid.variables['calibration_method'][idx] = method  # 0=Rayleigh, 1=Liquid_water_clouds
    if 'calibration_version' in ncid.variables:
        ncid.variables['calibration_version'][idx] = int(version)

    ncid.variables['laser_wavelength'][idx] = wavelength_nm

    # Write housekeeping data
    for key, var_name in [
        ('laser_life_time', 'laser_life_time'),
        ('status_detector', 'status_detector'),
        ('status_laser', 'status_laser'),
        ('temperature_optical_module', 'temperature_optical_module'),
        ('window_transmission', 'window_transmission'),
    ]:
        if key in housekeeping and not np.isnan(housekeeping[key]):
            ncid.variables[var_name][idx] = housekeeping[key]

    # Handle optical module ID specially (integer)
    if 'optical_module_id' in housekeeping:
        om_id = housekeeping['optical_module_id']
        if om_id is not None and not np.isnan(om_id):
            ncid.variables['optical_module_id'][idx] = int(om_id)

    ncid.close()
    return filepath


def strip_calibration_method(
    output_dir: Path,
    info: InstrumentInfo,
    method: int,
    t_start: Optional[float] = None,
    t_end: Optional[float] = None,
) -> int:
    """Remove rows tagged with ``method`` from this stream's yearly calibration file(s), rewriting each
    file without them. Used before an in-place re-run of a single method so that days which no longer
    yield a calibration do not keep a stale row (``write_calibration_result`` only overwrites on
    success, so a now-rejected day would otherwise retain its old value).

    The removal is restricted to rows whose ``time`` lies in ``[t_start, t_end]`` (days since epoch)
    when given, so re-running one date window never drops calibrations from other windows in the same
    yearly file. Returns the number of rows removed.

    All variables are indexed by the unlimited ``time`` dimension, so the surviving rows are copied
    back generically.
    """
    pattern = f"ALC_calibration_{info.wmo_id}_{info.identifier}*.nc"
    removed = 0
    for filepath in sorted(output_dir.glob(f"*/{pattern}")):
        with Dataset(filepath, 'r') as nc:
            nc.set_auto_mask(False)                      # keep raw fill values, copy them verbatim
            methods = np.asarray(nc.variables['calibration_method'][:])
            times = np.asarray(nc.variables['time'][:])
            remove = methods == method
            if t_start is not None:
                remove &= times >= t_start
            if t_end is not None:
                remove &= times <= t_end
            if not remove.any():
                continue
            keep = ~remove
            n_remove = int(remove.sum())
            kept = {v: np.asarray(nc.variables[v][:])[keep] for v in nc.variables}
        # rewrite the file with only the surviving rows (same structure, from info)
        ncid = _create_calibration_file(filepath, info)
        if any(n in kept for n in _V23_UNCERTAINTY_VARS):      # file had the v2.3 components
            _ensure_v23_uncertainty_variables(ncid, info.instrument_type.lidar_constant_units)
        for v, arr in kept.items():
            if arr.size:
                ncid.variables[v][0:arr.size] = arr
        ncid.close()
        removed += n_remove
    return removed
