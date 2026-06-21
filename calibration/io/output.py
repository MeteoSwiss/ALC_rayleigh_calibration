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

    if result.calibration_bottom_height is not None:
        ncid.variables['calibration_bottom_height'][idx] = result.calibration_bottom_height
    if result.calibration_top_height is not None:
        ncid.variables['calibration_top_height'][idx] = result.calibration_top_height

    ncid.variables['calibration_method'][idx] = method  # 0=Rayleigh, 1=Liquid_water_clouds

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
