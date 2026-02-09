"""
Data loading and preprocessing for ceilometer/lidar observations.

This module handles loading L1 NetCDF files and preprocessing the data
for Rayleigh calibration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, date as date_type
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset, num2date

from .config import InstrumentInfo, InstrumentType, CalibrationOptions


@dataclass
class CeilometerData:
    """Container for ceilometer observation data."""
    # Time information
    time: NDArray[np.float64]           # Days since epoch
    time_datetime: List[datetime]        # Datetime objects
    hours_since_start: NDArray[np.float64]

    # Spatial information
    range_alc: NDArray[np.float64]      # Range bins (m)
    altitude: float                      # Station altitude (m ASL)
    altitude_grid: NDArray[np.float64]  # Altitude of each range bin (m ASL)
    latitude: float                      # Station latitude (degrees north)
    longitude: float                     # Station longitude (degrees east)

    # Observations
    rcs: NDArray[np.float64]            # Range-corrected signal (time x range)
    cbh: NDArray[np.float64]            # Cloud base height (time x layers)

    # Housekeeping
    temperature_optical_module: NDArray[np.float64]
    window_transmission: NDArray[np.float64]
    status_laser: NDArray[np.float64]
    status_detector: NDArray[np.float64]
    laser_life_time: NDArray[np.float64]
    calibration_pulse: NDArray[np.float64]

    # Instrument metadata
    optical_module_id: Optional[int] = None
    instrument_serial_number: Optional[str] = None
    instrument_firmware_version: Optional[str] = None

    # Time metadata
    calendar: str = "gregorian"
    time_units: str = "days since 1970-01-01 00:00:00.000"


def build_file_paths(
    date_str: str,
    info: InstrumentInfo,
    options: CalibrationOptions,
) -> Tuple[Path, Path]:
    """
    Build file paths for current and previous day L1 files.

    Parameters
    ----------
    date_str : str
        Date string in YYYYMMDD format.
    info : InstrumentInfo
        Instrument configuration.
    options : CalibrationOptions
        Calibration options containing folder paths.

    Returns
    -------
    tuple of Path
        (current_day_file, previous_day_file)
    """
    current_date = datetime.strptime(date_str, '%Y%m%d')
    previous_date = current_date - timedelta(days=1)
    previous_str = previous_date.strftime("%Y%m%d")

    folder_current = (
        options.folder_root / info.wmo_id /
        date_str[:4] / date_str[4:6]
    )
    folder_previous = (
        options.folder_root / info.wmo_id /
        previous_str[:4] / previous_str[4:6]
    )

    file_current = folder_current / f"L1_{info.wmo_id}_{info.identifier}{date_str}.nc"
    file_previous = folder_previous / f"L1_{info.wmo_id}_{info.identifier}{previous_str}.nc"

    return file_current, file_previous


def _read_variable_safe(
    dataset: Dataset,
    var_name: str,
    default: Optional[float] = None,
) -> NDArray:
    """Safely read a variable from NetCDF, returning default if not found."""
    if var_name in dataset.variables:
        return dataset.variables[var_name][:]
    if default is not None:
        return np.array([default])
    raise KeyError(f"Variable '{var_name}' not found in dataset")


def _extract_instrument_metadata(
    dataset: Dataset,
    instrument_type: InstrumentType,
) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """
    Extract instrument metadata (serial number, optical module ID, firmware).

    Parameters
    ----------
    dataset : Dataset
        Open NetCDF dataset.
    instrument_type : InstrumentType
        Type of instrument.

    Returns
    -------
    tuple
        (optical_module_id, serial_number, firmware_version)
    """
    om_id = None
    serial = None
    firmware = None

    try:
        if instrument_type in {InstrumentType.CHM15k, InstrumentType.CHM8k}:
            # CHM instruments store these as global attributes with prefix
            serial_attr = getattr(dataset, 'instrument_serial_number', '')
            if serial_attr and len(serial_attr) > 3:
                serial = serial_attr[3:]  # Remove 'CHM' prefix

            om_attr = getattr(dataset, 'optical_module_id', '')
            if om_attr and len(om_attr) > 3:
                om_id = int(om_attr[3:])  # Remove 'TUB' prefix

            firmware = getattr(dataset, 'instrument_firmware_version', None)

        elif instrument_type in {InstrumentType.CL31, InstrumentType.CL51, InstrumentType.CL61}:
            serial = getattr(dataset, 'instrument_serial_number', None)
            # Vaisala instruments don't have optical_module_id

        elif instrument_type in {InstrumentType.MINI_MPL, InstrumentType.MPL}:
            serial = getattr(dataset, 'instrument_serial_number', None)
            fw = getattr(dataset, 'instrument_firmware_version', None)
            if fw and str(fw).lower() != 'unknown':
                firmware = str(fw)

    except (ValueError, AttributeError):
        pass  # Metadata extraction is non-critical

    return om_id, serial, firmware


def load_l1_data(
    file_list: List[Path],
    instrument_type: InstrumentType,
) -> Optional[CeilometerData]:
    """
    Load and concatenate L1 NetCDF files.

    Parameters
    ----------
    file_list : list of Path
        List of L1 file paths to load.
    instrument_type : InstrumentType
        Type of ceilometer/lidar.

    Returns
    -------
    CeilometerData or None
        Loaded data container, or None if loading fails.
    """
    if not file_list:
        return None

    # Initialize accumulators
    time_list = []
    rcs_list = []
    cbh_list = []
    temp_om_list = []
    window_trans_list = []
    status_laser_list = []
    status_detector_list = []
    laser_life_list = []
    cal_pulse_list = []

    range_alc = None
    altitude = None
    latitude = 0.0
    longitude = 0.0
    calendar = None
    time_units = None
    om_id = None
    serial = None
    firmware = None

    is_mini_mpl = instrument_type == InstrumentType.MINI_MPL

    for i, filepath in enumerate(file_list):
        if not filepath.exists():
            continue

        with Dataset(filepath, 'r') as data:
            # Time
            time_list.append(data.variables['time'][:])
            if i == 0:
                calendar = data.variables['time'].calendar
                time_units = data.variables['time'].units

            # Range, altitude, lat/lon (only from first file)
            if range_alc is None:
                range_alc = data.variables['range'][:]
                altitude = float(data.variables['station_altitude'][:])
                latitude = float(data.variables['station_latitude'][:])
                longitude = float(data.variables['station_longitude'][:])
            else:
                # Verify range hasn't changed
                if data.variables['range'].shape != range_alc.shape:
                    raise ValueError("Range shape changed between files")

            # Range-corrected signal
            rcs_list.append(data.variables['rcs_0'][:])

            # Cloud base height
            cbh_list.append(data.variables['cloud_base_height'][:])

            # Housekeeping variables
            temp_om_list.append(data.variables['temperature_optical_module'][:])

            if is_mini_mpl:
                # Mini-MPL doesn't have these variables
                window_trans_list.append(np.full(len(time_list[-1]), np.nan))
                status_laser_list.append(np.full(len(time_list[-1]), np.nan))
                status_detector_list.append(np.full(len(time_list[-1]), np.nan))
                laser_life_list.append(np.full(len(time_list[-1]), np.nan))
                cal_pulse_list.append(np.full(len(time_list[-1]), np.nan))
            else:
                window_trans_list.append(_read_variable_safe(data, 'window_transmission', np.nan))
                status_laser_list.append(_read_variable_safe(data, 'status_laser', np.nan))
                status_detector_list.append(_read_variable_safe(data, 'status_detector', np.nan))
                laser_life_list.append(_read_variable_safe(data, 'laser_life_time', np.nan))
                cal_pulse_list.append(_read_variable_safe(data, 'calibration_pulse', np.nan))

            # Extract metadata from first file
            if i == 0:
                om_id, serial, firmware = _extract_instrument_metadata(data, instrument_type)

    if not time_list:
        return None

    # Concatenate arrays
    time = np.concatenate(time_list)
    rcs = np.concatenate(rcs_list, axis=0)
    cbh = np.concatenate(cbh_list, axis=0)
    temp_om = np.concatenate(temp_om_list)
    window_trans = np.concatenate(window_trans_list)
    status_laser = np.concatenate(status_laser_list)
    status_detector = np.concatenate(status_detector_list)
    laser_life = np.concatenate(laser_life_list)
    cal_pulse = np.concatenate(cal_pulse_list)

    # Convert time to datetime
    try:
        time_datetime = list(num2date(time, time_units, calendar))
    except OverflowError:
        return None  # NaN in time array

    # Calculate hours since start
    hours_since_start = (time - time.min()) * 24

    # Calculate altitude grid
    altitude_grid = range_alc + altitude

    return CeilometerData(
        time=time,
        time_datetime=time_datetime,
        hours_since_start=hours_since_start,
        range_alc=range_alc,
        altitude=altitude,
        altitude_grid=altitude_grid,
        latitude=latitude,
        longitude=longitude,
        rcs=rcs,
        cbh=cbh,
        temperature_optical_module=temp_om,
        window_transmission=window_trans,
        status_laser=status_laser,
        status_detector=status_detector,
        laser_life_time=laser_life,
        calibration_pulse=cal_pulse,
        optical_module_id=om_id,
        instrument_serial_number=serial,
        instrument_firmware_version=firmware,
        calendar=calendar,
        time_units=time_units,
    )


def filter_time_range(
    data: CeilometerData,
    date_str: str,
    options: CalibrationOptions,
) -> CeilometerData:
    """
    Filter data to nighttime window around solar midnight.

    Solar time is computed from UTC using the station longitude:
        solar_hour = UTC_hour + longitude / 15.0

    For date YYYYMMDD (e.g., Feb 5), selects data from:
    - Previous solar day (Feb 4) where solar hour >= hour_min (evening)
    - Current solar day (Feb 5) where solar hour < hour_max (early morning)

    This captures the night from Feb 4 20:00 to Feb 5 04:00 (solar time).

    Parameters
    ----------
    data : CeilometerData
        Input data to filter.
    date_str : str
        Current date string (YYYYMMDD).
    options : CalibrationOptions
        Options containing time window parameters.

    Returns
    -------
    CeilometerData
        Filtered data.
    """
    current_date = datetime.strptime(date_str, '%Y%m%d').date()
    previous_date = current_date - timedelta(days=1)
    next_date = current_date + timedelta(days=1)

    # Solar time offset in hours (positive east of Greenwich)
    solar_offset_hours = data.longitude / 15.0

    # --- Vectorized computation (avoids slow Python loop) -----------------
    # Extract UTC hour/day from each profile
    n = len(data.time_datetime)
    utc_hours = np.empty(n)
    utc_days = np.empty(n, dtype='datetime64[D]')

    for i, dt in enumerate(data.time_datetime):
        utc_hours[i] = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
        utc_days[i] = np.datetime64(date_type(dt.year, dt.month, dt.day))

    # Solar hour and solar date
    solar_hours = utc_hours + solar_offset_hours
    solar_days = utc_days.copy()

    wrap_fwd = solar_hours >= 24
    wrap_bwd = solar_hours < 0
    solar_hours[wrap_fwd] -= 24
    solar_days[wrap_fwd] += np.timedelta64(1, 'D')
    solar_hours[wrap_bwd] += 24
    solar_days[wrap_bwd] -= np.timedelta64(1, 'D')

    prev_np = np.datetime64(previous_date)
    curr_np = np.datetime64(current_date)

    # Evening of previous solar day  OR  early morning of current solar day
    keep_mask = (
        ((solar_days == prev_np) & (solar_hours >= options.hour_min)) |
        ((solar_days == curr_np) & (solar_hours < options.hour_max))
    )

    # Apply mask
    return CeilometerData(
        time=data.time[keep_mask],
        time_datetime=[dt for dt, k in zip(data.time_datetime, keep_mask) if k],
        hours_since_start=data.hours_since_start[keep_mask],
        range_alc=data.range_alc,
        altitude=data.altitude,
        altitude_grid=data.altitude_grid,
        latitude=data.latitude,
        longitude=data.longitude,
        rcs=data.rcs[keep_mask],
        cbh=data.cbh[keep_mask],
        temperature_optical_module=data.temperature_optical_module[keep_mask],
        window_transmission=data.window_transmission[keep_mask],
        status_laser=data.status_laser[keep_mask],
        status_detector=data.status_detector[keep_mask],
        laser_life_time=data.laser_life_time[keep_mask],
        calibration_pulse=data.calibration_pulse[keep_mask],
        optical_module_id=data.optical_module_id,
        instrument_serial_number=data.instrument_serial_number,
        instrument_firmware_version=data.instrument_firmware_version,
        calendar=data.calendar,
        time_units=data.time_units,
    )


def filter_cloudy_profiles(
    data: CeilometerData,
    options: CalibrationOptions,
    no_cloud_value: float,
) -> Tuple[CeilometerData, bool, bool]:
    """
    Filter profiles affected by low clouds.

    Parameters
    ----------
    data : CeilometerData
        Input data.
    options : CalibrationOptions
        Calibration options.
    no_cloud_value : float
        Value indicating no cloud detected.

    Returns
    -------
    tuple
        (filtered_data, is_clear_night, is_partially_clear_night)
    """
    # Check if completely clear
    is_clear_night = np.all(data.cbh == no_cloud_value)

    if is_clear_night:
        return data, True, False

    # Calculate profiles per minute
    if len(data.hours_since_start) > 1:
        dt_hours = data.hours_since_start[1] - data.hours_since_start[0]
        profiles_per_min = int(np.round(1 / (dt_hours * 60)))
    else:
        profiles_per_min = 1

    min_profiles = int(options.min_time_range * 60 * profiles_per_min)

    # Find profiles with low clouds (below threshold and not "no cloud")
    has_low_cloud = np.logical_and(
        data.cbh[:, 0] != no_cloud_value,
        data.cbh[:, 0] <= options.z_low_cloud
    )

    n_clear = np.sum(~has_low_cloud)

    if n_clear < min_profiles:
        return data, False, False

    # Mark profiles contaminated by nearby clouds (15 min window)
    contamination_window = int(profiles_per_min * 15)
    contaminated = np.zeros(len(data.time), dtype=bool)

    for i in range(len(data.time)):
        start_idx = max(0, i - contamination_window + 1)
        end_idx = min(len(data.time), i + contamination_window)
        if np.any(has_low_cloud[start_idx:end_idx]):
            contaminated[i] = True

    n_remaining = np.sum(~contaminated)

    if n_remaining < min_profiles:
        return data, False, False

    # Remove contaminated profiles
    keep_mask = ~contaminated

    filtered_data = CeilometerData(
        time=data.time[keep_mask],
        time_datetime=[dt for dt, k in zip(data.time_datetime, keep_mask) if k],
        hours_since_start=data.hours_since_start[keep_mask],
        range_alc=data.range_alc,
        altitude=data.altitude,
        altitude_grid=data.altitude_grid,
        latitude=data.latitude,
        longitude=data.longitude,
        rcs=data.rcs[keep_mask].copy(),  # Copy to allow modification
        cbh=data.cbh[keep_mask],
        temperature_optical_module=data.temperature_optical_module[keep_mask],
        window_transmission=data.window_transmission[keep_mask],
        status_laser=data.status_laser[keep_mask],
        status_detector=data.status_detector[keep_mask],
        laser_life_time=data.laser_life_time[keep_mask],
        calibration_pulse=data.calibration_pulse[keep_mask],
        optical_module_id=data.optical_module_id,
        instrument_serial_number=data.instrument_serial_number,
        instrument_firmware_version=data.instrument_firmware_version,
        calendar=data.calendar,
        time_units=data.time_units,
    )

    # Mask signal above high clouds (500m below cloud base)
    n_partial = 0
    for i in range(len(filtered_data.time)):
        if filtered_data.cbh[i, 0] != no_cloud_value:
            cloud_height = filtered_data.cbh[i, 0]
            mask_above = filtered_data.range_alc >= (cloud_height - 500)
            filtered_data.rcs[i, mask_above] = np.nan
            n_partial += 1

    # Check if enough profiles remain after masking
    n_final = len(filtered_data.time) - n_partial
    if n_final < min_profiles:
        return data, False, False

    is_partially_clear = n_partial > 0
    return filtered_data, True, is_partially_clear
