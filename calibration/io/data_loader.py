"""
Data loading and preprocessing for ceilometer/lidar observations.

This module handles loading L1 NetCDF files and preprocessing the data
for Rayleigh calibration.
"""

from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timedelta, date as date_type
from functools import lru_cache
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset, num2date

from ..config import InstrumentInfo, InstrumentType, CalibrationOptions, DataLevel


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

    # Vaisala fog indicator. CL31/CL51/CL61 report a vertical_visibility (m) INSTEAD of a
    # cloud base when fog/precip obscures the beam — so a reported (finite, >0) value means
    # FOG and the profile must be excluded from the Rayleigh fit (it has no molecular column).
    # None when the variable is absent (e.g. CHM15k).
    vertical_visibility: Optional[NDArray[np.float64]] = None


def build_file_paths(
    date_str: str,
    info: InstrumentInfo,
    options: CalibrationOptions,
) -> List[Path]:
    """
    Build the list of input files needed to cover the night ending on *date_str*,
    for the configured data level (``options.data_level``).

    Layouts
    -------
    - ``L1``        : ``<root>/<WMO>/YYYY/MM/L1_<WMO>_<id><YYYYMMDD>.nc`` — previous + current day.
    - ``L2_daily``  : ``<root>/<WMO>/YYYY/MM/L2_<WMO>_<id><YYYYMMDD>.nc`` — previous + current day.
    - ``L2_monthly``: ``<root>/<WMO>/YYYY/L2_<WMO>_<id><YYYYMM>.nc`` — current month
      (plus the previous month only when the night starts on the 1st).

    Returns
    -------
    list of Path
        Candidate files in chronological order. (Existence is checked by the caller.)
    """
    current_date = datetime.strptime(date_str, "%Y%m%d")
    wmo, ident = info.wmo_id, info.identifier
    root = options.folder_root / wmo

    if options.data_level == DataLevel.L2_MONTHLY:
        months = [current_date]
        if current_date.day == 1:                      # night spills into previous month
            months.insert(0, current_date.replace(day=1) - timedelta(days=1))
        files = []
        for m in months:
            ym = m.strftime("%Y%m")
            files.append(root / ym[:4] / f"L2_{wmo}_{ident}{ym}.nc")
        return files

    if options.data_level == DataLevel.RAW:
        # Previous + current day. Canonical layout is ONE daily file
        # <root>/<WMO>/YYYYMMDD.nc (concatenated Cloudnet raw, or a CHM15k daily file);
        # the legacy per-day folder <root>/<WMO>/YYYYMMDD/*.nc (many short files) is still
        # accepted as a fallback for days not yet homogenised.
        previous_date = current_date - timedelta(days=1)
        files = []
        for d in (previous_date, current_date):
            stem = d.strftime("%Y%m%d")
            daily = root / f"{stem}.nc"
            folder = root / stem
            if daily.is_file():
                files.append(daily)
            elif folder.is_dir():
                files.extend(sorted(folder.glob("*.nc")))
            else:
                # flat layout: many timestamped files sharing <root> with the date in the
                # filename (e.g. Vaisala CL61 "PAY_CL6101_YYYYMMDD_HHMMSS.nc" / "live_YYYYMMDD_*.nc").
                files.extend(sorted(root.glob(f"*{stem}*.nc")))
        return files

    # Daily layouts (L1 or L2_daily): previous + current day
    prefix = "L1" if options.data_level == DataLevel.L1 else "L2"
    previous_date = current_date - timedelta(days=1)
    files = []
    for d in (previous_date, current_date):
        ymd = d.strftime("%Y%m%d")
        files.append(root / ymd[:4] / ymd[4:6] / f"{prefix}_{wmo}_{ident}{ymd}.nc")
    return files


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
    vert_vis_list = []

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

            # Housekeeping variables. Names differ by manufacturer: Lufft (CHM15k)
            # uses temperature_optical_module, Vaisala (CL61/CL51/CL31) uses
            # temperature_laser; some Vaisala L1 files also omit status_* /
            # laser_life_time. Missing HK is filled with NaN of the CORRECT length
            # (full time axis) so concatenation and the later keep-mask slicing
            # never fail. HK only feeds the output diagnostics, not the fit.
            n_t = len(time_list[-1])

            def _hk(dataset, *names):
                for nm in names:
                    if nm in dataset.variables:
                        return dataset.variables[nm][:]
                return np.full(n_t, np.nan)

            temp_om_list.append(_hk(data, 'temperature_optical_module', 'temperature_laser'))

            if is_mini_mpl:
                # Mini-MPL doesn't have these variables
                window_trans_list.append(np.full(n_t, np.nan))
                status_laser_list.append(np.full(n_t, np.nan))
                status_detector_list.append(np.full(n_t, np.nan))
                laser_life_list.append(np.full(n_t, np.nan))
                cal_pulse_list.append(np.full(n_t, np.nan))
            else:
                window_trans_list.append(_hk(data, 'window_transmission'))
                status_laser_list.append(_hk(data, 'status_laser', 'state_laser'))
                status_detector_list.append(_hk(data, 'status_detector', 'state_detector'))
                laser_life_list.append(_hk(data, 'laser_life_time'))
                cal_pulse_list.append(_hk(data, 'calibration_pulse'))

            # Fog indicator (Vaisala): NaN-filled if absent (CHM15k / Mini-MPL).
            vert_vis_list.append(_hk(data, 'vertical_visibility'))

            # Extract metadata from first file
            if i == 0:
                om_id, serial, firmware = _extract_instrument_metadata(data, instrument_type)

    if not time_list:
        return None

    # Concatenate arrays
    time = np.concatenate(time_list)
    rcs = np.concatenate(rcs_list, axis=0)
    cbh = np.concatenate(cbh_list, axis=0)
    # Normalize the L1 no-cloud sentinel to the instrument convention. E-PROFILE
    # L1 cloud_base_height encodes "no cloud" as -999.9 / -1000 (with a -999.9
    # _FillValue), whereas the cloud filter (and the L2 reader) test against the
    # instrument's no_cloud_value. Real cloud bases are positive, so any
    # non-finite, fill, or negative entry is mapped to the sentinel — otherwise
    # every L1 profile reads as cloudy and no night is ever "clear".
    no_cloud = instrument_type.no_cloud_value
    cbh = np.ma.filled(np.ma.masked_invalid(cbh.astype("f8")), np.nan)
    cbh[~np.isfinite(cbh) | (np.abs(cbh) > 1e30) | (cbh < 0)] = no_cloud
    temp_om = np.concatenate(temp_om_list)
    window_trans = np.concatenate(window_trans_list)
    status_laser = np.concatenate(status_laser_list)
    status_detector = np.concatenate(status_detector_list)
    laser_life = np.concatenate(laser_life_list)
    cal_pulse = np.concatenate(cal_pulse_list)
    vert_vis = np.concatenate(vert_vis_list)

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
        vertical_visibility=vert_vis,
        optical_module_id=om_id,
        instrument_serial_number=serial,
        instrument_firmware_version=firmware,
        calendar=calendar,
        time_units=time_units,
    )


# E-PROFILE L2 stores attenuated_backscatter_0 in 1E-6/(m*sr) units. The `units`
# attribute is unreliable (some files mislabel it "m^-1.sr^-1" while the values are
# still micro-scaled), so — matching the MATLAB reference loadL2Data.m which hardcodes
# `rcs = attBsc * 1e-6 * calConst` — we reconstruct with a FIXED 1e-6 factor.
L2_RCS_FACTOR = 1e-6


@lru_cache(maxsize=4)
def _read_l2_file(filepath_str: str, mtime: float, no_cloud: float):
    """
    Read and reconstruct one L2 file into (time, rcs, cbh, range_alc, station_alt,
    lat, lon, calendar, time_units), with the L1-equivalent rcs already built.

    Cached by (path, mtime): with L2_monthly, every night in a month resolves to the
    same monthly file, so without this cache the (tens-of-MB) file would be re-read
    ~30 times per month. The cache makes it read once. ``mtime`` keys on the file's
    modification time so a re-archived file is not served stale. Returned arrays are
    treated as read-only by the caller (it only concatenates copies), so sharing them
    across nights is safe.
    """
    with Dataset(filepath_str, "r") as d:
        tvar = d.variables["time"]
        time = np.asarray(tvar[:], dtype="f8")
        calendar = getattr(tvar, "calendar", "standard")
        time_units = tvar.units

        alt = np.asarray(d.variables["altitude"][:], dtype="f8")
        station_alt = float(d.variables["station_altitude"][:])
        full_range = alt - station_alt
        keep_range = full_range > 0          # avoid range=0 (divide-by-zero downstream)
        range_alc = full_range[keep_range]
        latitude = float(d.variables["station_latitude"][:])
        longitude = float(d.variables["station_longitude"][:])

        ab_var = d.variables["attenuated_backscatter_0"]
        atten = np.ma.filled(np.ma.masked_invalid(ab_var[:].astype("f8")), np.nan)
        if atten.shape == (len(alt), len(time)):   # (altitude, time) -> (time, range)
            atten = atten.T
        cc = np.asarray(d.variables["calibration_constant_0"][:], dtype="f8")
        # rcs = attBsc * 1e-6 * calConst  (fixed factor, MATLAB loadL2Data.m parity)
        rcs = (atten * cc[:, None] * L2_RCS_FACTOR)[:, keep_range]

        cbh = np.ma.filled(d.variables["cloud_base_height"][:].astype("f8"), np.nan)
        cbh[~np.isfinite(cbh) | (np.abs(cbh) > 1e30)] = no_cloud

        if "vertical_visibility" in d.variables:
            vv = np.ma.filled(np.ma.masked_invalid(d.variables["vertical_visibility"][:].astype("f8")), np.nan)
            vv = np.asarray(vv, dtype="f8").ravel()
            vert_vis = vv if vv.size == len(time) else np.full(len(time), np.nan)
        else:
            vert_vis = np.full(len(time), np.nan)

    return (time, rcs, cbh, vert_vis, range_alc, station_alt, latitude, longitude, calendar, time_units)


def _load_l2_data(
    file_list: List[Path],
    instrument_type: InstrumentType,
) -> Optional[CeilometerData]:
    """
    Load L2 files (daily or monthly) and reconstruct the L1-equivalent rcs.

    rcs[t, z] = attenuated_backscatter_0[z, t] * calibration_constant_0[t] * 1e-6
    (fixed factor, MATLAB loadL2Data.m parity). altitude (ASL) -> range =
    altitude - station_altitude; cloud no-cloud fill is mapped to the instrument's
    no-cloud sentinel; housekeeping is filled with NaN. Per-file reads are cached
    (see ``_read_l2_file``) so a monthly file is read once, not once per night.
    """
    if not file_list:
        return None

    no_cloud = instrument_type.no_cloud_value
    time_list, rcs_list, cbh_list, vv_list = [], [], [], []
    range_alc = altitude = None
    latitude = longitude = 0.0
    calendar = time_units = None

    for filepath in file_list:
        if not filepath.exists():
            continue
        (t, rcs_f, cbh_f, vv_f, rng, station_alt, lat, lon, cal, tu) = _read_l2_file(
            str(filepath), os.path.getmtime(filepath), no_cloud)
        time_list.append(t)
        rcs_list.append(rcs_f)
        cbh_list.append(cbh_f)
        vv_list.append(vv_f)
        if range_alc is None:
            range_alc, altitude = rng, station_alt
            latitude, longitude = lat, lon
            calendar, time_units = cal, tu

    if not time_list:
        return None

    time = np.concatenate(time_list)
    rcs = np.concatenate(rcs_list, axis=0)
    cbh = np.concatenate(cbh_list, axis=0)
    vert_vis = np.concatenate(vv_list)
    nan_hk = np.full(len(time), np.nan)

    try:
        time_datetime = list(num2date(time, time_units, calendar))
    except OverflowError:
        return None

    return CeilometerData(
        time=time,
        time_datetime=time_datetime,
        hours_since_start=(time - time.min()) * 24,
        range_alc=range_alc,
        altitude=altitude,
        altitude_grid=range_alc + altitude,
        latitude=latitude,
        longitude=longitude,
        rcs=rcs,
        cbh=cbh,
        temperature_optical_module=nan_hk.copy(),
        window_transmission=nan_hk.copy(),
        status_laser=nan_hk.copy(),
        status_detector=nan_hk.copy(),
        laser_life_time=nan_hk.copy(),
        calibration_pulse=nan_hk.copy(),
        vertical_visibility=vert_vis,
        calendar=calendar,
        time_units=time_units,
    )


def _sanitize_time_units(units: str) -> str:
    """Drop a trailing timezone offset some raw CHM15k files append (e.g.
    'seconds since 1904-01-01 00:00:00.000 00:00'), which cftime cannot parse."""
    u = str(units).strip()
    if "since" not in u.lower():
        return u
    head, _, rest = u.partition("since")
    parts = rest.strip().split()
    if not parts:
        return u
    date = parts[0]
    tim = parts[1].split(".")[0] if len(parts) > 1 else "00:00:00"
    return f"{head.strip()} since {date} {tim}"


def load_raw_data(
    file_list: List[Path],
    instrument_type: InstrumentType,
) -> Optional[CeilometerData]:
    """Load native RAW ceilometer NetCDF.

    Supports CHM15k (Lufft, variable ``beta_raw``, one daily file) and CL61 (Vaisala /
    Cloudnet, variable ``beta_att``, many ~5-min files per day) — concatenating all files in
    the list into one night. The range-corrected signal ``rcs`` is taken from ``beta_raw`` or
    ``beta_att`` (the molecular fit only needs a signal proportional to molecular in clean
    air; for CL61 the derived constant is then the correction to the factory β_att). Cloud
    base comes from ``cbh`` / ``cloud_base_heights``; housekeeping is NaN (diagnostics only).
    Returns a CeilometerData identical in shape to the L1/L2 readers.
    """
    files = [Path(f) for f in file_list if Path(f).exists()]
    if not files:
        return None
    time_list, rcs_list, cbh_list, vv_list = [], [], [], []
    range_alc = None
    altitude = 0.0
    latitude = 0.0
    longitude = 0.0
    time_units = None
    calendar = "standard"
    no_cloud = instrument_type.no_cloud_value

    for fp in files:
        try:
            d = Dataset(fp, "r")
        except OSError:
            continue
        with d:
            v = d.variables
            if "time" not in v or "range" not in v:
                continue
            sig = None
            for nm in ("beta_raw", "beta_att", "rcs_0"):
                if nm in v:
                    sig = np.asarray(v[nm][:], dtype="f8")
                    break
            if sig is None:
                continue
            this_range = np.asarray(v["range"][:], dtype="f8")
            if range_alc is None:
                range_alc = this_range
                tv = v["time"]
                time_units = _sanitize_time_units(getattr(tv, "units", "seconds since 1970-01-01 00:00:00"))
                calendar = getattr(tv, "calendar", "standard")

                def _scalar(*names, default=0.0):
                    for nm in names:
                        if nm in v:
                            arr = np.asarray(v[nm][:]).ravel()
                            if arr.size and np.isfinite(arr[0]):
                                return float(arr[0])
                    return default
                latitude = _scalar("latitude", "station_latitude")
                longitude = _scalar("longitude", "station_longitude")
                altitude = _scalar("altitude", "station_altitude", default=0.0)
            elif this_range.shape != range_alc.shape:
                continue   # different range grid — skip this sub-file
            time_list.append(np.asarray(v["time"][:], dtype="f8"))
            rcs_list.append(sig)
            cbh_v = None
            for nm in ("cbh", "cloud_base_heights", "cloud_base_height"):
                if nm in v:
                    cbh_v = np.asarray(v[nm][:], dtype="f8")
                    break
            if cbh_v is None:
                cbh_v = np.full((sig.shape[0], 1), no_cloud)
            cbh_list.append(cbh_v if cbh_v.ndim == 2 else cbh_v[:, None])
            if "vertical_visibility" in v:
                vvv = np.asarray(v["vertical_visibility"][:], dtype="f8").ravel()
                vv_list.append(vvv if vvv.size == sig.shape[0] else np.full(sig.shape[0], np.nan))
            else:
                vv_list.append(np.full(sig.shape[0], np.nan))

    if not time_list:
        return None
    time = np.concatenate(time_list)
    rcs = np.concatenate(rcs_list, axis=0)
    maxl = max(c.shape[1] for c in cbh_list)
    cbh = np.concatenate(
        [c if c.shape[1] == maxl else
         np.concatenate([c, np.full((c.shape[0], maxl - c.shape[1]), no_cloud)], axis=1)
         for c in cbh_list], axis=0)
    cbh = np.ma.filled(np.ma.masked_invalid(cbh.astype("f8")), np.nan)
    cbh[~np.isfinite(cbh) | (np.abs(cbh) > 1e30) | (cbh < 0)] = no_cloud

    vert_vis = np.ma.filled(np.ma.masked_invalid(np.concatenate(vv_list).astype("f8")), np.nan)
    order = np.argsort(time)              # short files may arrive unsorted
    time, rcs, cbh, vert_vis = time[order], rcs[order], cbh[order], vert_vis[order]
    try:
        time_datetime = list(num2date(time, time_units, calendar))
    except Exception:
        return None
    hours_since_start = (time - time.min()) / 3600.0   # raw time is in SECONDS
    nan = np.full(len(time), np.nan)
    return CeilometerData(
        time=time, time_datetime=time_datetime, hours_since_start=hours_since_start,
        range_alc=range_alc, altitude=altitude, altitude_grid=range_alc + altitude,
        latitude=latitude, longitude=longitude, rcs=rcs, cbh=cbh,
        temperature_optical_module=nan, window_transmission=nan, status_laser=nan,
        status_detector=nan, laser_life_time=nan, calibration_pulse=nan,
        vertical_visibility=vert_vis,
        calendar=calendar, time_units=time_units,
    )


def load_data(
    file_list: List[Path],
    instrument_type: InstrumentType,
    data_level: DataLevel = DataLevel.L1,
) -> Optional[CeilometerData]:
    """Dispatch to the L1 / RAW / L2 reader by data level."""
    if data_level == DataLevel.L1:
        return load_l1_data(file_list, instrument_type)
    if data_level == DataLevel.RAW:
        return load_raw_data(file_list, instrument_type)
    return _load_l2_data(file_list, instrument_type)


def _block_reduce_mean(arr: NDArray, factor: int, axis: int) -> NDArray:
    """NaN-aware block mean along ``axis`` by an integer ``factor``."""
    if factor <= 1:
        return arr
    n = arr.shape[axis]
    n_full = n // factor
    rem = n - n_full * factor
    arr = np.moveaxis(arr, axis, 0)
    out_blocks = []
    if n_full > 0:
        full = arr[: n_full * factor].reshape((n_full, factor) + arr.shape[1:])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_blocks.append(np.nanmean(full, axis=1))
    if rem > 0:
        tail = arr[n_full * factor:]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_blocks.append(np.nanmean(tail, axis=0, keepdims=True))
    if not out_blocks:
        return arr[:0]
    out = np.concatenate(out_blocks, axis=0)
    return np.moveaxis(out, 0, axis)


def _block_reduce_cloud_base(arr: NDArray, factor: int, axis: int) -> NDArray:
    """Block-reduce a cloud-base-height array by the LOWEST valid cloud base per block.

    The cloud base is the lowest cloud point, so a block of profiles is summarised by the
    *minimum* valid base, NOT the mean. Non-physical entries — the no-cloud sentinel, fill
    values, non-positive or absurdly large heights — are treated as NaN and ignored; a block
    with no valid cloud reduces to NaN (no cloud). Averaging instead (the previous behaviour)
    blends real heights with the no-cloud sentinel, fabricating phantom low clouds and dragging
    high cirrus below the low-cloud screen threshold -> spurious "not a clear night" rejections.
    """
    if factor <= 1:
        return arr
    a = np.moveaxis(np.asarray(arr, dtype="float64"), axis, 0)
    a = np.where(np.isfinite(a) & (a > 0.0) & (a < 20000.0), a, np.nan)
    n = a.shape[0]
    n_full = n // factor
    rem = n - n_full * factor
    out_blocks = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN block -> NaN (no cloud)
        if n_full > 0:
            full = a[: n_full * factor].reshape((n_full, factor) + a.shape[1:])
            out_blocks.append(np.nanmin(full, axis=1))
        if rem > 0:
            out_blocks.append(np.nanmin(a[n_full * factor:], axis=0, keepdims=True))
    if not out_blocks:
        return a[:0]
    return np.moveaxis(np.concatenate(out_blocks, axis=0), 0, axis)


def _datetime64_block_mean(time_dt: List[datetime], factor: int) -> List[datetime]:
    """Block-mean a datetime list and return Python datetimes at block centres."""
    if factor <= 1 or not time_dt:
        return list(time_dt)
    ints = np.asarray(time_dt, dtype="datetime64[ns]").astype("int64").astype("float64")
    reduced = _block_reduce_mean(ints, factor, axis=0)
    ns = reduced.round().astype("int64")
    return [datetime.utcfromtimestamp(int(val) / 1e9) for val in ns]


def average_ceilometer_data(
    data: CeilometerData,
    average_time_s: Optional[float] = None,
    average_range_m: Optional[float] = None,
) -> CeilometerData:
    """Downsample a CeilometerData object in time and range by block averaging.

    The helper is intentionally generic so Rayleigh and cloud calibration can use the
    same reduction logic. Time is averaged by block centres, range by block means.
    Non-signal 1D housekeeping arrays are time-averaged where possible; qualitative
    fields that are not meaningful to average are left untouched if they do not match
    the expected dimensionality.
    """
    if (average_time_s is None or average_time_s <= 0) and (average_range_m is None or average_range_m <= 0):
        return data

    out = CeilometerData(
        time=np.asarray(data.time, dtype="float64").copy(),
        time_datetime=list(data.time_datetime),
        hours_since_start=np.asarray(data.hours_since_start, dtype="float64").copy(),
        range_alc=np.asarray(data.range_alc, dtype="float64").copy(),
        altitude=float(data.altitude),
        altitude_grid=np.asarray(data.altitude_grid, dtype="float64").copy(),
        latitude=float(data.latitude),
        longitude=float(data.longitude),
        rcs=np.asarray(data.rcs, dtype="float64").copy(),
        cbh=np.asarray(data.cbh, dtype="float64").copy(),
        temperature_optical_module=np.asarray(data.temperature_optical_module, dtype="float64").copy(),
        window_transmission=np.asarray(data.window_transmission, dtype="float64").copy(),
        status_laser=np.asarray(data.status_laser, dtype="float64").copy(),
        status_detector=np.asarray(data.status_detector, dtype="float64").copy(),
        laser_life_time=np.asarray(data.laser_life_time, dtype="float64").copy(),
        calibration_pulse=np.asarray(data.calibration_pulse, dtype="float64").copy(),
        optical_module_id=data.optical_module_id,
        instrument_serial_number=data.instrument_serial_number,
        instrument_firmware_version=data.instrument_firmware_version,
        calendar=data.calendar,
        time_units=data.time_units,
        vertical_visibility=(None if data.vertical_visibility is None else np.asarray(data.vertical_visibility, dtype="float64").copy()),
    )

    if len(out.time_datetime) > 1 and average_time_s is not None and average_time_s > 0:
        dt_seconds = np.diff(np.asarray(out.time_datetime, dtype="datetime64[ns]").astype("int64")).astype("float64") / 1e9
        dt_seconds = dt_seconds[np.isfinite(dt_seconds) & (dt_seconds > 0)]
        if dt_seconds.size:
            time_factor = int(np.round(float(average_time_s) / float(np.median(dt_seconds))))
            if time_factor > 1:
                out.time = _block_reduce_mean(out.time, time_factor, axis=0)
                out.time_datetime = _datetime64_block_mean(out.time_datetime, time_factor)
                out.hours_since_start = _block_reduce_mean(out.hours_since_start, time_factor, axis=0)
                if out.cbh.size:
                    # cloud base = lowest point: reduce by min over valid bases, NOT mean
                    # (mean blends the no-cloud sentinel into phantom low clouds).
                    out.cbh = _block_reduce_cloud_base(out.cbh, time_factor, axis=0)
                if out.temperature_optical_module.size:
                    out.temperature_optical_module = _block_reduce_mean(out.temperature_optical_module, time_factor, axis=0)
                if out.window_transmission.size:
                    out.window_transmission = _block_reduce_mean(out.window_transmission, time_factor, axis=0)
                if out.status_laser.size:
                    out.status_laser = _block_reduce_mean(out.status_laser, time_factor, axis=0)
                if out.status_detector.size:
                    out.status_detector = _block_reduce_mean(out.status_detector, time_factor, axis=0)
                if out.laser_life_time.size:
                    out.laser_life_time = _block_reduce_mean(out.laser_life_time, time_factor, axis=0)
                if out.calibration_pulse.size:
                    out.calibration_pulse = _block_reduce_mean(out.calibration_pulse, time_factor, axis=0)
                if out.vertical_visibility is not None and out.vertical_visibility.size:
                    out.vertical_visibility = _block_reduce_mean(out.vertical_visibility, time_factor, axis=0)
                out.rcs = _block_reduce_mean(out.rcs, time_factor, axis=0)

    if len(out.range_alc) > 1 and average_range_m is not None and average_range_m > 0:
        dr = np.diff(out.range_alc)
        dr = dr[np.isfinite(dr) & (dr > 0)]
        if dr.size:
            range_factor = int(np.round(float(average_range_m) / float(np.median(dr))))
            if range_factor > 1:
                out.range_alc = _block_reduce_mean(out.range_alc, range_factor, axis=0)
                out.altitude_grid = out.range_alc + out.altitude
                out.rcs = _block_reduce_mean(out.rcs, range_factor, axis=1)

    return out


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

    use_sza = getattr(options, "use_sza_night", True) and \
        np.isfinite(data.latitude) and np.isfinite(data.longitude)

    if use_sza:
        # Darkness-adaptive window: keep the profiles that are actually dark
        # enough for the molecular fit (solar zenith angle > threshold), within
        # the night centred on the current date's solar midnight (previous-day
        # afternoon/evening + current-day morning, split at 12:00 solar). This
        # uses the full dark period (more in winter, less in summer) and adapts
        # to latitude, instead of a fixed solar-clock window.
        in_night = (
            ((solar_days == prev_np) & (solar_hours >= 12.0)) |
            ((solar_days == curr_np) & (solar_hours < 12.0))
        )
        # Solar zenith angle: hour angle from solar time, declination from doy.
        doy = current_date.timetuple().tm_yday
        decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
        lat = np.radians(data.latitude)
        hour_angle = np.radians((solar_hours - 12.0) * 15.0)
        cos_sza = (np.sin(lat) * np.sin(decl) +
                   np.cos(lat) * np.cos(decl) * np.cos(hour_angle))
        sza = np.degrees(np.arccos(np.clip(cos_sza, -1.0, 1.0)))
        keep_mask = in_night & (sza > getattr(options, "sza_night_threshold", 100.0))
    else:
        # Fallback: fixed solar-clock window (evening of previous solar day OR
        # early morning of current solar day).
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
        vertical_visibility=(data.vertical_visibility[keep_mask]
                             if data.vertical_visibility is not None else None),
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
    instrument_type: Optional["InstrumentType"] = None,
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
    # Fog: Vaisala CL31/CL51/CL61 report a vertical_visibility (m) INSTEAD of a cloud base
    # when fog/precip obscures the beam -> those profiles have no molecular column and must be
    # excluded from the Rayleigh fit (treated like a low cloud). (Edmonton 2026-02-01: ~1/3 of
    # the night was fog with cbh="no cloud", which previously slipped through.)
    # ONLY for the Vaisalas: the CHM15k/Mini-MPL report a vertical visibility ALONGSIDE clouds in
    # non-obscuring conditions, so treating their VV as fog would over-reject usable profiles.
    vv = getattr(data, "vertical_visibility", None)
    use_vv_fog = instrument_type is not None and instrument_type.reports_vv_obscuration
    if use_vv_fog and vv is not None:
        has_fog = np.isfinite(vv) & (vv > 0)
    else:
        has_fog = np.zeros(len(data.cbh), dtype=bool)

    # Check if completely clear (no cloud AND no fog). No-cloud is encoded either as the
    # instrument sentinel (no_cloud_value) or as NaN — the cloud-aware time binning leaves a
    # block with no valid cloud as NaN — so treat both as "no cloud".
    cbh_is_no_cloud = np.isnan(data.cbh) | (data.cbh == no_cloud_value)
    is_clear_night = bool(np.all(cbh_is_no_cloud) and not np.any(has_fog))

    if is_clear_night:
        return data, True, False

    # Profiles per minute. MUST stay a float: rounding to int collapses to 0 for coarse data
    # (e.g. L2 at ~5-min cadence -> 1/(0.083*60)=0.2 -> round 0), which would zero out
    # min_profiles and the contamination window and effectively DISABLE cloud screening — so a
    # cloudy L2 night would pass as "clear". Keeping the float makes the time-based thresholds
    # resolution-independent and identical to L1 at native resolution.
    if len(data.hours_since_start) > 1 and (data.hours_since_start[1] - data.hours_since_start[0]) > 0:
        dt_hours = data.hours_since_start[1] - data.hours_since_start[0]
        profiles_per_min = 1.0 / (dt_hours * 60.0)
    else:
        profiles_per_min = 1.0

    min_profiles = max(1, int(round(options.min_time_range * 60 * profiles_per_min)))

    # Find profiles with low clouds (below threshold and not "no cloud")
    has_low_cloud = np.logical_and(
        data.cbh[:, 0] != no_cloud_value,
        data.cbh[:, 0] <= options.z_low_cloud
    )
    has_low_cloud = has_low_cloud | has_fog   # fog obscures the column -> exclude like a low cloud

    n_clear = np.sum(~has_low_cloud)

    if n_clear < min_profiles:
        return data, False, False

    # Mark profiles contaminated by nearby clouds (15 min window); >=1 profile even on coarse data.
    # Vectorised dilation of has_low_cloud by +/-(window-1): a profile within 15 min of any low
    # cloud is contaminated. (The per-profile Python loop was O(N) -> too slow now the screen runs
    # at native resolution, before the L2-grid binning.)
    contamination_window = max(1, int(round(profiles_per_min * 15)))
    from scipy.ndimage import maximum_filter1d
    contaminated = maximum_filter1d(
        np.asarray(has_low_cloud, dtype=np.uint8),
        size=2 * contamination_window - 1, mode="constant", cval=0) > 0

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
        vertical_visibility=(data.vertical_visibility[keep_mask]
                             if data.vertical_visibility is not None else None),
        optical_module_id=data.optical_module_id,
        instrument_serial_number=data.instrument_serial_number,
        instrument_firmware_version=data.instrument_firmware_version,
        calendar=data.calendar,
        time_units=data.time_units,
    )

    # Mask signal above clouds (500 m below the cloud base) on every remaining profile -- the
    # attenuated signal above a cloud must never enter the fit. A profile is still USABLE for the
    # Rayleigh fit if its molecular window survives the mask: it is cloud-free, OR the cloud sits high
    # enough that the mask starts above the window top (cloud_base - 500 > range_end_m). Counting only
    # the fully-cloud-free profiles (the previous behaviour) wrongly rejected otherwise-fine nights
    # with persistent cirrus WELL ABOVE the 2-6 km window -- e.g. SOFIA 2025-01-28 had 349 fully-clear
    # profiles (< the 3 h minimum) but 2161 with the cloud base >= 6.5 km, i.e. a clean 2-6 km column.
    mask_margin = 500.0
    window_top = options.range_end_m
    n_cloud_masked = 0
    n_usable = 0
    for i in range(len(filtered_data.time)):
        if filtered_data.cbh[i, 0] != no_cloud_value:
            cloud_height = filtered_data.cbh[i, 0]
            filtered_data.rcs[i, filtered_data.range_alc >= (cloud_height - mask_margin)] = np.nan
            n_cloud_masked += 1
            if (cloud_height - mask_margin) > window_top:   # mask starts above the window -> intact
                n_usable += 1
        else:
            n_usable += 1

    # Enough profiles with a usable (clean) molecular window?
    if n_usable < min_profiles:
        return data, False, False

    # A night calibrated with some signal masked above clouds is a PARTIAL success (flag 0.5).
    is_partially_clear = n_cloud_masked > 0
    return filtered_data, True, is_partially_clear
