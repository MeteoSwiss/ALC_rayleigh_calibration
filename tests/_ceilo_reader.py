"""Legacy MATLAB-faithful ceilometer reader (L1 / L2 / Cloudnet-raw) + CeiloData averager.

Moved out of ``calibration.cloud.calibration`` (2026-07): the operational pipeline reads L1/L2
through the single shared loader (``calibration.io.data_loader`` -> ``_ceilo_from_ceilometerdata``)
and no longer needs this reader. It is retained here for the R&D validation scripts and any
Cloudnet-raw (``beta_raw``/``beta_att``) input, which the operational loader does not handle.

MATLAB is no longer the reference, so this is a reference/utility reader, not a parity target.
It reuses the beta-reconstruction primitives that stayed in the operational module.
"""
from __future__ import annotations

import warnings
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset

from calibration.cloud.calibration import (
    CeiloData, CloudCalConfig, set_defaults, _matlab_datenum,
    _beta_conversion_factor, _is_raw_signal, INSTRUMENT_CAL_DEFAULT,
)

def _ncread_matlab(nc: Dataset, name: str) -> NDArray:
    """Read a NetCDF variable in MATLAB ``ncread`` orientation (dimension order reversed).

    MATLAB stores arrays column-major and reverses the NetCDF dimension order on read, so
    a NetCDF variable with dims (altitude, time) is returned by ``ncread`` as
    [time, altitude]. The MATLAB ``read_ceilometer_data`` shape checks (transpose-or-not,
    keep-quality-flag-or-not) are written against that reversed order, so to match MATLAB
    bit-for-bit we present arrays in the same orientation, then apply the identical logic.
    """
    arr = np.asarray(nc.variables[name][:])
    if arr.ndim >= 2:
        arr = np.transpose(arr, axes=tuple(range(arr.ndim - 1, -1, -1)))
    return arr


def convert_time(time_raw: NDArray, time_units: str) -> Tuple[NDArray, NDArray]:
    """Port of ``convert_time``: NetCDF time -> (datetime64[ns], MATLAB datenum).

    Returns both representations; only the datenum is needed downstream (WV nearest-time
    interpolation and the YYYYMM derivation), but datetime64 is kept for clarity.
    """
    u = time_units.lower()
    epoch = None
    scale_s = None
    if "days since 1970-01-01" in u:
        epoch, scale_s = np.datetime64("1970-01-01T00:00:00", "ns"), 86400.0
    elif "seconds since 1970-01-01" in u:
        epoch, scale_s = np.datetime64("1970-01-01T00:00:00", "ns"), 1.0
    elif "hours since" in u:
        ref = time_units.split("hours since", 1)[1].strip()
        ref = ref.replace("T", " ").split(".")[0].strip()
        # try a few formats
        dt = None
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
            try:
                from datetime import datetime as _dt
                dt = _dt.strptime(ref.split()[0] + (" " + ref.split()[1] if len(ref.split()) > 1 else ""), fmt)
                break
            except (ValueError, IndexError):
                continue
        if dt is None:
            epoch = np.datetime64("1970-01-01T00:00:00", "ns")
        else:
            epoch = np.datetime64(dt).astype("datetime64[ns]")
        scale_s = 3600.0
    elif "days since 1904-01-01" in u:
        epoch, scale_s = np.datetime64("1904-01-01T00:00:00", "ns"), 86400.0
    else:
        epoch, scale_s = np.datetime64("1970-01-01T00:00:00", "ns"), 86400.0

    tr = np.asarray(time_raw, dtype="float64")
    ns = np.round(tr * scale_s * 1e9).astype("int64")
    dt64 = epoch + ns.astype("timedelta64[ns]")
    return dt64, _matlab_datenum(dt64)


# ===========================================================================
#  Optional pre-averaging (speed-up; NOT in the MATLAB reference)
# ===========================================================================
def _block_reduce_mean(arr: NDArray, factor: int, axis: int) -> NDArray:
    """NaN-aware block mean along ``axis`` by an integer ``factor``.

    The last partial block (if ``arr`` is not an exact multiple of ``factor``) is
    averaged over its remaining samples. All-NaN blocks return NaN.
    """
    if factor <= 1:
        return arr
    n = arr.shape[axis]
    n_full = n // factor
    rem = n - n_full * factor
    arr = np.moveaxis(arr, axis, 0)
    out_blocks = []
    if n_full > 0:
        full = arr[: n_full * factor]
        full = full.reshape((n_full, factor) + arr.shape[1:])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_blocks.append(np.nanmean(full, axis=1))
    if rem > 0:
        tail = arr[n_full * factor:]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_blocks.append(np.nanmean(tail, axis=0, keepdims=True))
    out = np.concatenate(out_blocks, axis=0)
    return np.moveaxis(out, 0, axis)


def _block_reduce_max(arr: NDArray, factor: int, axis: int) -> NDArray:
    """NaN-aware block max along ``axis`` by an integer ``factor`` (for flags)."""
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
            out_blocks.append(np.nanmax(full, axis=1))
    if rem > 0:
        tail = arr[n_full * factor:]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out_blocks.append(np.nanmax(tail, axis=0, keepdims=True))
    out = np.concatenate(out_blocks, axis=0)
    return np.moveaxis(out, 0, axis)


def _block_reduce_cloud_base(arr: NDArray, factor: int, axis: int) -> NDArray:
    """Block-reduce a cloud-base-height array by the LOWEST valid cloud base per block.

    The cloud base is the lowest cloud point, so a block is summarised by the *minimum* valid
    base, not the mean. Non-physical entries (no-cloud sentinel, fill, <=0 or >=20 km) are
    treated as NaN and ignored; a block with no valid cloud reduces to NaN. Averaging instead
    would blend real heights toward the no-cloud sentinel and mis-place the cloud base used by
    the 500-2400 m cloud-base filter.
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
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if n_full > 0:
            out_blocks.append(np.nanmin(a[: n_full * factor].reshape((n_full, factor) + a.shape[1:]), axis=1))
        if rem > 0:
            out_blocks.append(np.nanmin(a[n_full * factor:], axis=0, keepdims=True))
    return np.moveaxis(np.concatenate(out_blocks, axis=0), 0, axis)


def _matlab_unix_block_time(time_dt: NDArray, factor: int) -> NDArray:
    """Block-mean a datetime64 time vector, returning datetime64[ns] bin centres."""
    ints = time_dt.astype("datetime64[ns]").astype("int64").astype("float64")
    reduced = _block_reduce_mean(ints, factor, axis=0)
    return reduced.round().astype("int64").astype("datetime64[ns]")


def average_ceilo_data(data: CeiloData, config: CloudCalConfig) -> CeiloData:
    """Downsample ``data`` in time and/or range before processing (speed-up).

    Controlled by ``config.average_time_s`` and ``config.average_range_m``. Each is
    converted to an integer block factor from the native resolution (rounded, min 1);
    a factor of 1 leaves that axis untouched. Returns ``data`` unchanged when both
    factors are 1 (so the bit-for-bit MATLAB path is preserved when averaging is off).

    Averaged quantities:
      - beta            (range, time)  -> NaN-aware block mean on both axes
      - time/time_num   (time,)        -> block mean (bin centre)
      - range           (range,)       -> block mean; range_resol scaled
      - cbh             (time,)        -> block mean
      - quality_flag    (range, time)  -> block max (a block is flagged if ANY sample is)
      - window_transmission/laser_energy (time,) -> block mean
    """
    if data.time is None or np.size(data.time) == 0:
        return data

    # --- time factor from native cadence (median dt) ---
    t_factor = 1
    if config.average_time_s and config.average_time_s > 0 and data.time.size > 1:
        dt_s = np.median(np.diff(data.time.astype("datetime64[ns]").astype("int64"))) / 1e9
        if dt_s > 0:
            t_factor = max(1, int(round(config.average_time_s / dt_s)))

    # --- range factor from native gate spacing ---
    r_factor = 1
    if config.average_range_m and config.average_range_m > 0 and data.range.size > 1:
        dr = data.range_resol if data.range_resol > 0 else float(np.median(np.diff(data.range)))
        if dr > 0:
            r_factor = max(1, int(round(config.average_range_m / dr)))

    if t_factor == 1 and r_factor == 1:
        return data

    beta = data.beta  # (range, time)
    if r_factor > 1:
        beta = _block_reduce_mean(beta, r_factor, axis=0)
    if t_factor > 1:
        beta = _block_reduce_mean(beta, t_factor, axis=1)
    data.beta = np.ascontiguousarray(beta)

    # range axis
    if r_factor > 1:
        data.range = _block_reduce_mean(data.range, r_factor, axis=0)
        data.range_resol = float(data.range_resol * r_factor)

    # time axes + per-time vectors
    if t_factor > 1:
        tn = _matlab_unix_block_time(data.time, t_factor)
        data.time = tn
        data.time_num = _block_reduce_mean(np.asarray(data.time_num, dtype="float64"),
                                           t_factor, axis=0)
        if data.cbh is not None and data.cbh.size:
            # cloud base = lowest point: reduce by min over valid bases, not mean.
            data.cbh = _block_reduce_cloud_base(data.cbh, t_factor, axis=0)
        if data.window_transmission is not None and data.window_transmission.size:
            data.window_transmission = _block_reduce_mean(
                data.window_transmission, t_factor, axis=0)
        if data.laser_energy is not None and data.laser_energy.size:
            data.laser_energy = _block_reduce_mean(data.laser_energy, t_factor, axis=0)

    # quality flag (range, time): a reduced cell is "bad" if any contributing sample is
    if data.quality_flag is not None and data.quality_flag.size:
        qf = data.quality_flag.astype("float64")
        if r_factor > 1:
            qf = _block_reduce_max(qf, r_factor, axis=0)
        if t_factor > 1:
            qf = _block_reduce_max(qf, t_factor, axis=1)
        data.quality_flag = qf

    return data


# ===========================================================================
#  read_ceilometer_data
# ===========================================================================
def read_ceilometer_data(nc_file: str, config: CloudCalConfig) -> Tuple[CeiloData, int]:
    """Port of ``read_ceilometer_data`` (E-PROFILE L1/L2 + Cloudnet raw).

    Returns (data, status). status==0 on success. The instrument type read from the
    file's attributes updates ``config.instrument`` in place (as the MATLAB does).
    """
    try:
        with Dataset(nc_file, "r") as nc:
            var_names = set(nc.variables.keys())
            atts = {a: getattr(nc, a) for a in nc.ncattrs()}

            # --- instrument type from attributes (instrument_type, then title) ---
            found = False
            for key in ("instrument_type", "title"):
                if key in atts:
                    val = str(atts[key])
                    vlow = val.lower()
                    if "cl61" in vlow:
                        config.instrument = "CL61"; found = True
                    elif "cl51" in vlow:
                        config.instrument = "CL51"; found = True
                    elif "cl31" in vlow:
                        config.instrument = "CL31"; found = True
                    elif "chm" in vlow:
                        config.instrument = "CHM15k"; found = True
                    if found:
                        break

            # --- time ---
            if "time" not in var_names:
                raise ValueError("Time variable not found")
            time_raw = np.asarray(nc.variables["time"][:], dtype="float64")
            time_units = nc.variables["time"].units
            time_dt, time_num = convert_time(time_raw, time_units)
            n_time = time_dt.size

            # --- station altitude ---
            station_altitude = 0.0
            if "station_altitude" in var_names:
                sa = np.asarray(nc.variables["station_altitude"][:]).ravel()
                if sa.size:
                    station_altitude = float(sa[0])
            if (station_altitude == 0 and np.isfinite(config.station_latitude)
                    and not np.isnan(config.station_latitude)):
                pass  # MATLAB only falls back station_altitude from config.station_altitude
            if station_altitude == 0 and config.station_latitude is not None:
                # config has no station_altitude field in the MATLAB defaults that is
                # used here for altitude; MATLAB uses config.station_altitude (a separate
                # field). We do not have it; keep 0 (matches files that carry it).
                pass

            # --- station lat/lon ---
            # Prefer the canonical scalar station coordinate; fall back to a VALID per-profile value.
            # Reject fill values (in L2 the per-profile 'latitude'/'longitude' are often _FillValue
            # ~1e36) and out-of-range values, so the CAMS lookup never receives garbage coords -- that
            # silently picked the domain-edge cell (bad WV) and now would spuriously trip the
            # 'CAMS too far' (-10) guard.
            def _first_valid_coord(names, lim):
                for nm in names:
                    if nm in var_names:
                        v = np.asarray(nc.variables[nm][:], dtype="float64").ravel()
                        v = v[np.isfinite(v) & (np.abs(v) <= lim)]
                        if v.size:
                            return float(v[0])
                return float("nan")
            station_latitude = _first_valid_coord(("station_latitude", "latitude"), 90.0)
            station_longitude = _first_valid_coord(("station_longitude", "longitude"), 360.0)
            if np.isnan(station_latitude) and not np.isnan(config.station_latitude):
                station_latitude = config.station_latitude
            if np.isnan(station_longitude) and not np.isnan(config.station_longitude):
                station_longitude = config.station_longitude

            # --- range / height ---
            altitude_warning = False
            if "range" in var_names:
                rng = np.asarray(nc.variables["range"][:], dtype="float64")
            elif "height" in var_names:
                rng = np.asarray(nc.variables["height"][:], dtype="float64")
            elif "altitude" in var_names:
                alt = np.asarray(nc.variables["altitude"][:], dtype="float64")
                if station_altitude > 0 and np.min(alt) >= station_altitude:
                    rng = alt - station_altitude
                else:
                    if np.min(alt) < 100:
                        rng = alt
                    else:
                        altitude_warning = True
                        rng = alt - np.min(alt)
            else:
                raise ValueError("Range/height variable not found")

            rng = rng.copy()
            rng[rng < 0] = 0.0
            n_range = rng.size

            if "range_resol" in var_names:
                rr = np.asarray(nc.variables["range_resol"][:]).ravel()
                range_resol = float(rr[0]) if rr.size else float(rng[1] - rng[0])
            else:
                range_resol = float(rng[1] - rng[0])

            # --- beta ---
            # IMPORTANT (orientation): MATLAB ``ncread`` returns arrays with the NetCDF
            # dimension order REVERSED (column-major). To reproduce the MATLAB shape
            # checks verbatim (which decide whether to transpose, and whether the quality
            # flag is kept), we read every multidimensional variable in MATLAB order via
            # ``_ncread_matlab`` (i.e. reversed axes) and then apply the exact MATLAB
            # logic. For this CL61 file, attenuated_backscatter_0 is NetCDF
            # (altitude, time); MATLAB sees (time, altitude) and transposes to
            # (range, time) -- which we reproduce exactly below.
            beta_vars = ["attenuated_backscatter_0", "rcs_0", "beta", "beta_raw",
                         "attenuated_backscatter", "beta_att"]
            beta_var_name = ""
            beta_raw = None
            beta_units = None
            for nm in beta_vars:
                if nm in var_names:
                    beta_var_name = nm
                    raw = _ncread_matlab(nc, nm)
                    raw = np.ma.filled(np.ma.masked_invalid(raw.astype("float64")), np.nan)
                    beta_raw = raw
                    beta_units = getattr(nc.variables[nm], "units", None)
                    break
            if beta_var_name == "":
                raise ValueError("Backscatter variable not found")

            beta_factor = _beta_conversion_factor(config.instrument, beta_units)
            beta = beta_raw.astype("float64") * beta_factor
            # RAW range-corrected signals (L1 rcs_0, counts) are NOT physical backscatter: divide by
            # the calibration constant C so beta = rcs_0 / C is physical 1/(m*sr) (apparent lidar
            # ratio ~18 sr), the SAME scale as the already-calibrated L2 attbsc_0 (which is
            # rcs_0 / C * 1e6, recovered by the *1e-6 factor above). This makes RAW/L1 and L2 give
            # the same cloud coefficient.
            raw_ccal = None
            if _is_raw_signal(beta_units):
                raw_ccal = config.calibration_constant
                if raw_ccal is None:
                    raw_ccal = INSTRUMENT_CAL_DEFAULT.get(config.instrument, 1.0)
                if raw_ccal and np.isfinite(raw_ccal) and raw_ccal != 0:
                    beta = beta / raw_ccal

            # Orient to (range, time) -- exact MATLAB transpose logic on the MATLAB-order array
            if beta.ndim == 2 and beta.shape[0] == n_time and beta.shape[1] == n_range:
                beta = beta.T
            elif beta.ndim == 2 and beta.shape[0] != n_range:
                if beta.shape[1] == n_range:
                    beta = beta.T
                else:
                    raise ValueError(
                        f"Beta dimensions {beta.shape} do not match range {n_range} "
                        f"or time {n_time}")
            beta = np.ascontiguousarray(beta)

            # --- cloud base height (read in MATLAB order) ---
            cbh = None
            for nm in ("cloud_base_height", "cbh"):
                if nm in var_names:
                    cbh_raw = _ncread_matlab(nc, nm)
                    cbh_raw = np.ma.filled(
                        np.ma.masked_invalid(cbh_raw.astype("float64")), np.nan)
                    if cbh_raw.ndim == 2:
                        dims = cbh_raw.shape
                        if dims[0] == n_time:
                            cbh = cbh_raw[:, 0].astype("float64")
                        elif dims[1] == n_time:
                            cbh = cbh_raw[0, :].astype("float64")
                        else:
                            if dims[1] < 10:
                                cbh = cbh_raw[:, 0].astype("float64")
                            elif dims[0] < 10:
                                cbh = cbh_raw[0, :].astype("float64")
                    elif cbh_raw.ndim == 1 and cbh_raw.size == n_time:
                        cbh = cbh_raw.astype("float64")
                    break

            if cbh is not None:
                cbh = cbh.copy()
                cbh[(cbh < 0) | (cbh > 20000)] = np.nan
            else:
                cbh = np.full(n_time, np.nan)

            # --- quality flag (read in MATLAB order; only stored if [range x time]) ---
            # MATLAB: data.quality_flag = qf_raw ONLY if size(qf,1)==n_range AND
            # size(qf,2)==n_time, else []. For this CL61 file MATLAB sees the flag as
            # (time, range) -> condition fails -> [] -> the quality-flag filter is a
            # no-op (0 removed), which we must reproduce.
            quality_flag = None
            if "quality_flag" in var_names:
                qf = np.asarray(_ncread_matlab(nc, "quality_flag"))
                if qf.ndim == 2 and qf.shape[0] == n_range and qf.shape[1] == n_time:
                    quality_flag = qf
                else:
                    quality_flag = None

            # --- window transmission ---
            window_transmission = None
            if "window_transmission" in var_names:
                window_transmission = np.asarray(
                    nc.variables["window_transmission"][:], dtype="float64").ravel()

            # --- laser energy ---
            laser_energy = None
            if "laser_energy" in var_names:
                laser_energy = np.asarray(
                    nc.variables["laser_energy"][:], dtype="float64").ravel()
            elif "laser_pulse_energy" in var_names:
                laser_energy = np.asarray(
                    nc.variables["laser_pulse_energy"][:], dtype="float64").ravel()

            # --- operationally-applied lidar constant (E-PROFILE L2 'calibration_constant_0') ---
            # This is the Wiegner C_L = RCS/beta_att already baked into the file's
            # attenuated_backscatter; it lets the cloud method report an absolute
            # C_L = calibration_constant_applied / C (see module docstring). Robust median over
            # the finite, positive values; None when the variable is absent (e.g. raw L1).
            calibration_constant_applied = None
            for nm in ("calibration_constant_0", "calibration_constant"):
                if nm in var_names:
                    cc = np.asarray(nc.variables[nm][:], dtype="float64").ravel()
                    cc = cc[np.isfinite(cc) & (cc > 0)]
                    if cc.size:
                        calibration_constant_applied = float(np.median(cc))
                    break
            # RAW L1 carries no calibration_constant_0; the C we divided rcs_0 by IS the applied
            # constant, so report it (matches the L2 calibration_constant_0 -> same absolute C_L).
            if calibration_constant_applied is None and raw_ccal:
                calibration_constant_applied = float(raw_ccal)

        data = CeiloData(
            time=time_dt, time_num=time_num,
            station_altitude=station_altitude,
            station_latitude=station_latitude, station_longitude=station_longitude,
            range=rng, range_resol=range_resol, beta=beta, cbh=cbh,
            quality_flag=quality_flag, window_transmission=window_transmission,
            laser_energy=laser_energy, altitude_warning=altitude_warning,
            calibration_constant_applied=calibration_constant_applied,
        )
        return data, 0
    except Exception as exc:  # noqa: BLE001 - mirror MATLAB's try/catch -> status=1
        import warnings
        warnings.warn(f"Error reading NetCDF file: {exc}")
        return CeiloData(  # type: ignore[arg-type]
            time=np.array([]), time_num=np.array([]), station_altitude=0.0,
            station_latitude=float("nan"), station_longitude=float("nan"),
            range=np.array([]), range_resol=0.0, beta=np.zeros((0, 0)),
            cbh=np.array([]), quality_flag=None, window_transmission=None,
            laser_energy=None, altitude_warning=False), 1


def calibrate_file(config: CloudCalConfig):
    """Read a single file (L1 / L2 / Cloudnet-raw), block-average to ``config.average_*``, and run
    the liquid-cloud calibration — the historical ``liquid_cloud_calibration`` flow.

    Retained for Cloudnet-raw and research callers. Operational L1/L2 uses
    ``calibration.cloud.liquid_cloud_calibration`` (the single shared loader), which does not read
    Cloudnet-raw. Raises on read failure (and, via the calibration, on a missing WV correction).
    """
    from calibration.cloud.calibration import liquid_cloud_calibration_from_data
    config = set_defaults(config)
    data, status = read_ceilometer_data(config.nc_file, config)
    if status != 0:
        raise RuntimeError("Failed to read NetCDF file")
    data = average_ceilo_data(data, config)
    return liquid_cloud_calibration_from_data(data, config)
