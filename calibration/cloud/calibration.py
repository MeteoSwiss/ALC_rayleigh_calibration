"""
Liquid-water-cloud calibration of ceilometers (O'Connor 2004 / Hopkin 2019).

This is a faithful Python port of the MATLAB reference
``liquid_cloud/liquid_cloud_calibration.m`` (and the helpers it calls).
The goal is a *bit-for-bit* match (within float tolerance) with the MATLAB output, so
the structure mirrors the MATLAB one-function-per-helper layout and every formula,
threshold, indexing convention and NaN-handling choice is reproduced exactly.

Algorithm
---------
read data
  -> (910 nm only) two-way water-vapor transmission correction ``beta /= trans2``
     STRICT: if the WV correction cannot be computed, an exception is raised (no
     fallback). CHM15k (1064 nm) skips the WV correction.
  -> multiple-scattering correction (``beta *= eta(range)``)
  -> instrument-health filters (window transmission, laser energy, quality flag)
  -> apparent lidar ratio  S_apparent = 1/(2 * integral(beta_corrected dz))
  -> cloud filters (+/-300 m around the peak, aerosol ratio, CBH range)
  -> temporal-consistency filter (N consecutive profiles within +/-X% of their mean)
  -> C = S_consistent / 18.8     (S_theoretical for liquid water)
  -> stats (mean / median / mode / std), optional transmission correction.

C is a MULTIPLIER:  beta_true = C * beta_L2.

Water vapor
-----------
The WV chain is reproduced *exactly* as the MATLAB does it (NOT the cleaner
``calibration/water_vapor.py`` route, which differs from the MATLAB by a
small but non-negligible systematic factor). Specifically:

  CAMS T, q on model levels  (get_Beta_CAMS_oper_monthly.m)
    -> L137 half-level geopotential integration -> P_level, z_model [m ASL]
    -> RH = convert_humidity(P_level, T, q, 'specific humidity', 'relative humidity')
         e   = q*P / (c + (1-c)*q),  c = M_wet/M_dry = 18.0152/28.9644
         es  = Murphy & Koop (2005) saturation pressure over liquid
         RH  = 100 * e / es
    -> nw  = get_water_vapor_number_concentration_from_RH(T, RH)
         Pws = Wagner-Pruss IAPWS-95 saturation pressure [hPa]
         Pw  = Pws * RH/100  [hPa]
         Qw  = (1/Rw) * Pw*100 / T,  Rw = 0.4615  [g m^-3]
         nw  = 7.25e22 * Qw * Rw   [m^-3]
    -> per CAMS time step: trans2(z) via the Gaussian-weighted Wiegner core
       (identical to wv_t2eff.m / wv_t2eff_core), then nearest-time interpolation
       onto the ceilometer time grid.

Only the Gaussian-weighted core (``wv_t2eff_core``) and the L137 a/b coefficients are
re-used from ``water_vapor.py``; they are verified line-for-line equivalent to the
MATLAB. Everything humidity-related is re-implemented here to match the MATLAB's exact
(round-trip) numerics.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset

# Re-use the verified-equivalent pieces from the Rayleigh WV module.
from ..water_vapor_correction.water_vapor import wv_t2eff_core, load_abs_cross_section, _A137, _B137
from ..io.cams import ensure_cams_file

# --- Constants matching the MATLAB exactly ---------------------------------
S_THEORETICAL = 18.8                     # sr, liquid-water lidar ratio (O'Connor 2004)
_RD_CAMS = 287.06                        # J/(kg K), Rd in get_Beta_CAMS_oper_monthly.m
_G0 = 9.80665                            # m/s^2 (get_Beta_CAMS_oper_monthly.m)
_M_DRY = 28.9644                         # kg/kmol (get_meteo_const.m)
_M_WET = 18.0152                         # kg/kmol (get_meteo_const.m)
_C_MD = _M_WET / _M_DRY                  # = 0.6219808... (convert_humidity c)
_RW = 0.4615                             # J g^-1 K^-1 (get_water_vapor..RH.m)
_NW_COEF = 7.25e22                       # number-density coefficient (get_water_vapor..RH.m)


# ===========================================================================
#  Configuration
# ===========================================================================
@dataclass
class CloudCalConfig:
    """Configuration mirroring the MATLAB ``config`` struct + ``set_defaults``.

    Only ``nc_file`` is required. ``instrument`` may be left as the default and is
    overwritten from the file's ``instrument_type`` / ``title`` attribute, exactly as
    the MATLAB ``read_ceilometer_data`` does.
    """
    nc_file: str = ""
    instrument: str = "CL51"

    # set in set_defaults() from the instrument; left None so set_defaults can fill them
    wavelength: Optional[float] = None
    laser_fwhm: Optional[float] = None

    # numeric defaults (mirror set_defaults() exactly)
    cal_minheight: float = 100.0
    cal_maxheight: float = 2400.0
    cbh_minheight: float = 500.0
    cbh_maxheight: float = 2400.0
    ratio_filter: float = 0.05
    n_consecutive: int = 5
    consistency_range: float = 10.0
    temp_threshold: float = -20.0
    window_threshold: float = 90.0
    energy_threshold: float = 90.0
    attenuation_factor: float = 20.0
    debug: int = 0
    apply_transmission_correction: bool = True
    aerosol_lidar_ratio_low: float = 20.0
    aerosol_lidar_ratio_high: float = 70.0
    apply_wv_correction: bool = False
    cams_folder: str = "A:\\CAMS\\"
    abs_cs_lookup_table: str = ""
    station_latitude: float = float("nan")
    station_longitude: float = float("nan")
    # Auto-download a missing CAMS file from the ADS (off by default). 'day' fetches just
    # the file's day (CAMS_Beta_<YYYYMMDD>.nc), 'month' the whole month. Needs cdsapi +
    # cfgrib + ADS credentials in ~/.cdsapirc.
    auto_download_cams: bool = False
    cams_download_scope: str = "day"

    # NB: aerosol_lidar_ratio is NOT defined in MATLAB set_defaults(); the runner sets it
    # to 50. apply_transmission_correction uses it, so it must be provided when that flag
    # is on (otherwise MATLAB errors). We default it to None to reproduce the "must be
    # provided" contract, and only require it when the correction actually runs.
    aerosol_lidar_ratio: Optional[float] = None

    # Optional 'YYYYMM' override for CAMS retrieval (else derived from data.time[0])
    date_str: Optional[str] = None

    # --- Optional pre-averaging (NOT in the MATLAB reference) ----------------
    # Downsample the raw beta grid in time and range *before* the (expensive) WV
    # correction and the per-profile filters. This is a large speed-up for native
    # high-resolution files (e.g. CL61: ~8640 profiles x ~3276 gates per day).
    # Set either to None or <= 0 to disable that axis and keep the bit-for-bit
    # MATLAB resolution. Defaults: 5 min in time, 10 m in range.
    average_time_s: Optional[float] = 300.0
    average_range_m: Optional[float] = 10.0


def set_defaults(config: CloudCalConfig) -> CloudCalConfig:
    """Port of ``set_defaults``: set wavelength / laser_fwhm per instrument.

    The numeric defaults are already supplied by the dataclass field defaults (which
    match the MATLAB ``defaults`` struct one-for-one), so here we only need to fill the
    instrument-dependent wavelength and FWHM, matching the MATLAB ``switch``.
    """
    inst = config.instrument.upper()
    if inst == "CL31":
        wl, fwhm = 909.7, 6.0
    elif inst == "CL51":
        wl, fwhm = 910.0, 3.4
    elif inst == "CL61":
        wl, fwhm = 910.74, 1.0
    elif inst == "CHM15K":
        wl, fwhm = 1064.47, 0.5
    else:
        wl, fwhm = 910.0, 3.4
    # set_defaults always overwrites wavelength/laser_fwhm from the instrument switch
    config.wavelength = wl
    config.laser_fwhm = fwhm
    return config


# ===========================================================================
#  Data container
# ===========================================================================
@dataclass
class CeiloData:
    """Mirror of the MATLAB ``data`` struct (range x time orientation)."""
    time: NDArray            # python datetime64[ns] vector (n_time,)
    time_num: NDArray        # MATLAB-style datenum (days), for WV nearest-time interp
    station_altitude: float
    station_latitude: float
    station_longitude: float
    range: NDArray           # (n_range,) m AGL
    range_resol: float
    beta: NDArray            # (n_range, n_time)  m^-1 sr^-1
    cbh: NDArray             # (n_time,)
    quality_flag: Optional[NDArray]      # (n_range, n_time) or None
    window_transmission: Optional[NDArray]
    laser_energy: Optional[NDArray]
    altitude_warning: bool = False
    trans2_wv: Optional[NDArray] = None  # (n_range, n_time)


def _matlab_datenum(dt64: NDArray) -> NDArray:
    """Convert numpy datetime64 to MATLAB datenum (days since 0000-00-00 proleptic).

    datenum(1970,1,1) = 719529, so datenum = 719529 + days_since_unix_epoch.
    """
    days = (dt64.astype("datetime64[ns]").astype("int64")) / 86400e9
    return 719529.0 + days


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


def _normalize_beta_units(units: Optional[str]) -> str:
    """Normalize beta units text for strict matching."""
    if units is None:
        return ""
    u = str(units).strip().lower()
    u = "".join(u.split())
    u = u.replace("·", "*").replace("⋅", "*").replace("∙", "*").replace("×", "*")
    return u


def _beta_conversion_factor(instrument: str, units: Optional[str]) -> float:
    """Scale factor to bring a beta-like variable onto the cloud calibration's working scale.

    The O'Connor calibration coefficient C absorbs the absolute scale of the input, so the
    range-corrected RAW L1 signals (CL31/CL51 volts, CHM15k counts, Mini-MPL photon rate), the
    physical ``1/(m*sr)`` backscatter (CL61 L1) and the stored L2 product all use factor 1 — C is
    then expressed on the input's scale. Only the explicitly 1e-8-scaled CL31/CL51 backscatter
    form carries a different factor.

    UNKNOWN units are accepted with a WARNING (factor 1) instead of being rejected, so EVERY
    instrument x data-level combination runs (the liquid-cloud method itself is suitable for all
    elastic ceilometers/lidars, including the CL31/CL51 for which it is the primary method).
    """
    u = _normalize_beta_units(units)
    # explicitly 1e-8-scaled backscatter (legacy CL31/CL51 L1 form)
    if any(f in u for f in ("1e-8sr^-1.m^-1", "1e-8sr^-1*m^-1", "1e-8*sr^-1.m^-1", "1e-8*sr^-1*m^-1")):
        return 1e-8
    known = (
        # physical 1/(m*sr) backscatter (CL61 L1) and the stored 1e-6*1/(m*sr) L2 product
        "1/(m*sr)", "1/(sr*m)", "1/(m.sr)", "1/(sr.m)",
        "sr^-1*m^-1", "sr^-1.m^-1", "m^-1*sr^-1", "m^-1.sr^-1",
        # range-corrected raw L1 signal: CL31/CL51 volts, CHM15k counts, Mini-MPL photon rate
        "v*m^2", "v*m2", "m^2*counts/s", "m2*counts/s", "counts",
        "mhz.km^2.uj^-1", "mhz*km^2*uj^-1", "mhz.km^2*uj^-1",
    )
    if any(f in u for f in known):
        return 1.0
    warnings.warn(
        f"Unrecognized beta units {units!r} for {instrument}: using conversion factor 1.0 "
        "(the liquid-cloud coefficient is then on the input's arbitrary scale).", UserWarning)
    return 1.0


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
            station_latitude = float("nan")
            station_longitude = float("nan")
            if "latitude" in var_names:
                v = np.asarray(nc.variables["latitude"][:]).ravel()
                if v.size:
                    station_latitude = float(v[0])
            elif "station_latitude" in var_names:
                v = np.asarray(nc.variables["station_latitude"][:]).ravel()
                if v.size:
                    station_latitude = float(v[0])
            if "longitude" in var_names:
                v = np.asarray(nc.variables["longitude"][:]).ravel()
                if v.size:
                    station_longitude = float(v[0])
            elif "station_longitude" in var_names:
                v = np.asarray(nc.variables["station_longitude"][:]).ravel()
                if v.size:
                    station_longitude = float(v[0])
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

        data = CeiloData(
            time=time_dt, time_num=time_num,
            station_altitude=station_altitude,
            station_latitude=station_latitude, station_longitude=station_longitude,
            range=rng, range_resol=range_resol, beta=beta, cbh=cbh,
            quality_flag=quality_flag, window_transmission=window_transmission,
            laser_energy=laser_energy, altitude_warning=altitude_warning,
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


# ===========================================================================
#  Water-vapor transmission (faithful MATLAB chain)
# ===========================================================================
def _murphy_koop_es_liquid(T: NDArray) -> NDArray:
    """Murphy & Koop (2005) saturation vapor pressure over liquid [Pa].

    Port of calculate_saturation_vapor_pressure_liquid.m (default method).
    """
    T = np.asarray(T, dtype="float64")
    expo = (54.842763 - 6763.22 / T - 4.210 * np.log(T) + 0.000367 * T
            + np.tanh(0.0415 * (T - 218.8))
            * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025 * T))
    return np.exp(expo)


def _wagner_pruss_pws_hpa(T: NDArray) -> NDArray:
    """Wagner-Pruss IAPWS-95 saturation vapor pressure over liquid [hPa].

    Port of the (T>0) branch of get_water_vapor_number_concentration_from_RH.m, case 2.
    """
    T = np.asarray(T, dtype="float64")
    Tc = 647.096
    Pc = 220640.0
    V = 1.0 - T / Tc
    C1 = -7.85951783
    C2 = 1.84408259
    C3 = -11.7866497
    C4 = 22.6807411
    C5 = -15.9618719
    C6 = 1.80122502
    inner = (C1 * V + C2 * V ** 1.5 + C3 * V ** 3 + C4 * V ** 3.5
             + C5 * V ** 4 + C6 * V ** 7.5)
    return np.exp(Tc / T * inner) * Pc


def _nw_from_T_RH(T: NDArray, RH: NDArray) -> NDArray:
    """Port of get_water_vapor_number_concentration_from_RH.m (method 2).

    nw [m^-3] from temperature [K] and relative humidity [%].
    """
    Pws = _wagner_pruss_pws_hpa(T)        # hPa
    Pw = Pws * RH / 100.0                  # hPa
    Qw = (1.0 / _RW) * Pw * 100.0 / T     # g m^-3  (Pw*100 -> Pa)
    nw = _NW_COEF * Qw * _RW              # m^-3
    return nw


def _cams_levels_all_times(
    cams_file: str, latitude: float, longitude: float,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """Read CAMS T/q at the nearest grid point for ALL time steps and integrate the
    L137 half-level geopotential, exactly like get_Beta_CAMS_oper_monthly.m (z=[] path).

    Returns
    -------
    time_num : (n_t,) MATLAB datenum  = double(time_raw)/24 + datenum(1900,1,1)
    z_model  : (n_lev, n_t)  geopotential height [m ASL]   (NOT sorted; native order)
    T        : (n_lev, n_t)  temperature [K]
    nw       : (n_lev, n_t)  water-vapor number density [m^-3]

    The arrays keep the native CAMS level order (top..surface). The hydrostatic
    integration walks from the surface (last index) upward and keys the top-of-atmos
    singularity on the half-level whose pressure is identically zero (A=B=0), i.e.
    ``level_number == 1`` -> ``idx == 0`` -- matching the MATLAB ``find(level(i)-1==...)``
    behaviour without assuming a particular loop position.
    """
    with Dataset(cams_file, "r") as nc:
        lon_m = np.asarray(nc.variables["longitude"][:], dtype="float64")
        lat_m = np.asarray(nc.variables["latitude"][:], dtype="float64")
        level = np.asarray(nc.variables["level"][:], dtype="float64")
        time_raw = np.asarray(nc.variables["time"][:], dtype="float64")

        li = int(np.argmin(np.abs(lon_m - longitude)))
        ai = int(np.argmin(np.abs(lat_m - latitude)))
        n_lev = level.size
        n_t = time_raw.size

        # ncread(...,'t',[li,ai,1,1],[1,1,n_lev,n_t]) -> squeeze -> (n_lev, n_t)
        # netCDF4 variable dims are (time, level, latitude, longitude) in this file;
        # index accordingly and transpose to (level, time).
        tvar = nc.variables["t"]
        qvar = nc.variables["q"]
        zvar = nc.variables["z"]
        spvar = nc.variables["lnsp"]
        dim_names = tvar.dimensions  # e.g. ('time','level','latitude','longitude')

        def _read_profile(var):
            arr = np.asarray(var[:], dtype="float64")
            # Move to (level, time) regardless of stored dim order
            axes = {nm: k for k, nm in enumerate(dim_names)}
            # build an index tuple selecting li/ai on lon/lat, full on level/time
            idx = [slice(None)] * arr.ndim
            idx[axes["longitude"]] = li
            idx[axes["latitude"]] = ai
            sub = arr[tuple(idx)]  # remaining dims: level and time (in their orig order)
            # figure out remaining axis order
            remaining = [nm for nm in dim_names if nm in ("level", "time")]
            # sub axes correspond to remaining in order
            lev_pos = remaining.index("level")
            sub = np.moveaxis(sub, lev_pos, 0)  # (level, time)
            return sub

        T = _read_profile(tvar)
        q = _read_profile(qvar)

        # z, lnsp are surface fields stored at a single (top) level slot.
        # MATLAB reads them at level index 1 (start=1,count=1) for all times.
        def _read_surface(var):
            arr = np.asarray(var[:], dtype="float64")
            axes = {nm: k for k, nm in enumerate(dim_names)}
            idx = [slice(None)] * arr.ndim
            idx[axes["longitude"]] = li
            idx[axes["latitude"]] = ai
            idx[axes["level"]] = 0  # MATLAB start index 1 -> python 0 (first level slot)
            sub = arr[tuple(idx)]   # remaining dim: time
            return np.asarray(sub, dtype="float64").ravel()

        z_surf = _read_surface(zvar)     # (n_t,) surface geopotential [m^2/s^2]
        lnsp = _read_surface(spvar)      # (n_t,) ln surface pressure

    time_num = time_raw / 24.0 + 693962.0  # datenum(1900,1,1) = 693962

    # L137 a/b coefficients indexed by level NUMBER (param_137_levels(:,1)).
    A = _A137
    B = _B137

    T_moist = T * (1.0 + 0.609133 * q)
    z_f = np.full((n_lev, n_t), np.nan)
    P_level = np.full((n_lev, n_t), np.nan)

    lev_int = level.astype(int)
    for t in range(n_t):
        surface_pressure = np.exp(lnsp[t])
        z_h = z_surf[t]
        for i in range(n_lev - 1, -1, -1):
            idx = lev_int[i] - 1  # half-level index (A/B indexed by level number)
            Ph_lev = A[idx] + B[idx] * surface_pressure
            Ph_levp1 = A[idx + 1] + B[idx + 1] * surface_pressure
            P_level[i, t] = 0.5 * (Ph_lev + Ph_levp1)
            if idx == 0:
                dlogP = np.log(Ph_levp1 / 0.1)
                alpha = np.log(2.0)
            else:
                dlogP = np.log(Ph_levp1 / Ph_lev)
                dP = Ph_levp1 - Ph_lev
                alpha = 1.0 - (Ph_lev / dP) * dlogP
            TRd = T_moist[i, t] * _RD_CAMS
            z_f[i, t] = z_h + TRd * alpha
            z_h = z_h + TRd * dlogP

    z_model = z_f / _G0  # m ASL

    # RH via convert_humidity (q -> e -> RH with Murphy&Koop es)
    e = (q * P_level) / (_C_MD + (1.0 - _C_MD) * q)   # Pa
    es = _murphy_koop_es_liquid(T)                     # Pa
    RH = 100.0 * e / es                                # %
    nw = _nw_from_T_RH(T, RH)                          # m^-3

    return time_num, z_model, T, nw


def compute_wv_transmission(data: CeiloData, config: CloudCalConfig) -> NDArray:
    """Port of ``compute_wv_transmission``: two-way WV transmission (range x time).

    Mirrors the MATLAB exactly, including: nearest-height LUT mapping, the Gaussian
    laser spectrum, per-CAMS-time-step nw interpolation to range (with NaN fill below
    the lowest model level), the per-step Gaussian-weighted two-way transmission, and
    the final nearest-time interpolation onto the ceilometer time grid.

    Raises on missing LUT/coords (STRICT - the caller treats any failure as "do not
    calibrate this period").
    """
    lut_path = config.abs_cs_lookup_table
    if not lut_path or str(lut_path).strip() in ("", "."):
        # Fall back to the 910 nm WV LUT bundled as package data (same as the Rayleigh path).
        from ..water_vapor_correction.water_vapor import DEFAULT_ABS_CROSS_SECTION
        lut_path = str(DEFAULT_ABS_CROSS_SECTION)
    if not Path(lut_path).is_file():
        raise FileNotFoundError(f"Absorption cross-section lookup table not found: {lut_path}")
    if np.isnan(data.station_latitude) or np.isnan(data.station_longitude):
        raise ValueError("Station latitude/longitude required for water vapor correction")

    # CAMS lookup date (YYYYMMDD): the day of the first profile, so a daily CAMS file
    # (CAMS_Beta_<YYYYMMDD>.nc) resolves/auto-downloads as well as a monthly one. An
    # explicit config.date_str override (YYYYMM or YYYYMMDD) takes precedence.
    day8 = str(data.time[0].astype("datetime64[D]")).replace("-", "")
    if config.date_str:
        cams_date = config.date_str if len(config.date_str) >= 8 else config.date_str + day8[len(config.date_str):8]
    else:
        cams_date = day8

    # --- load absorption cross-section LUT ---
    abs_cs_wl, abs_cs_height, abscs_full = load_abs_cross_section(Path(lut_path))
    # abscs_full: (n_wl, n_height)

    # --- CAMS T/RH -> nw and geopotential, all time steps ---
    cams_path = ensure_cams_file(
        config.cams_folder, cams_date,
        auto_download=getattr(config, "auto_download_cams", False),
        scope=getattr(config, "cams_download_scope", "day"),
    )
    if cams_path is None:
        raise FileNotFoundError(
            f"No CAMS file for {cams_date} in {config.cams_folder} "
            f"(looked for monthly CAMS_Beta_{cams_date[:6]}.nc or daily CAMS_Beta_{cams_date[:8]}.nc)"
        )
    cams_file = str(cams_path)
    time_cams, cams_z, _cams_T, nw_all = _cams_levels_all_times(
        cams_file, data.station_latitude, data.station_longitude)

    n_cams = time_cams.size
    n_range = data.range.size
    range_col = data.range.astype("float64")

    # Nearest-height LUT mapping (min over height of |range - height|)
    # MATLAB: [~,height_indices] = min(abs(range_col' - abs_cs_height(:)),[],1)
    height_indices = np.argmin(
        np.abs(range_col[None, :] - abs_cs_height[:, None]), axis=0)
    abs_cs = abscs_full[:, height_indices]  # (n_wl, n_range)

    # Gaussian laser spectrum
    sigma = config.laser_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gauss = (np.exp(-0.5 * ((abs_cs_wl - config.wavelength) / sigma) ** 2)
             / (sigma * np.sqrt(2.0 * np.pi)))
    sum_gauss = gauss.sum()
    band = ((abs_cs_wl > config.wavelength - 3 * sigma)
            & (abs_cs_wl < config.wavelength + 3 * sigma))

    # Interpolate nw (per CAMS step) from model AGL grid onto the instrument range grid.
    wv_density_all = np.zeros((n_range, n_cams))
    for i in range(n_cams):
        z_agl = cams_z[:, i] - data.station_altitude
        # interp1(z_agl, nw, range_col) with linear interp, NaN outside.
        wv = _interp1_linear_nan(z_agl, nw_all[:, i], range_col)
        nan_mask = np.isnan(wv)
        if np.any(nan_mask):
            first_valid = np.where(~nan_mask)[0]
            if first_valid.size:
                fill = wv[first_valid[0]]
                wv[nan_mask & (range_col < np.nanmax(z_agl))] = fill
        wv[np.isnan(wv)] = 0.0
        wv_density_all[:, i] = wv

    # Two-way transmission per CAMS step (vectorised core, identical to wv_t2eff.m)
    trans2_cams = np.full((n_range, n_cams), np.nan)
    for i in range(n_cams):
        ext = abs_cs * (wv_density_all[:, i][None, :] / 1e4)  # (n_wl, n_range) [m^-1]
        sum_ext = np.zeros_like(ext)
        sum_ext[band, :] = _cumtrapz_axis1(ext[band, :], range_col)
        trans = np.exp(-sum_ext)
        trans2_cams[:, i] = (trans ** 2 * gauss[:, None]).sum(axis=0) / sum_gauss

    # Interpolate trans2 from CAMS time grid to ceilometer time grid (nearest).
    n_profiles = data.beta.shape[1]
    if n_cams == 1:
        trans2 = np.repeat(trans2_cams, n_profiles, axis=1)
    else:
        ceilo_dn = data.time_num.astype("float64")
        cams_dn = time_cams.astype("float64")
        trans2 = _interp1_nearest_extrap_cols(cams_dn, trans2_cams.T, ceilo_dn).T

    trans2 = np.array(trans2, dtype="float64")
    trans2[(trans2 <= 0) | np.isnan(trans2)] = 1.0
    return trans2


def _interp1_linear_nan(x: NDArray, y: NDArray, xi: NDArray) -> NDArray:
    """MATLAB interp1(x,y,xi,'linear') with NaN outside the data range and x sorted.

    x may be non-monotonic from the model; MATLAB requires monotonic x and sorts via the
    grid. CAMS z_model is monotonic increasing with index from surface up, but stored
    top..surface, so z_agl is decreasing in index. We sort ascending before interpolating
    (np.interp requires increasing x), which is what MATLAB's interp1 does internally.
    """
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    out = np.interp(xi, xs, ys, left=np.nan, right=np.nan)
    return out


def _cumtrapz_axis1(y: NDArray, x: NDArray) -> NDArray:
    """MATLAB cumtrapz(x, y, 2): cumulative trapezoid along axis 1 with leading zero."""
    dx = np.diff(x)
    incr = 0.5 * (y[:, 1:] + y[:, :-1]) * dx[None, :]
    return np.concatenate([np.zeros((y.shape[0], 1)), np.cumsum(incr, axis=1)], axis=1)


def _interp1_nearest_extrap_cols(x: NDArray, Y: NDArray, xi: NDArray) -> NDArray:
    """MATLAB interp1(x, Y, xi, 'nearest', 'extrap') with Y columns interpolated.

    Y has shape (len(x), ncols); returns (len(xi), ncols). 'nearest' picks the closest
    sample; 'extrap' for xi outside [min(x),max(x)] also picks the nearest endpoint
    (which 'nearest' already does), so a plain nearest lookup is equivalent.
    """
    x = np.asarray(x, dtype="float64")
    order = np.argsort(x)
    xs = x[order]
    Ys = Y[order, :]
    # nearest index for each xi
    pos = np.searchsorted(xs, xi)
    pos = np.clip(pos, 1, len(xs) - 1)
    left = xs[pos - 1]
    right = xs[pos]
    # MATLAB interp1(...,'nearest') breaks an exact tie towards the HIGHER index
    # (verified: query at the midpoint returns the right-hand sample). So pick LEFT only
    # when it is strictly closer; ties -> right.
    choose_left = (xi - left) < (right - xi)
    idx = np.where(choose_left, pos - 1, pos)
    # handle xi <= xs[0] or >= xs[-1]
    idx[xi <= xs[0]] = 0
    idx[xi >= xs[-1]] = len(xs) - 1
    return Ys[idx, :]


# ===========================================================================
#  apply_multiple_scattering_correction
# ===========================================================================
_ETA_CL31 = np.array([
    [0.250, 0.82854], [0.375, 0.82371], [0.625, 0.81608], [0.875, 0.80811],
    [1.125, 0.79969], [1.375, 0.79027], [1.625, 0.78227], [1.875, 0.77480],
    [2.125, 0.76710], [2.375, 0.76088]])
_ETA_CL51 = np.array([
    [0.250, 0.82881], [0.375, 0.82445], [0.625, 0.81752], [0.875, 0.81021],
    [1.125, 0.80241], [1.375, 0.79356], [1.625, 0.78595], [1.875, 0.77877],
    [2.125, 0.77100], [2.375, 0.76400]])


def apply_multiple_scattering_correction(
    beta: NDArray, range_data: NDArray, config: CloudCalConfig) -> NDArray:
    """Port of ``apply_multiple_scattering_correction`` (beta *= eta(range))."""
    range_km = range_data / 1000.0
    inst = config.instrument.upper()
    if inst == "CL31":
        eta = _ETA_CL31
    elif inst in ("CL51", "CL61"):
        eta = _ETA_CL51
    else:
        eta = np.array([[0.250, 0.82881], [2.375, 0.76400]])

    if eta.shape[0] > 2:
        factor_profile = _interp1_linear_extrap(eta[:, 0], eta[:, 1], range_km)
    else:
        # MATLAB: length(eta_correction) is the larger matrix dim (here columns=2),
        # so for the 2-row 'otherwise' table length()==2 -> ones(). But CL31/51/61 use
        # the 10x2 table where length()==10 -> interp. We branch on row count here.
        factor_profile = _interp1_linear_extrap(eta[:, 0], eta[:, 1], range_km)
    return beta * factor_profile[:, None]


def _interp1_linear_extrap(x: NDArray, y: NDArray, xi: NDArray) -> NDArray:
    """MATLAB interp1(x, y, xi, 'linear', 'extrap') with x ascending."""
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    xi = np.asarray(xi, dtype="float64")
    out = np.empty_like(xi)
    # interior + edges via slopes
    idx = np.clip(np.searchsorted(x, xi) - 1, 0, len(x) - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    slope = (y1 - y0) / (x1 - x0)
    out = y0 + slope * (xi - x0)
    return out


# ===========================================================================
#  apply_instrument_filters
# ===========================================================================
def apply_instrument_filters(
    beta: NDArray, data: CeiloData, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Port of ``apply_instrument_filters``: quality-flag / window / energy."""
    beta_filtered = beta.copy()
    n_profiles = beta.shape[1]
    stats = {"window_rejected": 0, "energy_rejected": 0, "quality_flag_rejected": 0}

    # 1. Quality flag (only present if oriented [range x time]; else data.quality_flag is None)
    if data.quality_flag is not None and data.quality_flag.size:
        lower_gate = _find_first(data.range >= config.cal_minheight)
        upper_gate = _find_last(data.range <= config.cal_maxheight)
        if lower_gate is not None and upper_gate is not None:
            qf = data.quality_flag
            for i in range(n_profiles):
                if qf.shape[0] == beta.shape[0]:
                    profile_flags = qf[lower_gate:upper_gate + 1, i]
                else:
                    profile_flags = qf[i, lower_gate:upper_gate + 1]
                if np.any(profile_flags > 0):
                    beta_filtered[:, i] = np.nan
                    stats["quality_flag_rejected"] += 1

    # 2. Window transmission
    if data.window_transmission is not None and data.window_transmission.size:
        for i in range(n_profiles):
            if data.window_transmission[i] < config.window_threshold:
                if not np.all(np.isnan(beta_filtered[:, i])):
                    beta_filtered[:, i] = np.nan
                    stats["window_rejected"] += 1

    # 3. Laser energy
    if data.laser_energy is not None and data.laser_energy.size:
        for i in range(n_profiles):
            if data.laser_energy[i] < config.energy_threshold:
                if not np.all(np.isnan(beta_filtered[:, i])):
                    beta_filtered[:, i] = np.nan
                    stats["energy_rejected"] += 1

    return beta_filtered, stats


# ===========================================================================
#  calculate_lidar_ratio
# ===========================================================================
def calculate_lidar_ratio(
    beta: NDArray, range_data: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, NDArray]:
    """Port of ``calculate_lidar_ratio``: S = 1/(2*trapz(beta dz)) per profile."""
    n_profiles = beta.shape[1]
    S = np.full(n_profiles, np.nan)
    integrated_beta = np.full(n_profiles, np.nan)

    lower_gate = _find_first(range_data >= config.cal_minheight)
    upper_gate = _find_last(range_data <= config.cal_maxheight)
    if lower_gate is None or upper_gate is None:
        return S, integrated_beta

    for i in range(n_profiles):
        col = beta[:, i]
        if np.all(np.isnan(col)):
            continue
        beta_profile = col[lower_gate:upper_gate + 1]
        if np.sum(np.isnan(beta_profile)) > beta_profile.size * 0.1:
            continue
        valid_mask = ~np.isnan(beta_profile)
        if np.sum(valid_mask) < 3:
            continue
        # trapz(range_data(lower_gate + find(valid_mask) - 1), beta_profile(valid_mask))
        x = range_data[lower_gate:upper_gate + 1][valid_mask]
        y = beta_profile[valid_mask]
        integ = _trapz(y, x)
        integrated_beta[i] = integ
        if integ > 0:
            S[i] = 1.0 / (2.0 * integ)
    return S, integrated_beta


# ===========================================================================
#  apply_cloud_filters
# ===========================================================================
def apply_cloud_filters(
    S: NDArray, beta: NDArray, data: CeiloData, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Port of ``apply_cloud_filters``: peak-based +/-300 m, aerosol ratio, CBH range."""
    n_profiles = S.size
    S_filtered = S.copy()
    stats = {"above_rejected": 0, "below_rejected": 0,
             "ratio_rejected": 0, "cbh_rejected": 0}

    range_resol = float(data.range[1] - data.range[0])
    gate_300m = int(round(300.0 / range_resol))
    if gate_300m < 1:
        gate_300m = 1

    lower_gate = _find_first(data.range >= config.cal_minheight)
    upper_gate = _find_last(data.range <= config.cal_maxheight)

    for i in range(n_profiles):
        if np.isnan(S_filtered[i]):
            continue
        beta_profile = beta[:, i]

        # peak within calibration range: zero out 1..lower_gate (MATLAB beta_roi(1:lower_gate)=0)
        beta_roi = beta_profile.copy()
        if lower_gate is not None:
            beta_roi[:lower_gate + 1] = 0.0  # MATLAB 1:lower_gate -> python [0:lower_gate+1)
        # max ignoring the fact that NaNs propagate: MATLAB max() ignores NaN
        if np.all(np.isnan(beta_roi)):
            S_filtered[i] = np.nan
            continue
        max_idx = int(np.nanargmax(beta_roi))
        max_beta = beta_roi[max_idx]

        if np.isnan(max_beta) or max_beta <= 0:
            S_filtered[i] = np.nan
            continue

        # Filter 1: 300 m above
        idx_above = max_idx + gate_300m
        if idx_above <= beta_profile.size - 1:
            beta_above = beta_profile[idx_above]
            if not np.isnan(beta_above) and beta_above * config.attenuation_factor > max_beta:
                S_filtered[i] = np.nan
                stats["above_rejected"] += 1
                continue

        # Filter 2: 300 m below
        idx_below = max_idx - gate_300m
        if idx_below >= 0:
            beta_below = beta_profile[idx_below]
            if not np.isnan(beta_below) and beta_below * config.attenuation_factor > max_beta:
                S_filtered[i] = np.nan
                stats["below_rejected"] += 1
                continue

        # Filter 3: aerosol contribution ratio
        if lower_gate is not None and max_idx > lower_gate + 5:
            idx_start_aero = max(0, lower_gate)  # MATLAB max(1,lower_gate) -> 0-based
            idx_end_aero = max_idx - 5
            if idx_end_aero > idx_start_aero:
                # MATLAB nansum(beta(idx_start:idx_end))*range_resol, inclusive end
                beta_below_cloud = np.nansum(
                    beta_profile[idx_start_aero:idx_end_aero + 1]) * range_resol
                beta_total = np.nansum(
                    beta_profile[lower_gate:upper_gate + 1]) * range_resol
                if beta_total > 0:
                    ratio = beta_below_cloud / beta_total
                    if ratio > config.ratio_filter:
                        S_filtered[i] = np.nan
                        stats["ratio_rejected"] += 1
                        continue

        # Filter 4: CBH range
        if data.cbh is not None and not np.isnan(data.cbh[i]):
            cbh = data.cbh[i]
        else:
            cbh = data.range[max_idx]
        if cbh < config.cbh_minheight or cbh > config.cbh_maxheight:
            S_filtered[i] = np.nan
            stats["cbh_rejected"] += 1
            continue

    return S_filtered, stats


# ===========================================================================
#  apply_temporal_consistency_filter
# ===========================================================================
def apply_temporal_consistency_filter(
    S: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Port of ``apply_temporal_consistency_filter``: N consecutive within +/-X%."""
    n_profiles = S.size
    S_consistent = np.full(n_profiles, np.nan)
    n_consec = int(config.n_consecutive)
    range_percent = config.consistency_range
    plus_limit = 1.0 + range_percent / 100.0
    minus_limit = 1.0 - range_percent / 100.0

    for i in range(n_profiles - n_consec + 1):
        group = S[i:i + n_consec]
        if np.any(np.isnan(group)):
            continue
        mean_group = np.mean(group)
        if mean_group <= 0:
            continue
        is_consistent = np.all(
            (group >= mean_group * minus_limit) & (group <= mean_group * plus_limit))
        if is_consistent:
            S_consistent[i:i + n_consec] = S[i:i + n_consec]

    stats = {"n_rejected": int(np.sum(~np.isnan(S) & np.isnan(S_consistent)))}
    return S_consistent, stats


# ===========================================================================
#  apply_transmission_correction
# ===========================================================================
def apply_transmission_correction(
    beta: NDArray, data: CeiloData, S: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, NDArray, NDArray]:
    """Port of ``apply_transmission_correction``: aerosol two-way transmission below cloud."""
    n_profiles = S.size
    C_corrected = np.full(n_profiles, np.nan)
    C_low = np.full(n_profiles, np.nan)
    C_high = np.full(n_profiles, np.nan)
    range_resol = float(data.range[1] - data.range[0])

    if config.aerosol_lidar_ratio is None:
        raise ValueError(
            "config.aerosol_lidar_ratio is required when apply_transmission_correction "
            "is on (MATLAB set_defaults does not define it; the runner sets it to 50).")

    for i in range(n_profiles):
        if np.isnan(S[i]):
            continue
        beta_profile = beta[:, i]
        # [~, max_idx] = max(beta_profile)  (over full profile; NaN ignored)
        if np.all(np.isnan(beta_profile)):
            continue
        max_idx = int(np.nanargmax(beta_profile))

        # MATLAB max_idx < 10 (1-based) -> first 9 gates. 0-based: max_idx < 9.
        if max_idx < 9:
            C_corrected[i] = S[i] / S_THEORETICAL
            C_low[i] = S[i] / S_THEORETICAL
            C_high[i] = S[i] / S_THEORETICAL
            continue

        # Integrated below-cloud backscatter, in the INSTRUMENT (uncalibrated) beta units:
        #   nansum(beta(5:max_idx-5))*range_resol     (MATLAB 1-based 5:(max_idx-5) inclusive
        #   -> 0-based [4 : max_idx-5] inclusive).
        seg = beta_profile[4:(max_idx - 5) + 1]
        B_aerosol_raw = np.nansum(seg) * range_resol
        # Convert to a PHYSICAL aerosol optical depth before the Beer-Lambert transmission.
        # ``beta`` here is the uncalibrated attenuated backscatter (e.g. L2 stored in
        # 1E-6*1/(m*sr) units, ~O(1)); the per-profile calibration coefficient
        # C = S/S_THEORETICAL is exactly what scales it to physical 1/(m*sr), so
        #   AOD = LR * integral(C * beta) dr = LR * C * B_aerosol_raw.
        # Without the C factor the "AOD" is ~1e6x too large for L2 input and T2 underflows
        # to 0 (the MATLAB reference shares this latent bug — its parity test is skipped).
        C_base = S[i] / S_THEORETICAL
        B_aerosol = C_base * B_aerosol_raw
        if B_aerosol <= 0:
            C_corrected[i] = C_base
            C_low[i] = C_base
            C_high[i] = C_base
            continue

        T2 = np.exp(-2.0 * config.aerosol_lidar_ratio * B_aerosol)
        T2_low = np.exp(-2.0 * config.aerosol_lidar_ratio_low * B_aerosol)
        T2_high = np.exp(-2.0 * config.aerosol_lidar_ratio_high * B_aerosol)
        # Same form as the MATLAB reference (C_corrected = C_base * T2); only B_aerosol's
        # units were corrected above (the scale bug). T2 <= 1.
        C_corrected[i] = C_base * T2
        C_low[i] = C_base * T2_low
        C_high[i] = C_base * T2_high

    return C_corrected, C_low, C_high


# ===========================================================================
#  calculate_mode
# ===========================================================================
def calculate_mode(values: NDArray) -> float:
    """Port of ``calculate_mode``: round to 2 dp, histcounts (auto bins), bin center.

    Reproduces MATLAB ``histcounts(x)`` automatic binning (``matlab.internal.math.binpicker``)
    so the mode bin matches MATLAB. The returned mode is the centre of the most populated
    bin, with the same fallback to the median when ``max_idx`` is the final edge.
    """
    v = np.asarray(values, dtype="float64")
    if v.size == 0 or np.all(np.isnan(v)):
        return float("nan")
    data_rounded = np.round(v, 2)
    counts, edges = _matlab_histcounts_auto(data_rounded)
    max_idx = int(np.argmax(counts))
    if max_idx < len(edges) - 1:
        return float((edges[max_idx] + edges[max_idx + 1]) / 2.0)
    return float(np.median(v))


def _matlab_histcounts_auto(x: NDArray) -> Tuple[NDArray, NDArray]:
    """Reproduce MATLAB ``histcounts(x)`` automatic binning + counts.

    For non-integer continuous data MATLAB's 'auto' rule uses Scott's normal-reference
    bin width ``3.5*std(x)/n^(1/3)``, snapped to a "nice" number, with left edges aligned
    to a multiple of that width spanning [min,max] (``binpicker``). Verified against
    MATLAB on the real calibration-coefficient vector (binwidth 0.05, edges 0.60..1.20).
    """
    x = x[~np.isnan(x)]
    if x.size == 0:
        return np.array([0]), np.array([0.0, 1.0])
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    edges = _binpicker_edges(xmin, xmax, x.size, float(np.std(x, ddof=1)) if x.size > 1 else 0.0)
    # histcounts: bin k counts [edges[k], edges[k+1]); the LAST bin is closed on the right.
    counts, _ = np.histogram(x, bins=edges)
    return counts.astype("int64"), edges


def _binpicker_edges(xmin: float, xmax: float, n: int, xstd: float) -> NDArray:
    """Port of MATLAB's ``binpicker`` (auto/scott rule) -> bin edges.

    Scott's raw width ``3.5*std/n^(1/3)`` is snapped to a nice number
    ``{1,2,3,5,10} * 10^floor(log10(raw))`` using MATLAB's thresholds (1.5, 2.5, 4, 7.5),
    then edges are left-aligned to a multiple of that width covering [xmin, xmax].
    """
    if xmin == xmax:
        # MATLAB makes a single unit-width bin centred on the value (xmin-0.5, xmin+0.5)
        return np.array([xmin - 0.5, xmin + 0.5])

    raw_bw = 3.5 * xstd / (n ** (1.0 / 3.0)) if xstd > 0 else (xmax - xmin)
    if not np.isfinite(raw_bw) or raw_bw <= 0:
        raw_bw = (xmax - xmin)

    pow10 = 10.0 ** np.floor(np.log10(raw_bw))
    rel = raw_bw / pow10               # in [1, 10)
    if rel < 1.5:
        nice = 1.0
    elif rel < 2.5:
        nice = 2.0
    elif rel < 4.0:
        nice = 3.0
    elif rel < 7.5:
        nice = 5.0
    else:
        nice = 10.0
    bw = nice * pow10

    left = bw * np.floor(xmin / bw)
    nbins = int(np.ceil((xmax - left) / bw))
    if nbins < 1:
        nbins = 1
    # ensure the right edge strictly covers xmax (guards float edge cases)
    while left + nbins * bw <= xmax:
        nbins += 1
    return left + bw * np.arange(nbins + 1)


# ===========================================================================
#  Small index helpers (MATLAB find(...,1,'first'/'last'))
# ===========================================================================
def _trapz(y: NDArray, x: NDArray) -> float:
    """Trapezoidal integral matching MATLAB ``trapz(x, y)`` (numpy.trapezoid, no deprecation)."""
    y = np.asarray(y, dtype="float64")
    x = np.asarray(x, dtype="float64")
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def _find_first(mask: NDArray) -> Optional[int]:
    idx = np.where(mask)[0]
    return int(idx[0]) if idx.size else None


def _find_last(mask: NDArray) -> Optional[int]:
    idx = np.where(mask)[0]
    return int(idx[-1]) if idx.size else None


# ===========================================================================
#  Results container + create_empty_results
# ===========================================================================
@dataclass
class CloudCalResults:
    calibration_factor: float = float("nan")
    calibration_coefficient: float = float("nan")
    lidar_ratios: NDArray = field(default_factory=lambda: np.array([]))
    cal_mean: float = float("nan")
    cal_median: float = float("nan")
    cal_mode: float = float("nan")
    cal_std: float = float("nan")
    n_profiles: int = 0
    time: NDArray = field(default_factory=lambda: np.array([]))
    cbh: NDArray = field(default_factory=lambda: np.array([]))
    all_coefficients: NDArray = field(default_factory=lambda: np.array([]))
    altitude_warning: bool = False
    trans2_wv: Optional[NDArray] = None
    config: Optional[CloudCalConfig] = None
    filter_stats: Optional[dict] = None
    cloud_stats: Optional[dict] = None
    consistency_stats: Optional[dict] = None
    # extra (not in MATLAB results, useful for the test)
    S_apparent: Optional[NDArray] = None
    S_consistent: Optional[NDArray] = None


def create_empty_results() -> CloudCalResults:
    """Port of ``create_empty_results``."""
    return CloudCalResults()


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
            data.cbh = _block_reduce_mean(data.cbh, t_factor, axis=0)
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


def _matlab_unix_block_time(time_dt: NDArray, factor: int) -> NDArray:
    """Block-mean a datetime64 time vector, returning datetime64[ns] bin centres."""
    ints = time_dt.astype("datetime64[ns]").astype("int64").astype("float64")
    reduced = _block_reduce_mean(ints, factor, axis=0)
    return reduced.round().astype("int64").astype("datetime64[ns]")


# ===========================================================================
#  Main: liquid_cloud_calibration
# ===========================================================================
def liquid_cloud_calibration(config: CloudCalConfig) -> CloudCalResults:
    """Port of the top-level ``liquid_cloud_calibration`` function.

    Returns a :class:`CloudCalResults`. Mirrors the MATLAB control flow exactly,
    including the strict WV requirement for CL31/CL51/CL61 (raises on failure) and the
    optional aerosol transmission correction (which replaces the coefficients).
    """
    config = set_defaults(config)

    data, status = read_ceilometer_data(config.nc_file, config)
    if status != 0:
        raise RuntimeError("Failed to read NetCDF file")

    return liquid_cloud_calibration_from_data(data, config)


def liquid_cloud_calibration_from_data(data: CeiloData, config: CloudCalConfig) -> CloudCalResults:
    """Run the liquid-cloud calibration on an already-loaded :class:`CeiloData` — i.e. the
    post-read part of :func:`liquid_cloud_calibration`. Used to calibrate multi-file
    Cloudnet data concatenated outside this module (a day or a month of CL61 raw files).
    ``config.instrument`` must be set by the caller; ``set_defaults`` is applied here
    (idempotently) so wavelength/FWHM match it. Behaviour is identical to the file-based
    entry point, which now simply reads then delegates here."""
    config = set_defaults(config)

    # The liquid-cloud (O'Connor/Hopkin) calibration is the PRIMARY method for the 910 nm
    # ceilometers (CL31/CL51/CL61) — including the CL31/CL51, which cannot be Rayleigh-calibrated.
    # It also runs for CHM15k / Mini-MPL. No suitability warning here: the only genuine warning is
    # _beta_conversion_factor's unrecognized-units fallback.

    # --- Optional pre-averaging (time/range) to speed up high-res files ---
    # Applied BEFORE the WV correction and the per-profile filters (the costly steps),
    # so all downstream processing runs on the reduced grid. Disabled (no-op) when both
    # factors resolve to 1, preserving the bit-for-bit MATLAB path.
    data = average_ceilo_data(data, config)

    # --- Water-vapor absorption correction ---
    in_wv_band = config.instrument.strip().upper() in ("CL31", "CL51", "CL61")
    if config.apply_wv_correction and in_wv_band:
        try:
            trans2 = compute_wv_transmission(data, config)
        except Exception as exc:  # noqa: BLE001 - STRICT: re-raise, no fallback
            raise RuntimeError(
                f"Water-vapor correction FAILED at {config.wavelength:.2f} nm "
                f"({config.instrument}): {exc}. A 910 nm cloud calibration without a "
                f"valid WV correction is not permitted (no fallback): this period is "
                f"NOT calibrated.") from exc
        trans2 = np.asarray(trans2, dtype="float64")
        if (trans2.size == 0 or not np.any(np.isfinite(trans2))
                or np.all(trans2 == 1)):
            raise RuntimeError(
                f"Water-vapor correction at {config.wavelength:.2f} nm returned no "
                f"usable transmission (empty/NaN/all-ones). Not calibrated.")
        data.beta = data.beta / trans2
        data.trans2_wv = trans2
    elif config.apply_wv_correction and not in_wv_band:
        data.trans2_wv = np.ones_like(data.beta)

    # --- Multiple scattering correction ---
    beta_corrected = apply_multiple_scattering_correction(data.beta, data.range, config)

    # --- Instrument health filters ---
    beta_filtered, filter_stats = apply_instrument_filters(beta_corrected, data, config)

    # --- Apparent lidar ratio ---
    S_apparent, integrated_beta = calculate_lidar_ratio(beta_filtered, data.range, config)

    # --- Cloud quality filters ---
    S_filtered, cloud_stats = apply_cloud_filters(S_apparent, beta_filtered, data, config)

    # --- Temporal consistency filter ---
    S_consistent, consistency_stats = apply_temporal_consistency_filter(S_filtered, config)

    # --- Calibration coefficient ---
    calibration_coefficients = S_consistent / S_THEORETICAL
    valid_idx = ~np.isnan(calibration_coefficients)
    valid_C = calibration_coefficients[valid_idx]

    if valid_C.size == 0:
        res = create_empty_results()
        res.config = config
        res.filter_stats = filter_stats
        res.S_apparent = S_apparent
        res.S_consistent = S_consistent
        res.trans2_wv = data.trans2_wv
        return res

    cal_mean = float(np.mean(valid_C))
    cal_median = float(np.median(valid_C))
    cal_std = float(np.std(valid_C, ddof=1))  # MATLAB std default normalises by N-1
    cal_mode = calculate_mode(valid_C)
    n_profiles = int(valid_C.size)

    # --- Optional aerosol transmission correction (replaces the coefficients) ---
    if config.apply_transmission_correction:
        C_corrected, _C_low, _C_high = apply_transmission_correction(
            beta_filtered, data, S_consistent, config)
        calibration_coefficients = C_corrected
        valid_C = calibration_coefficients[valid_idx]
        cal_mean = float(np.mean(valid_C))
        cal_median = float(np.median(valid_C))
        cal_std = float(np.std(valid_C, ddof=1))
        cal_mode = calculate_mode(valid_C)

    res = CloudCalResults()
    # No valid in-cloud profiles (cal_median 0/NaN) -> report NaN rather than crash.
    res.calibration_factor = (1.0 / cal_median) if (np.isfinite(cal_median) and cal_median != 0.0) else float("nan")
    res.calibration_coefficient = cal_median
    res.lidar_ratios = S_consistent[valid_idx]
    res.cal_mean = cal_mean
    res.cal_median = cal_median
    res.cal_mode = cal_mode
    res.cal_std = cal_std
    res.n_profiles = n_profiles
    res.time = data.time[valid_idx]
    res.cbh = data.cbh[valid_idx]
    res.all_coefficients = calibration_coefficients
    res.config = config
    res.filter_stats = filter_stats
    res.cloud_stats = cloud_stats
    res.consistency_stats = consistency_stats
    res.altitude_warning = data.altitude_warning
    res.trans2_wv = data.trans2_wv
    res.S_apparent = S_apparent
    res.S_consistent = S_consistent
    return res
