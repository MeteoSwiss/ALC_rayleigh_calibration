"""Water-vapour two-way transmission for the liquid-cloud calibration (910 nm band).

Split out of ``cloud/calibration.py`` and co-located here with the other water-vapour code
(2026-07). It is the cloud method's OWN transmission chain: the humidity physics (Murphy-Koop /
Wagner-Pruss / number density), the CAMS & ERA5 model-level readers, and ``compute_wv_transmission``,
built on the shared Gaussian core (``wv_t2eff_core``) + L137 coefficients in ``water_vapor.py``.
Functions duck-type ``data``/``config`` (type hints only, stringised by ``from __future__``), so this
module does not import back into the cloud package. Re-exported from ``cloud.calibration`` for
back-compat. (It historically diverged from the shared route for MATLAB parity; now that MATLAB is
no longer the reference the two water-vapour paths could be reconciled — a future cleanup.)
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset

from .water_vapor import (
    wv_t2eff_core, load_abs_cross_section, in_water_vapor_band, _A137, _B137)
from ..io.cams import ensure_cams_file

_RD_CAMS = 287.06                        # J/(kg K), Rd in get_Beta_CAMS_oper_monthly.m


_G0 = 9.80665                            # m/s^2 (get_Beta_CAMS_oper_monthly.m)


_M_DRY = 28.9644                         # kg/kmol (get_meteo_const.m)


_M_WET = 18.0152                         # kg/kmol (get_meteo_const.m)


_C_MD = _M_WET / _M_DRY                  # = 0.6219808... (convert_humidity c)


_RW = 0.4615                             # J g^-1 K^-1 (get_water_vapor..RH.m)


_NW_COEF = 7.25e22                       # number-density coefficient (get_water_vapor..RH.m)


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


@lru_cache(maxsize=32)
def _cams_levels_all_times(
    cams_file: str, latitude: float, longitude: float,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """Read CAMS T/q at the nearest grid point for ALL time steps and integrate the
    L137 half-level geopotential, exactly like get_Beta_CAMS_oper_monthly.m (z=[] path).

    Cached by (file, lat, lon): a station processes many days from the same monthly CAMS
    file at one location, so without the cache the (tens-of-MB) file is re-read and the
    hydrostatic integration redone for every single day (~8 s each). The returned arrays
    are treated as read-only by ``compute_wv_transmission`` (it only interpolates from
    them), so sharing them across days is safe.

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
        if "station" in nc.dimensions:
            # Per-station MARS format (mars_wv_station_v1): t/q are (time, station, level),
            # lnsp (time, station), z (station) static orography. Pick the nearest station
            # (an exact match) and read its column; the L137 integration below is identical.
            lat_m = np.asarray(nc.variables["lat"][:], dtype="float64")
            lon_m = np.asarray(nc.variables["lon"][:], dtype="float64")
            level = np.asarray(nc.variables["level"][:], dtype="float64")
            t_cf = np.asarray(nc.variables["time"][:], dtype="float64")   # seconds since 1970-01-01
            dlon = ((lon_m - longitude + 180.0) % 360.0) - 180.0
            si = int(np.argmin((lat_m - latitude) ** 2 + dlon ** 2))
            n_lev = level.size
            n_t = t_cf.size
            T = np.asarray(nc.variables["t"][:, si, :], dtype="float64").T           # (level, time)
            q = np.asarray(nc.variables["q"][:, si, :], dtype="float64").T
            lnsp = np.asarray(nc.variables["lnsp"][:, si], dtype="float64").ravel()  # (time,)
            z_surf = np.full(n_t, float(np.asarray(nc.variables["z"][si])), dtype="float64")  # static
            time_num = t_cf / 86400.0 + 719529.0   # datenum(1970,1,1) = 719529
        else:
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
                # Slice the NetCDF variable DIRECTLY so only this station's (level, time) column
                # is read from disk. Never materialise the full 4-D (time, level, lat, lon) array
                # via var[:] -- that is multi-GB for a monthly CAMS file and was the cause of the
                # ~10 GB/worker RSS growth (freed numpy memory is not returned to the OS).
                axes = {nm: k for k, nm in enumerate(dim_names)}
                idx = [slice(None)] * len(dim_names)
                idx[axes["longitude"]] = li
                idx[axes["latitude"]] = ai
                sub = np.asarray(var[tuple(idx)], dtype="float64")  # remaining dims: level, time (orig order)
                remaining = [nm for nm in dim_names if nm in ("level", "time")]
                lev_pos = remaining.index("level")
                sub = np.moveaxis(sub, lev_pos, 0)  # (level, time)
                return sub

            T = _read_profile(tvar)
            q = _read_profile(qvar)

            # z, lnsp are surface fields stored at a single (top) level slot.
            # MATLAB reads them at level index 1 (start=1,count=1) for all times.
            def _read_surface(var):
                axes = {nm: k for k, nm in enumerate(dim_names)}
                idx = [slice(None)] * len(dim_names)
                idx[axes["longitude"]] = li
                idx[axes["latitude"]] = ai
                idx[axes["level"]] = 0  # MATLAB start index 1 -> python 0 (first level slot)
                sub = np.asarray(var[tuple(idx)], dtype="float64")  # read only the column; remaining dim: time
                return sub.ravel()

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


_ERA5_CACHE_DATA: dict = {}


def _era5_levels_all_times(era5_cache: str, latitude: float, longitude: float):
    """Per-step ERA5 humidity profile at the station nearest (lat, lon), from a prefetched cache.

    Mirrors the return contract of :func:`_cams_levels_all_times` so ``compute_wv_transmission`` is
    source-agnostic:

        time_num : (n_t,)          MATLAB datenum (days; datenum(1970,1,1)=719529)
        z_asl    : (n_lev, n_t)    geopotential height [m ASL]
        T        : (n_lev, n_t)    temperature [K]
        nw       : (n_lev, n_t)    water-vapour number density [m^-3]

    The cache (built offline by ``scripts/prefetch_era5_edh.py`` from the DestinE Earth Data Hub)
    is a NetCDF with dims (station, time, level), coords ``lat``/``lon`` per station, ``time``
    (datetime64), ``level`` (pressure [hPa]), and variables ``q`` [kg/kg], ``t`` [K],
    ``z`` [geopotential m^2/s^2]. Cached in-process per file (one station stream reads it once).
    ERA5-only policy: a missing cache, a too-far nearest station (> ~0.75 deg), or an all-NaN
    profile RAISES, so the caller flags the night uncalibrated (flag -4) — never a WV-free run.
    """
    from ..water_vapor_correction.water_vapor import KB, EPS, G0
    if not era5_cache or not Path(era5_cache).is_file():
        raise FileNotFoundError(
            f"ERA5 water-vapour cache not found: {era5_cache!r} (wv_source='era5' needs the "
            f"prefetched cache from scripts/prefetch_era5_edh.py)")
    key = str(era5_cache)
    cached = _ERA5_CACHE_DATA.get(key)
    if cached is None:
        import xarray as xr
        with xr.open_dataset(era5_cache) as ds:
            cached = (
                np.asarray(ds["lat"].values, float), np.asarray(ds["lon"].values, float),
                np.asarray(ds["level"].values, float), np.asarray(ds["time"].values),
                np.asarray(ds["q"].values, float), np.asarray(ds["t"].values, float),
                np.asarray(ds["z"].values, float))
        _ERA5_CACHE_DATA[key] = cached
    slat, slon, level, times, q, T, z = cached
    dlon = ((slon - longitude + 180.0) % 360.0) - 180.0
    dist = np.sqrt((slat - latitude) ** 2 + dlon ** 2)
    si = int(np.argmin(dist))
    if dist[si] > 0.75:
        raise ValueError(
            f"No ERA5 cache profile near station ({latitude:.2f},{longitude:.2f}); "
            f"nearest cached point {dist[si]:.2f} deg away")
    qs, Ts, zs = q[si], T[si], z[si]              # each (time, level)
    if not np.any(np.isfinite(qs) & np.isfinite(Ts) & np.isfinite(zs)):
        raise ValueError(
            f"ERA5 cache profile all-NaN at station ({latitude:.2f},{longitude:.2f})")
    P = level[None, :] * 100.0                     # (1, level) Pa
    Pw = qs * P / (EPS + (1.0 - EPS) * qs)         # water-vapour partial pressure [Pa]
    nw = Pw / (KB * Ts)                            # number density [m^-3]
    height = zs / G0                               # geopotential height [m ASL]
    dn = times.astype("datetime64[s]").astype("float64") / 86400.0 + 719529.0
    return dn, height.T, Ts.T, nw.T                # (n_t,), (n_lev,n_t) x3


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

    # --- Humidity source (T/q/z -> nw and geopotential height), all time steps ---
    # The rest of this routine is source-agnostic: it only needs the per-step model times
    # (MATLAB datenum), geopotential HEIGHT [m ASL] and water-vapour number density [m^-3].
    if str(getattr(config, "wv_source", "cams")).lower() == "era5":
        # Prefetched ERA5 profile cache (all stations, pressure levels; built offline from the
        # DestinE Earth Data Hub). ERA5-only: a station/time with no usable profile raises ->
        # flag -4 (not calibrated), exactly like a missing CAMS WV correction.
        time_cams, cams_z, _cams_T, nw_all = _era5_levels_all_times(
            config.era5_cache, data.station_latitude, data.station_longitude)
    else:
        _cams_kw = dict(
            auto_download=getattr(config, "auto_download_cams", False),
            scope=getattr(config, "cams_download_scope", "day"),
            latitude=data.station_latitude, longitude=data.station_longitude,
        )
        cams_path = ensure_cams_file(config.cams_folder, cams_date, **_cams_kw)
        if cams_path is None and getattr(config, "cams_folder_fallback", ""):
            # Month absent from the primary archive (e.g. the 0.4 deg set) -> fall back to the coarser
            # 1 deg archive (also L137 model levels), so the month still calibrates instead of -4.
            cams_path = ensure_cams_file(config.cams_folder_fallback, cams_date, **_cams_kw)
        if cams_path is None:
            raise FileNotFoundError(
                f"No CAMS file for {cams_date} in {config.cams_folder} "
                f"(nor fallback {config.cams_folder_fallback or '-'}); looked for monthly "
                f"CAMS_Beta_{cams_date[:6]}.nc or daily CAMS_Beta_{cams_date[:8]}.nc")
        cams_file = str(cams_path)
        from ..water_vapor_correction.water_vapor import cams_point_too_far
        if cams_point_too_far(cams_file, data.station_latitude, data.station_longitude):
            # Regional CAMS has no data at this station (e.g. New Zealand): the nearest grid point is
            # at the domain edge, thousands of km away -> a bogus WV correction. Fail loudly so the
            # runner emits flag -10 ('Closest CAMS data too far') instead of a false success.
            raise ValueError(
                f"Closest CAMS data too far from station "
                f"({data.station_latitude:.2f},{data.station_longitude:.2f}); station outside CAMS domain")
        time_cams, cams_z, _cams_T, nw_all = _cams_levels_all_times(
            cams_file, data.station_latitude, data.station_longitude)

    # A monthly CAMS file holds ~248 (3-hourly) steps, but one daily ceilometer file spans
    # ~1 day. The nearest-time interpolation at the end only ever selects the CAMS steps that
    # bracket the data, so restrict to that window first (the per-CAMS-step WV integration is
    # the dominant cost). Using the same float representation as the final interp keeps it
    # exact regardless of the time units, and the +/-1 step pad preserves nearest-extrapolation.
    if time_cams.size > 1:
        cdn = time_cams.astype("float64")
        edn = data.time_num.astype("float64")
        t0, t1 = float(edn.min()), float(edn.max())
        inside = np.where((cdn >= t0) & (cdn <= t1))[0]
        if inside.size:
            lo, hi = int(inside[0]), int(inside[-1])
        else:  # data falls between two CAMS steps -> keep the single nearest step
            lo = hi = int(np.argmin(np.abs(cdn - t0)))
        lo = max(lo - 1, 0)
        hi = min(hi + 1, time_cams.size - 1)
        sl = slice(lo, hi + 1)
        time_cams, cams_z, nw_all = time_cams[sl], cams_z[:, sl], nw_all[:, sl]

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

    # Two-way transmission for ALL CAMS steps at once (matrix form of wv_t2eff.m's per-step
    # loop). Build the in-band extinction (n_band, n_range, n_cams), cumulative-trapezoid the
    # optical depth along range, then Gaussian-average exp(-2*tau). Out-of-band wavelengths
    # carry zero optical depth (transmission 1) and enter only via the Gaussian normalization,
    # i.e. as the constant gauss[~band].sum() term -- exactly as the loop's sum over all wl.
    ab = abs_cs[band, :]                                          # (n_band, n_range)
    ext = ab[:, :, None] * (wv_density_all[None, :, :] / 1e4)     # (n_band, n_range, n_cams) [m^-1]
    dx = np.diff(range_col)
    incr = 0.5 * (ext[:, 1:, :] + ext[:, :-1, :]) * dx[None, :, None]
    tau = np.concatenate([np.zeros((ext.shape[0], 1, n_cams)),
                          np.cumsum(incr, axis=1)], axis=1)        # (n_band, n_range, n_cams)
    weighted = np.exp(-2.0 * tau) * gauss[band][:, None, None]     # trans^2 * gauss, in band
    trans2_cams = (weighted.sum(axis=0) + float(gauss[~band].sum())) / sum_gauss  # (n_range, n_cams)

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
