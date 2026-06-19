"""
Water-vapor two-way transmission for the Rayleigh calibration of 910 nm ALCs.

At 905-911 nm, water-vapor absorption attenuates the lidar signal, including the
molecular-reference region used for the Rayleigh calibration. Without removing it
the derived lidar constant is biased low. This module computes the spectrally
averaged two-way water-vapor transmission T2_wv(range) so the molecular model
used in the lidar-constant fit can be corrected. CHM15k (1064 nm) is outside the
absorption band and needs no correction.

Faithful Python port of the validated MATLAB routines:
  wv_t2eff.m, compute_wv_transmission.m, get_water_vapor_number_concentration_from_RH.m
(Wiegner & Gasteiger 2015). Humidity is taken from CAMS specific humidity q.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

# Physical constants
KB = 1.380649e-23     # Boltzmann constant [J/K]
RD = 287.058          # specific gas constant, dry air [J/(kg K)]
G0 = 9.80665          # standard gravity [m/s^2]
EPS = 0.621981        # Mw/Md (ratio of molar masses water/dry air)

# Per-instrument laser spectrum (Payerne Qmini campaign 2026-06-02).
# central wavelength [nm], spectral FWHM [nm]. Used for the Gaussian weighting.
LASER_SPECTRUM = {
    "CL31":   (909.7, 6.0),
    "CL51":   (910.0, 3.4),
    "CL61":   (910.74, 1.0),
    "CHM15k": (1064.47, 0.5),
    "CHM8k":  (1064.0, 0.5),
}


def in_water_vapor_band(wavelength_nm: float) -> bool:
    """True if the wavelength is inside the 900-920 nm water-vapor absorption band."""
    return 900.0 <= wavelength_nm <= 920.0


def laser_spectrum_for(instrument_type: str, fallback_nm: float) -> Tuple[float, float]:
    """Return (lambda0_nm, fwhm_nm) for an instrument type."""
    if instrument_type in LASER_SPECTRUM:
        return LASER_SPECTRUM[instrument_type]
    return (fallback_nm, 3.4)


# 910 nm water-vapour absorption LUT shipped as package data: the 903-918 nm band of
# abs_cross_647_full_levels_1000.nc (HITRAN/MT-CKD), abscs_ave only. Used when no explicit
# LUT path is configured, so the WV correction works out of the box for 910 nm instruments.
DEFAULT_ABS_CROSS_SECTION = (
    Path(__file__).resolve().parent.parent / "data" / "abs_cross_wv_910nm.nc"
)


def load_abs_cross_section(lut_path: Optional[Path] = None) -> Tuple[NDArray, NDArray, NDArray]:
    """
    Load the HITRAN/MT-CKD absorption cross-section LUT.

    With no path (or a missing / empty one) the bundled 910 nm band LUT
    (:data:`DEFAULT_ABS_CROSS_SECTION`) is used, so the water-vapour correction
    runs without any external file.

    Returns
    -------
    wl_nm   : (n_wl,) wavelength grid [nm]
    height_m: (n_height,) height grid [m]
    abscs   : (n_wl, n_height) absorption cross-section [cm^2]
    """
    if lut_path is None or str(lut_path).strip() in ("", ".") or not Path(lut_path).exists():
        lut_path = DEFAULT_ABS_CROSS_SECTION
    from netCDF4 import Dataset
    with Dataset(lut_path) as nc:
        wl_nm = np.asarray(nc.variables["lambda"][:], dtype=float)
        height_m = np.asarray(nc.variables["height_in_km"][:], dtype=float) * 1000.0
        abscs = np.asarray(nc.variables["abscs_ave"][:], dtype=float)
    # Ensure abscs is [n_wl, n_height]
    if abscs.shape[0] != wl_nm.size and abscs.shape[1] == wl_nm.size:
        abscs = abscs.T
    return wl_nm, height_m, abscs


# ECMWF L137 half-level coefficients (a [Pa], b [-]) for levels 0..137.
# Pressure at half level k: Ph(k) = A[k] + B[k] * surface_pressure.
# (Port of param_137_levels in get_Beta_CAMS_oper_monthly.m.)
_A137 = np.array([
    0.0, 2.000365, 3.102241, 4.666084, 6.827977, 9.746966, 13.605424, 18.608931,
    24.985718, 32.98571, 42.879242, 54.955463, 69.520576, 86.895882, 107.415741,
    131.425507, 159.279404, 191.338562, 227.968948, 269.539581, 316.420746,
    368.982361, 427.592499, 492.616028, 564.413452, 643.339905, 729.744141,
    823.967834, 926.34491, 1037.201172, 1156.853638, 1285.610352, 1423.770142,
    1571.622925, 1729.448975, 1897.519287, 2076.095947, 2265.431641, 2465.770508,
    2677.348145, 2900.391357, 3135.119385, 3381.743652, 3640.468262, 3911.490479,
    4194.930664, 4490.817383, 4799.149414, 5119.89502, 5452.990723, 5798.344727,
    6156.074219, 6526.946777, 6911.870605, 7311.869141, 7727.412109, 8159.354004,
    8608.525391, 9076.400391, 9562.682617, 10065.97852, 10584.63184, 11116.66211,
    11660.06738, 12211.54785, 12766.87305, 13324.66895, 13881.33106, 14432.13965,
    14975.61523, 15508.25684, 16026.11523, 16527.32227, 17008.78906, 17467.61328,
    17901.62109, 18308.43359, 18685.71875, 19031.28906, 19343.51172, 19620.04297,
    19859.39063, 20059.93164, 20219.66406, 20337.86328, 20412.30859, 20442.07813,
    20425.71875, 20361.81641, 20249.51172, 20087.08594, 19874.02539, 19608.57227,
    19290.22656, 18917.46094, 18489.70703, 18006.92578, 17471.83984, 16888.6875,
    16262.04688, 15596.69531, 14898.45313, 14173.32422, 13427.76953, 12668.25781,
    11901.33984, 11133.30469, 10370.17578, 9617.515625, 8880.453125, 8163.375,
    7470.34375, 6804.421875, 6168.53125, 5564.382813, 4993.796875, 4457.375,
    3955.960938, 3489.234375, 3057.265625, 2659.140625, 2294.242188, 1961.5,
    1659.476563, 1387.546875, 1143.25, 926.507813, 734.992188, 568.0625,
    424.414063, 302.476563, 202.484375, 122.101563, 62.78125, 22.835938,
    3.757813, 0.0, 0.0])
_B137 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7e-06, 2.4e-05, 5.9e-05,
    0.000112, 0.000199, 0.00034, 0.000562, 0.00089, 0.001353, 0.001992, 0.002857,
    0.003971, 0.005378, 0.007133, 0.009261, 0.011806, 0.014816, 0.018318,
    0.022355, 0.026964, 0.032176, 0.038026, 0.044548, 0.051773, 0.059728,
    0.068448, 0.077958, 0.088286, 0.099462, 0.111505, 0.124448, 0.138313,
    0.153125, 0.16891, 0.185689, 0.203491, 0.222333, 0.242244, 0.263242,
    0.285354, 0.308598, 0.332939, 0.358254, 0.384363, 0.411125, 0.438391,
    0.466003, 0.4938, 0.521619, 0.549301, 0.576692, 0.603648, 0.630036,
    0.655736, 0.680643, 0.704669, 0.727739, 0.749797, 0.770798, 0.790717,
    0.809536, 0.827256, 0.843881, 0.859432, 0.873929, 0.887408, 0.8999,
    0.911448, 0.922096, 0.931881, 0.94086, 0.949064, 0.95655, 0.963352,
    0.969513, 0.975078, 0.980072, 0.984542, 0.9885, 0.991984, 0.995003,
    0.99763, 1.0])
_RD_CAMS = 287.06  # J/(kg K), as in get_Beta_CAMS_oper_monthly.m


def _cams_levels(
    cams_file: Path,
    latitude: float,
    longitude: float,
    t_start: np.datetime64,
    t_end: np.datetime64,
) -> Optional[Tuple[NDArray, NDArray, NDArray, NDArray]]:
    """
    Core CAMS model-level read + hydrostatic integration over a time window at the
    nearest grid point. Shared by the water-vapor correction and the optional CAMS
    molecular profile.

    Faithful port of get_Beta_CAMS_oper_monthly.m: model-level pressure from the
    ECMWF L137 a/b coefficients (CAMS z/lnsp are SURFACE fields, not profiles),
    geopotential integrated hydrostatically from the surface, then
        Pw   = q*P / (eps + (1-eps)*q)   (water-vapor partial pressure)
        n_wv = Pw / (kB * T)             (number density) [m^-3]

    Returns
    -------
    (H, T, P_level, n_wv) each sorted ascending in altitude, or None if no data:
        H        : geopotential height [m ASL]
        T        : temperature [K]
        P_level  : full model-level pressure [Pa]
        n_wv     : water-vapor number density [m^-3]
    """
    import xarray as xr

    ds = xr.open_dataset(cams_file)
    try:
        sub = ds.sel(latitude=latitude, longitude=longitude, method="nearest")
        tmask = (sub.time.values >= t_start) & (sub.time.values <= t_end)
        if not np.any(tmask):
            it = int(np.abs(sub.time.values - (t_start + (t_end - t_start) / 2)).argmin())
            tmask = np.zeros(sub.time.size, dtype=bool)
            tmask[it] = True
        sub = sub.isel(time=np.where(tmask)[0])

        level = np.asarray(ds["level"].values, dtype=int)                 # model level numbers
        T = np.asarray(sub["t"].mean("time").values, dtype=float)         # [level] K
        q = np.asarray(sub["q"].mean("time").values, dtype=float)         # [level] kg/kg
        z_raw = np.asarray(sub["z"].mean("time").values, dtype=float)     # surface geopotential (one finite level)
        lnsp = np.asarray(sub["lnsp"].mean("time").values, dtype=float)   # ln surface pressure (one finite level)
    finally:
        ds.close()

    # Sort levels ascending (1=top .. 137=surface). The hydrostatic integration
    # below starts at the surface (last index) and walks up, and the top-of-atmos
    # special case is keyed off the half-level pressure being zero (idx==0); both
    # require ascending order. CAMS files are already ascending, so this is a no-op
    # in practice but makes the routine robust to a reordered level axis.
    sort_idx = np.argsort(level)
    level = level[sort_idx]
    T = T[sort_idx]
    q = q[sort_idx]

    # z and lnsp are surface fields stored at a single level slot -> take the finite value
    z_surf = z_raw[np.isfinite(z_raw)]
    lnsp_s = lnsp[np.isfinite(lnsp)]
    if z_surf.size == 0 or lnsp_s.size == 0:
        return None
    surface_pressure = float(np.exp(lnsp_s[0]))
    z_h = float(z_surf[0])                       # surface geopotential [m^2/s^2]

    T_moist = T * (1.0 + 0.609133 * q)
    nlev = len(level)
    z_f = np.full(nlev, np.nan)
    P_level = np.full(nlev, np.nan)

    # Integrate from the surface (highest level number / last index) upward.
    for i in range(nlev - 1, -1, -1):
        idx = int(level[i]) - 1                  # half-level (level_num-1); A/B indexed by level number
        Ph_lev = _A137[idx] + _B137[idx] * surface_pressure
        Ph_levp1 = _A137[idx + 1] + _B137[idx + 1] * surface_pressure
        P_level[i] = 0.5 * (Ph_lev + Ph_levp1)
        if idx == 0:
            # Top of atmosphere: upper half-level pressure is exactly 0 (A[0]=B[0]=0),
            # so log(Ph_levp1/Ph_lev) would diverge -> ECMWF replacement (0.1 Pa, ln2).
            dlogP = np.log(Ph_levp1 / 0.1)
            alpha = np.log(2.0)
        else:
            dlogP = np.log(Ph_levp1 / Ph_lev)
            alpha = 1.0 - (Ph_lev / (Ph_levp1 - Ph_lev)) * dlogP
        TRd = T_moist[i] * _RD_CAMS
        z_f[i] = z_h + TRd * alpha
        z_h = z_h + TRd * dlogP

    H = z_f / G0                                  # geopotential height [m ASL]
    Pw = q * P_level / (EPS + (1.0 - EPS) * q)    # water-vapor partial pressure [Pa]
    n_wv = Pw / (KB * T)                          # [m^-3]

    order = np.argsort(H)
    return H[order], T[order], P_level[order], n_wv[order]


def cams_water_vapor_profile(
    cams_file: Path,
    latitude: float,
    longitude: float,
    t_start: np.datetime64,
    t_end: np.datetime64,
) -> Optional[Tuple[NDArray, NDArray]]:
    """
    Mean CAMS water-vapor number-density profile over a time window at the nearest
    grid point (thin wrapper over _cams_levels; kept for API stability).

    Returns
    -------
    (altitude_m_asl, n_wv) sorted ascending in altitude, or None if no data.
    """
    levels = _cams_levels(cams_file, latitude, longitude, t_start, t_end)
    if levels is None:
        return None
    h, _t, _p, n_wv = levels
    return h, n_wv


def cams_temperature_pressure_profile(
    cams_file: Path,
    latitude: float,
    longitude: float,
    t_start: np.datetime64,
    t_end: np.datetime64,
) -> Optional[Tuple[NDArray, NDArray, NDArray]]:
    """
    Mean CAMS temperature/pressure profile over a time window at the nearest grid
    point, for building the molecular reference (beta_mol proportional to P/T) from
    CAMS instead of the US Standard 1976 atmosphere (thin wrapper over _cams_levels).

    Returns
    -------
    (altitude_m_asl, temperature_K, pressure_Pa) sorted ascending in altitude,
    or None if no data.
    """
    levels = _cams_levels(cams_file, latitude, longitude, t_start, t_end)
    if levels is None:
        return None
    h, t, p, _n_wv = levels
    return h, t, p


def two_way_wv_transmission(
    altitude_grid_m: NDArray,
    station_altitude_m: float,
    n_wv_alt_m: NDArray,
    n_wv: NDArray,
    lut_path: Path,
    wavelength_nm: float,
    fwhm_nm: float,
) -> NDArray:
    """
    Spectrally averaged two-way water-vapor transmission on the instrument grid.

    Port of wv_t2eff.m / compute_wv_transmission.m.

    Parameters
    ----------
    altitude_grid_m   : (n_range,) instrument altitude grid [m ASL]
    station_altitude_m: station altitude [m ASL]
    n_wv_alt_m, n_wv  : water-vapor number-density profile (altitude [m ASL], value [m^-3])
    lut_path          : HITRAN cross-section LUT path
    wavelength_nm, fwhm_nm : laser central wavelength and FWHM

    Returns
    -------
    t2 : (n_range,) two-way transmission in (0, 1]
    """
    wl, lut_height_m, abscs = load_abs_cross_section(lut_path)
    range_m = altitude_grid_m - station_altitude_m            # AGL

    # Map LUT cross-sections to each range gate by nearest height (AGL grid vs LUT height)
    hidx = np.abs(range_m[None, :] - lut_height_m[:, None]).argmin(axis=0)
    abscs_col = abscs[:, hidx]                                # [n_wl, n_range] cm^2

    # Water vapor on the instrument range grid (fill below lowest level with first valid)
    wv = np.interp(altitude_grid_m, n_wv_alt_m, n_wv, left=np.nan, right=0.0)
    if np.any(np.isnan(wv)):
        first = np.where(~np.isnan(wv))[0]
        if first.size:
            wv[: first[0]] = wv[first[0]]
        wv = np.nan_to_num(wv, nan=0.0)

    return wv_t2eff_core(wavelength_nm, fwhm_nm, wl, abscs_col, range_m, wv)


def wv_t2eff_core(
    lambda0_nm: float,
    fwhm_nm: float,
    wl_nm: NDArray,
    abscs_col_cm2: NDArray,
    range_m: NDArray,
    n_wv: NDArray,
) -> NDArray:
    """
    Gaussian-weighted two-way water-vapor transmission (faithful port of wv_t2eff.m).

    Parameters
    ----------
    lambda0_nm, fwhm_nm : laser central wavelength and FWHM [nm]
    wl_nm        : (n_wl,) LUT wavelength grid [nm]
    abscs_col_cm2: (n_wl, n_range) absorption cross-section mapped to range [cm^2]
    range_m      : (n_range,) range gates [m]
    n_wv         : (n_range,) water-vapor number density [m^-3]

    Returns
    -------
    t2 : (n_range,) two-way transmission in (0, 1]
    """
    wl = np.asarray(wl_nm, dtype=float)
    sigma = fwhm_nm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gauss = np.exp(-0.5 * ((wl - lambda0_nm) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    band = (wl > lambda0_nm - 3 * sigma) & (wl < lambda0_nm + 3 * sigma)

    ext = np.asarray(abscs_col_cm2, float) * (np.asarray(n_wv, float)[None, :] / 1e4)  # cm^2->m^2 * m^-3 = m^-1
    sum_ext = np.zeros_like(ext)
    sum_ext[band, :] = _cumtrapz(ext[band, :], np.asarray(range_m, float), axis=1)
    trans = np.exp(-sum_ext)
    t2 = (trans ** 2 * gauss[:, None]).sum(axis=0) / gauss.sum()
    t2[(t2 <= 0) | ~np.isfinite(t2)] = 1.0
    return t2


def _cumtrapz(y: NDArray, x: NDArray, axis: int = -1) -> NDArray:
    """Cumulative trapezoidal integral with a leading zero (like MATLAB cumtrapz)."""
    dx = np.diff(x)
    y = np.moveaxis(y, axis, -1)
    incr = 0.5 * (y[..., 1:] + y[..., :-1]) * dx
    out = np.concatenate([np.zeros(y.shape[:-1] + (1,)), np.cumsum(incr, axis=-1)], axis=-1)
    return np.moveaxis(out, -1, axis)
