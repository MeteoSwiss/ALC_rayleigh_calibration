"""Read CAMS ground-referenced aerosol backscatter at a station.

The CAMS_Beta_*.nc files carry ``aerbackscatgnd355/532/1064`` (attenuated
backscatter due to aerosol, from the ground), in m^-1 sr^-1, on the model
levels. We reuse :func:`calibration.cloud.calibration._cams_levels_all_times`
for the model-level geopotential height and time axis (it is LRU-cached, so the
extra call is free) and read only the two aerosol-backscatter columns here.

Wavelength handling (matches E_PROFILE_ALC_Monthly_OB.m):
* 532 nm  (Mini-MPL)      -> aerbackscatgnd532 directly
* 1064 nm (CHM15k)        -> aerbackscatgnd1064 directly
* 910 nm  (CL31/CL51/CL61)-> Angstrom interpolation between 532 and 1064:
      ang  = -log(b532 / b1064) / log(532 / 1064)
      b910 = b1064 * (910 / 1064) ** (-ang)
  The water-vapour absorption at 910 nm is applied separately to the
  *observation* (see :mod:`calibration.omb.omb`), not to the CAMS aerosol
  backscatter (which is a modelled, absorption-free quantity).
"""

from __future__ import annotations

import numpy as np
from netCDF4 import Dataset
from numpy.typing import NDArray

from ..cloud.calibration import _cams_levels_all_times

__all__ = ["cams_aerosol_backscatter", "angstrom_interpolate"]

_WL_532 = 532.0
_WL_1064 = 1064.0


def angstrom_interpolate(
    b532: NDArray, b1064: NDArray, wavelength_nm: float
) -> NDArray:
    """Aerosol backscatter at ``wavelength_nm`` from the 532/1064 nm pair via the
    Angstrom exponent. Gates where either anchor is non-positive (top of
    atmosphere) return NaN."""
    b532 = np.asarray(b532, dtype=float)
    b1064 = np.asarray(b1064, dtype=float)
    out = np.full(b532.shape, np.nan)
    ok = (b532 > 0) & (b1064 > 0)
    if np.any(ok):
        ang = -np.log(b532[ok] / b1064[ok]) / np.log(_WL_532 / _WL_1064)
        out[ok] = b1064[ok] * (wavelength_nm / _WL_1064) ** (-ang)
    return out


def _read_aer_column(nc: Dataset, name: str, li: int, ai: int) -> NDArray:
    """Read aerbackscatgnd<name> at one grid point as (level, time).

    Slices the NetCDF variable directly (never materialises the 4-D array),
    mirroring ``_cams_levels_all_times._read_profile``.
    """
    var = nc.variables[name]
    dim_names = var.dimensions  # ('time','level','latitude','longitude')
    axes = {nm: k for k, nm in enumerate(dim_names)}
    idx = [slice(None)] * len(dim_names)
    idx[axes["longitude"]] = li
    idx[axes["latitude"]] = ai
    sub = np.asarray(var[tuple(idx)], dtype="float64")  # remaining: level, time (orig order)
    remaining = [nm for nm in dim_names if nm in ("level", "time")]
    lev_pos = remaining.index("level")
    return np.moveaxis(sub, lev_pos, 0)  # (level, time)


def cams_aerosol_backscatter(
    cams_file: str, latitude: float, longitude: float, wavelength_nm: float
) -> tuple[NDArray, NDArray, NDArray]:
    """CAMS aerosol backscatter at the nearest grid point, on the model levels.

    Returns
    -------
    time_num : (n_t,)        MATLAB datenum of the CAMS steps
    z_model  : (n_lev, n_t)  geopotential height [m ASL] (native CAMS order)
    beta_aer : (n_lev, n_t)  aerosol backscatter [m^-1 sr^-1] at wavelength_nm
    """
    # Heights + time come from the cached molecular reader (same grid point).
    time_num, z_model, _T, _nw = _cams_levels_all_times(cams_file, latitude, longitude)

    with Dataset(cams_file, "r") as nc:
        lon_m = np.asarray(nc.variables["longitude"][:], dtype="float64")
        lat_m = np.asarray(nc.variables["latitude"][:], dtype="float64")
        li = int(np.argmin(np.abs(lon_m - longitude)))
        ai = int(np.argmin(np.abs(lat_m - latitude)))
        wl = round(float(wavelength_nm))
        if wl == 532:
            beta_aer = _read_aer_column(nc, "aerbackscatgnd532", li, ai)
        elif wl == 1064:
            beta_aer = _read_aer_column(nc, "aerbackscatgnd1064", li, ai)
        elif wl == 355:
            beta_aer = _read_aer_column(nc, "aerbackscatgnd355", li, ai)
        else:
            b532 = _read_aer_column(nc, "aerbackscatgnd532", li, ai)
            b1064 = _read_aer_column(nc, "aerbackscatgnd1064", li, ai)
            beta_aer = angstrom_interpolate(b532, b1064, float(wavelength_nm))

    return time_num, z_model, beta_aer
