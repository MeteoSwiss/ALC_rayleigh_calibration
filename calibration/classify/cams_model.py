"""Build a ceiloclass ``Model`` (temperature on the ceilometer grid) from an
E-PROFILE ``CAMS_Beta_<month>.nc`` file — the *operational* temperature source,
replacing a Cloudnet ECMWF model file.

In these CAMS files ``z``/``lnsp`` are *surface* fields (geopotential per level is
mostly fill), so per-level pressure and geometric height are reconstructed from the
ECMWF L137 a/b coefficients + surface pressure. We reuse the calibration pipeline's
own shared CAMS reader (:func:`calibration.water_vapor_correction.water_vapor.cams_temperature_pressure_profile`)
— the same function the WV correction uses — rather than re-deriving the coefficients.

We take its time-mean (altitude, T, P) profile over the night and interpolate T onto the
ceilometer range grid (aligned to a.s.l. via the site altitude). Dry-bulb by default, which
is all aerosol/cloud/ice *detection* needs; pressure is returned too, so wet-bulb (a sharper
ice/liquid split) is a small extension.
"""

from __future__ import annotations

from os import PathLike

import numpy as np

from ceiloclass.model import Model

from ..water_vapor_correction.water_vapor import _matlab_datenum_days, cams_levels_all_times

_ISA_LAPSE = 6.5e-3  # K/m, standard tropospheric lapse rate (T rises downward)


def _temperature_on_heights(
    obs_h_asl: np.ndarray, h: np.ndarray, t: np.ndarray
) -> np.ndarray:
    """Temperature at obs heights: interpolate within CAMS, lapse-extrapolate below.

    Above the lowest CAMS level the profile is interpolated directly. *Below* it — which
    happens whenever the station sits under its (coarse) CAMS grid-cell orography —
    temperature is extrapolated with the standard-atmosphere lapse rate anchored at the
    lowest CAMS value, rather than a naive two-point linear extrapolation whose slope is set
    by the ~20 m-spaced, often surface-inverted lowest model levels and can run to absurd
    values (e.g. -110 degC at the surface for a station ~900 m below its grid cell).
    """
    out = np.interp(obs_h_asl, h, t)  # clamps outside [h[0], h[-1]]
    below = obs_h_asl < h[0]
    out[below] = t[0] + _ISA_LAPSE * (h[0] - obs_h_asl[below])
    return out


def cams_to_model(
    cams_file: str | PathLike,
    lat: float,
    lon: float,
    obs_time: np.ndarray,
    obs_range: np.ndarray,
    altitude: float,
) -> Model:
    """Time-varying temperature ``Model`` on the (obs_time, obs_range) grid from a CAMS file.

    Each CAMS time step is integrated to geometric height independently (via
    :func:`cams_levels_all_times`), mapped onto the observation height grid, then interpolated in
    time to every observation profile. This preserves CAMS's native (3-hourly) temporal resolution
    so the melting layer *moves* through the day -- a day-averaged profile both freezes the 0 degC
    isotherm and, on days with a surface-pressure swing, distorts it (model-level averaging invents a
    warm layer; see :func:`..water_vapor_correction.water_vapor._cams_levels`).

    Args:
        cams_file: ``CAMS_Beta_<month>.nc`` (or regional variant).
        lat, lon: Station coordinates; the nearest CAMS grid point is used.
        obs_time: Ceilometer time (datetimes).
        obs_range: Ceilometer range (m).
        altitude: Site altitude (m a.s.l.), to align range (a.g.l.) with CAMS a.s.l.

    Returns:
        A ceiloclass ``Model`` (dry-bulb temperature + extrapolation flag).
    """
    time_num, z_t, temp_t, _n_wv = cams_levels_all_times(str(cams_file), float(lat), float(lon))
    if time_num.size == 0:
        msg = f"No CAMS data at ({lat},{lon}) in {cams_file}"
        raise ValueError(msg)

    obs_h_asl = np.asarray(obs_range, dtype=float) + float(altitude)
    # 1) each CAMS time -> temperature on the obs height grid (a.s.l.), lapse-extrapolated below
    t_ct = np.empty((time_num.size, obs_h_asl.size))
    top_asl = np.empty(time_num.size)
    for j in range(time_num.size):
        zj = np.asarray(z_t[:, j], dtype=float)
        order = np.argsort(zj)
        t_ct[j] = _temperature_on_heights(obs_h_asl, zj[order], np.asarray(temp_t[order, j]))
        top_asl[j] = zj[order][-1]
    # 2) interpolate across CAMS time to each observation profile (moving melting layer)
    obs_num = _matlab_datenum_days(np.asarray(obs_time, dtype="datetime64[ns]"))
    tw = np.empty((obs_num.size, obs_h_asl.size))
    for i in range(obs_h_asl.size):
        tw[:, i] = np.interp(obs_num, time_num, t_ct[:, i])
    top_at_obs = np.interp(obs_num, time_num, top_asl)
    extrapolated = obs_h_asl[None, :] > top_at_obs[:, None]
    return Model(tw, extrapolated)
