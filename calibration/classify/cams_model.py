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
from pathlib import Path

import numpy as np

from ceiloclass.model import Model

from ..water_vapor_correction.water_vapor import cams_temperature_pressure_profile

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
    """Temperature ``Model`` on the (obs_time, obs_range) grid from a CAMS file.

    Args:
        cams_file: ``CAMS_Beta_<month>.nc`` (or regional variant).
        lat, lon: Station coordinates; the nearest CAMS grid point is used.
        obs_time: Ceilometer time (datetimes); its span sets the CAMS averaging window.
        obs_range: Ceilometer range (m).
        altitude: Site altitude (m a.s.l.), to align range (a.g.l.) with CAMS a.s.l.

    Returns:
        A ceiloclass ``Model`` (dry-bulb temperature + extrapolation flag).
    """
    t_start = np.datetime64(min(obs_time))
    t_end = np.datetime64(max(obs_time))
    profile = cams_temperature_pressure_profile(
        Path(cams_file), float(lat), float(lon), t_start, t_end
    )
    if profile is None:
        msg = f"No CAMS data at ({lat},{lon}) in {cams_file}"
        raise ValueError(msg)
    h_asl, temp, _pressure = profile  # sorted ascending in altitude

    obs_h_asl = np.asarray(obs_range, dtype=float) + float(altitude)
    t_on_range = _temperature_on_heights(obs_h_asl, np.asarray(h_asl), np.asarray(temp))
    # Temperature varies little over a night: use the mean profile for every step.
    tw = np.broadcast_to(t_on_range, (len(obs_time), obs_h_asl.size)).copy()
    extrapolated = np.broadcast_to(
        obs_h_asl > h_asl[-1], (len(obs_time), obs_h_asl.size)
    ).copy()
    return Model(tw, extrapolated)
