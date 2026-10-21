"""Polar-summer (low-sun) fallback for the night selection of the Rayleigh calibration (v2.3).

Why
---
The Rayleigh calibration is night-only: with ``use_sza_night=1`` it keeps the profiles whose solar
zenith angle (SZA) exceeds ``sza_night_threshold`` (100 deg) inside the night window centred on solar
midnight (``calibration.io.data_loader.filter_time_range``). The largest SZA a night reaches is the
one at solar midnight, ``180 - |latitude + declination|``. North of ~56.6 N (and south of ~56.6 S in
austral summer) it stays below 100 deg around the solstice, so no profile qualifies and the station
gets no Rayleigh calibration for weeks to months (Hopen 76.5 N: early April to mid-September).

What this module does (only when ``options.polar_night_fallback`` is True; default False)
------------------------------------------------------------------------------------------
1. Trigger. The fallback replaces the operational selection only when the night window offers no
   more than ``polar_trigger_dark_hours`` hours with SZA > ``sza_night_threshold``. The dark
   duration is computed from the solar GEOMETRY (same declination and hour-angle formula as
   ``filter_time_range``), not from the profiles present, so a data gap during a real night never
   pulls daylight profiles in.
2. Selection = the lowest sun. The profiles closest to solar midnight are kept: SZA above a floor
   ``max(min(SZA_max - polar_sza_delta, SZA(polar_min_hours)), polar_min_sza)``, where
   ``SZA(h)`` is the SZA ``h/2`` hours either side of solar midnight. The window is therefore always
   contiguous, centred on solar midnight, at least ``polar_min_hours`` long, and never admits a sun
   higher than ``90 - polar_min_sza`` deg. It is applied through the operational
   ``filter_time_range`` itself (with the floor as threshold), so the SZA formula and the profile
   bookkeeping are exactly the operational ones.
3. Geometric rejection. If no ``polar_min_hours``-long window with every profile at SZA >=
   ``polar_min_sza`` exists (sun too high), the night is rejected with ``LOW_SUN_REJECT_FLAG``.
4. Solar-background QC. By day the CHM15k under-subtracts the solar background: an additive
   residual that is constant in raw signal, i.e. constant in ``signal = rcs / r**2``, and biases the
   pointwise C_L by ``delta / (slope * p_mol(z))`` (report 19, q15e: +6 %-pt at 4 km at noon,
   growing as r**2/beta_mol). The residual is measured where the atmosphere is molecular and the
   relative effect largest, the far-range band ``polar_residual_band_m`` (8-12 km): the median of
   ``signal - slope * p_mol`` there, converted into the implied bias at the window centre. On dark
   nights it is ~0 by construction; a low-sun night whose implied bias exceeds
   ``polar_max_residual_pct`` is rejected (``LOW_SUN_REJECT_FLAG``).
5. Flag. A low-sun calibration that passes every gate is flagged ``LOW_SUN_FLAG`` (0.25) instead of
   1 / 0.5, with an explicit message, so the series and the dashboard can tell it apart.

Nothing here runs when the option is off: ``calibrate_rayleigh`` only imports this module behind the
option, so eprof_v2.2 results stay bit-identical.
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Tuple

import numpy as np

from ..config import CalibrationResult
from ..flags import LOW_SUN_FLAG, LOW_SUN_REJECT_FLAG
from ..io.data_loader import filter_time_range

logger = logging.getLogger(__name__)

#: Far-range band (m above the instrument) of the residual-baseline test. 8-12 km: above the aerosol
#: and the calibration window, below the 12-15 km bins where the night profile collapses (standard-
#: atmosphere / tropopause mismatch, q15e) and where the CHM15k signal is pure noise.
DEFAULT_RESIDUAL_BAND_M = (8000.0, 12000.0)


# ---------------------------------------------------------------------------------------------
# Solar geometry (identical formulas to calibration.io.data_loader.filter_time_range)
# ---------------------------------------------------------------------------------------------
def declination_rad(date_str: str) -> float:
    """Solar declination of the night key date, exactly as ``filter_time_range`` computes it."""
    doy = datetime.strptime(date_str, "%Y%m%d").timetuple().tm_yday
    return float(np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25))


def sza_at_hour_angle(lat_deg: float, decl_rad: float, hour_angle_deg) -> np.ndarray:
    """Solar zenith angle (deg) for an hour angle from solar NOON (deg; 180 = solar midnight)."""
    lat = np.radians(lat_deg)
    cos_sza = (np.sin(lat) * np.sin(decl_rad)
               + np.cos(lat) * np.cos(decl_rad) * np.cos(np.radians(hour_angle_deg)))
    return np.degrees(np.arccos(np.clip(cos_sza, -1.0, 1.0)))


def hours_above_sza(lat_deg: float, decl_rad: float, sza_deg: float) -> float:
    """Hours of the 24-h night window during which SZA > ``sza_deg`` (a window centred on midnight).

    SZA > s  <=>  cos(H) < q  with  q = (cos s - sin(lat) sin(decl)) / (cos(lat) cos(decl)),
    i.e. |H - 180| < 180 - arccos(q) for the hour angle H from solar noon.
    """
    lat = np.radians(lat_deg)
    den = np.cos(lat) * np.cos(decl_rad)
    num = np.cos(np.radians(sza_deg)) - np.sin(lat) * np.sin(decl_rad)
    if abs(den) < 1e-12:                       # at the pole the SZA does not change through the day
        return 24.0 if num > 0 else 0.0
    q = num / den
    if q <= -1.0:
        return 0.0
    if q >= 1.0:
        return 24.0
    return float(2.0 * (180.0 - np.degrees(np.arccos(q))) / 15.0)


@dataclass
class NightGeometry:
    """Solar geometry of one night window (night key date ``date_str``)."""
    sza_max: float        # deg, SZA at solar midnight = the lowest sun of the night
    dark_hours: float     # h with SZA > sza_night_threshold
    decl_rad: float


def night_geometry(date_str: str, lat_deg: float, threshold_deg: float) -> NightGeometry:
    """Lowest sun and dark duration of the night ending on ``date_str`` at latitude ``lat_deg``."""
    decl = declination_rad(date_str)
    sza_max = float(sza_at_hour_angle(lat_deg, decl, 180.0))
    return NightGeometry(sza_max=sza_max, dark_hours=hours_above_sza(lat_deg, decl, threshold_deg),
                         decl_rad=decl)


def low_sun_floor(geo: NightGeometry, lat_deg: float, sza_delta: float, min_hours: float,
                  min_sza: float) -> float:
    """SZA floor of the low-sun window: the larger of the 'within delta of the lowest sun' and the
    'min_hours centred on solar midnight' windows, never admitting SZA < min_sza."""
    floor_delta = geo.sza_max - max(float(sza_delta), 0.0)
    half_h = max(float(min_hours), 0.0) / 2.0
    floor_hours = float(sza_at_hour_angle(lat_deg, geo.decl_rad, 180.0 - 15.0 * half_h))
    return max(min(floor_delta, floor_hours), float(min_sza))


# ---------------------------------------------------------------------------------------------
# Fallback selection
# ---------------------------------------------------------------------------------------------
@dataclass
class LowSunNight:
    """What the fallback did on one night (attached to the result message and to fit_inputs_out)."""
    sza_max: float
    dark_hours: float
    sza_floor: float = np.nan
    window_hours: float = np.nan        # geometric duration of the window SZA > sza_floor
    n_profiles: int = 0                 # native profiles selected (before the cloud screen)
    residual_pct: float = np.nan        # implied window bias of the far-range residual (QC)
    x_far: float = np.nan               # molecular-normalised signal in the far-range band
    rejected: bool = False
    reason: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def rejection_result(self) -> CalibrationResult:
        return CalibrationResult(lidar_constant=-1, flag=LOW_SUN_REJECT_FLAG, uncertainty=0,
                                 message=f"Low-sun night rejected: {self.reason}")

    def describe(self) -> str:
        txt = (f"Low-sun (polar summer): SZA > {self.sza_floor:.1f} deg over "
               f"{self.window_hours:.1f} h (night max {self.sza_max:.1f}, dark "
               f"{self.dark_hours:.1f} h)")
        if np.isfinite(self.residual_pct):
            txt += f", far-range residual {self.residual_pct:+.1f} %"
        return txt

    def mark_result(self, result: CalibrationResult) -> CalibrationResult:
        """Turn a successful calibration into a flagged low-sun one (in place)."""
        partial = " partial," if result.flag == 0.5 else ""
        result.flag = LOW_SUN_FLAG
        result.message = f"{self.describe()};{partial} {result.message}"
        return result


def apply_low_sun_fallback(data_all, data_night, date_str: str, options
                           ) -> Tuple[object, Optional[LowSunNight]]:
    """Replace the operational night selection by the low-sun window when there is no real night.

    ``data_all`` is the loaded data BEFORE the time filter, ``data_night`` what the operational
    ``filter_time_range`` returned. Returns ``(data, low_sun)``: ``low_sun`` is None when the fallback
    does not apply (option off, clock window, coordinates missing, or a real night exists) and the
    operational ``data_night`` is returned untouched. When it applies, ``data`` holds the low-sun
    profiles; when the sun is too high ``low_sun.rejected`` is set and the caller returns
    ``low_sun.rejection_result`` (``data_night`` is passed through unchanged in that case).
    """
    if not getattr(options, "polar_night_fallback", False):
        return data_night, None
    lat, lon = float(data_all.latitude), float(data_all.longitude)
    if not (getattr(options, "use_sza_night", True) and np.isfinite(lat) and np.isfinite(lon)):
        return data_night, None                # clock-window selection: no SZA, no fallback
    thr = float(getattr(options, "sza_night_threshold", 100.0))
    geo = night_geometry(date_str, lat, thr)
    trigger = float(getattr(options, "polar_trigger_dark_hours", 0.0))
    if geo.dark_hours > trigger:
        return data_night, None                # a real night exists: operational selection

    min_sza = float(getattr(options, "polar_min_sza", 0.0))
    min_hours = float(getattr(options, "polar_min_hours", 0.0))
    floor = low_sun_floor(geo, lat, getattr(options, "polar_sza_delta", 0.0), min_hours, min_sza)
    window_h = hours_above_sza(lat, geo.decl_rad, floor)
    ls = LowSunNight(sza_max=geo.sza_max, dark_hours=geo.dark_hours, sza_floor=floor,
                     window_hours=window_h)
    if geo.sza_max <= min_sza or window_h + 1e-9 < min_hours:
        ls.rejected = True
        ls.reason = (f"sun too high (night max SZA {geo.sza_max:.1f} deg; no {min_hours:.1f}-h window "
                     f"with SZA >= {min_sza:.1f} deg)")
        logger.warning(ls.reason)
        return data_night, ls

    selected = filter_time_range(data_all, date_str, _with_threshold(options, floor))
    ls.n_profiles = len(selected.time)
    logger.info("Low-sun fallback: dark %.2f h <= %.2f h -> SZA > %.2f deg (night max %.2f), "
                "%.1f h window, %d profiles", geo.dark_hours, trigger, floor, geo.sza_max,
                window_h, ls.n_profiles)
    return selected, ls


def _with_threshold(options, threshold: float):
    o = copy.copy(options)
    o.sza_night_threshold = float(threshold)
    return o


# ---------------------------------------------------------------------------------------------
# Solar-background QC: far-range residual baseline
# ---------------------------------------------------------------------------------------------
def far_range_residual(signal, p_mol, range_alc, slope: float, window_center_m: float,
                       band_m=DEFAULT_RESIDUAL_BAND_M) -> Tuple[float, float]:
    """Implied C_L bias (%) at the window centre of an additive residual measured at far range.

    The under-subtracted solar background is constant in ``signal = rcs / r**2``:
    ``signal = slope * p_mol + delta``. ``delta`` is estimated as the MEDIAN of
    ``signal - slope * p_mol`` over the far-range band (robust to a thin cirrus layer covering part
    of the band and to noise spikes); the pointwise C_L at the window centre is then biased by
    ``delta / (slope * p_mol(z_c))``. Returns ``(bias_pct, x_far)`` with ``x_far`` the band's
    molecular-normalised signal ``1 + delta / (slope * median p_mol)``; NaNs when not measurable.
    """
    sig = np.asarray(signal, float)
    pm = np.asarray(p_mol, float)
    rng = np.asarray(range_alc, float)
    if not (np.isfinite(slope) and slope > 0 and np.isfinite(window_center_m)):
        return np.nan, np.nan
    m = (rng >= band_m[0]) & (rng <= band_m[1]) & np.isfinite(sig) & np.isfinite(pm) & (pm > 0)
    if m.sum() < 10:
        return np.nan, np.nan
    delta = float(np.median(sig[m] - slope * pm[m]))
    ok = np.isfinite(pm) & (pm > 0)
    pm_c = float(np.interp(window_center_m, rng[ok], pm[ok]))
    if not pm_c > 0:
        return np.nan, np.nan
    return 100.0 * delta / (slope * pm_c), 1.0 + delta / (slope * float(np.median(pm[m])))


def low_sun_residual_qc(low_sun: LowSunNight, signal, p_mol, range_alc, fit_result, options):
    """Background QC of a low-sun night: ``(None, None)`` if it passes, else ``(flag, message)``."""
    band = tuple(getattr(options, "polar_residual_band_m", DEFAULT_RESIDUAL_BAND_M))
    bias, x_far = far_range_residual(signal, p_mol, range_alc, float(fit_result.slope),
                                     float(fit_result.center_range_m), band)
    low_sun.residual_pct, low_sun.x_far = bias, x_far
    thr = float(getattr(options, "polar_max_residual_pct", np.inf))
    if np.isfinite(bias) and abs(bias) > thr:
        low_sun.rejected = True
        low_sun.reason = (f"far-range residual baseline implies {bias:+.1f} % at the window "
                          f"(|limit| {thr:.1f} %; X({band[0] / 1e3:.0f}-{band[1] / 1e3:.0f} km) "
                          f"= {x_far:.2f})")
        logger.warning(low_sun.reason)
        return LOW_SUN_REJECT_FLAG, f"Low-sun night rejected: {low_sun.reason}"
    return None, None
