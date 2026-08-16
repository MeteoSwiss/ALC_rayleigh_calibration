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
  -> C = S_consistent / 18.8     (S_theoretical for liquid water; O'Connor 2004 multiplier)
  -> stats (mean / median / mode / std), optional transmission correction.

Coefficient convention (Wiegner & Geiss 2012)
---------------------------------------------
The reported calibration constant is the Wiegner lidar constant ``C_L = RCS / beta_att``
(the SAME definition and units as the Rayleigh product), so the cloud and Rayleigh
constants are directly comparable on one axis. The O'Connor cloud coefficient ``C`` is a
MULTIPLIER on the file's PHYSICAL attenuated backscatter:  ``beta_true = C * beta_file``.
The method first brings the input onto a physical 1/(m*sr) scale so the apparent lidar ratio
is ~18 sr: the L2 product (stored "1E-6*1/(m*sr)" == Mm^-1 sr^-1) is multiplied by 1e-6, and a
RAW range-corrected signal (L1 rcs_0 in V*m^2, counts) is divided by the calibration constant C
(``beta = rcs_0 / C``; ``config.calibration_constant`` or the per-instrument default). RAW/L1
and L2 therefore feed the SAME physical beta and return the SAME C. ``C ~ 1`` when the applied
constant is the true one (CL61: S_apparent ~18-22 sr, C ~1.1); a FIXED nominal placeholder
(CL31/CL51 calibration_constant_0 = 1e8) leaves the apparent lidar ratio high (~78 sr) so C ~4
(= 78/18.8) is the genuine correction. Because the file's ``beta_att`` was produced with the
applied constant ``calibration_constant_0`` (= RCS/beta_att = the operationally-applied C_L),
the absolute Wiegner constant is

      C_L = calibration_constant_0 / C .

:class:`CloudCalResults` therefore exposes, in order of preference:
  - ``lidar_constant``           = C_L            (headline, Wiegner; NaN if no applied const)
  - ``calibration_factor``       = 1 / C          (dimensionless Wiegner-sense inverse)
  - ``calibration_coefficient``  = C              (O'Connor multiplier; internal/diagnostic)

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
from functools import lru_cache
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import numpy as np
from numpy.typing import NDArray
from netCDF4 import Dataset

# Re-use the verified-equivalent pieces from the Rayleigh WV module.
from ..water_vapor_correction.water_vapor import (
    wv_t2eff_core, load_abs_cross_section, in_water_vapor_band, _A137, _B137)
from ..io.cams import ensure_cams_file

# --- split-out submodules (re-exported so `from ...calibration import X` keeps working) ---
# The cloud water-vapour chain lives with the other WV code in ../water_vapor_correction/.
from ..water_vapor_correction.cloud_water_vapor import (  # noqa: F401,E402
    compute_wv_transmission, _murphy_koop_es_liquid, _wagner_pruss_pws_hpa, _nw_from_T_RH,
    _era5_levels_all_times, _interp1_linear_nan, _interp1_nearest_extrap_cols)
from ._filters import (  # noqa: F401,E402
    S_THEORETICAL, apply_multiple_scattering_correction, apply_instrument_filters,
    calculate_lidar_ratio, apply_cloud_filters, apply_temporal_consistency_filter,
    apply_transmission_correction, _trapz, _find_first, _find_last, _interp1_linear_extrap,
    _ETA_CL31, _ETA_CL51, _ETA_CL61, _ETA_CHM15K, _ETA_MINIMPL, _ETA_MPL)
# NOTE for eta-swap experiments: apply_multiple_scattering_correction reads the eta tables from
# calibration.cloud._filters (their real home) — re-exporting them here is read-only; to override
# the tables at runtime, patch calibration.cloud._filters._ETA_*, not this module's names.


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
    # MATLAB legacy gate (kept for set_defaults() parity only; superseded by
    # window_correction_threshold below and no longer used by the filters).
    window_threshold: float = 90.0
    # Window-transmission handling: the reported window transmission is an arbitrary
    # manufacturer-scaled diagnostic, and the true window attenuation is absorbed by the
    # calibration coefficient itself — so the signal is NOT corrected for it (the legacy
    # beta / (T/100)^2 correction would double-count). Profiles are only REJECTED when
    # T < window_correction_threshold (window so degraded the signal is untrustworthy).
    window_correction_threshold: float = 50.0
    energy_threshold: float = 90.0
    attenuation_factor: float = 20.0
    debug: int = 0
    apply_transmission_correction: bool = True
    aerosol_lidar_ratio_low: float = 20.0
    aerosol_lidar_ratio_high: float = 70.0
    apply_wv_correction: bool = True
    cams_folder: str = "A:\\CAMS\\"
    abs_cs_lookup_table: str = ""
    station_latitude: float = float("nan")
    station_longitude: float = float("nan")
    # Auto-download a missing CAMS file from the ADS (off by default). 'day' fetches just
    # the file's day (CAMS_Beta_<YYYYMMDD>.nc), 'month' the whole month. Needs cdsapi +
    # cfgrib + ADS credentials in ~/.cdsapirc.
    auto_download_cams: bool = False
    cams_download_scope: str = "day"
    # Fallback CAMS folder tried when the primary cams_folder lacks the month's file (e.g. a 0.4 deg
    # archive backed up by a 1 deg archive for months the 0.4 deg download has not covered). Both are
    # ECMWF/IFS MODEL levels (L137), so the water-vapour vertical resolution is preserved either way.
    cams_folder_fallback: str = ""

    # Humidity source for the water-vapour two-way transmission (910 nm only).
    #   'cams'  -> CAMS MODEL levels (L137) at the nearest grid point (OPERATIONAL DEFAULT; reads
    #              cams_folder, with cams_folder_fallback for missing months). Dense in the boundary
    #              layer, where most of the sub-cloud water vapour sits.
    #   'era5'  -> a PREFETCHED ERA5 profile cache (era5_cache), a NetCDF of all stations' q/t/z on
    #              ERA5 PRESSURE levels (built by scripts/prefetch_era5_edh.py from the DestinE Earth
    #              Data Hub). RESEARCH option only (ERA5-vs-CAMS WV comparison): the Hub's ERA5 is a
    #              19-level pressure subset (~4 levels below 3 km) -> coarser boundary layer than the
    #              CAMS model levels, so it is NOT used operationally. ERA5-only: a station/time with
    #              no cache profile is NOT calibrated (raises -> flag -4), like a missing CAMS WV.
    wv_source: str = "cams"
    era5_cache: str = ""

    # NB: aerosol_lidar_ratio is NOT defined in MATLAB set_defaults(); the runner sets it
    # to 50. apply_transmission_correction uses it, so it must be provided when that flag
    # is on (otherwise MATLAB errors). We default it to None to reproduce the "must be
    # provided" contract, and only require it when the correction actually runs.
    aerosol_lidar_ratio: Optional[float] = None

    # Optional 'YYYYMM' override for CAMS retrieval (else derived from data.time[0])
    date_str: Optional[str] = None

    # Calibration constant C used to turn a RAW range-corrected signal (L1 rcs_0 in V*m^2,
    # CHM counts, Mini-MPL photon rate) into physical attenuated backscatter beta = rcs_0 / C
    # (the L1->L2 conversion is attbsc_0[Mm^-1 sr^-1] = rcs_0 / C * 1e6, so physical beta =
    # attbsc_0 * 1e-6 = rcs_0 / C). This brings the apparent lidar ratio onto the ~18 sr scale
    # so RAW/L1 and the already-calibrated L2 give the same cloud coefficient. None -> use the
    # per-instrument default (INSTRUMENT_CAL_DEFAULT). Ignored for already-physical inputs
    # (L2 attbsc_0, CL61 L1 beta_att).
    calibration_constant: Optional[float] = None

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
    # Le couple (lambda0, FWHM) qui pilote la correction vapeur d'eau vient d'UNE seule table,
    # `water_vapor.LASER_SPECTRUM` (mesures Qmini de Payerne, 2026-06-02) -- il etait duplique
    # ici en dur. Cela compte : la position de raie du CL61 est en cours d'arbitrage (spec
    # constructeur 910,55 vs mesure 910,74 ; l'absorption effective differe d'un facteur ~3
    # entre les deux), et le jour ou elle changera il ne doit y avoir qu'un seul endroit a
    # editer. NB : deux autres longueurs d'onde CL61 coexistent volontairement dans le code et
    # ne doivent PAS etre alignees sur celle-ci : 910,0 nm pour le calcul moleculaire
    # (`InstrumentType.wavelength_nm` ; effet 0,33 % sur C_L, insensible) et 910,55 nm pour la
    # generation des tables eta de diffusion multiple (`validation/multiple_scattering_eta.py`,
    # insensible aussi).
    from ..water_vapor_correction.water_vapor import LASER_SPECTRUM

    inst = config.instrument.upper()
    _canon = {k.upper(): k for k in LASER_SPECTRUM}
    if inst in _canon:
        wl, fwhm = LASER_SPECTRUM[_canon[inst]]
    elif inst in ("MINI-MPL", "MINIMPL", "MINI_MPL", "MPL"):
        # Mini-MPL is a 532 nm system (M. Hervo): far outside the 910 nm water-vapor
        # band, so it must NEVER get the WV correction. Giving it the correct wavelength
        # makes the wavelength-based WV gate (below) skip it by physics.
        wl, fwhm = 532.0, 0.5
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
    # Operationally-applied Wiegner lidar constant C_L = RCS/beta_att that produced the
    # file's attenuated_backscatter (E-PROFILE L2 'calibration_constant_0'). Used to turn
    # the O'Connor multiplier C into an absolute C_L = calibration_constant_applied / C.
    # None when the file carries no such constant (e.g. raw L1 in counts).
    calibration_constant_applied: Optional[float] = None


def build_cloud_input(cd, config: "CloudCalConfig", rcs_units,
                               cal_const_applied=None) -> "CeiloData":
    """Build the cloud ``CeiloData`` from a ``CeilometerData`` grid (the single cloud reader path).

    Reconstructs physical attenuated backscatter oriented (range, time):
      * **L1** (``cal_const_applied`` None): ``beta = rcs * unit_factor``, and for RAW signals
        (L1 rcs_0) divide by the calibration constant ``C`` (``config.calibration_constant`` or
        the per-instrument default) so beta is physical 1/(m*sr).
      * **L2** (``cal_const_applied`` set): the loader's rcs already folds in
        ``calibration_constant_0`` (rcs = attbsc_0 * cc * 1e-6), so physical ``beta = rcs / cc``.
        The applied constant is reported so the cloud method can quote an absolute C_L.

    ``cd`` is duck-typed (``.rcs``/``.cbh``/``.range_alc``/``.time_datetime``/``.altitude``/
    ``.latitude``/``.longitude`` plus optional ``.window_transmission``/``.laser_energy``) to avoid
    a cloud <-> io.data_loader import cycle. Works on any grid — the coarse working grid the
    read-once passes share, or a directly-loaded L1/L2 file.
    """
    rcs = np.ma.filled(np.ma.masked_invalid(np.asarray(cd.rcs, dtype="float64")), np.nan)
    if cal_const_applied is not None and np.isfinite(cal_const_applied) and cal_const_applied != 0:
        beta = rcs / cal_const_applied                       # L2: recover physical beta
        applied = float(cal_const_applied)
    else:
        factor = _beta_conversion_factor(config.instrument, rcs_units)
        beta = rcs * factor
        applied = None
        if _is_raw_signal(rcs_units):
            raw_ccal = config.calibration_constant
            if raw_ccal is None:
                raw_ccal = INSTRUMENT_CAL_DEFAULT.get(config.instrument, 1.0)
            if raw_ccal and np.isfinite(raw_ccal) and raw_ccal != 0:
                beta = beta / raw_ccal
                applied = float(raw_ccal)        # the applied C (as read_ceilometer_data reported)
    beta = np.ascontiguousarray(beta.T)  # (range, time)

    # cbh: lowest layer, clipped.
    cbh = np.asarray(cd.cbh, dtype="float64")
    if cbh.ndim == 2:
        cbh = cbh[:, 0]
    cbh = cbh.copy()
    cbh[(cbh < 0) | (cbh > 20000)] = np.nan

    rng = np.asarray(cd.range_alc, dtype="float64")
    range_resol = float(rng[1] - rng[0]) if rng.size > 1 else float("nan")
    time_dt = np.asarray(cd.time_datetime, dtype="datetime64[ns]")
    wt = getattr(cd, "window_transmission", None)
    # Laser pulse energy drives the energy_rejected filter (Vaisala laser_energy; NaN for CHM15k,
    # which the < threshold test keeps) so the gate is applied, never silently disabled.
    le = getattr(cd, "laser_energy", None)
    return CeiloData(
        time=time_dt,
        time_num=_matlab_datenum(time_dt),
        station_altitude=float(cd.altitude),
        station_latitude=float(cd.latitude),
        station_longitude=float(cd.longitude),
        range=rng,
        range_resol=range_resol,
        beta=beta,
        cbh=cbh,
        quality_flag=None,
        window_transmission=(None if wt is None else np.asarray(wt, dtype="float64").ravel()),
        laser_energy=(None if le is None else np.asarray(le, dtype="float64").ravel()),
        trans2_wv=None,
        calibration_constant_applied=applied,
    )


def build_cloud_input_from_day(idd, config: "CloudCalConfig") -> "CeiloData":
    """Thin wrapper: build cloud ``CeiloData`` from a shared ``InstrumentDayData``'s coarse
    working grid (the read-once path used by the daily runner and the file entry).

    ``idd`` duck-typed (``.working`` / ``.rcs_units`` / ``.instrument_type``) to avoid an import
    cycle. Uses the coarse ``working`` grid (30 s/10 m) — the same data the other passes use —
    so cloud no longer averages separately.
    """
    config.instrument = idd.instrument_type.value
    return build_cloud_input(
        idd.working, config, idd.rcs_units,
        getattr(idd.working, "calibration_constant_applied", None))


def _matlab_datenum(dt64: NDArray) -> NDArray:
    """Convert numpy datetime64 to MATLAB datenum (days since 0000-00-00 proleptic).

    datenum(1970,1,1) = 719529, so datenum = 719529 + days_since_unix_epoch.
    """
    days = (dt64.astype("datetime64[ns]").astype("int64")) / 86400e9
    return 719529.0 + days






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
    # E-PROFILE L2 attenuated backscatter: stored as "1E-6 * 1/(m*sr)" == Mm^-1 sr^-1. Multiply
    # by 1e-6 to recover physical 1/(m*sr) so the apparent lidar ratio is ~18 sr (NOT ~1e-5).
    # MUST precede the bare "1/(m*sr)" raw form below, which it contains as a substring.
    if any(f in u for f in ("1e-6*1/(m*sr)", "1e-6*1/(sr*m)", "1e-6*1/(m.sr)", "1e-6*1/(sr.m)",
                            "1e-6sr^-1*m^-1", "1e-6sr^-1.m^-1", "1e-6m^-1*sr^-1", "1e-6m^-1.sr^-1")):
        return 1e-6
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


# Per-instrument default calibration constant C used to turn a RAW range-corrected signal into
# physical attenuated backscatter (beta = rcs_0 / C). From the L1->L2 code (section 2): these are
# the fallbacks used when no instrument-specific lidar_constant_0 is available.
INSTRUMENT_CAL_DEFAULT = {"CL31": 1e8, "CL51": 1e8, "CL61": 1.0, "CHM15k": 3e11, "Mini-MPL": 5e5}


def _is_raw_signal(units: Optional[str]) -> bool:
    """True for a RAW range-corrected signal (V*m^2, counts, MHz*km^2/uJ) that must be divided by
    the calibration constant C to become physical backscatter; False for already-physical
    backscatter (the 1/(m*sr) forms, including the 1e-6-scaled L2 product and CL61 L1 beta_att)."""
    u = _normalize_beta_units(units)
    return any(f in u for f in ("v*m^2", "v*m2", "m^2*counts/s", "m2*counts/s", "counts",
                                "mhz.km^2.uj^-1", "mhz*km^2*uj^-1", "mhz.km^2*uj^-1"))




# ===========================================================================
#  Water-vapor transmission (faithful MATLAB chain)
# ===========================================================================




















# ===========================================================================
#  apply_multiple_scattering_correction
# ===========================================================================
# Multiple-scattering factor eta(cloud-base height) for the O'Connor/Hopkin method, computed with
# the Photon Variance-Covariance model of Hogan (2006) -- a line-for-line port of the reference
# multiscatter 1.2.11 small_angle.c ('original' algorithm, incl. Eloranta's exact double scattering
# and within-gate multiple scattering; validated digit-level against the reference binary) -- for
# each instrument's receiver FOV / divergence / wavelength at droplet equivalent-area RADIUS
# a_G = 5.5 um (11 um diameter) and in-cloud extinction alpha = 10 /km.
#
# a_G = 5.5 um is the Cloudnet-MEASURED effective radius of the clouds the O'Connor method actually
# calibrates on (drizzle-free homogeneous warm stratocumulus -- ~2 um smaller than the general cloud
# population), and it is what the beta-vs-CHM15k closure selects as best across CL61/CL51/CL31 (report
# sections 9-10). It SUPERSEDES the legacy Hopkin ladder / fitted a_G = 8 um: at low cloud base eta
# rises to ~0.95 (was ~0.83), removing an over-correction of the most common (low-cloud) calibration
# scenes; the correction is slightly stronger aloft. Each Vaisala type now has its OWN table -- CL31's
# wider 0.83 mrad FOV collects more forward-scattered light than the 0.56 mrad CL51/CL61, so it gets a
# stronger correction (it no longer borrows the CL51 table). CHM15k / Mini-MPL / MPL use the same
# a_G = 5.5 um for consistency (they are never cloud-calibrated operationally -- their reference is the
# Rayleigh method and cloud calibration on them warns about detector saturation).
# Derivation, validation and figures: validation/multiple_scattering_eta.py and
# doc/reports/multiple_scattering_check.md.






# ===========================================================================
#  apply_instrument_filters
# ===========================================================================


# ===========================================================================
#  calculate_lidar_ratio
# ===========================================================================


# ===========================================================================
#  apply_cloud_filters
# ===========================================================================


# ===========================================================================
#  apply_temporal_consistency_filter
# ===========================================================================


# ===========================================================================
#  apply_transmission_correction
# ===========================================================================


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






# ===========================================================================
#  Results container + create_empty_results
# ===========================================================================
@dataclass
class CloudCalResults:
    # Wiegner lidar constant C_L = RCS / beta_att (same definition/units as the Rayleigh
    # product). For the cloud method C_L = calibration_constant_applied / C, where C is the
    # O'Connor multiplier below. NaN when the file carries no applied constant (e.g. raw L1).
    lidar_constant: float = float("nan")          # C_L  (headline, Wiegner)
    # 1/C: dimensionless Wiegner-sense correction factor (multiply the file's C_L by this).
    calibration_factor: float = float("nan")      # 1 / C  (Wiegner-sense inverse)
    # C (O'Connor 2004 multiplier): beta_true = C * beta_file; ~1 when the file is well calibrated.
    calibration_coefficient: float = float("nan")  # C  (O'Connor multiplier, diagnostic)
    # Operationally-applied constant from the file (calibration_constant_0); used to map the
    # per-profile O'Connor C to per-profile Wiegner C_L = calibration_constant_applied / C.
    # NaN when the file carries no applied constant (e.g. raw L1).
    calibration_constant_applied: float = float("nan")
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
#  Main: liquid_cloud_calibration
# ===========================================================================
def _detect_cloud_input(nc_file: str, config: "CloudCalConfig"):
    """Resolve (InstrumentType, DataLevel, rcs_units) for a cloud input file.

    Instrument comes from ``config.instrument`` if set, else the file's
    ``instrument_type``/``title`` attribute. Level is inferred from the variables present:
    ``rcs_0`` -> L1, ``attenuated_backscatter_0`` -> L2. Cloudnet-raw (``beta_raw``/
    ``beta_att``) is not handled operationally (use ``tests/_ceilo_reader.py``). Also updates
    ``config.instrument`` in place.
    """
    from ..config import InstrumentType, DataLevel
    with Dataset(nc_file, "r") as nc:
        var_names = set(nc.variables.keys())
        inst = config.instrument
        if not inst:
            atts = " ".join(str(getattr(nc, a, "")) for a in ("instrument_type", "title")).lower()
            for key, name in (("cl61", "CL61"), ("cl51", "CL51"), ("cl31", "CL31"), ("chm", "CHM15k")):
                if key in atts:
                    inst = name
                    break
        if not inst:
            raise RuntimeError(f"Cannot determine instrument type for {nc_file}")
        config.instrument = inst
        if "rcs_0" in var_names:
            level = DataLevel.L1
            rcs_units = getattr(nc.variables["rcs_0"], "units", None)
        elif "attenuated_backscatter_0" in var_names:
            level = DataLevel.L2_DAILY
            rcs_units = None            # L2 uses the calibration-constant path, not units
        elif "beta_raw" in var_names or "beta_att" in var_names:
            raise NotImplementedError(
                "Cloudnet-raw input is research-only; use tests/_ceilo_reader.read_ceilometer_data")
        else:
            raise RuntimeError(f"No rcs_0 or attenuated_backscatter_0 in {nc_file}")
    return InstrumentType(inst), level, rcs_units


def liquid_cloud_calibration(config: CloudCalConfig) -> CloudCalResults:
    """Calibrate a single L1 or L2 file via the shared loader (``load_data``).

    Reads once into a ``CeilometerData``, block-averages to ``config.average_time_s`` /
    ``average_range_m``, reconstructs the cloud ``CeiloData`` (:func:`build_cloud_input`),
    and runs :func:`liquid_cloud_calibration_from_data`. L1 and L2 both go through the ONE reader;
    the strict WV requirement (CL31/CL51/CL61) still raises on failure. Cloudnet-raw is not handled
    here (research only).
    """
    from ..io.data_loader import load_data, average_ceilometer_data
    config = set_defaults(config)

    itype, level, rcs_units = _detect_cloud_input(config.nc_file, config)
    cd = load_data([Path(config.nc_file)], itype, level)
    if cd is None:
        raise RuntimeError("Failed to read NetCDF file")
    cd = average_ceilometer_data(cd, config.average_time_s, config.average_range_m)
    data = build_cloud_input(cd, config, rcs_units, cd.calibration_constant_applied)

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
    # The photon-counting instruments (CHM15k/CHM8k, Mini-MPL, MPL) SATURATE in the strong
    # liquid-cloud return, so their integrated backscatter is biased low and the derived
    # coefficient is unreliable — warn whenever the cloud method is run on them (their
    # reference method is the Rayleigh calibration).
    if str(getattr(config, "instrument", "")).upper() in (
            "CHM15K", "CHM8K", "MINI-MPL", "MINIMPL", "MPL"):
        warnings.warn(
            f"{config.instrument}: photon-counting detector saturates in liquid clouds; "
            "the cloud-calibration coefficient is unreliable (biased by detector "
            "non-linearity). Use the Rayleigh calibration for this instrument type.",
            UserWarning)

    # NOTE: no averaging here. The input arrives pre-coarsened — the daily runner shares the
    # loader's working grid (30 s/10 m) and the file entry above averages via the loader — so
    # cloud has no separate averaging mechanism (the old average_ceilo_data path is retired).

    # --- Water-vapor absorption correction ---
    # Gate on the laser WAVELENGTH, not the instrument name: the WV correction applies
    # only inside the 900-920 nm H2O absorption band (CL31/CL51/CL61). 1064 nm (CHM15k)
    # and 532 nm (Mini-MPL) are outside it and must never be WV-corrected — enforced here
    # by physics rather than a hard-coded instrument list (mirrors the Rayleigh path).
    in_wv_band = in_water_vapor_band(config.wavelength)
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
        # The failed-calibration path is exactly where the dominant-rejection flag
        # (-22..-26) is derived from these stats — they must survive the early return.
        res.cloud_stats = cloud_stats
        res.consistency_stats = consistency_stats
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
    # Absolute Wiegner lidar constant C_L = calibration_constant_applied / C, comparable to the
    # Rayleigh product. Raw-beta L1 (notably CL61) carries NO applied constant, so C_L would be NaN
    # and the dashboard would drop an otherwise-good calibration. Fall back to the instrument default
    # (INSTRUMENT_CAL_DEFAULT; 1.0 for CL61) so a headline C_L is still reported -- a normalised
    # correction relative to the L1 beta, consistent with the CL61 theoretical value of 1.0.
    cc_applied = getattr(data, "calibration_constant_applied", None)
    if cc_applied is None or not np.isfinite(cc_applied):
        cc_applied = INSTRUMENT_CAL_DEFAULT.get(str(getattr(config, "instrument", "")), float("nan"))
    if cc_applied is not None and np.isfinite(cc_applied):
        res.calibration_constant_applied = float(cc_applied)
        if np.isfinite(cal_median) and cal_median != 0.0:
            res.lidar_constant = float(cc_applied / cal_median)
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
