"""
Uncertainty budget of the nightly Rayleigh lidar constant (E-PROFILE calibration v2.3).

Why a new budget
----------------
The published ``lidar_constant_uncertainty`` (v2.2) is 2 x the robust sigma of a 125-member
ensemble: 5 aerosol lidar ratios x 5 window shifts (+-200 m) x 5 time samples. On 147 CHM15k units
its median 1-sigma is 2.7 % while the observed night-to-night scatter is 8.9 % (pull std 3.3, 51 %
of the nights outside 2 sigma; Payerne 6.1 % reported vs 22 % observed). Three things are missing:

* the window-shift dimension sees only the LOCAL gradient of C_L with height (u ~ 0.3 x g for a
  gradient g in %/km), whereas the molecular window moves 1-2 km between nights. Electronic profile
  distortion, a faint aerosol layer or a non-molecular window all make C_L depend on the window
  height, and nothing in the profile can reliably tell these apart -- a large uncertainty is the
  honest answer;
* the molecular model (US Standard 1976 instead of the real T/p) is ~2 % off at 2-6 km, seasonally;
* the photon noise is only seen through 70 % random time subsets of the same profiles.

The v2.3 budget (``uncertainty_model = "v2.3"``) keeps the reported constant unchanged (same
ensemble median) and adds four terms in quadrature to the ensemble term:

    sigma^2 = sigma_ensemble^2 + sigma_vertical^2 + sigma_molecular^2 + sigma_aerosol^2
              + sigma_noise^2
    lidar_constant_uncertainty = 2 sigma * C_L              (coverage factor k = 2)

* ``sigma_vertical`` -- vertical consistency. C_L is recomputed with the operational Klett/C_L
  machinery (``_compute_cl_for_perturbation``: Klett from the ground with the reference in the
  sub-window, aerosol LR as configured) in contiguous 500 m sub-windows tiling the window, up to
  ``extend_above_m`` above it while the photon SNR allows, and ``extend_below_m`` under it but never
  under ``min_altitude_m`` (4 km: only free troposphere is probed under a window; under a low window
  the layer below is the boundary-layer aerosol). A noise-weighted straight line through ln C_L(z) gives the
  gradient g = d ln C_L / dz; the term is ``|g| x window_displacement_m``, i.e. how much the
  constant changes if the window moves by a representative inter-night displacement (1 km). The
  residual scatter of the sub-window constants about that line that the noise does not explain (a
  layer, a kink) is added in quadrature. The term then saturates at ``vertical_cap_rel`` (15 %,
  1 sigma): the independent review (2026-10) found the observed error of the steepest nights 2.3 x
  below |g| x D -- the window search keeps the cleanest part of the column.
* ``sigma_molecular`` -- US-1976 vs the real atmosphere. Default: a climatology of the CAMS-vs-
  US-1976 ln-bias fitted on 22 618 CHM15k nights of 127 units at 37-55 N (report 19, q15d), as a
  function of latitude, day of year, window height ABOVE SEA LEVEL (the fit's convention) and
  station altitude, plus its synoptic residual. With ``molecular_source = "cams"`` the reference is
  already CAMS and only a floor remains.
* ``sigma_aerosol`` -- a fixed allowance for a vertically uniform aerosol in the window, which no
  single-profile test can see (2 % at 1064 nm, report 19 / q15d; scaled as (lambda/1064)^3).
* ``sigma_noise`` -- the photon noise of the night-mean constant in the window. Default
  (``noise_estimator = "window_mean"``): the noise of the per-profile window mean of signal/p_mol,
  measured by first differences on the binned profile stack (window_mean_noise_rel), which keeps the
  noise common to the gates of a profile. Pre-review ("per_gate"): the per-gate noise of eprof_v2.2
  (temporal first differences on the native profiles), corrected for its time-binning error
  (``time_binning_factor``) and combined as independent gates -- 1.6 x too small on the CHM15k and
  ~4 x on the CL61 against split-night constants.

The pre-review v2.3 behaviour is reproduced by ``uncertainty_params = {"noise_estimator":
"per_gate", "molecular_height_reference": "agl", "vertical_cap_rel": None}``.

The QC gate ``u > C`` (flag -6, practically never fires) is replaced in v2.3 mode by a RELATIVE
limit ``u / C > max_rel_uncertainty`` (2-sigma, default from the validation, see doc/v2.3).

Nothing here runs unless ``options.uncertainty_model == "v2.3"``; with the default "ensemble" the
calibration and its NetCDF output are bit-identical to v2.2.
"""

from __future__ import annotations

import dataclasses
import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

UNCERTAINTY_MODELS = ("ensemble", "v2.3")

# Tunables of the v2.3 budget; overridable key-by-key through options.uncertainty_params.
# The values are the ones validated in doc/v2.3/uncertainty.md (30 CHM15k, 3 CL61, 3 Mini-MPL units,
# 2025-01..2026-07; time-series and same-night twin pulls), as revised by the independent review
# (2026-10; 37 fresh units, section 8): noise_estimator, molecular_height_reference, vertical_cap_rel.
# The pre-review budget: {"noise_estimator": "per_gate", "molecular_height_reference": "agl",
# "vertical_cap_rel": None}.
DEFAULT_PARAMS = dict(
    # --- vertical consistency -------------------------------------------------------------
    subwindow_m=500.0,             # width of one sub-window (contiguous, aligned on the window bottom)
    # How far under the window bottom the sub-windows start, and the floor they never go below.
    # Under a LOW window the layer below is the boundary-layer aerosol the window search avoided:
    # probing it measures that aerosol, not the window (1 km below from 1.5 km: median |g| 13.7 %/km
    # instead of 3.4, twin pulls 0.25-0.36 = over-covered x3). Under a HIGH window it is free
    # troposphere, where a strong gradient is a real warning (Payerne 2026-09-07: C_L x2 1 km below
    # the window). Hence 1 km below, but never under 4 km AGL (doc/v2.3/uncertainty.md section 3.4).
    extend_below_m=1000.0,
    min_altitude_m=4000.0,
    extend_above_m=1000.0,         # how far above the window top they may go ...
    snr_min=3.0,                   # ... while the mean per-gate SNR of the sub-window >= this
    max_altitude_m=10000.0,        # hard upper limit (m AGL)
    min_gates=5,                   # a sub-window needs this many finite gates
    noise_floor_rel=0.01,          # 1-sigma floor added to every sub-window's noise (weights)
    window_displacement_m=1000.0,  # representative inter-night window displacement D (1 sigma)
    include_residual=True,         # add the noise-unexplained scatter about the line
    # Saturation of the vertical term (1 sigma, relative): s -> s / sqrt(1 + (s / cap)^2); None = off
    # (pre-review). The window search keeps the cleanest part of the column, so the window constant's
    # error does not grow without bound with the gradient measured ACROSS the probed range: the
    # steepest decile scattered 2.3 x less than |g| x D on the time-series truth (an upper bound).
    # 0.15 was chosen on the 36 design units and validated on 37 fresh units (top-decile pull 0.46 ->
    # 1.02, decile spread -64 %; review 2026-10, doc/v2.3/uncertainty.md section 8).
    vertical_cap_rel=0.15,
    # --- molecular model ----------------------------------------------------------------------
    molecular="climatology",       # "climatology" | "none"
    molecular_floor_rel=0.005,     # CAMS T/p error floor when molecular_source == "cams"
    # --- residual, vertically uniform aerosol in the window ---------------------------------------
    # Invisible to any single-profile test (R(z) stays flat). 1 sigma at 1064 nm from report 19 /
    # q15d (CAMS in-window aerosol, scaled by its observed transfer to C_L: 2.2 % night-to-night),
    # scaled to other wavelengths as (lambda / 1064 nm)^3 (beta_mol ~ lambda^-4, beta_aer ~ lambda^-1).
    aerosol_floor_rel=0.02,
    # Height of the window centre the climatology is evaluated at: "asl" is the convention the
    # coefficients were fitted with (validation/v23_uncertainty/u23r_molclim.py reproduces them);
    # "agl" is the pre-review behaviour (kept to reproduce it).
    molecular_height_reference="asl",
    # --- photon noise ---------------------------------------------------------------------------
    noise_time_aggregation_penalty=1.0,   # extra factor on the propagated noise (1 = white noise)
    # "per_gate": per-gate noise (sigma_signal / sqrt(k)) combined as independent gates (pre-review;
    # 1.6 x too small on the CHM15k, 3.4 x on the CL61 against split-night constants);
    # "window_mean": the noise of the per-profile window mean measured on the binned stack
    # (window_mean_noise_rel), which keeps range-correlated noise; per-gate where not measurable.
    noise_estimator="window_mean",
    # --- QC -----------------------------------------------------------------------------------
    # Flag -6 when u / C_L (2 sigma) exceeds this. 0.6 = a 1-sigma uncertainty of 30 %, ~4 x the
    # median night (7.7 %): beyond it kalman_best_estimate stops down-weighting (its weight ratio is
    # clipped to x4 in sigma). With the pre-review budget it rejected 7.7 % of the nights, 3.4 x
    # enriched in nights off by > 25 %, but 40-50 % of them good nights, and with the Kalman weighting
    # on the rejection made co-located best estimates LESS consistent (review 2026-10). With the
    # saturated vertical term it fires on < 0.5 % of the nights (extreme noise / ensemble spread).
    max_rel_uncertainty=0.6,
)

# ---------------------------------------------------------------------------------------------
# Molecular-model climatology (CAMS vs US Standard 1976), ln-bias of C_L in percent.
# Least-squares fit on report 19 / q15d per-night biases (validation/chm15k_stability/
# q15d_atmosphere.py: b_mol = ln slope of the CAMS molecular signal against the US-1976 one inside
# each night's own window), 22 618 nights, 127 CHM15k units, 2025-01..2026-07, latitudes 36-60 N:
#   b[%] = c0 + c1 cos(w d) + c2 sin(w d) + c3 (lat-47) + c4 (lat-47) cos(w d)
#          + c5 (zc-4.5) + c6 (zc-4.5) cos(w d) + c7 h_stn
# with d = day of year, w = 2 pi / 365.25, zc = window centre (km AGL), h_stn = station altitude (km).
# Residual (synoptic) standard deviation 1.35 %. Overall RMS of the realised bias 1.98 %.
_MOL_COEF = (-0.421, 1.516, 0.799, 0.080, -0.017, -0.493, -0.210, 0.008)
_MOL_RESID_PCT = 1.35
_MOL_LAT_RANGE = (35.0, 62.0)      # latitude span of the fit; outside it the edge value is used


def molecular_bias_climatology(latitude: float, day_of_year: float, window_centre_km: float,
                               station_altitude_km: float = 0.0) -> float:
    """Expected ln-bias of C_L (fraction, not %) from using US-1976 instead of the real T/p.

    Southern-hemisphere stations use the mirrored season (day + half a year) and |latitude|; the
    latitude is clipped to the span the fit was made on.
    """
    lat = float(latitude)
    doy = float(day_of_year)
    if lat < 0:
        lat = -lat
        doy = doy + 365.25 / 2.0
    lat = float(np.clip(lat, *_MOL_LAT_RANGE))
    w = 2.0 * np.pi / 365.25
    c = np.cos(w * doy)
    s = np.sin(w * doy)
    zc = float(window_centre_km) if np.isfinite(window_centre_km) else 4.5
    k = _MOL_COEF
    pct = (k[0] + k[1] * c + k[2] * s + k[3] * (lat - 47.0) + k[4] * (lat - 47.0) * c
           + k[5] * (zc - 4.5) + k[6] * (zc - 4.5) * c + k[7] * float(station_altitude_km))
    return pct / 100.0


@dataclass
class UncertaintyBudget:
    """The v2.3 uncertainty budget of one night. All ``u_*`` are absolute 2-sigma values (units of
    C_L); ``sigma_*_rel`` the corresponding relative 1-sigma values."""
    model: str
    lidar_constant: float
    u_total: float
    u_ensemble: float
    u_vertical: float
    u_molecular: float
    u_noise: float
    sigma_ensemble_rel: float
    sigma_vertical_rel: float
    sigma_molecular_rel: float
    sigma_noise_rel: float
    u_aerosol: float = 0.0
    sigma_aerosol_rel: float = 0.0
    gradient_per_km: float = np.nan          # d ln C_L / dz (1/km) from the sub-windows
    gradient_sigma_per_km: float = np.nan    # its 1-sigma standard error (noise only)
    residual_rel: float = np.nan             # noise-unexplained scatter about the line (1 sigma)
    chi2_red: float = np.nan                 # reduced chi-square of the line fit
    molecular_bias_rel: float = np.nan       # expected US-1976 ln-bias (climatology), signed
    n_subwindows: int = 0
    z_lo_m: float = np.nan                   # sub-window range used (m AGL)
    z_hi_m: float = np.nan
    exceeds_limit: bool = False              # u / C_L > max_rel_uncertainty (QC flag -6)
    max_rel_uncertainty: float = np.nan
    sub_centre_m: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    sub_cl: NDArray[np.float64] = field(default_factory=lambda: np.array([]))
    sub_sigma_rel: NDArray[np.float64] = field(default_factory=lambda: np.array([]))

    @property
    def rel_total(self) -> float:
        """u_total / C_L (2 sigma)."""
        return self.u_total / self.lidar_constant if self.lidar_constant > 0 else np.nan

    def components(self) -> dict:
        """Scalar components for the calibration NetCDF / CSV (absolute 2-sigma + diagnostics)."""
        return dict(
            lidar_constant_uncertainty_ensemble=self.u_ensemble,
            lidar_constant_uncertainty_vertical=self.u_vertical,
            lidar_constant_uncertainty_molecular=self.u_molecular,
            lidar_constant_uncertainty_aerosol=self.u_aerosol,
            lidar_constant_uncertainty_noise=self.u_noise,
            lidar_constant_height_gradient=self.gradient_per_km,
        )


def resolve_params(options) -> dict:
    """DEFAULT_PARAMS overridden by ``options.uncertainty_params`` (unknown keys are ignored)."""
    p = dict(DEFAULT_PARAMS)
    over = getattr(options, "uncertainty_params", None) or {}
    for k, v in dict(over).items():
        if k in p:
            p[k] = v
        else:
            logger.warning("uncertainty_params: unknown key %r ignored", k)
    return p


def _expected_signal(signal, p_mol, ref_mask):
    """Molecular signal expected at every gate, c_hat * p_mol, with c_hat the median of
    signal / p_mol in the reference (operational) window. The noise is referred to THIS rather than
    to the measured signal: at a per-gate SNR below ~1 the measured signal crosses zero and
    sigma / |signal| explodes, whereas the constant is built from a window MEAN of the signal."""
    s = np.asarray(signal, float)
    p = np.asarray(p_mol, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = s[ref_mask] / p[ref_mask]
    ratio = ratio[np.isfinite(ratio)]
    if ratio.size < 3:
        return None
    c_hat = float(np.median(ratio))
    if not np.isfinite(c_hat) or c_hat <= 0:
        return None
    expected = c_hat * p
    expected[~np.isfinite(expected) | (expected <= 0)] = np.nan
    return expected


def gate_relative_noise(signal: NDArray[np.float64],
                        sigma_signal: Optional[NDArray[np.float64]],
                        p_mol: Optional[NDArray[np.float64]] = None,
                        ref_mask: Optional[NDArray[np.bool_]] = None
                        ) -> Optional[NDArray[np.float64]]:
    """Per-gate relative 1-sigma noise of the night-mean signal, sigma / (c_hat p_mol), or None.

    Without ``p_mol``/``ref_mask`` the measured |signal| is used as the reference (only safe at a
    high per-gate SNR)."""
    if sigma_signal is None:
        return None
    s = np.asarray(signal, float)
    sg = np.asarray(sigma_signal, float)
    if sg.shape != s.shape:
        return None
    ref = None
    if p_mol is not None and ref_mask is not None:
        ref = _expected_signal(s, p_mol, ref_mask)
    if ref is None:
        ref = np.abs(s)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = sg / ref
    e[~np.isfinite(e) | (e <= 0)] = np.nan
    return e


def time_binning_factor(sigma_native: Optional[dict], hours) -> float:
    """Native profiles per fit-grid profile, k = dt_binned / dt_native (>= 1).

    calibrate_rayleigh's ``sigma_signal`` divides the noise of ONE NATIVE profile (temporal first
    differences, 12-15 s on the CHM15k) by sqrt(number of 300-s BINNED profiles), so it is too large
    by sqrt(k) (~4.5 on the CHM15k; measured 4.1-4.9 against the scatter of consecutive binned
    profiles, validation/v23_uncertainty/u23_noise_check.py). The v2.2 gates keep using it unchanged
    (bit-identical); the v2.3 budget divides it by sqrt(k)."""
    try:
        dt_n = float((sigma_native or {}).get("dt_s", np.nan))
        h = np.sort(np.asarray(hours, float))
        if h.size < 3 or not np.isfinite(dt_n) or dt_n <= 0:
            return 1.0
        dt_b = float(np.median(np.diff(h))) * 3600.0
        k = dt_b / dt_n
        return float(k) if np.isfinite(k) and k > 1.0 else 1.0
    except Exception:                       # noqa: BLE001 - a diagnostic must never fail a night
        return 1.0


def window_noise_rel(e: Optional[NDArray[np.float64]], mask: NDArray[np.bool_]) -> float:
    """Relative 1-sigma noise of the constant from one window: C_L is the window mean of
    rcs / beta_mol times the Klett transmission term (the gate-to-gate noise cancels in the Klett
    normalisation, see doc/v2.3/uncertainty.md), so its noise is that of a mean over the gates."""
    if e is None:
        return np.nan
    ee = e[mask]
    ee = ee[np.isfinite(ee)]
    if ee.size == 0:
        return np.nan
    return float(np.sqrt(np.sum(ee ** 2)) / ee.size)


def window_mean_noise_rel(signal_stack: Optional[NDArray[np.float64]],
                          p_mol: NDArray[np.float64], mask: NDArray[np.bool_],
                          median_aggregation: bool = False, min_rows: int = 8) -> float:
    """Relative 1-sigma noise of the constant of one (sub-)window, measured on the binned stack.

    C_L is, to first order, the window mean of rcs / beta_mol of the night-mean profile, i.e. the
    time mean of the per-profile window means of signal / p_mol. Forming that window mean FIRST keeps
    every noise component that is common to the gates of a profile (a per-profile background
    subtraction error, range-correlated detector noise, fast atmospheric fluctuation), which the
    per-gate estimate (window_noise_rel: independent gates, sqrt(N) reduction) averages away. The
    robust first difference of consecutive profiles / sqrt(2) is the noise of one binned profile,
    / sqrt(n) that of the night mean (x 1.2533 for a median aggregation). Consecutive differences
    also carry the atmospheric change over one bin (5 min): an upper bound, as for
    molecular_methods._binned_mean_noise. NaN when fewer than ``min_rows`` profiles are usable.
    (Review 2026-10, doc/v2.3/uncertainty.md section 8: the per-gate estimate is 1.6 x too small on
    the CHM15k and 2-6 x on the CL61 against split-night constants.)
    """
    if signal_stack is None:
        return np.nan
    st = np.asarray(signal_stack, float)
    pm = np.asarray(p_mol, float)
    m = np.asarray(mask, bool) & np.isfinite(pm) & (pm > 0)
    if st.ndim != 2 or st.shape[1] != pm.size or st.shape[0] < min_rows or m.sum() < 3:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = np.nanmean(st[:, m] / pm[m][None, :], axis=1)
    r = r[np.isfinite(r)]
    if r.size < min_rows:
        return np.nan
    d = np.diff(r)
    s1 = 1.4826 * float(np.median(np.abs(d - np.median(d)))) / np.sqrt(2.0)
    level = float(np.mean(r))
    if not np.isfinite(level) or level <= 0:
        return np.nan
    return s1 / np.sqrt(r.size) * (1.2533 if median_aggregation else 1.0) / level


def _fallback_gate_noise(signal, p_mol, range_alc, start_m, end_m):
    """Gate noise when no measured noise profile exists (e.g. too few profiles): a constant
    ABSOLUTE noise (background-limited lidar) estimated from the robust in-window scatter of the
    signal about c_hat p_mol -- conservative, since it also contains atmospheric structure --
    referred to the expected molecular signal, so the relative noise grows with height."""
    m = (range_alc >= start_m) & (range_alc <= end_m) & np.isfinite(signal) & np.isfinite(p_mol)
    if m.sum() < 5:
        return None
    expected = _expected_signal(signal, p_mol, m)
    if expected is None:
        return None
    res = (signal - expected)[m]
    res = res[np.isfinite(res)]
    if res.size < 5:
        return None
    sigma_abs = 1.4826 * np.median(np.abs(res - np.median(res)))
    sigma_abs = max(sigma_abs, 1e-4 * float(np.nanmedian(expected[m])))
    with np.errstate(divide="ignore", invalid="ignore"):
        e = sigma_abs / expected
    e[~np.isfinite(e) | (e <= 0)] = np.nan
    return e


def subwindow_constants(cl_in_window: Callable[[float, float], Optional[float]],
                        range_alc: NDArray[np.float64],
                        gate_noise: Optional[NDArray[np.float64]],
                        valid: NDArray[np.bool_],
                        window_start_m: float, window_end_m: float,
                        params: dict,
                        sub_noise: Optional[Callable[[NDArray[np.bool_]], float]] = None):
    """C_L in contiguous sub-windows tiling the usable molecular range.

    The sub-windows are aligned on the operational window bottom (so they tile the window itself)
    and extend ``extend_below_m`` under it (not below ``min_altitude_m``) and above the window top
    while the mean per-gate SNR of the sub-window stays >= ``snr_min``. ``sub_noise(mask)``, when
    given, returns the relative noise of a sub-window's constant (the per-gate estimate is used where
    it returns NaN); the SNR stop always uses the per-gate noise.

    Returns (centres_m, cl, sigma_rel, z_lo, z_hi) as numpy arrays / floats.
    """
    width = float(params["subwindow_m"])
    zb, zt = float(window_start_m), float(window_end_m)
    z_lo = max(zb - float(params["extend_below_m"]), float(params["min_altitude_m"]))
    z_lo = min(z_lo, zb)                              # never start above the window itself
    z_max = min(float(params["max_altitude_m"]), float(range_alc[-1]),
                max(zt + float(params.get("extend_above_m", np.inf)), zt + width))
    snr_min = float(params["snr_min"])
    min_gates = int(params["min_gates"])

    k_lo = int(np.floor((z_lo - zb) / width + 1e-9))  # negative: sub-windows under the window
    k = k_lo
    centres, cls, sig = [], [], []
    while True:
        a = zb + k * width
        b = a + width
        if b > z_max + 1e-6:
            break
        if a < z_lo - 1e-6:                            # partial first sub-window: clip to z_lo
            a = z_lo
        above_top = a >= zt - 1e-6
        m = (range_alc >= a) & (range_alc <= b) & valid
        if m.sum() >= min_gates:
            if above_top and gate_noise is not None:
                # SNR test only ABOVE the operational window: the window itself is always covered
                ee = gate_noise[m]
                ee = ee[np.isfinite(ee)]
                mean_snr = float(np.mean(1.0 / ee)) if ee.size else 0.0
                if mean_snr < snr_min:
                    break
            c = cl_in_window(a, b)
            if c is not None and np.isfinite(c) and c > 0:
                centres.append(0.5 * (a + b))
                cls.append(float(c))
                s_k = sub_noise(m) if sub_noise is not None else np.nan
                if not np.isfinite(s_k):
                    s_k = window_noise_rel(gate_noise, m) if gate_noise is not None else np.nan
                sig.append(s_k)
        elif above_top:
            break                                       # ran out of data above the window
        if above_top and gate_noise is None and b >= zt + width - 1e-6:
            break                                       # no noise profile: one sub-window above only
        k += 1
    z_hi = (centres[-1] + 0.5 * width) if centres else np.nan
    return (np.asarray(centres, float), np.asarray(cls, float), np.asarray(sig, float),
            float(z_lo), float(z_hi))


def fit_log_gradient(centres_m, cl, sigma_rel, window_centre_m, noise_floor_rel=0.01):
    """Noise-weighted straight line through ln C_L(z).

    Returns dict(g, g_sigma, chi2_red, resid_rel, n): g = d ln C / dz in 1/km, g_sigma its standard
    error from the sub-window noise, chi2_red the reduced chi-square of the line, resid_rel the
    1-sigma scatter about the line in excess of the noise (0 when the noise explains it)."""
    out = dict(g=np.nan, g_sigma=np.nan, chi2_red=np.nan, resid_rel=np.nan, n=0)
    z = np.asarray(centres_m, float)
    y = np.log(np.asarray(cl, float))
    s = np.asarray(sigma_rel, float)
    ok = np.isfinite(z) & np.isfinite(y)
    z, y, s = z[ok], y[ok], s[ok]
    n = int(z.size)
    out["n"] = n
    if n < 3:
        return out
    s = np.where(np.isfinite(s) & (s > 0), s, np.nanmedian(s[np.isfinite(s) & (s > 0)])
                 if np.any(np.isfinite(s) & (s > 0)) else 0.05)
    s = np.sqrt(s ** 2 + float(noise_floor_rel) ** 2)
    w = 1.0 / s ** 2
    x = (z - float(window_centre_m)) / 1000.0           # km from the window centre
    xm = np.sum(w * x) / np.sum(w)
    ym = np.sum(w * y) / np.sum(w)
    sxx = np.sum(w * (x - xm) ** 2)
    if not np.isfinite(sxx) or sxx <= 0:
        return out
    g = float(np.sum(w * (x - xm) * (y - ym)) / sxx)
    res = y - (ym + g * (x - xm))
    chi2 = float(np.sum(w * res ** 2))
    dof = n - 2
    chi2_red = chi2 / dof if dof > 0 else np.nan
    # Excess scatter about the line: weighted mean square of the residuals minus what the noise
    # explains, per sub-window. A layer or a kink shows up here; pure noise gives 0.
    msr = float(np.sum(w * res ** 2) / np.sum(w))
    noise_ms = float(np.sum(w * s ** 2) / np.sum(w)) * dof / n
    resid = float(np.sqrt(max(msr - noise_ms, 0.0)))
    out.update(g=g, g_sigma=float(1.0 / np.sqrt(sxx)), chi2_red=chi2_red, resid_rel=resid)
    return out


def uncertainty_budget(*, cl: float, u_ensemble: float, rcs_mean: NDArray[np.float64],
                       range_alc: NDArray[np.float64], p_mol: NDArray[np.float64],
                       fit_result, sigma_signal: Optional[NDArray[np.float64]],
                       cl_in_window: Callable[[float, float], Optional[float]],
                       options, latitude: float, station_altitude: float,
                       date_str: str, wavelength_nm: float = 1064.0,
                       signal_stack: Optional[NDArray[np.float64]] = None,
                       median_aggregation: bool = False) -> UncertaintyBudget:
    """Compute the v2.3 uncertainty budget for one calibrated night (see module docstring).

    ``cl_in_window(start_m, end_m)`` must return the lidar constant the OPERATIONAL machinery gives
    when the molecular window is [start_m, end_m] (nominal lidar ratio, all profiles), or None.
    ``signal_stack`` (binned, screened profiles / r^2, the stack rcs_mean is the time mean of) is
    read only with ``noise_estimator = "window_mean"``.
    Never raises: a failing component degrades to a conservative fallback and is logged.
    """
    p = resolve_params(options)
    rng = np.asarray(range_alc, float)
    rcs = np.asarray(rcs_mean, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        signal = rcs / rng ** 2
    zb, zt = float(fit_result.range_start_m), float(fit_result.range_end_m)
    zc = 0.5 * (zb + zt)

    # (a) ensemble -------------------------------------------------------------------------------
    sig_ens = float(u_ensemble) / (2.0 * cl) if cl > 0 else np.nan

    # (d) photon noise in the operational window ---------------------------------------------
    in_window = (rng >= zb) & (rng <= zt) & np.isfinite(rcs)
    pm = np.asarray(p_mol, float)
    gate_noise = gate_relative_noise(signal, sigma_signal, pm, in_window)
    if gate_noise is None:
        gate_noise = _fallback_gate_noise(signal, pm, rng, zb, zt)
    if gate_noise is not None:
        gate_noise = gate_noise * float(p["noise_time_aggregation_penalty"])
    sig_noise = window_noise_rel(gate_noise, in_window)
    sub_noise = None
    if str(p.get("noise_estimator", "per_gate")).lower() == "window_mean" and signal_stack is not None:
        # The noise of the per-profile WINDOW MEAN, measured on the stack: keeps the noise common to
        # the gates of a profile that the independent-gate sum averages away (review 2026-10).
        penalty = float(p["noise_time_aggregation_penalty"])

        def sub_noise(mask):
            return window_mean_noise_rel(signal_stack, pm, mask, median_aggregation) * penalty
        s_wm = sub_noise(in_window)
        if np.isfinite(s_wm):
            sig_noise = s_wm
    if not np.isfinite(sig_noise):
        sig_noise = 0.0

    # (b) vertical consistency ------------------------------------------------------------------
    valid = np.isfinite(rcs)
    centres, sub_cl, sub_sig, z_lo, z_hi = (np.array([]),) * 3 + (np.nan, np.nan)
    # Nominal sub-window width first; a short usable range (low window against min_altitude_m and
    # no SNR above it) is retried with half-width sub-windows so that a line can still be fitted.
    for width in (float(p["subwindow_m"]), 0.5 * float(p["subwindow_m"])):
        try:
            centres, sub_cl, sub_sig, z_lo, z_hi = subwindow_constants(
                cl_in_window, rng, gate_noise, valid, zb, zt, dict(p, subwindow_m=width),
                sub_noise=sub_noise)
        except Exception:                   # noqa: BLE001 - a diagnostic must never lose a night
            logger.debug("sub-window constants failed", exc_info=True)
            centres, sub_cl, sub_sig, z_lo, z_hi = (np.array([]),) * 3 + (np.nan, np.nan)
        if sub_cl.size >= 3:
            break
    fit = fit_log_gradient(centres, sub_cl, sub_sig, zc, p["noise_floor_rel"])
    D_km = float(p["window_displacement_m"]) / 1000.0
    if np.isfinite(fit["g"]):
        sig_vert = abs(fit["g"]) * D_km
        if p["include_residual"] and np.isfinite(fit["resid_rel"]):
            sig_vert = float(np.hypot(sig_vert, fit["resid_rel"]))
    elif sub_cl.size == 2:
        # two sub-windows: the gradient is their difference, its noise is not separable
        g2 = (np.log(sub_cl[1]) - np.log(sub_cl[0])) / ((centres[1] - centres[0]) / 1000.0)
        fit["g"] = float(g2)
        sig_vert = abs(g2) * D_km
    else:
        # Not measurable at all. Conservative fallback: attribute the WHOLE ensemble spread to the
        # +-200 m window shifts (robust sigma ~0.3 x g for a gradient g, see the module docstring),
        # i.e. g ~ 2 sigma_ens / 0.3 per km, and scale it to the displacement D.
        sig_vert = 2.0 * sig_ens / 0.3 * D_km if np.isfinite(sig_ens) else 0.0
    cap = p.get("vertical_cap_rel")
    if cap is not None and np.isfinite(float(cap)) and float(cap) > 0:
        sig_vert = float(sig_vert / np.sqrt(1.0 + (sig_vert / float(cap)) ** 2))

    # (c) molecular model -------------------------------------------------------------------------
    mol_bias = np.nan
    mode = str(p["molecular"]).lower()
    if mode == "none":
        sig_mol = 0.0
    elif str(getattr(options, "molecular_source", "standard")) == "cams":
        sig_mol = float(p["molecular_floor_rel"])
    else:
        try:
            doy = datetime.strptime(str(date_str)[:8], "%Y%m%d").timetuple().tm_yday
        except ValueError:
            doy = 182
        # The coefficients were fitted with the window centre ABOVE SEA LEVEL (q15d windows, m ASL);
        # "agl" reproduces the pre-review evaluation (window centre above ground).
        zc_km = zc / 1000.0
        if str(p.get("molecular_height_reference", "agl")).lower() == "asl":
            zc_km += float(station_altitude) / 1000.0
        mol_bias = molecular_bias_climatology(latitude, doy, zc_km,
                                              float(station_altitude) / 1000.0)
        sig_mol = float(np.hypot(mol_bias, _MOL_RESID_PCT / 100.0))

    # (e) residual, vertically uniform aerosol -------------------------------------------------
    wl = float(wavelength_nm) if np.isfinite(wavelength_nm) and wavelength_nm > 0 else 1064.0
    sig_aer = float(p["aerosol_floor_rel"]) * (wl / 1064.0) ** 3

    sig_total = float(np.sqrt(np.nansum(np.square([sig_ens, sig_vert, sig_mol, sig_aer,
                                                    sig_noise]))))
    u_total = 2.0 * sig_total * cl
    limit = float(p["max_rel_uncertainty"])
    budget = UncertaintyBudget(
        model="v2.3", lidar_constant=float(cl), u_total=float(u_total),
        u_ensemble=float(u_ensemble), u_vertical=2.0 * sig_vert * cl,
        u_molecular=2.0 * sig_mol * cl, u_noise=2.0 * sig_noise * cl,
        u_aerosol=2.0 * sig_aer * cl, sigma_aerosol_rel=sig_aer,
        sigma_ensemble_rel=sig_ens, sigma_vertical_rel=float(sig_vert),
        sigma_molecular_rel=float(sig_mol), sigma_noise_rel=float(sig_noise),
        gradient_per_km=float(fit["g"]), gradient_sigma_per_km=float(fit["g_sigma"]),
        residual_rel=float(fit["resid_rel"]), chi2_red=float(fit["chi2_red"]),
        molecular_bias_rel=float(mol_bias), n_subwindows=int(sub_cl.size),
        z_lo_m=float(z_lo), z_hi_m=float(z_hi),
        exceeds_limit=bool(np.isfinite(u_total) and u_total > limit * cl),
        max_rel_uncertainty=limit,
        sub_centre_m=centres, sub_cl=sub_cl, sub_sigma_rel=sub_sig,
    )
    logger.info(
        "v2.3 uncertainty: total 2s %.1f%% = ens %.1f, vert %.1f (g %+.1f%%/km over %d sub-windows "
        "%.0f-%.0f m), mol %.1f, aer %.1f, noise %.1f [%% 2s]",
        100 * budget.rel_total, 200 * sig_ens, 200 * sig_vert, 100 * fit["g"],
        budget.n_subwindows, z_lo, z_hi, 200 * sig_mol, 200 * sig_aer, 200 * sig_noise)
    return budget


def replace_window(fit_result, start_m: float, end_m: float):
    """A copy of a RayleighFitResult with another molecular window (for the sub-window runs)."""
    return dataclasses.replace(fit_result, range_start_m=float(start_m), range_end_m=float(end_m))


def operational_cl_in_window(rcs_mean, range_alc, beta_mol, fit_result, options):
    """``cl_in_window(start_m, end_m)`` built on the OPERATIONAL machinery
    (calibration._compute_cl_for_perturbation: Klett from the ground with the reference in the
    window, nominal aerosol lidar ratio, no shift, all profiles)."""
    from .calibration import _compute_cl_for_perturbation   # lazy: avoids a circular import

    def cl_in_window(start_m, end_m):
        pr = _compute_cl_for_perturbation(
            rcs_mean=rcs_mean, range_alc=range_alc, beta_mol=beta_mol,
            fit_result=replace_window(fit_result, start_m, end_m),
            lidar_ratio_aerosol=options.lidar_ratio_aerosol, altitude_shift_m=0,
            subtract_background=options.subtract_background,
            consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
            sign_error_v10=options.sign_error_v10,
        )
        return None if pr is None else float(pr.lidar_constant)

    return cl_in_window


def budget_for_night(*, cl, u_ensemble, rcs_mean, range_alc, mol_props, fit_result,
                     sigma_signal, options, info, date_str, sigma_native=None,
                     hours=None, signal_stack=None) -> UncertaintyBudget:
    """Entry point used by calibrate_rayleigh (the whole v2.3 hook is this one call).

    ``sigma_native`` / ``hours`` (native noise dict, fit-grid profile times in hours) correct the
    time-binning error of ``sigma_signal`` (see time_binning_factor). ``signal_stack`` is the
    binned, screened profile stack rcs_mean averages (for noise_estimator = "window_mean")."""
    try:
        if sigma_signal is not None:
            sigma_signal = np.asarray(sigma_signal, float) / np.sqrt(
                time_binning_factor(sigma_native, hours))
        return uncertainty_budget(
            cl=cl, u_ensemble=u_ensemble, rcs_mean=rcs_mean, range_alc=range_alc,
            p_mol=mol_props.p_mol, fit_result=fit_result, sigma_signal=sigma_signal,
            cl_in_window=operational_cl_in_window(rcs_mean, range_alc, mol_props.beta_mol,
                                                  fit_result, options),
            options=options, latitude=info.latitude, station_altitude=info.altitude,
            date_str=date_str, wavelength_nm=float(info.instrument_type.wavelength_nm),
            signal_stack=signal_stack,
            median_aggregation=getattr(options, "time_aggregation", "mean") == "median")
    except Exception:                       # noqa: BLE001 - never lose a calibrated night
        logger.warning("v2.3 uncertainty budget failed; reporting the ensemble term only",
                       exc_info=True)
        sig_ens = float(u_ensemble) / (2.0 * cl) if cl > 0 else np.nan
        limit = float(resolve_params(options)["max_rel_uncertainty"])
        return UncertaintyBudget(
            model="v2.3-fallback", lidar_constant=float(cl), u_total=float(u_ensemble),
            u_ensemble=float(u_ensemble), u_vertical=np.nan, u_molecular=np.nan, u_noise=np.nan,
            u_aerosol=np.nan, sigma_aerosol_rel=np.nan,
            sigma_ensemble_rel=sig_ens, sigma_vertical_rel=np.nan, sigma_molecular_rel=np.nan,
            sigma_noise_rel=np.nan, exceeds_limit=bool(u_ensemble > limit * cl),
            max_rel_uncertainty=limit)
