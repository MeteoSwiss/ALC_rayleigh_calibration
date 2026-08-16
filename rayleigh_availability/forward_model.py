# -*- coding: utf-8 -*-
"""Forward model of a ceilometer night + closure test of the SHIPPED Rayleigh retrieval.

WHY THIS MODULE EXISTS
----------------------
The observational study (``rayleigh_availability/``) shows that the retrieved lidar constant
is NOT independent of the altitude of the molecular fit window: v2.2-recovered nights, which
fit 1-2 km higher than v2.0-kept nights, return a constant that differs by -22.5 % (Payerne
CHM15k) to +32 % (Aosta), and the sign varies by site. From observations alone we cannot say
whether that is

  (a) an ATMOSPHERIC effect (the recovered nights really are aerosol-loaded, so a fit at 4.5 km
      sees a different -- and wrong -- reference than a fit at 3 km), or
  (b) an ESTIMATOR effect (something in the shipped chain -- window search, Klett, the optical
      depth integration, the median over the window -- makes C_L depend on where you look).

The only way to separate them is a forward model: build a synthetic profile from a KNOWN lidar
constant ``C_true``, push it through the ACTUAL production code, and measure the error exactly.
That is the closure the observations cannot give.

This module is the foundation. It deliberately contains NO mechanism study: it builds the
atmosphere, builds the signal, runs the shipped retrieval, and proves that with no perturbation
the retrieval returns ``C_true``. Mechanism modules import it.

THE FORWARD EQUATION
--------------------
For a single-scattering, elastic-backscatter lidar the received power at range ``z`` is

    P(z) = C_true * O(z) * beta_tot(z) * exp(-2 * INT_0^z ext_tot(z') dz') / z^2      [1]

and the quantity the E-PROFILE pipeline calls ``rcs`` (range-corrected signal, ``rcs_0`` in L1) is

    rcs(z) = P(z) * z^2 = C_true * O(z) * beta_tot(z) * T2(z)                          [2]

with

    beta_tot = beta_mol + beta_aer                       [m^-1 sr^-1]
    ext_tot  = beta_mol * S_mol + ext_aer                [m^-1]
    S_mol    = 8*pi/3 sr   (Rayleigh extinction-to-backscatter ratio, exact)
    ext_aer  = beta_aer * S_aer,  S_aer the aerosol lidar ratio [sr]
    T2(z)    = exp(-2 * INT_0^z ext_tot dz')             [-]  two-way transmission
    O(z)     = overlap function [-], 1 in the 2-6 km fit band for every ALC considered here

``C_true`` carries whatever units make [2] dimensionally consistent with the L1 ``rcs_0`` of the
instrument; for a Payerne CHM15k it is ~6.5e11 (see ``PAYERNE_CHM15K``). Everything else is SI:
range in m, backscatter in m^-1 sr^-1, extinction in m^-1, lidar ratios in sr.

WHAT THE RETRIEVAL DOES (and which shipped functions we call)
-------------------------------------------------------------
``retrieve()`` reproduces the production chain exactly by CALLING it, never by reimplementing it:

  * ``calibration.rayleigh.atmosphere.load_standard_atmosphere``      (T/p, US Std 1976)
  * ``calibration.rayleigh.atmosphere.calculate_molecular_properties`` (beta_mol, alpha_mol, p_mol)
  * ``calibration.rayleigh.rayleigh_fit.find_optimal_molecular_window`` (window search + gates)
  * ``calibration.rayleigh.calibration._compute_cl_for_perturbation``  (Klett + C_L, one config)
        which itself calls ``atmosphere.klett_inversion`` and
        ``rayleigh_fit.calculate_lidar_constant``.

``_compute_cl_for_perturbation`` is private, and importing it is deliberate: it IS the production
kernel that turns one (window, lidar ratio, altitude shift) into one constant. Re-writing it here
would test our copy instead of the shipped estimator, which is the whole point of the exercise.
The only step of ``calibrate_rayleigh`` we do NOT run is the front end (file reading, night/cloud
screening, L1->L2 binning) -- we hand the retrieval the binned night-mean profile it would have
produced, on the real instrument grid (see ``default_grid``).

The production best estimate is the MEDIAN over 5 lidar ratios x 5 window shifts x 5 time subsets
(``calibration.LR_DELTAS``, ``ALT_SHIFTS_M``). ``retrieve(..., perturbations=True)`` reproduces the
LR x shift part (there is only one profile, so no time subsets); ``perturbations=False`` (default)
returns the nominal configuration, which is what a per-window ladder needs.

CLOSURE
-------
``run_closure_tests()`` runs the three tests the study requires and writes
``doc/reports/figs_altitude_audit/fwd_closure.png``. Run it with::

    python rayleigh_availability/forward_model.py

Reference: Wiegner & Geiss, AMT 5, 1953 (2012) for the calibration; Fernald, Appl. Opt. 23, 652
(1984) / Klett, Appl. Opt. 24, 1638 (1985) with the sign correction of Speidel & Vogelmann,
Appl. Opt. 62, 861 (2023) for the inversion.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.stats import linregress

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from calibration.config import CalibrationOptions                      # noqa: E402
from calibration.rayleigh.atmosphere import (                          # noqa: E402
    DEFAULT_STANDARD_ATMOSPHERE,
    MOLECULAR_LIDAR_RATIO,
    calculate_molecular_properties,
    load_cams_atmosphere,
    load_standard_atmosphere,
)
from calibration.rayleigh.calibration import (                         # noqa: E402
    ALT_SHIFTS_M,
    LR_DELTAS,
    _compute_cl_for_perturbation,
)
from calibration.rayleigh.rayleigh_fit import (                        # noqa: E402
    RayleighFitResult,
    find_optimal_molecular_window,
)

__all__ = [
    "PAYERNE_CHM15K",
    "default_grid",
    "make_atmosphere",
    "aerosol_layer",
    "exponential_aerosol",
    "forward_signal",
    "retrieve",
    "scan_windows",
    "run_closure_tests",
]


# ---------------------------------------------------------------------------
# Reference instrument: Payerne CHM15k (0-20000-0-06610, identifier A)
# ---------------------------------------------------------------------------
# Every number below was MEASURED on real L1 nights through the production front end
# (calibrate_rayleigh(..., fit_inputs_out=...)), not assumed:
#
#   grid            512 gates, dz = 29.97 m, first gate 22.4775 m  -> the L2 grid the pipeline
#                   bins L1 onto (300 s x 30 m); top gate 15337 m.
#   C_TRUE          6.5e11  (recorded constants: 6.45e11 on 2025-01-15, 6.48e11 on 2025-03-06,
#                   6.91e11 on 2026-05-22)
#   SIGMA_NIGHT_*   1-sigma photon noise of the NIGHT-MEAN range-normalised signal (rcs/z^2) on
#                   the fit grid. Measured two independent ways on 2025-01-15 (135 binned
#                   profiles = 2700 native 15 s profiles):
#                     - first_difference_sigma on the native profiles, propagated to the fit grid:
#                       8.57e-3 per native gate at 3 km -> 1.17e-4 for the night mean;
#                     - direct scatter of the night-mean signal over 12-15 km (pure noise, no
#                       atmospheric signal there): 1.14e-4 (std) / 1.14e-4 (first difference).
#                   The two agree to 3 %, so 1.2e-4 is the honest clean-night value.
#                   A noisy summer night (2026-05-22, 78 binned profiles) gives 5.0e-4.
#
#   NOTE / KNOWN BUG. The pipeline's own ``sigma_signal`` (calibration.py::_sigma_on_fit_grid) is
#   5.30e-4 on that same night, i.e. sqrt(300/15) = 4.472x TOO LARGE, because the propagation
#   divides by sqrt(n_binned_profiles) instead of sqrt(n_native_profiles). Measured here as
#   5.3042e-4 / 1.1861e-4 = 4.47. It affects the v2.2 chi2red GATES only, never the constant.
PAYERNE_CHM15K = dict(
    wavelength_nm=1064.0,
    station_alt_m=490.0,
    n_gates=512,
    dz_m=29.97,
    z_first_m=22.4775,
    C_true=6.5e11,
    sigma_night_clean=1.2e-4,      # signal units, night mean, 300 s x 30 m grid
    sigma_night_noisy=5.0e-4,
    sigma_native_profile=8.6e-3,   # signal units, one native 15 s x 15 m profile
    n_profiles_binned=135,
    n_profiles_native=2700,
    additive_residual_signal=2.5e-5,   # mean signal left at 12-15 km on 2025-01-15
)


def default_grid(n_gates: int = 512, dz_m: float = 29.97,
                 z_first_m: float = 22.4775) -> NDArray[np.float64]:
    """Return the real Payerne CHM15k post-binning range grid (m AGL, gate centres).

    This is the grid ``average_ceilometer_data`` produces for an L1 CHM15k night
    (300 s x 30 m L2 grid): 512 gates from 22.4775 m to 15337.1 m in steps of 29.97 m.
    Using the true grid rather than a round 30 m one matters because the window search
    converts metres to bin indices by flooring, so gate offsets shift window edges.
    """
    return z_first_m + dz_m * np.arange(int(n_gates), dtype=float)


# ---------------------------------------------------------------------------
# 1. Atmosphere
# ---------------------------------------------------------------------------
def make_atmosphere(
    z: NDArray[np.float64],
    station_alt_m: float = 490.0,
    t0: float = 288.15,
    p0: float = 1013.25,
    cams_file: Optional[Path] = None,
    *,
    wavelength_nm: float = 1064.0,
    latitude: Optional[float] = None,
    longitude: Optional[float] = None,
    t_start: Optional[np.datetime64] = None,
    t_end: Optional[np.datetime64] = None,
    std_atm_file: Optional[Path] = None,
) -> dict:
    """Build the molecular reference atmosphere on the instrument range grid.

    Uses the REPO's own molecular routine, so the ``p_mol`` returned here is bit-for-bit the
    profile the production window search consumes:

      * temperature / pressure: ``calibration.rayleigh.atmosphere.load_standard_atmosphere``
        (US Standard Atmosphere 1976, the operational default ``molecular_source='standard'``),
        or ``load_cams_atmosphere`` when ``cams_file`` is given;
      * scattering: ``calibration.rayleigh.atmosphere.calculate_molecular_properties``
        (Bucholtz 1995 cross-section, depolarisation 0.0301, ``alpha_mol = beta_mol * 8*pi/3``).

    Parameters
    ----------
    z : ndarray
        Range from the instrument, m AGL (gate centres). Use :func:`default_grid`.
    station_alt_m : float
        Station altitude, m ASL. The US-Std table is sampled at ``station_alt_m + z``.
    t0, p0 : float
        SEA-LEVEL reference temperature (K) and pressure (hPa). The US-Std profile is
        multiplied by ``t0 / 288.15`` and ``p0 / 1013.25`` respectively, so the defaults are
        an exact no-op. This is the only knob for "a warmer / lower-pressure night"; the
        SHAPE of the profile is always US-Std (or CAMS).
    cams_file : Path, optional
        A ``CAMS_Beta_*.nc`` file. When given, ``latitude``, ``longitude``, ``t_start`` and
        ``t_end`` are required and the CAMS T/p replaces the standard atmosphere.
    wavelength_nm : float
        Laser wavelength (1064 for CHM15k, 910 for CL31/CL51/CL61, 532 for MPL).

    Returns
    -------
    dict
        ``z`` (m AGL), ``altitude_asl`` (m), ``temperature`` (K), ``pressure`` (Pa),
        ``beta_mol`` (m^-1 sr^-1), ``alpha_mol`` (m^-1), ``p_mol`` (the range-normalised
        molecular power the window search regresses against), ``transmission_mol`` (two-way,
        molecular only), ``beta_att_mol``, plus ``wavelength_nm``, ``station_alt_m``,
        ``lidar_ratio_mol`` (= 8*pi/3 sr).
    """
    z = np.asarray(z, dtype=float)
    altitude_asl = station_alt_m + z

    if cams_file is not None:
        if latitude is None or longitude is None or t_start is None or t_end is None:
            raise ValueError("cams_file requires latitude, longitude, t_start and t_end")
        prof = load_cams_atmosphere(Path(cams_file), latitude, longitude,
                                    t_start, t_end, altitude_asl)
        if prof is None:
            raise ValueError(f"no usable CAMS T/p in {cams_file}")
    else:
        prof = load_standard_atmosphere(std_atm_file or DEFAULT_STANDARD_ATMOSPHERE, altitude_asl)

    # Optional uniform rescaling of the T/p profile (no-op with the default sea-level values).
    temperature = prof.temperature * (float(t0) / 288.15)
    pressure = prof.pressure * (float(p0) / 1013.25)

    mol = calculate_molecular_properties(temperature, pressure, z, wavelength_nm * 1e-9)

    return dict(
        z=z,
        altitude_asl=altitude_asl,
        temperature=temperature,
        pressure=pressure,
        beta_mol=mol.beta_mol,
        alpha_mol=mol.alpha_mol,
        p_mol=mol.p_mol,
        transmission_mol=mol.transmission,
        beta_att_mol=mol.beta_att_mol,
        wavelength_nm=float(wavelength_nm),
        station_alt_m=float(station_alt_m),
        lidar_ratio_mol=float(MOLECULAR_LIDAR_RATIO),
    )


# ---------------------------------------------------------------------------
# 2. Aerosol
# ---------------------------------------------------------------------------
def aerosol_layer(
    z: NDArray[np.float64],
    base_m: float,
    top_m: float,
    beta_peak: float,
    lidar_ratio_sr: float,
    *,
    taper_fraction: float = 0.25,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """A single aerosol layer: smooth top-hat between ``base_m`` and ``top_m``.

    The layer is flat at ``beta_peak`` in the middle and tapers to zero with a raised cosine
    (Hann) over the outer ``taper_fraction`` of its thickness at each edge. A hard-edged
    top-hat would put a step in ``beta_tot`` right where the window search evaluates a
    Rayleigh-shape residual, which produces gate behaviour that is an artefact of the step
    rather than of the aerosol; a real advected layer has metre-to-hundred-metre edges.

    Parameters
    ----------
    z : ndarray
        Range grid, m AGL.
    base_m, top_m : float
        Layer bottom and top, m AGL.
    beta_peak : float
        Peak aerosol backscatter, m^-1 sr^-1. For scale: a Payerne residual free-tropospheric
        haze is ~1e-7, a clear-night molecular beta_mol at 3 km / 1064 nm is ~6.7e-8, so
        beta_peak = 1e-8 is a scattering ratio of ~1.15 at 3 km.
    lidar_ratio_sr : float
        Aerosol extinction-to-backscatter ratio S_aer, sr (operational default 52).

    Returns
    -------
    (beta_aer, ext_aer)
        Both on ``z``: backscatter in m^-1 sr^-1 and extinction ``= beta_aer * S_aer`` in m^-1.
    """
    z = np.asarray(z, dtype=float)
    beta_aer = np.zeros_like(z)
    if not (top_m > base_m) or beta_peak == 0.0:
        return beta_aer, beta_aer.copy()

    thickness = float(top_m) - float(base_m)
    taper = max(float(taper_fraction), 0.0) * thickness
    inside = (z >= base_m) & (z <= top_m)
    shape = np.zeros_like(z)
    shape[inside] = 1.0
    if taper > 0:
        d_bot = (z - base_m) / taper
        d_top = (top_m - z) / taper
        edge = np.minimum(d_bot, d_top)
        ramp = inside & (edge < 1.0)
        # Hann half-window: 0 at the edge, 1 one taper length inside.
        shape[ramp] = 0.5 * (1.0 - np.cos(np.pi * np.clip(edge[ramp], 0.0, 1.0)))

    beta_aer = float(beta_peak) * shape
    ext_aer = beta_aer * float(lidar_ratio_sr)
    return beta_aer, ext_aer


def exponential_aerosol(
    z: NDArray[np.float64],
    beta_surface: float,
    scale_height_m: float,
    lidar_ratio_sr: float,
    *,
    top_m: Optional[float] = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Exponentially decaying haze, ``beta_aer(z) = beta_surface * exp(-z / scale_height_m)``.

    This is the shape a "hazy but cloud-free" night has in the free troposphere, and it is the
    one that can produce a smooth, monotonic signal/p_mol slope like the -19 %/km measured
    between 2 and 6 km on the v2.2-recovered nights. Optionally truncated above ``top_m``.
    Units as in :func:`aerosol_layer`.
    """
    z = np.asarray(z, dtype=float)
    beta_aer = float(beta_surface) * np.exp(-z / float(scale_height_m))
    if top_m is not None:
        beta_aer = np.where(z <= float(top_m), beta_aer, 0.0)
    return beta_aer, beta_aer * float(lidar_ratio_sr)


# ---------------------------------------------------------------------------
# 3. Forward signal
# ---------------------------------------------------------------------------
def forward_signal(
    z: NDArray[np.float64],
    C_true: float,
    beta_mol: NDArray[np.float64],
    beta_aer: Optional[NDArray[np.float64]] = None,
    ext_aer: Optional[NDArray[np.float64]] = None,
    lidar_ratio_mol: float = float(MOLECULAR_LIDAR_RATIO),
    overlap: Optional[NDArray[np.float64]] = None,
    background: float = 0.0,
    additive_residual: float = 0.0,
    noise_sigma=None,
    rng: Optional[np.random.Generator] = None,
) -> NDArray[np.float64]:
    """Synthesise ``rcs`` -- exactly the quantity the pipeline calls ``rcs`` (L1 ``rcs_0``).

    Implements equations [1]-[2] of the module docstring::

        ext_tot = beta_mol * lidar_ratio_mol + ext_aer
        T2(z)   = exp(-2 * INT_0^z ext_tot dz')
        rcs(z)  = C_true * O(z) * (beta_mol + beta_aer) * T2(z)
                  + background * z^2 + additive_residual + noise * z^2

    Integration convention. ``INT_0^z`` is a cumulative TRAPEZOID over the range grid, plus a
    leading rectangle ``ext_tot[0] * z[0]`` for the un-sampled segment between the instrument and
    the first gate centre. The pipeline's own optical depth
    (``rayleigh_fit.calculate_lidar_constant``) is the same cumulative trapezoid but starts its
    integral AT the first gate, i.e. it silently drops that leading rectangle. On a CHM15k grid
    (z[0] = 22.5 m) that is 2 * 7.9e-7 * 22.5 = 3.6e-5, so it costs 0.0036 % of C_L -- included
    here on purpose so the closure residual reports it rather than hiding it.

    Parameters
    ----------
    z : ndarray
        Range grid, m AGL. Must be strictly increasing.
    C_true : float
        The imposed lidar constant. Units are whatever makes rcs match the instrument's L1
        ``rcs_0``; ~6.5e11 for a Payerne CHM15k.
    beta_mol, beta_aer : ndarray
        Molecular / aerosol backscatter, m^-1 sr^-1. ``beta_aer=None`` means a pure molecular
        atmosphere.
    ext_aer : ndarray, optional
        Aerosol extinction, m^-1. Defaults to zeros; note it is NOT derived from ``beta_aer``
        here -- pass the pair returned by :func:`aerosol_layer` so the lidar ratio is explicit.
    lidar_ratio_mol : float
        Molecular extinction-to-backscatter ratio, sr. The physical value is 8*pi/3 = 8.3776
        and is what the retrieval assumes; exposed so a mechanism study can perturb it.
    overlap : ndarray, optional
        Overlap function O(z), dimensionless, defaults to 1 everywhere. Above ~1.5 km every
        instrument here is in full overlap, so it is irrelevant to the fit band -- but it does
        change the optical-depth integral the retrieval starts at the ground.
    background : float
        Constant offset of the range-normalised signal (rcs/z^2), i.e. an un-subtracted detector
        background or stray light. Adds ``background * z^2`` to rcs. Same units as ``signal``.
    additive_residual : float
        Constant offset of ``rcs`` ITSELF (flat in range-corrected units), e.g. an afterpulse
        residual. Adds ``additive_residual`` to rcs.
    noise_sigma : float or ndarray, optional
        1-sigma photon noise OF THE SIGNAL (rcs/z^2), per gate; scalar or one value per gate.
        Gaussian noise ``N(0, noise_sigma) * z^2`` is added to rcs. This is the same quantity as
        the pipeline's ``sigma_signal``; on a clean Payerne CHM15k night-mean profile it is
        1.2e-4 and essentially flat with range (background-limited regime).
    rng : numpy Generator, optional
        Source of the noise. Pass a seeded generator for reproducibility.

    Returns
    -------
    ndarray
        ``rcs`` on ``z``.
    """
    z = np.asarray(z, dtype=float)
    beta_mol = np.asarray(beta_mol, dtype=float)
    beta_aer = np.zeros_like(z) if beta_aer is None else np.asarray(beta_aer, dtype=float)
    ext_aer = np.zeros_like(z) if ext_aer is None else np.asarray(ext_aer, dtype=float)
    overlap = np.ones_like(z) if overlap is None else np.asarray(overlap, dtype=float)

    beta_tot = beta_mol + beta_aer
    ext_tot = beta_mol * float(lidar_ratio_mol) + ext_aer

    # Cumulative optical depth from the GROUND (see the integration-convention note above).
    optical_depth = np.empty_like(z)
    optical_depth[0] = ext_tot[0] * z[0]
    optical_depth[1:] = optical_depth[0] + np.cumsum(
        0.5 * np.diff(z) * (ext_tot[:-1] + ext_tot[1:]))
    two_way_transmission = np.exp(-2.0 * optical_depth)

    rcs = float(C_true) * overlap * beta_tot * two_way_transmission

    if background:
        rcs = rcs + float(background) * z ** 2
    if additive_residual:
        rcs = rcs + float(additive_residual)
    if noise_sigma is not None:
        generator = rng if rng is not None else np.random.default_rng()
        sigma = np.broadcast_to(np.asarray(noise_sigma, dtype=float), z.shape)
        rcs = rcs + generator.normal(0.0, 1.0, size=z.shape) * sigma * z ** 2
    return rcs


# ---------------------------------------------------------------------------
# 4. Retrieval (the SHIPPED chain)
# ---------------------------------------------------------------------------
def default_options(molecular_method: str = "eprof_v2") -> CalibrationOptions:
    """Operational options (``options.json`` at the repo root) with the given window method.

    Guarantees the same ``lidar_ratio_aerosol`` (52 sr), ``subtract_background`` (False),
    ``consider_points_lower_than_molecular`` (True), ``half_length_options_m`` and search band
    (2000-6000 m) as production. A FRESH object every call -- deliberately not cached, so a
    caller that tweaks one gate cannot silently change every later retrieval.
    """
    options = CalibrationOptions.from_json(REPO / "options.json")
    options.molecular_method = str(molecular_method)
    options.plot_main = options.plot_all = False
    return options


def _forced_fit_result(
    signal: NDArray[np.float64],
    p_mol: NDArray[np.float64],
    z: NDArray[np.float64],
    bottom_m: float,
    top_m: float,
) -> RayleighFitResult:
    """Build the ``RayleighFitResult`` for a window imposed by hand.

    The free fit ``signal = a * p_mol + b`` is computed with the same ``scipy.stats.linregress``
    call ``molecular_methods.compute_window_grid`` uses, over the same in-window points, so a
    forced window is byte-comparable with one the search would have returned. ``relative_error``
    is the production metric |a - median(signal/p_mol)| / median * 100.
    """
    mask = (z >= bottom_m) & (z <= top_m) & np.isfinite(signal) & np.isfinite(p_mol)
    if mask.sum() < 5:
        raise ValueError(f"forced window {bottom_m:.0f}-{top_m:.0f} m has < 5 usable gates")
    x, y = p_mol[mask], signal[mask]
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    median_ratio = float(np.median(y / x))
    rel_error = abs((slope - median_ratio) / median_ratio * 100.0) if median_ratio else np.inf
    centre = 0.5 * (float(bottom_m) + float(top_m))
    half = 0.5 * (float(top_m) - float(bottom_m))
    return RayleighFitResult(
        slope=float(slope), intercept=float(intercept), r_squared=float(r_value ** 2),
        std_error=float(std_err), p_value=float(p_value),
        center_range_m=centre, half_length_m=half,
        range_start_m=float(bottom_m), range_end_m=float(top_m),
        altitude_start=0.0, altitude_end=0.0, relative_error=float(rel_error),
    )


def retrieve(
    z: NDArray[np.float64],
    rcs: NDArray[np.float64],
    atmosphere: dict,
    options: Optional[CalibrationOptions] = None,
    window: Optional[tuple] = None,
    *,
    perturbations: bool = False,
    sigma_signal: Optional[NDArray[np.float64]] = None,
    signal_stack: Optional[NDArray[np.float64]] = None,
) -> dict:
    """Run the SHIPPED Rayleigh retrieval on one synthetic night-mean profile.

    Pipeline functions called, in order (nothing is reimplemented):

      1. ``signal = rcs / z**2``                              -- as ``calibrate_rayleigh`` step 6;
      2. ``rayleigh_fit.find_optimal_molecular_window``       -- only when ``window is None``;
         otherwise the window is imposed and its free fit is computed with the same
         ``linregress`` call the production grid search uses (see ``_forced_fit_result``);
      3. ``calibration._compute_cl_for_perturbation``         -- which runs
         ``atmosphere.klett_inversion`` then ``rayleigh_fit.calculate_lidar_constant``
         (median of ``rcs / beta_tot * exp(2*OD)`` over the window, ``subtract_background``
         taken from ``options``).

    Parameters
    ----------
    z, rcs : ndarray
        Range grid (m AGL) and the synthetic range-corrected signal from :func:`forward_signal`.
    atmosphere : dict
        Output of :func:`make_atmosphere` (must be on the same ``z``).
    options : CalibrationOptions, optional
        Defaults to :func:`default_options` = operational ``options.json`` with ``eprof_v2``.
    window : (bottom_m, top_m), optional
        Force the molecular window (m AGL) instead of searching. This is what makes the
        C_L-versus-height ladder possible.
    perturbations : bool
        False (default) -> the nominal configuration only (LR = ``options.lidar_ratio_aerosol``,
        altitude shift 0). True -> reproduce production's LR x altitude-shift grid
        (``LR_DELTAS`` x ``ALT_SHIFTS_M``) and also return their median and robust spread.

    Returns
    -------
    dict
        ``ok`` (bool), ``message``, ``C_L`` (nominal), ``C_L_median`` and ``C_L_uncertainty``
        (only with ``perturbations=True``, else NaN), ``window_bottom`` / ``window_top`` /
        ``window_centre`` (m AGL), ``*_asl`` equivalents, ``slope``, ``intercept``,
        ``r_squared``, ``rel_error``, ``n_gates``, and the diagnostic profiles ``cl_profile``,
        ``beta_tot``, ``ext_tot``, ``signal``, plus ``fit_result``.
    """
    z = np.asarray(z, dtype=float)
    rcs = np.asarray(rcs, dtype=float)
    options = options if options is not None else default_options()
    p_mol = np.asarray(atmosphere["p_mol"], dtype=float)
    beta_mol = np.asarray(atmosphere["beta_mol"], dtype=float)
    if p_mol.shape != z.shape:
        raise ValueError("atmosphere is on a different grid than z")

    signal = rcs / (z ** 2)

    failed = dict(ok=False, C_L=np.nan, C_L_median=np.nan, C_L_uncertainty=np.nan,
                  window_bottom=np.nan, window_top=np.nan, window_centre=np.nan,
                  window_bottom_asl=np.nan, window_top_asl=np.nan,
                  slope=np.nan, intercept=np.nan, r_squared=np.nan, rel_error=np.nan,
                  n_gates=0, cl_profile=None, beta_tot=None, ext_tot=None,
                  signal=signal, fit_result=None)

    if window is None:
        fit_result = find_optimal_molecular_window(
            signal=signal, p_mol=p_mol, range_alc=z,
            half_length_options_m=options.half_length_options_m,
            range_start_m=options.range_start_m, range_end_m=options.range_end_m,
            increment_bins=options.fit_range_increment_bins,
            min_window_start_m=options.min_window_start_m,
            min_r2=options.min_window_r2, max_rel_error=options.max_window_rel_error,
            method=getattr(options, "molecular_method", "eprof_v2"),
            signal_stack=signal_stack,
            method_params=getattr(options, "molecular_params", None) or None,
            sigma_signal=sigma_signal,
        )
        if not np.isfinite(fit_result.slope):
            return dict(failed, message="no molecular window passed the validity gates")
    else:
        try:
            fit_result = _forced_fit_result(signal, p_mol, z, float(window[0]), float(window[1]))
        except ValueError as exc:
            return dict(failed, message=str(exc))

    lidar_ratios = [options.lidar_ratio_aerosol]
    shifts = [0.0]
    if perturbations:
        lidar_ratios = [options.lidar_ratio_aerosol + d for d in LR_DELTAS
                        if options.lidar_ratio_aerosol + d > 0]
        shifts = list(ALT_SHIFTS_M)

    values = []
    nominal = None
    for lidar_ratio in lidar_ratios:
        for shift in shifts:
            result = _compute_cl_for_perturbation(
                rcs_mean=rcs, range_alc=z, beta_mol=beta_mol, fit_result=fit_result,
                lidar_ratio_aerosol=lidar_ratio, altitude_shift_m=shift,
                subtract_background=options.subtract_background,
                consider_points_lower_than_molecular=options.consider_points_lower_than_molecular,
                sign_error_v10=getattr(options, "sign_error_v10", False),
                return_diagnostics=(lidar_ratio == options.lidar_ratio_aerosol and shift == 0.0),
            )
            if result is None:
                continue
            values.append(result.lidar_constant)
            if lidar_ratio == options.lidar_ratio_aerosol and shift == 0.0:
                nominal = result

    if nominal is None:
        return dict(failed, message="the nominal Klett / lidar-constant configuration failed")

    cl_median, cl_uncertainty = np.nan, np.nan
    if perturbations and values:
        arr = np.asarray(values, float)
        cl_median = float(np.median(arr))
        q25, q75 = np.percentile(arr, [25, 75])
        sigma_iqr = (q75 - q25) / 1.349
        sigma_mad = float(np.median(np.abs(arr - cl_median))) / 0.6745
        cl_uncertainty = 2.0 * max(sigma_iqr, sigma_mad)

    in_window = (z >= fit_result.range_start_m) & (z <= fit_result.range_end_m)
    station_alt = float(atmosphere.get("station_alt_m", 0.0))
    return dict(
        ok=True, message="ok",
        C_L=float(nominal.lidar_constant),
        C_L_median=cl_median, C_L_uncertainty=cl_uncertainty,
        window_bottom=float(fit_result.range_start_m),
        window_top=float(fit_result.range_end_m),
        window_centre=float(fit_result.center_range_m),
        window_bottom_asl=float(fit_result.range_start_m) + station_alt,
        window_top_asl=float(fit_result.range_end_m) + station_alt,
        slope=float(fit_result.slope), intercept=float(fit_result.intercept),
        r_squared=float(fit_result.r_squared), rel_error=float(fit_result.relative_error),
        n_gates=int(in_window.sum()),
        cl_profile=nominal.cl_profile, beta_tot=nominal.beta_tot, ext_tot=nominal.ext_tot,
        signal=signal, fit_result=fit_result,
    )


def scan_windows(
    z: NDArray[np.float64],
    rcs: NDArray[np.float64],
    atmosphere: dict,
    centres_m: Sequence[float],
    half_width_m: float,
    *,
    options: Optional[CalibrationOptions] = None,
    perturbations: bool = False,
) -> list:
    """The C_L(window centre) ladder -- the object the whole altitude study is about.

    For each centre in ``centres_m`` the molecular window is FORCED to
    ``[centre - half_width_m, centre + half_width_m]`` (m AGL) and the shipped retrieval is run.
    A perfectly altitude-independent calibration returns the same constant at every centre;
    the measured network gradient is ~7 %/km in |dC_L/dz| and up to -15 %/km at Payerne.

    Returns
    -------
    list of (centre_m, C_L)
        ``C_L`` is NaN where the retrieval failed at that centre.
    """
    ladder = []
    for centre in centres_m:
        result = retrieve(z, rcs, atmosphere,
                          options=options,
                          window=(float(centre) - float(half_width_m),
                                  float(centre) + float(half_width_m)),
                          perturbations=perturbations)
        value = (result["C_L_median"] if (perturbations and result["ok"])
                 else (result["C_L"] if result["ok"] else np.nan))
        ladder.append((float(centre), float(value)))
    return ladder


# ---------------------------------------------------------------------------
# 5. Closure tests
# ---------------------------------------------------------------------------
def _ladder_slope_pct_per_km(centres_m, values, reference) -> float:
    """Least-squares slope of (C_L/reference - 1) versus window centre, in % per km."""
    centres = np.asarray(centres_m, float)
    relative = np.asarray(values, float) / float(reference) - 1.0
    good = np.isfinite(relative)
    if good.sum() < 3:
        return np.nan
    slope = np.polyfit(centres[good] / 1000.0, relative[good] * 100.0, 1)[0]
    return float(slope)


def run_closure_tests(
    fig_path: Optional[Path] = None,
    n_seeds: int = 200,
    centres_m: Sequence[float] = tuple(np.arange(2500.0, 6501.0, 500.0)),
    half_width_m: float = 490.0,
    verbose: bool = True,
) -> dict:
    """Run the three closure tests and (optionally) write the closure figure.

    Test 1 -- PURE MOLECULAR, noiseless, overlap = 1, background = 0. The ladder must be flat and
             equal to ``C_true``. Any deviation is an ESTIMATOR defect, not atmosphere.
    Test 2 -- same, plus realistic photon noise (sigma from a real Payerne CHM15k night, see
             ``PAYERNE_CHM15K``), ``n_seeds`` realisations. Bias (mean over seeds) and scatter
             (std over seeds) are reported SEPARATELY per window height.
    Test 3 -- linearity: doubling ``C_true`` must double every retrieved C_L exactly.

    Returns a dict of numbers (all relative deviations in %, slopes in %/km).
    """
    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atmosphere = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                                 wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    beta_mol = atmosphere["beta_mol"]
    centres = np.asarray(centres_m, float)

    # ---- Test 1: noiseless pure molecular -------------------------------------------------
    rcs_clean = forward_signal(z, C_true, beta_mol)
    ladder_clean = scan_windows(z, rcs_clean, atmosphere, centres, half_width_m)
    clean_values = np.array([v for _, v in ladder_clean], float)
    clean_dev_pct = (clean_values / C_true - 1.0) * 100.0
    clean_max_abs = float(np.nanmax(np.abs(clean_dev_pct)))
    clean_slope = _ladder_slope_pct_per_km(centres, clean_values, C_true)

    # The free window search on the same profile (what production would have chosen).
    searched = retrieve(z, rcs_clean, atmosphere)
    searched_dev_pct = ((searched["C_L"] / C_true - 1.0) * 100.0) if searched["ok"] else np.nan

    # ---- Test 2: realistic photon noise ---------------------------------------------------
    sigma = PAYERNE_CHM15K["sigma_night_clean"]
    rng = np.random.default_rng(20260815)
    noisy = np.full((int(n_seeds), centres.size), np.nan)
    for seed in range(int(n_seeds)):
        rcs_noisy = forward_signal(z, C_true, beta_mol, noise_sigma=sigma, rng=rng)
        noisy[seed, :] = [v for _, v in scan_windows(z, rcs_noisy, atmosphere,
                                                     centres, half_width_m)]
    noisy_bias_pct = (np.nanmean(noisy, axis=0) / C_true - 1.0) * 100.0
    noisy_scatter_pct = np.nanstd(noisy, axis=0) / C_true * 100.0
    noisy_slope = _ladder_slope_pct_per_km(centres, np.nanmean(noisy, axis=0), C_true)
    # Standard error of the mean bias, so "is the bias significant?" is answerable.
    noisy_bias_sem_pct = noisy_scatter_pct / np.sqrt(np.sum(np.isfinite(noisy), axis=0))

    # ---- Test 3: linearity in C_true ------------------------------------------------------
    rcs_double = forward_signal(z, 2.0 * C_true, beta_mol)
    ladder_double = scan_windows(z, rcs_double, atmosphere, centres, half_width_m)
    double_values = np.array([v for _, v in ladder_double], float)
    with np.errstate(invalid="ignore"):
        linearity_ratio = double_values / clean_values
    linearity_max_err = float(np.nanmax(np.abs(linearity_ratio - 2.0)))

    results = dict(
        C_true=C_true,
        centres_m=centres.tolist(),
        half_width_m=float(half_width_m),
        window_gates=int(retrieve(z, rcs_clean, atmosphere,
                                  window=(centres[0] - half_width_m,
                                          centres[0] + half_width_m))["n_gates"]),
        clean_dev_pct=clean_dev_pct.tolist(),
        clean_max_abs_dev_pct=clean_max_abs,
        clean_slope_pct_per_km=clean_slope,
        searched_window=(searched["window_bottom"], searched["window_top"]) if searched["ok"] else None,
        searched_dev_pct=float(searched_dev_pct),
        noise_sigma_signal=float(sigma),
        n_seeds=int(n_seeds),
        noisy_bias_pct=noisy_bias_pct.tolist(),
        noisy_bias_sem_pct=noisy_bias_sem_pct.tolist(),
        noisy_scatter_pct=noisy_scatter_pct.tolist(),
        noisy_slope_pct_per_km=noisy_slope,
        linearity_ratio=linearity_ratio.tolist(),
        linearity_max_abs_err=linearity_max_err,
    )

    if verbose:
        print(f"C_true = {C_true:.4e}   grid {z.size} gates, dz = {np.diff(z)[0]:.2f} m, "
              f"window half-width {half_width_m:.0f} m ({results['window_gates']} gates)")
        print("\nTEST 1 - pure molecular, noiseless, O=1, background=0")
        for centre, dev in zip(centres, clean_dev_pct):
            print(f"   centre {centre:6.0f} m AGL   C_L/C_true - 1 = {dev:+9.5f} %")
        print(f"   max |deviation| = {clean_max_abs:.5f} %   residual slope = "
              f"{clean_slope:+.5f} %/km")
        if searched["ok"]:
            print(f"   free window search picked {searched['window_bottom']:.0f}-"
                  f"{searched['window_top']:.0f} m AGL -> {searched_dev_pct:+.5f} %")
        print(f"\nTEST 2 - photon noise sigma_signal = {sigma:.2e} (flat), {n_seeds} seeds")
        for centre, bias, sem, scatter in zip(centres, noisy_bias_pct, noisy_bias_sem_pct,
                                              noisy_scatter_pct):
            print(f"   centre {centre:6.0f} m AGL   bias = {bias:+7.3f} +- {sem:.3f} %   "
                  f"scatter = {scatter:6.3f} %")
        print(f"   mean-ladder residual slope = {noisy_slope:+.4f} %/km")
        print("\nTEST 3 - C_true doubled")
        print(f"   max |C_L(2C)/C_L(C) - 2| = {linearity_max_err:.3e}")

    if fig_path is not None:
        _closure_figure(Path(fig_path), z, atmosphere, rcs_clean, sigma, C_true, centres,
                        clean_values, noisy, half_width_m)
        results["fig_path"] = str(fig_path)

    return results


def _closure_figure(fig_path: Path, z, atmosphere, rcs_clean, sigma, C_true, centres,
                    clean_values, noisy, half_width_m) -> None:
    """Left: the synthetic profiles. Middle/right: the C_L ladder, noiseless then with noise.

    Two ladder panels because the two effects live three orders of magnitude apart: the
    noiseless closure is a few 1e-3 % and would be an invisible vertical line on the +-3 %
    axis the noise scatter needs. Altitude is on Y in all three panels.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    signal_clean = rcs_clean / z ** 2
    rng = np.random.default_rng(1)
    signal_noisy = forward_signal(z, C_true, atmosphere["beta_mol"],
                                  noise_sigma=sigma, rng=rng) / z ** 2

    clean_pct = (clean_values / C_true - 1.0) * 100.0
    bias_pct = (np.nanmean(noisy, axis=0) / C_true - 1.0) * 100.0
    scatter_pct = np.nanstd(noisy, axis=0) / C_true * 100.0
    sem_pct = scatter_pct / np.sqrt(np.sum(np.isfinite(noisy), axis=0))
    clean_slope = _ladder_slope_pct_per_km(centres, clean_values, C_true)
    noisy_slope = _ladder_slope_pct_per_km(centres, np.nanmean(noisy, axis=0), C_true)
    y_km = centres / 1000.0

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.4))

    # --- panel 1: the synthetic profiles ---------------------------------------------------
    ax = axes[0]
    ax.plot(signal_noisy, z / 1000.0, color="0.75", lw=0.8,
            label=f"synthetic night mean + noise ($\\sigma$={sigma:.1e})")
    ax.plot(signal_clean, z / 1000.0, color="C0", lw=1.8, label="synthetic, noiseless")
    ax.plot(C_true * atmosphere["p_mol"], z / 1000.0, color="C3", lw=1.0, ls="--",
            label="$C_{true}\\cdot p_{mol}$ (pipeline molecular reference)")
    ax.axhspan(2.0, 7.0, color="C2", alpha=0.07)
    ax.set_xscale("log")
    ax.set_xlim(1e-5, 5e0)
    ax.set_ylim(0, 10)
    ax.set_xlabel("range-normalised signal  rcs / $z^2$   [signal units]")
    ax.set_ylabel("range above the instrument  [km AGL]")
    ax.set_title("Synthetic pure-molecular night\n(Payerne CHM15k grid, 1064 nm)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)

    # --- panel 2: noiseless closure ladder, zoomed --------------------------------------------
    ax = axes[1]
    ax.axvline(0.0, color="C3", lw=1.6, ls="--", label="$C_{true}$ (reference)")
    ax.plot(clean_pct, y_km, "o-", color="C0", ms=6, lw=1.6, label="retrieved $C_L$")
    ax.set_xlim(-0.012, 0.012)
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title("TEST 1 - noiseless closure\n(note the $\\pm$0.012 % scale)", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(min(y_km) - 1.15, max(y_km) + 0.35)
    ax.text(0.03, 0.02,
            f"max |dev| = {np.nanmax(np.abs(clean_pct)):.4f} %\n"
            f"slope = {clean_slope:+.5f} %/km\n"
            f"(the $-$0.0034 % offset is the optical depth\n"
            f"below the first gate, z[0] = {z[0]:.1f} m)",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    # --- panel 3: with photon noise ----------------------------------------------------------
    ax = axes[2]
    ax.axvline(0.0, color="C3", lw=1.6, ls="--", label="$C_{true}$ (reference)")
    ax.errorbar(bias_pct, y_km, xerr=scatter_pct, fmt="none", ecolor="0.65", capsize=3, lw=1.0,
                label="per-night scatter, 1$\\sigma$")
    ax.errorbar(bias_pct, y_km, xerr=sem_pct, fmt="s", color="0.25", ms=5, capsize=2, lw=1.4,
                label=f"mean bias $\\pm$ s.e.m. ({noisy.shape[0]} seeds)")
    ax.set_xlabel("$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre  [km AGL]")
    ax.set_title(f"TEST 2 - with photon noise\n($\\sigma$ = {sigma:.1e}, clean Payerne night)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(min(y_km) - 1.15, max(y_km) + 0.35)
    ax.text(0.03, 0.02,
            f"bias slope = {noisy_slope:+.3f} %/km\n"
            f"scatter {scatter_pct[0]:.2f} % at {y_km[0]:.1f} km "
            f"-> {scatter_pct[-1]:.2f} % at {y_km[-1]:.1f} km\n"
            f"observed Payerne gradient: $-$15 %/km",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    fig.suptitle(f"Rayleigh forward-model closure: the shipped retrieval on a known $C_{{true}}$ "
                 f"= {C_true:.2e}   |   forced windows of $\\pm${half_width_m:.0f} m",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    out = REPO / "doc" / "reports" / "figs_altitude_audit" / "fwd_closure.png"
    run_closure_tests(fig_path=out)
    print(f"\nfigure -> {out}")
