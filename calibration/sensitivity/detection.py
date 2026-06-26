"""Aerosol detection thresholds from a measured noise floor.

Ports the detection-threshold block of ``dark_measurement_cl61_chm_cl31.m``
(and the identical block in the ambient script).

Hypothesis: the per-gate temporal noise sigma_dark(r) of attenuated
backscatter (Mm^-1 sr^-1) is the limiting factor. The minimum *detectable*
attenuated backscatter at range r for a required SNR and averaging time tau is

    beta_att_min(r) = SNR * sigma(r, tau),   sigma(r, tau) = sigma0(r) * sqrt(dt/tau)

(white-noise tau^-1/2 scaling, valid below the Allan drift knee). To recover
the *true* aerosol backscatter the instrument can detect we divide by the
two-way molecular transmission T_mol^2 (aerosol two-way transmission neglected
at threshold - the thin-layer limit):

    beta_min(r) = beta_att_min(r) / T_mol^2(r)
    alpha_min(r) = LR * beta_min(r)                       [Mm^-1]
    M_min(r)     = alpha_min(r) / MEC                     [-> ug/m^3]

LR  = lidar ratio [sr];  MEC = mass extinction coefficient [m^2/g].
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "VolcanicScenario",
    "VOLCANIC_SCENARIOS",
    "ICAO_LEVELS_UG_M3",
    "ICAO_LABELS",
    "H_SCALE_DEFAULT",
    "ALPHA_MOL0",
    "alpha_mol0_for_wavelength",
    "two_way_molecular_transmission",
    "scale_sigma_to_tau",
    "min_detectable_backscatter",
    "min_detectable_extinction",
    "min_detectable_mass",
    "detection_altitude",
]


@dataclass(frozen=True)
class VolcanicScenario:
    """Optical/microphysical assumptions for volcanic ash."""

    name: str
    LR: float   # lidar ratio [sr]
    MEC: float  # mass extinction coefficient [m^2/g]


# EARLINET-style fine volcanic ash (Ansmann et al. 2011; Eyjafjallajokull 2010).
VOLCANIC_SCENARIOS: tuple[VolcanicScenario, ...] = (
    VolcanicScenario("Ash: LR=60 sr, MEC=0.60 m^2/g", 60.0, 0.60),
    VolcanicScenario("Ash: LR=50 sr, MEC=0.80 m^2/g", 50.0, 0.80),
)

# ICAO EUR/NAT Volcanic Ash Contingency Plan (EUR Doc 019 / NAT Doc 006,
# 2010 edition; London VAAC ash-concentration charts):
#   low    : 0.2e-3 .. 2e-3 g/m^3  (200 .. 2000 ug/m^3)
#   medium : 2e-3   .. 4e-3 g/m^3  (2000 .. 4000 ug/m^3)
#   high   : > 4e-3 g/m^3          (> 4000 ug/m^3)
ICAO_LEVELS_UG_M3: tuple[float, ...] = (200.0, 2000.0, 4000.0)
ICAO_LABELS: tuple[str, ...] = (
    "200 ug/m3 (low edge)",
    "2000 ug/m3 (low/medium)",
    "4000 ug/m3 (medium/high)",
)

H_SCALE_DEFAULT = 8000.0  # atmospheric scale height [m]

# Standard-atmosphere molecular extinction at the surface, per wavelength
# [m^-1]. 910 and 1064 nm are the MATLAB values; 532 nm is scaled by the
# Rayleigh lambda^-4 law from the 1064 nm anchor (5e-7 * (1064/532)^4).
ALPHA_MOL0: Dict[int, float] = {
    532: 5.0e-7 * (1064.0 / 532.0) ** 4,  # ~8.0e-6 m^-1
    910: 1.0e-6,
    1064: 5.0e-7,
}


def alpha_mol0_for_wavelength(wavelength_nm: float) -> float:
    """Surface molecular extinction for the nearest tabulated wavelength, or a
    lambda^-4 scaling from 1064 nm for anything else."""
    wl = int(round(wavelength_nm))
    if wl in ALPHA_MOL0:
        return ALPHA_MOL0[wl]
    nearest = min(ALPHA_MOL0, key=lambda k: abs(k - wl))
    if abs(nearest - wl) <= 60:  # 910<->910.5, 1064<->1064.5, 532
        return ALPHA_MOL0[nearest]
    return 5.0e-7 * (1064.0 / wavelength_nm) ** 4


def two_way_molecular_transmission(
    r: NDArray, alpha_mol0: float, h_scale: float = H_SCALE_DEFAULT
) -> NDArray:
    """T_mol^2(r) = exp(-2 * integral_0^r alpha_mol dr') with an exponential
    Rayleigh profile alpha_mol(z) = alpha_mol0 * exp(-z/h_scale)."""
    r = np.asarray(r, dtype=float)
    return np.exp(-2.0 * alpha_mol0 * h_scale * (1.0 - np.exp(-r / h_scale)))


def scale_sigma_to_tau(sigma0: NDArray, dt: float, tau: float) -> NDArray:
    """White-noise tau^-1/2 scaling: sigma(tau) = sigma0 * sqrt(dt/max(tau,dt)).

    ``sigma0`` is the per-gate noise at the native sampling time ``dt``.
    """
    return np.asarray(sigma0, dtype=float) * np.sqrt(dt / max(tau, dt))


def min_detectable_backscatter(
    sigma0: NDArray,
    r: NDArray,
    dt: float,
    tau: float,
    snr: float,
    alpha_mol0: float,
    h_scale: float = H_SCALE_DEFAULT,
) -> NDArray:
    """Minimum detectable *true aerosol* backscatter [Mm^-1 sr^-1]."""
    sig_tau = scale_sigma_to_tau(sigma0, dt, tau)
    tmol2 = two_way_molecular_transmission(r, alpha_mol0, h_scale)
    return snr * sig_tau / tmol2


def min_detectable_extinction(beta_min: NDArray, lr: float) -> NDArray:
    """alpha_min = LR * beta_min  [Mm^-1] (beta_min in Mm^-1 sr^-1)."""
    return lr * np.asarray(beta_min, dtype=float)


def min_detectable_mass(beta_min: NDArray, lr: float, mec: float) -> NDArray:
    """Minimum detectable mass concentration [ug/m^3].

    alpha_min[Mm^-1] -> *1e-6 -> m^-1; M[g/m^3] = alpha/MEC; *1e6 -> ug/m^3.
    """
    alpha_min_per_m = lr * np.asarray(beta_min, dtype=float) * 1e-6
    return (alpha_min_per_m / mec) * 1e6


def detection_altitude(
    m_min_ug: NDArray, r: NDArray, level_ug: float
) -> float:
    """Highest altitude [m] at which the minimum detectable mass is still
    below ``level_ug`` (ug/m^3). Returns NaN if the threshold is never met.

    This is the network-map headline scalar: how high can the instrument see a
    given ICAO ash concentration.
    """
    m_min_ug = np.asarray(m_min_ug, dtype=float)
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(m_min_ug) & (m_min_ug <= level_ug)
    if not np.any(ok):
        return np.nan
    return float(np.nanmax(r[ok]))
