"""
Atmospheric physics calculations for Rayleigh calibration.

This module contains functions for calculating molecular backscatter,
extinction, and related atmospheric parameters using the Bucholtz (1995)
formulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d


# Physical constants
STANDARD_TEMPERATURE = 288.15  # K (15°C)
STANDARD_PRESSURE = 101325.0   # Pa
DEPOLARIZATION_RATIO = 0.0301
MOLECULAR_DENSITY_STP = 2.547e25  # m^-3
GAS_CONSTANT = 8.314418  # J/(mol·K)
MOLECULAR_LIDAR_RATIO = 8 * np.pi / 3  # sr


@dataclass
class AtmosphericProfile:
    """Temperature and pressure profiles on instrument range grid."""
    temperature: NDArray[np.float64]  # K
    pressure: NDArray[np.float64]     # Pa
    altitude: NDArray[np.float64]     # m ASL


@dataclass
class MolecularProperties:
    """Molecular scattering properties calculated from atmospheric profiles."""
    beta_mol: NDArray[np.float64]       # Molecular backscatter coefficient (m^-1 sr^-1)
    alpha_mol: NDArray[np.float64]      # Molecular extinction coefficient (m^-1)
    beta_att_mol: NDArray[np.float64]   # Attenuated molecular backscatter
    p_mol: NDArray[np.float64]          # Molecular power (range-normalized)
    transmission: NDArray[np.float64]   # Two-way transmission


def calculate_refractive_index(wavelength_m: float) -> float:
    """
    Calculate the refractive index of air at a given wavelength.

    Uses the Edlén formula for the refractive index of standard air.

    Parameters
    ----------
    wavelength_m : float
        Wavelength in meters.

    Returns
    -------
    float
        Refractive index of air.
    """
    wavelength_um = wavelength_m * 1e6  # Convert to micrometers
    sigma2 = (1 / wavelength_um) ** 2   # Wavenumber squared

    # Edlén formula coefficients
    m = (5791817 / (238.0185 - sigma2) + 167909 / (57.362 - sigma2)) * 1e-8 + 1
    return m


def calculate_rayleigh_phase_function(depol_ratio: float = DEPOLARIZATION_RATIO) -> float:
    """
    Calculate the Rayleigh phase function at 180° (backscatter direction).

    Parameters
    ----------
    depol_ratio : float
        Depolarization ratio (default: 0.0301 for air).

    Returns
    -------
    float
        Phase function value at 180°.
    """
    gamma = depol_ratio / (2 - depol_ratio)
    # Phase function at theta = pi (backscatter)
    cos_theta_sq = np.cos(np.pi) ** 2
    p_ray = 3 / 4 / (1 + 2 * gamma) * (1 + 3 * gamma + (1 - gamma) * cos_theta_sq)
    return p_ray


def calculate_molecular_properties(
    temperature: NDArray[np.float64],
    pressure: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    wavelength_m: float,
) -> MolecularProperties:
    """
    Calculate molecular scattering properties using Bucholtz (1995) formulation.

    This is a vectorized implementation that avoids explicit Python loops
    for better performance.

    Parameters
    ----------
    temperature : ndarray
        Temperature profile in Kelvin.
    pressure : ndarray
        Pressure profile in Pascals.
    range_alc : ndarray
        Range bins from instrument in meters.
    wavelength_m : float
        Laser wavelength in meters.

    Returns
    -------
    MolecularProperties
        Dataclass containing all molecular scattering properties.
    """
    # Calculate refractive index and phase function
    m = calculate_refractive_index(wavelength_m)
    p_ray = calculate_rayleigh_phase_function()

    # Bucholtz (1995) formula for total molecular scattering cross-section
    # Vectorized calculation
    numerator = 24 * np.pi**3 * ((m**2 - 1)**2) * (6 + 3 * DEPOLARIZATION_RATIO)
    denominator = (
        wavelength_m**4 * MOLECULAR_DENSITY_STP**2 *
        (m**2 + 2)**2 * (6 - 7 * DEPOLARIZATION_RATIO)
    )

    # Total molecular backscatter coefficient (before phase function normalization)
    beta_mol_total = (
        numerator / denominator *
        MOLECULAR_DENSITY_STP *
        STANDARD_TEMPERATURE * pressure / STANDARD_PRESSURE / temperature
    )

    # Backscatter coefficient (normalized by phase function)
    beta_mol = beta_mol_total / (4 * np.pi) * p_ray

    # Extinction coefficient (related to backscatter by molecular lidar ratio)
    alpha_mol = beta_mol * MOLECULAR_LIDAR_RATIO

    # Calculate two-way transmission using cumulative sum (vectorized)
    dz = np.abs(range_alc[1] - range_alc[0]) if len(range_alc) > 1 else 1.0

    # Cumulative optical depth (two-way)
    optical_depth = np.cumsum(alpha_mol) * dz
    # Shift to get integral from 0 to z (not including current bin)
    optical_depth = np.insert(optical_depth[:-1], 0, 0)
    transmission = np.exp(-2 * optical_depth)

    # Attenuated molecular backscatter
    beta_att_mol = beta_mol * transmission

    # Molecular power (range-normalized signal)
    p_mol = beta_mol * transmission / (range_alc**2 + 1e-10)  # Avoid division by zero

    return MolecularProperties(
        beta_mol=beta_mol,
        alpha_mol=alpha_mol,
        beta_att_mol=beta_att_mol,
        p_mol=p_mol,
        transmission=transmission,
    )


def load_standard_atmosphere(
    filepath: Path,
    altitude_grid: NDArray[np.float64],
) -> AtmosphericProfile:
    """
    Load US Standard Atmosphere 1976 and interpolate to instrument grid.

    Parameters
    ----------
    filepath : Path
        Path to standard atmosphere CSV file.
    altitude_grid : ndarray
        Target altitude grid in meters ASL.

    Returns
    -------
    AtmosphericProfile
        Temperature and pressure interpolated to altitude grid.
    """
    data = np.genfromtxt(filepath, delimiter=',', names=True)

    # Ensure altitude_grid is a regular numpy array (not masked)
    if np.ma.isMaskedArray(altitude_grid):
        if np.any(np.ma.getmaskarray(altitude_grid)):
            raise ValueError("Cannot interpolate with masked altitude values")
        altitude_grid = altitude_grid.data

    temperature = interp1d(
        data['Altitude'], data['Temperature'],
        bounds_error=False, fill_value=np.nan
    )(altitude_grid)

    pressure = interp1d(
        data['Altitude'], data['Pressure'],
        bounds_error=False, fill_value=np.nan
    )(altitude_grid)

    return AtmosphericProfile(
        temperature=temperature,
        pressure=pressure,
        altitude=altitude_grid,
    )


def load_ecmwf_profile(
    filepath: Path,
    latitude: float,
    longitude: float,
    altitude_grid: NDArray[np.float64],
) -> Optional[AtmosphericProfile]:
    """
    Load ECMWF MACC reanalysis data and interpolate to instrument grid.

    Parameters
    ----------
    filepath : Path
        Path to ECMWF NetCDF file.
    latitude : float
        Station latitude.
    longitude : float
        Station longitude.
    altitude_grid : ndarray
        Target altitude grid in meters ASL.

    Returns
    -------
    AtmosphericProfile or None
        Temperature and pressure interpolated to altitude grid,
        or None if file not found.
    """
    from netCDF4 import Dataset

    if not filepath.exists():
        return None

    # Ensure altitude_grid is a regular numpy array
    if np.ma.isMaskedArray(altitude_grid):
        if np.any(np.ma.getmaskarray(altitude_grid)):
            raise ValueError("Cannot interpolate with masked altitude values")
        altitude_grid = altitude_grid.data

    with Dataset(filepath, 'r') as data:
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]

        # Find nearest grid point
        idx_lon = np.abs(lon - longitude).argmin()
        idx_lat = np.abs(lat - latitude).argmin()

        T_raw = data.variables['t'][:, :, idx_lat, idx_lon]
        level = data.variables['level'][:]
        z_ecmwf = data.variables['z'][:, :, idx_lat, idx_lon] / 9.80655  # Geopotential to height

    # Use first time step, interpolate to instrument grid
    # (original code used index [1, :] but that seems like a bug - using [0, :])
    T = interp1d(z_ecmwf[0, :], T_raw[0, :], bounds_error=False, fill_value=np.nan)(altitude_grid)
    P = interp1d(z_ecmwf[0, :], level, bounds_error=False, fill_value=np.nan)(altitude_grid) * 100  # hPa to Pa

    # Fill NaN at bottom with nearest valid value
    valid_idx = np.where(~np.isnan(T))[0]
    if len(valid_idx) > 0:
        first_valid = valid_idx[0]
        T[:first_valid] = T[first_valid]
        P[:first_valid] = P[first_valid]

    return AtmosphericProfile(
        temperature=T,
        pressure=P,
        altitude=altitude_grid,
    )


def klett_inversion(
    beta_att: NDArray[np.float64],
    beta_mol: NDArray[np.float64],
    range_alc: NDArray[np.float64],
    reference_index: int,
    lidar_ratio_aerosol: float,
    reference_value: float,
    i_start: int,
    i_end: int,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """
    Perform Klett inversion to retrieve aerosol backscatter and extinction.

    This implementation is optimized but maintains the original algorithm
    for scientific reproducibility.

    Parameters
    ----------
    beta_att : ndarray
        Attenuated backscatter profile.
    beta_mol : ndarray
        Molecular backscatter profile.
    range_alc : ndarray
        Range bins in meters.
    reference_index : int
        Index of reference altitude (in molecular region).
    lidar_ratio_aerosol : float
        Assumed aerosol lidar ratio in sr.
    reference_value : float
        Reference total backscatter at reference altitude.
    i_start : int
        Start index for inversion.
    i_end : int
        End index for inversion.

    Returns
    -------
    tuple of ndarray
        (beta_aerosol, beta_total, extinction_aerosol)
    """
    n_bins = len(beta_att)
    dz = np.abs(range_alc[1] - range_alc[0]) if len(range_alc) > 1 else 1.0

    beta_aer = np.zeros(n_bins)
    beta_tot = np.zeros(n_bins)
    ext_aer = np.zeros(n_bins)

    lr_diff = lidar_ratio_aerosol - MOLECULAR_LIDAR_RATIO

    # --- Vectorised Klett inversion over [i_start, i_end) ------------------
    # Pre-compute the cumulative integral of beta_mol * lr_diff * dz from
    # index 0 upward.  The integral from bin k to reference_index is then
    #   cum[reference_index] - cum[k].
    mol_dz = beta_mol * lr_diff * dz          # (n_bins,)
    cum = np.zeros(n_bins + 1)                # cum[k] = sum(mol_dz[0:k])
    np.cumsum(mol_dz, out=cum[1:])

    # qt[R] = sum(beta_mol[R:ref] * lr_diff * dz) = cum[ref] - cum[R]
    R_idx = np.arange(i_start, i_end)
    qt = cum[reference_index] - cum[R_idx]     # (n_R,)
    numerator = beta_att[R_idx] * np.exp(-2 * qt)  # (n_R,)

    # T_factor[r] = -2 * (cum[ref] - cum[r])  for each bin r
    T_factor_all = -2.0 * (cum[reference_index] - cum[i_start:reference_index])
    weighted_all = beta_att[i_start:reference_index] * np.exp(T_factor_all) * lidar_ratio_aerosol * dz

    # denominator_sum[R] = sum over r in [R, reference_index) of weighted_all[r]
    # Use reverse cumulative sum so that rev_cumsum[k] = sum(weighted_all[k:end])
    rev_cumsum = np.cumsum(weighted_all[::-1])[::-1]  # length = ref - i_start
    # For R_idx, the offset into weighted_all is R - i_start
    denom_sum = np.empty(len(R_idx))
    offset = R_idx - i_start
    valid_off = offset < len(rev_cumsum)
    denom_sum[valid_off] = rev_cumsum[offset[valid_off]]
    denom_sum[~valid_off] = 0.0

    denominator = reference_value + 2 * denom_sum

    beta_aer[R_idx] = numerator / denominator - beta_mol[R_idx]
    beta_tot[R_idx] = beta_aer[R_idx] + beta_mol[R_idx]
    ext_aer[R_idx] = np.maximum(0, beta_aer[R_idx]) * lidar_ratio_aerosol

    # Fill below start index with first valid value
    if i_start > 0:
        beta_aer[:i_start] = beta_aer[i_start]
        beta_tot[:i_start] = beta_tot[i_start]
        ext_aer[:i_start] = ext_aer[i_start]

    return beta_aer, beta_tot, ext_aer
