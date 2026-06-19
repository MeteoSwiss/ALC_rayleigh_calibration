"""Core Rayleigh (molecular) calibration.

The molecular-calibration pipeline (:func:`calibrate_rayleigh`), the molecular-window
grid search, the atmosphere/molecular-properties model, and the selectable
molecular-window detection methods.
"""

from .calibration import calibrate_rayleigh
from .rayleigh_fit import (
    find_optimal_molecular_window,
    calculate_lidar_constant,
    validate_calibration,
    RayleighFitResult,
)
from .atmosphere import (
    DEFAULT_STANDARD_ATMOSPHERE,
    load_standard_atmosphere,
    load_ecmwf_profile,
    load_cams_atmosphere,
    calculate_molecular_properties,
    klett_inversion,
    MolecularProperties,
    AtmosphericProfile,
    MOLECULAR_LIDAR_RATIO,
)

__all__ = [
    "calibrate_rayleigh",
    "find_optimal_molecular_window",
    "calculate_lidar_constant",
    "validate_calibration",
    "RayleighFitResult",
    "DEFAULT_STANDARD_ATMOSPHERE",
    "load_standard_atmosphere",
    "load_ecmwf_profile",
    "load_cams_atmosphere",
    "calculate_molecular_properties",
    "klett_inversion",
    "MolecularProperties",
    "AtmosphericProfile",
    "MOLECULAR_LIDAR_RATIO",
]
