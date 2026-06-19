"""
Rayleigh Calibration Package for Automated Lidars and Ceilometers.

This package provides tools for performing Rayleigh (molecular) calibration
of ceilometers and lidars, following the methodology described in various
E-PROFILE publications.

Main Features
-------------
- Support for multiple instrument types (CHM15k, CL51, CL61, Mini-MPL)
- Automatic detection of optimal molecular scattering window
- Klett inversion for extinction retrieval
- Quality control and validation
- NetCDF output compatible with E-PROFILE standards

Example Usage
-------------
>>> from rayleigh_calibration import calibrate_rayleigh, load_instruments, CalibrationOptions
>>> 
>>> # Load instrument configuration
>>> instruments = load_instruments("instruments.json")
>>> options = CalibrationOptions.from_json("options.json")
>>> 
>>> # Run calibration for a single instrument
>>> result = calibrate_rayleigh("20240115", instruments[0], options)
>>> print(f"Lidar constant: {result.lidar_constant:.4e}")

Authors
-------
Original code: hem (2015)
Modernized version: Claude/Anthropic (2024)
"""

__version__ = "2.0.0"
__author__ = "E-PROFILE"

from .config import (
    InstrumentType,
    DataLevel,
    InstrumentInfo,
    CalibrationOptions,
    CalibrationResult,
    load_instruments,
)

from .rayleigh.calibration import calibrate_rayleigh

from .rayleigh.atmosphere import (
    DEFAULT_STANDARD_ATMOSPHERE,
    calculate_molecular_properties,
    load_standard_atmosphere,
    load_ecmwf_profile,
    MolecularProperties,
    AtmosphericProfile,
)

from .io.data_loader import (
    CeilometerData,
    load_l1_data,
    load_data,
    build_file_paths,
)

from .rayleigh.rayleigh_fit import (
    find_optimal_molecular_window,
    RayleighFitResult,
)


__all__ = [
    # Main function
    "calibrate_rayleigh",
    # Configuration
    "InstrumentType",
    "DataLevel",
    "InstrumentInfo",
    "CalibrationOptions",
    "CalibrationResult",
    "load_instruments",
    # Atmosphere
    "calculate_molecular_properties",
    "load_standard_atmosphere",
    "DEFAULT_STANDARD_ATMOSPHERE",
    "load_ecmwf_profile",
    "MolecularProperties",
    "AtmosphericProfile",
    # Data
    "CeilometerData",
    "load_l1_data",
    "load_data",
    "build_file_paths",
    # Fitting
    "find_optimal_molecular_window",
    "RayleighFitResult",
]
