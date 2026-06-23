"""Liquid-water-cloud calibration (O'Connor / Hopkin method).

Python port of the MATLAB ``liquid_cloud_calibration`` implementation, bit-for-bit
validated against the reference. Reads E-PROFILE L2 or Cloudnet CL61 raw and applies
the mandatory water-vapour correction.
"""

from .calibration import (
    liquid_cloud_calibration,
    liquid_cloud_calibration_from_data,
    CloudCalConfig,
    CloudCalResults,
)

__all__ = [
    "liquid_cloud_calibration",
    "liquid_cloud_calibration_from_data",
    "CloudCalConfig",
    "CloudCalResults",
]
