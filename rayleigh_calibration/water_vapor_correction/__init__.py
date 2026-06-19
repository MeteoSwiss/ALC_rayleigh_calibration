"""Water-vapour absorption correction for 910 nm ALCs.

Spectral two-way water-vapour transmission from CAMS monthly means and a HITRAN
absorption cross-section LUT; a port of the MATLAB ``wv_t2eff`` /
``compute_wv_transmission`` routines.
"""

from .water_vapor import (
    in_water_vapor_band,
    laser_spectrum_for,
    cams_water_vapor_profile,
    two_way_wv_transmission,
    wv_t2eff_core,
    load_abs_cross_section,
)

__all__ = [
    "in_water_vapor_band",
    "laser_spectrum_for",
    "cams_water_vapor_profile",
    "two_way_wv_transmission",
    "wv_t2eff_core",
    "load_abs_cross_section",
]
