"""Observation-minus-Background (OmB) for E-PROFILE ALCs against CAMS.

Compares cloud-free, time-averaged attenuated backscatter from the ceilometer
network against the CAMS aerosol-backscatter forecast, reproducing the MATLAB
``E_PROFILE_ALC_Monthly_OB.m`` methodology in Python and reusing the existing
calibration infrastructure (CAMS level reader, water-vapour correction).

* :mod:`calibration.omb.cams_aerosol` - read CAMS ground-referenced aerosol
  backscatter at a station, at the instrument wavelength (532 / 1064 nm native,
  910 nm via the Angstrom exponent between 532 and 1064 nm).
* :mod:`calibration.omb.omb` - average + cloud-screen the observation, apply the
  910 nm water-vapour correction, interpolate onto the CAMS grid and form the
  bias (observation minus background) for one or more calibration sources
  (operational L2 constant and / or our Kalman best-estimate C_L).
"""

from . import cams_aerosol, omb  # noqa: F401
