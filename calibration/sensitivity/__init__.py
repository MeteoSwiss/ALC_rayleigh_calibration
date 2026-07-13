"""Instrument-sensitivity / noise-floor analysis for E-PROFILE ALCs.

This package ports the noise-characterisation and aerosol-detection-threshold
methodology from the MATLAB scripts

    dark_measurement_cl61_chm_cl31.m   (hood-on dark noise, full range)
    ambient_noise_cl61_chm_cl31.m      (clear-sky ambient noise, estimator b)

into reusable, instrument-agnostic kernels:

* :mod:`calibration.sensitivity.noise`     - per-gate noise estimators
  (temporal first difference = estimator "b"), robust statistics, altitude
  binning, the solar day/night split and the overlapping Allan deviation with
  a chi-squared confidence interval.
* :mod:`calibration.sensitivity.detection` - the minimum-detectable
  backscatter / extinction / mass-concentration thresholds at a given SNR and
  averaging time, the two-way molecular transmission and the ICAO volcanic-ash
  contamination levels.

The :mod:`validation.dark_measurement_payerne` script exercises these kernels
on Payerne L1 data so that the same code is validated before it feeds the
network-wide sensitivity product on the monitoring dashboard.
"""

from . import noise, detection  # noqa: F401
