#!/usr/bin/env python3
"""Test how processing choices reduce the CL61 Rayleigh per-night scatter.

The Allan analysis showed the per-night Rayleigh constant is heavily noise-dominated
(sigma~39 % at 1 day). Here we test the *processing* levers on the same Payerne CL61
2026 nights and measure the resulting night-to-night CV and the number of usable nights:

  baseline      : options.json defaults (fit 2-6 km, quality 15, >=3 h night)
  wide_vertical : extend the molecular fit window (2-9 km) + larger half-lengths
                  -> more vertical averaging per night
  strict_quality: tighten the method-agreement gate (quality 15 -> 8)
                  -> noise filtering (keep only clean nights)
  long_night    : require >=5 h of night -> more temporal averaging per night
"""
from __future__ import annotations
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import calendar as cal
import logging, warnings
from pathlib import Path
import numpy as np
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.CRITICAL)

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "C", 46.8137, 6.9425, 491.0
BASE = Path("A:/E-PROFILE_L2_monthly")
info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
dates = [f"2026{m:02d}{d:02d}" for m in range(1, 6) for d in range(1, cal.monthrange(2026, m)[1] + 1)]


def make_opts(**kw):
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = BASE
    o.data_level = DataLevel.L2_MONTHLY
    o.plot_all = o.plot_main = False
    o.folder_output = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/cl61_param_test")  # writable
    for k, v in kw.items():
        setattr(o, k, v)
    return o


variants = {
    "baseline":       make_opts(),
    "wide_vertical":  make_opts(range_end_m=9000.0,
                                half_length_options_m=(250, 490, 730, 970, 1450, 1930, 2400, 3000)),
    "strict_quality": make_opts(threshold_quality=8.0),
    "long_night":     make_opts(min_time_range=5),
}

print(f"CL61 Payerne Rayleigh processing-parameter test ({len(dates)} candidate nights, 2026)\n")
print(f"{'variant':<15}{'n_valid':>8}{'CV%':>8}{'median CL':>11}")
for name, o in variants.items():
    cls = []
    for ds in dates:
        try:
            r = calibrate_rayleigh(ds, info, o)
            if r.flag in (1, 1.0, 0.5) and r.lidar_constant > 0:
                cls.append(r.lidar_constant)
        except Exception:
            pass
    cls = np.array(cls)
    if cls.size:
        cv = 100 * np.std(cls) / np.mean(cls)
        print(f"{name:<15}{cls.size:>8}{cv:>8.1f}{np.median(cls):>11.4f}")
    else:
        print(f"{name:<15}{0:>8}{'--':>8}")
print("\nDONE")
