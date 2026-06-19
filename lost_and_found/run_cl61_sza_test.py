#!/usr/bin/env python3
"""Compare the new SZA-based (darkness-adaptive) night selection against the old
fixed solar-clock window, on Payerne CL61 2026. Lower per-night CV and/or more
usable nights = the SZA window is the better averaging choice."""
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
info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
BASE = Path("A:/E-PROFILE_L2_monthly")
dates = [f"2026{m:02d}{d:02d}" for m in range(1, 6) for d in range(1, cal.monthrange(2026, m)[1] + 1)]


def make_opts(use_sza, thresh=100.0):
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = BASE
    o.data_level = DataLevel.L2_MONTHLY
    o.plot_all = o.plot_main = False
    o.folder_output = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/cl61_sza_test")
    o.use_sza_night = use_sza
    o.sza_night_threshold = thresh
    return o


variants = [("clock 20->04 (old)", make_opts(False)),
            ("SZA > 96 (civil)", make_opts(True, 96)),
            ("SZA > 100 (nautical, new default)", make_opts(True, 100)),
            ("SZA > 102", make_opts(True, 102)),
            ("SZA > 108 (astronomical)", make_opts(True, 108))]

print("Payerne CL61 2026 — night selection: SZA-based vs clock\n")
print(f"{'method':<36}{'n_valid':>8}{'CV%':>8}{'median CL':>11}")
for name, o in variants:
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
        print(f"{name:<36}{cls.size:>8}{100*np.std(cls)/np.mean(cls):>8.1f}{np.median(cls):>11.4f}")
    else:
        print(f"{name:<36}{0:>8}{'--':>8}")
print("\nDONE")
