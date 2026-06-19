#!/usr/bin/env python3
"""Rayleigh is night-only: the nighttime averaging window sets the per-night noise.
Current pipeline = FIXED solar-clock window hour_min=20 -> hour_max=4 (8 h).

(1) How much DARK time is actually available at Payerne by month (solar zenith angle
    > 100 deg = nautical night)? -> what the fixed 8 h window leaves unused.
(2) Does WIDENING the window reduce the per-night Rayleigh scatter / recover nights?
    Re-run Payerne CL61 2026 with 8 h / 12 h / 14 h / 16 h windows; report CV and
    the number of usable nights.
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
import datetime as dt
import logging, warnings
from pathlib import Path
import numpy as np
from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from rayleigh_calibration.config import InstrumentType

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.CRITICAL)
WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "C", 46.8137, 6.9425, 491.0
info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                      instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
BASE = Path("A:/E-PROFILE_L2_monthly")


def solar_decl(doy):
    return np.radians(-23.44) * np.cos(2 * np.pi * (doy + 10) / 365.25)


def dark_hours(lat_deg, doy, sza_thresh=100.0):
    """Hours per day with solar zenith angle > sza_thresh (nautical night ~100 deg)."""
    lat = np.radians(lat_deg)
    d = solar_decl(doy)
    cz = np.cos(np.radians(sza_thresh))
    # cos(SZA) = sin(lat)sin(dec) + cos(lat)cos(dec)cos(H); dark when cos(SZA) < cz
    x = (cz - np.sin(lat) * np.sin(d)) / (np.cos(lat) * np.cos(d))
    if x <= -1:
        return 0.0          # never dark enough (white night)
    if x >= 1:
        return 24.0         # always dark (polar night)
    H = np.degrees(np.arccos(x))
    return 24.0 * (180.0 - H) / 180.0   # fraction of day darker than threshold


print("=== (1) available nautical-dark hours at Payerne (SZA>100 deg) vs the fixed 8 h window ===")
print(f"{'month':>6}{'dark h (SZA>100)':>18}{'currently used':>16}")
for mo in range(1, 13):
    doy = dt.date(2026, mo, 15).timetuple().tm_yday
    dh = dark_hours(LAT, doy)
    print(f"{mo:>6}{dh:>16.1f} h{min(8.0, dh):>14.1f} h")

print("\n=== (2) per-night Rayleigh CV vs night-window width (Payerne CL61, 2026) ===")
dates = [f"2026{m:02d}{d:02d}" for m in range(1, 6) for d in range(1, cal.monthrange(2026, m)[1] + 1)]


def make_opts(hmin, hmax):
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = BASE
    o.data_level = DataLevel.L2_MONTHLY
    o.plot_all = o.plot_main = False
    o.folder_output = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/cl61_nightwin_test")
    o.hour_min = hmin
    o.hour_max = hmax
    return o


windows = [("8 h (current 20->04)", 20, 4),
           ("12 h (18->06)", 18, 6),
           ("14 h (17->07)", 17, 7),
           ("16 h (16->08)", 16, 8)]
print(f"{'window':<22}{'n_valid':>8}{'CV%':>8}{'median CL':>11}")
for name, hmin, hmax in windows:
    o = make_opts(hmin, hmax)
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
        print(f"{name:<22}{cls.size:>8}{100*np.std(cls)/np.mean(cls):>8.1f}{np.median(cls):>11.4f}")
    else:
        print(f"{name:<22}{0:>8}{'--':>8}")
print("\nDONE")
