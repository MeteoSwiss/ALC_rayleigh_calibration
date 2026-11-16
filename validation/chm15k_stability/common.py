# -*- coding: utf-8 -*-
"""Shared paths and loaders for the CHM15k calibration-stability study (report 19).

The study asks whether the CHM15k Rayleigh lidar constant C_L depends on the instrument temperature
and whether it has a diurnal / short-term pattern. All analysis scripts read the 10-min housekeeping
+ binned-signal cache built by `hk_cache.py` and the operational nightly constants published by
E-PROFILE (dashboard NetCDF, one file per stream and year).
"""
from __future__ import annotations

import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202610_CHM15k_stability")
CACHE = OUT / "hk"
FIGS = REPO / "doc/reports/figs_chm15k_stability"
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
NC_DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/calib_dashboard/nc")
CENSUS = REPO / "validation/scope_l1_2026_census.json"
YEARS = (2025, 2026)

# 10-min blocks: 144 per day. Range bins: 100 m to 1 km (overlap / near range), 250 m above.
BLOCKS_PER_DAY = 144
RANGE_EDGES = np.r_[np.arange(0, 1000, 100), np.arange(1000, 8001, 250)].astype(float)
RANGE_MID = 0.5 * (RANGE_EDGES[:-1] + RANGE_EDGES[1:])

# Housekeeping variables kept per 10-min block (median). Temperatures are converted to degC.
HK_VARS = ("temp_int", "temp_ext", "temperature_detector", "temperature_optical_module",
           "calibration_pulse", "window_transmission", "status_laser", "status_detector",
           "laser_life_time", "laser_pulses", "bckgrd_rcs_0", "stddev",
           # Vaisala (CL31 / CL61) housekeeping; NaN in CHM15k files (and absent from caches built before)
           "temperature_laser", "laser_energy")
TEMP_VARS = ("temp_int", "temp_ext", "temperature_detector", "temperature_optical_module", "temperature_laser")

# Same-wavelength co-located CHM15k pairs (aerosol cancels in the ratio, so day and night work).
PAIR_SITES = {
    "0-20000-0-06240": ("A", "B", "C", "D"),   # Amsterdam Schiphol (units a few km apart)
    "0-20000-0-12843": ("A", "B"),             # Budapest
    "0-20008-0-INO": ("A", "B"),               # Magurele / Magurele RADO
}

EPOCH = datetime(1970, 1, 1)


def chm15k_streams():
    """All CHM15k streams of the L1 census as (wmo, ident, site, lat, lon)."""
    out = []
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["type"].upper().startswith("CHM"):
            out.append((s["wmo"], s["ident"], s["site"], s["lat"], s["lon"]))
    return out


def stream_key(wmo, ident):
    return f"{wmo}_{ident}"


def cache_file(wmo, ident):
    return CACHE / f"{stream_key(wmo, ident)}.npz"


def load_cache(wmo, ident):
    """Return the cache as a dict of arrays, or None when the stream was not cached."""
    f = cache_file(wmo, ident)
    if not f.exists():
        return None
    with np.load(f, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def load_constants(wmo, ident):
    """Operational nightly Rayleigh constants (calibration_method 0) of one stream, 2025-2026.

    Returns a dict of arrays: t, t0, t1 (days since 1970), c (lidar constant), u (uncertainty),
    zb, zt (calibration window, m ASL), module (optical_module_id string), version.
    """
    import netCDF4
    keys = ("t", "t0", "t1", "c", "u", "zb", "zt", "module")
    acc = {k: [] for k in keys}
    for y in YEARS:
        f = NC_DIR / f"ALC_calibration_{wmo}_{ident}{y}.nc"
        if not f.exists():
            continue
        with netCDF4.Dataset(f) as ds:
            m = np.ma.filled(ds["calibration_method"][:], -1) == 0
            if not m.any():
                continue
            acc["t"].append(np.ma.filled(ds["time"][:], np.nan)[m])
            acc["t0"].append(np.ma.filled(ds["start_time"][:], np.nan)[m])
            acc["t1"].append(np.ma.filled(ds["end_time"][:], np.nan)[m])
            acc["c"].append(np.ma.filled(ds["lidar_constant"][:].astype("f8"), np.nan)[m])
            acc["u"].append(np.ma.filled(ds["lidar_constant_uncertainty"][:].astype("f8"), np.nan)[m])
            acc["zb"].append(np.ma.filled(ds["calibration_bottom_height"][:].astype("f8"), np.nan)[m])
            acc["zt"].append(np.ma.filled(ds["calibration_top_height"][:].astype("f8"), np.nan)[m])
            mod = ds["optical_module_id"][:]
            mod = np.array([str(x) for x in np.asarray(mod)])[m]
            acc["module"].append(mod)
    if not acc["t"]:
        return None
    out = {k: np.concatenate(v) for k, v in acc.items()}
    order = np.argsort(out["t"])
    return {k: v[order] for k, v in out.items()}


def days_to_datetime(t):
    return np.array([EPOCH + timedelta(days=float(x)) for x in np.atleast_1d(t)])


def solar_hour(t_days, lon):
    """Local mean solar time in hours for UTC days-since-1970."""
    return np.mod((np.asarray(t_days) % 1.0) * 24.0 + lon / 15.0, 24.0)
