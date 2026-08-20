# -*- coding: utf-8 -*-
"""Shared plumbing for the remote-dark campaign (M1 anchored fit + M3 modulation route).

Everything here is deliberately thin: paths, one L1 day reader that also returns the housekeeping
modulators M3 regresses on, and the solar elevation (imported from the v4 estimator so the two
campaigns can never disagree on what "night" or "twilight" means).

Outputs live OUTSIDE the repo, next to the other campaign data:
    C:/DATA/Projects/202606_E-PROFILE_calibration/remote_dark/
"""
from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from rayleigh_availability.dark_from_clearsky import solar_elevation_deg  # noqa: E402  (shared def)

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
OUT_DIR = DATA / "remote_dark"
V4_DIR = DATA / "rayleigh_availability" / "dark_clearsky"          # v4 night caches + estimates
HOOD_NPZ = DATA / "rayleigh_availability" / "dark_profiles_payerne.npz"   # the truth

L1_ROOT = Path("D:/E-PROFILE_L1_2026")

#: Payerne, the validation site. lat/lon feed the solar elevation; ident -> instrument type.
PAYERNE = {"wmo": "0-20000-0-06610", "lat": 46.813, "lon": 6.943,
           "types": {"A": "CHM15k", "B": "CL31", "C": "CL61"}}

#: Housekeeping variables per profile, in preference order per field. Read tolerantly: a missing
#: variable yields NaN, never an exception — M3 treats "no modulator" as "this component is not
#: retrievable here", which is an answer, not an error.
HK_VARS = {
    "bckgrd": ["bckgrd_rcs_0"],
    "t_int": ["temp_int", "temperature_laser", "temperature_optical_module"],
    "laser": ["laser_energy", "laser_life_time"],
    "window": ["window_transmission"],
}


def ensure_out(sub: str = "") -> Path:
    p = OUT_DIR / sub if sub else OUT_DIR
    p.mkdir(parents=True, exist_ok=True)
    return p


def l1_file(wmo: str, ident: str, day: datetime) -> Path:
    """Payerne uses the year/month layout; network streams are flat. Try both."""
    a = L1_ROOT / wmo / f"{day:%Y}" / f"{day:%m}" / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
    if a.exists():
        return a
    return L1_ROOT / wmo / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"


def read_day(wmo: str, ident: str, day: datetime):
    """One L1 day -> dict(times, rng, rcs, hk) or None.

    rcs stays in rcs_0 units (range-corrected counts): that is the unit the whole dark campaign
    works in, and the unit in which the hood truth is stored. Gates at z <= 0 are NaN'd (the CL61
    first gate is 0 m and would blow up every /z^2 view downstream).
    """
    f = l1_file(wmo, ident, day)
    if not f.exists():
        return None
    try:
        with netCDF4.Dataset(f) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f8"), np.nan)
            hk = {}
            for field, cands in HK_VARS.items():
                v = None
                for nm in cands:
                    if nm in ds.variables:
                        a = np.ma.filled(ds.variables[nm][:].astype("f8"), np.nan).ravel()
                        if a.size == t.size and np.isfinite(a).any():
                            v = a
                            break
                hk[field] = v if v is not None else np.full(t.size, np.nan)
    except Exception:                                                       # noqa: BLE001
        return None
    if rcs.ndim != 2 or rcs.shape[1] != rng.size:
        return None
    times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
    rcs = np.where(rng[None, :] > 0, rcs, np.nan)
    return {"times": times, "rng": rng, "rcs": rcs, "hk": hk, "file": f}


def days_between(d0: datetime, d1: datetime):
    d = d0
    while d <= d1:
        yield d
        d += timedelta(days=1)


__all__ = ["DATA", "OUT_DIR", "V4_DIR", "HOOD_NPZ", "L1_ROOT", "PAYERNE",
           "ensure_out", "l1_file", "read_day", "days_between", "solar_elevation_deg"]
