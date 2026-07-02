"""Aosta CL61 (0-380-5-1, ident B): does the dark-offset explanation apply, and did something
change around March 2026? Monthly offset proxy = median over nights of the nightly-mean (21-04 UT)
beta_att in the 9-13 km band (atmosphere there ~ +0.02*C_L: tiny; changes trace the instrument
offset). Plus: distinct global attributes (firmware/software/serial) of the L1 files per month."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

WMO, IDENT = "0-380-5-1", "B"
L1 = Path("D:/E-PROFILE_L1_2026")
MONTHS = ["202511", "202512", "202601", "202602", "202603", "202604", "202605", "202606"]
ATTRS = ("software_version", "firmware_version", "history", "source", "instrument_serial_number",
         "instrument_firmware_version", "comment", "title")

attr_seen = {}
print("month    nights  offset proxy 9-13 km [Mm-1 sr-1]   3-6 km residual-ish")
for ym in MONTHS:
    vals36, vals913 = [], []
    for day in range(2, 29, 4):                       # ~7 nights/month
        ds = f"{ym}{day:02d}"
        f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            for a in ATTRS:
                if hasattr(nc, a):
                    attr_seen.setdefault(a, {}).setdefault(str(getattr(nc, a))[:90], []).append(ym)
            tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
            tv = np.asarray(nc.variables["time"][:], "f8")
            frac = (tv - np.floor(tv)) if "day" in tu else (tv % 86400) / 86400.0
            ni = (frac >= 21 / 24) | (frac <= 4 / 24)
            if ni.sum() < 100:
                continue
            rng = np.asarray(nc.variables["range"][:], "f8")
            x = np.asarray(nc.variables["rcs_0"][:], "f8")
            if x.shape[0] != tv.size:
                x = x.T
            m = np.nanmean(x[ni], axis=0) * 1e6
        z1 = (rng >= 9000) & (rng <= 13000); z2 = (rng >= 3000) & (rng <= 6000)
        if np.nanmax(m[(rng >= 6000) & (rng <= 13000)]) > 1.0:   # cirrus guard
            continue
        vals913.append(float(np.nanmean(m[z1]))); vals36.append(float(np.nanmean(m[z2])))
    if vals913:
        print(f"{ym}   {len(vals913):4d}    {np.median(vals913):+8.4f}                     {np.median(vals36):+8.4f}")
print("\ndistinct file attributes and the months they appear in:")
for a, d in attr_seen.items():
    if len(d) > 1:
        for v, ms in d.items():
            print(f"  {a}: '{v}' -> {sorted(set(ms))}")
    else:
        v = next(iter(d))
        print(f"  {a}: constant ('{v[:60]}...')" if len(v) > 60 else f"  {a}: constant ('{v}')")
print("AOSTA_HISTORY_DONE")
