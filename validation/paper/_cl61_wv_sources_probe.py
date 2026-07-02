"""WV-source sensitivity for the CL61 Rayleigh calibration (Payerne):
T2_wv at each night's molecular fit window from
  (a) CAMS 1 deg monthly  (operational baseline),
  (b) Payerne radiosonde (00 UT, D:/Soundings/sounding_pay_2026.csv, rh/T -> n_wv),
  (c) CAMS 0.4 deg monthly (D:/CAMS_Monthly_04 — only months 202501/202502/202506/202507 exist).
The Rayleigh C_L scales as 1/T2_wv(window) -> dC_L(%) = 100*(T2_a/T2_x - 1)."""
import sys, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import pandas as pd

from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, DEFAULT_ABS_CROSS_SECTION)
from calibration.io.cams import ensure_cams_file
from calibration.cloud.calibration import _nw_from_T_RH

LAT, LON, ALT = 46.813, 6.943, 491.0
LAM0, FWHM = 910.74, 1.0
NIGHTS = [("20260316", 2630, 6010), ("20260328", 2659, 4599), ("20260402", 3000, 5000),
          ("20260407", 3000, 5000), ("20260410", 3000, 5000), ("20260417", 3000, 5000),
          ("20260422", 3000, 5000), ("20260424", 3000, 5000), ("20260428", 3000, 5000),
          ("20260521", 3000, 5000), ("20260602", 3000, 5000)]
SND = pd.read_csv("D:/Soundings/sounding_pay_2026.csv")
SND["dt"] = pd.to_datetime(SND["t"] - 719529, unit="D")   # MATLAB datenum -> pandas
# recent launches are stored as per-day files (the yearly compilation stops mid-March)
def sounding_day(ds):
    f = Path(f"D:/Soundings/sounding_pay_{ds}.csv")
    if f.is_file():
        s = pd.read_csv(f)
        if "RH" in s.columns:            # per-day files name the column RH (yearly: rh)
            s = s.rename(columns={"RH": "rh"})
        s["dt"] = pd.to_datetime(s["t"] - 719529, unit="D")
        return s
    day = pd.Timestamp(f"{ds[:4]}-{ds[4:6]}-{ds[6:8]}")
    return SND[(SND["dt"] >= day - pd.Timedelta(hours=2)) & (SND["dt"] <= day + pd.Timedelta(hours=28))]

zgrid = np.arange(0, 8001, 30.0) + ALT

def t2_from_cams(folder, ds):
    cams = ensure_cams_file(Path(folder), ds, auto_download=False)
    if cams is None:
        return None
    m0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-01")
    prof = cams_water_vapor_profile(cams, LAT, LON, m0, m0 + np.timedelta64(27, "D"))
    if prof is None:
        return None
    h, n = prof
    return np.asarray(two_way_wv_transmission(zgrid, ALT, h, n, DEFAULT_ABS_CROSS_SECTION, LAM0, FWHM), "f8")

def t2_from_sounding(ds):
    s = sounding_day(ds)
    s = s.dropna(subset=["z", "T", "rh"])
    s = s[(s["z"] >= ALT - 20) & (s["z"] <= 12000)].sort_values("z")
    # keep the FIRST launch of the day (00 UT) — the calibration nights end ~04 UT
    if len(s) < 20:
        return None
    t0 = s["dt"].min()
    s = s[s["dt"] <= t0 + pd.Timedelta(hours=3)]
    if len(s) < 20:
        return None
    nw = _nw_from_T_RH(s["T"].to_numpy(), s["rh"].to_numpy())   # rh already in percent
    return np.asarray(two_way_wv_transmission(zgrid, ALT, s["z"].to_numpy(), nw,
                                              DEFAULT_ABS_CROSS_SECTION, LAM0, FWHM), "f8")

print("night      window[km]   T2(1deg)  T2(sonde)  T2(0.4deg)   dCL sonde   dCL 0.4deg")
d_sonde, d_04 = [], []
for ds, zb, zt in NIGHTS:
    zm = (zgrid - ALT >= zb) & (zgrid - ALT <= zt)
    t2a = t2_from_cams("D:/CAMS", ds)
    t2s = t2_from_sounding(ds)
    t204 = t2_from_cams("D:/CAMS_Monthly_04", ds[:6] + "01")
    va = float(np.nanmean(t2a[zm])) if t2a is not None else np.nan
    vs = float(np.nanmean(t2s[zm])) if t2s is not None else np.nan
    v4 = float(np.nanmean(t204[zm])) if t204 is not None else np.nan
    ds_pct = 100 * (va / vs - 1) if np.isfinite(vs) else np.nan
    d4_pct = 100 * (va / v4 - 1) if np.isfinite(v4) else np.nan
    if np.isfinite(ds_pct):
        d_sonde.append(ds_pct)
    if np.isfinite(d4_pct):
        d_04.append(d4_pct)
    print(f"{ds}  {zb/1000:.1f}-{zt/1000:.1f}      {va:.3f}     {vs if np.isfinite(vs) else float('nan'):.3f}      "
          f"{v4 if np.isfinite(v4) else float('nan'):.3f}      {ds_pct:+6.1f}%    {d4_pct:+6.1f}%")
print(f"\nmedian dC_L if the sounding were truth : {np.median(d_sonde):+.1f}%  (n={len(d_sonde)})")
print(f"median dC_L 1deg vs 0.4deg (same month) : {np.median(d_04):+.1f}%  (n={len(d_04)})" if d_04 else
      "0.4deg: no overlapping months among the calibration nights (only 202501/02, 202506/07 downloaded)")
# generic resolution check on the months that DO exist at both resolutions
print("\ngeneric 1deg-vs-0.4deg check (fixed 3-5 km window):")
zm = (zgrid - ALT >= 3000) & (zgrid - ALT <= 5000)
for ds in ("20250115", "20250215", "20260615"):
    a = t2_from_cams("D:/CAMS", ds); b = t2_from_cams("D:/CAMS_Monthly_04", ds[:6] + "01")
    if a is None or b is None:
        print(f"  {ds[:6]}: missing"); continue
    print(f"  {ds[:6]}: T2(1deg)={np.nanmean(a[zm]):.3f}  T2(0.4deg)={np.nanmean(b[zm]):.3f}  "
          f"dC_L={100*(np.nanmean(a[zm])/np.nanmean(b[zm])-1):+.1f}%")
print("WV_SOURCES_DONE")
