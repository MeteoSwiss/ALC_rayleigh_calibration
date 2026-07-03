"""CL61 Rayleigh recalibration with RADIOSONDE water vapour (proper evaluation).

Three full eprof_v2 calibrations per night:
  (a) operational: monthly-CAMS WV, original L1;
  (b) sounding WV (Payerne 00 UT launch, rh/T -> n_wv), original L1;
  (c) sounding WV + measured covered-telescope offset b_dark(z) removed (best physics).
The sounding replaces the CAMS *WV profile* via injection into the calibration module
(everything else — molecular model, gates, Klett, window search — identical).
Metrics: per-night C_L, median gap to the cloud constant 1.425, night-to-night scatter.
"""
import sys, csv, shutil, tempfile, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset

from calibration.config import CalibrationOptions, InstrumentInfo, InstrumentType, DataLevel
import calibration.rayleigh.calibration as RC
from calibration.cloud.calibration import _nw_from_T_RH
from validation.paper.calib_benchmark import RAYLEIGH_METHOD, L1_ROOT, CAMS, REPO

WMO, IDENT, LAT, LON, ALT = "0-20000-0-06610", "C", 46.813, 6.943, 491.0
NIGHTS = ["20260225", "20260316", "20260328", "20260402", "20260407", "20260410", "20260411",
          "20260417", "20260420", "20260422", "20260423", "20260424", "20260428", "20260521", "20260602"]
CL_CLOUD = 1.4252

# ---- sounding provider (per-day files first, yearly compilation fallback) ----
_YR = pd.read_csv("D:/Soundings/sounding_pay_2026.csv")
_YR["dt"] = pd.to_datetime(_YR["t"] - 719529, unit="D")

def sounding_profile(ds):
    f = Path(f"D:/Soundings/sounding_pay_{ds}.csv")
    if f.is_file():
        s = pd.read_csv(f).rename(columns={"RH": "rh"})
        s["dt"] = pd.to_datetime(s["t"] - 719529, unit="D")
    else:
        day = pd.Timestamp(f"{ds[:4]}-{ds[4:6]}-{ds[6:8]}")
        s = _YR[(_YR["dt"] >= day - pd.Timedelta(hours=2)) & (_YR["dt"] <= day + pd.Timedelta(hours=28))]
    s = s.dropna(subset=["z", "T", "rh"])
    s = s[(s["z"] >= ALT - 20) & (s["z"] <= 12000)].sort_values("z")
    if len(s) < 20:
        return None
    t0 = s["dt"].min()
    s = s[s["dt"] <= t0 + pd.Timedelta(hours=3)]        # the 00 UT launch
    if len(s) < 20:
        return None
    nw = _nw_from_T_RH(s["T"].to_numpy(), s["rh"].to_numpy())   # rh in percent
    return s["z"].to_numpy(), nw

_orig_wvprof = RC.cams_water_vapor_profile
def _sounding_wvprof(cams_file, lat, lon, t_start, t_end):
    ds = str(np.datetime64(t_start, "D")).replace("-", "")
    p = sounding_profile(ds)
    return p if p is not None else _orig_wvprof(cams_file, lat, lon, t_start, t_end)

# ---- offset-corrected archive (rebuild) ----
dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
rng_d, b_dark = dk["rng"], dk["b_dark"]
b_s = dk["b_smooth"]   # robust per-gate median (330 m running median), see _cl61_dark_windows
CORR = Path(tempfile.mkdtemp(prefix="cl61_corr2_"))
for ds in NIGHTS:
    src = L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
    if not src.exists():
        continue
    dst = CORR / WMO / ds[:4] / ds[4:6] / src.name
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    with Dataset(dst, "r+") as nc:
        r = np.asarray(nc.variables["range"][:], "f8")
        corr = np.interp(r, rng_d, b_s) * 1e-6
        v = nc.variables["rcs_0"]; x = np.asarray(v[:], "f8")
        v[:] = (x - corr[None, :]) if x.shape[-1] == r.size else (x.T - corr[None, :]).T

def run(root, use_sounding):
    RC.cams_water_vapor_profile = _sounding_wvprof if use_sounding else _orig_wvprof
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = root; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = True
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                          instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
    out = {}
    for ds in NIGHTS:
        try:
            r = RC.calibrate_rayleigh(ds, info, o)
        except Exception:
            continue
        if r.flag in (1, 0.5) and np.isfinite(r.lidar_constant) and r.lidar_constant > 0:
            out[ds] = float(r.lidar_constant)
    RC.cams_water_vapor_profile = _orig_wvprof
    return out

cfgs = {"a_cams": run(L1_ROOT, False), "b_sonde": run(L1_ROOT, True), "c_sonde_offset": run(CORR, True)}
print("night      (a) CAMS-WV   (b) sounding-WV   (c) sounding+offset")
allds = sorted(set().union(*[set(v) for v in cfgs.values()]))
for ds in allds:
    row = "  ".join(f"{cfgs[c].get(ds, float('nan')):12.3f}" for c in ("a_cams", "b_sonde", "c_sonde_offset"))
    print(f"{ds}  {row}")
print()
for c, lbl in (("a_cams", "(a) CAMS-WV, original      "),
               ("b_sonde", "(b) sounding-WV, original  "),
               ("c_sonde_offset", "(c) sounding-WV + offset   ")):
    v = np.array(list(cfgs[c].values()))
    v = v[v < 3 * np.median(v)] if v.size else v
    if v.size:
        print(f"{lbl}: n={v.size}  median C_L={np.median(v):.3f}  scatter(robust)={1.4826*np.median(np.abs(v-np.median(v))):.3f} "
              f"({100*1.4826*np.median(np.abs(v-np.median(v)))/np.median(v):.1f}%)  gap to cloud: {100*(np.median(v)/CL_CLOUD-1):+.1f}%")
# paired (b)-(a) on common nights
com = sorted(set(cfgs["a_cams"]) & set(cfgs["b_sonde"]))
if com:
    d = np.array([100 * (cfgs["b_sonde"][x] / cfgs["a_cams"][x] - 1) for x in com])
    print(f"\npaired sounding-vs-CAMS on {len(com)} common nights: median dC_L = {np.median(d):+.1f}%  (predicted -3.6%)")
print("SOUNDING_CAL_DONE")
