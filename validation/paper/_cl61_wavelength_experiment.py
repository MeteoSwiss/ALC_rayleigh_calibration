"""CL61 Rayleigh-calibration wavelength experiment (Payerne):
recalibrate every successful night with the MEASURED laser line (910.74 nm, Qmini 2026-06-02)
vs the MANUFACTURER nominal (910.55 nm), same FWHM 1.0 nm. The cloud method of the same
instrument gives C_L = 1.425 — whichever line brings the Rayleigh C_L closer wins.
Also: pure two-way WV transmission at the fit-window altitudes for both lines -> the expected
C_L scaling, to check the recalibration result is understood physically."""
import sys, csv, tempfile, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np

from calibration.config import CalibrationOptions, InstrumentInfo, InstrumentType, DataLevel
from calibration.rayleigh.calibration import calibrate_rayleigh
from calibration.water_vapor_correction import water_vapor as WVMOD
from calibration.io.cams import ensure_cams_file
from validation.paper.calib_benchmark import RAYLEIGH_METHOD, L1_ROOT, CAMS, REPO

WMO, IDENT = "0-20000-0-06610", "C"
LAT, LON, ALT = 46.813, 6.943, 491.0
CAL_CSV = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/fullcal_l1_2026") / f"{WMO}_{IDENT}" / f"{WMO}_{IDENT}_cal.csv"

# successful nights + their fit windows from the operational calout
nights = []
for r in csv.DictReader(open(CAL_CSV, encoding="utf-8")):
    if r["method"] == "rayleigh" and str(r["flag"]) in ("1", "1.0", "0.5"):
        nights.append((r["date"], float(r["cal_value"]),
                       float(r["bottom_height"] or 0), float(r["top_height"] or 0)))
print(f"{len(nights)} successful nights; operational C_L values: "
      f"median={np.median([n[1] for n in nights]):.3f}")

def run_all(lam0, fwhm, tag):
    WVMOD.LASER_SPECTRUM["CL61"] = (lam0, fwhm)
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.molecular_method = RAYLEIGH_METHOD
    o.folder_root = L1_ROOT; o.data_level = DataLevel.L1
    o.cams_folder = Path(CAMS); o.abs_cs_lookup_table = Path("")
    o.apply_wv_correction = True
    o.folder_output = Path(tempfile.mkdtemp()); o.plot_main = o.plot_all = False
    info = InstrumentInfo(site_name=WMO, wmo_id=WMO, identifier=IDENT,
                          instrument_type=InstrumentType.CL61, latitude=LAT, longitude=LON, altitude=ALT)
    out = {}
    for ds, cop, bh, th in nights:
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception as e:
            print(f"  {ds} {tag}: EXC {e}"); continue
        if r.flag in (1, 0.5) and np.isfinite(r.lidar_constant) and r.lidar_constant > 0:
            out[ds] = float(r.lidar_constant)
    return out

res_meas = run_all(910.74, 1.0, "measured")
res_manu = run_all(910.55, 1.0, "manufacturer")
WVMOD.LASER_SPECTRUM["CL61"] = (910.74, 1.0)   # restore

common = sorted(set(res_meas) & set(res_manu))
cm = np.array([res_meas[d] for d in common]); cu = np.array([res_manu[d] for d in common])
# drop the far outliers (>3x median) for the robust summary, report both
ok = (cm < 3 * np.median(cm)) & (cu < 3 * np.median(cu))
print(f"\ncommon nights recalibrated: {len(common)} ({ok.sum()} after outlier drop)")
print(f"C_L measured   (910.74/1.0): median={np.median(cm[ok]):.3f}")
print(f"C_L manufactur (910.55/1.0): median={np.median(cu[ok]):.3f}")
print(f"per-night ratio manu/meas: median={np.median(cu[ok]/cm[ok]):.3f}")
print(f"cloud-method reference C_L = 1.425;  gap measured={100*(np.median(cm[ok])/1.425-1):+.1f}%  "
      f"manufacturer={100*(np.median(cu[ok])/1.425-1):+.1f}%")

# pure-transmission check at the fit windows (one CAMS month, Payerne)
cams = ensure_cams_file(Path(CAMS), "20260401", auto_download=False)
from calibration.water_vapor_correction.water_vapor import cams_water_vapor_profile, two_way_wv_transmission
prof = cams_water_vapor_profile(cams, LAT, LON, np.datetime64("2026-04-01"), np.datetime64("2026-05-01"))
h, n = prof
zgrid = np.arange(0, 8001, 30.0) + ALT
t2_meas = two_way_wv_transmission(zgrid, ALT, h, n, WVMOD.DEFAULT_ABS_CROSS_SECTION, 910.74, 1.0)
t2_manu = two_way_wv_transmission(zgrid, ALT, h, n, WVMOD.DEFAULT_ABS_CROSS_SECTION, 910.55, 1.0)
for z in (2500, 3500, 4500):
    j = np.argmin(np.abs(zgrid - ALT - z))
    print(f"T2_wv @{z} m AGL: measured={t2_meas[j]:.3f} manufacturer={t2_manu[j]:.3f} "
          f"-> C_L ratio manu/meas = {t2_meas[j]/t2_manu[j]:.3f}")
print("WAVELENGTH_EXPERIMENT_DONE")
