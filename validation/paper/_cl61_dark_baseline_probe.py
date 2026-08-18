"""CL61 weak-signal / baseline probe (task: confirm or kill the 'weak-signal artefact' hypothesis).

On the successful Rayleigh-calibration nights, compare the NIGHTLY-MEAN CL61 attenuated
backscatter against the attenuated molecular model beta_mol(910)*T2 over 2-14 km:
 - 8-14 km the atmosphere is essentially molecular-only and tiny -> any systematic departure of
   the nightly mean from the molecular curve IS a baseline/afterpulse residual (random noise
   averages down by sqrt(N) ~ 25 over a night).
 - A baseline of ~ -0.016 Mm^-1 sr^-1 at the 2.6-6 km fit windows would explain the -12 % C_L.
Uses the same WV transmission as the calibration (910.74/1.0, monthly CAMS).
"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from calibration.water_vapor_correction.water_vapor import (
    cams_water_vapor_profile, two_way_wv_transmission, DEFAULT_ABS_CROSS_SECTION)
from calibration.io.cams import ensure_cams_file

WMO, IDENT, ALT, LAT, LON = "0-20000-0-06610", "C", 491.0, 46.813, 6.943
L1 = Path("D:/E-PROFILE_L1_2026")
STD = Path("calibration/data/standard_atmosphere_US_1976_50km.csv")
NIGHTS = ["20260316", "20260328", "20260402", "20260407", "20260410", "20260411",
          "20260417", "20260420", "20260422", "20260423", "20260424", "20260428"]
CL_CLOUD = 1.4252   # the cloud-method (CHM15k-anchored) constant
CL_RAY = 1.251      # the Rayleigh-method constant

rows, mean_profiles = [], []
rng_ref = None
for ds in NIGHTS:
    f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
    if not f.exists():
        continue
    with Dataset(f) as nc:
        t = np.asarray(nc.variables["time"][:], "f8")
        tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
        rng = np.asarray(nc.variables["range"][:], "f8")
        rcs = np.asarray(nc.variables["rcs_0"][:], "f8")
        if rcs.shape[0] != t.size:
            rcs = rcs.T
    # night window 21-04 UT (index by fraction of day)
    frac = (t - np.floor(t)) if "day" in tu else (t % 86400) / 86400.0
    ni = (frac >= 21 / 24) | (frac <= 4 / 24)
    if ni.sum() < 100:
        continue
    prof = np.nanmean(rcs[ni], axis=0) * 1e6          # -> Mm^-1 sr^-1
    n = int(ni.sum())
    if rng_ref is None:
        rng_ref = rng
    mean_profiles.append((ds, prof, n))

print(f"nights loaded: {len(mean_profiles)}")

# attenuated molecular model at 910 nm (+ WV transmission from the same month's CAMS)
atm = load_standard_atmosphere(STD, rng_ref + ALT)
mol = calculate_molecular_properties(atm.temperature, atm.pressure, rng_ref, 910.0e-9)
bmol_att = mol.beta_mol * mol.transmission * 1e6      # Mm^-1 sr^-1 (two-way molecular)
cams = ensure_cams_file(Path("D:/CAMS"), "20260401", auto_download=False)
h, nwv = cams_water_vapor_profile(cams, LAT, LON, np.datetime64("2026-04-01"), np.datetime64("2026-05-01"))
t2wv = two_way_wv_transmission(rng_ref + ALT, ALT, h, nwv, DEFAULT_ABS_CROSS_SECTION, 910.74, 1.0)
model = bmol_att * t2wv                                # what C_L * model should match

# per-night residual baseline in altitude bands (measured_mean - CL_cloud*model)
print("\nresidual baseline b = mean(beta_meas) - C_L(cloud)*beta_mol*T2  [Mm^-1 sr^-1]")
print("night      N     3-6 km      8-11 km    11-14 km   (need ~ -0.016*C_L at 3-6 km, i.e. -0.023, to explain -12%)")
b36_all, b811_all, b1114_all = [], [], []
for ds, prof, n in mean_profiles:
    res = prof - CL_CLOUD * model
    def bandmean(lo, hi):
        m = (rng_ref >= lo) & (rng_ref <= hi)
        return float(np.nanmean(res[m]))
    b36, b811, b1114 = bandmean(3000, 6000), bandmean(8000, 11000), bandmean(11000, 14000)
    b36_all.append(b36); b811_all.append(b811); b1114_all.append(b1114)
    print(f"{ds}  {n:4d}  {b36:+.4f}    {b811:+.4f}    {b1114:+.4f}")
print(f"\nMEDIAN     -    {np.median(b36_all):+.4f}    {np.median(b811_all):+.4f}    {np.median(b1114_all):+.4f}")
print(f"required baseline for the -12% C_L gap at the fit windows: {-0.12*CL_CLOUD*np.nanmean(model[(rng_ref>=3000)&(rng_ref<=6000)]):+.4f}")

# figure: nightly-mean profiles vs the two calibrated molecular models
fig, ax = plt.subplots(1, 2, figsize=(15, 6))
for ds, prof, n in mean_profiles:
    ax[0].plot(prof, rng_ref, lw=0.8, alpha=0.5)
ax[0].plot(CL_CLOUD * model, rng_ref, "k-", lw=2.2, label=f"C_L(cloud)={CL_CLOUD:.2f} x mol x T2")
ax[0].plot(CL_RAY * model, rng_ref, "b--", lw=2.2, label=f"C_L(Rayleigh)={CL_RAY:.2f} x mol x T2")
ax[0].set_xlim(-0.05, 0.6); ax[0].set_ylim(0, 14000); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=9)
ax[0].set_xlabel(r"nightly-mean $\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]"); ax[0].set_ylabel("range AGL [m]")
ax[0].set_title("(a) CL61 nightly means (21-04 UT) vs calibrated molecular models")
med = np.nanmedian(np.array([p for _, p, _ in mean_profiles]), axis=0)
ax[1].plot(med - CL_CLOUD * model, rng_ref, "k-", lw=1.8, label="residual vs C_L(cloud)")
ax[1].plot(med - CL_RAY * model, rng_ref, "b--", lw=1.8, label="residual vs C_L(Rayleigh)")
ax[1].axvline(0, color="0.5", lw=0.8)
ax[1].set_xlim(-0.08, 0.08); ax[1].set_ylim(0, 14000); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=9)
ax[1].set_xlabel(r"median residual [Mm$^{-1}$sr$^{-1}$]")
ax[1].set_title("(b) median nightly-mean minus model: a flat offset = baseline artefact")
fig.suptitle("CL61 Payerne baseline probe on the Rayleigh-calibration nights", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_dark_baseline.png")
fig.savefig(out, dpi=150)
print("saved", out)
print("DARK_PROBE_DONE")
