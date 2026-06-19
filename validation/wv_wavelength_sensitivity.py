"""wv_wavelength_sensitivity.py

Impact of the laser-wavelength configuration on the water-vapor (WV) correction of
the Payerne ALCs. H2O absorption varies sharply across 905-915 nm, so the assumed
emission wavelength lambda0 and spectral width (FWHM) change the two-way WV
transmission T2_wv and hence the calibrated/validated attenuated backscatter.

Compares, for the Payerne CL61 and CL31:
  - manufacturer CL61 wavelength (910.55 nm) vs the Qmini-measured 910.74 nm;
  - CL61 FWHM sensitivity (the measured line is narrow, spectrometer-limited);
  - CL31 spectral breadth (~5-7 nm FWHM) and inter-acquisition wandering
    (peak 909.0-910.1 nm) — propagated as a T2_wv spread.

Metric: median T2_wv over the 500-3000 m comparison band (drives the validation
bias; beta_corr = beta / T2_wv) and over the 2-6 km Rayleigh reference (drives the
lidar constant). Impact on beta = T2_baseline / T2_config - 1.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from calibration.water_vapor_correction.water_vapor import (  # noqa: E402
    cams_water_vapor_profile, two_way_wv_transmission, load_abs_cross_section,
)

LAT, LON, ALT = 46.8137, 6.9425, 491.0
CAMS_DIR = Path("D:/CAMS")
LUT = Path(r"C:/Users/hervo/OneDrive/Documents/MATLAB/MDA/monitoring_alc_monthly/abs_cross_647_full_levels_1000.nc")
OUTFIG = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/wv_wavelength_sensitivity.png")

range_alc = np.arange(0.0, 12000.0, 30.0)
alt_grid = ALT + range_alc
band = (range_alc >= 500) & (range_alc <= 3000)   # comparison band (validation bias)
ref = (range_alc >= 2000) & (range_alc <= 6000)   # Rayleigh fit window (constant)


def night_window(ym, day=15):
    night = np.datetime64(f"{ym[:4]}-{ym[4:6]}-{day:02d}")
    return night - np.timedelta64(4, "h"), night + np.timedelta64(4, "h")


def t2_for(h_wv, n_wv, l0, fwhm):
    return two_way_wv_transmission(alt_grid, ALT, h_wv, n_wv, LUT, l0, fwhm)


# ---- Main month: March 2026 ----
YM = "202603"
ts, te = night_window(YM)
h_wv, n_wv = cams_water_vapor_profile(CAMS_DIR / f"CAMS_Beta_{YM}.nc", LAT, LON, ts, te)

configs = [
    ("CL61 manuf. 910.55 (FWHM 1.0)", 910.55, 1.0, "CL61"),
    ("CL61 meas.  910.74 (FWHM 1.0)", 910.74, 1.0, "CL61"),
    ("CL61 910.74 FWHM 0.1", 910.74, 0.1, "CL61f"),
    ("CL61 910.74 FWHM 0.5", 910.74, 0.5, "CL61f"),
    ("CL61 910.74 FWHM 1.5", 910.74, 1.5, "CL61f"),
    ("CL31 meas.  909.7 (FWHM 6.0)", 909.7, 6.0, "CL31"),
    ("CL31 909.0 FWHM 6.0", 909.0, 6.0, "CL31w"),
    ("CL31 910.1 FWHM 6.0", 910.1, 6.0, "CL31w"),
    ("CL31 909.7 FWHM 5.0", 909.7, 5.0, "CL31w"),
    ("CL31 909.7 FWHM 7.0", 909.7, 7.0, "CL31w"),
]
res = {}
for label, l0, fw, grp in configs:
    t2 = t2_for(h_wv, n_wv, l0, fw)
    res[label] = dict(l0=l0, fw=fw, grp=grp, t2=t2,
                      band=float(np.nanmedian(t2[band])),
                      ref=float(np.nanmedian(t2[ref])))

base = res["CL61 meas.  910.74 (FWHM 1.0)"]["band"]
print(f"Payerne {YM}: median T2_wv over 500-3000 m (baseline = CL61 measured 910.74)")
print(f"{'config':32s} {'T2(500-3000m)':>13s} {'T2(2-6km)':>10s} {'beta impact':>12s}")
for label in res:
    r = res[label]
    impact = base / r["band"] - 1.0   # beta change vs baseline (beta ~ 1/T2)
    print(f"{label:32s} {r['band']:13.4f} {r['ref']:10.4f} {impact*100:+11.2f}%")

# ---- CL61 manufacturer vs measured, per month Feb-May ----
print("\nCL61 manufacturer (910.55) vs measured (910.74), FWHM 1.0 - per month")
print(f"{'month':7s} {'T2 manuf':>9s} {'T2 meas':>9s} {'beta impact (meas vs manuf)':>27s}")
for ym in ("202602", "202603", "202604", "202605"):
    f = CAMS_DIR / f"CAMS_Beta_{ym}.nc"
    if not f.exists():
        continue
    a, b = night_window(ym)
    hh, nn = cams_water_vapor_profile(f, LAT, LON, a, b)
    t2_man = float(np.nanmedian(t2_for(hh, nn, 910.55, 1.0)[band]))
    t2_mea = float(np.nanmedian(t2_for(hh, nn, 910.74, 1.0)[band]))
    impact = t2_man / t2_mea - 1.0   # beta(meas)/beta(manuf) - 1
    print(f"{ym:7s} {t2_man:9.4f} {t2_mea:9.4f} {impact*100:+26.2f}%")

# ---- CL31 spread (wandering + FWHM) over the March profile ----
cl31 = [res[k]["band"] for k in res if res[k]["grp"] in ("CL31", "CL31w")]
print(f"\nCL31 T2(500-3000m) spread over lambda0 909.0-910.1 nm x FWHM 5-7 nm: "
      f"{min(cl31):.4f} - {max(cl31):.4f}  (range {(max(cl31)/min(cl31)-1)*100:.1f}% in beta)")

# ---- CL61: does a narrow FWHM amplify the +-0.10 nm center uncertainty? ----
print("\nCL61 center sensitivity to the +-0.10 nm measurement uncertainty, March, 500-3000 m")
for fw in (0.1, 1.0):
    v = {dl: float(np.nanmedian(t2_for(h_wv, n_wv, 910.74 + dl, fw)[band]))
         for dl in (-0.10, 0.0, +0.10)}
    swing = (max(v.values()) / min(v.values()) - 1) * 100
    print(f"  FWHM {fw:3.1f}: T2(-0.1)={v[-0.10]:.4f}  T2(0)={v[0.0]:.4f}  "
          f"T2(+0.1)={v[0.10]:.4f}   -> beta swing {swing:.2f}%")

# =====================================================================
# Figure (landscape)
# =====================================================================
wl, lut_h_m, abscs = load_abs_cross_section(LUT)
hidx = int(np.argmin(np.abs(lut_h_m - 3000.0)))   # ~3 km cross-section
sigma_3km = abscs[:, hidx]
sel = (wl >= 905) & (wl <= 915)

fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))

# (a) absorption spectrum + laser spectra
axa = ax[0]
axa.semilogy(wl[sel], sigma_3km[sel], color="0.35", lw=0.8, label="H$_2$O abs. x-section (~3 km)")
axa.set_xlim(905, 915)
axa.set_xlabel("wavelength [nm]")
axa.set_ylabel("H$_2$O cross-section [cm$^2$]  (log)")
axa.set_title("(a) H$_2$O absorption band + laser spectra")
axt = axa.twinx()
lasers = [("CL61 manuf. 910.55", 910.55, 1.0, "tab:red", "--"),
          ("CL61 meas. 910.74", 910.74, 1.0, "tab:blue", "-"),
          ("CL31 909.7 (FWHM 6)", 909.7, 6.0, "tab:green", "-")]
for lab, l0, fw, c, ls in lasers:
    s = fw / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((wl[sel] - l0) / s) ** 2)
    axt.plot(wl[sel], g, color=c, ls=ls, lw=1.8, label=lab)
axt.axvline(910.55, color="tab:red", lw=0.6, alpha=0.5)
axt.axvline(910.74, color="tab:blue", lw=0.6, alpha=0.5)
axt.set_ylabel("laser spectrum (normalised)")
axt.set_ylim(0, 1.25)
h1, l1 = axa.get_legend_handles_labels()
h2, l2 = axt.get_legend_handles_labels()
axt.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")

# (b) T2_wv profiles
axb = ax[1]
prof = [("CL61 manuf. 910.55 (FWHM 1.0)", "tab:red", "--"),
        ("CL61 meas.  910.74 (FWHM 1.0)", "tab:blue", "-"),
        ("CL31 meas.  909.7 (FWHM 6.0)", "tab:green", "-")]
for key, c, ls in prof:
    axb.plot(res[key]["t2"], range_alc, color=c, ls=ls, lw=1.8,
             label=key.replace("  ", " "))
# CL61 FWHM band (0.1 vs 1.5) as a shaded envelope
t2lo = res["CL61 910.74 FWHM 0.1"]["t2"]
t2hi = res["CL61 910.74 FWHM 1.5"]["t2"]
axb.fill_betweenx(range_alc, t2lo, t2hi, color="tab:blue", alpha=0.12,
                  label="CL61 FWHM 0.1-1.5 nm")
axb.axhspan(500, 3000, color="0.85", alpha=0.4, zorder=0)
axb.set_xlabel("two-way WV transmission T$^2_{wv}$")
axb.set_ylabel("range AGL [m]")
axb.set_ylim(0, 6000)
axb.set_title(f"(b) T$^2_{{wv}}$ profile (Payerne {YM})")
axb.legend(fontsize=7, loc="lower left")
axb.grid(alpha=0.3)

# (c) median T2 over the comparison band per config
axc = ax[2]
labels = [k for k in res]
vals = [res[k]["band"] for k in labels]
cols = {"CL61": "tab:blue", "CL61f": "tab:cyan", "CL31": "tab:green", "CL31w": "tab:olive"}
bcol = [cols[res[k]["grp"]] for k in labels]
short = [k.replace("CL61 ", "CL61 ").replace("CL31 ", "CL31 ") for k in labels]
y = np.arange(len(labels))
axc.barh(y, vals, color=bcol)
axc.set_yticks(y)
axc.set_yticklabels(short, fontsize=7)
axc.invert_yaxis()
axc.set_xlabel("median T$^2_{wv}$ over 500-3000 m")
axc.set_xlim(min(vals) * 0.97, max(vals) * 1.01)
axc.set_title("(c) WV transmission by laser config")
axc.axvline(base, color="tab:blue", lw=1.0, ls=":", label="CL61 measured")
for yi, v in zip(y, vals):
    axc.text(v, yi, f" {v:.3f}", va="center", fontsize=6.5)
axc.grid(axis="x", alpha=0.3)

fig.suptitle("Payerne ALC — water-vapor correction vs laser-wavelength configuration",
             fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])
OUTFIG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTFIG, dpi=200)
print(f"\nSaved figure: {OUTFIG}")
print("WV_WL_DONE")
