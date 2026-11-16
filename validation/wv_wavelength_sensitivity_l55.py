"""wv_wavelength_sensitivity_l55.py -- report 03 section 3 redone with the correct CL61 line.

The original study (validation/wv_wavelength_sensitivity.py, June 2026) compared 910.55 vs 910.74 nm
at an assumed FWHM of 1.0 nm and concluded "<0.5 %, either centre is fine, keep ~1 nm averaging".
The correct CL61 line is 910.55 nm / FWHM 0.188 nm (sigma 0.08 nm; Vaisala datasheet B211780EN,
DA10 guide M212895EN-E). This script keeps the original setup -- Payerne, CAMS 1 deg monthly
(D:/CAMS; here 0.4 deg by default, --cams1deg for the July grid), the same HITRAN LUT, night of the 15th +-4 h, 500-3000 m and 2-6 km medians -- and adds:
  * the correct configuration 910.55 / 0.188 and its ratio to the retired 910.74 / 1.0;
  * the FWHM ladder at lambda0 = 910.55 (0.188, 0.5, 1.0, 1.5);
  * the +-0.10 nm centre sensitivity at FWHM 0.188 around 910.55;
  * the per-month table Feb-May 2026 for the three CL61 models.
Nothing of the original script is modified. Run: python -m validation.wv_wavelength_sensitivity_l55
Output: stdout tables + doc/reports/figs_wv_resolution/wv_wavelength_sensitivity_l55.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration.water_vapor_correction.water_vapor import (  # noqa: E402
    cams_water_vapor_profile, two_way_wv_transmission, load_abs_cross_section)

LAT, LON, ALT = 46.8137, 6.9425, 491.0
# CAMS: 0.4 deg monthly by default (the owner's rule: never the 1 deg folders for a number labelled
# correct; at Payerne the 1 deg grid point is ~894 m too high, PWV -26 %). `--cams1deg` replicates the
# July setup (D:/CAMS) for comparison only.
CAMS1 = "--cams1deg" in sys.argv
CAMS_DIR = Path("D:/CAMS") if CAMS1 else Path("D:/CAMS_Monthly_04")
LUT = Path(r"C:/Users/hervo/OneDrive/Documents/MATLAB/MDA/monitoring_alc_monthly/abs_cross_647_full_levels_1000.nc")
OUTFIG = REPO / "doc" / "reports" / "figs_wv_resolution" / (
    "wv_wavelength_sensitivity_l55_cams1deg.png" if CAMS1 else "wv_wavelength_sensitivity_l55.png")

RETIRED = (910.74, 1.0)
CORRECT = (910.55, 0.188)
range_alc = np.arange(0.0, 12000.0, 30.0)
alt_grid = ALT + range_alc
band = (range_alc >= 500) & (range_alc <= 3000)
ref = (range_alc >= 2000) & (range_alc <= 6000)


def night_window(ym, day=15):
    night = np.datetime64(f"{ym[:4]}-{ym[4:6]}-{day:02d}")
    return night - np.timedelta64(4, "h"), night + np.timedelta64(4, "h")


def t2_for(h_wv, n_wv, l0, fwhm):
    return two_way_wv_transmission(alt_grid, ALT, h_wv, n_wv, LUT, l0, fwhm)


def absorption(t2, mask):
    """Median two-way absorption 1 - T2 over a range mask."""
    return 1.0 - float(np.nanmedian(t2[mask]))


def main():
    YM = "202603"
    a, b = night_window(YM)
    h_wv, n_wv = cams_water_vapor_profile(CAMS_DIR / f"CAMS_Beta_{YM}.nc", LAT, LON, a, b)
    configs = [
        ("CL61 910.74 / 1.0 (retired model)", *RETIRED),
        ("CL61 910.55 / 1.0 (old 'manufacturer' test)", 910.55, 1.0),
        ("CL61 910.55 / 0.188 (correct)", *CORRECT),
        ("CL61 910.55 / 0.5", 910.55, 0.5),
        ("CL61 910.55 / 1.5", 910.55, 1.5),
        ("CL61 910.74 / 0.188", 910.74, 0.188),
    ]
    t2s = {}
    base = None
    print(f"Payerne {YM} ({CAMS_DIR}, night of the 15th), baseline = retired 910.74/1.0")
    print(f"{'config':45s} {'T2(0.5-3km)':>11s} {'T2(2-6km)':>10s} {'1-T2 (0.5-3km)':>14s} "
          f"{'beta vs retired':>15s}")
    for lab, l0, fw in configs:
        t2 = t2_for(h_wv, n_wv, l0, fw)
        t2s[lab] = t2
        tb, tr = float(np.nanmedian(t2[band])), float(np.nanmedian(t2[ref]))
        if base is None:
            base = (tb, tr, t2)
        print(f"{lab:45s} {tb:11.4f} {tr:10.4f} {100*(1-tb):13.2f}% {100*(base[0]/tb-1):+14.2f}%")
    t_ret, t_cor = t2s[configs[0][0]], t2s[configs[2][0]]
    for name, m in (("0.5-3 km", band), ("2-6 km (Rayleigh window)", ref)):
        r = absorption(t_ret, m) / absorption(t_cor, m)
        print(f"absorption ratio retired/correct over {name}: x{r:.1f}")

    print("\nPer month, 500-3000 m: T2 for the three CL61 models and the beta change retired->correct")
    print(f"{'month':7s} {'910.74/1.0':>11s} {'910.55/1.0':>11s} {'910.55/0.188':>13s} "
          f"{'beta(correct)/beta(retired)-1':>30s} {'abs. ratio':>10s}")
    months = {}
    for ym in ("202602", "202603", "202604", "202605"):
        f = CAMS_DIR / f"CAMS_Beta_{ym}.nc"
        if not f.exists():
            continue
        a, b = night_window(ym)
        hh, nn = cams_water_vapor_profile(f, LAT, LON, a, b)
        v = [float(np.nanmedian(t2_for(hh, nn, l0, fw)[band]))
             for l0, fw in (RETIRED, (910.55, 1.0), CORRECT)]
        months[ym] = v
        print(f"{ym:7s} {v[0]:11.4f} {v[1]:11.4f} {v[2]:13.4f} {100*(v[0]/v[2]-1):+29.2f}% "
              f"{(1-v[0])/(1-v[2]):10.1f}")

    print("\nCentre sensitivity (+-0.10 nm), March, 500-3000 m")
    for l0, fw in ((910.55, 0.188), (910.74, 1.0)):
        v = {dl: float(np.nanmedian(t2_for(h_wv, n_wv, l0 + dl, fw)[band])) for dl in (-0.10, 0.0, 0.10)}
        swing = (max(v.values()) / min(v.values()) - 1) * 100
        print(f"  {l0}/{fw}: T2(-0.1)={v[-0.10]:.4f} T2(0)={v[0.0]:.4f} T2(+0.1)={v[0.10]:.4f} "
              f"-> beta swing {swing:.2f}%")

    # figure (landscape)
    wl, lut_h_m, abscs = load_abs_cross_section(LUT)
    hidx = int(np.argmin(np.abs(lut_h_m - 3000.0)))
    sel = (wl >= 909.0) & (wl <= 912.5)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    ax[0].semilogy(wl[sel], abscs[sel, hidx], color="0.35", lw=0.8, label="H$_2$O cross-section (~3 km)")
    ax[0].set_xlabel("wavelength [nm]")
    ax[0].set_ylabel("H$_2$O cross-section [cm$^2$] (log)")
    ax[0].set_title("(a) H$_2$O lines and the two CL61 laser models")
    axt = ax[0].twinx()
    for (l0, fw), c, lab in ((RETIRED, "tab:blue", "910.74 / 1.0 nm (retired model)"),
                             (CORRECT, "tab:red", "910.55 / 0.188 nm (correct)")):
        s = fw / (2 * np.sqrt(2 * np.log(2)))
        axt.plot(wl[sel], np.exp(-0.5 * ((wl[sel] - l0) / s) ** 2), color=c, lw=1.8, label=lab)
    axt.set_ylabel("normalised laser spectrum")
    axt.legend(loc="upper left", fontsize=8)
    for lab, c in ((configs[0][0], "tab:blue"), (configs[2][0], "tab:red"),
                   (configs[1][0], "tab:gray")):
        ax[1].plot(t2s[lab], range_alc / 1000.0, color=c, label=lab)
    ax[1].axhspan(0.5, 3.0, color="0.9", zorder=0)
    ax[1].set_xlabel("two-way WV transmission T$^2_{wv}$")
    ax[1].set_ylabel("range [km]")
    ax[1].set_ylim(0, 8)
    ax[1].set_title(f"(b) Payerne {YM}: T$^2_{{wv}}$ vs range")
    ax[1].legend(fontsize=8)
    yms = sorted(months)
    x = np.arange(len(yms))
    for j, (lab, c) in enumerate((("910.74 / 1.0 (retired)", "tab:blue"),
                                  ("910.55 / 1.0", "tab:gray"),
                                  ("910.55 / 0.188 (correct)", "tab:red"))):
        ax[2].bar(x + (j - 1) * 0.27, [100 * (1 - months[m][j]) for m in yms], 0.27, color=c, label=lab)
    ax[2].set_xticks(x)
    ax[2].set_xticklabels(yms)
    ax[2].set_ylabel("median two-way absorption 1 - T$^2_{wv}$, 0.5-3 km [%]")
    ax[2].set_title(f"(c) absorption per month ({'CAMS 1°' if CAMS1 else 'CAMS 0.4°'})")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    OUTFIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTFIG, dpi=140)
    print(f"\nfigure -> {OUTFIG}")


if __name__ == "__main__":
    main()
