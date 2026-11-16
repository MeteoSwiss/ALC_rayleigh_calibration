"""Slide 6 of the Annual Meeting 2026 deck, v2.3 edition ("Water vapour at 910 nm"): the measured emission
spectra of the Payerne ALCs (paper Fig. A1, the existing slide image) with a panel ABOVE it showing the
water-vapour absorption on the same wavelength axis.

The spectra image is kept pixel for pixel (its data are not in the repository); the new panel is drawn in
a matplotlib axes placed exactly over the pixel span of the image's panel (a) axes (904-917 nm, measured on the
image: x = 53.5 .. 582.5 of 897 px), so the two share one wavelength axis. Absorption = two-way absorption over
a 1 km path at 10 g m-3 of water vapour, 1 - exp(-2 sigma n L), from the operational HITRAN/MT-CKD cross-section
table calibration/data/abs_cross_wv_910nm.nc (surface layer, cm^2 per molecule, 0.0083 nm resolution). The
laser lines (centre +- FWHM/2) of the CL31, CL51 and CL61 are shaded as in the deck's text. Above panel (b)
(CHM15k, 1060-1069 nm) only a note: the table does not cover 1064 nm, where the water-vapour absorption is
negligible (< 0.1 %).

Run from the repo root:  python presentation/slide6_wv_absorption_v23.py <path of the spectra png>
Output: presentation/figs_annual2026/fig_wv_absorption_spectra_v23_slide.png
"""
import sys
from pathlib import Path

import netCDF4
import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_wv_absorption_spectra_v23_slide.png"
# axes boxes of the spectra image, pixels (left, right) measured on the 897 x 369 px image
PANEL_A = (53.5, 582.5, 904.0, 917.0)      # x0, x1, lambda at x0, lambda at x1
PANEL_B = (624.0, 888.5, 1060.0, 1069.0)
RHO_V = 10.0        # g m-3 of water vapour (a humid European day)
PATH_M = 1000.0     # one-way path, m
LINES = [("CL31", 909.7, 6.0, "#1f77b4"), ("CL51", 910.0, 3.4, "#2ca02c"), ("CL61", 910.55, 0.188, "#000000")]


def main(png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    img = np.asarray(Image.open(png).convert("RGB"))
    H, W = img.shape[:2]
    ds = netCDF4.Dataset(LUT)
    lam = np.asarray(ds["lambda"][:], float)
    sig_cm2 = np.asarray(ds["abscs_ave"][:], float)[0]           # surface layer, cm^2 per molecule
    ds.close()
    n_wv = RHO_V / 18.015 * 6.02214e23                            # molecules m-3
    tau2 = 2.0 * sig_cm2 * 1e-4 * n_wv * PATH_M                   # two-way optical depth over PATH_M
    absorption_pct = 100.0 * (1.0 - np.exp(-tau2))

    top_h, gap, head = 200, 6, 44                                 # px: panel height, gap, two title rows
    Htot = H + gap + top_h + head
    dpi = 100
    fig = plt.figure(figsize=(W / dpi, Htot / dpi), dpi=dpi)
    ax_img = fig.add_axes([0, 0, 1, H / Htot])
    ax_img.imshow(img, aspect="auto", interpolation="none")
    ax_img.set_axis_off()
    x0, x1, l0, l1 = PANEL_A
    ax = fig.add_axes([x0 / W, (H + gap) / Htot, (x1 - x0) / W, top_h / Htot])
    ax.plot(lam, absorption_pct, color="#d62728", lw=0.9)
    ax.set_xlim(l0, l1)
    ax.set_yscale("log")
    ax.set_ylim(0.05, 100)
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"], fontsize=9)
    ax.set_ylabel("absorption [%]", fontsize=9)
    ax.tick_params(axis="x", labelbottom=False, direction="in")
    ax.grid(alpha=0.3)
    for name, c, fwhm, col in LINES:
        ax.axvspan(c - fwhm / 2, c + fwhm / 2, color=col, alpha=0.10 if name != "CL61" else 0.35, lw=0)
    ax.text(909.7, 0.062, "CL31", ha="center", va="bottom", fontsize=8, color="#1f77b4", fontweight="bold")
    ax.text(910.0, 0.11, "CL51", ha="center", va="bottom", fontsize=8, color="#2ca02c", fontweight="bold")
    ax.text(910.75, 0.062, "CL61", ha="left", va="bottom", fontsize=8, color="#000000", fontweight="bold")
    # line-averaged absorption, written above the panel with the title
    avg = {}
    for name, c, fwhm, col in LINES:
        sel = (lam >= c - fwhm / 2) & (lam <= c + fwhm / 2)
        avg[name] = float(np.mean(absorption_pct[sel]))
    ax.text(0.0, 1.12, "(c) H$_2$O absorption (HITRAN), two-way over 1 km at 10 g m$^{-3}$", transform=ax.transAxes,
            va="bottom", ha="left", fontsize=9, fontweight="bold")
    ax.text(0.0, 1.03, "line average: " + ",  ".join(f"{k} {v:.1f} %" for k, v in avg.items()),
            transform=ax.transAxes, va="bottom", ha="left", fontsize=8.5)
    xb0, xb1, lb0, lb1 = PANEL_B
    bx = fig.add_axes([xb0 / W, (H + gap) / Htot, (xb1 - xb0) / W, top_h / Htot])
    bx.set_xlim(lb0, lb1)
    bx.set_yscale("log")
    bx.set_ylim(0.05, 100)
    bx.set_yticks([0.1, 1, 10, 100])
    bx.set_yticklabels(["0.1", "1", "10", "100"], fontsize=9)
    bx.tick_params(axis="x", labelbottom=False, direction="in")
    bx.axhline(0.1, color="#d62728", lw=1.0, ls="--")
    bx.text(0.5, 0.55, "1064 nm: no significant\nH$_2$O absorption (< 0.1 %)", transform=bx.transAxes,
            ha="center", va="center", fontsize=9)
    bx.grid(alpha=0.3)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=dpi)
    print("wrote", OUT_PNG, "size", W, Htot)
    print("line-averaged two-way absorption over 1 km at 10 g/m3:", {k: round(v, 2) for k, v in avg.items()},
          "peak", round(float(absorption_pct.max()), 1), "%")


if __name__ == "__main__":
    main(sys.argv[1])
