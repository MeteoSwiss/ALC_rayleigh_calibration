"""Generate landscape figures for the calibration status presentation.

All numbers are taken from the project validation reports:
  - doc/reports/l1_2026_variability_report.md
  - doc/reports/network_v2_vs_v11_report.md
  - validation/compare_daily.py (sign-error fingerprint)
Figures are written to presentation/figs/.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# ---------------------------------------------------------------- palette
INK     = "#0E1B33"   # darkest navy
NAVY    = "#13294B"   # primary dark
SKY     = "#2E7DA1"   # supporting teal-blue
GREEN   = "#2FBF71"   # 532 nm / Mini-MPL / "best" accent
AMBER   = "#E9A23B"   # highlight / caution
RED     = "#D64550"   # error / alert
LIGHT   = "#F5F7FB"   # light background
GRID    = "#D9E0EA"
MUTED   = "#5A6B82"
INKTXT  = "#1A2A44"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INKTXT,
    "text.color": INKTXT,
    "xtick.color": INKTXT,
    "ytick.color": INKTXT,
    "axes.linewidth": 0.8,
    "figure.dpi": 200,
    "savefig.dpi": 200,
})

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def save(fig, name):
    path = os.path.join(FIGS, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white", pad_inches=0.12)
    plt.close(fig)
    print("wrote", path)


# ====================================================================== 1
# Method precision ranking (sigma_SD by molecular-window method)
# mean over 24 instruments, L1 2026  (l1_2026_variability_report.md table 3)
def fig_method_ranking():
    methods = [
        ("E-PROF v2", 9.5, GREEN, "recommended"),
        ("E-PROF v1.0 (sign error)", 10.3, RED, "current operational"),
        ("EARLINET / SCC", 10.5, SKY, ""),
        ("E-PROF v0.25 (MATLAB)", 11.4, NAVY, ""),
        ("E-PROF v1.1 (sign corr.)", 11.6, NAVY, ""),
        ("E-PROF v1.2 (improved)", 12.9, NAVY, ""),
        ("Bellini / ALICENET", 14.1, "#9AA7B8", ""),
    ]
    methods = methods[::-1]  # best at top
    labels = [m[0] for m in methods]
    vals   = [m[1] for m in methods]
    cols   = [m[2] for m in methods]
    tags   = [m[3] for m in methods]
    y = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    bars = ax.barh(y, vals, color=cols, height=0.62, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=12.5)
    ax.set_xlim(0, 18.0)
    ax.set_xlabel(r"$\sigma_{SD}$  —  night-to-night precision (% of median $C_L$)   ·   lower is better",
                  fontsize=12.5)
    style(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    for yi, v, tag in zip(y, vals, tags):
        ax.text(v + 0.25, yi, f"{v:.1f}", va="center", ha="left",
                fontsize=12, fontweight="bold", color=INKTXT)
        if tag:
            ax.text(v + 1.05, yi, f"◀ {tag}", va="center", ha="left",
                    fontsize=10.5, fontstyle="italic",
                    color=GREEN if tag == "recommended" else RED, fontweight="bold")
    ax.text(0, 1.135, "Night-to-night precision by molecular-window method",
            transform=ax.transAxes, fontsize=15, fontweight="bold", color=INKTXT)
    ax.text(0, 1.05, "E-PROFILE L1 2026  ·  mean of 24 instruments (10 CHM15k, 4 Mini-MPL, 10 CL61)",
            transform=ax.transAxes, fontsize=11, color=MUTED)
    save(fig, "fig_method_ranking.png")


# ====================================================================== 2
# E-PROF v2 vs v1.1 across the network (L2 2026)
# network_v2_vs_v11_report.md
def fig_v2_vs_v11():
    types = ["CHM15k\n(1064 nm)", "CL61\n(910 nm)", "Mini-MPL\n(532 nm)"]
    yield_v2  = [87.5, 75.0, 100.0]
    yield_v11 = [73.9, 69.8, 100.0]
    sd_v2     = [9.5, 7.0, 6.4]
    sd_v11    = [11.7, 9.7, 6.8]

    x = np.arange(len(types))
    w = 0.38
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 4.9))

    # --- yield
    axL.bar(x - w/2, yield_v11, w, label="E-PROF v1.1", color=NAVY)
    axL.bar(x + w/2, yield_v2,  w, label="E-PROF v2",   color=GREEN)
    axL.set_ylim(0, 116)
    axL.set_ylabel("valid-on-clear yield (%)", fontsize=12)
    axL.set_title("Yield  —  higher is better", fontsize=13, fontweight="bold", loc="left")
    for xi, (a, b) in enumerate(zip(yield_v11, yield_v2)):
        axL.text(xi - w/2, a + 1.5, f"{a:.0f}", ha="center", fontsize=10.5, color=NAVY, fontweight="bold")
        axL.text(xi + w/2, b + 1.5, f"{b:.0f}", ha="center", fontsize=10.5, color="#1E7A48", fontweight="bold")
    axL.set_xticks(x); axL.set_xticklabels(types, fontsize=11)
    style(axL); axL.yaxis.grid(True, color=GRID, linewidth=0.8)

    # --- sigma_SD
    axR.bar(x - w/2, sd_v11, w, label="E-PROF v1.1", color=NAVY)
    axR.bar(x + w/2, sd_v2,  w, label="E-PROF v2",   color=GREEN)
    axR.set_ylim(0, 14)
    axR.set_ylabel(r"$\sigma_{SD}$ (%)", fontsize=12)
    axR.set_title(r"Precision $\sigma_{SD}$  —  lower is better", fontsize=13, fontweight="bold", loc="left")
    for xi, (a, b) in enumerate(zip(sd_v11, sd_v2)):
        axR.text(xi - w/2, a + 0.15, f"{a:.1f}", ha="center", fontsize=10.5, color=NAVY, fontweight="bold")
        axR.text(xi + w/2, b + 0.15, f"{b:.1f}", ha="center", fontsize=10.5, color="#1E7A48", fontweight="bold")
    axR.set_xticks(x); axR.set_xticklabels(types, fontsize=11)
    style(axR); axR.yaxis.grid(True, color=GRID, linewidth=0.8)
    axR.legend(frameon=False, fontsize=11, loc="upper right", ncol=1)

    fig.suptitle("E-PROF v2 vs v1.1 across the E-PROFILE network (L2 2026)",
                 fontsize=15, fontweight="bold", color=INKTXT, x=0.5, y=1.02)
    fig.text(0.5, -0.02, "L1 shows the same pattern.  v2 lifts CHM15k yield by ~+10 pp and cuts "
             r"$\sigma_{SD}$ by ~2 pp; Mini-MPL is already precise.",
             ha="center", fontsize=10.5, color=MUTED)
    save(fig, "fig_v2_vs_v11.png")


# ====================================================================== 3
# Mini-MPL per-unit precision (l1_2026_variability_report.md, table 4)
def fig_minimpl_precision():
    units = ["Mini-MPL-Mobile", "Aléria", "Lille", "Brest"]
    sd    = [6.8, 6.8, 10.2, 20.5]
    valid = [71, 81, 59, 39]
    cols  = [GREEN, GREEN, GREEN, RED]
    x = np.arange(len(units))

    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    bars = ax.bar(x, sd, width=0.6, color=cols, edgecolor="white", linewidth=0.8)
    ax.axhline(8.5, color=NAVY, ls="--", lw=1.6)
    ax.text(3.46, 8.5 + 0.4, "group median 8.5%", ha="right", va="bottom",
            fontsize=11, color=NAVY, fontstyle="italic", fontweight="bold")
    ax.axhspan(0, 5, color=SKY, alpha=0.10)
    ax.text(0.02, 4.35, "~5% irreducible floor", fontsize=10.5, color=NAVY, fontstyle="italic", fontweight="bold")

    for xi, (v, vd) in enumerate(zip(sd, valid)):
        ax.text(xi, v + 0.35, f"{v:.1f}%", ha="center", fontsize=13, fontweight="bold", color=INKTXT)
        ax.text(xi, 0.6, f"{vd}% valid", ha="center", fontsize=10.5, color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(units, fontsize=12.5)
    ax.set_ylim(0, 23)
    ax.set_ylabel(r"$\sigma_{SD}$  (% of median $C_L$)", fontsize=12.5)
    style(ax); ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.text(0, 1.115, "Mini-MPL Rayleigh calibration precision per unit  (E-PROF v2, L1 2026)",
            transform=ax.transAxes, fontsize=14, fontweight="bold", color=INKTXT)
    ax.text(0, 1.04, "Three of four units reach CHM15k-class precision; Brest is background-contaminated.",
            transform=ax.transAxes, fontsize=11, color=MUTED)
    save(fig, "fig_minimpl_precision.png")


# ====================================================================== 4
# Sign-fix impact: wavelength dependence + stability unchanged
def fig_signfix_impact():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 4.9),
                                   gridspec_kw={"width_ratios": [1.25, 1.0]})

    # --- left: |C_L| change when sign corrected, by wavelength
    labels = ["CHM15k\n1064 nm", "CL61\n910 nm", "Mini-MPL\n532 nm"]
    centre = [0.5, 0.5, 13.0]
    cols   = [NAVY, SKY, GREEN]
    x = np.arange(3)
    axL.bar(x, centre, width=0.6, color=cols, edgecolor="white", linewidth=0.8)
    # range whisker for the 532 nm bar (8-18%)
    axL.plot([2, 2], [8, 18], color=INK, lw=2.2)
    axL.plot([1.93, 2.07], [8, 8], color=INK, lw=2.2)
    axL.plot([1.93, 2.07], [18, 18], color=INK, lw=2.2)
    axL.text(0, 1.1, "~0.5%", ha="center", fontsize=11.5, fontweight="bold", color=INKTXT)
    axL.text(1, 1.1, "~0.5%", ha="center", fontsize=11.5, fontweight="bold", color=INKTXT)
    axL.text(2, 18.7, "+8 to +18%", ha="center", fontsize=12.5, fontweight="bold", color="#1E7A48")
    axL.set_xticks(x); axL.set_xticklabels(labels, fontsize=11.5)
    axL.set_ylim(0, 21)
    axL.set_ylabel(r"absolute change in $C_L$ (%)", fontsize=12)
    style(axL); axL.yaxis.grid(True, color=GRID, linewidth=0.8)
    axL.set_title("Absolute $C_L$ shift when the sign is corrected",
                  fontsize=13, fontweight="bold", loc="left")
    axL.text(0.0, -0.30, "Effect scales with aerosol optical depth → negligible in the\nclear-air "
             "windows at 1064/910 nm, large at 532 nm.",
             transform=axL.transAxes, fontsize=10, color=MUTED)

    # --- right: stability unchanged
    bars_x = [0, 1]
    sd = [10.3, 11.6]
    axR.bar(bars_x, sd, width=0.55, color=[AMBER, GREEN], edgecolor="white", linewidth=0.8)
    for xi, v in zip(bars_x, sd):
        axR.text(xi, v + 0.2, f"{v:.1f}%", ha="center", fontsize=13, fontweight="bold", color=INKTXT)
    axR.set_xticks(bars_x)
    axR.set_xticklabels(["v1.0\n(sign error)", "v1.1\n(corrected)"], fontsize=11.5)
    axR.set_ylim(0, 14)
    axR.set_ylabel(r"$\sigma_{SD}$ of $C_L$ (%)", fontsize=12)
    style(axR); axR.yaxis.grid(True, color=GRID, linewidth=0.8)
    axR.set_title(r"$C_L$ night-to-night stability", fontsize=13, fontweight="bold", loc="left")
    axR.annotate("≈ unchanged", xy=(0.5, 12.2), ha="center", fontsize=12,
                 color=NAVY, fontstyle="italic", fontweight="bold")
    axR.text(0.0, -0.30, "The sign error corrupts the downstream $\\beta_{att}$ /\nAOD product, "
             "not the calibration constant $C_L$.",
             transform=axR.transAxes, fontsize=10, color=MUTED)

    fig.suptitle("Impact of the Klett sign correction (E-PROF v1.0 → v1.1)",
                 fontsize=15, fontweight="bold", color=INKTXT, x=0.5, y=1.03)
    save(fig, "fig_signfix_impact.png")


# ====================================================================== 5
# Rayleigh / molecular-window method schematic (illustrative)
def fig_rayleigh_schematic():
    rng = np.random.default_rng(7)
    z = np.linspace(0.05, 8.0, 400)               # km
    H = 8.0
    beta_mol = np.exp(-z / H)                      # molecular ~ density
    aerosol = 1.6 * np.exp(-((z - 0.9) / 0.7) ** 2) + 0.25 * np.exp(-z / 1.2)
    signal = beta_mol * (1 + aerosol)             # range-corrected signal (a.u.)
    signal_noisy = signal * (1 + 0.04 * rng.standard_normal(z.size))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 5.0),
                                   gridspec_kw={"width_ratios": [1.05, 1.0]})

    # left: profiles, altitude on y
    axL.plot(signal_noisy, z, color=SKY, lw=1.8, label="measured signal  $RCS/r^2$")
    axL.plot(beta_mol * signal[np.argmin(np.abs(z - 3.5))] / beta_mol[np.argmin(np.abs(z - 3.5))],
             z, color=NAVY, lw=2.0, ls="--", label=r"molecular reference  $\beta_{mol}$")
    z0, z1 = 2.5, 5.0
    axL.axhspan(z0, z1, color=GREEN, alpha=0.16)
    axL.text(0.96, (z0 + z1) / 2, "molecular\nwindow\n(2–6 km)", ha="right", va="center",
             fontsize=11, color="#1E7A48", fontweight="bold")
    axL.annotate("aerosol layer", xy=(signal_noisy[np.argmin(np.abs(z-0.9))], 0.9),
                 xytext=(1.1, 1.7), fontsize=10.5, color=AMBER, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.5))
    axL.set_xlim(0, 1.25); axL.set_ylim(0, 8)
    axL.set_xlabel("range-corrected signal (a.u.)", fontsize=12)
    axL.set_ylabel("altitude (km AGL)", fontsize=12)
    style(axL); axL.grid(True, color=GRID, linewidth=0.7)
    axL.legend(frameon=False, fontsize=10.5, loc="upper right")
    axL.set_title("1 · Find a clean molecular window", fontsize=13, fontweight="bold", loc="left")

    # right: linear fit signal vs molecular in the window
    mask = (z >= z0) & (z <= z1)
    xm = beta_mol[mask]
    ym = signal_noisy[mask]
    A = np.polyfit(xm, ym, 1)
    xx = np.linspace(xm.min(), xm.max(), 50)
    axR.scatter(xm, ym, s=22, color=GREEN, edgecolor="white", linewidth=0.5, zorder=3,
                label="window points")
    axR.plot(xx, np.polyval(A, xx), color=NAVY, lw=2.4, zorder=2,
             label=r"fit:  signal $= C_L\,\beta_{mol}$")
    axR.set_xlabel(r"molecular backscatter  $\beta_{mol}$  (a.u.)", fontsize=12)
    axR.set_ylabel("signal in window (a.u.)", fontsize=12)
    style(axR); axR.grid(True, color=GRID, linewidth=0.7)
    axR.legend(frameon=False, fontsize=10.5, loc="upper left")
    axR.set_title(r"2 · Slope $\Rightarrow$ lidar constant $C_L$", fontsize=13, fontweight="bold", loc="left")
    axR.text(0.98, 0.06, "+ Klett inversion below the window\n+ sensitivity test (LR ±20 sr, alt ±200 m, time subsets)",
             transform=axR.transAxes, ha="right", va="bottom", fontsize=9.8, color=MUTED,
             bbox=dict(boxstyle="round,pad=0.4", fc=LIGHT, ec=GRID))

    fig.suptitle("Rayleigh (molecular) calibration — how it works  (schematic)",
                 fontsize=15, fontweight="bold", color=INKTXT, x=0.5, y=1.02)
    save(fig, "fig_rayleigh_schematic.png")


# ====================================================================== 6
# Sign error: aerosol backscatter goes negative (illustrative)
def fig_signerror_negative():
    z = np.linspace(0.1, 5.0, 300)               # km
    # physical aerosol backscatter (v1.1, +2 sign): positive, decaying with a layer
    good = 0.9 * np.exp(-z / 1.1) + 0.5 * np.exp(-((z - 0.8) / 0.5) ** 2) + 0.04
    # buggy (v1.0, -2 sign): exponential factor with wrong sign drags it negative aloft
    bug = good * np.exp(-(z) * 0.9) - 0.55 * (1 - np.exp(-z / 1.3))

    fig, ax = plt.subplots(figsize=(11.2, 5.0))
    ax.axhline(0, color=MUTED, lw=1.2)
    ax.fill_between(z, bug, 0, where=(bug < 0), color=RED, alpha=0.16, interpolate=True)
    ax.plot(z, good, color=GREEN, lw=3.0, label=r"v1.1  (+2 sign, corrected) — physical")
    ax.plot(z, bug,  color=RED,   lw=3.0, ls="--", label=r"v1.0  ($-2$ sign error) — unphysical")
    ax.text(3.4, -0.32, "aerosol backscatter\ndriven negative", color=RED, fontsize=12,
            fontweight="bold", ha="center")
    ax.set_xlabel("altitude (km AGL)", fontsize=12.5)
    ax.set_ylabel(r"retrieved aerosol backscatter $\beta_{aer}$ (a.u.)", fontsize=12.5)
    ax.set_xlim(0, 5); ax.set_ylim(-0.7, 1.5)
    style(ax); ax.grid(True, color=GRID, linewidth=0.7)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    ax.set_title(r"The Klett sign error: $\exp(-2\,q_t)$ instead of $\exp(+2\,q_t)$  (schematic)",
                 fontsize=14.5, fontweight="bold", color=INKTXT, loc="left", pad=10)
    save(fig, "fig_signerror_negative.png")


if __name__ == "__main__":
    fig_method_ranking()
    fig_v2_vs_v11()
    fig_minimpl_precision()
    fig_signfix_impact()
    fig_rayleigh_schematic()
    fig_signerror_negative()
    print("done")
