"""High-level flowchart of the liquid-cloud calibration procedure: paper Figure A3-1.

No input data (the text is hard-coded). Writes fig_cloud_flowchart.png into doc/reports/figs_multiple_scattering/
and, if given, into the folder passed as first argument. The paper image is byte-identical to the output of
the run with matplotlib 3.10.3 (2026-07); newer matplotlib versions differ by a few pixels.

Run:  python validation/_cloud_flowchart.py [<extra output folder>]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO_FIGS = Path(__file__).resolve().parents[1] / "doc" / "reports" / "figs_multiple_scattering"
EXTRA_OUT = [Path(a) for a in sys.argv[1:2]]   # optional extra output folder (was .../ms_fov)

COL = dict(data="#d9e8f5", corr="#fde9c9", test="#e7f3d9", calc="#f5d9d9", out="#d9f0ec")


def box(ax, x, y, w, h, text, color, fs=9.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.05",
                 fc=color, ec="#555555", lw=1.2, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight=("bold" if bold else "normal"), zorder=3)


def arrow(ax, x0, y0, x1, y1, style="-|>", color="#333333"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                 mutation_scale=16, lw=1.5, color=color, zorder=1))


fig, ax = plt.subplots(figsize=(13.5, 9))
ax.set_xlim(0, 12); ax.set_ylim(0, 15); ax.axis("off")

box(ax, 6, 14.2, 6.4, 0.9, "Attenuated backscatter $\\beta'$(z,t)   (L1 rcs$_0$ or L2)", COL["data"], 11, True)

# corrections column
ax.text(6, 13.2, "Profile corrections", ha="center", fontsize=10.5, style="italic", color="#8a5a00")
box(ax, 2.3, 12.2, 3.6, 0.95, "Multiple scattering\n$\\beta' \\leftarrow \\eta(z)\\,\\beta'$   (Hogan 2006)", COL["corr"], 9)
box(ax, 6.0, 12.2, 3.6, 0.95, "Water-vapour (910 nm)\n$\\beta' \\leftarrow \\beta' / T^2_{WV}(z)$  (CAMS+HITRAN)", COL["corr"], 9)
box(ax, 9.7, 12.2, 3.6, 0.95, "Instrument health\nquality flag / laser energy /\nwindow transmission gate", COL["corr"], 8.3)
for x in (2.3, 6.0, 9.7):
    arrow(ax, 6, 13.75, x, 12.68)
    arrow(ax, x, 11.72, 6, 11.1)

# cloud selection tests
box(ax, 6, 10.6, 7.4, 0.55, "CLOUD SELECTION  —  a profile is kept only if ALL pass:", "#eeeeee", 10, True)
tests = [
    "Totally attenuating: signal 300 m ABOVE peak $\\geq$ 20$\\times$ weaker",
    "Isolated base: signal 300 m BELOW peak $\\geq$ 20$\\times$ weaker",
    "Low aerosol: sub-cloud backscatter < 5 % of column (0.1-2.4 km)",
    "Cloud base 0.5-2.4 km AGL   (overlap / integration limits)",
    "Temporal consistency: 5 consecutive S within $\\pm$10 % of mean",
]
y = 9.75
for t in tests:
    box(ax, 6, y, 7.4, 0.62, t, COL["test"], 8.8)
    y -= 0.75
arrow(ax, 6, 10.33, 6, 10.06)

# integration + aerosol transmission + coefficient
arrow(ax, 6, y + 0.38, 6, y - 0.05)
box(ax, 6, y - 0.55, 6.6, 0.9,
    "Integrate  $B=\\int_{0.1\\,km}^{z_{cb}} \\beta'\\,dz$   +  sub-cloud aerosol two-way\n"
    "transmission correction ($S_{aer}$=50 sr, 20-70 sr bounds)", COL["calc"], 9)
y2 = y - 1.55
arrow(ax, 6, y - 1.0, 6, y2 + 0.05)
box(ax, 6, y2 - 0.5, 6.6, 0.95,
    "Per-profile coefficient   $C = S_{apparent}/S_{liq}$,   $S_{liq}=18.8$ sr\n"
    "$S_{apparent} = 1/(2B)$    (O'Connor 2004; Hopkin 2019)", COL["calc"], 9.5)
y3 = y2 - 1.5
arrow(ax, 6, y2 - 1.0, 6, y3 + 0.05)
box(ax, 6, y3 - 0.5, 6.6, 0.9,
    "Daily median C  $\\rightarrow$  Kalman-smoothed calibration constant\n"
    "$\\beta_{true} = C\\,\\beta'$   (relative gain of the instrument)", COL["out"], 9.5, True)

fig.suptitle("Liquid-cloud calibration procedure (O'Connor/Hopkin, as implemented)",
             fontweight="bold", fontsize=13, y=0.97)
fig.tight_layout()
for d in [REPO_FIGS] + EXTRA_OUT:
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / "fig_cloud_flowchart.png", dpi=160, bbox_inches="tight")
print("  fig_cloud_flowchart.png")
