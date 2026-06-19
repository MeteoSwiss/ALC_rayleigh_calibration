"""Draw a flowchart of the 'optimal' molecular-window detection method."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/molecular_methods")
OUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(15, 8.5))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

BLUE = "#dbe7f3"; BLUEE = "#1f77b4"
RED = "#f6d6d6"; REDE = "#d62728"
GREEN = "#d7ecd9"; GREENE = "#2ca02c"
GREY = "#ececec"; GREYE = "#7f7f7f"
ORANGE = "#fde3cf"; ORANGEE = "#ff7f0e"


def box(cx, cy, w, h, text, fc, ec, fs=10, weight="normal"):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0.6,rounding_size=2",
                                fc=fc, ec=ec, lw=1.6, zorder=2))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, zorder=3, weight=weight)


def arrow(x1, y1, x2, y2, color="#333333", text=None, tx=0, ty=0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=18, lw=1.8, color=color, zorder=1))
    if text:
        ax.text((x1 + x2) / 2 + tx, (y1 + y2) / 2 + ty, text, ha="center", va="center",
                fontsize=9, color=color, style="italic")


# --- nodes ---
box(35, 95, 52, 7, "INPUT  •  night profile stack  signal(t, z)   +   molecular profile  p_mol(z)",
    BLUE, BLUEE, 10, "bold")

box(35, 84, 52, 8,
    "TIME-RESOLVED AEROSOL/CLOUD FLAGGING  (the distinctive step)\n"
    "per altitude: flag cells  signal/p_mol  >  median + 4·MAD\n"
    "→ removes aerosol present only part of the night, a passing cloud, BL aerosol",
    RED, REDE, 9.5)

box(35, 73, 52, 6.5, "TIME-CLEANED MEAN PROFILE\nmean over UN-flagged cells per altitude", GREEN, GREENE, 10)

box(35, 62.5, 52, 6.5, "GRID SEARCH over windows (centre × half-length)\nfit  signal = a·p_mol + b  in each window", GREY, GREYE, 10)

box(35, 44, 60, 16,
    "ELIGIBILITY GATES  (a window must pass ALL)\n"
    "• start ≥ 2 km AGL            (above the boundary-layer aerosol)\n"
    "• R² ≥ 0.5 , slope a > 0 , |b| < a   (real, positive molecular fit)\n"
    "• shape residual ≤ 12 %        (Rayleigh-shape match)\n"
    "• scattering ratio R ≤ 1.1     (aerosol-free; Wiegner & Geiß)\n"
    "• in-window SNR ok             (low ratio scatter)\n"
    "• |slope − median ratio| ≤ 15 %  (proportional; = pipeline QC)\n"
    "• temporal CV ≤ 0.5            (steady in time = molecular, not aerosol)",
    BLUE, BLUEE, 9.5)

box(86, 44, 22, 8, "no window passes\n→ NO CALIBRATION\n(flag −2; Kalman skips)", RED, REDE, 9.5, "bold")

box(35, 27, 60, 8,
    "COMPOSITE QUALITY SCORE  (maximise)\n"
    "Q = R² − 0.25·|R−1| − 0.20·resid − 0.10·SNR − 0.35·tCV − 0.20·rel + 0.10·n",
    GREEN, GREENE, 9.5)

box(35, 16, 60, 5.5, "SELECT the max-Q window  =  molecular reference window", GREY, GREYE, 10, "bold")

box(35, 6, 60, 6,
    "C_L = median( signal / p_mol )  in the window\n→ Klett β_att  +  downstream pipeline QC", ORANGE, ORANGEE, 10)

# --- arrows ---
arrow(35, 91.5, 35, 88.2)
arrow(35, 80.0, 35, 76.3)
arrow(35, 69.7, 35, 65.9)
arrow(35, 59.2, 35, 52.2)
arrow(35, 35.8, 35, 31.2)
arrow(35, 23.0, 35, 18.9)
arrow(35, 13.2, 35, 9.2)
arrow(65, 44, 75, 44, REDE)

ax.text(50, 99.2, "The 'optimal' molecular-window detection method", ha="center", va="center",
        fontsize=14, weight="bold")

fig.savefig(OUT / "optimal_flowchart.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("saved optimal_flowchart.png")
