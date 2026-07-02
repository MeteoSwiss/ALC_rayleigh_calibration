"""Draw the two analysis flowcharts for the paper-validation report (landscape PNGs):
  fig_flow_intercompare.png : multi-ceilometer station intercomparison pipeline
  fig_flow_earlinet.png     : ceilometer vs EARLINET research-lidar pipeline
Same visual style as validation/make_optimal_flowchart.py."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")

BLUE = "#dbe7f3"; BLUEE = "#1f77b4"
RED = "#f6d6d6"; REDE = "#d62728"
GREEN = "#d7ecd9"; GREENE = "#2ca02c"
GREY = "#ececec"; GREYE = "#7f7f7f"
ORANGE = "#fde3cf"; ORANGEE = "#ff7f0e"
PURPLE = "#e6def0"; PURPLEE = "#9467bd"


def _mk(figsize=(16, 9)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    return fig, ax


def box(ax, cx, cy, w, h, text, fc, ec, fs=9.5, weight="normal"):
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0.5,rounding_size=1.6",
                                fc=fc, ec=ec, lw=1.6, zorder=2))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, zorder=3, weight=weight)


def arrow(ax, x1, y1, x2, y2, color="#333333", text=None, tx=1.2, ty=0, fs=8.5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16,
                                 color=color, lw=1.5, zorder=1))
    if text:
        ax.text((x1 + x2) / 2 + tx, (y1 + y2) / 2 + ty, text, fontsize=fs, color=color, ha="left")


def flow_intercompare():
    fig, ax = _mk()
    # inputs row
    box(ax, 12, 90, 20, 9, "E-PROFILE L2 monthly\n$\\beta_{att}$ 30 s native\n(per channel)", BLUE, BLUEE, weight="bold")
    box(ax, 38, 90, 20, 9, "Calibration series: lidar constant C$_L$\nnightly Rayleigh (eprof_v2) /\ndaily cloud (O'Connor) + Kalman", ORANGE, ORANGEE)
    box(ax, 64, 90, 18, 9, "CAMS monthly\nwater vapour", PURPLE, PURPLEE)
    box(ax, 87, 90, 18, 9, "US standard\natmosphere (molecular)", PURPLE, PURPLEE)
    # calibrate
    box(ax, 25, 74, 30, 8, "Apply calibration (linear interp of daily Kalman C$_L$)\nRayleigh: $\\beta$·(calc/C$_L$)  cloud: $\\beta$·(default/C$_L$)\nCL61: default for both (calc = vendor metadata)", GREY, GREYE)
    arrow(ax, 12, 85.5, 20, 78); arrow(ax, 38, 85.5, 30, 78)
    # WV + wavelength
    box(ax, 25, 62, 30, 8, "910 nm channels: divide by two-way WV\ntransmission T$^2_{wv}$ (per month, CAMS profile)", GREY, GREYE)
    arrow(ax, 64, 85.5, 33, 66); arrow(ax, 25, 70, 25, 66)
    box(ax, 25, 50, 30, 8, "Normalise wavelength to reference\nAngstrom $\\alpha$=1; Mini-MPL: molaer\n($\\beta_{aer}$/2 + $\\beta_{mol}$(1064) exact)", GREY, GREYE)
    arrow(ax, 87, 85.5, 38, 54, tx=0); arrow(ax, 25, 58, 25, 54)
    # screening
    box(ax, 72, 68, 34, 10, "SCREEN (science stream)\nquality_flag > 0 · clouds (any CBH 0-20 km)\nfog / vertical visibility · ±15 min expansion", RED, REDE, weight="bold")
    arrow(ax, 40, 50, 58, 63.5)
    # retime + SNR
    box(ax, 72, 52, 34, 10, "HOURLY MEDIAN (60 min $\\geq$ 30 min)\nbin kept only if $\\geq$ 30 min of samples\nremove gates with window SNR < 3\n(SNR = med / ($\\sigma_{rob}/\\sqrt{n}$))", RED, REDE, weight="bold")
    arrow(ax, 72, 63, 72, 57.5)
    # sync + grid
    box(ax, 72, 37, 34, 8, "UNION time grid across channels\ncommon altitude grid (bin average)", GREY, GREYE)
    arrow(ax, 72, 47, 72, 41)
    # outputs
    box(ax, 30, 20, 32, 12, "MEDIAN PROFILES (panel a)\nonly hours where ALL channels report\n(common-time sample, N annotated)\nmedian ± IQR, noise-floor truncated", GREEN, GREENE, weight="bold")
    box(ax, 72, 20, 34, 12, "PAIRWISE STATS vs reference (500-3000 m)\nrelbias (mean) · med relbias · r (linear)\nlog r (positive pairs) · N\nscatter (panel b) + difference hist (panel c)", GREEN, GREENE, weight="bold")
    arrow(ax, 62, 33, 40, 26.5); arrow(ax, 72, 33, 72, 26.5)
    ax.text(50, 4, "display stream: quality_flag masking only (curtain panels d-g show clouds/fog for context)",
            ha="center", fontsize=9, style="italic", color=GREYE)
    fig.suptitle("Multi-ceilometer station intercomparison — analysis flow", fontweight="bold", fontsize=13)
    fig.savefig(OUT / "fig_flow_intercompare.png", dpi=170, bbox_inches="tight"); plt.close(fig)
    print("fig_flow_intercompare.png")


def flow_earlinet():
    fig, ax = _mk()
    # EARLINET branch (left)
    box(ax, 16, 90, 26, 9, "EARLINET SCC L2 b1064 files\nparticle backscatter (m$^{-1}$sr$^{-1}$)\ncloud-screened by SCC", BLUE, BLUEE, weight="bold")
    box(ax, 16, 76, 26, 9, "dedup (start,end): keep highest\n(qc level, version) · window midpoint\n(time_bounds) as profile time", GREY, GREYE)
    arrow(ax, 16, 85.5, 16, 80.5)
    box(ax, 16, 62, 26, 9, "interp to 15 m grid\nNaN below instrument overlap\n(no constant fill)", RED, REDE, weight="bold")
    arrow(ax, 16, 71.5, 16, 66.5)
    box(ax, 16, 46, 26, 12, "$\\beta_{att}$ = ($\\beta_p$+$\\beta_{mol}$)·T$^2$\nLR: per-scene SCC value (else 50 sr)\nOD: lowest trusted extinction extended\nto ground; molecular from std atm", GREY, GREYE)
    arrow(ax, 16, 57.5, 16, 52)
    # CHM branch (right)
    box(ax, 62, 90, 26, 9, "E-PROFILE L2 CHM15k\n$\\beta_{att}$ 30 s native", BLUE, BLUEE, weight="bold")
    box(ax, 88, 90, 18, 9, "Rayleigh calibration\n(eprof_v2 + Kalman)", ORANGE, ORANGEE)
    box(ax, 75, 76, 30, 8, "calibrate: $\\beta$·(calc/C$_L$)", GREY, GREYE)
    arrow(ax, 62, 85.5, 68, 80); arrow(ax, 88, 85.5, 82, 80)
    box(ax, 75, 63, 30, 9, "SCREEN: quality_flag · clouds (CBH)\nfog/VV · ±15 min expansion", RED, REDE, weight="bold")
    arrow(ax, 75, 72, 75, 67.5)
    # matching (center)
    box(ax, 50, 33, 44, 13, "MATCH per EARLINET profile\nCHM median within ±30 min ($\\geq$ 30 min of samples required)\nremove gates with window SNR < 3\ninterp to EARLINET 15 m grid · skip all-NaN pairs", RED, REDE, weight="bold")
    arrow(ax, 16, 40, 34, 39.5); arrow(ax, 75, 58.5, 62, 39.5)
    # outputs
    box(ax, 24, 12, 34, 10, "median matched profile ± IQR (panel a)\ncurtains by date (panels c, d)", GREEN, GREENE, weight="bold")
    box(ax, 72, 12, 36, 10, "density scatter 500-5000 m (panel b)\nrelbias · med relbias · r · log r · N matched", GREEN, GREENE, weight="bold")
    arrow(ax, 40, 26.5, 32, 17.5); arrow(ax, 60, 26.5, 68, 17.5)
    fig.suptitle("Ceilometer vs EARLINET research lidar — analysis flow", fontweight="bold", fontsize=13)
    fig.savefig(OUT / "fig_flow_earlinet.png", dpi=170, bbox_inches="tight"); plt.close(fig)
    print("fig_flow_earlinet.png")


if __name__ == "__main__":
    flow_intercompare()
    flow_earlinet()
