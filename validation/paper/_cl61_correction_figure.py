"""Figure for report section 1c: per-night CL61 Rayleigh C_L before/after removing the measured
dark offset (values from _cl61_offset_correction_experiment run, 2026-07-02)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nights = ["20260225", "20260316", "20260328", "20260407", "20260410", "20260411", "20260417",
          "20260420", "20260421", "20260422", "20260423", "20260424", "20260521", "20260602"]
orig = {"20260316": 1.268, "20260328": 1.147, "20260410": 1.341, "20260411": 1.047,
        "20260417": 1.009, "20260420": 1.324, "20260602": 1.282}
corr = {"20260225": 1.139, "20260316": 1.318, "20260407": 1.202, "20260411": 1.286,
        "20260417": 1.227, "20260422": 1.270, "20260423": 1.269, "20260424": 1.043, "20260521": 1.536}
ns = sorted(set(orig) | set(corr))
x = np.arange(len(ns))
fig, ax = plt.subplots(figsize=(12, 5.5))
for i, n in enumerate(ns):
    a, b = orig.get(n), corr.get(n)
    if a is not None:
        ax.plot(i, a, "o", color="0.45", ms=8, zorder=3)
    if b is not None:
        ax.plot(i, b, "o", color="#1f77b4", ms=8, zorder=3)
    if a is not None and b is not None:
        ax.annotate("", xy=(i, b), xytext=(i, a),
                    arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.6))
ax.axhline(1.4252, color="#404040", lw=2, ls="--", label="C$_L$(cloud) = 1.425 (CHM15k-anchored)")
ax.axhline(1.251, color="0.45", lw=1.2, ls=":", label="operational Rayleigh Kalman (1.251)")
ax.plot([], [], "o", color="0.45", label="original")
ax.plot([], [], "o", color="#1f77b4", label="dark-offset corrected")
ax.set_xticks(x); ax.set_xticklabels([n[4:6] + "-" + n[6:] for n in ns], rotation=45, fontsize=8)
ax.set_ylabel("lidar constant C$_L$"); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower right")
ax.set_title("Payerne CL61 Rayleigh calibration per night: removing the measured dark offset "
             "(arrows: paired nights, +4/+22/+23 %)")
fig.tight_layout()
fig.savefig("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_offset_correction.png", dpi=150)
print("saved")
