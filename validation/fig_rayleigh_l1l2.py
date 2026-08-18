"""
fig_rayleigh_l1l2.py — figure for the L1<->L2 Rayleigh root-cause finding. Two panels:
(A) demonstration day (CHM15k 06610, 2026-02-25, full D: files): lidar constant C_L by method and
    data level/grid, showing v1.1 works on both, v2/earlinet reject NATIVE L1 but agree with L2 once
    L1 is binned to the L2 grid (30 m x 300 s); v1.2 rejects this day at both levels.
(B) 5-stream CHM15k check (rayleigh_l1_grid_check.py): valid fraction per combo.
Values are the measured results quoted in network_v2_vs_v11_report.md.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_l1l2")
OUT.mkdir(parents=True, exist_ok=True)

# (A) demonstration day: C_L in 1e11 units; np.nan = rejected (flag != 1)
methods = ["v1.1", "v2", "earlinet", "v1.2"]
variants = ["L1 native", "L1 binned\n(30 m x 300 s)", "L2 (30 m x 300 s)"]
# rows = method, cols = variant
CL = np.array([
    [5.17, np.nan, 5.04],   # v1.1 (native L1 works directly; "binned" not needed -> n/a shown as gap)
    [np.nan, 5.00, 4.92],   # v2: native rejected, binned matches L2
    [np.nan, 5.11, 5.09],   # earlinet: native rejected, binned matches L2
    [np.nan, np.nan, np.nan],  # v1.2: rejected at both levels this day
])
# v1.1 native shown under the "L1 native" column; place it there instead of binned:
CL_v11_native = 5.17

# (B) 5-stream check valid %
combos = ["L1 v1.1\nnative", "L1 v2\nnative", "L1 v2\nbinned", "L2 v2"]
validpct = [13, 9, 15, 10]
medCL = [1.873, 1.919, 2.018, 1.999]  # 1e11

fig, (axA, axB) = plt.subplots(1, 2, figsize=(18, 6.4))

# Panel A
x = np.arange(len(methods)); w = 0.26
colors = ["#7f7f7f", "#1f77b4", "#2ca02c"]
for k, var in enumerate(variants):
    vals = CL[:, k].copy()
    if k == 0:
        vals[0] = CL_v11_native      # v1.1 native value in the native column
    bars = axA.bar(x + (k - 1) * w, np.nan_to_num(vals, nan=0.0), w, label=var, color=colors[k])
    for b, v in zip(bars, vals):
        if np.isfinite(v):
            axA.annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=8)
        else:
            axA.annotate("rejected", (b.get_x() + b.get_width() / 2, 0.1), ha="center", va="bottom",
                         fontsize=7, color="#d62728", rotation=90)
axA.set_xticks(x); axA.set_xticklabels(methods)
axA.set_ylabel(r"lidar constant $C_L$  (10$^{11}$, a.u.)")
axA.set_title("CHM15k 06610, 2026-02-25 (full day): $C_L$ by method and grid\n"
              "same data — L1 binned to the L2 grid agrees with L2", fontsize=11)
axA.set_ylim(0, 6.2); axA.legend(fontsize=9); axA.grid(axis="y", alpha=0.3)

# Panel B
xb = np.arange(len(combos))
bars = axB.bar(xb, validpct, color=["#7f7f7f", "#d62728", "#1f77b4", "#2ca02c"])
for b, v, c in zip(bars, validpct, medCL):
    axB.annotate(f"{v}%\nC_L {c:.2f}", (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=8)
axB.set_xticks(xb); axB.set_xticklabels(combos)
axB.set_ylabel("valid Rayleigh nights (% of sampled days)")
axB.set_title("5 CHM15k streams (sparse winter sample): native-L1 v2 is suppressed;\n"
              "binning L1 to the L2 grid recovers it; v1.1 is grid-robust", fontsize=11)
axB.set_ylim(0, 20); axB.grid(axis="y", alpha=0.3)

fig.suptitle("L1 vs L2 Rayleigh: same data — a method/grid interaction, not a data difference", fontsize=14)
fig.tight_layout(); fig.savefig(OUT / "fig_rayleigh_l1l2.png", dpi=130); plt.close(fig)
print("saved", OUT / "fig_rayleigh_l1l2.png")
