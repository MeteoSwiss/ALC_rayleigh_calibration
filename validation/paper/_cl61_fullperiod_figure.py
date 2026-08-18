# -*- coding: utf-8 -*-
"""Full-period Payerne CL61 Rayleigh recalibration figure (2026-02-24..06-30): per-day C_L WITH
and WITHOUT the hood-offset correction, vs the CHM15k-anchored cloud constant. The correction both
UNLOCKS nights (the offset was spoiling the fit) and shifts the constant toward the cloud value."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CL_CLOUD = 1.4252
CSV = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_fullperiod_recal.csv"
df = pd.read_csv(CSV)
df["dt"] = pd.to_datetime(df.date.astype(str), format="%Y%m%d")

def good(flagcol, clcol, flags):
    m = df[flagcol].isin(flags) & df[clcol].notna() & (df[clcol] > 0) & (df[clcol] < 3)
    return df[m]

o1 = good("flag_orig", "CL_orig", [1]); c1 = good("flag_corr", "CL_corr", [1])
o0 = good("flag_orig", "CL_orig", [1, 0.5, 0]); c0 = good("flag_corr", "CL_corr", [1, 0.5, 0])
paired = df[df.flag_orig.isin([1, 0.5, 0]) & df.flag_corr.isin([1, 0.5, 0]) & df.CL_orig.notna() & df.CL_corr.notna()
            & (df.CL_orig > 0) & (df.CL_corr > 0)]

fig = plt.figure(figsize=(16, 6.6))
gs = fig.add_gridspec(1, 3, width_ratios=[2.4, 1, 1], wspace=0.28)

ax = fig.add_subplot(gs[0, 0])
ax.axhline(CL_CLOUD, color="#2ca02c", lw=2, ls="--", label=f"cloud constant (CHM15k-anchored) = {CL_CLOUD}")
ax.plot(o0.dt, o0.CL_orig, "o", ms=5, mfc="none", mec="0.5", label=f"original (flag≥0, n={len(o0)})")
ax.plot(c0.dt, c0.CL_corr, "o", ms=5, color="#1f77b4", alpha=0.55, label=f"offset-corrected (flag≥0, n={len(c0)})")
ax.plot(o1.dt, o1.CL_orig, "s", ms=8, mfc="none", mec="k", mew=1.4, label=f"original flag=1 (n={len(o1)})")
ax.plot(c1.dt, c1.CL_corr, "s", ms=8, color="#08519c", label=f"corrected flag=1 (n={len(c1)})")
ax.axhline(o1.CL_orig.median(), color="0.5", lw=1.2, ls=":")
ax.axhline(c1.CL_corr.median(), color="#08519c", lw=1.2, ls=":")
ax.set_ylabel(r"CL61 Rayleigh lidar constant $C_L$"); ax.set_xlabel("date (2026)")
ax.set_ylim(0.8, 1.8); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left", ncol=1)
ax.set_title("(a) full-period per-night $C_L$: original vs offset-corrected")

# paired arrows on common nights
axp = fig.add_subplot(gs[0, 1])
for _, r in paired.iterrows():
    axp.plot([0, 1], [r.CL_orig, r.CL_corr], "-", color="0.7", lw=0.8, zorder=1)
axp.plot(np.zeros(len(paired)), paired.CL_orig, "o", color="0.5", ms=6)
axp.plot(np.ones(len(paired)), paired.CL_corr, "o", color="#1f77b4", ms=6)
axp.axhline(CL_CLOUD, color="#2ca02c", lw=2, ls="--")
sh = 100 * (paired.CL_corr.values / paired.CL_orig.values - 1)
axp.set_xticks([0, 1]); axp.set_xticklabels(["original", "corrected"]); axp.set_xlim(-0.3, 1.3)
axp.set_ylim(0.8, 1.8); axp.grid(alpha=0.3); axp.set_ylabel(r"$C_L$")
axp.set_title(f"(b) paired nights (n={len(paired)})\nmedian shift {np.median(sh):+.0f}%")

axh = fig.add_subplot(gs[0, 2])
bins = np.linspace(0.8, 1.8, 22)
axh.hist(o0.CL_orig, bins=bins, orientation="horizontal", color="0.6", alpha=0.6, label="original")
axh.hist(c0.CL_corr, bins=bins, orientation="horizontal", color="#1f77b4", alpha=0.6, label="corrected")
axh.axhline(CL_CLOUD, color="#2ca02c", lw=2, ls="--")
axh.axhline(o0.CL_orig.median(), color="0.4", lw=1.4, ls=":")
axh.axhline(c0.CL_corr.median(), color="#08519c", lw=1.4, ls=":")
axh.set_ylim(0.8, 1.8); axh.set_xlabel("count (flag≥0)"); axh.grid(alpha=0.3); axh.legend(fontsize=8, loc="upper right")
axh.set_title("(c) distribution")

g_o = 100 * (o0.CL_orig.median() / CL_CLOUD - 1); g_c = 100 * (c0.CL_corr.median() / CL_CLOUD - 1)
fig.suptitle(f"Payerne CL61 Rayleigh — full record ({df.dt.min():%Y-%m-%d} to {df.dt.max():%Y-%m-%d}), "
             f"with/without hood-offset correction\n"
             f"correction UNLOCKS nights (flag=1: {len(o1)}→{len(c1)}; any: {len(o0)}→{len(c0)}); "
             f"paired shift {np.median(sh):+.0f}%; median gap to cloud {g_o:+.0f}%→{g_c:+.0f}%",
             fontweight="bold", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.92))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_fullperiod_recal.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
print(f"orig flag=1 n={len(o1)} med={o1.CL_orig.median():.3f} | corr flag=1 n={len(c1)} med={c1.CL_corr.median():.3f}")
print(f"any-flag: orig n={len(o0)} med={o0.CL_orig.median():.3f} gap {g_o:+.1f}% | corr n={len(c0)} med={c0.CL_corr.median():.3f} gap {g_c:+.1f}%")
print(f"paired n={len(paired)} median shift {np.median(sh):+.1f}%")
