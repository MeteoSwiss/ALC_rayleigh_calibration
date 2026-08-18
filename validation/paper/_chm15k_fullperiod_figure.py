# -*- coding: utf-8 -*-
"""Full-period Payerne CHM15k Rayleigh recalibration figure (counterpart of the CL61 one). The
CHM15k has no cloud-method anchor, so the reference is its own NATIVE median: a MATERIAL offset
would move the corrected constant away from native by a lot and unlock many nights (as the CL61
does); here the correction is small and within scatter -> the offset is not material."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/chm15k_fullperiod_recal.csv"
S = 1e-11                                   # plot C_L in 1e11 units
df = pd.read_csv(CSV); df["dt"] = pd.to_datetime(df.date.astype(str), format="%Y%m%d")

def good(fc, cc, fl):
    m = df[fc].isin(fl) & df[cc].notna() & (df[cc] > 0)
    return df[m]
o1, c1 = good("flag_orig", "CL_orig", [1]), good("flag_corr", "CL_corr", [1])
o0, c0 = good("flag_orig", "CL_orig", [1, 0.5, 0]), good("flag_corr", "CL_corr", [1, 0.5, 0])
paired = df[df.flag_orig.isin([1, 0.5, 0]) & df.flag_corr.isin([1, 0.5, 0]) & df.CL_orig.notna() & df.CL_corr.notna()
            & (df.CL_orig > 0) & (df.CL_corr > 0)]
REF = o0.CL_orig.median()                   # native median = the trusted reference (EARLINET-validated)
shift = 100 * (c0.CL_corr.median() / REF - 1)
psh = 100 * (paired.CL_corr.values / paired.CL_orig.values - 1)

fig = plt.figure(figsize=(16, 6.6))
gs = fig.add_gridspec(1, 3, width_ratios=[2.4, 1, 1], wspace=0.28)
ax = fig.add_subplot(gs[0, 0])
ax.axhline(REF * S, color="#d62728", lw=2, ls="--", label=f"native median (reference) = {REF*S:.2f}e11")
ax.plot(o0.dt, o0.CL_orig * S, "o", ms=5, mfc="none", mec="0.5", label=f"original (flag≥0, n={len(o0)})")
ax.plot(c0.dt, c0.CL_corr * S, "o", ms=5, color="#1f77b4", alpha=0.55, label=f"offset-corrected (flag≥0, n={len(c0)})")
ax.plot(o1.dt, o1.CL_orig * S, "s", ms=8, mfc="none", mec="k", mew=1.4, label=f"original flag=1 (n={len(o1)})")
ax.plot(c1.dt, c1.CL_corr * S, "s", ms=8, color="#08519c", label=f"corrected flag=1 (n={len(c1)})")
ax.axhline(c0.CL_corr.median() * S, color="#08519c", lw=1.2, ls=":")
ax.set_ylabel(r"CHM15k Rayleigh constant $C_L$  [$\times10^{11}$]"); ax.set_xlabel("date (2026)")
ax.set_ylim(4.5, 9.0); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
ax.set_title("(a) full-period per-night $C_L$: original vs offset-corrected")

axp = fig.add_subplot(gs[0, 1])
for _, r in paired.iterrows():
    axp.plot([0, 1], [r.CL_orig * S, r.CL_corr * S], "-", color="0.7", lw=0.8, zorder=1)
axp.plot(np.zeros(len(paired)), paired.CL_orig * S, "o", color="0.5", ms=6)
axp.plot(np.ones(len(paired)), paired.CL_corr * S, "o", color="#1f77b4", ms=6)
axp.axhline(REF * S, color="#d62728", lw=2, ls="--")
axp.set_xticks([0, 1]); axp.set_xticklabels(["original", "corrected"]); axp.set_xlim(-0.3, 1.3)
axp.set_ylim(4.5, 9.0); axp.grid(alpha=0.3); axp.set_ylabel(r"$C_L$ [$\times10^{11}$]")
axp.set_title(f"(b) paired nights (n={len(paired)})\nmedian shift {np.median(psh):+.1f}%")

axh = fig.add_subplot(gs[0, 2])
bins = np.linspace(4.5, 9.0, 22)
axh.hist(o0.CL_orig * S, bins=bins, orientation="horizontal", color="0.6", alpha=0.6, label="original")
axh.hist(c0.CL_corr * S, bins=bins, orientation="horizontal", color="#1f77b4", alpha=0.6, label="corrected")
axh.axhline(REF * S, color="#d62728", lw=2, ls="--")
axh.set_ylim(4.5, 9.0); axh.set_xlabel("count (flag≥0)"); axh.grid(alpha=0.3); axh.legend(fontsize=8, loc="upper right")
axh.set_title("(c) distribution")

fig.suptitle(f"Payerne CHM15k Rayleigh — full record ({df.dt.min():%Y-%m-%d} to {df.dt.max():%Y-%m-%d}), "
             f"with/without hood-offset correction\npaired shift {np.median(psh):+.1f}%, flag=1 {len(o1)}→{len(c1)} "
             f"(magnitude LIKE the CL61) — but the intercomparison moves CHM15k −7.1% AWAY from the cloud/EARLINET "
             f"anchor (CL61 moves +13.5%→+1.0% TOWARD it): correct the CL61, keep the CHM15k native",
             fontweight="bold", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.90))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_chm15k_fullperiod_recal.png"
fig.savefig(OUT, dpi=150); print("saved", OUT)
print(f"native med {REF*S:.3f}e11, corrected med {c0.CL_corr.median()*S:.3f}e11, shift {shift:+.1f}%, "
      f"paired {np.median(psh):+.1f}%, flag=1 {len(o1)}->{len(c1)}")
