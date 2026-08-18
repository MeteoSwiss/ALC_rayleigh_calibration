# -*- coding: utf-8 -*-
"""The correction profiles b_dark(z) actually subtracted from L1 rcs_0 for the CL61 and the CHM15k,
i.e. the physical-model offset that the recalibration removes. Shown in the range-corrected rcs_0
units (what is subtracted, = P_model x z^2) with the pooled hood median for context, plus the
non-range-corrected P view where the model is defined."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
import numpy as np
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/"
cl61 = np.load(D + "cl61_b_dark.npz")
chm = np.load(D + "chm15k_b_dark.npz")

fig, ax = plt.subplots(1, 2, figsize=(15, 6.4))

# --- CL61 ---
r = cl61["rng"]; zkm = r / 1000.0
Pm = cl61["P_med"]; bphys = cl61["b_phys"]           # b_phys in Mm-1 sr-1 (range-corrected)
dr = np.median(np.diff(r)); kk = int(round(600 / dr)); kk += 1 - kk % 2
a = ax[0]; a.axhline(0, color="0.6", lw=0.8)
a.plot(medfilt(np.nan_to_num(Pm * zkm ** 2), kk), r, "-", color="0.55", lw=1.2, label="pooled hood median")
a.plot(bphys, r, "-", color="#d62728", lw=2.4, label="correction b(z) subtracted (AC-coupling model)")
a.set_xlim(-0.05, 0.02); a.set_ylim(0, 15000); a.grid(alpha=0.3)
a.set_xlabel(r"$b(z)$ subtracted from rcs$_0$  [Mm$^{-1}$sr$^{-1}$]"); a.set_ylabel("range [m]")
a.set_title("CL61 correction (910 nm): fast lobe + slow AC undershoot"); a.legend(fontsize=9, loc="lower left")
a.axhspan(3000, 5000, color="orange", alpha=0.10)

# --- CHM15k ---
r2 = chm["rng"]; zkm2 = r2 / 1000.0
Pm2 = chm["P_med"]; bphys2 = chm["b_phys"]           # counts/s.m2 (range-corrected)
dr2 = np.median(np.diff(r2)); k2 = int(round(600 / dr2)); k2 += 1 - k2 % 2
b = ax[1]; b.axhline(0, color="0.6", lw=0.8)
S = 1e-3
b.plot(medfilt(np.nan_to_num(Pm2 * zkm2 ** 2), k2) * S, r2, "-", color="0.55", lw=1.2, label="pooled hood median")
b.plot(bphys2 * S, r2, "-", color="#d62728", lw=2.4, label="correction b(z) (single-exp over-subtraction model)")
b.set_ylim(0, 15000); b.grid(alpha=0.3)
b.set_xlabel(r"$b(z)$ subtracted from rcs$_0$  [10$^3$ counts s$^{-1}$ m$^2$]"); b.set_ylabel("range [m]")
b.set_title("CHM15k correction (1064 nm): single negative relaxation"); b.legend(fontsize=9, loc="lower right")
b.axhspan(3000, 5000, color="orange", alpha=0.10)

fig.suptitle("Offset correction profiles b(z) subtracted from L1 rcs$_0$ — physical models "
             "(orange band = Rayleigh fit window 3-5 km)", fontweight="bold", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = D + "fig_hood_correction_profiles.png"
fig.savefig(OUT, dpi=150); print("saved", OUT)
