"""CL61 (ident C) nightly-mean profiles before/after dark-offset correction, vs both calibrated
molecular models. Clean nights; shaded = typical Rayleigh fit window."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from validation.paper._cl61_dark_baseline_probe import mean_profiles, rng_ref, model, CL_CLOUD, CL_RAY

dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
k = max(1, int(round(300 / np.median(np.diff(dk["rng"])))))
b_s = np.convolve(np.nan_to_num(dk["b_dark"]), np.ones(k) / k, mode="same")
b_s[dk["rng"] < 300] = 0.0
b = np.interp(rng_ref, dk["rng"], b_s)
CLEAN = ("20260402", "20260422", "20260424")
fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
for ax, ds in zip(axes, CLEAN):
    prof = next(p for d, p, n in mean_profiles if d == ds)
    ax.plot(prof, rng_ref, color="0.55", lw=1.4, label="original nightly mean")
    ax.plot(prof - b, rng_ref, color="#1f77b4", lw=1.6, label="offset-corrected")
    ax.plot(CL_CLOUD * model, rng_ref, "--", color="#d62728", lw=2.0,
            label=r"C$_L$(cloud) $\times\ \beta_{mol}$T$^2_{mol}$T$^2_{wv}$")
    ax.plot(CL_RAY * model, rng_ref, ":", color="0.3", lw=1.5, label=r"C$_L$(Rayleigh) $\times$ model")
    ax.axhspan(3000, 5000, color="orange", alpha=0.10)
    ax.set_xlim(-0.02, 0.25); ax.set_ylim(0, 9000); ax.grid(alpha=0.3)
    ax.set_title(f"night {ds[:4]}-{ds[4:6]}-{ds[6:]}")
    ax.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]")
axes[0].set_ylabel("range AGL [m]")
axes[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Payerne CL61 (ident C) nightly-mean profiles before/after dark-offset correction - "
             "shaded: typical Rayleigh fit window", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_profiles_beforeafter.png", dpi=150)
print("saved")
