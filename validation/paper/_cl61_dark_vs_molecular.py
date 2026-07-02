"""Uncorrected dark-signal histograms in beta_att space vs the molecular return: the direct
estimate of the offset's impact on the molecular (Rayleigh) calibration per window altitude."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
from validation.paper._cl61_dark_windows import load, WINDOWS, COLS
from datetime import datetime

ALT, CL = 491.0, 1.4252
dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
rng, b_dark = dk["rng"], dk["b_dark"]
atm = load_standard_atmosphere(Path("calibration/data/standard_atmosphere_US_1976_50km.csv"), rng + ALT)
mol = calculate_molecular_properties(atm.temperature, atm.pressure, rng, 910.0e-9)
molsig = CL * mol.beta_mol * mol.transmission * 1e6 * 0.78     # x typical T2_wv -> expected signal

print("z[km]  b_dark    molecular signal   b/mol = beta bias = dC_L if window here")
for z in (2000, 3000, 4000, 5000, 6000):
    j = np.argmin(np.abs(rng - z))
    bs = float(np.nanmedian(b_dark[max(0, j-5):j+5]))
    print(f"{z/1000:4.0f}   {bs:+.4f}   {molsig[j]:.4f}            {100*bs/molsig[j]:+6.1f}%")

fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for i, (a, b) in enumerate(WINDOWS):
    X, r = load(datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M"))
    for k, (lo, hi) in enumerate([(2000, 4000), (4000, 6000)]):
        zb = (r >= lo) & (r <= hi)
        v = (X[:, zb] * 1e6).ravel(); v = v[np.isfinite(v)]
        p = np.percentile(v, [1, 99])
        ax[k].hist(v, 120, range=tuple(p), histtype="step", density=True, color=COLS[i], lw=1.2,
                   label=f"{a[5:10]} mean {np.mean(v):+.3f}")
for k, (lo, hi) in enumerate([(2000, 4000), (4000, 6000)]):
    j = (rng >= lo) & (rng <= hi)
    ms = float(np.nanmean(molsig[j]))
    ax[k].axvline(0, color="k", lw=1)
    ax[k].axvline(ms, color="#d62728", lw=2, ls="--", label=f"molecular signal {ms:.3f}")
    ax[k].axvline(float(np.nanmean(b_dark[j])), color="0.3", lw=2, ls=":", label=f"dark offset {np.nanmean(b_dark[j]):+.3f}")
    ax[k].set_title(f"uncorrected dark $\\beta_{{att}}$, {lo/1000:.0f}-{hi/1000:.0f} km vs the molecular return")
    ax[k].set_xlabel("$\\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]"); ax[k].grid(alpha=0.3); ax[k].legend(fontsize=8)
fig.suptitle("CL61 dark signal (no correction) against the Rayleigh-calibration target signal", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_dark_vs_molecular.png", dpi=150)
print("saved fig_cl61_dark_vs_molecular.png")
