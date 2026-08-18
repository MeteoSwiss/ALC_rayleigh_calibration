# -*- coding: utf-8 -*-
"""P-space (beta/z^2, homoscedastic) background estimation + noise-vs-range test +
3-hourly noise/offset during the 25.5-h hood window coloured by laser temperature."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from validation.paper._cl61_dark_windows import load, WINDOWS

# ---- pooled hood data ----
pools = []
for a, b in WINDOWS:
    X, rng = load(datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M"))
    pools.append(X)
X = np.vstack(pools) * 1e6                     # beta_att, Mm-1 sr-1
zkm = rng / 1000.0
P = X / zkm[None, :] ** 2                      # NON-range-corrected shape (homoscedastic noise)

# per-gate robust stats in P space
P_med = np.nanmedian(P, axis=0)
P_mad = 1.4826 * np.nanmedian(np.abs(P - P_med[None, :]), axis=0)
k = int(round(330 / np.median(np.diff(rng)))); k += (k + 1) % 2
P_smooth = medfilt(np.nan_to_num(P_med), k)
b_new = P_smooth * zkm ** 2                    # back to beta space
b_new[rng < 300] = 0.0
dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
np.savez("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz",
         rng=rng, b_dark=dk["b_dark"], b_smooth=b_new, b_mad=dk["b_mad"],
         P_med=P_med, P_mad=P_mad, n_profiles=X.shape[0])
print("b_dark (P-space smoothed, x z^2):  3 km %+.4f   5 km %+.4f   8 km %+.4f   12 km %+.4f"
      % tuple(b_new[np.argmin(np.abs(rng - z))] for z in (3000, 5000, 8000, 12000)))
# is sigma_P flat with range? (the beta-space fan = z^2 amplification, or real?)
for z in (2000, 5000, 10000, 14000):
    j = np.argmin(np.abs(rng - z))
    print(f"sigma_P @{z/1000:4.0f} km = {np.nanmedian(P_mad[max(0,j-20):j+20]):.5f}  "
          f"-> sigma_beta = {np.nanmedian(P_mad[max(0,j-20):j+20])*zkm[j]**2:.4f}")

# ---- 3-hourly noise + temperature during the 25.5-h window ----
t1, t2 = datetime(2026, 5, 26, 11, 45), datetime(2026, 5, 27, 13, 15)
tt, TL = [], []
for ds in ("20260526", "20260527"):
    f = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610") / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_C{ds}.nc"
    with Dataset(f) as nc:
        tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
        tv = np.asarray(nc.variables["time"][:], "f8")
        t = [datetime(1970, 1, 1) + timedelta(seconds=x) for x in (tv * 86400.0 if "day" in tu else tv)]
        Tl = np.asarray(nc.variables["temperature_laser"][:], "f8").ravel()
        for x, T in zip(t, Tl):
            if t1 <= x <= t2:
                tt.append(x); TL.append(T)
tt = np.array(tt); TL = np.array(TL)
Xw = pools[1] * 1e6
Pw = Xw / zkm[None, :] ** 2
assert len(tt) == Xw.shape[0]
blocks = []
b0 = t1
while b0 < t2:
    sel = (tt >= b0) & (tt < b0 + timedelta(hours=3))
    if sel.sum() > 60:
        med = np.nanmedian(Pw[sel], axis=0)
        sig = 1.4826 * np.nanmedian(np.abs(Pw[sel] - med[None, :]), axis=0)
        blocks.append((b0, float(np.nanmedian(TL[sel])), med, sig))
    b0 += timedelta(hours=3)

fig, ax = plt.subplots(1, 3, figsize=(17, 6))
Ts = [b[1] for b in blocks]
norm = plt.Normalize(min(Ts), max(Ts))
for b0, T, med, sig in blocks:
    c = cm.plasma(norm(T))
    ax[0].plot(medfilt(sig, k), rng, lw=1.4, color=c)
    ax[1].plot(medfilt(med, k), rng, lw=1.4, color=c, label=b0.strftime("%d %Hh"))
sm = cm.ScalarMappable(norm=norm, cmap="plasma")
fig.colorbar(sm, ax=ax[1]).set_label("median laser temperature [degC]")
ax[0].set_xlabel("sigma_P (robust) [Mm-1 sr-1 km-2]"); ax[0].set_title("(a) 3-hourly NOISE in P space")
ax[1].set_xlabel("median P [Mm-1 sr-1 km-2]"); ax[1].set_title("(b) 3-hourly OFFSET in P space")
ax[1].set_xlim(-0.004, 0.004); ax[0].set_xlim(0, 0.02)
for a_ in ax[:2]:
    a_.set_ylim(0, 15000); a_.grid(alpha=0.3); a_.set_ylabel("range [m]")
zb = (rng >= 3000) & (rng <= 6000)
ax[2].scatter([b[1] for b in blocks], [float(np.nanmean(b[2][zb] * zkm[zb] ** 2)) for b in blocks],
              c=[b[1] for b in blocks], cmap="plasma", s=60, edgecolors="k")
ax[2].set_xlabel("median laser temperature [degC]"); ax[2].set_ylabel("offset b(3-6 km) [Mm-1 sr-1]")
ax[2].set_title("(c) 3-6 km offset vs laser temperature"); ax[2].grid(alpha=0.3)
fig.suptitle("CL61 25.5-h hood test: 3-hourly noise and offset in NON-range-corrected space, coloured by laser temperature",
             fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_hood_3hourly.png"
fig.savefig(out, dpi=150)
print("saved", out)
