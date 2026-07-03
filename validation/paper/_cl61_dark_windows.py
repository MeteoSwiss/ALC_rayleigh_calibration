"""Payerne CL61 covered-telescope (hood-on) dark measurements — full characterisation.
Three windows (2026-05-12 09:35-14:50, 05-26 11:45-13:15, 06-09 09:20-11:50); the 2026-05-12
01:35-09:35 period is DISCARDED (pre-cover: contains fog/cloud returns). Outputs: window-mean
beta_att, smoothed b_dark(z), non-range-corrected mean P(z)=beta/z^2 and per-band HISTOGRAMS of
the per-sample P (distribution center = the offset; shape = detector noise)."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WMO, IDENT = "0-20000-0-06610", "C"
L1 = Path("D:/E-PROFILE_L1_2026")
WINDOWS = [("2026-05-12 09:35", "2026-05-12 14:50"),
           ("2026-05-26 11:45", "2026-05-27 13:15"),
           ("2026-06-09 09:20", "2026-06-09 11:50")]
COLS = ["#1f77b4", "#ff7f0e", "#2ca02c"]

def load(t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    P, rng = [], None
    for ds in sorted(days):
        f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
            tv = np.asarray(nc.variables["time"][:], "f8")
            tt = np.array([datetime(1970, 1, 1) + timedelta(seconds=x)
                           for x in (tv * 86400.0 if "day" in tu else tv)])
            r = np.asarray(nc.variables["range"][:], "f8")
            x = np.asarray(nc.variables["rcs_0"][:], "f8")
            if x.shape[0] != tt.size:
                x = x.T
        s = (tt >= t1) & (tt <= t2)
        if s.any():
            P.append(x[s]); rng = r
    return (np.vstack(P), rng) if P else (None, None)

fig = plt.figure(figsize=(19, 10))
gs = fig.add_gridspec(2, 3)
axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
BANDS = [(2000, 4000), (4000, 8000), (8000, 14000)]
axH = [fig.add_subplot(gs[1, j]) for j in range(3)]
means, pools, rng = [], [], None
print("window                 N     dur[h]  b(3km)    b(5km)    b(8km)   [Mm-1 sr-1]")
for i, (a, b) in enumerate(WINDOWS):
    t1, t2 = datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M")
    X, rng = load(t1, t2)
    pools.append(X)
    m = np.nanmean(X, axis=0) * 1e6
    means.append(m)
    at = lambda z: float(m[np.argmin(np.abs(rng - z))])
    print(f"{a}..{b[11:]}  {X.shape[0]:5d}  {(t2-t1).total_seconds()/3600:5.1f}  {at(3000):+8.4f}  {at(5000):+8.4f}  {at(8000):+8.4f}")
    lab = f"{a[5:16]} (n={X.shape[0]})"
    zkm = rng / 1000.0
    axA.plot(m, rng, ".", ms=2.0, color=COLS[i], label=lab)
    axB.plot(m, rng, ".", ms=2.0, color=COLS[i], alpha=0.6)
    axC.plot(m / zkm**2, rng, ".", ms=2.0, color=COLS[i], alpha=0.6)
    Pns = X * 1e6 / zkm[None, :]**2
    for j, (lo, hi) in enumerate(BANDS):
        zb = (rng >= lo) & (rng <= hi)
        v = Pns[:, zb].ravel(); v = v[np.isfinite(v)]
        p1, p99 = np.percentile(v, [0.5, 99.5])
        axH[j].hist(v, bins=120, range=(p1, p99), histtype="step", density=True,
                    color=COLS[i], lw=1.3, label=f"{a[5:10]}: mean {np.mean(v):+.4f}")

# ROBUST estimator: per-gate MEDIAN over ALL pooled hood profiles (outlier/spike-immune),
# then a RUNNING MEDIAN over ~330 m (robust smooth), zero below 300 m (hood near-field).
from scipy.signal import medfilt
P = np.vstack(pools) * 1e6
b_med = np.nanmedian(P, axis=0)
k = int(round(330 / np.median(np.diff(rng))))
k += (k + 1) % 2
b_smooth = medfilt(np.nan_to_num(b_med), k)
b_smooth[rng < 300] = 0.0
b_mad = 1.4826 * np.nanmedian(np.abs(P - b_med[None, :]), axis=0)
np.savez("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz",
         rng=rng, b_dark=b_med, b_smooth=b_smooth, b_mad=b_mad, n_profiles=P.shape[0])
for z in (3000, 5000, 8000, 12000):
    j = np.argmin(np.abs(rng - z))
    print(f"robust b_dark @{z/1000:.0f} km: {b_smooth[j]:+.4f}  (MAD/sqrt(n) = {b_mad[j]/np.sqrt(P.shape[0]):.4f})")

axB.plot(b_smooth, rng, "k-", lw=2.4,
         label=f"per-gate median of {P.shape[0]} profiles, 330 m running median")
axC.plot(b_smooth / (rng / 1000.0) ** 2, rng, "k-", lw=2.2, label="robust median / z$^2$")
axA.set_xlim(-0.15, 0.3); axA.set_ylim(0, 15000); axA.legend(fontsize=8, markerscale=4)
axA.set_title("(a) window-mean $\\beta_{att}$ (covered telescope)")
axB.set_xlim(-0.09, 0.05); axB.set_ylim(0, 15000); axB.axvline(0, color="0.5", lw=0.8); axB.legend(fontsize=8)
axB.set_title("(b) zoom: the altitude-dependent offset and its ROBUST estimate")
axC.set_xlim(-0.02, 0.02); axC.set_ylim(0, 15000); axC.axvline(0, color="0.5", lw=0.8); axC.legend(fontsize=8)
axC.set_title("(c) non-range-corrected $\\beta_{att}/z^2$ + robust median")
for a_ in (axA, axB, axC):
    a_.grid(alpha=0.3); a_.set_ylabel("range [m]")
    a_.set_xlabel("[Mm$^{-1}$ sr$^{-1}$]")
axC.set_xlabel("[Mm$^{-1}$ sr$^{-1}$ km$^{-2}$]")
for j, (lo, hi) in enumerate(BANDS):
    axH[j].axvline(0, color="k", lw=0.9)
    axH[j].set_title(f"(d{j+1}) per-sample $\\beta_{{att}}/z^2$, {lo/1000:.0f}-{hi/1000:.0f} km")
    axH[j].set_xlabel("[Mm$^{-1}$ sr$^{-1}$ km$^{-2}$]"); axH[j].grid(alpha=0.3); axH[j].legend(fontsize=8)
fig.suptitle("Payerne CL61 hood-on dark measurements: robust offset profile and signal statistics",
             fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_dark_windows.png", dpi=150)
print("saved fig + robust cl61_b_dark.npz")
print("DARK_WINDOWS_DONE")
