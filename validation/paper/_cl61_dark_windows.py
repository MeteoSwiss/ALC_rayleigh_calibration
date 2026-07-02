"""Real hood-on dark windows for the Payerne CL61 (dates from dark_measurement_cl61_chm_cl31.m).
If the operational unit was covered, the L1 archive during these windows contains PURE dark
signal: its MEAN profile is the direct measurement of the baseline b_dark(z) in beta_att units,
and b_dark(z)/z^2 shows whether the artefact is flat in the NON-range-corrected signal
(electronic offset) or in beta (processing residual)."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
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
           ("2026-06-09 09:20", "2026-06-09 11:50"),
           ("2026-05-12 01:35", "2026-05-12 09:35")]   # 4th: pre-cover part of the long window


def load_window(t1, t2):
    days = {t1.strftime("%Y%m%d")}
    d = t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    profs, rng = [], None
    for ds in sorted(days):
        f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            tu = getattr(nc.variables["time"], "units", "days since 1970-01-01")
            tv = np.asarray(nc.variables["time"][:], "f8")
            base = datetime(1970, 1, 1)
            tt = np.array([base + timedelta(seconds=x) for x in (tv * 86400.0 if "day" in tu else tv)])
            r = np.asarray(nc.variables["range"][:], "f8")
            rcs = np.asarray(nc.variables["rcs_0"][:], "f8")
            if rcs.shape[0] != tt.size:
                rcs = rcs.T
        sel = (tt >= t1) & (tt <= t2)
        if sel.any():
            profs.append(rcs[sel]); rng = r
    if not profs:
        return None, None
    return np.vstack(profs), rng


fig, ax = plt.subplots(1, 3, figsize=(17, 6))
print("window                                  N     mean beta_att [Mm^-1 sr^-1] at")
print("                                              1 km      3 km      5 km      8 km     12 km    max(0-3km)")
results = []
for i, (a, b) in enumerate(WINDOWS):
    t1 = datetime.strptime(a, "%Y-%m-%d %H:%M"); t2 = datetime.strptime(b, "%Y-%m-%d %H:%M")
    X, rng = load_window(t1, t2)
    if X is None:
        print(f"{a} .. {b}: no L1 data"); continue
    mean = np.nanmean(X, axis=0) * 1e6
    def at(z):
        return float(mean[np.argmin(np.abs(rng - z))])
    mx = float(np.nanmax(mean[rng <= 3000]))
    # covered telescope: the MEAN at 1 km is ~0 (real daytime atmosphere: 0.1-0.5); spike maxima
    # (electronics/cross-talk) do not defeat the mean criterion
    dark = abs(at(1000)) < 0.05
    results.append((a, b, X.shape[0], mean, rng, dark))
    print(f"{a}..{b.split(' ')[1]}  {X.shape[0]:5d}  {at(1000):+8.4f}  {at(3000):+8.4f}  {at(5000):+8.4f}  "
          f"{at(8000):+8.4f}  {at(12000):+8.4f}  {mx:8.3f}  {'<- DARK' if dark else '<- atmosphere'}")
    lab = f"{a[5:16]} (n={X.shape[0]})"
    ax[0].plot(mean, rng, lw=1.4, label=lab)
    ax[1].plot(mean, rng, lw=1.4)
    # non-range-corrected shape: P(z) = beta/z^2 (arbitrary scale) — flat P = electronic offset
    P = mean / (rng / 1000.0) ** 2
    ax[2].plot(P, rng, lw=1.4)
ax[0].set_xlim(-0.1, 1.0); ax[0].set_ylim(0, 15000); ax[0].set_title("(a) window-mean $\\beta_{att}$")
ax[1].set_xlim(-0.08, 0.08); ax[1].set_ylim(0, 15000); ax[1].set_title("(b) zoom around zero")
ax[1].axvline(0, color="k", lw=0.8)
ax[2].set_xlim(-0.05, 0.05); ax[2].set_ylim(0, 15000)
ax[2].set_title("(c) $\\beta_{att}/z^2$ (non-range-corrected shape)"); ax[2].axvline(0, color="k", lw=0.8)
for a_ in ax:
    a_.grid(alpha=0.3); a_.set_ylabel("range [m]")
    a_.set_xlabel("[Mm$^{-1}$ sr$^{-1}$]" if a_ is not ax[2] else "[Mm$^{-1}$ sr$^{-1}$ km$^{-2}$]")
ax[0].legend(fontsize=8)
fig.suptitle("Payerne CL61 — L1 archive during the hood-on dark windows", fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_dark_windows.png")
fig.savefig(out, dpi=150)
# save the mean dark profile of the DARK windows for the correction experiment
darks = [m for (_, _, _, m, r, d) in results if d]
if darks:
    b_dark = np.nanmedian(np.vstack(darks), axis=0)
    np.savez(Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz"),
             rng=results[0][4], b_dark=b_dark)
    print(f"\nsaved median dark profile from {len(darks)} covered window(s) -> cl61_b_dark.npz")
print("saved", out)
print("DARK_WINDOWS_DONE")
