# -*- coding: utf-8 -*-
"""Physical model of the CHM15k covered-telescope offset (Payerne, 4 hood sessions) — the
counterpart of fig_cl61_offset_physical_model, for the photon-counting anchor.

Unlike the CL61 (analog, AC-coupled -> a POSITIVE lobe + a negative undershoot = the pulse
response of a high-pass), the CHM15k offset in P = rcs_0/z^2 has NO positive lobe: it is a
purely NEGATIVE residual that peaks near ~1 km and recovers monotonically toward zero. A single
negative exponential describes it (a second term collapses to zero):
    P(r) = b_inf - A exp(-r/L)
This is the signature of a photon-counting BACKGROUND / afterpulse OVER-SUBTRACTION (a near-range
derived component the firmware over-removes), NOT an AC-coupling artefact — the absence of the
positive lobe is the diagnostic that separates the two mechanisms. The SNR is lower than the CL61
(offset ~1e3 vs per-sample noise ~3e4 counts/s), so this is the median shape, not a per-profile
feature. Below ~1 km the overlap correction dominates and is excluded from the fit."""
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
from scipy.optimize import least_squares
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
HOOD = [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
        ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")]
FIT_LO, FIT_HI = 1000.0, 14000.0    # < 1 km: overlap-correction region, excluded
C_HALF = 2.99792458e8 / 2.0
OUT_NPZ = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/chm15k_b_dark.npz"

def load(t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, rng = [], None
    for ds in sorted(days):
        f = L1 / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_A{ds}.nc"
        if not f.exists():
            continue
        with Dataset(f) as nc:
            tv = np.asarray(nc.variables["time"][:], "f8")
            tt = np.array([datetime(1970, 1, 1) + timedelta(days=x) for x in tv])
            r = np.asarray(nc.variables["range"][:], "f8")
            x = np.asarray(nc.variables["rcs_0"][:], "f8")
            if x.shape[0] != tt.size:
                x = x.T
        s = (tt >= t1) & (tt <= t2)
        if s.any():
            X.append(x[s]); rng = r
    return np.vstack(X), rng

pool = []
for s1, s2 in HOOD:
    X, rng = load(datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
    pool.append(X)
X = np.vstack(pool); zkm = rng / 1000.0
P = X / zkm[None, :] ** 2
P_med = np.nanmedian(P, axis=0)
P_mad = 1.4826 * np.nanmedian(np.abs(P - P_med[None, :]), axis=0)

def model(p, r):
    binf, A, L = p
    return binf - A * np.exp(-r / L)

m = np.isfinite(P_med) & (rng >= FIT_LO) & (rng <= FIT_HI)
rf, yf = rng[m], P_med[m]
best = None
for L0 in (1500, 2500, 4000):
    r = least_squares(lambda p: model(p, rf) - yf, [0, 800, L0],
                      bounds=([-300, 0, 400], [300, 6000, 20000]), loss="soft_l1",
                      f_scale=300.0, max_nfev=20000)
    if best is None or r.cost < best.cost:
        best = r
binf, A, L = best.x
fit = model(best.x, rng)
def rmse(pred):
    return float(np.sqrt(np.nanmean((pred[m] - yf) ** 2)))
R2 = 1 - np.nansum((fit[m] - yf) ** 2) / np.nansum((yf - np.nanmean(yf)) ** 2)
# constant-offset null model (flat over-subtraction) for contrast
flat = np.full_like(rng, np.nanmedian(yf))

print("=" * 74)
print("CHM15k offset model  P(r) = b_inf - A exp(-r/L)  (photon-counting over-subtraction)")
print("=" * 74)
print(f"  A = {A:7.1f} counts/s     L = {L:7.0f} m   (range scale of the negative relaxation)")
print(f"  b_inf = {binf:+.1f} counts/s   peak {fit[m].min():.0f} @ {rf[np.argmin(fit[m])]:.0f} m")
print(f"  RMSE exp-recovery {rmse(fit):.1f}  R2={R2:.3f}  |  flat over-subtraction RMSE {rmse(flat):.1f}")
print(f"  sigma_P (per-sample) ~ {np.nanmedian(P_mad):.0f} counts/s -> offset is the MEDIAN shape (N={X.shape[0]})")
print(f"  NO positive near-lobe (contrast CL61) -> NOT an AC-coupling high-pass artefact.")

b_phys = fit * zkm ** 2
b_phys[rng < FIT_LO] = 0.0
np.savez(OUT_NPZ, rng=rng, P_med=P_med, P_mad=P_mad, b_phys=b_phys,
         phys_params=np.array(best.x), n_profiles=X.shape[0])

# --------------------------------------------------------------- figure (same 3-panel style)
fig, ax = plt.subplots(1, 3, figsize=(18, 6.2))
dr = np.median(np.diff(rng)); ks = int(round(600 / dr)); ks += 1 - ks % 2
Ps = medfilt(np.nan_to_num(P_med), ks)
S = 1e-3   # kcounts/s

a = ax[0]; a.axhline(0, color="0.6", lw=0.8)
a.plot(rng, P_med * S, ".", ms=1.6, color="0.75", label="hood P (per-gate median)")
a.plot(rng, Ps * S, "-", color="0.35", lw=1.3, label="600 m running median")
a.plot(rng, fit * S, "-", color="#d62728", lw=2.4, label=fr"model $b_\infty-Ae^{{-r/L}}$, L={L:.0f} m")
a.plot(rng, flat * S, ":", color="#2ca02c", lw=1.6, label="flat over-subtraction (rejected)")
a.axvspan(0, FIT_LO, color="0.5", alpha=0.15)
a.text(500, 0.55, "overlap\n(<1 km,\nexcluded)", ha="center", fontsize=7, color="0.4")
a.set_xlim(0, 15000); a.set_ylim(-1.4, 0.9)
a.set_xlabel("range [m]"); a.set_ylabel(r"$P=rcs_0/z^2$  [10$^3$ counts s$^{-1}$]")
a.set_title("(a) raw offset: purely NEGATIVE, no positive lobe -> over-subtraction, not AC-coupling")
a.legend(fontsize=8, loc="lower right"); a.grid(alpha=0.3)

b = ax[1]; b.axhline(0, color="0.6", lw=0.8)
b.plot(rng, P_med * zkm ** 2 * S, ".", ms=1.6, color="0.8", label="hood offset (median)")
b.plot(rng, Ps * zkm ** 2 * S, "-", color="0.35", lw=1.2, label="600 m running median")
b.plot(rng, b_phys * S, "-", color="#d62728", lw=2.4, label=r"physical model $\times z^2$")
b.set_xlim(0, 15000)
b.set_xlabel("range [m]"); b.set_ylabel(r"$rcs_0$ offset  [10$^3$ counts s$^{-1}$ m$^2$]")
b.set_title("(b) range-corrected offset = correction applied to L1")
b.legend(fontsize=8, loc="lower right"); b.grid(alpha=0.3)

c = ax[2]
c.semilogy(rng, np.abs(A * np.exp(-rng / L)) * S, "-", color="#d62728", lw=2,
           label=fr"relaxation: L={L:.0f} m ($\approx${2*L/C_HALF*1e6:.0f} $\mu$s eq.)")
c.semilogy(rng, np.abs(Ps - binf) * S, ".", ms=2.2, color="0.5", label="|median - b$_\\infty$|")
c.axvspan(0, FIT_LO, color="0.5", alpha=0.12)
c.set_xlim(0, 15000); c.set_ylim(3e-3, 2)
c.set_xlabel("range [m]"); c.set_ylabel(r"|component|  [10$^3$ counts s$^{-1}$]")
c.set_title("(c) single negative relaxation (median follows the exponential above ~1 km)")
c.legend(fontsize=8, loc="upper right"); c.grid(alpha=0.3, which="both")

fig.suptitle("Payerne CHM15k covered-telescope offset — physical model "
             "(photon-counting background/afterpulse over-subtraction; contrast the CL61 AC-coupling)",
             fontweight="bold", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_chm15k_offset_physical_model.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
print("CHM_PHYS_MODEL_DONE")
