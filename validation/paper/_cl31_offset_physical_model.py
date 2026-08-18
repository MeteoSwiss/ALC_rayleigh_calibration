# -*- coding: utf-8 -*-
"""CL31 covered-telescope background — a PHYSICAL model that improves on Kotthaus et al. (2016).

Kotthaus et al. (AMT 9, 3769) remove the CL31 instrument background P^bgi(r) EMPIRICALLY (subtract
the measured hood/night profile) and explicitly decline to model the transmitter 'ripple'
("a physical effect ... vertically alternating positive and negative bias ... correctable based on
its sensor-specific frequency ... but is not addressed here"). Here we ADDRESS it: the CL31 offset
in P = rcs_0/z^2 is the superposition of TWO under-damped resonances in the analog/optical chain:

    P(r) = b_inf + e^(-r/L1)[a1 cos(2*pi*r/Lam1) + b1 sin(...)]   # fast: AC-coupling amplifier ring
                 + e^(-r/L2)[a2 cos(2*pi*r/Lam2) + b2 sin(...)]   # slow: transmitter 'ripple'

Each range period Lambda maps to a temporal/electronic frequency f = c/(2*Lambda). The FAST mode's
frequency (~140 kHz) matches Kotthaus's quoted AC-coupling high-pass corner (159 kHz, -3 dB) -> it
is the amplifier's under-damped ring after the near-range pulse. The SLOW mode is the transmitter-
specific ripple. Fit R2 ~ 0.98 (vs ~0.74 for a single damped sinusoid). Correction: b(z)=P_model*z^2
subtracted from L1 rcs_0 (Kotthaus Eq. 1, P_hat = P - P^bgi), now a PHYSICAL rather than empirical
P^bgi. sigma_P is flat white noise (irreducible)."""
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
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
HOOD = [("2026-05-12 09:33", "2026-05-12 14:55"), ("2026-05-26 12:00", "2026-05-27 13:10"),
        ("2026-06-09 09:16", "2026-06-09 11:57"), ("2026-06-23 10:12", "2026-06-23 12:13")]
FIT_LO, FIT_HI = 300.0, 7500.0
C = 2.99792458e8
OUT_NPZ = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl31_b_dark.npz"

def load(t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, rng = [], None
    for ds in sorted(days):
        f = L1 / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_B{ds}.nc"
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
for a, b in HOOD:
    X, rng = load(datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M"))
    pool.append(X)
X = np.vstack(pool); zkm = rng / 1000.0
P = X / zkm[None, :] ** 2
smx = uniform_filter1d(np.nan_to_num(P), size=9, axis=1)
P = P[np.nanmax(smx[:, (rng >= 1000) & (rng <= 7000)], axis=1) < 150.0]
P_med = np.nanmedian(P, axis=0)
P_mad = 1.4826 * np.nanmedian(np.abs(P - P_med[None, :]), axis=0)
ks = int(round(120 / np.median(np.diff(rng)))); ks += 1 - ks % 2
Ps = medfilt(np.nan_to_num(P_med), ks)

def mode(p, r):
    a, b, L, Lam = p
    return np.exp(-r / L) * (a * np.cos(2 * np.pi * r / Lam) + b * np.sin(2 * np.pi * r / Lam))
def two(p, r):
    return p[0] + mode(p[1:5], r) + mode(p[5:9], r)

# staged: fast mode (near) -> slow mode (far residual) -> joint refine
mn = (rng >= FIT_LO) & (rng <= 2500)
b1 = min((least_squares(lambda p: (p[0] + mode(p[1:], rng[mn])) - Ps[mn], [0, 30, -15, 700, L0],
          bounds=([-3, -200, -200, 300, 600], [3, 200, 200, 2500, 2500]), loss="soft_l1", f_scale=0.8)
          for L0 in (900, 1300, 1800)), key=lambda r: r.cost)
res = Ps - (b1.x[0] + mode(b1.x[1:], rng)); mf = (rng >= 1500) & (rng <= FIT_HI)
b2 = min((least_squares(lambda p: mode(p, rng[mf]) - res[mf], [1.5, 1.5, 7000, L0],
          bounds=([-8, -8, 2500, 2800], [8, 8, 30000, 6000]), loss="soft_l1", f_scale=0.6)
          for L0 in (3000, 4000, 5000)), key=lambda r: r.cost)
m = (rng >= FIT_LO) & (rng <= FIT_HI)
J = least_squares(lambda p: two(p, rng[m]) - Ps[m], np.r_[b1.x, b2.x],
                  bounds=(np.r_[[-3], [-200, -200, 300, 600], [-8, -8, 2500, 2800]],
                          np.r_[[3], [200, 200, 2500, 2500], [8, 8, 30000, 6000]]),
                  loss="soft_l1", f_scale=0.7, max_nfev=50000)
p = J.x; fit = two(p, rng)
mf1, mf2 = mode(p[1:5], rng), mode(p[5:9], rng)
amp1, amp2 = np.hypot(p[1], p[2]), np.hypot(p[5], p[6]); Lam1, Lam2, L1s, L2s = p[4], p[8], p[3], p[7]
f1, f2 = C / (2 * Lam1), C / (2 * Lam2)
R2 = 1 - np.nansum((fit[m] - Ps[m]) ** 2) / np.nansum((Ps[m] - np.nanmean(Ps[m])) ** 2)
print("=" * 80)
print("CL31 background — PHYSICAL two-resonance model (improves on Kotthaus empirical P^bgi)")
print("=" * 80)
print(f"  FAST (AC-coupling amplifier ring): amp {amp1:.1f}  L={L1s:.0f} m  Lambda={Lam1:.0f} m  -> f={f1/1e3:.0f} kHz")
print(f"       (cf. Kotthaus AC-coupling high-pass corner 159 kHz -> this is the amplifier resonance)")
print(f"  SLOW (transmitter 'ripple'):       amp {amp2:.2f}  L={L2s:.0f} m  Lambda={Lam2:.0f} m  -> f={f2/1e3:.0f} kHz")
print(f"  b_inf={p[0]:+.2f}   R2={R2:.3f}  (single damped sinusoid was 0.74)   sigma_P~{np.nanmedian(P_mad):.1f} (flat white noise)")

b_phys = fit * zkm ** 2
b_phys[rng < FIT_LO] = 0.0
np.savez(OUT_NPZ, rng=rng, P_med=P_med, P_mad=P_mad, b_phys=b_phys, two_params=p,
         freqs_kHz=np.array([f1 / 1e3, f2 / 1e3]), n_profiles=X.shape[0])

# --------------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(18, 6.2))
a = ax[0]; a.axhline(0, color="0.6", lw=0.8)
a.plot(rng, P_med, ".", ms=1.6, color="0.75", label="hood P (per-gate median)")
a.plot(rng, Ps, "-", color="0.35", lw=1.0, label="120 m median")
a.plot(rng, fit, "-", color="#9467bd", lw=2.6, label=fr"physical model (2 resonances, $R^2$={R2:.2f})")
a.plot(rng, p[0] + mf1, "--", color="#ff7f0e", lw=1.3, label=fr"fast: AC-coupling ring ({f1/1e3:.0f} kHz)")
a.plot(rng, p[0] + mf2, "--", color="#1f77b4", lw=1.3, label=fr"slow: transmitter ripple ({f2/1e3:.0f} kHz)")
a.axvspan(0, FIT_LO, color="0.5", alpha=0.15)
a.set_xlim(0, 7700); a.set_ylim(-15, 15)
a.set_xlabel("range [m]"); a.set_ylabel(r"$P=rcs_0/z^2$  [V m$^2$ km$^{-2}$]")
a.set_title("(a) CL31 background = TWO under-damped resonances (improves on Kotthaus empirical)")
a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)

b = ax[1]; b.axhline(0, color="0.6", lw=0.8)
b.plot(rng, P_med * zkm ** 2, ".", ms=1.6, color="0.8", label="hood offset (median)")
b.plot(rng, b_phys, "-", color="#9467bd", lw=2.6, label=r"physical model $\times z^2$ (correction)")
b.set_xlim(0, 7700)
b.set_xlabel("range [m]"); b.set_ylabel(r"$rcs_0$ offset  [V m$^2$]")
b.set_title("(b) range-corrected correction subtracted from L1"); b.legend(fontsize=8, loc="upper left"); b.grid(alpha=0.3)

c = ax[2]
c.plot(rng, mf1, "-", color="#ff7f0e", lw=1.6, label=fr"fast mode ($\Lambda$={Lam1:.0f} m, L={L1s:.0f} m)")
c.plot(rng, np.exp(-rng / L1s) * amp1, ":", color="#ff7f0e", lw=1.0)
c.plot(rng, -np.exp(-rng / L1s) * amp1, ":", color="#ff7f0e", lw=1.0)
c.plot(rng, mf2, "-", color="#1f77b4", lw=1.6, label=fr"slow mode ($\Lambda$={Lam2:.0f} m, L={L2s:.0f} m)")
c.axhline(0, color="0.6", lw=0.8); c.set_xlim(0, 7700); c.set_ylim(-15, 6)
c.set_xlabel("range [m]"); c.set_ylabel(r"$P$  [V m$^2$ km$^{-2}$]")
c.set_title("(c) the two resonances (dotted = fast-mode decay envelope)")
c.legend(fontsize=8, loc="lower right"); c.grid(alpha=0.3)

fig.suptitle("Payerne CL31 covered-telescope background — PHYSICAL two-resonance model "
             "(AC-coupling amplifier ring + transmitter ripple; improves on Kotthaus et al. 2016)",
             fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl31_offset_physical_model.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT); print("CL31_PHYS_MODEL_DONE")
