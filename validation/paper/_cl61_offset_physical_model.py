# -*- coding: utf-8 -*-
"""Physical model of the CL61 shot-synchronous electronic offset (hood-on, Payerne).

The offset is measured in NON-range-corrected space P(r) = beta_att / z^2 (homoscedastic
detector noise), pooled over the three covered-telescope windows (cl61_b_dark.npz -> P_med).

Observed shape of P(r):
  * strong POSITIVE near-range lobe (fast decay, 0-1.7 km),
  * zero crossing ~1.7 km, NEGATIVE undershoot minimum ~2.5 km,
  * slow (quasi-linear) recovery back toward zero through 3-10 km,
  * flat ~0 baseline beyond ~11 km.

This is the impulse/pulse response of an AC-coupled (high-pass, "CR-differentiator") analog
chain to the shot-synchronous near-range transient (internal reflection + electrical pickup +
detector recovery). A high-pass blocks DC, so its response to a fast input pulse is the pulse
itself MINUS a slow undershoot whose area exactly cancels it (zero net area = DC blocked):

    P(r) = A_p * exp(-r/L_p)        (fast POSITIVE lobe: near-range transient / afterpulse)
         - A_u * exp(-r/L_u)        (slow NEGATIVE undershoot: AC-coupling recovery, tau_u = RC)
         + b_inf                    (residual DC floor, ~0)

with A_p, A_u >= 0 and L_p << L_u. Range<->time<->RC: r = c*t/2, so a scale L maps to a
time constant tau = 2L/c = L / (c/2), c/2 = 1.499e8 m/s. The DC-blocking (AC) constraint is
testable: Integral(P dr) = A_p*L_p - A_u*L_u  -> 0 for a pure high-pass.

The slow undershoot term, Taylor-expanded for r << L_u, is  -A_u + (A_u/L_u) r  -- i.e. the
earlier empirical LINEAR RAMP (linfit_a ~ A_u/L_u, linfit_c ~ b_inf - A_u). So the linear-ramp
correction used so far is exactly the small-argument approximation of this physical model.
"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NPZ = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz"
C_HALF = 2.99792458e8 / 2.0            # m/s -> range/time mapping
FIT_LO, FIT_HI = 350.0, 15000.0        # m (below 350 m: hood near-field, clipped)

d = np.load(NPZ)
rng = d["rng"].astype(float)
P = d["P_med"].astype(float)           # raw-space per-gate median offset [Mm-1 sr-1 km-2]
Pmad = d["P_mad"].astype(float)
nprof = int(d["n_profiles"])
zkm = rng / 1000.0

# ---- physical model: high-pass pulse response ----
def model(p, r):
    Ap, Lp, Au, Lu, binf = p
    return Ap * np.exp(-r / Lp) - Au * np.exp(-r / Lu) + binf

m = np.isfinite(P) & (rng >= FIT_LO) & (rng <= FIT_HI)
rf, Pf = rng[m], P[m]

def resid(p):
    return model(p, rf) - Pf

# multi-start robust NLS (soft_l1 so the near-range spike does not dominate)
lb = [0.0,   80.0,  0.0,  2000.0, -1e-4]
ub = [1.0,  1500.0, 1.0, 15000.0,  1e-4]
best = None
for Lp0 in (130, 200, 350):
    for Lu0 in (3000, 5000, 8000):
        p0 = [0.02, Lp0, 1.5e-3, Lu0, 0.0]
        try:
            r = least_squares(resid, p0, bounds=(lb, ub), loss="soft_l1",
                              f_scale=5e-4, max_nfev=20000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
Ap, Lp, Au, Lu, binf = best.x
fit = model(best.x, rng)

# goodness vs the two simpler models, on the fit region
def rmse(pred):
    return float(np.sqrt(np.nanmean((pred[m] - Pf) ** 2)))
ss_res = np.nansum((fit[m] - Pf) ** 2)
ss_tot = np.nansum((Pf - np.nanmean(Pf)) ** 2)
R2 = 1 - ss_res / ss_tot
# linear ramp (stored)
lin = float(d["linfit_a"]) * rng + float(d["linfit_c"])
# single negative exponential (undershoot only, no positive lobe)
def resid1(p):
    Au1, Lu1, b1 = p
    return (-Au1 * np.exp(-rf / p[1]) + b1) - Pf
r1 = least_squares(resid1, [1.5e-3, 5000, 0.0],
                   bounds=([0, 2000, -1e-3], [1, 15000, 1e-3]), loss="soft_l1", f_scale=5e-4)
sing = -r1.x[0] * np.exp(-rng / r1.x[1]) + r1.x[2]

# time constants + AC (DC-blocking) area test
tau_p, tau_u = Lp / C_HALF * 1e6, Lu / C_HALF * 1e6      # microseconds
area_pos, area_neg = Ap * Lp, Au * Lu                     # integral of each exp term
xcross = -Lu * np.log((binf) / Au) if (0 < binf < Au) else np.nan

print("=" * 74)
print("PHYSICAL MODEL  P(r) = Ap*exp(-r/Lp) - Au*exp(-r/Lu) + b_inf   (AC high-pass)")
print("=" * 74)
print(f"  A_p (positive lobe)  = {Ap:+.4e}  Mm-1 sr-1 km-2")
print(f"  L_p (fast scale)     = {Lp:8.1f} m   -> tau_p = {tau_p:6.2f} us   (near-range transient)")
print(f"  A_u (undershoot)     = {Au:+.4e}  Mm-1 sr-1 km-2")
print(f"  L_u (slow scale)     = {Lu:8.1f} m   -> tau_u = {tau_u:6.2f} us   (AC-coupling RC)")
print(f"  b_inf (DC floor)     = {binf:+.3e}")
print(f"  undershoot minimum   ~ {rng[np.argmin(fit)]:.0f} m,  P_min = {fit.min():+.2e}")
print("-" * 74)
print(f"  DC-blocking (AC) test:  A_p*L_p = {area_pos:.3e}   A_u*L_u = {area_neg:.3e}"
      f"   ratio = {area_pos/area_neg:5.2f}   (=1 for a pure high-pass)")
print(f"  slow-term linearisation: A_u/L_u = {Au/Lu:.3e}  vs stored linfit_a = {float(d['linfit_a']):.3e}")
print("-" * 74)
print(f"  RMSE  physical(5p) = {rmse(fit):.3e}   R2 = {R2:.3f}")
print(f"  RMSE  single-exp   = {rmse(sing):.3e}")
print(f"  RMSE  linear ramp  = {rmse(lin):.3e}")

# save physical-model offset for recalibration (beta-space, zero the hood near-field)
b_phys = fit * zkm ** 2
b_phys[rng < FIT_LO] = 0.0
out = {k: d[k] for k in d.files}
out.update(dict(b_phys=b_phys, phys_params=np.array(best.x),
                phys_names=np.array(["Ap", "Lp", "Au", "Lu", "b_inf"])))
np.savez(NPZ, **out)
print("saved b_phys + phys_params into", NPZ.split('/')[-1])

# ----------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(18, 6.2))
dr = np.median(np.diff(rng)); ks = int(round(600 / dr)); ks += 1 - ks % 2
Ps = medfilt(np.nan_to_num(P), ks)

# (a) raw space P(r): data + model + components + regimes
a = ax[0]
a.axhline(0, color="0.6", lw=0.8)
a.plot(rng, P * 1e3, ".", ms=1.6, color="0.75", label="hood P (per-gate median)")
a.plot(rng, Ps * 1e3, "-", color="0.35", lw=1.3, label="600 m running median")
a.plot(rng, fit * 1e3, "-", color="#d62728", lw=2.4, label="physical model (5 par)")
a.plot(rng, (Ap * np.exp(-rng / Lp) + binf) * 1e3, "--", color="#ff7f0e", lw=1.5,
       label=r"fast +lobe  $A_p e^{-r/L_p}$")
a.plot(rng, (-Au * np.exp(-rng / Lu) + binf) * 1e3, "--", color="#1f77b4", lw=1.5,
       label=r"slow undershoot  $-A_u e^{-r/L_u}$")
a.axvspan(FIT_LO, 1700, color="#ff7f0e", alpha=0.06)
a.axvspan(1700, 10500, color="#1f77b4", alpha=0.06)
a.axvspan(10500, 15000, color="0.5", alpha=0.06)
a.set_xlim(0, 15000); a.set_ylim(-2.2, 4.5)
a.set_xlabel("range [m]"); a.set_ylabel(r"$P=\beta_{att}/z^2$  [Mm$^{-1}$sr$^{-1}$km$^{-2}$]  $\times10^{3}$")
a.set_title("(a) raw (non-range-corrected) offset: AC high-pass pulse response")
a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
a.text(850, -1.9, "fast\ndecay", ha="center", fontsize=8, color="#b5651d")
a.text(6000, -1.9, "slow (quasi-linear) undershoot recovery", ha="center", fontsize=8, color="#1f5fa8")
a.text(12700, -1.9, "flat", ha="center", fontsize=8, color="0.4")
a.annotate("internal pulse → +16\n(off-scale, <350 m,\nexcluded from fit)", xy=(300, 4.35),
           xytext=(1500, 3.3), fontsize=7, color="0.4",
           arrowprops=dict(arrowstyle="->", color="0.55", lw=0.8))

# (b) beta space (x z^2): the actual correction, model vs empirical linear ramp
b = ax[1]
b.axhline(0, color="0.6", lw=0.8)
b.plot(rng, d["b_dark"], ".", ms=1.6, color="0.8", label=r"hood $\beta_{att}$ (median)")
b.plot(rng, Ps * zkm ** 2, "-", color="0.35", lw=1.2, label="600 m running median")
b.plot(rng, b_phys, "-", color="#d62728", lw=2.4, label=r"physical model $\times z^2$")
b.plot(rng, np.minimum(lin, 0) * zkm ** 2 * ((rng >= 300) & (rng <= 12000)),
       ":", color="#2ca02c", lw=1.8, label=r"linear-ramp (clamped) $\times z^2$")
b.set_xlim(0, 15000); b.set_ylim(-0.06, 0.03)
b.set_xlabel("range [m]"); b.set_ylabel(r"$\beta_{att}$ offset  [Mm$^{-1}$sr$^{-1}$]")
b.set_title("(b) range-corrected offset = correction applied to L1")
b.legend(fontsize=8, loc="lower left"); b.grid(alpha=0.3)

# (c) two exponential relaxations: each component line vs (data - the OTHER component).
# Subtracting the slow undershoot from the median leaves the fast lobe (orange dots) -> it
# DOES follow the fast line; subtracting the fast lobe leaves the slow undershoot (blue dots).
c = ax[2]
fastc = Ap * np.exp(-rng / Lp); slowc = Au * np.exp(-rng / Lu)
c.semilogy(rng, fastc * 1e3, "-", color="#ff7f0e", lw=2,
           label=fr"fast lobe: $\tau_p$={tau_p:.1f} $\mu$s ($L_p$={Lp:.0f} m)")
c.semilogy(rng, np.abs(Ps + slowc - binf) * 1e3, ".", ms=2.6, color="#ff7f0e", alpha=0.45,
           label="median − slow  (isolates fast)")
c.semilogy(rng, slowc * 1e3, "-", color="#1f77b4", lw=2,
           label=fr"slow undershoot: $\tau_u$={tau_u:.1f} $\mu$s ($L_u$={Lu:.0f} m)")
c.semilogy(rng, np.abs(Ps - fastc - binf) * 1e3, ".", ms=2.6, color="#1f77b4", alpha=0.45,
           label="median − fast  (isolates slow)")
c.axvspan(0, 500, color="0.5", alpha=0.12)
c.text(250, 12, "internal\npulse", ha="center", fontsize=7, color="0.35")
c.set_xlim(0, 15000); c.set_ylim(1e-3, 30)
c.set_xlabel("range [m]"); c.set_ylabel(r"|value|  [Mm$^{-1}$sr$^{-1}$km$^{-2}$]  $\times10^{3}$")
c.set_title("(c) two relaxations: each followed by data minus the other component")
c.legend(fontsize=7.5, loc="upper right"); c.grid(alpha=0.3, which="both")

fig.suptitle("Payerne CL61 electronic offset — physical model (AC-coupling high-pass pulse response)",
             fontweight="bold", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_offset_physical_model.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
print("PHYS_MODEL_DONE")
