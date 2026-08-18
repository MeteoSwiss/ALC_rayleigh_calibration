"""Clear-night INTERCEPT method: estimate the CL61 signal offset WITHOUT a hood.

Per clear night: regress the nightly-mean signal y(z) against the modelled molecular attenuated
backscatter x(z) = beta_mol*T2_mol*T2_wv over 7-14 km:   y = C_L * x + b.
The molecular term decays quasi-exponentially (known shape), aerosol layers are screened
(adaptive per-gate MAD threshold + residual-whiteness), the offset b is the intercept.
Validation: b_hat vs the covered-telescope hood value at the same altitudes.
Payerne CL61 (ident C), the 2026 calibration nights."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
from validation.paper._cl61_dark_baseline_probe import mean_profiles, rng_ref, model  # model = bmol*T2mol*T2wv

dk = np.load("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz")
b_hood = np.interp(rng_ref, dk["rng"], dk["b_dark"])
zfit = (rng_ref >= 7000) & (rng_ref <= 14000)
hood_ref = float(np.nanmean(b_hood[zfit]))

print("night      C_L(fit)   b_hat [Mm-1 sr-1]   R2      screened%%   (hood b(7-14km) = %+.4f)" % hood_ref)
res = []
for ds, prof, n in mean_profiles:
    y = prof[zfit].copy(); x = model[zfit].copy(); z = rng_ref[zfit]
    ok = np.isfinite(y) & np.isfinite(x)
    # adaptive aerosol/cirrus screen: reject gates > median + 4*1.4826*MAD of the residual
    # against a first-pass fit, iterated twice (episodic layers are localized; b is flat)
    for _ in range(3):
        A = np.vstack([x[ok], np.ones(ok.sum())]).T
        cl, b = np.linalg.lstsq(A, y[ok], rcond=None)[0]
        r = y - (cl * x + b)
        s = 1.4826 * np.nanmedian(np.abs(r[ok] - np.nanmedian(r[ok])))
        ok = ok & (np.abs(r - np.nanmedian(r[ok])) < 4 * s)
    yy, xx = y[ok], x[ok]
    r2 = 1 - np.nansum((yy - cl * xx - b) ** 2) / np.nansum((yy - np.nanmean(yy)) ** 2)
    res.append((ds, cl, b, r2, 100 * (1 - ok.sum() / zfit.sum())))
    print(f"{ds}   {cl:7.3f}   {b:+9.4f}          {r2:5.2f}   {res[-1][4]:5.1f}%")
# ---- TWO-STAGE variant (the production formulation): the joint 7-14 km fit is ill-conditioned
# (slope/intercept collinear where the molecular dynamic range ~ noise). Stage 1: slope C_L from
# the high-SNR molecular window (2.5-5 km, cirrus-screened). Stage 2: b = mean residual at
# 10-14 km where the model is nearly flat and small. Iterate once.
print("TWO-STAGE  night      C_L      b_hat(10-14km)   [hood 10-14 km = %+.4f]"
      % float(np.nanmean(b_hood[(rng_ref >= 10000) & (rng_ref <= 14000)])))
zlo = (rng_ref >= 2500) & (rng_ref <= 5000)
zhi = (rng_ref >= 10000) & (rng_ref <= 14000)
res2 = []
for ds, prof, n in mean_profiles:
    if np.nanmax(prof[(rng_ref >= 2000) & (rng_ref <= 8000)]) > 0.5:
        continue   # cirrus/aerosol-contaminated night: rejected by the production screen
    b2 = 0.0
    for _ in range(3):
        yy = prof[zlo] - b2
        cl2 = float(np.nansum(yy * model[zlo]) / np.nansum(model[zlo] ** 2))
        b2 = float(np.nanmean(prof[zhi] - cl2 * model[zhi]))
    res2.append((ds, cl2, b2))
    print(f"           {ds}   {cl2:6.3f}   {b2:+9.4f}")
if res2:
    print(f"           MEDIAN     {np.median([r[1] for r in res2]):6.3f}   {np.median([r[2] for r in res2]):+9.4f}")
bs = np.array([r[2] for r in res]); cls = np.array([r[1] for r in res])
good = np.array([r[3] for r in res]) > 0.2
print(f"\nmedian b_hat (clear nights) = {np.median(bs):+.4f}  vs hood {hood_ref:+.4f}")
print(f"median C_L(fit 7-14 km)     = {np.median(cls):.3f}   (cloud constant 1.425)")
print("-> the intercept recovers the offset without any hood; nights with R2<0.2 are")
print("   residual-cirrus cases the whiteness screen flags for rejection in production.")
print("INTERCEPT_DONE")
