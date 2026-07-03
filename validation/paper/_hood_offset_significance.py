# -*- coding: utf-8 -*-
"""Is the hood offset at the Rayleigh fit window STATISTICALLY significant, or an averaging
artefact? For each instrument, per-profile 3-5 km band-mean of P = rcs_0/z^2, then the median over
profiles with a bootstrap 95% CI and the median standard error. The decisive number is the offset
in units of its own uncertainty (sigma). Also the clear-night molecular reference the same way, so
the '-25%' fractional bias can be error-propagated. hervo63 challenge: CHM15k offset looks lost in
noise and the CHM15k calibration validates against EARLINET/CL61-cloud."""
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

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
# ident, C_L (to express P in beta/r^2), hood windows, band
INST = {
    "CHM15k": ("A", 5e11, [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
                            ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")]),
    "CL61":   ("C", 1.0,  [("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-26 11:45", "2026-05-27 13:15"),
                            ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-23 10:00", "2026-06-23 12:25")]),
}
NIGHTS = [("2026-04-02 20:00", "2026-04-03 03:30"), ("2026-04-22 20:00", "2026-04-23 03:30"),
          ("2026-04-24 20:00", "2026-04-25 03:30")]
BAND = (3000, 5000)
rng_boot = np.random.RandomState  # not used (Math.random banned in workflows, fine in plain python)

def load(ident, t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, rng = [], None
    for ds in sorted(days):
        f = L1 / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_{ident}{ds}.nc"
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
    return (np.vstack(X), rng) if X else (None, None)

def band_vals(ident, wins, CL):
    """Per-profile 3-5 km band-mean of beta/r^2 = (rcs_0/CL)/z^2, pooled over windows."""
    out = []
    for s1, s2 in wins:
        X, rng = load(ident, datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
        if X is None:
            continue
        zkm = rng / 1000.0
        P = (X / CL) / zkm[None, :] ** 2         # beta/r^2 in the MATLAB sense (per-profile)
        b = (rng >= BAND[0]) & (rng <= BAND[1])
        out.append(np.nanmean(P[:, b], axis=1))
    return np.concatenate(out) if out else np.array([])

def stats(v):
    v = v[np.isfinite(v)]
    n = v.size; med = np.median(v)
    se = 1.2533 * np.std(v) / np.sqrt(n)                    # standard error of the median
    boot = np.array([np.median(v[np.random.randint(0, n, n)]) for _ in range(3000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return n, med, se, lo, hi

print(f"3-5 km band-mean of beta/r^2 (per-profile), pooled 4 hood sessions")
print(f"{'inst':8s} {'N':>6s} {'median':>12s} {'SE(med)':>10s} {'sigma':>6s}   95% CI            molecular(clear)   FRAC BIAS")
res = {}
for inst, (ident, CL, wins) in INST.items():
    v = band_vals(ident, wins, CL)
    n, med, se, lo, hi = stats(v)
    # clear-night molecular, same band/units
    mv = band_vals(ident, NIGHTS, CL)
    mn, mmed, mse, mlo, mhi = stats(mv)
    sig = med / se
    fb = 100 * med / mmed
    fb_err = abs(fb) * np.sqrt((se / med) ** 2 + (mse / mmed) ** 2)
    res[inst] = (med, se, sig, mmed, mse, fb, fb_err)
    print(f"{inst:8s} {n:6d} {med:+12.3e} {se:10.2e} {sig:6.1f}   [{lo:+.2e},{hi:+.2e}]  "
          f"mol {mmed:+.2e}+-{mse:.1e}  {fb:+.0f}% +- {fb_err:.0f}%")
print()
print("Interpretation:")
for inst, (med, se, sig, mmed, mse, fb, fb_err) in res.items():
    verdict = ("SIGNIFICANT" if abs(sig) >= 3 else "MARGINAL" if abs(sig) >= 2 else "NOT SIGNIFICANT (lost in noise)")
    molsig = "molecular>3sigma" if abs(mmed / mse) >= 3 else "MOLECULAR ALSO NOISE-DOMINATED"
    print(f"  {inst}: offset {sig:+.1f}sigma -> {verdict};  reference {molsig} ({mmed/mse:+.1f}sigma)")
print("SIGNIF_DONE")
