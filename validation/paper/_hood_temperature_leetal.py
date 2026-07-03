# -*- coding: utf-8 -*-
"""Temperature analysis in the style of Le et al. (2026), fig. panels (f)/(g): for each instrument
the INSTRUMENTAL BIAS mu(beta/r^2) and the NOISE variance sigma^2(beta/r^2) as a function of range,
one curve per internal-temperature bin (coloured by temp_int). Pooled over the four 2026 hood
sessions. beta/r^2 = (rcs_0/C_L)/z^2 (CL61 C_L=1; CHM15k C_L=5e11). mu is the per-gate MEAN over
the profiles in the bin (the DC bias); sigma^2 the per-gate variance (the detector noise)."""
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

L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
K2C = 273.15
INST = {"CL61": ("C", 1.0, [("2026-05-12 09:17", "2026-05-12 14:58"), ("2026-05-26 11:45", "2026-05-27 13:15"),
                             ("2026-06-09 09:10", "2026-06-09 11:51"), ("2026-06-23 10:00", "2026-06-23 12:25")]),
        "chm15k": ("A", 5e11, [("2026-05-12 09:24", "2026-05-12 14:53"), ("2026-05-26 12:00", "2026-05-27 13:15"),
                               ("2026-06-09 09:16", "2026-06-09 11:55"), ("2026-06-23 12:32", "2026-06-23 15:09")])}

def load(ident, t1, t2):
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    X, Ti, rng = [], [], None
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
            ti = np.asarray(nc.variables["temp_int"][:], "f8").ravel() - K2C
        s = (tt >= t1) & (tt <= t2)
        if s.any():
            X.append(x[s]); Ti.append(ti[s]); rng = r
    return (np.vstack(X), np.concatenate(Ti), rng) if X else (None, None, None)

fig, ax = plt.subplots(2, 2, figsize=(15, 10))
NB = 6
for row, (inst, (ident, CL, wins)) in enumerate(INST.items()):
    Xs, Ts, rng = [], [], None
    for s1, s2 in wins:
        X, Ti, rng = load(ident, datetime.strptime(s1, "%Y-%m-%d %H:%M"), datetime.strptime(s2, "%Y-%m-%d %H:%M"))
        if X is not None:
            Xs.append(X); Ts.append(Ti)
    X = np.vstack(Xs); Ti = np.concatenate(Ts); zkm = rng / 1000.0
    # contamination screen: reject profiles with a COHERENT real return (hood on/off transitions,
    # cloud/aerosol leaking in before full coverage). Two criteria, both robust to photon-counting
    # single-gate spikes via a 5-gate boxcar: (i) a persistent haze -> median(1.5-4 km) > 0.1;
    # (ii) a LOCALISED cloud anywhere 0.8-4 km -> smoothed peak > 5 (dark noise is < 2, clouds >> 5;
    # this catches the 9 Jun 09:10 hood-on transition, a 412 Mm-1sr-1 cloud at 1.35 km that the
    # band-median missed because the layer sits below the 1.5 km band).
    from scipy.ndimage import uniform_filter1d
    betaMm = (X / CL) * 1e6
    sm = uniform_filter1d(np.nan_to_num(betaMm), size=5, axis=1)
    keep = (np.nanmedian(betaMm[:, (rng >= 1500) & (rng <= 4000)], axis=1) < 0.1) \
        & (np.nanmax(sm[:, (rng >= 800) & (rng <= 4000)], axis=1) < 5.0)
    nrej = (~keep).sum()
    X, Ti = X[keep], Ti[keep]
    print(f"  {inst}: screened {nrej} contaminated profiles ({100*nrej/keep.size:.1f}%)")
    P = (X / CL) / zkm[None, :] ** 2                       # beta/r^2 (per-profile), a.u.
    edges = np.quantile(Ti, np.linspace(0, 1, NB + 1))
    norm = plt.Normalize(np.median(Ti[Ti <= edges[1]]), np.median(Ti[Ti >= edges[-2]]))
    dr = np.median(np.diff(rng)); kk = int(round(200 / dr)); kk += 1 - kk % 2
    amu, asg = ax[row, 0], ax[row, 1]
    for i in range(NB):
        sel = (Ti >= edges[i]) & (Ti < edges[i + 1] if i < NB - 1 else Ti <= edges[i + 1])
        if sel.sum() < 100:
            continue
        Tc = float(np.median(Ti[sel])); c = cm.coolwarm(norm(Tc))
        mu = np.nanmean(P[sel], axis=0)                    # instrumental bias
        sg = np.nanvar(P[sel], axis=0)                     # noise variance
        amu.plot(medfilt(np.nan_to_num(mu), kk), rng, lw=1.4, color=c)
        asg.plot(medfilt(np.nan_to_num(sg), kk), rng, lw=1.4, color=c)
    sm = cm.ScalarMappable(norm=norm, cmap="coolwarm")
    for a in (amu, asg):
        a.set_ylim(0, 15000); a.grid(alpha=0.3); a.set_ylabel("Range [m]")
        fig.colorbar(sm, ax=a).set_label("temp_int [°C]")
    amu.axvline(0, color="0.5", lw=0.8); amu.axhspan(3000, 5000, color="orange", alpha=0.10)
    # zoom mu to the free-troposphere structure (near-field spike goes off-scale)
    mu_all = np.array([np.nanmean(P[(Ti >= edges[i]) & (Ti < edges[i + 1])][:, (rng >= 2000) & (rng <= 12000)])
                       for i in range(NB)])
    lim = 4 * np.nanmax(np.abs([np.nanmean(P[(Ti >= edges[i]) & (Ti < edges[i + 1])], axis=0)[(rng >= 2500) & (rng <= 6000)]
                                for i in range(NB)]))
    amu.set_xlim(-lim, lim)
    # zoom sigma^2 to the flat homoscedastic level (near-field variance goes off-scale)
    sg_all = np.nanvar(P, axis=0)
    ref = np.nanmedian(sg_all[(rng >= 5000) & (rng <= 12000)])
    asg.set_xlim(0.3 * ref, 2.2 * ref); asg.axhspan(3000, 5000, color="orange", alpha=0.10)
    amu.set_xlabel(r"$\mu_{\beta/r^2}$  (instrumental bias) [a.u.]")
    amu.set_title(f"({'f' if row == 0 else 'h'}) {inst}: bias  $\\mu(\\beta/r^2)$ vs range, per temp bin")
    asg.set_xlabel(r"$\sigma^2_{\beta/r^2}$  (noise) [a.u.]")
    asg.set_title(f"({'g' if row == 0 else 'i'}) {inst}: noise  $\\sigma^2(\\beta/r^2)$ vs range, per temp bin")

fig.suptitle("Payerne hood dark: instrumental bias (μ) and noise (σ²) vs range, coloured by internal "
             "temperature — Le et al. (2026) style\nCL61 (top): bias grows toward the surface & deepens when cold; "
             "CHM15k (bottom): weaker bias. σ² is ~flat (homoscedastic) and rises when warm (more dark current).",
             fontweight="bold", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_hood_temperature_leetal.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
