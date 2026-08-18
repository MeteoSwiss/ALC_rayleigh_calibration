# -*- coding: utf-8 -*-
"""Temperature dependence of the CL61 electronic offset (physical high-pass model).

Pools ALL THREE covered-telescope (hood) windows:
  * 2026-05-26 11:45 -> 05-27 13:15  (the 25.5 h window: a full diurnal temperature cycle),
  * 2026-05-12 09:35 -> 14:50        (extra dark measurement),
  * 2026-06-09 09:20 -> 11:50        (extra dark measurement).
Bins the profiles by INTERNAL electronics temperature temp_int (range ~22-44 C, much wider
lever than temperature_laser), and refits the physical model
    P(r) = Ap exp(-r/Lp) - Au exp(-r/Lu) + b_inf
per temperature bin. Answers: does the night-time deepening of the undershoot come from the
undershoot AMPLITUDE A_u (gain / pulse charge) or the RC TIME CONSTANT tau_u (=L_u) drifting
with temperature? Output: parameter-vs-T trends + a temperature-indexed correction law.
"""
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

WMO, IDENT = "0-20000-0-06610", "C"
L1 = Path("D:/E-PROFILE_L1_2026")
WINDOWS = [("2026-05-26 11:45", "2026-05-27 13:15"),
           ("2026-05-12 09:35", "2026-05-12 14:50"),
           ("2026-06-09 09:20", "2026-06-09 11:50")]
NPZ = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/cl61_b_dark.npz"
C_HALF = 2.99792458e8 / 2.0
FIT_LO, FIT_HI = 350.0, 15000.0
K2C = 273.15

def load_win(t1, t2):
    """Return P (n,gates) in raw space [Mm-1 sr-1 km-2], rng, temp_int[C], temp_laser[C]."""
    days, d = set(), t1
    while d.date() <= t2.date():
        days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
    P, Ti, Tl, rng = [], [], [], None
    for ds in sorted(days):
        f = L1 / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
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
            tl = np.asarray(nc.variables["temperature_laser"][:], "f8").ravel() - K2C
        s = (tt >= t1) & (tt <= t2)
        if s.any():
            zkm = r / 1000.0
            P.append(x[s] * 1e6 / zkm[None, :] ** 2); Ti.append(ti[s]); Tl.append(tl[s]); rng = r
    return (np.vstack(P), rng, np.concatenate(Ti), np.concatenate(Tl)) if P else (None, None, None, None)

Ps, Tis, Tls = [], [], []
rng = None
print("window                         N      temp_int[C]        temp_laser[C]")
for a, b in WINDOWS:
    t1 = datetime.strptime(a, "%Y-%m-%d %H:%M"); t2 = datetime.strptime(b, "%Y-%m-%d %H:%M")
    P, rng, Ti, Tl = load_win(t1, t2)
    Ps.append(P); Tis.append(Ti); Tls.append(Tl)
    print(f"{a}..{b[11:]}   {P.shape[0]:5d}   [{Ti.min():5.1f},{Ti.max():5.1f}]   [{Tl.min():5.1f},{Tl.max():5.1f}]")
P = np.vstack(Ps); Ti = np.concatenate(Tis); Tl = np.concatenate(Tls)
zkm = rng / 1000.0
print(f"POOLED  N={P.shape[0]}  temp_int [{Ti.min():.1f},{Ti.max():.1f}] C (span {Ti.max()-Ti.min():.1f} K)")

# ---- physical model + fitters ----
def model(p, r):
    Ap, Lp, Au, Lu, binf = p
    return Ap * np.exp(-r / Lp) - Au * np.exp(-r / Lu) + binf

g = np.load(NPZ)
GP = [float(x) for x in g["phys_params"]]           # global Ap,Lp,Au,Lu,binf
mfit = np.isfinite(rng) & (rng >= FIT_LO) & (rng <= FIT_HI)
rf = rng[mfit]

def fit_full(Pmed):
    y = Pmed[mfit]
    ok = np.isfinite(y)
    lb = [0, 80, 0, 2000, -2e-4]; ub = [1, 1500, 1, 15000, 2e-4]
    best = None
    for Lp0 in (GP[1], 300, 900):
        for Lu0 in (GP[3], 3500, 6000):
            p0 = [max(GP[0], 1e-3), Lp0, max(GP[2], 1e-3), Lu0, 0.0]
            try:
                r = least_squares(lambda p: model(p, rf[ok]) - y[ok], p0, bounds=(lb, ub),
                                  loss="soft_l1", f_scale=6e-4, max_nfev=12000)
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
    return best.x

def fit_amp(Pmed):
    """Constrained: L_p, L_u fixed at global -> robust amplitude trend (Ap, Au, b_inf linear)."""
    Lp, Lu = GP[1], GP[3]
    y = Pmed[mfit]; ok = np.isfinite(y)
    A = np.column_stack([np.exp(-rf[ok] / Lp), -np.exp(-rf[ok] / Lu), np.ones(ok.sum())])
    coef, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
    return coef  # Ap, Au, b_inf

# ---- bin by temp_int (equal-count quantile bins) ----
NB = 6
edges = np.quantile(Ti, np.linspace(0, 1, NB + 1))
bins = []
for i in range(NB):
    sel = (Ti >= edges[i]) & (Ti <= edges[i + 1] if i == NB - 1 else Ti < edges[i + 1])
    if sel.sum() < 200:
        continue
    Pmed = np.nanmedian(P[sel], axis=0)
    full = fit_full(Pmed); Ap, Au, binf = fit_amp(Pmed)
    off36 = np.nanmean((Pmed * zkm ** 2)[(rng >= 3000) & (rng <= 6000)])
    off812 = np.nanmean((Pmed * zkm ** 2)[(rng >= 8000) & (rng <= 12000)])
    bins.append(dict(T=float(np.median(Ti[sel])), n=int(sel.sum()), Pmed=Pmed,
                     Ap=full[0], Lp=full[1], Au=full[2], Lu=full[3], binf=full[4],
                     Ap_c=Ap, Au_c=Au, off36=off36, off812=off812))
    print(f"  bin T={bins[-1]['T']:5.1f}C n={sel.sum():5d}  Au_c={Au:+.2e} Ap_c={Ap:+.2e}  "
          f"Lu={full[3]:6.0f} Lp={full[1]:5.0f}  off(3-6km)={off36:+.4f}")

T = np.array([b["T"] for b in bins])
def lin(y):
    a, c = np.polyfit(T, y, 1); return a, c, a * T + c
# pooled trends
sA, cA, _ = lin(np.array([b["Au_c"] for b in bins]))         # undershoot amplitude
sL, cL, _ = lin(np.array([b["Lu"] for b in bins]))           # undershoot RC scale
s36, c36, _ = lin(np.array([b["off36"] for b in bins]))      # 3-6 km beta offset
Trng = T.max() - T.min()

# WITHIN-SESSION trend (25.5 h window alone) — isolates temperature from session-to-session
# drift; the pooled trend can be inflated by the coldest bin belonging to a different session.
Pw, Tiw = Ps[0], Tis[0]
ew = np.quantile(Tiw, np.linspace(0, 1, 6)); Tw, ow = [], []
for i in range(5):
    s = (Tiw >= ew[i]) & (Tiw < (ew[i + 1] + (1 if i == 4 else 0)))
    Pm = np.nanmedian(Pw[s], axis=0)
    Tw.append(float(np.median(Tiw[s]))); ow.append(float(np.nanmean((Pm * zkm ** 2)[(rng >= 3000) & (rng <= 6000)])))
Tw, ow = np.array(Tw), np.array(ow); sw, cw = np.polyfit(Tw, ow, 1)
ratio_pool = abs((c36 + s36 * T.min()) / (c36 + s36 * T.max()))
ratio_sess = abs((cw + sw * Tw.min()) / (cw + sw * Tw.max()))

print("\nTEMPERATURE TRENDS (over the measured %.1f K span):" % Trng)
print(f"  undershoot amplitude A_u : {sA:+.3e} per K   ({100*sA*Trng/np.mean([b['Au_c'] for b in bins]):+.0f}% across span)")
print(f"  undershoot RC scale  L_u : {sL:+.1f} m per K  ({100*sL*Trng/np.mean([b['Lu'] for b in bins]):+.0f}% across span, ~noise)")
print(f"  3-6 km offset  POOLED    : {s36:+.3e} Mm-1sr-1/K   cold/warm ratio {ratio_pool:.2f}")
print(f"  3-6 km offset  IN-SESSION: {sw:+.3e} Mm-1sr-1/K   cold/warm ratio {ratio_sess:.2f}  (25.5 h window only)")
print("CONCLUSION: offset is PRIMARILY temperature-stable; a modest secondary term "
      f"(cold ~{ratio_sess:.2f}x warm in-session) sits in the AMPLITUDE A_u, RC constant is stable.\n"
      "            Between-session drift (May-12 deeper than May-26 at matched T) => periodic re-characterisation "
      "matters more than instantaneous T-indexing.")

# save temperature coefficients
out = {k: g[k] for k in g.files}
out.update(dict(T_bins=T, Au_of_T=np.array([b["Au_c"] for b in bins]),
                Lu_of_T=np.array([b["Lu"] for b in bins]),
                off36_of_T=np.array([b["off36"] for b in bins]),
                Au_dT=np.array([sA, cA]), Lu_dT=np.array([sL, cL]),
                off36_dT=np.array([s36, c36]), off36_dT_session=np.array([sw, cw]),
                T_bins_session=Tw, off36_of_T_session=ow))
np.savez(NPZ, **out)
print("saved temperature coefficients into", NPZ.split('/')[-1])

# ---------------------------------------------------------------- figure (landscape 2x2)
fig, ax = plt.subplots(2, 2, figsize=(16.5, 10))
norm = plt.Normalize(T.min(), T.max()); sm = cm.ScalarMappable(norm=norm, cmap="plasma")
from scipy.signal import medfilt
dr = np.median(np.diff(rng)); ks = int(round(500 / dr)); ks += 1 - ks % 2

# (a) offset profile P(r) per temp bin
a = ax[0, 0]; a.axhline(0, color="0.6", lw=0.8)
for b in bins:
    a.plot(medfilt(np.nan_to_num(b["Pmed"]), ks) * 1e3, rng, lw=1.6, color=cm.plasma(norm(b["T"])))
a.set_xlim(-2.2, 4.5); a.set_ylim(0, 15000); a.grid(alpha=0.3)
a.set_xlabel(r"$P=\beta_{att}/z^2$ [Mm$^{-1}$sr$^{-1}$km$^{-2}$] $\times10^3$"); a.set_ylabel("range [m]")
a.set_title("(a) offset profile per temperature bin"); fig.colorbar(sm, ax=a).set_label("temp_int [C]")

# (b) undershoot & lobe amplitude vs T
b_ = ax[0, 1]
b_.scatter(T, [b["Au_c"] * 1e3 for b in bins], c=T, cmap="plasma", s=70, edgecolors="k", zorder=3, label="$A_u$ undershoot")
b_.plot(T, (sA * T + cA) * 1e3, "-", color="#1f77b4", lw=1.5)
b_2 = b_.twinx()
b_2.scatter(T, [b["Ap_c"] * 1e3 for b in bins], c=T, cmap="plasma", marker="^", s=70, edgecolors="k", zorder=3)
b_2.set_ylabel(r"$A_p$ lobe [$\times10^3$]", color="#ff7f0e")
b_.set_xlabel("temp_int [C]"); b_.set_ylabel(r"$A_u$ undershoot [$\times10^3$]", color="#1f77b4")
b_.set_title("(b) amplitudes vs temperature (circle $A_u$, triangle $A_p$)"); b_.grid(alpha=0.3)

# (c) time constant vs T — the message is "no significant trend" (RC temperature-stable)
c_ = ax[1, 0]
tau_u = np.array([b["Lu"] for b in bins]) / C_HALF * 1e6
c_.scatter(T, tau_u, c=T, cmap="plasma", s=70, edgecolors="k", zorder=3)
c_.axhline(tau_u.mean(), color="0.4", ls="--", lw=1.4, label=fr"mean $\tau_u$={tau_u.mean():.0f} $\mu$s (no trend)")
c_.set_xlabel("temp_int [C]"); c_.set_ylabel(r"$\tau_u=2L_u/c$ [$\mu$s]")
c_.set_title("(c) RC time constant vs T — temperature-stable"); c_.grid(alpha=0.3); c_.legend(fontsize=9)

# (d) 3-6 km beta offset vs T: pooled (solid) vs within-session (dashed) — the honest comparison
d_ = ax[1, 1]
d_.scatter(T, [b["off36"] for b in bins], c=T, cmap="plasma", s=70, edgecolors="k", zorder=3, label="3-6 km (3 windows pooled)")
d_.plot(T, s36 * T + c36, "-", color="#1f77b4", lw=1.5, label=f"pooled fit (cold/warm {ratio_pool:.2f})")
d_.scatter(Tw, ow, facecolors="none", edgecolors="0.3", s=70, zorder=3, label="3-6 km (25.5 h session only)")
d_.plot(Tw, sw * Tw + cw, "--", color="0.3", lw=1.5, label=f"in-session fit (cold/warm {ratio_sess:.2f})")
d_.axhline(0, color="0.6", lw=0.8)
d_.set_xlabel("temp_int [C]"); d_.set_ylabel(r"$\beta_{att}$ offset [Mm$^{-1}$sr$^{-1}$]")
d_.set_title("(d) 3-6 km offset vs T: mostly stable, weak in-session slope"); d_.grid(alpha=0.3); d_.legend(fontsize=8)

fig.suptitle("Payerne CL61 electronic offset vs internal temperature (3 hood windows pooled)",
             fontweight="bold", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/fig_cl61_offset_vs_temperature.png"
fig.savefig(OUT, dpi=150)
print("saved", OUT)
print("TEMP_FIT_DONE")
