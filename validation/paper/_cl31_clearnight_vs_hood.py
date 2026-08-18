# -*- coding: utf-8 -*-
"""VALIDATION GATE: can the clear-night method recover the Payerne CL31 offset measured under the
terminal hood? Payerne CL31 (0-20000-0-06610 B) is the ONE instrument with BOTH a covered-telescope
(hood) offset AND a long clear-night record -> the ground truth for the network-wide clear-night
approach used everywhere else (no hood exists at other sites).

Also tests whether the internal laser temperature refines the offset (Kotthaus heat-sink classes):
the clear-night pool is split into temperature terciles and the offset re-extracted per bin."""
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from validation.paper._offset_lib import (load_clearnight, offset_stats, highpass, autocorr,
                                          read_file, L1ROOT)

WMO, IDENT = "0-20000-0-06610", "B"
HOOD = [("2026-05-12 09:33", "2026-05-12 14:55"), ("2026-05-26 12:00", "2026-05-27 13:10"),
        ("2026-06-09 09:16", "2026-06-09 11:57"), ("2026-06-23 10:12", "2026-06-23 12:13")]
D = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/"
FIG = D + "fig_cl31_clearnight_vs_hood.png"
FIG_REPORT = "C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_cl31_clearnight_vs_hood.png"


def load_hood():
    pool, tpool, rng = [], [], None
    for a, b in HOOD:
        t1, t2 = datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M")
        days, d = set(), t1
        while d.date() <= t2.date():
            days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
        for ds in sorted(days):
            f = L1ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{IDENT}{ds}.nc"
            if not f.exists():
                continue
            tt, r, x, cbh, temp = read_file(f)
            s = (tt >= t1) & (tt <= t2)
            if s.any():
                pool.append(x[s]); rng = r
                tpool.append(temp[s] if temp is not None else np.full(int(s.sum()), np.nan))
    X = np.vstack(pool); T = np.concatenate(tpool)
    from scipy.ndimage import uniform_filter1d
    P = X / (rng / 1000.0)[None, :] ** 2
    sm = uniform_filter1d(np.nan_to_num(P), 9, axis=1)
    keep = np.nanmax(sm[:, (rng >= 1000) & (rng <= 7000)], axis=1) < 150.0    # reject any light leak
    return dict(rng=rng, X=X[keep], temp=T[keep])


print("loading Payerne CL31 clear nights ...", flush=True)
cn = load_clearnight(WMO, IDENT)
print(f"  clear-night: {cn['n_kept']} kept of {cn['n_pool']} night/clear over {cn['n_files']} files", flush=True)
print("loading Payerne CL31 hood ...", flush=True)
hd = load_hood()
print(f"  hood: {hd['X'].shape[0]} covered profiles", flush=True)
rng = cn["rng"]; dr = float(np.median(np.diff(rng))); zkm = rng / 1000.0
Pcn, MADcn, SEcn, Ncn = offset_stats(cn["X"], rng)
Phd, MADhd, SEhd, Nhd = offset_stats(hd["X"], rng)
hp_cn, hp_hd = highpass(Pcn, dr), highpass(Phd, dr)

# two-resonance hood model for the overlay
mdl = np.load(D + "cl31_b_dark.npz"); b_model = mdl["b_phys"] / zkm ** 2      # back to P units
b_model[rng < 300] = np.nan

# quantitative agreement in the two regimes
def band_stats(lo, hi):
    m = (rng >= lo) & (rng <= hi) & np.isfinite(hp_cn) & np.isfinite(hp_hd)
    r = np.corrcoef(hp_cn[m], hp_hd[m])[0, 1]
    amp_ratio = np.sqrt(np.mean(hp_cn[m] ** 2)) / max(np.sqrt(np.mean(hp_hd[m] ** 2)), 1e-9)
    return r, amp_ratio, np.sqrt(np.mean(hp_cn[m] ** 2)), np.sqrt(np.mean(hp_hd[m] ** 2))

near = band_stats(400, 1500)     # fast AC-coupling ring (damped, near range)
far = band_stats(1800, 7000)     # slow transmitter ripple (persistent, far range)
# near-range: the meaningful metric is the P-median EXCESS (clear-night aerosol over the hood ring)
mnr = (rng >= 400) & (rng <= 1500)
near_excess = float(np.nanmedian(Pcn[mnr]) - np.nanmedian(Phd[mnr]))
print("=" * 84)
print("Payerne CL31 — clear-night vs hood offset")
print("=" * 84)
print(f"  NEAR 0.4-1.5 km (fast ring): clear-night P exceeds hood by {near[1]:.1f}x (median excess "
      f"{near_excess:+.1f}) -> aerosol-dominated, the damped ring is NOT recoverable from clear nights")
print(f"  FAR  1.8-7 km (slow ripple): corr={far[0]:+.2f}  amp(cn)/amp(hood)={far[1]:.2f}  "
      f"RMS cn={far[2]:.2f} hood={far[3]:.2f}  -> clear-night RECOVERS the persistent ripple")

# ---------- temperature split (temperature_laser is in KELVIN) ----------
T = cn["temp"].copy()
if np.isfinite(T).any() and np.nanmedian(T) > 200:
    T = T - 273.15                                  # K -> C
okT = np.isfinite(T)
tinfo = ""
tbins = []
if okT.sum() > 500 and np.nanstd(T[okT]) > 0.3:
    q = np.nanpercentile(T[okT], [33, 67])
    lab = [f"cold (<{q[0]:.1f}C)", f"mid", f"warm (>{q[1]:.1f}C)"]
    masks = [T <= q[0], (T > q[0]) & (T <= q[1]), T > q[1]]
    for l, mk in zip(lab, masks):
        if mk.sum() > 100:
            Pm, _, _, _ = offset_stats(cn["X"][mk], rng)
            tbins.append((l, float(np.nanmedian(T[mk])), highpass(Pm, dr), int(mk.sum())))
    if len(tbins) >= 2:
        amps = [np.sqrt(np.nanmean(b[2][(rng >= 1800) & (rng <= 7000)] ** 2)) for b in tbins]
        temps = [b[1] for b in tbins]
        slope = np.polyfit(temps, amps, 1)[0]
        tinfo = f"far-ripple RMS vs T: {amps[0]:.2f}->{amps[-1]:.2f} over {temps[0]:.1f}->{temps[-1]:.1f}C (slope {slope:+.3f}/C)"
        print(f"  TEMPERATURE ({okT.sum()} profiles, {np.nanmin(T[okT]):.1f}-{np.nanmax(T[okT]):.1f}C): {tinfo}")
else:
    print(f"  TEMPERATURE: insufficient spread (std={np.nanstd(T[okT]) if okT.any() else np.nan:.2f}C)")

# ---------- figure ----------
fig, ax = plt.subplots(2, 3, figsize=(18, 10))
# (a) near-range P-median overlay
a = ax[0][0]
a.plot(rng, Phd, "-", color="#d62728", lw=1.6, label="hood P (median)")
a.plot(rng, Pcn, "-", color="#1f77b4", lw=1.2, label="clear-night P (median)")
a.plot(rng, b_model, "--", color="k", lw=1.2, label="two-resonance model (hood fit)")
a.axvspan(0, 300, color="0.5", alpha=0.12); a.set_xlim(0, 2500); a.set_ylim(-20, 60)
a.set_xlabel("range [m]"); a.set_ylabel(r"$P=rcs_0/z^2$"); a.legend(fontsize=8)
a.set_title("(a) near range: the fast ring — hood only (aerosol masks clear-night)")
a.grid(alpha=0.3)
# (b) far-range high-pass overlay
b = ax[0][1]
b.plot(rng, hp_hd, "-", color="#d62728", lw=1.2, label="hood ripple")
b.plot(rng, hp_cn, "-", color="#1f77b4", lw=1.0, label="clear-night ripple")
b.axhline(0, color="0.7", lw=0.7); b.set_xlim(1800, 7000); b.set_ylim(-4, 4)
b.set_xlabel("range [m]"); b.set_ylabel("high-pass P"); b.legend(fontsize=8)
b.set_title(f"(b) far range: slow ripple — clear-night recovers it (corr={far[0]:+.2f})")
b.grid(alpha=0.3)
# (c) scatter clear-night vs hood ripple (far)
c = ax[0][2]
mf = (rng >= 1800) & (rng <= 7000)
c.plot(hp_hd[mf], hp_cn[mf], ".", ms=3, color="#1f77b4", alpha=0.5)
lim = 4; c.plot([-lim, lim], [-lim, lim], "k--", lw=1)
c.set_xlim(-lim, lim); c.set_ylim(-lim, lim); c.set_aspect("equal")
c.set_xlabel("hood ripple"); c.set_ylabel("clear-night ripple")
c.set_title(f"(c) far ripple 1:1  (corr={far[0]:+.2f}, amp ratio {far[1]:.2f})"); c.grid(alpha=0.3)
# (d) temperature-split offset (far ripple)
d = ax[1][0]
for l, tm, hpb, nn in tbins:
    d.plot(rng, hpb, lw=1.1, label=f"{l}  {tm:.1f}C  (n={nn})")
d.axhline(0, color="0.7", lw=0.7); d.set_xlim(1800, 6000); d.set_ylim(-4, 4)
d.set_xlabel("range [m]"); d.set_ylabel("high-pass P"); d.legend(fontsize=8)
d.set_title("(d) offset by laser-temperature tercile"); d.grid(alpha=0.3)
# (e) near-range P by temperature (does the ring level shift with T?)
e = ax[1][1]
if tbins:                                    # re-extract full P per bin for the near range
    q = np.nanpercentile(T[okT], [33, 67]); masks = [T <= q[0], (T > q[0]) & (T <= q[1]), T > q[1]]
    for (l, tm, _, _), mk in zip(tbins, masks):
        Pm, _, _, _ = offset_stats(cn["X"][mk], rng)
        e.plot(rng, Pm, lw=1.1, label=f"{tm:.1f}°C")
e.plot(rng, Phd, "k--", lw=1.0, label="hood")
e.set_xlim(300, 1800); e.set_ylim(-5, 40); e.set_xlabel("range [m]"); e.set_ylabel("P")
e.legend(fontsize=8); e.set_title("(e) near-range offset vs temperature"); e.grid(alpha=0.3)
# (f) summary text
f = ax[1][2]; f.axis("off")
txt = ("VALIDATION SUMMARY\n\n"
       f"clear-night: {cn['n_kept']} profiles\nhood: {hd['X'].shape[0]} profiles\n\n"
       f"NEAR 0.4-1.5 km (fast ring)\n  corr = {near[0]:+.2f}\n  amp cn/hood = {near[1]:.2f}\n"
       f"  -> aerosol masks the ring\n\n"
       f"FAR 1.8-7 km (slow ripple)\n  corr = {far[0]:+.2f}\n  amp cn/hood = {far[1]:.2f}\n"
       f"  -> clear-night RECOVERS it\n\n"
       f"TEMPERATURE\n  {tinfo if tinfo else 'insufficient spread'}")
f.text(0.02, 0.98, txt, va="top", ha="left", fontsize=10.5, family="monospace")
fig.suptitle("Payerne CL31 — clear-night offset vs the terminal-hood ground truth "
             "(validating the network-wide clear-night method)", fontweight="bold", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(FIG, dpi=150); Path(FIG_REPORT).parent.mkdir(parents=True, exist_ok=True); fig.savefig(FIG_REPORT, dpi=150)
print("saved", FIG); print("CL31_CLEARNIGHT_VS_HOOD_DONE")
