# -*- coding: utf-8 -*-
"""Apply the CL31 offset model to the Uccle CL51 (0-20000-0-06447 ident A) — and find it does NOT
transfer. Uccle has no terminal-hood campaign (only Payerne did), so the offset is extracted from
clear-night per-gate medians of P = rcs_0/z^2 (nighttime, cloud-free, cleanest-50% aerosol), where
the FIXED electronic pattern survives the median while the variable atmosphere averages to a smooth
baseline. A gaussian high-pass (sigma~350 m) removes that smooth molecular/aerosol baseline and
isolates the electronic ripple.

Result: the CL51 offset is dominated by a **fixed 40 m (4-range-gate) ripple**, undamped to >10 km
(autocorrelation ~0.99 at 40/80/120 m; RMS flat at ~1.55 from 2-10 km). 40 m = 4x the 10 m gate ->
f = c/2Lambda = 3.75 MHz = f_sample/4 (the 15 MHz range-gate clock): a 4-way ADC-interleave /
digitizer fixed-pattern ripple. This is the "sensor-specific frequency" transmitter/receiver ripple
Kotthaus et al. (2016) describe but decline to model. It is ABSENT as a coherent feature in the
Payerne CL31 (its 40 m content is incoherent noise, RMS 0.56, autocorr <0.2): the two-resonance CL31
model (amplifier ring 1053 m + ripple 5080 m) does not describe the CL51. Each Vaisala unit imprints
its OWN fixed additive pattern; the universal principle is the fixed offset (Kotthaus P^bgi),
not a universal frequency.

Caveat: without a hood, the near range (<1.8 km) is dominated by Uccle's real maritime boundary-layer
aerosol, so only the oscillatory ripple (>1.8 km) is cleanly recoverable, not the smooth near-range
electronic offset. b_phys therefore carries the ripple only (zeroed below 1.8 km)."""
import sys, warnings, glob
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
from scipy.ndimage import gaussian_filter1d, uniform_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 2.99792458e8
UCCLE_A = Path("D:/E-PROFILE_L1_2026/0-20000-0-06447/2026")   # CL51
PAY_B = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")           # CL31 hood (for the contrast)
MONTHS = ("03", "04", "05", "06")
HP_SIGMA_M = 350.0
CLEAN_LO, CLEAN_HI = 1800.0, 12000.0    # window where atmosphere is smooth -> ripple is clean
D = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/"
OUT_NPZ = D + "cl51_b_dark.npz"
OUT_FIG = D + "fig_cl51_uccle_offset_model.png"
FIG_REPORT = "C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_cl51_uccle_offset_model.png"

HOOD_CL31 = [("2026-05-12 09:33", "2026-05-12 14:55"), ("2026-05-26 12:00", "2026-05-27 13:10"),
             ("2026-06-09 09:16", "2026-06-09 11:57"), ("2026-06-23 10:12", "2026-06-23 12:13")]


def _read(f):
    with Dataset(f) as nc:
        tv = np.asarray(nc.variables["time"][:], "f8")
        tt = np.array([datetime(1970, 1, 1) + timedelta(days=x) for x in tv])
        r = np.asarray(nc.variables["range"][:], "f8")
        x = np.asarray(nc.variables["rcs_0"][:], "f8")
        if x.shape[0] != tt.size:
            x = x.T
        cbh = None
        if "cloud_base_height" in nc.variables:
            cbh = np.asarray(nc.variables["cloud_base_height"][:], "f8")
            if cbh.ndim == 1:
                cbh = cbh[:, None]
            if cbh.shape[0] != tt.size:
                cbh = cbh.T
    return tt, r, x, cbh


def load_uccle_clearnight():
    """Pooled nighttime (22-03 UTC), cloud-free (all CBH layers == fill), cleanest-50% aerosol."""
    pool, rng = [], None
    for mm in MONTHS:
        for f in sorted(glob.glob(str(UCCLE_A / mm / "L1_0-20000-0-06447_A*.nc"))):
            tt, r, x, cbh = _read(f)
            hr = np.array([d.hour + d.minute / 60 for d in tt])
            night = (hr >= 22) | (hr <= 3)
            nocloud = np.nanmax(cbh, axis=1) < 0 if cbh is not None else np.ones(tt.size, bool)
            sel = night & nocloud
            if sel.any():
                pool.append(x[sel]); rng = r
    X = np.vstack(pool); zkm = rng / 1000.0
    P = X / zkm[None, :] ** 2
    sm = uniform_filter1d(np.nan_to_num(P), 9, axis=1)
    lowaer = np.nanmedian(sm[:, (rng >= 1500) & (rng <= 4000)], axis=1)
    keep = lowaer < np.nanpercentile(lowaer, 50)
    return X[keep], rng, P[keep]


def load_cl31_hood():
    pool, rng = [], None
    for a, b in HOOD_CL31:
        t1, t2 = datetime.strptime(a, "%Y-%m-%d %H:%M"), datetime.strptime(b, "%Y-%m-%d %H:%M")
        days, d = set(), t1
        while d.date() <= t2.date():
            days.add(d.strftime("%Y%m%d")); d += timedelta(days=1)
        for ds in sorted(days):
            f = PAY_B / ds[:4] / ds[4:6] / f"L1_0-20000-0-06610_B{ds}.nc"
            if not f.exists():
                continue
            tt, r, x, _ = _read(f)
            s = (tt >= t1) & (tt <= t2)
            if s.any():
                pool.append(x[s]); rng = r
    X = np.vstack(pool); zkm = rng / 1000.0
    P = X / zkm[None, :] ** 2
    sm = uniform_filter1d(np.nan_to_num(P), 9, axis=1)
    P = P[np.nanmax(sm[:, (rng >= 1000) & (rng <= 7000)], axis=1) < 150.0]
    return rng, P


def highpass(Pmed, dr):
    return Pmed - gaussian_filter1d(np.nan_to_num(Pmed), HP_SIGMA_M / dr)


def autocorr(y, dr):
    y = np.nan_to_num(y - np.nanmean(y))
    ac = np.correlate(y, y, "full")[len(y) - 1:]
    return np.arange(len(ac)) * dr, ac / ac[0]


def rms_by_km(hp, rng):
    # from 2 km up: below that Uccle's real boundary-layer aerosol contaminates the clear-night high-pass
    los = np.arange(2000, 10000, 1000)
    out = []
    for lo in los:
        seg = hp[(rng >= lo) & (rng < lo + 1000)]
        out.append(np.sqrt(np.nanmean(seg ** 2)) if np.isfinite(seg).any() else np.nan)
    return los / 1000.0, np.array(out)


# ---------- load ----------
print("loading Uccle CL51 clear nights ...", flush=True)
Xc, rng, Pc = load_uccle_clearnight()
dr = float(np.median(np.diff(rng))); zkm = rng / 1000.0
Pmed = np.nanmedian(Pc, axis=0)
Pmad = 1.4826 * np.nanmedian(np.abs(Pc - Pmed[None, :]), axis=0)
hp = highpass(Pmed, dr)
print(f"  CL51: {Xc.shape[0]} clear-night profiles, dr={dr:.0f} m", flush=True)
print("loading Payerne CL31 hood (contrast) ...", flush=True)
rng31, P31 = load_cl31_hood(); dr31 = float(np.median(np.diff(rng31)))
Pmed31 = np.nanmedian(P31, axis=0); hp31 = highpass(Pmed31, dr31)

# ---------- characterise the CL51 ripple ----------
lag, ac = autocorr(hp[(rng >= CLEAN_LO) & (rng <= CLEAN_HI)], dr)
pk, _ = find_peaks(ac, height=0.3)
Lam = lag[pk[0]] if len(pk) else np.nan          # fundamental period = first strong autocorr peak
f_ripple = C / (2 * Lam) if np.isfinite(Lam) else np.nan
f_sample = C / (2 * dr)                            # range-gate clock
# fit the fixed period-4 pattern: DC + 40 m (fundamental) + 20 m (Nyquist harmonic), undamped
m = (rng >= CLEAN_LO) & (rng <= CLEAN_HI); r = rng[m]; y = np.nan_to_num(hp[m])
A = np.c_[np.ones_like(r), np.cos(2 * np.pi * r / Lam), np.sin(2 * np.pi * r / Lam),
          np.cos(2 * np.pi * r / (Lam / 2)), np.sin(2 * np.pi * r / (Lam / 2))]
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
fit = A @ coef
amp_fund = np.hypot(coef[1], coef[2]); amp_harm = np.hypot(coef[3], coef[4])
R2 = 1 - np.sum((fit - y) ** 2) / np.sum((y - y.mean()) ** 2)
kmc51, rms51 = rms_by_km(hp, rng); kmc31, rms31 = rms_by_km(hp31, rng31)
lag31, ac31 = autocorr(hp31[(rng31 >= 2000) & (rng31 <= 7000)], dr31)

print("=" * 84)
print("Uccle CL51 offset — the CL31 two-resonance model does NOT transfer")
print("=" * 84)
print(f"  dominant ripple: Lambda={Lam:.0f} m ({Lam/dr:.0f} range gates) -> f=c/2L={f_ripple/1e6:.2f} MHz")
print(f"                   = f_sample/{f_sample/f_ripple:.0f}  (range-gate clock f_sample={f_sample/1e6:.1f} MHz)")
print(f"  amplitude: fundamental {amp_fund:.2f}, Nyquist-harmonic {amp_harm:.2f} (P=rcs_0/z^2 units)")
print(f"  undamped: RMS {np.nanmin(rms51):.2f}-{np.nanmax(rms51):.2f} across 2-10 km (flat); fit R2={R2:.3f}")
print(f"  CL31 contrast: NO coherent {Lam:.0f} m ripple (autocorr peak {np.max(ac31[3:]):.2f}, RMS {np.nanmedian(rms31):.2f} noise-like)")

# ---------- correction (ripple only; near range needs a hood) ----------
b_phys = np.zeros_like(rng)
mm = rng >= CLEAN_LO
Afull = np.c_[np.ones(mm.sum()), np.cos(2 * np.pi * rng[mm] / Lam), np.sin(2 * np.pi * rng[mm] / Lam),
              np.cos(2 * np.pi * rng[mm] / (Lam / 2)), np.sin(2 * np.pi * rng[mm] / (Lam / 2))]
ripple = (Afull @ coef) - coef[0]                 # ripple about the local mean (drop the DC baseline)
b_phys[mm] = ripple * zkm[mm] ** 2                 # -> rcs_0 units
np.savez(OUT_NPZ, rng=rng, P_med=Pmed, P_mad=Pmad, b_phys=b_phys, Lambda=Lam,
         f_ripple=f_ripple, amp_fund=amp_fund, amp_harm=amp_harm, R2=R2, n_profiles=Xc.shape[0],
         coef=coef, clean_lo=CLEAN_LO)

# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
# (a) native-resolution ripple + fit
a = ax[0]; w = (rng >= 2500) & (rng <= 3100)
a.plot(rng[w], hp[w], ".-", ms=3, lw=0.8, color="#1f77b4", label="CL51 high-pass residual (10 m)")
a.plot(rng[w], (np.c_[np.ones(w.sum()), np.cos(2*np.pi*rng[w]/Lam), np.sin(2*np.pi*rng[w]/Lam),
        np.cos(2*np.pi*rng[w]/(Lam/2)), np.sin(2*np.pi*rng[w]/(Lam/2))] @ coef),
       "-", color="#d62728", lw=1.8, label=fr"period-4 fit ($\Lambda$={Lam:.0f} m)")
a.plot(rng31[(rng31 >= 2500) & (rng31 <= 3100)], hp31[(rng31 >= 2500) & (rng31 <= 3100)],
       ".-", ms=2, lw=0.6, color="0.6", label="CL31 (no coherent ripple)")
a.axhline(0, color="0.7", lw=0.7); a.set_xlim(2500, 3100)
a.set_xlabel("range [m]"); a.set_ylabel(r"high-pass $P$  [V m$^2$ km$^{-2}$]")
a.set_title(fr"(a) CL51: fixed 40 m ripple = $f_{{sample}}/4$ ADC interleave ({f_ripple/1e6:.2f} MHz)")
a.legend(fontsize=8, loc="upper right"); a.grid(alpha=0.3)
# (b) undamped RMS vs range
b = ax[1]
b.plot(kmc51, rms51, "o-", color="#1f77b4", lw=1.8, label="CL51 (Uccle, clear nights)")
b.plot(kmc31, rms31, "s-", color="0.5", lw=1.5, label="CL31 (Payerne, hood)")
b.set_xlabel("range [km]"); b.set_ylabel("ripple RMS  [V m$^2$ km$^{-2}$]")
b.set_title("(b) CL51 ripple is UNDAMPED to 10 km; CL31 has none")
b.legend(fontsize=9); b.grid(alpha=0.3)
b.set_ylim(0, np.nanmax(np.r_[rms51, rms31]) * 1.3)
# (c) autocorrelation
c = ax[2]
c.plot(lag[:40], ac[:40], "o-", color="#1f77b4", lw=1.6, ms=4, label="CL51")
c.plot(lag31[:40], ac31[:40], "s-", color="0.5", lw=1.3, ms=3, label="CL31")
for k in (1, 2, 3, 4):
    c.axvline(k * Lam, color="#d62728", ls=":", lw=0.9)
c.axhline(0, color="0.7", lw=0.7); c.set_xlim(0, 12 * dr)
c.set_xlabel("lag [m]"); c.set_ylabel("autocorrelation")
c.set_title(fr"(c) sharp {Lam:.0f} m periodicity (CL51) vs noise (CL31)")
c.legend(fontsize=9); c.grid(alpha=0.3)
fig.suptitle("Uccle CL51 electronic offset — a fixed 4-gate (40 m, $f_s$/4) digitizer ripple, "
             "distinct from the Payerne CL31 (the two-resonance model does not transfer)",
             fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT_FIG, dpi=150)
Path(FIG_REPORT).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG_REPORT, dpi=150)
print("saved", OUT_FIG)
print("saved", FIG_REPORT)
print("CL51_UCCLE_MODEL_DONE")
