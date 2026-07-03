# -*- coding: utf-8 -*-
"""Post-process the network offset scan (_network_offset_scan.py npz cache) into the FINAL operational
coefficient table. The first-pass scan fit a FIXED 4-gate (f_sample/4) ripple; the network shows the
ripple period is UNIT-specific (40/60/80 m ...), and an exact period matters (a global sinusoid at
40 m dephases over 10 km if the true period is 40.4 m). Here we refine the period per instrument by a
fine least-squares scan over a short, phase-stable window, then report the ripple amplitude, coherence,
undamped extent and a correctability verdict. Temperature slope is merged from the scan CSV.

Deliverable: network_offset_coeffs_final.csv/json -- per instrument the coefficients to de-ripple the
aerosol backscatter (subtract b_phys, cached per instrument in network_offset/<key>.npz)."""
import sys, glob, csv, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from validation.paper._offset_lib import autocorr, C_LIGHT

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUTDIR = D / "network_offset"
FIT_LO, FIT_HI = 2000.0, 6000.0          # short, phase-stable window for the ripple period/amplitude
AC_CORRECTABLE = 0.50                      # coherence threshold to call the ripple operationally correctable
REL_CORRECTABLE = 1.0                      # and rel-amplitude threshold (% of mid-range signal)


def refine_ripple(hp, rng, dr, lam0):
    """Fine period scan around the coarse autocorr period, then linear LSQ amplitude/phase."""
    m = (rng >= FIT_LO) & (rng <= FIT_HI); r = rng[m]; y = np.nan_to_num(hp[m])
    if r.size < 50 or not np.isfinite(lam0):
        return None
    lo, hi = max(2 * dr, lam0 * 0.6), lam0 * 1.6
    grid = np.arange(lo, hi, 0.2)
    def cost(L):
        A = np.c_[np.ones_like(r), np.cos(2 * np.pi * r / L), np.sin(2 * np.pi * r / L)]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        return np.sum((A @ c - y) ** 2)
    Lam = float(grid[int(np.argmin([cost(L) for L in grid]))])
    A = np.c_[np.ones_like(r), np.cos(2 * np.pi * r / Lam), np.sin(2 * np.pi * r / Lam),
              np.cos(2 * np.pi * r / (Lam / 2)), np.sin(2 * np.pi * r / (Lam / 2))]
    c, *_ = np.linalg.lstsq(A, y, rcond=None); fit = A @ c
    R2 = 1 - np.sum((fit - y) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-30)
    return dict(Lam=Lam, amp=float(np.hypot(c[1], c[2])), amp_harm=float(np.hypot(c[3], c[4])),
                R2=float(R2), f_MHz=C_LIGHT / (2 * Lam) / 1e6, gates=Lam / dr)


# merge temperature slope from the first-pass scan CSV
tslope = {}
scan_csv = D / "network_offset_coeffs.csv"
if scan_csv.exists():
    for row in csv.DictReader(open(scan_csv, encoding="utf-8")):
        tslope[(row["wmo"], row["ident"], row["itype"])] = (row["temp_slope_per_C"], row["temp_lo_C"], row["temp_hi_C"])

rows = []
for f in sorted(glob.glob(str(OUTDIR / "*.npz"))):
    z = np.load(f, allow_pickle=True)
    name = Path(f).stem                      # <wmo>_<ident>_<itype>
    itype = name.split("_")[-1]; ident = name.split("_")[-2]; wmo = "_".join(name.split("_")[:-2])
    site = str(z["site"])
    rng = z["rng"]; hp = z["hp"]; Pmed = z["Pmed"]; dr = float(np.median(np.diff(rng)))
    Lam_dom = float(z["Lam_dom"]) if np.isfinite(z["Lam_dom"]) else np.nan
    sig_ref = float(np.nanmedian(np.abs(Pmed[(rng >= 800) & (rng <= 1500)])))
    sig_ref = sig_ref if sig_ref > 0 else np.nan
    # coherence at the coarse period
    lag, ac = autocorr(hp[(rng >= FIT_LO) & (rng <= FIT_HI)], dr)
    coh = float(np.interp(Lam_dom, lag, ac)) if np.isfinite(Lam_dom) else np.nan
    ref = refine_ripple(hp, rng, dr, Lam_dom) if np.isfinite(Lam_dom) else None
    # undamped extent: RMS(6-10 km) / RMS(2-4 km)
    def rms(a, b):
        s = hp[(rng >= a) & (rng <= b)]
        return float(np.sqrt(np.nanmean(s ** 2))) if np.isfinite(s).any() else np.nan
    undamp = rms(6000, 10000) / rms(2000, 4000) if rms(2000, 4000) else np.nan
    rel_amp = 100 * ref["amp"] / sig_ref if (ref and np.isfinite(sig_ref)) else 0.0
    relrms = 100 * rms(2000, 8000) / sig_ref if np.isfinite(sig_ref) else np.nan
    correctable = bool(np.isfinite(coh) and coh >= AC_CORRECTABLE and rel_amp >= REL_CORRECTABLE)
    ts = tslope.get((wmo, ident, itype), ("nan", "nan", "nan"))
    rows.append(dict(itype=itype, site=site, wmo=wmo, ident=ident, dr=dr, sig_ref=sig_ref,
                     Lam=ref["Lam"] if ref else np.nan, gates=ref["gates"] if ref else np.nan,
                     f_MHz=ref["f_MHz"] if ref else np.nan, amp=ref["amp"] if ref else np.nan,
                     rel_amp_pct=rel_amp, R2=ref["R2"] if ref else np.nan, coherence=coh,
                     relrms_pct=relrms, undamp=undamp, tslope=float(ts[0]) if ts[0] not in ("nan", "") else np.nan,
                     correctable=correctable, _rng=rng, _hp=hp))

rows.sort(key=lambda r: ({"CL51": 0, "CL31": 1, "CL61": 2}[r["itype"]], -np.nan_to_num(r["rel_amp_pct"])))
COLS = ["itype", "site", "wmo", "ident", "range_gate_m", "ripple_period_m", "ripple_gates",
        "ripple_f_MHz", "ripple_rel_amp_pct", "ripple_R2", "coherence", "relrms_pct", "undamped_ratio",
        "temp_slope_per_C", "correctable"]
with open(D / "network_offset_coeffs_final.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(COLS)
    for r in rows:
        w.writerow([r["itype"], r["site"], r["wmo"], r["ident"], f"{r['dr']:.0f}", f"{r['Lam']:.1f}",
                    f"{r['gates']:.2f}", f"{r['f_MHz']:.3f}", f"{r['rel_amp_pct']:.2f}", f"{r['R2']:.2f}",
                    f"{r['coherence']:.2f}", f"{r['relrms_pct']:.2f}", f"{r['undamp']:.2f}",
                    f"{r['tslope']:.4g}", int(r["correctable"])])
Path(D / "network_offset_coeffs_final.json").write_text(json.dumps(
    [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("_")}
     for r in rows], indent=1, default=str))

# ---------- summary ----------
nc = {t: sum(1 for r in rows if r["itype"] == t) for t in ("CL51", "CL31", "CL61")}
ncorr = {t: sum(1 for r in rows if r["itype"] == t and r["correctable"]) for t in ("CL51", "CL31", "CL61")}
print("=" * 92)
print("NETWORK OFFSET COEFFICIENTS (clear-night, refined per-unit period)")
print("=" * 92)
for r in rows:
    flag = "CORRECTABLE" if r["correctable"] else "-"
    print(f"  {r['itype']:5} {r['site'][:22]:22} Λ={r['Lam']:6.1f}m ({r['gates']:4.1f} gates, {r['f_MHz']:5.2f}MHz)"
          f"  amp={r['rel_amp_pct']:4.1f}%  coh={r['coherence']:+.2f}  undamp={r['undamp']:.2f}"
          f"  Tslp={r['tslope']:+.3f}  {flag}")
print(f"\n  correctable ripple: CL51 {ncorr['CL51']}/{nc['CL51']}, CL31 {ncorr['CL31']}/{nc['CL31']}, "
      f"CL61 {ncorr['CL61']}/{nc['CL61']}")

# ---------- figure ----------
COL = {"CL61": "#d62728", "CL51": "#1f77b4", "CL31": "#2ca02c"}
fig, ax = plt.subplots(2, 2, figsize=(16, 10))
a = ax[0][0]; x = 0; ticks = []; labs = []
for it in ("CL51", "CL31", "CL61"):
    for r in [r for r in rows if r["itype"] == it]:
        col = COL[it] if r["correctable"] else "0.75"
        a.bar(x, r["rel_amp_pct"], color=col, edgecolor=COL[it], lw=1.2); ticks.append(x); labs.append(r["site"][:11]); x += 1
    x += 1
a.axhline(REL_CORRECTABLE, color="k", ls=":", lw=1, label=f"correctable threshold {REL_CORRECTABLE}%")
a.set_xticks(ticks); a.set_xticklabels(labs, rotation=90, fontsize=6.5)
a.set_ylabel("ripple amplitude [% of mid-range signal]")
a.set_title("(a) recoverable ripple per instrument (solid = correctable, grey = below threshold)")
a.grid(alpha=0.3, axis="y"); a.legend(fontsize=8)
b = ax[0][1]
for it in ("CL51", "CL31", "CL61"):
    grp = [r for r in rows if r["itype"] == it and np.isfinite(r["Lam"])]
    b.scatter([r["gates"] for r in grp], [r["coherence"] for r in grp], c=COL[it], s=45, label=it)
for g in (4, 6, 8):
    b.axvline(g, color="0.8", ls=":", lw=0.8)
b.set_xlabel("ripple period [range gates]"); b.set_ylabel("coherence (autocorr at Λ)")
b.set_title("(b) ripple period clusters at integer gates (f$_s$/N digitizer interleave)")
b.legend(fontsize=9); b.grid(alpha=0.3)
c = ax[1][0]
for it in ("CL51", "CL31", "CL61"):
    grp = sorted([r for r in rows if r["itype"] == it], key=lambda r: -np.nan_to_num(r["rel_amp_pct"]))
    if grp and np.isfinite(grp[0]["sig_ref"]):
        r = grp[0]; m = (r["_rng"] >= 2500) & (r["_rng"] <= 3100)
        c.plot(r["_rng"][m], 100 * r["_hp"][m] / r["sig_ref"], color=COL[it], lw=1.1,
               label=f"{it}: {r['site'][:12]} (Λ={r['Lam']:.0f}m)")
c.axhline(0, color="0.7", lw=0.7); c.set_xlim(2500, 3100)
c.set_xlabel("range [m]"); c.set_ylabel("relative ripple [%]")
c.set_title("(c) strongest recoverable ripple per type"); c.legend(fontsize=8); c.grid(alpha=0.3)
d = ax[1][1]; x = 0; ticks = []; labs = []
for it in ("CL51", "CL31", "CL61"):
    for r in [r for r in rows if r["itype"] == it and np.isfinite(r["tslope"])]:
        d.bar(x, r["tslope"], color=COL[it]); ticks.append(x); labs.append(r["site"][:11]); x += 1
    x += 1
d.axhline(0, color="0.5", lw=0.8); d.set_xticks(ticks); d.set_xticklabels(labs, rotation=90, fontsize=6.5)
d.set_ylabel("ripple RMS slope [per °C]"); d.set_title("(d) internal-temperature dependence (mostly weak)")
d.grid(alpha=0.3, axis="y")
fig.suptitle("Network electronic-offset coefficients (clear-night, per-unit refined period): the ripple "
             "is a UNIT-specific f$_s$/N digitizer pattern — correctable for some CL51/CL31, absent in CL61",
             fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
FIG = D / "fig_network_offset_coeffs.png"
FIGR = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_network_offset_coeffs.png")
fig.savefig(FIG, dpi=150); fig.savefig(FIGR, dpi=150)
print("saved", FIG); print("NETWORK_COEFFS_DONE")
