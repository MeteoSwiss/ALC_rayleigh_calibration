# -*- coding: utf-8 -*-
"""Network-scale electronic-offset scan: extract the clear-night offset for all usable CL61 (10),
10 CL51 and 10 CL31 (+ the hood-anchored Payerne CL31), fit a compact ripple model per instrument,
test the internal-temperature dependence, and emit an operational coefficient table.

Per instrument we report (validated on Payerne CL31 in _cl31_clearnight_vs_hood.py — clear-night
recovers the PERSISTENT ripple, not the near-range damped ring):
  * digitizer ripple at f_sample/4 (period = 4 range gates): the Vaisala ADC-interleave fixed pattern
  * the dominant coherent autocorrelation period (if any) and the ripple RMS in 2-8 km
  * a temperature slope (ripple RMS per deg C of internal laser temperature)
  * b_phys(range): the empirical correction (ripple x z^2, in rcs_0 units), zeroed below the
    aerosol-contaminated near range -> subtract from L1 rcs_0 to de-ripple the aerosol backscatter.

Outputs: network_offset_coeffs.csv (+ .json), per-instrument npz, and a summary figure."""
import sys, json, csv, warnings
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
from validation.paper._offset_lib import (load_clearnight, offset_stats, highpass, autocorr,
                                          dominant_period, fit_ripple, C_LIGHT)

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUTDIR = D / "network_offset"; OUTDIR.mkdir(exist_ok=True)
STRIDE = 2                       # every 2nd file: ~60 clear-night days is plenty for a fixed pattern
CLEAN_LO = 1800.0                # below this, boundary-layer aerosol contaminates the clear-night median
RMS_LO, RMS_HI = 2000.0, 8000.0

SAMPLE = json.loads(Path(r"C:/Users/hervo/AppData/Local/Temp/claude/"
                         r"C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/"
                         r"7771e6db-207a-4e4d-a515-e8c3bb4cd7c9/scratchpad/network_sample.json").read_text())
# add the hood-anchored Payerne CL31 for cross-check
SAMPLE["CL31"].append(["0-20000-0-06610", "B", 0, "PAYERNE(hood)"])

TEMP_VAR = {"CL61": "temperature_laser", "CL51": "temperature_laser", "CL31": "temperature_laser"}


def analyse(itype, wmo, ident, site):
    cn = load_clearnight(wmo, ident, temp_var=TEMP_VAR[itype], file_stride=STRIDE)
    if cn is None or cn["X"].shape[0] < 200:
        return None
    rng = cn["rng"]; dr = float(np.median(np.diff(rng)))
    Pmed, MAD, SE, N = offset_stats(cn["X"], rng)
    hp = highpass(Pmed, dr)
    zkm = rng / 1000.0
    # dimensionless reference signal (mid-range, present on clear nights) so CL61 (physical beta ~1e-7)
    # and CL31/CL51 (volts ~1-50) are comparable; rel_ripple = ripple as a fraction of it
    sig_ref = float(np.nanmedian(np.abs(Pmed[(rng >= 800) & (rng <= 1500)])))
    # (1) digitizer ripple at f_sample/4  (period = 4 gates)
    Lam_dig = 4 * dr
    dig = fit_ripple(hp, rng, Lam_dig, CLEAN_LO, min(rng.max(), 12000.0))
    # (2) dominant coherent period (autocorr) + ripple RMS
    Lam_dom, ac_dom, rms = dominant_period(hp, rng, dr, RMS_LO, RMS_HI, pmin=2 * dr, pmax=1500.0)
    # (3) temperature slope of the ripple RMS
    T = cn["temp"].copy()
    if np.isfinite(T).any() and np.nanmedian(T) > 200:
        T = T - 273.15
    okT = np.isfinite(T); tslope = np.nan; trange = (np.nan, np.nan); tamps = None
    if okT.sum() > 500 and np.nanstd(T[okT]) > 0.5:
        q = np.nanpercentile(T[okT], [25, 50, 75])
        masks = [T <= q[0], (T > q[0]) & (T <= q[2]), T > q[2]]
        pts = []
        for mk in masks:
            if mk.sum() > 100:
                Pm, *_ = offset_stats(cn["X"][mk], rng)
                rr = np.sqrt(np.nanmean(highpass(Pm, dr)[(rng >= RMS_LO) & (rng <= RMS_HI)] ** 2))
                pts.append((float(np.nanmedian(T[mk])), float(rr)))
        if len(pts) >= 2:
            tt, aa = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
            tslope = float(np.polyfit(tt, aa, 1)[0]); trange = (float(tt.min()), float(tt.max()))
            tamps = (tt, aa)
    # (4) correction b_phys = ripple x z^2, zeroed in the contaminated near range
    b_phys = hp.copy(); b_phys[rng < CLEAN_LO] = 0.0; b_phys = b_phys * zkm ** 2
    np.savez(OUTDIR / f"{wmo}_{ident}_{itype}.npz", rng=rng, Pmed=Pmed, MAD=MAD, hp=hp,
             b_phys=b_phys, dig_coef=dig["coef"], Lam_dig=Lam_dig, Lam_dom=Lam_dom, site=site)
    sr = sig_ref if (sig_ref and np.isfinite(sig_ref) and sig_ref > 0) else np.nan
    return dict(itype=itype, wmo=wmo, ident=ident, site=site, n=cn["n_kept"], dr=dr,
                dig_amp=dig["amp_fund"], dig_R2=dig["R2"], dig_fMHz=dig["f_MHz"],
                Lam_dom=Lam_dom, ac_dom=ac_dom, ripple_rms=rms, sig_ref=sr,
                dig_rel_pct=100 * dig["amp_fund"] / sr, rel_ripple_pct=100 * rms / sr,
                tslope=tslope, tlo=trange[0], thi=trange[1],
                _rng=rng, _hp=hp, _tamps=tamps)


rows = []
for itype in ("CL61", "CL51", "CL31"):
    for wmo, ident, _nd, site in SAMPLE[itype]:
        print(f"[{itype}] {wmo}_{ident} {site} ...", flush=True)
        try:
            r = analyse(itype, wmo, ident, site)
        except Exception as exc:
            print(f"   FAILED: {exc}"); continue
        if r is None:
            print("   skipped (too few clear nights)"); continue
        rows.append(r)
        print(f"   n={r['n']} dr={r['dr']:.0f}m  fs/4 ripple {r['dig_rel_pct']:.1f}% (R2={r['dig_R2']:.2f})  "
              f"dom Lambda={r['Lam_dom']:.0f}m (ac={r['ac_dom']:.2f})  relrms={r['rel_ripple_pct']:.1f}%  "
              f"Tslope={r['tslope']:+.3f}/C", flush=True)

# ---------- coefficient table ----------
cols = ["itype", "site", "wmo", "ident", "n_profiles", "range_gate_m", "digitizer_ripple_amp",
        "digitizer_ripple_pct", "digitizer_ripple_R2", "digitizer_f_MHz", "dominant_period_m",
        "dominant_autocorr", "ripple_rms_2_8km", "rel_ripple_pct", "temp_slope_per_C",
        "temp_lo_C", "temp_hi_C"]
with open(D / "network_offset_coeffs.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(cols)
    for r in rows:
        w.writerow([r["itype"], r["site"], r["wmo"], r["ident"], r["n"], f"{r['dr']:.0f}",
                    f"{r['dig_amp']:.4g}", f"{r['dig_rel_pct']:.2f}", f"{r['dig_R2']:.3f}",
                    f"{r['dig_fMHz']:.3f}", f"{r['Lam_dom']:.0f}", f"{r['ac_dom']:.3f}",
                    f"{r['ripple_rms']:.4g}", f"{r['rel_ripple_pct']:.2f}",
                    f"{r['tslope']:.4g}", f"{r['tlo']:.1f}", f"{r['thi']:.1f}"])
Path(D / "network_offset_coeffs.json").write_text(json.dumps(
    [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("_")}
     for r in rows], indent=1))
print(f"\nwrote {len(rows)} instruments -> network_offset_coeffs.csv/json")

# ---------- summary figure ----------
COL = {"CL61": "#d62728", "CL51": "#1f77b4", "CL31": "#2ca02c"}
fig, ax = plt.subplots(2, 2, figsize=(16, 10))
# (a) digitizer (fs/4) ripple amplitude per instrument, grouped by type
a = ax[0][0]; x = 0; ticks = []; labs = []
for itype in ("CL61", "CL51", "CL31"):
    grp = [r for r in rows if r["itype"] == itype]
    for r in grp:
        a.bar(x, r["dig_rel_pct"], color=COL[itype]); ticks.append(x); labs.append(r["site"][:10]); x += 1
    x += 1
a.set_xticks(ticks); a.set_xticklabels(labs, rotation=90, fontsize=6.5)
a.set_ylabel(r"$f_s/4$ ripple [% of mid-range signal]")
a.set_title("(a) digitizer (f$_s$/4) ripple — strong in CL51, weak/absent in CL61 & CL31")
a.grid(alpha=0.3, axis="y")
for it in ("CL61", "CL51", "CL31"):
    a.bar(np.nan, np.nan, color=COL[it], label=it)
a.legend(fontsize=9)
# (b) ripple strength vs coherence (autocorrelation)
b = ax[0][1]
for itype in ("CL61", "CL51", "CL31"):
    grp = [r for r in rows if r["itype"] == itype]
    b.scatter([r["dig_rel_pct"] for r in grp], [r["ac_dom"] for r in grp], c=COL[itype], label=itype, s=45)
b.set_xlabel(r"$f_s/4$ ripple [%]"); b.set_ylabel("dominant autocorrelation")
b.set_title("(b) ripple strength vs coherence"); b.legend(fontsize=9); b.grid(alpha=0.3)
# (c) example ripple profiles, normalised (one per type, the strongest)
c = ax[1][0]
for itype in ("CL61", "CL51", "CL31"):
    grp = sorted([r for r in rows if r["itype"] == itype], key=lambda r: -r["dig_rel_pct"])
    if grp:
        r = grp[0]; m = (r["_rng"] >= 2500) & (r["_rng"] <= 3100)
        c.plot(r["_rng"][m], 100 * r["_hp"][m] / r["sig_ref"], color=COL[itype], lw=1.1,
               label=f"{itype}: {r['site'][:12]}")
c.axhline(0, color="0.7", lw=0.7); c.set_xlim(2500, 3100)
c.set_xlabel("range [m]"); c.set_ylabel("relative ripple [%]")
c.set_title("(c) example recoverable ripple (strongest per type)"); c.legend(fontsize=8); c.grid(alpha=0.3)
# (d) temperature slope per instrument
d = ax[1][1]; x = 0; ticks = []; labs = []
for itype in ("CL61", "CL51", "CL31"):
    for r in [r for r in rows if r["itype"] == itype]:
        if np.isfinite(r["tslope"]):
            d.bar(x, r["tslope"], color=COL[itype]); ticks.append(x); labs.append(r["site"][:10]); x += 1
    x += 1
d.axhline(0, color="0.5", lw=0.8); d.set_xticks(ticks); d.set_xticklabels(labs, rotation=90, fontsize=6.5)
d.set_ylabel("ripple RMS slope per °C"); d.set_title("(d) internal-temperature dependence of the ripple")
d.grid(alpha=0.3, axis="y")
fig.suptitle("Network electronic-offset scan — 10 CL61 + 10 CL51 + 11 CL31 (clear-night method): "
             "the f$_s$/4 digitizer ripple is a CL51 signature; CL61/CL31 far-ripple is weak",
             fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
FIG = D / "fig_network_offset_scan.png"
FIG_REPORT = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_network_offset_scan.png")
fig.savefig(FIG, dpi=150); fig.savefig(FIG_REPORT, dpi=150)
print("saved", FIG); print("NETWORK_OFFSET_SCAN_DONE")
