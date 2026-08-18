# -*- coding: utf-8 -*-
"""Post-process the network offset scan (npz cache) into the FINAL operational coefficient table +
summary figure. The scan (v2) decides correctability by SPLIT-HALF REPRODUCIBILITY (a fixed
instrumental pattern reproduces across independent halves of the nights; atmosphere/noise does not)
and stores a FULL-RANGE b_phys (gate-fold for phase-locked digitizer ripples; empirical free-
troposphere pattern for reproducible slow distortions). Here we refine the ripple period for
reporting, add the undamped ratio, merge the temperature slope, and emit
network_offset_coeffs_final.csv/json + fig_network_offset_coeffs.png."""
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
from validation.paper._offset_lib import C_LIGHT

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
NPZ = D / "network_offset"
FIT_LO, FIT_HI = 2000.0, 6000.0


def refine_period(hp, rng, lam0, dr):
    if not np.isfinite(lam0):
        return lam0
    m = (rng >= FIT_LO) & (rng <= FIT_HI); r = rng[m]; y = np.nan_to_num(hp[m])
    if r.size < 40:
        return lam0
    grid = np.arange(max(2 * dr, lam0 * 0.6), lam0 * 1.6, 0.2)
    def cost(L):
        A = np.c_[np.ones_like(r), np.cos(2 * np.pi * r / L), np.sin(2 * np.pi * r / L)]
        c, *_ = np.linalg.lstsq(A, y, rcond=None); return np.sum((A @ c - y) ** 2)
    return float(grid[int(np.argmin([cost(L) for L in grid]))])


# temperature slope from the scan CSV
tslope = {}
sc = D / "network_offset_scan.csv"
if sc.exists():
    for row in csv.DictReader(open(sc, encoding="utf-8")):
        tslope[(row["wmo"], row["ident"], row["itype"])] = row.get("temp_slope_per_C", "nan")

rows = []
for f in sorted(glob.glob(str(NPZ / "*.npz"))):
    z = np.load(f, allow_pickle=True)
    name = Path(f).stem
    itype = name.split("_")[-1]; ident = name.split("_")[-2]; wmo = "_".join(name.split("_")[:-2])
    rng = z["rng"]; hp = z["hp"]; dr = float(np.median(np.diff(rng)))
    Lam0 = float(z["Lam_dom"]) if np.isfinite(z["Lam_dom"]) else np.nan
    Lam = refine_period(hp, rng, Lam0, dr)
    repro = float(z["repro"]); rel = float(z["rel"])
    gate = bool(z["is_gate_ripple"]); N = int(z["N"]); sig = float(z["sig_ref"])

    def rms(a, b):
        s = hp[(rng >= a) & (rng <= b)]
        return float(np.sqrt(np.nanmean(s ** 2))) if np.isfinite(s).any() else np.nan
    undamp = rms(6000, 10000) / rms(2000, 4000) if rms(2000, 4000) else np.nan
    # FINAL correctability: a real fixed offset the clear-night method can recover, requiring
    #  (i)  REPRODUCIBLE across independent halves of the nights (repro>=0.6) -> a fixed pattern, not noise;
    #  (ii) SIGNIFICANT (>=1.5% of the mid-range signal) -> worth correcting;
    #  (iii) NOT decreasing-with-range like atmosphere (undamp>=0.5) -> an electronic offset stays ~flat
    #        in P=rcs_0/z^2, whereas a persistent aerosol layer decays (undamp->0). This rejects e.g.
    #        Athens (undamp 0.14, a near-range atmospheric feature) while keeping the digitizer ripples.
    corr = bool(repro >= 0.6 and rel >= 1.5 and np.isfinite(undamp) and undamp >= 0.5)
    mode = f"gate-fold N={N}" if (gate and corr) else ("empirical-FT" if corr else "none")
    ts = tslope.get((wmo, ident, itype), "nan")
    rows.append(dict(itype=itype, site=str(z["site"]), wmo=wmo, ident=ident, dr=dr,
                     Lam=Lam, gates=Lam / dr if np.isfinite(Lam) else np.nan,
                     f_MHz=C_LIGHT / (2 * Lam) / 1e6 if np.isfinite(Lam) else np.nan,
                     rel=rel, repro=repro, undamp=undamp, mode=mode, correctable=corr,
                     tslope=float(ts) if ts not in ("nan", "") else np.nan, sig=sig,
                     _rng=rng, _hp=hp, _hpa=z["hpa"], _hpb=z["hpb"], _bphys=z["b_phys"]))

rows.sort(key=lambda r: ({"CL51": 0, "CL31": 1, "CL61": 2}[r["itype"]], -int(r["correctable"]), -np.nan_to_num(r["rel"])))
COLS = ["itype", "site", "wmo", "ident", "range_gate_m", "ripple_period_m", "ripple_gates",
        "ripple_f_MHz", "rel_amp_pct", "split_half_repro", "undamped_ratio", "correction_mode",
        "temp_slope_per_C", "correctable"]
with open(D / "network_offset_coeffs_final.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(COLS)
    for r in rows:
        w.writerow([r["itype"], r["site"], r["wmo"], r["ident"], f"{r['dr']:.0f}", f"{r['Lam']:.1f}",
                    f"{r['gates']:.2f}", f"{r['f_MHz']:.3f}", f"{r['rel']:.2f}", f"{r['repro']:.3f}",
                    f"{r['undamp']:.2f}", r["mode"], f"{r['tslope']:.4g}", int(r["correctable"])])
Path(D / "network_offset_coeffs_final.json").write_text(json.dumps(
    [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("_")}
     for r in rows], indent=1, default=str))

nc = {t: (sum(1 for r in rows if r["itype"] == t and r["correctable"]),
          sum(1 for r in rows if r["itype"] == t)) for t in ("CL51", "CL31", "CL61")}
print("=" * 96)
print("NETWORK OFFSET COEFFICIENTS (split-half reproducibility; full-range correction)")
print("=" * 96)
for r in rows:
    print(f"  {r['itype']:5} {r['site'][:22]:22} Λ={r['Lam']:6.1f}m amp={r['rel']:5.1f}% repro={r['repro']:+.2f}"
          f"  undamp={r['undamp']:4.1f}  [{r['mode']:14}]  {'✔ CORRECTABLE' if r['correctable'] else '-'}")
print(f"\n  correctable: CL51 {nc['CL51'][0]}/{nc['CL51'][1]}, CL31 {nc['CL31'][0]}/{nc['CL31'][1]}, "
      f"CL61 {nc['CL61'][0]}/{nc['CL61'][1]}")

# ---------- summary figure ----------
COL = {"CL61": "#d62728", "CL51": "#1f77b4", "CL31": "#2ca02c"}
fig, ax = plt.subplots(2, 2, figsize=(16, 10))
a = ax[0][0]; x = 0; ticks = []; labs = []
for it in ("CL51", "CL31", "CL61"):
    for r in [r for r in rows if r["itype"] == it]:
        col = COL[it] if r["correctable"] else "0.78"
        a.bar(x, r["rel"], color=col, edgecolor=COL[it], lw=1.1); ticks.append(x); labs.append(r["site"][:11]); x += 1
    x += 1
a.set_xticks(ticks); a.set_xticklabels(labs, rotation=90, fontsize=6.3)
a.set_ylabel("offset amplitude [% of mid-range signal]")
a.set_title("(a) recoverable offset amplitude (solid = correctable)"); a.grid(alpha=0.3, axis="y")
b = ax[0][1]
for it in ("CL51", "CL31", "CL61"):
    grp = [r for r in rows if r["itype"] == it]
    b.scatter([r["repro"] for r in grp], [r["rel"] for r in grp], c=COL[it], s=45, label=it)
b.axvline(0.6, color="k", ls=":", lw=1); b.axhline(1.5, color="k", ls=":", lw=1)
b.set_xlabel("split-half reproducibility"); b.set_ylabel("offset amplitude [%]")
b.set_title("(b) correctable = reproducible (≥0.6) AND ≥1.5% AND not decaying-with-range")
b.legend(fontsize=9); b.grid(alpha=0.3)
c = ax[1][0]
for it in ("CL51", "CL31", "CL61"):
    grp = sorted([r for r in rows if r["itype"] == it and r["correctable"]], key=lambda r: -r["rel"])
    if grp and np.isfinite(grp[0]["sig"]):
        r = grp[0]; m = (r["_rng"] >= 2500) & (r["_rng"] <= 3100)
        c.plot(r["_rng"][m], 100 * r["_hp"][m] / r["sig"], color=COL[it], lw=1.1, label=f"{it}: {r['site'][:12]}")
c.axhline(0, color="0.7", lw=0.7); c.set_xlim(2500, 3100)
c.set_xlabel("range [m]"); c.set_ylabel("relative offset [%]")
c.set_title("(c) strongest correctable offset per type"); c.legend(fontsize=8); c.grid(alpha=0.3)
d = ax[1][1]; x = 0; ticks = []; labs = []
for it in ("CL51", "CL31", "CL61"):
    for r in [r for r in rows if r["itype"] == it and np.isfinite(r["tslope"])]:
        d.bar(x, r["tslope"], color=COL[it]); ticks.append(x); labs.append(r["site"][:11]); x += 1
    x += 1
d.axhline(0, color="0.5", lw=0.8); d.set_xticks(ticks); d.set_xticklabels(labs, rotation=90, fontsize=6.3)
d.set_ylabel("offset RMS slope [per °C]"); d.set_title("(d) internal-temperature dependence (mostly weak)")
d.grid(alpha=0.3, axis="y")
fig.suptitle("Network electronic-offset coefficients — correctability by split-half reproducibility, "
             "full-range correction (gate-fold ripple / empirical distortion)", fontweight="bold", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
FIGR = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report/fig_network_offset_coeffs.png")
fig.savefig(D / "fig_network_offset_coeffs.png", dpi=150); fig.savefig(FIGR, dpi=150)
print("saved fig_network_offset_coeffs.png"); print("NETWORK_COEFFS_DONE")
