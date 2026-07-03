# -*- coding: utf-8 -*-
"""Per-station offset diagnostics for ALL scanned instruments, so the `correctable` verdict is
auditable. One panel per station overlays the clear-night high-pass offset from the TWO INDEPENDENT
HALVES of the nights (even / odd profiles): if they agree, the offset is a real FIXED instrumental
pattern (correctable); if they diverge, it is atmosphere/noise (not). Title carries Λ, amplitude,
split-half reproducibility and the verdict; panels are framed green (correctable) or grey (not). An
inset shows the FULL 0–8 km offset so large, non-periodic distortions are visible too."""
import sys, glob, warnings
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

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
NPZ = D / "network_offset"
REPORT_FIG = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report")
ZOOM_LO, ZOOM_HI = 2000.0, 3600.0


import csv
# FINAL verdict (with the undamped/instrumental gate) from the coefficient table
_final = {}
_cf = D / "network_offset_coeffs_final.csv"
if _cf.exists():
    for r in csv.DictReader(open(_cf, encoding="utf-8")):
        _final[(r["wmo"], r["ident"], r["itype"])] = (r["correctable"] == "1", r["correction_mode"])


def load(itype):
    out = []
    for f in sorted(glob.glob(str(NPZ / f"*_{itype}.npz"))):
        z = np.load(f, allow_pickle=True)
        name = Path(f).stem
        wmo = "_".join(name.split("_")[:-2]); ident = name.split("_")[-2]
        corr, cmode = _final.get((wmo, ident, itype), (bool(z["correctable"]), ""))
        out.append(dict(site=str(z["site"]), rng=z["rng"], hp=z["hp"], hpa=z["hpa"], hpb=z["hpb"],
                        repro=float(z["repro"]), rel=float(z["rel"]), Lam=float(z["Lam_dom"]),
                        N=int(z["N"]), corr=corr, gate="gate-fold" in cmode))
    out.sort(key=lambda d: (-int(d["corr"]), -np.nan_to_num(d["rel"])))
    return out


for itype in ("CL51", "CL31", "CL61"):
    rows = load(itype)
    n = len(rows); ncol = 3; nrow = int(np.ceil(n / ncol))
    fig, ax = plt.subplots(nrow, ncol, figsize=(16, 3.0 * nrow), squeeze=False)
    for k, d in enumerate(rows):
        a = ax[k // ncol][k % ncol]
        rng = d["rng"]; zm = (rng >= ZOOM_LO) & (rng <= ZOOM_HI)
        a.plot(rng[zm], d["hpa"][zm], "-", lw=0.8, color="#1f77b4", label="half A")
        a.plot(rng[zm], d["hpb"][zm], "-", lw=0.8, color="#ff7f0e", label="half B")
        a.axhline(0, color="0.7", lw=0.6)
        col = "#2ca02c" if d["corr"] else "0.6"
        for s in a.spines.values():
            s.set_color(col); s.set_linewidth(2.4 if d["corr"] else 1.0)
        mode = f"fold N={d['N']}" if d["gate"] else ("emp-FT" if d["corr"] else "")
        a.set_title(f"{d['site'][:20]}  Λ={d['Lam']:.0f}m amp={d['rel']:.1f}% repro={d['repro']:+.2f}"
                    + (f"  ✔{mode}" if d["corr"] else ""), fontsize=8.2,
                    color=("#1a7a1a" if d["corr"] else "0.25"))
        a.set_xlim(ZOOM_LO, ZOOM_HI); a.tick_params(labelsize=7)
        if k == 0:
            a.legend(fontsize=7, loc="upper right")
        ins = a.inset_axes([0.62, 0.62, 0.36, 0.36])
        fm = rng <= 8000
        ins.plot(d["hp"][fm], rng[fm] / 1000, lw=0.5, color="0.4")
        ins.axvline(0, color="0.7", lw=0.5); ins.set_xticks([]); ins.tick_params(labelsize=5.5); ins.set_ylim(0, 8)
    for k in range(n, nrow * ncol):
        ax[k // ncol][k % ncol].axis("off")
    ncorr = sum(1 for d in rows if d["corr"])
    fig.suptitle(f"{itype} offset diagnostics — {ncorr}/{n} correctable. Two colours = two independent "
                 f"halves of the nights; they COINCIDE for a real fixed pattern (green), DIVERGE for "
                 f"noise/atmosphere (grey). Native high-pass {ZOOM_LO/1000:.1f}–{ZOOM_HI/1000:.1f} km; inset = 0–8 km.",
                 fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = D / f"fig_offset_diagnostics_{itype}.png"
    fig.savefig(out, dpi=140); fig.savefig(REPORT_FIG / out.name, dpi=140)
    print(f"saved {out.name}  ({ncorr}/{n} correctable)", flush=True)
print("OFFSET_DIAGNOSTICS_DONE")
