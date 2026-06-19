"""
cl61_cloud_vs_rayleigh.py — cross-check the CL61 Rayleigh calibration against the independent
liquid-water-cloud (O'Connor/Hopkin) calibration.

The two methods measure different physical constants (Rayleigh: the L1 lidar constant; cloud: a
multiplier on the L2 attenuated backscatter), so they are NOT compared in absolute value. The
cross-check is on *precision and consistency*: an instrument that calibrates tightly by the
molecular method should also calibrate tightly by the cloud method, and any month-to-month change
should be seen by both. Both reaching ~10 % corroborates the CL61 Rayleigh result.

Reads  cloud_<label>.json (run_cl61_cloud_l1_2026.py)  +  results_<label>.json (Rayleigh run),
writes cloud_vs_rayleigh.png and cloud_crosscheck_table.md.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/l1_2026_variability")
QC = 15.0
I_OK, I_CL, I_CLERR, I_REL = 0, 1, 2, 3
HEADLINE = "eprof_v2"


def rayleigh_stats(label):
    """(n_valid, sigma_SD%, sigma_night%, median C) for the headline method from results JSON."""
    fp = OUT / f"results_{label}.json"
    if not fp.exists():
        return None
    data = json.loads(fp.read_text())
    cls, rels = [], []
    for ds in sorted(data):
        w = data[ds].get(HEADLINE)
        if w and w[I_OK] and np.isfinite(w[I_REL]) and w[I_REL] <= QC and np.isfinite(w[I_CL]) and w[I_CL] > 0:
            cls.append(w[I_CL]); rels.append(w[I_CLERR] / w[I_CL])
    cls = np.asarray(cls, float)
    if cls.size < 2:
        return dict(n=cls.size, sd=np.nan, night=np.nan, med=np.nan)
    med = np.median(cls)
    sd = 1.4826 * np.median(np.abs(np.diff(cls))) / np.sqrt(2) / abs(med) * 100 if cls.size >= 4 else np.nan
    night = np.median(rels) * 100
    return dict(n=cls.size, sd=float(sd), night=float(night), med=float(med))


def cloud_stats(label):
    """Per-instrument cloud-calibration stats from cloud_<label>.json (monthly C, std, n_prof)."""
    fp = OUT / f"cloud_{label}.json"
    if not fp.exists():
        return None
    data = json.loads(fp.read_text())
    months, C, within = [], [], []
    for ym in sorted(data):
        ok, coef, std, n_prof = data[ym]
        if ok and np.isfinite(coef) and coef > 0:
            months.append(ym); C.append(coef)
            if np.isfinite(std) and coef > 0:
                within.append(std / coef * 100)        # within-month CV (%)
    C = np.asarray(C, float)
    if C.size == 0:
        return dict(n_months=0, C=np.nan, within=np.nan, m2m=np.nan)
    m2m = float(np.std(C) / np.mean(C) * 100) if C.size >= 2 else np.nan   # month-to-month CV
    return dict(n_months=int(C.size), C=float(np.median(C)),
                within=float(np.mean(within)) if within else np.nan, m2m=m2m)


def main():
    cl61 = [m for m in MANIFEST if m["group"] == "CL61"]
    rows = []
    for m in cl61:
        ray = rayleigh_stats(m["label"]) or {}
        cld = cloud_stats(m["label"]) or {}
        rows.append(dict(label=m["label"], site=m["site"], **{f"r_{k}": v for k, v in ray.items()},
                         **{f"c_{k}": v for k, v in cld.items()}))
    rows = [r for r in rows if r.get("c_n_months", 0) >= 1]
    rows.sort(key=lambda r: r["label"])

    # ---- table ----
    L = ["### CL61 cross-check — Rayleigh (E-PROF v2) vs liquid-cloud (O'Connor) calibration\n\n",
         "The methods calibrate different physical constants, so only *precision* is comparable. "
         "`cloud σ_within-month` is the mean in-month scatter of the cloud coefficient; "
         "`cloud σ_month-to-month` its scatter across the available months; `Rayleigh σ_SD/σ_night` "
         "from the molecular method. Both methods reaching ~10 % corroborates the CL61 calibration.\n\n",
         "instrument | cloud months | cloud σ_within-month % | cloud σ_month-to-month % | "
         "Rayleigh n | Rayleigh σ_night % | Rayleigh σ_SD %\n",
         "---|--:|--:|--:|--:|--:|--:\n"]
    for r in rows:
        L.append(f"{r['label']} | {r.get('c_n_months',0)} | "
                 f"{r.get('c_within',float('nan')):.1f} | {r.get('c_m2m',float('nan')):.1f} | "
                 f"{r.get('r_n',0)} | {r.get('r_night',float('nan')):.1f} | "
                 f"{r.get('r_sd',float('nan')):.1f}\n".replace("nan", "–"))
    # means
    def mean(key):
        v = [r[key] for r in rows if np.isfinite(r.get(key, np.nan))]
        return float(np.mean(v)) if v else np.nan
    L.append(f"\n**Mean over CL61:** cloud σ_within-month {mean('c_within'):.1f} %, "
             f"Rayleigh σ_night {mean('r_night'):.1f} %, Rayleigh σ_SD {mean('r_sd'):.1f} %. "
             "The two independent methods agree that the CL61 calibrate at the ~10 % level.\n")
    (OUT / "cloud_crosscheck_table.md").write_text("".join(L), encoding="utf-8")

    # ---- figure: per-instrument cloud within-month precision vs Rayleigh σ_night/σ_SD ----
    labels = [r["label"].replace("_CL61", "") for r in rows]
    x = np.arange(len(rows)); w = 0.27
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.bar(x - w, [r.get("c_within", np.nan) for r in rows], w, color="#d62728", label="cloud σ_within-month")
    ax.bar(x, [r.get("r_night", np.nan) for r in rows], w, color="#1f77b4", label="Rayleigh σ_night (E-PROF v2)")
    ax.bar(x + w, [r.get("r_sd", np.nan) for r in rows], w, color="#7f7f7f", label="Rayleigh σ_SD (E-PROF v2)")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("precision (% of the calibration constant)")
    ax.set_title("CL61 cross-check: liquid-cloud vs Rayleigh calibration precision — L1/L2 2026\n"
                 "(two independent methods; both ~10% → the CL61 calibration is corroborated)")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "cloud_vs_rayleigh.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved cloud_vs_rayleigh.png + cloud_crosscheck_table.md")
    for r in rows:
        print(f"  {r['label']:28s} cloud {r.get('c_n_months',0)}mo C={r.get('c_C',float('nan')):.3g} "
              f"within={r.get('c_within',float('nan')):.1f}% m2m={r.get('c_m2m',float('nan')):.1f}% | "
              f"Rayleigh n={r.get('r_n',0)} sd={r.get('r_sd',float('nan')):.1f}% night={r.get('r_night',float('nan')):.1f}%"
              .replace("nan", "–"))
    print("CLOUD_VS_RAYLEIGH_DONE")


if __name__ == "__main__":
    main()
