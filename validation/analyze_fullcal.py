#!/usr/bin/env python3
"""First-pass analysis of the all-station L2-monthly calibration (in progress).

Reads fullcal_all/summary.csv (deduped by WMO, last wins), prints per-type statistics
and per-station flag breakdown (sampled from the station CSVs), and writes plots to
fullcal_all/analysis/.
"""
import csv
import glob
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all")
OUT = BASE / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

# --- dedupe summary (last row per wmo wins) ---
rows = {}
for r in csv.DictReader(open(BASE / "summary.csv")):
    rows[r["wmo"]] = r
recs = []
for r in rows.values():
    try:
        recs.append(dict(wmo=r["wmo"], itype=r["itype"], n=int(r["n_dates"]),
                         ok=int(r["n_success"]), med=float(r["median_cl"])))
    except ValueError:
        pass

print(f"Stations analysed: {len(recs)}")
by_type = defaultdict(list)
for r in recs:
    by_type[r["itype"]].append(r)

print("\n=== Per-type summary ===")
print(f"{'type':10s}{'stns':>5}{'nights':>9}{'ok':>8}{'ok%':>6}   median-CL across stations [min, median, max]")
for t, rs in sorted(by_type.items()):
    nights = sum(r["n"] for r in rs)
    ok = sum(r["ok"] for r in rs)
    meds = np.array([r["med"] for r in rs if r["ok"] > 0 and np.isfinite(r["med"])])
    mstr = (f"[{meds.min():.2e}, {np.median(meds):.2e}, {meds.max():.2e}]"
            if len(meds) else "(none)")
    print(f"{t:10s}{len(rs):>5}{nights:>9}{ok:>8}{100*ok/max(nights,1):>5.1f}%   {mstr}")

# --- flag breakdown across all station CSVs done so far ---
flagmap = {1: "success", 0.5: "partial", 0: "no_data", -1: "not_clear",
           -2: "not_proportional", -3: "method_disagree", -4: "no_model",
           -5: "rcs_nan/pert_fail", -6: "unc>value", -7: "neg_slope", -8: "fit|b|>a", -99: "driver_error"}
flags = Counter()
for f in glob.glob(str(BASE / "*/*_cl.csv")):
    for r in csv.DictReader(open(f)):
        try:
            flags[float(r["flag"])] += 1
        except ValueError:
            pass
total = sum(flags.values())
print(f"\n=== Flag breakdown ({total} night-calibrations) ===")
for fl, n in sorted(flags.items(), key=lambda x: -x[1]):
    print(f"  {flagmap.get(fl, fl):20s} {n:7d}  {100*n/total:4.1f}%")

# --- plots ---
chm = [r for r in by_type.get("CHM15k", []) if r["ok"] > 0 and np.isfinite(r["med"])]
if chm:
    meds = np.array([r["med"] for r in chm])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    ax[0].hist(meds, bins=np.linspace(0, 6e11, 40), color="tab:blue", alpha=0.8)
    ax[0].axvline(np.median(meds), color="k", ls="--", label=f"median {np.median(meds):.2e}")
    ax[0].set_xlabel("per-station median C_L"); ax[0].set_ylabel("# stations")
    ax[0].set_title(f"CHM15k median C_L ({len(chm)} stations)"); ax[0].legend()
    rates = np.array([100 * r["ok"] / r["n"] for r in chm])
    ax[1].hist(rates, bins=np.linspace(0, 60, 30), color="tab:green", alpha=0.8)
    ax[1].axvline(np.median(rates), color="k", ls="--", label=f"median {np.median(rates):.0f}%")
    ax[1].set_xlabel("success rate [%]"); ax[1].set_ylabel("# stations")
    ax[1].set_title("CHM15k per-station success rate"); ax[1].legend()
    fig.tight_layout(); fig.savefig(OUT / "chm15k_distributions.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\nsaved {OUT/'chm15k_distributions.png'}")

# notable CHM15k outliers
if chm:
    chm_sorted = sorted(chm, key=lambda r: r["med"])
    print("\nCHM15k lowest median C_L:", [(r["wmo"], f"{r['med']:.2e}") for r in chm_sorted[:3]])
    print("CHM15k highest median C_L:", [(r["wmo"], f"{r['med']:.2e}") for r in chm_sorted[-3:]])
