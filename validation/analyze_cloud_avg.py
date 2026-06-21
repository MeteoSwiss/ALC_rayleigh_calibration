"""
analyze_cloud_avg.py — optimal pre-averaging for the liquid-cloud calibration. From run_cloud_avg_sweep
(native L1, K6 gates, average_time_s swept), compute per averaging level: valid % (valid days / processed
days) and sigma_SD (robust successive-difference precision of the nightly coefficient), per instrument
type. Plots valid% and sigma_SD vs average_time_s and prints the optimum.
"""
from __future__ import annotations
import glob, json, sys
from collections import defaultdict
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
AD = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_avg")
META = {m["label"]: m for m in json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())}
LEVELS = [30, 60, 120, 300, 600, 900, 1200]
TYPES = ["CL31", "CL51", "CL61"]
TCOLOR = {"CL31": "#1f77b4", "CL51": "#ff7f0e", "CL61": "#2ca02c"}


def sd(c):
    c = np.asarray(c, float)
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) / abs(np.median(c)) * 100) \
        if c.size >= 4 and np.median(c) > 0 else float("nan")


def main():
    # per (type, level): list over streams of (valid%, sigma_SD)
    agg = {t: {a: {"v": [], "s": []} for a in LEVELS} for t in TYPES}
    nstream = defaultdict(int)
    for fp in glob.glob(str(AD / "avg_*.json")):
        lab = Path(fp).stem[len("avg_"):]
        t = META.get(lab, {}).get("group")
        if t not in TYPES:
            continue
        nstream[t] += 1
        days = json.loads(Path(fp).read_text())
        order = sorted(days)
        n = len(order)
        if not n:
            continue
        for a in LEVELS:
            key = str(a) if str(a) in days[order[0]] else a
            cls = [days[ds][key][0] for ds in order
                   if key in days[ds] and np.isfinite(days[ds][key][0]) and days[ds][key][0] > 0 and days[ds][key][1] >= 1]
            agg[t][a]["v"].append(100.0 * len(cls) / n)
            agg[t][a]["s"].append(sd(cls))

    def med(x):
        x = [v for v in x if np.isfinite(v)]
        return float(np.median(x)) if x else float("nan")

    summary = {t: {a: {"valid": med(agg[t][a]["v"]), "sigma_sd": med(agg[t][a]["s"]),
                       "n": nstream[t]} for a in LEVELS} for t in TYPES}
    (AD / "cloud_avg_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    for t in TYPES:
        v = [summary[t][a]["valid"] for a in LEVELS]
        s = [summary[t][a]["sigma_sd"] for a in LEVELS]
        axes[0].plot(LEVELS, v, "o-", color=TCOLOR[t], label=f"{t} (n={nstream[t]})")
        axes[1].plot(LEVELS, s, "o-", color=TCOLOR[t], label=t)
    axes[0].axvline(300, color="0.6", ls="--", lw=1, label="current default 300 s")
    axes[1].axvline(300, color="0.6", ls="--", lw=1)
    axes[0].set_xlabel("average_time_s (s)"); axes[0].set_ylabel("valid cloud calibrations (% of days)")
    axes[0].set_title("Yield vs averaging"); axes[0].set_xscale("log"); axes[0].grid(alpha=0.3); axes[0].legend()
    axes[1].set_xlabel("average_time_s (s)"); axes[1].set_ylabel("σ_SD (% of median C_L)")
    axes[1].set_title("Short-term variability vs averaging"); axes[1].set_xscale("log"); axes[1].grid(alpha=0.3); axes[1].legend()
    fig.suptitle("Optimal pre-averaging for liquid-cloud calibration (native L1, K6 gates)", fontsize=14)
    fig.tight_layout(); fig.savefig(AD / "fig_cloud_avg.png", dpi=130); plt.close(fig)

    print("avg_time_s | " + " | ".join(f"{t} valid/σ" for t in TYPES))
    for a in LEVELS:
        cells = [f"{summary[t][a]['valid']:.0f}/{summary[t][a]['sigma_sd']:.1f}".replace("nan", "-") for t in TYPES]
        print(f"{a:9d}s | " + " | ".join(cells))
    for t in TYPES:
        best = max(LEVELS, key=lambda a: (summary[t][a]["valid"] if np.isfinite(summary[t][a]["valid"]) else -1))
        print(f"  {t}: max-valid at average_time_s = {best}s ({summary[t][best]['valid']:.0f}%, σ_SD {summary[t][best]['sigma_sd']:.1f})")
    print("saved cloud_avg_summary.json + fig_cloud_avg.png")


if __name__ == "__main__":
    main()
