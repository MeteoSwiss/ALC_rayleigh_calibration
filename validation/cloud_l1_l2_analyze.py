"""
cloud_l1_l2_analyze.py — summarize cloud_l1_l2_native.py. Per (type, variant, config): valid %
(valid days / days the variant was attempted) and sigma_SD (robust successive-difference precision,
% of median; convention-independent, labelled C_L). Variants: L1_native, L1_300s, L2. Configs:
K0/K6/K7. Produces a grouped-bar figure (yield) + sigma_SD table and prints the key contrasts:
 - added value of the config change (K0 -> K6/K7)
 - added value of native L1 vs L2.
"""
from __future__ import annotations
import glob, json, sys
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
LL = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_l1l2")
META = {m["label"]: m for m in json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())}
TYPES = ["CL31", "CL51", "CL61"]
VARIANTS = ["L1_native", "L1_300s", "L2"]
CONFIGS = ["K0", "K6", "K7"]
CCOLOR = {"K0": "#7f7f7f", "K6": "#1f77b4", "K7": "#d62728"}


def sd(c):
    c = np.asarray(c, float)
    if c.size < 4:
        return np.nan
    m = np.median(c)
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) / abs(m) * 100) if m > 0 else np.nan


def main():
    # agg[type][variant][config] = {"v":[valid% per stream], "s":[sigma_SD per stream]}
    agg = {t: {v: {c: {"v": [], "s": []} for c in CONFIGS} for v in VARIANTS} for t in TYPES}
    nstream = {t: 0 for t in TYPES}
    for fp in glob.glob(str(LL / "ll_*.json")):
        lab = Path(fp).stem[len("ll_"):]
        t = META.get(lab, {}).get("group")
        if t not in TYPES:
            continue
        days = json.loads(Path(fp).read_text())
        if not days:
            continue
        nstream[t] += 1
        order = sorted(days)
        for v in VARIANTS:
            dv = [d for d in order if v in days[d]]
            n_proc = len(dv)
            if not n_proc:
                for c in CONFIGS:
                    agg[t][v][c]["v"].append(np.nan); agg[t][v][c]["s"].append(np.nan)
                continue
            for c in CONFIGS:
                vals = [days[d][v][c][0] for d in dv
                        if np.isfinite(days[d][v][c][0]) and days[d][v][c][0] > 0 and days[d][v][c][1] >= 1]
                agg[t][v][c]["v"].append(100.0 * len(vals) / n_proc)
                agg[t][v][c]["s"].append(sd(vals))

    def med(x):
        x = [v for v in x if np.isfinite(v)]
        return float(np.median(x)) if x else float("nan")

    summary = {t: {v: {c: {"valid": med(agg[t][v][c]["v"]), "sigma_sd": med(agg[t][v][c]["s"])}
                       for c in CONFIGS} for v in VARIANTS} for t in TYPES}
    summary["_n_streams"] = nstream
    (LL / "cloud_l1l2_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")

    # --- figure: 1 row, 3 cols (per type); grouped bars variant x config (valid%) ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.2), sharey=True)
    x = np.arange(len(VARIANTS)); w = 0.26
    for ax, t in zip(axes, TYPES):
        for k, c in enumerate(CONFIGS):
            vals = [summary[t][v][c]["valid"] for v in VARIANTS]
            bars = ax.bar(x + (k - 1) * w, vals, w, label=c, color=CCOLOR[c])
            for b, vv, v in zip(bars, vals, VARIANTS):
                s = summary[t][v][c]["sigma_sd"]
                if np.isfinite(vv):
                    ax.annotate(f"{vv:.0f}\nσ{s:.0f}" if np.isfinite(s) else f"{vv:.0f}",
                                (b.get_x() + b.get_width() / 2, vv), ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels(VARIANTS, fontsize=9)
        ax.set_title(f"{t} (n={nstream[t]})", fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("valid cloud calibrations (% of attempted days)\nbar label: valid% / σ_SD%")
    axes[0].legend(title="config")
    fig.suptitle("Liquid-cloud calibration: L1 native vs L1@300 s vs L2, baseline K0 vs K6/K7 "
                 f"({sum(nstream.values())} streams)", fontsize=14)
    fig.tight_layout(); fig.savefig(LL / "fig_cloud_l1l2.png", dpi=130); plt.close(fig)

    # --- text ---
    print(f"streams: " + ", ".join(f"{t}={nstream[t]}" for t in TYPES))
    for t in TYPES:
        print(f"\n=== {t} :: valid% / sigma_SD%  (median over streams) ===")
        print(f"{'variant':10s} " + " ".join(f"{c:>12s}" for c in CONFIGS))
        for v in VARIANTS:
            cells = [f"{summary[t][v][c]['valid']:4.0f}/{summary[t][v][c]['sigma_sd']:4.1f}".replace("nan", " - ")
                     for c in CONFIGS]
            print(f"{v:10s} " + " ".join(f"{x:>12s}" for x in cells))
    print("\nsaved cloud_l1l2_summary.json + fig_cloud_l1l2.png")


if __name__ == "__main__":
    main()
