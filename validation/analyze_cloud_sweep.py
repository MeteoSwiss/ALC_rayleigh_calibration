"""
analyze_cloud_sweep.py — summarize the liquid-cloud calibration config sweep (run_cloud_sweep.py).
Per (level, type, config): number of valid cloud calibrations (days with cal_median>0 and >=1
qualifying in-cloud profile), the valid fraction (valid days / processed days), and sigma_SD (robust
successive-difference precision of the nightly cal_median, % of median). Also L1-vs-L2 comparison.

Usage: python analyze_cloud_sweep.py <phase>     (phase 1 = 10/type, 2 = full)
Outputs (figs_paper_validation/cloud_sweep/): cloud_summary_<phase>.json,
fig_cloud_pareto_<phase>.png, fig_cloud_validbars_<phase>.png.
"""
from __future__ import annotations
import glob
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
DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_sweep")
META = {m["label"]: m for m in json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())}
TYPES = ["CL31", "CL51", "CL61"]
LEVELS = ["L1", "L2"]
CONFIGS = ["K0_baseline", "K1_consec3", "K2_consist20", "K3_ratio0.1", "K4_cbh3000",
           "K5_tempcold", "K6_balanced", "K7_aggressive"]
CLABEL = {"K0_baseline": "K0 baseline", "K1_consec3": "K1 n_consec=3", "K2_consist20": "K2 consist≤20%",
          "K3_ratio0.1": "K3 ratio≤0.10", "K4_cbh3000": "K4 cbh/cal≤3 km", "K5_tempcold": "K5 temp≥-25°C",
          "K6_balanced": "K6 balanced", "K7_aggressive": "K7 aggressive"}
COLORS = {c: plt.cm.tab10(i) for i, c in enumerate(CONFIGS)}
I_CL, I_N = 0, 1


def sigma_sd(c):
    c = np.asarray(c, float)
    if c.size < 4:
        return np.nan
    m = np.median(c)
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) / abs(m) * 100) if m > 0 else np.nan


def load(phase):
    data = {}
    for fp in glob.glob(str(DIR / f"cloud_{phase}_*.json")):
        name = Path(fp).stem[len(f"cloud_{phase}_"):]      # <level>_<label>
        level, label = name.split("_", 1)
        try:
            data[(level, label)] = json.loads(Path(fp).read_text())
        except Exception:
            pass
    return data


def per_stream(data):
    rows = []
    for (level, label), days in data.items():
        g = META.get(label, {}).get("group")
        if g is None or not days:
            continue
        order = sorted(days)
        n_proc = len(order)
        rec = dict(level=level, label=label, group=g, n_proc=n_proc)
        for c in CONFIGS:
            cls = [days[d][c][I_CL] for d in order
                   if np.isfinite(days[d][c][I_CL]) and days[d][c][I_CL] > 0 and days[d][c][I_N] >= 1]
            rec[f"valid_{c}"] = len(cls)
            rec[f"vfrac_{c}"] = 100.0 * len(cls) / n_proc if n_proc else np.nan
            rec[f"sd_{c}"] = sigma_sd(cls)
        rows.append(rec)
    return rows


def aggregate(rows):
    agg = {lvl: {t: {} for t in TYPES} for lvl in LEVELS}
    for lvl in LEVELS:
        for t in TYPES:
            sub = [r for r in rows if r["level"] == lvl and r["group"] == t]
            for c in CONFIGS:
                vf = [r[f"vfrac_{c}"] for r in sub if np.isfinite(r[f"vfrac_{c}"])]
                sd = [r[f"sd_{c}"] for r in sub if np.isfinite(r[f"sd_{c}"])]
                vd = [r[f"valid_{c}"] for r in sub]
                agg[lvl][t][c] = dict(
                    vfrac=float(np.median(vf)) if vf else np.nan,
                    sigma_sd=float(np.median(sd)) if sd else np.nan,
                    valid_days_total=int(np.sum(vd)),
                    valid_days_median=float(np.median(vd)) if vd else 0.0,
                    n_streams=len(sub), n_sd=len(sd))
    return agg


def fig_pareto(agg, path, phase):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))
    for ax, t in zip(axes, TYPES):
        for lvl, mk in zip(LEVELS, ["o", "s"]):
            for c in CONFIGS:
                a = agg[lvl][t][c]
                if not (np.isfinite(a["vfrac"]) and np.isfinite(a["sigma_sd"])):
                    continue
                ax.scatter(a["vfrac"], a["sigma_sd"], s=150, marker=mk, color=COLORS[c],
                           edgecolor="k", linewidth=0.6, zorder=3)
                ax.annotate(c.split("_")[0], (a["vfrac"], a["sigma_sd"]), fontsize=8, ha="center", va="center")
        ax.set_title(t, fontsize=13, fontweight="bold")
        ax.set_xlabel("valid cloud calibrations (% of processed days)")
        ax.grid(alpha=0.3)
        ax.annotate("best ↘\n(more valid, lower σ)", xy=(0.97, 0.05), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=8, color="green",
                    bbox=dict(boxstyle="round", fc="honeydew", ec="green", alpha=0.7))
    axes[0].set_ylabel("σ_SD (% of median C) — short-term variability")
    h = [plt.Line2D([], [], marker="o", color="0.5", ls="", label="L1"),
         plt.Line2D([], [], marker="s", color="0.5", ls="", label="L2")]
    h += [plt.Line2D([], [], marker="o", color=COLORS[c], ls="", label=CLABEL[c]) for c in CONFIGS]
    axes[2].legend(handles=h, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.suptitle(f"Cloud-calibration configs — yield vs variability (phase {phase}; median over streams)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 0.93, 1]); fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_validbars(agg, path, phase):
    fig, axes = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    x = np.arange(len(CONFIGS)); w = 0.26
    for ax, lvl in zip(axes, LEVELS):
        for k, t in enumerate(TYPES):
            ax.bar(x + (k - 1) * w, [agg[lvl][t][c]["vfrac"] for c in CONFIGS], w, label=t)
        ax.set_xticks(x); ax.set_xticklabels([c.split("_")[0] for c in CONFIGS])
        ax.set_title(lvl, fontsize=13, fontweight="bold"); ax.grid(axis="y", alpha=0.3)
        ax.set_xlabel("cloud config")
    axes[0].set_ylabel("valid cloud calibrations (% of processed days)\nmedian over streams")
    axes[0].legend(title="type")
    fig.suptitle(f"Valid cloud-calibration fraction by config (phase {phase})", fontsize=14)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "1"
    data = load(phase)
    rows = per_stream(data)
    agg = aggregate(rows)
    (DIR / f"cloud_summary_{phase}.json").write_text(
        json.dumps({"aggregate": agg, "per_stream": rows}, indent=1), encoding="utf-8")
    fig_pareto(agg, DIR / f"fig_cloud_pareto_{phase}.png", phase)
    fig_validbars(agg, DIR / f"fig_cloud_validbars_{phase}.png", phase)
    print(f"phase {phase}: {len(data)} (level,stream) files, {len(rows)} usable streams")
    for lvl in LEVELS:
        print(f"\n=== {lvl} :: valid% / σ_SD%  (median over streams) ===")
        print(f"{'config':14s} " + " ".join(f"{t:>14s}" for t in TYPES))
        for c in CONFIGS:
            cells = []
            for t in TYPES:
                a = agg[lvl][t][c]
                cells.append(f"{a['vfrac']:4.0f}/{a['sigma_sd']:4.1f}".replace("nan", " - "))
            print(f"{c:14s} " + " ".join(f"{x:>14s}" for x in cells))
    print(f"\nsaved cloud_summary_{phase}.json + 2 figures")


if __name__ == "__main__":
    main()
