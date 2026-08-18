"""
analyze_v2_sweep.py — turn run_v2_sweep.py output into (a) the clear-night failure diagnostic
(which v2 gate is the binding constraint per instrument type / level) and (b) the config trade-off
between valid-calibration fraction and short-term variability (sigma_SD).

sigma_SD is the robust successive-difference (von Neumann) precision used throughout this project:
    sigma_SD = 1.4826 * median(|diff(C)|) / sqrt(2) / |median(C)| * 100   [% of median C]

Outputs (under figs_paper_validation/v2_sweep/):
    v2_sweep_summary.json      machine-readable summary for the report + verification agents
    fig_v2_pareto.png          valid% vs sigma_SD per config, one panel per instrument type
    fig_v2_validbars.png       valid% per config, grouped by type, L1 vs L2
    fig_v2_bottleneck.png      leave-one-gate-out recovery of failed clear nights, per type/level
"""
from __future__ import annotations
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/v2_sweep")
MANIFEST = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
LABEL_GROUP = {i["label"]: i["group"] for i in MANIFEST}
TYPES = ["CHM15k", "Mini-MPL", "CL61"]
LEVELS = ["L1", "L2"]
I_OK, I_CL = 0, 1
QC = 15.0

OPT_CONFIGS = ["C0_baseline", "C1_tcv0.8", "C2_tcv1.2", "C3_shape",
               "C4_r2_0.35", "C5_start1200", "C6_balanced", "C7_aggressive", "C8_recommended"]
CONFIG_LABEL = {
    "C0_baseline": "C0 baseline (production v2)",
    "C1_tcv0.8": "C1 temporal_cv≤0.8",
    "C2_tcv1.2": "C2 temporal_cv≤1.2",
    "C3_shape": "C3 looser shape/ratio",
    "C4_r2_0.35": "C4 R²≥0.35",
    "C5_start1200": "C5 start≥1.2 km",
    "C6_balanced": "C6 balanced relax",
    "C7_aggressive": "C7 aggressive relax",
    "C8_recommended": "C8 RECOMMENDED (scatter≤1.15 + C6)",
}
LOO_GATES = {
    "LOO_temporal": "temporal_cv", "LOO_r2": "R²", "LOO_residual": "Rayleigh residual",
    "LOO_scattering": "scattering ratio", "LOO_ratiostd": "in-window ratio std",
    "LOO_start": "window start height",
}


def sigma_sd(cls):
    """Robust successive-difference precision, % of median; NaN if < 4 points."""
    cls = np.asarray(cls, float)
    if cls.size < 4:
        return np.nan
    med = np.median(cls)
    if med <= 0:
        return np.nan
    return float(1.4826 * np.median(np.abs(np.diff(cls))) / np.sqrt(2) / abs(med) * 100)


def load():
    """{(level,label): {ds: {config: [ok,cl,...]}}} for every sweep_*.json present."""
    data = {}
    for fp in glob.glob(str(DIR / "sweep_*.json")):
        name = Path(fp).stem            # sweep_L1_Payerne_CHM15k
        _, level, label = name.split("_", 2)
        data[(level, label)] = json.loads(Path(fp).read_text())
    return data


def per_instrument(data):
    """rows: one per (level,label,config) with valid_pct + sigma_SD."""
    rows = []
    for (level, label), nights in data.items():
        group = LABEL_GROUP.get(label)
        if group is None:
            continue
        order = sorted(nights)
        n_fit = len(order)
        for cfg in OPT_CONFIGS:
            cls = [nights[ds][cfg][I_CL] for ds in order if nights[ds][cfg][I_OK]]
            rows.append(dict(level=level, label=label, group=group, config=cfg,
                             n_fit=n_fit, n_valid=len(cls),
                             valid_pct=100.0 * len(cls) / n_fit if n_fit else np.nan,
                             sigma_sd=sigma_sd(cls)))
    return rows


def aggregate(rows):
    """median valid% and median sigma_SD across instruments, per (level,type,config)."""
    agg = {lvl: {t: {} for t in TYPES} for lvl in LEVELS}
    for lvl in LEVELS:
        for t in TYPES:
            for cfg in OPT_CONFIGS:
                sub = [r for r in rows if r["level"] == lvl and r["group"] == t and r["config"] == cfg]
                vp = [r["valid_pct"] for r in sub if np.isfinite(r["valid_pct"])]
                sd = [r["sigma_sd"] for r in sub if np.isfinite(r["sigma_sd"])]
                agg[lvl][t][cfg] = dict(
                    valid_pct=float(np.median(vp)) if vp else np.nan,
                    sigma_sd=float(np.median(sd)) if sd else np.nan,
                    n_inst=len(vp), n_inst_sd=len(sd),
                    fit_nights=int(np.sum([r["n_fit"] for r in sub])))
    return agg


def bottleneck(data):
    """Among baseline-FAILED clear nights, fraction recovered by relaxing ONE gate. Per level,type."""
    out = {lvl: {t: {} for t in TYPES} for lvl in LEVELS}
    for lvl in LEVELS:
        for t in TYPES:
            failed = 0
            rec = defaultdict(int)
            for (level, label), nights in data.items():
                if level != lvl or LABEL_GROUP.get(label) != t:
                    continue
                for ds, cfgs in nights.items():
                    if cfgs["C0_baseline"][I_OK]:
                        continue
                    failed += 1
                    for g in LOO_GATES:
                        if cfgs[g][I_OK]:
                            rec[g] += 1
            out[lvl][t] = dict(n_failed=failed,
                               recovery={LOO_GATES[g]: (100.0 * rec[g] / failed if failed else 0.0)
                                         for g in LOO_GATES})
    return out


# ───────────────────────────── figures ──────────────────────────────────────
COLORS = {c: plt.cm.tab10(i) for i, c in enumerate(OPT_CONFIGS)}


def fig_pareto(agg, path):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))
    for ax, t in zip(axes, TYPES):
        for lvl, mk in zip(LEVELS, ["o", "s"]):
            for cfg in OPT_CONFIGS:
                a = agg[lvl][t][cfg]
                if not (np.isfinite(a["valid_pct"]) and np.isfinite(a["sigma_sd"])):
                    continue
                ax.scatter(a["valid_pct"], a["sigma_sd"], s=150, marker=mk,
                           color=COLORS[cfg], edgecolor="k", linewidth=0.6, zorder=3,
                           label=None)
                ax.annotate(cfg.split("_")[0], (a["valid_pct"], a["sigma_sd"]),
                            fontsize=8, ha="center", va="center", zorder=4)
        ax.set_title(f"{t}", fontsize=13, fontweight="bold")
        ax.set_xlabel("valid calibrations on clear nights (%)")
        ax.grid(alpha=0.3)
        ax.annotate("best ↘\n(more valid, lower σ)", xy=(0.97, 0.05), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=8, color="green",
                    bbox=dict(boxstyle="round", fc="honeydew", ec="green", alpha=0.7))
    axes[0].set_ylabel("σ_SD  (% of median C_L)  — short-term variability")
    h = [plt.Line2D([], [], marker="o", color="0.5", ls="", label="L1"),
         plt.Line2D([], [], marker="s", color="0.5", ls="", label="L2")]
    h += [plt.Line2D([], [], marker="o", color=COLORS[c], ls="", label=CONFIG_LABEL[c])
          for c in OPT_CONFIGS]
    axes[2].legend(handles=h, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.suptitle("v2 gate configurations — yield vs short-term variability (median over instruments)",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 0.93, 1])
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_validbars(agg, path):
    fig, axes = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    x = np.arange(len(OPT_CONFIGS))
    w = 0.26
    for ax, lvl in zip(axes, LEVELS):
        for k, t in enumerate(TYPES):
            vals = [agg[lvl][t][c]["valid_pct"] for c in OPT_CONFIGS]
            ax.bar(x + (k - 1) * w, vals, w, label=t)
        ax.set_xticks(x)
        ax.set_xticklabels([c.split("_")[0] for c in OPT_CONFIGS], rotation=0)
        ax.set_title(f"{lvl}", fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        ax.set_xlabel("v2 configuration")
    axes[0].set_ylabel("valid calibrations on clear nights (%)\n(median over instruments)")
    axes[0].legend(title="instrument", fontsize=10)
    fig.suptitle("Valid-calibration fraction on clear nights by v2 configuration", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_bottleneck(bn, path):
    gates = list(LOO_GATES.values())
    fig, axes = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    x = np.arange(len(gates))
    w = 0.26
    for ax, lvl in zip(axes, LEVELS):
        for k, t in enumerate(TYPES):
            vals = [bn[lvl][t]["recovery"][g] for g in gates]
            ax.bar(x + (k - 1) * w, vals, w, label=f"{t} (n_fail={bn[lvl][t]['n_failed']})")
        ax.set_xticks(x)
        ax.set_xticklabels(gates, rotation=20, ha="right")
        ax.set_title(f"{lvl}", fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        ax.set_xlabel("gate relaxed (leave-one-out)")
    axes[0].set_ylabel("baseline-failed clear nights\nrecovered by relaxing ONLY this gate (%)")
    for ax in axes:
        ax.legend(fontsize=9)
    fig.suptitle("Why clear nights fail v2 — single-gate recovery of failed nights "
                 "(higher = more binding constraint)", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    data = load()
    print(f"loaded {len(data)} (level,instrument) sweep files")
    rows = per_instrument(data)
    agg = aggregate(rows)
    bn = bottleneck(data)
    summary = dict(aggregate=agg, bottleneck=bn,
                   per_instrument=rows, config_label=CONFIG_LABEL)
    (DIR / "v2_sweep_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    fig_pareto(agg, DIR / "fig_v2_pareto.png")
    fig_validbars(agg, DIR / "fig_v2_validbars.png")
    fig_bottleneck(bn, DIR / "fig_v2_bottleneck.png")

    # console digest
    for lvl in LEVELS:
        print(f"\n=== {lvl} :: valid% / sigma_SD%  (median over instruments) ===")
        print(f"{'config':14s} " + " ".join(f"{t:>16s}" for t in TYPES))
        for cfg in OPT_CONFIGS:
            cells = []
            for t in TYPES:
                a = agg[lvl][t][cfg]
                cells.append(f"{a['valid_pct']:5.1f}/{a['sigma_sd']:4.1f}".replace("nan", " - "))
            print(f"{cfg:14s} " + " ".join(f"{c:>16s}" for c in cells))
    print("\n=== binding gate (single-gate recovery of failed clear nights) ===")
    for lvl in LEVELS:
        for t in TYPES:
            r = bn[lvl][t]["recovery"]
            top = sorted(r.items(), key=lambda kv: -kv[1])[:2]
            print(f"  {lvl} {t:9s} n_fail={bn[lvl][t]['n_failed']:4d} -> "
                  + ", ".join(f"{g}:{v:.0f}%" for g, v in top))
    print("\nsaved v2_sweep_summary.json + 3 figures to", DIR)


if __name__ == "__main__":
    main()
