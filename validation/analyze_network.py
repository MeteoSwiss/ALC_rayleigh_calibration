"""
analyze_network.py — summarize the network-wide v2(C8)-vs-v1.1 comparison (run_network_v2_vs_v11.py).
Per stream and per level: valid% (valid calibrations / clear nights) and sigma_SD (robust
successive-difference precision, % of median C). Aggregated per instrument type, with PAIRED
per-stream deltas (v2 - v1.1). Landscape figures + a summary JSON + a markdown report fragment.
"""
from __future__ import annotations
import glob
import json
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/network_v2_v11")
MANIFEST = json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
LABEL_GROUP = {m["label"]: m["group"] for m in MANIFEST}
TYPES = ["CHM15k", "Mini-MPL", "CL61"]
LEVELS = ["L1", "L2"]
METHODS = ["v2", "v1.1"]
COLOR = {"v2": "#1f77b4", "v1.1": "#ff7f0e"}
TCOLOR = {"CHM15k": "#1f77b4", "Mini-MPL": "#2ca02c", "CL61": "#d62728"}
I_OK, I_CL = 0, 1


def sigma_sd(cls):
    cls = np.asarray(cls, float)
    if cls.size < 4:
        return np.nan
    med = np.median(cls)
    if med <= 0:
        return np.nan
    return float(1.4826 * np.median(np.abs(np.diff(cls))) / np.sqrt(2) / abs(med) * 100)


def load():
    data = {}
    for fp in glob.glob(str(DIR / "net_*.json")):
        name = Path(fp).stem            # net_L1_<wmo>_<id>
        _, level, label = name.split("_", 2)
        try:
            data[(level, label)] = json.loads(Path(fp).read_text())
        except Exception:
            pass
    return data


def per_stream(data):
    rows = []
    for (level, label), nights in data.items():
        group = LABEL_GROUP.get(label)
        if group is None or not nights:
            continue
        order = sorted(nights)
        n_fit = len(order)
        rec = dict(level=level, label=label, group=group, n_fit=n_fit)
        for m in METHODS:
            cls = [nights[ds][m][I_CL] for ds in order if nights[ds][m][I_OK]]
            rec[f"valid_{m}"] = 100.0 * len(cls) / n_fit if n_fit else np.nan
            rec[f"sd_{m}"] = sigma_sd(cls)
            rec[f"n_{m}"] = len(cls)
        rows.append(rec)
    return rows


def aggregate(rows):
    agg = {lvl: {t: {} for t in TYPES} for lvl in LEVELS}
    for lvl in LEVELS:
        for t in TYPES:
            sub = [r for r in rows if r["level"] == lvl and r["group"] == t and r["n_fit"] >= 4]
            d = {"n_streams": len(sub)}
            for m in METHODS:
                vp = [r[f"valid_{m}"] for r in sub if np.isfinite(r[f"valid_{m}"])]
                sd = [r[f"sd_{m}"] for r in sub if np.isfinite(r[f"sd_{m}"])]
                d[f"valid_{m}"] = float(np.median(vp)) if vp else np.nan
                d[f"sd_{m}"] = float(np.median(sd)) if sd else np.nan
            # paired deltas (same stream, both methods)
            dv = [r["valid_v2"] - r["valid_v1.1"] for r in sub
                  if np.isfinite(r["valid_v2"]) and np.isfinite(r["valid_v1.1"])]
            dsd = [r["sd_v2"] - r["sd_v1.1"] for r in sub
                   if np.isfinite(r["sd_v2"]) and np.isfinite(r["sd_v1.1"])]
            d["d_valid_median"] = float(np.median(dv)) if dv else np.nan
            d["d_valid_winfrac"] = float(np.mean([x > 0 for x in dv])) if dv else np.nan
            d["d_sd_median"] = float(np.median(dsd)) if dsd else np.nan
            d["d_sd_better_frac"] = float(np.mean([x < 0 for x in dsd])) if dsd else np.nan
            agg[lvl][t] = d
    return agg


def fig_bars(agg, path):
    fig, axes = plt.subplots(2, 2, figsize=(20, 11))
    x = np.arange(len(TYPES)); w = 0.36
    for col, lvl in enumerate(LEVELS):
        a_valid = axes[0][col]; a_sd = axes[1][col]
        for k, m in enumerate(METHODS):
            v = [agg[lvl][t][f"valid_{m}"] for t in TYPES]
            s = [agg[lvl][t][f"sd_{m}"] for t in TYPES]
            a_valid.bar(x + (k - 0.5) * w, v, w, label=m, color=COLOR[m])
            a_sd.bar(x + (k - 0.5) * w, s, w, label=m, color=COLOR[m])
        for ax, ttl in ((a_valid, f"{lvl} — valid% on clear nights"), (a_sd, f"{lvl} — σ_SD (%)")):
            ax.set_xticks(x); ax.set_xticklabels([f"{t}\n(n={agg[lvl][t]['n_streams']})" for t in TYPES])
            ax.set_title(ttl, fontweight="bold"); ax.grid(axis="y", alpha=0.3); ax.legend()
        a_valid.set_ylabel("valid calibrations (%)  median over streams")
        a_sd.set_ylabel("σ_SD (% of median C)  median over streams")
    fig.suptitle("Network-wide: optimized E-PROF v2 (C8) vs E-PROF v1.1 — 164 streams, 2026",
                 fontsize=15)
    fig.tight_layout()
    fig.savefig(path, dpi=130); plt.close(fig)


def fig_paired(rows, path):
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    for col, lvl in enumerate(LEVELS):
        for row_i, (key, lab, lim) in enumerate(
                [("valid", "valid% on clear nights", (0, 100)), ("sd", "σ_SD (%)", (0, 25))]):
            ax = axes[row_i][col]
            for t in TYPES:
                sub = [r for r in rows if r["level"] == lvl and r["group"] == t]
                xs = [r[f"{key}_v1.1"] for r in sub]
                ys = [r[f"{key}_v2"] for r in sub]
                pts = [(a, b) for a, b in zip(xs, ys) if np.isfinite(a) and np.isfinite(b)]
                if pts:
                    ax.scatter(*zip(*pts), s=30, alpha=0.6, color=TCOLOR[t], label=t, edgecolor="k", linewidth=0.3)
            ax.plot(lim, lim, "k--", lw=1, alpha=0.6)
            ax.set_xlim(lim); ax.set_ylim(lim)
            ax.set_xlabel(f"E-PROF v1.1 — {lab}"); ax.set_ylabel(f"optimized v2 (C8) — {lab}")
            better = "above" if key == "valid" else "below"
            ax.set_title(f"{lvl} — {lab}  (v2 {better} dashed = better)", fontweight="bold")
            ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Paired per-stream comparison: optimized v2 vs v1.1 (each point = one instrument)",
                 fontsize=15)
    fig.tight_layout()
    fig.savefig(path, dpi=130); plt.close(fig)


def main():
    data = load()
    rows = per_stream(data)
    agg = aggregate(rows)
    (DIR / "network_summary.json").write_text(
        json.dumps({"aggregate": agg, "per_stream": rows}, indent=1), encoding="utf-8")
    fig_bars(agg, DIR / "fig_net_bars.png")
    fig_paired(rows, DIR / "fig_net_paired.png")

    print(f"loaded {len(data)} (level,stream) files; {len(rows)} usable streams")
    for lvl in LEVELS:
        print(f"\n=== {lvl}  (median over streams) ===")
        print(f"{'type':10s} {'n':>4s} | {'valid v2':>9s} {'valid v1.1':>10s} {'Δvalid':>7s} {'win%':>5s} | "
              f"{'σSD v2':>7s} {'σSD v1.1':>9s} {'ΔσSD':>6s}")
        for t in TYPES:
            a = agg[lvl][t]
            print(f"{t:10s} {a['n_streams']:4d} | {a['valid_v2']:9.1f} {a['valid_v1.1']:10.1f} "
                  f"{a['d_valid_median']:+7.1f} {100*a['d_valid_winfrac']:4.0f}% | "
                  f"{a['sd_v2']:7.1f} {a['sd_v1.1']:9.1f} {a['d_sd_median']:+6.1f}".replace("nan", " - "))
    print("\nsaved network_summary.json + fig_net_bars.png + fig_net_paired.png")


if __name__ == "__main__":
    main()
