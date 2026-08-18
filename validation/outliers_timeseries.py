"""
outliers_timeseries.py — calibration-constant time series + per-calibration outlier rate, from the
network run (run_network_v2_vs_v11.py). For each stream/level/method we take the date-ordered valid
nightly lidar constants and flag OUTLIERS in a drift-aware, robust way:

    residual_i = C_i - rolling_median(C, 9 nights)
    sigma_rob  = 1.4826 * median(|residual - median(residual)|)
    outlier    = |residual - median(residual)| > 3 * sigma_rob

(detrending removes slow seasonal drift so only genuine jumps/spikes count). outlier% = flagged/valid.

Outputs (figs_paper_validation/network_v2_v11/):
    outlier_summary.json                      per-stream outlier% (v2 & v1.1, L1 & L2) + type medians
    fig_ts_CL61.png / _MiniMPL.png / _CHM15k.png   nightly-C time series (L2, v2) with outliers in red
    fig_outlier_overview.png                  outlier% distribution by type, v2 vs v1.1, both levels
"""
from __future__ import annotations
import glob
import json
import sys
from datetime import datetime
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
DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/network_v2_v11")
MANIFEST = json.loads((REPO / "validation" / "scope_network_2026.json").read_text())
META = {m["label"]: m for m in MANIFEST}
TYPES = ["CHM15k", "Mini-MPL", "CL61"]
TCOLOR = {"CHM15k": "#1f77b4", "Mini-MPL": "#2ca02c", "CL61": "#d62728"}
METHODS = ["v2", "v1.1"]
I_OK, I_CL = 0, 1
K = 3.0
WIN = 9


def rolling_median(x, win):
    n = len(x); h = win // 2
    return np.array([np.median(x[max(0, i - h):min(n, i + h + 1)]) for i in range(n)])


def flag_outliers(cls):
    cls = np.asarray(cls, float)
    n = cls.size
    if n < 5:
        return np.zeros(n, bool), np.nan
    resid = cls - rolling_median(cls, WIN)
    c = resid - np.median(resid)
    s = 1.4826 * np.median(np.abs(c))
    if not (s > 0):
        return np.zeros(n, bool), 0.0
    out = np.abs(c) > K * s
    return out, 100.0 * out.sum() / n


def series(nights, method):
    order = sorted(nights)
    ds, cls = [], []
    for d in order:
        rec = nights[d][method]
        if rec[I_OK] and np.isfinite(rec[I_CL]) and rec[I_CL] > 0:
            ds.append(datetime.strptime(d, "%Y%m%d")); cls.append(rec[I_CL])
    return ds, np.asarray(cls, float)


def load():
    data = {}
    for fp in glob.glob(str(DIR / "net_*.json")):
        _, level, label = Path(fp).stem.split("_", 2)
        try:
            data[(level, label)] = json.loads(Path(fp).read_text())
        except Exception:
            pass
    return data


def panel(ax, ds, cls, out, title):
    if cls.size:
        med = np.median(cls)
        resid = cls - rolling_median(cls, WIN)
        s = 1.4826 * np.median(np.abs(resid - np.median(resid)))
        ax.fill_between(ds, med - K * s, med + K * s, color="0.85", zorder=0, label="±3σ (robust)")
        ax.axhline(med, color="green", lw=1.0, zorder=1)
        ax.plot(ds, cls, "o", ms=3, color="#1f77b4", alpha=0.7, zorder=2)
        if out.any():
            ax.plot(np.array(ds)[out], cls[out], "o", ms=5, color="red", zorder=3, label="outlier")
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30); lab.set_ha("right")


def ts_grid(data, group, path, level="L2", ncol=4, pick=None):
    labels = [lab for (lvl, lab) in data if lvl == level and META.get(lab, {}).get("group") == group]
    recs = []
    for lab in labels:
        ds, cls = series(data[(level, lab)], "v2")
        out, pct = flag_outliers(cls)
        recs.append((lab, ds, cls, out, pct))
    recs = [r for r in recs if r[2].size >= 5]
    if pick is not None and len(recs) > pick:
        recs.sort(key=lambda r: -r[4])                  # show the worst offenders first
        recs = recs[:pick]
    recs.sort(key=lambda r: META[r[0]].get("site", r[0]))
    n = len(recs)
    if n == 0:
        return
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.8 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for k, (lab, ds, cls, out, pct) in enumerate(recs):
        ax = axes[k // ncol][k % ncol]; ax.axis("on")
        site = META[lab].get("site", lab)[:18]
        panel(ax, ds, cls, out, f"{site} [{lab}]\noutliers {pct:.1f}% (n={cls.size})")
    axes[0][0].legend(fontsize=7, loc="best")
    fig.suptitle(f"{group} — nightly lidar constant ({level}, optimized v2), outliers in red", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(path, dpi=120); plt.close(fig)


def overview(rows, path):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    # (a) box/strip of v2 outlier% by type, per level
    ax = axes[0]
    pos = 0; ticks = []; lbls = []
    for t in TYPES:
        for lvl, off, col in [("L1", -0.18, "#9ecae1"), ("L2", 0.18, "#3182bd")]:
            vals = [r["out_v2"] for r in rows if r["group"] == t and r["level"] == lvl and np.isfinite(r["out_v2"])]
            if vals:
                ax.scatter(np.full(len(vals), pos + off) + np.random.uniform(-0.05, 0.05, len(vals)),
                           vals, s=18, alpha=0.6, color=col)
                ax.scatter([pos + off], [np.median(vals)], marker="_", s=600, color="k", zorder=5)
        ticks.append(pos); lbls.append(t); pos += 1
    ax.set_xticks(ticks); ax.set_xticklabels(lbls)
    ax.set_ylabel("outlier rate (% of valid nights)")
    ax.set_title("v2 outlier rate by type (L1 light / L2 dark; bar = median)", fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    # (b) v2 vs v1.1 outlier% paired (L2)
    ax = axes[1]
    for t in TYPES:
        sub = [r for r in rows if r["group"] == t and r["level"] == "L2"]
        xy = [(r["out_v11"], r["out_v2"]) for r in sub if np.isfinite(r["out_v11"]) and np.isfinite(r["out_v2"])]
        if xy:
            ax.scatter(*zip(*xy), s=28, alpha=0.6, color=TCOLOR[t], label=t, edgecolor="k", linewidth=0.3)
    lim = [0, max(1, max([r["out_v2"] for r in rows if np.isfinite(r["out_v2"])] + [r["out_v11"] for r in rows if np.isfinite(r["out_v11"])]))]
    ax.plot(lim, lim, "k--", lw=1, alpha=0.6); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("v1.1 outlier %"); ax.set_ylabel("v2 outlier %")
    ax.set_title("L2: outlier rate v2 vs v1.1 (below line = v2 cleaner)", fontsize=11)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    # (c) sorted per-stream v2 outlier% (L2)
    ax = axes[2]
    sub = sorted([r for r in rows if r["level"] == "L2" and np.isfinite(r["out_v2"])], key=lambda r: r["out_v2"])
    ax.bar(range(len(sub)), [r["out_v2"] for r in sub],
           color=[TCOLOR[r["group"]] for r in sub])
    ax.set_xlabel("stream (sorted)"); ax.set_ylabel("v2 outlier % (L2)")
    ax.set_title("Per-stream v2 outlier rate, L2 (sorted)", fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=TCOLOR[t], label=t) for t in TYPES], fontsize=8)
    fig.suptitle("Calibration outlier rates — network, optimized v2 vs v1.1 (drift-aware robust, k=3)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=130); plt.close(fig)


def main():
    np.random.seed(0)
    data = load()
    rows = []
    for (level, label), nights in data.items():
        g = META.get(label, {}).get("group")
        if g is None:
            continue
        rec = dict(level=level, label=label, group=g, site=META[label].get("site", label))
        for m in METHODS:
            _, cls = series(nights, m)
            _, pct = flag_outliers(cls)
            rec[f"n_{m}"] = int(cls.size)
            rec[f"out_{'v2' if m=='v2' else 'v11'}"] = pct
        rows.append(rec)
    (DIR / "outlier_summary.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")

    ts_grid(data, "CL61", DIR / "fig_ts_CL61.png", ncol=4)
    ts_grid(data, "Mini-MPL", DIR / "fig_ts_MiniMPL.png", ncol=3)
    ts_grid(data, "CHM15k", DIR / "fig_ts_CHM15k.png", ncol=4, pick=12)
    overview(rows, DIR / "fig_outlier_overview.png")

    print(f"{len(rows)} (level,stream) rows")
    print(f"\n{'type':10s} {'lvl':3s} | {'med out% v2':>11s} {'med out% v1.1':>13s} | {'max v2':>7s}")
    digest = {}
    for t in TYPES:
        for lvl in ["L1", "L2"]:
            sub = [r for r in rows if r["group"] == t and r["level"] == lvl]
            o2 = [r["out_v2"] for r in sub if np.isfinite(r["out_v2"])]
            o1 = [r["out_v11"] for r in sub if np.isfinite(r["out_v11"])]
            if not o2:
                continue
            print(f"{t:10s} {lvl:3s} | {np.median(o2):11.1f} {np.median(o1):13.1f} | {max(o2):7.1f}")
    print("\nsaved outlier_summary.json + 4 figures")


if __name__ == "__main__":
    main()
