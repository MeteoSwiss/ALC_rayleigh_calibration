"""Plot C_L time series (baseline 300s/K0 vs P5 native/K7+fallback) per station, so the night-to-night
scatter behind sigma_SD can be eyeballed. One landscape figure per type. Uses cloud_yield_experiment
output. Run: python validation/cloud_yield_plot_timeseries.py
"""
import glob
import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CYD = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_yield")
REPORT = Path(__file__).resolve().parents[1] / "doc" / "reports" / "cloud_yield_l1"
META = {m["label"]: m for m in json.load(open(Path(__file__).resolve().parents[1] / "validation" / "scope_cloud_2026.json"))}
CONST = {"CL31": 1e8, "CL51": 1e8, "CL61": 1.0}          # P5 uses the CL61 fallback
CONST_BASE = {"CL31": 1e8, "CL51": 1e8, "CL61": np.nan}  # baseline: CL61 has no constant


def series(days, variant, config, const):
    xs, ys = [], []
    for ds in sorted(days):
        cell = days[ds].get(variant, {}).get(config)
        if cell and cell[1] >= 1 and np.isfinite(cell[0]) and cell[0] > 0 and np.isfinite(const) and const > 0:
            xs.append(datetime.strptime(ds, "%Y%m%d"))
            ys.append(const / cell[0])
    return xs, np.array(ys)


def sigma_sd(y):
    y = np.asarray(y, float)
    if y.size < 4:
        return np.nan
    m = np.median(y)
    return 1.4826 * np.median(np.abs(np.diff(y))) / np.sqrt(2) / abs(m) * 100 if m > 0 else np.nan


def plot_type(t, labels, ncol=4):
    n = len(labels)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.0 * nrow), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for i, lab in enumerate(labels):
        ax = axes[i // ncol][i % ncol]
        ax.axis("on")
        days = json.loads((CYD / f"cy_{lab}.json").read_text())
        xb, yb = series(days, "base_300s", "K0", CONST_BASE[t])
        xp, yp = series(days, "native", "K7", CONST[t])
        if len(yb):
            ax.plot(xb, yb, "o-", ms=4, lw=0.8, color="#7f7f7f", label=f"baseline (n={len(yb)}, σ{sigma_sd(yb):.0f})")
            ax.axhline(np.median(yb), color="#7f7f7f", ls=":", lw=0.8)
        if len(yp):
            ax.plot(xp, yp, "s-", ms=4, lw=0.8, color="#d62728", alpha=0.8,
                    label=f"P5 combo (n={len(yp)}, σ{sigma_sd(yp):.0f})")
            md = np.median(yp)
            ax.axhline(md, color="#d62728", ls="--", lw=0.8)
            ax.axhspan(md * (1 - sigma_sd(yp) / 100), md * (1 + sigma_sd(yp) / 100), color="#d62728", alpha=0.07)
        ax.set_title(f"{lab}", fontsize=9)
        ax.legend(fontsize=7, loc="best")
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2)
        for lb in ax.get_xticklabels():
            lb.set_rotation(30); lb.set_ha("right")
    fig.suptitle(f"{t} — C_L time series: baseline (300 s / K0) vs P5 (native / K7 + CL61 fallback). "
                 f"Dashed = median, shaded = ±σ_SD.", fontsize=13)
    fig.tight_layout()
    out = REPORT / f"cl_timeseries_{t}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)


def main():
    by = {"CL31": [], "CL51": [], "CL61": []}
    for fp in glob.glob(str(CYD / "cy_*.json")):
        lab = Path(fp).stem[3:]
        t = META.get(lab, {}).get("group")
        if t in by:
            by[t].append(lab)
    plot_type("CL61", sorted(by["CL61"]), ncol=4)                 # all 11
    plot_type("CL31", sorted(by["CL31"])[:8], ncol=4)
    plot_type("CL51", sorted(by["CL51"])[:8], ncol=4)


if __name__ == "__main__":
    main()
