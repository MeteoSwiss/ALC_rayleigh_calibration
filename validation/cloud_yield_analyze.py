"""cloud_yield_analyze.py -- summarize cloud_yield_experiment.py into per-(type, proposition) results.

A calibration is DASHBOARD-VALID when it yields a positive lidar constant C_L = const_type / coef
(const: CL31/CL51 = 1e8; CL61 = NaN at baseline, or 1.0 with the fallback). Each proposition selects
an averaging variant + a gate config + a CL61 C_L rule:

  baseline   base_300s  K0   CL61 const = NaN   (the current dashboard)
  P1_avg     native     K0   CL61 const = NaN   (adapt averaging)
  P2_consec  base_300s  K1   CL61 const = NaN   (n_consecutive 5->3)
  P3_gates   base_300s  K7   CL61 const = NaN   (ease cloud gates)
  P4_cl61    base_300s  K0   CL61 const = 1.0   (CL61 applied-constant fallback)
  P5_combo   native     K7   CL61 const = 1.0   (combined)

Per (type, proposition): valid% = valid days / days-with-data (median over streams), and sigma_SD =
robust successive-difference precision of C_L (% of median). Writes a JSON summary + a bar figure.
Run:  python validation/cloud_yield_analyze.py
"""
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
CYD = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloud_yield")
META = {m["label"]: m for m in json.loads((REPO / "validation" / "scope_cloud_2026.json").read_text())}
TYPES = ["CL31", "CL51", "CL61"]
CONST = {"CL31": 1e8, "CL51": 1e8, "CL61": np.nan}            # baseline C_L constant per type
CONST_FB = {"CL31": 1e8, "CL51": 1e8, "CL61": 1.0}            # with the CL61 fallback

# proposition -> (variant, config, const-table)
PROPS = {
    "baseline":  ("base_300s", "K0", CONST),
    "P1_avg":    ("native",    "K0", CONST),
    "P2_consec": ("base_300s", "K1", CONST),
    "P3_gates":  ("base_300s", "K7", CONST),
    "P4_cl61":   ("base_300s", "K0", CONST_FB),
    "P5_combo":  ("native",    "K7", CONST_FB),
}


def sigma_sd(cl):
    """Robust successive-difference sigma of C_L over valid days, % of median (convention-free)."""
    c = np.asarray([x for x in cl if np.isfinite(x) and x > 0], float)
    if c.size < 4:
        return np.nan
    m = np.median(c)
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) / abs(m) * 100) if m > 0 else np.nan


def cl_series(days, variant, config, const):
    """C_L per day for a (variant, config) under a type constant; NaN where not valid."""
    out = []
    for ds in sorted(days):
        cell = days[ds].get(variant, {}).get(config)
        if not cell:
            out.append(np.nan); continue
        coef, n = cell
        cl = (const / coef) if (n >= 1 and np.isfinite(coef) and coef > 0 and np.isfinite(const) and const > 0) else np.nan
        out.append(cl)
    return out


def main():
    # per type: list over streams of {prop: (valid%, sigma_SD)}
    per = {t: [] for t in TYPES}
    nstream = {t: 0 for t in TYPES}
    for fp in glob.glob(str(CYD / "cy_*.json")):
        lab = Path(fp).stem[3:]
        t = META.get(lab, {}).get("group")
        if t not in TYPES:
            continue
        days = json.loads(Path(fp).read_text())
        if not days:
            continue
        nstream[t] += 1
        ndata = len(days)
        rec = {}
        for pname, (variant, config, ctab) in PROPS.items():
            cls = cl_series(days, variant, config, ctab[t])
            nvalid = sum(1 for x in cls if np.isfinite(x) and x > 0)
            rec[pname] = (100.0 * nvalid / ndata if ndata else np.nan, sigma_sd(cls))
        per[t].append(rec)

    def med(xs):
        xs = [x for x in xs if np.isfinite(x)]
        return float(np.median(xs)) if xs else float("nan")

    summary = {}
    for t in TYPES:
        summary[t] = {"n_streams": nstream[t], "props": {}}
        for pname in PROPS:
            summary[t]["props"][pname] = {
                "valid_pct": med([r[pname][0] for r in per[t]]),
                "sigma_sd": med([r[pname][1] for r in per[t]]),
            }
    (CYD / "cloud_yield_summary.json").write_text(json.dumps(summary, indent=2))

    # print table
    print(f"{'type':6s} {'n':>3s}  " + "  ".join(f"{p:>10s}" for p in PROPS))
    for t in TYPES:
        vp = "  ".join(f"{summary[t]['props'][p]['valid_pct']:6.0f}%   " for p in PROPS)
        print(f"{t:6s} {nstream[t]:3d}  valid%: {vp}")
        ss = "  ".join(f"{summary[t]['props'][p]['sigma_sd']:6.1f}    " for p in PROPS)
        print(f"{'':6s} {'':3s}  sigma : {ss}")

    # figure: grouped bars (valid% per proposition, one panel per type) + sigma annotation
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    x = np.arange(len(PROPS))
    colors = ["#7f7f7f", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#d62728"]
    for ax, t in zip(axes, TYPES):
        vals = [summary[t]["props"][p]["valid_pct"] for p in PROPS]
        sig = [summary[t]["props"][p]["sigma_sd"] for p in PROPS]
        ax.bar(x, vals, color=colors)
        for xi, v, s in zip(x, vals, sig):
            ax.text(xi, (v if np.isfinite(v) else 0) + 1,
                    f"{v:.0f}%\nσ{s:.0f}" if np.isfinite(v) else "0%", ha="center", fontsize=9)
        ax.set_xticks(x); ax.set_xticklabels(list(PROPS), rotation=35, ha="right", fontsize=9)
        ax.set_title(f"{t}  (n={nstream[t]} streams)")
        ax.set_ylabel("dashboard-valid % (median over streams)")
        ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Liquid-cloud calibration yield (DASHBOARD-VALID C_L) — 5 propositions vs baseline (L1)",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(CYD / "cloud_yield_figure.png", dpi=150)
    print("wrote", CYD / "cloud_yield_summary.json", "and cloud_yield_figure.png")


if __name__ == "__main__":
    main()
