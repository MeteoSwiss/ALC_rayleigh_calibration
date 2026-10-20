# -*- coding: utf-8 -*-
"""Paper Section 6.1: cloud versus Rayleigh calibration on every instrument type of the network.

Both methods were FORCED on every selected instrument (run_method_test.ps1: operational v2.2
settings, Feb-Jun 2026), although the operational routing excludes the Rayleigh method for the
CL31/CL51 and the cloud method for the photon-counting CHM15k/Mini-MPL. Per instrument:
  * success rate of each method (% of days with a valid calibration, flag 1 or 0.5);
  * ratio C_L,Rayleigh / C_L,cloud: every nightly constant of one method divided by the median of
    the other method's constants within +/-15 days (Rayleigh night / cloud for the Vaisala types,
    cloud day / Rayleigh for the photon-counting types, so that the NON-operational method is the
    one tested); shown as C_Rayleigh / C_cloud in both cases;
  * day-to-day scatter of each method (1.4826 x MAD of log C, %).

Run:  python validation/paper/fig_method_comparison.py
"""
from __future__ import annotations
import csv
import json
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RUN = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/method_test_2026")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08")
STREAMS = json.loads((REPO / "validation/paper/method_test_streams.json").read_text())
N_DAYS = (date(2026, 6, 29) - date(2026, 2, 1)).days + 1
ORDER = ["CL31", "CL51", "CL61", "CHM15k", "Mini-MPL"]
TCOL = {"CL31": "#1b9e77", "CL51": "#66a61e", "CL61": "#7570b3", "CHM15k": "#d95f02", "Mini-MPL": "#e7298a"}
PHOTON_COUNTING = {"CHM15k", "Mini-MPL"}
PAIR_DAYS = 15


def load(key):
    f = RUN / key / f"{key}_cal.csv"
    out = {"rayleigh": {}, "cloud": {}}
    files = [f] if f.exists() else sorted(RUN.glob(f"_parts/{key}__*/{key}/{key}_cal.csv"))   # unmerged run
    for r in (row for p in files for row in csv.DictReader(open(p, encoding="utf-8"))):
        try:
            flag, val = float(r["flag"]), float(r["cal_value"])
        except ValueError:
            continue
        if flag in (1.0, 0.5) and val > 0 and r["method"] in out:
            d = date(int(r["date"][:4]), int(r["date"][4:6]), int(r["date"][6:8]))
            out[r["method"]][d] = val
    return out


def paired_ratio(test, ref):
    """Each test value / median of the reference values within +/-PAIR_DAYS (else whole period)."""
    if not test or not ref:
        return np.array([])
    rd = np.array([d.toordinal() for d in ref])
    rv = np.array(list(ref.values()))
    out = []
    for d, v in test.items():
        m = np.abs(rd - d.toordinal()) <= PAIR_DAYS
        out.append(v / np.median(rv[m] if m.sum() >= 3 else rv))
    return np.array(out)


def scatter(vals):
    if len(vals) < 5:
        return np.nan
    y = np.log(np.array(list(vals)))
    return 100 * 1.4826 * np.median(np.abs(y - np.median(y)))


def main():
    rows = []
    for s in STREAMS:
        c = load(s["key"])
        if s["type"] in PHOTON_COUNTING:
            rr = 1.0 / paired_ratio(c["cloud"], c["rayleigh"])       # C_Rayleigh / C_cloud(day)
        else:
            rr = paired_ratio(c["rayleigh"], c["cloud"])
        rows.append(dict(key=s["key"], type=s["type"], site=s["site"].replace("_", " ").title()[:16],
                         n_ray=len(c["rayleigh"]), n_cld=len(c["cloud"]),
                         p_ray=100 * len(c["rayleigh"]) / N_DAYS, p_cld=100 * len(c["cloud"]) / N_DAYS,
                         sc_ray=scatter(c["rayleigh"].values()), sc_cld=scatter(c["cloud"].values()),
                         ratio=rr))
    rows.sort(key=lambda r: (ORDER.index(r["type"]), r["site"]))
    summary = {}
    for r in rows:
        q = r["ratio"]
        r["r_med"] = float(np.median(q)) if q.size else np.nan
        r["r_q1"], r["r_q3"] = (float(np.percentile(q, 25)), float(np.percentile(q, 75))) if q.size else (np.nan, np.nan)
        summary[r["key"]] = {k: (float(v) if isinstance(v, (float, np.floating)) else v)
                             for k, v in r.items() if k != "ratio"}
        print(f"{r['type']:8s} {r['site']:16s} Rayleigh {r['n_ray']:3d} ({r['p_ray']:4.0f} %, scatter "
              f"{r['sc_ray']:5.1f} %) | cloud {r['n_cld']:3d} ({r['p_cld']:4.0f} %, scatter {r['sc_cld']:5.1f} %) | "
              f"C_Ray/C_cloud {r['r_med']:.2f} [{r['r_q1']:.2f}-{r['r_q3']:.2f}] n={q.size}")
    for t in ORDER:
        sub = [r for r in rows if r["type"] == t]
        meds = [r["r_med"] for r in sub if np.isfinite(r["r_med"])]
        summary[f"type_{t}"] = dict(
            n_instruments=len(sub), p_ray_median=float(np.median([r["p_ray"] for r in sub])),
            p_cld_median=float(np.median([r["p_cld"] for r in sub])),
            ratio_median_of_medians=float(np.median(meds)) if meds else None,
            ratio_range=[float(min(meds)), float(max(meds))] if meds else None,
            sc_ray_median=float(np.nanmedian([r["sc_ray"] for r in sub])),
            sc_cld_median=float(np.nanmedian([r["sc_cld"] for r in sub])))
        print(t, summary[f"type_{t}"])
    (OUT / "method_comparison_summary.json").write_text(json.dumps(summary, indent=1))
    figure(rows)


def figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(len(rows))
    fig, axs = plt.subplots(3, 1, figsize=(18, 11), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.5, 1]))
    # (a) success rates
    ax = axs[0]
    ax.bar(x - 0.2, [r["p_ray"] for r in rows], 0.4, color="#4c72b0", label="Rayleigh")
    ax.bar(x + 0.2, [r["p_cld"] for r in rows], 0.4, color="#999999", label="liquid cloud")
    ax.set_ylabel("days with a valid\ncalibration [%]")
    ax.set_title("(a) Success rate of each method, 1 February - 29 June 2026", fontsize=11, loc="left", pad=22)
    ax.legend(fontsize=9, ncol=2, loc="upper left")
    # (b) ratio
    ax = axs[1]
    for i, r in enumerate(rows):
        q = r["ratio"]
        col = TCOL[r["type"]]
        if q.size == 0:
            ax.text(i, 1.0, "no valid pair", ha="center", va="center", fontsize=7.5, color="0.4", rotation=90)
            continue
        p10, p25, p50, p75, p90 = np.percentile(q, [10, 25, 50, 75, 90])
        ax.vlines(i, p10, p90, color=col, lw=1.2)
        ax.add_patch(plt.Rectangle((i - 0.3, p25), 0.6, p75 - p25, facecolor=col, alpha=0.45, edgecolor=col))
        ax.hlines(p50, i - 0.3, i + 0.3, color="k", lw=2)
        ax.text(i, 0.13, f"{q.size}", ha="center", fontsize=7.5, color="0.3")
    ax.axhline(1, color="k", lw=0.9)
    ax.axhspan(0.9, 1.1, color="0.9", zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 10)
    ax.set_ylabel(r"C$_{L,Rayleigh}$ / C$_{L,cloud}$")
    ax.set_title("(b) Ratio of the two calibrations (box: IQR, line: median, whiskers: 10-90 %; "
                 "number of calibrations at the bottom; grey band: +/-10 %)", fontsize=11, loc="left")
    # (c) day-to-day scatter
    ax = axs[2]
    ax.plot(x - 0.1, [r["sc_ray"] for r in rows], "o", color="#4c72b0", label="Rayleigh")
    ax.plot(x + 0.1, [r["sc_cld"] for r in rows], "s", color="#999999", label="liquid cloud")
    ax.set_ylabel("day-to-day scatter\n(1.48 MAD of ln C) [%]")
    ax.set_title("(c) Stability of each method", fontsize=11, loc="left")
    ax.legend(fontsize=9, ncol=2, loc="upper left")
    ax.set_xticks(x)
    ax.set_xticklabels([r["site"] for r in rows], rotation=60, ha="right", fontsize=8.5)
    for ax in axs:
        ax.grid(alpha=0.3, axis="y")
        start = 0
        for t in ORDER:
            n = sum(1 for r in rows if r["type"] == t)
            if n:
                ax.axvspan(start - 0.5, start + n - 0.5, color=TCOL[t], alpha=0.07, lw=0)
                if ax is axs[0]:
                    ax.text(start + n / 2 - 0.5, 1.03, t, transform=ax.get_xaxis_transform(), ha="center",
                            fontsize=11, fontweight="bold", color=TCOL[t])
            start += n
    for lab, r in zip(axs[2].get_xticklabels(), rows):
        lab.set_color(TCOL[r["type"]])
    fig.tight_layout()
    fig.savefig(OUT / "fig_v08_methods.png", dpi=200)
    print(f"-> {OUT / 'fig_v08_methods.png'}")


if __name__ == "__main__":
    main()
