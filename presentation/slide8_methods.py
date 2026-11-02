"""Slide 8 of the E-PROFILE Annual Meeting 2026 deck: cloud versus Rayleigh calibration on every
instrument type (paper Fig. 11, validation/paper/fig_method_comparison.py) redrawn for a slide: same
data and statistics (reuses load / paired_ratio from the paper script), only panels (a) success rate
and (b) C_Rayleigh / C_cloud, with larger fonts, short titles and the ratio axis zoomed to 0.2-10.

Run from the repo root:  python presentation/slide8_methods.py
Output: presentation/figs_annual2026/fig_methods_slide.png (+ per-type summary printed as JSON)
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation.paper import fig_method_comparison as M  # noqa: E402

OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_methods_slide.png"

# short station names so the rotated labels fit on the slide
SHORT = {"Athens Observato": "Athens", "Keflavikurflugvo": "Keflavik", "Wien Hohe Warte": "Wien",
         "Amsterdam Ap Sch": "Amsterdam", "Mini Mpl Mobile": "Mobile", "Saint-Christophe": "St-Christophe",
         "Sion Epfl Alpole": "Sion"}


def rows_():
    rows = []
    for s in M.STREAMS:
        c = M.load(s["key"])
        if s["type"] in M.PHOTON_COUNTING:
            rr = 1.0 / M.paired_ratio(c["cloud"], c["rayleigh"])
        else:
            rr = M.paired_ratio(c["rayleigh"], c["cloud"])
        rows.append(dict(key=s["key"], type=s["type"], site=SHORT.get(s["site"].replace("_", " ").title()[:16], s["site"].replace("_", " ").title()[:13]),
                         p_ray=100 * len(c["rayleigh"]) / M.N_DAYS, p_cld=100 * len(c["cloud"]) / M.N_DAYS,
                         ratio=rr))
    rows.sort(key=lambda r: (M.ORDER.index(r["type"]), r["site"]))
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 19, "axes.titlesize": 20, "axes.labelsize": 20,
                         "xtick.labelsize": 17, "ytick.labelsize": 18, "legend.fontsize": 18})
    rows = rows_()
    x = np.arange(len(rows))
    fig, axs = plt.subplots(2, 1, figsize=(17.28, 6.4), sharex=True,
                            gridspec_kw=dict(height_ratios=[1, 1.35], hspace=0.12))
    fig.subplots_adjust(left=0.07, right=0.995, top=0.925, bottom=0.235)
    ax = axs[0]
    ax.bar(x - 0.2, [r["p_ray"] for r in rows], 0.4, color="#4c72b0", label="Rayleigh")
    ax.bar(x + 0.2, [r["p_cld"] for r in rows], 0.4, color="#999999", label="Liquid cloud")
    ax.set_ylabel("Days [%]")
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 40, 80])
    ax.legend(ncol=2, loc="upper right", frameon=False, borderaxespad=0.1, handlelength=1.2)
    ax = axs[1]
    summ = {}
    for i, r in enumerate(rows):
        q = r["ratio"]
        col = M.TCOL[r["type"]]
        if q.size == 0:
            ax.text(i, 1.0, "no pair", ha="center", va="center", fontsize=14, color="0.4", rotation=90)
            continue
        p10, p25, p50, p75, p90 = np.percentile(q, [10, 25, 50, 75, 90])
        ax.vlines(i, p10, p90, color=col, lw=1.6)
        ax.add_patch(plt.Rectangle((i - 0.32, p25), 0.64, p75 - p25, facecolor=col, alpha=0.5, edgecolor=col))
        ax.hlines(p50, i - 0.32, i + 0.32, color="k", lw=2.4)
    for t in M.ORDER:
        meds = [float(np.median(r["ratio"])) for r in rows if r["type"] == t and r["ratio"].size]
        summ[t] = dict(n=sum(1 for r in rows if r["type"] == t),
                       ratio_median_of_medians=round(float(np.median(meds)), 2) if meds else None,
                       ratio_range=[round(min(meds), 2), round(max(meds), 2)] if meds else None,
                       p_ray_median=round(float(np.median([r["p_ray"] for r in rows if r["type"] == t])), 1),
                       p_cld_median=round(float(np.median([r["p_cld"] for r in rows if r["type"] == t])), 1))
    ax.axhline(1, color="k", lw=1.0)
    ax.axhspan(0.9, 1.1, color="0.85", zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 30)
    ax.set_yticks([0.2, 0.5, 1, 2, 5, 10, 20])
    ax.set_yticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
    ax.set_ylabel(r"C$_{Rayleigh}$ / C$_{cloud}$")
    ax.text(0.995, 0.03, "box: IQR · whiskers: 10–90 % · grey: ±10 %",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=16, color="0.3")
    ax.set_xticks(x)
    ax.set_xticklabels([r["site"] for r in rows], rotation=45, ha="right", rotation_mode="anchor")
    for ax in axs:
        ax.grid(alpha=0.3, axis="y")
        ax.set_xlim(-0.6, len(rows) - 0.4)
        start = 0
        for t in M.ORDER:
            n = sum(1 for r in rows if r["type"] == t)
            if n:
                ax.axvspan(start - 0.5, start + n - 0.5, color=M.TCOL[t], alpha=0.08, lw=0)
                if ax is axs[0]:
                    ax.text(start + n / 2 - 0.5, 1.03, t, transform=ax.get_xaxis_transform(), ha="center",
                            fontsize=22, fontweight="bold", color=M.TCOL[t])
            start += n
    for lab, r in zip(axs[1].get_xticklabels(), rows):
        lab.set_color(M.TCOL[r["type"]])
    fig.align_ylabels(axs)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print("wrote", OUT_PNG)
    print("STATS", json.dumps(summ))


if __name__ == "__main__":
    main()
