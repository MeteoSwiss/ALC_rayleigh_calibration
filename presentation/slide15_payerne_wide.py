"""Slide 15 of the E-PROFILE Annual Meeting 2026 deck: the paper's Payerne figure (fig_v09_payerne,
validation/paper/run_paper_v09.py::fig_station) re-laid out as ONE wide row of four panels with
larger fonts for a 16:9 slide. Same cache, same statistics, same panel functions -- only the layout,
the titles and the font sizes change.

Run from the repo root:  python presentation/slide15_payerne_wide.py
Output: presentation/figs_annual2026/fig14_payerne_wide.png
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from validation.paper import run_paper_v09 as P  # noqa: E402
from validation.paper import figures as FIG  # noqa: E402

OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig14_payerne_wide.png"


def main(key="payerne"):
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15,
                         "xtick.labelsize": 14, "ytick.labelsize": 14})
    st = P.STATIONS[key]
    c = np.load(P.OUT / f"cache_{key}.npz")
    z, beta, ref = c["z"], c["beta"].astype("f8"), int(c["ref"])
    labels, chans = list(c["labels"]), st["channels"]
    cols = FIG.channel_colors(chans)
    comps = [k for k in range(len(chans)) if k != ref]
    minn = max(30, int(0.02 * beta.shape[1]))
    band = P.ZBAND_STATION
    common = np.all(np.isfinite(beta), axis=0)
    n_common_h = int(np.any(common[:, (z >= band[0]) & (z <= band[1])], axis=1).sum())
    prof = [(np.where(common, beta[k], np.nan), cols[k], labels[k] + (" (ref.)" if k == ref else ""))
            for k in range(len(chans))]
    layers = [(100, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 5000)]
    stats, rd = {}, []
    for k in comps:
        s = P.pair_stats(beta[k], beta[ref], z, band, layers)
        stats[labels[k]] = s
        rd.append((P.rel_diff(beta[k], beta[ref]), cols[k], labels[k], s["med"]))
    zmax = 4000
    fig, axs = plt.subplots(1, 4, figsize=(27, 8.2),
                            gridspec_kw=dict(width_ratios=[1, 1, 1.12, 1.12], wspace=0.32),
                            constrained_layout=False)
    fig.subplots_adjust(left=0.04, right=0.985, top=0.86, bottom=0.13)
    P.panel_profiles(axs[0], z, prof, band, zmax, minn,
                     f"(a) Median profiles (± IQR)\n{n_common_h} common hours")
    P.panel_reldiff(axs[1], z, rd, band, zmax, minn, f"(b) Difference to {labels[ref]}")
    zb = (z >= band[0]) & (z <= band[1])
    for ax, k, lt in zip(axs[2:], comps, "cd"):
        s = stats[labels[k]]
        P.panel_binscatter(ax, beta[ref][:, zb].ravel(), beta[k][:, zb].ravel(), cols[k],
                           f"{labels[ref]}  " + P.BLAB, f"{labels[k]}  " + P.BLAB,
                           f"({lt}) {labels[k]} vs {labels[ref]}, {int(band[0])}–{int(band[1])} m\n"
                           f"median {s['med']:+.1f} %, log r = {s['r_log']:.2f}")
    for ax in axs:
        leg = ax.get_legend()
        if leg is not None:
            plt.setp(leg.get_texts(), fontsize=13)
    axs[1].set_ylabel("")
    fig.suptitle(f"{st['name']}: {P._period(c['time'])}, operational calibration",
                 fontsize=18, fontweight="bold", y=0.985)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)
    print("wrote", OUT_PNG, {k: round(v["med"], 1) for k, v in stats.items()}, n_common_h)


if __name__ == "__main__":
    main()
