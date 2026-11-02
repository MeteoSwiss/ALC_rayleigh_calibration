"""Slide 5 of the E-PROFILE Annual Meeting 2026 deck: the Cloudnet droplet effective-radius climatology
(paper Fig. A4-1, validation/_cloudnet_droplet_climatology.py) redrawn for a slide: same data, same
samples (all liquid pixels, cloud base, the valid-calibration scenes) and the same 0.75 um bins, but
one shared legend, no suptitle, larger fonts, and the 5.5 um radius of the PVC tables marked.

The expensive part (reading every Cloudnet `der` file on D:) is done once and cached as fine
histograms in presentation/figs_annual2026/droplet_hist_cache.npz; rerun with --refresh to rebuild.

Run from the repo root:  python presentation/slide5_droplets.py [--refresh]
Output: presentation/figs_annual2026/fig_droplets_slide.png (+ medians printed as JSON)
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation import _cloudnet_droplet_climatology as D  # noqa: E402

FIGDIR = REPO / "presentation" / "figs_annual2026"
CACHE = FIGDIR / "droplet_hist_cache.npz"
OUT_PNG = FIGDIR / "fig_droplets_slide.png"
FINE = np.arange(0.0, 21.75 + 1e-9, 0.25)      # 3 fine bins = one 0.75 um bin of the paper figure
SAMPLES = ("all", "base", "cal")


def build_cache():
    out = {}
    for s in D.SITES:
        per = dict(all=[], base=[], cal=D._calib_scene_re(s))
        files = sorted(D.CLOUDNET.glob(f"*_{s}_der.nc"))
        print(f"{s}: {len(files)} daily files", flush=True)
        with ProcessPoolExecutor(max_workers=D.WORKERS) as ex:
            for res in ex.map(D.read_one, files, chunksize=4):
                if res is None:
                    continue
                a, b = res
                if a.size:
                    per["all"].append(a)
                if b.size:
                    per["base"].append(b)
        for lab in ("all", "base"):
            per[lab] = np.concatenate(per[lab]) if per[lab] else np.array([])
        for lab in SAMPLES:
            x = per[lab][np.isfinite(per[lab])]
            xin = x[x < 22]                          # the paper figure's display sample
            out[f"{s}_{lab}_counts"] = np.histogram(xin, bins=FINE)[0]
            out[f"{s}_{lab}_med22"] = np.median(xin) if xin.size else np.nan
            out[f"{s}_{lab}_med"] = np.median(x) if x.size else np.nan
            out[f"{s}_{lab}_n"] = xin.size
    np.savez_compressed(CACHE, fine_edges=FINE, **out)
    return dict(np.load(CACHE))


def main():
    c = build_cache() if ("--refresh" in sys.argv or not CACHE.exists()) else dict(np.load(CACHE))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 18, "axes.labelsize": 16,
                         "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14})
    edges = c["fine_edges"][::3]                     # 0.75 um bins, 0 .. 21.75
    width = np.diff(edges)
    style = {"all": dict(fill=True, color="#6BAED6", alpha=0.55, label="All liquid-cloud pixels"),
             "base": dict(fill=True, color="#9E9E9E", alpha=0.55, label="Cloud base"),
             "cal": dict(fill=False, color="#D8232A", lw=3.2, label="Calibration scenes")}
    disp = {"payerne": "Payerne", "palaiseau": "Palaiseau", "lindenberg": "Lindenberg"}
    fig, axes = plt.subplots(1, 3, figsize=(11.28, 6.0), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.80, bottom=0.14, wspace=0.08)
    stats = {}
    for ax, s in zip(axes, D.SITES):
        ax.axvspan(4.0, 10.0, color="green", alpha=0.10, lw=0, label="Le et al. 2026: 4–10 µm")
        for lab in SAMPLES:
            k = np.asarray(c[f"{s}_{lab}_counts"], float)
            k = k[: 3 * (k.size // 3)].reshape(-1, 3).sum(axis=1)
            if k.sum() == 0:
                continue
            ax.stairs(k / (k.sum() * width), edges, **style[lab])
        ax.axvline(5.5, color="k", ls="--", lw=2.0, label="PVC tables: 5.5 µm")
        med = float(c[f"{s}_cal_med22"])
        ax.text(0.97, 0.95, f"calibration scenes\nmedian {med:.1f} µm", transform=ax.transAxes,
                ha="right", va="top", fontsize=15, color="#D8232A")
        ax.set_title(disp[s])
        ax.set_xlim(0, 20)
        ax.set_xticks([0, 5, 10, 15, 20])
        ax.grid(alpha=0.3)
        stats[s] = {lab: dict(median=round(float(c[f"{s}_{lab}_med22"]), 2),
                              median_all_values=round(float(c[f"{s}_{lab}_med"]), 2),
                              n=int(c[f"{s}_{lab}_n"])) for lab in SAMPLES}
    axes[0].set_ylabel("Probability density")
    fig.supxlabel("Droplet effective radius [µm]", fontsize=16, y=0.015)
    h, l = axes[0].get_legend_handles_labels()
    order = [l.index(x) for x in ("All liquid-cloud pixels", "Cloud base", "Calibration scenes",
                                   "Le et al. 2026: 4–10 µm", "PVC tables: 5.5 µm") if x in l]
    fig.legend([h[i] for i in order], [l[i] for i in order], loc="upper center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, 1.0), columnspacing=1.6)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=200)
    print("wrote", OUT_PNG)
    print("STATS", json.dumps(stats))


if __name__ == "__main__":
    main()
