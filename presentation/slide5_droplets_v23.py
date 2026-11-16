"""Slide 5 of the Annual Meeting 2026 deck, v2.3 edition ("Multiple scattering per instrument"): the droplet
size behind the generation-2 multiple-scattering tables (a_G = 3.3 um), paper Fig. A4-2 of manuscript v1.0
(validation/paper/fig_v10_droplet_size.py) redrawn with slide-sized fonts. Same data: the Cloudnet equivalent-
area radius a_G = 0.832 r_e in the lowest 100 m above the liquid base at every profile a cloud calibration used
(C:/DATA/Projects/202611_eta_study/step1d/scenes_msfov.csv) and the site medians of the Cloudnet sites of the
E-PROFILE domain (C:/DATA/Projects/202611_eta_study/cloudnet_sites/site_table.csv).

Run from the repo root:  python presentation/slide5_droplets_v23.py
Output: presentation/figs_annual2026/fig_droplets_v23_slide.png (+ a STATS line)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation.paper.fig_v10_droplet_size import A_G, ETA, UNITS  # noqa: E402

OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_droplets_v23_slide.png"


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 15, "axes.labelsize": 15,
                         "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12.5})
    sc = pd.read_csv(ETA / "step1d" / "scenes_msfov.csv", usecols=["key", "date", "n_st1_100", "aG_um_100"])
    sc = sc[(sc.n_st1_100 >= 1) & np.isfinite(sc.aG_um_100)]
    pooled = sc.aG_um_100.values
    st = pd.read_csv(ETA / "cloudnet_sites" / "site_table.csv")
    st = st[(st.eprofile_box.astype(str).str.lower() == "true") & (st.n_days >= 20)].sort_values("aG_med").reset_index(drop=True)
    med_sites = float(np.median(st.aG_med))

    # the slide shows the figure in a 1120 x 596 px box next to the text column: draw it at that aspect ratio
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.2, 5.96), gridspec_kw=dict(width_ratios=[1, 1.15]))
    fig.subplots_adjust(left=0.085, right=0.99, top=0.925, bottom=0.15, wspace=0.62)
    bins = np.arange(0, 10.01, 0.25)
    short = {"Payerne CL61": "Payerne CL61", "Payerne CL31": "Payerne CL31", "Palaiseau CL31": "Palaiseau CL31",
             "Lindenberg CL61": "Lindenberg CL61"}
    for key, lab, col in UNITS:
        v = sc.aG_um_100[sc.key == key].values
        if v.size == 0:
            continue
        ax.hist(np.clip(v, bins[0], bins[-1] - 1e-6), bins, histtype="step", lw=2.4, color=col, density=True,
                label=f"{short.get(lab, lab)}: {np.median(v):.1f} µm")
    ax.axvline(A_G, color="k", ls="--", lw=2, label=f"adopted a$_G$ = {A_G} µm")
    ax.set_xlabel("a$_G$ in the lowest 100 m [µm]")
    ax.set_ylabel("fraction of calibration profiles")
    ax.set_xlim(0, 10)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, handlelength=1.4, fontsize=11)
    ax.set_title("Calibration profiles, 4 units")
    # (b) the Cloudnet sites of the domain, every other label to keep the axis readable
    y = np.arange(len(st))
    bx.barh(y, st.aG_med, color="#9ecae1", edgecolor="0.3", lw=0.5)
    bx.errorbar(st.aG_med, y, xerr=[st.aG_med - st.aG_p25, st.aG_p75 - st.aG_med], fmt="none", ecolor="0.25",
                elinewidth=1.0, capsize=2)
    bx.axvline(A_G, color="k", ls="--", lw=2, label=f"adopted a$_G$ = {A_G} µm")
    bx.axvline(med_sites, color="tab:red", ls=":", lw=2.2, label=f"median of the sites: {med_sites:.1f} µm")
    bx.set_yticks(y)
    bx.set_yticklabels([str(n) for n in st.name], fontsize=9.5)
    bx.set_xlabel("site median a$_G$ [µm], bars: IQR")
    bx.set_xlim(0, 8)
    bx.set_ylim(-0.7, len(st) - 0.3)
    bx.grid(alpha=0.3, axis="x")
    bx.legend(loc="lower right", frameon=True, framealpha=0.95, fontsize=11)
    bx.set_title(f"{len(st)} Cloudnet sites, 2025–2026")
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    stats = dict(n_profiles=int(pooled.size), n_days=int(sc.date.nunique()), aG_median_all=round(float(np.median(pooled)), 2),
                 aG_p10=round(float(np.percentile(pooled, 10)), 1), aG_p90=round(float(np.percentile(pooled, 90)), 1),
                 n_sites=int(len(st)), site_median=round(med_sites, 2),
                 site_q25=round(float(np.percentile(st.aG_med, 25)), 2), site_q75=round(float(np.percentile(st.aG_med, 75)), 2),
                 units={lab: dict(n=int((sc.key == key).sum()), median=round(float(np.median(sc.aG_um_100[sc.key == key])), 2))
                        for key, lab, _ in UNITS if (sc.key == key).any()})
    print("wrote", OUT_PNG)
    print("STATS", json.dumps(stats))


if __name__ == "__main__":
    main()
