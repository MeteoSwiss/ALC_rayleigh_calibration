"""Backup slide of the Annual Meeting 2026 deck, v2.3 edition ("Choosing the in-cloud extinction"): the alpha
ladder of the PVC eta tables at the adopted droplet radius a_G = 3.3 um, redrawn with slide-sized fonts from
the outputs of validation/eta_study/alpha_ladder.py (C:/DATA/Projects/202611_eta_study/alpha_ladder/: 43 Vaisala
units, 1.0 M stored in-cloud profiles rescaled exactly per profile; doc/reports/eta_switch_decision_2026-11-02.md
section 8). Three panels: the C_L change against the 2026-07 tables at cloud base 1 km, the CBH slope of ln C_L
(unit x month fixed effects, 95 % CI) and the CL61 cloud / Rayleigh ratio (Payerne with the hood dark v3, and the
median of the dark-free CL61 units).

Run from the repo root:  python presentation/slide9_alpha_ladder_v23.py
Output: presentation/figs_annual2026/fig_alpha_ladder_v23_slide.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
LADDER = Path("C:/DATA/Projects/202611_eta_study/alpha_ladder")
OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_alpha_ladder_v23_slide.png"
ALPHAS = [10, 20, 30]
TCOL = {"CL31": "#ff7f0e", "CL51": "#9467bd", "CL61": "#2ca02c"}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 17, "axes.titlesize": 18, "axes.labelsize": 17,
                         "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15})
    lev = pd.read_csv(LADDER / "level_by_cbh.csv")
    slo = pd.read_csv(LADDER / "cbh_slope.csv")
    rat = pd.read_csv(LADDER / "cl61_ratios.csv")
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17.28, 5.9))
    fig.subplots_adjust(left=0.055, right=0.995, top=0.9, bottom=0.15, wspace=0.3)
    for t, col in TCOL.items():
        y = [float(lev[(lev.itype == t) & (lev.bin == "1 km") & (lev.variant == f"ag3.3_a{a}")].change_pct.iloc[0]) for a in ALPHAS]
        a1.plot(ALPHAS, y, "o-", color=col, lw=2.4, ms=9, label=t)
        s = [slo[(slo.itype == t) & (slo.variant == f"ag3.3_a{a}")].iloc[0] for a in ALPHAS]
        a2.errorbar(np.array(ALPHAS) + {"CL31": -0.5, "CL51": 0, "CL61": 0.5}[t], [r.slope for r in s],
                    yerr=[[r.slope - r.lo for r in s], [r.hi - r.slope for r in s]], fmt="o-", color=col, lw=2.4,
                    ms=9, capsize=4, label=t)
    a1.axhline(0, color="k", lw=1)
    a1.axvline(20, color="k", ls=":", lw=1.5)
    a1.set_ylabel("C$_L$ change vs the July 2026 tables [%]")
    a1.set_title("Level at cloud base 1 km")
    a1.legend(frameon=False)
    a2.axhline(0, color="k", lw=1)
    a2.axvline(20, color="k", ls=":", lw=1.5)
    a2.set_ylabel("slope of ln C$_L$ vs cloud base [%/km]")
    a2.set_title("Cloud-base dependence (a check)")
    a2.legend(frameon=False)
    pay = rat[(rat.key == "0-20000-0-06610_C") & (rat.reference == "Rayleigh+dark v3")].iloc[0]
    free = rat[rat.reference == "Rayleigh no dark"]
    a3.axhspan(0.95, 1.05, color="0.85", lw=0, label="±5 %")
    a3.plot(ALPHAS, [pay[f"ag3.3_a{a}"] for a in ALPHAS], "s-", color="#1f77b4", lw=2.4, ms=9, label="Payerne, Rayleigh + hood dark")
    a3.plot(ALPHAS, [float(np.median(free[f"ag3.3_a{a}"])) for a in ALPHAS], "o-", color="0.25", lw=2.4, ms=9,
            label=f"median of {len(free)} dark-free CL61")
    a3.axhline(1, color="k", lw=1)
    a3.axvline(20, color="k", ls=":", lw=1.5)
    a3.set_ylabel("CL61 cloud C$_L$ / Rayleigh C$_L$")
    a3.set_title("CL61 cloud against Rayleigh")
    a3.legend(frameon=True, framealpha=0.95, loc="upper right")
    for ax in (a1, a2, a3):
        ax.set_xticks(ALPHAS)
        ax.set_xlabel("in-cloud extinction α of the table [km⁻¹]")
        ax.grid(alpha=0.3)
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print("wrote", OUT_PNG)
    print("level 1 km:", {t: [float(lev[(lev.itype == t) & (lev.bin == "1 km") & (lev.variant == f"ag3.3_a{a}")].change_pct.iloc[0]) for a in ALPHAS] for t in TCOL})
    print("slopes a20:", {t: float(slo[(slo.itype == t) & (slo.variant == "ag3.3_a20")].slope.iloc[0]) for t in TCOL})
    print("CL61 ratio a20: Payerne+dark", float(pay["ag3.3_a20"]), "dark-free median", float(np.median(free["ag3.3_a20"])), "n", len(free))


if __name__ == "__main__":
    main()
