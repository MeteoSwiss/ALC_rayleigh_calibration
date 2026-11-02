# -*- coding: utf-8 -*-
"""Figure of the alpha ladder (alpha_ladder.py outputs): C_L level vs the 2026-07 tables, CBH slope and the
CL61 cloud/Rayleigh ratio, against alpha = 10 / 20 / 30 /km for a_G 3.3 and 3.7 um.

Run:  python validation/eta_study/fig_alpha_ladder.py [--out doc/reports/figs_04v2]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
LAD = Path("C:/DATA/Projects/202611_eta_study/alpha_ladder")
TYPE_COL = {"CL31": "#ff7f0e", "CL51": "#9467bd", "CL61": "#2ca02c"}
ALPHAS = (10, 20, 30)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "doc" / "reports" / "figs_04v2"))
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    lev = pd.read_csv(LAD / "level_by_cbh.csv")
    slo = pd.read_csv(LAD / "cbh_slope.csv")
    rat = pd.read_csv(LAD / "cl61_ratios.csv")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))
    for ag, ls, lab in ((3.3, "-", "a_G 3.3 um"), (3.7, "--", "a_G 3.7 um")):
        ax = axes[0]
        for it in ("CL31", "CL51", "CL61"):
            y = [float(lev[(lev.itype == it) & (lev.bin == "1 km") & (lev.variant == f"ag{ag:.1f}_a{a}")].change_pct.iloc[0]) for a in ALPHAS]
            ax.plot(ALPHAS, y, ls, marker="o", color=TYPE_COL[it], lw=2, ms=7, label=f"{it}, {lab}")
        ax = axes[1]
        for it in ("CL31", "CL51", "CL61"):
            s = slo[(slo.itype == it)]
            y = [float(s[s.variant == f"ag{ag:.1f}_a{a}"].slope.iloc[0]) for a in ALPHAS]
            lo = [float(s[s.variant == f"ag{ag:.1f}_a{a}"].lo.iloc[0]) for a in ALPHAS]
            hi = [float(s[s.variant == f"ag{ag:.1f}_a{a}"].hi.iloc[0]) for a in ALPHAS]
            ax.errorbar(np.array(ALPHAS) + (0.4 if ag == 3.7 else -0.4), y, yerr=[np.subtract(y, lo), np.subtract(hi, y)],
                        fmt=ls, marker="o", color=TYPE_COL[it], lw=2, ms=7, capsize=3, label=f"{it}, {lab}")
        ax = axes[2]
        pay = rat[(rat.key == "0-20000-0-06610_C") & (rat.reference == "Rayleigh+dark v3")].iloc[0]
        oth = rat[(rat.reference == "Rayleigh no dark") & (rat.key != "0-20000-0-06610_C")]
        ax.plot(ALPHAS, [pay[f"ag{ag:.1f}_a{a}"] for a in ALPHAS], ls, marker="s", color="#1f77b4", lw=2, ms=7,
                label=f"Payerne CL61 vs Rayleigh + hood dark v3, {lab}")
        ax.plot(ALPHAS, [oth[f"ag{ag:.1f}_a{a}"].median() for a in ALPHAS], ls, marker="o", color="#404040", lw=2, ms=7,
                label=f"median of 10 dark-free CL61 units, {lab}")
    axes[0].axhline(0, color="grey", lw=1)
    axes[0].set_ylabel("C_L change vs the 2026-07 tables at CBH 1 km [%]")
    axes[0].set_title("level: 43 study units, median per-profile ratio", fontsize=11)
    axes[1].axhline(0, color="grey", lw=1)
    axes[1].set_ylabel("slope of ln C_L vs CBH [%/km], unit x month FE")
    axes[1].set_title("CBH dependence (a check, not a fit)", fontsize=11)
    axes[2].axhline(1, color="grey", lw=1)
    axes[2].axhspan(0.95, 1.05, color="grey", alpha=0.12, lw=0)
    axes[2].set_ylabel("cloud C_L / Rayleigh C_L (ratio of medians)")
    axes[2].set_title("CL61 cloud against its own Rayleigh constant", fontsize=11)
    for ax in axes:
        ax.set_xlabel("in-cloud extinction alpha of the eta table [1/km]")
        ax.set_xticks(ALPHAS)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    for ax, txt in ((axes[0], "operational 2026-07: a_G 5.5 um, alpha 10"), (axes[1], "measured tail decay: 28-31 /km (lower bound); O'Connor 2004 assumed 15-20")):
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=8, va="top", color="#555555")
    fig.suptitle("Alpha ladder of the PVC eta tables: exact per-profile rescaling of the stored in-cloud windows (1.0 M profiles)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "fig_alpha_ladder.png", dpi=140)
    plt.close(fig)
    print("->", out / "fig_alpha_ladder.png")


if __name__ == "__main__":
    main()
