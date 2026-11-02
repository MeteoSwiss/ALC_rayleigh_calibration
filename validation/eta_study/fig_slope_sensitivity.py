"""Two summary figures of the eta study from existing CSVs (no computation):
  fig_cbh_slope_variants.png   de-seasoned CBH slope of ln C_L per eta variant and type (slope/network/cbh_slope_network.csv)
  fig_sensitivity.png          sensitivities (not fits) of slope and 1-km level per input (sensitivity/sensitivity_table.csv)

Run:  python validation/eta_study/fig_slope_sensitivity.py [--out DIR]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

BASE = Path("C:/DATA/Projects/202611_eta_study")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(BASE / "fr_figs"))
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)

    s = pd.read_csv(BASE / "slope" / "network" / "cbh_slope_network.csv")
    s = s[s.group.isin(["CL31 all", "CL51 all", "CL61 all"])]
    names = {"a": "operational (PVC 5.5 um, 10 /km)", "b": "legacy Hopkin", "c": "proposed (PVC 3.3 um, alpha_k)"}
    cols = {"a": "k", "b": "C1", "c": "C2"}
    fig, ax = plt.subplots(figsize=(12, 5.2))
    x = np.arange(3)
    for i, v in enumerate(("a", "b", "c")):
        g = s[s.variant == v].set_index("group").loc[["CL31 all", "CL51 all", "CL61 all"]]
        err = np.vstack([g.unit_month_FE - g.unit_month_lo, g.unit_month_hi - g.unit_month_FE])
        ax.bar(x + (i - 1) * 0.27, g.unit_month_FE, 0.25, yerr=err, color=cols[v], alpha=0.8, capsize=4, label=names[v])
        for xx, val in zip(x + (i - 1) * 0.27, g.unit_month_FE):
            ax.text(xx, val + (0.3 if val >= 0 else -0.6), f"{val:+.1f}", ha="center", fontsize=9)
    ax.axhline(0, color="grey", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t}\n({int(n)} units)" for t, n in zip(("CL31", "CL51", "CL61"),
                        s[s.variant == "a"].set_index("group").loc[["CL31 all", "CL51 all", "CL61 all"]].n_units)])
    ax.set_ylabel("slope of ln C_L vs cloud-base height [%/km]")
    ax.set_title("CBH slope of the cloud lidar constant per eta variant (per profile, unit x month fixed effects,\n"
                 "95 % CI clustered by unit-day; a check, not a fit)", fontsize=11)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_cbh_slope_variants.png", dpi=140)
    plt.close(fig)

    t = pd.read_csv(BASE / "sensitivity" / "sensitivity_table.csv")
    t = t[~t.input.isin(["operational", "central"])]
    inputs = list(dict.fromkeys(t.input))
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True)
    y = np.arange(len(inputs))
    for i, (it, c) in enumerate((("CL31", "C0"), ("CL51", "C2"), ("CL61", "C3"))):
        g = t[t.itype == it].set_index("input").reindex(inputs)
        axes[0].barh(y + (i - 1) * 0.27, g.d_slope_pct_km, 0.25, color=c, label=it)
        axes[1].barh(y + (i - 1) * 0.27, g.d_level_1km_pct, 0.25, color=c, label=it)
    for ax, lab in zip(axes, ("change of CBH slope [%/km]", "change of C_L at 1 km [%]")):
        ax.axvline(0, color="grey", lw=1)
        ax.set_xlabel(lab)
        ax.grid(alpha=0.3, axis="x")
        ax.legend(fontsize=9)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(inputs)
    fig.suptitle("Sensitivities (not fits) of the cloud calibration to the eta inputs, around a_G 3.3 um / network alpha_k",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / "fig_sensitivity.png", dpi=140)
    plt.close(fig)
    print("written", out, flush=True)


if __name__ == "__main__":
    main()
