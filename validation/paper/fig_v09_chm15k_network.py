"""Section 5.1.1 figure of the paper (v0.9): network ranking of the Rayleigh lidar constants C_L of the CHM15k
(operational / dashboard version: eprof_v2.2), the counterpart of the liquid-cloud ranking of Figure 5.

Source: the network v2.2 run calout_v22_04 (433 streams, 2025-01-01 .. 2026-08-13), verified identical to the
dashboard nightly constants for the CHM15k (median ratio 1.00000). Per instrument: median of the successful
nightly constants (flags 1 and 0.5, at least MIN_NIGHTS nights), bar = interquartile range over the nights,
which includes the seasonal cycle of the constant. The CHM15k units discussed in the paper are labelled.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v09")
DEFAULT_CL = 3e11              # E-PROFILE default constant of the CHM15k (Section 4.3)
UNIT = 1e11
MIN_NIGHTS = 10
PAPER_UNITS = {                # CHM15k discussed in the paper -> short label
    "0-20000-0-06610_A": "Payerne",
    "0-20000-0-06736_A": "Kleine Scheidegg",
    "0-20000-0-10962_0": "Hohenpeissenberg",
    "0-20000-0-10471_0": "Leipzig",
    "0-20008-0-INO_A": "Măgurele A",
    "0-20008-0-INO_B": "Măgurele B",
    "0-20000-0-06240_A": "Schiphol A",
    "0-20000-0-06240_B": "Schiphol B",
    "0-20000-0-06240_C": "Schiphol C",
    "0-20000-0-06240_D": "Schiphol D",
    "0-250-1001-07151_B": "Palaiseau",
}


def collect():
    items = json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))
    rows = []
    for e in items:
        if e.get("type") != "CHM15k":
            continue
        key = f"{e['wmo']}_{e['ident']}"
        f = RUN / key / f"{key}_cal.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        d = d[(d.method == "rayleigh") & d.flag.isin([1, 0.5]) & (d.cal_value > 0)]
        c = d.cal_value.values / UNIT
        if c.size >= MIN_NIGHTS:
            rows.append(dict(key=key, site=e.get("site", ""), n=c.size, med=np.median(c),
                             q25=np.percentile(c, 25), q75=np.percentile(c, 75)))
    return pd.DataFrame(rows).sort_values("med").reset_index(drop=True)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = collect()
    df.to_csv(OUT / "fig_v09_chm15k_network.csv", index=False)
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(17, 6.2), dpi=200)
    ax.bar(x, df.med, width=0.85, color="#c0392b", label=f"CHM15k ({len(df)})")
    ax.errorbar(x, df.med, yerr=[df.med - df.q25, df.q75 - df.med], fmt="none", ecolor="0.25",
                elinewidth=0.6, capsize=0)
    net_med = float(np.median(df.med))
    ax.axhline(DEFAULT_CL / UNIT, color="red", ls="--", lw=1.4, label="default constant (3 × 10$^{11}$)")
    ax.axhline(net_med, color="blue", ls="--", lw=1.4, label=f"network median ({net_med:.2f} × 10$^{{11}}$)")
    top = max(float(df.q75.quantile(0.98)) * 1.15, DEFAULT_CL / UNIT * 1.3)
    ax.set_ylim(0, top)
    ax.set_xlim(-1, len(df))
    last_i, lift = -99, 0.0
    for i, r in df.iterrows():
        if r.key in PAPER_UNITS:
            lift = (top * 0.17 if lift == 0.0 else 0.0) if i - last_i <= 3 else 0.0   # de-overlap neighbours
            y = min(r.q75, top * 0.62) + top * 0.015
            ax.plot(i, y, marker="v", color="k", ms=4)
            if lift:
                ax.plot([i, i], [y + top * 0.01, y + lift], color="0.5", lw=0.5)
            ax.text(i, y + top * 0.015 + lift, PAPER_UNITS[r.key], rotation=90, ha="center", va="bottom",
                    fontsize=8)
            last_i = i
    n_off = int((df.q75 > top).sum())
    if n_off:
        ax.text(0.93, 0.97, f"{n_off} bar{'s' if n_off > 1 else ''} clipped at the top", transform=ax.transAxes,
                ha="right", va="top",
                fontsize=9, color="0.3", bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))
    ax.set_xlabel("instrument (sorted by lidar constant)")
    ax.set_ylabel("Rayleigh lidar constant $C_L$  (10$^{11}$)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_v09_chm15k_network.png")
    print(f"saved {OUT / 'fig_v09_chm15k_network.png'}")
    q = df.med.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(2).to_dict()
    print(f"instruments {len(df)}, median nights {int(df.n.median())}, median C_L quantiles (1e11) {q}")
    print(f"range {df.med.min():.2f}-{df.med.max():.2f}; below default 3e11: {(df.med < 3).mean():.2f}; "
          f"within x2 of the median: {((df.med > net_med / 2) & (df.med < net_med * 2)).mean():.2f}; "
          f"median IQR/median {np.median((df.q75 - df.q25) / df.med):.2f}; |C_L/default-1|>0.2: "
          f"{(abs(df.med / (DEFAULT_CL / UNIT) - 1) > 0.2).mean():.2f}; outside x2 of default: "
          f"{((df.med < DEFAULT_CL / UNIT / 2) | (df.med > DEFAULT_CL / UNIT * 2)).mean():.3f}")
    for k, lab in PAPER_UNITS.items():
        m = df[df.key == k]
        if len(m):
            r = m.iloc[0]
            print(f"  {lab:18s} rank {int(m.index[0]) + 1}/{len(df)}  C_L {r.med:.2f}e11 (IQR {r.q25:.2f}-{r.q75:.2f}), n={r.n}")
        else:
            print(f"  {lab:18s} not in the ranking")


if __name__ == "__main__":
    main()
