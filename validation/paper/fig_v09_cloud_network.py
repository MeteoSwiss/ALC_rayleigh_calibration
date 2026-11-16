"""Figure 5 of the paper (v0.9): network ranking of the liquid-cloud LIDAR CONSTANTS C_L of the operational
(dashboard) version, C_L = C_0 / C (Wiegner convention, beta_att = rcs_0 / C_L).

C = C_0 / C_L is the O'Connor coefficient relative to the constant C_0 already applied to the
Level-1 signal (the E-PROFILE default 1e8 for the CL31/CL51, 1 for the CL61 whose L1 signal is
calibrated by the manufacturer): C > 1 means that the input signal under-estimates the backscatter.

Sources (both equal to the dashboard nightly values):
  * CL31/CL51: the network v2.2 run calout_v22_04 (433 streams, 2025-01-01 .. 2026-08-13; verified
    identical to the dashboard nightly constants, median ratio 1.00000);
  * CL61: the dashboard series fetched by fetch_dashboard_calib.py (calib_dashboard_cl61/), because
    calout_v22_04 predates the switch of the CL61 emission spectrum to 910.55 nm / 0.188 nm.
Per instrument: median of the daily coefficients (>= MIN_DAYS days), bar = interquartile range.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04")
CL61_DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v09/calib_dashboard_cl61")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v09")
C0 = {"CL31": 1e8, "CL51": 1e8, "CL61": 1.0}
MIN_DAYS = 10
COLORS = {"CL31": "#1f77b4", "CL51": "#ff7f0e", "CL61": "#2ca02c"}


def census_types():
    items = json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))
    return {f"{e['wmo']}_{e['ident']}": (e.get("type"), e.get("site", "")) for e in items}


def collect():
    types = census_types()
    rows = []
    for key, (itype, site) in types.items():
        if itype not in ("CL31", "CL51"):
            continue
        f = RUN / key / f"{key}_cal.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        d = d[(d.method == "cloud") & d.flag.isin([1, 0.5]) & (d.cal_value > 0)]
        c = d.cal_value.values / C0[itype]            # C_L in units of the default constant C_0
        if c.size >= MIN_DAYS:
            rows.append(dict(key=key, site=site, itype=itype, n=c.size, med=np.median(c),
                             q25=np.percentile(c, 25), q75=np.percentile(c, 75),
                             first=str(d.date.min()), last=str(d.date.max())))
    for f in sorted(CL61_DIR.glob("*_cloud.csv")):
        key = f.name[:-len("_cloud.csv")]
        d = pd.read_csv(f).dropna(subset=["C_daily"])
        d = d[d.C_daily > 0]
        c = d.C_daily.values / C0["CL61"]            # C_L of the CL61 (C_0 = 1)
        if c.size >= MIN_DAYS:
            rows.append(dict(key=key, site=types.get(key, ("", ""))[1], itype="CL61", n=c.size,
                             med=np.median(c), q25=np.percentile(c, 25), q75=np.percentile(c, 75),
                             first=str(d.time.min()), last=str(d.time.max())))
    return pd.DataFrame(rows).sort_values("med").reset_index(drop=True)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = collect()
    df.to_csv(OUT / "fig_v09_cloud_network.csv", index=False)
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(17, 6.2), dpi=200)
    for t in ("CL31", "CL51", "CL61"):
        m = (df.itype == t).values
        ax.bar(x[m], df.med[m], width=0.85, color=COLORS[t], label=f"{t} ({int(m.sum())})")
        ax.errorbar(x[m], df.med[m], yerr=[df.med[m] - df.q25[m], df.q75[m] - df.med[m]],
                    fmt="none", ecolor="0.25", elinewidth=0.6, capsize=0)
    net_med = float(np.median(df.med))
    ax.axhline(1.0, color="red", ls="--", lw=1.4, label="default constant (10$^8$ CL31/CL51, 1 CL61)")
    ax.axhline(net_med, color="blue", ls="--", lw=1.4, label=f"network median ({net_med:.2f})")
    ax.set_xlim(-1, len(df))
    top = 1.6
    ax.set_ylim(0, top)
    ax.set_xlabel("instrument (sorted by lidar constant)")
    ax.set_ylabel("$C_L$ of the CL31 and CL51  (10$^8$)")
    ax2 = ax.twinx()
    ax2.set_ylim(0, top)
    ax2.set_ylabel("$C_L$ of the CL61")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    n_off = int((df.q75 > top).sum())
    if n_off:
        ax.text(0.93, 0.97, f"{n_off} bars clipped at the top", transform=ax.transAxes, ha="right", va="top",
                fontsize=9, color="0.3", bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))
    fig.tight_layout()
    fig.savefig(OUT / "fig_v09_cloud_network.png")
    print(f"saved {OUT / 'fig_v09_cloud_network.png'}")
    by = df.groupby("itype").med.describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]).round(3)
    print(by.to_string())
    vai = df[df.itype.isin(["CL31", "CL51"])]
    cl61 = df[df.itype == "CL61"]
    print(f"instruments {len(df)}; network median C_L/C_0 {net_med:.3f}; CL31/CL51 below the default: "
          f"{(vai.med < 1).mean():.3f} ({int((vai.med < 1).sum())}/{len(vai)}); CL31/CL51 C_L/1e8 range "
          f"{vai.med.min():.2f}-{vai.med.max():.2f}; CL61 range {cl61.med.min():.3f}-{cl61.med.max():.3f}, median "
          f"{cl61.med.median():.3f}; median n days {int(df.n.median())}")

if __name__ == "__main__":
    main()
