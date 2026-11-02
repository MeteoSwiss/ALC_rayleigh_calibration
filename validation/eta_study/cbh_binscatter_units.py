# -*- coding: utf-8 -*-
"""Binned scatter of the cloud lidar constant against cloud-base height, per instrument type and per
unit, for the three eta variants of the study (a = operational PVC 5.5 um / 10 /km, b = legacy Hopkin,
c = proposed PVC a_G 3.3 um / network alpha_k). Answers "does the proposed table remove the CBH
dependence of C_L, unit by unit?".

Input: the DAILY medians written by cbh_slope_windowed.py --stats (one row per unit-day: median CBH and
median ln C_L per variant, from the stored in-cloud windows of every accepted profile), for the 270
Vaisala units: slope/pc_study (43), slope/pc_network (8), slope/linard (219).

Per unit, ln C_L of each variant is centred on the unit's own median (the constant drops out, only the
CBH shape is compared), then binned in 0.25 km of CBH. The per-unit slope is the Theil-Sen slope of the
daily points (the "daily" statistic of report 04 v2 section 5.1).

Outputs (doc/reports/figs_04v2/ unless --out):
  fig_cbh_binscatter_types.png        pooled per type, three variants, binned median and IQR
  fig_cbh_binscatter_CL61_units.png   one panel per CL61 unit (+ the distorted Payerne CL31)
  fig_cbh_binscatter_CL51_units.png   one panel per CL51 unit
  fig_cbh_binscatter_CL31_units_p<N>.png   CL31 units, 42 per page
  fig_cbh_slope_per_unit_hist.png     distribution of the per-unit slopes per type and variant
  C:/DATA/Projects/202611_eta_study/slope/network/cbh_slope_per_unit.csv

Run:  python validation/eta_study/cbh_binscatter_units.py [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
STUDY = Path("C:/DATA/Projects/202611_eta_study")
STATS_DIRS = [STUDY / "slope" / "pc_study", STUDY / "slope" / "pc_network", STUDY / "slope" / "linard"]
VARS = {"a": ("operational (PVC 5.5 um, 10 /km)", "#404040", "-"),
        "b": ("legacy Hopkin", "#ff7f0e", "--"),
        "c": ("proposed (PVC 3.3 um, alpha_k)", "#2ca02c", "-")}
EDGES = np.arange(0.0, 3.01, 0.25)
MID = 0.5 * (EDGES[:-1] + EDGES[1:])
MIN_DAYS = 5
PAYERNE_CL31 = "0-20000-0-06610_B"


def load_daily():
    frames = []
    for d in STATS_DIRS:
        f = d / "daily_lnC.csv.gz"
        if f.exists():
            frames.append(pd.read_csv(f))
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates(["key", "date"], keep="first")          # pc_study wins over linard
    for v in VARS:
        df[f"d_{v}"] = df[f"lnC_{v}"] - df.groupby("key")[f"lnC_{v}"].transform("median")
    return df


def site_names():
    names = {}
    f = STUDY / "cscs_unit_list.csv"
    if f.exists():
        with open(f, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                names[r["key"]] = r["site"]
    try:
        for r in json.loads((REPO / "validation" / "scope_l1_2026_census.json").read_text(encoding="utf-8")):
            names.setdefault(f"{r['wmo']}_{r['ident']}", r["site"])
    except Exception:
        pass
    return names


def binned(x, y):
    med = np.full(MID.size, np.nan)
    q1 = np.full(MID.size, np.nan)
    q3 = np.full(MID.size, np.nan)
    n = np.zeros(MID.size, int)
    for i in range(MID.size):
        m = (x >= EDGES[i]) & (x < EDGES[i + 1]) & np.isfinite(y)
        n[i] = m.sum()
        if n[i] >= 3:
            med[i], q1[i], q3[i] = np.median(y[m]), np.percentile(y[m], 25), np.percentile(y[m], 75)
    return med, q1, q3, n


def unit_slopes(df):
    from scipy.stats import theilslopes
    rows = []
    for key, g in df.groupby("key"):
        if len(g) < MIN_DAYS or g.cbh.nunique() < 3:
            continue
        rec = dict(key=key, itype=g.itype.iloc[0], n_days=len(g), cbh_med_km=float(g.cbh.median()))
        for v in VARS:
            sl = theilslopes(g[f"d_{v}"].to_numpy(), g.cbh.to_numpy())
            rec[f"slope_{v}"] = 100 * sl[0]
            rec[f"slope_{v}_lo"] = 100 * sl[2]
            rec[f"slope_{v}_hi"] = 100 * sl[3]
        rows.append(rec)
    return pd.DataFrame(rows)


def fig_types(df, slopes, out, plt):
    fe = pd.read_csv(STUDY / "slope" / "network" / "cbh_slope_network.csv")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), sharey=True)
    for ax, it in zip(axes, ("CL31", "CL51", "CL61")):
        g = df[df.itype == it]
        for v, (lab, c, ls) in VARS.items():
            med, q1, q3, n = binned(g.cbh.to_numpy(), g[f"d_{v}"].to_numpy())
            ok = n >= 20
            ax.fill_between(MID[ok], 100 * q1[ok], 100 * q3[ok], color=c, alpha=0.10, lw=0)
            r = fe[(fe.group == f"{it} all") & (fe.variant == v)]
            s_fe = f"{float(r.unit_month_FE.iloc[0]):+.1f} [{float(r.unit_month_lo.iloc[0]):+.1f}, {float(r.unit_month_hi.iloc[0]):+.1f}] %/km" if len(r) else ""
            ax.plot(MID[ok], 100 * med[ok], ls, color=c, lw=2.2, marker="o", ms=5, label=f"{lab}: {s_fe}")
        ax.axhline(0, color="grey", lw=1)
        ax.set_xlabel("cloud-base height [km AGL], daily median")
        ax.set_title(f"{it}: {g.key.nunique()} units, {len(g)} unit-days", fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower left", title="per-profile slope, unit x month fixed effects", title_fontsize=8)
        ax.set_xlim(0, 3)
    axes[0].set_ylabel("ln C_L - unit median [%]  (binned median, shaded IQR)")
    axes[0].set_ylim(-30, 30)
    fig.suptitle("Cloud lidar constant against cloud-base height, all Vaisala units of the network, per eta variant "
                 "(a check, not a fit: the inputs were fixed from Cloudnet and the in-cloud decay)", fontsize=11.5)
    fig.tight_layout()
    fig.savefig(out / "fig_cbh_binscatter_types.png", dpi=140)
    plt.close(fig)


def _panel(ax, g, sl, title, points=True, small=False):
    if points:
        ax.plot(g.cbh, 100 * g.d_a, ".", color="#bbbbbb", ms=3, alpha=0.6)
    for v in ("b", "a", "c"):
        lab, c, ls = VARS[v]
        med, q1, q3, n = binned(g.cbh.to_numpy(), g[f"d_{v}"].to_numpy())
        ok = n >= 3
        ax.plot(MID[ok], 100 * med[ok], ls, color=c, lw=1.2 if small else 2, marker="o", ms=2 if small else 4)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlim(0, 3)
    ax.set_ylim(-40, 40)
    ax.grid(alpha=0.25)
    if sl is not None:
        txt = f"{title}\nn={int(sl.n_days)} d  a {sl.slope_a:+.1f}  c {sl.slope_c:+.1f} %/km"
    else:
        txt = title
    ax.set_title(txt, fontsize=7 if small else 9)
    ax.tick_params(labelsize=6 if small else 8)


def fig_units(df, slopes, names, out, plt, itype, keys, fname, ncols, nrows, points, suptitle):
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax in axes:
        ax.set_visible(False)
    for ax, key in zip(axes, keys):
        ax.set_visible(True)
        g = df[df.key == key]
        sl = slopes[slopes.key == key]
        sl = sl.iloc[0] if len(sl) else None
        tag = names.get(key, key)[:22]
        if key == PAYERNE_CL31:
            tag += " (distorted CL31)"
        _panel(ax, g, sl, f"{tag} [{g.itype.iloc[0]}]", points=points, small=(nrows * ncols > 16))
    for ax in axes[-ncols:]:
        if ax.get_visible():
            ax.set_xlabel("CBH [km]", fontsize=8)
    fig.suptitle(suptitle + "   (black = operational, green = proposed, orange dashed = legacy Hopkin; "
                 "y = ln C_L - unit median [%])", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out / fname, dpi=120)
    plt.close(fig)
    print("->", out / fname)


def fig_hist(slopes, out, plt):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), sharex=True)
    edges = np.arange(-20, 20.1, 1.0)
    for ax, it in zip(axes, ("CL31", "CL51", "CL61")):
        s = slopes[slopes.itype == it]
        for v, (lab, c, ls) in VARS.items():
            x = s[f"slope_{v}"].to_numpy()
            ax.hist(x, edges, histtype="step", color=c, lw=2, ls=ls,
                    label=f"{lab}: median {np.median(x):+.1f}, IQR {np.percentile(x, 25):+.1f}..{np.percentile(x, 75):+.1f} %/km")
        ax.axvline(0, color="grey", lw=1)
        ax.set_title(f"{it}: {len(s)} units with >= {MIN_DAYS} days", fontsize=11)
        ax.set_xlabel("per-unit Theil-Sen slope of daily ln C_L vs CBH [%/km]")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("number of units")
    fig.suptitle("Distribution of the per-unit CBH slopes per eta variant", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "fig_cbh_slope_per_unit_hist.png", dpi=140)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "doc" / "reports" / "figs_04v2"))
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    df = load_daily()
    names = site_names()
    slopes = unit_slopes(df)
    slopes["site"] = [names.get(k, "") for k in slopes.key]
    slopes.to_csv(STUDY / "slope" / "network" / "cbh_slope_per_unit.csv", index=False, float_format="%.3f")
    print(f"{df.key.nunique()} units, {len(df)} unit-days; slopes for {len(slopes)} units")
    for it in ("CL31", "CL51", "CL61"):
        s = slopes[slopes.itype == it]
        print(f"  {it}: n {len(s)}  median slope a {s.slope_a.median():+.2f}  b {s.slope_b.median():+.2f}  c {s.slope_c.median():+.2f} %/km"
              f" | units within +-2 %/km: a {(s.slope_a.abs() <= 2).sum()}  c {(s.slope_c.abs() <= 2).sum()}"
              f" | |slope| median a {s.slope_a.abs().median():.2f}  c {s.slope_c.abs().median():.2f}")
    fig_types(df, slopes, out, plt)
    fig_hist(slopes, out, plt)
    keys61 = sorted(df[df.itype == "CL61"].key.unique()) + [PAYERNE_CL31]
    fig_units(df, slopes, names, out, plt, "CL61", keys61, "fig_cbh_binscatter_CL61_units.png", 4, 3, True,
              "CL61 units (and the distorted Payerne CL31): daily cloud C_L vs CBH")
    keys51 = sorted(df[df.itype == "CL51"].key.unique())
    fig_units(df, slopes, names, out, plt, "CL51", keys51, "fig_cbh_binscatter_CL51_units.png", 8, 6, False,
              "CL51 units: daily cloud C_L vs CBH")
    keys31 = sorted(df[df.itype == "CL31"].key.unique())
    per = 42
    for p in range(0, len(keys31), per):
        fig_units(df, slopes, names, out, plt, "CL31", keys31[p:p + per], f"fig_cbh_binscatter_CL31_units_p{p // per + 1}.png",
                  7, 6, False, f"CL31 units {p + 1}-{min(p + per, len(keys31))} of {len(keys31)}: daily cloud C_L vs CBH")


if __name__ == "__main__":
    main()
