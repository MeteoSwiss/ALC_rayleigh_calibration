"""Network-wide alpha_k: cloud or receiver?

Merges the per-unit summaries of the three runs (collect_unit_summary.py output):
  PC study run    C:/DATA/Projects/202611_eta_study/summary/
  PC network run  C:/DATA/Projects/202611_eta_study/network/pc/summary/
  Linard array    C:/DATA/Projects/202611_eta_study/network/linard/summary/  (copied back)
checks them against network/split.csv (every unit exactly once), then asks whether alpha_k
follows the CLIMATE (latitude, region, season) or the INSTRUMENT (type, gate spacing, firmware):
  - group medians and spreads per type / gate / firmware / region,
  - a rank-based variance split: share of the between-unit variance of alpha_k explained by each
    factor alone (eta^2 of a one-way ANOVA on the unit medians), and by latitude (Spearman),
  - within-type seasonality from the unit monthly medians.

Run:  python validation/eta_study/network_alpha_k_analysis.py
Outputs: C:/DATA/Projects/202611_eta_study/network/ (unit_summary_all.csv, factor_table.csv,
         fig_network_alpha_k.png, summary.txt)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

BASE = Path("C:/DATA/Projects/202611_eta_study")
SOURCES = {"PC-study": BASE / "summary", "PC-network": BASE / "network" / "pc" / "summary",
           "Linard": BASE / "network" / "linard" / "summary"}
OUT = BASE / "network"
REGION = {  # coarse climate regions from the L1 site_location country
    "ICELAND": "subarctic maritime", "NORWAY": "subarctic maritime", "SWEDEN": "northern Europe",
    "FINLAND": "northern Europe", "DENMARK": "northern Europe", "UK": "NW Atlantic", "IRELAND": "NW Atlantic",
    "FRANCE": "western Europe", "BELGIUM": "western Europe", "NETHERLANDS": "western Europe",
    "LUXEMBOURG": "western Europe", "GERMANY": "central Europe", "SWITZERLAND": "central Europe",
    "AUSTRIA": "central Europe", "CZECH_REPUBLIC": "central Europe", "CZECHIA": "central Europe",
    "SLOVAKIA": "central Europe", "POLAND": "central Europe", "HUNGARY": "central Europe",
    "SLOVENIA": "central Europe", "SPAIN": "Iberia/Mediterranean", "PORTUGAL": "Iberia/Mediterranean",
    "ITALY": "Iberia/Mediterranean", "GREECE": "Iberia/Mediterranean", "CYPRUS": "Iberia/Mediterranean",
    "CROATIA": "central Europe", "ROMANIA": "central Europe", "BULGARIA": "central Europe",
    "ESTONIA": "northern Europe", "LATVIA": "northern Europe", "LITHUANIA": "northern Europe",
}


def eta2(values, groups):
    """Share of the variance of `values` explained by the grouping (one-way ANOVA eta^2)."""
    v = np.asarray(values, float)
    g = np.asarray(groups)
    ok = np.isfinite(v)
    v, g = v[ok], g[ok]
    if v.size < 3:
        return np.nan
    ss_tot = ((v - v.mean()) ** 2).sum()
    ss_b = sum(((v[g == k].mean() - v.mean()) ** 2) * (g == k).sum() for k in np.unique(g))
    return float(ss_b / ss_tot) if ss_tot > 0 else np.nan


def main():
    import pandas as pd
    from scipy.stats import spearmanr
    log = open(OUT / "summary.txt", "w", encoding="utf-8")

    def say(s=""):
        print(s, flush=True)
        log.write(s + "\n")

    parts, months = [], []
    for name, d in SOURCES.items():
        f = d / "unit_summary.csv"
        if f.exists():
            t = pd.read_csv(f)
            t["run"] = name
            parts.append(t)
            m = d / "unit_month_alpha_k.csv"
            if m.exists():
                months.append(pd.read_csv(m))
        else:
            say(f"missing summary: {f}")
    df = pd.concat(parts, ignore_index=True)
    # a run lists every unit it was GIVEN; units moved elsewhere show status "no output" there --
    # keep only the run that actually produced each unit, then check against split.csv
    df = df[df.status == "ok"].copy()
    split = pd.read_csv(OUT / "split.csv")
    dup = df.key[df.key.duplicated()].tolist()
    gaps = sorted(set(split.key) - set(df.key))
    say(f"units: {len(df)} with output, split.csv {len(split)}; duplicates {dup}; missing {len(gaps)} {gaps[:10]}")
    df = df.drop_duplicates("key")
    df["region"] = df.country.map(REGION).fillna("other")
    df["fw_main"] = df.firmware.astype(str).str.split(";").str[-1]
    df.to_csv(OUT / "unit_summary_all.csv", index=False, float_format="%.4g")

    say("\n== alpha_k per type (unit medians) ==")
    say(df.groupby("itype").alpha_k_med.describe(percentiles=[0.1, 0.5, 0.9]).round(2).to_string())
    for col in ("gate_m", "fw_main", "region"):
        say(f"\n== alpha_k by {col} (per type; unit medians) ==")
        t = df.groupby(["itype", col]).alpha_k_med.agg(["size", "median", "min", "max"]).round(2)
        say(t.to_string())

    say("\n== variance of the unit-median alpha_k explained by each factor (eta^2, within type) ==")
    rows = []
    for it, g in df.groupby("itype"):
        r = dict(itype=it, n_units=len(g), alpha_k_med=float(g.alpha_k_med.median()),
                 sd_units=float(g.alpha_k_med.std()))
        for col in ("gate_m", "fw_main", "region", "country"):
            r[f"eta2_{col}"] = eta2(g.alpha_k_med, g[col])
            r[f"n_levels_{col}"] = int(g[col].nunique())
        r["spearman_lat"] = float(spearmanr(g.lat, g.alpha_k_med).statistic) if len(g) > 4 else np.nan
        r["spearman_alt"] = float(spearmanr(g.alt, g.alpha_k_med).statistic) if len(g) > 4 else np.nan
        r["spearman_cbh"] = float(spearmanr(g.cbh_med, g.alpha_k_med).statistic) if len(g) > 4 else np.nan
        rows.append(r)
    allr = dict(itype="all", n_units=len(df), alpha_k_med=float(df.alpha_k_med.median()),
                sd_units=float(df.alpha_k_med.std()))
    for col in ("itype", "gate_m", "fw_main", "region", "country"):
        allr[f"eta2_{col}"] = eta2(df.alpha_k_med, df[col])
        allr[f"n_levels_{col}"] = int(df[col].nunique())
    allr["spearman_lat"] = float(spearmanr(df.lat, df.alpha_k_med).statistic)
    allr["spearman_alt"] = float(spearmanr(df.alt, df.alpha_k_med).statistic)
    allr["spearman_cbh"] = float(spearmanr(df.cbh_med, df.alpha_k_med).statistic)
    rows.append(allr)
    ft = pd.DataFrame(rows)
    ft.to_csv(OUT / "factor_table.csv", index=False, float_format="%.3f")
    say(ft.round(3).to_string(index=False))
    say("(eta^2 with many levels is inflated by chance: compare with n_levels; country has the most levels)")

    # region effect controlled for type: residual from the type median
    df["resid_type"] = df.alpha_k_med - df.groupby("itype").alpha_k_med.transform("median")
    say("\n== region effect after removing the type median (km-1) ==")
    say(df.groupby("region").resid_type.agg(["size", "median", "min", "max"]).round(2).to_string())
    say(f"eta^2(region | type-residual) = {eta2(df.resid_type, df.region):.3f};  "
        f"eta^2(firmware | type-residual) = {eta2(df.resid_type, df.fw_main):.3f};  "
        f"eta^2(gate | type-residual) = {eta2(df.resid_type, df.gate_m):.3f}")

    # confounding: region vs instrument configuration (does the low-alpha group follow the SITE or the setup?)
    say("\n== confounding checks (unit-median alpha_k, type residual) ==")
    c31 = df[df.itype == "CL31"]
    say("CL31 region x gate (median residual, n):\n" + c31.pivot_table(
        index="region", columns="gate_m", values="resid_type", aggfunc=lambda x: f"{np.median(x):+.2f} ({len(x)})").to_string())
    for fw in ("205", "321", "174", "300"):
        g = c31[c31.fw_main == fw]
        if len(g) >= 5:
            say(f"CL31 firmware {fw} only ({len(g)} units): region medians " + ", ".join(
                f"{r} {v:+.2f} (n={n})" for r, v, n in zip(*[g.groupby('region').resid_type.median().index,
                                                               g.groupby('region').resid_type.median().values,
                                                               g.groupby('region').size().values])))
    g = c31[c31.gate_m == 10.0]
    say(f"CL31 10 m gates only: eta2(region) {eta2(g.resid_type, g.region):.3f}, eta2(firmware) {eta2(g.resid_type, g.fw_main):.3f}")
    g = c31[c31.region.isin(['central Europe', 'western Europe'])]
    say(f"CL31 central+western Europe only ({len(g)} units): eta2(firmware) {eta2(g.resid_type, g.fw_main):.3f}, "
        f"eta2(gate) {eta2(g.resid_type, g.gate_m):.3f}; firmware medians " +
        ", ".join(f"{f} {v:+.2f} (n={n})" for f, v, n in zip(g.groupby('fw_main').resid_type.median().index,
                                                             g.groupby('fw_main').resid_type.median().values,
                                                             g.groupby('fw_main').size().values)))
    low = df[df.region.isin(["subarctic maritime", "northern Europe"])]
    say(f"low-alpha group (subarctic maritime + northern Europe): {len(low)} units, types {low.itype.value_counts().to_dict()}, "
        f"firmware {low.fw_main.value_counts().to_dict()}, gates {low.gate_m.value_counts().to_dict()}")
    if months:
        mm = pd.concat(months, ignore_index=True)
        mm = mm[mm.n >= 200]
        mm["mon"] = mm.month.astype(str).str[-2:].astype(int)
        mm["anom"] = mm.alpha_k_med / mm.groupby("key").alpha_k_med.transform("median") - 1
        say("\n== seasonality: median monthly anomaly of alpha_k (unit month / unit median - 1), per type ==")
        say(mm.pivot_table(index="mon", columns="itype", values="anom", aggfunc="median").mul(100).round(1).to_string())
    figure(df, say)
    log.close()


def figure(df, say):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {"CL31": "C0", "CL51": "C2", "CL61": "C3"}
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    ax = axes[0]
    sc = ax.scatter(df.lon, df.lat, c=df.alpha_k_med, cmap="viridis", s=35, vmin=18, vmax=36,
                    edgecolor=[cols.get(t, "k") for t in df.itype], linewidth=1.2)
    fig.colorbar(sc, ax=ax, label="unit median alpha_k [km$^{-1}$]")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("alpha_k per unit (edge colour = type)")
    ax.grid(alpha=0.3)
    ax = axes[1]
    for it, g in df.groupby("itype"):
        ax.plot(g.alpha_k_med, g.lat, "o", color=cols.get(it, "k"), label=f"{it} (n={len(g)})", alpha=0.8)
    ax.set_xlabel("unit median alpha_k [km$^{-1}$]")
    ax.set_ylabel("latitude")
    ax.set_title("alpha_k vs latitude")
    ax.grid(alpha=0.3)
    ax.legend()
    ax = axes[2]
    order = df.groupby("fw_main").alpha_k_med.median().sort_values().index.tolist()
    data = [df.loc[df.fw_main == f, "alpha_k_med"].values for f in order]
    ax.boxplot(data, vert=False, labels=[f"{f} (n={len(d)})" for f, d in zip(order, data)])
    ax.set_xlabel("unit median alpha_k [km$^{-1}$]")
    ax.set_title("alpha_k by firmware version")
    ax.grid(alpha=0.3)
    fig.suptitle("Network alpha_k (tail slope): climate or instrument?", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_network_alpha_k.png", dpi=130)
    plt.close(fig)
    say("figure written")


if __name__ == "__main__":
    main()
