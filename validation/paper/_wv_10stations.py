"""10-station horizontal-resolution impact of the coarse CAMS on the WV cloud-calibration correction.

Resolution ladder, per station (5 flat + 5 mountain/valley), winter (Jan) + summer (Jul) 2025:
  - CAMS 1 deg    (~111 km, legacy)          D:/CAMS
  - CAMS 0.4 deg  (~44 km, OPERATIONAL)      D:/CAMS_Monthly_04
  - ERA5 0.25 deg (~28 km, independent ref)  OUT/era5_{season}.nc
Diagnostics: model orography error (cell surface - station), operational WV correction T2_wv(2 km),
and the calibration error dC(%) = 100*(T2_source(cbh)/T2_ERA5(cbh) - 1). ERA5 is the finer independent
reference (for deep valleys even ERA5 smooths the terrain -> see the CERRA add-on).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.paper.wv_resolution_lib import (  # noqa: E402
    STATIONS, OUT, REPO_FIGS, cams_month_file, cams1_month_file, nwv_cams_fast, nwv_era5,
    cams_surface_alt, t2_of, t2_at, iwv_mm, dC_pct)

SEASONS = {"winter": "202501", "summer": "202507"}
DAYS = range(10, 21)
CBH = 2000.0
ERA5 = {s: OUT / f"era5_{s}.nc" for s in SEASONS}


def collect():
    rows = []
    for st in STATIONS:
        for season, ym in SEASONS.items():
            f04 = cams_month_file(ym + "15")
            f1 = cams1_month_file(ym + "15")
            era5f = ERA5[season] if ERA5[season].is_file() and ERA5[season].stat().st_size > 1000 else None
            if not f04.is_file():
                continue
            oro04 = cams_surface_alt(f04, st["lat"], st["lon"]) - st["alt"]
            oro1 = (cams_surface_alt(f1, st["lat"], st["lon"]) - st["alt"]) if f1.is_file() else np.nan
            for d in DAYS:
                for hr in (0, 12):
                    when = np.datetime64(f"{ym[:4]}-{ym[4:6]}-{d:02d}T{hr:02d}:00:00")
                    c04 = nwv_cams_fast(f04, st["lat"], st["lon"], when)
                    t2_04 = t2_at(*t2_of(*c04, st["alt"], st["inst"]), st["alt"], CBH)
                    rec = dict(name=st["name"], inst=st["inst"], terrain=st["terrain"], alt=st["alt"],
                               season=season, oro04=oro04, oro1=oro1, t2_04=t2_04, iwv=iwv_mm(*c04))
                    if f1.is_file():
                        c1 = nwv_cams_fast(f1, st["lat"], st["lon"], when)
                        rec["t2_1d"] = t2_at(*t2_of(*c1, st["alt"], st["inst"]), st["alt"], CBH)
                    if era5f is not None:
                        er = nwv_era5(era5f, st["lat"], st["lon"], when)
                        t2_e = t2_at(*t2_of(*er, st["alt"], st["inst"]), st["alt"], CBH)
                        rec["t2_era5"] = t2_e
                        rec["dC_04"] = dC_pct(t2_04, t2_e)                       # operational vs ERA5
                        rec["dC_1d"] = dC_pct(rec.get("t2_1d", np.nan), t2_e)    # 1 deg vs ERA5
                    if "t2_1d" in rec:
                        rec["dC_1d_04"] = dC_pct(rec["t2_1d"], t2_04)            # 1 deg vs 0.4 deg
                    rows.append(rec)
    return pd.DataFrame(rows)


def summary(df):
    print(f"\n10-STATION WV horizontal-resolution impact (cloud base {CBH/1000:.0f} km; dC = calib. error %)")
    print(f"  {'station':12s} {'terr':4s} {'alt':>5s} | {'oro err 1deg/0.4':>16s} | {'WVcorr%':>7s} | "
          f"{'dC 0.4-ERA5 w/s':>16s} | {'dC 1deg-ERA5 w/s':>16s}")
    for _, s0 in df.drop_duplicates("name").iterrows():
        sub = df[df.name == s0["name"]]
        wv = 100 * (1 - sub["t2_04"].median())

        def med(col, season):
            if col not in sub:
                return np.nan
            x = sub[sub.season == season][col].dropna()
            return np.median(x) if len(x) else np.nan
        print(f"  {s0['name']:12s} {s0['terrain'][:4]:4s} {s0['alt']:5.0f} | "
              f"{s0['oro1']:+7.0f}/{s0['oro04']:+7.0f} | {wv:6.1f} | "
              f"{med('dC_04','winter'):+6.1f}/{med('dC_04','summer'):+6.1f}  | "
              f"{med('dC_1d','winter'):+6.1f}/{med('dC_1d','summer'):+6.1f}")


def figs(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    order = df.drop_duplicates("name").sort_values("oro04")
    names = order["name"].tolist()
    col = {"flat": "C0", "mountain": "C3"}
    has_e = "dC_04" in df.columns and df["dC_04"].notna().any()

    def med(n, season, c):
        if c not in df:
            return np.nan
        x = df[(df.name == n) & (df.season == season)][c].dropna()
        return np.median(x) if len(x) else np.nan

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))
    y = np.arange(len(names))
    # (a) orography error 1deg vs 0.4deg
    ax = axes[0]
    o1 = [order[order.name == n]["oro1"].iloc[0] for n in names]
    o4 = [order[order.name == n]["oro04"].iloc[0] for n in names]
    ax.barh(y + 0.2, o1, height=0.38, color="0.6", label="CAMS 1$\\degree$")
    ax.barh(y - 0.2, o4, height=0.38, color="C2", label="CAMS 0.4$\\degree$ (oper.)")
    ax.set_yticks(y); ax.set_yticklabels([f"{n} ({order[order.name==n]['alt'].iloc[0]:.0f} m)" for n in names], fontsize=8.5)
    ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("model orography error [m]  (cell $-$ station)")
    ax.set_title("(a) grid-cell altitude mismatch"); ax.grid(alpha=0.3, axis="x"); ax.legend(fontsize=8, loc="lower right")

    # (b) dC vs ERA5 (summer) for 1deg and 0.4deg
    ax = axes[1]
    if has_e:
        d04 = [med(n, "summer", "dC_04") for n in names]
        d1 = [med(n, "summer", "dC_1d") for n in names]
        ax.barh(y - 0.2, d04, height=0.38, color="C2", label="CAMS 0.4$\\degree$ $-$ ERA5")
        ax.barh(y + 0.2, d1, height=0.38, color="0.6", label="CAMS 1$\\degree$ $-$ ERA5")
        ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8.5)
        ax.axvline(0, color="k", lw=0.8); ax.set_xlabel("$\\Delta C$ at 2 km [%]  (summer, vs ERA5)")
        ax.set_title("(b) calibration error vs independent ERA5"); ax.legend(fontsize=8); ax.grid(alpha=0.3, axis="x")
    else:
        ax.text(0.5, 0.5, "ERA5 pending", ha="center", transform=ax.transAxes)

    # (c) |dC 0.4-ERA5| vs |orography error|
    ax = axes[2]
    if has_e:
        for terr in ("flat", "mountain"):
            sub = df[df.terrain == terr]
            for nm, g in sub.groupby("name"):
                xx = abs(g["oro04"].iloc[0]); yy = np.nanmedian(np.abs(g["dC_04"]))
                ax.scatter(xx, yy, s=70, color=col[terr], alpha=0.85, marker=("o" if terr == "flat" else "^"))
                ax.annotate(nm, (xx, yy), fontsize=7, xytext=(4, 2), textcoords="offset points")
        ax.set_xlabel("|CAMS 0.4$\\degree$ orography error| [m]"); ax.set_ylabel("|$\\Delta C$| CAMS 0.4$\\degree$$-$ERA5 at 2 km [%]")
        ax.set_title("(c) impact grows with orography error"); ax.grid(alpha=0.3)
        ax.legend(handles=[Patch(color="C0", label="flat"), Patch(color="C3", label="mountain/valley")], fontsize=9)
    fig.suptitle("Coarse-CAMS horizontal-resolution impact on the water-vapour correction — 10 stations (2025)",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for dd in (OUT, REPO_FIGS):
        dd.mkdir(parents=True, exist_ok=True); fig.savefig(dd / "fig_wv_10stations_impact.png", dpi=150)
    plt.close(fig)
    print("  fig_wv_10stations_impact.png")

    # ---- profile case study: flat (Payerne) vs valley (Aosta), summer ----
    when = np.datetime64("2025-07-15T12:00:00")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, nm in zip(axes, ("Payerne", "Aosta")):
        st = next(s for s in STATIONS if s["name"] == nm)
        f04 = cams_month_file("20250715"); f1 = cams1_month_file("20250715")
        srcs = [("CAMS 1$\\degree$", nwv_cams_fast(f1, st["lat"], st["lon"], when), "0.5", ":"),
                ("CAMS 0.4$\\degree$ (oper.)", nwv_cams_fast(f04, st["lat"], st["lon"], when), "C3", "--")]
        if ERA5["summer"].is_file():
            srcs.append(("ERA5 0.25$\\degree$", nwv_era5(ERA5["summer"], st["lat"], st["lon"], when), "C0", "-"))
        for lab, pr, c, ls in srcs:
            zg, t2 = t2_of(*pr, st["alt"], st["inst"])
            ax.plot(t2, (zg - st["alt"]) / 1000.0, ls, color=c, lw=1.9, label=lab)
        oro = cams_surface_alt(f04, st["lat"], st["lon"]) - st["alt"]
        ax.set_ylim(0, 5); ax.set_xlabel("$T^2_{WV}$ (two-way)"); ax.set_ylabel("height AGL [km]")
        ax.set_title(f"{nm} ({st['alt']:.0f} m, {st['inst']}) — CAMS 0.4$\\degree$ oro err {oro:+.0f} m")
        ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower left")
    fig.suptitle("WV transmission by model resolution — flat plateau vs deep Alpine valley (2025-07-15 12 UT)",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for dd in (OUT, REPO_FIGS):
        fig.savefig(dd / "fig_wv_10stations_profiles.png", dpi=150)
    plt.close(fig)
    print("  fig_wv_10stations_profiles.png")


def main():
    df = collect()
    df.to_csv(OUT / "wv_10stations.csv", index=False)
    summary(df)
    figs(df)
    print("WV_10STATIONS_DONE")


if __name__ == "__main__":
    main()
