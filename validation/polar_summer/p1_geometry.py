# -*- coding: utf-8 -*-
"""p1 - polar summer: which stations lose their Rayleigh night, when, and how much sun is left.

1. Solar geometry (operational formula of calibration.io.data_loader.filter_time_range: declination
   of the night key date, hour angle from the station longitude, no equation of time), for every
   E-PROFILE L1 census stream and every night key date of 2025:
     * the LOWEST sun of the night = SZA at solar midnight = 180 - |lat + decl|;
     * the dark duration (hours with SZA > threshold) for thresholds 100, 98, 96, 94, 92, 90;
     * the number of nights per year without a single qualifying profile, per threshold;
     * the low-sun window the v2.3 fallback would use (calibration.rayleigh.polar, defaults).
2. Measured solar background: CHM15k bckgrd_rcs_0 (10-min HK cache of report 19) versus SZA at the
   polar stations, relative to each station's dark median (SZA >= 105, clear blocks).

CPU only (no L1 read). Outputs in C:/DATA/Projects/202610_v23/polar/:
  p1_stations.csv, p1_stations.md, p1_geometry.png, p1_background.csv, p1_background.png
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))

from calibration.rayleigh import polar  # noqa: E402

OUT = Path("C:/DATA/Projects/202610_v23/polar")
CENSUS = REPO / "validation/scope_l1_2026_census.json"
THRESHOLDS = (100.0, 98.0, 96.0, 94.0, 92.0, 90.0)
RAYLEIGH_TYPES = {"CHM15k", "CL61", "Mini-MPL"}       # scripts/run_network_calibration.RAYLEIGH_TYPES
DEFAULT = dict(delta=0.0, min_hours=6.0, min_sza=78.0)    # CalibrationOptions defaults (v2.3)
EPOCH = datetime(1970, 1, 1)


def nights_2025():
    d0 = datetime(2025, 1, 1)
    return [(d0 + timedelta(days=k)).strftime("%Y%m%d") for k in range(365)]


def station_table():
    census = json.load(open(CENSUS, encoding="utf-8"))
    dates = nights_2025()
    rows = []
    curves = {}
    for s in census:
        lat = float(s["lat"])
        geo = [polar.night_geometry(d, lat, 100.0) for d in dates]
        szamax = np.array([g.sza_max for g in geo])
        row = dict(key=f"{s['wmo']}_{s['ident']}", site=s["site"], type=s["type"], lat=lat,
                   lon=float(s["lon"]), rayleigh=s["type"] in RAYLEIGH_TYPES,
                   sza_max_min=float(szamax.min()))
        for thr in THRESHOLDS:
            none = szamax <= thr
            row[f"n_none_{int(thr)}"] = int(none.sum())
            if none.any():
                idx = np.where(none)[0]
                row[f"first_{int(thr)}"] = dates[idx[0]]
                row[f"last_{int(thr)}"] = dates[idx[-1]]
        dark = np.array([polar.hours_above_sza(lat, g.decl_rad, 100.0) for g in geo])
        row["n_dark_lt1h"] = int((dark < 1.0).sum())
        row["n_dark_lt3h"] = int((dark < 3.0).sum())
        # fallback window at the summer solstice of the station's hemisphere
        sol = "20250621" if lat >= 0 else "20251221"
        g = polar.night_geometry(sol, lat, 100.0)
        fl = polar.low_sun_floor(g, lat, DEFAULT["delta"], DEFAULT["min_hours"], DEFAULT["min_sza"])
        row.update(sza_max_solstice=g.sza_max, fallback_floor_solstice=fl,
                   fallback_window_h_solstice=polar.hours_above_sza(lat, g.decl_rad, fl))
        rows.append(row)
        if row["n_none_100"] > 0 and s["type"] in RAYLEIGH_TYPES:
            curves[row["key"]] = (s["site"], s["type"], lat, szamax)
    return rows, curves


def write_tables(rows):
    import pandas as pd
    df = pd.DataFrame(rows).sort_values("lat", ascending=False)
    df.to_csv(OUT / "p1_stations.csv", index=False)
    aff = df[df.n_none_100 > 0]
    lines = ["| Station | Type | Lat | Lowest sun at solstice (SZA) | Nights/yr without SZA>100 | "
             "period | >98 | >96 | >94 | >92 | >90 | nights with < 1 h / < 3 h dark | fallback window at solstice |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in aff.iterrows():
        per = f"{r.first_100[4:6]}-{r.first_100[6:]} .. {r.last_100[4:6]}-{r.last_100[6:]}"
        rl = "" if r.rayleigh else " (no Rayleigh)"
        lines.append(
            f"| {r.site.title().replace('_', ' ')} | {r.type}{rl} | {r.lat:.1f} | {r.sza_max_solstice:.1f} | "
            f"{r.n_none_100} | {per} | {r.n_none_98} | {r.n_none_96} | {r.n_none_94} | {r.n_none_92} | "
            f"{r.n_none_90} | {r.n_dark_lt1h} / {r.n_dark_lt3h} | SZA > {r.fallback_floor_solstice:.1f}, "
            f"{r.fallback_window_h_solstice:.1f} h |")
    (OUT / "p1_stations.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    south = df[df.lat < 0]
    print(f"{len(df)} streams; {len(aff)} with nights without SZA>100 "
          f"({int(aff.rayleigh.sum())} Rayleigh-calibrated types); southern: "
          + ", ".join(f"{r.site} {r.lat:.1f} lowest-sun {r.sza_max_min:.1f}" for _, r in south.iterrows()))
    print(aff[["site", "type", "lat", "sza_max_solstice", "n_none_100", "n_none_98", "n_none_96",
               "n_none_94", "n_none_92", "n_none_90", "n_dark_lt1h", "n_dark_lt3h",
               "fallback_floor_solstice", "fallback_window_h_solstice"]].to_string())
    return df


def background_vs_sza():
    """CHM15k bckgrd_rcs_0 (HK cache, clear 10-min blocks) / dark median versus SZA."""
    import pandas as pd
    from common import load_cache
    census = {(s["wmo"], s["ident"]): s for s in json.load(open(CENSUS, encoding="utf-8"))}
    units = [k for k in census if census[k]["type"] == "CHM15k" and census[k]["lat"] >= 56.6]
    edges = np.array([60, 70, 75, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 105, 180.0])
    rows = []
    for wmo, ident in units:
        s = census[(wmo, ident)]
        c = load_cache(wmo, ident)
        if c is None:
            continue
        names = [str(x) for x in c["hk_vars"]]
        bg = c["hk"][:, :, names.index("bckgrd_rcs_0")].astype(float).ravel()
        clear = ((c["n_clear"] > 0) & (c["n_clear"] >= 0.8 * c["n_all"])).ravel()
        t = (np.asarray(c["day"], float)[:, None] + (np.arange(144)[None, :] + 0.5) / 144.0).ravel()
        tsol = t + s["lon"] / 360.0
        night = np.floor(tsol + 0.5)
        sh = np.mod((tsol - np.floor(tsol)) * 24.0, 24.0)
        doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in night])
        decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
        la = np.radians(s["lat"])
        cz = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(np.radians((sh - 12) * 15))
        sza = np.degrees(np.arccos(np.clip(cz, -1, 1)))
        ok = np.isfinite(bg) & (bg > 0) & clear
        dark = ok & (sza >= 105)
        if dark.sum() < 100:
            continue
        ref = float(np.median(bg[dark]))
        for a, b in zip(edges[:-1], edges[1:]):
            m = ok & (sza >= a) & (sza < b)
            if m.sum() >= 20:
                r = bg[m] / ref
                rows.append(dict(key=f"{wmo}_{ident}", site=s["site"], lat=s["lat"], sza_lo=a,
                                 sza_hi=b, ratio_median=float(np.median(r)),
                                 ratio_p25=float(np.percentile(r, 25)),
                                 ratio_p75=float(np.percentile(r, 75)),
                                 solar_abs=float(np.median(bg[m]) - ref), bg_dark=ref, n=int(m.sum())))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "p1_background.csv", index=False)
    return df


def figures(df, curves, bgdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig, ax = plt.subplots(1, 2, figsize=(17, 6.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    doy = np.arange(1, 366)
    lats = [v[2] for v in curves.values()]
    norm = matplotlib.colors.Normalize(vmin=min(lats) - 1, vmax=max(lats) + 1)
    for key, (site, typ, lat, szamax) in sorted(curves.items(), key=lambda kv: -kv[1][2]):
        a.plot(doy, szamax, color=cm.viridis_r(norm(lat)), lw=1.6,
               ls="-" if typ == "CHM15k" else "--",
               label=f"{site.title().replace('_', ' ')} {lat:.1f} N ({typ})")
    for thr in THRESHOLDS:
        a.axhline(thr, color="tab:red" if thr == 100 else "0.6", lw=1.2 if thr == 100 else 0.7,
                  ls="-" if thr == 100 else ":")
    a.text(3, 100.4, "operational threshold SZA 100", color="tab:red", fontsize=8)
    a.axhline(DEFAULT["min_sza"], color="k", lw=0.8, ls="--")
    a.text(250, DEFAULT["min_sza"] + 0.4, f"fallback floor SZA {DEFAULT['min_sza']:.0f}", fontsize=8)
    a.set_xlim(1, 365)
    a.set_ylim(76, 112)
    a.set_xlabel("day of year (night key date, 2025)")
    a.set_ylabel("lowest sun of the night: SZA at solar midnight [deg]")
    a.set_title("Lowest sun reached each night (operational SZA formula)")
    a.legend(fontsize=7.5, ncol=1, loc="lower left")
    a.grid(alpha=0.3)

    a = ax[1]
    lat_grid = np.arange(50.0, 80.01, 0.25)
    dates = nights_2025()
    for thr, col in zip(THRESHOLDS, ("tab:red", "tab:orange", "gold", "tab:green", "tab:blue",
                                     "tab:purple")):
        n = []
        for la in lat_grid:
            sm = np.array([polar.night_geometry(d, la, 100.0).sza_max for d in dates[::1]])
            n.append(int((sm <= thr).sum()))
        a.plot(lat_grid, n, color=col, lw=1.8, label=f"no profile with SZA > {thr:.0f}")
    aff = df[(df.n_none_100 > 0) & df.rayleigh].sort_values("lat")
    a.plot(aff.lat, aff.n_none_100, "o", color="tab:red", ms=6, mec="k", zorder=5)
    # hand-placed label offsets (points) for the clustered stations
    offs = {"OSLO": (-30, 6), "LERWICK": (6, -12), "FLESLAND": (6, 2), "BIRKENES": (6, -11),
            "STORNOWAY": (6, 3), "ANDOYA": (-38, 4), "JAN_MAYEN": (6, -10)}
    for _, r in aff.iterrows():
        a.annotate(r.site.title().replace("_", " "), (r.lat, r.n_none_100), fontsize=7,
                   xytext=offs.get(r.site, (6, -10)), textcoords="offset points")
    a.set_xlabel("latitude [deg N] (same in the south, six months later)")
    a.set_ylabel("nights per year without a qualifying profile")
    a.set_title("Nights lost to the SZA threshold")
    a.legend(fontsize=8, loc="upper left")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "p1_geometry.png", dpi=130)
    plt.close(fig)

    if bgdf is not None and len(bgdf):
        fig, ax = plt.subplots(1, 2, figsize=(16, 5.6))
        a = ax[0]
        for key, g in bgdf.groupby("key"):
            lat = g.lat.iloc[0]
            mid = 0.5 * (g.sza_lo + np.minimum(g.sza_hi, 110))
            a.plot(mid, g.ratio_median, "o-", ms=3, color=cm.viridis_r(norm(lat)),
                   label=f"{g.site.iloc[0].title().replace('_', ' ')} {lat:.1f} N")
        a.axvline(100, color="tab:red", lw=1)
        a.axvline(DEFAULT["min_sza"], color="k", ls="--", lw=0.8)
        a.axhline(1, color="0.5", lw=0.8)
        a.set_xlim(62, 110)
        a.set_xlabel("solar zenith angle [deg]")
        a.set_ylabel("background / dark-night background (SZA >= 105)")
        a.set_title("CHM15k solar background (bckgrd_rcs_0), clear 10-min blocks")
        a.legend(fontsize=7.5, ncol=2)
        a.grid(alpha=0.3)
        a = ax[1]
        for key, g in bgdf.groupby("key"):
            lat = g.lat.iloc[0]
            mid = 0.5 * (g.sza_lo + np.minimum(g.sza_hi, 110))
            a.semilogy(mid, np.maximum(g.solar_abs, 1e-6), "o-", ms=3, color=cm.viridis_r(norm(lat)))
        a.axvline(100, color="tab:red", lw=1)
        a.axvline(DEFAULT["min_sza"], color="k", ls="--", lw=0.8)
        a.axhline(1e-6, color="0.6", lw=0.8, ls=":")
        a.text(63, 1.25e-6, "<= 1e-6: not measurable", fontsize=7, color="0.4")
        a.set_xlim(62, 110)
        a.set_xlabel("solar zenith angle [deg]")
        a.set_ylabel("solar part of the background [photons/shot]")
        a.set_title("absolute solar background (background - dark median)")
        a.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(OUT / "p1_background.png", dpi=130)
        plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, curves = station_table()
    df = write_tables(rows)
    bg = background_vs_sza()
    print(bg.groupby(["sza_lo"])["ratio_median"].median().to_string())
    figures(df, curves, bg)


if __name__ == "__main__":
    main()
