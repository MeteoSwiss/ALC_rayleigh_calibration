"""Payerne anchor for the CAMS-resolution / WV study (build on _cl61_wv_sources_probe.py).

Compares the water-vapour two-way transmission T2_wv(z) and its impact on the cloud calibration
(C ~ T2_wv(cbh)) from three sources with DECREASING resolution error:
  - radiosonde   (in-situ vertical truth, 00/12 UT)                        <- truth
  - IFS @ point  (Cloudnet ..._ecmwf.nc, native ~9 km, HOURLY, 137 lev)   <- isolates horizontal+temporal
  - CAMS 0.4 deg (operational, 3-hourly, 101 lev, smoothed orography)     <- operational baseline

Winter (Jan 2025) vs summer (Jul 2025). Outputs a summary table + three landscape figures.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validation.paper.wv_resolution_lib import (  # noqa: E402
    PAYERNE, OUT, REPO_FIGS, cams_month_file, ifs_file_payerne, nwv_cams_fast, nwv_ifs,
    nwv_ifs_all_hours, nwv_sounding, cams_surface_alt, cams_grid, t2_of, t2_at, iwv_mm, dC_pct)

CBH = [500.0, 1000.0, 1500.0, 2000.0, 3000.0]          # cloud base heights [m AGL]
SEASONS = {"winter": "202501", "summer": "202507"}
P = PAYERNE


def collect():
    """Loop days x {00,12 UT}; record T2@cbh for CAMS/IFS/sonde + IWV, per season."""
    rows = []
    for season, ym in SEASONS.items():
        camsf = cams_month_file(ym + "15")
        ndays = pd.Period(f"{ym[:4]}-{ym[4:6]}").days_in_month
        for d in range(1, ndays + 1):
            day8 = f"{ym}{d:02d}"
            ifsf = ifs_file_payerne(day8)
            if not ifsf.is_file() or not camsf.is_file():
                continue
            for hr in (0, 12):
                when = np.datetime64(f"{ym[:4]}-{ym[4:6]}-{d:02d}T{hr:02d}:00:00")
                snd = nwv_sounding(day8, hr, P["alt"])
                if snd is None:
                    continue
                cam = nwv_cams_fast(camsf, P["lat"], P["lon"], when)
                ifs = nwv_ifs(ifsf, when)
                rec = dict(season=season, day=day8, hr=hr,
                           iwv_snd=iwv_mm(*snd), iwv_ifs=iwv_mm(*ifs), iwv_cam=iwv_mm(*cam))
                zg_s, t2_s = t2_of(*snd, P["alt"], P["inst"])
                zg_i, t2_i = t2_of(*ifs, P["alt"], P["inst"])
                zg_c, t2_c = t2_of(*cam, P["alt"], P["inst"])
                for cb in CBH:
                    s = t2_at(zg_s, t2_s, P["alt"], cb)
                    i = t2_at(zg_i, t2_i, P["alt"], cb)
                    c = t2_at(zg_c, t2_c, P["alt"], cb)
                    rec[f"t2s_{cb:.0f}"] = s
                    rec[f"t2i_{cb:.0f}"] = i
                    rec[f"t2c_{cb:.0f}"] = c
                    rec[f"dC_cam_{cb:.0f}"] = dC_pct(c, s)   # CAMS vs sonde truth
                    rec[f"dC_ifs_{cb:.0f}"] = dC_pct(i, s)   # IFS  vs sonde truth
                    rec[f"dC_ci_{cb:.0f}"] = dC_pct(c, i)    # CAMS vs IFS (pure horizontal)
                rows.append(rec)
    return pd.DataFrame(rows)


def summary(df):
    off, sp = cams_grid(cams_month_file("20250115"), P["lat"], P["lon"])
    csurf = cams_surface_alt(cams_month_file("20250115"), P["lat"], P["lon"])
    print(f"\nPAYERNE  station alt={P['alt']:.0f} m  |  CAMS 0.4deg cell alt={csurf:.0f} m "
          f"(orography error {csurf-P['alt']:+.0f} m; nearest-point offset {off:.2f} deg)")
    print(f"samples: winter n={ (df.season=='winter').sum() }, summer n={ (df.season=='summer').sum() }")
    print("\n  T2_wv error -> cloud-calibration error dC (%), median [P10,P90], truth = radiosonde")
    print(f"  {'season':7s} {'cbh':>5s} | {'IWV(mm)':>8s} | {'dC CAMS-sonde':>22s} | {'dC IFS-sonde':>20s} | {'dC CAMS-IFS':>18s}")
    for season in SEASONS:
        g = df[df.season == season]
        iwv = g["iwv_snd"].median()
        for cb in (1000.0, 2000.0, 3000.0):
            def q(col):
                x = g[col].to_numpy(); x = x[np.isfinite(x)]
                return np.median(x), np.percentile(x, 10), np.percentile(x, 90)
            mc, lc, hc = q(f"dC_cam_{cb:.0f}")
            mi, li, hi = q(f"dC_ifs_{cb:.0f}")
            mci, lci, hci = q(f"dC_ci_{cb:.0f}")
            print(f"  {season:7s} {cb/1000:.1f}km | {iwv:8.1f} | "
                  f"{mc:+5.1f} [{lc:+4.1f},{hc:+4.1f}]   | {mi:+5.1f} [{li:+4.1f},{hi:+4.1f}] | "
                  f"{mci:+5.1f} [{lci:+4.1f},{hci:+4.1f}]")


def profile_fig(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5))
    days = {"winter": "20250115", "summer": "20250715"}
    for r, (season, day8) in enumerate(days.items()):
        when = np.datetime64(f"{day8[:4]}-{day8[4:6]}-{day8[6:8]}T12:00:00")
        snd = nwv_sounding(day8, 12, P["alt"]); ifs = nwv_ifs(ifs_file_payerne(day8), when)
        cam = nwv_cams_fast(cams_month_file(day8), P["lat"], P["lon"], when)
        srcs = [("radiosonde (truth)", snd, "k", "-"), ("IFS @ point (hourly)", ifs, "C0", "-"),
                ("CAMS 0.4$\\degree$ (3-hourly)", cam, "C3", "--")]
        # (col0) n_wv profile
        ax = axes[r, 0]
        for lab, pr, c, ls in srcs:
            ax.plot(np.asarray(pr[1]) / 1e22, (np.asarray(pr[0]) - P["alt"]) / 1000.0, ls, color=c, lw=1.7, label=lab)
        ax.set_ylim(0, 5); ax.set_xlim(left=0)
        ax.set_xlabel("$n_{H_2O}$ [$10^{22}$ m$^{-3}$]"); ax.set_ylabel("height AGL [km]")
        ax.set_title(f"({'ad'[r]}) {season} 2025-{day8[4:6]}-{day8[6:8]} 12 UT  —  humidity")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
        # (col1) T2_wv profile
        ax = axes[r, 1]
        for lab, pr, c, ls in srcs:
            zg, t2 = t2_of(*pr, P["alt"], P["inst"])
            ax.plot(t2, (zg - P["alt"]) / 1000.0, ls, color=c, lw=1.7, label=lab)
        ax.set_ylim(0, 5); ax.set_xlabel("$T^2_{WV}$ (two-way)"); ax.set_ylabel("height AGL [km]")
        ax.set_title(f"({'be'[r]}) {season}  —  WV transmission  (IWV={iwv_mm(*snd):.0f} mm)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
        # (col2) dC vs cloud height
        ax = axes[r, 2]
        zg_s, t2_s = t2_of(*snd, P["alt"], P["inst"])
        zg_i, t2_i = t2_of(*ifs, P["alt"], P["inst"])
        zg_c, t2_c = t2_of(*cam, P["alt"], P["inst"])
        hh = np.linspace(0.3, 4.0, 40)
        dcc = [dC_pct(t2_at(zg_c, t2_c, P["alt"], h*1000), t2_at(zg_s, t2_s, P["alt"], h*1000)) for h in hh]
        dci = [dC_pct(t2_at(zg_i, t2_i, P["alt"], h*1000), t2_at(zg_s, t2_s, P["alt"], h*1000)) for h in hh]
        ax.plot(dcc, hh, "C3--", lw=1.8, label="CAMS $-$ sonde")
        ax.plot(dci, hh, "C0-", lw=1.8, label="IFS $-$ sonde")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_ylim(0.3, 4.0); ax.set_xlabel("cloud-calibration error $\\Delta C$ [%]"); ax.set_ylabel("cloud base height [km]")
        ax.set_title(f"({'cf'[r]}) {season}  —  impact on C"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Payerne — water-vapour correction vs profile source and season (radiosonde = truth)",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for d in (OUT, REPO_FIGS):
        d.mkdir(parents=True, exist_ok=True); fig.savefig(d / "fig_wv_payerne_profiles.png", dpi=155)
    plt.close(fig)
    print("  fig_wv_payerne_profiles.png")


def seasonal_fig(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    # (a) dC distributions at 2 km, CAMS-sonde and IFS-sonde, winter vs summer
    ax = axes[0]
    data, labels, colors = [], [], []
    for season in ("winter", "summer"):
        g = df[df.season == season]
        data += [g["dC_cam_2000"].dropna(), g["dC_ifs_2000"].dropna()]
        labels += [f"{season}\nCAMS", f"{season}\nIFS"]
        colors += ["C3", "C0"]
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False, widths=0.6)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("$\\Delta C$ at 2 km cloud base [%]")
    ax.set_title("(a) calibration error by source and season"); ax.grid(alpha=0.3, axis="y")
    # (b) T2@2km vs IWV: the correction magnitude grows with humidity
    ax = axes[1]
    for season, mk in (("winter", "o"), ("summer", "s")):
        g = df[df.season == season]
        ax.scatter(g["iwv_snd"], 100 * (1 - g["t2s_2000"]), s=18, marker=mk, alpha=0.6, label=f"{season} (sonde)")
    ax.set_xlabel("integrated water vapour IWV [mm]"); ax.set_ylabel("WV correction at 2 km  100·(1$-T^2$) [%]")
    ax.set_title("(b) correction magnitude vs humidity"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    # (c) |dC| CAMS-sonde vs IWV: error grows with humidity
    ax = axes[2]
    for season, mk, c in (("winter", "o", "C0"), ("summer", "s", "C3")):
        g = df[df.season == season]
        ax.scatter(g["iwv_snd"], g["dC_cam_2000"], s=18, marker=mk, color=c, alpha=0.6, label=f"{season}")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("integrated water vapour IWV [mm]"); ax.set_ylabel("$\\Delta C$ CAMS$-$sonde at 2 km [%]")
    ax.set_title("(c) CAMS error grows with humidity"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle("Payerne — seasonal water-vapour correction error (CAMS 0.4$\\degree$ operational vs radiosonde truth)",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for d in (OUT, REPO_FIGS):
        fig.savefig(d / "fig_wv_payerne_seasonal.png", dpi=155)
    plt.close(fig)
    print("  fig_wv_payerne_seasonal.png")


def temporal_fig():
    """1h vs 3h: from hourly IFS, the error of 3-hourly sampling (linear interp) on T2@2km -> dC."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    # (a) representative summer day diurnal T2@2km: 1h vs 3h-subsampled
    day8 = "20250715"
    times, profs = nwv_ifs_all_hours(ifs_file_payerne(day8))
    hrs = (times - times.astype("datetime64[D]")[0]).astype("timedelta64[m]").astype(float) / 60.0
    t2_1h = np.array([t2_at(*t2_of(*pr, P["alt"], P["inst"]), P["alt"], 2000.0) for pr in profs])
    sub = np.arange(0, len(hrs), 3)                          # 0,3,6,... UT (3-hourly, like CAMS)
    t2_3h = np.interp(hrs, hrs[sub], t2_1h[sub])
    ax = axes[0]
    ax.plot(hrs, t2_1h, "C0-o", ms=3, lw=1.4, label="hourly (IFS truth)")
    ax.plot(hrs, t2_3h, "C3--", lw=1.6, label="3-hourly sampled + interp")
    ax.plot(hrs[sub], t2_1h[sub], "C3s", ms=7, mfc="none", label="3-hourly nodes")
    ax.set_xlabel("hour (UT)"); ax.set_ylabel("$T^2_{WV}$ at 2 km"); ax.set_xlim(0, 24)
    ax.set_title(f"(a) diurnal WV transmission — summer {day8[4:6]}-{day8[6:8]}"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    # (b),(c) distribution of the 3h-sampling error over each season
    for k, (season, ym) in enumerate(SEASONS.items()):
        errs = []
        ndays = pd.Period(f"{ym[:4]}-{ym[4:6]}").days_in_month
        for d in range(1, ndays + 1):
            f = ifs_file_payerne(f"{ym}{d:02d}")
            if not f.is_file():
                continue
            tt, pr = nwv_ifs_all_hours(f)
            if tt.size < 24:
                continue
            hh = (tt - tt.astype("datetime64[D]")[0]).astype("timedelta64[m]").astype(float) / 60.0
            t1 = np.array([t2_at(*t2_of(*p, P["alt"], P["inst"]), P["alt"], 2000.0) for p in pr])
            ss = np.arange(0, len(hh), 3)
            t3 = np.interp(hh, hh[ss], t1[ss])
            errs.append(100 * (t3 / t1 - 1))                # dC (%) from 3h sampling
        errs = np.concatenate(errs) if errs else np.array([0.0])
        ax = axes[1 + k]
        ax.hist(errs, bins=np.linspace(-3, 3, 60), color=("C0" if season == "winter" else "C3"), alpha=0.7)
        rms = np.sqrt(np.nanmean(errs**2)); p95 = np.nanpercentile(np.abs(errs), 95)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("$\\Delta C$ from 3-hourly sampling [%]"); ax.set_ylabel("count")
        ax.set_title(f"({'bc'[k]}) {season}: RMS={rms:.2f}%, P95={p95:.2f}%"); ax.grid(alpha=0.3)
    fig.suptitle("Payerne — benefit of 1-hourly vs 3-hourly WV sampling (hourly IFS as truth, error at 2 km cloud base)",
                 fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for d in (OUT, REPO_FIGS):
        fig.savefig(d / "fig_wv_payerne_temporal.png", dpi=155)
    plt.close(fig)
    print("  fig_wv_payerne_temporal.png")


def main():
    OUT.mkdir(parents=True, exist_ok=True); REPO_FIGS.mkdir(parents=True, exist_ok=True)
    df = collect()
    df.to_csv(OUT / "wv_payerne_sources.csv", index=False)
    summary(df)
    profile_fig(df)
    seasonal_fig(df)
    temporal_fig()
    print("PAYERNE_WV_DONE")


if __name__ == "__main__":
    main()
