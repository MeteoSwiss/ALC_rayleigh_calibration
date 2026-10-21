# -*- coding: utf-8 -*-
"""Simplified, slot-sized data figures for the CHM15k temperature deck, from the report-19 result files.

Run with the project venv:  .venv/Scripts/python.exe presentation/chm15k_temperature/make_data_figures.py
"""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
SKILL = Path("C:/Users/hervo/.claude/plugins/cache/presentation/meteoswiss-presentation/0.1.0/"
             "skills/meteoswiss-presentation/scripts")
sys.path.insert(0, str(SKILL))
sys.path.insert(0, str(REPO / "validation/chm15k_stability"))
from illustration import C, new_figure, save  # noqa: E402
from common import OUT, load_cache, load_constants  # noqa: E402
from q1_temperature import ANOM_WINDOW, night_hk, robust_mask, running_median  # noqa: E402

FIG = Path(__file__).resolve().parent / "fig"


def rows(name):
    return list(csv.DictReader(open(OUT / name, encoding="utf-8")))


def season_vs_anomaly(wmo="0-20000-0-10393", ident="0", label="Lindenberg CHM15k"):
    """Slide 3: raw C vs T (season mixed in) and the warm/cold-spell anomaly view, one unit."""
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    T = night_hk(cache, cons["t0"], cons["t1"])[:, list(cache["hk_vars"]).index("temp_int")]
    lc = np.log(cons["c"])
    m = np.isfinite(T) & np.isfinite(lc)
    m[m] = robust_mask(lc[m])
    t, lc, T = cons["t"][m], lc[m], T[m]
    o = np.argsort(t)
    t, lc, T = t[o], lc[o], T[o]
    la = lc - running_median(t, lc, ANOM_WINDOW)
    ta = T - running_median(t, T, ANOM_WINDOW)
    doy = (t - 0.0) % 365.25
    fig, axs = new_figure("figure", ncols=2)
    ax = axs[0]
    sc = ax.scatter(T, 100 * (lc - np.median(lc)), c=doy, cmap="twilight", s=12, vmin=0, vmax=365)
    ax.set_xlabel("internal temperature [°C]")
    ax.set_ylabel("C$_L$ − median [%]")
    ax.set_title("raw (colour = season)", fontsize=13)
    ok = np.isfinite(la) & np.isfinite(ta)
    ax = axs[1]
    ax.scatter(ta[ok], 100 * la[ok], s=12, color=C["blue_dark"])
    p = np.polyfit(ta[ok], 100 * la[ok], 1)
    xx = np.linspace(np.min(ta[ok]), np.max(ta[ok]), 10)
    ax.plot(xx, np.polyval(p, xx), "-", color=C["red"], lw=2.5, label=f"{p[0]:+.1f} %/K")
    ax.set_xlabel("temperature anomaly vs 31-day median [K]")
    ax.set_ylabel("C$_L$ anomaly [%]")
    ax.set_title("warm / cold spells only", fontsize=13)
    ax.legend(frameon=False, loc="upper left")
    save(fig, FIG / "fig_season_vs_anomaly.png")


def pairs():
    """Slide 5: diurnal cycle of the signal ratio of three same-wavelength pairs."""
    z = np.load(OUT / "q2b_frames.npz")
    band = "1.5-2.5 km"
    picks = [("AMSTERDAM_AP_SCHIP B/C", "Schiphol B / C", C["red"]),
             ("BUDAPEST A/B", "Budapest A / B", C["brown"]),
             ("MAGURELE A/B", "Măgurele A / B", C["blue_dark"])]
    fig, ax = new_figure("figure")
    ax.axvspan(0, 4, color=C["grey_light"], lw=0)
    ax.axvspan(22, 24, color=C["grey_light"], lw=0)
    for key, lab, col in picks:
        sh = z[f"{key}|sh"]
        rel = z[f"{key}|rel_{band}"]
        hb = np.floor(sh).astype(int)
        med = [np.nanmedian(rel[hb == h]) if np.isfinite(rel[hb == h]).sum() > 10 else np.nan for h in range(24)]
        ax.plot(np.arange(24) + 0.5, 100 * np.array(med), "-", color=col, lw=2.8, label=lab)
    ax.axhline(0, color=C["grey"], lw=0.8)
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 3))
    ax.set_xlabel("local solar hour")
    ax.set_ylabel("signal ratio vs night [%]")
    ax.legend(frameon=False, loc="upper left")
    save(fig, FIG / "fig_pairs.png")


def prediction():
    """Slide 6: daytime anomaly predicted by the night-time coefficient vs observed."""
    r = rows("q3_synthesis.csv")
    fig, ax = new_figure("figure")
    lim = (-20, 34)
    ax.plot(lim, lim, "--", color=C["grey"], lw=1)
    for kind, mk, col, lab in (("pair", "o", C["blue_dark"], "co-located CHM15k pairs"),
                               ("CHM/CL61", "s", C["brown"], "CHM15k vs CL61")):
        sel = [x for x in r if x["kind"] == kind and "Payerne" not in x["label"]]
        ax.scatter([float(x["pred"]) for x in sel], [float(x["obs"]) for x in sel], s=70, marker=mk, color=col,
                   label=lab, zorder=3)
    pay = [x for x in r if "Payerne" in x["label"]][0]
    ax.scatter(float(pay["pred"]), float(pay["obs"]), s=70, marker="s", facecolor="white", edgecolor=C["red"],
               lw=2, zorder=3, label="Payerne (abnormal receiver)")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("predicted from the night-time coefficient [%]")
    ax.set_ylabel("observed day − night [%]")
    ax.annotate("Payerne", (float(pay["pred"]), float(pay["obs"])), xytext=(8, -4), textcoords="offset points",
                fontsize=11, color=C["red"])
    ax.legend(frameon=False, loc="upper left", fontsize=11)
    save(fig, FIG / "fig_prediction.png")


def excess():
    """Slide 7: distribution over the network of the noon excess of the effective constant."""
    r = rows("q3_daytime_excess.csv")
    fig, ax = new_figure("figure")
    bins = np.arange(-6, 32, 2)
    for key, col, lab in (("excess", C["blue_dark"], "all seasons"), ("excess_JJA", C["red"], "summer")):
        x = np.array([float(v[key]) for v in r])
        x = x[np.isfinite(x)]
        ax.hist(np.clip(x, -6, 30), bins=bins, histtype="step", lw=2.5, color=col,
                label=f"{lab}: median {np.median(x):+.1f} %, > 10 %: {100 * np.mean(x > 10):.0f} % of units")
    ax.axvline(0, color=C["grey"], lw=0.8)
    ax.set_xlabel("noon effective constant vs night calibration [%]")
    ax.set_ylabel("number of CHM15k")
    ax.legend(frameon=False, loc="upper right", fontsize=11)
    save(fig, FIG / "fig_excess.png")


def validation():
    """Slide 8: day − night bias against independent references, before / after correction."""
    r = rows("q9_validation_tcorr.csv")
    comps = list(dict.fromkeys(x["comparison"] for x in r))
    names = []
    for c in comps:
        site, rest = c.split(":", 1)
        site = site.replace("Bucharest-Măgurele", "Măgurele").replace("Leipzig 0", "Leipzig")
        ref = "CL61" if "CL61" in rest else "EARLINET"
        names.append(f"{site}\nvs {ref}")
    fig, ax = new_figure("figure")
    x = np.arange(len(comps))
    for var, col, dx, lab in (("operational", C["grey"], -0.2, "operational (night constant)"),
                              ("T-corrected", C["red"], 0.2, "temperature-corrected")):
        v = [float(next(y["day_minus_night"] for y in r if y["comparison"] == c and y["variant"] == var)) for c in comps]
        ax.bar(x + dx, v, width=0.38, color=col, label=lab)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=11)
    ax.set_ylabel("day − night bias [%-pt]")
    ax.legend(frameon=False, loc="upper right", fontsize=11)
    save(fig, FIG / "fig_validation.png")


def generation():
    """Generation slide: kappa, night background and regulated set-points vs module production year."""
    import matplotlib.pyplot as plt
    r = rows("q14_generation.csv")
    yr = np.array([int(x["year"]) for x in r])
    years = [y for y in range(2011, 2025) if np.sum(yr == y) >= 2]
    fig0, _ = new_figure("figure")
    plt.close(fig0)
    fig, axs = plt.subplots(1, 3, figsize=(8.8, 2.95))
    jit = np.random.default_rng(0).uniform(-0.25, 0.25, len(yr))

    def med_line(ax, x):
        ax.plot(years, [np.nanmedian(x[yr == y]) for y in years], "o-", color=C["red"], lw=2, ms=4)
    k = np.array([float(x["kappa_in"]) for x in r])
    ax = axs[0]
    ax.plot(yr + jit, k, ".", color=C["grey"], alpha=0.5, ms=4)
    med_line(ax, k)
    ax.set_ylim(-1, 4.5)
    ax.set_title("κ [%/K]", fontsize=12)
    bg = np.array([float(x["background"]) for x in r])
    ax = axs[1]
    ax.plot(yr + jit, bg, ".", color=C["grey"], alpha=0.5, ms=4)
    med_line(ax, bg)
    ax.set_yscale("log")
    ax.set_ylim(5e-4, 5e-2)
    ax.set_title("night background (dark)", fontsize=12)
    td = np.array([float(x["T_det"]) for x in r])
    tl = np.array([float(x["T_lom"]) for x in r])
    ax = axs[2]
    share = [100 * np.mean(td[yr == y] > 27) for y in years]
    ax.bar(years, share, color=C["blue_light"], edgecolor=C["blue_dark"], label="detector at 29.85 °C [% units]")
    ax.set_ylim(0, 100)
    ax.set_title("regulated set-points", fontsize=12)
    ax2 = ax.twinx()
    ax2.plot(years, [np.nanmedian(tl[yr == y]) for y in years], "s-", color=C["brown"], lw=2, ms=4,
             label="optical module [°C]")
    ax2.set_ylim(33.5, 35.2)
    ax2.tick_params(labelsize=10)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.5, loc="upper left")
    for a in axs:
        a.set_xticks([2012, 2016, 2020, 2024])
        a.tick_params(labelsize=10)
        a.axvline(2015.5, color=C["grey"], ls=":", lw=1)
        a.axvline(2019.5, color=C["grey"], ls=":", lw=1)
    axs[1].set_xlabel("optical-module production year (TUB serial)", fontsize=11)
    save(fig, FIG / "fig_generation.png")


if __name__ == "__main__":
    generation()
    season_vs_anomaly()
    pairs()
    prediction()
    excess()
    validation()
    print("ok")
