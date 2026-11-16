# -*- coding: utf-8 -*-
"""Paper Kalman figure (Kleine Scheidegg CHM15k), operational vs temperature-corrected nightly constants.

Top panel = the paper figure (fig_v09_kalman_kleine_scheidegg). Bottom panel = the same operational
improve_alc Kalman filter run on the nightly constants normalised to a reference internal temperature:

    C_norm = C_night * exp(-kappa_u (T_cal - T_ref))

with T_cal the median temp_int over each night's calibration window (10-min HK cache) and kappa_u the
unit's intra-night coefficient (report 19, q2a; the conservative, deployable estimate).

Run with the project venv:  .venv/Scripts/python.exe presentation/chm15k_temperature/make_kalman_tcorr.py
"""
import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation/chm15k_stability"))
from validation.paper.kalman_improve_alc import fitfunc_sin, kalman_improve_alc  # noqa: E402
from common import OUT, load_cache, load_constants  # noqa: E402
from q1_temperature import night_hk  # noqa: E402

CD = "C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/calib_dashboard"
# station from the command line: make_kalman_tcorr.py <wmo> <ident> <name> <tag>  (default Kleine Scheidegg)
WMO, IDENT, NAME, TAG = (sys.argv[1:5] if len(sys.argv) >= 5 else ("0-20000-0-06736", "A", "Kleine Scheidegg", "kse"))
FIG = Path(__file__).resolve().parent / "fig"


def night_temperature():
    """Calibration-window median temp_int per calibration date."""
    cons, cache = load_constants(WMO, IDENT), load_cache(WMO, IDENT)
    T = night_hk(cache, cons["t0"], cons["t1"])[:, list(cache["hk_vars"]).index("temp_int")]
    dates = pd.to_datetime(cons["t"], unit="D").normalize()
    return pd.Series(T, index=dates).groupby(level=0).first()


def kappa():
    for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8")):
        if r["stream"] == f"{WMO}_{IDENT}":
            return float(r["b_T_ctrl"])
    raise KeyError(WMO)


def panel(ax, t, v, title, show_legend):
    r = kalman_improve_alc(list(t), v, "CHM15k")
    out = np.array([pd.Timestamp(x).date() in set(r["outlier_dates"]) for x in t])
    g = pd.to_datetime(r["grid"])
    st, sd = r["state"], r["std"]
    ax.fill_between(g, (st - sd) / 1e11, (st + sd) / 1e11, color="tab:green", alpha=0.25, lw=0,
                    label="best estimate ±1σ")
    ax.plot(g, st / 1e11, color="tab:green", lw=1.6, label="Kalman best estimate")
    amp = np.nan
    if r["model"] == "annual cycle":
        ax.plot(g, fitfunc_sin(np.array([x.timetuple().tm_yday for x in g]), *r["params"]) / 1e11, ":", color="k",
                lw=1.2, label="annual-cycle prediction model")
        amp = 100 * abs(r["params"][0]) / r["params"][2]
    ax.plot(t[~out], v[~out] / 1e11, "o", color="tab:red", ms=3.5, label=f"nightly constant ({int((~out).sum())} nights)")
    ax.plot(t[out], v[out] / 1e11, "x", color="0.45", ms=6, mew=1.2, label=f"outlier, IQR criterion ({int(out.sum())})")
    rms = 100 * np.sqrt(r["var_meas"]) / np.median(v[~out])
    ax.set_title(f"{title}   —   annual-cycle amplitude ±{amp:.0f} %, scatter around the model {rms:.0f} %, "
                 f"Kalman 1σ {100 * np.median(sd / st):.1f} %", fontsize=10, loc="left")
    ax.set_ylabel("C$_L$ (10$^{11}$ a.u.)")
    ax.grid(alpha=0.3)
    if show_legend:
        ax.legend(loc="upper left", ncol=3, fontsize=8.5, framealpha=0.9)
    return amp, rms


def main():
    d = pd.read_csv(f"{CD}/{WMO}_{IDENT}_rayleigh.csv", parse_dates=["time"])
    n = d.C_daily.notna()
    t = np.array(list(d.time[n]))
    v = d.C_daily[n].values
    Tn = night_temperature().reindex(pd.to_datetime(d.time[n]).dt.normalize()).values
    k = kappa()
    t_ref = float(np.nanmedian(Tn))
    ok = np.isfinite(Tn)
    vc = v * np.exp(-k / 100.0 * (Tn - t_ref))
    fig, axs = plt.subplots(2, 1, figsize=(16, 7.6), dpi=200, sharex=True, sharey=True)
    # both panels on the same nights (the HK cache, hence T_cal, ends in July 2026)
    a0, r0 = panel(axs[0], t[ok], v[ok], "Operational nightly constants (paper figure, same nights as below)", True)
    a1, r1 = panel(axs[1], t[ok], vc[ok],
                   f"Normalised to {t_ref:.0f} °C internal temperature, κ = {k:+.2f} %/K", False)
    axs[1].set_xlim(pd.Timestamp("2024-12-01"), pd.Timestamp("2026-07-31"))
    axs[1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    last = pd.Timestamp(t[ok][-1]).strftime("%d %B %Y")
    fig.suptitle(f"CHM15k {NAME} (WIGOS {WMO}, instrument {IDENT}), 1 January 2025 – {last}",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / f"fig_kalman_tcorr_{TAG}.png")
    print(f"kappa {k:+.2f} %/K, T_ref {t_ref:.1f} degC, nights {ok.sum()}/{len(v)}")
    print(f"annual amplitude {a0:.1f} -> {a1:.1f} %, scatter around model {r0:.1f} -> {r1:.1f} %")


def slide_version():
    """Slot-sized single panel for the deck: both Kalman best estimates and annual-cycle models overlaid."""
    sys.path.insert(0, "C:/Users/hervo/.claude/plugins/cache/presentation/meteoswiss-presentation/0.1.0/"
                       "skills/meteoswiss-presentation/scripts")
    from illustration import C, new_figure, save
    d = pd.read_csv(f"{CD}/{WMO}_{IDENT}_rayleigh.csv", parse_dates=["time"])
    n = d.C_daily.notna()
    t = np.array(list(d.time[n]))
    v = d.C_daily[n].values
    Tn = night_temperature().reindex(pd.to_datetime(d.time[n]).dt.normalize()).values
    ok = np.isfinite(Tn)
    k = kappa()
    t_ref = float(np.nanmedian(Tn))
    t, v, Tn = t[ok], v[ok], Tn[ok]
    vc = v * np.exp(-k / 100.0 * (Tn - t_ref))
    fig, ax = new_figure("figure")
    for vals, col, lab, z in ((v, C["grey"], "operational", 2), (vc, C["red"], "temperature-corrected", 3)):
        r = kalman_improve_alc(list(t), vals, "CHM15k")
        g = pd.to_datetime(r["grid"])
        amp = 100 * abs(r["params"][0]) / r["params"][2]
        out = np.array([pd.Timestamp(x).date() in set(r["outlier_dates"]) for x in t])
        ax.plot(t[~out], vals[~out] / 1e11, "o", color=col, ms=2.2, alpha=0.35, zorder=z)
        ax.plot(g, r["state"] / 1e11, "-", color=col, lw=2.0, zorder=z + 1, label=f"{lab}: seasonal ±{amp:.0f} %")
        ax.plot(g, fitfunc_sin(np.array([x.timetuple().tm_yday for x in g]), *r["params"]) / 1e11, ":", color=col,
                lw=1.5, zorder=z + 1)
    ax.set_ylabel("C$_L$ (10$^{11}$ a.u.)")
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.legend(frameon=False, loc="upper left", fontsize=11)
    save(fig, FIG / f"fig_kalman_slide_{TAG}.png")


if __name__ == "__main__":
    main()
    slide_version()
