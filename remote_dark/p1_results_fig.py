# -*- coding: utf-8 -*-
"""Fig 23 — P1 archive sweep results on one page.

Left: the network timeline of all flagged full-day candidates, classified by the pattern that
identifies their nature without opening a single file: RUNS of >=3 consecutive days on one
stream (cover or outage — hardware doesn't heal in a day), days SYNCHRONISED across >=5
streams (weather — snow on windows), and isolated singles/pairs (mostly alpine snow).
Payerne's hood-session days are the positive control (green stars).

Right: the two case verifications — Berus (the best true-cover candidate, with its
housekeeping verdict in the title) and the Payerne control day, each against a normal day.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from remote_dark.common import ensure_out, read_day, PAYERNE

FIG = ensure_out("figs")
SCR = Path(r"C:/Users/hervo/AppData/Local/Temp/claude"
           r"/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration"
           r"/265f36ab-1a40-4592-bae6-6c3550f2fac1/scratchpad")
CAND = Path(__file__).parent / "p1_cover_candidates.csv"


def med_profile(rng, rcs):
    z = np.where(rng > 0, rng, np.nan)
    return np.nanmedian(rcs / z[None, :] ** 2, axis=0)


def day_from_nc(path):
    with netCDF4.Dataset(path) as ds:
        rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
        rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f4"), np.nan)
        hk = {}
        for v in ("status_detector", "status_laser", "window_transmission"):
            if v in ds.variables:
                hk[v] = float(np.nanmedian(
                    np.ma.filled(ds.variables[v][:].astype("f8"), np.nan)))
    return rng, rcs, hk


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    rows = list(csv.DictReader(open(CAND)))
    by_stream = defaultdict(list)
    by_date = defaultdict(set)
    for r in rows:
        key = (r["wmo"], r["ident"], r["site"])
        by_stream[key].append(r["date"])
        by_date[r["date"]].add(key)

    # classify each candidate day
    kinds = {}
    for key, ds in by_stream.items():
        days = sorted(datetime.strptime(d, "%Y%m%d") for d in ds)
        in_run = set()
        i = 0
        while i < len(days):
            j = i
            while j + 1 < len(days) and (days[j + 1] - days[j]).days <= 1:
                j += 1
            if j - i + 1 >= 3:
                in_run.update(days[i:j + 1])
            i = j + 1
        for d in days:
            ds8 = f"{d:%Y%m%d}"
            if key[2] == "PAYERNE":
                kinds[(key, ds8)] = "control"
            elif d in in_run:
                kinds[(key, ds8)] = "run"
            elif len(by_date[ds8]) >= 5:
                kinds[(key, ds8)] = "sync"
            else:
                kinds[(key, ds8)] = "single"

    order = {"run": 0, "control": 1, "single": 2, "sync": 3}
    streams = sorted(by_stream, key=lambda k: (
        min(order[kinds[(k, d)]] for d in by_stream[k]), k[2]))
    colors = {"run": "crimson", "sync": "steelblue", "single": "0.55",
              "control": "seagreen"}

    fig = plt.figure(figsize=(16.5, 9.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.9, 1.0], hspace=0.28, wspace=0.16)
    ax = fig.add_subplot(gs[:, 0])
    for iy, key in enumerate(streams):
        for d in by_stream[key]:
            k = kinds[(key, d)]
            t = mdates.date2num(datetime.strptime(d, "%Y%m%d"))
            ax.plot(t, iy, "*" if k == "control" else "s",
                    color=colors[k], ms=9 if k == "control" else 4,
                    mec="none" if k != "control" else "k")
    ax.set_yticks(range(len(streams)))
    ax.set_yticklabels([k[2][:20] for k in streams], fontsize=5.5)
    ax.set_ylim(-1, len(streams))
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(alpha=0.2, axis="x")
    for lab, colr in (("run ≥3 d (cover/outage)", "crimson"),
                      ("synchronised ≥5 stations (weather)", "steelblue"),
                      ("isolated (mostly snow)", "0.55"),
                      ("Payerne hood control", "seagreen")):
        ax.plot([], [], "s" if "hood" not in lab else "*", color=colr, label=lab,
                ms=5 if "hood" not in lab else 9)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(f"226 flagged full-day candidates, {len(streams)} streams — "
                 "the temporal pattern is the classifier", fontsize=10)

    # --- Berus verification ------------------------------------------------------------------
    axb = fig.add_subplot(gs[0, 1])
    rng_c, rcs_c, hk_c = day_from_nc(SCR / "L1_0-20000-0-10704_020260207.nc")
    rng_n, rcs_n, _ = day_from_nc(SCR / "L1_0-20000-0-10704_020260220.nc")
    for rng, rcs, colr, lab in ((rng_n, rcs_n, "0.6", "normal day 2026-02-20"),
                                (rng_c, rcs_c, "crimson", "candidate 2026-02-07")):
        med = med_profile(rng, rcs)
        m = (rng >= 100) & (rng <= 15000)
        axb.plot(med[m], rng[m] / 1e3, "-", color=colr, lw=1.5, label=lab)
    axb.axvline(0, color="k", lw=0.5)
    med_c = med_profile(rng_c, rcs_c)
    span = np.nanpercentile(np.abs(med_c[(rng_c >= 300) & (rng_c <= 10000)]), 98)
    axb.set_xlim(-8 * span, 8 * span)
    hks = "  ".join(f"{k.split('_')[-1]}={v:.0f}" for k, v in hk_c.items())
    axb.set_title(f"BERUS candidate — HK: {hks}", fontsize=9)
    axb.set_xlabel("P-view (rcs/z²)")
    axb.set_ylabel("altitude (km)")
    axb.legend(fontsize=7)
    axb.grid(alpha=0.25)
    print("BERUS 2026-02-07 HK:", hk_c)

    # --- Payerne control day -----------------------------------------------------------------
    axp = fig.add_subplot(gs[1, 1])
    for date, colr, lab in ((datetime(2026, 5, 20), "0.6", "normal day 2026-05-20"),
                            (datetime(2026, 5, 26), "seagreen", "hood day 2026-05-26")):
        d = read_day(PAYERNE["wmo"], "A", date)
        med = med_profile(d["rng"], d["rcs"])
        m = (d["rng"] >= 100) & (d["rng"] <= 15000)
        axp.plot(med[m], d["rng"][m] / 1e3, "-", color=colr, lw=1.5, label=lab)
        if "hood" in lab:
            span = np.nanpercentile(np.abs(med[(d["rng"] >= 300) & (d["rng"] <= 10000)]), 98)
            axp.set_xlim(-8 * span, 8 * span)
    axp.axvline(0, color="k", lw=0.5)
    axp.set_title("PAYERNE positive control (known hood session)", fontsize=9)
    axp.set_xlabel("P-view (rcs/z²)")
    axp.set_ylabel("altitude (km)")
    axp.legend(fontsize=7)
    axp.grid(alpha=0.25)

    fig.suptitle("P1 archive sweep — every detected case (left) and the case verifications "
                 "(right): profile collapse + housekeeping decide cover vs fault vs weather",
                 fontsize=11)
    f = FIG / "fig23_p1_results.png"
    fig.savefig(f, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
