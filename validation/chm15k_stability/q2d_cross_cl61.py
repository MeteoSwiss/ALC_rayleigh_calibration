# -*- coding: utf-8 -*-
"""Q2d of report 19: absolute day/night test of the CHM15k against a co-located CL61.

Same-wavelength CHM15k pairs only see DIFFERENTIAL effects: a sensitivity loss common to both units
(e.g. the daytime drop of the internal calibration pulse with solar background, seen at every unit)
cancels in their ratio. The CL61 uses an analog APD receiver without photon counting, so it gives
an independent reference by day. The ratio CHM15k/CL61 is affected by the 1064/910 nm wavelength
difference (aerosol Angstrom exponent, water vapour at 910 nm), which can change between the
nocturnal residual layer and the daytime mixed layer by a few percent; the test is therefore
sensitive to effects of ~10 % and more, which is the size of the calibration-pulse drop.

Per site: diurnal composite of ln(rcs_CHM / rcs_CL61) relative to the night, and its regression on
the CHM15k's own calibration-pulse anomaly and on its solar background. The firmware divides beta_raw
by p_calc, so: slope 0 = the daytime p_calc drop is a real detector sensitivity loss, correctly
compensated; slope -1 = the drop is an artefact of the reference and the daytime signal is
over-corrected (inflated).

Needs the CL61 streams in the cache:  python validation/chm15k_stability/q2d_cross_cl61.py --build
Run:                                   python validation/chm15k_stability/q2d_cross_cl61.py
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, OUT, chm15k_streams, load_cache  # noqa: E402
from q2b_pairs import BANDS, MAIN_BAND, boot_median_ci, hourly, is_night, ols_slope, pair_frame, relative  # noqa: E402

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")

# (wmo, CHM15k ident, CL61 ident, label)
SITES = (("0-20000-0-06610", "A", "C", "Payerne"),
         ("0-20000-0-10393", "0", "C", "Lindenberg"),
         ("0-20000-0-03808", "A", "C", "Camborne"),
         ("0-380-5-1", "0", "B", "Aosta"))


def build():
    from concurrent.futures import ProcessPoolExecutor
    from hk_cache import _worker
    jobs = [(w, cl, False) for w, _, cl, _ in SITES]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for key, status, nd in ex.map(_worker, jobs):
            print(f"  {key}: {nd} d {status}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    a = ap.parse_args()
    if a.build:
        build()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lon_of = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    rows, frames = [], []
    for wmo, chm, cl, label in SITES:
        cc, cv = load_cache(wmo, chm), load_cache(wmo, cl)
        if cc is None or cv is None:
            print(f"  {label}: missing cache ({chm}: {cc is not None}, {cl}: {cv is not None})")
            continue
        fr = pair_frame(hourly(cc), hourly(cv), lon_of[(wmo, chm)])
        months = np.array([(datetime(1970, 1, 1) + timedelta(days=float(t))).month for t in fr["t"]])
        with np.errstate(invalid="ignore", divide="ignore"):
            fr["lnp"] = np.log(fr["calibration_pulse_A"])
            fr["lnbg"] = np.log(fr["bckgrd_rcs_0_A"])
        rel_p = relative(fr, "lnp")
        day = (fr["sh"] >= 10) & (fr["sh"] < 16)
        row = dict(site=label, n_days=len(np.unique(np.floor(fr["t"]))))
        for band in BANDS:
            rel = relative(fr, band)
            fr["rel_" + band] = rel
            m, lo, hi = boot_median_ci(rel[day])
            row[f"dmn_{band}"] = 100 * m
            row[f"dmn_ci_{band}"] = f"[{100 * lo:+.1f},{100 * hi:+.1f}]"
            for s, mm in (("JJA", (6, 7, 8)), ("DJF", (12, 1, 2))):
                sel = day & np.isin(months, mm)
                row[f"dmn_{s}_{band}"] = 100 * np.nanmedian(rel[sel]) if np.isfinite(rel[sel]).sum() > 20 else np.nan
        rel = fr["rel_" + MAIN_BAND]
        b, se, n = ols_slope(rel_p, rel)
        row.update(b_lnp=b, se_lnp=se, n_lnp=n, sd_relp_day=100 * float(np.nanstd(rel_p[day])),
                   calpulse_day=100 * float(np.nanmedian(rel_p[day])))
        fr["rel_p"] = rel_p
        rows.append(row)
        frames.append((label, fr))
        print(f"  {label:10s} days {row['n_days']:4d}  day-night CHM/CL61 " +
              "  ".join(f"{bd} {row[f'dmn_{bd}']:+5.1f}% {row[f'dmn_ci_{bd}']}" for bd in BANDS) +
              f"  JJA {row[f'dmn_JJA_{MAIN_BAND}']:+.1f} DJF {row[f'dmn_DJF_{MAIN_BAND}']:+.1f}  "
              f"CHM cal-pulse day {row['calpulse_day']:+.1f} %  dlnR/dln p_cal {b:+.2f}+/-{se:.2f}", flush=True)
    if not rows:
        return
    with open(OUT / "q2d_cross_cl61.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    fig, axs = plt.subplots(2, len(frames), figsize=(5.2 * len(frames), 8), squeeze=False)
    for j, (label, fr) in enumerate(frames):
        ax = axs[0, j]
        hb = np.floor(fr["sh"]).astype(int)
        for band, col in zip(BANDS, ("#9467bd", "#1f77b4", "#17becf")):
            rel = fr["rel_" + band]
            med = [np.nanmedian(rel[hb == h]) if np.isfinite(rel[hb == h]).sum() > 10 else np.nan for h in range(24)]
            ax.plot(np.arange(24) + 0.5, 100 * np.array(med), "-", color=col, lw=1.6, label=f"CHM/CL61 {band}")
        rp = fr["rel_p"]
        med = [np.nanmedian(rp[hb == h]) if np.isfinite(rp[hb == h]).sum() > 10 else np.nan for h in range(24)]
        ax.plot(np.arange(24) + 0.5, 100 * np.array(med), "k--", lw=1.4, label="CHM15k calibration pulse")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(0, 24)
        ax.set_xticks(range(0, 25, 3))
        ax.set_ylim(-45, 25)
        ax.set_xlabel("local solar hour")
        ax.set_ylabel("relative to night [%]")
        ax.set_title(f"{label}: CHM15k / CL61 signal ratio", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
        ax = axs[1, j]
        rel = fr["rel_" + MAIN_BAND]
        ok = np.isfinite(rel) & np.isfinite(rp)
        ax.plot(100 * rp[ok], 100 * rel[ok], ".", ms=1.5, alpha=0.2)
        if ok.sum() > 50:
            edges = np.percentile(rp[ok], np.linspace(1, 99, 12))
            mid = [100 * np.median(rp[ok & (rp >= a) & (rp < b)]) for a, b in zip(edges[:-1], edges[1:])]
            med = [100 * np.median(rel[ok & (rp >= a) & (rp < b)]) for a, b in zip(edges[:-1], edges[1:])]
            p = np.polyfit(rp[ok], rel[ok], 1)
            ax.plot(mid, med, "o-", color="#d62728", lw=2, label=f"slope {p[0]:+.2f}")
        ax.plot([-50, 20], [50, -20], "k--", lw=0.8, label="slope -1: p_cal artefact, over-corrected")
        ax.set_xlim(-50, 20)
        ax.set_ylim(-50, 30)
        ax.set_xlabel("CHM15k ln p$_{cal}$ - night ref [%]")
        ax.set_ylabel(f"ln(CHM/CL61) {MAIN_BAND} - night ref [%]")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q2d_cross_cl61.png", dpi=130)
    print(f"-> {FIGS / 'q2d_cross_cl61.png'}")
    np.savez_compressed(OUT / "q2d_frames.npz", **{f"{lab}|{k}": v for lab, fr in frames for k, v in fr.items()})


if __name__ == "__main__":
    main()
