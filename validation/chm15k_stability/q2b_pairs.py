# -*- coding: utf-8 -*-
"""Q2b of report 19: diurnal / short-term stability from co-located same-wavelength CHM15k pairs.

Two CHM15k at the same site look at (nearly) the same air at the same wavelength, so in the ratio
of their clear-sky signals the atmosphere cancels - aerosol included, by day as well as by night.
What is left is the ratio of their lidar constants (plus noise and horizontal inhomogeneity for units
a few km apart). If the constant of either unit moves with time of day, temperature, solar
background or loss of thermal regulation, the ratio moves.

For each pair and hourly clear block:
  lnR(z-band) = ln( <rcs_A> / <rcs_B> ) in bands 1.0-1.5, 1.5-2.5, 2.5-4.0 km
  lnR_rel     = lnR - (median of night-time lnR within +-15 days)   (removes the slow level/swaps)
Analyses: diurnal composite of lnR_rel by local solar hour and season; regression on the
instrument-temperature difference (night and day separately), on the solar background, and on the
optical-module / detector excursion above set-point (regulation loss, Rosu et al. 2026).

Run:  python validation/chm15k_stability/q2b_pairs.py
"""
from __future__ import annotations

import csv
import itertools
import sys
from datetime import datetime, timedelta
from pathlib import Path

import warnings

import numpy as np

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, FIGS, OUT, PAIR_SITES, RANGE_MID, chm15k_streams, load_cache,  # noqa: E402
                    solar_hour)

BANDS = {"1.0-1.5 km": (1000, 1500), "1.5-2.5 km": (1500, 2500), "2.5-4.0 km": (2500, 4000)}
MAIN_BAND = "1.5-2.5 km"
MIN_CLEAR_FRAC = 0.9
NIGHT_HOURS = (22.0, 4.0)       # local solar time used as the night reference
REF_HALF_WINDOW = 15.0          # days; only has to follow slow drifts (cloudy sites have few clear nights)
SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}


def hourly(cache):
    """Aggregate the 10-min cache to hourly means over clear blocks. Returns dict keyed by hour index."""
    nb = cache["rcs"].shape[-1]
    nd = len(cache["day"])
    rcs = cache["rcs"].reshape(nd, 24, 6, nb)
    clear = (cache["n_clear"] >= MIN_CLEAR_FRAC * np.maximum(cache["n_all"], 1)) & (cache["n_all"] > 0)
    clear = clear.reshape(nd, 24, 6)
    hk = cache["hk"].reshape(nd, 24, 6, -1)
    rcs = np.where(clear[..., None], rcs, np.nan)
    with np.errstate(invalid="ignore"):
        rcs_h = np.nanmean(rcs, axis=2)                       # (nd, 24, nb)
        hk_h = np.nanmedian(hk, axis=2)                       # (nd, 24, nvar)
    nclear = clear.sum(axis=2)
    rcs_h[nclear < 5] = np.nan                                # hour must be >= 50 of 60 min clear
    # Integer hour index (hours since 1970) of each hour START; exact, no rounding ambiguity.
    t = (np.round(cache["day"]).astype(np.int64)[:, None] * 24 + np.arange(24)[None, :]).ravel()
    return dict(t=t, rcs=rcs_h.reshape(-1, nb), hk=hk_h.reshape(-1, hk.shape[-1]),
                names=list(cache["hk_vars"]), module=np.repeat(cache["module"], 24))


def band_mean(rcs, band):
    sel = (RANGE_MID >= band[0]) & (RANGE_MID <= band[1])
    with np.errstate(invalid="ignore"):
        return np.nanmean(rcs[:, sel], axis=1)


def pair_frame(ha, hb, lon):
    """Aligned hourly arrays for one pair."""
    common, ia, ib = np.intersect1d(ha["t"], hb["t"], return_indices=True)
    t_days = (common + 0.5) / 24.0          # hour centre
    out = dict(t=t_days, sh=solar_hour(t_days, lon),
               mod=np.char.add(np.char.add(ha["module"][ia], "|"), hb["module"][ib]))
    for name, band in BANDS.items():
        a = band_mean(ha["rcs"][ia], band)
        b = band_mean(hb["rcs"][ib], band)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[name] = np.where((a > 0) & (b > 0), np.log(a / b), np.nan)
    for v in ("temp_int", "temp_ext", "temperature_optical_module", "temperature_detector", "bckgrd_rcs_0",
              "calibration_pulse", "window_transmission", "temperature_laser", "laser_energy"):
        # caches built before a variable was added simply lack it -> NaN
        for side, h, idx in (("_A", ha, ia), ("_B", hb, ib)):
            out[v + side] = h["hk"][idx, h["names"].index(v)] if v in h["names"] else np.full(len(idx), np.nan)
    return out


def is_night(sh):
    return (sh >= NIGHT_HOURS[0]) | (sh < NIGHT_HOURS[1])


def relative(fr, band):
    """lnR minus the median night-time lnR within +-REF_HALF_WINDOW days, same module pair only.

    The reference never crosses an optical-module swap, so a step in either unit cannot leak into
    the daytime anomaly.
    """
    x = fr[band]
    ref = np.full_like(x, np.nan)
    for mod in np.unique(fr["mod"]):
        same = fr["mod"] == mod
        night = is_night(fr["sh"]) & np.isfinite(x) & same
        tn, xn = fr["t"][night], x[night]
        idx = np.where(same)[0]
        lo = np.searchsorted(tn, fr["t"][idx] - REF_HALF_WINDOW)
        hi = np.searchsorted(tn, fr["t"][idx] + REF_HALF_WINDOW)
        for k, i in enumerate(idx):
            if hi[k] - lo[k] >= 4:
                ref[i] = np.median(xn[lo[k]:hi[k]])
    return x - ref


def excursion(v, pct=20):
    floor = np.nanpercentile(v, pct)
    return np.clip(v - floor, 0, None)


def ols_slope(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30 or np.std(x[ok]) == 0:
        return np.nan, np.nan, int(ok.sum())
    p, cov = np.polyfit(x[ok], y[ok], 1, cov=True)
    return p[0], np.sqrt(cov[0, 0]), int(ok.sum())


def boot_median_ci(x, n=300, rng=np.random.default_rng(0)):
    x = x[np.isfinite(x)]
    if x.size < 10:
        return np.nan, np.nan, np.nan
    meds = np.median(rng.choice(x, size=(n, x.size), replace=True), axis=1)
    return np.median(x), np.percentile(meds, 2.5), np.percentile(meds, 97.5)


def analyse_pair(fr, label):
    rows = {}
    months = np.array([(datetime(1970, 1, 1) + timedelta(days=float(t))).month for t in fr["t"]])
    for band in BANDS:
        rel = relative(fr, band)
        fr["rel_" + band] = rel
        night = is_night(fr["sh"])
        day = (fr["sh"] >= 10) & (fr["sh"] < 16)
        m_d, lo_d, hi_d = boot_median_ci(rel[day])
        dT = fr["temp_int_A"] - fr["temp_int_B"]
        dT_rel = dT - np.nanmedian(dT)
        b_n, se_n, n_n = ols_slope(dT_rel[night], rel[night])
        b_d, se_d, n_d = ols_slope(dT_rel[day], rel[day])
        b_a, se_a, n_a = ols_slope(dT_rel, rel)
        rows[band] = dict(pair=label, band=band, n_hours=int(np.isfinite(rel).sum()),
                          day_minus_night=100 * m_d, ci_lo=100 * lo_d, ci_hi=100 * hi_d,
                          sd_night=100 * np.nanstd(rel[night]), sd_day=100 * np.nanstd(rel[day]),
                          dT_sd=float(np.nanstd(dT_rel)),
                          b_dT_all=100 * b_a, se_dT_all=100 * se_a,
                          b_dT_night=100 * b_n, se_dT_night=100 * se_n,
                          b_dT_day=100 * b_d, se_dT_day=100 * se_d)
        for s, mm in SEASONS.items():
            sel = day & np.isin(months, mm)
            rows[band]["dmn_" + s] = 100 * np.nanmedian(rel[sel]) if np.isfinite(rel[sel]).sum() > 20 else np.nan
    # Regulation loss: day hours with the optical module / detector of either unit above set-point.
    rel = fr["rel_" + MAIN_BAND]
    day = (fr["sh"] >= 10) & (fr["sh"] < 18)
    exc = np.fmax(excursion(fr["temperature_optical_module_A"]), excursion(fr["temperature_optical_module_B"]))
    hot = day & (exc > 1.0)
    cool = day & (exc <= 0.2)
    rows[MAIN_BAND].update(n_hot=int(np.isfinite(rel[hot]).sum()),
                           hot_minus_cool=100 * (np.nanmedian(rel[hot]) - np.nanmedian(rel[cool]))
                           if np.isfinite(rel[hot]).sum() > 10 else np.nan)
    # Internal reference: the calibration pulse drops by day with the solar background (detector
    # sensitivity?). The firmware divides beta_raw by p_calc (beta_raw = ((P/lp)-b)/(cs o(r) p_calc) r^2).
    # Slope 0: the drop is a real sensitivity loss that the division compensates; slope -1: the drop
    # is an artefact of the reference and the division over-corrects the daytime signal.
    # Both terms are taken relative to their own night reference.
    with np.errstate(invalid="ignore", divide="ignore"):
        fr["lnp"] = np.log(fr["calibration_pulse_A"] / fr["calibration_pulse_B"])
    dlnp = relative(fr, "lnp")
    fr["rel_lnp"] = dlnp
    for sel_name, sel in (("all", np.ones_like(day)), ("day", day)):
        b, se, n = ols_slope(dlnp[sel], rel[sel])
        rows[MAIN_BAND].update(**{f"b_lnp_{sel_name}": b, f"se_lnp_{sel_name}": se, f"n_lnp_{sel_name}": n})
    rows[MAIN_BAND]["sd_dlnp_day"] = 100 * float(np.nanstd(dlnp[day]))
    # Solar background: daytime ratio vs background difference (photon-counting dead time).
    dbg = fr["bckgrd_rcs_0_A"] - fr["bckgrd_rcs_0_B"]
    bg = 0.5 * (fr["bckgrd_rcs_0_A"] + fr["bckgrd_rcs_0_B"])
    b, se, _ = ols_slope(bg[day] / np.nanmedian(bg[day]), rel[day])
    rows[MAIN_BAND].update(b_bg=100 * b, se_bg=100 * se, corr_dbg=float(np.corrcoef(
        *[v[np.isfinite(dbg) & np.isfinite(rel) & day] for v in (dbg, rel)])[0, 1]) if np.sum(
        np.isfinite(dbg) & np.isfinite(rel) & day) > 30 else np.nan)
    return rows


def figures(frames):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIGS.mkdir(parents=True, exist_ok=True)
    n = len(frames)
    # Fig A: diurnal composites per pair (main band, all seasons + JJA/DJF)
    ncol = min(n, 4)
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.9 * nrow), squeeze=False)
    for ax, (label, fr) in zip(axs.ravel(), frames):
        months = np.array([(datetime(1970, 1, 1) + timedelta(days=float(t))).month for t in fr["t"]])
        for band, col in zip(BANDS, ("#9467bd", "#1f77b4", "#17becf")):
            rel = fr["rel_" + band]
            hb = np.floor(fr["sh"]).astype(int)
            med = [np.nanmedian(rel[hb == h]) if np.isfinite(rel[hb == h]).sum() > 10 else np.nan for h in range(24)]
            ax.plot(np.arange(24) + 0.5, 100 * np.array(med), "-", color=col, lw=1.4 if band == MAIN_BAND else 0.9,
                    label=band)
        rel = fr["rel_" + MAIN_BAND]
        hb = np.floor(fr["sh"]).astype(int)
        for mm, col, lab in (((6, 7, 8), "#d62728", "JJA"), ((12, 1, 2), "#2ca02c", "DJF")):
            med = [np.nanmedian(rel[(hb == h) & np.isin(months, mm)]) if np.isfinite(
                rel[(hb == h) & np.isin(months, mm)]).sum() > 8 else np.nan for h in range(24)]
            ax.plot(np.arange(24) + 0.5, 100 * np.array(med), "--", color=col, lw=1, label=f"{MAIN_BAND} {lab}")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylim(-15, 15)
        ax.set_xlim(0, 24)
        ax.set_xticks(range(0, 25, 3))
        ax.set_xlabel("local solar hour")
        ax.set_ylabel("ln(R) - night ref [%]")
        ax.set_title(label, fontsize=9)
        ax.grid(alpha=0.3)
    axs.ravel()[0].legend(fontsize=7, ncol=2)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Same-wavelength CHM15k pairs: diurnal cycle of the signal ratio (atmosphere cancels)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "q2b_pairs_diurnal.png", dpi=130)
    plt.close(fig)

    # Fig B: ratio vs temperature difference, night and day, main band
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.9 * nrow), squeeze=False)
    for ax, (label, fr) in zip(axs.ravel(), frames):
        rel = fr["rel_" + MAIN_BAND]
        dT = fr["temp_int_A"] - fr["temp_int_B"]
        dT = dT - np.nanmedian(dT)
        night = is_night(fr["sh"])
        day = (fr["sh"] >= 10) & (fr["sh"] < 16)
        for sel, col, lab in ((night, "#1f77b4", "night"), (day, "#ff7f0e", "day 10-16")):
            ok = sel & np.isfinite(rel) & np.isfinite(dT)
            ax.plot(dT[ok], 100 * rel[ok], ".", ms=1.5, alpha=0.25, color=col)
            if ok.sum() > 30:
                edges = np.percentile(dT[ok], np.linspace(1, 99, 11))
                mid = [np.median(dT[ok & (dT >= a) & (dT < b)]) for a, b in zip(edges[:-1], edges[1:])]
                med = [100 * np.median(rel[ok & (dT >= a) & (dT < b)]) for a, b in zip(edges[:-1], edges[1:])]
                p = np.polyfit(dT[ok], rel[ok], 1)
                ax.plot(mid, med, "o-", color=col, lw=2, label=f"{lab}: {100 * p[0]:+.2f} %/K")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylim(-25, 25)
        ax.set_xlabel("T$_{int,A}$ - T$_{int,B}$ (centred) [K]")
        ax.set_ylabel("ln(R) - night ref [%]")
        ax.set_title(label, fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGS / "q2b_pairs_vs_dT.png", dpi=130)
    plt.close(fig)

    # Fig C: ratio vs the calibration-pulse ratio (is rcs_0 normalised by the internal reference?)
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.9 * nrow), squeeze=False)
    for ax, (label, fr) in zip(axs.ravel(), frames):
        rel = fr["rel_" + MAIN_BAND]
        x = fr["rel_lnp"]
        ok = np.isfinite(rel) & np.isfinite(x)
        ax.plot(100 * x[ok], 100 * rel[ok], ".", ms=1.5, alpha=0.25, color="#1f77b4")
        if ok.sum() > 50:
            edges = np.percentile(x[ok], np.linspace(1, 99, 13))
            mid = [100 * np.median(x[ok & (x >= a) & (x < b)]) for a, b in zip(edges[:-1], edges[1:])]
            med = [100 * np.median(rel[ok & (x >= a) & (x < b)]) for a, b in zip(edges[:-1], edges[1:])]
            p = np.polyfit(x[ok], rel[ok], 1)
            ax.plot(mid, med, "o-", color="#d62728", lw=2, label=f"slope {p[0]:+.2f}")
        lim = 30
        ax.plot([-lim, lim], [lim, -lim], "k--", lw=0.8, label="slope -1: p_cal artefact, over-corrected")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("ln(p$_{cal,A}$/p$_{cal,B}$) - night ref [%]")
        ax.set_ylabel("ln(R) - night ref [%]")
        ax.set_title(label, fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGS / "q2b_pairs_vs_calpulse.png", dpi=130)
    plt.close(fig)


def main():
    lon_of = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    site_of = {(w, i): s for w, i, s, la, lo in chm15k_streams()}
    rows, frames = [], []
    for wmo, units in PAIR_SITES.items():
        hours = {}
        for u in units:
            c = load_cache(wmo, u)
            if c is not None:
                hours[u] = hourly(c)
        for a, b in itertools.combinations(sorted(hours), 2):
            label = f"{site_of[(wmo, a)][:18]} {a}/{b}"
            fr = pair_frame(hours[a], hours[b], lon_of[(wmo, a)])
            res = analyse_pair(fr, label)
            frames.append((label, fr))
            for band, r in res.items():
                rows.append(r)
            r = res[MAIN_BAND]
            print(f"  {label:28s} hours {r['n_hours']:5d}  day-night {r['day_minus_night']:+5.1f} % "
                  f"[{r['ci_lo']:+.1f},{r['ci_hi']:+.1f}]  sd night {r['sd_night']:.1f} day {r['sd_day']:.1f}  "
                  f"dT sd {r['dT_sd']:.1f} K  b_dT night {r['b_dT_night']:+.2f}+/-{r['se_dT_night']:.2f} "
                  f"day {r['b_dT_day']:+.2f}+/-{r['se_dT_day']:.2f} %/K  hot-cool {r.get('hot_minus_cool', np.nan):+.1f} % "
                  f"(n {r.get('n_hot', 0)})  bg {r.get('b_bg', np.nan):+.2f}  "
                  f"dlnR/dln(p_A/p_B) all {r['b_lnp_all']:+.2f}+/-{r['se_lnp_all']:.2f} "
                  f"day {r['b_lnp_day']:+.2f}+/-{r['se_lnp_day']:.2f} (sd dlnp day {r['sd_dlnp_day']:.1f} %)",
                  flush=True)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("pair", "band"), k))
    with open(OUT / "q2b_pairs.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    figures(frames)
    np.savez_compressed(OUT / "q2b_frames.npz", **{f"{lab}|{k}": v for lab, fr in frames for k, v in fr.items()})
    print("Q2B_DONE")


if __name__ == "__main__":
    main()
