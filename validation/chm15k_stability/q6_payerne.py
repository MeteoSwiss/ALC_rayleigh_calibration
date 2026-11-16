# -*- coding: utf-8 -*-
"""Why does Payerne not fit? (report 19, part 6)

At Payerne the night-time temperature coefficient of the CHM15k (+2.0 %/K) and its daytime warming
(+10 K) predict a +20 % day-night CHM15k/CL61 ratio; +2.7 % is observed. Tests:

A. Triangle CHM15k (A) / CL31 (B) / CL61 (C): which instrument carries the daytime cycle?
B. Diurnal cycle of each instrument's own housekeeping (CHM15k T_int, T_lom; Vaisala internal and laser
   temperature, laser energy) at Payerne and at the other CL61 sites.
C. Daytime coefficient of the CHM15k measured directly: day-to-day variations of the CHM15k/reference ratio
   at a FIXED solar hour against the CHM15k T_int anomaly at that hour (removes the mean diurnal cycle, so
   aerosol-diurnal and Angstrom effects drop out).
D. History of the Payerne coefficient (11-yr archive eras, network 2025-26, intra-night).

Run:  python validation/chm15k_stability/q6_payerne.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, OUT, chm15k_streams, load_cache  # noqa: E402
from q2b_pairs import boot_median_ci, hourly, pair_frame, relative  # noqa: E402

warnings.filterwarnings("ignore")
BANDS = ("1.5-2.5 km", "2.5-4.0 km")
DAY = (10.0, 16.0)
SITES = {"Payerne": ("0-20000-0-06610", "A", "C", "B"), "Lindenberg": ("0-20000-0-10393", "0", "C", None),
         "Camborne": ("0-20000-0-03808", "A", "C", None), "Aosta": ("0-380-5-1", "0", "B", None)}


def composite(x, sh):
    hb = np.floor(sh).astype(int)
    return np.array([np.nanmedian(x[hb == h]) if np.isfinite(x[hb == h]).sum() > 8 else np.nan for h in range(24)])


def anom_same_hour(fr, x, half=15.0):
    out = np.full_like(x, np.nan)
    hb = np.floor(fr["sh"]).astype(int)
    for h in range(24):
        idx = np.where((hb == h) & np.isfinite(x))[0]
        t = fr["t"][idx]
        for k, i in enumerate(idx):
            m = np.abs(t - t[k]) <= half
            if m.sum() >= 5:
                out[i] = x[i] - np.median(x[idx[m]])
    return out


def fit(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 40 or np.std(x[ok]) == 0:
        return np.nan, np.nan, int(ok.sum())
    p, c = np.polyfit(x[ok], y[ok], 1, cov=True)
    return 100 * p[0], 100 * np.sqrt(c[0, 0]), int(ok.sum())


def rel_key(fr, key):
    fr["_x"] = fr[key]
    return relative(fr, "_x")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lon = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    hrs = {}

    def H(w, u):
        if (w, u) not in hrs:
            c = load_cache(w, u)
            hrs[(w, u)] = hourly(c) if c is not None else None
        return hrs[(w, u)]

    rows, comps, hkcomps = [], {}, {}
    wP = SITES["Payerne"][0]
    # A. Payerne triangle, full overlap of each pair and the common CL61 period
    t61 = H(wP, "C")["t"]
    for (a, b, lab) in (("A", "C", "CHM15k/CL61"), ("A", "B", "CHM15k/CL31"), ("C", "B", "CL61/CL31")):
        fr = pair_frame(H(wP, a), H(wP, b), lon[(wP, "A")])
        common = np.isin(np.floor(fr["t"] * 24).astype(np.int64), t61)
        for period, sel_p in (("all", np.ones(len(fr["t"]), bool)), ("CL61 period", common)):
            for band in BANDS:
                rel = relative(fr, band)
                rel = np.where(sel_p, rel, np.nan)
                day = (fr["sh"] >= DAY[0]) & (fr["sh"] < DAY[1])
                m, lo, hi = boot_median_ci(rel[day])
                rows.append(dict(site="Payerne", ratio=lab, period=period, band=band, n=int(np.isfinite(rel).sum()),
                                 day_minus_night=100 * m, ci_lo=100 * lo, ci_hi=100 * hi))
                if band == BANDS[0]:
                    comps[(lab, period)] = 100 * composite(rel, fr["sh"])
                print(f"  Payerne {lab:12s} {period:12s} {band}: day-night {100 * m:+6.1f} % [{100 * lo:+.1f},{100 * hi:+.1f}]"
                      f"  (hours {int(np.isfinite(rel).sum())})", flush=True)

    # B. housekeeping diurnal cycles (relative to the night) of each instrument
    for site, (w, chm, cl61, cl31) in SITES.items():
        fr = pair_frame(H(w, chm), H(w, cl61), lon[(w, chm)])
        for key, lab in (("temp_int_A", "CHM15k T_int"), ("temperature_optical_module_A", "CHM15k T_lom"),
                         ("temp_int_B", "CL61 T_int"), ("temperature_laser_B", "CL61 T_laser"),
                         ("laser_energy_B", "CL61 laser energy")):
            hkcomps[(site, lab)] = composite(rel_key(fr, key), fr["sh"])
        if cl31:
            fr3 = pair_frame(H(w, chm), H(w, cl31), lon[(w, chm)])
            for key, lab in (("temperature_laser_B", "CL31 T_laser"), ("laser_energy_B", "CL31 laser energy")):
                hkcomps[(site, lab)] = composite(rel_key(fr3, key), fr3["sh"])
        noon = slice(11, 15)
        print(f"  {site:10s} noon vs night: " + "  ".join(
            f"{lab} {np.nanmean(v[noon]):+.1f}" for (s, lab), v in hkcomps.items() if s == site), flush=True)

    # C. daytime coefficient of the CHM15k from day-to-day variations at a fixed solar hour
    kappa_rows = []
    for site, (w, chm, cl61, cl31) in SITES.items():
        for ref, lab in ((cl61, "CL61"), (cl31, "CL31")):
            if ref is None:
                continue
            fr = pair_frame(H(w, chm), H(w, ref), lon[(w, chm)])
            sel = (fr["sh"] >= 9) & (fr["sh"] < 17) & np.isfinite(fr[BANDS[0]])
            x = anom_same_hour(fr, np.where(sel, fr["temp_int_A"], np.nan))
            y = anom_same_hour(fr, np.where(sel, fr[BANDS[0]], np.nan))
            b, se, n = fit(x, y)
            # the reference's own temperature, same treatment (does the Vaisala co-vary?)
            xr = anom_same_hour(fr, np.where(sel, fr["temp_int_B"], np.nan)) if lab == "CL61" else None
            br = fit(xr, y) if xr is not None else (np.nan, np.nan, 0)
            r_tt = float(np.corrcoef(x[np.isfinite(x) & np.isfinite(xr)], xr[np.isfinite(x) & np.isfinite(xr)])[0, 1]) \
                if xr is not None and np.sum(np.isfinite(x) & np.isfinite(xr)) > 30 else np.nan
            kappa_rows.append(dict(site=site, ref=lab, b_day=b, se_day=se, n=n, sd_dT=float(np.nanstd(x)),
                                  b_vs_ref_T=br[0], se_vs_ref_T=br[1], r_Tchm_Tref=r_tt))
            print(f"  {site:10s} CHM15k/{lab}: daytime same-hour anomaly slope vs CHM T_int {b:+.2f}+/-{se:.2f} %/K "
                  f"(n {n}, dT sd {np.nanstd(x):.2f} K); vs ref T_int {br[0]:+.2f}+/-{br[1]:.2f}; r(T_chm,T_ref) {r_tt:+.2f}",
                  flush=True)

    # D. Payerne coefficient history
    hist = []
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["tvar"] == "temp_int" and (r["stream"] == "pay" or r["stream"].startswith("0-20000-0-06610_A")):
            hist.append((r["arm"], r["era"], int(r["n"]), float(r["b2"]), float(r["se2"]), float(r["std_logc"])))
    for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8")):
        if r["stream"] == "0-20000-0-06610_A":
            print(f"  Payerne intra-night: b_T_ctrl {float(r['b_T_ctrl']):+.2f}+/-{float(r['se_T_ctrl']):.2f}, "
                  f"hourly scatter {float(r['sd_obs']):.1f} % vs noise {float(r['sd_noise']):.1f} % (network median ~4.3)")
    for h in hist:
        print(f"  Payerne M2 arm {h[0]} era {h[1]}: n={h[2]} b2 {h[3]:+.2f}+/-{h[4]:.2f} %/K, record scatter {h[5]:.0f} %")

    for name, rr in (("q6_payerne_triangle.csv", rows), ("q6_daytime_kappa.csv", kappa_rows)):
        with open(OUT / name, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rr[0].keys()))
            w.writeheader()
            w.writerows(rr)

    # Figure
    fig, axs = plt.subplots(1, 3, figsize=(21, 5.6))
    x = np.arange(24) + 0.5
    ax = axs[0]
    for (lab, period), col in zip([k for k in comps if k[1] == "CL61 period"], ("#d62728", "#1f77b4", "#2ca02c")):
        ax.plot(x, comps[(lab, period)], "-", color=col, lw=2, label=f"{lab} (common CL61 period)")
    ax.plot(x, comps[("CHM15k/CL31", "all")], "--", color="#1f77b4", lw=1, label="CHM15k/CL31 (2025-26)")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 3))
    ax.set_ylim(-20, 30)
    ax.set_xlabel("local solar hour")
    ax.set_ylabel(f"ln(ratio) - night ref, {BANDS[0]} [%]")
    ax.set_title("(a) Payerne triangle: which instrument moves by day?", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    for ax, labs, title, unit in ((axs[1], ("CHM15k T_int", "CHM15k T_lom", "CL61 T_int", "CL61 T_laser", "CL31 T_laser"),
                                   "(b) instrument temperatures vs night", "K"),
                                  (axs[2], ("CL61 laser energy", "CL31 laser energy"),
                                   "(c) Vaisala laser energy vs night", "%-pt")):
        for site, ls in zip(SITES, ("-", "--", ":", "-.")):
            for lab, col in zip(labs, ("#d62728", "#ff7f0e", "#1f77b4", "#17becf", "#2ca02c")):
                if (site, lab) in hkcomps and np.isfinite(hkcomps[(site, lab)]).any():
                    ax.plot(x, hkcomps[(site, lab)], ls, color=col, lw=2 if site == "Payerne" else 1,
                            label=f"{site} {lab}")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(0, 24)
        ax.set_xticks(range(0, 25, 3))
        ax.set_xlabel("local solar hour")
        ax.set_ylabel(f"relative to night [{unit}]")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q6_payerne.png", dpi=130)
    print(f"-> {FIGS / 'q6_payerne.png'}")


if __name__ == "__main__":
    main()
