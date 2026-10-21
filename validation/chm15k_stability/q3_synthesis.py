# -*- coding: utf-8 -*-
"""Synthesis of report 19: does one per-unit temperature coefficient explain the daytime behaviour?

1. Per-unit consistency of the temperature coefficient between two independent estimators:
   the synoptic-anomaly model M2 on the nightly Rayleigh constants (q1, night-to-night, weeks) and
   the intra-night slope with elapsed time held fixed (q2a, hours).
2. Prediction of the daytime signal anomaly from that night-time coefficient alone:
     predicted day-night of ln(signal) = kappa_T(unit) x (T_int(day) - T_int(night ref))
   compared with
     * the observed day-night ratio of same-wavelength pairs  (pred_A - pred_B),
     * the observed day-night CHM15k / CL61 ratio              (pred_CHM; CL61 assumed stable).
   If the points fall on 1:1 the diurnal pattern is the temperature coefficient seen by day;
   a residual is the part driven by something else (solar background / detector).

Run after q1, q2a:  python validation/chm15k_stability/q3_synthesis.py
"""
from __future__ import annotations

import csv
import itertools
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, OUT, PAIR_SITES, chm15k_streams, load_cache  # noqa: E402
from q2b_pairs import MAIN_BAND, hourly, is_night, pair_frame, relative  # noqa: E402
from q2d_cross_cl61 import SITES as CL61_SITES  # noqa: E402

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")
DAY = (10.0, 16.0)


def read_csv(path):
    return list(csv.DictReader(open(path, encoding="utf-8")))


def kappas_intranight():
    out = {}
    for r in read_csv(OUT / "q2a_intranight.csv"):
        out[r["stream"]] = (float(r["b_T_ctrl"]), float(r["se_T_ctrl"]))
    return out


def kappas_m2():
    """Per stream: M2 slope of the era with the most nights (temp_int, arm A)."""
    best = {}
    for r in read_csv(OUT / "q1_slopes.csv"):
        if r["arm"] != "A" or r["tvar"] != "temp_int" or r["b2"] in ("", "nan"):
            continue
        n = int(r["n"])
        if r["stream"] not in best or n > best[r["stream"]][2]:
            best[r["stream"]] = (float(r["b2"]), float(r["se2"]), n)
    return best


def day_temperature_rise(fr, key, valid):
    """Median of T - night reference over daytime hours where the ratio is valid."""
    fr["_T"] = fr[key]
    rel_t = relative(fr, "_T")
    day = (fr["sh"] >= DAY[0]) & (fr["sh"] < DAY[1]) & valid
    return float(np.nanmedian(rel_t[day])), rel_t


def hourly_composite(x, sh):
    hb = np.floor(sh).astype(int)
    return np.array([np.nanmedian(x[hb == h]) if np.isfinite(x[hb == h]).sum() > 10 else np.nan for h in range(24)])


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bi = kappas_intranight()
    bm = kappas_m2()
    lon_of = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    site_of = {(w, i): s for w, i, s, la, lo in chm15k_streams()}

    # 1. estimator consistency
    common = sorted(set(bi) & set(bm))
    x = np.array([bm[k][0] for k in common])
    xe = np.array([bm[k][1] for k in common])
    y = np.array([bi[k][0] for k in common])
    ye = np.array([bi[k][1] for k in common])
    ok = np.isfinite(x) & np.isfinite(y)
    r_est = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 3 else np.nan
    print(f"Per-unit kappa: M2 (night-to-night) vs intra-night: N={ok.sum()}, r={r_est:+.2f}, "
          f"median M2 {np.median(x[ok]):+.2f}, median intra-night {np.median(y[ok]):+.2f} %/K")

    # 2. daytime prediction
    rows, curves = [], []
    for wmo, units in PAIR_SITES.items():
        hrs = {u: hourly(load_cache(wmo, u)) for u in units if load_cache(wmo, u) is not None}
        for a, b in itertools.combinations(sorted(hrs), 2):
            ka, kb = f"{wmo}_{a}", f"{wmo}_{b}"
            if ka not in bi or kb not in bi:
                continue
            fr = pair_frame(hrs[a], hrs[b], lon_of[(wmo, a)])
            rel = relative(fr, MAIN_BAND)
            valid = np.isfinite(rel)
            dTa, relTa = day_temperature_rise(fr, "temp_int_A", valid)
            dTb, relTb = day_temperature_rise(fr, "temp_int_B", valid)
            day = (fr["sh"] >= DAY[0]) & (fr["sh"] < DAY[1])
            obs = 100 * float(np.nanmedian(rel[day]))
            pred = bi[ka][0] * dTa - bi[kb][0] * dTb
            pe = np.hypot(bi[ka][1] * dTa, bi[kb][1] * dTb)
            label = f"{site_of[(wmo, a)][:10]} {a}/{b}"
            rows.append(dict(kind="pair", label=label, obs=obs, pred=pred, pred_se=pe, dT_A=dTa, dT_B=dTb,
                             kappa_A=bi[ka][0], kappa_B=bi[kb][0]))
            # predicted curve on the SAME clear hours as the observed one
            relTa = np.where(valid, relTa, np.nan)
            relTb = np.where(valid, relTb, np.nan)
            curves.append((label, 100 * hourly_composite(rel, fr["sh"]),
                           bi[ka][0] * hourly_composite(relTa, fr["sh"]) - bi[kb][0] * hourly_composite(relTb, fr["sh"])))
    for wmo, chm, cl, name in CL61_SITES:
        kc = f"{wmo}_{chm}"
        cc, cv = load_cache(wmo, chm), load_cache(wmo, cl)
        if kc not in bi or cc is None or cv is None:
            continue
        fr = pair_frame(hourly(cc), hourly(cv), lon_of[(wmo, chm)])
        rel = relative(fr, MAIN_BAND)
        dT, relT = day_temperature_rise(fr, "temp_int_A", np.isfinite(rel))
        day = (fr["sh"] >= DAY[0]) & (fr["sh"] < DAY[1])
        obs = 100 * float(np.nanmedian(rel[day]))
        rows.append(dict(kind="CHM/CL61", label=f"{name} CHM/CL61", obs=obs, pred=bi[kc][0] * dT,
                         pred_se=bi[kc][1] * dT, dT_A=dT, dT_B=np.nan, kappa_A=bi[kc][0], kappa_B=np.nan))
        relT = np.where(np.isfinite(rel), relT, np.nan)
        curves.append((f"{name} CHM/CL61", 100 * hourly_composite(rel, fr["sh"]),
                       bi[kc][0] * hourly_composite(relT, fr["sh"])))
    for r in rows:
        print(f"  {r['kind']:9s} {r['label']:26s} kappa {r['kappa_A']:+.2f}/{r['kappa_B']:+.2f} %/K  "
              f"dT_day {r['dT_A']:+.1f}/{r['dT_B']:+.1f} K  predicted {r['pred']:+6.1f} +/- {r['pred_se']:.1f} %  "
              f"observed {r['obs']:+6.1f} %")
    with open(OUT / "q3_synthesis.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    pairs = [r for r in rows if r["kind"] == "pair"]
    po = np.array([r["obs"] for r in pairs])
    pp = np.array([r["pred"] for r in pairs])
    r_pair = float(np.corrcoef(po, pp)[0, 1])
    slope_pair = float(np.sum(po * pp) / np.sum(pp ** 2))          # through-origin regression
    print(f"Pairs: r(observed, predicted) = {r_pair:+.2f}, observed = {slope_pair:.2f} x predicted (through origin)")

    # 3. What a temperature normalisation would buy: robust scatter of the nightly constant
    #    (successive-difference, i.e. short-term) and of the whole record, before and after
    #    C_norm = C x exp(-kappa (T_int - T_ref)) with the unit's own intra-night kappa.
    from q1_temperature import night_hk, robust_mask
    from common import load_constants
    gains = []
    for key, (kappa, se) in bi.items():
        wmo, ident = key.rsplit("_", 1)
        cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
        if cons is None or cache is None or not np.isfinite(kappa):
            continue
        T = night_hk(cache, cons["t0"], cons["t1"])[:, list(cache["hk_vars"]).index("temp_int")]
        logc = np.log(np.where(cons["c"] > 0, cons["c"], np.nan))
        main = cons["module"] == max(set(cons["module"]), key=list(cons["module"]).count)
        m = main & np.isfinite(T) & np.isfinite(logc)
        m[m] = robust_mask(logc[m])
        if m.sum() < 30:
            continue
        lc, tt = logc[m], T[m]
        ln = lc - kappa / 100 * (tt - np.median(tt))

        def rsd(v):
            return 100 * 1.4826 * np.median(np.abs(v - np.median(v)))

        def ssd(v):          # successive-difference scatter (short-term), per single night
            return 100 * 1.4826 * np.median(np.abs(np.diff(v))) / np.sqrt(2)
        gains.append(dict(stream=key, n=int(m.sum()), kappa=kappa, sd_raw=rsd(lc), sd_norm=rsd(ln),
                          ssd_raw=ssd(lc), ssd_norm=ssd(ln)))
    if gains:
        g = {k: np.array([r[k] for r in gains]) for k in ("sd_raw", "sd_norm", "ssd_raw", "ssd_norm")}
        print(f"T normalisation (unit kappa, N={len(gains)}): record scatter median {np.median(g['sd_raw']):.1f} -> "
              f"{np.median(g['sd_norm']):.1f} %, night-to-night {np.median(g['ssd_raw']):.1f} -> "
              f"{np.median(g['ssd_norm']):.1f} %; units improved: record {np.sum(g['sd_norm'] < g['sd_raw'])}/"
              f"{len(gains)}, night-to-night {np.sum(g['ssd_norm'] < g['ssd_raw'])}/{len(gains)}")
        with open(OUT / "q3_normalisation_gain.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(gains[0].keys()))
            w.writeheader()
            w.writerows(gains)

    # figure
    FIGS.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(21, 8.5))
    gs = fig.add_gridspec(2, 6)
    ax = fig.add_subplot(gs[:, 0:2])
    ax.errorbar(x[ok], y[ok], xerr=xe[ok], yerr=ye[ok], fmt="o", ms=4, lw=0.6, alpha=0.7, color="#1f77b4")
    lim = [-4, 8]
    ax.plot(lim, lim, "k--", lw=0.8, label="1:1")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("M2: night-to-night synoptic anomaly [%/K]")
    ax.set_ylabel("intra-night, elapsed time held fixed [%/K]")
    ax.set_title(f"(a) temperature coefficient per unit, two estimators\nN={ok.sum()}, r={r_est:+.2f}", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax = fig.add_subplot(gs[:, 2:4])
    for kind, mk, col in (("pair", "o", "#d62728"), ("CHM/CL61", "s", "#2ca02c")):
        sel = [r for r in rows if r["kind"] == kind]
        ax.errorbar([r["pred"] for r in sel], [r["obs"] for r in sel], xerr=[r["pred_se"] for r in sel],
                    fmt=mk, ms=6, color=col, label=kind)
        for r in sel:
            ax.annotate(r["label"].replace("AMSTERDAM_", "AMS "), (r["pred"], r["obs"]), fontsize=7,
                        xytext=(4, 2), textcoords="offset points")
    lim = [-25, 35]
    ax.plot(lim, lim, "k--", lw=0.8, label="1:1")
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"predicted from night-time $\kappa$ x daytime warming [%]")
    ax.set_ylabel("observed day - night signal ratio (10-16 h solar) [%]")
    ax.set_title(f"(b) daytime anomaly predicted by the temperature coefficient\npairs r={r_pair:+.2f}, "
                 f"obs = {slope_pair:.2f} x pred", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    preferred = ("AMSTERDAM_ B/C", "BUDAPEST A/B", "Lindenberg CHM/CL61", "Aosta CHM/CL61")
    show = [c for want in preferred for c in curves if c[0] == want] or curves[:4]
    for j, (label, o, p) in enumerate(show[:4]):
        ax = fig.add_subplot(gs[j // 2, 4 + j % 2])
        ax.plot(np.arange(24) + 0.5, o, "-", color="#1f77b4", lw=1.8, label="observed")
        ax.plot(np.arange(24) + 0.5, p, "--", color="#d62728", lw=1.8, label=r"predicted ($\kappa$ x dT)")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(0, 24)
        ax.set_xticks(range(0, 25, 6))
        ax.set_title(label, fontsize=9)
        ax.set_xlabel("local solar hour")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=7)
            ax.set_ylabel("relative to night [%]")
    fig.tight_layout()
    fig.savefig(FIGS / "q3_synthesis.png", dpi=130)
    print(f"-> {FIGS / 'q3_synthesis.png'}")


def generation_and_daytime_figure():
    """Temperature coefficient by optical-module generation (TUBxx = production year) and the
    resulting daytime excess of the effective constant over the night-time Rayleigh value."""
    import collections
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bi = kappas_intranight()
    qc = {r["stream"]: r for r in read_csv(OUT / "q2c_hk_diurnal.csv")}
    gen = collections.defaultdict(list)
    excess = {"all seasons": [], "JJA": []}
    rows = []
    for key, (b, se) in bi.items():
        if not np.isfinite(b):
            continue
        c = load_cache(*key.rsplit("_", 1))
        mods = [m for m in c["module"] if m]
        g = max(set(mods), key=mods.count)[:5] if mods else ""
        if g:
            gen[g].append(b)
        if key in qc:
            ea = b * float(qc[key]["temp_int_noon"])
            ej = b * float(qc[key]["temp_int_noon_JJA"])
            excess["all seasons"].append(ea)
            excess["JJA"].append(ej)
            rows.append(dict(stream=key, generation=g, kappa=b, noon_dT=float(qc[key]["temp_int_noon"]),
                             noon_dT_JJA=float(qc[key]["temp_int_noon_JJA"]), excess=ea, excess_JJA=ej))
    with open(OUT / "q3_daytime_excess.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    fig, axs = plt.subplots(1, 2, figsize=(19, 5.5), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]
    keys = sorted(gen)
    ax.boxplot([gen[k] for k in keys], showfliers=False)
    for i, k in enumerate(keys, 1):
        ax.plot(np.full(len(gen[k]), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(gen[k])), gen[k], ".",
                color="#1f77b4", alpha=0.6)
    ax.set_xticks(range(1, len(keys) + 1))
    ax.set_xticklabels([f"{k}\n(n={len(gen[k])})" for k in keys], fontsize=8)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_ylim(-2, 5)
    ax.set_xlabel("optical-module generation (TUBxx = production year 20xx)")
    ax.set_ylabel("intra-night d ln C / dT$_{int}$ [%/K]")
    ax.set_title("(a) temperature coefficient by optical-module generation", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1]
    bins = np.linspace(-10, 30, 41)
    for (name, x), col in zip(excess.items(), ("#1f77b4", "#d62728")):
        x = np.array(x)
        x = x[np.isfinite(x)]
        ax.hist(np.clip(x, -10, 30), bins=bins, histtype="step", lw=1.8, color=col,
                label=f"{name}: median {np.median(x):+.1f} %, >10 %: {100 * np.mean(x > 10):.0f} % of units")
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel(r"noon effective constant vs night calibration [%]  ($\kappa$ x noon T$_{int}$ rise)")
    ax.set_ylabel("number of units")
    ax.set_title("(b) daytime excess not seen by the night-only Rayleigh calibration", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "q3_generation_daytime.png", dpi=130)
    print(f"-> {FIGS / 'q3_generation_daytime.png'}")


if __name__ == "__main__":
    main()
    generation_and_daytime_figure()
