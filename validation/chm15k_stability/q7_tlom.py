# -*- coding: utf-8 -*-
"""Does the laser-optical-module temperature (temp_lom = temperature_optical_module) matter? (report 19, part 7)

temp_lom is regulated (~33.85 degC) and reported in whole degrees, so it carries information only when
  (i) it flips by 1 K at night (regulation dither), or
  (ii) it leaves regulation on hot afternoons (Rosu et al. 2026: up to 40-48 degC).
Tests:
  N. night, every unit: 10-min clear-sky signal in the Rayleigh window, intra-night anomalies, joint
     regression on dT_int, dT_lom and elapsed time -> partial slope of T_lom with T_int held fixed.
  D. day, same-wavelength pairs and CHM15k/CL61 sites: same-hour anomalies (mean diurnal cycle removed)
     of the ratio regressed jointly on the T_int and T_lom anomalies (pairs: A-B differences), plus the
     residual after the T_int model in hours with a T_lom excursion > 1 K vs none.

Run:  python validation/chm15k_stability/q7_tlom.py
"""
from __future__ import annotations

import csv
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CENSUS, FIGS, OUT, PAIR_SITES, RANGE_MID, chm15k_streams, load_cache, load_constants, stream_key  # noqa: E402
from q2a_intranight import MIN_CLEAR_FRAC, flatten  # noqa: E402
from q2b_pairs import MAIN_BAND, hourly, pair_frame  # noqa: E402
from q6_payerne import SITES as CL61_SITES, anom_same_hour  # noqa: E402

warnings.filterwarnings("ignore")


def ols(X, y):
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[ok], y[ok]
    if len(y) < 60:
        return None
    med = np.median(y)
    keep = np.abs(y - med) <= 4 * 1.4826 * np.median(np.abs(y - med)) + 1e-9
    X, y = X[keep], y[keep]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    cov = (res @ res / max(len(y) - X.shape[1], 1)) * np.linalg.pinv(X.T @ X)
    return coef, np.sqrt(np.diag(cov)), len(y)


def night_test():
    alts = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    rows = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
        if cons is None or cache is None:
            continue
        names = list(cache["hk_vars"])
        fl = flatten(cache)
        cos_t = np.cos(np.deg2rad(float(cache["tilt"]) if np.isfinite(cache["tilt"]) else 0.0))
        Y, Ti, Tl, H = [], [], [], []
        for k in range(len(cons["t"])):
            i0 = np.searchsorted(fl["t"], cons["t0"][k])
            i1 = np.searchsorted(fl["t"], cons["t1"][k], side="right")
            if i1 - i0 < 24:
                continue
            win = (RANGE_MID >= (cons["zb"][k] - alts[(wmo, ident)]) / cos_t) & \
                  (RANGE_MID <= (cons["zt"][k] - alts[(wmo, ident)]) / cos_t)
            if not win.any():
                continue
            s = np.nanmean(fl["rcs"][i0:i1][:, win], axis=1)
            clear = fl["n_clear"][i0:i1] >= MIN_CLEAR_FRAC * np.maximum(fl["n_all"][i0:i1], 1)
            s = np.where(clear & (s > 0), s, np.nan)
            if np.isfinite(s).sum() < 12:
                continue
            ti = fl["hk"][i0:i1, names.index("temp_int")]
            tl = fl["hk"][i0:i1, names.index("temperature_optical_module")]
            Y.append(np.log(s / np.nanmedian(s)))
            Ti.append(ti - np.nanmedian(ti))
            Tl.append(tl - np.nanmedian(tl))
            H.append((fl["t"][i0:i1] - fl["t"][i0]) * 24)
        if len(Y) < 25:
            continue
        y, ti, tl, h = (np.concatenate(v) for v in (Y, Ti, Tl, H))
        frac_lom = float(np.mean(np.abs(tl[np.isfinite(tl)]) > 0.5))
        r = ols(np.column_stack([np.ones_like(y), ti, tl, h]), y)
        if r is None:
            continue
        coef, se, n = r
        rows.append(dict(stream=stream_key(wmo, ident), site=site, n=n, frac_lom_moving=frac_lom,
                         b_Tint=100 * coef[1], se_Tint=100 * se[1], b_Tlom=100 * coef[2], se_Tlom=100 * se[2]))
    return rows


def day_test():
    lon = {(w, i): lo for w, i, s, la, lo in chm15k_streams()}
    rows = []
    cases = []
    for wmo, units in PAIR_SITES.items():
        hrs = {u: hourly(load_cache(wmo, u)) for u in units}
        for a, b in itertools.combinations(sorted(hrs), 2):
            cases.append((f"{wmo}_{a}/{b}", pair_frame(hrs[a], hrs[b], lon[(wmo, a)]), True))
    for site, (w, chm, cl61, _) in CL61_SITES.items():
        cases.append((f"{site} CHM15k/CL61", pair_frame(hourly(load_cache(w, chm)), hourly(load_cache(w, cl61)),
                                                         lon[(w, chm)]), False))
    for label, fr, is_pair in cases:
        sel = (fr["sh"] >= 9) & (fr["sh"] < 17) & np.isfinite(fr[MAIN_BAND])
        if is_pair:
            ti = fr["temp_int_A"] - fr["temp_int_B"]
            tl = fr["temperature_optical_module_A"] - fr["temperature_optical_module_B"]
            exc = np.fmax(fr["temperature_optical_module_A"], fr["temperature_optical_module_B"])
        else:
            ti = fr["temp_int_A"]
            tl = fr["temperature_optical_module_A"]
            exc = fr["temperature_optical_module_A"]
        exc = exc - np.nanpercentile(exc, 20)
        y = anom_same_hour(fr, np.where(sel, fr[MAIN_BAND], np.nan))
        ta = anom_same_hour(fr, np.where(sel, ti, np.nan))
        la = anom_same_hour(fr, np.where(sel, tl, np.nan))
        r1 = ols(np.column_stack([np.ones_like(y), ta]), y)
        r2 = ols(np.column_stack([np.ones_like(y), ta, la]), y)
        row = dict(case=label, n=r2[2] if r2 else 0, sd_lom_anom=float(np.nanstd(la)),
                   frac_hot=float(np.nanmean(exc[sel] > 1.0)))
        if r1:
            row.update(b_Tint_only=100 * r1[0][1], se_Tint_only=100 * r1[1][1])
        if r2:
            row.update(b_Tint=100 * r2[0][1], se_Tint=100 * r2[1][1], b_Tlom=100 * r2[0][2], se_Tlom=100 * r2[1][2])
        # residual of the T_int-only model in hot (regulation lost) vs regulated hours
        if r1:
            res = y - (r1[0][0] + r1[0][1] * ta)
            hot, cool = sel & (exc > 1.0) & np.isfinite(res), sel & (exc <= 0.2) & np.isfinite(res)
            row.update(n_hot=int(hot.sum()),
                       resid_hot_minus_cool=100 * (np.median(res[hot]) - np.median(res[cool])) if hot.sum() > 15 else np.nan)
        rows.append(row)
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    print("N. night (regulation dither):")
    nrows = night_test()
    bl = np.array([r["b_Tlom"] for r in nrows])
    sl = np.array([r["se_Tlom"] for r in nrows])
    bi = np.array([r["b_Tint"] for r in nrows])
    fm = np.array([r["frac_lom_moving"] for r in nrows])
    ok = np.isfinite(bl) & np.isfinite(sl) & (sl > 0) & (fm > 0.02)
    w = 1 / sl[ok] ** 2
    print(f"  units {len(nrows)}, with T_lom moving at night (>2 % of blocks) {ok.sum()}; "
          f"T_lom partial slope median {np.median(bl[ok]):+.2f} %/K, pooled {np.sum(w * bl[ok]) / np.sum(w):+.2f} "
          f"+/- {1 / np.sqrt(np.sum(w)):.2f}, sig + {np.sum(bl[ok] > 3 * sl[ok])} / - {np.sum(bl[ok] < -3 * sl[ok])}; "
          f"T_int partial median {np.median(bi[ok]):+.2f} %/K")
    print("D. day (pairs, CHM15k/CL61), same-hour anomalies:")
    drows = day_test()
    for r in drows:
        print(f"  {r['case']:30s} n {r['n']:5d}  T_int only {r.get('b_Tint_only', np.nan):+.2f}+/-{r.get('se_Tint_only', np.nan):.2f}"
              f"  joint T_int {r.get('b_Tint', np.nan):+.2f}+/-{r.get('se_Tint', np.nan):.2f}  T_lom {r.get('b_Tlom', np.nan):+.2f}"
              f"+/-{r.get('se_Tlom', np.nan):.2f} %/K  (T_lom anom sd {r['sd_lom_anom']:.2f} K, hot hours {r['frac_hot'] * 100:.0f} %, "
              f"resid hot-cool {r.get('resid_hot_minus_cool', np.nan):+.1f} %, n_hot {r.get('n_hot', 0)})")
    for name, rr in (("q7_tlom_night.csv", nrows), ("q7_tlom_day.csv", drows)):
        keys = sorted({k for r in rr for k in r}, key=lambda k: (k not in ("stream", "site", "case"), k))
        with open(OUT / name, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys)
            wr.writeheader()
            wr.writerows(rr)
    fig, axs = plt.subplots(1, 2, figsize=(18, 5.5))
    ax = axs[0]
    ax.errorbar(bi[ok], bl[ok], yerr=sl[ok], fmt="o", ms=3, lw=0.5, alpha=0.6)
    ax.axhline(0, color="grey", lw=0.8)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("partial slope on T_int [%/K]")
    ax.set_ylabel("partial slope on T_lom [%/K]")
    ax.set_ylim(-15, 15)
    ax.set_title(f"(a) night, intra-night joint fit: T_lom vs T_int (N={ok.sum()} units with T_lom dither)", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1]
    labs = [r["case"].replace("0-20000-0-", "").replace("0-20008-0-", "") for r in drows]
    xx = np.arange(len(drows))
    for key, se, col, lab, dx in (("b_Tint", "se_Tint", "#1f77b4", "T_int (joint)", -0.15),
                                  ("b_Tlom", "se_Tlom", "#d62728", "T_lom (joint)", 0.15)):
        ax.errorbar(xx + dx, [r.get(key, np.nan) for r in drows], yerr=[r.get(se, np.nan) for r in drows], fmt="o",
                    color=col, label=lab)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xticks(xx)
    ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("daytime same-hour slope [%/K]")
    ax.set_ylim(-10, 12)
    ax.set_title("(b) day: joint slopes on T_int and T_lom (pairs: A-B differences)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q7_tlom.png", dpi=130)
    print(f"-> {FIGS / 'q7_tlom.png'}")


if __name__ == "__main__":
    main()
