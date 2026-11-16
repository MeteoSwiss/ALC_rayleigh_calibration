# -*- coding: utf-8 -*-
"""Can aerosol explain the CHM15k temperature coefficient? AERONET control (report 19, part 5).

Warm spells are often hazier. Aerosol leaking into the "molecular" Rayleigh window inflates the
retrieved constant, so a warm-hazy covariance would masquerade as a temperature coefficient. The
cache-based proxy (low-level aerosol contrast) already left the slope unchanged; here the
independent, column-integrated AERONET AOD is used. It also sees elevated layers (dust, smoke)
that a boundary-layer proxy misses.

For every CHM15k stream with an AERONET site within 15 km (Level 1.5, daily averages, 2025-2026):
  night AOD = mean of the daily AOD on the evening date and the morning date of the calibration night;
  AOD at 1064 nm from AOD_1020 and the 440-870 Angstrom exponent.
Models, on 31-day anomalies (M2 of q1):
  la ~ ta                         (reference, AERONET-matched nights only)
  la ~ ta + aod_a                 (aerosol load)
  la ~ ta + aod_a + ang_a         (+ particle size)
and M2 restricted to the cleanest third of nights (AOD1064 below the unit's 33rd percentile).

Run:  python validation/chm15k_stability/q5_aeronet.py --aeronet-dir <dir with AOD15_daily_<site>.txt and aeronet_match.json>
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BLOCKS_PER_DAY, FIGS, OUT, load_cache, load_constants  # noqa: E402
from q1_temperature import ANOM_WINDOW, night_hk, robust_mask, running_median  # noqa: E402

warnings.filterwarnings("ignore")
EPOCH = datetime(1970, 1, 1)


def read_aeronet_daily(path):
    """{date: (aod1064, aod500, angstrom_440_870)} from an AERONET v3 daily-average download."""
    out = {}
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    hdr = next((i for i, l in enumerate(lines) if l.startswith("AERONET_Site,")), None)
    if hdr is None:
        return out
    cols = lines[hdr].split(",")
    ix = {c: k for k, c in enumerate(cols)}
    for l in lines[hdr + 1:]:
        p = l.split(",")
        if len(p) < len(cols) - 2:
            continue
        try:
            d = datetime.strptime(p[ix["Date(dd:mm:yyyy)"]], "%d:%m:%Y").date()
            a1020 = float(p[ix["AOD_1020nm"]])
            a500 = float(p[ix["AOD_500nm"]])
            ang = float(p[ix["440-870_Angstrom_Exponent"]])
        except (ValueError, KeyError):
            continue
        if a1020 < -900 or ang < -900:
            continue
        out[d] = (a1020 * (1020.0 / 1064.0) ** ang, a500 if a500 > -900 else np.nan, ang)
    return out


def ols(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    cov = (res @ res / max(len(y) - X.shape[1], 1)) * np.linalg.inv(X.T @ X)
    return coef, np.sqrt(np.diag(cov))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aeronet-dir", required=True)
    a = ap.parse_args()
    d = Path(a.aeronet_dir)
    match = json.load(open(d.parent / "aeronet_match.json", encoding="utf-8"))
    rows, pooled = [], {"ta": [], "aod": [], "la": []}
    for key, site, ae_site, dist in match:
        f = d / f"AOD15_daily_{ae_site}.txt"
        if not f.exists():
            continue
        aer = read_aeronet_daily(f)
        if len(aer) < 50:
            continue
        wmo, ident = key.rsplit("_", 1)
        cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
        if cons is None or cache is None:
            continue
        T = night_hk(cache, cons["t0"], cons["t1"])[:, list(cache["hk_vars"]).index("temp_int")]
        logc = np.log(np.where(cons["c"] > 0, cons["c"], np.nan))
        main_mod = max(set(cons["module"]), key=list(cons["module"]).count)
        m = (cons["module"] == main_mod) & np.isfinite(T) & np.isfinite(logc)
        m[m] = robust_mask(logc[m])
        t = cons["t"][m]
        lc, tt = logc[m], T[m]
        aod = np.full(len(t), np.nan)
        ang = np.full(len(t), np.nan)
        for k, (t0, t1) in enumerate(zip(cons["t0"][m], cons["t1"][m])):
            vals = [aer[x] for x in {(EPOCH + timedelta(days=float(t0))).date(),
                                     (EPOCH + timedelta(days=float(t1))).date()} if x in aer]
            if vals:
                aod[k] = np.mean([v[0] for v in vals])
                ang[k] = np.mean([v[2] for v in vals])
        order = np.argsort(t)
        t, lc, tt, aod, ang = t[order], lc[order], tt[order], aod[order], ang[order]
        la = lc - running_median(t, lc, ANOM_WINDOW)
        ta = tt - running_median(t, tt, ANOM_WINDOW)
        aa = aod - running_median(t, aod, ANOM_WINDOW)
        ga = ang - running_median(t, ang, ANOM_WINDOW)
        ok = np.isfinite(la) & np.isfinite(ta) & np.isfinite(aa) & np.isfinite(ga)
        if ok.sum() < 30 or np.std(ta[ok]) < 0.3:
            continue
        one = np.ones(ok.sum())
        c0, s0 = ols(np.column_stack([one, ta[ok]]), la[ok])
        c1, s1 = ols(np.column_stack([one, ta[ok], aa[ok]]), la[ok])
        c2, s2 = ols(np.column_stack([one, ta[ok], aa[ok], ga[ok]]), la[ok])
        clean = ok & (aod <= np.nanpercentile(aod[ok], 33))
        if clean.sum() >= 15 and np.std(ta[clean]) > 0.3:
            cc, sc = ols(np.column_stack([np.ones(clean.sum()), ta[clean]]), la[clean])
        else:
            cc, sc = (np.nan, np.nan), (np.nan, np.nan)
        row = dict(stream=key, site=site, aeronet=ae_site, dist_km=dist, n=int(ok.sum()),
                   aod_median=float(np.nanmedian(aod[ok])), r_ta_aod=float(np.corrcoef(ta[ok], aa[ok])[0, 1]),
                   b_ref=100 * c0[1], se_ref=100 * s0[1], b_aod=100 * c1[1], se_aod=100 * s1[1],
                   g_aod_per_0p1=100 * c1[2] * 0.1, b_aod_ang=100 * c2[1], se_aod_ang=100 * s2[1],
                   n_clean=int(clean.sum()), b_clean=100 * cc[1], se_clean=100 * sc[1])
        rows.append(row)
        for k_, v in (("ta", ta[ok]), ("aod", aa[ok]), ("la", la[ok])):
            pooled[k_].append(v)
        print(f"  {key:26s} {ae_site:20s} n {row['n']:4d} AOD {row['aod_median']:.3f}  r(dT,dAOD) {row['r_ta_aod']:+.2f}  "
              f"kappa {row['b_ref']:+.2f}+/-{row['se_ref']:.2f} -> +AOD {row['b_aod']:+.2f} -> +AOD+Ang {row['b_aod_ang']:+.2f}"
              f"+/-{row['se_aod_ang']:.2f}  clean-third {row['b_clean']:+.2f}+/-{row['se_clean']:.2f} (n {row['n_clean']})  "
              f"dlnC per +0.1 AOD {row['g_aod_per_0p1']:+.1f} %", flush=True)
    with open(OUT / "q5_aeronet.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    def pool(k, s):
        b = np.array([r[k] for r in rows], float)
        e = np.array([r[s] for r in rows], float)
        ok = np.isfinite(b) & np.isfinite(e) & (e > 0)
        wts = 1 / e[ok] ** 2
        return np.sum(wts * b[ok]) / np.sum(wts), 1 / np.sqrt(np.sum(wts)), np.median(b[ok]), ok.sum()
    print("\nPooled (inverse-variance) / median over units:")
    for k, s, lab in (("b_ref", "se_ref", "M2, AERONET-matched nights"), ("b_aod", "se_aod", "M2 + AOD"),
                      ("b_aod_ang", "se_aod_ang", "M2 + AOD + Angstrom"), ("b_clean", "se_clean", "M2, cleanest third")):
        p, e, md, n = pool(k, s)
        print(f"  {lab:28s} pooled {p:+.2f} +/- {e:.2f} %/K, median {md:+.2f}, N={n}")
    r = np.array([x["r_ta_aod"] for x in rows])
    print(f"  r(T anomaly, AOD anomaly): median {np.median(r):+.2f}  [p10 {np.percentile(r, 10):+.2f}, p90 {np.percentile(r, 90):+.2f}]")
    g = np.array([x["g_aod_per_0p1"] for x in rows])
    print(f"  C response to +0.1 AOD1064: median {np.median(g):+.1f} %")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.5))
    ax = axs[0]
    b0 = np.array([x["b_ref"] for x in rows])
    b2 = np.array([x["b_aod_ang"] for x in rows])
    e2 = np.array([x["se_aod_ang"] for x in rows])
    ax.errorbar(b0, b2, yerr=e2, fmt="o", ms=4, lw=0.6)
    ax.plot([-4, 8], [-4, 8], "k--", lw=0.8, label="1:1")
    ax.set_xlabel("M2 slope, no aerosol term [%/K]")
    ax.set_ylabel("M2 slope with AERONET AOD + Angstrom [%/K]")
    ax.set_title(f"(a) temperature slope with / without AERONET control (N={len(rows)})", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    bc = np.array([x["b_clean"] for x in rows])
    ec = np.array([x["se_clean"] for x in rows])
    ax.errorbar(b0, bc, yerr=ec, fmt="o", ms=4, lw=0.6, color="#2ca02c")
    ax.plot([-4, 8], [-4, 8], "k--", lw=0.8)
    ax.set_xlabel("M2 slope, all nights [%/K]")
    ax.set_ylabel("M2 slope, cleanest third of nights [%/K]")
    ax.set_title("(b) slope on the cleanest (lowest-AOD) nights", fontsize=10)
    ax.set_ylim(-10, 15)
    ax.grid(alpha=0.3)
    ax = axs[2]
    ta, aa = np.concatenate(pooled["ta"]), np.concatenate(pooled["aod"])
    ax.plot(ta, aa, ".", ms=2, alpha=0.3)
    ax.set_xlabel("T$_{int}$ anomaly [K]")
    ax.set_ylabel("AOD$_{1064}$ anomaly (AERONET)")
    ax.set_title(f"(c) warm spells vs aerosol: pooled r = {np.corrcoef(ta, aa)[0, 1]:+.2f}", fontsize=10)
    ax.set_ylim(-0.15, 0.3)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q5_aeronet.png", dpi=130)
    print(f"-> {FIGS / 'q5_aeronet.png'}")


if __name__ == "__main__":
    main()
