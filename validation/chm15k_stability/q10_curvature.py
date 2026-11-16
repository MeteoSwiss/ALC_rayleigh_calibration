# -*- coding: utf-8 -*-
"""Curvature of the temperature response (report 19, part 10): laser energy vs laser-line/filter mismatch.

A narrowband interference filter whose passband drifts across a fixed laser line gives a transmission
that is a smooth peaked function of temperature -> curvature over a wide temperature range, with the
sign of the local slope set by which flank the line sits on. A pulse-energy change (E ~ 1/sigma) is
close to linear over ~20 K. Per unit-era, on 31-day anomalies around the seasonal state (M2 of q1):

    la = b1 * ta + b2 * (T - T_mean) * ta        (slope varying linearly with the absolute temperature)

and directly on the season-controlled model M1 + quadratic term:

    ln C = a + b T + c (T - T_mean)^2 + annual harmonic + trend

b2 / c significantly negative everywhere -> the slope decreases at high T (approaching a passband peak);
mixed signs -> unit-specific flank positions.

Run:  python validation/chm15k_stability/q10_curvature.py   (after q1)
"""
from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FIGS, OUT  # noqa: E402
import q1_temperature as Q  # noqa: E402

warnings.filterwarnings("ignore")


def curvature(series):
    rows = []
    for key, (t, logc, T, la, ta) in series.items():
        if key[2] != "temp_int" or len(t) < 80:
            continue
        Tm = np.mean(T)
        one = np.ones(len(t))
        doy = 2 * np.pi * (t % 365.25) / 365.25
        yrs = (t - t.mean()) / 365.25
        X = np.column_stack([one, T - Tm, (T - Tm) ** 2, yrs, np.sin(doy), np.cos(doy)])
        coef, se, _ = Q.ols(X, logc)
        ok = np.isfinite(la) & np.isfinite(ta)
        X2 = np.column_stack([np.ones(ok.sum()), ta[ok], (T[ok] - Tm) * ta[ok]])
        c2, s2, _ = Q.ols(X2, la[ok])
        rows.append(dict(stream=key[0], era=key[1], n=len(t), T_range=float(np.percentile(T, 95) - np.percentile(T, 5)),
                         b_lin=100 * coef[1], c_quad=100 * coef[2], se_quad=100 * se[2],
                         slope_anom=100 * c2[1], dslope_dT=100 * c2[2], se_dslope=100 * s2[2]))
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, sa = Q.arm_a()
    _, sb = Q.arm_b()
    out = []
    for arm, s in (("A", sa), ("B", sb)):
        rows = curvature(s)
        for r in rows:
            r["arm"] = arm
        out += rows
        c = np.array([r["c_quad"] for r in rows])
        e = np.array([r["se_quad"] for r in rows])
        d = np.array([r["dslope_dT"] for r in rows])
        de = np.array([r["se_dslope"] for r in rows])
        w = 1 / e ** 2
        wd = 1 / de ** 2
        print(f"arm {arm}: N={len(rows)}  quadratic term c (%/K^2): median {np.median(c):+.3f}, pooled "
              f"{np.sum(w * c) / np.sum(w):+.3f}+/-{1 / np.sqrt(np.sum(w)):.3f}, sig + {np.sum(c > 3 * e)} / - {np.sum(c < -3 * e)}")
        print(f"        anomaly slope change per K of absolute T (%/K^2): median {np.median(d):+.3f}, pooled "
              f"{np.sum(wd * d) / np.sum(wd):+.3f}+/-{1 / np.sqrt(np.sum(wd)):.3f}, sig + {np.sum(d > 3 * de)} / - {np.sum(d < -3 * de)}")
    with open(OUT / "q10_curvature.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    fig, axs = plt.subplots(1, 2, figsize=(17, 5))
    for ax, arm in zip(axs, ("A", "B")):
        rr = [r for r in out if r["arm"] == arm]
        ax.errorbar([r["b_lin"] for r in rr], [r["c_quad"] for r in rr], yerr=[r["se_quad"] for r in rr], fmt="o",
                    ms=3, lw=0.5, alpha=0.7)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("linear slope at the mean temperature [%/K]")
        ax.set_ylabel("quadratic term [%/K$^2$]")
        ax.set_title(f"curvature of ln C vs T$_{{int}}$ (season + trend controlled), arm {arm}", fontsize=10)
        ax.set_ylim(-0.5, 0.5)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "q10_curvature.png", dpi=130)
    print(f"-> {FIGS / 'q10_curvature.png'}")


if __name__ == "__main__":
    main()
