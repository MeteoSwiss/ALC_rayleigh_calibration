# -*- coding: utf-8 -*-
"""Physical model of the CHM15k temperature coefficient: receiver filter drifting against the laser line.

Model (report 19, part 11):
  * The laser line lambda_L is fixed (crystal in the regulated laser-optical module).
  * The receiver interference filter sits at a temperature T_f that follows the housing temperature
    T_int with a first-order lag tau:  dT_f/dt = (T_int - T_f) / tau.
  * Its centre wavelength drifts linearly, lambda_c(T_f) = lambda_c0 + alpha (T_f - T0), alpha > 0
    (dielectric filters shift to longer wavelength when warmed).
  * Passband ~ Gaussian of FWHM W:  Tr = exp(-4 ln2 ((lambda_L - lambda_c)/W)^2).
  Then ln C (proportional to ln Tr) is exactly QUADRATIC in T_f:
      d ln Tr / dT          = 8 ln2 (alpha/W) (Delta/W),   Delta = lambda_L - lambda_c
      (1/2) d2 ln Tr / dT2  = -4 ln2 (alpha/W)^2          (the same for every detuning)
  so from each unit's fitted linear slope b and quadratic term c (q10, season + trend controlled):
      alpha/W = sqrt(-c / (4 ln2)),  Delta0/W = b / (8 ln2 alpha/W),
      transmission loss at the mean temperature = 1 - exp(-4 ln2 (Delta0/W)^2),
      best-match (peak) temperature T* = T_mean + b / (-2 c).
  Only ratios to W are identifiable; W itself needs the filter specification (Lufft).

Test of the lag: intra-night fits of the Rayleigh-window signal against T_f(tau) for tau = 0 ... 12 h;
the best tau minimises the residual variance, and kappa(tau) should rise towards the multi-day value.

Run:  python validation/chm15k_stability/q11_filter_model.py   (after q1, q2a, q10)
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CENSUS, FIGS, OUT, RANGE_MID, chm15k_streams, load_cache, load_constants  # noqa: E402
from q2a_intranight import MIN_CLEAR_FRAC, flatten  # noqa: E402

warnings.filterwarnings("ignore")
LN2 = np.log(2.0)
LOM_SETPOINT = 33.85      # degC, regulated optical-module temperature reported at night
TAUS_H = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]


def curvature_inversion():
    tmean = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int":
            tmean[(r["stream"], r["era"])] = float(r["t_mean"])
    rows = []
    for r in csv.DictReader(open(OUT / "q10_curvature.csv", encoding="utf-8")):
        if r["arm"] != "A":
            continue
        b = float(r["b_lin"]) / 100.0
        c = float(r["c_quad"]) / 100.0
        se = float(r["se_quad"]) / 100.0
        key = (r["stream"], r["era"])
        if key not in tmean or not np.isfinite(c):
            continue
        row = dict(stream=r["stream"], era=r["era"], b=100 * b, c=100 * c, se_c=100 * se, T_mean=tmean[key],
                   c_significant=bool(c < -3 * se))
        if c < 0:
            aw = np.sqrt(-c / (4 * LN2))
            d0 = b / (8 * LN2 * aw)
            row.update(alpha_over_W=aw, delta0_over_W=d0, loss_at_Tmean=1 - np.exp(-4 * LN2 * d0 ** 2),
                       T_peak=tmean[key] + b / (-2 * c))
        rows.append(row)
    return rows


def ema(t_days, x, tau_h):
    """First-order low-pass of x on an irregular time axis (gaps > 6 h restart the filter)."""
    if tau_h <= 0:
        return x.copy()
    y = np.full_like(x, np.nan)
    prev_t, prev_y = None, np.nan
    for i in range(len(x)):
        xi = x[i]
        if not np.isfinite(xi):
            continue
        if prev_t is None or (t_days[i] - prev_t) * 24 > 6 or not np.isfinite(prev_y):
            prev_y = xi
        else:
            a = 1 - np.exp(-(t_days[i] - prev_t) * 24 / tau_h)
            prev_y = prev_y + a * (xi - prev_y)
        prev_t = t_days[i]
        y[i] = prev_y
    return y


def lag_scan(wmo, ident, alt):
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    names = list(cache["hk_vars"])
    fl = flatten(cache)
    T = fl["hk"][:, names.index("temp_int")]
    Tf = {tau: ema(fl["t"], T, tau) for tau in TAUS_H}
    cos_t = np.cos(np.deg2rad(float(cache["tilt"]) if np.isfinite(cache["tilt"]) else 0.0))
    Y, H, TT = [], [], {tau: [] for tau in TAUS_H}
    for k in range(len(cons["t"])):
        i0 = np.searchsorted(fl["t"], cons["t0"][k])
        i1 = np.searchsorted(fl["t"], cons["t1"][k], side="right")
        if i1 - i0 < 24:
            continue
        win = (RANGE_MID >= (cons["zb"][k] - alt) / cos_t) & (RANGE_MID <= (cons["zt"][k] - alt) / cos_t)
        if not win.any():
            continue
        s = np.nanmean(fl["rcs"][i0:i1][:, win], axis=1)
        clear = fl["n_clear"][i0:i1] >= MIN_CLEAR_FRAC * np.maximum(fl["n_all"][i0:i1], 1)
        s = np.where(clear & (s > 0), s, np.nan)
        if np.isfinite(s).sum() < 12:
            continue
        # hourly means to tame photon noise
        nh = (i1 - i0) // 6
        sh = np.array([np.nanmean(s[j * 6:(j + 1) * 6]) for j in range(nh)])
        if np.isfinite(sh).sum() < 4 or np.nanmedian(sh) <= 0:
            continue
        Y.append(np.log(sh / np.nanmedian(sh)))
        tt = np.array([np.mean(fl["t"][i0 + j * 6:i0 + (j + 1) * 6]) for j in range(nh)])
        H.append((tt - tt[0]) * 24)
        for tau in TAUS_H:
            th = np.array([np.nanmean(Tf[tau][i0 + j * 6:i0 + (j + 1) * 6]) for j in range(nh)])
            TT[tau].append(th - np.nanmedian(th))
    if len(Y) < 25:
        return None
    y, h = np.concatenate(Y), np.concatenate(H)
    out = {}
    for tau in TAUS_H:
        x = np.concatenate(TT[tau])
        ok = np.isfinite(y) & np.isfinite(x)
        med = np.median(y[ok])
        ok &= np.abs(y - med) <= 4 * 1.4826 * np.median(np.abs(y[ok] - med))
        X = np.column_stack([np.ones(ok.sum()), x[ok], h[ok]])
        coef, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
        res = y[ok] - X @ coef
        out[tau] = (100 * coef[1], float(np.mean(res ** 2)), int(ok.sum()))
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = curvature_inversion()
    sig = [r for r in rows if r["c_significant"]]
    neg = [r for r in rows if "T_peak" in r]
    with open(OUT / "q11_filter_inversion.csv", "w", newline="", encoding="utf-8") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    def summ(rr, k):
        x = np.array([r[k] for r in rr], float)
        x = x[np.isfinite(x)]
        return np.median(x), np.percentile(x, 25), np.percentile(x, 75), len(x)
    print(f"curvature: {len(rows)} unit-eras, c<0 in {len(neg)}, significantly negative in {len(sig)}")
    for lab, rr in (("all c<0", neg), ("significant c<0", sig)):
        for k in ("alpha_over_W", "delta0_over_W", "loss_at_Tmean", "T_peak", "T_mean"):
            m, lo, hi, n = summ(rr, k)
            print(f"  {lab:16s} {k:14s} median {m:8.4f}  IQR [{lo:8.4f}, {hi:8.4f}]  n {n}")
    # implied filter width for literature drift rates
    aw = summ(sig, "alpha_over_W")[0]
    for alpha in (0.002, 0.005, 0.010, 0.020):
        print(f"  alpha = {1000 * alpha:4.0f} pm/K  ->  implied FWHM W = {alpha / aw:.2f} nm")

    # lag scan on all units
    alts = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    K = {tau: [] for tau in TAUS_H}
    R = {tau: [] for tau in TAUS_H}
    best = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        res = lag_scan(wmo, ident, alts[(wmo, ident)])
        if res is None:
            continue
        v0 = res[0.0][1]
        for tau in TAUS_H:
            K[tau].append(res[tau][0])
            R[tau].append(res[tau][1] / v0)
        best.append(min(TAUS_H, key=lambda tau: res[tau][1]))
    print(f"lag scan on {len(best)} units: best tau median {np.median(best):.1f} h, IQR "
          f"[{np.percentile(best, 25):.1f}, {np.percentile(best, 75):.1f}]")
    for tau in TAUS_H:
        print(f"  tau {tau:4.1f} h: kappa median {np.median(K[tau]):+.2f} %/K, residual variance / tau=0: "
              f"{np.median(R[tau]):.4f}")
    with open(OUT / "q11_lag_scan.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["tau_h", "kappa_median", "kappa_p25", "kappa_p75", "resvar_ratio_median"])
        for tau in TAUS_H:
            w.writerow([tau, np.median(K[tau]), np.percentile(K[tau], 25), np.percentile(K[tau], 75),
                        np.median(R[tau])])

    # figure
    fig, axs = plt.subplots(1, 3, figsize=(20, 5.4))
    ax = axs[0]
    tp = np.array([r["T_peak"] for r in sig])
    tm = np.array([r["T_mean"] for r in sig])
    bins = np.arange(10, 61, 2.5)
    ax.hist(np.clip(tm, 10, 60), bins=bins, histtype="step", lw=2, color="#7f7f7f", label="mean internal T of the unit")
    ax.hist(np.clip(tp, 10, 60), bins=bins, histtype="stepfilled", alpha=0.6, color="#d62728",
            label="best-match (peak) temperature T*")
    ax.axvline(LOM_SETPOINT, color="k", ls="--", lw=1.5, label=f"optical-module set-point {LOM_SETPOINT} °C")
    ax.set_xlabel("temperature [°C]")
    ax.set_ylabel("number of unit-eras (significant curvature)")
    ax.set_title(f"(a) where the filter would match the laser line (N={len(sig)})", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    d0 = np.array([r["delta0_over_W"] for r in neg])
    ax.hist(np.clip(d0, -0.6, 0.6), bins=np.linspace(-0.6, 0.6, 37), color="#1f77b4")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("detuning Δ0 / FWHM at the mean temperature (laser − filter centre)")
    ax.set_ylabel("number of unit-eras")
    loss = np.array([r["loss_at_Tmean"] for r in neg])
    ax.set_title(f"(b) filter detuning: median transmission loss {100 * np.median(loss):.0f} %", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[2]
    km = [np.median(K[t]) for t in TAUS_H]
    rv = [np.median(R[t]) for t in TAUS_H]
    ax.plot(TAUS_H, km, "o-", color="#d62728", lw=2, label="κ (intra-night fit on lagged T)")
    ax.set_xlabel("filter thermal lag τ [h]")
    ax.set_ylabel("κ [%/K]", color="#d62728")
    ax2 = ax.twinx()
    ax2.plot(TAUS_H, rv, "s--", color="#1f77b4", lw=1.5, label="residual variance (rel. to τ = 0)")
    ax2.set_ylabel("residual variance / (τ = 0)", color="#1f77b4")
    ax.set_title(f"(c) thermal lag: best τ median {np.median(best):.1f} h", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "q11_filter_model.png", dpi=130)
    print(f"-> {FIGS / 'q11_filter_model.png'}")


if __name__ == "__main__":
    main()
