# -*- coding: utf-8 -*-
"""How much of the seasonal cycle of the CHM15k constant does the temperature normalisation remove?
(report 19, part 13)

For every CHM15k (main optical module, 2025-2026 nightly constants), the paper's Kalman prediction model
(annual sinusoid, kalman_improve_alc.fit_model after the global-IQR outlier screen) is fitted to
  * the operational nightly constants,
  * the constants normalised with the unit's intra-night kappa and with its multi-day (M2) kappa,
and the annual amplitude and scatter around the model are compared. The expected removable part is
kappa x (annual amplitude of the calibration-night internal temperature), in phase with C.

Run:  python validation/chm15k_stability/q13_kalman_network.py
"""
from __future__ import annotations

import csv
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from validation.paper.kalman_improve_alc import fit_model, flag_outliers  # noqa: E402
from common import FIGS, OUT, chm15k_streams, load_cache, load_constants, stream_key  # noqa: E402
from q1_temperature import night_hk  # noqa: E402

warnings.filterwarnings("ignore")
EPOCH = datetime(1970, 1, 1)


def annual(times, values):
    """(amplitude % of offset, peak day of year, scatter %) of the Kalman annual-cycle model."""
    v = np.asarray(values, float)
    out = flag_outliers(v)
    t = [x for x, o in zip(times, out) if not o]
    y = v[~out]
    if len(y) < 60:
        return np.nan, np.nan, np.nan
    predict, params, model = fit_model(t, y)
    if model != "annual cycle":
        return np.nan, np.nan, np.nan
    amp, phase, off = params
    res = y - np.array([predict(x) for x in t])
    # peak day of amp*sin(w doy + phase): w doy + phase = pi/2 (amp>0) or 3pi/2 (amp<0)
    w = 2 * np.pi / 365
    peak = ((np.pi / 2 if amp > 0 else 3 * np.pi / 2) - phase) / w % 365
    return 100 * abs(amp) / off, peak, 100 * np.sqrt(np.mean(res ** 2)) / np.median(y)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    k_in = {r["stream"]: float(r["b_T_ctrl"]) for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8"))}
    k_m2 = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int" and r["b2"] not in ("", "nan"):
            k_m2[(r["stream"], r["era"])] = float(r["b2"])
    rows = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = stream_key(wmo, ident)
        cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
        if cons is None or cache is None or key not in k_in:
            continue
        mod = max(set(cons["module"]), key=list(cons["module"]).count)
        T = night_hk(cache, cons["t0"], cons["t1"])[:, list(cache["hk_vars"]).index("temp_int")]
        m = (cons["module"] == mod) & np.isfinite(T) & (cons["c"] > 0)
        if m.sum() < 100:
            continue
        times = [EPOCH + timedelta(days=float(x)) for x in cons["t"][m]]
        c, T = cons["c"][m], T[m]
        tref = float(np.median(T))
        a_op, p_op, s_op = annual(times, c)
        a_in, _, s_in = annual(times, c * np.exp(-k_in[key] / 100 * (T - tref)))
        km2 = k_m2.get((key, mod), np.nan)
        a_m2, _, s_m2 = annual(times, c * np.exp(-km2 / 100 * (T - tref))) if np.isfinite(km2) else (np.nan,) * 3
        # annual cycle of the calibration-night internal temperature (additive, K)
        from validation.paper.kalman_improve_alc import fitfunc_sin
        from scipy.optimize import curve_fit
        doy = np.array([x.timetuple().tm_yday for x in times], float)
        try:
            pT, _ = curve_fit(fitfunc_sin, doy, T, [5, 3 * np.pi / 2, tref], maxfev=20000)
            aT = abs(pT[0])
        except RuntimeError:
            aT = np.nan
        rows.append(dict(stream=key, site=site, module=mod, n=int(m.sum()), kappa_in=k_in[key], kappa_m2=km2,
                         T_amp=aT, amp_op=a_op, peak_op=p_op, scat_op=s_op, amp_in=a_in, scat_in=s_in,
                         amp_m2=a_m2, scat_m2=s_m2, expected_removable=k_in[key] * aT))
    with open(OUT / "q13_kalman_network.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    A = {k: np.array([r[k] for r in rows], float) for k in rows[0] if k not in ("stream", "site", "module")}
    ok = np.isfinite(A["amp_op"]) & np.isfinite(A["amp_in"])
    print(f"units {ok.sum()}: seasonal amplitude median op {np.median(A['amp_op'][ok]):.1f} % -> intra-night kappa "
          f"{np.median(A['amp_in'][ok]):.1f} %, reduced in {np.mean(A['amp_in'][ok] < A['amp_op'][ok]) * 100:.0f} % of units")
    ok2 = ok & np.isfinite(A["amp_m2"])
    print(f"  with multi-day kappa: {np.median(A['amp_m2'][ok2]):.1f} %, reduced in "
          f"{np.mean(A['amp_m2'][ok2] < A['amp_op'][ok2]) * 100:.0f} %")
    print(f"  scatter median op {np.median(A['scat_op'][ok]):.1f} -> {np.median(A['scat_in'][ok]):.1f} (in) / "
          f"{np.nanmedian(A['scat_m2'][ok2]):.1f} (m2)")
    print(f"  night T_int annual amplitude median {np.nanmedian(A['T_amp']):.1f} K; expected removable "
          f"{np.nanmedian(A['expected_removable']):.1f} %")
    frac = 1 - A["amp_m2"] / A["amp_op"]
    print(f"  fraction of the seasonal amplitude removed (M2 kappa): median {np.nanmedian(frac[ok2]) * 100:.0f} %")
    # ranking: strong kappa and large thermal swing
    top = sorted([r for r in rows if np.isfinite(r["amp_m2"]) and np.isfinite(r["amp_op"])],
                 key=lambda r: -(r["kappa_m2"] * r["T_amp"]))[:12]
    for r in top:
        print(f"  {r['stream']:28s} {r['site'][:18]:18s} {r['module']} kappa in {r['kappa_in']:+.2f} m2 {r['kappa_m2']:+.2f}"
              f" T_amp {r['T_amp']:4.1f} K  amp {r['amp_op']:5.1f} -> {r['amp_in']:5.1f} (in) / {r['amp_m2']:5.1f} (m2)"
              f"  scatter {r['scat_op']:4.1f} -> {r['scat_m2']:4.1f}")
    fig, axs = plt.subplots(1, 2, figsize=(16, 5.5))
    ax = axs[0]
    ax.scatter(A["amp_op"][ok2], A["amp_m2"][ok2], c=A["kappa_m2"][ok2], cmap="viridis", vmin=0, vmax=4, s=25)
    lim = [0, np.nanpercentile(A["amp_op"][ok2], 98) * 1.05]
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel("seasonal amplitude, operational [%]")
    ax.set_ylabel("seasonal amplitude, temperature-normalised [%]")
    ax.set_title("Kalman annual-cycle amplitude per unit (colour = κ multi-day)", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1]
    ax.scatter(A["expected_removable"][ok2], (A["amp_op"] - A["amp_m2"])[ok2], s=25, color="#1f77b4")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("κ (intra-night) × night T$_{int}$ annual amplitude [%]")
    ax.set_ylabel("seasonal amplitude removed [%-pt]")
    ax.set_title("removal grows with κ × the unit's seasonal temperature swing", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "q13_kalman_network.png", dpi=130)
    print(f"-> {FIGS / 'q13_kalman_network.png'}")


if __name__ == "__main__":
    main()
