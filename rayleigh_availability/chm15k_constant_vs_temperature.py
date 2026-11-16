# -*- coding: utf-8 -*-
"""Payerne CHM15k: does the Rayleigh constant follow the instrument temperature (internal, external,
detector)?  Operational v2.2 constants vs the dark + background corrected re-run, 2025-2026.

For every calibration (dashboard NetCDF: start_time / end_time of the night window), the median of
temp_int, temp_ext and temperature_detector over that window is read from the L1 files. The constant is
regressed on each temperature (log C vs T -> slope in %/K), and the residual monthly scatter after
removing the temperature dependence is compared with the raw scatter.

Run:  python rayleigh_availability/chm15k_constant_vs_temperature.py
"""
from __future__ import annotations
import csv
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
NC = D / "figs_paper_validation/paper_v08/calib_dashboard/nc"
CORR = D / "diag_v22_darkbg_v2/0-20000-0-06610_A/0-20000-0-06610_A_cal.csv"
L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
FIG = REPO / "doc/reports/figs_remote_dark/chm15k_constant_vs_temperature.png"
TVARS = ("temp_int", "temp_ext", "temperature_detector")
_cache = {}


def l1_day(day):
    if day not in _cache:
        f = L1 / f"{day:%Y}" / f"{day:%m}" / f"L1_0-20000-0-06610_A{day:%Y%m%d}.nc"
        if not f.exists():
            _cache[day] = None
        else:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                _cache[day] = (t, {v: np.ma.filled(ds[v][:].astype("f8"), np.nan) - 273.15 for v in TVARS})
    return _cache[day]


def temps(t0, t1):
    """Median of each temperature over [t0, t1] (days since 1970)."""
    out = {v: [] for v in TVARS}
    d0 = datetime(1970, 1, 1) + timedelta(days=float(t0))
    d1 = datetime(1970, 1, 1) + timedelta(days=float(t1))
    day = d0.date()
    while day <= d1.date():
        got = l1_day(day)
        if got:
            t, x = got
            m = (t >= t0) & (t <= t1)
            for v in TVARS:
                out[v].append(x[v][m])
        day += timedelta(days=1)
    return {v: float(np.nanmedian(np.concatenate(out[v]))) if out[v] and np.concatenate(out[v]).size else np.nan
            for v in TVARS}


def load():
    rows = []
    for y in (2025, 2026):
        with netCDF4.Dataset(NC / f"ALC_calibration_0-20000-0-06610_A{y}.nc") as ds:
            m = np.ma.filled(ds["calibration_method"][:], -1) == 0
            t = np.ma.filled(ds["time"][:], np.nan)[m]
            s = np.ma.filled(ds["start_time"][:], np.nan)[m]
            e = np.ma.filled(ds["end_time"][:], np.nan)[m]
            c = np.ma.filled(ds["lidar_constant"][:], np.nan)[m]
        for ti, si, ei, ci in zip(t, s, e, c):
            d = (datetime(1970, 1, 1) + timedelta(days=float(ti))).strftime("%Y-%m-%d")
            rows.append(dict(date=d, c_op=ci, **temps(si, ei)))
    corr = {}
    for r in csv.DictReader(open(CORR, encoding="utf-8")):
        if r["method"] == "rayleigh" and r["flag"] in ("1", "0.5"):
            corr[f"{r['date'][:4]}-{r['date'][4:6]}-{r['date'][6:8]}"] = float(r["cal_value"])
    for r in rows:
        r["c_corr"] = corr.get(r["date"], np.nan)
    return rows


def fit(c, T):
    ok = np.isfinite(c) & np.isfinite(T) & (c > 0)
    y = np.log(c[ok])
    A = np.vstack([T[ok], np.ones(ok.sum())]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    se = np.sqrt(np.sum(resid ** 2) / (ok.sum() - 2) / np.sum((T[ok] - T[ok].mean()) ** 2))
    r = np.corrcoef(T[ok], y)[0, 1]
    return dict(n=int(ok.sum()), slope_pct=100 * coef[0], slope_se=100 * se, r=r,
                std_raw=100 * np.std(y), std_resid=100 * np.std(resid), coef=coef, ok=ok)


def main():
    rows = load()
    T = {v: np.array([r[v] for r in rows]) for v in TVARS}
    res = {}
    for cname in ("c_op", "c_corr"):
        c = np.array([r[cname] for r in rows])
        print(f"== {'operational' if cname == 'c_op' else 'dark + background corrected'} constant")
        for v in TVARS:
            f = fit(c, T[v])
            res[(cname, v)] = f
            print(f"   vs {v:21s}: n={f['n']:3d}  r = {f['r']:+.2f}  slope {f['slope_pct']:+.2f} +/- "
                  f"{f['slope_se']:.2f} %/K   scatter {f['std_raw']:.1f} % -> {f['std_resid']:.1f} % after T removal")
    for v in TVARS:
        x = T[v][np.isfinite(T[v])]
        print(f"   {v}: range {np.percentile(x, 5):.1f} .. {np.percentile(x, 95):.1f} degC (5-95 %)")
    figure(rows, T, res)


def figure(rows, T, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    dates = np.array([datetime.strptime(r["date"], "%Y-%m-%d") for r in rows])
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.5), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    ax = axs[0]
    for cname, col, lab in (("c_op", "#d62728", "operational"), ("c_corr", "#2ca02c", "dark + background corrected")):
        c = np.array([r[cname] for r in rows])
        ax.plot(dates, c / np.nanmedian(c), "o", ms=3.5, color=col, alpha=0.7, label=lab)
    ax.set_ylabel("C$_L$ / median")
    ax2 = ax.twinx()
    ax2.plot(dates, T["temp_int"], "-", color="#ff7f0e", lw=1, alpha=0.8, label="temp_int")
    ax2.set_ylabel("internal temperature during the calibration [\u00b0C]", color="#ff7f0e")
    ax.set_title("Payerne CHM15k nightly Rayleigh constants and internal temperature", fontsize=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    for ax, cname, title in ((axs[1], "c_op", "operational"), (axs[2], "c_corr", "dark + background corrected")):
        c = np.array([r[cname] for r in rows])
        f = res[(cname, "temp_int")]
        x = T["temp_int"]
        ok = f["ok"]
        sc = ax.scatter(x[ok], c[ok] / np.nanmedian(c[ok]), c=[d.month for d in dates[ok]], cmap="twilight",
                        s=14, vmin=1, vmax=12)
        xx = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 50)
        ax.plot(xx, np.exp(f["coef"][0] * xx + f["coef"][1]) / np.nanmedian(c[ok]), "k--", lw=1.2,
                label=f"{f['slope_pct']:+.2f} \u00b1 {f['slope_se']:.2f} %/K, r = {f['r']:+.2f}")
        ax.set_xlabel("temp_int [\u00b0C]")
        ax.set_ylabel("C$_L$ / median")
        ax.set_title(f"{title} constant vs internal temperature", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    plt.colorbar(sc, ax=axs[2], label="month")
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=130)
    print(f"-> {FIG}")


if __name__ == "__main__":
    main()
