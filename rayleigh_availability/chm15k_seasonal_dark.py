# -*- coding: utf-8 -*-
"""Is the seasonal cycle of the Payerne CHM15k Rayleigh constant partly an artefact of the baseline
(dark + background bias), through the seasonal change of the fit-window altitude?

Nightly data, 2025-01 .. 2026-10:
  * operational v2.2 constants with their fit windows and optical-module temperature (dashboard
    NetCDF ALC_calibration_<stream><year>.nc);
  * the same calibration re-run with the measured dark + far-field background corrected
    (diag_v22_darkbg_v2 <stream>_cal.csv, bottom/top heights per night).
Per night the expected baseline bias is f = mean over the fit window of (dark + background) / molecular
(attenuated, 1064 nm, US standard atmosphere), so that C_op ~ C_true * (1 + f).

Run:  python rayleigh_availability/chm15k_seasonal_dark.py
"""
from __future__ import annotations
import csv
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from validation.paper import intercompare as IC   # noqa: E402

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
NC = D / "figs_paper_validation/paper_v08/calib_dashboard/nc"
CORR = D / "diag_v22_darkbg_v2/0-20000-0-06610_A/0-20000-0-06610_A_cal.csv"
DARK = D / "rayleigh_availability/dark_profiles_payerne_v2.npz"
BIAS = D / "rayleigh_availability/background_bias_operational.json"
FIG = REPO / "doc/reports/figs_remote_dark/chm15k_seasonal_dark.png"
ALT = 491.0


def operational():
    out = {}
    for y in (2025, 2026):
        with netCDF4.Dataset(NC / f"ALC_calibration_0-20000-0-06610_A{y}.nc") as ds:
            m = np.ma.filled(ds["calibration_method"][:], -1) == 0
            t = np.ma.filled(ds["time"][:], np.nan)[m]
            c = np.ma.filled(ds["lidar_constant"][:], np.nan)[m]
            b = np.ma.filled(ds["calibration_bottom_height"][:], np.nan)[m] - ALT
            tp = np.ma.filled(ds["calibration_top_height"][:], np.nan)[m] - ALT
            to = np.ma.filled(ds["temperature_optical_module"][:], np.nan)[m]
        for ti, ci, bi, pi, oi in zip(t, c, b, tp, to):
            d = (datetime(1970, 1, 1) + timedelta(days=float(ti))).strftime("%Y-%m-%d")
            out[d] = dict(c=ci, bot=bi, top=pi, temp=oi)
    return out


def corrected():
    out = {}
    for r in csv.DictReader(open(CORR, encoding="utf-8")):
        if r["method"] == "rayleigh" and r["flag"] in ("1", "0.5"):
            d = f"{r['date'][:4]}-{r['date'][4:6]}-{r['date'][6:8]}"
            out[d] = dict(c=float(r["cal_value"]), bot=float(r["bottom_height"]) - ALT,
                          top=float(r["top_height"]) - ALT)
    return out


def baseline_fraction(bot, top, c_true):
    """Mean of (dark + background) / molecular over [bot, top] m a.g.l., for a true constant c_true."""
    z = np.load(DARK)
    r, dark = z["A_range"], z["A_b_rcs"]
    cbg = json.loads(BIAS.read_text())["A"]["bias_med"]
    m = (r >= bot) & (r <= top)
    if m.sum() < 3:
        return np.nan
    mol = IC._molecular_beta(r[m], ALT, 1064.47) * 1e-6 * c_true          # molecular in rcs_0 units
    base = dark[m] + cbg * (r[m] / 1000.0) ** 2
    return float(np.mean(base / mol))


def monthly(dates, vals):
    mon = np.array([d[:7] for d in dates])
    keys = sorted(set(mon))
    return keys, np.array([np.nanmedian(np.asarray(vals)[mon == k]) for k in keys]), \
        np.array([int((mon == k).sum()) for k in keys])


def main():
    op, co = operational(), corrected()
    dd = sorted(op)
    c_op = np.array([op[d]["c"] for d in dd])
    mid = np.array([(op[d]["bot"] + op[d]["top"]) / 2 for d in dd])
    temp = np.array([op[d]["temp"] for d in dd])
    ok = np.isfinite(c_op) & np.isfinite(mid) & np.isfinite(temp)
    print(f"operational nights: {ok.sum()}")
    print(f"  corr(C_op, fit-window mid-height) = {np.corrcoef(c_op[ok], mid[ok])[0, 1]:+.2f}")
    print(f"  corr(C_op, optical-module temperature) = {np.corrcoef(c_op[ok], temp[ok])[0, 1]:+.2f}")
    print(f"  corr(fit-window mid-height, temperature) = {np.corrcoef(mid[ok], temp[ok])[0, 1]:+.2f}")
    # predicted baseline bias at each operational night's window (true constant ~ corrected Kalman level)
    c_true_level = np.nanmedian([v["c"] for v in co.values()])
    f = np.array([baseline_fraction(op[d]["bot"], op[d]["top"], c_true_level) for d in dd])
    c_pred_true = c_op / (1 + f)                         # operational constant with the bias removed
    ok2 = ok & np.isfinite(f)
    print(f"  predicted baseline bias f at the fit window: median {100 * np.nanmedian(f):+.1f} %, "
          f"range (10-90 %) {100 * np.nanpercentile(f, 10):+.1f} .. {100 * np.nanpercentile(f, 90):+.1f} %")
    print(f"  corr(f, fit-window mid-height) = {np.corrcoef(f[ok2], mid[ok2])[0, 1]:+.2f}")
    # seasonal amplitude: monthly medians, 2025-01..2026-09
    res = {}
    for name, dates, vals in (("operational", dd, c_op),
                              ("operational / (1+f) (predicted)", dd, c_pred_true),
                              ("re-run dark + background", sorted(co), [co[d]["c"] for d in sorted(co)])):
        k, v, n = monthly(dates, vals)
        good = n >= 3
        amp = (np.nanmax(v[good]) - np.nanmin(v[good])) / np.nanmedian(v[good])
        cv = np.nanstd(v[good]) / np.nanmedian(v[good])
        res[name] = (k, v, n)
        print(f"  {name:34s}: monthly medians peak-to-peak {100 * amp:5.1f} %, std {100 * cv:4.1f} % "
              f"({good.sum()} months with >= 3 nights)")
    k, v, n = monthly(dd, mid)
    res["mid"] = (k, v, n)
    k, v, n = monthly(dd, temp)
    res["temp"] = (k, v, n)
    figure(res)


def figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(17, 5.5))
    ax = axs[0]
    cols = {"operational": "#d62728", "operational / (1+f) (predicted)": "#9467bd",
            "re-run dark + background": "#2ca02c"}
    for name, col in cols.items():
        k, v, n = res[name]
        good = n >= 3
        x = [datetime.strptime(s, "%Y-%m") for s in np.array(k)[good]]
        ax.plot(x, v[good] / np.nanmedian(v[good]), "o-", color=col, label=name)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_ylabel("monthly median C$_L$ / its own median")
    ax.set_title("Payerne CHM15k Rayleigh constant: seasonal cycle with and without baseline correction",
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    ax = axs[1]
    k, v, n = res["mid"]
    good = n >= 3
    x = [datetime.strptime(s, "%Y-%m") for s in np.array(k)[good]]
    ax.plot(x, v[good] / 1000, "o-", color="#1f77b4", label="fit-window mid-height [km]")
    ax.set_ylabel("fit-window mid-height a.g.l. [km]", color="#1f77b4")
    ax2 = ax.twinx()
    k, v, n = res["temp"]
    ax2.plot(x, v[good] - 273.15 if np.nanmedian(v) > 200 else v[good], "s--", color="#ff7f0e",
             label="optical-module temperature")
    ax2.set_ylabel("optical-module temperature [\u00b0C]", color="#ff7f0e")
    ax.set_title("Fit-window height and optical-module temperature (monthly medians)", fontsize=10)
    ax.grid(alpha=0.3)
    fig.autofmt_xdate()
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=130)
    print(f"-> {FIG}")


if __name__ == "__main__":
    main()
