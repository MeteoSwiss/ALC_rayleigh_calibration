# -*- coding: utf-8 -*-
"""Network-wide estimate of the operational far-field background bias for all CHM15k and CL61 with local
L1 data: is Payerne exceptional?

Same test as background_bias_sites.py (clear nights, sun < -10 deg, no cloud base; measured minus
modelled molecular P at 12-15 km; expressed as % of the molecular signal in the 2.5-6.5 km Rayleigh fit
range), on every second night of March-June 2026. Lidar constants from the live dashboard (CHM15k
Rayleigh Kalman, CL61 cloud Kalman, Rayleigh if no cloud series). Only the far-field gates are read.

Run:  python rayleigh_availability/background_bias_network.py
Out:  <DATA>/rayleigh_availability/background_bias_network.json + doc/reports/figs_remote_dark/
      background_bias_network.png
"""
from __future__ import annotations
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
from validation.paper import intercompare as IC                      # noqa: E402
from validation.paper import fetch_dashboard_calib as FD              # noqa: E402
from background_bias_sites import sun_elevation, FAR, FIT, WL         # noqa: E402

L1 = Path("D:/E-PROFILE_L1_2026")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/calib_dashboard")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/background_bias_network.json")
FIG = REPO / "doc" / "reports" / "figs_remote_dark" / "background_bias_network.png"
D0, D1, STEP = datetime(2026, 3, 1), datetime(2026, 6, 30), 2


def streams():
    c = json.loads((REPO / "validation" / "scope_l1_2026_census.json").read_text())
    return [r for r in c if r["type"] in ("CHM15k", "CL61")]


def fetch(r):
    key = f"{r['wmo']}_{r['ident']}"
    need = CALIB / f"{key}_{'cloud' if r['type'] == 'CL61' else 'rayleigh'}.csv"
    if need.exists() and need.stat().st_mtime > datetime(2026, 10, 9).timestamp():
        return key, "cached"
    try:
        FD.write_stream(key, CALIB)
        return key, "ok"
    except Exception as e:                                   # station page missing / no Kalman
        return key, f"fail {e!r}"


def kalman(key, itype):
    for m in (("cloud", "rayleigh") if itype == "CL61" else ("rayleigh",)):
        f = CALIB / f"{key}_{m}.csv"
        if f.exists():
            rows = [x for x in csv.DictReader(open(f, encoding="utf-8")) if x["C_kalman"] not in ("", "nan")]
            if rows:
                return m, (np.array([np.datetime64(x["time"]) for x in rows]).astype("datetime64[D]").astype("f8"),
                           np.array([float(x["C_kalman"]) for x in rows]))
    return None, None


def analyse(r):
    key, itype = f"{r['wmo']}_{r['ident']}", r["type"]
    method, kal = kalman(key, itype)
    base = dict(site=r.get("site", ""), stream=key, type=itype, lat=r["lat"], lon=r["lon"], n=0)
    if kal is None:
        return {**base, "skip": "no dashboard constant"}
    kd, kv = kal
    vals, rmax, alt = [], None, None
    day = D0
    while day <= D1:
        f = L1 / r["wmo"] / f"{day:%Y}" / f"{day:%m}" / f"L1_{r['wmo']}_{r['ident']}{day:%Y%m%d}.nc"
        day += timedelta(days=STEP)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                rmax = float(np.nanmax(rng))
                if rmax < FAR[1] - 500:
                    return {**base, "skip": f"range {rmax:.0f} m < far window"}
                far = np.where((rng >= FAR[0]) & (rng <= FAR[1]))[0]
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                x = ds["rcs_0"]
                xf = np.ma.filled((x[:, far[0]:far[-1] + 1] if x.dimensions[0] == "time"
                                   else x[far[0]:far[-1] + 1, :].T).astype("f8"), np.nan)
                cbh = np.ma.filled(ds["cloud_base_height"][:].astype("f8"), np.nan)
                alt = float(np.ravel(ds["station_altitude"][:])[0])
        except Exception:
            continue
        tt = np.array([datetime(1970, 1, 1) + timedelta(days=float(v)) for v in t])
        cb = cbh if cbh.ndim == 2 else cbh[:, None]
        if cb.shape[0] != tt.size:
            cb = cb.T
        m = (sun_elevation(tt, r["lat"], r["lon"]) < -10) & ~np.any(np.isfinite(cb) & (cb > 0), axis=1)
        if m.sum() < 60:
            continue
        rk = rng[far] / 1000.0
        pmeas = float(np.nanmean(xf[m] / rk[None, :] ** 2))
        cl = float(np.interp(np.datetime64(tt[0].date()).astype("datetime64[D]").astype("f8"), kd, kv))
        pexp = float(np.mean(cl * IC._molecular_beta(rng[far], alt, WL[itype]) * 1e-6 / rk ** 2))
        zf = np.arange(FIT[0], FIT[1] + 1, 30.0)
        pfit = float(np.mean(cl * IC._molecular_beta(zf, alt, WL[itype]) * 1e-6 / (zf / 1000) ** 2))
        if np.isfinite(pmeas) and pfit > 0:
            vals.append(100 * (pmeas - pexp) / pfit)
    if len(vals) < 8:
        return {**base, "skip": f"only {len(vals)} clear nights", "n": len(vals)}
    v = np.array(vals)
    return {**base, "n": len(v), "method": method, "pct_med": float(np.median(v)),
            "pct_se": float(1.2533 * np.std(v) / np.sqrt(len(v))),
            "pct_q25": float(np.percentile(v, 25)), "pct_q75": float(np.percentile(v, 75))}


def main():
    ss = streams()
    with ThreadPoolExecutor(max_workers=8) as ex:
        fetched = list(ex.map(fetch, ss))
    print("dashboard fetch:", sum(1 for _, s in fetched if s != "cached" and not s.startswith("fail")), "new,",
          sum(1 for _, s in fetched if s.startswith("fail")), "failed")
    with ProcessPoolExecutor(max_workers=24) as ex:
        res = list(ex.map(analyse, ss))
    OUT.write_text(json.dumps(res, indent=1))
    ok = [x for x in res if "pct_med" in x]
    for t in ("CHM15k", "CL61"):
        v = np.array([x["pct_med"] for x in ok if x["type"] == t])
        pay = [x for x in ok if x["type"] == t and x["stream"].startswith("0-20000-0-06610")]
        rank = (np.sum(v <= pay[0]["pct_med"]) if pay else None)
        print(f"{t}: {len(v)} units; median {np.median(v):+.1f} %, IQR {np.percentile(v, 25):+.1f} .. "
              f"{np.percentile(v, 75):+.1f} %, 5-95 % {np.percentile(v, 5):+.1f} .. {np.percentile(v, 95):+.1f} %"
              + (f"; Payerne {pay[0]['pct_med']:+.1f} % (rank {rank} of {len(v)} from most negative)" if pay else ""))
    skipped = [x for x in res if "skip" in x]
    print(f"skipped {len(skipped)}:", {k: sum(1 for x in skipped if x['skip'].split()[0] == k)
                                       for k in {x['skip'].split()[0] for x in skipped}})
    figure(ok)


def figure(ok):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(18, 6.5), gridspec_kw=dict(width_ratios=[3, 1]))
    col = {"CHM15k": "#d62728", "CL61": "#1f77b4"}
    for t, ax in zip(("CHM15k", "CL61"), axs):
        xs = sorted([x for x in ok if x["type"] == t], key=lambda x: x["pct_med"])
        y = [x["pct_med"] for x in xs]
        e = [x["pct_se"] for x in xs]
        k = np.arange(len(xs))
        ax.errorbar(k, y, yerr=e, fmt="o", ms=3.5, color=col[t], ecolor="0.6", lw=0.8, capsize=0)
        for i, x in enumerate(xs):
            if x["stream"].startswith("0-20000-0-06610"):
                ax.plot(i, x["pct_med"], "*", ms=16, color="k", label="Payerne")
            if any(s in x["stream"] for s in ("0-20000-0-10393", "0-20000-0-03808", "0-380-5-1")):
                ax.plot(i, x["pct_med"], "D", ms=7, mfc="none", mec="k", label=None)
        ax.axhline(0, color="k", lw=0.8)
        ax.axhline(np.median(y), color=col[t], lw=1, ls="--", label=f"network median {np.median(y):+.1f} %")
        ax.set_xlabel(f"{t} units, sorted ({len(xs)})")
        ax.set_ylabel("background bias [% of molecular signal in 2.5-6.5 km]")
        ax.set_title(f"{t}: far-field background bias on clear nights (Mar-Jun 2026)\n"
                     "(diamonds: Lindenberg, Camborne, Aosta)", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        lo, hi = np.percentile(y, 2), np.percentile(y, 98)
        ax.set_ylim(min(lo, -8) - 2, max(hi, 4) + 2)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=130)
    print(f"-> {FIG}")


if __name__ == "__main__":
    main()
