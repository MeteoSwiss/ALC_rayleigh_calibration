# -*- coding: utf-8 -*-
"""Operational far-field background bias of co-located CHM15k + CL61 instruments (Payerne, Lindenberg,
Camborne, Aosta).

On clear nights (sun > 10 deg below the horizon, no cloud base in the profile) the non-range-corrected
signal P = rcs_0 / r_km^2 at 12-15 km should equal the molecular signal there (US standard atmosphere at
the station altitude, two-way molecular transmission, operational lidar constant of the day: CHM15k
Rayleigh, CL61 cloud). bias = measured - expected, per night; reported as the median with its standard
error and as % of the molecular signal in the 2.5-6.5 km Rayleigh fit range (the bias acts as c*r^2 on the
range-corrected signal, i.e. as that fraction of the molecular reference).

Run:  python rayleigh_availability/background_bias_sites.py
"""
from __future__ import annotations
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validation.paper import intercompare as IC   # noqa: E402

L1 = Path("D:/E-PROFILE_L1_2026")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/calib_dashboard")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/background_bias_sites.json")
FAR = (12000.0, 15000.0)
FIT = (2500.0, 6500.0)
WL = {"CHM15k": 1064.47, "CL61": 910.55}
STREAMS = [  # site, wmo, ident, type, calibration method used for the expected molecular signal
    ("Payerne", "0-20000-0-06610", "A", "CHM15k", "rayleigh"),
    ("Payerne", "0-20000-0-06610", "C", "CL61", "cloud"),
    ("Lindenberg", "0-20000-0-10393", "0", "CHM15k", "rayleigh"),
    ("Lindenberg", "0-20000-0-10393", "C", "CL61", "cloud"),
    ("Camborne", "0-20000-0-03808", "A", "CHM15k", "rayleigh"),
    ("Camborne", "0-20000-0-03808", "C", "CL61", "cloud"),
    ("Aosta", "0-380-5-1", "0", "CHM15k", "rayleigh"),
    ("Aosta", "0-380-5-1", "B", "CL61", "cloud"),
]
D0, D1 = datetime(2025, 1, 1), datetime(2026, 10, 8)


def sun_elevation(t, lat, lon):
    """Approximate solar elevation [deg] (NOAA simplified), t in UTC datetimes."""
    doy = np.array([x.timetuple().tm_yday for x in t])
    hr = np.array([x.hour + x.minute / 60 + x.second / 3600 for x in t])
    g = 2 * np.pi / 365 * (doy - 1 + (hr - 12) / 24)
    dec = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
           + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    ha = np.radians((hr * 60 + eqt + 4 * lon) / 4 - 180)
    la = np.radians(lat)
    return np.degrees(np.arcsin(np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)))


def kalman(wmo, ident, method):
    f = CALIB / f"{wmo}_{ident}_{method}.csv"
    rows = [r for r in csv.DictReader(open(f, encoding="utf-8")) if r["C_kalman"] not in ("", "nan")]
    return (np.array([np.datetime64(r["time"]) for r in rows]).astype("datetime64[D]").astype("f8"),
            np.array([float(r["C_kalman"]) for r in rows]))


def analyse(stream):
    site, wmo, ident, itype, method = stream
    kd, kv = kalman(wmo, ident, method)
    out, meta = [], None
    day = D0
    while day <= D1:
        f = L1 / wmo / f"{day:%Y}" / f"{day:%m}" / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
        day += timedelta(days=1)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
                r = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
                x = np.ma.filled(ds["rcs_0"][:].astype("f8"), np.nan)
                cbh = np.ma.filled(ds["cloud_base_height"][:].astype("f8"), np.nan)
                if meta is None:
                    meta = dict(lat=float(np.ravel(ds["station_latitude"][:])[0]),
                                lon=float(np.ravel(ds["station_longitude"][:])[0]),
                                alt=float(np.ravel(ds["station_altitude"][:])[0]), rmax=float(np.nanmax(r)))
        except Exception:
            continue
        if x.shape[0] != t.size:
            x = x.T
        if r.max() < FAR[1] - 500:
            continue
        tt = np.array([datetime(1970, 1, 1) + timedelta(days=float(v)) for v in t])
        dark = sun_elevation(tt, meta["lat"], meta["lon"]) < -10
        cb = cbh if cbh.ndim == 2 else cbh[:, None]
        clear = ~np.any(np.isfinite(cb) & (cb > 0), axis=1)
        m = dark & clear
        if m.sum() < 60:
            continue
        far = (r >= FAR[0]) & (r <= FAR[1])
        rk = r[far] / 1000.0
        pmeas = float(np.nanmean(x[m][:, far] / rk[None, :] ** 2))
        cl = float(np.interp(np.datetime64(tt[0].date()).astype("datetime64[D]").astype("f8"), kd, kv))
        bmol = IC._molecular_beta(r[far], meta["alt"], WL[itype]) * 1e-6
        pexp = float(np.mean(cl * bmol / rk ** 2))
        zf = np.arange(FIT[0], FIT[1] + 1, 30.0)
        pfit = float(np.mean(cl * IC._molecular_beta(zf, meta["alt"], WL[itype]) * 1e-6 / (zf / 1000) ** 2))
        out.append((tt[0].strftime("%Y-%m-%d"), pmeas, pexp, pfit, int(m.sum())))
    if not out:
        return dict(site=site, stream=f"{wmo}_{ident}", type=itype, n=0)
    pm, pe, pf = (np.array([o[i] for o in out]) for i in (1, 2, 3))
    bias = pm - pe
    pct = 100 * bias / pf
    n = len(bias)
    return dict(site=site, stream=f"{wmo}_{ident}", type=itype, n=n, first=out[0][0], last=out[-1][0],
                station_alt=meta["alt"], bias_med=float(np.median(bias)),
                bias_se=float(1.2533 * np.std(bias) / np.sqrt(n)), exp_med=float(np.median(pe)),
                meas_med=float(np.median(pm)), pct_fit_med=float(np.median(pct)),
                pct_fit_se=float(1.2533 * np.std(pct) / np.sqrt(n)),
                pct_fit_q25=float(np.percentile(pct, 25)), pct_fit_q75=float(np.percentile(pct, 75)),
                ratio_meas_to_mol=float(np.median(pm) / np.median(pe)))


def main():
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(analyse, STREAMS))
    OUT.write_text(json.dumps(res, indent=1))
    print(f"{'site':11s} {'instr':7s} {'nights':>6s}  {'period':23s} {'meas/mol far':>12s}  "
          f"{'bias % of fit-range molecular':>30s}")
    for s in res:
        if not s["n"]:
            print(f"{s['site']:11s} {s['type']:7s}      0")
            continue
        print(f"{s['site']:11s} {s['type']:7s} {s['n']:6d}  {s['first']}..{s['last']}  "
              f"{s['ratio_meas_to_mol']:+12.2f}   {s['pct_fit_med']:+6.1f} +/- {s['pct_fit_se']:.1f} % "
              f"(IQR {s['pct_fit_q25']:+.1f} .. {s['pct_fit_q75']:+.1f})")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
