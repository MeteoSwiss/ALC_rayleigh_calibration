#!/usr/bin/env python3
"""
Full Rayleigh calibration for one representative instrument per type, using the
current (fixed) package code and the `data_level` switch.

  chm15k  : Payerne 0-20000-0-06610   L1          D:/E-PROFILE_L1
  cl61    : 0-203-10-LNG              L2_monthly  A:/E-PROFILE_L2_monthly
  minimpl : Toulouse 0-20000-0-07617  L2_monthly  A:/E-PROFILE_L2_monthly

Parallelises across years (one worker per year). Writes <inst>_cl.csv and a CL
time-series PNG under C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal/<inst>/.

Usage:  python run_calibration.py --inst chm15k|cl61|minimpl [--workers N]
"""
from __future__ import annotations

import argparse
import calendar as _cal
import csv
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal")
A_MONTHLY = Path("A:/E-PROFILE_L2_monthly")

STATIONS = {
    "chm15k":  dict(root=Path("D:/E-PROFILE_L1"), wmo="0-20000-0-06610", ident="A",
                    itype=InstrumentType.CHM15k, level=DataLevel.L1,
                    lat=46.81, lon=6.94, alt=490.0, site="PAYERNE CHM15k"),
    "cl61":    dict(root=A_MONTHLY, wmo="0-203-10-LNG", ident="A",
                    itype=InstrumentType.CL61, level=DataLevel.L2_MONTHLY,
                    lat=48.70, lon=16.93, alt=150.0, site="LNG CL61"),
    "minimpl": dict(root=A_MONTHLY, wmo="0-20000-0-07617", ident="A",
                    itype=InstrumentType.MINI_MPL, level=DataLevel.L2_MONTHLY,
                    lat=43.578, lon=1.374, alt=154.64, site="TOULOUSE Mini-MPL"),
}


def _info(inst):
    s = STATIONS[inst]
    return InstrumentInfo(site_name=s["site"], wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=s["itype"], latitude=s["lat"],
                          longitude=s["lon"], altitude=s["alt"])


def _options(inst):
    o = CalibrationOptions.from_json(Path("options.json"))
    s = STATIONS[inst]
    o.folder_root = s["root"]
    o.data_level = s["level"]
    o.folder_output = OUT / inst
    o.plot_all = o.plot_main = False
    return o


def _discover_dates(inst):
    s = STATIONS[inst]
    d = s["root"] / s["wmo"]
    dates = set()
    if s["level"] == DataLevel.L2_MONTHLY:
        for f in d.glob(f"*/L2_{s['wmo']}_{s['ident']}*.nc"):
            ym = "".join(c for c in f.stem.split("_")[-1] if c.isdigit())[-6:]
            if len(ym) == 6:
                y, mo = int(ym[:4]), int(ym[4:6])
                for day in range(1, _cal.monthrange(y, mo)[1] + 1):
                    dates.add(f"{y:04d}{mo:02d}{day:02d}")
    else:
        prefix = "L1" if s["level"] == DataLevel.L1 else "L2"
        for f in d.glob(f"*/*/{prefix}_{s['wmo']}_{s['ident']}*.nc"):
            digits = "".join(c for c in f.stem.split("_")[-1] if c.isdigit())
            if len(digits) >= 8:
                dates.add(digits[-8:])
    return sorted(dates)


def _process_year(payload):
    year, dates, inst = payload
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    info, options = _info(inst), _options(inst)
    rows = []
    for ds in dates:
        try:
            r = calibrate_rayleigh(ds, info, options)
            rows.append(dict(date=ds, flag=r.flag, flag_meaning=r.flag_meaning,
                             lidar_constant=r.lidar_constant, uncertainty=r.uncertainty,
                             bottom_height=r.calibration_bottom_height,
                             top_height=r.calibration_top_height, message=r.message))
        except Exception as exc:
            rows.append(dict(date=ds, flag=-99, flag_meaning="driver_error",
                             lidar_constant=-1, uncertainty=0, bottom_height=None,
                             top_height=None, message=f"{type(exc).__name__}: {exc}"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inst", required=True, choices=list(STATIONS))
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    inst = args.inst

    dates = _discover_dates(inst)
    by_year = {}
    for d in dates:
        by_year.setdefault(d[:4], []).append(d)
    payloads = [(y, ds, inst) for y, ds in sorted(by_year.items())]
    nw = min(args.workers, len(payloads)) or 1
    print(f"[{inst}] {len(dates)} dates / {len(payloads)} years / {nw} workers "
          f"/ {STATIONS[inst]['level'].value}", flush=True)

    all_rows = []
    with ProcessPoolExecutor(max_workers=nw) as pool:
        futs = {pool.submit(_process_year, p): p[0] for p in payloads}
        for fut in as_completed(futs):
            rows = fut.result(); all_rows.extend(rows)
            n_ok = sum(1 for r in rows if r["flag"] in (1, 0.5))
            print(f"  year {futs[fut]} done: {n_ok}/{len(rows)} ok", flush=True)
    all_rows.sort(key=lambda r: r["date"])

    (OUT / inst).mkdir(parents=True, exist_ok=True)
    csv_path = OUT / inst / f"{inst}_cl.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader(); w.writerows(all_rows)
    ok = [r for r in all_rows if r["flag"] in (1, 0.5)]
    print(f"[{inst}] {len(ok)}/{len(all_rows)} successful -> {csv_path}", flush=True)

    # CL time-series plot
    if ok:
        dts = [datetime.strptime(r["date"], "%Y%m%d") for r in ok]
        cl = np.array([r["lidar_constant"] for r in ok])
        unc = np.array([r["uncertainty"] for r in ok])
        fig, ax = plt.subplots(figsize=(14, 5))
        ax.fill_between(dts, cl - unc, cl + unc, color="tab:blue", alpha=0.15)
        ax.plot(dts, cl, "o-", color="tab:blue", ms=3, lw=0.7)
        ax.set_ylabel("Lidar constant"); ax.grid(True, alpha=0.3)
        ax.set_title(f"{STATIONS[inst]['site']} ({STATIONS[inst]['wmo']}) — Rayleigh lidar "
                     f"constant\n{len(ok)} successful nights, {STATIONS[inst]['level'].value}, "
                     f"median {np.median(cl):.3e}")
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        fig.autofmt_xdate(); fig.tight_layout()
        fig.savefig(OUT / inst / f"{inst}_cl_timeseries.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"[{inst}] median CL = {np.median(cl):.4e}", flush=True)


if __name__ == "__main__":
    main()
