#!/usr/bin/env python3
"""
Calibrate ALL CHM15k / CL61 / Mini-MPL stations that have data in A:/E-PROFILE_L2_monthly,
reading directly from L2-monthly (data_level=L2_monthly) with the current (fixed) code.

Resumable: a station whose CSV already exists is skipped, so the job can be re-launched
after an interruption and simply continues. Per-station results go to
  C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all/<WMO>/<WMO>_cl.csv
and a rolling summary is appended to fullcal_all/summary.csv.

Usage:  python run_all_l2monthly.py [--workers 6] [--types CHM15k,CL61,Mini-MPL]
"""
from __future__ import annotations

# Keep CPU low: force single-threaded BLAS/numpy BEFORE numpy is imported (applies to
# the main process and to every spawned worker, which re-import this module on Windows).
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import calendar as _cal
import csv
import json
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType

BASE = Path("A:/E-PROFILE_L2_monthly")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all")
ITYPE = {"CHM15k": InstrumentType.CHM15k, "CL61": InstrumentType.CL61, "Mini-MPL": InstrumentType.MINI_MPL}


def _key(station):
    """Per-instrument key '<WMO>_<identifier>' (a WMO may host several instruments)."""
    return f"{station['wmo']}_{station['identifier']}"


def _dates_for(station):
    d = BASE / station["wmo"]
    dates = set()
    for f in d.glob(f"*/L2_{station['wmo']}_{station['identifier']}*.nc"):
        ym = "".join(c for c in f.stem.split("_")[-1] if c.isdigit())[-6:]
        if len(ym) == 6:
            y, mo = int(ym[:4]), int(ym[4:6])
            for day in range(1, _cal.monthrange(y, mo)[1] + 1):
                dates.add(f"{y:04d}{mo:02d}{day:02d}")
    return sorted(dates)


_FIELDS = ["date", "flag", "lidar_constant", "uncertainty",
           "bottom_height", "top_height", "message"]


def _partial_path(station, year):
    """Per-year checkpoint CSV, so a process killed mid-station resumes by year."""
    return OUT / _key(station) / "_partial" / f"{year}.csv"


def _process_year(payload):
    year, dates, station = payload
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    info = InstrumentInfo(site_name=station["wmo"], wmo_id=station["wmo"],
                          identifier=station["identifier"], instrument_type=ITYPE[station["itype"]],
                          latitude=station["lat"], longitude=station["lon"], altitude=station["alt"])
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = BASE
    o.data_level = DataLevel.L2_MONTHLY
    o.folder_output = OUT / _key(station)
    o.plot_all = o.plot_main = False
    rows = []
    for ds in dates:
        try:
            r = calibrate_rayleigh(ds, info, o)
            rows.append(dict(date=ds, flag=r.flag, lidar_constant=r.lidar_constant,
                             uncertainty=r.uncertainty, bottom_height=r.calibration_bottom_height,
                             top_height=r.calibration_top_height, message=r.message))
        except Exception as exc:
            rows.append(dict(date=ds, flag=-99, lidar_constant=-1, uncertainty=0,
                             bottom_height=None, top_height=None, message=f"{type(exc).__name__}: {exc}"))
    # Checkpoint this year atomically (write temp, then replace) so a kill never
    # leaves a half-written partial that would be trusted on resume.
    pp = _partial_path(station, year)
    pp.parent.mkdir(parents=True, exist_ok=True)
    tmp = pp.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_FIELDS)
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, pp)
    return year


def _calibrate_station(station, workers):
    dates = _dates_for(station)
    by_year = {}
    for d in dates:
        by_year.setdefault(d[:4], []).append(d)
    # Skip years already checkpointed (resume after a mid-station kill).
    todo_years = [(y, ds, station) for y, ds in sorted(by_year.items())
                  if not _partial_path(station, y).exists()]
    if todo_years:
        nw = min(workers, len(todo_years)) or 1
        with ProcessPoolExecutor(max_workers=nw) as pool:
            futs = [pool.submit(_process_year, p) for p in todo_years]
            for fut in as_completed(futs):
                fut.result()
    # Assemble the full station record from all year partials.
    all_rows = []
    for y in sorted(by_year):
        pp = _partial_path(station, y)
        with open(pp, newline="", encoding="utf-8") as f:
            all_rows.extend(csv.DictReader(f))
    all_rows.sort(key=lambda r: r["date"])
    return all_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--types", default="CHM15k,CL61,Mini-MPL")
    args = ap.parse_args()
    wanted = set(args.types.split(","))

    manifest = json.load(open("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json"))
    stations = [s for s in manifest if s["itype"] in wanted]
    OUT.mkdir(parents=True, exist_ok=True)
    summary_path = OUT / "summary.csv"
    if not summary_path.exists():
        with open(summary_path, "w", newline="") as f:
            csv.writer(f).writerow(["key", "wmo", "identifier", "itype", "n_dates", "n_success", "median_cl"])

    # Output is keyed per instrument '<WMO>_<identifier>' (a WMO may host several).
    done = {p.parent.name for p in OUT.glob("*/*_cl.csv")}
    # Smallest instruments first: many quick ones complete per ~30-min process lifetime
    # (banking progress), and each large multi-year one gets a fresh cycle to itself.
    todo = sorted([s for s in stations if _key(s) not in done],
                  key=lambda s: s.get("n_months", 0))
    print(f"{len(stations)} target instruments; {len(done)} already done; {len(todo)} to do", flush=True)

    for k, s in enumerate(todo, 1):
        key = _key(s)
        rows = _calibrate_station(s, args.workers)
        # rows come back from the year-partial CSVs as strings -> parse for the summary.
        cl = [float(r["lidar_constant"]) for r in rows
              if str(r["flag"]) in ("1", "1.0", "0.5")]
        sdir = OUT / key
        sdir.mkdir(parents=True, exist_ok=True)
        with open(sdir / f"{key}_cl.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=_FIELDS)
            w.writeheader(); w.writerows(rows)
        # Instrument is fully banked; drop its per-year checkpoints.
        for pp in (sdir / "_partial").glob("*.csv"):
            pp.unlink()
        med = float(np.median(cl)) if cl else float("nan")
        with open(summary_path, "a", newline="") as f:
            csv.writer(f).writerow([key, s["wmo"], s["identifier"], s["itype"],
                                    len(rows), len(cl), f"{med:.4e}"])
        print(f"[{k}/{len(todo)}] {key} ({s['itype']}): {len(cl)}/{len(rows)} ok, med={med:.3e}", flush=True)

    print("ALL STATIONS DONE", flush=True)
    Path(OUT / "ALL_DONE.flag").write_text("done")


if __name__ == "__main__":
    main()
