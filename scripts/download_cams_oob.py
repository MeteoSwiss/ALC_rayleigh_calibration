#!/usr/bin/env python3
"""Download per-station LEAN CAMS boxes for out-of-domain E-PROFILE stations.

The regional CAMS archive (Europe/N-Atlantic) does not cover far-flung affiliates
(Canada, New Zealand, ...). A 910 nm ceilometer there cannot do its water-vapor
correction and is flagged -10 ("Closest CAMS data too far"). This script gives each
such station its OWN small CAMS box, downloading only the calibration fields
(t/q/z/lnsp -- NO aerosol backscatter), so the station calibrates but produces no OmB.

It reads ``validation/out_of_domain_stations.json`` (wmo, lat/lon, box [N,W,S,E], and a
start..end YYYYMM month range), then for each station writes one monthly file per month::

    <out>/<wmo>/CAMS_Beta_YYYYMM.nc

The calibration finds these via per-station routing (``calibration/io/cams.py``
``find_cams_file(station_id=...)``): a file under ``cams_folder/<wmo>/`` wins over the
shared archive, and is treated as authoritative (no re-download from the European domain).

Why chunked + parallel: the ADS per-request cost limit (variables x levels x steps x
dates) caps a single request at only ~3 days even for a tiny box, and each request waits
~minutes in the ADS queue. So a month is split into small day-chunks that are fetched
CONCURRENTLY (one ADS job each), then merged into the monthly file. The GRIB messages are
self-describing, so concatenated chunks form one valid multi-message file.

Usage
-----
    python scripts/download_cams_oob.py                       # all stations in the registry
    python scripts/download_cams_oob.py --only 0-20008-0-LAU  # one station
    python scripts/download_cams_oob.py --months 202602       # just these month(s) (comma-sep)
    python scripts/download_cams_oob.py --out D:/CAMS         # CAMS root (per-station subdirs)
    python scripts/download_cams_oob.py --chunk-days 3 --workers 6
    python scripts/download_cams_oob.py --force               # re-download existing monthly files

Resumable: a monthly file that already exists is skipped (unless --force). Needs cdsapi +
cfgrib + ADS credentials (see calibration/io/download_cams_beta.py).
"""
from __future__ import annotations

import argparse
import calendar
import json
import os
import shutil
import sys
import tempfile
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import xarray as xr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration.io import download_cams_beta as D  # noqa: E402

DEFAULT_REGISTRY = REPO / "validation" / "out_of_domain_stations.json"


def month_list(start: str, end: str) -> list[str]:
    """Inclusive list of 'YYYYMM' from start to end."""
    y, m = int(start[:4]), int(start[4:6])
    ey, em = int(end[:4]), int(end[4:6])
    out = []
    while (y, m) <= (ey, em):
        out.append(f"{y:04d}{m:02d}")
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return out


def month_days(ym: str, last_allowed: date) -> list[str]:
    """'YYYY-MM-DD' for every day of the month, capped at last_allowed."""
    y, m = int(ym[:4]), int(ym[4:6])
    n = calendar.monthrange(y, m)[1]
    out = []
    for d in range(1, n + 1):
        day = date(y, m, d)
        if day <= last_allowed:
            out.append(day.strftime("%Y-%m-%d"))
    return out


def _retrieve_chunk(days, variables, area, part_path):
    """One ADS request for *days* -> *part_path* GRIB. On a cost-limit rejection, split
    variables (then days) and recurse, concatenating the sub-parts. A fresh client per
    call keeps this thread-safe (cdsapi wraps a non-thread-safe requests.Session)."""
    client = D._ads_client()
    try:
        client.retrieve(
            D.DATASET,
            {
                "date": list(days),
                "time": D.RUN_TIME,
                "leadtime_hour": D.LEADTIME,
                "type": ["forecast"],
                "variable": list(variables),
                "model_level": D.MODEL_LEVELS,
                "area": list(area),
                "data_format": "grib",
                "download_format": "unarchived",
            },
            part_path,
        )
        return
    except Exception as exc:  # noqa: BLE001
        if not D._request_too_large(exc):
            raise
        if len(variables) > 1:
            mid = len(variables) // 2
            subs = [variables[:mid], variables[mid:]]
            split_axis, items = "vars", subs
        elif len(days) > 1:
            mid = len(days) // 2
            items = [days[:mid], days[mid:]]
            split_axis = "days"
        else:
            raise
        sub_parts = []
        for k, it in enumerate(items):
            sp = f"{part_path}.{split_axis}{k}"
            if split_axis == "vars":
                _retrieve_chunk(days, it, area, sp)
            else:
                _retrieve_chunk(it, variables, area, sp)
            sub_parts.append(sp)
        with open(part_path, "wb") as out:
            for sp in sub_parts:
                with open(sp, "rb") as f:
                    shutil.copyfileobj(f, out)
                os.remove(sp)


def download_month(ym, area, out_path, chunk_days=3, workers=6):
    """Download one month for *area* into *out_path* (CAMS_Beta_YYYYMM.nc) via concurrent
    day-chunks merged to one GRIB, then GRIB->netCDF->build_output (lean t/q/z/lnsp)."""
    yesterday = datetime.now(timezone.utc).date() - timedelta(days=1)   # UTC day convention
    days = month_days(ym, yesterday)
    if not days:
        print(f"  [skip] {ym}: no days on/before {yesterday}")
        return False
    chunks = [days[i:i + chunk_days] for i in range(0, len(days), chunk_days)]
    workdir = tempfile.mkdtemp(prefix=f"cams_oob_{ym}_")
    parts = [os.path.join(workdir, f"chunk_{k:02d}.grib") for k in range(len(chunks))]
    print(f"  [get ] {Path(out_path).name}  {days[0]}..{days[-1]}  "
          f"{len(chunks)} chunks x<={chunk_days}d, {workers} concurrent", flush=True)
    try:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_retrieve_chunk, c, D.CALIBRATION_VARIABLES, area, p): k
                    for k, (c, p) in enumerate(zip(chunks, parts))}
            for fut in as_completed(futs):
                fut.result()  # re-raise any chunk failure
        grib_path = os.path.join(workdir, "month.grib")
        with open(grib_path, "wb") as out:
            for p in parts:
                with open(p, "rb") as f:
                    shutil.copyfileobj(f, out)
        raw_nc = os.path.join(workdir, "raw.nc")
        D.grib_to_netcdf(grib_path, raw_nc)
        with xr.open_dataset(raw_nc) as ds:
            D.build_output(ds, str(out_path))
        return True
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    ap.add_argument("--out", default=os.environ.get("ALC_CAMS_DIR", "D:/CAMS"),
                    help="CAMS root folder (per-station <wmo>/ subdirs are created)")
    ap.add_argument("--only", default=None, help="download just this WMO id")
    ap.add_argument("--months", default=None,
                    help="comma-separated YYYYMM to restrict to (default: each station's range)")
    ap.add_argument("--chunk-days", type=int, default=3, help="days per ADS request")
    ap.add_argument("--workers", type=int, default=6, help="concurrent ADS requests")
    ap.add_argument("--force", action="store_true", help="re-download existing monthly files")
    args = ap.parse_args()

    reg = json.loads(Path(args.registry).read_text())
    stations = reg["stations"]
    if args.only:
        stations = [s for s in stations if s["wmo"] == args.only]
        if not stations:
            sys.exit(f"no station {args.only!r} in {args.registry}")
    only_months = set(args.months.split(",")) if args.months else None

    root = Path(args.out)
    n_done = n_skip = n_fail = 0
    for st in stations:
        wmo = st["wmo"]
        area = [float(x) for x in st["box"]]
        out_dir = root / wmo
        out_dir.mkdir(parents=True, exist_ok=True)
        months = month_list(st["start"], st["end"])
        if only_months:
            months = [m for m in months if m in only_months]
        print(f"\n[oob] {wmo} ({st.get('site', '')}) box(N,W,S,E)={area}  "
              f"{months[0] if months else '-'}..{months[-1] if months else '-'} "
              f"({len(months)} months) -> {out_dir}")
        for ym in months:
            out_path = out_dir / f"CAMS_Beta_{ym}.nc"
            if out_path.exists() and not args.force:
                print(f"  [skip] {out_path.name} (exists)")
                n_skip += 1
                continue
            try:
                ok = download_month(ym, area, out_path,
                                    chunk_days=args.chunk_days, workers=args.workers)
                if ok:
                    print(f"  [ok  ] {out_path.name}  ({out_path.stat().st_size/1e6:.1f} MB)")
                    n_done += 1
            except Exception as exc:  # noqa: BLE001 - keep going; resumable
                print(f"  [FAIL] {ym}: {exc}")
                traceback.print_exc()
                n_fail += 1
    print(f"\n[oob] done: {n_done} downloaded, {n_skip} skipped, {n_fail} failed")


if __name__ == "__main__":
    main()
