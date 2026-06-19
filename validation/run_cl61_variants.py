#!/usr/bin/env python3
"""Controlled re-run of the CL61 Rayleigh calibration in a sensitivity matrix.

Variants (year 2026, the period the L1 archive covers, for a fair comparison):
    L2_WVon  L2_WVoff   -- L2 monthly  (A:/E-PROFILE_L2_monthly)
    L1_WVon  L1_WVoff   -- L1 daily    (D:/E-PROFILE_L1_2026)

Kalman is NOT a variant here: it is a downstream smoothing of the SAME daily
CSV, so "no Kalman" = daily median of the CSV and "Kalman" = the filtered series.

L1 vs L2 instrument identifiers differ (E-PROFILE assigns letters per level), so
each CL61's L1 identifier (the 910.55 nm channel) is resolved separately; the
output is still keyed by the manifest key '<WMO>_<id>' so it lines up with the
cloud calibration. Payerne CL61 (0-20000-0-06610_C) has no network L1 (its only
910 nm L1 channel is the operational CL31) and is skipped in the L1 variants.

Output: C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/cl61_verify/<variant>/<key>/<key>_cl.csv
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import calendar as _cal
import csv
import glob
import json
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from rayleigh_calibration.config import InstrumentType

YEAR = 2026
L2BASE = Path("A:/E-PROFILE_L2_monthly")
L1BASE = Path("D:/E-PROFILE_L1_2026")
OUTROOT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/cl61_verify")
MANIFEST = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json")
OPTS = Path("options.json")
PAYERNE = "0-20000-0-06610_C"   # research CL61, no network L1

VARIANTS = {
    "L2_WVon":  dict(level="L2_MONTHLY", base=L2BASE, wv=True),
    "L2_WVoff": dict(level="L2_MONTHLY", base=L2BASE, wv=False),
    "L1_WVon":  dict(level="L1",         base=L1BASE, wv=True),
    "L1_WVoff": dict(level="L1",         base=L1BASE, wv=False),
}
_FIELDS = ["date", "flag", "lidar_constant", "uncertainty",
           "bottom_height", "top_height", "message"]


def cl61_stations():
    man = json.load(open(MANIFEST))
    return [m for m in man if m["itype"].upper() == "CL61"]


def l1_identifier(wmo):
    """L1 identifier whose global attribute instrument_type == 'CL61' (the type
    attribute is authoritative; a co-located CL51 is also ~910 nm so wavelength
    alone cannot tell them apart)."""
    from netCDF4 import Dataset
    fs = glob.glob(f"{L1BASE}/{wmo}/**/L1_{wmo}_*.nc", recursive=True)
    best = {}
    for f in fs:
        tail = os.path.basename(f).split("_")[-1]      # <id><YYYYMMDD>.nc
        ident = tail.split(str(YEAR))[0]               # everything before the date
        best.setdefault(ident, f)
    for ident, f in sorted(best.items()):
        try:
            d = Dataset(f)
            itype = str(getattr(d, "instrument_type", "")).strip()
            d.close()
        except Exception:
            continue
        if itype.upper() == "CL61":
            return ident
    return None


def l2_dates(wmo, ident):
    out = set()
    for f in glob.glob(f"{L2BASE}/{wmo}/**/L2_{wmo}_{ident}*.nc", recursive=True):
        ym = "".join(c for c in os.path.basename(f).split("_")[-1] if c.isdigit())[-6:]
        if len(ym) == 6 and int(ym[:4]) == YEAR:
            y, mo = int(ym[:4]), int(ym[4:6])
            for day in range(1, _cal.monthrange(y, mo)[1] + 1):
                out.add(f"{y:04d}{mo:02d}{day:02d}")
    return sorted(out)


def l1_dates(wmo, l1id):
    out = set()
    for f in glob.glob(f"{L1BASE}/{wmo}/**/L1_{wmo}_{l1id}*.nc", recursive=True):
        digits = "".join(c for c in os.path.basename(f) if c.isdigit())
        ymd = digits[-8:]
        if len(ymd) == 8 and ymd.startswith(str(YEAR)):
            out.add(ymd)
    return sorted(out)


def _worker(task):
    key, wmo, calib_id, itype, lat, lon, alt, ds, level, base, wv, outdir = task
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    info = InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier=calib_id,
                          instrument_type=InstrumentType.CL61,
                          latitude=lat, longitude=lon, altitude=alt)
    o = CalibrationOptions.from_json(OPTS)
    o.folder_root = Path(base)
    o.data_level = DataLevel.L1 if level == "L1" else DataLevel.L2_MONTHLY
    o.apply_wv_correction = bool(wv)
    o.plot_all = o.plot_main = False
    o.folder_output = Path(outdir)
    try:
        r = calibrate_rayleigh(ds, info, o)
        return key, dict(date=ds, flag=r.flag, lidar_constant=r.lidar_constant,
                         uncertainty=r.uncertainty, bottom_height=r.calibration_bottom_height,
                         top_height=r.calibration_top_height, message=r.message)
    except Exception as exc:
        return key, dict(date=ds, flag=-99, lidar_constant=-1, uncertainty=0,
                         bottom_height=None, top_height=None, message=f"{type(exc).__name__}: {exc}")


def build_tasks(variant):
    v = VARIANTS[variant]
    outbase = OUTROOT / variant
    tasks, plan = [], []
    for s in cl61_stations():
        key = f"{s['wmo']}_{s['identifier']}"
        if v["level"] == "L1":
            if key == PAYERNE:
                continue
            l1id = l1_identifier(s["wmo"])
            if l1id is None:
                print(f"  [{variant}] no 910.55 nm L1 channel for {key}; skip", flush=True)
                continue
            calib_id, dates = l1id, l1_dates(s["wmo"], l1id)
        else:
            calib_id, dates = s["identifier"], l2_dates(s["wmo"], s["identifier"])
        outdir = outbase / key
        if (outdir / f"{key}_cl.csv").exists():
            print(f"  [{variant}] {key} already done", flush=True)
            continue
        plan.append(f"{key} (calib_id={calib_id}, {len(dates)} dates)")
        for ds in dates:
            tasks.append((key, s["wmo"], calib_id, s["itype"], s["lat"], s["lon"],
                          s["alt"], ds, v["level"], str(v["base"]), v["wv"], str(outdir)))
    print(f"[{variant}] {len(plan)} instruments, {len(tasks)} night-calibrations:")
    for p in plan:
        print("   ", p)
    return tasks, outbase


def main():
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True, choices=list(VARIANTS))
    ap.add_argument("--workers", type=int, default=20)
    args = ap.parse_args()

    tasks, outbase = build_tasks(args.variant)
    if not tasks:
        print(f"[{args.variant}] nothing to do"); return
    outbase.mkdir(parents=True, exist_ok=True)

    results = {}
    nw = min(args.workers, len(tasks)) or 1
    with ProcessPoolExecutor(max_workers=nw) as pool:
        futs = [pool.submit(_worker, t) for t in tasks]
        done = 0
        for fut in as_completed(futs):
            key, row = fut.result()
            results.setdefault(key, []).append(row)
            done += 1
            if done % 50 == 0:
                print(f"  [{args.variant}] {done}/{len(tasks)}", flush=True)

    for key, rows in results.items():
        rows.sort(key=lambda r: r["date"])
        sdir = outbase / key
        sdir.mkdir(parents=True, exist_ok=True)
        with open(sdir / f"{key}_cl.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=_FIELDS)
            w.writeheader(); w.writerows(rows)
        cl = [float(r["lidar_constant"]) for r in rows if str(r["flag"]) in ("1", "1.0", "0.5")]
        med = float(np.median(cl)) if cl else float("nan")
        print(f"  [{args.variant}] {key}: {len(cl)}/{len(rows)} nights ok, med={med:.4g}", flush=True)

    print(f"VARIANT_DONE {args.variant}", flush=True)


if __name__ == "__main__":
    main()
