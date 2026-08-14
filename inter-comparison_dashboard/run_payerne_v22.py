# -*- coding: utf-8 -*-
"""Rayleigh constants for the Payerne CHM15k and CL61 under eprof_v2 and eprof_v2.2.

The CHM15k comes from the availability corpus, but the Payerne CL61 is not in it (it was used there
only as a cloud-calibrated REFERENCE), and the dashboard needs its Rayleigh series too -- the CL61
carries both a cloud and a Rayleigh calibration, and v2.2 changes the Rayleigh one.

Writes {date: [flag, C_L, unc, bottom, top, message]} per (method, ident), the same record format
the availability study uses, so l1_l2_calib can consume either.

Run:  python inter-comparison_dashboard/run_payerne_v22.py
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel  # noqa: E402
from calibration.config import InstrumentType                                   # noqa: E402

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/calib_raw")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
WMO = "0-20000-0-06610"
STA = dict(latitude=46.813690185546875, longitude=6.942546844482422, altitude=490.0)
# CL61 L1 at Payerne starts 2026-02-24; the CHM15k covers the whole study window.
UNITS = {"C": ("CL61", "20260224", "20260814"), "A": ("CHM15k", "20250101", "20260814")}
METHODS = {"v2.0": ("eprof_v2", {}), "v2.2": ("eprof_v2.2", dict(max_chi2red=2.5))}


def dates(first, last):
    d, d1 = datetime.strptime(first, "%Y%m%d"), datetime.strptime(last, "%Y%m%d")
    while d <= d1:
        yield d.strftime("%Y%m%d")
        d += timedelta(days=1)


def _f(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def run_one(job):
    tag, ident = job
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    itype, first, last = UNITS[ident]
    meth, params = METHODS[tag]
    out_f = OUT / f"payerne_{ident}_{tag}.json"
    rec = json.loads(out_f.read_text()) if out_f.exists() and "--force" not in sys.argv else {}
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.molecular_method = meth
    o.molecular_params = dict(params)
    o.plot_main = o.plot_all = False
    o.folder_output = OUT / "tmp"
    info = InstrumentInfo(site_name=f"PAYERNE_{itype}_{ident}", wmo_id=WMO, identifier=ident,
                          instrument_type=InstrumentType(itype), **STA)
    for ds in dates(first, last):
        if ds in rec:
            continue
        if not (L1_ROOT / WMO / ds[:4] / ds[4:6] / f"L1_{WMO}_{ident}{ds}.nc").exists():
            continue
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception as exc:
            rec[ds] = [-99, None, None, None, None, f"EXC {type(exc).__name__}"[:80]]
            continue
        if r is None:
            continue
        rec[ds] = [_f(r.flag), _f(r.lidar_constant), _f(r.uncertainty),
                   _f(r.calibration_bottom_height), _f(r.calibration_top_height),
                   str(r.message or "")[:100]]
    OUT.mkdir(parents=True, exist_ok=True)
    out_f.write_text(json.dumps(rec), encoding="utf-8")
    n_ok = sum(1 for v in rec.values() if v[0] in (1.0, 0.5))
    return tag, ident, itype, len(rec), n_ok


def main():
    jobs = [(t, i) for t in METHODS for i in UNITS]
    workers = int(os.environ.get("RA_WORKERS", "4"))
    print(f"Payerne Rayleigh: {len(jobs)} jobs over {workers} workers -> {OUT}", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(run_one, j) for j in jobs]):
            tag, ident, itype, n, n_ok = fut.result()
            print(f"  {tag:5s} {ident} ({itype:6s}) {n:4d} nights, {n_ok:3d} valid", flush=True)
    print("PAYERNE_V22_DONE", flush=True)


if __name__ == "__main__":
    main()
