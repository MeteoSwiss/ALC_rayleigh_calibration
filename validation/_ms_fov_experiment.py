"""Tests 1-3 of the multiple-scattering validation plan, on co-located different-FOV
ceilometers at Payerne / Lindenberg / Palaiseau (each with a Cloudnet der product).

For every instrument stream, over 2025-01..2026-05, runs the cloud calibration in three
eta modes (noeta / legacy / pvc) and saves PER IN-CLOUD PROFILE: (time, cloud height,
apparent lidar ratio in each mode). From these:
  - raw apparent etaS = S(noeta)  -> model-independent (Tests 1, 3, and the theoretical fit)
  - calibrated C_L per mode = cc_applied * S_THEO / S(mode)  -> C_L-vs-altitude flatness test
  - concurrent cross-instrument matching by time -> the FOV-ratio killer test (Test 1)

Instrument FOV / wavelength (half-angle receiver FOV):
  CHM15k 0.23 mrad / 1064 nm   (narrowest -> least MS)
  CL61   0.56 mrad / 910 nm
  CL31   0.83 mrad / 910 nm    (widest -> most MS)

Caveats baked in:
  - Lindenberg CL61 had a ~2x hardware step-change at 2025-05: EXCLUDE < 2025-05-01.
  - Payerne CL61 L1 exists only for 2026 (auto).
  - CL61 carries a near-range electronic undershoot (uncorrectable; see memory
    network-offset-models) that biases low-cloud integrals -> CHM15k is the CLEAN reference.
"""
from __future__ import annotations

import copy
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import date, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

L1_ROOT = Path("D:/E-PROFILE_L1_2026")
CAMS = Path("D:/CAMS")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/ms_fov")
START, END = date(2025, 1, 1), date(2026, 5, 31)
WORKERS = 30   # home gaming PC has 32 cores
S_THEO = 18.8
FLAT = np.array([[0.25, 1.0], [4.0, 1.0]])

# stream, instrument type, receiver FOV [mrad], wavelength [nm], optional earliest date
STREAMS = [
    # Payerne
    ("0-20000-0-06610_C", "CL61", 0.56, 910.0, None),
    ("0-20000-0-06610_B", "CL31", 0.83, 910.0, None),
    ("0-20000-0-06610_A", "CHM15k", 0.23, 1064.0, None),
    # Lindenberg (CL61 post step-change only)
    ("0-20000-0-10393_C", "CL61", 0.56, 910.0, "20250501"),
    ("0-20000-0-10393_0", "CHM15k", 0.23, 1064.0, None),
    # Palaiseau / SIRTA
    ("0-250-1001-07151_A", "CL31", 0.83, 910.0, None),
    ("0-250-1001-07151_B", "CHM15k", 0.23, 1064.0, None),
]
SITE_OF_STREAM = {
    "0-20000-0-06610_C": "Payerne", "0-20000-0-06610_B": "Payerne", "0-20000-0-06610_A": "Payerne",
    "0-20000-0-10393_C": "Lindenberg", "0-20000-0-10393_0": "Lindenberg",
    "0-250-1001-07151_A": "Palaiseau", "0-250-1001-07151_B": "Palaiseau",
}
COORD = {"Payerne": (46.81, 6.94), "Lindenberg": (52.21, 14.12), "Palaiseau": (48.72, 2.21)}

# Pristine eta tables captured at import (before any worker mutates the module globals).
import calibration.cloud.calibration as _CCP
import calibration.cloud._filters as _CF  # eta tables live here post-split
_LEG = dict(CL31=_CF._ETA_CL31.copy(), CL51=_CF._ETA_CL51.copy())         # Hopkin ladders
_PVC_CHM = _CF._ETA_CHM15K.copy()
_PVC_MINI = _CF._ETA_MINIMPL.copy()
_PVC_MPL = _CF._ETA_MPL.copy()


def _pvc_cl():
    from validation.multiple_scattering_eta import eta_table, CBH_GRID_KM
    return (np.column_stack([CBH_GRID_KM, eta_table("CL31", 8.0, alpha=10e-3)]),
            np.column_stack([CBH_GRID_KM, eta_table("CL51", 8.0, alpha=10e-3)]))


def _set_mode(CC, mode, pvc31, pvc51):
    if mode == "noeta":
        _CF._ETA_CL31 = FLAT.copy(); _CF._ETA_CL51 = FLAT.copy()
        _CF._ETA_CHM15K = FLAT.copy(); _CF._ETA_MINIMPL = FLAT.copy(); _CF._ETA_MPL = FLAT.copy()
    elif mode == "legacy":
        # Hopkin ladders for the Vaisalas; NO correction for the photon counters (MATLAB otherwise)
        _CF._ETA_CL31 = _LEG["CL31"].copy(); _CF._ETA_CL51 = _LEG["CL51"].copy()
        _CF._ETA_CHM15K = FLAT.copy(); _CF._ETA_MINIMPL = FLAT.copy(); _CF._ETA_MPL = FLAT.copy()
    else:  # pvc
        _CF._ETA_CL31 = np.asarray(pvc31); _CF._ETA_CL51 = np.asarray(pvc51)
        _CF._ETA_CHM15K = _PVC_CHM.copy(); _CF._ETA_MINIMPL = _PVC_MINI.copy(); _CF._ETA_MPL = _PVC_MPL.copy()


def _l1(stream, ds):
    wmo, ident = stream.rsplit("_", 1)
    return L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"


def worker(job):
    stream, itype, ds, pvc31, pvc51 = job
    warnings.filterwarnings("ignore")
    import calibration.cloud.calibration as CC
    from calibration.cloud.calibration import (
        CloudCalConfig, set_defaults, liquid_cloud_calibration_from_data)
    from tests._ceilo_reader import read_ceilometer_data
    fp = _l1(stream, ds)
    site = SITE_OF_STREAM[stream]
    lat, lon = COORD[site]
    if not fp.exists():
        return dict(key=stream, itype=itype, date=ds, empty="nofile")
    try:
        base = set_defaults(CloudCalConfig(
            nc_file=str(fp), instrument=itype, apply_wv_correction=(itype != "CHM15k"),
            apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
            cams_folder=str(CAMS), abs_cs_lookup_table="",
            station_latitude=lat, station_longitude=lon,
            average_time_s=30.0, average_range_m=10.0))
        data0, status = read_ceilometer_data(base.nc_file, base)
        if status != 0 or data0 is None:
            return dict(key=stream, itype=itype, date=ds, empty="nodata")
        out = dict(key=stream, itype=itype, site=site, date=ds)
        cc = None
        for mode in ("noeta", "legacy", "pvc"):
            _set_mode(CC, mode, pvc31, pvc51)
            res = liquid_cloud_calibration_from_data(copy.deepcopy(data0), base)
            n = int(getattr(res, "n_profiles", 0))
            S = np.asarray(getattr(res, "lidar_ratios", []), dtype=float)
            cbh = np.asarray(getattr(res, "cbh", []), dtype=float)
            tt = np.asarray(getattr(res, "time", []))
            k = min(S.size, cbh.size, tt.size)
            ok = (S[:k] > 0) & np.isfinite(S[:k]) & np.isfinite(cbh[:k])
            if mode == "noeta":
                # epoch seconds for cross-instrument concurrent matching
                if tt.size and np.issubdtype(tt.dtype, np.datetime64):
                    ts = tt[:k].astype("datetime64[s]").astype(np.int64)
                else:
                    ts = np.round(np.asarray(tt[:k], dtype=float)).astype(np.int64)
                out["time"] = ts[ok].tolist()
                out["cbh"] = cbh[:k][ok].tolist()
                out["S_noeta"] = S[:k][ok].tolist()
                out["n"] = int(ok.sum())
                cc = getattr(res, "calibration_constant_applied", None)
            else:
                out[f"S_{mode}"] = S[:k][ok].tolist()
        out["cc_applied"] = float(cc) if cc is not None and np.isfinite(cc) else None
        return out
    except Exception as exc:  # noqa: BLE001
        return dict(key=stream, itype=itype, date=ds, err=f"{type(exc).__name__}: {exc}"[:140])


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    pvc31, pvc51 = _pvc_cl()
    np.savez(OUT / "tables.npz", pvc31=pvc31, pvc51=pvc51,
             leg31=_LEG["CL31"], leg51=_LEG["CL51"],
             pvc_chm=_PVC_CHM, pvc_mini=_PVC_MINI)
    days = [(START + timedelta(days=k)).strftime("%Y%m%d")
            for k in range((END - START).days + 1)]
    jobs = []
    for stream, itype, fov, wl, earliest in STREAMS:
        for ds in days:
            if earliest and ds < earliest:
                continue
            jobs.append((stream, itype, ds, pvc31, pvc51))
    res_path = OUT / "results.jsonl"
    done = set()
    if res_path.exists():
        for l in open(res_path, encoding="utf-8"):
            r = json.loads(l)
            done.add((r["key"], r["date"]))
    jobs = [j for j in jobs if (j[0], j[2]) not in done]
    print(f"{len(jobs)} jobs to run ({len(done)} already done)", flush=True)
    n = ok = 0
    with open(res_path, "a", encoding="utf-8") as f:
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for r in ex.map(worker, jobs, chunksize=6):
                n += 1
                if r is not None:
                    f.write(json.dumps(r) + "\n"); ok += 1
                if n % 300 == 0:
                    print(f"  {n}/{len(jobs)} ({ok} written)", flush=True)
    print(f"done: {n} jobs, {ok} rows -> {res_path}", flush=True)


if __name__ == "__main__":
    run()
