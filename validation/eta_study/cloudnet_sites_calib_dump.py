"""Per-profile cloud-calibration scenes for the Vaisala units next to Cloudnet `der` sites without a dump.

The step-1d scene lists (old/ms_fov/results.jsonl, cloud_profile_dump*/<key>_profiles.npz) cover only
Payerne B/C, Palaiseau A and Lindenberg C. The other E-PROFILE Vaisala ceilometers at (or near) a
Cloudnet site with the `der` product have no per-profile list, only per-day results
(rayleigh_availability/cbh_native/cbh_scenes.csv). This driver produces that list for them by running
rayleigh_availability/cloud_profile_dump.process_day (the operational liquid-cloud calibration of
scripts/run_network_calibration.py::_do_cloud, replicated verbatim there) on the unit's known
cloud-calibration days that also have a Cloudnet der file.

Only process_day is used: cloud_profile_dump.main() would also redraw a committed report figure.
CAMS: the module's CAMS_FOLDERS (0.4 deg monthlies D:/CAMS_Monthly_04, then the dailies).
ALC_ETA_TABLES=pvc55 is set here (inherited by the workers): the eta tables of the existing dumps
(2026-07 generation), so the accepted-profile selection is made the same way as for the units of
step 1d. The relative temporal-consistency filter makes the selection nearly independent of eta.

Output (same fields as cloud_profile_dump): C:/DATA/Projects/202611_eta_study/cloudnet_sites/
calib_dump/<key>_profiles.npz with time_ns, date, cbh_m (m AGL), c_oconnor, s_apparent, meta.

Run:  python validation/eta_study/cloudnet_sites_calib_dump.py [--workers 8]
Resumable: a unit whose npz exists is skipped (delete it to redo).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
os.environ["ALC_ETA_TABLES"] = "pvc55"

import argparse
import csv
import json
import subprocess
import sys
import time as time_mod
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import cloud_profile_dump as CPD  # noqa: E402  (also puts REPO on sys.path)

CLOUDNET = Path("D:/Cloudnet_Droplet_radius")
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
SCENES_CSV = REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv"
OUT = Path("C:/DATA/Projects/202611_eta_study/cloudnet_sites/calib_dump")
FILELIST = OUT.parent / "der_filelist.csv"

# unit key -> (Cloudnet site, distance unit-site in km): co-located (<= 1 km) or nearby (<= 11 km)
UNITS = {
    "0-826-300-3_A": ("chilbolton", 0.1),       # CL51 Chilbolton
    "0-20008-0-CEO_A": ("laquila", 0.1),        # CL51 L'Aquila
    "0-20008-0-NOR_A": ("norunda", 0.1),        # CL51 Norunda
    "0-20000-0-03749_A": ("chilbolton", 9.2),   # CL31 Middle Wallop (nearby)
    "0-20000-0-17601_A": ("limassol", 10.9),    # CL31 Akrotiri (nearby)
}


def der_days(site: str) -> set[str]:
    """YYYYMMDD of every der file of a Cloudnet site: the portal list written by
    cloudnet_der_download.py (so the dump can run while the download is still going), plus the disk."""
    days = set()
    if FILELIST.exists():
        with open(FILELIST, newline="", encoding="utf-8") as f:
            days = {r["date"].replace("-", "") for r in csv.DictReader(f) if r["site"] == site}
    days |= {p.name[:8] for p in (CLOUDNET / site).glob(f"*_{site}_der.nc")}
    days |= {p.name[:8] for p in CLOUDNET.glob(f"*_{site}_der.nc")}
    return days


def run_day(task):
    """Worker entry point (ALC_ETA_TABLES above is re-applied when the worker imports this module)."""
    return CPD.process_day(task)


def code_state() -> str:
    """git HEAD + whether the calibration code differs from it (another session may be editing)."""
    try:
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=30).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "calibration"],
                               capture_output=True, text=True, timeout=30).stdout.split("\n")
        dirty = [d[3:] for d in dirty if d.strip()]
        return f"{head}" + (f" + uncommitted: {', '.join(dirty)}" if dirty else "")
    except Exception as exc:  # noqa: BLE001
        return f"unknown ({type(exc).__name__})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--units", default=",".join(UNITS))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    census = {f"{c['wmo']}_{c['ident']}": c for c in json.load(open(CENSUS, encoding="utf-8"))}
    cal_days = {}
    with open(SCENES_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cal_days.setdefault(r["key"], set()).add(r["date"])

    tasks = []
    for key in a.units.split(","):
        if (OUT / f"{key}_profiles.npz").exists():
            print(f"{key}: npz exists - skipped", flush=True)
            continue
        c = census[key]
        days = sorted(cal_days.get(key, set()) & der_days(UNITS[key][0]))
        print(f"{key} {c['type']} {c['site']}: {len(cal_days.get(key, ()))} cloud-calibration days, "
              f"{len(days)} with a Cloudnet der file ({UNITS[key][0]})", flush=True)
        wmo, ident = key.rsplit("_", 1)
        tasks += [(wmo, ident, c["type"], c["lat"], c["lon"], ds) for ds in days]
    if not tasks:
        return
    state = code_state()
    print(f"{len(tasks)} unit-days on {a.workers} workers; ALC_ETA_TABLES={os.environ['ALC_ETA_TABLES']}; "
          f"CAMS {CPD.CAMS_FOLDERS}; code {state}", flush=True)

    acc, fails = {}, {}
    t0 = time_mod.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = [pool.submit(run_day, t) for t in tasks]
        for k, fut in enumerate(as_completed(futs), 1):
            key, ds, payload, msg = fut.result()
            if payload is None:
                fails.setdefault(key, {})
                fails[key][msg] = fails[key].get(msg, 0) + 1
            else:
                acc.setdefault(key, []).append(payload)
            if k % 100 == 0 or k == len(tasks):
                print(f"  {k}/{len(tasks)} ({time_mod.time() - t0:.0f} s)", flush=True)
    for key in {t[0] + "_" + t[1] for t in tasks}:
        payloads = acc.get(key, [])
        if not payloads:
            print(f"{key}: no accepted profile ({fails.get(key)})", flush=True)
            continue
        arrays = {f: np.concatenate([p[f] for p in payloads])
                  for f in ("time_ns", "cbh_m", "c_oconnor", "s_apparent", "date")}
        order = np.argsort(arrays["time_ns"])
        arrays = {f: v[order] for f, v in arrays.items()}
        c = census[key]
        meta = dict(key=key, type=c["type"], site=c["site"], lat=c["lat"], lon=c["lon"],
                    cloudnet_site=UNITS[key][0], distance_km=UNITS[key][1], n_days_ok=len(payloads),
                    fails=fails.get(key, {}), eta_tables=os.environ["ALC_ETA_TABLES"],
                    cams_folders=CPD.CAMS_FOLDERS, code=state)
        np.savez_compressed(OUT / f"{key}_profiles.npz", **arrays, meta=json.dumps(meta))
        print(f"{key}: {arrays['cbh_m'].size} profiles on {len(payloads)} days; failures "
              f"{fails.get(key, {})}", flush=True)


if __name__ == "__main__":
    main()
