#!/usr/bin/env python3
"""Targeted corrected-Rayleigh calibration for specific instruments.

Reuses run_all_l2monthly's per-station machinery but only for the keys passed on
the command line, so we can fill in the few instruments needed for the paper
validation without re-running the whole network.

Usage:
    python run_targeted_rayleigh.py [--workers N] KEY [KEY ...]
        KEY = "<WMO>_<identifier>"   e.g.  0-20000-0-06610_A

Writes C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all/<KEY>/<KEY>_cl.csv in the
same format as run_all_l2monthly.py. Resumable: a KEY whose _cl.csv already
exists is skipped.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import run_all_l2monthly as R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("keys", nargs="+", help="instrument keys '<WMO>_<identifier>'")
    args = ap.parse_args()

    manifest = json.load(open("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json"))
    by_key = {f"{s['wmo']}_{s['identifier']}": s for s in manifest}

    R.OUT.mkdir(parents=True, exist_ok=True)
    for key in args.keys:
        if key not in by_key:
            print(f"[skip] {key} not in manifest", flush=True)
            continue
        station = by_key[key]
        sdir = R.OUT / key
        if (sdir / f"{key}_cl.csv").exists():
            print(f"[done] {key} already calibrated", flush=True)
            continue

        print(f"[run ] {key} ({station['itype']}) ...", flush=True)
        rows = R._calibrate_station(station, args.workers)

        cl = [float(r["lidar_constant"]) for r in rows
              if str(r["flag"]) in ("1", "1.0", "0.5")]
        sdir.mkdir(parents=True, exist_ok=True)
        with open(sdir / f"{key}_cl.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=R._FIELDS)
            w.writeheader()
            w.writerows(rows)
        # Drop per-year checkpoints now the instrument is fully banked.
        for pp in (sdir / "_partial").glob("*.csv"):
            pp.unlink()
        med = float(np.median(cl)) if cl else float("nan")
        print(f"[ok  ] {key}: {len(cl)}/{len(rows)} nights ok, med={med:.4e}", flush=True)

    print("TARGETED DONE", flush=True)


if __name__ == "__main__":
    main()
