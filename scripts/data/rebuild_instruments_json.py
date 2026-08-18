"""
rebuild_instruments_json.py — rebuild instruments.json from the ACTUAL fleet present in the
E-PROFILE archive for a given month (default: May 2025), reading the instrument metadata
straight from the NetCDF files. The shipped instruments.json predates the CL61 fleet and
carries several fields the calibration code never uses; this regenerates a clean, current list.

Per (WMO, identifier) stream — a single station can host several instruments with several
identifiers (e.g. Payerne A=CHM15k, B=CL31, C=CL61), each handled as a separate entry — we read:
    Type      <- instrument_type           (global attribute)
    Latitude  <- station_latitude          (variable)
    Longitude <- station_longitude         (variable)
    Altitude  <- station_altitude          (variable)
    Serial    <- instrument_serial_number  (global attribute)
    SiteName  <- site_location / title      (global attribute), preferring the curated old value

Fields KEPT (in the file and/or used by the code: see calibration/config.py, main.py):
    WMO, Identifier, Type, SiteName, Latitude, Longitude, Altitude, Serial, Calibrated
Fields DROPPED (not in the files AND never used by the code):
    Reference, FLength, NWS, Status
`Calibrated` is not in the files but IS used (main.py: skip uncalibrated instruments), so it is
preserved from the old instruments.json where available, else defaulted to "1".

Usage:  python scripts/data/rebuild_instruments_json.py [YYYYMM] [--archive D:/E-PROFILE_L2_2025] [--write]
Without --write it only previews (writes instruments.rebuilt.json for inspection).
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import shutil
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

from netCDF4 import Dataset

REPO = Path(__file__).resolve().parents[2]
OLD_JSON = REPO / "instruments.json"

KEEP_TYPES = {"CHM15k", "CHM8k", "CL31", "CL51", "CL61", "Mini-MPL", "MPL"}


def parse_name(fname: str):
    """Handle monthly 'L2_<WMO>_<id>202505.nc' and daily 'L2_<WMO>_<id>20250501.nc'."""
    core = fname[3:-3]                      # strip 'L2_' / 'L1_' and '.nc'
    if core[-8:].isdigit():
        date, wmo_ident = core[-8:], core[:-8]      # daily YYYYMMDD
    elif core[-6:].isdigit():
        date, wmo_ident = core[-6:], core[:-6]      # monthly YYYYMM
    else:
        return None
    if "_" not in wmo_ident:
        return None
    wmo, ident = wmo_ident.rsplit("_", 1)
    return wmo, ident, date


def scan_month(archive: Path, ym: str):
    """{(wmo, ident): one representative file path} for month ym (YYYYMM). Supports both the
    monthly layout (<WMO>/<YYYY>/L2_..<YYYYMM>.nc) and the daily layout
    (<WMO>/<YYYY>/<MM>/L2_..<YYYYMMDD>.nc)."""
    one = {}
    yyyy, mm = ym[:4], ym[4:6]
    for wmo_dir in sorted(archive.iterdir()):
        if not wmo_dir.is_dir():
            continue
        for pat in (str(wmo_dir / yyyy / f"L?_*{ym}.nc"),
                    str(wmo_dir / yyyy / mm / f"L?_*{ym}*.nc")):
            for f in glob.glob(pat):
                p = parse_name(os.path.basename(f))
                if p is None:
                    continue
                wmo, ident, _ = p
                one.setdefault((wmo, ident), f)
    return one


def read_meta(path: str):
    with Dataset(path, "r") as d:
        def attr(name, default=""):
            return str(getattr(d, name)) if name in d.ncattrs() else default

        def var0(name, default=None):
            if name in d.variables:
                try:
                    import numpy as np
                    return float(np.asarray(d.variables[name][:]).ravel()[0])
                except Exception:
                    return default
            return default

        itype = attr("instrument_type", "?")
        serial = attr("instrument_serial_number", "")
        site = attr("site_location", "") or attr("title", "").split(" ")[0]
        lat = var0("station_latitude")
        lon = var0("station_longitude")
        alt = var0("station_altitude", 0.0)
    return itype, serial, site, lat, lon, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("month", nargs="?", default="202505", help="YYYYMM (default 202505)")
    ap.add_argument("--archive", default="A:/E-PROFILE_L2_monthly",
                    help="archive root (default A:/E-PROFILE_L2_monthly — full fleet, monthly files)")
    ap.add_argument("--write", action="store_true", help="overwrite instruments.json (backs up first)")
    args = ap.parse_args()

    archive = Path(args.archive)
    old = {(e["WMO"], e.get("Identifier", "A")): e for e in json.load(open(OLD_JSON))}
    streams = scan_month(archive, args.month)
    print(f"{args.month}: {len(streams)} (WMO, identifier) streams in {archive}")

    entries = []
    census = defaultdict(int)
    skipped = 0
    multi = defaultdict(list)
    for (wmo, ident), f in sorted(streams.items()):
        try:
            itype, serial, site, lat, lon, alt = read_meta(f)
        except Exception as e:
            print(f"  WARN {wmo}_{ident}: {e}")
            skipped += 1
            continue
        if itype not in KEEP_TYPES or lat is None or lon is None:
            skipped += 1
            continue
        prev = old.get((wmo, ident), {})
        # Prefer the curated SiteName / Calibrated from the old list where present.
        site_name = prev.get("SiteName") or site or wmo
        calibrated = prev.get("Calibrated", "1")
        entries.append({
            "WMO": wmo,
            "Identifier": ident,
            "Type": itype,
            "SiteName": site_name,
            "Latitude": round(lat, 5),
            "Longitude": round(lon, 5),
            "Altitude": round(alt or 0.0, 1),
            "Serial": serial,
            "Calibrated": calibrated,
        })
        census[itype] += 1
        multi[wmo].append(ident)

    entries.sort(key=lambda e: (e["WMO"], e["Identifier"]))
    multi_stations = {w: ids for w, ids in multi.items() if len(ids) > 1}

    print(f"\nBuilt {len(entries)} instrument entries ({skipped} skipped: unknown type / no coords)")
    print("Type census:", dict(sorted(census.items(), key=lambda kv: -kv[1])))
    print(f"Multi-instrument stations: {len(multi_stations)} "
          f"(e.g. {dict(list(multi_stations.items())[:4])})")
    n_new = sum(1 for e in entries if (e["WMO"], e["Identifier"]) not in old)
    print(f"New since old list: {n_new};  old had {len(old)} entries")

    preview = REPO / "instruments.rebuilt.json"
    preview.write_text(json.dumps(entries, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote preview -> {preview.name}")

    if args.write:
        backup = REPO / f"instruments.backup_{date.today().isoformat()}.json"
        shutil.copy2(OLD_JSON, backup)
        OLD_JSON.write_text(json.dumps(entries, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"Backed up old -> {backup.name};  overwrote instruments.json")
    else:
        print("(dry run — pass --write to overwrite instruments.json)")


if __name__ == "__main__":
    main()
