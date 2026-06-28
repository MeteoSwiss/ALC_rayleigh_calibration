#!/usr/bin/env python3
"""Refresh the station census (ALC_CENSUS) from the live L1 archive.

Scans the L1 archive for every ceilometer/lidar stream (a station may host several
streams under different identifiers) and MERGES what it finds into the existing census
JSON, so the daily calibration automatically picks up newly-installed stations without a
manual edit. Existing entries are updated in place (coverage ``last``/``n_days`` extended,
metadata refreshed); brand-new streams are appended; entries already in the census are
never dropped even if a stream stops reporting.

Archive layout (E-PROFILE):
    <L1_ROOT>/<wmo>/<year>/<month>/L1_<wmo>_<ident><YYYYMMDD>.nc

Paths come from the ALC_* env vars (ops/config.sh): ALC_L1_ROOT, ALC_CENSUS. Run stand-alone
or let ops_daily.py call it at the start of each daily run.

    python scripts/refresh_census.py            # refresh ALC_CENSUS in place
    python scripts/refresh_census.py --year 2026
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
from netCDF4 import Dataset

REPO = Path(__file__).resolve().parents[1]
L1_ROOT = Path(os.environ.get("ALC_L1_ROOT", "D:/E-PROFILE_L1_2026"))
CENSUS = Path(os.environ.get("ALC_CENSUS", str(REPO / "validation" / "scope_l1_2026_census.json")))

# Instrument types the calibration knows how to process (matches run_all_l1_2026.ITYPE).
TARGET_TYPES = {"CL31", "CL51", "CL61", "CHM15k", "Mini-MPL"}
DATE_RE = re.compile(r"(\d{8})\.nc$")


def _scalar(d, *names, default=0.0):
    """First finite scalar among the named variables, else *default*."""
    for n in names:
        if n in d.variables:
            a = np.asarray(d.variables[n][:]).ravel()
            if a.size and np.isfinite(a[0]):
                return float(a[0])
    return default


def _key(wmo: str, ident: str) -> str:
    return f"{wmo}_{ident}"


def scan_archive(root: Path, year: str) -> dict:
    """Return {key: stream_dict} for every TARGET stream found under *root* for *year*."""
    found: dict[str, dict] = {}
    if not root.exists():
        print(f"[refresh] L1 root does not exist: {root}", flush=True)
        return found

    for wmo in sorted(os.listdir(root)):
        files = glob.glob(str(root / wmo / year / "*" / f"L1_{wmo}_*.nc"))
        if not files:
            continue
        by_id: dict[str, list] = defaultdict(list)
        for f in files:
            rest = os.path.basename(f)[len(f"L1_{wmo}_"):-3]   # <ident><YYYYMMDD>
            if len(rest) <= 8:
                continue
            by_id[rest[:-8]].append(f)

        for ident, fs in by_id.items():
            fs = sorted(fs)
            dates = sorted(m.group(1) for m in (DATE_RE.search(os.path.basename(f)) for f in fs) if m)
            if not dates:
                continue
            try:
                with Dataset(fs[-1]) as d:                       # newest file: freshest metadata
                    itype = str(getattr(d, "instrument_type", "?"))
                    if itype not in TARGET_TYPES:
                        continue
                    lat = _scalar(d, "station_latitude", "latitude")
                    lon = _scalar(d, "station_longitude", "longitude")
                    alt = _scalar(d, "station_altitude", "altitude")
                    site = str(getattr(d, "site_location",
                                       getattr(d, "wigos_station_id", wmo))).split(",")[0]
            except Exception as exc:  # noqa: BLE001 - a single unreadable file must not abort the scan
                print(f"[refresh] skip {wmo}_{ident}: {type(exc).__name__}: {exc}", flush=True)
                continue
            found[_key(wmo, ident)] = dict(
                wmo=wmo, ident=ident, type=itype, site=site,
                lat=lat, lon=lon, alt=alt,
                first=dates[0], last=dates[-1], n_days=len(dates),
            )
    return found


def merge(existing: list, found: dict) -> tuple[list, int, int]:
    """Merge scanned streams into the existing census list.

    New streams are appended; existing ones have coverage/metadata refreshed (and ``first``
    kept at the earliest ever seen). Returns (merged_list, n_new, n_updated)."""
    by_key = {_key(s["wmo"], s["ident"]): s for s in existing}
    n_new = n_updated = 0

    for key, info in found.items():
        if key in by_key:
            cur = by_key[key]
            changed = False
            # extend coverage (census first/last span the union of what we have ever seen)
            if info["first"] < cur.get("first", info["first"]):
                cur["first"] = info["first"]; changed = True
            if info["last"] > cur.get("last", ""):
                cur["last"] = info["last"]; changed = True
            for f in ("type", "site", "lat", "lon", "alt", "n_days"):
                if cur.get(f) != info[f]:
                    cur[f] = info[f]; changed = True
            n_updated += int(changed)
        else:
            by_key[key] = info
            n_new += 1

    merged = sorted(by_key.values(), key=lambda s: _key(s["wmo"], s["ident"]))
    return merged, n_new, n_updated


def _write_atomic(path: Path, data: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1)
        os.chmod(tmp, 0o644)   # mkstemp is 0600; keep the census group/other-readable (dashboard reads it)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def refresh(root: Path = L1_ROOT, census_path: Path = CENSUS, year: str = "2026") -> int:
    """Scan *root* and merge into *census_path* in place. Returns the number of new streams."""
    try:
        existing = json.loads(census_path.read_text(encoding="utf-8")) if census_path.exists() else []
    except Exception as exc:  # noqa: BLE001
        print(f"[refresh] could not read census {census_path}: {exc}; starting fresh", flush=True)
        existing = []

    found = scan_archive(root, year)
    if not found:
        print("[refresh] no streams found in archive -> census left unchanged", flush=True)
        return 0

    merged, n_new, n_updated = merge(existing, found)
    _write_atomic(census_path, merged)
    print(f"[refresh] census {census_path.name}: {len(merged)} streams "
          f"(+{n_new} new, {n_updated} updated)", flush=True)
    if n_new:
        new_keys = [k for k in (_key(s["wmo"], s["ident"]) for s in merged)
                    if k not in {_key(s["wmo"], s["ident"]) for s in existing}]
        print(f"[refresh] new streams: {', '.join(sorted(new_keys))}", flush=True)
    return n_new


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", default="2026", help="archive year subfolder to scan (default 2026)")
    ap.add_argument("--l1-root", default=str(L1_ROOT), help="override L1 archive root")
    ap.add_argument("--census", default=str(CENSUS), help="override census JSON path")
    args = ap.parse_args()
    refresh(Path(args.l1_root), Path(args.census), args.year)
    return 0


if __name__ == "__main__":
    sys.exit(main())
