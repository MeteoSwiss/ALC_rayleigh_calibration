#!/usr/bin/env python3
"""Fill the empty ``version`` cell of the cloud rows in a v2.3 tree's <key>_cal.csv with 1100.

Up to the v2.3 release run the runner left the CSV version column empty for cloud rows (the NetCDF
always carried the code); the runner writes it since commit 2e42786. Only EMPTY cells of cloud rows are
filled, and only with the code the stream's NetCDF already carries for its cloud rows (1100, eta tables
v2), so the CSV and the NetCDF say the same thing. Rows that already carry a version are left alone.

    python scripts/fill_cloud_version.py <tree> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path

CLOUD_CODE = "1100"


def fill(cal: Path, dry_run: bool) -> int:
    with cal.open(newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        fields = rd.fieldnames or []
        rows = list(rd)
    if "version" not in fields:
        return 0
    n = 0
    for r in rows:
        if r.get("method") == "cloud" and not str(r.get("version", "")).strip():
            r["version"] = CLOUD_CODE
            n += 1
    if n and not dry_run:
        fd, tmp = tempfile.mkstemp(dir=cal.parent, prefix=".fill_", suffix=".csv")
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        os.chmod(tmp, cal.stat().st_mode & 0o777)
        os.replace(tmp, cal)
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    total = streams = 0
    for cal in sorted(a.tree.glob("*/*_cal.csv")):
        n = fill(cal, a.dry_run)
        if n:
            streams += 1
            total += n
    print(f"{'would fill' if a.dry_run else 'filled'} {total} cloud rows with version {CLOUD_CODE} "
          f"in {streams} streams under {a.tree}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
