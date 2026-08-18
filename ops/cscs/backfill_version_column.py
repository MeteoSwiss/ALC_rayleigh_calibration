"""Stamp the `version` column into a calibration archive written before that column existed.

Why this is a backfill and not a re-run
--------------------------------------
The v2.2 release tree was produced by `ALC_MOLECULAR_METHOD=eprof_v2.2`, but by code that predates
the `version` field in CSV_FIELDS. The numbers are v2.2; only the provenance stamp is missing, and
without it the dashboard cannot name the vintage and silently falls back to an unlabelled series --
the exact failure that once presented a v2.2 archive as "v2.0".

Re-running the network would produce byte-identical constants and cost hours, so this writes the
column the current runner would have written:

  * rayleigh rows -> VERSION_CODES["eprof_v2.2"] = 220
  * cloud rows    -> VERSION_CODES["cloud_oconnor"] = 1000

The rayleigh value is exact rather than assumed: `version_code(o.molecular_method)` records the
method that was REQUESTED, and v2.2 falling back to the v2 gates for a noise-starved night does not
change that name (molecular_methods.py) -- so 220 is what the runner writes for every row of a v2.2
run, fallback or not.

Safety
------
* refuses to touch a file whose header already has `version` (idempotent);
* refuses a method it does not recognise, rather than guessing;
* writes to a temp file in the same directory and renames -> a crash cannot truncate an archive;
* --dry-run reports what it would do and changes nothing.

    python backfill_version_column.py <tree> [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

RAYLEIGH_CODE = "220"          # VERSION_CODES["eprof_v2.2"]
CLOUD_CODE = "1000"            # VERSION_CODES["cloud_oconnor"]
FIELDS = ["date", "method", "flag", "cal_value", "uncertainty",
          "n_profiles", "bottom_height", "top_height", "version", "message"]


def code_for(method: str) -> str:
    m = (method or "").strip().lower()
    if m == "rayleigh":
        return RAYLEIGH_CODE
    if m in ("cloud", "cloud_oconnor", "liquid-cloud", "liquid_cloud"):
        return CLOUD_CODE
    raise ValueError(f"unrecognised method {method!r} — refusing to guess its version")


def backfill(path: Path, dry: bool) -> tuple[str, int]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
        f.seek(0)
        header = f.readline().strip().split(",")
    if "version" in header:
        return "already stamped", len(rows)
    if not rows:
        return "empty", 0
    out = []
    for r in rows:
        r = {k: (r.get(k) or "") for k in FIELDS if k != "version"} | {
            "version": code_for(r.get("method", ""))}
        out.append(r)
    if dry:
        return "would stamp", len(out)
    tmp = path.with_suffix(".csv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    os.replace(tmp, path)
    return "stamped", len(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tree", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if not a.tree.is_dir():
        print(f"no such tree: {a.tree}", file=sys.stderr)
        return 2
    tally: dict[str, int] = {}
    rows = 0
    bad = []
    for d in sorted(p for p in a.tree.iterdir() if p.is_dir()):
        f = d / f"{d.name}_cal.csv"
        if not f.exists():
            continue
        try:
            what, n = backfill(f, a.dry_run)
        except (ValueError, OSError) as exc:
            bad.append(f"{d.name}: {exc}")
            continue
        tally[what] = tally.get(what, 0) + 1
        rows += n
    for k, v in sorted(tally.items()):
        print(f"  {k:16s} {v:4d} files")
    print(f"  rows seen        {rows}")
    for b in bad:
        print(f"  FAILED {b}", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
