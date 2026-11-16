"""Rebuild the Appendix 2 read cache after the identifier fix without re-reading the whole supplement archive.

`fig_appendix2_availability.py --rescan` re-reads every supplement daily file (about 13,000 Level-1 files of
~25 MB on F:, several hours). Only the NUMERIC identifiers (DWD CHM15k "0") were missing from the previous
scan, so this helper:
  * rescans the main hourly archive completely (35,800 small files, minutes);
  * keeps the supplement records of the letter identifiers from the previous cache (same files, same
    reader), minus the station-months that the main archive now covers;
  * reads the supplement days of the numeric identifiers that the main archive does not cover;
  * recomputes the supplement report from a fresh listing of the supplement archive.
The previous cache is kept as monthly_records_<date>.json.gz. Afterwards run
`python -m validation.paper.fig_appendix2_availability` (without --rescan) to build figure, stats and text.
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

from validation.paper import fig_appendix2_availability as A

WORKERS = 12


def supplement_from_days(by_day: dict) -> list:
    """Monthly supplement records from {(wigos, ident, day): path} (same logic as A.supplement_records)."""
    keys = list(by_day)
    if not keys:
        return []
    with ProcessPoolExecutor(WORKERS) as pool:
        results = list(pool.map(A.read_daily_hours, [by_day[k] for k in keys], chunksize=16))
    monthly = defaultdict(lambda: {"hours": set(), "meta": None})
    for (wigos, ident, day), res in zip(keys, results):
        rec = monthly[(wigos, ident, day[:6])]
        rec["hours"].update(res["hours"])
        if rec["meta"] is None or rec["meta"]["type"] == "Other":
            rec["meta"] = res
    out = []
    for (wigos, ident, ym), rec in monthly.items():
        if not rec["hours"]:
            continue
        hours = sorted(rec["hours"])
        out.append({"wigos": wigos, "ident": ident, "ym": ym, "source": "supplement", "hours": len(hours),
                    "tmin": hours[0] / 24, "tmax": (hours[-1] + 1) / 24, "type": rec["meta"]["type"],
                    "site": rec["meta"]["site"], "lat": rec["meta"]["lat"]})
    return out


def main() -> None:
    cache = A.DEFAULT_OUT / "cache" / "monthly_records.json.gz"
    with gzip.open(cache, "rt") as fh:
        old = json.load(fh)
    backup = cache.with_name("monthly_records_20261009.json.gz")
    if not backup.exists():
        shutil.copyfile(cache, backup)

    main_files = A.list_main_archive()
    print(f"main archive: {len(main_files)} hourly-monthly files", flush=True)
    with ProcessPoolExecutor(WORKERS) as pool:
        results = list(pool.map(A.read_hourly_file, [f[3] for f in main_files], chunksize=32))
    records = [{"wigos": w, "ident": i, "ym": ym, "source": "main", **res}
               for (w, i, ym, _p), res in zip(main_files, results) if res["hours"] > 0]
    main_sm = {(r["wigos"], r["ym"]) for r in records}
    n_num = sum(1 for r in records if r["ident"].isdigit())
    print(f"main records {len(records)} (numeric identifiers {n_num})", flush=True)

    daily = A.list_supplement(WORKERS)
    print(f"supplement listing: {len(daily)} daily files", flush=True)
    old_sup = [r for r in old["records"] if r["source"] == "supplement"
               and not r["ident"].isdigit() and (r["wigos"], r["ym"]) not in main_sm]
    extra_num = [d for d in daily if d[2].isdigit() and (d[1], d[3][:6]) not in main_sm
                 and not d[4].startswith("zip:")]
    by_day = {}
    for level, wigos, ident, day, path in sorted(extra_num):
        by_day.setdefault((wigos, ident, day), path)
    print(f"numeric-identifier supplement days to read: {len(by_day)}", flush=True)
    new_sup = supplement_from_days(by_day)
    sup = old_sup + new_sup

    added = sorted({r["ym"] for r in sup})
    report = {"root": A.SUPPLEMENT_ROOT, "n_daily_files_listed": len(daily),
              "n_daily_files_in_zip": sum(1 for d in daily if d[4].startswith("zip:")),
              "n_station_months_listed": len({(d[1], d[3][:6]) for d in daily}),
              "n_station_months_added": len({(r["wigos"], r["ym"]) for r in sup}),
              "n_stations_with_added_months": len({r["wigos"] for r in sup}),
              "added_month_range": [added[0], added[-1]] if added else None,
              "added_station_months_by_year": dict(sorted(Counter(
                  r["ym"][:4] for r in {(r["wigos"], r["ym"]): r for r in sup}.values()).items())),
              "note": "letter-identifier supplement records reused from the 2026-10-09 scan; numeric identifiers read on "
                      "2026-10-10 (appendix2_incremental_cache.py)"}
    with gzip.open(cache, "wt") as fh:
        json.dump({"records": records + sup, "supplement_report": report}, fh)
    print(f"cache written: {len(records)} main + {len(sup)} supplement records "
          f"({len(old_sup)} reused, {len(new_sup)} new)", flush=True)
    print("INCREMENTAL_DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
