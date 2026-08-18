"""
download_lindenberg_missing.py
==============================

Download missing daily CL61 raw files for Lindenberg from the Cloudnet open API
into the layout expected by ``run_lindenberg_cl61_cal.py``:

    A:\\CL61_Cloudnet\\Lindenberg\\YYYYMMDD.nc

The Cloudnet per-file download URL embeds a per-file UUID (e.g.
``/api/download/raw/<uuid>/20241002_lindenberg_cl61_t1034980.nc``), so we cannot
just edit the date in a single URL. Instead we query the raw-files API per day to
get each file's real ``downloadUrl``, pick the **daily concatenated** file (the big
``YYYYMMDD_lindenberg_cl61_t*.nc``, not the 288 small ``live_*`` files), download it
and save it renamed to ``YYYYMMDD.nc``.

By default it scans the calibration date range (2024-05-01 .. 2026-05-31), finds the
days with no ``YYYYMMDD.nc`` present, and downloads only those. You can also pass
explicit dates.

Usage
-----
    # download every missing day in the default range
    .venv\\Scripts\\python.exe download_lindenberg_missing.py

    # download specific days
    .venv\\Scripts\\python.exe download_lindenberg_missing.py 20241002 20241008

    # custom inclusive range
    .venv\\Scripts\\python.exe download_lindenberg_missing.py --range 2024-10-01 2024-12-31

Downloads run concurrently (``--workers``, default 6 — network/disk bound, not CPU).
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta
from pathlib import Path

# --- Cloudnet API ----------------------------------------------------------
API = "https://cloudnet.fmi.fi/api/raw-files"
SITE = "lindenberg"
# CL61 instrument PID at Lindenberg (same as download_cloudnet_cl61.py).
INSTRUMENT_PID = "https://hdl.handle.net/21.12132/3.695573e5981845d9"

# --- Local layout (must match run_lindenberg_cl61_cal.py) ------------------
DEST_DIR = Path("A:/CL61_Cloudnet/Lindenberg")

DATE_START = date(2024, 5, 1)
DATE_END = date(2026, 5, 31)   # inclusive

# A real daily file is ~454 MB; treat anything well below that as incomplete.
MIN_VALID_BYTES = 100_000_000


def daterange(d0: date, d1: date):
    d = d0
    while d <= d1:
        yield d
        d += timedelta(days=1)


def present_days() -> set[str]:
    """YYYYMMDD stems already downloaded (non-empty, plausibly complete)."""
    days = set()
    if DEST_DIR.is_dir():
        for p in DEST_DIR.glob("????????.nc"):
            try:
                if p.stat().st_size >= MIN_VALID_BYTES:
                    days.add(p.stem)
            except OSError:
                pass
    return days


def query_day(ymd: str) -> list[dict]:
    """Query the raw-files API for one YYYYMMDD day."""
    iso = f"{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}"
    q = {"site": SITE, "instrumentPid": INSTRUMENT_PID, "dateFrom": iso, "dateTo": iso}
    url = API + "?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def pick_daily_file(files: list[dict], ymd: str) -> dict | None:
    """Choose the big daily concatenated file for ``ymd`` (not a ``live_*`` file).

    The daily file is named ``YYYYMMDD_lindenberg_cl61_t*.nc``. If several non-live
    candidates exist, take the largest (the full-day concat).
    """
    candidates = [
        f for f in files
        if not f["filename"].lower().startswith("live")
        and str(f.get("measurementDate", "")).replace("-", "")[:8] == ymd
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda f: int(f.get("size", 0)))


def download_one(ymd: str) -> tuple[str, str, int]:
    """Download the daily file for ``ymd`` -> ``DEST_DIR/YYYYMMDD.nc``.

    Returns (ymd, status, bytes) where status is 'ok' | 'skip' | 'missing' | 'error:<msg>'.
    """
    dest = DEST_DIR / f"{ymd}.nc"
    if dest.exists() and dest.stat().st_size >= MIN_VALID_BYTES:
        return ymd, "skip", dest.stat().st_size

    try:
        files = query_day(ymd)
    except Exception as exc:  # noqa: BLE001
        return ymd, f"error:query {exc}", 0

    daily = pick_daily_file(files, ymd)
    if daily is None:
        return ymd, "missing", 0

    tmp = dest.with_suffix(".nc.part")
    try:
        urllib.request.urlretrieve(daily["downloadUrl"], tmp)
        size = tmp.stat().st_size
        if size < MIN_VALID_BYTES:
            tmp.unlink(missing_ok=True)
            return ymd, f"error:too small ({size} B)", 0
        tmp.replace(dest)   # atomic rename only after a complete download
        return ymd, "ok", size
    except Exception as exc:  # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return ymd, f"error:download {exc}", 0


def parse_args(argv: list[str]):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dates", nargs="*",
                   help="Explicit YYYYMMDD days to fetch. If omitted, fetch all missing "
                        "days in the range.")
    p.add_argument("--range", nargs=2, metavar=("FROM", "TO"),
                   help="Inclusive date range (YYYY-MM-DD YYYY-MM-DD) to scan for missing days.")
    p.add_argument("--workers", type=int, default=6,
                   help="Concurrent downloads (default 6; network/disk bound).")
    p.add_argument("--force", action="store_true",
                   help="Re-download even if a valid local file already exists.")
    return p.parse_args(argv)


def main(argv: list[str]) -> None:
    args = parse_args(argv)
    DEST_DIR.mkdir(parents=True, exist_ok=True)

    if args.dates:
        targets = sorted(set(args.dates))
        for d in targets:
            if len(d) != 8 or not d.isdigit():
                print(f"Ignoring invalid date '{d}' (want YYYYMMDD).")
        targets = [d for d in targets if len(d) == 8 and d.isdigit()]
    else:
        if args.range:
            d0 = datetime.strptime(args.range[0], "%Y-%m-%d").date()
            d1 = datetime.strptime(args.range[1], "%Y-%m-%d").date()
        else:
            d0, d1 = DATE_START, DATE_END
        all_days = [d.strftime("%Y%m%d") for d in daterange(d0, d1)]
        have = set() if args.force else present_days()
        targets = [d for d in all_days if d not in have]

    if not targets:
        print("Nothing to download — all days present.")
        return

    print(f"Lindenberg CL61: {len(targets)} day(s) to fetch -> {DEST_DIR}")
    print(f"  {args.workers} concurrent downloads")

    n_ok = n_skip = n_missing = n_err = 0
    nbytes = 0
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(download_one, d): d for d in targets}
        for fut in as_completed(futures):
            ymd, status, size = fut.result()
            done += 1
            if status == "ok":
                n_ok += 1
                nbytes += size
                print(f"  [{done}/{len(targets)}] {ymd}: downloaded {size/1e6:.0f} MB")
            elif status == "skip":
                n_skip += 1
            elif status == "missing":
                n_missing += 1
                print(f"  [{done}/{len(targets)}] {ymd}: NOT on Cloudnet (no daily file)")
            else:
                n_err += 1
                print(f"  [{done}/{len(targets)}] {ymd}: {status}")

    print(f"\nDone. downloaded={n_ok} ({nbytes/1e9:.1f} GB)  skipped={n_skip}  "
          f"not-on-server={n_missing}  errors={n_err}")


if __name__ == "__main__":
    main(sys.argv[1:])
