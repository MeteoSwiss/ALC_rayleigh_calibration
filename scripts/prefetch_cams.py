#!/usr/bin/env python3
"""Pre-download every CAMS box for the day's calibration (cron ~06:00, ahead of 15:00).

The daily calibration (ops_daily.py) pre-fetches the Europe box; the small regional boxes
that serve the far-flung affiliate stations would otherwise be fetched lazily during the
calibration. This script downloads them ALL up front — the default Europe box plus each
regional box in ``download_cams_beta.CAMS_REGIONS`` — so nothing waits on an ADS round-trip
when the 15:00 calibration starts.

Target day = D - ALC_DAY_LAG (default yesterday) + the ALC_BACKFILL_DAYS preceding days,
mirroring ops_daily's self-healing window. Idempotent (ensure_cams_file skips a box whose
file already exists) and best-effort (a box that is not yet on the ADS, or fails, is logged
and retried on the next run / by the 15:00 calibration).

    python scripts/prefetch_cams.py             # the daily window (D-LAG + backfill)
    python scripts/prefetch_cams.py 20260601    # one specific day
"""
from __future__ import annotations

import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from calibration.io.cams import ensure_cams_file            # noqa: E402
from calibration.io.download_cams_beta import CAMS_REGIONS   # noqa: E402

CAMS_DIR = os.environ.get("ALC_CAMS_DIR", str(REPO / "_cams"))
DAY_LAG = int(os.environ.get("ALC_DAY_LAG") or "1")
BACKFILL_DAYS = int(os.environ.get("ALC_BACKFILL_DAYS") or "5")


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z] {msg}", flush=True)


def target_days(today_utc: datetime) -> list:
    base = today_utc - timedelta(days=DAY_LAG)
    days = [(base - timedelta(days=k)).strftime("%Y%m%d") for k in range(BACKFILL_DAYS, -1, -1)]
    return sorted(set(days))


def main() -> int:
    days = [sys.argv[1]] if len(sys.argv) > 1 else target_days(datetime.now(timezone.utc))
    regions = ["europe"] + list(CAMS_REGIONS)
    log(f"prefetch CAMS | dir={CAMS_DIR} | days={','.join(days)} | regions={','.join(regions)}")
    ready = missing = 0
    for ds in days:
        for region in regions:
            try:
                path = ensure_cams_file(CAMS_DIR, ds, auto_download=True, scope="day",
                                        require_backscatter=True, region=region)
            except Exception as exc:  # noqa: BLE001 - one box failing must not abort the rest
                path = None
                log(f"  {ds} {region}: ERROR {type(exc).__name__}: {exc}")
            if path is not None:
                ready += 1
                log(f"  {ds} {region}: ready ({Path(path).name})")
            else:
                missing += 1
                log(f"  {ds} {region}: NOT available yet (will retry next run / at 15:00)")
    log(f"prefetch done: {ready} ready, {missing} missing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
