#!/usr/bin/env python3
"""Daily operational driver for the ALC calibration pipeline (run from ops/run_daily.sh via cron).

For each target day -- the day D-LAG (default yesterday, since CAMS lands the next day) plus a small
backfill of the last few days that have not been processed yet (self-healing after an outage) -- it:

  1. ensures the CAMS file for that day (downloads it from the ADS if missing; retried -- this is the
     slow, network-dependent step, deliberately isolated so a failure is one clear log line);
  2. runs the calibration for that day across the whole network: Rayleigh + liquid-cloud + Kalman
     (scripts/run_all_l1_2026.py, which now MERGES into the per-stream CSVs instead of overwriting);

then, once, if anything changed:

  3. updates the dashboard incrementally (build_dashboard.py --changed-only -> only the changed station
     pages re-render, the summary always does);
  4. writes a heartbeat (ALC_DASHBOARD_DIR/.last_success) so a monitor can tell the pipeline is alive.

All paths come from the ALC_* env vars (ops/config.sh). Dates are UTC. A day is recorded as processed
once its CAMS was present and the calibration ran, so it is not redone; CAMS-missing days are retried
on the next run. Exit code is non-zero only on a hard failure (dashboard build failed, or there was
work to do but no day could be processed).

Manual use:
    python ops/ops_daily.py                 # normal daily run (D-LAG + backfill)
    python ops/ops_daily.py --day 20260401  # (re)process one specific day
    python ops/ops_daily.py --no-dashboard  # calibration only
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from calibration.io.cams import ensure_cams_file   # noqa: E402


def _envp(name: str, default: str) -> Path:
    return Path(os.environ.get(name) or default)


CAMS_DIR = _envp("ALC_CAMS_DIR", "D:/CAMS")
FULLCAL_DIR = _envp("ALC_FULLCAL_DIR", str(REPO / "_fullcal"))
DASHBOARD_DIR = _envp("ALC_DASHBOARD_DIR", str(REPO / "_dashboard"))
OPCOEFF_CSV = os.environ.get("ALC_OPCOEFF_CSV") or ""
L2_DIR = os.environ.get("ALC_L2_DIR") or ""
DAY_LAG = int(os.environ.get("ALC_DAY_LAG") or "1")
BACKFILL_DAYS = int(os.environ.get("ALC_BACKFILL_DAYS") or "5")
WORKERS = str(os.environ.get("ALC_WORKERS") or "6")
PY = sys.executable
PUBLISH = (os.environ.get("ALC_PUBLISH") or "0") == "1"
PUBLISH_SH = REPO / "ops" / "publish.sh"

PROCESSED = DASHBOARD_DIR / ".processed_days"     # one YYYYMMDD per line: days CAMS was present + run
HEARTBEAT = DASHBOARD_DIR / ".last_success"


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z] {msg}", flush=True)


def load_processed() -> set:
    try:
        return set(PROCESSED.read_text(encoding="utf-8").split())
    except OSError:
        return set()


def mark_processed(ds: str) -> None:
    done = load_processed()
    done.add(ds)
    PROCESSED.parent.mkdir(parents=True, exist_ok=True)
    PROCESSED.write_text("\n".join(sorted(done)) + "\n", encoding="utf-8")


def target_days(today_utc: datetime) -> list:
    """The day D-LAG plus the BACKFILL_DAYS days before it (oldest first)."""
    base = today_utc - timedelta(days=DAY_LAG)
    days = [(base - timedelta(days=k)).strftime("%Y%m%d") for k in range(BACKFILL_DAYS, -1, -1)]
    return sorted(set(days))


def fetch_cams(ds: str, retries: int = 3) -> bool:
    for attempt in range(1, retries + 1):
        try:
            path = ensure_cams_file(CAMS_DIR, ds, auto_download=True, scope="day")
            if path is not None:
                return True
            log(f"  CAMS {ds}: not found and download returned nothing (attempt {attempt}/{retries})")
        except Exception as exc:  # noqa: BLE001 - network/credential/grib failures must not crash the run
            log(f"  CAMS {ds}: attempt {attempt}/{retries} failed: {type(exc).__name__}: {exc}")
        if attempt < retries:
            time.sleep(min(60 * attempt, 300))
    return False


def calibrate(ds: str) -> bool:
    cmd = [PY, str(REPO / "scripts" / "run_all_l1_2026.py"),
           "--start", ds, "--end", ds, "--per-type", "0", "--ignore-coverage",
           "--workers", WORKERS, "--methods", "rayleigh,cloud", "--force",
           "--sens", "--omb"]
    log(f"  calibrate {ds}: {' '.join(cmd[1:])}")
    return subprocess.run(cmd, cwd=str(REPO)).returncode == 0


def update_opcoeff(ds: str) -> None:
    """Append the day's operational L2 calibration_constant_0 to the opcoeff CSV (the
    dashboard's %-of-operational ratio). Resumable + best-effort -- never fails the run."""
    if not (L2_DIR and OPCOEFF_CSV):
        return
    cmd = [PY, str(REPO / "scripts" / "extract_l2_opcoeff.py"), L2_DIR, OPCOEFF_CSV,
           "--start", ds, "--end", ds]
    log(f"  opcoeff {ds}: extract_l2_opcoeff --start {ds} --end {ds}")
    try:
        subprocess.run(cmd, cwd=str(REPO), timeout=1800)
    except Exception as exc:  # noqa: BLE001
        log(f"  opcoeff {ds}: failed: {type(exc).__name__}: {exc}")


def update_dashboard() -> bool:
    cmd = [PY, str(REPO / "scripts" / "build_dashboard.py"),
           "--fullcal", str(FULLCAL_DIR), "--out", str(DASHBOARD_DIR), "--changed-only"]
    if L2_DIR:
        cmd += ["--l2dir", L2_DIR]
    if OPCOEFF_CSV:
        cmd += ["--opcoeff", OPCOEFF_CSV]
    log(f"  dashboard: {' '.join(cmd[1:])}")
    return subprocess.run(cmd, cwd=str(REPO)).returncode == 0


def publish() -> bool:
    """Publish the freshly-built site to the EWC (images -> S3 bucket, HTML -> web VM) via
    ops/publish.sh. Best-effort: a publish failure is logged and flagged in the heartbeat but does NOT
    fail the run -- the calibration and the local dashboard already succeeded."""
    if not PUBLISH_SH.exists():
        log("  publish: ops/publish.sh missing -> skip")
        return False
    log(f"  publish: bash {PUBLISH_SH}")
    return subprocess.run(["bash", str(PUBLISH_SH)], cwd=str(REPO)).returncode == 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--day", default=None, help="process this single day YYYYMMDD (overrides D-LAG + backfill)")
    ap.add_argument("--no-dashboard", action="store_true", help="run the calibration only, skip the dashboard")
    ap.add_argument("--no-publish", action="store_true", help="skip the EWC publish step even if ALC_PUBLISH=1")
    ap.add_argument("--force-all", action="store_true", help="ignore the processed-days record (reprocess the window)")
    args = ap.parse_args()

    log(f"ALC daily pipeline starting | CAMS={CAMS_DIR} | fullcal={FULLCAL_DIR} | dashboard={DASHBOARD_DIR}")

    if args.day:
        days = [args.day]
    else:
        processed = set() if args.force_all else load_processed()
        days = [d for d in target_days(datetime.now(timezone.utc)) if d not in processed]

    if not days:
        log("nothing to do (all target days already processed)")
        HEARTBEAT.parent.mkdir(parents=True, exist_ok=True)
        HEARTBEAT.write_text(f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z idle\n", encoding="utf-8")
        return 0

    log(f"target days: {', '.join(days)}")
    ran, cams_failures = [], []
    for ds in days:
        log(f"day {ds}: fetching CAMS ...")
        if not fetch_cams(ds):
            log(f"day {ds}: CAMS unavailable -> skip (will retry next run)")
            cams_failures.append(ds)
            continue
        ok = calibrate(ds)
        update_opcoeff(ds)
        mark_processed(ds)                       # CAMS was present + we ran; don't redo even if some streams failed
        ran.append(ds)
        log(f"day {ds}: calibration {'ok' if ok else 'completed WITH ERRORS (see per-stream output)'}")

    dash_ok = True
    published = None
    if ran and not args.no_dashboard:
        log("updating dashboard (incremental) ...")
        dash_ok = update_dashboard()
        log(f"dashboard update {'ok' if dash_ok else 'FAILED'}")
        if dash_ok and PUBLISH and not args.no_publish:
            log("publishing to the EWC ...")
            published = publish()
            log(f"publish {'ok' if published else 'FAILED (non-fatal)'}")
    elif not ran:
        log("no day processed -> dashboard not updated")

    HEARTBEAT.parent.mkdir(parents=True, exist_ok=True)
    HEARTBEAT.write_text(
        f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z ran={','.join(ran) or '-'} "
        f"cams_missing={','.join(cams_failures) or '-'} dashboard={'ok' if dash_ok else 'fail'}"
        f"{'' if published is None else (' publish=' + ('ok' if published else 'fail'))}\n",
        encoding="utf-8")

    # hard failure: dashboard build failed, or there was work but nothing could be processed
    if not dash_ok:
        log("FAILURE: dashboard build failed")
        return 1
    if days and not ran:
        log("FAILURE: had target days but none could be processed (CAMS unavailable for all)")
        return 1
    log("ALC daily pipeline done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
