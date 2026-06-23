#!/usr/bin/env bash
# Daily ALC calibration pipeline -- cron entry point. Thin glue around ops/ops_daily.py:
# sources the path config, activates the venv, holds a single-instance lock, logs, alerts on failure.
#
# Install (run `crontab -e`) -- e.g. 06:30 UTC, once yesterday's L1 files and CAMS have landed:
#   30 6 * * *  /PATH/TO/ALC_rayleigh_calibration/ops/run_daily.sh
#
# NOT `set -e`: we capture ops_daily.py's exit code and act on it ourselves.
set -uo pipefail

OPS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$OPS_DIR/config.sh"

LOG_DIR="${ALC_LOG_DIR:-$OPS_DIR/logs}"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/daily_$(date -u +%Y%m%d_%H%M%S).log"

# single-instance lock: if a previous run is still going, exit quietly (don't pile up)
exec 9>"${ALC_LOCK:-$LOG_DIR/.daily.lock}"
if ! flock -n 9; then
  echo "[$(date -u +%FT%TZ)] another run holds the lock -> exit" | tee -a "$LOG"
  exit 0
fi

# python environment (skip if running the system python)
if [ -f "${ALC_VENV:-}/bin/activate" ]; then
  # shellcheck source=/dev/null
  source "$ALC_VENV/bin/activate"
fi

echo "[$(date -u +%FT%TZ)] === ALC daily pipeline ===" | tee -a "$LOG"
python "$OPS_DIR/ops_daily.py" "$@" 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
echo "[$(date -u +%FT%TZ)] === exit $rc ===" | tee -a "$LOG"

# keep only the last 30 daily logs
ls -1t "$LOG_DIR"/daily_*.log 2>/dev/null | tail -n +31 | xargs -r rm -f

# email the log tail on failure (best effort; needs a working `mail`)
if [ "$rc" -ne 0 ] && [ -n "${ALC_ALERT_EMAIL:-}" ]; then
  tail -n 50 "$LOG" | mail -s "ALC daily pipeline FAILED (rc=$rc)" "$ALC_ALERT_EMAIL" 2>/dev/null || true
fi
exit "$rc"
