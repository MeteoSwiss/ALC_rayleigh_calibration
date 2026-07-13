#!/usr/bin/env bash
# Daily ALC calibration pipeline -- cron entry point. Thin glue around ops/ops_daily.py:
# sources the path config, activates the venv, holds a single-instance lock, logs, and prints a
# concise run summary to stdout so cron mails it (the full ~30 MB run stays in the log file).
#
# Install (run `crontab -e`) -- e.g. 15:00 local, once yesterday's L1 files and CAMS have landed:
#   0 15 * * *  /PATH/TO/ALC_rayleigh_calibration/ops/run_daily.sh
# Cron mails a job's stdout to MAILTO (or the crond user's local mail if MAILTO is unset). Set
#   MAILTO=you@example.org
# at the top of the crontab to choose the recipient of the summary below.
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
  echo "[$(date -u +%FT%TZ)] another run holds the lock -> exit"
  exit 0
fi

# python environment (skip if running the system python)
if [ -f "${ALC_VENV:-}/bin/activate" ]; then
  # shellcheck source=/dev/null
  source "$ALC_VENV/bin/activate"
fi

# The full pipeline output (~30 MB: 433 per-stream logs + publish) goes to the LOG FILE ONLY.
# It is deliberately NOT sent to stdout: cron would try to mail all of it and the mailer drops a
# message that large -- which is exactly why the daily run produced no email. stdout carries only
# the short summary printed at the end.
echo "[$(date -u +%FT%TZ)] === ALC daily pipeline ===" >> "$LOG"
python "$OPS_DIR/ops_daily.py" "$@" >> "$LOG" 2>&1
rc=$?
echo "[$(date -u +%FT%TZ)] === exit $rc ===" >> "$LOG"

# keep only the last 30 daily logs
ls -1t "$LOG_DIR"/daily_*.log 2>/dev/null | tail -n +31 | xargs -r rm -f

# --- concise run summary -> stdout (cron mails it) + appended to the log --------------------------
{
  echo "ALC daily calibration — exit=$rc — $(date -u +%FT%TZ)"
  [ -f "${ALC_DASHBOARD_DIR:-}/.last_success" ] && echo "heartbeat: $(cat "$ALC_DASHBOARD_DIR/.last_success")"
  echo "log: $LOG"
  echo
  echo "pipeline:"
  grep -aE '^\[20[0-9-]+[ T]' "$LOG" \
    | grep -aiE 'starting|target day|fetching CAMS|CAMS.*(missing|FAILED)|calibration (ok|FAILED)|dashboard update|publishing|pipeline done' \
    | sed -E 's/ \| CAMS=.*//; s/^/  /'
  proc=$(grep -acE '^\[[0-9]+/[0-9]+\]' "$LOG")
  ok=$(grep -acE '^\[[0-9]+/[0-9]+\].*: ok' "$LOG")
  echo "  calibration streams: $ok/$proc ok"
  echo
  echo "publish:"
  grep -aE '^publish:' "$LOG" | tail -8 | sed 's/^/  /'
  # real errors (capped); the benign per-stream 'method err' and the publish summary lines are excluded
  errs=$(grep -aiE 'Traceback|CRITICAL|: ERROR|ERROR:' "$LOG" \
         | grep -avE 'method err|publish: completed WITH ERRORS|publish: ERROR:' | head -10)
  [ -n "$errs" ] && { echo; echo "errors (first 10 — see log for detail):"; printf '%s\n' "$errs" | sed 's/^/  /'; }
} | tee -a "$LOG"

exit "$rc"
