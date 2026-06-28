#!/usr/bin/env bash
# Pre-download the day's CAMS boxes (Europe + the small regional boxes) ahead of the 15:00
# calibration, so no ADS round-trip blocks the daily run. Thin glue around
# scripts/prefetch_cams.py: source config, activate venv, single-instance lock, log.
#
# Install (crontab -e) -- a few minutes past 06:00, once the day's CAMS forecast has landed:
#   2 6 * * *  /data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code/ops/prefetch_cams.sh
set -uo pipefail

OPS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$OPS_DIR/config.sh"

LOG_DIR="${ALC_LOG_DIR:-$OPS_DIR/logs}"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/prefetch_cams_$(date -u +%Y%m%d_%H%M%S).log"

# single-instance lock (don't overlap with a previous slow download)
exec 9>"${ALC_PREFETCH_LOCK:-$LOG_DIR/.prefetch_cams.lock}"
if ! flock -n 9; then
  echo "[$(date -u +%FT%TZ)] another prefetch holds the lock -> exit" | tee -a "$LOG"
  exit 0
fi

if [ -f "${ALC_VENV:-}/bin/activate" ]; then
  # shellcheck source=/dev/null
  source "$ALC_VENV/bin/activate"
fi

echo "[$(date -u +%FT%TZ)] === CAMS prefetch ===" | tee -a "$LOG"
python "$OPS_DIR/../scripts/prefetch_cams.py" "$@" 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
echo "[$(date -u +%FT%TZ)] === exit $rc ===" | tee -a "$LOG"

# keep only the last 30 prefetch logs
ls -1t "$LOG_DIR"/prefetch_cams_*.log 2>/dev/null | tail -n +31 | xargs -r rm -f
exit "$rc"
