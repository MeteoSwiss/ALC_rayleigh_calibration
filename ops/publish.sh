#!/usr/bin/env bash
# Publish the built dashboard to the European Weather Cloud (or any S3 + web host):
#   * the bulky image trees (diag/ ombsens/ flagex/)  -> a public S3 bucket (rclone or aws),
#   * the static HTML + assets                          -> the web server's docroot (rsync over ssh).
# Driven by the ALC_* vars in ops/config.sh. Each leg is skipped (with a note) when its target is not
# configured. Symlinked diagnostics are dereferenced on upload so the real bytes land in the bucket.
#
# On failure each leg reports WHICH leg, its exit code (with the rsync exit-code meaning), and the
# tail of its actual error output -- so a red publish is diagnosable from the log, not a bare rc.
#
# Behind a proxy (e.g. a MeteoSwiss server): the S3 client (aws/rclone) honours http_proxy/https_proxy
# from the environment; the VM rsync uses $ALC_VM_SSH, where you put the key and a ProxyCommand. See
# ops/config.sh + ops/README.md.
#
# Manual use:  ops/publish.sh        # publish whatever is configured
# From cron:   called by ops/ops_daily.py after a successful build when ALC_PUBLISH=1.
set -uo pipefail

OPS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
# Config file is overridable via ALC_CONFIG (e.g. a home/laptop config with no proxy and local
# paths); defaults to the committed server config.sh so cron behaviour is unchanged.
source "${ALC_CONFIG:-$OPS_DIR/config.sh}"

SITE="${ALC_DASHBOARD_DIR:?ALC_DASHBOARD_DIR not set}"
if [ ! -d "$SITE" ]; then
  echo "publish: site dir '$SITE' does not exist -> nothing to publish" >&2
  exit 1
fi

S3_TOOL="${ALC_S3_TOOL:-rclone}"     # rclone | aws
VM_SSH="${ALC_VM_SSH:-ssh}"          # ssh command for the rsync (add -i KEY and a ProxyCommand for a proxy)
rc=0
FAILED=()                            # human-readable list of legs that failed (for the summary line)

# Human meaning of an rsync exit code (man rsync, EXIT VALUES) -- so the log says WHAT the code means.
rsync_meaning() {
  case "$1" in
    1)   echo "syntax or usage error" ;;
    2)   echo "protocol incompatibility" ;;
    3)   echo "errors selecting input/output files or dirs" ;;
    4)   echo "requested action not supported" ;;
    5)   echo "error starting client-server protocol (often the remote shell/rsync)" ;;
    10)  echo "error in socket I/O" ;;
    11)  echo "error in file I/O" ;;
    12)  echo "error in rsync protocol data stream" ;;
    13)  echo "errors with program diagnostics" ;;
    22)  echo "error allocating core memory buffers" ;;
    23)  echo "partial transfer due to error (some files were not transferred)" ;;
    24)  echo "partial transfer due to vanished source files" ;;
    30)  echo "timeout in data send/receive" ;;
    35)  echo "timeout waiting for daemon connection" ;;
    255) echo "ssh/connection failure" ;;
    *)   echo "see 'man rsync' EXIT VALUES" ;;
  esac
}

# sync one local dir to <bucket>/<sub>, dereferencing symlinks, with the configured client.
# aws uses --only-show-errors so routine per-file uploads don't bury real errors in the log.
s3_sync() {
  local src="$1" sub="$2"
  case "$S3_TOOL" in
    aws)
      aws ${ALC_AWS_PROFILE:+--profile "$ALC_AWS_PROFILE"} \
          ${ALC_S3_ENDPOINT:+--endpoint-url "$ALC_S3_ENDPOINT"} \
          s3 sync "$src" "s3://$ALC_S3_BUCKET/$sub" --follow-symlinks --only-show-errors ;;
    *)
      rclone sync --copy-links --transfers 16 --checkers 32 \
          "$src" "$ALC_S3_REMOTE:$ALC_S3_BUCKET/$sub" ;;
  esac
}

# --- 1. images -> public S3 bucket ---------------------------------------------------------------
if [ -n "${ALC_S3_BUCKET:-}" ] && { [ -n "${ALC_S3_REMOTE:-}" ] || [ "$S3_TOOL" = "aws" ]; }; then
  # 'data' = the daily panel's per-night JSON payloads + station index. They are bucket cargo
  # like the images: tens of MB per stream, immutable per night -- the ssh HTML leg must never
  # carry them. The page fetches same-origin first and falls back to the bucket URL.
  for d in diag ombsens flagex nc data; do
    [ -d "$SITE/$d" ] || continue
    echo "publish: $S3_TOOL sync $d/ -> $ALC_S3_BUCKET/$d"
    err=$(s3_sync "$SITE/$d" "$d" 2>&1 1>/dev/null); lrc=$?
    if [ "$lrc" -ne 0 ]; then
      rc=$lrc
      # A dangling symlink surfaces as "... File does not exist" and forces aws to exit 2; count those
      # separately from real errors (auth, network, permission) so the log distinguishes the two.
      skips=$(printf '%s\n' "$err" | grep -ci 'does not exist' || true)
      detail=$(printf '%s\n' "$err" \
        | grep -viE 'does not exist|Skipping file' \
        | grep -iE 'error|denied|forbidden|fatal|expired|unable|could not|resolve|timed out|connection' \
        | head -4)
      m="S3 sync '$d' failed (exit $lrc)"
      [ "${skips:-0}" -gt 0 ] && m="$m -- $skips missing-source file(s) skipped (dangling links?)"
      echo "publish: ERROR: $m" >&2
      [ -n "$detail" ] && printf '%s\n' "$detail" | sed 's/^/publish:   aws: /' >&2
      FAILED+=("s3:$d(exit $lrc)")
    fi
  done
else
  echo "publish: S3 target not configured (need ALC_S3_BUCKET + ALC_S3_REMOTE, or ALC_S3_TOOL=aws) -> skip images"
fi

# --- 1b. images-in-bucket: prune now-published local diag PNGs to reclaim disk -----------------
# Diagnostic PNGs (cloud/rayleigh) under plots/ AND the classification curtains under classification/
# are both staged into diag/ and synced above; prune both local copies (the classification NetCDF is
# kept -- it is the per-day marker the dashboard index reads).
if [ "${ALC_IMAGES_IN_BUCKET:-0}" = "1" ] && [ "${rc:-0}" -eq 0 ] && [ -n "${ALC_FULLCAL_DIR:-}" ]; then
  n=$(find "$ALC_FULLCAL_DIR" \( -name '*_diag_compact.png' -o -name '*_classification.png' \) -type f 2>/dev/null | wc -l)
  if [ "$n" -gt 0 ]; then
    echo "publish: images-in-bucket -> pruning $n local diag/classification PNGs already synced to the bucket"
    find "$ALC_FULLCAL_DIR" \( -name '*_diag_compact.png' -o -name '*_classification.png' \) -type f -delete
    find "$ALC_FULLCAL_DIR" \( -path '*/plots/*' -o -path '*/classification/*' \) -type d -empty -delete 2>/dev/null || true
  fi
fi

# --- 2. HTML + assets -> web server docroot ------------------------------------------------------
# Exclude the image trees (they live in the bucket), the SQLite index, and the build dotfiles.
# --chmod makes files world-readable so nginx (www-data) can serve them regardless of the build umask.
# NOTE ON --delete: rsync does NOT delete excluded paths on the receiver (that would need
# --delete-excluded), so every --exclude here also PROTECTS that path in the docroot. That is what
# keeps 'intercomparison/' alive: the inter-comparison dashboard is built on a workstation from
# local L1/CAMS archives (inter-comparison_dashboard/build_v3.py + render_v3.py) and uploaded
# straight to the VM, so it never exists in $SITE and an unprotected --delete would wipe it on the
# next daily publish.
if [ -n "${ALC_VM_RSYNC_TARGET:-}" ]; then
  echo "publish: rsync HTML/assets -> $ALC_VM_RSYNC_TARGET"
  err=$(rsync -az --delete --chmod=D755,F644 -e "$VM_SSH" \
    --exclude 'diag/' --exclude 'ombsens/' --exclude 'flagex/' --exclude 'nc/' --exclude 'data/' --exclude 'fullcal_l1_2026/' \
    --exclude 'intercomparison/' \
    --exclude 'calib_index.sqlite' --exclude '.last_build' \
    --exclude '.processed_days' --exclude '.last_success' --exclude '.git*' \
    "$SITE/" "$ALC_VM_RSYNC_TARGET/" 2>&1 1>/dev/null); lrc=$?
  if [ "$lrc" -ne 0 ]; then
    rc=$lrc
    echo "publish: ERROR: rsync HTML -> $ALC_VM_RSYNC_TARGET failed (exit $lrc: $(rsync_meaning "$lrc"))" >&2
    [ -n "$err" ] && printf '%s\n' "$err" | tail -8 | sed 's/^/publish:   rsync: /' >&2
    FAILED+=("rsync(exit $lrc)")
  fi
else
  echo "publish: ALC_VM_RSYNC_TARGET not set -> skipping HTML upload"
fi

if [ "$rc" -ne 0 ]; then
  echo "publish: completed WITH ERRORS (rc=$rc) -- failed leg(s): ${FAILED[*]}" >&2
else
  echo "publish: OK -- all configured legs succeeded"
fi
exit "$rc"
