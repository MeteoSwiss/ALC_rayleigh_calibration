#!/usr/bin/env bash
# Publish the built dashboard to the European Weather Cloud (or any S3 + web host):
#   * the bulky image trees (diag/ ombsens/ flagex/)  -> a public S3 bucket (rclone or aws),
#   * the static HTML + assets                          -> the web server's docroot (rsync over ssh).
# Driven by the ALC_* vars in ops/config.sh. Each leg is skipped (with a note) when its target is not
# configured. Symlinked diagnostics are dereferenced on upload so the real bytes land in the bucket.
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
source "$OPS_DIR/config.sh"

SITE="${ALC_DASHBOARD_DIR:?ALC_DASHBOARD_DIR not set}"
if [ ! -d "$SITE" ]; then
  echo "publish: site dir '$SITE' does not exist -> nothing to publish" >&2
  exit 1
fi

S3_TOOL="${ALC_S3_TOOL:-rclone}"     # rclone | aws
VM_SSH="${ALC_VM_SSH:-ssh}"          # ssh command for the rsync (add -i KEY and a ProxyCommand for a proxy)
rc=0

# sync one local dir to <bucket>/<sub>, dereferencing symlinks, with the configured client.
s3_sync() {
  local src="$1" sub="$2"
  case "$S3_TOOL" in
    aws)
      aws ${ALC_AWS_PROFILE:+--profile "$ALC_AWS_PROFILE"} \
          ${ALC_S3_ENDPOINT:+--endpoint-url "$ALC_S3_ENDPOINT"} \
          s3 sync "$src" "s3://$ALC_S3_BUCKET/$sub" --follow-symlinks --no-progress ;;
    *)
      rclone sync --copy-links --transfers 16 --checkers 32 \
          "$src" "$ALC_S3_REMOTE:$ALC_S3_BUCKET/$sub" ;;
  esac
}

# --- 1. images -> public S3 bucket ---------------------------------------------------------------
if [ -n "${ALC_S3_BUCKET:-}" ] && { [ -n "${ALC_S3_REMOTE:-}" ] || [ "$S3_TOOL" = "aws" ]; }; then
  for d in diag ombsens flagex; do
    [ -d "$SITE/$d" ] || continue
    echo "publish: $S3_TOOL sync $d/ -> $ALC_S3_BUCKET/$d"
    s3_sync "$SITE/$d" "$d" || rc=$?
  done
else
  echo "publish: S3 target not configured (need ALC_S3_BUCKET + ALC_S3_REMOTE, or ALC_S3_TOOL=aws) -> skip images"
fi

# --- 1b. images-in-bucket: prune now-published local diag PNGs to reclaim disk -----------------
if [ "${ALC_IMAGES_IN_BUCKET:-0}" = "1" ] && [ "${rc:-0}" -eq 0 ] && [ -n "${ALC_FULLCAL_DIR:-}" ]; then
  n=$(find "$ALC_FULLCAL_DIR" -name '*_diag_compact.png' -type f 2>/dev/null | wc -l)
  if [ "$n" -gt 0 ]; then
    echo "publish: images-in-bucket -> pruning $n local diag PNGs already synced to the bucket"
    find "$ALC_FULLCAL_DIR" -name '*_diag_compact.png' -type f -delete
    find "$ALC_FULLCAL_DIR" -path '*/plots/*' -type d -empty -delete 2>/dev/null || true
  fi
fi

# --- 2. HTML + assets -> web server docroot ------------------------------------------------------
# Exclude the image trees (they live in the bucket), the SQLite index, and the build dotfiles.
# --chmod makes files world-readable so nginx (www-data) can serve them regardless of the build umask.
if [ -n "${ALC_VM_RSYNC_TARGET:-}" ]; then
  echo "publish: rsync HTML/assets -> $ALC_VM_RSYNC_TARGET"
  rsync -az --delete --chmod=D755,F644 -e "$VM_SSH" \
    --exclude 'diag/' --exclude 'ombsens/' --exclude 'flagex/' --exclude 'fullcal_l1_2026/' \
    --exclude 'calib_index.sqlite' --exclude '.last_build' \
    --exclude '.processed_days' --exclude '.last_success' --exclude '.git*' \
    "$SITE/" "$ALC_VM_RSYNC_TARGET/" || rc=$?
else
  echo "publish: ALC_VM_RSYNC_TARGET not set -> skipping HTML upload"
fi

if [ "$rc" -ne 0 ]; then
  echo "publish: completed WITH ERRORS (rc=$rc)" >&2
fi
exit "$rc"
