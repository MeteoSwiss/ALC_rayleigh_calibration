#!/usr/bin/env bash
# Publish the built dashboard to the European Weather Cloud (or any S3 + web host):
#   * the bulky image trees (diag/ ombsens/ flagex/)  -> a public S3 bucket (rclone),
#   * the static HTML + assets                          -> the web server's docroot (rsync over ssh).
# Driven by the ALC_* vars in ops/config.sh. Each leg is skipped (with a note) when its target is not
# configured, so a partial setup still works. Symlinked diagnostics are dereferenced on upload
# (--copy-links) so the real bytes land in the bucket. The bucket base URL must also be baked into the
# HTML at build time via ALC_IMG_BASE_URL so the pages point their images at the bucket.
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

rc=0

# --- 1. images -> public S3 bucket ---------------------------------------------------------------
if [ -n "${ALC_S3_REMOTE:-}" ] && [ -n "${ALC_S3_BUCKET:-}" ]; then
  for d in diag ombsens flagex; do
    if [ -d "$SITE/$d" ]; then
      echo "publish: rclone sync $d/ -> $ALC_S3_REMOTE:$ALC_S3_BUCKET/$d"
      rclone sync --copy-links --transfers 16 --checkers 32 \
        "$SITE/$d" "$ALC_S3_REMOTE:$ALC_S3_BUCKET/$d" || rc=$?
    fi
  done
else
  echo "publish: ALC_S3_REMOTE/ALC_S3_BUCKET not set -> skipping image upload"
fi

# --- 2. HTML + assets -> web server docroot ------------------------------------------------------
# Exclude the image trees (they live in the bucket), the SQLite index, and the build dotfiles.
if [ -n "${ALC_VM_RSYNC_TARGET:-}" ]; then
  echo "publish: rsync HTML/assets -> $ALC_VM_RSYNC_TARGET"
  # --chmod makes files world-readable on the web server (nginx runs as www-data); without it the
  # build host's umask can land non-readable files -> nginx 403/404.
  rsync -az --delete --chmod=D755,F644 \
    --exclude 'diag/' --exclude 'ombsens/' --exclude 'flagex/' \
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
