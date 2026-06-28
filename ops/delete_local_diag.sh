#!/bin/bash
cd /data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code || exit 1
source ops/config.sh >/dev/null 2>&1
EP="$ALC_S3_ENDPOINT"; BK="s3://eprofile-alc-dashboard/diag"
del_pngs=0; freed=0; skipped=0
echo "START $(date -u +%FT%TZ)  before: $(df -h /data/zue | tail -1 | awk '{print $4" free ("$5")"}')"
for d in "$ALC_FULLCAL_DIR"/*/; do
  k=$(basename "$d"); [ -d "$d/plots" ] || continue
  L=$(find "$d/plots" -name "*_diag_compact.png" -type f 2>/dev/null | wc -l)
  [ "$L" -gt 0 ] || continue
  B=$(aws --profile ewc --endpoint-url "$EP" s3 ls "$BK/$k/" 2>/dev/null | grep -c '\.png')
  if [ "$B" -ge "$L" ]; then
    find "$d/plots" -name "*_diag_compact.png" -type f -delete
    find "$d/plots" -type d -empty -delete 2>/dev/null
    freed=$((freed+1)); del_pngs=$((del_pngs+L))
  else
    echo "SKIP $k local=$L bucket=$B"; skipped=$((skipped+1))
  fi
done
echo "IMGDEL_DONE deleted_pngs=$del_pngs freed_stations=$freed skipped=$skipped"
echo "after: $(df -h /data/zue | tail -1 | awk '{print $4" free ("$5")"}')  $(date -u +%FT%TZ)"
