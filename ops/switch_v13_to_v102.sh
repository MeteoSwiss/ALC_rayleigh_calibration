#!/usr/bin/env bash
# One-shot coordinated switch: rename the v13 Rayleigh-test data dir to v1.0.2 and repoint
# everything that references it (ops/config.sh ALC_V13_DIR + the dev rayleigh_calibration
# options.json folder_output that the 18:18 cron writes through), then rebuild + publish.
#
# MUST run when the dev `rayleigh_calibration.main` recalc is NOT running -- it writes into
# the dir, and renaming it mid-run would split / corrupt its output. The guard below refuses
# to proceed while that process is alive. Safe to re-run (idempotent).
set -uo pipefail
OLD=/data/zue/E_PROFILE/ALC/Calibration/rayleigh-test-v13
NEW=/data/zue/E_PROFILE/ALC/Calibration/rayleigh-v1.0.2
REPO=/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code
DEVOPT=/proj/pay/E-PROFILE/Calibration_codes/dev/rayleigh_calibration/options.json

if pgrep -f 'rayleigh_calibration\.main' >/dev/null; then
  echo "REFUSING: rayleigh_calibration.main is still running (it writes $OLD)."
  echo "          Re-run this script once that recalc has finished."
  exit 3
fi

# 1. rename the data dir
if [ -d "$OLD" ] && [ ! -d "$NEW" ]; then
  mv "$OLD" "$NEW" && echo "[1] renamed dir -> $NEW"
elif [ -d "$NEW" ]; then
  echo "[1] dir already renamed ($NEW exists)"
else
  echo "[1] WARN: neither $OLD nor $NEW found as expected"
fi

# 2. ops/config.sh ALC_V13_DIR -> NEW path
if grep -q 'rayleigh-test-v13' "$REPO/ops/config.sh"; then
  sed -i 's#rayleigh-test-v13#rayleigh-v1.0.2#g' "$REPO/ops/config.sh" && echo "[2] config.sh ALC_V13_DIR updated"
else
  echo "[2] config.sh already updated"
fi

# 3. dev rayleigh_calibration options.json folder_output -> NEW path (the 18:18 cron writes here)
if grep -q 'rayleigh-test-v13' "$DEVOPT"; then
  sed -i 's#rayleigh-test-v13#rayleigh-v1.0.2#g' "$DEVOPT" && echo "[3] dev options.json folder_output updated"
else
  echo "[3] dev options.json already updated"
fi

# 4. rebuild dashboard + publish (overlay now read from the renamed dir, labelled v1.0.2)
echo "[4] rebuild + publish ..."
cd "$REPO"
# shellcheck source=/dev/null
source ops/config.sh
# shellcheck source=/dev/null
[ -f "$ALC_VENV/bin/activate" ] && source "$ALC_VENV/bin/activate"
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR" \
  ${ALC_L2_DIR:+--l2dir "$ALC_L2_DIR"} ${ALC_OPCOEFF_CSV:+--opcoeff "$ALC_OPCOEFF_CSV"}
bash ops/publish.sh
echo "[done] v13 -> v1.0.2 switch complete."
