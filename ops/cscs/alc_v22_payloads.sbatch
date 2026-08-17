#!/bin/bash
#SBATCH --job-name=alc_v22_pay
#SBATCH --partition=postproc
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_pay_%j.log
#
# Step 3 of the release: the per-day PAYLOADS that feed the interactive daily-calibration panel,
# for the full 2025-01-01 -> 2026-08-13 archive at full resolution.
#
# Why this is a third job and not a flag on the other two:
#  * it does not write into the calibration archive at all -- it only reads it and writes JSON into
#    a separate tree, so it can be killed and restarted without touching a single constant;
#  * it is the only step whose output is sized in tens of GB, so it must be restartable on its own.
#
# ⚠ The two release switches below are NOT optional here even though this "is not a calibration":
# the payload generator RE-RUNS the retrieval in-process to capture the curtain and the diagnostic
# panels. Without them the panels would be eprof_v2 + the superseded CL61 spectrum sitting under
# v2.2 numbers -- a mixed-vintage page, in the one place no gate looks.
#
# Balfrin house rules: CPU partition only (postproc), one node of the 15-node pool, no --account,
# run from /scratch with the code read-only in $HOME.
set -u

case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac

mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc

CODE=$HOME/alc_v22_code
export PYTHONPATH=$CODE
export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_CAMS_DIR='/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS'
export ALC_FULLCAL_DIR=/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel
export ALC_MOLECULAR_METHOD=eprof_v2.2
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=28800

OUT=/scratch/mch/mhrvo/dash_payloads
START=20250101
END=20260813

if [ ! -d "$ALC_FULLCAL_DIR" ]; then
  echo "REFUSING: $ALC_FULLCAL_DIR does not exist" >&2
  exit 1
fi

# Stream list = every directory of the release tree. That tree was itself built from the union of
# the live census and the pre-existing archive, so taking it here inherits that union -- a station
# commissioned since the census snapshot is not dropped.
KEYS=()
for d in "$ALC_FULLCAL_DIR"/*/; do
  k=$(basename "$d")
  [ -f "$d/${k}_cal.csv" ] || continue
  KEYS+=(--key "$k")
done
echo "streams: $(( ${#KEYS[@]} / 2 ))   window: $START..$END   out: $OUT"
echo "method:  $ALC_MOLECULAR_METHOD   wv: $ALC_WV_SPECTRUM"

python "$CODE/scripts/build_station_dashboard.py" \
  "${KEYS[@]}" \
  --start "$START" --end "$END" \
  --cal-dir "$ALC_FULLCAL_DIR" --l1-root "$ALC_L1_ROOT" \
  --out "$OUT" --payloads --no-pages --workers 200
rc=$?
echo "payloads rc=$rc"

echo -n "stations with payloads: "; ls -d "$OUT"/data/*/ 2>/dev/null | wc -l
echo -n "payload files:          "; find "$OUT/data" -name '*.json' 2>/dev/null | wc -l
echo -n "total size:             "; du -sh "$OUT/data" 2>/dev/null | cut -f1
echo "ALC_V22_PAYLOADS_DONE rc=$rc"
exit $rc
