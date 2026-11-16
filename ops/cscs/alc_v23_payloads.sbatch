#!/bin/bash
#SBATCH --job-name=alc_v23_pay
#SBATCH --partition=postproc
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_pay_%A_%a.log
#
# v2.3.0 cutover: the per-day PAYLOADS of the interactive daily-calibration panel for the whole
# operational v2.3 tree (the release tree with the fix-run streams and the gap fill), as
# alc_v22_payloads.sbatch did for v2.2. The generator RE-RUNS the retrieval in-process to capture
# the curtain and the diagnostic panels, so it gets exactly the release switches: without them the
# panel would show v2.2 diagnostics under v2.3 numbers. It only READS the calibration tree.
# ALC_CLASSIFY stays unset, like the release calibration (no classification screen in the history).
#     sbatch --array=0-3 --export=ALL,END=20261115 ops/cscs/alc_v23_payloads.sbatch
# An array task takes every Nth stream (N = the array size) with its own work dir; the payloads land
# in per-station folders of one OUT tree. CPU partition only, one node per task (4 of Linard's 14).
set -u
case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac
mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc
CODE=${CODE:-/scratch/mch/mhrvo/alc_v23_final_code}
export PYTHONPATH=$CODE
export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_CAMS_DIR='/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS'
export ALC_FULLCAL_DIR=${TREE:-/scratch/mch/mhrvo/E_PROFILE_calout_v23_ops}
export ALC_MOLECULAR_METHOD=eprof_v2.3
export ALC_UNCERTAINTY_MODEL=v2.3
export ALC_KALMAN_UNCERTAINTY_WEIGHTING=1
export ALC_POLAR_NIGHT_FALLBACK=1
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=28800
OUT=${OUT:-/scratch/mch/mhrvo/dash_payloads_v23}
START=${START:-20250101}
: "${END:?set END=YYYYMMDD (the last day of the tree)}"
[ -d "$ALC_FULLCAL_DIR" ] || { echo "REFUSING: $ALC_FULLCAL_DIR does not exist" >&2; exit 1; }
IDX=${SLURM_ARRAY_TASK_ID:-0}
NT=${SLURM_ARRAY_TASK_COUNT:-1}
KEYS=()
i=0
for d in "$ALC_FULLCAL_DIR"/*/; do
  k=$(basename "$d")
  [ -f "$d/${k}_cal.csv" ] || continue
  [ $(( i % NT )) -eq "$IDX" ] && KEYS+=(--key "$k")
  i=$((i + 1))
done
echo "task $IDX/$NT  streams: $(( ${#KEYS[@]} / 2 ))   window: $START..$END   out: $OUT   code $(cat $CODE/GIT_HEAD 2>/dev/null)"
python "$CODE/scripts/build_station_dashboard.py" "${KEYS[@]}" \
  --start "$START" --end "$END" \
  --cal-dir "$ALC_FULLCAL_DIR" --l1-root "$ALC_L1_ROOT" \
  --out "$OUT" --work "$OUT/_work_$IDX" --payloads --no-pages --workers 240
rc=$?
echo -n "stations with payloads: "; ls -d "$OUT"/data/*/ 2>/dev/null | wc -l
echo -n "payload files:          "; find "$OUT/data" -name '*.json' 2>/dev/null | wc -l
echo -n "total size:             "; du -sh "$OUT/data" 2>/dev/null | cut -f1
echo "ALC_V23_PAYLOADS_DONE rc=$rc"
exit $rc
