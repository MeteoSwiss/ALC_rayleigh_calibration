#!/bin/bash
#SBATCH --job-name=alc_v23_gap
#SBATCH --partition=postproc
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_gap_%A_%a.log
#
# v2.3.0 cutover gap fill: the days between the end of the release window (2026-11-01) and the switch,
# per stream, exactly as the daily run does them (ops/ops_daily.py: --methods rayleigh,cloud --force --sens
# --omb --ignore-coverage), merged into the operational v2.3 tree (the release tree with the fix-run
# streams in place). Classification is not run: its history is copied from the v2.2 tree at the switch.
#     sbatch --array=0-433%96 --export=ALL,GAP_START=20261102,GAP_END=20261115 ops/cscs/alc_v23_gap_array.sbatch
# A later top-up uses the next GAP_START. CPU partitions only.
set -u
case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac
mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc
for v in OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS; do export $v=1; done

CODE=${CODE:-/scratch/mch/mhrvo/alc_v23_final_code}
export PYTHONPATH=$CODE
export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_CAMS_DIR='/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS'
export ALC_FULLCAL_DIR=${OUT_DIR:-/scratch/mch/mhrvo/E_PROFILE_calout_v23_ops}
export ALC_MOLECULAR_METHOD=eprof_v2.3
export ALC_UNCERTAINTY_MODEL=v2.3
export ALC_KALMAN_UNCERTAINTY_WEIGHTING=1
export ALC_POLAR_NIGHT_FALLBACK=1
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=0
export PLOTS=0
: "${GAP_START:?set GAP_START=YYYYMMDD}" "${GAP_END:?set GAP_END=YYYYMMDD}"
KEYS=${KEYS_FILE:-/scratch/mch/mhrvo/keys_v23_release.txt}
IDX=${SLURM_ARRAY_TASK_ID:-0}
KEY=$(sed -n "$((IDX + 1))p" "$KEYS" | tr -d '\r')
[ -z "$KEY" ] && { echo "no key for index $IDX"; exit 1; }
echo "task $IDX  stream $KEY  gap $GAP_START..$GAP_END  code $(cat $CODE/GIT_HEAD 2>/dev/null || echo '?')"
python $CODE/scripts/run_network_calibration.py \
  --stream "$KEY" --start $GAP_START --end $GAP_END \
  --methods rayleigh,cloud --force --sens --omb --ignore-coverage --workers 1
rc=$?
echo "stream $KEY gap rc=$rc"
exit $rc
