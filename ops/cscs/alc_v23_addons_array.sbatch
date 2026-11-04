#!/bin/bash
#SBATCH --job-name=alc_v23_adda
#SBATCH --partition=postproc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_adda_%A_%a.log
#
# OmB + sensitivity add-ons for the v2.3 release, one stream per array task, on top of the archive written by
# alc_v23_release_array.sbatch. Submit after that array has finished (or with --dependency=afterany:<id>):
#     sbatch --array=0-$((N-1))%64 ops/cscs/alc_v23_addons_array.sbatch
# --mem=24G: OmB loads a whole month of native L1 (CL61 ~9-12 GB). --no-cal only adds files (caches the
# daily runner's regression guard needs).
set -u
case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac
mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc
for v in OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS; do export $v=1; done
CODE=$HOME/alc_v23_code
export PYTHONPATH=$CODE
export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_CAMS_DIR='/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS'
export ALC_FULLCAL_DIR=/scratch/mch/mhrvo/E_PROFILE_calout_v23_rel
export ALC_MOLECULAR_METHOD=eprof_v2.3
export ALC_UNCERTAINTY_MODEL=v2.3
export ALC_KALMAN_UNCERTAINTY_WEIGHTING=1
export ALC_POLAR_NIGHT_FALLBACK=1
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=0
export PLOTS=0
KEYS=${KEYS_FILE:-/scratch/mch/mhrvo/keys_v23_release.txt}   # KEYS_FILE: a subset keylist (sbatch --export=ALL,KEYS_FILE=...)
# ADDONS_FORCE: "--force" (default) recomputes every day; export ADDONS_FORCE="" for a rerun that keeps the
# per-day caches a timed-out task already wrote (sbatch --time=20:00:00 for the long Dutch/DWD streams).
IDX=${SLURM_ARRAY_TASK_ID:-0}
KEY=$(sed -n "$((IDX + 1))p" "$KEYS" | tr -d '\r')
[ -z "$KEY" ] && { echo "no key for index $IDX"; exit 1; }
[ -d "$ALC_FULLCAL_DIR/$KEY" ] || { echo "no archive for $KEY -- nothing to add"; exit 0; }
echo "task $IDX  stream $KEY"
python $CODE/scripts/run_network_calibration.py \
  --stream "$KEY" --start 20250101 --end 20261101 \
  --no-cal --sens --omb --ignore-coverage ${ADDONS_FORCE---force} --workers 1
rc=$?
echo "stream $KEY addons rc=$rc"
exit $rc
