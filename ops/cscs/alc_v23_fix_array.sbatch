#!/bin/bash
#SBATCH --job-name=alc_v23_fix
#SBATCH --partition=postproc
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_fix_%A_%a.log
#
# v2.3.0 fix run (2026-11-16): recompute the streams that the unit-swap and Kalman changes of commit
# b9212e3 move, over the release window, calibration then add-ons, into a SEPARATE tree that replaces
# those streams in the release tree afterwards:
#   * the four CL31 -> CL51 swaps (per-day Vaisala type + Kalman restart): Plzen 0-20000-0-11450_A,
#     Tusimice 0-20000-0-11438_A, Libus 0-20000-0-11520_A, Usti nad Labem 0-20000-0-11502_A;
#   * Payerne CL31 0-20000-0-06610_B (Kalman restart at the 2026-07-07 optical block);
#   * Oslo CHM15k 0-20000-0-01492_A (low-sun nights out of its Kalman).
# The OmB / sensitivity per-day caches carry the Kalman C_L, so the add-ons are forced too.
#     sbatch --array=0-5 --export=ALL,KEYS_FILE=/scratch/mch/mhrvo/keys_v23_fix.txt ops/cscs/alc_v23_fix_array.sbatch
# Switches identical to alc_v23_release_array.sbatch. CPU partitions only.
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
export ALC_FULLCAL_DIR=${OUT_DIR:-/scratch/mch/mhrvo/E_PROFILE_calout_v23_fix}
export ALC_MOLECULAR_METHOD=eprof_v2.3
export ALC_UNCERTAINTY_MODEL=v2.3
export ALC_KALMAN_UNCERTAINTY_WEIGHTING=1
export ALC_POLAR_NIGHT_FALLBACK=1
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=0
export PLOTS=0
START=20250101
END=20261101
KEYS=${KEYS_FILE:-/scratch/mch/mhrvo/keys_v23_fix.txt}
IDX=${SLURM_ARRAY_TASK_ID:-0}
KEY=$(sed -n "$((IDX + 1))p" "$KEYS" | tr -d '\r')
[ -z "$KEY" ] && { echo "no key for index $IDX"; exit 1; }
mkdir -p "$ALC_FULLCAL_DIR"
echo "task $IDX  stream $KEY  window $START..$END  code $(cat $CODE/GIT_HEAD 2>/dev/null || echo '?')"
python $CODE/scripts/run_network_calibration.py \
  --stream "$KEY" --start $START --end $END \
  --methods rayleigh,cloud --force --ignore-coverage --workers 1
rc=$?
echo "stream $KEY cal rc=$rc"
[ $rc -ne 0 ] && exit $rc
python $CODE/scripts/run_network_calibration.py \
  --stream "$KEY" --start $START --end $END \
  --no-cal --sens --omb --ignore-coverage --force --workers 1
rc=$?
echo "stream $KEY addons rc=$rc"
exit $rc
