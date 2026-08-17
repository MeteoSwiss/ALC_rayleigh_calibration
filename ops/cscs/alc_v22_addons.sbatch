#!/bin/bash
#SBATCH --job-name=alc_v22_add
#SBATCH --partition=postproc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_add_%j.log
#
# Step 2 of the release recompute: the OmB and sensitivity products, added ON TOP of the archive
# written by alc_v22_release.sbatch. Run it only after that job reports ALC_V22_RELEASE_DONE rc=0.
#
# Why a separate job rather than --sens --omb on the first one:
#  * these read CAMS per day and are the slow half; a time-limit kill would otherwise take the
#    CALIBRATION archive down with them, and this way each half is restartable on its own;
#  * --no-cal makes this step consume the Kalman constants already on disk, so it cannot change a
#    single calibration value -- it only adds files.
#
# It is also the step that writes _omb_cache.npz / _sens_cache.npz. Those caches are NOT optional:
# the daily runner's regression guard compares the cache against the CSV and SKIPS the product
# rather than rebuild it, so an archive transferred without them freezes OmB and sensitivity
# permanently on the operational server.
#
# Balfrin house rules: CPU partition only (postproc), one node of the 15-node pool (~7 %, inside the
# 50 % rule), no --account (no accounting on balfrin), run from /scratch with the code read-only in
# $HOME. See ops/cscs/alc_v22_release.sbatch for the full rationale.
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
export STREAM_TIMEOUT=28800
export PLOTS=0

if [ ! -d "$ALC_FULLCAL_DIR" ]; then
  echo "REFUSING: $ALC_FULLCAL_DIR does not exist -- run alc_v22_release.sbatch first" >&2
  exit 1
fi

echo "tree:    $ALC_FULLCAL_DIR"
echo "streams: $(ls -d $ALC_FULLCAL_DIR/*/ 2>/dev/null | wc -l)"

python $CODE/scripts/run_network_calibration.py \
  --start 20250101 --end 20260813 \
  --no-cal --sens --omb \
  --per-type 0 --ignore-coverage --force \
  --workers 200
rc=$?
echo "addons rc=$rc"

echo -n "streams with _sens_cache.npz: "; ls $ALC_FULLCAL_DIR/*/_sens_cache.npz 2>/dev/null | wc -l
echo -n "streams with _omb_cache.npz:  "; ls $ALC_FULLCAL_DIR/*/_omb_cache.npz 2>/dev/null | wc -l
echo -n "streams with _sens.csv:       "; ls $ALC_FULLCAL_DIR/*/*_sens.csv 2>/dev/null | wc -l
echo -n "streams with _omb.csv:        "; ls $ALC_FULLCAL_DIR/*/*_omb.csv 2>/dev/null | wc -l
echo "ALC_V22_ADDONS_DONE rc=$rc"
exit $rc
