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
# The SAME two release switches as alc_v22_release.sbatch, and they are NOT redundant under
# --no-cal: that flag stops the stored constants being recomputed, but the sensitivity product
# re-runs the retrieval under perturbed lidar ratio / altitude shifts, and OmB re-derives the
# 910 nm water-vapour correction. Without these the add-ons would be computed with eprof_v2 and the
# superseded CL61 spectrum while the archive beside them is v2.2 -- a mixed-vintage tree in the one
# place the vintage gate does not look.
export ALC_MOLECULAR_METHOD=eprof_v2.2
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export STREAM_TIMEOUT=28800
# Root cause of the 14/433 first run (job 5122633) -- NOT concurrency: a day whose L1 slice holds
# ZERO profiles reached numpy's char.replace, whose empty-array buffersizes.max() raises ValueError;
# one such day aborted the stream's whole OmB pass AND (single shared try) its sensitivity pass.
# Fixed at three layers in run_network_calibration.py: empty-day skip in the OmB loop, an empty-time
# guard in _const_per_profile, and omb/sens in SEPARATE trys with full tracebacks printed.
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
