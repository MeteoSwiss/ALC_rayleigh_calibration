#!/bin/bash
#SBATCH --job-name=alc_v23_add
#SBATCH --partition=postproc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_add_%j.log
#
# Step 2 of the v2.3 release recompute: OmB and sensitivity products on top of the archive written by
# alc_v23_release.sbatch (submit with --dependency=afterok:<release job id>). --no-cal consumes the
# constants already on disk and only adds files (the _omb_cache.npz / _sens_cache.npz the daily runner's
# regression guard needs). Same switches as the release job: the sensitivity product re-runs the retrieval
# under perturbations and OmB re-derives the 910 nm water-vapour correction. See alc_v22_addons.sbatch for
# the full rationale and the empty-day fixes.
set -u
case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac
mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc
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
export STREAM_TIMEOUT=28800
export PLOTS=0
if [ ! -d "$ALC_FULLCAL_DIR" ]; then
  echo "REFUSING: $ALC_FULLCAL_DIR does not exist -- run alc_v23_release.sbatch first" >&2
  exit 1
fi
echo "tree:    $ALC_FULLCAL_DIR"
echo "streams: $(ls -d $ALC_FULLCAL_DIR/*/ 2>/dev/null | wc -l)"
python $CODE/scripts/run_network_calibration.py \
  --start 20250101 --end 20261101 \
  --no-cal --sens --omb \
  --per-type 0 --ignore-coverage --force \
  --workers 200
rc=$?
echo "addons rc=$rc"
echo -n "streams with _sens_cache.npz: "; ls $ALC_FULLCAL_DIR/*/_sens_cache.npz 2>/dev/null | wc -l
echo -n "streams with _omb_cache.npz:  "; ls $ALC_FULLCAL_DIR/*/_omb_cache.npz 2>/dev/null | wc -l
echo "ALC_V23_ADDONS_DONE rc=$rc"
exit $rc
