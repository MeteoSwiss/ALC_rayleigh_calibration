#!/bin/bash
#SBATCH --job-name=alc_v22_sens
#SBATCH --partition=postproc
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_sens_%j.log
#
# Re-run of the SENSITIVITY half only, after the empty-day fix.
#
# The first add-ons pass produced _sens_cache.npz for 15 of 434 streams. Not a crash in the log and
# not a skip: a day whose L1 slice holds ZERO profiles reached sensitivity_over_period, which
# indexed time[0] -> IndexError. The runner forces logging to CRITICAL inside _process_stream, so
# that exception went to the worker's stdout and never appeared in the batch log -- it presented as
# "sensitivity silently produces nothing". Short windows survived because they happened to contain
# no such day; the full 2025-2026 window did not.
#
# --no-cal --sens: consumes the Kalman constants already on disk and writes only the sensitivity
# product, so it cannot change a single calibration value. Run it AFTER the omb pass finishes
# (--dependency=afterok) so two jobs never write the same stream directory.
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
export PLOTS=0

echo -n "streams with _sens_cache BEFORE: "; ls $ALC_FULLCAL_DIR/*/_sens_cache.npz 2>/dev/null | wc -l

python "$CODE/scripts/run_network_calibration.py" \
  --start 20250101 --end 20260813 \
  --no-cal --sens \
  --per-type 0 --ignore-coverage --force \
  --workers 200
rc=$?

echo -n "streams with _sens_cache AFTER:  "; ls $ALC_FULLCAL_DIR/*/_sens_cache.npz 2>/dev/null | wc -l
echo -n "streams with _sens.csv:          "; ls $ALC_FULLCAL_DIR/*/*_sens.csv 2>/dev/null | wc -l
echo "ALC_V22_SENS_DONE rc=$rc"
exit $rc
