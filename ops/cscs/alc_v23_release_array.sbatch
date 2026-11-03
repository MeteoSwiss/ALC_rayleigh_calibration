#!/bin/bash
#SBATCH --job-name=alc_v23_arr
#SBATCH --partition=postproc
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_arr_%A_%a.log
#
# Release recompute of calibration v2.3 as a PER-STREAM JOB ARRAY (the queue-friendly twin of
# alc_v23_release.sbatch, for a busy cluster: Balfrin's CPU partitions were down on 2026-11-02, Linard had
# 12 of 14 nodes allocated). One task = one stream over the whole window, written straight into the
# operational layout <key>/<year>/ALC_calibration_*.nc (no chunking, no promotion step), so tasks never
# share a file. Submit with
#     sbatch --array=0-$((N-1))%96 ops/cscs/alc_v23_release_array.sbatch
# where N is the stream count printed by `python ... make_keylist` below (434 for the 2026 census);
# widen or narrow live with `scontrol update JobId=<id> ArrayTaskThrottle=<n>`. House rules: CPU partitions
# only, at most ~6 of the 14 Linard nodes (96 tasks x 4 cores x 16 GB = 6 nodes by memory), code read-only
# in $HOME, run from /scratch. Switches identical to alc_v23_release.sbatch.
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
START=20250101
END=20261101
KEYS=${KEYS_FILE:-/scratch/mch/mhrvo/keys_v23_release.txt}   # KEYS_FILE: a subset keylist (sbatch --export=ALL,KEYS_FILE=...)

if [ ! -f "$KEYS" ]; then
  # the stream list, from the census shipped with the code (LF line endings: a CRLF keylist breaks the task
  # lookup silently -- see the Linard notes)
  python - <<PY > "$KEYS"
import json
c = json.load(open('$CODE/validation/scope_l1_2026_census.json'))
rows = c if isinstance(c, list) else c.get('streams', [])
print('\n'.join(sorted({r['wmo'] + '_' + r['ident'] for r in rows})))
PY
fi
N=$(wc -l < "$KEYS")
IDX=${SLURM_ARRAY_TASK_ID:-0}
KEY=$(sed -n "$((IDX + 1))p" "$KEYS" | tr -d '\r')
if [ -z "$KEY" ]; then echo "no key for index $IDX (N=$N)"; exit 1; fi
mkdir -p "$ALC_FULLCAL_DIR"
echo "task $IDX/$N  stream $KEY  window $START..$END  code $(cat $CODE/GIT_HEAD 2>/dev/null || echo '?')"
python $CODE/scripts/run_network_calibration.py \
  --stream "$KEY" --start $START --end $END \
  --methods rayleigh,cloud --force --ignore-coverage \
  --workers 1
rc=$?
echo "stream $KEY rc=$rc"
if [ "$IDX" = "0" ]; then
  cat > "$ALC_FULLCAL_DIR/RELEASE_PROVENANCE.txt" <<EOF
release: calibration v2.3 (Rayleigh eprof_v2.3 + uncertainty v2.3 + Kalman weighting + polar low-sun fallback; cloud eta tables v2, version 1100)
code: $CODE $(cat $CODE/GIT_HEAD 2>/dev/null || echo unknown)
job: array ${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID} on $(hostname), $(date -u +%Y-%m-%dT%H:%M:%SZ), $N streams
window: $START..$END
switches: ALC_MOLECULAR_METHOD=$ALC_MOLECULAR_METHOD ALC_UNCERTAINTY_MODEL=$ALC_UNCERTAINTY_MODEL ALC_KALMAN_UNCERTAINTY_WEIGHTING=$ALC_KALMAN_UNCERTAINTY_WEIGHTING ALC_POLAR_NIGHT_FALLBACK=$ALC_POLAR_NIGHT_FALLBACK ALC_WV_SPECTRUM=$ALC_WV_SPECTRUM ALC_ETA_TABLES=${ALC_ETA_TABLES:-unset(v2)} ALC_DARK_PROFILE=${ALC_DARK_PROFILE:-unset}
cams: $ALC_CAMS_DIR
EOF
fi
exit $rc
