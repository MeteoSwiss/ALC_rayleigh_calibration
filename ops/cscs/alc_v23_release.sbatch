#!/bin/bash
#SBATCH --job-name=alc_v23_rel
#SBATCH --partition=postproc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v23_rel_%j.log
#
# Release recompute of the whole network for calibration v2.3 (doc/v2.3/network_rerun_plan_2026-11.md):
#   * Rayleigh eprof_v2.3 (measured tier-2 noise + vertical-consistency test on the CHM15k, noise only
#     and a 5 km tier-2 window cap on the CL61, v2.2 elsewhere; calibration_version 230 / 231 low-sun)
#   * uncertainty budget v2.3 with the Kalman weighted by it
#   * polar-summer low-sun fallback (flag 0.25, version 231)
#   * liquid-cloud eta tables v2 (PVC a_G 3.3 um, alpha 30 /km, S 18.8 sr; cloud version 1100)
#   * Payerne CHM15k A on tier 1 only (options.json station_molecular_params)
# The CL61 water-vapour line 910.55 / 0.188 nm is the code default since 2026-11-02 (exported below
# anyway, so a stale code copy cannot silently fall back to the retired model).
#
# Same house rules as alc_v22_release.sbatch: CPU partition only (the guard refuses anything else),
# one postproc node (~7 % of the pool), code read-only from $HOME via PYTHONPATH, run from /scratch,
# never two jobs writing one output tree, --chunks 1 so the yearly NetCDFs land in the operational
# layout (one chunk per stream, parallelism across streams).
#
# Code: rsync the operational subset of the repository to $HOME/alc_v23_code before submitting
# (calibration monitoring ops scripts validation rayleigh_availability tests options.json
# pyproject.toml, no __pycache__), and write the commit hash to $HOME/alc_v23_code/GIT_HEAD.
# Data: L1 /scratch/mch/mhrvo/E-PROFILE_L1 (network to 2026-08-13 at the time of the first run,
# Payerne to 2026-11-01), CAMS 0.4 deg monthlies to 2026-06 + daily files to 2026-11-01.
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
export ALC_MOLECULAR_METHOD=eprof_v2.3
export ALC_UNCERTAINTY_MODEL=v2.3
export ALC_KALMAN_UNCERTAINTY_WEIGHTING=1
export ALC_POLAR_NIGHT_FALLBACK=1
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
# Deliberately NOT set: ALC_DARK_PROFILE (network homogeneity, operator decision 2026-08-17) and
# ALC_ETA_TABLES (unset = the v2 tables; pvc55 / hopkin_legacy are comparison runs only).
export STREAM_TIMEOUT=28800
export PLOTS=0

OUT=/scratch/mch/mhrvo/E_PROFILE_calout_v23_rel
START=20250101
END=20261101

STREAMS=$(python - <<PY
import json
c = json.load(open('$CODE/validation/scope_l1_2026_census.json'))
rows = c if isinstance(c, list) else c.get('streams', [])
print(','.join(sorted({r['wmo'] + '_' + r['ident'] for r in rows})))
PY
)
echo "streams: $(echo $STREAMS | tr ',' '\n' | wc -l)"
echo "out:     $OUT   window $START..$END"
echo "method:  $ALC_MOLECULAR_METHOD  uncertainty $ALC_UNCERTAINTY_MODEL  kalman-weighting $ALC_KALMAN_UNCERTAINTY_WEIGHTING  polar $ALC_POLAR_NIGHT_FALLBACK"
echo "code:    $CODE ($(cat $CODE/GIT_HEAD 2>/dev/null || echo 'no GIT_HEAD'))"

python $CODE/scripts/run_streams_parallel.py \
  --streams "$STREAMS" \
  --start $START --end $END \
  --chunks 1 --methods rayleigh,cloud \
  --cams '/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS' \
  --workers 230 \
  --out "$OUT"
rc=$?
echo "run_streams_parallel rc=$rc"

# Promote the yearly NetCDFs into the OPERATIONAL layout <key>/<year>/ALC_calibration_*.nc (one chunk
# per stream, so this is a pure move; a stream with several chunks is refused loudly).
promoted=0; skipped=0
for sdir in "$OUT"/*/; do
  pn="${sdir}parts_nc"
  [ -d "$pn" ] || continue
  n=$(ls -d "$pn"/chunk* 2>/dev/null | wc -l)
  if [ "$n" -ne 1 ]; then
    echo "SKIP $(basename "$sdir"): $n chunks under parts_nc -- cannot promote safely"
    skipped=$((skipped + 1)); continue
  fi
  for y in "$pn"/chunk0/*/; do
    [ -d "$y" ] || continue
    mv -n "$y" "$sdir" && promoted=$((promoted + 1))
  done
  rmdir "$pn/chunk0" "$pn" 2>/dev/null || true
done
echo "netcdf promoted: $promoted year-dirs, $skipped streams skipped"
[ "$skipped" -gt 0 ] && rc=1

# Provenance sidecar beside the archive.
cat > "$OUT/RELEASE_PROVENANCE.txt" <<EOF
release: calibration v2.3 (Rayleigh eprof_v2.3 + uncertainty v2.3 + Kalman weighting + polar low-sun fallback; cloud eta tables v2, version 1100)
code: $CODE $(cat $CODE/GIT_HEAD 2>/dev/null || echo unknown)
job: $SLURM_JOB_ID on $(hostname), $(date -u +%Y-%m-%dT%H:%M:%SZ)
window: $START..$END   streams: $(echo $STREAMS | tr ',' '\n' | wc -l)
switches: ALC_MOLECULAR_METHOD=$ALC_MOLECULAR_METHOD ALC_UNCERTAINTY_MODEL=$ALC_UNCERTAINTY_MODEL ALC_KALMAN_UNCERTAINTY_WEIGHTING=$ALC_KALMAN_UNCERTAINTY_WEIGHTING ALC_POLAR_NIGHT_FALLBACK=$ALC_POLAR_NIGHT_FALLBACK ALC_WV_SPECTRUM=$ALC_WV_SPECTRUM ALC_ETA_TABLES=${ALC_ETA_TABLES:-unset(v2)} ALC_DARK_PROFILE=${ALC_DARK_PROFILE:-unset}
cams: /scratch/mch/mhrvo/CAMS_Monthly_04 ; /scratch/mch/mhrvo/CAMS
rc: $rc
EOF
echo "ALC_V23_RELEASE_DONE rc=$rc"
exit $rc
