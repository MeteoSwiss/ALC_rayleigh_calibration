#!/bin/bash
#SBATCH --job-name=alc_v22_rel
#SBATCH --partition=postproc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_rel_%j.log
#
# Release recompute of the whole network with eprof_v2.2 AND the CL61 constructor water-vapour
# spectrum -- the archive the operational cutover flips onto (doc/OPERATIONS_v22_cutover.md).
#
# Balfrin house rules honoured here (CSCS "SLURM partitions: Rules and Guidelines"):
#   * postproc is a CPU partition (data analysis). NEVER a GPU partition -- debug, short,
#     short-shared, normal, normal-shared, lowprio, preemptible and production are all GPU, and a
#     CPU-only job there wastes GPUs and money. The guard below refuses to run anywhere else.
#   * One node out of the 15-node CPU pool = ~7 %, far inside the "no more than 50 % of the machine
#     for more than 15 minutes" rule, so no Teams announcement is needed.
#   * No --account / -A: there is no accounting on balfrin.
#   * The job runs from /scratch (--chdir) and reads the code read-only from $HOME via PYTHONPATH:
#     compute jobs must not run from the home directory, which is shared with production.
#   * The whole node is requested because the run peaks around 147 GiB (measured, job 5110089).
#     Balfrin has no memory control, so cores are the proxy for memory: 147/512*256 ~ 74 cores would
#     be the floor, and this job genuinely uses the cores anyway (230 workers).
set -u

case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac

mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc

CODE=$HOME/alc_v22_code          # code ONLY (~8.5 MB): calibration, monitoring, scripts, validation
export PYTHONPATH=$CODE

export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_MOLECULAR_METHOD=eprof_v2.2
# THE release change: the CL61 emission line measured by the manufacturer (Le & O'Connor response to
# RC1, corroborated by the DA10 guide M212895EN-E and by HITRAN line-by-line). The old 910.74/1.0
# model over-corrects water vapour ~9x on this type. CL31/CL51 keep their own entries; the CHM15k at
# 1064 nm is outside the absorption band and unaffected.
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
# Deliberately NOT set: ALC_DARK_PROFILE (operator decision -- keep the network homogeneous).
# The default per-stream timeout kills a long stream AND DISCARDS its output, silently truncating
# the archive; a full-history stream needs hours.
export STREAM_TIMEOUT=28800
export PLOTS=0

OUT=/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel

STREAMS=$(python - <<PY
import json
c = json.load(open('$CODE/validation/scope_l1_2026_census.json'))
rows = c if isinstance(c, list) else c.get('streams', [])
print(','.join(sorted({r['wmo'] + '_' + r['ident'] for r in rows})))
PY
)
echo "streams: $(echo $STREAMS | tr ',' '\n' | wc -l)"
echo "out:     $OUT"
echo "method:  $ALC_MOLECULAR_METHOD   wv_spectrum: $ALC_WV_SPECTRUM"

# --chunks 1 ON PURPOSE. Chunking splits a stream's DATE RANGE across workers, and
# run_streams_parallel deliberately does NOT merge the yearly NetCDF (see its header): each chunk
# leaves its own under parts_nc/chunk<N>/. That is fine for a study, but the operational tree needs
# <key>/<year>/ALC_calibration_<key><year>.nc -- it is the E-PROFILE exchange product and the
# dashboard's download link. With one chunk per stream the NetCDF is written straight into the
# operational layout, and there is still no concurrent writer on a stream (parallelism is ACROSS
# streams), which is the corruption the chunking exists to avoid. It costs wall-clock: the critical
# path becomes the longest single stream's full history, hence 12 h rather than 6 h.
python $CODE/scripts/run_streams_parallel.py \
  --streams "$STREAMS" \
  --start 20250101 --end 20260813 \
  --chunks 1 --methods rayleigh,cloud \
  --cams '/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS' \
  --workers 230 \
  --out "$OUT"
rc=$?
echo "run_streams_parallel rc=$rc"

# Promote the yearly NetCDFs into the OPERATIONAL layout <key>/<year>/ALC_calibration_*.nc.
# run_streams_parallel always parks them under parts_nc/chunk<N>/ and never merges them (its own
# header says so), which is fine for a study but not for the tree the ops server serves and the
# hub reads. With --chunks 1 there is exactly ONE chunk, so its year directories ARE the complete
# files and this is a pure move (they are hard links to real files, not symlinks into _parts).
# If a second chunk ever appears the year exists in several parts and picking chunk0 would silently
# truncate the archive -- so refuse that stream loudly instead.
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

# Provenance beside the archive: which code, which switches, which window. Deliberately a sidecar
# file rather than a new NetCDF version code -- the yearly NetCDF is the exchange product read by
# consumers outside this repository, and a changed version code there could silently drop every
# calibration from a downstream filter.
cat > "$OUT/RUN_PROVENANCE.json" <<PY
{
  "job_id": "${SLURM_JOB_ID:-}",
  "partition": "${SLURM_JOB_PARTITION:-}",
  "molecular_method": "eprof_v2.2",
  "wv_spectrum": {"CL61": [910.55, 0.188]},
  "dark_profile": null,
  "window": {"start": "20250101", "end": "20260813"},
  "cams": "CAMS_Monthly_04 (0.4 deg) + CAMS dailies",
  "code": "$CODE",
  "git_describe": "$(cat $CODE/.release_commit 2>/dev/null || echo unknown)"
}
PY
echo "ALC_V22_RELEASE_DONE rc=$rc"
exit $rc
