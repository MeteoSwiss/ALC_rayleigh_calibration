#!/bin/bash
#SBATCH --job-name=alc_v22_fin
#SBATCH --partition=postproc
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_fin_%j.log
#
# Step 4: turn the finished archive + payloads into a published site. Submitted with
#   --dependency=afterok:<addons>:<payloads>
# so it starts only once BOTH producers have succeeded, and runs unattended overnight.
#
#   1. stamp the `version` column into the archive (it predates the column; see the script)
#   2. full site build, bucket mode, into the payload tree so pages sit beside their data/
#   3. the site test tier, as a gate that is REPORTED but does not block the data upload
#   4. push data/ to the object store
#   5. tar the HTML for the operator's relay to the web VM
#
# Deliberately NOT done here:
#   * ombsens/ upload -- those keys are live NOW under the v2.0 site, and replacing them before the
#     v2.2 HTML lands would leave the public page showing v2.2 OmB plots beside v2.0 constants. It
#     goes with the HTML flip instead (publish_v22_payloads_to_ewc.sh ombsens).
#   * the HTML upload -- the web VM's ssh key is on the operator's workstation, not here.
#
# postproc compute nodes DO have internet here (verified: HTTP 200 to the ECMWF endpoint), which is
# what lets the upload run inside a batch job rather than on the login node.
set -u

case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac

source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc

CODE=$HOME/alc_v22_code
export PYTHONPATH=$CODE
export PATH="$PATH:/store_new/mch/msclim/share/software/bin"
TREE=/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel
SITE=/scratch/mch/mhrvo/dash_payloads
export ALC_IMG_BASE_URL=https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/

echo "=============== 1. version stamp"
python "$CODE/ops/cscs/backfill_version_column.py" "$TREE" || exit 1

echo "=============== 2. site build (full, no --changed-only)"
t0=$SECONDS
python "$CODE/scripts/build_dashboard.py" --fullcal "$TREE" --out "$SITE" --workers 96
rc_build=$?
echo "build rc=$rc_build in $(( (SECONDS - t0) / 60 )) min"
echo -n "station pages: "; ls "$SITE"/stations/*.html 2>/dev/null | wc -l
[ $rc_build -ne 0 ] && { echo "ALC_V22_FINALIZE_DONE rc=$rc_build (build failed)"; exit $rc_build; }

echo "=============== 3. site test tier (reported, non-blocking for the data leg)"
ALC_SITE_DIR="$SITE" python -m pytest "$CODE/tests/test_dashboard_site.py" -q --no-header 2>&1 | tail -15
rc_test=${PIPESTATUS[0]}
echo "site tests rc=$rc_test"

echo "=============== 4. payload upload -> object store"
aws configure set default.s3.max_concurrent_requests 64 --profile ewc 2>/dev/null || true
aws configure set default.s3.max_queue_size 20000 --profile ewc 2>/dev/null || true
t0=$SECONDS
aws --profile ewc --endpoint-url https://object-store.os-api.cci2.ecmwf.int \
    s3 sync "$SITE/data" s3://eprofile-alc-dashboard/data --only-show-errors --size-only
rc_s3=$?
echo "s3 data rc=$rc_s3 in $(( (SECONDS - t0) / 60 )) min"
echo -n "local payloads:  "; find "$SITE/data" -type f | wc -l
echo -n "remote payloads: "; aws --profile ewc --endpoint-url https://object-store.os-api.cci2.ecmwf.int \
    s3 ls s3://eprofile-alc-dashboard/data/ --recursive --summarize 2>/dev/null \
    | awk '/Total Objects:/{print $3}'

echo "=============== 5. HTML bundle for the operator's relay"
tar -C "$SITE" -czf /scratch/mch/mhrvo/alc_v22_html.tgz \
    --exclude=data --exclude=diag --exclude=ombsens --exclude=nc --exclude=flagex . 2>/dev/null
echo -n "html bundle: "; du -h /scratch/mch/mhrvo/alc_v22_html.tgz 2>/dev/null | cut -f1

echo "ALC_V22_FINALIZE_DONE build=$rc_build tests=$rc_test s3=$rc_s3"
exit 0
