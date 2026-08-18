#!/bin/bash
#SBATCH --job-name=alc_v22_up
#SBATCH --partition=postproc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_up_%j.log
#
# The payload upload, pulled OUT of the finalize job so it can run in parallel with the site test
# tier rather than behind it.
#
# Safe to start early and safe to run twice: the payload tree is written once by the payload job
# (COMPLETED before this is submitted) and never touched again, and `aws s3 sync --size-only` is
# idempotent -- finalize's own sync afterwards then finds the objects already present and becomes a
# near no-op. The test tier walks tens of thousands of 400 kB payloads, so serialising a 102 GB
# upload behind it would have put the release into the middle of the working day for no reason.
set -u

case "${SLURM_JOB_PARTITION:-unset}" in
  postproc|pp-short|pp-serial|pp-long) ;;
  *) echo "REFUSING: '${SLURM_JOB_PARTITION:-unset}' is not a CPU partition" >&2; exit 1 ;;
esac

export PATH="$PATH:/store_new/mch/msclim/share/software/bin"
SITE=/scratch/mch/mhrvo/dash_payloads
BUCKET=eprofile-alc-dashboard
ENDPOINT=https://object-store.os-api.cci2.ecmwf.int

aws configure set default.s3.max_concurrent_requests 64 --profile ewc 2>/dev/null || true
aws configure set default.s3.max_queue_size 20000 --profile ewc 2>/dev/null || true

echo -n "local payloads: "; find "$SITE/data" -type f | wc -l
echo -n "local size:     "; du -sh "$SITE/data" | cut -f1
t0=$SECONDS
aws --profile ewc --endpoint-url "$ENDPOINT" \
    s3 sync "$SITE/data" "s3://$BUCKET/data" --only-show-errors --size-only
rc=$?
echo "s3 data rc=$rc in $(( (SECONDS - t0) / 60 )) min"
echo -n "remote payloads: "
aws --profile ewc --endpoint-url "$ENDPOINT" s3 ls "s3://$BUCKET/data/" --recursive --summarize \
    2>/dev/null | awk '/Total Objects:/{print $3}'
echo "ALC_V22_UPLOAD_DONE rc=$rc"
exit $rc
