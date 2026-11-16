#!/bin/bash
# Dry run of the eta-study array on ONE unit-month (Payerne CL61, 2026-04), to compare with the
# home-PC result (1 819 profiles, identical C to the operational dump). CPU-only partition
# pp-short (1 h; GRES null on Linard). Writes to $ETA_BASE/dryrun, not to the real output dirs.
#SBATCH --job-name=eta_dryrun
#SBATCH --account=msmda
#SBATCH --partition=pp-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/mch/mhrvo/eta_study/logs/dryrun_%j.out
set -euo pipefail
source /scratch/mch/mhrvo/eta_study/code/validation/eta_study/cscs/env_eta.sh
export ETA_OUT=$ETA_BASE/dryrun/incloud
export ETA_BASE_FIT=$ETA_BASE/dryrun
KEY=0-20000-0-06610_C
echo "dry run $KEY 202604 on $(hostname) $(date -Is)"
SECONDS=0
"$PY" -u -W ignore validation/eta_study/extract_incloud_profiles.py --unit "$KEY" \
    --months 202604 --workers "$SLURM_CPUS_PER_TASK"
echo "extraction wall ${SECONDS} s"
SECONDS=0
ETA_BASE=$ETA_BASE_FIT "$PY" -u -W ignore validation/eta_study/fit_alpha_shape.py --unit "$KEY"
echo "fit wall ${SECONDS} s"
sstat -j "$SLURM_JOB_ID" --format=MaxRSS 2>/dev/null || true
echo "done $(date -Is)"
