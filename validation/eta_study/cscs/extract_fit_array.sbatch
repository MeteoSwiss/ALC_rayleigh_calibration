#!/bin/bash
# Full-network eta-study array on CSCS Linard: one task per unit = in-cloud profile extraction
# (validation/eta_study/extract_incloud_profiles.py --unit) then the step-2 shape fit
# (fit_alpha_shape.py --unit). Keylist: $ETA_KEYLIST (270 units, largest first).
#
# CPU ONLY. Linard partitions without GPUs (sinfo GRES = null, checked 2026-11-02):
#   postproc (default, 1 d), pp-serial (5 d), pp-short (1 h), pp-long (5 d), xfer.
# NEVER use normal, debug, short, short-shared, lowprio, normal-half, normal-shared or
# preemptible (gpu:4 nodes). Never request --gres=gpu or --gpus.
#
# Submit (only after the go):  sbatch validation/eta_study/cscs/extract_fit_array.sbatch
# The array range must match the keylist length: wc -l $ETA_KEYLIST  -> --array=0-(N-1)
#SBATCH --job-name=eta_incloud
#SBATCH --account=msmda
#SBATCH --partition=postproc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# dry run 335708 (CL61, 1 active worker): MaxRSS 1.33 GB -> 8 workers x ~1 GB x 1.5
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --array=0-269%40
#SBATCH --output=/scratch/mch/mhrvo/eta_study/logs/%A_%a.out
set -euo pipefail
source /scratch/mch/mhrvo/eta_study/code/validation/eta_study/cscs/env_eta.sh
# tr: a keylist written on Windows carries CRLF -> the key would end in '\r' (first launch: 80 tasks
# failed with KeyError '...\r')
KEY=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$ETA_KEYLIST" | tr -d '\r')
echo "task $SLURM_ARRAY_TASK_ID -> $KEY on $(hostname) $(date -Is)"
"$PY" -u -W ignore validation/eta_study/extract_incloud_profiles.py --unit "$KEY" --workers "$SLURM_CPUS_PER_TASK"
"$PY" -u -W ignore validation/eta_study/fit_alpha_shape.py --unit "$KEY"
echo "done $(date -Is)"
