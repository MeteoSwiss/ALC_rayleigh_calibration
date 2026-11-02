#!/bin/bash
# Per-unit summary of the Linard eta array (validation/eta_study/collect_unit_summary.py): small CSVs
# in $ETA_BASE/summary/ -- the only files copied back to the PC. CPU-only partition (pp-short).
# Chain after the arrays:  sbatch --dependency=afterany:<array1>:<array2> collect_summary.sbatch
#SBATCH --job-name=eta_collect
#SBATCH --account=msmda
#SBATCH --partition=pp-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --chdir=/scratch/mch/mhrvo/eta_study
#SBATCH --output=/scratch/mch/mhrvo/eta_study/logs/collect_%j.out
set -euo pipefail
source /scratch/mch/mhrvo/eta_study/code/validation/eta_study/cscs/env_eta.sh
"$PY" -u -W ignore validation/eta_study/collect_unit_summary.py
ls -la "$ETA_BASE/summary"
