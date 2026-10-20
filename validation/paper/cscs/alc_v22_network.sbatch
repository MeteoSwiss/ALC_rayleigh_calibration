#!/bin/bash
#SBATCH --job-name=alc_v22_net
#SBATCH --account=msmda
#SBATCH --partition=postproc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --chdir=/scratch/mch/mhrvo
#SBATCH --output=/scratch/mch/mhrvo/logs/alc_v22_net_%j.log

# Network v2.2 rerun ("calout_v22_04") used by the paper: Rayleigh (eprof_v2.2, chi2red <= 2.5) and
# liquid-cloud (PVC eta tables) calibration of every stream of the 2026 L1 census, 2025-01-01..2026-08-13,
# on the 0.4-deg CAMS monthlies (1-deg monthlies as per-month fallback).
# Paper use: Figure 4 (CHM15k ranking, fig_v09_chm15k_network.py) and Figure 6 (CL31/CL51 cloud ranking,
# fig_v09_cloud_network.py). Verified identical to the dashboard nightly constants (median ratio 1.00000).
#
# Submitted on balfrin (CSCS) on 2026-08-15 from ~/alc_v22_network.sbatch, one full CPU postproc node.
# Code: a snapshot of this repository in $HOME/alc_v22_cscs (PYTHONPATH), synced before submission.
# Output: /scratch/mch/mhrvo/E_PROFILE_calout_v22_04, copied to
#         C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04 (<key>/<key>_cal.csv per stream).
# Recovered from the session log of 2026-08-15 (the file existed only in the balfrin home directory).
set -u
mkdir -p /scratch/mch/mhrvo/logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate alc
export PYTHONPATH=$HOME/alc_v22_cscs

export ALC_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1
export ALC_MOLECULAR_METHOD=eprof_v2.2
export ALC_MOLECULAR_PARAMS='{"max_chi2red": 2.5}'
export PLOTS=0

STREAMS=$(python - <<'PY'
import json
c = json.load(open('/users/mhrvo/alc_v22_cscs/validation/scope_l1_2026_census.json'))
rows = c if isinstance(c, list) else c.get('streams', [])
print(','.join(sorted({r['wmo'] + '_' + r['ident'] for r in rows})))
PY
)
echo "streams: $(echo $STREAMS | tr ',' '\n' | wc -l)"

python $HOME/alc_v22_cscs/scripts/run_streams_parallel.py \
  --streams "$STREAMS" \
  --start 20250101 --end 20260813 --chunks 2 --methods rayleigh,cloud \
  --cams '/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS' \
  --workers 230 \
  --out /scratch/mch/mhrvo/E_PROFILE_calout_v22_04
