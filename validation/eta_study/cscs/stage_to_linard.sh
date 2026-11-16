#!/bin/bash
# Driver (run from the home PC, Git Bash): stage the eta-study code + forward LUT + keylist on
# Linard /scratch. ~115 MB, nothing else is copied (L1 and 0.4-deg CAMS are already on /scratch).
#   bash validation/eta_study/cscs/stage_to_linard.sh            # stage only
#   then on Linard:  sbatch $ETA_CODE/validation/eta_study/cscs/dryrun_pp_short.sbatch
#   full array ONLY after the go:  sbatch $ETA_CODE/validation/eta_study/cscs/extract_fit_array.sbatch
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
LUT=${LUT:-C:/DATA/Projects/202611_eta_study/forward_lut_ext.npz}
REMOTE=/scratch/mch/mhrvo/eta_study
KEYS=$(mktemp)
cd "$REPO"
python -W ignore validation/eta_study/extract_incloud_profiles.py --network --list --keylist "$KEYS" | tail -1
echo "keylist: $(wc -l < "$KEYS") units"
ssh -o BatchMode=yes linard "mkdir -p $REMOTE/code $REMOTE/logs $REMOTE/incloud"
# code: the calibration package (incl. its data LUTs), the eta-study scripts, the PVC module,
# the census and the per-day cloud-calibration table (day lists)
tar czf - --exclude=__pycache__ calibration validation/eta_study validation/multiple_scattering_eta.py \
    validation/scope_l1_2026_census.json rayleigh_availability/cbh_native/cbh_scenes.csv |
    ssh -o BatchMode=yes linard "tar xzf - -C $REMOTE/code"
scp -q -o BatchMode=yes "$LUT" "linard:$REMOTE/forward_lut_ext.npz"
scp -q -o BatchMode=yes "$KEYS" "linard:$REMOTE/keys_network.txt"
rm -f "$KEYS"
ssh -o BatchMode=yes linard "du -sh $REMOTE/code $REMOTE/forward_lut_ext.npz; wc -l $REMOTE/keys_network.txt"
