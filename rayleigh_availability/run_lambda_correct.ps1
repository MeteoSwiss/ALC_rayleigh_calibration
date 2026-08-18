# Calibration CL61 avec le spectre d'emission CORRECT, tel que les auteurs de
# Le & O'Connor (reponse au relecteur RC1, juin 2026) le tiennent de Vaisala :
#   lambda0 = 910,55 nm, sigma = 0,08 nm  ->  FWHM = 2,3548 x 0,08 = 0,188 nm
# (leur valeur CL51 companion, 910 +- 1,44 nm, redonne exactement le FWHM 3,4 nm de
#  Wiegner & Gasteiger 2015 : la convention est bien sigma, pas FWHM.)
# Notre table operationnelle (910,74 / 1,0) sur-corrige d'un facteur ~9.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$env:ALC_WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'
$base = "C:/DATA/Projects/202606_E-PROFILE_calibration"

"=== calibration CL61 spectre correct (910,55 / FWHM 0,188) ==="
python scripts/run_streams_parallel.py --streams 0-20000-0-06610_C,0-20000-0-10393_C `
    --start 20250101 --end 20260813 --methods rayleigh,cloud `
    --out "$base/diag_v22_l55s008" --cams "A:/CAMS_Monthly_04;D:/CAMS_daily" --chunks 6 --workers 8
"=== fin calibration rc=$LASTEXITCODE ==="

"=== extraction par profil (heatmap Hopkin) pour la meme configuration ==="
$env:ALC_DUMP_TAG = "l55s008"
$env:ALC_DUMP_UNITS = "0-20000-0-06610_C,0-20000-0-06610_B,0-20000-0-10393_C"
python rayleigh_availability/cloud_profile_dump.py --workers 8
"=== fin dump rc=$LASTEXITCODE ==="
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_DUMP_TAG,Env:ALC_DUMP_UNITS -ErrorAction SilentlyContinue
"LAMBDA_CORRECT_DONE"
