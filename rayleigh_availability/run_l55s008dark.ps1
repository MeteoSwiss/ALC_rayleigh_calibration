# Calibration Rayleigh CL61 Payerne COMPLETEMENT corrigee : spectre constructeur
# (910,55 / sigma 0,08 -> FWHM 0,188) + dark mesure sous capot soustrait.
# C'est l'etat de reference du canal Rayleigh CL61 — les deux corrections etaient jusqu'ici
# produites separement (diag_v22_l55s008 sans dark ; diag_v22_dark avec l'ancien modele WV).
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$env:ALC_WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'
$env:ALC_DARK_PROFILE = "C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_profiles_payerne.npz"
python scripts/run_streams_parallel.py --streams 0-20000-0-06610_C --start 20250101 --end 20260813 `
    --methods rayleigh --out "C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_l55s008dark" `
    --cams "A:/CAMS_Monthly_04;D:/CAMS_daily" --chunks 6 --workers 12
"rc=$LASTEXITCODE"
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue
"L55S008DARK_DONE"
