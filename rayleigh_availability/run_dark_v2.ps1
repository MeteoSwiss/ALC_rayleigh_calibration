# Dark-corrected Rayleigh reruns for the paper appendix, extended to 2026-10-08 with the v2 dark profile
# (10 sessions, covered intervals verified in the L1 data: dark_windows_detect.py + build_dark_v2.py).
# Identical to the operational v2.2 except the dark subtraction; CL61 with the operational spectrum.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$env:ALC_WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'
$env:ALC_DARK_PROFILE = "C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_profiles_payerne_v2.npz"
python scripts/run_streams_parallel.py --streams 0-20000-0-06610_A,0-20000-0-06610_C --start 20250101 --end 20261008 `
    --methods rayleigh --out "C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_dark_v2" `
    --cams "D:/CAMS_Monthly_04;D:/CAMS_daily" --chunks 10 --workers 20
"rc=$LASTEXITCODE"
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue
"DARK_V2_DONE"
