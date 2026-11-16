# No-water-vapour calibration run diag_v22_nowv (paper Section 5.1.1, validation/paper/wv_magnitude_cl61.py).
# Identical to the operational v2.2 except ALC_WV_DISABLE=1. Payerne CL31 + CL61 and Lindenberg CL61,
# 2025-01-01 .. 2026-08-13. Launched on 2026-08-16 as two inline commands (system python 3.12); recovered here.
# CAMS: the 0.4-deg monthlies were on A: then; they now live on D:/CAMS_Monthly_04.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_WV_DISABLE = "1"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$OUT = "C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_nowv"
$CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"     # was 'A:/CAMS_Monthly_04;D:/CAMS_daily'
python -u scripts/run_streams_parallel.py --streams 0-20000-0-06610_B,0-20000-0-06610_C --start 20250101 --end 20260813 --methods rayleigh,cloud --out $OUT --cams $CAMS --chunks 6 --workers 20
python -u scripts/run_streams_parallel.py --streams 0-20000-0-10393_C --start 20250101 --end 20260813 --methods rayleigh,cloud --out $OUT --cams $CAMS --chunks 6 --workers 18
