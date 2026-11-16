# Paper Section 6.1: cloud vs Rayleigh on every instrument type, BOTH methods forced on all types
# (operational routing excludes Rayleigh for CL31/CL51 and cloud for CHM15k/Mini-MPL; this run shows why).
# Feb-Jun 2026, operational v2.2 settings, CL61 with the operational (910.55 nm / 0.188 nm) spectrum, no dark.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$env:ALC_WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'
$env:ALC_RAYLEIGH_TYPES_EXTRA = "CL31,CL51"
$env:ALC_CLOUD_TYPES_EXTRA = "CHM15k,Mini-MPL"
python scripts/run_streams_parallel.py --streams 0-20000-0-06610_B,0-20000-0-03808_B,0-250-1001-07151_A,0-20000-0-08224_A,0-20000-0-11538_A,0-20000-0-02453_A,0-20000-0-16714_A,0-20000-0-06447_A,0-20000-0-04018_A,0-20000-0-11035_A,0-20000-0-11518_A,0-196-0-CYPN_A,0-20000-0-06610_A,0-20000-0-10393_0,0-20000-0-03808_A,0-20000-0-06240_A,0-250-1001-07151_B,0-20008-0-INO_A,0-20008-0-UGR_A,0-20000-0-07014_A,0-20000-0-07110_A,0-20000-0-07617_A,0-20000-0-07774_A,0-20000-0-03808_C,0-20000-0-06418_B,0-20000-0-06447_B,0-20000-0-06610_C,0-20000-0-10393_C,0-20000-0-11538_B,0-20008-0-BIR_A,0-20008-0-LAU_A,0-203-10-LNG_A,0-380-5-1_B,0-756-4-EERLCL61_A --start 20260201 --end 20260629 `
    --methods rayleigh,cloud --out "C:/DATA/Projects/202606_E-PROFILE_calibration/method_test_2026" `
    --cams "D:/CAMS_Monthly_04;D:/CAMS_daily" --chunks 5 --workers 30
"rc=$LASTEXITCODE"
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_RAYLEIGH_TYPES_EXTRA,Env:ALC_CLOUD_TYPES_EXTRA -ErrorAction SilentlyContinue
"METHOD_TEST_DONE"
