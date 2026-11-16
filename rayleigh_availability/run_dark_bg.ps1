# Appendix 6 variants: (1) measured hood dark + operational far-field background bias, (2) background bias
# only (hood-free). Background bias = median over clear nights of [measured - modelled molecular] P at
# 12-15 km (rayleigh_availability/background_bias_operational.json), added back as c*r^2.
# Identical to the operational v2.2 otherwise; CL61 with the operational spectrum.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$env:ALC_WV_SPECTRUM = '{"CL61": [910.55, 0.188]}'
$d = "C:/DATA/Projects/202606_E-PROFILE_calibration"
foreach ($j in @(@{npz="$d/rayleigh_availability/dark_profiles_payerne_v2_bg.npz"; out="$d/diag_v22_darkbg_v2"},
                 @{npz="$d/rayleigh_availability/background_only_payerne.npz"; out="$d/diag_v22_bgonly"})) {
    $env:ALC_DARK_PROFILE = $j.npz
    "=== $($j.out) ==="
    python scripts/run_streams_parallel.py --streams 0-20000-0-06610_A,0-20000-0-06610_C --start 20250101 --end 20261008 `
        --methods rayleigh --out $j.out --cams "D:/CAMS_Monthly_04;D:/CAMS_daily" --chunks 10 --workers 20
    "rc=$LASTEXITCODE"
}
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue
"DARK_BG_DONE"
