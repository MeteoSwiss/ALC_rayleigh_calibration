# Runs apparies avec/sans dark estime -- Amsterdam (4 CHM15k) et Lindenberg (CHM15k+CL61).
# Strictement identiques sauf ALC_DARK_PROFILE. Methode eprof_v2.2, CAMS 0.4 deg, Rayleigh seul.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$cams = "A:/CAMS_Monthly_04;D:/CAMS_daily"
$base = "C:/DATA/Projects/202606_E-PROFILE_calibration/dark_test"
$dk = "C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_clearsky"
New-Item -ItemType Directory -Force $base | Out-Null

$jobs = @(
    @{streams="0-20000-0-06240_A,0-20000-0-06240_B,0-20000-0-06240_C,0-20000-0-06240_D"; out="$base/amst_nodark"; dark=""},
    @{streams="0-20000-0-06240_A,0-20000-0-06240_B,0-20000-0-06240_C,0-20000-0-06240_D"; out="$base/amst_dark";   dark="$dk/dark_est_0-20000-0-06240.npz"},
    @{streams="0-20000-0-10393_0,0-20000-0-10393_C"; out="$base/lind_nodark"; dark=""},
    @{streams="0-20000-0-10393_0,0-20000-0-10393_C"; out="$base/lind_dark";   dark="$dk/dark_est_0-20000-0-10393.npz"}
)
foreach ($j in $jobs) {
    if ($j.dark) { $env:ALC_DARK_PROFILE = $j.dark } else { Remove-Item Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue }
    "=== $($j.out) (dark='$($j.dark)') ==="
    python scripts/run_streams_parallel.py --streams $j.streams --start 20250101 --end 20260813 `
        --methods rayleigh --out $j.out --cams $cams --chunks 6 --workers 24
    "=== fin $($j.out) rc=$LASTEXITCODE ==="
}
"TOUT_TERMINE"
