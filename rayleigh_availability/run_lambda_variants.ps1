# Calibrations CL61 sous les differentes hypotheses de spectre laser, pour comparaison dashboard.
# Echelle physique complete, de "toute l'absorption" a "aucune" :
#   910,74 / 1,0 nm  = mesure Qmini (OPERATIONNEL, deja calcule dans calout_v22_04)
#   910,55 / 1,0 nm  = spec constructeur, meme largeur   -> isole l'effet de POSITION
#   910,55 / 0,1 nm  = spec + "very narrow bandwidth"    -> l'hypothese de mitigation par conception
#   sans WV          = limite f = 0                      (deja calcule dans diag_v22_nowv)
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$cams = "A:/CAMS_Monthly_04;D:/CAMS_daily"
$base = "C:/DATA/Projects/202606_E-PROFILE_calibration"
$streams = "0-20000-0-06610_C,0-20000-0-10393_C"     # les deux CL61 des sites du dashboard

$variants = @(
    @{tag = "l55w10"; spec = '{"CL61": [910.55, 1.0]}'},
    @{tag = "l55w01"; spec = '{"CL61": [910.55, 0.1]}'}
)
foreach ($v in $variants) {
    $env:ALC_WV_SPECTRUM = $v.spec
    $out = "$base/diag_v22_$($v.tag)"
    "=== $out  (ALC_WV_SPECTRUM=$($v.spec)) ==="
    python scripts/run_streams_parallel.py --streams $streams --start 20250101 --end 20260813 `
        --methods rayleigh,cloud --out $out --cams $cams --chunks 6 --workers 22
    "=== fin $out rc=$LASTEXITCODE ==="
}
Remove-Item Env:ALC_WV_SPECTRUM -ErrorAction SilentlyContinue
"LAMBDA_VARIANTS_DONE"
