# Recalibration nuage du CL31 de Payerne (B) sous les spectres alternatifs :
#   constructeur : lambda0 = 910,0 nm (datasheet Vaisala CL31 "910 +- 10 nm at 25 C",
#                  la valeur que Wiegner & Gasteiger 2015 utilisent aussi), largeur mesuree 6,0
#   Wiegner      : 910,0 / FWHM 3,4 (la gaussienne 3-4 nm de Wiegner 2015)
# Le modele actuel (909,7 / 6,0, mesure Qmini) reste la reference operationnelle.
# Sensibilite LUT attendue : +0,5 % et +3,6 % d'absorption effective -> effets faibles, mais
# c'est precisement ce que le dashboard doit montrer.
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$cams = "A:/CAMS_Monthly_04;D:/CAMS_daily"
$base = "C:/DATA/Projects/202606_E-PROFILE_calibration"

$variants = @(
    @{tag = "cl31l910";  spec = '{"CL31": [910.0, 6.0]}'},
    @{tag = "cl31wieg";  spec = '{"CL31": [910.0, 3.4]}'}
)
foreach ($v in $variants) {
    $env:ALC_WV_SPECTRUM = $v.spec
    "=== calibration $($v.tag) (ALC_WV_SPECTRUM=$($v.spec)) ==="
    python scripts/run_streams_parallel.py --streams 0-20000-0-06610_B --start 20250101 --end 20260813 `
        --methods cloud --out "$base/diag_v22_$($v.tag)" --cams $cams --chunks 6 --workers 12
    "=== fin calibration rc=$LASTEXITCODE ==="
    $env:ALC_DUMP_TAG = $v.tag
    $env:ALC_DUMP_UNITS = "0-20000-0-06610_B"
    python rayleigh_availability/cloud_profile_dump.py --workers 10
    "=== fin dump rc=$LASTEXITCODE ==="
}
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_DUMP_TAG,Env:ALC_DUMP_UNITS -ErrorAction SilentlyContinue
"CL31_LAMBDA_DONE"
