# Extractions par profil (paires C_i, CBH_i) pour la heatmap Hopkin du dashboard, une par
# CONFIGURATION NUAGE. Seules les unites des sites du dashboard sont necessaires ; le script
# saute les npz deja presents, donc on restreint via ALC_DUMP_UNITS.
#   nominal  : 910,74 / 1,0  (operationnel)
#   l55w10   : 910,55 / 1,0  (spec constructeur, meme largeur)
#   l55w01   : 910,55 / 0,1  (spec + "very narrow")
#   nowv     : sans correction WV
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$env:ALC_DUMP_UNITS = "0-20000-0-06610_C,0-20000-0-06610_B,0-20000-0-10393_C"

$cfgs = @(
    @{tag = "nominal"; spec = ""; nowv = "0"},
    @{tag = "l55w10";  spec = '{"CL61": [910.55, 1.0]}'; nowv = "0"},
    @{tag = "l55w01";  spec = '{"CL61": [910.55, 0.1]}'; nowv = "0"},
    @{tag = "nowv";    spec = ""; nowv = "1"}
)
foreach ($c in $cfgs) {
    $env:ALC_DUMP_TAG = $c.tag
    if ($c.spec) { $env:ALC_WV_SPECTRUM = $c.spec } else { Remove-Item Env:ALC_WV_SPECTRUM -ErrorAction SilentlyContinue }
    if ($c.nowv -eq "1") { $env:ALC_WV_DISABLE = "1" } else { Remove-Item Env:ALC_WV_DISABLE -ErrorAction SilentlyContinue }
    "=== dump configuration $($c.tag) (spec='$($c.spec)' nowv=$($c.nowv)) ==="
    python rayleigh_availability/cloud_profile_dump.py --workers 10
    "=== fin $($c.tag) rc=$LASTEXITCODE ==="
}
Remove-Item Env:ALC_WV_SPECTRUM,Env:ALC_WV_DISABLE,Env:ALC_DUMP_TAG,Env:ALC_DUMP_UNITS -ErrorAction SilentlyContinue
"HOPKIN_VARIANTS_DONE"
