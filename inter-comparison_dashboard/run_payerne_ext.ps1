# Payerne inter-comparison dashboard: calibration reruns that carry every run-based variant of the
# page past the end of its original run, to the end of the extended window (2026-10-31).
# Settings = those of the original runs (eprof_v2.2; the WV-spectrum / no-WV switch of each ladder
# rung; the hood dark v3 for the dark pair).  variants_v3.EXT_TAG maps each output to its variant;
# the builder appends only the nights dated after the original series' last one.  "v2.0" is not here: it is
# the published v2.0 NetCDFs, which cannot be rerun.
# Needs the October CAMS daily files in D:/CAMS_daily (synced from balfrin:/scratch/mch/mhrvo/CAMS).
$ErrorActionPreference = "Continue"
Set-Location C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration
$D = "C:/DATA/Projects/202606_E-PROFILE_calibration"
$OUT = "$D/dashboard_ext_20261031"
$cams = "D:/CAMS_Monthly_04;D:/CAMS_daily"
$env:ALC_MOLECULAR_METHOD = "eprof_v2.2"
$env:ALC_L1_ROOT = "D:/E-PROFILE_L1_2026"
$A = "0-20000-0-06610_A"; $B = "0-20000-0-06610_B"; $C = "0-20000-0-06610_C"

function Run($tag, $streams, $methods, $start, $spec, $nowv, $dark) {
    if ($spec) { $env:ALC_WV_SPECTRUM = $spec } else { Remove-Item Env:ALC_WV_SPECTRUM -ErrorAction SilentlyContinue }
    if ($nowv) { $env:ALC_WV_DISABLE = "1" } else { Remove-Item Env:ALC_WV_DISABLE -ErrorAction SilentlyContinue }
    if ($dark) { $env:ALC_DARK_PROFILE = "$D/rayleigh_availability/dark_profiles_payerne_v3.npz" }
    else { Remove-Item Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue }
    "=== $tag : $streams [$methods] from $start (spec='$spec' nowv=$nowv dark=$dark) ==="
    python scripts/run_streams_parallel.py --streams $streams --start $start --end 20261031 `
        --methods $methods --out "$OUT/$tag" --cams $cams --chunks 8 --workers 8
    "=== $tag rc=$LASTEXITCODE ==="
}

$op = '{"CL61": [910.55, 0.188]}'                  # operational spectrum (CL31 keeps 909.7/6.0)
# from the day after the original runs end (diag_v22_*, the v2.2 json corpus: 2026-08-13)
Run "l55s008"  $A "rayleigh"       "20260814" $op $false $false   # CHM15k v2.2 (no WV term at 1064 nm)
Run "l55s008"  $C "rayleigh,cloud" "20260814" $op $false $false   # CL61 * (both methods)
Run "l55s008"  $B "cloud"          "20260814" $op $false $false   # CL31 cloud, its operational 909.7/6.0
Run "nominal"  $C "rayleigh,cloud" "20260814" '{"CL61": [910.74, 1.0]}' $false $false   # retired CL61 model
Run "l55w10"   $C "rayleigh,cloud" "20260814" '{"CL61": [910.55, 1.0]}' $false $false
Run "l55w01"   $C "rayleigh,cloud" "20260814" '{"CL61": [910.55, 0.1]}' $false $false
Run "nowv"     $C "rayleigh,cloud" "20260814" "" $true $false     # no WV correction
Run "nowv"     $B "cloud"          "20260814" "" $true $false
Run "cl31l910" $B "cloud"          "20260814" '{"CL31": [910.0, 6.0]}' $false $false   # datasheet nominal
Run "cl31wieg" $B "cloud"          "20260814" '{"CL31": [910.0, 3.4]}' $false $false   # Wiegner 2015
# hood dark v3: from the day after diag_v22_dark_v3 ends (2026-10-08)
Run "dark_v3"  "$A,$C" "rayleigh"  "20261009" $op $false $true

Remove-Item Env:ALC_WV_SPECTRUM, Env:ALC_WV_DISABLE, Env:ALC_DARK_PROFILE -ErrorAction SilentlyContinue
"PAYERNE_EXT_DONE"
