# Runs the resumable 170-station L2 calibration until the ALL_DONE flag appears.
# Intended to be launched by the Windows Task Scheduler so it is independent of
# the Claude session (which kills session-spawned processes after ~30 min).
$ErrorActionPreference = "Continue"
$repo = "C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration"
$out  = "D:\E-PROFILE_calibration_rayleigh\fullcal_all"
$flag = Join-Path $out "ALL_DONE.flag"
Set-Location $repo
$i = 0
while (-not (Test-Path $flag)) {
    $i++
    "$(Get-Date -Format o)  loop iteration $i starting" | Out-File -Append "$out\sched.log"
    python run_all_l2monthly.py --workers 6 *>> "$out\sched.log"
    Start-Sleep -Seconds 5
}
"$(Get-Date -Format o)  ALL DONE flag present - exiting loop" | Out-File -Append "$out\sched.log"
