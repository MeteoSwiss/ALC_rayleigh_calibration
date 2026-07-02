# CSCS runbook — CAMS 0.4° + OmB + sensitivity (2025–2026)

End-to-end steps to produce the Observation-minus-Background (OmB) and instrument-
sensitivity products for the E-PROFILE network on CSCS (balfrin), and rebuild the
dashboard. Everything is **resumable** — re-run any step after an interruption.

> **Why this can't be launched from the dev machine:** the ADS download and the CSCS
> jobs must run on CSCS (internet/ADS auth + the cluster filesystem). The dev machine
> has no CSCS session. Run the steps below on CSCS.

## 0. Environment (once per shell / in `ops/config.sh`)
```bash
conda activate alc
export ALC_L1_ROOT=/capstor/scratch/<user>/E-PROFILE_L1          # L1 archive (<wmo>/<YYYY>/<MM>/)
export ALC_L2_DIR=/capstor/scratch/<user>/E-PROFILE_L2           # L2 (for the OmB operational overlay)
export ALC_CAMS_DIR=/capstor/scratch/<user>/CAMS                 # CAMS_Beta_*.nc (this runbook fills it)
export ALC_FULLCAL_DIR=/capstor/scratch/<user>/fullcal_l1_2026   # calibration + OmB/sens outputs
export ALC_CENSUS=$REPO/validation/scope_l1_2026_census.json
export ADS_API_KEY=<your-ADS-key>            # or have ~/.cdsapirc with the key
export HTTPS_PROXY=http://proxy.cscs.ch:8080 # only if CSCS requires a proxy for outbound HTTPS
```

## 1. CAMS 0.4° download (aerosol + T/RH)  — run on a LOGIN / DATA node
The ADS download is network-bound and needs internet, which CSCS **compute nodes lack**.
Run it on a login/data-mover node inside `tmux`/`screen` (it is resumable, so a dropped
session is fine). Native 0.4° grid (`REGRID_TO_1DEG=False`), Europe+Arctic box
`AREA=[80,-30,27,45]` covering 421/427 stations.
```bash
tmux new -s cams
python scripts/download_cams_cscs.py --start 202501 --end 202612
```
- ~5–8 GB per monthly file (0.4°), ~100–200 GB for 24 months; ADS is queued → can take
  **days** of wall-clock. The per-month skip (file present AND has backscatter) makes
  re-runs safe.
- The 6 non-European affiliates (Canada×3, Bonaire, NZ×2) are outside the box → no OmB
  (sensitivity still works for them; it needs no CAMS aerosol). For those, either accept
  sensitivity-only or run a separate small-box download (`ALC_CAMS_AREA="N,W,S,E"`).

## 2. Calibration (only if not already done for 2025–2026)
The OmB/sens pass REUSES the existing per-stream Kalman, so this is only needed if the
calibration has not been run for the window. SLURM array as usual (~3–5 min/month).
```bash
python scripts/run_all_l1_2026.py --start 20250101 --end 20261231 \
    --per-type 0 --workers <N> --ignore-coverage
```
This writes `<key>_cal.csv`, `<key>_kalman.csv`, `<key>_hk.csv` per stream.

## 3. OmB + sensitivity (reuse the Kalman) — needs step 1 done
```bash
python scripts/run_all_l1_2026.py --no-cal --omb --sens \
    --start 20250101 --end 20261231 --per-type 0 --workers <N> --ignore-coverage
```
- Reads `<key>_kalman.csv` (no recalibration), writes per stream:
  `<key>_omb.png` + `<key>_omb.csv` and `<key>_sens.png` + `<key>_sens.csv`.
- A stream with no Kalman C_L, or no CAMS-with-backscatter for the month, is **skipped**
  (no fabricated numbers). Resumable on the `_sens.csv` marker.
- **Run month-by-month** (e.g. `--start 20250101 --end 20250131`, then 20250201…): OmB
  resolves the CAMS file for the window's month, so a single multi-month window only
  compares the end month. Sensitivity is unaffected by this, but month-aligned runs keep
  both correct and the figures show one month of "evolution over time".
- Memory: OmB loads a full month of native L1 per stream (CL61 ~2 GB → float64 copies in
  `compute_omb` ~9 GB). Size the SLURM task memory accordingly, or lower `--workers`.

## 4. Real-time / daily (D-1)
```bash
Y=$(date -u -d 'yesterday' +%Y%m%d)
python scripts/run_all_l1_2026.py --start $Y --end $Y --ignore-coverage   # cal + Kalman update
python scripts/run_all_l1_2026.py --no-cal --omb --sens --start $Y --end $Y --ignore-coverage
```
(OmB needs that day's CAMS forecast, which publishes ~next day — hence D-1.)

## 5. Dashboard
```bash
python scripts/build_dashboard.py --fullcal $ALC_FULLCAL_DIR \
    --manifest $ALC_CENSUS --l2dir $ALC_L2_DIR --start 20250101 --end 20261231
```
Picks up `<key>_omb.csv`/`<key>_sens.csv` for the two new summary maps (mean OmB bias,
ICAO detection altitude; markers symboled by instrument type) and embeds
`<key>_omb.png`/`<key>_sens.png` on each station page.
