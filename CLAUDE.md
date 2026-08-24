# CLAUDE.md — E-PROFILE ALC calibration (durable operator context)

This file is **committed to the repo** so Claude Code auto-loads it on *every* machine that
checks this repo out (it is the cross-machine equivalent of the per-machine assistant memory).
Keep it in sync with `doc/OPERATIONS.md` (the long-form ops guide) — this is the condensed
"what you must not forget" sheet.

## System & where it lives

The E-PROFILE ALC Rayleigh/cloud calibration pipeline **and** its dashboard run on
**`zueub434.meteoswiss.ch`** under `/data/zue/E_PROFILE/ALC/Calibration/`:

- code: `ALC_calibration_v2.0_code/` (this repo)
- calibration output (fullcal): `ALC_calibration_v2.0/` — per-station `<key>/` dirs (~436)
- dashboard build dir: `dashboard/`
- L1 input: `/data/zue/E_PROFILE/ALC/L1_FILES/<wmo>/<year>/<month>/L1_<wmo>_<ident><YYYYMMDD>.nc`
- CAMS cache: `/data/zue/E_PROFILE/ALC/CAMS/CAMS_Beta_<YYYYMMDD>.nc`

Relocated here in 2026-06 off `/mnt/amaroc_data/alc_calib`, a shared NAS stuck at its **inode
ceiling** (`df -i` ~100% while TBs of bytes are free → errors show as "No space left on device").
`/mnt/amaroc_data/...` is the OLD location and must not be written to. The real fix is an admin
`maxfiles` raise on that volume — the camera thinning cron can't outpace the cameras.

Always orient first: `source ops/config.sh` (exports all `ALC_*` paths + S3 creds + `ALC_VENV`),
then `source "$ALC_VENV/bin/activate"`.

Live dashboard: <https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/>

## Daily flow

`cron 0 15 * * *` → `ops/run_daily.sh` → `ops/ops_daily.py`:

1. **refresh census** (`scripts/refresh_census.py`) — scan the L1 archive and merge new stations
   into `validation/scope_l1_2026_census.json` (new streams appended, existing never dropped) so a
   newly-installed station is calibrated the same day.
2. **fetch CAMS** for D-1 (ADS download, retried).
3. **calibrate** D-1 across the network: `scripts/run_network_calibration.py --sens --omb` (Rayleigh +
   liquid-cloud + Kalman; per-day caches in `calibration/incremental.py` + a regression guard so a
   missing cache never overwrites a rich history; MERGES into per-stream CSVs, no overwrite).
4. **update_opcoeff** (`extract_l2_opcoeff.py` → `operational_coefficients.csv`).
5. **build dashboard** (`build_dashboard.py --changed-only`, bucket-mode: images served from the EWC
   S3 bucket, HTML references bucket URLs).
6. **publish** (`ops/publish.sh`: images→S3 bucket `eprofile-alc-dashboard`, HTML→web VM
   `hem@136.156.139.31:/var/www/alc`, then prunes local diag PNGs once on the bucket).

Target days = D-`ALC_DAY_LAG`(=1) + the last `ALC_BACKFILL_DAYS`(=5) unprocessed days (self-healing).
Everything is driven by `ALC_*` env vars in `ops/config.sh` — the single file you edit. Dates UTC.

## Operational gotchas

- **Become `rem`**: only `sudo su - rem` is NOPASSWD. `sudo -u rem`, `sudo -n`, and
  `sudo su - rem -c '...'` all demand a password → pass commands via **STDIN**
  (`echo 'cmd' | sudo su - rem`), never `-c`.
- `rem`'s login shell has **`noclobber`** → `cat > existing_file` fails; `rm -f file` first.
- `ops/publish.sh` HTML→VM rsync intermittently returns **rc=2** yet the HTML usually still lands;
  if the live site lags, just re-run `bash ops/publish.sh`.
- Internet egress (pip, CAMS/ADS API) needs `export https_proxy=http://proxy.meteoswiss.ch:8080`
  (and `http_proxy`); `config.sh` already sets these.
- E-PROFILE **L1 for a day lands the next morning** (~03:30Z) → a daily run must fire after that
  (15:00 cron is fine; a pre-dawn run finds no data for "yesterday").
- **balfrin (CSCS)**: from zueub434 `ssh -n -o BatchMode=yes balfrin`. Compute nodes have no
  internet; SLURM can park failed array tasks as held ("launch failed requeued held") — release
  with `scontrol release <jobid>`.
  - **⚠️ CPU nodes ONLY — never GPU. GPU runs cost a fortune (a recent batch was > $10,000).**
    Always submit to a **CPU partition**: `pp-short` (1 h; CI/pre-post), `pp-serial` (120 h, 1 core;
    verification), `postproc` (24 h; analysis), `pp-long` (120 h; long analysis),
    `pp-production`/`pp-prodntc`/`pp-dispntc` (production, restricted users).
    **Do NOT use** the GPU partitions `debug`, `short`, `short-shared`, `normal`, `normal-shared`,
    `lowprio`, `preemptible`, `production`. If a job seems to need a GPU, **stop and ask first**.

## Rayleigh overlay (dashboard)

Station time-series charts show one extra Rayleigh overlay, **hidden by default**
(`visible="legendonly"`, appears only when clicked in the legend), reloaded fresh each build via
`monitoring/render.py::_load_oldray`:

- **v1.0** (operational) from `ALC_OLDRAY_DIR=/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh`
  (year subdirs `2025/`,`2026/`; the RAW `ALC_calibration_<key><YYYY>.nc`, NOT the `kalman/` subdir).

The v1.0.2 / "v13" test overlay (`ALC_V13_DIR`) was retired 2026-07 and removed from code and
docs; hem's `18 18 * * *` cron that produced it should be deleted from the server crontab.

## CAMS domains (regional boxes)

The daily CAMS download is regional. The default **Europe+Arctic** box (`ALC_CAMS_AREA`, N80/W-30/S27/E45,
0.4°) covers ~427 census stations and keeps the legacy `CAMS_Beta_<date>.nc` name. A handful of
affiliates fall outside it and are served by their own **small** boxes (file
`CAMS_Beta_<region>_<date>.nc`), routed by station lat/lon in `calibration.io.download_cams_beta`:
`namerica_west` (Edmonton), `ontario` (Western/London ON), `caribbean` (Bonaire), `newzealand`
(Lauder + Auckland). A cron at `02 06 * * *` (`ops/prefetch_cams.sh`) downloads all the day's boxes
ahead of the 15:00 calibration. The water-vapour correction is **mandatory** for 910 nm instruments —
a degraded no-WV mode is rejected (it worsens results); a 910 nm night without usable CAMS is flagged,
never calibrated WV-free.

## Multiple-scattering & water-vapour source (cloud calibration)

- **MS correction (2026-07)**: the O'Connor/Hopkin η(cloud-base) tables are now the **PVC (Hogan 2006)
  tables at droplet radius a_G = 5.5 µm** (the Cloudnet-measured calibration-scene size, 11 µm
  diameter; alpha = 10 /km), replacing the legacy Hopkin / fitted 8 µm ladder. Each Vaisala type has
  its **own** table — CL31's wider 0.83 mrad FOV gets a stronger correction; CL61 has its own table
  (no longer borrows CL51). At low cloud base η ≈ 0.95 (was 0.83) → ~9 % lower C for low clouds (the
  bulk of calibration scenes). Derivation/validation: `validation/multiple_scattering_eta.py`,
  `doc/reports/multiple_scattering_check.md`.
- **WV humidity source**: operational default is **CAMS model levels (L137)** (`wv_source='cams'`,
  dense in the boundary layer). An ERA5 path exists (`wv_source='era5'` + a prefetched Earth Data Hub
  cache from `scripts/prefetch_era5_edh.py`) but is **research-only** — the Hub's ERA5 is a 19-level
  pressure subset (~4 levels below 3 km → coarser BL than CAMS). ERA5 model levels aren't efficiently
  available (EDH lacks them; ARCO-ERA5 model levels are chunked whole-field → multi-TB for per-station
  time-series; ECMWF Polytope `class=ea` blocked). So CAMS model levels stay operational.
- **OmB + WV CAMS resolution**: monthly **0.4°** (`ALC_CAMS_DIR`) with a per-month **1° fallback**
  (`ALC_CAMS_DIR_FALLBACK`) for months the 0.4° download has not covered — both are L137 model levels,
  so vertical resolution is preserved. Applies to the cloud WV correction and OmB alike.

## Candidate-method & CAMS plumbing (2026-08)

- `cams_folder` (options.json / `ALC_CAMS_DIR`) accepts a **`;`-separated folder list**, searched in
  order (monthly then daily per folder) — e.g. `"A:/CAMS_Monthly_04;D:/CAMS_daily"`. ⚠️ Never point
  it at a folder holding **1° monthlies** (`D:/CAMS_run_v20`, `D:/CAMS`): the 1° grid-point orography
  at Payerne is 894 m too high → PWV −26 % → WV column truncated.
- Runner env overrides (leave unset → options.json wins, operational cron unaffected):
  `ALC_MOLECULAR_METHOD` / `ALC_MOLECULAR_PARAMS` (JSON) — run a candidate algorithm network-wide;
  `ALC_DARK_PROFILE` — measured dark-baseline npz (see below); `ALC_L1_ROOT`, `ALC_CAMS_DIR`.
- `scripts/run_streams_parallel.py` — few-streams × long-window runs on many cores (chunks the date
  range, per-chunk output dirs — **mandatory**, concurrent writers corrupt the shared per-ident
  NetCDF — then merges; verified bit-identical to a sequential run).
- **Measured dark baseline** (`dark_profile_file` option / `ALC_DARK_PROFILE`): subtracts the
  covered-telescope b(z) profile (rcs_0 units, `rayleigh_availability/dark_profiles.py`) before the
  Rayleigh fit. At Payerne: +24 % on the CHM15k constant, +24 nights/yr, gradient −10.5 → −2.5 %/km.
  Do **NOT** confuse with `subtract_background` (fitted intercept = atmosphere in disguise, keep 0).
- Campaign reports index: `doc/reports/README.md` §2026-08.

## CL61 water-vapour spectrum & the v3 dashboard (2026-08-16)

- **The CL61 emission line is (910.55 nm, σ 0.08 → FWHM 0.188 nm)** — Le & O'Connor response to
  referee RC1 (Vaisala pers. comm.), corroborated by the PUBLIC DA10 DIAL guide M212895EN-E
  ("910.55 and 910.99 nm", offline channel "with low water vapor absorption") and Mariani 2021
  (DA10 transmitter 0.19–0.21 nm). The old model (910.74/1.0, Qmini centroid bias + wrong width)
  **over-corrects WV ×9** → it CREATED the CL61 cloud dC/dCBH slope (+9 %/km) and ~−12 % of level.
  Fully-corrected CL61 (★ spectrum + measured dark) lands within ~2 % of the CHM15k reference.
- Env switches (default unset = operational unchanged): **`ALC_WV_SPECTRUM`** JSON overrides the
  laser spectrum for BOTH methods (e.g. `{"CL61": [910.55, 0.188]}`; single source of truth
  `water_vapor.LASER_SPECTRUM`, cloud `set_defaults` reads it too); **`ALC_WV_DISABLE=1`** turns
  the WV correction off. Runs per hypothesis: `diag_v22_l55s008` (★), `_l55w10`, `_l55w01`,
  `_nowv`, `_l55s008dark` (★+dark, the CL61 Rayleigh reference), CL31 ladder `_cl31l910`/`_cl31wieg`.
- Three CL61 wavelengths coexist ON PURPOSE: 910.74→WV correction table (update to 910.55/0.188
  pending final arbitration), 910.0→molecular Rayleigh reference (0.33 % effect), 910.55→η tables.
- **Dashboard v3** = `inter-comparison_dashboard/index.html` (one page, 3 sites, per-instrument
  method/variant/dark controls, all-variants ladder table, per-config Hopkin, binned-PWV panel);
  architecture + editing guide in `inter-comparison_dashboard/README.md`; **`variants_v3.py` is
  the file to edit when a new run lands**. Anti-fossil rule: state-dependent numbers in page
  prose are injected at build (`prose_tokens`), never hard-coded.
- **Noise filter (2026-08-24)**: « Filtre bruit » = SNR≥3 admission masks precomputed by
  `nf_v3.py` from the NATIVE L1 files (`l1_l2_io.read_l1_native` — the `_streams_*` cache is
  HOURLY, never use it for noise statistics), windows 5 min (default)/30/60/3 h, four modes
  (par instrument / intersection / moyenne-d'abord / masque scène, referee = the site CHM15k,
  `nf_ref` in variants_v3). SNR filtering conditions the sample on signal: at Payerne 5 min the
  CL31 keeps 25 % of band hours and the measured sampling bias (CHM15k median under the CL31
  mask) is **+26 %** — shown live on the page; the scene/intersection/aggregate modes are the
  unbiased comparisons. Masks frozen at site-default constants; L2 inherits L1 masks;
  curtains/PWV stay static. Verification: `check_v3.py` (gate 0.1 %).
- Cloud calibration is dark-immune (+0.08 % CL31, ~0.001 % CL61, measured) and its fixed
  100–2400 m window makes a static dark CBH-flat by construction — never rerun cloud for dark.
- L2 as distributed: CHM15k **and CL61** = operational Rayleigh v1.0 (the CL61's v1.0
  calibration went live recently: constant 1.0 default on 48 d then 2.1257 from mid-2026-06 on
  the paired window — the transition is visible in the L2 row tooltip); CL31 = UNCALIBRATED
  (flat default 1e8 — the gap the cloud calibration fills).

## Working preferences (user hervo63)

- **Ask before relocating/moving data**; when proposing to write somewhere, state *what* / *how much*
  / *where* up front. In-place edits and computation are fine without asking.
- `/home/pay/...` is **horribly slow** — prefer `/tmp` for transient/scratch files.
- Wants honest trade-off analysis ending in a **clear recommendation**, not a menu of options.
- Comfortable in French (often writes FR, mixes EN).
- **The user pushes git themselves** — prepare commits, do not `git push`.
