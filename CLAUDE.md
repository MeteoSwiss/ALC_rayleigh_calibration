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
3. **calibrate** D-1 across the network: `scripts/run_all_l1_2026.py --sens --omb` (Rayleigh +
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

## Rayleigh overlays (dashboard)

Station time-series charts show two extra Rayleigh overlays, both **hidden by default**
(`visible="legendonly"`, appear only when clicked in the legend), reloaded fresh each build via
`monitoring/render.py::_load_oldray`:

- **v1.0** (operational) from `ALC_OLDRAY_DIR=/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh`
  (year subdirs `2025/`,`2026/`; the RAW `ALC_calibration_<key><YYYY>.nc`, NOT the `kalman/` subdir).
- **v1.0.2** (was "v13", a rewrite of v1.0) from `ALC_V13_DIR`. Produced by hem's cron
  `18 18 * * *` running `rayleigh_calibration.main` (repo
  `/proj/pay/E-PROFILE/Calibration_codes/dev/rayleigh_calibration`). Needs `netCDF4` in pyenv
  3.12.12.

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

## Working preferences (user hervo63)

- **Ask before relocating/moving data**; when proposing to write somewhere, state *what* / *how much*
  / *where* up front. In-place edits and computation are fine without asking.
- `/home/pay/...` is **horribly slow** — prefer `/tmp` for transient/scratch files.
- Wants honest trade-off analysis ending in a **clear recommendation**, not a menu of options.
- Comfortable in French (often writes FR, mixes EN).
- **The user pushes git themselves** — prepare commits, do not `git push`.
