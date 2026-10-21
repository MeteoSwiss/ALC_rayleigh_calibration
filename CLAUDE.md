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
- **The 910.55/0.188 spectrum is now OPERATIONAL (verified 2026-08-31).** The arbitration is
  settled: the per-night NetCDFs published by E-PROFILE (EWC bucket `nc/<wmo>_<ident>/`,
  `calibration_version` 220 Rayleigh / 1000 cloud) reproduce our `diag_v22_l55s008` reruns
  **night-for-night** — cloud 70/70 nights within 2 % at a median ratio of 1.0000, Rayleigh 37/40 —
  and match NO nominal-spectrum run (0/65). Consequences: ★ variants = the operational model, the
  910.74 ones = retired; the operational Rayleigh constant is **dark-free** (it matches
  `v2.2_l55s008`, not `_l55s008dark`); CHM15k (1064 nm, no WV term) and CL31 (wide 6 nm line,
  insensitive: all its spectra agree within 0.5 %) are unaffected — only the CL61 moved.
  ⚠ **The local `NC_DIR` copy (`C:/Users/hervo/Downloads`) predates the switch** (nights to
  2026-08-06, no `calibration_version`), so the dashboard's `v2.0`/`cloudWV` variants still show
  the retired configuration until those files are refreshed from the bucket.
- Two other CL61 wavelengths coexist ON PURPOSE: 910.0→molecular Rayleigh reference (0.33 %
  effect), 910.55→η tables.
- **Dashboard v3** = `inter-comparison_dashboard/index.html` (one page, 3 sites, per-instrument
  method/variant/dark controls, all-variants ladder table, per-config Hopkin, binned-PWV panel);
  architecture + editing guide in `inter-comparison_dashboard/README.md`; **`variants_v3.py` is
  the file to edit when a new run lands**. Anti-fossil rule: state-dependent numbers in page
  prose are injected at build (`prose_tokens`), never hard-coded.
- **Noise filter (2026-08-24)**: « Filtre bruit » = SNR≥3 admission masks precomputed by
  `nf_v3.py` from the NATIVE L1 files (`l1_l2_io.read_l1_native` — the `_streams_*` cache is
  HOURLY, never use it for noise statistics), windows 5 min (default)/30/60/3 h, four modes
  (par instrument / intersection / moyenne-d'abord / masque scène, referee = the site CHM15k,
  `nf_ref` in variants_v3; scene thresholds = ABSOLUTE β_att 0.1/0.25/0.5 Mm⁻¹sr⁻¹ measured
  AND detected by the referee — operator decisions 2026-08-24). SNR filtering conditions the
  sample on signal: at Payerne 5 min the
  CL31 keeps 25 % of band hours and the measured sampling bias (CHM15k median under the CL31
  mask) is **+26 %** — shown live on the page; the scene/intersection/aggregate modes are the
  unbiased comparisons. Masks frozen at site-default constants; L2 inherits L1 masks;
  curtains/PWV stay static. Verification: `check_v3.py` (gate 0.1 %).
- **CHM15k overlap option (2026-08-25)**: dashboard checkbox « recouvrement T (Hervo 2016) »,
  CHM15k rows only — the temperature-dependent overlap residual, rcs₀·(1 + (a(z)·T_int + b(z))/100),
  kernel `validation/paper/overlap.py` (validated vs MATLAB), models one NetCDF per instrument in
  `ALC_OVERLAP_DIR` (default `D:/TEMP_MODELS/202606`; all 6 CHM15k of the 3 sites covered).
  Applied after the dark, before the C_L division, and to **L2 as well** (dimensionless factor,
  same optics — unlike the rcs₀-unit dark). It is a **no-op above ~700–840 m** (a and b vanish at
  full overlap), so it reshapes the near range only and moves the 500–3000 m band median by
  +0.01 %: the ladder/statistics and the dark-free constants pairing stay coherent. Verified in
  `check_v3.py` against the reference kernel (gate 1e-5; observed 2.1e-06 = a/b storage precision).
- Cloud calibration is dark-immune (+0.08 % CL31, ~0.001 % CL61, measured) and its fixed
  100–2400 m window makes a static dark CBH-flat by construction — never rerun cloud for dark.
- L2 as distributed: CHM15k **and CL61** = operational Rayleigh v1.0 (the CL61's v1.0
  calibration went live recently: constant 1.0 default on 48 d then 2.1257 from mid-2026-06 on
  the paired window — the transition is visible in the L2 row tooltip); CL31 = UNCALIBRATED
  (flat default 1e8 — the gap the cloud calibration fills).

## CHM15k temperature coefficient & diurnal bias (2026-10-20, report 19)

- Notation: **κ** = temperature coefficient d ln C_L/dT_int (%/K). Never call it β (β = backscatter).
- The CHM15k Rayleigh constant rises with **`temp_int`** by **+1.3…+2.1 %/K**. Estimators: synoptic
  anomalies (season removed), intra-night with elapsed time held fixed, the 11-yr archive, and
  same-wavelength pairs. It is unit-specific (0…+3.5) and grows with module generation (TUB≤15 ≈ +1,
  TUB20–24 ≈ +2.5 %/K). Never regress C on raw `temp_int`: it is ~70 % season. Use anomalies.
- Calibration is night-only, so **daytime L2 is biased high** by κ×(noon T_int rise): median +4.8 %,
  +7.5 % in summer, >10 % for ⅓ of units. The night κ predicts the day (Schiphol/Budapest/Măgurele pairs
  r = 0.98).
- **Mechanism** = multiplicative optical-throughput change tied to the HOUSING temperature (no lag, τ ≈ 0).
  Leading hypothesis: laser line vs the **1 nm** receiver filter (Lufft datasheet; the 0.1 nm on the Köln wiki is
  unsupported). Co-equal rival: far-field alignment drift (FOV 0.45 mrad ≈ divergence). The manual puts the filter
  INSIDE the heater-only laser-optical module (LOM = transmit+receive unit), so the filter theory needs a heat
  leak; only Lufft can decide (11 questions, report 19 §7.8). Working assumption (operator, 2026-10-20):
  **W = 1 nm** → the fit needs a filter drift **α ≈ 20 pm/K** (IQR 15–24) = exactly Jenoptik's published
  JENOFILT coefficient (ΔCWL = (T−22)×0.02 nm; Jenoptik built the CHM15k before Lufft) → a conventional
  (cemented/evaporated) filter; hard-coated ones (2–5 pm/K) would give only 0.2–0.5 %/K. Laser energy unlikely (PRF −0.06 %/K). Excluded:
  detector (Peltier), `p_calc` reference, additive dark, overlap, aerosol (AERONET: slope +2.36 → +2.01).
  `status_laser` is stuck at 100 %. `temp_lom` matters only when the heater-only LOM overheats on hot
  afternoons (unit-specific sign).
- κ normalisation before the Kalman shrinks the seasonal amplitude: network median ±16 → ±10 %, 35 % removed
  (Budapest B ±20 → ±7 %; Kleine Scheidegg weak, ±18 → ±14 %: old module, small swing). κ by module generation:
  ≤2015 +0.95, 2016–19 +1.65, ≥2020 +2.0 %/K. Breaks: 2015 (dark background ÷4), ~2020 (detector set-point
  24.85 → 29.85 °C, LOM 33.85 → 34.35–34.85 °C). Pulse count cannot proxy pulse energy (±0.5 %, −0.07 %/K).
- **Annual cycle attribution (q15, §7.13):** κ×T_int explains only 25–45 % of the ±16 % annual cycle. The rest
  (~±9 %) is a height-independent GAIN, common to co-located CHM15k (pairs cancel it), absent in the CL61, peaking
  near the SOLSTICE (DOY ~150–177, vs T_int peak ~206). q16 (§7.14): it is a REAL CHM15k-specific gain (±17 %),
  every year since 2014 in the independent MATLAB/NCEP archive, same in every fixed window 3–8 km, not twilight,
  not insolation/humidity/LR/noise/code/raw2l1 (bit-identical beta_raw)/firmware; absent in the CL61. Apply κ only to
  SUB-DAILY T deviations; keep an annual term. Candidates: sun on the biaxial head (alignment) or firmware → Lufft Q12.
  Code bug: v2.2 fallback tier nights +9 % (33 % of June nights).
  Excluded (each ≤ ~1.5 %): aerosol in/below window, molecular model (−1.6 %, anti-phase), H2O at 1064 (0.1 %),
  laser/PRF, p_calc (division works), window transmission, noise/dark/afterpulse/overlap, window selection.
  `laser_life_time` wraps at 65 536 h (16-bit). Daytime residual solar baseline +6 %-pt at 4 km (L2 aloft only).
- The simple exp model beats the v4 housekeeping ML (bagged tree): 13.0 vs 14.2 % operationally, 20.2 vs
  22.3 % HK-only. The T-correction cuts the day−night bias vs CL61/EARLINET from 7–23 to ≤ 4.5 pts at 6 of
  7 comparisons.
- **Payerne CHM15k (TUB140016) is an outlier.** Its receiver is noisy (night bg 5–25× the network), so its
  night κ is unreliable, and its LOM loses regulation on 37 % of daytime hours with sign −2 %/K. Never
  apply its night κ. Screen archive outliers PER MODULE (serlom), never per station.
- `calibration_pulse` drops −15 % (−24 % summer) by day with solar background. That is a real detector
  loss; the firmware divides it out (`beta_raw = ((P/lp)−b)/(cs·o(r)·p_calc)·r²`). Not an artefact.
- CHM15k L1 carries `optical_module_id` (TUB…) per file: use it to split eras. Schiphol A–D, Budapest A/B
  and Măgurele A/B share one station position in L1 and L2.
- Cache: `C:/DATA/Projects/202610_CHM15k_stability/hk/` (10-min HK + binned clear rcs, all 154 CHM15k +
  4 CL61). Code: `validation/chm15k_stability/`.

## Working preferences (user hervo63)

- **Ask before relocating/moving data**; when proposing to write somewhere, state *what* / *how much*
  / *where* up front. In-place edits and computation are fine without asking.
- `/home/pay/...` is **horribly slow** — prefer `/tmp` for transient/scratch files.
- Wants honest trade-off analysis ending in a **clear recommendation**, not a menu of options.
- Comfortable in French (often writes FR, mixes EN).
- **The user pushes git themselves** — prepare commits, do not `git push`.
