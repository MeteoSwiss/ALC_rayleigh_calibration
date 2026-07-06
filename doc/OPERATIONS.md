# E-PROFILE ALC calibration — Operations & architecture (host `zueub434`)

Operator guide for the **live** E-PROFILE Automatic-Lidar-and-Ceilometer (ALC) calibration +
dashboard system as it runs on **`zueub434.meteoswiss.ch`** after the **2026-06 migration** off the
inode-exhausted NAS (`/mnt/amaroc_data`) to `/data/zue`.

> Scope. This is the *operations* doc: where things live on this host, what the daily pipeline does,
> how to run/recover it, and the host-specific gotchas. For the science/method see `README.md` (top
> level) and `doc/README.md`; for the one-time EWC publish setup see `doc/reports/ewc_dashboard_deployment.md`
> and `ops/README.md`. Paths and values below were read from the live `ops/config.sh` on 2026-06-28.

---

## 1. Overview

Each day the pipeline calibrates the **previous day** of E-PROFILE ALC data across the whole network
(Rayleigh / molecular + liquid-cloud + Kalman smoothing), plus two windowed products (instrument
**sensitivity** and **OmB** = observation-minus-CAMS-background), then refreshes a static **dashboard**
and publishes it to the **European Weather Cloud** (EWC). The aerosol/thermodynamic background comes
from a daily **CAMS** download.

```
cron 0 15 * * * -> ops/run_daily.sh -> ops/ops_daily.py
  (15:00, live)      (glue: config.sh,      0. refresh census            (scripts/refresh_census.py, merge-only)
                     venv, flock, log,      1. fetch CAMS for D-1        (ADS download, retried; all regional boxes)
                     alert mail)            2. calibrate D-1 network     (scripts/run_all_l1_2026.py --sens --omb)
                                            3. update_opcoeff            (extract_l2_opcoeff.py -> operational_coefficients.csv)
                                            4. dashboard (incremental)   (scripts/build_dashboard.py --changed-only)
                                            5. publish to EWC            (ops/publish.sh: images->S3, HTML->web VM)
                                            6. heartbeat                 (.last_success)

Live site:  https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/
```

Everything is driven by `ALC_*` environment variables in **`ops/config.sh`** — the single file you
edit. Dates are UTC. Repo root on this host:
`/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code`.

---

## 2. `/data/zue` layout (from `ops/config.sh`)

The whole calibration + dashboard tree was relocated under
`/data/zue/E_PROFILE/ALC/Calibration/`. Key exported paths (source `ops/config.sh` to load them):

| Env var | Value on `zueub434` | What it is |
|---|---|---|
| `ALC_REPO` | `/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code` | this repo (code) |
| `ALC_VENV` | `$ALC_REPO/.venv` | Python venv (cfgrib, eccodes, plotly, netCDF4, pandas, cdsapi) |
| `ALC_FULLCAL_DIR` | `/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0` | calibration output: per-station `<key>/` dirs (~436), each with `<key>_cal.csv`, `_kalman.csv`, `_sens.csv`, `_omb.csv`, `_hk.csv`, npz caches, `plots/` |
| `ALC_DASHBOARD_DIR` | `/data/zue/E_PROFILE/ALC/Calibration/dashboard` | built static site (`index.html`, `stations/`, `diag/`, `ombsens/`, `calib_index.sqlite`, dotfiles) |
| `ALC_L1_ROOT` | `/data/zue/E_PROFILE/ALC/L1_FILES` | E-PROFILE L1 input: `<wmo>/<year>/<month>/L1_<wmo>_<ident><YYYYMMDD>.nc` |
| `ALC_L2_DIR` | `/data/zue/E_PROFILE/ALC/L2_FILES` | L2 archive (station name/country + operational `calibration_constant_0`) |
| `ALC_CAMS_DIR` | `/data/zue/E_PROFILE/ALC/CAMS` | CAMS cache: `CAMS_Beta_<YYYYMMDD>.nc`, fetched daily |
| `ALC_CAMS_DIR_FALLBACK` | *(unset)* | Optional coarser CAMS archive (e.g. 1°) tried per-month when the primary (0.4°) lacks a month — used by the cloud WV correction **and** OmB; both L137 model levels |
| `ALC_OPCOEFF_CSV` | `/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0/operational_coefficients.csv` | operational L2 constants for the dashboard %-of-operational ratio |
| `ALC_OLDRAY_DIR` | `/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh` | **v1.0** (old operational) Rayleigh overlay source (raw `2025/` + `2026/` yearly NetCDFs) |
| `ALC_CENSUS` | `$ALC_REPO/validation/scope_l1_2026_census.json` | station census (wmo/ident/type/lat/lon) = dashboard manifest |
| `PATH` (prepend) | `/data/zue/E_PROFILE/ALC/Calibration/tools/bin` | AWS CLI v2 installed off the slow NFS home, for `publish.sh` |

`ALC_IMAGES_IN_BUCKET=1` is set (bucket mode — see section 5).

> Note. Only the *new* tree under `/data/zue/E_PROFILE/ALC/Calibration/` is authoritative.
> `/mnt/amaroc_data/...` is the **old** NAS location and must not be written to (it is at its inode
> ceiling — the reason for the move). See the verification caveat in section 11.

---

## 3. Daily pipeline, step by step

Entry point `ops/run_daily.sh` (cron glue): sources `ops/config.sh`, activates the venv, takes a
single-instance `flock` (overlapping runs exit quietly), logs to `ops/logs/daily_<ts>.log` (last 30
kept), and on a non-zero exit mails the last 50 log lines to `ALC_ALERT_EMAIL` (`hem@meteoswiss.ch`).
It then runs `ops/ops_daily.py`, which for each **target day** does:

**Target days** = the day `D - ALC_DAY_LAG` (lag `1` -> yesterday, because CAMS lands the next morning)
plus the previous `ALC_BACKFILL_DAYS` (`5`) not-yet-processed days — a self-healing backfill after an
outage. Processed days are recorded in `$ALC_DASHBOARD_DIR/.processed_days`; CAMS-missing days are
retried on the next run.

0. **Refresh the census** (once per run, before the day loop): `scripts/refresh_census.py` scans
   `$ALC_L1_ROOT` and merges what it finds into `$ALC_CENSUS` — new streams appended, existing ones
   never dropped (coverage/metadata refreshed, atomic write) — so a newly installed station is
   calibrated the same day. Best-effort: a census failure is logged and never blocks the run.
1. **Fetch CAMS** (`calibration.io.cams.ensure_cams_file`, 3 retries, isolated so a network failure is
   one log line). Downloads `CAMS_Beta_<D>.nc` covering the night (`D-1`..`D`) from the ADS if missing.
2. **Calibrate the day** across the network:
   `scripts/run_all_l1_2026.py --start D --end D --per-type 0 --ignore-coverage --workers 6
   --methods rayleigh,cloud --force --sens --omb`.
   Runs `ALC_WORKERS=6` parallel instrument streams (per-stream subprocess timeout
   `STREAM_TIMEOUT=1800 s`). It **merges** into the per-stream CSVs (it no longer overwrites them).
   `--sens`/`--omb` additionally update the sensitivity and OmB products via per-day caches (section 4).
3. **`update_opcoeff`** (`scripts/extract_l2_opcoeff.py $ALC_L2_DIR $ALC_OPCOEFF_CSV --start D --end D`):
   appends the day's operational L2 `calibration_constant_0` to `$ALC_OPCOEFF_CSV` (the dashboard's
   %-of-operational ratio). Best-effort — never fails the run.
4. **`mark_processed(D)`** — the day is recorded as done (CAMS was present + we ran) even if some
   individual streams errored, so it is not redone.

Then, **once**, if any day was processed:

5. **Dashboard** (`scripts/build_dashboard.py --fullcal $ALC_FULLCAL_DIR --out $ALC_DASHBOARD_DIR
   --changed-only [--l2dir ...] [--opcoeff ...]`): re-renders only station pages whose `<key>_cal.csv`
   changed since the last build (marker `$ALC_DASHBOARD_DIR/.last_build`); the summary always rebuilds.
   The v1.0 Rayleigh overlay and `--opcoeff` ratios are wired in here.
6. **Publish** (`ops/publish.sh`) when `ALC_PUBLISH=1` and the build succeeded — see section 5.
   Publish failure is **non-fatal** (logged + flagged in the heartbeat).
7. **Heartbeat** `$ALC_DASHBOARD_DIR/.last_success` (UTC) records `ran=...`, `cams_missing=...`,
   `dashboard=ok|fail`, `publish=ok|fail`. Alert if its age exceeds ~1 day.

`ops_daily.py` exits non-zero **only** on a hard failure: the dashboard build failed, or there were
target days but none could be processed (CAMS unavailable for all).

**Manual operation**
```bash
cd /data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code && source ops/config.sh
ops/run_daily.sh                          # exactly what cron runs (with lock/log/alert)
python ops/ops_daily.py --day 20260401    # (re)process one specific day
python ops/ops_daily.py --no-dashboard    # calibration only
python ops/ops_daily.py --no-publish      # build dashboard, skip the EWC publish
python ops/ops_daily.py --force-all       # ignore the processed-days record (reprocess the window)
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR"   # full rebuild
```

---

## 4. Cache + regression-guard design

**Sensitivity** and **OmB** are *windowed* products: their headline figures aggregate over many days.
Re-reading the whole L1 archive every day is exactly what we avoid. Each station keeps a compact
per-day cache next to its outputs (`calibration/incremental.py`):

- `<key>/_sens_cache.npz` — one `bmin_night`/`bmin_day` column (min detectable backscatter vs altitude)
  per day. Aggregation = per-altitude nanmedian; **exact** vs a single multi-day pass.
- `<key>/_omb_cache.npz` — per-CAMS-time columns of bias + interpolated obs/CAMS profiles; aggregation
  concatenates along the CAMS-time axis and recomputes the scalar/profile stats exactly as
  `compute_omb` does.

A daily run computes only yesterday's column(s) from yesterday's L1, appends to the cache
(de-duplicating by date / CAMS timestamp so a forced re-run **replaces** rather than doubles), then
re-derives the full-period snapshot (`<key>_sens.csv` / `<key>_omb.csv` + PNG) from the cache.

**Regression guard** — `scripts/run_all_l1_2026.py::_cache_coverage_regression(key, sdir, cache_name,
output_name)`. The historic caches were built **offline** (chunked monthly jobs) over the full
2025–2026 window. A daily single-day run for a station **whose cache file is missing** would otherwise
create a brand-new 1-day cache and overwrite the rich historic `<key>_sens.csv`/`_omb.csv` with a
1-day value — silently regressing the panel. The guard returns **True -> SKIP** (leave the output
untouched, print `REGRESSION-GUARD: ...`) iff **both**: the product cache file does **not** exist, AND
a **non-empty** historic output already exists (>=1 data row beyond the CSV header). A genuinely new
station (no cache, no output) proceeds normally and starts accumulating; an existing cache is appended
to as usual. No calibration math is touched — it only gates the CSV write. Invoked at the top of both
`_do_omb` and `_do_sens`.

> Practical consequence: if you wipe a station's `_sens_cache.npz`/`_omb_cache.npz` but keep its CSVs,
> the daily run will **skip** that product (logging the guard) rather than clobber history. To rebuild
> a station's windowed product, regenerate its cache offline.

---

## 5. Bucket mode, publish targets, live URL

**Bucket mode** (`ALC_IMAGES_IN_BUCKET=1`). The bulky per-night diagnostic and OmB/sensitivity images
live in the EWC **S3 bucket** `eprofile-alc-dashboard`; the HTML references their **bucket URLs**
(absolute, under `ALC_IMG_BASE_URL`). Because the images are in the bucket, the local PNGs are
disposable:

- `monitoring/index.py::_scan_diagnostics` (+ `_diag_rows_from_csv`) lists the diagnostic set from each
  station's `<key>_cal.csv` (1 row per rayleigh/cloud calibration = 1 PNG), **not** from local PNGs, so
  the diag viewer survives even after the local images are deleted. `render.py` likewise gates the
  ombsens/diag panels on the CSV/data and emits bucket URLs (`config.IMAGES_IN_BUCKET`,
  `config.IMG_BASE_URL`).
- `ops/publish.sh` **block 1b**: after a successful S3 sync, prunes local `*_diag_compact.png` under
  `$ALC_FULLCAL_DIR` (and removes the now-empty `plots/` dirs) to reclaim disk. (`ops/delete_local_diag.sh`
  is the standalone, safer variant: it deletes per station only after confirming the bucket holds at
  least as many PNGs as local.)

**Publish targets** (`ops/publish.sh`, driven by `ALC_*`):

| Leg | Source | Target |
|---|---|---|
| Images | `$ALC_DASHBOARD_DIR/{diag,ombsens,flagex}/` | S3 bucket `eprofile-alc-dashboard` (`aws s3 sync --follow-symlinks`, endpoint `object-store.os-api.cci2.ecmwf.int`, profile `ewc`) |
| HTML + assets | `$ALC_DASHBOARD_DIR/` (excluding image trees, `calib_index.sqlite`, dotfiles) | web VM `hem@136.156.139.31:/var/www/alc` (`rsync -az --delete --chmod=D755,F644` over `$ALC_VM_SSH`) |

Image base URL baked into the HTML:
`ALC_IMG_BASE_URL=https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/`.
**Live site:** <https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/>.
(The web VM serves the HTML; the browser fetches images straight from the bucket — no CORS needed for
`<img>`. Full EWC deployment record: `doc/reports/ewc_dashboard_deployment.md`.)

---

## 6. Dashboard panels

Built by `scripts/build_dashboard.py` -> `monitoring/{index,render,charts}.py`, Jinja templates in
`monitoring/templates/` (`summary.html`, `station.html`, `flags.html`).

- **Summary (`index.html`)** — network maps (theoretical vs operational constant, %-of-operational
  ratio, OmB map), a sortable per-(station,method) table with value sparklines, country filters.
- **Station page (`stations/<key>.html`)**, per method (rayleigh / cloud):
  - **Calibration time series** (`charts.series_timeseries`): our calibrated constant +/- uncertainty,
    the operational constant line (from `--opcoeff`), and on the **Rayleigh** series one overlay —
    **v1.0** (black `x`, from `ALC_OLDRAY_DIR`), `visible="legendonly"` — **hidden until you click it
    in the legend**. The overlay is read fresh each build from yearly
    `ALC_calibration_<key><YYYY>.nc` files (recursive glob over `2025/`+`2026/`), taking
    `lidar_constant` where `calibration_method == 0` (`render.py::_load_oldray`). (The v13/v1.0.2
    test overlay was retired and removed from the code in 2026-07.)
  - **Per-calibration diagnostic viewer** (keyboard-navigable: left/right valid cals, Ctrl+left/right
    all days, up/down change station), with a QC-flag widget (flags persisted in the browser; export
    via the widget).
  - **OmB vs CAMS** panel and **instrument sensitivity / ICAO detection-limit** panel (images from the
    bucket in bucket mode).

---

## 7. CAMS dependency + regional boxes (out-of-Europe stations)

The Rayleigh (molecular + water-vapor) and liquid-cloud + OmB calibrations need a daily CAMS
model-level file. Download config (`calibration/io/download_cams_beta.py`):

- **Default area** `ALC_CAMS_AREA = "80,-30,27,45"` (N,W,S,E) — the Europe+Arctic box (lat
  **27–80 N**, lon **-30–45 E**), **0.4 deg** grid (the ADS pre-interpolates to 0.4 deg and ignores
  the `grid` keyword). North edge 80 N includes Hopen (76.5 N). ~**540 MB/day**. Keeps the legacy
  un-suffixed name `CAMS_Beta_<YYYYMMDD>.nc`.
- **Variables**: `aerbackscatgnd532` / `aerbackscatgnd1064` (aerosol attenuated backscatter, for OmB) +
  `t`, `q`, `z`, `lnsp` (temperature, specific humidity, geopotential, ln surface pressure — for the
  molecular/water-vapor correction). Model levels `1` + `38..137`; forecast lead times 3,6,...,24 h.
- A file lacking aerosol backscatter is re-downloaded (when `auto_download`) so it is OmB-usable;
  partial/corrupt downloads are deleted so a later run retries cleanly.

**Regional boxes (implemented 2026-06/07).** The census stations outside the Europe box are served by
small dedicated boxes (`CAMS_REGIONS` in `download_cams_beta.py`, each padded ~1.5 deg around the
cluster, a few MB/day): `namerica_west` (Edmonton), `ontario` (Western/London ON), `caribbean`
(Bonaire), `newzealand` (Lauder + Auckland). Each station is routed by lat/lon via `region_for()`
(the SMALLEST containing small box wins; else the Europe box; else no domain). Regional files are
suffixed `CAMS_Beta_<region>_<YYYYMMDD>.nc` and are both produced (prefetch + auto-download) and
consumed (`calibration.io.cams`) under that name. The `cams_point_too_far` guard (flag **-10**)
ensures a station can never silently pick up data from a wrong box's edge cell. The small boxes are
**lean** (`t/q/z/lnsp` only, no aerosol backscatter): out-of-Europe stations get the full Rayleigh +
cloud **calibration** but **no OmB** product. A station outside every box stays uncalibrated
(flag -10) until a region is added to `CAMS_REGIONS` — one dict entry; the prefetch cron and the
calibration routing pick it up automatically.

The **910 nm WV rule** is strict everywhere: a night whose CAMS file is missing, unusable or too far
is flagged (-4 / -10) and **never calibrated WV-free** (enforced in `calibration/rayleigh/calibration.py`
and `calibration/cloud/calibration.py`; hardened 2026-07 — the former silent WV-free fallback on an
unusable WV profile is fixed).

A prefetch cron (`2 6 * * *`, `ops/prefetch_cams.sh` → `scripts/prefetch_cams.py`) downloads ALL the
day's boxes shortly after the CAMS forecast lands, so no ADS round-trip blocks the 15:00 run.

---

## 8. Cron jobs (user `hem` on `zueub434`)

`crontab -l` as `hem` (ALC-relevant rows; the full crontab also has REM radiometer/webcam jobs):

| When (server time) | Job | Status |
|---|---|---|
| `0 15 * * *` | `.../ALC_calibration_v2.0_code/ops/run_daily.sh` (this ALC daily pipeline) | **active** (re-enabled after the 2026-06 migration checks) |
| `2 6 * * *` | `.../ALC_calibration_v2.0_code/ops/prefetch_cams.sh` (pre-download the day's CAMS boxes) | **active** |

> The former `18 18 * * *` v13 test-Rayleigh cron (repo `/proj/pay/E-PROFILE/Calibration_codes/dev/rayleigh_calibration`)
> fed the v1.0.2 dashboard overlay, which was retired and removed from the code in 2026-07 —
> **delete that crontab line**; nothing consumes its output anymore.

Camera-image thinning runs under user **`rem`** (`sudo /bin/su - rem`, then `crontab -l`):

| When | Job | Purpose |
|---|---|---|
| `30 3 * * *` | `/home/pay/users/rem/scripts/thin_old_cameras.sh` | thins Kloten/Geneva camera images older than 100 days under `/mnt/amaroc_data/PROD`: zips all images per day, keeps 1/view/5 min decompressed, deletes the rest. **Frees inodes on the NAS** (the operational reason the ALC tree was moved off it). |

**If the daily cron ever has to be re-enabled** (e.g. after a future migration), confirm first:
(a) `ops/config.sh` points at the `/data/zue` tree and a fresh shell exports the right
`ALC_DASHBOARD_DIR` (see the section 11 caveat), (b) the venv activates and ADS credentials work,
and (c) the schedule is after yesterday's L1 + CAMS are available (L1 is delivered the next
morning; 15:00 is fine). Validate with a manual `python ops/ops_daily.py --day <yesterday>` first.

---

## 9. Recovery & monitoring quick reference

- **Did it run?** `tail $ALC_DASHBOARD_DIR/.last_success` (heartbeat) and the newest
  `ops/logs/daily_*.log`. `.processed_days` lists the days already done.
- **Missed a day / outage:** the next run backfills the last `ALC_BACKFILL_DAYS` (5) automatically; or
  force one with `python ops/ops_daily.py --day YYYYMMDD`.
- **Reprocess a window:** `python ops/ops_daily.py --force-all` (ignores `.processed_days`).
- **Rebuild the whole dashboard:** drop `$ALC_DASHBOARD_DIR/.last_build` (forces a full render) or run
  `build_dashboard.py` without `--changed-only`.
- **Disk pressure:** bucket mode auto-prunes diag PNGs after publish; `ops/delete_local_diag.sh` prunes
  more conservatively (bucket-verified).

---

## 10. Gotchas (host-specific)

- **`sudo` is stdin-only to `rem`.** `hem` may run exactly `sudo /bin/su - rem` (NOPASSWD) and nothing
  else (`sudo -l`). There is **no** `sudo -u rem <cmd>`; pipe commands in:
  `echo "crontab -l" | sudo /bin/su - rem`.
- **`rem` has `set -o noclobber`** in its `~/.bashrc`. Inside a `su - rem` shell, `>` will **refuse to
  overwrite** an existing file (`cannot overwrite existing file`). Use `>>` (append) or `>|` (force) —
  or write to a fresh path.
- **`publish.sh` rc=2.** `publish.sh` propagates the `aws`/`rsync` exit code; **rsync exit 2 = protocol
  / stream error** surfaces as `publish: completed WITH ERRORS (rc=2)`. It is **non-fatal** to the
  daily run (the calibration + local dashboard already succeeded) — investigate the VM SSH / proxy, but
  it does not fail the pipeline.
- **Proxy.** The box is behind a proxy. `ops/config.sh` exports
  `https_proxy=http://proxy.meteoswiss.ch:8080` (and `http_proxy`, `no_proxy=localhost,127.0.0.1`) for
  the ADS download and the `aws`/`rclone` image upload. The **HTML rsync** uses a **different** path:
  `ALC_VM_SSH` carries `-i ~/.ssh/EWC` and a SOCKS `ProxyCommand` via
  `connect -S proxy.meteoswiss.ch:1080` (port **1080**, SOCKS — not 8080). Don't conflate the two ports.
- **AWS CLI off the NFS home.** `aws` v2 lives in `/data/zue/E_PROFILE/ALC/Calibration/tools/bin`
  (prepended to `PATH` in `config.sh`) because the NFS home is slow.
- **L1 latency.** Yesterday's L1 files are **delivered the next morning** — hence `ALC_DAY_LAG=1` and a
  15:00 schedule. Running too early means no L1 (and the day's CAMS is not available until the next day
  either).
- **`config.sh.*` backups.** `ops/` keeps timestamped `config.sh.pre_*` snapshots (pre_relocate,
  pre_bucketflag, pre_v13, ...) — handy to diff what changed in the migration; only `config.sh` is live.

---

## 11. Verification caveats

- The **newest** daily log (`ops/logs/daily_*.log`, a manual run at 15:25Z 2026-06-28) printed
  `dashboard=/mnt/amaroc_data/alc_calib` — the **old NAS** path, not the relocated
  `/data/zue/.../dashboard`. That run used a **stale environment** (old `ALC_DASHBOARD_DIR`), so it is
  not representative of the committed `config.sh`. Confirm that the shell/cron environment that will
  run `run_daily.sh` sources the **current** `ops/config.sh` (which correctly exports the `/data/zue`
  dashboard) before re-enabling cron. The live `dashboard/` tree (built 2026-06-28 19:30) is the
  correct, current one; `.processed_days` there lists 2026-06-24..27.
- Everything else in this doc was verified against the live code/config on 2026-06-28: paths in
  `ops/config.sh`; the pipeline in `ops/ops_daily.py` + `run_daily.sh` + `publish.sh`; the regression
  guard and caches in `run_all_l1_2026.py` + `calibration/incremental.py`; bucket mode in
  `monitoring/{index,render,config}.py`; the v1.0 overlay in `charts.py` + `render.py`; the CAMS
  boxes/vars/routing in `download_cams_beta.py`; the cron table (`hem` + `rem`);
  and the `rem` `noclobber` / `sudo` constraints. (Section 7 + the cron table re-verified 2026-07-06
  against `download_cams_beta.py`, `ops/ops_daily.py` and `ops/prefetch_cams.sh`.)
