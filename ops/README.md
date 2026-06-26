# Operational daily pipeline (Linux server)

Each day this calibrates the previous day's E-PROFILE ALC data (**Rayleigh + liquid-cloud + Kalman**),
downloads the **CAMS** aerosol file it needs, and **incrementally** refreshes the static dashboard.

```
cron ──> ops/run_daily.sh ──> ops/ops_daily.py
   (glue: config, venv,          1. fetch CAMS for D-1   (ADS download, retried)
    flock, log, alert)           2. calibrate D-1         (scripts/run_all_l1_2026.py, all streams)
                                 3. dashboard             (scripts/build_dashboard.py --changed-only)
                                 4. heartbeat             (.last_success)
nginx/apache ──> serves $ALC_DASHBOARD_DIR  (HTTP, separate path)
```

Everything is driven by `ALC_*` env vars set in **`ops/config.sh`** — the only file you edit to migrate.

## One-time setup

1. **Clone + venv + deps**
   ```bash
   git clone <repo> ~/ALC_rayleigh_calibration && cd ~/ALC_rayleigh_calibration
   python3 -m venv .venv && . .venv/bin/activate
   pip install -r requirements.txt
   ```
   `cfgrib` converts the CAMS GRIB to netCDF in pure Python — **no conda, no system eccodes, no
   `grib_to_netcdf` CLI** needed.

2. **CAMS / ADS credentials** — put your ADS key in `~/.cdsapirc` and accept the CAMS licences once.
   Smoke-test the download:
   ```bash
   . ops/config.sh
   python -c "from calibration.io.cams import ensure_cams_file as e; print(e('$ALC_CAMS_DIR','20260401',auto_download=True))"
   ```

3. **Edit `ops/config.sh`** — paths (L1, CAMS, fullcal, dashboard), venv, proxy (if the box needs one
   for the ADS), workers, optional alert email. This is the only file to change.

4. **First full dashboard build** (one-time, ~minutes; daily updates after are incremental):
   ```bash
   . ops/config.sh
   python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR"
   ```
   No calibration output yet? Seed a few days first: `python ops/ops_daily.py --day YYYYMMDD` (repeat),
   or a range with `scripts/run_all_l1_2026.py --start … --end … --per-type 0 --methods rayleigh,cloud`.

5. **Serve it** — point your web server at `$ALC_DASHBOARD_DIR`, e.g. nginx
   `location /alc/ { alias /…/dashboard/; index index.html; }`. Serving over **HTTP** (not `file://`)
   is what makes the QC-flag persistence reliable.

6. **Schedule** — `crontab -e` (no sudo needed):
   ```
   30 6 * * *  /home/you/ALC_rayleigh_calibration/ops/run_daily.sh
   ```
   06:30 is an example — choose a time after the previous day's L1 files **and** CAMS are available.

## Behaviour
- **`ALC_DAY_LAG=1`** → process yesterday (CAMS lands the next day).
- **`ALC_BACKFILL_DAYS=5`** → also retry the last 5 not-yet-processed days, so a missed night (server
  down, late CAMS) self-heals. Done days are recorded in `$ALC_DASHBOARD_DIR/.processed_days`;
  CAMS-missing days are retried until CAMS appears.
- **Lock** — `flock` prevents overlapping runs. **Logs** — `ops/logs/daily_<ts>.log` (last 30 kept);
  failures mail `$ALC_ALERT_EMAIL` if set. **Heartbeat** — `$ALC_DASHBOARD_DIR/.last_success` (UTC);
  alert if its age exceeds a day to catch a silently-stopped cron.

## Manual operations
```bash
ops/run_daily.sh                          # exactly what cron runs
python ops/ops_daily.py --day 20260401    # (re)process one specific day
python ops/ops_daily.py --no-dashboard    # calibration only
python ops/ops_daily.py --no-publish      # build the dashboard but skip the EWC publish step
python ops/ops_daily.py --force-all       # ignore the processed-days record
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR"  # full rebuild
ops/publish.sh                            # publish the current build to the EWC (images->bucket, HTML->VM)
```

## Publish online (European Weather Cloud)

Optional: after each successful build, push the site online. Set in `ops/config.sh`: `ALC_PUBLISH=1`
plus the bucket/VM targets. The split is **images → public S3 bucket**, **HTML → web VM (nginx)**; the
bucket base URL is baked into the HTML via `ALC_IMG_BASE_URL` so pages point their images at the bucket.

```
ops/ops_daily.py ──(after a successful build)──> ops/publish.sh
                                                  ├─ rclone sync diag/ ombsens/ flagex/ -> $ALC_S3_REMOTE:$ALC_S3_BUCKET
                                                  └─ rsync HTML+assets (no images)      -> $ALC_VM_RSYNC_TARGET (nginx docroot)
```

**One-time setup**
1. **Bucket** `eprofile-alc-dashboard` (created on cci2). Apply the public-read policy (anonymous
   `s3:GetObject`, no listing) once: `s3cmd setpolicy ops/ewc_bucket_policy.json s3://eprofile-alc-dashboard`
   (or via awscli `s3api put-bucket-policy`).
2. **rclone remote** pointing at `https://object-store.os-api.cci2.ecmwf.int` (provider `Other`, your
   Cypher S3 keys); set `ALC_S3_REMOTE` / `ALC_S3_BUCKET`.
3. **Web VM** with nginx serving the `$ALC_VM_RSYNC_TARGET` path; authorise the ops host's SSH key.
4. Set `ALC_IMG_BASE_URL=https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/`.

**Cutover** — switching from relative to absolute image URLs needs ONE full rebuild so every page is
rewritten, then a publish:
```
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR" \
       --img-base-url "$ALC_IMG_BASE_URL"
ops/publish.sh
```
After that the daily `--changed-only` build + publish is incremental (changed pages rsync'd, new images
synced; the first bulk image upload is tens of GB so consider seeding it from a fast machine).

**Behind a proxy** (e.g. MeteoSwiss server 434): the S3 client (`aws`/`rclone`) picks up
`https_proxy`/`http_proxy` from the environment for the image upload; for the HTML rsync set `ALC_VM_SSH`
to an ssh command carrying the key and a ProxyCommand, e.g.
`export ALC_VM_SSH='ssh -i ~/.ssh/EWC -o ProxyCommand="connect -S proxy.meteoswiss.ch:1080 %h %p"'`.
Pick the upload client with `ALC_S3_TOOL` (`rclone` with `ALC_S3_REMOTE`, or `aws` with `ALC_AWS_PROFILE`
+ `ALC_S3_ENDPOINT`). If the host can reach S3 over HTTPS but cannot SSH out to the VM at all, switch HTML
delivery to a pull model instead: `rclone sync` the HTML to an `html/` prefix in the bucket and run
`rclone sync <remote>:<bucket>/html /var/www/alc` from a cron on the VM (VM↔bucket is inside EWC, fast).

Publish failures are **non-fatal** (logged + flagged in `.last_success` as `publish=fail`) so they never
mask a good calibration.
