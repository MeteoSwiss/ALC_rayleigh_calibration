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
python ops/ops_daily.py --force-all       # ignore the processed-days record
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out "$ALC_DASHBOARD_DIR"  # full rebuild
```
