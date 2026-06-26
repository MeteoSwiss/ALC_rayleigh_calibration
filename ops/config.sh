# Operational configuration for the daily ALC calibration pipeline.
# Sourced by ops/run_daily.sh. Edit these for your server -- they override the in-code defaults via the
# ALC_* environment variables, so migrating to Linux is a one-file change (no code edits).

# --- repo + python environment ------------------------------------------------------------------
export ALC_REPO="${ALC_REPO:-$HOME/ALC_rayleigh_calibration}"        # this repo's location on the server
export ALC_VENV="${ALC_VENV:-$ALC_REPO/.venv}"                       # venv with deps (cfgrib, eccodes, plotly, netCDF4, pandas, cdsapi)

# --- data inputs --------------------------------------------------------------------------------
export ALC_L1_ROOT="/data/zue/E_PROFILE/ALC/L1_FILES"               # E-PROFILE L1: <wmo>/<year>/<month>/L1_<wmo>_<ident><YYYYMMDD>.nc
export ALC_CAMS_DIR="/data/zue/E_PROFILE/ALC/CAMS"                  # CAMS cache: CAMS_Beta_<YYYYMMDD>.nc (fetched daily)
export ALC_CENSUS="$ALC_REPO/validation/scope_l1_2026_census.json" # station census (wmo/ident/type/lat/lon) -- also used as the dashboard manifest

# --- calibration output (per-stream <key>_cal.csv + yearly NetCDFs) -----------------------------
export ALC_FULLCAL_DIR="/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/fullcal"

# --- dashboard (static site, served by your web server at a SEPARATE path) ----------------------
export ALC_DASHBOARD_DIR="/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/dashboard"
export ALC_L2_DIR=""                                                # optional: L2 archive for station name/country (blank = skip)
export ALC_OPCOEFF_CSV=""                                           # optional: operational-constant CSV for the comparison maps (blank = skip)

# --- publish to the European Weather Cloud (optional; all blank/0 = don't publish) ---------------
# Push the built site online: bulky images (diag/ombsens/flagex) -> a public S3 bucket, and the static
# HTML+assets -> a web VM's docroot (nginx). See ops/publish.sh + ops/README.md. The bucket base URL
# must match ALC_IMG_BASE_URL so the HTML (built by build_dashboard.py) points its images at the bucket.
export ALC_PUBLISH="${ALC_PUBLISH:-0}"                              # 1 = publish after each successful dashboard build
export ALC_IMG_BASE_URL="${ALC_IMG_BASE_URL:-}"                    # public bucket base URL baked into the HTML, e.g.
                                                                   #   https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/
export ALC_S3_REMOTE="${ALC_S3_REMOTE:-}"                          # rclone remote name (~/.config/rclone/rclone.conf)
export ALC_S3_BUCKET="${ALC_S3_BUCKET:-}"                          # bucket name, e.g. eprofile-alc-dashboard
export ALC_VM_RSYNC_TARGET="${ALC_VM_RSYNC_TARGET:-}"             # ssh host:path of the web docroot, e.g. ewc-alc:/var/www/alc

# --- pipeline behaviour -------------------------------------------------------------------------
export ALC_DAY_LAG="${ALC_DAY_LAG:-1}"                              # process day D-LAG (CAMS available next day -> 1)
export ALC_BACKFILL_DAYS="${ALC_BACKFILL_DAYS:-5}"                  # also (re)try the last N still-missing days -> self-healing after an outage
export ALC_WORKERS="${ALC_WORKERS:-6}"                             # parallel streams (slow no-sudo server -> keep modest)
export STREAM_TIMEOUT="${STREAM_TIMEOUT:-1800}"                     # per-stream subprocess timeout (seconds)
export PLOTS="${PLOTS:-1}"                                          # 1 = render diagnostic PNGs (needed for the per-calibration viewer)

# --- network (no-sudo server behind a proxy): uncomment + set if CAMS download needs it ----------
# export http_proxy="http://proxy.example:8080"; export https_proxy="$http_proxy"; export no_proxy="localhost,127.0.0.1"

# --- failure notification (optional) ------------------------------------------------------------
export ALC_ALERT_EMAIL=""                                          # blank = no mail; else run_daily.sh mails this address on failure
