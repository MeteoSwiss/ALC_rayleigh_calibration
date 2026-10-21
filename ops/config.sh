# Operational configuration for the daily ALC calibration pipeline.
# Sourced by ops/run_daily.sh. Edit these for your server -- they override the in-code defaults via the
# ALC_* environment variables, so migrating to Linux is a one-file change (no code edits).

# --- repo + python environment ------------------------------------------------------------------
export ALC_REPO="${ALC_REPO:-/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code}"        # this repo's location on the server
export ALC_VENV="${ALC_VENV:-$ALC_REPO/.venv313}"                       # venv with deps (cfgrib, eccodes, plotly, netCDF4, pandas, cdsapi)
export PATH="/data/zue/E_PROFILE/ALC/Calibration/tools/bin:$PATH"   # AWS CLI v2 (installed off the slow NFS home) for ops/publish.sh

# --- data inputs --------------------------------------------------------------------------------
export ALC_L1_ROOT="/data/zue/E_PROFILE/ALC/L1_FILES"               # E-PROFILE L1: <wmo>/<year>/<month>/L1_<wmo>_<ident><YYYYMMDD>.nc
export ALC_CAMS_DIR="/data/zue/E_PROFILE/ALC/CAMS"                  # CAMS cache: CAMS_Beta_<YYYYMMDD>.nc (fetched daily)
export ALC_CENSUS="$ALC_REPO/validation/scope_l1_2026_census.json" # station census (wmo/ident/type/lat/lon) -- also used as the dashboard manifest

# --- calibration output (per-stream <key>_cal.csv + yearly NetCDFs) -----------------------------
export ALC_FULLCAL_DIR="/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0"

# ==================================================================================================
# v2.2 CUTOVER BLOCK -- uncomment ALL FOUR LINES TOGETHER, never individually.
# ==================================================================================================
# This is the flip. Procedure, prerequisites and rollback: doc/OPERATIONS_v22_cutover.md.
#
# Why one block: the daily runner MERGES into the per-stream CSVs, so v2.2 code writing into the
# v2.0 tree interleaves the two vintages irreversibly and destroys the archive you would roll back
# to. The method, the water-vapour spectrum and the output tree must therefore change in the SAME
# edit -- which is also why they ship commented out: deploying the code must not flip anything.
#
# Do NOT uncomment before scripts/check_v22_archive.py passes against the new tree.
# Rollback = re-comment these lines (or restore ops/config.sh.pre_v22_<date>), then rebuild the
# dashboard WITHOUT --changed-only and publish.
#
export ALC_FULLCAL_DIR="/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.2"
export ALC_MOLECULAR_METHOD="eprof_v2.2"
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export ALC_OPCOEFF_CSV="/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.2/operational_coefficients.csv"
#
# Deliberately NOT set: ALC_DARK_PROFILE. The measured hood dark exists only at Payerne, and the
# operator's decision (2026-08-17) is to keep the network homogeneous; Payerne CHM15k therefore
# keeps its known ~-24 % near-top-of-window bias (phase4_network_validation.md §4.4).
# ==================================================================================================

# --- dashboard (static site, served by your web server at a SEPARATE path) ----------------------
export ALC_DASHBOARD_DIR="/data/zue/E_PROFILE/ALC/Calibration/dashboard"
export ALC_L2_DIR="/data/zue/E_PROFILE/ALC/L2_FILES"                                                # optional: L2 archive for station name/country (blank = skip)
export ALC_OPCOEFF_CSV="${ALC_OPCOEFF_CSV:-/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0/operational_coefficients.csv}"                       # optional: operational-constant CSV for the comparison maps (blank = skip). Guarded: the v2.2 cutover block above must win.
export ALC_OLDRAY_DIR="/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh"   # old operational Rayleigh (v1) overlay
export ALC_CEDA_LINKS="$ALC_REPO/validation/ceda_links.json"     # {key: CEDA-L2 URL} for the per-page CEDA link (committed; blank = skip). Refresh occasionally: python scripts/build_ceda_links.py --out "$ALC_CEDA_LINKS"

# --- publish to the European Weather Cloud (optional; all blank/0 = don't publish) ---------------
# Push the built site online: bulky images (diag/ombsens/flagex) -> a public S3 bucket, and the static
# HTML+assets -> a web VM's docroot (nginx). See ops/publish.sh + ops/README.md. The bucket base URL
# must match ALC_IMG_BASE_URL so the HTML (built by build_dashboard.py) points its images at the bucket.
export ALC_PUBLISH="${ALC_PUBLISH:-1}"                              # 1 = publish after each successful dashboard build
export ALC_IMG_BASE_URL="${ALC_IMG_BASE_URL:-https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/}"                    # public bucket base URL baked into the HTML, e.g.
                                                                   #   https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/
export ALC_S3_BUCKET="${ALC_S3_BUCKET:-eprofile-alc-dashboard}"                          # bucket name, e.g. eprofile-alc-dashboard
export ALC_S3_TOOL="${ALC_S3_TOOL:-aws}"                        # image-upload client: rclone | aws
export ALC_S3_REMOTE="${ALC_S3_REMOTE:-}"                          # rclone remote name (ALC_S3_TOOL=rclone), ~/.config/rclone/rclone.conf
export ALC_AWS_PROFILE="${ALC_AWS_PROFILE:-ewc}"                      # aws profile (ALC_S3_TOOL=aws), e.g. ewc
export ALC_S3_ENDPOINT="${ALC_S3_ENDPOINT:-https://object-store.os-api.cci2.ecmwf.int}"                      # aws endpoint (ALC_S3_TOOL=aws), e.g. https://object-store.os-api.cci2.ecmwf.int
export ALC_VM_RSYNC_TARGET="${ALC_VM_RSYNC_TARGET:-hem@136.156.139.31:/var/www/alc}"             # ssh host:path of the web docroot, e.g. hem@136.156.139.31:/var/www/alc
export ALC_VM_SSH="${ALC_VM_SSH:-ssh -i ~/.ssh/EWC -o ProxyCommand=\"nc -x proxy.meteoswiss.ch:1080 -X 5 %h %p\"}"                             # ssh command for the rsync; SOCKS5 via nc (the 'connect' helper does not exist on zueub439)
# Behind a proxy (operational host is zueub439; zueub434 was the pre-2026-07-14 host) -- aws/rclone honour http_proxy/https_proxy; the VM rsync
# uses ALC_VM_SSH. Example:
#   export https_proxy="http://proxy.meteoswiss.ch:<port>"; export http_proxy="$https_proxy"
#   export no_proxy="localhost,127.0.0.1"
#   export ALC_VM_SSH='ssh -i ~/.ssh/EWC -o ProxyCommand="connect -S proxy.meteoswiss.ch:1080 %h %p"'
#   # rclone over SOCKS only (no HTTP proxy): export ALL_PROXY="socks5://proxy.meteoswiss.ch:1080"

# --- pipeline behaviour -------------------------------------------------------------------------
export ALC_DAY_LAG="${ALC_DAY_LAG:-1}"                              # process day D-LAG (CAMS available next day -> 1)
export ALC_BACKFILL_DAYS="${ALC_BACKFILL_DAYS:-5}"                  # also (re)try the last N still-missing days -> self-healing after an outage
export ALC_WORKERS="${ALC_WORKERS:-6}"                             # parallel streams (slow no-sudo server -> keep modest)
export STREAM_TIMEOUT="${STREAM_TIMEOUT:-1800}"                     # per-stream subprocess timeout (seconds)
# PLOTS was retired 2026-08-18: the per-calibration diagnostic PNGs are no longer produced (the
# interactive daily panel renders those nights from JSON payloads); the classification curtain
# PNG is now always emitted by the classify step.
export ALC_CLASSIFY="${ALC_CLASSIFY:-1}"                           # 1 = also run Cloudnet target classification (needs ceiloclass in $ALC_VENV; set 1 after install)

# --- network (no-sudo server behind a proxy): uncomment + set if CAMS download needs it ----------
export https_proxy="http://proxy.meteoswiss.ch:8080"; export http_proxy="$https_proxy"; export no_proxy="localhost,127.0.0.1"

# --- failure notification (optional) ------------------------------------------------------------
export ALC_ALERT_EMAIL="hem@meteoswiss.ch"                                          # blank = no mail; else run_daily.sh mails this address on failure

# images live in the EWC bucket; build lists diag/ombsens from CSVs, local PNGs deletable
export ALC_IMAGES_IN_BUCKET=1

# TLS CA bundle for cron egress: venv certifi lacks the ADS chain CA in the bare cron env.
export REQUESTS_CA_BUNDLE="/etc/ssl/certs/ca-certificates.crt"

# --- scratch space: keep ALL temp files off the tiny system /tmp --------------------------------
# /tmp on the ops host is a small separate volume (~3 GB) that fills up and aborts the CAMS
# download with "[Errno 28] No space left on device". Point every tool's scratch at the big /data
# share instead. CAMS_TMPDIR is read by download_cams_beta.py; TMPDIR/TMP/TEMP cover everything
# else (cfgrib/eccodes, xarray, aws cli, ...).
export ALC_TMPDIR="${ALC_TMPDIR:-/data/zue/E_PROFILE/ALC/Calibration/tmp}"
mkdir -p "$ALC_TMPDIR" 2>/dev/null || true
export TMPDIR="$ALC_TMPDIR"; export TMP="$ALC_TMPDIR"; export TEMP="$ALC_TMPDIR"
export CAMS_TMPDIR="$ALC_TMPDIR"

# Curated per-flag example images for flags.html. Without this the DAILY rebuild silently strips
# the examples the operator curated (build_dashboard only kept them when --flagex was passed).
# export ALC_FLAGEX_DIR="$ALC_BASE/flag_examples"

# Daily-calibration panel payloads (interactive station pages). OFF until the panel is rolled out:
# when =1, ops_daily regenerates data/<key>/<date>_<method>.json for panel-enabled stations after
# each calibrated day, and publish.sh ships data/ on the S3 bucket leg.
export ALC_PANEL_PAYLOADS=1
