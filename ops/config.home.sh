#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------------
# HOME-PC (Windows + WSL) config for building & publishing the ALC dashboard to the EWC.
# It replaces the server config.sh assumptions: NO MeteoSwiss proxy, local /mnt paths, aws
# 'ewc' profile. Machine-specific -> keep it local (add to .gitignore); do not deploy on the server.
#
# Publish leg (runs in WSL, sourced by publish.sh via the ALC_CONFIG override):
#     wsl -- bash -lc "cd /mnt/c/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration \
#                      && ALC_CONFIG=ops/config.home.sh bash ops/publish.sh"
#
# The BUILD + opcoeff + (optional) v1 reprocess run with the WINDOWS venv (.venv/Scripts/python.exe)
# using Windows paths (C:/... , D:/...). Only this publish leg uses the /mnt/... WSL paths below.
# ---------------------------------------------------------------------------------------------

# --- where the built site + reprocessed fullcal live (WSL view; C:\ == /mnt/c) ---------------
export ALC_DASHBOARD_DIR="${ALC_DASHBOARD_DIR:-/mnt/c/DATA/Projects/202606_E-PROFILE_calibration/dashboard_home}"
export ALC_FULLCAL_DIR="${ALC_FULLCAL_DIR:-/mnt/c/DATA/Projects/202606_E-PROFILE_calibration/fullcal_home}"

# --- images -> EWC public S3 bucket (aws v2; add the 'ewc' profile to WSL ~/.aws first) -------
export ALC_S3_TOOL=aws
export ALC_AWS_PROFILE=ewc
export ALC_S3_ENDPOINT="https://object-store.os-api.cci2.ecmwf.int"
export ALC_S3_BUCKET="eprofile-alc-dashboard"
export ALC_IMG_BASE_URL="https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/"
export ALC_IMAGES_IN_BUCKET=1     # HTML references bucket URLs; local diag PNGs pruned after upload
                                  #   set 0 to KEEP local PNGs (home has disk); then upload still happens

# --- HTML -> EWC web VM (nginx docroot). Home reaches 136.156.139.31 directly -> NO proxy. -----
export ALC_VM_RSYNC_TARGET="hem@136.156.139.31:/var/www/alc"
export ALC_VM_SSH='ssh -i ~/.ssh/EWC'

# --- explicitly clear any MeteoSwiss proxy (home has open egress) ------------------------------
unset https_proxy http_proxy ALL_PROXY
