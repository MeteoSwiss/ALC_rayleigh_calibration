#!/bin/bash
# Environment of the eta-study array on CSCS Linard (sourced by every sbatch in this folder).
# Code is STAGED under $ETA_BASE/code (the Linard $HOME is over quota, so no git clone there);
# inputs are the existing /scratch mirrors -- nothing is copied by the jobs themselves.
export ETA_BASE=/scratch/mch/mhrvo/eta_study
export ETA_CODE=$ETA_BASE/code
export ETA_OUT=$ETA_BASE/incloud                     # per-unit extraction (NPZ) + units/*.json
export ETA_LUT=$ETA_BASE/forward_lut_ext.npz         # PVC forward table (window base-100..+380 m)
export ETA_L1_ROOT=/scratch/mch/mhrvo/E-PROFILE_L1   # E-PROFILE L1 mirror (<wmo>/<yyyy>/<mm>/L1_*.nc)
# 0.4 deg monthly CAMS ONLY. Never /scratch/mch/mhrvo/CAMS_Monthly (1 deg, orography bias): the
# extraction script refuses a 1-deg folder name.
export ETA_CAMS=/scratch/mch/mhrvo/CAMS_Monthly_04
# keylist: all 270 network units, or a pre-set ETA_KEYLIST (e.g. keys_linard.txt = the units the
# home PC has not already processed)
export ETA_KEYLIST=${ETA_KEYLIST:-$ETA_BASE/keys_network.txt}
export PY=$HOME/miniforge3/envs/alc/bin/python       # conda env 'alc' (numpy/scipy/netCDF4/xarray)
export PYTHONPATH=$ETA_CODE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ETA_CODE"
