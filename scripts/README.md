# scripts/

Runnable entry points: production / batch calibration runners and data-preparation
tools. **Run them from the repository root** — they read `options.json` /
`instruments.json` from the working directory. Figures and result archives are written
to `C:\DATA\Projects\202606_E-PROFILE_calibration` (configured near the top of each
script).

## Calibration runners
| Script | Purpose |
|---|---|
| `run_all_l2monthly.py` | Calibrate every CHM15k / CL61 / Mini-MPL station in the E-PROFILE L2-monthly archive (resumable, parallel) |
| `run_calibration.py` | Calibrate one representative instrument per type (CHM15k / CL61 / Mini-MPL) |
| `run_targeted_rayleigh.py` | Calibrate specific WMO + identifier keys |
| `calibrate_cloudnet_cl61.py` | Rayleigh + liquid-cloud calibration of Cloudnet CL61 raw (Lindenberg, Hyytiälä) |
| `run_lindenberg_cl61_cal.py` | Lindenberg CL61 Rayleigh + cloud calibration with WV correction and Kalman smoothing |
| `run_all_loop.ps1` | PowerShell wrapper that re-runs `run_all_l2monthly.py` until the archive is complete |

## data/ — acquisition & preparation
| Script | Purpose |
|---|---|
| `download_cloudnet_cl61.py`, `download_cloudnet_longrun.py`, `download_lindenberg_missing.py` | Fetch CL61 raw NetCDF from the ACTRIS-Cloudnet API |
| `homogenize_cl61_daily.py` | Consolidate per-minute Cloudnet CL61 day-folders into one daily file |
| `convert_minimpl_l2.py` | Reconstruct daily L1 RCS files for the Toulouse Mini-MPL from L2-monthly |
| `build_l2_manifest.py` | Scan the L2-monthly archive → `stations_l2_manifest.json` (project data dir) |
| `build_cop_lookup.py` | Extract operational calibration constants → `cop_lookup.json` (project data dir) |
| `get_aeronet.py` | Fetch AERONET AOD (interpolated to 532 nm) |
| `export_matlab_daily.m` | (MATLAB) export per-instrument MATLAB Rayleigh results to CSV |

## doc_tools/
| Script | Purpose |
|---|---|
| `md_to_pdf.py` | Convert a Markdown report to PDF |
| `embed_images.py` | Inline local images as base64 data URIs for a portable Markdown report |
