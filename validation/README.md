# validation/

Comparison, sensitivity and reproduction studies behind the E-PROFILE ALC paper.
**Run from the repository root**; figures / tables are written under
`C:\DATA\Projects\202606_E-PROFILE_calibration\figs_paper_validation\…`. The written-up
results live in `doc/reports/`.

## MATLAB ↔ Python parity & method comparison
| Script | Purpose |
|---|---|
| `compare_matlab_python.py` | Compare MATLAB vs Python Rayleigh coefficients across stations |
| `compare_daily.py` | Day-by-day MATLAB vs Python comparison on the same nights |
| `compare_molecular_methods.py` | Compare the molecular-window methods on a night; defines display labels/colors |
| `validate_cams_molecular.py` | Validate the CAMS-molecular option vs US Standard 1976 + the MATLAB formula |
| `reproduce_report.py` | Reproduce the MeteoFrance Klett-method report (multiple code versions, via git worktrees) |

## Long-run / multi-site studies
| Script | Purpose |
|---|---|
| `longrun_methods.py` | Full-archive comparison of the molecular-window methods (CHM15k + Mini-MPL) |
| `precision_longrun.py` | Drift-insensitive precision metrics from long-run checkpoints |
| `reaggregate_longrun.py` | Robust (MAD) re-aggregation + ranking figure |
| `run_cl61_variants.py` | CL61 sensitivity matrix (L1/L2 × WV on/off) |
| `run_eprof_v10.py` | Pre-sign-error E-PROF v1.0 baseline for historical comparison |
| `compare_lindenberg_chm_cl61.py` | §7.6 Lindenberg CHM15k vs Cloudnet CL61 β_att (WV + wavelength correction) |
| `analyze_fullcal.py` | Per-station calibration summary statistics / flag breakdown |

## Diagnostics & figures
| Script | Purpose |
|---|---|
| `run_rayleigh_diag_payerne.py` | Per-night Rayleigh diagnostics for Payerne CL61 |
| `make_optimal_flowchart.py` | Flowchart figure of the `optimal` molecular-window method |
| `plot_cl_timeseries_site.py` | Lidar-constant time-series comparison plots for a site |
| `wv_wavelength_sensitivity.py` | WV-correction sensitivity to laser wavelength / FWHM |
