# Paper code: "The E-Profile Automatic Lidar and Ceilometer network"

This folder holds the code behind the figures, tables and numbers of the manuscript
(`doc/E-PROFILE ALC Paper_v0.8.docx`, the version shared for review; v0.9 is an identical copy).
This README maps every figure and table to the code, the input data and the command that produced it.
It also states what could not be reproduced. Provenance was checked on 2026-10-20 against the images embedded
in the manuscript and the Claude Code session logs, which reach back to 2026-07-14.

## Conventions

- **Data roots.**
  - `DATA` = `C:/DATA/Projects/202606_E-PROFILE_calibration`.
  - `PV` = `DATA/figs_paper_validation`. The figures of manuscript v0.8 are in `PV/paper_v08` and those of v0.9
    in `PV/paper_v09`.
  - E-Profile L1 data are in `D:/E-PROFILE_L1_2026`, CAMS in `D:/CAMS_Monthly_04` (0.4°) and `D:/CAMS_daily`.
- **Run Python from the repository root.** Use `python -m validation.paper.<script>` unless stated otherwise.
- **Python environment matters for byte-identical images.** The v0.8 figures (Figures 11, A4-1, A5-2 and A6-2)
  were rendered with the system Python 3.12 and matplotlib 3.10.3. With that interpreter their re-runs are
  byte-identical. The repository `.venv` has matplotlib 3.11.0, which gives the same content with a few percent
  of pixels different. The v0.9 figures were made with the `.venv`.
- **MATLAB figures.** They were made with R2025b. The scripts are archived in `matlab/`, with a provenance
  header in each file.

## Figures and tables

Status key:
- **byte-identical**: re-run on 2026-10-20 and the PNG is identical to the embedded image or its source file.
- **numbers reproduced**: the quoted values were recomputed exactly.
- **traced**: the code is known but no re-run was possible.
- **not reproducible**: code or data are missing.

| Item | Code | Command / notes | Inputs | Output | Status |
|---|---|---|---|---|---|
| Fig. 1 left (network map) | MeteoSwiss `mapmetnet` (github.com/MeteoSwiss/mapmetnet), stamped on the image | Script and station list not found on this PC | Station list on 2026-04-02 | — | **not reproducible** (ask the person who ran mapmetnet) |
| Fig. 1 right (stations per month) | `matlab/count_stations_by_month.m` | Run in MATLAB with no arguments; writes `station_counts_by_month.png` into the current folder. The paper image is a cropped screenshot. | `A:/E-PROFILE_L2_monthly`, Nov 2016 – Dec 2025 | `station_counts_by_month.png` | traced |
| Fig. 2 (processing overview) | Hand-drawn diagram | No editable source (`.drawio`, `.vsdx`) found | — | — | **not reproducible** (source file needed from its author) |
| Fig. 3 (Rayleigh example, Hohenpeissenberg) | Dashboard daily panel: `monitoring/panel.py`, payloads from `scripts/build_station_dashboard.py --payloads` | On zueub434: `python scripts/build_station_dashboard.py --key 0-20000-0-10962_0 --cal-dir $ALC_FULLCAL_DIR --status-dir $ALC_FULLCAL_DIR --start 20260912 --end 20260912 --l1-root $ALC_L1_ROOT --out <dir> --payloads`, then export the PNG from the station page (Plotly modebar) | L1 `0-20000-0-10962`, night 2026-09-12/13, CAMS | browser download | traced |
| Fig. 4 (CHM15k Rayleigh C_L ranking) | `fig_v09_chm15k_network.py` | `python -m validation.paper.fig_v09_chm15k_network` | `DATA/calout_v22_04` (network run, see Upstream runs) | `PV/paper_v09/fig_v09_chm15k_network.{png,csv}` | traced |
| Fig. 5 (daily cloud calibration, Palaiseau CL31) | `fig_v09_cloud_example.py` (uses `validation/_showcase_cloud_figure.py`) | `python -m validation.paper.fig_v09_cloud_example` | L1 Palaiseau 2026-06-12, CAMS, Cloudnet `der` (`D:/Cloudnet_Droplet_radius`) | `PV/paper_v09/fig_v09_cloud_example.png` | traced; the headline n and C are asserted against the operational pipeline |
| Fig. 6 (cloud C_L ranking) | `fig_v09_cloud_network.py` | `python -m validation.paper.fig_v09_cloud_network` | `DATA/calout_v22_04` (CL31/CL51); `PV/paper_v09/calib_dashboard_cl61` (CL61, from `fetch_dashboard_calib.py`) | `PV/paper_v09/fig_v09_cloud_network.{png,csv}` | traced |
| Fig. 7 (Kalman, Kleine Scheidegg) | `fig_v09_kalman_kleine_scheidegg.py` | `python -m validation.paper.fig_v09_kalman_kleine_scheidegg`; prints the Section 5.1.3 statistics | `PV/paper_v08/calib_dashboard/0-20000-0-06736_A_rayleigh.csv` | `PV/paper_v09/fig_v09_kalman_kleine_scheidegg.png` | byte-identical |
| Fig. 8 (network CHM15k overlap) | `validation/overlap_clearsky/fig_paper_network.py` | `cd validation/overlap_clearsky && python fig_paper_network.py <out.png>`; prints kept/rejected counts and z90 | Overlap-estimation results of `validation/overlap_clearsky`, temperature models `D:/TEMP_MODELS/202606` | `PV/paper_v09/fig_v09_overlap_network.png` | traced |
| Fig. 9 (monthly monitoring) | MDA package `monitoring_alc_monthly`: `E_PROFILE_ALC_Monthly_monitoring.m` → `plot_errors_warnings.m`, `plot_window_transmission.m`, `plot_laser_energy.m` | Operational monthly report (MeteoSwiss GitLab, MDAR); not part of this repository | Monitoring files `M:/zue-data-zue/E_PROFILE/ALC/Monitoring/`, a December, year unknown | monthly Word report | traced (outside this repository) |
| Fig. 10 (comparison with CAMS) | Same package, `E_PROFILE_ALC_Monthly_OB.m` | As Fig. 9 | L3 monthly file Hoyerswerda `0-20000-0-10495`, CAMS forecast; year unknown | monthly report | traced (outside this repository) |
| Fig. 11 + Section 6.1 numbers (cloud vs Rayleigh on every type) | `fig_method_comparison.py` on the forced-method run `DATA/method_test_2026` | Run: `powershell -ExecutionPolicy Bypass -File validation\paper\run_method_test.ps1` (34 streams in `method_test_streams.json`, Lauder later dropped: no CAMS). Figure: `python -W ignore validation/paper/fig_method_comparison.py` (system Python). Extra counts: `python -m validation.paper.method_test_flag_counts` | L1, CAMS; needs `ALC_RAYLEIGH_TYPES_EXTRA`/`ALC_CLOUD_TYPES_EXTRA` in `scripts/run_network_calibration.py` | `PV/paper_v08/fig_v08_methods.png`, `method_comparison_summary.json` | byte-identical; numbers reproduced |
| Figs. 12–16 + Section 6.2/6.3 numbers (EARLINET, Payerne, Amsterdam, Palaiseau) | `run_paper_v09.py` (builds on `run_paper_validation.py`, `intercompare.py`, `earlinet.py`, `figures.py`) | `ALC_VAL_SITE_WORKERS=8 python -m validation.paper.run_paper_v09` (log `PV/paper_v09/run_v09.log`); `--plot-only` re-plots from the caches | `PV/paper_v09/calib_improve_alc`, L1, EARLINET files, CAMS | `PV/paper_v09/fig_v09_{leipzig,magurele,payerne,amsterdam,palaiseau}.png`, `summary_v09.json` (every Section 6 number) | traced |
| Fig. 17 (CBH, Payerne) | `matlab/read_cloud_CEDA_fig17_20260130.m` (station-browser section) | MATLAB; the paper image is a screenshot of the figure | `MATLAB/ALC/eprofile-data.mat` (2026-01-29, from `readEProfileData` on `D:/E-PROFILE_L2`) | screenshot | numbers reproduced (876,241 observations) |
| Fig. 18 (network CBH climatology) | `matlab/read_cloud_CEDA.m`, section "Diurnal and Monthly Climatology - CBH" | MATLAB | Not resolved: no MAT file on disk reproduces the station set of the figure | `C:/DATA/Projects/E-PROFILE_Paper/` (missing) | traced; **input not reproducible** (regenerate from `eprofile-cloud-data-2016-2025.mat`) |
| Fig. A1 + Table A1 (emission spectra, Qmini) | not found | Values first entered the code in commit `1da14fb` (2026-06-16); see `doc/reports/03_water_vapour_correction.md` §2.5 | Qmini exports of 2026-06-02 (not found) | — | **not reproducible** (raw spectra and plotting code missing) |
| Fig. A2-1 + Appendix 2 text | `fig_appendix2_availability.py` (+ `appendix2_incremental_cache.py`) | `python -m validation.paper.appendix2_incremental_cache` (cache rebuild after the identifier fix), then `python -m validation.paper.fig_appendix2_availability`; the text is generated with every number injected | `A:/E-PROFILE_L2_hourly_monthly_all` (+ L1 daily supplement before Nov 2016) | `PV/paper_v08/appendix2/{fig_A2_data_availability.png, appendix2_stats.json, appendix2_text.md}` | traced |
| Fig. A3-1 (cloud flowchart) | `validation/_cloud_flowchart.py` | `python validation/_cloud_flowchart.py [<extra folder>]` | none | `doc/reports/figs_multiple_scattering/fig_cloud_flowchart.png` | byte-identical with matplotlib 3.10.3 |
| Table A4-1 + η sensitivity text | Operational tables `_ETA_*` in `calibration/cloud/_filters.py` (= `validation/multiple_scattering_eta.py::eta_table(type, 5.5, alpha=10e-3)`) | `python -m validation.paper.eta_alpha_sensitivity` | none | stdout | numbers reproduced |
| Fig. A4-1 (droplet effective radius) | `validation/_cloudnet_droplet_climatology.py --paper` | `python validation/_cloudnet_droplet_climatology.py --paper` | Cloudnet `der` files `D:/Cloudnet_Droplet_radius`; calibration-scene radii `DATA/old/ms_fov/match_cache.npz` (from `validation/_ms_cloudnet_tests.py` on `results.jsonl` of `validation/_ms_fov_experiment.py`) | `doc/reports/figs_multiple_scattering/fig_droplet_climatology_paper.png` | byte-identical |
| Fig. A5-1 (dark of 12 May 2026) | `matlab/dark_measurement_cl61_chm_cl31.m`, "Figure 2" section | MATLAB | L1 Payerne 2026-05-12 | — | **not reproducible exactly** (the saved script differs from the version that made the figure; see its header) |
| Fig. A5-2 + Appendix 5 numbers | `rayleigh_availability/dark_evolution.py` | `python rayleigh_availability/dark_evolution.py` (system Python) | Session windows `DATA/rayleigh_availability/dark_windows_detected.json` (`dark_windows_detect.py`), L1 Payerne, `PV/paper_v08/calib_dashboard` | `PV/paper_v08/fig_v08_dark_evolution.png`, `dark_evolution_stats.json` | byte-identical; numbers reproduced |
| Fig. A6-1 + Appendix 6 Payerne numbers | `run_paper_v09.py` (appendix part) | as Figs. 12–16 | Dark reruns `DATA/diag_v22_dark_v2` (`rayleigh_availability/run_dark_v2.ps1`, dark profile from `build_dark_v2.py`) | `PV/paper_v09/fig_v09_payerne_cl61dark.png`, `summary_v09.json` | traced |
| Fig. A6-2 + Appendix 6 network numbers | `remote_dark/cl61_cloud_dark_network.py` | Full run: `python -m remote_dark.cl61_cloud_dark_network`; the paper copy: `python -W ignore -m remote_dark.cl61_cloud_dark_network --replot` (dpi 200), then copy `doc/reports/figs_remote_dark/cl61_cloud_dark_network.png` to `PV/paper_v08/fig_v08_cl61_cloud_dark.png` | L1 CL61 2025-10-01 – 2026-10-08, census, per-unit caches `DATA/remote_dark/cl61_cloud_dark/*.npz`, `DATA/method_test_2026`, hood dark `dark_profiles_payerne_v3.npz` (`build_dark_v3.py`) | PNG + `DATA/remote_dark/cl61_cloud_dark/summary.json` | byte-identical; numbers reproduced |

## Numbers in the text without a figure

| Text | Code | Inputs |
|---|---|---|
| Section 5.1.1: uncorrected water vapour biases the CL61 constant low by 3–4 % | `wv_magnitude_cl61.py` | Runs `DATA/diag_v22_l55s008` (`rayleigh_availability/run_lambda_correct.ps1`) and `DATA/diag_v22_nowv` (`rayleigh_availability/run_diag_v22_nowv.ps1`) |
| Section 5.1.3: Kleine Scheidegg nights, outliers, 2A/C0 = 36 %, σ_m = 11 % | `fig_v09_kalman_kleine_scheidegg.py` | as Fig. 7 |
| Section 6.1: Mini-MPL 68–91 %, Lanzhot laser energy | `method_test_flag_counts.py` | `DATA/method_test_2026` |
| Sections 6.2, 6.3, Appendix 6 Payerne | `run_paper_v09.py` → `summary_v09.json` | as Figs. 12–16 |
| Appendix 2 | `fig_appendix2_availability.py` → `appendix2_stats.json`, `appendix2_text.md` | as Fig. A2-1 |

## Upstream runs (inputs of several figures)

| Run | Produced by | Used by |
|---|---|---|
| `DATA/calout_v22_04`: network v2.2 calibration, all census streams, 2025-01-01 – 2026-08-13 | `cscs/alc_v22_network.sbatch` on balfrin (CPU `postproc` node), output copied from `/scratch/mch/mhrvo/E_PROFILE_calout_v22_04` | Figs. 4, 6 |
| `PV/paper_v08/calib_dashboard`: nightly constants of the live dashboard, fetched 2026-10-09 | `fetch_dashboard_calib.py` (`python -m validation.paper.fetch_dashboard_calib <out> <WMO_ident> ...`) and, for the network CHM15k/CL61, `rayleigh_availability/background_bias_network.py` (calls `fetch_dashboard_calib.write_stream`) | Fig. 7, A5-2, `calib_improve_alc` |
| `PV/paper_v09/calib_improve_alc`: dashboard nightly values smoothed with the operational improve_alc Kalman | `build_calib_improve_alc.py` (filter: `kalman_improve_alc.py`, primitives from improve_alc_calib) | Figs. 12–16, A6-1 |
| `PV/paper_v09/calib_dashboard_cl61` | `python -m validation.paper.fetch_dashboard_calib <out> 0-20000-0-03808_C 0-20000-0-06418_B 0-20000-0-06447_B 0-20000-0-06610_C 0-20000-0-10393_C 0-20000-0-11538_B 0-20000-0-14015_B 0-20008-0-BIR_A 0-20008-0-EDT_B 0-20008-0-LAU_A 0-203-10-LNG_A 0-380-5-1_B 0-756-4-EERLCL61_A` | Fig. 6 |
| `DATA/method_test_2026`: both methods forced on 34 instruments, Feb–Jun 2026 | `run_method_test.ps1` | Fig. 11, A6-2 |
| `DATA/diag_v22_dark_v2`: Payerne Rayleigh rerun with the measured dark | `rayleigh_availability/run_dark_v2.ps1`; dark profile from `dark_windows_detect.py` + `build_dark_v2.py` | Fig. A6-1 |
| `DATA/diag_v22_l55s008`, `DATA/diag_v22_nowv` | `run_lambda_correct.ps1`, `run_diag_v22_nowv.ps1` | Section 5.1.1 |
| `DATA/old/ms_fov`: multiple-scattering / Cloudnet matching experiment | `validation/_ms_fov_experiment.py` → `validation/_ms_cloudnet_tests.py` | Fig. A4-1 |

## Rebuild order (v0.9 chain)

1. Network run `calout_v22_04` (balfrin) and the dashboard fetches. The fetches give the state of 2026-10-09;
   a new fetch returns the dashboard as it is on the day you run it.
2. `python -m validation.paper.build_calib_improve_alc`
3. `powershell -File rayleigh_availability/run_dark_v2.ps1` (Appendix 6 rerun)
4. `ALC_VAL_SITE_WORKERS=8 python -m validation.paper.run_paper_v09`
5. Figure scripts: `fig_v09_chm15k_network`, `fig_v09_cloud_example`, `fig_v09_cloud_network`,
   `fig_v09_kalman_kleine_scheidegg`, `fig_paper_network.py`, `fig_appendix2_availability`.
6. Number checks: `wv_magnitude_cl61`, `eta_alpha_sensitivity`, `method_test_flag_counts`.

## Not reproducible, or reproducible only outside this repository

- **Fig. 1 left:** the mapmetnet call and station list are missing.
- **Fig. 2:** no editable source of the diagram.
- **Figs. 9 and 10:** operational MDA monitoring package, outside this repository; the month and year shown are
  unknown.
- **Fig. 18:** the input MAT file of the figure was overwritten or made elsewhere.
- **Fig. A1 and Table A1:** the Qmini raw spectra and the plotting code are missing.
- **Fig. A5-1:** the saved MATLAB script differs from the version that made the figure.

## Text mismatches found during the check (for the authors)

- Fig. 3 caption says the grey points are "excluded by the cloud screen"; the panel legend says
  "excluded — signal outlier, no low cloud".
- Section 5.1.1: "median over 2025–2026 at Payerne and Lindenberg". The Payerne nights cover Mar–Aug 2026 only.
- Section 6.1: "eight units within ±10 % (0.91–1.10)": Aosta is 1.102. Mini-MPL "68–91 %" uses days with data
  as denominator, the other percentages use the 149 calendar days.
- Appendix 5: the CL61 range "−19 to −6 %" omits the +2.4 % session of 29 Sep; the CL31 did not decrease
  steadily (309, 109, 168, 125, 80, 9 %); the light-leak test fits 150–3000 m, not 0.5–1.2 km.
- Appendix 6: "20–25 % of the molecular" is Payerne −25.1 % and Uccle −19.6 %; Edmonton has no value.
- Table A4-1: the beam-divergence column lists manufacturer values; the model uses one divergence per type
  (`multiple_scattering_eta.py`).
