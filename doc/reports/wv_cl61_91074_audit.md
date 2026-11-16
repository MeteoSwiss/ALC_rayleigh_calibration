# Audit: CL61 numbers computed with the retired WV model 910.74 nm / 1.0 nm

**STATUS: DONE (2026-11-02 ~16:50).** Flags, the inventory, the reruns at 910.55 nm on CAMS 0.4° for reports 03,
04 §6–§9, 05 (CL61 rows), 07 §4.7, the CBH studies and the WV-resolution study, and the report updates are
complete. §6 lists what was not rerun and why.

**The correct CL61 line is 910.55 nm / FWHM 0.188 nm (σ 0.08 nm)** (Vaisala datasheet B211780EN, DA10 guide
M212895EN-E). The retired Qmini model (910.74 nm / 1.0 nm) over-corrected the CL61 water-vapour (WV) absorption
6–7× (absorption ratio ×6.4–7.3 at Payerne, CAMS 0.4°). It raised the CL61 constants (Rayleigh ×1/0.81,
cloud ×1/0.86) and created the CL61 cloud dC/dCBH slope (+9.05 → +0.79 %/km at the correct line). CL31
(909.7/6.0), CL51 (910.0/3.4) and CHM15k are not affected.

Conventions in this file and in the report updates (each update block starts with `<!-- l55-update -->` and
**Update 2026-11-02**):
- **"910.74 nm (retired model)"** = the original (July/August) numbers. They mostly used CAMS 1°.
- **"910.55 nm (correct)"** = reruns at 910.55/0.188 on **CAMS 0.4°** ("D:/CAMS_Monthly_04;D:/CAMS_daily"),
  today's code. The 1° folders are never used for numbers labelled correct: at Payerne the 1° grid point is
  ≈ 894 m too high (PWV −26 %).
- Owner decision 2026-11-02 15:00: no further 910.74 controls. Where a 910.74 control on 0.4° (`l74ctl04`)
  already existed (report 04 §6 eta_impact, §7 ms_seasonal), it is shown as a middle column to quantify the
  "1° → 0.4° CAMS + code changes since July" part of the difference. Elsewhere a note states that this part is
  included in the difference.
- `DATA` = `C:/DATA/Projects/202606_E-PROFILE_calibration`.

## 1. How the affected items were dated

- `water_vapor.LASER_SPECTRUM["CL61"]` was (910.74, 1.0) from commit 1da14fb (2026-06-16) until 2026-11-02. The
  (910.55, 0.188) entry in the working copy is uncommitted.
- Cloud path: until a67f5d2 (2026-08-16) the values were hard-coded in `set_defaults` as 910.74/1.0, and no env
  override was possible.
- `validation/paper/intercompare.py` `WV_PARAMS["CL61"]` was 910.74/1.0 until a25977e (2026-10-20). It does not
  read the env var.
- MATLAB (`Documents/MATLAB/ALC/liquid_cloud_calibration/liquid_cloud_calibration.m:317`,
  `paper_val_process.m:387`) has 910.74/1.0 hard-coded since d345333 (2026-06-16).
- `validation/paper/wv_resolution_lib.py` read the raw table (bypassing the env var). It is fixed in the working
  copy (uncommitted).
- `calout_v22_04` (network v2.2, 2026-08-15/16, `validation/paper/cscs/alc_v22_network.sbatch`) does not export
  `ALC_WV_SPECTRUM`, so its CL61 rows used 910.74. The same holds for `diag_v22`, `diag_v22_04` and
  `cloud_profile_dump` (2026-08-16 11:09).
- The July report-05 numbers used CAMS 1° monthly for most months (report 05 §1.2: "monthly 1° (June 2026: daily
  0.4°)"; in July `D:/CAMS_Monthly_04` held only 2025-01/02/06/07).
- Correct-spectrum runs that already existed: `DATA/diag_v22_l55s008`, `_l55s008dark`,
  `cloud_profile_dump_l55s008`, `method_test_2026`, `diag_v22_dark_v2/v3`, and every ops/cscs sbatch since the
  cutover.

## 2. Inventory, results and verdicts

Grades: **H** = level or CBH-slope number, verdict at risk. **M** = level enters partly. **L** = yield,
precision or OmB. Every row has its flag line (82 flags in 16 reports, plus the Q-table of
cbh_slope_heterogeneity) and an update block.

| # | Report / section | 910.74 nm (retired) → **910.55 nm (correct, 0.4°)** | Script / data | Sens. | Verdict |
|---|---|---|---|---|---|
| 1 | 04 §6 | CL61 \|r(ln C,CBH)\| legacy→PVC 0.421→0.215 → **0.202→0.024**; PVC slope −5.9 → **+0.65 %/km** (C); band ratios unchanged | `_eta_impact_experiment_l55.py`, `_eta_impact_tables_l55.py`; `DATA/eta_impact_l55` (+ `_l74ctl04`) | H | PVC kept, cleaner; CL31 now also favours PVC on the documented estimator |
| 2 | 04 §7 | slope test Payerne 0.242/0.057 → **0.022/0.212**; Lindenberg 0.169/0.107 → **0.101/0.039**; Aosta 0.552/0.388 → **0.394/0.160**; cloud/ray PVC 0.93/0.98/0.87 → **1.080/0.994/0.910** | `_ms_seasonal_experiment_l55.py`, `_ms_seasonal_fix_errors_l55.py`; `DATA/ms_seasonal_l55` | H | **changed**: Payerne inconclusive (thin Rayleigh sample); Lindenberg and Aosta favour PVC |
| 3 | 04 §8 | C_L-alt Payerne CL61 0.20→0.01 → **0.065→0.265**; Lindenberg 0.53→0.33 → **0.362→0.081**; Test 2 Payerne 0.198→0.065 → **0.071/0.078**; Lindenberg 0.533→0.234 → **0.361→0.105** | `_ms_fov_experiment_l55.py`, `_ms_fov_analyze_l55.py`, `_ms_cloudnet_tests_l55.py`; `DATA/ms_fov_l55` | H | **changed** at Payerne (legacy/tie); Lindenberg PVC stronger |
| 4 | 04 §9 | β vs CHM15k, PVC 5.5 µm: Payerne +12.6 → **+15.6 %**, Lindenberg +16.8 → **+14.4 %** | `_ms_beta_intercompare_l55.py` | H | unchanged (5.5 µm best); Palaiseau CL31 anomaly (§7) |
| 4b | 04 consolidated verdict | — | — | — | note added: "every test" no longer holds at Payerne |
| 5 | 03 §2.2, §2.5 | config (910.74, 1.0) → (910.55, 0.188); FWHM documented | text | — | **changed** |
| 6 | 03 §3.1–3.3, §3.6 | "910.55 vs 910.74 < 0.5 %" → **β −11.5…−15.8 %, absorption ×6.4–7.3**; ±0.1 nm swing 27 % → **3.5 %** | `validation/wv_wavelength_sensitivity_l55.py` (0.4°; `--cams1deg` replicates July) | H | **changed** |
| 7 | 03 §4 | WV in CL61 Rayleigh: +16 % → **+3.9 % (Payerne), +2.8 % (Lindenberg)**; same 7 nights 1.227 → **1.006** (0.4°) | `validation/paper/wv_magnitude_cl61.py`, `_cl61_wavelength_experiment_l55.py` | H | **changed** (still mandatory, 3–4 % term) |
| 8 | 03 §5 | CL61 WV correction at 2 km 17–21 % → **2.6–3.5 %**; CAMS−sonde error ±2 % → **±0.5 %** | `validation/paper/_wv_resolution_l55.py` | M | direction unchanged, stronger |
| 9 | 07 §1, §2.2, §2.3, §2.5 | WV toggle Rayleigh +20.8 % / cloud −12.8 % → **+3.9 / +2.8 %**; ray/cloud ×≈0.94 | `wv_magnitude_cl61.py`; re-expression #20 | H | **changed** ("< 1 % convergence" was the WV error) |
| 10 | 07 §3.1–3.4, §3.6 | WV ≈ 18 pts → **≈ 3–4 pts** | MATLAB (not rerun) | H | **changed** |
| 11 | 07 §4, §4.3, §4.4, §4.7–4.9 | 910.74 vs 910.55 @FWHM 1.0 "+0.5 % immaterial" → **−17 % (0.4°; −18 % 1°)**; T² 0.78 → 0.96 | `validation/paper/_cl61_wavelength_experiment_l55.py` | H | **changed** |
| 12 | 07 §7 | ray/cloud Payerne/Lindenberg/Aosta/Camborne → v22_04 re-expressed **0.805/0.973/1.059/0.899** | #20 | H | **changed** (Aosta cycle not a WV-lever effect) |
| 13 | 07 §8.4 | CL61 cloud − CHM: +0.1/−1.7/+12.0/−2.9 → **−2.1/−4.1/+8.2/−7.2 %** | report-05 rerun (#14) | M | conversion verdict unchanged |
| 14 | 05 §2.3, §3.3–3.5, §4.1, §5, §6 | CL61 table (see the 05 §2.3 update); WV in-band +15–35 % → **+2.5–4.5 %**; Payerne cloud/ray split 0.7 → **7.1 pts**; Uccle ΔWV +1.1 → **−13.4 %** | `validation/paper/_paper_validation_cl61_l55.py`; `DATA/fullcal_cl61_l55`, `DATA/figs_paper_validation/paper_python_l55` | H | **changed** (§3.3, §3.5, §4.1); §3.4 unchanged |
| 15 | 08 §3.4, §3.6 | levels ×0.82 / ×0.85; Payerne CL61 −2.1 / +5.0 % vs CHM | estimates + #14 | H | conclusions unchanged |
| 16 | 01 §6.4, §6.5 | CL61 σ | not rerun | L | unchanged (σ can only decrease) |
| 17 | 02 §4.1, §4.2, §5.2, §6.3, §7 | CL61 yield, σ_SD | not rerun | L/M | adopted config unchanged |
| 18 | 09 §2.4 | WV toggle → **+3.9 / +2.8 %** | `wv_magnitude_cl61.py` | H | **changed** |
| 19 | 06 CL61 | OmB within a few hundredths | estimate + #14 | L | unchanged |
| 20 | phase4 §1, §3, §5.2 | CL61 CBH residual **+9.05 ± 1.46 → +0.79 ± 1.68 %/km**; ray/cloud **1.019 → 0.937** | `rayleigh_availability/cbh_native/rescale_cl61_l55.py`, `rescale_cl61_rayleigh_l55.py`, `analyze_cbh_l55.py`; `DATA/cbh_native_l55` | H | **changed** (residual gone; action 3 closed) |
| 21 | cbh_slope_heterogeneity Q-table, §1.3, §1.4, §2, §4, §5.4 | type effect +9 → **+0.79 %/km**; mechanism = spectral WV | idem | H | **changed** |
| 22 | cbh_hopkin_study §1.2, §1.4–1.6, §1bis, §2–§4 | 9/9 positive → **5/9, none > 2σ** | idem | H | **changed** |
| 23 | altitude_independence_audit summary, §2.2–2.4, §3.2, §4–§6 | WV term +7.42/+10.3 → **+1.8 %/km**; residual +7.69 → **+0.79** | idem | H | **changed** (η/WV not degenerate; A2 moot) |
| 24 | dark_aeronet_sonde_audit §2.4, §4, §5.2, §6 | CL61 spectral WV term ≈ 6× smaller | estimates (sonde scripts not rerun) | H | **changed** (WV-or-η answered: WV) |
| 25 | cl61_cloud_vs_rayleigh_origin §1, §2.3, §3, §4, §5 | Payerne ray/cloud 0.862 → **0.805**; "WV common-mode" refuted; λ_mol → 910.55 | #20 | H | **changed** |
| 26 | dark_clearsky_method §3, §4 | CL61 θ with the retired molecular WV model | not rerun | L | expected small |
| 27 | paper_v08_review | statements about the committed default at review time | — | — | not errors |

Validation of the re-expression (#20): rescaled vs full 910.55 reruns (`diag_v22_l55s008`) give cloud 0.989
(Payerne C, 70 d) and 0.987 (Lindenberg C, 289 d), and Rayleigh 0.979 / 0.984.

## 3. Judged NOT errors

- `cbh_slope_heterogeneity` §2.3bis–§2.3quater (deliberate spectrum arbitration).
- The runs `diag_v22_l55s008/_l55w10/_l55w01/_nowv/_nominal/_cl31l910/_cl31wieg` and the dashboard v3 ladder.
- The paper v0.9 chain (`run_paper_v08/09.py`, `run_method_test.ps1`, `run_dark_v2/v3.ps1`,
  `wv_magnitude_cl61.py`, `fig_v09_cloud_network.py`).
- `paper_v08_review` (it describes the configuration as committed at the time).
- `13_v22_release_dashboard_ewc_plan.md`.
- The other CL61 wavelengths: 910.0 (molecular, 0.33 %), 910.55 (η tables), and the Ångström ratio (0.02 %).
- Reports without CL61 WV-dependent numbers: 10, 11, 12, 14, 15, 16 (the remote-dark code has no WV term), 17,
  18, 19 (+plan), altitude_forward_model_study (it uses 910.55), rayleigh_availability (CHM15k),
  dark_network_diagnostic, cloudnet_classification_plan.
- 03 §3.1 and 07 §4.7 were deliberate comparisons, but at FWHM 1.0. Their conclusions were wrong and are
  updated.

## 4. Code-drift findings (they affect any rerun of the July experiment scripts)

1. **Averaging.** Since 1727870 (2026-07-10), `liquid_cloud_calibration_from_data` no longer averages to
   30 s / 10 m. The July scripts then run the CL61 on native 4.8 m gates: eta_impact gives 1 valid CL61 day
   instead of 121, ms_seasonal gives Payerne NaN, and Aosta cloud/ray comes out at 0.32. The wrappers re-apply
   `tests/_ceilo_reader.average_ceilo_data`. With it, today's code reproduces the July 1° values to ≈ 1e-4
   (12 Payerne/Camborne CL61 day × mode pairs).
2. **η tables.** In July `_ETA_CL31/_ETA_CL51` were the legacy Hopkin ladders and the CL61 borrowed the CL51
   table. The wrappers restore both: the ladders are copied verbatim from git 83f451f^.
3. **`options.json`** `cams_folder` now points at a non-existent `A:/CAMS_Monthly_04/`. The wrappers set the
   CAMS folder explicitly.

## 5. Reproduction (CPU only; outputs only to new `_l55` folders/files)

All report-04 wrappers take `--tag l55` (`l74ctl04` = the 910.74 control on 0.4°, `l74ctl` = 1°, invalid).
Shared helper: `validation/_wv_l55_common.py` (spectrum env, CAMS 0.4°, July η tables, averaging, figure suffix,
seeding).

| item | command | output |
|---|---|---|
| 04 §6 | `python -m validation._eta_impact_experiment_l55 --tag l55`; tables `python -m validation._eta_impact_tables_l55 <results.jsonl>` | `DATA/eta_impact_l55/`, `figs_multiple_scattering/fig_eta_impact_*_l55.png` |
| 04 §7 | `python -m validation._ms_seasonal_experiment_l55 --tag l55`, then `python -m validation._ms_seasonal_fix_errors_l55 --tag l55` | `DATA/ms_seasonal_l55/`, `fig_ms_*_l55.png`, `fig_S5_*_l55.png` |
| 04 §8 | `_ms_fov_experiment_l55`, `_ms_fov_analyze_l55`, `_ms_cloudnet_tests_l55` (`--tag l55`) | `DATA/ms_fov_l55/`, `fig_test*_l55.png`, `fig_cl_vs_altitude_l55.png` |
| 04 §9 | `python -m validation._ms_beta_intercompare_l55 --tag l55` | `DATA/ms_fov_l55/beta_intercompare.json`, `fig_beta_vs_chm15k_l55.png` |
| 05 | `python -m validation.paper._paper_validation_cl61_l55 --tag l55 --step all` | `DATA/fullcal_cl61_l55/`, `DATA/figs_paper_validation/paper_python_l55/`, `figs_l1_validation/fig_*_l55.png` |
| 03 §3 | `python -m validation.wv_wavelength_sensitivity_l55 [--cams1deg]` | `figs_wv_resolution/wv_wavelength_sensitivity_l55[_cams1deg].png` |
| 03 §4 etc. | `python -m validation.paper.wv_magnitude_cl61` | `DATA/wv_magnitude_cl61_l55.log` |
| 03 §5 | `python -m validation.paper._wv_resolution_l55 --tag l55` | `DATA/wv_resolution_l55/`, `figs_wv_resolution/fig_wv_*_l55.png` |
| 07 §4.7 | `python -m validation.paper._cl61_wavelength_experiment_l55 [--cams1deg]` | `DATA/cl61_wavelength_l55[_cams1deg]/` |
| CBH reports | `python rayleigh_availability/cbh_native/rescale_cl61_l55.py`; `..._rayleigh_l55.py`; `analyze_cbh_l55.py` | `DATA/cbh_native_l55/` |
| NetCDF check | `python -m validation._check_netcdf_outputs_l55 <folders>` | `DATA/l55_netcdf_check.log` |

## 6. Not rerun (with reason), and open items

- **MATLAB chains** (07 §2/§3, the MATLAB part of 03 §4): outside this repository, CL61 spectrum hard-coded.
  The Python reruns above replace their headline numbers.
- `_cl61_network_coherence.py`, `_cl61_aeronet_closure.py`, `_cl61_dark_baseline_probe.py`,
  `_cl61_wv_sources_probe.py`, `_cl61_methodology_experiment.py` (07 §4.3/4.4/4.8/4.9, §7, §8): they hard-code
  910.74 in function calls. The sections carry estimates or the corresponding rerun numbers (#12, #13).
- The sonde scripts of dark_aeronet_sonde_audit, `dark_from_clearsky.py`, and 01/02/06/08 (precision, yield, or
  levels that scale).
- The 1bis per-profile heatmaps and §2/§3 of cbh_hopkin_study were not regenerated; the stream-level slopes
  replace them.
- **Open anomaly (outside scope):** 04 §9 Palaiseau CL31 β bias vs CHM15k is −55 % in today's pipeline, against
  +9.5 % in July. The CL31 cloud constants are identical to 0.3 %, so this is not WV- or calibration-related. It
  comes from today's `intercompare` chain or the 0.4° CAMS for that site and needs a separate look.
- **Open:** `_ms_seasonal_experiment_l55 --tag l55` exited with rc=127 during its figure stage, after all rows
  were written (results.jsonl valid: 4644 lines, 0 unparsable). The serial-fix step re-ran the full §7
  analysis and figures successfully, so nothing is missing, but the cause of the 127 was not identified.
- The Fig. A4-1 droplet radii (`old/ms_fov/match_cache.npz`) come from 910.74-era scene selection. A new
  `DATA/ms_fov_l55/match_cache.npz` exists but was not compared.

## 7. Log

- 12:05 checkpoint (inventory only).
- 13:00 flags: 82 flag lines in 16 reports.
- 13:10 `wv_wavelength_sensitivity_l55`, first on 1° (renamed `_cams1deg`), then on 0.4° (10 s each);
  `wv_magnitude_cl61`.
- 13:15 `rescale_cl61_l55.py` (781 CL61 stream-days, 48 s); `analyze_cbh_l55.py`; `rescale_cl61_rayleigh_l55.py`
  (837 nights, 13 streams).
- 13:20 `_wv_resolution_l55` with tags l55 and l74ctl (≈ 5.5 min each; its 1° column is part of the study design).
- 13:21–13:46 first report-04 chain (l74ctl, 1°): **invalid** (averaging drift, §4). Dirs `*_l74ctl` are marked
  `INVALID_README.txt`; the invalid figures this run had just written into figs_multiple_scattering were deleted.
- 13:50 `_cl61_wavelength_experiment_l55` on 0.4° (≈ 3 min); the earlier 1° run is kept as `_cams1deg`.
- 13:47–14:57 chain 2, l74ctl04 (0.4°): eta_impact 348 s; ms_seasonal 3010 s; ms_fov 14:43–14:57.
  **Incident:** the first chain's bash loop had survived TaskStop. It ran l74ctl fov_analyze, cloudnet and beta
  on invalid inputs and started eta_impact l55 at 14:40, putting about 33 workers on the PC for a few minutes. It
  was killed; its eta_impact l55 (valid code path) completed and was reused. The 4 invalid `*_l74ctl.png`
  figures it wrote were deleted.
- NetCDF clashes: 16 ms_seasonal rows failed on shared-NetCDF writes, 15 of them CHM15k rows that also failed in
  the July cache. The seasonal wrapper now writes per-process `_w<pid>/_ray_out`. `_ms_seasonal_fix_errors_l55`
  recomputed all of them serially: l74ctl04 16/16 and l55 14/14 recovered. The l74ctl04 fix process exited at
  ≈ 15:47.
- Skipped: Lindenberg CL61 cloud 2025-01-01. It needs CAMS 2024-12-31, outside the 0.4° archive, and was not run
  on 1°.
- 14:57–15:20 chain 3: eta_impact l55 analysis, ms_seasonal l55 (1201 s, rc=127, see §6).
- 15:20 owner rules on PC use. Measured total CPU: 66 % before and ≈ 95 % after launch (shared with another
  session). Our share was capped at ≈ 11 cores: ms_fov l55 with 7 workers, report-05 recalibration 2 streams at a
  time, and 2 serial fixes; BLAS threads = 1.
- 15:20–16:38 parallel phase: ms_fov l55 1405 s; fov_analyze 565 s; cloudnet 203 s; beta 1592 s; report-05
  calibrate 2795 s, calib 2 s, validate 1826 s.
- 16:38 NetCDF integrity: 98 files under ms_seasonal_l74ctl04, ms_seasonal_l55, fullcal_cl61_l55 and
  ms_seasonal_l74ctl; **0 corrupt**. No July output was touched.
- 16:40–16:50 report updates for 04 §6–§9 (plus the consolidated-verdict note) and 05 §2.3–§6; earlier estimates
  in 07 §8.4, 08 §3.6 and 06 replaced by the measured report-05 shifts.
