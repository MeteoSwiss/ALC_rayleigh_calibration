# Documentation

Reference literature, method descriptions, and the analysis / validation reports
behind the E-PROFILE automatic-lidar-and-ceilometer (ALC) calibration work.

```
doc/
├── *.pdf, *.txt          reference papers (+ plain-text extracts) and method descriptions
├── WATER_VAPOR_AUDIT.md  audit of the water-vapour port (Python vs MATLAB / atmoslib)
└── reports/              generated analysis & validation reports (E-PROFILE ALC paper)
```

> **Note.** Large binaries — figures and result archives — live **outside** the repo
> under `C:\DATA\Projects\202606_E-PROFILE_calibration`. The reports below reference
> figures stored there.

## Reference papers

Filenames carry the first author and/or the journal DOI slug; open the PDF for the
authoritative citation. A `.txt` extract is available for the files marked ✓.

| File | Reference | Topic | .txt |
|---|---|---|:--:|
| `ao-62-4-861.pdf` | Speidel & Vogelmann (2023), *Applied Optics* **62**(4), 861 | Corrected Klett–Fernald algorithm; backscatter-retrieval **sign error** & sensitivity | ✓ |
| `Wiegner - Water-vapor-correction-amt-8-3971-2015.pdf` | Wiegner & Gasteiger (2015), *Atmos. Meas. Tech.* **8**, 3971 | Spectral water-vapour absorption correction (WAPL) — basis of `water_vapor_correction/` | |
| `Wiegner- Ceilinex-amt-12-471-2019.pdf` | Wiegner et al. (2019), *Atmos. Meas. Tech.* **12**, 471 | CeilInex ceilometer inter-comparison; attenuated-backscatter / water-vapour validation | |
| `Wiegner-amt-7-1979-2014.pdf` | Wiegner et al. (2014), *Atmos. Meas. Tech.* **7**, 1979 | Aerosol lidar / ceilometer retrieval | |
| `Kotthaus-amt-9-3769-2016.pdf` | Kotthaus et al. (2016), *Atmos. Meas. Tech.* **9**, 3769 | Recommendations for processing Vaisala CL31 attenuated-backscatter profiles | |
| `Investigating the seasonal fluctuations of the CHM15K Ceilometer calibration constant.pdf` | — | CHM15k calibration-constant seasonal drift/stability | |
| `Looschelders-Haefelin-Meteorological Applications - 2025 - Looschelders - Inter‐Instrument Variability of Vaisala CL61 …pdf` | Looschelders et al. (2025), *Meteorol. Appl.* | CL61 inter-instrument variability | |
| `LE-OConnor_egusphere-2025-6331.pdf` (+ `-supplement.pdf`) | O'Connor et al. (2025), EGUsphere preprint 2025-6331 | Liquid-cloud calibration reference (O'Connor cloud method) | |
| `Laffineur_egusphere-2026-948.pdf` | Laffineur et al. (2026), EGUsphere preprint 2026-948 | recent calibration / validation work (see file) | |

## Method-description documents

| File | What it is | .txt |
|---|---|:--:|
| `E-Profile_calibration_method_description.pdf` | E-PROFILE Rayleigh-calibration algorithm: L1 loading, cloud filtering, molecular-window fit, Klett inversion, lidar-constant, NetCDF output | ✓ |
| `June09_0940_Vogelmann_Klett.pdf` | Speidel & Vogelmann talk — "Is your aerosol backscatter retrieval afflicted by a sign error?" (history + corrected Klett form) | ✓ |
| `Report_Klett_method_implementation_E_Profile-2.pdf` | Implementation report: negative-extinction / AOD discrepancy in the E-PROFILE Klett inversion; sign-error hypothesis & consequences | ✓ |

## Generated reports (`reports/`)

Analysis & validation reports for the E-PROFILE ALC paper (M. Hervo, MeteoSwiss, 2026).
`reports/README.md` is the original MATLAB-package guide (kept for reference).

### Foundations / robust ensemble calibration
| Report | Summary |
|---|---|
| `calibration_coefficient_convention.md` | **Convention (read first).** Single source of truth for the coefficient — the Wiegner & Geiß (2012) lidar constant `C_L = RCS/β_att` — and how it is named, defined and reported across the code, output files, figures and the other reports |
| `molecular_window_detection_methods_report.md` | Core methodology: makes molecular-window detection pluggable (7 strategies); recommends `improved` as default, `optimal` for aerosol-rich scenes |
| `IMPLEMENTATION_SUMMARY.md` | New robust-ensemble calibration modules (multi-window / multi-LR uncertainty) |
| `INTEGRATION_COMPLETE.md` | Status: robust ensemble calibration integrated into the main workflow |
| `ROBUST_CALIBRATION_README.md` | Rationale + usage for GUM-compliant ensemble uncertainty estimation |
| `DOCUMENTATION_VERIFICATION.md` | Check that implementation matches documentation |

### Method comparison (molecular-window detection)
| Report | Summary |
|---|---|
| `method_comparison_multisite.md` | 7 strategies over 35 nights at Payerne / Amsterdam / EDT, ranked by stability and yield |
| `molecular_methods_longrun_report.md` | Full-archive comparison of 7 methods across 14 instruments (~5 months each); robust CV |
| `molecular_methods_longrun_report_embedded.md` | As above, with figures embedded inline |
| `molecular_window_detection_methods_report_embedded.md` | Methodology report with figures embedded inline |
| `precision_longrun.md` | Drift-insensitive precision (14 sites): separates noise from drift; `optimal` most precise |
| `ranking_robust_longrun.md` | Per-instrument MAD-based robust CV of the 7 methods (CHM15k, Mini-MPL tables) |

### E-PROF v2 optimization & network deployment
| Report | Summary |
|---|---|
| `v2_optimization_report.md` | Why clear nights fail the molecular fit, and a tuned **E-PROF v2**: the **scattering-ratio gate** (`max_scattering_ratio = 1.10`) is the dominant binding constraint (leave-one-gate-out recovers 77 % of failed L1-CL61 nights); 24 instruments × L1 + L2 |
| `network_v2_vs_v11_report.md` | Network-wide validation of optimized **v2 (config C8**, the `eprof_v2` default) vs v1.1 across **148 CHM15k / 11 CL61 / 5 Mini-MPL** streams, both levels, every clear night of 2026 |
| `cloud_optimization_report.md` | Liquid-cloud calibration gate-config sweep (CL31 / CL51 / CL61): 8 configurations of the O'Connor/Hopkin method scored on valid-calibration count and short-term variability (σ_SD) |
| `rayleigh_network_diagnosis_report.md` | Rayleigh calibration across the whole network (153 CHM15k + Mini-MPL streams, 2026): C_L time series, problematic-station diagnosis from L1 housekeeping, and L1-vs-L2 outcome |
| `calibration_outliers_report.md` | Network C_L time series and per-calibration **outlier rate** (2026): drift-aware robust flagging; optimized v2 (C8) vs v1.1, L1 + L2 |
| `calibration_v2_camsfar_modifications_2026-06-26.md` | June 2026 changes to the Rayleigh + cloud calibration and dashboard (branch `wv-correction`): the evidence behind each change, per-instrument-type validation, and the network-wide redeployment |

### Calibration stability & variability
| Report | Summary |
|---|---|
| `l1_2026_variability_report.md` (+ `_embedded`) | Per-instrument night-to-night variability from **L1 2026** for 10 CHM15k + 4 Mini-MPL + **10 CL61** (the long-run study had no CL61). Methods renamed to **E-PROF versions** (v0.25/v1.0/v1.1/v1.2/v2); **E-PROF v2** most precise (σ_SD 9.5 %), CL61 median σ_SD 7.7 %. Includes **E-PROF v1.0 (sign error)** — same C_L stability as v1.1 — and an independent **liquid-cloud cross-check** of the CL61 (agrees at ~10 %, both flag Zeebrugge). **`calipso` dropped** (no stratospheric molecular reference for a ground-up ALC) |
| `calibration_stability_report.md` | Long-term stability drivers; instrumental drift dominates (shown via WV-insensitive 1064 nm CHM15k) |
| `calibration_short_term_variability_report.md` | Daily/per-night scatter is mostly measurement noise; irreducible instrumental floor 8–20 % |
| `cl61_calibration_verification_report.md` | CL61 network (9 instruments): Rayleigh vs liquid-cloud, ±WV, ±Kalman, on L1 & L2 |
| `ambient_noise_report_20260529-30.md` | Ambient (no-hood) noise characterization for CL31 / CL61 / CHM15k at Payerne |

### CL61 investigations (offset & Rayleigh failures)
| Report | Summary |
|---|---|
| `cl61_rayleigh_investigation.md` | Root cause of C_L(Rayleigh) ≠ C_L(cloud) for the CL61: a systematic **background over-subtraction in the vendor β_att** (Payerne +13.5 % Rayleigh vs -0.6 % cloud) |
| `cl61_chm15k_offset_correction.md` | CL61 & CHM15k **electronic-offset** characterisation and correction from covered-telescope ("hood") dark measurements at Payerne; companion to the root-cause note |
| `payerne_cl61_detection_report.md` | Does cloud/fog detection explain the CL61 Rayleigh failures? **No** — the CL61 does not over-detect cloud; net yield equals the CHM15k (6/96 nights); the failures are genuinely cloudy nights (CHM15k is the stricter instrument on fog) |
| `dark_measurement_payerne.md` | Python port of the MATLAB detector-noise / dark-measurement analysis (CL31 / CL61 / CHM15k) on operational Payerne L1: noise floor and detection thresholds (uncovered telescopes → low-altitude noise is an upper bound) |

### Attenuated-backscatter validation (paper)
| Report | Summary |
|---|---|
| `paper_validation_report.md` | Main paper validation write-up — material, methods and results (2026-07-02); operational Python calibration + validation, figures in `figs_paper_report/` |
| `paper_python_validation.md` | Operational Python β_att validation at the benchmark stations: the dashboard `fullcal_l1_2026` calibration (Rayleigh `eprof_v2` + O'Connor cloud, Kalman-smoothed) reproducing the MATLAB figure layouts, with MATLAB kept as a legacy reference |
| `l1_vs_l2_validation.md` | L1-vs-L2 calibration validation (`eprof_v2`): native L1 binned to the L2 grid calibrates **identically** to L2 (Rayleigh L1/L2 median 0.985); cloud coefficients physical O(1) |
| `l1_validation_cscs.md` | L1 β_att validation applying the **CSCS Kalman C_L** directly to the native L1 `rcs_0` (no L2, no vendor files), then WV / wavelength / screening / gridding vs the CHM15k Rayleigh reference |
| `omb_payerne.md` | Observation-minus-Background spot-check at Payerne (L1 vs CAMS, 2026-05): operational L2 vs Kalman C_L; 910 nm via Ångström interpolation + WV correction, CHM15k native |

### Water-vapour & wavelength sensitivity
| Report | Summary |
|---|---|
| `attbsc_validation_technical.md` | Validation of calibrated attenuated backscatter via six independent strategies |
| `attbsc_wv_literature_review.md` | WV-absorption-correction literature review; anchors on Wiegner & Gasteiger (2015) WAPL |
| `payerne_cl61_wv_sensitivity.md` | WV correction brings CL61 to within +1.7 % of CHM15k (vs +19.3 % uncorrected) |
| `payerne_cl61_calibration_sensitivity.md` | CL61 Rayleigh vs liquid-cloud sensitivity at Payerne vs colocated CHM15k |
| `wv_fwhm_literature_review.md` | Laser-emission FWHM review (CL31/CL51/CL61): manufacturer vs measured (Qmini) |
| `wv_wavelength_sensitivity.md` | Sensitivity of the WV correction to laser wavelength (910.55 vs 910.74 nm) and FWHM |

### Operations & deployment
| Report | Summary |
|---|---|
| `cscs_omb_sens_runbook.md` | CSCS (balfrin) runbook to produce the CAMS 0.4° + OmB + sensitivity products and rebuild the dashboard; every step resumable (must run on CSCS — needs ADS auth + the cluster filesystem) |
| `ewc_dashboard_deployment.md` | European Weather Cloud deployment of the monitoring dashboard (publish-only): per-night diagnostic images to a public S3 bucket, static HTML on a small VM referencing them by URL (live 2026-06-26) |

Also in `doc/`: **`WATER_VAPOR_AUDIT.md`** — code audit of the Python water-vapour
two-way-transmission implementation vs MATLAB and ACTRIS-Cloudnet `atmoslib`; verdict:
the Python port is correct (and slightly more accurate; one latent bug fixed in Python).
