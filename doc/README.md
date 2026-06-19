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

Journal-coded filenames are identified best-effort from the DOI slug; open the PDF for
the authoritative citation. A `.txt` extract is available for the files marked ✓.

| File | Reference | Topic | .txt |
|---|---|---|:--:|
| `ao-62-4-861.pdf` | Speidel & Vogelmann (2023), *Applied Optics* **62**(4), 861 | Corrected Klett–Fernald algorithm; backscatter-retrieval **sign error** & sensitivity | ✓ |
| `amt-8-3971-2015.pdf` | Wiegner & Gasteiger (2015), *Atmos. Meas. Tech.* **8**, 3971 | Spectral water-vapour absorption correction (WAPL) — basis of `water_vapor_correction/` | |
| `amt-12-471-2019.pdf` | *Atmos. Meas. Tech.* **12**, 471 (2019) | Ceilometer attenuated-backscatter / water-vapour validation | |
| `amt-7-1979-2014.pdf` | *Atmos. Meas. Tech.* **7**, 1979 (2014) | Aerosol lidar / ceilometer retrieval | |
| `Investigating the seasonal fluctuations of the CHM15K Ceilometer calibration constant.pdf` | — | CHM15k calibration-constant seasonal drift/stability | |
| `Meteorological Applications - 2025 - Looschelders - Inter‐Instrument Variability of Vaisala CL61 …pdf` | Looschelders et al. (2025), *Meteorol. Appl.* | CL61 inter-instrument variability | |
| `egusphere-2025-6331.pdf` (+ `-supplement.pdf`) | EGUsphere preprint 2025-6331 | calibration / inter-instrument study (see file) | |
| `egusphere-2026-948.pdf` | EGUsphere preprint 2026-948 | recent calibration / validation work (see file) | |

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

### Calibration stability & variability
| Report | Summary |
|---|---|
| `l1_2026_variability_report.md` (+ `_embedded`) | Per-instrument night-to-night variability from **L1 2026** for 10 CHM15k + 4 Mini-MPL + **10 CL61** (the long-run study had no CL61). Methods renamed to **E-PROF versions** (v0.25/v1.0/v1.1/v1.2/v2); **E-PROF v2** most precise (σ_SD 9.5 %), CL61 median σ_SD 7.7 %. Includes **E-PROF v1.0 (sign error)** — same C_L stability as v1.1 — and an independent **liquid-cloud cross-check** of the CL61 (agrees at ~10 %, both flag Zeebrugge). **`calipso` dropped** (no stratospheric molecular reference for a ground-up ALC) |
| `calibration_stability_report.md` | Long-term stability drivers; instrumental drift dominates (shown via WV-insensitive 1064 nm CHM15k) |
| `calibration_short_term_variability_report.md` | Daily/per-night scatter is mostly measurement noise; irreducible instrumental floor 8–20 % |
| `cl61_calibration_verification_report.md` | CL61 network (9 instruments): Rayleigh vs liquid-cloud, ±WV, ±Kalman, on L1 & L2 |
| `ambient_noise_report_20260529-30.md` | Ambient (no-hood) noise characterization for CL31 / CL61 / CHM15k at Payerne |

### Water-vapour & wavelength sensitivity
| Report | Summary |
|---|---|
| `attbsc_validation_technical.md` | Validation of calibrated attenuated backscatter via six independent strategies |
| `attbsc_wv_literature_review.md` | WV-absorption-correction literature review; anchors on Wiegner & Gasteiger (2015) WAPL |
| `payerne_cl61_wv_sensitivity.md` | WV correction brings CL61 to within +1.7 % of CHM15k (vs +19.3 % uncorrected) |
| `payerne_cl61_calibration_sensitivity.md` | CL61 Rayleigh vs liquid-cloud sensitivity at Payerne vs colocated CHM15k |
| `wv_fwhm_literature_review.md` | Laser-emission FWHM review (CL31/CL51/CL61): manufacturer vs measured (Qmini) |
| `wv_wavelength_sensitivity.md` | Sensitivity of the WV correction to laser wavelength (910.55 vs 910.74 nm) and FWHM |

Also in `doc/`: **`WATER_VAPOR_AUDIT.md`** — code audit of the Python water-vapour
two-way-transmission implementation vs MATLAB and ACTRIS-Cloudnet `atmoslib`; verdict:
the Python port is correct (and slightly more accurate; one latent bug fixed in Python).
