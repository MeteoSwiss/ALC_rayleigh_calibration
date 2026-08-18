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

Analysis & validation reports for the E-PROFILE ALC paper (M. Hervo, MeteoSwiss). **Consolidated
2026-07-10** from ~50 working notes into the **10 thematic reports** below (duplication removed,
superseded results demoted to labelled history, claims spot-checked against the code). See
[`reports/README.md`](reports/README.md) for the full index, conventions and current defaults.

| # | Report | Scope |
|---|--------|-------|
| 1 | [Rayleigh (molecular) calibration](reports/01_rayleigh_calibration.md) | Molecular-window methods; the `eprof_v2` (C8) default; method comparison / ranking / precision; v2 optimization; network diagnosis; night-to-night variability |
| 2 | [Liquid-cloud (O'Connor) calibration](reports/02_cloud_calibration.md) | Cloud method; the fixed 100-2400 m integration gate; gate-config sweep; network yield |
| 3 | [Water-vapour correction](reports/03_water_vapour_correction.md) | 910 nm WV correction: literature, per-type λ₀ / FWHM, wavelength- and CAMS-resolution sensitivity |
| 4 | [Multiple scattering](reports/04_multiple_scattering.md) | PVC (Hogan 2006) η tables at a_G = 5.5 µm, one per instrument type |
| 5 | [Attenuated-backscatter validation](reports/05_attbsc_validation.md) | Uniform-L1 paper validation; L1-vs-L2; CSCS L1; noise-filter sensitivity |
| 6 | [OmB — Observation-minus-Background](reports/06_omb.md) | Payerne OmB vs CAMS; operational L2 constant vs Kalman `C_L` |
| 7 | [CL61 calibration deep-dive](reports/07_cl61_calibration.md) | Rayleigh-vs-cloud reconciliation; network verification; sensitivity; detection; 910↔1064 nm conversion |
| 8 | [Overlap, near-range tilt & electronic offset](reports/08_overlap_nearrange_offset.md) | Overlap-from-noise; the resolved Payerne tilt; offset / dark characterisation; ambient noise |
| 9 | [Calibration stability, conventions & monitoring](reports/09_calibration_stability_monitoring.md) | The Wiegner `C_L` convention; stability drivers; short-term variability; outliers; v2-vs-v1.1; June-2026 changelog |
| 10 | [Operations, deployment & pipeline architecture](reports/10_operations_deployment.md) | The read-once pipeline; CSCS OmB+sensitivity runbook; EWC dashboard; ceiloclass plan |


Also in `doc/`: **`WATER_VAPOR_AUDIT.md`** — code audit of the Python water-vapour
two-way-transmission implementation vs MATLAB and ACTRIS-Cloudnet `atmoslib`; verdict:
the Python port is correct (and slightly more accurate; one latent bug fixed in Python).
