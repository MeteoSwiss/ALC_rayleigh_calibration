# E-PROFILE ALC calibration — technical reports

*Index. Consolidated 2026-07-10 (M. Hervo, MeteoSwiss). The ~50 working notes that had
accumulated here were merged into the **10 thematic reports** below — duplication removed,
superseded results demoted to labelled history, every load-bearing claim spot-checked against the
current code. The original notes were removed (recoverable from git history); their content lives on
in the consolidated reports, whose provenance lines list exactly which notes each one absorbed.*

## The reports

| # | Report | Scope |
|---|--------|-------|
| 1 | [Rayleigh (molecular) calibration](01_rayleigh_calibration.md) | Molecular-window detection; the pluggable method set and the `eprof_v2` (C8) default; robust/ensemble uncertainty (history); multi-site & full-archive method comparison, ranking and precision; v2 optimization; network diagnosis; night-to-night variability. |
| 2 | [Liquid-cloud (O'Connor) calibration](02_cloud_calibration.md) | The liquid-water-cloud method; the fixed 100–2400 m integration gate (CBH-robust); gate-config sweep; network yield; coupling to multiple scattering and the read-once cloud path. |
| 3 | [Water-vapour correction](03_water_vapour_correction.md) | The ~910 nm WV absorption correction (mandatory; CAMS L137): literature basis, per-type laser λ₀ + FWHM, wavelength-config sensitivity, CL61 sensitivity, and CAMS spatial/temporal resolution sensitivity. |
| 4 | [Multiple scattering](04_multiple_scattering.md) | The η(cloud-base) correction of the cloud calibration — the current PVC (Hogan 2006) tables at a_G = 5.5 µm, one per instrument type, and the derivation that replaced the legacy 8 µm ladder. |
| 5 | [Attenuated-backscatter validation](05_attbsc_validation.md) | Methodology & conventions; the uniform-L1 paper validation (benchmark + recent stations); L1-vs-L2; L1 validation on CSCS; noise-filter sensitivity. |
| 6 | [OmB — Observation-minus-Background](06_omb.md) | Per-station OmB vs the CAMS aerosol forecast (Payerne spot-check); operational L2 constant vs Kalman best-estimate `C_L`. |
| 7 | [CL61 calibration deep-dive](07_cl61_calibration.md) | Why `C_L`(Rayleigh) once differed from `C_L`(cloud) for CL61 and how the per-type η fix reconciled them; network verification; sensitivity; cloud/fog detection; 910↔1064 nm wavelength conversion. |
| 8 | [Overlap, near-range tilt & electronic offset](08_overlap_nearrange_offset.md) | Hood-free overlap reconstruction from clear-sky noise; the resolved Payerne CL61–CHM15k near-range tilt; single-site & network electronic-offset/dark characterisation; ambient noise & detection thresholds. |
| 9 | [Calibration stability, conventions & monitoring](09_calibration_stability_monitoring.md) | The Wiegner `C_L = RCS/β_att` convention (canonical); stability drivers; short-term-variability diagnosis; per-calibration outlier rate; network v2(C8)-vs-v1.1; June-2026 changelog. |
| 10 | [Operations, deployment & pipeline architecture](10_operations_deployment.md) | The read-once / shared-30 s×10 m-grid pipeline (implemented); the daily flow; the CSCS OmB + sensitivity runbook; the EWC dashboard deployment; the ceiloclass integration plan. |

## Conventions and current defaults (as of 2026-07-10)

- **Calibration coefficient:** the Wiegner lidar constant **`C_L = RCS/β_att`** everywhere (report 9 §1).
- **Rayleigh molecular-window method:** repo default **`eprof_v2`, config C8**. `calipso` is retired.
  (The *deployed* E-PROFILE operational network still runs the older **E-PROF v1.0** — distinct from
  the repo default, and still carrying the historical Klett sign error.)
- **Multiple scattering:** PVC (Hogan 2006) η tables at **a_G = 5.5 µm**, one per instrument type.
- **Water vapour:** **mandatory** for 910 nm (CL31/CL51/CL61); source = **CAMS model levels (L137)**,
  monthly 0.4° with a 1° fallback; a 910 nm night without usable CAMS is flagged, never calibrated WV-free.
- **Cloud calibration:** liquid-water O'Connor method over a **fixed 100–2400 m** integration gate.
- **Pipeline:** read-once / share-many — each instrument-day loads L1 and CAMS once, coarsens to a
  shared **30 s × 10 m** working grid, and fans out to classification + Rayleigh + cloud + housekeeping
  + OmB + sensitivity.
- **Operational system:** `zueub434.meteoswiss.ch`, `/data/zue/E_PROFILE/ALC/Calibration/`.

## Figures

Figures are stored alongside the reports in the `figs_*/` (and a few named) sub-directories and are
embedded inline in each report. `figs_extracted/` holds PNGs that were previously inlined as base64
in the source notes. A small number of figures produced on the R&D branch (`figs_paper_validation/`)
and MATLAB-side (`figs_ceilo_ambient/`) were never committed to this repo; where a report needed one,
its caption is retained with a note and a pointer to the script that regenerates it.
