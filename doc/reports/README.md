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

## Campaign reports — altitude independence & dark noise (2026-08)

*Six linked studies (2026-08-14/16) triggered by the eprof_v2.2 availability work: "is the
calibration constant altitude-independent, and why do cloud and Rayleigh disagree on the CL61?"
Every headline figure in these reports passed an independent adversarial re-derivation; refuted
claims are labelled as such inside each report.*

| Report | One-line verdict |
|--------|------------------|
| [Altitude-independence audit (observations)](altitude_independence_audit.md) | Rayleigh `C_L` is height-independent on the windows the gates actually allow (3–5 % spread, covered by the published uncertainty) but falls −7.5 %/km below 3.75 km on clear nights; cloud `C` rose +6–14 %/km with CBH under the legacy η tables — PVC cures CL31, half-cures CL51, not CL61. Within-night vs across-night gradients have OPPOSITE signs (endogeneity trap). |
| [Forward-model study (theory, closed simulators)](altitude_forward_model_study.md) | The shipped estimator is exact (closure 0.005 %, 1e-5 %/km); the observed gradient is ~95 % atmosphere / ~5 % estimator; overlap and photon noise are quantitatively ruled out; `subtract_background=True` subtracts atmosphere in disguise — never enable it; the published uncertainty is blind to the backscatter term (haze nights weigh ~6× too much in the Kalman). |
| [Dark campaign, saturation, AERONET, radiosondes](dark_aeronet_sonde_audit.md) | The covered-telescope campaign: CHM15k carries a constant −17 % pedestal (no proven drift); CL61 carries a range-growing baseline explaining 66 % of its within-night gradient; CL61 cloud-CBH saturation ruled out; WV exonerated within-night (sonde swap); AERONET closes the "AOD ×36" question — same column, different scale height, S = 52 sr validated; CAMS-1° PWV −26 % at Payerne (grid-point orography). |
| [CL61 cloud-vs-Rayleigh origin](cl61_cloud_vs_rayleigh_origin.md) | The disagreement closes: dark (−10 %) + in-window aerosol (+3 %) + λ_mol (−0.3 %) inside the cloud method's own floor (S_c ± 4.3 %). η tables PROVEN correct (no static table can produce the CBH residual). New instrument findings: Payerne CL61 diode thermal regulation lost 2026-06-11; the FOV ± convention supports 0.56 mrad half-angle. |
| [Phase-4 network validation (CSCS run, 433 streams)](phase4_network_validation.md) | v2.2 validated for deployment: availability CHM15k 86.6 %, corpus gain ×1.49, superset 99.5 %, recovered-night offset +0.37 % (time-paired). CL61 CBH residual NOT cured by native PVC + 0.4° CAMS (+9.05 ± 1.49 %/km, 11/11 streams). Payerne dark correction +24 % / +24 nights/yr confirmed on an independent run; ~3 % of CHM15k streams (5) carry a strong dark; the far-field proxy does NOT generalise (null result). |
| [Dark from clear-sky nights — method & network scan (FR)](dark_clearsky_method.md) | Per-gate slope+intercept across nights (CAMS molecular + aerosol regressors) separates electronics from molecular EXCEPT for the molecular-shaped component (proven invisible by injection). Usable score: θ on the Payerne hood template, affine-calibrated per stream; validated at Payerne (A θ≈1, C θ_corr=0.96). Messina confirmed (θ≈+5), Montsec's FE gradient is NOT window-range dark; blind below 2 km. |
| [dC/dCBH heterogeneity across units (cloud calibration, FR)](cbh_slope_heterogeneity.md) | Per-unit CBH slopes are REAL (CL31 τ ≈ 2.5–3.2 %/km, split-half r = 0.56): CL31 network flatness is a cancellation between genuinely positive and negative units, CL51 is shifted (+4.5) AND heterogeneous, CL61 is homogeneous at +8.3 (type effect). No measurable L1 factor explains the intra-type spread (FOV/divergence untestable; firmware confounded with country); flat-slope units are NOT better calibrated (raw \|slope\|–scatter link is mechanical). Recommendation: no network dC/dCBH correction, no η retuning (static-table ceiling +8.5 < 9.05 %/km); expose the per-flux slope as a QC diagnostic. ⚠ CORRECTED 2026-08-16: the "CL61 = additive baseline" route is REFUTED for the CLOUD slope (hood-measured closure −0.01 %/km vs +9) — the baseline only explains the RAYLEIGH side; the cloud mechanism stays open (prime lead: spectral λ₀/WV). |
| [Network dark diagnostic, hardware timeline & change detection (FR)](dark_network_diagnostic.md) | The FE indicator retrieves all 4 known darks in the top 6 of 153 streams (8 in the strong-dark zone; NEW candidates: Bern, Twenthe; 9 transients). L1 attrs give 20 CHM15k module swaps (~7.6 %/stream/yr) + a 2026-06 firmware rollout confound. VERDICT: the dark diagnostic canNOT detect hardware changes on its own (recall ≤ 33 %, precision ~26 %, 2026 common-mode false alarms) — detect via L1 serial/module attrs, qualify via dark, restart the Kalman at each swap. θ per era PROVES the dark follows the optical module both ways (arrives at Payerne with TUB140016, leaves HKZA with TUB150043); 3 swaps change sensitivity 26–39 % with no dark. Watch: rolling FE (network-median-corrected) < −6 %/km on 2 consecutive windows + θ confirmation; hood priorities Messina, Bern (θ first), Twenthe, Payerne CL61. |
| [Hopkin heatmaps, trust metric, slope stability & correction test (FR)](cbh_hopkin_study.md) | Hopkin-fig-6-style C_L(CBH) heatmaps for 241 units (gallery off-repo): pooled CL31 passes Hopkin's flatness test (band means 98–100 %), CL51/CL61 tilt. Trust metric \|slope\| × IQR(CBH): 1/266 units > 10 % impact, 92 % < 5 % — the long-term mean is protected by the CBH climatology (seasonal wobble median 1–3 %), the individual scene is not (±8–13 % worst-case). Slopes are stable unit properties (ratio ~1.2) modulated by a common ±2–5 %/km seasonal cycle; 30 genuine drifters. Out-of-sample correction test (raw/shrunk/type-mean): CL31 no gain, CL51 modest −4 % scatter with the fixed type slope, CL61 gains masked-symptom only — correct NOTHING downstream; publish slope+impact+stability as dashboard QC; inspect CL51 11487_A (3 independent flags). |

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
