# Rayleigh (molecular) calibration — methods, optimization and network performance

*Consolidated 2026-07-10 (M. Hervo, MeteoSwiss E-PROFILE ALC). Merges: molecular_window_detection_methods_report.md, molecular_methods_longrun_report.md, method_comparison_multisite.md, v2_optimization_report.md, ranking_robust_longrun.md, precision_longrun.md, rayleigh_network_diagnosis_report.md, l1_2026_variability_report.md, ROBUST_CALIBRATION_README.md, IMPLEMENTATION_SUMMARY.md, INTEGRATION_COMPLETE.md, DOCUMENTATION_VERIFICATION.md.*

Rayleigh (molecular) calibration is the core absolute-calibration path for elastic ceilometers and
lidars in the E-PROFILE ALC network (CL31/CL51/CL61 Vaisala, CHM15k Lufft/Jenoptik, Mini-MPL). An
elastic instrument measures only *total* attenuated backscatter, so to obtain an absolute constant it
must find an altitude window it can treat as **aerosol-free** and force the range-corrected signal to
match the molecular (Rayleigh) profile computed from T/p. The quality of the whole calibration is
hostage to that window choice. This document describes the pluggable molecular-window detector, the
tuned repo default (`eprof_v2` / config **C8**), the ensemble-uncertainty method (historical MATLAB
scaffolding), the multi-site and full-archive method comparisons, the v2 gate optimization, the
network diagnosis, and the night-to-night variability across instrument types.

Calibration coefficients follow the single **Wiegner lidar constant** convention `C_L = RCS/β_att`
across code, outputs, figures and reports.

---

## Table of contents

1. [Method — molecular-window detection and the pluggable method set](#1-method--molecular-window-detection-and-the-pluggable-method-set)
2. [Robust / ensemble uncertainty (historical MATLAB scaffolding)](#2-robust--ensemble-uncertainty-historical-matlab-scaffolding)
3. [Multi-site & long-run method comparison](#3-multi-site--long-run-method-comparison)
4. [v2 optimization — why clear nights fail, and the tuned C8](#4-v2-optimization--why-clear-nights-fail-and-the-tuned-c8)
5. [Network diagnosis — problematic-station causes, L1 vs L2](#5-network-diagnosis--problematic-station-causes-l1-vs-l2)
6. [Night-to-night variability — L1 2026 (CHM15k / Mini-MPL / CL61)](#6-night-to-night-variability--l1-2026-chm15k--mini-mpl--cl61)

---

## 1. Method — molecular-window detection and the pluggable method set

### 1.1 Principle

The molecular detector runs a **grid search** over candidate windows (centre × half-length) on the
night-mean range-normalized signal `signal = RCS/r²`. For each candidate window it computes the
free-fit slope/intercept/R², the forced-intercept-0 slope/RMSE/R², the median signal/molecular ratio
(the `C_L` proxy), the in-window SNR, the Rayleigh-shape residual, the scattering ratio, and — when
the per-profile stack is supplied — the **temporal** variability of the window-mean ratio across the
night. Different **selectors** then apply different eligibility masks and selection rules to that one
grid. The lidar constant is `C_L = median(signal/p_mol)` in the chosen window, passed to the Klett
β_att step and the downstream proportionality QC.

All logic lives in `calibration/rayleigh/molecular_methods.py`:
`compute_window_grid(signal, p_mol, range, …, signal_stack=…)` runs the grid once and returns every
candidate's statistics; the selectors consume it; `select_molecular_window(method, …)` is the
dispatcher. `rayleigh_fit.py` wires it into `calibrate_rayleigh`, and `config.py` / `options.json`
expose `molecular_method`. The `min_window_*` / `max_window_*` gates in `options.json` tune the
`eprof_v1.2` path only; every other method uses its own tuned defaults.

### 1.2 The pluggable method set (E-PROF version keys)

Methods are named by their **E-PROF calibration version**. The current code registers
`METHODS = ("eprof_v1.1", "eprof_v1.2", "eprof_v0.25", "earlinet", "eprof_v2", "eprof_v2p", "bellini")`;
legacy keys (`main`/`improved`/`matlab`/`optimal`/`eprof_v10`) are still accepted as **input aliases**.

| key (options.json) | label | fit | eligibility (gate) | selection rule | basis |
|---|---|---|---|---|---|
| `eprof_v1.0` | E-PROF v1.0 (sign error) | free | none | min Σ\|b\| centre, then max R² | legacy `main` window **+ pre-fix Klett sign error**; historical operational baseline |
| `eprof_v1.1` | E-PROF v1.1 (sign cor) | free intercept | none (any finite fit) | min Σ\|intercept\| centre, then max R² | legacy `main` window, sign-corrected — **degenerate** (no R² floor / above-aerosol gate) |
| `eprof_v0.25` | E-PROF v0.25 (MATLAB) | **intercept forced = 0** | R²₀≥0.5, slope>0; centres ≤5 km | min Σ RMSE centre, then max R² | MATLAB `Auto_Calib_25` (Hervo & Poltera 2014) |
| `eprof_v1.2` | E-PROF v1.2 (improved) | free intercept | start ≥2 km, R²≥0.5, slope>0, \|b\|<a, slope≈median ratio | **max R²** among eligible | production fix; Mattis 2016 shape gate + above-aerosol |
| `earlinet` | EARLINET/SCC | forced 0 | Rayleigh-shape residual ≤10 %, SNR gate, scattering ratio ≤1.1 | **lowest** qualifying window | EARLINET SCC (Mattis 2016 "take the lowest"), Freudenthaler 2018 |
| `eprof_v2` | **E-PROF v2 (default)** | free | all physical gates **+ temporal-variability ≤ threshold** | max composite quality (R² + shape + purity + SNR + **temporal steadiness**) | "optimal" best-of; the temporal idea is novel here. Config **C8** (§4). |
| `bellini` | Bellini/ALICENET | free | 3–7 km, width 600–3000 m; residuals not autocorrelated (Breusch–Godfrey); slope>0, intercept≈0; border-residual sign<0; E_CL≤40 % | max **M_Ray = (adjR²+(1−\|b\|))/std(b)** | ALICENET (Bellini et al. 2024, AMT 17, 6119) |

**`calipso` (CALIOP "highest clean layer") is RETIRED and removed from the method set.** A ground-up
ALC has no stratospheric pure-Rayleigh reference, so normalizing to the *highest* clean layer chases
noise. Do not use it. *(Superseded 2026-06: earlier drafts listed `calipso` as a live seventh method;
in the long-run tables below it appears only as historical evidence for why it fails.)*

An additional **`eprof_v2p` ("purity")** selector exists in the current code as a Stage-0 pilot for a
molecular-window-selection refinement: it uses the **same gates** as `eprof_v2` but, among eligible
windows, **minimises aerosol** (curvature rel-error dominant, then scattering ratio) instead of
maximising R², so a tight but aerosol-tainted low window cannot win on R² alone. Its weights are
provisional (to be optimised on the network sweep); it is **not** the default.

**Klett sign error (status).** The historical Klett sign error is **fixed in the active codebases**.
Only the **deployed operational E-PROF v1.0** still carries it (see §1.4). Importantly, the sign error
corrupts the downstream **attenuated-backscatter / AOD product**, **not** the Rayleigh calibration
constant `C_L` (quantified in §6): v1.0 and the sign-corrected versions share essentially the same
`C_L` variability.

### 1.3 The `eprof_v2` ("optimal") temporal-variability aerosol rejection

All single-profile methods collapse the night to **one mean profile** and must guess whether a
smooth, linear-looking layer is molecular or aerosol. `eprof_v2` additionally uses the **full profile
time series**: molecular scattering is **steady in time** (it only tracks T/p), whereas **aerosol
advects and fluctuates**. For each candidate window it computes the temporal coefficient of variation
of the window-mean signal/molecular ratio across the night (averaging over the window's range bins
first suppresses photon noise, leaving mostly atmospheric variability). Temporally variable windows
are **rejected as aerosol** and steadier windows are rewarded. This catches aerosol that *looks*
linear in the mean (high R²) but betrays itself by fluctuating — information every single-profile
method discards.

**Time-resolved layer flagging.** Beyond the per-window temporal CV, `eprof_v2` flags individual
*time-altitude cells* contaminated by aerosol/cloud: per altitude it takes the temporal median and
MAD of signal/molecular and flags cells exceeding `median + 4·MAD` (an upper-tail outlier test, so
clean molecular noise is left intact and the cleaned mean stays unbiased — an earlier percentile
threshold biased it low and had to be replaced). It then fits on the *time-cleaned* mean, using the
clean part of an otherwise-contaminated night (aerosol only at the start, a cloud only at the end).

**Flowchart of `eprof_v2`** (the red box is the distinctive time-resolved flagging step):

*(Figure — flowchart of the optimal molecular-window detection method. Source figure
`molecular_methods/optimal_flowchart.png` not committed to this repo.)*

Step by step: (1) take the night's per-profile range-normalized signal `signal(t,z)` and the
molecular profile `p_mol(z)`; (2) **flag and remove** aerosol/cloud cells per altitude
(`signal/p_mol > median + 4·MAD`); (3) average the un-flagged cells into a **time-cleaned mean**
profile; (4) grid-search windows (centre × half-length) and fit `signal = a·p_mol + b` in each; (5)
keep only windows passing **all** eligibility gates — start above the boundary-layer aerosol, R²≥0.40,
slope>0, |b|<a, Rayleigh-shape residual ≤16 %, scattering ratio ≤1.15, adequate in-window SNR,
slope ≈ median ratio (rel-error ≤15 %), and temporal CV ≤0.8 (C8 gates, §4); (6) if none pass, emit
**no calibration** (flag −2) and let the Kalman skip the night; (7) otherwise choose the window
maximising the composite quality `Q = R² − w_ratio·|R−1| − w_resid·resid − w_snr·SNR − w_tvar·tCV −
w_rel·rel + w_npts·n`; (8) `C_L = median(signal/p_mol)` in that window.

### 1.4 Repo default vs deployed operational network

- **Repo default = `eprof_v2`, config C8** (retuned 2026-06-20; `config.py` `molecular_method="eprof_v2"`).
  All "v2 (C8) vs v1.1/v1.2" comparisons in this document reflect the current repo default.
- **Deployed E-PROFILE operational network still runs E-PROF v1.0** — the older path that still
  contains the historical Klett sign error. Do **not** conflate the two. The v1.0 result is retained
  here as a citable baseline (it isolates the sign error's effect to β_att, not `C_L`; §6).

### 1.5 What EarthCARE (ATLID) does — and why an elastic ALC still needs a window search

ATLID is a **355 nm High-Spectral-Resolution Lidar (HSRL)** with three channels (co-polar Mie,
co-polar Rayleigh/molecular, cross-polar). A Fabry-Pérot **high-spectral-resolution etalon** separates
the spectrally broad molecular return (→ Rayleigh channel) from the narrow particulate peak (→ Mie
channel). Because it **measures the molecular return directly at every range gate**, ATLID does **not**
need to search for an aerosol-free window — the molecular channel *is* the calibration reference
throughout the profile (Donovan et al. 2024; Wehr et al. 2023). It still uses a high-altitude
(~30–40 km) pure-Rayleigh band, but for a *different* purpose: to determine spectral **cross-talk
coefficients** (χ = Mie-into-Rayleigh, ε = Rayleigh-into-Mie; a height-dependent channel matrix) and
inter-channel gain — not to fix an absolute backscatter constant.

**What an elastic ceilometer can borrow:** the discipline of *defining* clean air (volcanic-aerosol
awareness, averaging to beat noise, treating the clean-air assumption as testable — the temporal flag
is one such test); **per-shot background by interpolating a pre-/post-echo estimate**; a proper
T/p-dependent molecular term; and surface/cloud echoes as auxiliary references (parallels the
liquid-cloud calibration). **What cannot transfer:** the direct continuous molecular reference (needs
spectral separation), the χ/ε cross-talk correction, and extinction without a lidar-ratio assumption
— exactly what a single-channel elastic system lacks, which is *why* the molecular-window search
remains necessary. *(Sources: Donovan et al. 2024 AMT 17 5301; Wehr et al. 2023 AMT 16 3581; Eisinger
et al. 2024 AMT 17 839; Irbah et al. 2023 AMT 16 3631.)*

### 1.6 Literature basis for the selection criteria

- **Mattis, D'Amico, Baars et al. 2016** (EARLINET SCC, AMT 9, 3009): the minimum-signal search is
  **degenerate** — it "would find a minimum also in the case that there are fewer particles … large
  errors." Remedy: SNR/std gate + Rayleigh-shape test + **take the lowest qualifying window**. → basis
  of `earlinet`, and the reason `main`/`eprof_v1.1` fails.
- **Freudenthaler et al. 2018** (EARLINET QA): Rayleigh-fit **relative residual ≤ 1 %** + deviation
  plot. → the shape-residual gate.
- **Wiegner & Geiß 2012** (AMT 5, 1953): scattering ratio **R ≤ 1.1** as an acceptance criterion and
  systematic-error term; slope matching. → the scattering-ratio gate.
- **Baars et al. 2016** (PollyNET, ACP 16, 5111): automated Rayleigh-shape test.
- **Bellini et al. 2024** (ALICENET, AMT 17, 6119, Suppl. S3): two-step E-PROFILE-based Rayleigh fit
  in 3–7 km with a **Breusch–Godfrey residual-autocorrelation test**, the **M_Ray** window metric, and
  an **E_CL** relative-uncertainty gate. → basis of `bellini`.

### 1.7 Detailed examples — profiles, selected windows, time-resolved flagging

Each method runs on the **identical prepared profile** per night. In the profile figures the brackets
are each method's selected window; in the pcolor figures the colour is the signal/molecular ratio over
the night (molecular = vertically and temporally uniform; aerosol/cloud = enhanced and variable), the
dashed lines are the selected window centres, and the **hatched cells are those `eprof_v2` flagged as
aerosol/cloud and excluded** from its fit. The bottom row shows the resulting `C_L` ± uncertainty per
method, so the agreement of the constants is read directly beneath the windows that produced them.

**Payerne CL61 (910 nm)** — 2026-03-12 (aerosol-laden), -16 and -28 (cleaner); -16 has aerosol only at
the *start* of the night and -28 a thin cloud near the *end*, both appearing as hatched flagged cells:

*(Figures — Payerne CL61 vertical profiles with each method's selected window bracketed, and the
signal/molecular ratio time-height with window centres and the optimal-flagged cells. Source figures
`molecular_methods/profiles_Payerne_CL61.png`, `pcolor_Payerne_CL61.png` not committed to this repo.)*

**Amsterdam CHM15k (1064 nm)** — high-SNR reference; the signal tracks the molecular line cleanly from
~2 km to 7 km on clear nights.
*(Figures `molecular_methods/profiles_Amsterdam_CHM15k.png`, `pcolor_Amsterdam_CHM15k.png` not
committed to this repo.)*

**EDT CL61 (910 nm, Edmonton 53.5°N)** — high-latitude site (also exercises the darkness-adaptive
night window and the WV correction).
*(Figures `molecular_methods/profiles_EDT_CL61.png`, `pcolor_EDT_CL61.png` not committed to this
repo.)*

What the examples show: (1) on an **aerosol night** no method calibrates — the gated methods select no
eligible window; the ungated legacy path picks a degenerate high window (R²≈0.06, negative scattering
ratio, negative constant) that the pipeline's proportionality QC then rejects. (2) On **clear nights**
the gated methods agree on the constant to a few %; `earlinet` favours the lowest qualifying window,
`eprof_v2`/`eprof_v1.2`/`eprof_v0.25` land mid-profile. (3) `eprof_v2`'s time-resolved flagging
excises the start-of-night aerosol, the end-of-night thin cloud, and the persistent boundary-layer
aerosol, keeping the clean part of a partly-contaminated night usable.

---

## 2. Robust / ensemble uncertainty (historical MATLAB scaffolding)

> **Status — historical (MATLAB-era, Feb 2026).** The multi-window × multi-lidar-ratio GUM ensemble
> described below was implemented as a set of MATLAB `.m` functions
> (`findOptimalMolecularWindowEnsemble.m`, `calcLidarConstantRobust.m`, `validateEnsembleResults.m`,
> `calibrateRayleighWithDataRobust.m`) in Feb 2026. **The current Python pipeline does NOT use it.**
> Verified against `calibration/rayleigh/calibration.py`: the operational path selects a **single**
> molecular window and derives the uncertainty as a **robust 2σ** = 2 × max(half-IQR, MAD) of the
> in-window signal/molecular spread — not a window×LR ensemble budget. The ensemble method is kept
> here for its scientific rationale and as a design reference for a future Python port; treat the
> `.m` filenames, `use_robust_calibration` JSON flag, and per-instrument physical-CL bounds as
> MATLAB-era, not current behaviour.

The ensemble approach aimed at a scientifically rigorous, GUM-compliant uncertainty by (1) using
**multiple molecular windows** instead of one "best" window, (2) testing **multiple aerosol lidar
ratios** to quantify systematic LR uncertainty, and (3) combining Type A (statistical) and Type B
(systematic) components with a k=2 expanded uncertainty.

**Multi-window ensemble.** Select the top N non-overlapping windows (default N=5–7,
`min_window_separation` ≈ 300 m, `min_window_thickness` ≥ 200 m); average over independent regions to
reduce random uncertainty and detect atmospheric inhomogeneity. (Ref: Leblanc et al. 2016, NDACC
standards.)

**Lidar-ratio sensitivity.** Instead of a fixed LR (52 sr), test an ensemble — either an explicit
range (e.g. `[30,35,…,70]`) or an uncertainty-based band (`LRaer=52` ± `lidar_ratio_uncertainty=10`
→ `[32,42,52,62,72]`). Aerosol LR spans 20–30 sr (marine), 40–60 sr (continental/dust), 50–80 sr
(urban/pollution); the extinction correction is LR-dependent and this term was often the **dominant,
previously-unaccounted** uncertainty source. (Ref: Wandinger et al. 2016.)

**Uncertainty budget (GUM).**

```
u_random     = sqrt(u_statistical² + u_windows² + u_fit²)
u_systematic = sqrt(u_lidarRatio² + u_atmosphere² + u_overlap² + u_background²)
u_combined   = sqrt(u_random² + u_systematic²)
u_expanded   = 2 × u_combined          (k=2, ~95 % confidence)
```

Defaults: atmosphere (T/P model) 3 %, overlap 0 % (above full overlap), background 2 %. (Ref:
ISO/IEC Guide 98-3:2008.)

**Multi-criteria window quality score** (used by the ensemble window finder):

| criterion | description | threshold |
|---|---|---|
| R² fit | goodness of linear regression | > 0.95 |
| physical validity | \|intercept\| << \|slope\| | < 0.5 |
| statistical stability | coefficient of variation | < 15 % |
| precision | relative fit uncertainty | < 15 % |
| thickness | window size | > 200 m |
| smoothness | aerosol screening | < 10 % roughness |
| relative error | method agreement | < 15 % |

`Q = (R²)^0.25 × (intercept)^0.20 × (stability)^0.15 × (precision)^0.15 × (thickness)^0.10 ×
(smoothness)^0.05 × (rel_error)^0.10`

**Quality validation** ran 9 automated tests (ensemble CV, inter-window spread, LR sensitivity,
modified-Z outliers, u/CL ratio, physical-range, ensemble size, window R², relative uncertainty) and
graded EXCELLENT / GOOD / ACCEPTABLE / MARGINAL / FAILED. Instrument-specific physical-CL bounds were
used (CHM15k/CHM8k 1e10–1e13; Vaisala CL61/CL31/CL51 0.1–100; Mini-MPL 1e4–1e7; generic 1e-4–1e4).

Standard vs robust, as designed:

| aspect | standard | robust (MATLAB, historical) |
|---|---|---|
| windows | 1 (best) | 5–10 (ensemble) |
| lidar ratios | 1 (fixed) | 5–10 (range) |
| uncertainty components | 2 (std + fit) | 7 (full budget) |
| quality checks | basic | 9 comprehensive tests |
| uncertainty estimation | ~2σ approximation | GUM-compliant |
| typical rel. uncertainty | 10–15 % | 5–10 % (more realistic) |
| computation time | 1× | 3–5× |

*References (ensemble method):* Wiegner & Geiß 2012 (AMT 5, 1953); ISO/IEC Guide 98-3:2008 (GUM);
Leblanc et al. 2016 (AMT 9, 4029); Freudenthaler et al. 2018 (AMT 11, 4723); Wandinger et al. 2016
(AMT 9, 1001).

---

## 3. Multi-site & long-run method comparison

Three complementary evaluations: a small cross-site spot-check (Payerne / Amsterdam / EDT), the
full-archive long-run over 14 instruments, and drift-insensitive precision metrics. **All three favour
`eprof_v2` (C8) / `earlinet` for precision, with the caveat that "best" flips with the metric.**

### 3.1 Cross-site spot-check (Payerne, Amsterdam, EDT)

*Regenerated 2026-06-21 with the current method set (`eprof_v2` optimized C8; `calipso` retired). L2
monthly data; 35 nights (every 2nd day, Feb+Mar 2026) per instrument. `n_cal` = nights calibrating
through the full pipeline (rel_error ≤ 15 %); `CL_CV` = night-to-night `C_L` scatter.*

> **Small-sample caveat.** Only 35 nights are sampled and most are cloudy, so `n_cal` is small (0–8
> per instrument) and `CL_CV` is noisy (e.g. v1.1's 71.8 % mean is driven by a couple of EDT/Payerne
> outliers). Treat this as an illustrative spot-check; the statistically robust evidence is the
> full-network comparison (164 streams) and the long-run ranking below, both favouring **v2 (C8)** for
> precision and yield.

![Usable-night fraction (left) and night-to-night lidar-constant CV (right) per method across the key stations.](figs_extracted/method_comparison_multisite_01.png)

inst | type | method | nights | n_cal | CL_CV% | med_R2 | med_tcv | med_rel%
---|---|---|---|---|---|---|---|---
Payerne_CHM15k | CHM15k | eprof_v1.1 | 11 | 8 | 77.6 | 0.654 | 1.14 | 3.5
Payerne_CHM15k | CHM15k | eprof_v1.2 | 11 | 4 | 10.1 | 0.963 | 0.54 | 11.0
Payerne_CHM15k | CHM15k | eprof_v0.25 | 11 | 3 | 8.1 | 0.894 | 0.80 | 3.1
Payerne_CHM15k | CHM15k | earlinet | 11 | 2 | - | 0.904 | 0.33 | 11.3
Payerne_CHM15k | CHM15k | eprof_v2 | 11 | 5 | 10.4 | 0.957 | 0.33 | 10.9
Payerne_CHM15k | CHM15k | bellini | 11 | 1 | - | 0.784 | 0.68 | 2.1
Payerne_CL31 | CL31 | eprof_v1.1 | 13 | 3 | 25.5 | 0.391 | 0.90 | 7.5
Payerne_CL31 | CL31 | eprof_v1.2 | 13 | 0 | - | - | - | -
Payerne_CL31 | CL31 | eprof_v0.25 | 13 | 1 | - | 0.545 | 7.99 | 10.3
Payerne_CL31 | CL31 | earlinet | 13 | 0 | - | - | - | -
Payerne_CL31 | CL31 | eprof_v2 | 13 | 0 | - | - | - | -
Payerne_CL31 | CL31 | bellini | 13 | 0 | - | - | - | -
Payerne_CL61 | CL61 | eprof_v1.1 | 11 | 6 | 6.8 | 0.981 | 0.24 | 8.2
Payerne_CL61 | CL61 | eprof_v1.2 | 11 | 3 | 6.7 | 0.991 | 0.21 | 4.1
Payerne_CL61 | CL61 | eprof_v0.25 | 11 | 4 | 7.2 | 0.981 | 0.29 | 8.3
Payerne_CL61 | CL61 | earlinet | 11 | 3 | 9.0 | 0.691 | 0.23 | 5.6
Payerne_CL61 | CL61 | eprof_v2 | 11 | 4 | 9.3 | 0.956 | 0.30 | 8.9
Payerne_CL61 | CL61 | bellini | 11 | 5 | 12.9 | 0.756 | 0.42 | 9.3
Amsterdam_CHM15k | CHM15k | eprof_v1.1 | 4 | 3 | 18.0 | 0.932 | 0.31 | 4.8
Amsterdam_CHM15k | CHM15k | eprof_v1.2 | 4 | 3 | 8.2 | 0.987 | 0.08 | 4.5
Amsterdam_CHM15k | CHM15k | eprof_v0.25 | 4 | 3 | 12.9 | 0.963 | 0.28 | 5.2
Amsterdam_CHM15k | CHM15k | earlinet | 4 | 3 | 3.5 | 0.962 | 0.08 | 3.6
Amsterdam_CHM15k | CHM15k | eprof_v2 | 4 | 4 | 10.7 | 0.978 | 0.14 | 7.2
Amsterdam_CHM15k | CHM15k | bellini | 4 | 3 | 13.7 | 0.949 | 0.23 | 3.4
EDT_CL51 | CL51 | eprof_v1.1 | 6 | 3 | 231.3 | 0.983 | 0.79 | 5.5
EDT_CL51 | CL51 | eprof_v1.2 | 6 | 3 | 22.9 | 0.899 | 0.17 | 12.8
EDT_CL51 | CL51 | eprof_v0.25 | 6 | 3 | 16.5 | 0.979 | 0.33 | 1.6
EDT_CL51 | CL51 | earlinet | 6 | 0 | - | - | - | -
EDT_CL51 | CL51 | eprof_v2 | 6 | 0 | - | - | - | -
EDT_CL51 | CL51 | bellini | 6 | 0 | - | - | - | -
EDT_CL61 | CL61 | eprof_v1.1 | 5 | 2 | - | 0.888 | 0.76 | 5.8
EDT_CL61 | CL61 | eprof_v1.2 | 5 | 1 | - | 0.985 | 0.08 | 10.6
EDT_CL61 | CL61 | eprof_v0.25 | 5 | 1 | - | 0.938 | 0.53 | 0.2
EDT_CL61 | CL61 | earlinet | 5 | 0 | - | - | - | -
EDT_CL61 | CL61 | eprof_v2 | 5 | 0 | - | - | - | -
EDT_CL61 | CL61 | bellini | 5 | 1 | - | 0.842 | 0.76 | 14.8

![Per-method nightly lidar-constant time series across the sampled nights (one panel per instrument).](figs_extracted/method_comparison_multisite_02.png)

**Ranking (mean over instruments):**

method | calibrated-fraction | mean CL_CV% | mean temporal_cv | score
---|---|---|---|---
eprof_v1.2 | 0.35 | 12.0 | 0.22 | 0.286
eprof_v2 | 0.30 | 10.1 | 0.26 | 0.241
earlinet | 0.20 | 6.2 | 0.21 | 0.165
bellini | 0.25 | 13.3 | 0.52 | 0.164
eprof_v0.25 | 0.36 | 11.2 | 1.70 | 0.133
eprof_v1.1 | 0.53 | 71.8 | 0.69 | 0.087

On this small multi-site sample `eprof_v1.2` scores highest (best usable-night fraction at competitive
stability), with **`eprof_v2` (C8) a close second** (0.241 vs 0.286) that wins on the much larger
network and long-run samples. `earlinet` is the most *stable* where it calibrates (mean CL_CV 6.2 %)
but lowest-yield (0.20). v1.1's apparent high yield comes with a runaway CL_CV (71.8 %) from a few
unstable nights — the precision penalty the gated v2/v1.2 methods were designed to remove.

Representative detail figures:

![Payerne CL61 — per-method molecular windows on representative nights, and the signal/molecular time-height with selected centres.](figs_extracted/method_comparison_multisite_03.png)
![Amsterdam CHM15k — per-method molecular windows.](figs_extracted/method_comparison_multisite_04.png)
![EDT CL61 — per-method molecular windows.](figs_extracted/method_comparison_multisite_05.png)

### 3.2 Full-archive long-run (14 instruments)

The methods were run over the **entire E-PROFILE L2-monthly archive** for **14 instruments** — 10
CHM15k (Payerne, Lindenberg, Aosta, Palaiseau, Granada, Magurele, Bergen, Oslo, Hamburg,
Hohenpeissenberg) + 4 Mini-MPL (Brest, Toulouse, Corsica, SIRTA) — sampling 5 nights/month over each
site's full period (~80–113 months; ~400–520 fit-nights/site). This stresses every method on years of
varied atmospheres. *(This archive predates the CL61 rollout, so CL61 is covered separately in §6.
Method names in the long-run tables retain the pre-rename keys `main`/`improved`/`matlab`/`optimal`/
`calipso`/`eprof_v10`, which map to the current keys via the aliases in §1.2; `calipso` figures are
historical evidence for its retirement.)*

Over hundreds of nights the **classic std/mean CV explodes (100–1500 %) for the higher-yield methods**
— each admits a few nights with a spurious window and an extreme constant. The **robust (MAD-based) CV
collapses to 31–42 % for all methods** (the inflation was rare outliers the Kalman would smooth). The
two views together give the real ranking:

| method (current key) | calibrated-frac | robust CV % | std CV % (outlier-sensitive) | med R² | med temporal_cv |
|---|---|---|---|---|---|
| `eprof_v2` (optimal) | 0.45 | **31** | 128 | **0.97** | 0.15 |
| `earlinet` | 0.45 | **32** | 41 | 0.90 | **0.13** |
| `bellini` | 0.42 | 33 | 43 | 0.87 | 0.21 |
| `eprof_v0.25` (matlab) | 0.52 | 35 | 380 | 0.94 | 0.22 |
| `eprof_v1.2` (improved) | 0.46 | 38 | 293 | 0.97 | 0.16 |
| `eprof_v1.1` (main) | 0.57 | 37 | 445 | 0.91 | 0.25 |
| `calipso` *(retired)* | 0.79 | 42 | 385 | 0.65 | 0.27 |

*(Figures — usable nights + robust night-to-night CV per method (full archive, 14 sites); and the
calibration-constant time series per method per site. Source figures
`molecular_methods_longrun/summary_robust_longrun.png`, `timeseries_longrun.png` not committed to this
repo.)*

**Robust-CV ranking (mean over instruments, yield-weighted score):**

method | calibrated-fraction | mean rob_CV% | mean temporal_cv | score
---|---|---|---|---
calipso *(retired)* | 0.79 | 42.5 | 0.27 | 0.293
optimal (`eprof_v2`) | 0.45 | 31.3 | 0.15 | 0.253
earlinet | 0.45 | 32.1 | 0.13 | 0.248
main (`eprof_v1.1`) | 0.57 | 37.4 | 0.25 | 0.242
matlab (`eprof_v0.25`) | 0.52 | 35.2 | 0.22 | 0.242
bellini | 0.42 | 32.8 | 0.21 | 0.212
improved (`eprof_v1.2`) | 0.46 | 37.7 | 0.16 | 0.211
eprof_v10 (`eprof_v1.0`) | 0.46 | 33.2 | - | 0.140

Note that a *yield-weighted* robust score ranks `calipso` first purely on raw night count — a warning
that yield-weighting rewards the degenerate high-yield behaviour. On **quality + stability** the winner
is `eprof_v2` (or `earlinet`). This is exactly why the retired `calipso` must not be used: high count,
lowest R² (0.65, noisy high windows), highest robust scatter.

Per-instrument findings: the **CHM15k seasonal cycle** (internal-temperature / laser-ageing) is
visible in the stable methods' time series — exactly what a calibration monitor should track; the
unstable methods bury it under scatter. **Mini-MPL (532 nm) is the easy case** — all methods reach
R²≈0.99 and low temporal CV; the clean high-SNR molecular column means even `main`/`calipso` behave,
though `eprof_v2` still has the lowest CV. Full per-instrument numbers are in §3.4.

### 3.3 Precision — the right way to evaluate (drift-insensitive)

**CV is the wrong metric**: it mixes measurement *precision* with the real *seasonal + laser-ageing
drift* a calibration is supposed to track, so it penalises a precise instrument for having a seasonal
cycle. We use metrics that remove slow drift:

- **σ_SD** — *successive-difference* precision (von Neumann): robust scatter of |CLᵢ₊₁−CLᵢ| between
  time-ordered consecutive calibrations ÷ √2. Slow drift cancels in the difference ⇒ short-term noise
  only. **Headline precision metric.**
- **σ_detrend** — robust scatter of CL minus its ~2-month rolling median (seasonal + trend removed).
- **σ_within-month** — robust scatter of CL pooled within each calendar month (season ≈ const).
- **σ_night** — average single-night spread (in-window std of signal/molecular ÷ CL).
- **valid %** — fraction of sampled nights yielding a valid calibration.

**Per method (mean over the 14 instruments):**

| method (current key) | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_month % | CV % (ref) |
|---|---|---|---|---|---|---|
| `eprof_v2` (optimal) | 45 | 8.6 | **13.6** | 12.5 | 8.4 | 128 |
| `earlinet` | 45 | 4.2 | **15.9** | 13.6 | 10.0 | 41 |
| `eprof_v1.0` (eprof_v10) | 46 | 4.6 | **16.5** | 13.4 | 10.3 | 412 |
| `bellini` | 42 | 10.7 | **18.3** | 15.4 | 11.0 | 43 |
| `eprof_v0.25` (matlab) | 52 | 16.7 | **18.5** | 15.3 | 12.2 | 380 |
| `eprof_v1.2` (improved) | 46 | 9.4 | **19.2** | 16.1 | 12.3 | 293 |
| `eprof_v1.1` (main) | 57 | 20.0 | **20.4** | 16.0 | 12.1 | 444 |
| `calipso` *(retired)* | 79 | 14.9 | **24.9** | 19.9 | 16.4 | 385 |

*(Figure — yield and drift-insensitive precision metrics per method. Source figure
`molecular_methods_longrun/precision_longrun.png` not committed to this repo.)*

**This resolves the earlier metric ambiguity.** The decisive facts:

1. **CV (128–444 %) ≫ σ_SD (14–25 %) for every method** — *most of the CV was real drift, not noise.*
   The genuine per-night precision is ~14–25 %, far better than CV implied. (The two "low-CV" methods —
   earlinet 41 %, bellini 43 % — only looked good because they reject the drift-revealing nights;
   their σ_SD is mid-pack.)
2. **`eprof_v2` is the most precise method** by all three drift-insensitive measures (σ_SD 13.6 %,
   σ_detrend 12.5 %, σ_within-month 8.4 %) and has a low within-night spread (8.6 %).
3. **`earlinet` is second on σ_SD with the tightest within-night spread (4.2 %)** — its narrow,
   low-altitude windows are internally very consistent.
4. **`calipso` is the least precise (σ_SD 24.9 %)** despite the highest yield (79 %) — quantity at the
   cost of precision (its noisy high windows). Reinforces its retirement.
5. **The drift CV reveals is physical and useful** — the seasonal/laser-ageing cycle is exactly what
   the Kalman should track; a precision metric must not conflate it with noise.

### 3.4 Full per-instrument long-run tables

**Robust aggregation** (`rob_CV` = 1.4826·MAD/median; `std_CV` = classic std/mean, reference):

inst | type | method | nights | n_cal | rob_CV% | std_CV% | med_R2 | med_tcv
---|---|---|---|---|---|---|---|---
Aosta_CHM15k | CHM15k | eprof_v10 | 29 | 13 | 16.1 | 15 | - | -
Aosta_CHM15k | CHM15k | main | 29 | 17 | 14.8 | 13 | 0.894 | 0.44
Aosta_CHM15k | CHM15k | improved | 29 | 13 | 14.7 | 11 | 0.947 | 0.29
Aosta_CHM15k | CHM15k | optimal | 29 | 10 | 10.2 | 9 | 0.987 | 0.23
Aosta_CHM15k | CHM15k | matlab | 29 | 13 | 13.7 | 10 | 0.947 | 0.36
Aosta_CHM15k | CHM15k | calipso | 29 | 23 | 21.3 | 327 | 0.686 | 0.41
Aosta_CHM15k | CHM15k | earlinet | 29 | 10 | 8.3 | 9 | 0.842 | 0.24
Aosta_CHM15k | CHM15k | bellini | 29 | 15 | 15.4 | 28 | 0.829 | 0.32
Bergen_CHM15k | CHM15k | eprof_v10 | 411 | 89 | 35.9 | 596 | - | -
Bergen_CHM15k | CHM15k | main | 411 | 139 | 48.6 | 592 | 0.885 | 0.39
Bergen_CHM15k | CHM15k | improved | 411 | 100 | 48.8 | 358 | 0.961 | 0.21
Bergen_CHM15k | CHM15k | optimal | 411 | 86 | 53.1 | 39 | 0.966 | 0.18
Bergen_CHM15k | CHM15k | matlab | 411 | 130 | 50.3 | 590 | 0.898 | 0.37
Bergen_CHM15k | CHM15k | calipso | 411 | 231 | 49.3 | 460 | 0.634 | 0.35
Bergen_CHM15k | CHM15k | earlinet | 411 | 74 | 53.6 | 39 | 0.913 | 0.14
Bergen_CHM15k | CHM15k | bellini | 411 | 65 | 35.0 | 38 | 0.756 | 0.30
Brest-MPL_MPL | Mini-MPL | eprof_v10 | 436 | 197 | 64.5 | 66 | - | -
Brest-MPL_MPL | Mini-MPL | main | 436 | 268 | 75.0 | 77 | 0.998 | 0.10
Brest-MPL_MPL | Mini-MPL | improved | 436 | 245 | 77.9 | 61 | 0.999 | 0.11
Brest-MPL_MPL | Mini-MPL | optimal | 436 | 195 | 63.3 | 46 | 0.999 | 0.06
Brest-MPL_MPL | Mini-MPL | matlab | 436 | 257 | 72.0 | 74 | 0.997 | 0.09
Brest-MPL_MPL | Mini-MPL | calipso | 436 | 308 | 90.3 | 71 | 0.995 | 0.24
Brest-MPL_MPL | Mini-MPL | earlinet | 436 | 246 | 68.0 | 56 | 0.990 | 0.09
Brest-MPL_MPL | Mini-MPL | bellini | 436 | 224 | 64.7 | 58 | 0.993 | 0.10
Corsica-MPL_MPL | Mini-MPL | eprof_v10 | 431 | 317 | 35.4 | 41 | - | -
Corsica-MPL_MPL | Mini-MPL | main | 431 | 353 | 37.7 | 47 | 0.999 | 0.05
Corsica-MPL_MPL | Mini-MPL | improved | 431 | 324 | 34.6 | 41 | 1.000 | 0.05
Corsica-MPL_MPL | Mini-MPL | optimal | 431 | 294 | 30.3 | 31 | 0.999 | 0.04
Corsica-MPL_MPL | Mini-MPL | matlab | 431 | 336 | 36.6 | 46 | 0.999 | 0.05
Corsica-MPL_MPL | Mini-MPL | calipso | 431 | 370 | 42.0 | 49 | 0.999 | 0.06
Corsica-MPL_MPL | Mini-MPL | earlinet | 431 | 328 | 33.4 | 42 | 0.994 | 0.05
Corsica-MPL_MPL | Mini-MPL | bellini | 431 | 261 | 34.1 | 40 | 0.999 | 0.05
Granada_CHM15k | CHM15k | eprof_v10 | 465 | 275 | 17.6 | 123 | - | -
Granada_CHM15k | CHM15k | main | 465 | 306 | 16.5 | 496 | 0.953 | 0.30
Granada_CHM15k | CHM15k | improved | 465 | 257 | 14.8 | 151 | 0.988 | 0.17
Granada_CHM15k | CHM15k | optimal | 465 | 263 | 13.1 | 87 | 0.980 | 0.19
Granada_CHM15k | CHM15k | matlab | 465 | 288 | 15.9 | 469 | 0.971 | 0.25
Granada_CHM15k | CHM15k | calipso | 465 | 437 | 19.7 | 189 | 0.730 | 0.37
Granada_CHM15k | CHM15k | earlinet | 465 | 253 | 13.2 | 87 | 0.917 | 0.14
Granada_CHM15k | CHM15k | bellini | 465 | 279 | 15.2 | 84 | 0.899 | 0.24
Hamburg_CHM15k | CHM15k | eprof_v10 | 495 | 168 | 22.9 | 486 | - | -
Hamburg_CHM15k | CHM15k | main | 495 | 216 | 24.4 | 492 | 0.935 | 0.26
Hamburg_CHM15k | CHM15k | improved | 495 | 175 | 21.7 | 294 | 0.976 | 0.15
Hamburg_CHM15k | CHM15k | optimal | 495 | 199 | 16.3 | 18 | 0.969 | 0.16
Hamburg_CHM15k | CHM15k | matlab | 495 | 200 | 22.8 | 500 | 0.955 | 0.21
Hamburg_CHM15k | CHM15k | calipso | 495 | 375 | 23.6 | 403 | 0.667 | 0.27
Hamburg_CHM15k | CHM15k | earlinet | 495 | 172 | 14.9 | 17 | 0.906 | 0.13
Hamburg_CHM15k | CHM15k | bellini | 495 | 166 | 24.7 | 33 | 0.862 | 0.24
Hohenpeiss_CHM15k | CHM15k | eprof_v10 | 482 | 193 | 28.5 | 572 | - | -
Hohenpeiss_CHM15k | CHM15k | main | 482 | 246 | 29.0 | 572 | 0.913 | 0.32
Hohenpeiss_CHM15k | CHM15k | improved | 482 | 178 | 32.3 | 121 | 0.979 | 0.18
Hohenpeiss_CHM15k | CHM15k | optimal | 482 | 198 | 26.5 | 54 | 0.977 | 0.17
Hohenpeiss_CHM15k | CHM15k | matlab | 482 | 218 | 27.9 | 613 | 0.947 | 0.28
Hohenpeiss_CHM15k | CHM15k | calipso | 482 | 388 | 35.5 | 549 | 0.672 | 0.33
Hohenpeiss_CHM15k | CHM15k | earlinet | 482 | 170 | 21.7 | 25 | 0.921 | 0.13
Hohenpeiss_CHM15k | CHM15k | bellini | 482 | 164 | 26.6 | 40 | 0.887 | 0.28
Lindenberg_CHM15k | CHM15k | eprof_v10 | 495 | 209 | 29.9 | 681 | - | -
Lindenberg_CHM15k | CHM15k | main | 495 | 258 | 31.7 | 648 | 0.901 | 0.33
Lindenberg_CHM15k | CHM15k | improved | 495 | 181 | 28.3 | 244 | 0.965 | 0.20
Lindenberg_CHM15k | CHM15k | optimal | 495 | 182 | 29.9 | 87 | 0.963 | 0.18
Lindenberg_CHM15k | CHM15k | matlab | 495 | 223 | 32.3 | 627 | 0.942 | 0.28
Lindenberg_CHM15k | CHM15k | calipso | 495 | 400 | 37.9 | 298 | 0.639 | 0.32
Lindenberg_CHM15k | CHM15k | earlinet | 495 | 158 | 30.0 | 30 | 0.903 | 0.15
Lindenberg_CHM15k | CHM15k | bellini | 495 | 168 | 32.6 | 37 | 0.815 | 0.28
Magurele_CHM15k | CHM15k | eprof_v10 | 368 | 172 | 22.8 | 590 | - | -
Magurele_CHM15k | CHM15k | main | 368 | 191 | 24.5 | 566 | 0.956 | 0.23
Magurele_CHM15k | CHM15k | improved | 368 | 151 | 28.9 | 137 | 0.985 | 0.15
Magurele_CHM15k | CHM15k | optimal | 368 | 204 | 21.7 | 47 | 0.970 | 0.16
Magurele_CHM15k | CHM15k | matlab | 368 | 178 | 24.9 | 590 | 0.968 | 0.18
Magurele_CHM15k | CHM15k | calipso | 368 | 319 | 28.7 | 375 | 0.769 | 0.25
Magurele_CHM15k | CHM15k | earlinet | 368 | 188 | 22.1 | 50 | 0.904 | 0.13
Magurele_CHM15k | CHM15k | bellini | 368 | 169 | 22.2 | 24 | 0.893 | 0.20
Oslo_CHM15k | CHM15k | eprof_v10 | 400 | 150 | 28.6 | 623 | - | -
Oslo_CHM15k | CHM15k | main | 400 | 197 | 40.5 | 563 | 0.885 | 0.32
Oslo_CHM15k | CHM15k | improved | 400 | 141 | 34.7 | 809 | 0.952 | 0.21
Oslo_CHM15k | CHM15k | optimal | 400 | 164 | 21.3 | 43 | 0.963 | 0.19
Oslo_CHM15k | CHM15k | matlab | 400 | 169 | 27.9 | 541 | 0.944 | 0.28
Oslo_CHM15k | CHM15k | calipso | 400 | 279 | 40.0 | 387 | 0.647 | 0.30
Oslo_CHM15k | CHM15k | earlinet | 400 | 146 | 18.9 | 43 | 0.878 | 0.17
Oslo_CHM15k | CHM15k | bellini | 400 | 113 | 14.9 | 35 | 0.807 | 0.28
Palaiseau_CHM15k | CHM15k | eprof_v10 | 524 | 232 | 23.0 | 605 | - | -
Palaiseau_CHM15k | CHM15k | main | 524 | 292 | 25.8 | 574 | 0.896 | 0.33
Palaiseau_CHM15k | CHM15k | improved | 524 | 209 | 26.7 | 616 | 0.965 | 0.20
Palaiseau_CHM15k | CHM15k | optimal | 524 | 189 | 22.0 | 34 | 0.961 | 0.19
Palaiseau_CHM15k | CHM15k | matlab | 524 | 267 | 22.7 | 572 | 0.937 | 0.29
Palaiseau_CHM15k | CHM15k | calipso | 524 | 417 | 26.5 | 623 | 0.605 | 0.31
Palaiseau_CHM15k | CHM15k | earlinet | 524 | 167 | 22.0 | 32 | 0.895 | 0.16
Palaiseau_CHM15k | CHM15k | bellini | 524 | 155 | 23.8 | 32 | 0.830 | 0.27
Payerne_CHM15k | CHM15k | eprof_v10 | 518 | 163 | 43.6 | 1267 | - | -
Payerne_CHM15k | CHM15k | main | 518 | 215 | 47.6 | 1455 | 0.908 | 0.34
Payerne_CHM15k | CHM15k | improved | 518 | 135 | 56.5 | 1149 | 0.975 | 0.17
Payerne_CHM15k | CHM15k | optimal | 518 | 148 | 40.3 | 1207 | 0.967 | 0.18
Payerne_CHM15k | CHM15k | matlab | 518 | 184 | 44.2 | 574 | 0.940 | 0.24
Payerne_CHM15k | CHM15k | calipso | 518 | 389 | 62.5 | 1534 | 0.650 | 0.33
Payerne_CHM15k | CHM15k | earlinet | 518 | 133 | 41.2 | 38 | 0.887 | 0.15
Payerne_CHM15k | CHM15k | bellini | 518 | 129 | 42.3 | 54 | 0.888 | 0.20
SIRTA-MPL_MPL | Mini-MPL | eprof_v10 | 432 | 274 | 49.0 | 46 | - | -
SIRTA-MPL_MPL | Mini-MPL | main | 432 | 307 | 51.8 | 53 | 0.999 | 0.10
SIRTA-MPL_MPL | Mini-MPL | improved | 432 | 300 | 53.0 | 48 | 0.999 | 0.10
SIRTA-MPL_MPL | Mini-MPL | optimal | 432 | 281 | 45.7 | 40 | 0.999 | 0.08
SIRTA-MPL_MPL | Mini-MPL | matlab | 432 | 300 | 49.8 | 51 | 0.998 | 0.10
SIRTA-MPL_MPL | Mini-MPL | calipso | 432 | 372 | 61.4 | 63 | 0.997 | 0.13
SIRTA-MPL_MPL | Mini-MPL | earlinet | 432 | 306 | 52.0 | 47 | 0.983 | 0.10
SIRTA-MPL_MPL | Mini-MPL | bellini | 432 | 260 | 53.8 | 49 | 0.997 | 0.11
Toulouse-MPL_MPL | Mini-MPL | eprof_v10 | 436 | 285 | 47.5 | 61 | - | -
Toulouse-MPL_MPL | Mini-MPL | main | 436 | 338 | 56.0 | 66 | 0.999 | 0.07
Toulouse-MPL_MPL | Mini-MPL | improved | 436 | 305 | 54.1 | 58 | 0.999 | 0.06
Toulouse-MPL_MPL | Mini-MPL | optimal | 436 | 290 | 44.1 | 49 | 0.999 | 0.05
Toulouse-MPL_MPL | Mini-MPL | matlab | 436 | 322 | 52.7 | 63 | 0.998 | 0.06
Toulouse-MPL_MPL | Mini-MPL | calipso | 436 | 369 | 56.2 | 67 | 0.996 | 0.08
Toulouse-MPL_MPL | Mini-MPL | earlinet | 436 | 318 | 50.3 | 55 | 0.990 | 0.06
Toulouse-MPL_MPL | Mini-MPL | bellini | 436 | 255 | 53.5 | 56 | 0.997 | 0.06

**Drift-insensitive precision, per instrument × method** (all scatters % of median CL):

inst | type | method | valid_% | sigma_night% | sigma_SD% | sigma_dt% | sigma_im% | CV%
---|---|---|---|---|---|---|---|---
Aosta_CHM15k | CHM15k | eprof_v10 | 45 | 4.8 | 10.6 | 10.0 | 10.2 | 15
Aosta_CHM15k | CHM15k | main | 59 | 28.3 | 14.9 | 9.9 | 4.9 | 13
Aosta_CHM15k | CHM15k | improved | 45 | 18.4 | 9.9 | 13.7 | 10.0 | 11
Aosta_CHM15k | CHM15k | optimal | 34 | 9.7 | 8.6 | 9.6 | 7.2 | 9
Aosta_CHM15k | CHM15k | matlab | 45 | 19.1 | 10.1 | 8.1 | 14.4 | 10
Aosta_CHM15k | CHM15k | calipso | 79 | 24.4 | 11.7 | 16.9 | 14.8 | 327
Aosta_CHM15k | CHM15k | earlinet | 34 | 5.7 | 12.0 | 5.2 | 4.9 | 9
Aosta_CHM15k | CHM15k | bellini | 52 | 13.5 | 14.4 | 6.0 | 4.3 | 28
Bergen_CHM15k | CHM15k | eprof_v10 | 22 | 6.1 | 14.3 | 11.9 | 8.7 | 596
Bergen_CHM15k | CHM15k | main | 33 | 34.5 | 19.8 | 13.8 | 9.0 | 592
Bergen_CHM15k | CHM15k | improved | 24 | 14.0 | 12.9 | 11.9 | 8.2 | 358
Bergen_CHM15k | CHM15k | optimal | 21 | 10.9 | 9.3 | 9.8 | 5.8 | 39
Bergen_CHM15k | CHM15k | matlab | 32 | 33.2 | 12.8 | 10.0 | 8.3 | 590
Bergen_CHM15k | CHM15k | calipso | 56 | 21.2 | 22.9 | 14.9 | 13.0 | 460
Bergen_CHM15k | CHM15k | earlinet | 18 | 5.7 | 11.8 | 8.7 | 5.3 | 39
Bergen_CHM15k | CHM15k | bellini | 16 | 16.9 | 17.8 | 12.9 | 6.6 | 38
Brest-MPL_MPL | Mini-MPL | eprof_v10 | 45 | 3.5 | 28.0 | 22.2 | 16.9 | 66
Brest-MPL_MPL | Mini-MPL | main | 61 | 2.9 | 46.9 | 37.9 | 27.7 | 77
Brest-MPL_MPL | Mini-MPL | improved | 56 | 2.6 | 41.9 | 38.3 | 26.7 | 61
Brest-MPL_MPL | Mini-MPL | optimal | 45 | 2.6 | 22.2 | 25.5 | 14.4 | 46
Brest-MPL_MPL | Mini-MPL | matlab | 59 | 3.0 | 38.5 | 33.4 | 24.1 | 74
Brest-MPL_MPL | Mini-MPL | calipso | 71 | 4.0 | 56.4 | 50.5 | 39.2 | 71
Brest-MPL_MPL | Mini-MPL | earlinet | 56 | 1.6 | 34.3 | 31.9 | 24.1 | 56
Brest-MPL_MPL | Mini-MPL | bellini | 51 | 2.2 | 39.8 | 34.2 | 26.2 | 58
Corsica-MPL_MPL | Mini-MPL | eprof_v10 | 74 | 1.2 | 17.8 | 14.5 | 10.4 | 41
Corsica-MPL_MPL | Mini-MPL | main | 81 | 1.3 | 23.0 | 16.9 | 13.1 | 47
Corsica-MPL_MPL | Mini-MPL | improved | 75 | 1.2 | 20.0 | 14.9 | 14.1 | 41
Corsica-MPL_MPL | Mini-MPL | optimal | 68 | 1.7 | 15.3 | 13.0 | 9.4 | 31
Corsica-MPL_MPL | Mini-MPL | matlab | 78 | 1.3 | 23.7 | 18.1 | 14.3 | 46
Corsica-MPL_MPL | Mini-MPL | calipso | 86 | 2.1 | 21.7 | 19.8 | 13.4 | 49
Corsica-MPL_MPL | Mini-MPL | earlinet | 76 | 1.2 | 18.7 | 14.8 | 11.9 | 42
Corsica-MPL_MPL | Mini-MPL | bellini | 61 | 0.7 | 18.9 | 15.2 | 11.8 | 40
Granada_CHM15k | CHM15k | eprof_v10 | 59 | 5.0 | 11.3 | 10.1 | 7.6 | 123
Granada_CHM15k | CHM15k | main | 66 | 20.0 | 12.3 | 10.2 | 8.5 | 496
Granada_CHM15k | CHM15k | improved | 55 | 9.9 | 13.9 | 10.9 | 9.2 | 151
Granada_CHM15k | CHM15k | optimal | 57 | 11.3 | 10.9 | 9.4 | 6.2 | 87
Granada_CHM15k | CHM15k | matlab | 62 | 15.9 | 12.4 | 10.6 | 8.7 | 469
Granada_CHM15k | CHM15k | calipso | 94 | 19.1 | 17.9 | 13.4 | 11.7 | 189
Granada_CHM15k | CHM15k | earlinet | 54 | 4.8 | 11.1 | 9.0 | 7.1 | 87
Granada_CHM15k | CHM15k | bellini | 60 | 12.4 | 12.6 | 11.1 | 8.0 | 84
Hamburg_CHM15k | CHM15k | eprof_v10 | 34 | 5.7 | 13.3 | 11.3 | 7.6 | 486
Hamburg_CHM15k | CHM15k | main | 44 | 24.1 | 15.0 | 11.8 | 8.4 | 492
Hamburg_CHM15k | CHM15k | improved | 35 | 12.4 | 18.7 | 12.8 | 10.5 | 294
Hamburg_CHM15k | CHM15k | optimal | 40 | 11.2 | 10.4 | 9.6 | 6.6 | 18
Hamburg_CHM15k | CHM15k | matlab | 40 | 20.6 | 15.0 | 12.4 | 8.4 | 500
Hamburg_CHM15k | CHM15k | calipso | 76 | 18.5 | 18.2 | 13.6 | 10.4 | 403
Hamburg_CHM15k | CHM15k | earlinet | 35 | 4.9 | 11.1 | 9.5 | 7.1 | 17
Hamburg_CHM15k | CHM15k | bellini | 34 | 13.2 | 16.3 | 12.9 | 10.2 | 33
Hohenpeiss_CHM15k | CHM15k | eprof_v10 | 40 | 5.1 | 15.3 | 12.6 | 8.7 | 572
Hohenpeiss_CHM15k | CHM15k | main | 51 | 25.4 | 13.5 | 12.2 | 8.1 | 572
Hohenpeiss_CHM15k | CHM15k | improved | 37 | 10.9 | 12.9 | 8.3 | 7.1 | 121
Hohenpeiss_CHM15k | CHM15k | optimal | 41 | 11.2 | 9.9 | 8.2 | 5.2 | 54
Hohenpeiss_CHM15k | CHM15k | matlab | 45 | 21.0 | 14.1 | 10.1 | 7.2 | 613
Hohenpeiss_CHM15k | CHM15k | calipso | 80 | 19.0 | 22.3 | 13.5 | 13.2 | 549
Hohenpeiss_CHM15k | CHM15k | earlinet | 35 | 4.5 | 9.2 | 8.1 | 5.1 | 25
Hohenpeiss_CHM15k | CHM15k | bellini | 34 | 14.4 | 10.8 | 8.7 | 7.3 | 40
Lindenberg_CHM15k | CHM15k | eprof_v10 | 42 | 5.9 | 11.2 | 9.0 | 7.6 | 681
Lindenberg_CHM15k | CHM15k | main | 52 | 28.4 | 14.2 | 10.8 | 8.9 | 648
Lindenberg_CHM15k | CHM15k | improved | 37 | 12.8 | 15.5 | 12.1 | 8.6 | 244
Lindenberg_CHM15k | CHM15k | optimal | 37 | 11.3 | 11.0 | 8.1 | 5.0 | 87
Lindenberg_CHM15k | CHM15k | matlab | 45 | 25.9 | 14.3 | 11.8 | 9.6 | 627
Lindenberg_CHM15k | CHM15k | calipso | 81 | 20.6 | 22.6 | 14.8 | 13.6 | 298
Lindenberg_CHM15k | CHM15k | earlinet | 32 | 5.6 | 10.8 | 8.7 | 5.6 | 30
Lindenberg_CHM15k | CHM15k | bellini | 34 | 16.0 | 14.7 | 13.3 | 9.3 | 37
Magurele_CHM15k | CHM15k | eprof_v10 | 47 | 5.1 | 14.7 | 10.3 | 9.1 | 590
Magurele_CHM15k | CHM15k | main | 52 | 16.4 | 16.5 | 10.3 | 9.2 | 566
Magurele_CHM15k | CHM15k | improved | 41 | 8.3 | 19.9 | 14.9 | 10.5 | 137
Magurele_CHM15k | CHM15k | optimal | 55 | 10.6 | 11.9 | 10.6 | 7.3 | 47
Magurele_CHM15k | CHM15k | matlab | 48 | 13.1 | 15.6 | 11.8 | 9.6 | 590
Magurele_CHM15k | CHM15k | calipso | 87 | 14.7 | 18.8 | 14.3 | 11.1 | 375
Magurele_CHM15k | CHM15k | earlinet | 51 | 4.7 | 12.2 | 11.1 | 7.6 | 50
Magurele_CHM15k | CHM15k | bellini | 46 | 11.4 | 12.6 | 8.8 | 6.9 | 24
Oslo_CHM15k | CHM15k | eprof_v10 | 38 | 5.1 | 12.8 | 10.1 | 7.6 | 623
Oslo_CHM15k | CHM15k | main | 48 | 31.1 | 14.6 | 10.4 | 6.9 | 563
Oslo_CHM15k | CHM15k | improved | 35 | 13.2 | 9.3 | 8.0 | 6.4 | 809
Oslo_CHM15k | CHM15k | optimal | 41 | 12.0 | 8.2 | 9.0 | 4.0 | 43
Oslo_CHM15k | CHM15k | matlab | 42 | 26.5 | 13.5 | 10.0 | 6.0 | 541
Oslo_CHM15k | CHM15k | calipso | 70 | 18.7 | 15.5 | 11.6 | 9.5 | 387
Oslo_CHM15k | CHM15k | earlinet | 36 | 6.2 | 8.5 | 7.7 | 4.4 | 43
Oslo_CHM15k | CHM15k | bellini | 28 | 18.3 | 13.1 | 11.8 | 6.9 | 35
Palaiseau_CHM15k | CHM15k | eprof_v10 | 44 | 7.1 | 11.1 | 9.0 | 6.7 | 605
Palaiseau_CHM15k | CHM15k | main | 56 | 31.1 | 13.1 | 12.1 | 8.0 | 574
Palaiseau_CHM15k | CHM15k | improved | 40 | 14.5 | 13.0 | 10.5 | 7.8 | 616
Palaiseau_CHM15k | CHM15k | optimal | 36 | 11.8 | 12.5 | 9.8 | 8.2 | 34
Palaiseau_CHM15k | CHM15k | matlab | 51 | 26.2 | 13.1 | 11.7 | 9.5 | 572
Palaiseau_CHM15k | CHM15k | calipso | 80 | 21.6 | 18.5 | 13.6 | 11.7 | 623
Palaiseau_CHM15k | CHM15k | earlinet | 32 | 6.2 | 12.9 | 10.9 | 7.9 | 32
Palaiseau_CHM15k | CHM15k | bellini | 30 | 16.8 | 12.8 | 12.7 | 7.9 | 32
Payerne_CHM15k | CHM15k | eprof_v10 | 31 | 5.7 | 17.6 | 15.8 | 10.1 | 1267
Payerne_CHM15k | CHM15k | main | 42 | 31.7 | 18.2 | 12.3 | 9.9 | 1455
Payerne_CHM15k | CHM15k | improved | 26 | 9.1 | 25.7 | 20.4 | 11.5 | 1149
Payerne_CHM15k | CHM15k | optimal | 29 | 11.6 | 13.8 | 12.1 | 9.6 | 1207
Payerne_CHM15k | CHM15k | matlab | 36 | 24.6 | 17.0 | 14.1 | 9.0 | 574
Payerne_CHM15k | CHM15k | calipso | 75 | 19.0 | 31.1 | 20.1 | 17.8 | 1534
Payerne_CHM15k | CHM15k | earlinet | 26 | 5.3 | 12.4 | 11.3 | 8.3 | 38
Payerne_CHM15k | CHM15k | bellini | 25 | 11.9 | 14.8 | 14.1 | 9.0 | 54
SIRTA-MPL_MPL | Mini-MPL | eprof_v10 | 63 | 2.6 | 24.8 | 20.3 | 15.9 | 46
SIRTA-MPL_MPL | Mini-MPL | main | 71 | 1.9 | 27.2 | 26.3 | 22.1 | 53
SIRTA-MPL_MPL | Mini-MPL | improved | 69 | 1.7 | 26.9 | 25.7 | 20.5 | 48
SIRTA-MPL_MPL | Mini-MPL | optimal | 65 | 2.2 | 22.8 | 19.4 | 15.3 | 40
SIRTA-MPL_MPL | Mini-MPL | matlab | 69 | 1.8 | 27.1 | 27.3 | 21.0 | 51
SIRTA-MPL_MPL | Mini-MPL | calipso | 86 | 3.1 | 36.2 | 32.0 | 25.8 | 63
SIRTA-MPL_MPL | Mini-MPL | earlinet | 71 | 1.7 | 27.6 | 22.2 | 20.4 | 47
SIRTA-MPL_MPL | Mini-MPL | bellini | 60 | 1.4 | 28.8 | 27.6 | 20.5 | 49
Toulouse-MPL_MPL | Mini-MPL | eprof_v10 | 65 | 2.1 | 28.5 | 20.9 | 16.5 | 61
Toulouse-MPL_MPL | Mini-MPL | main | 77 | 2.2 | 36.7 | 29.0 | 24.4 | 66
Toulouse-MPL_MPL | Mini-MPL | improved | 70 | 2.2 | 29.0 | 22.8 | 20.5 | 58
Toulouse-MPL_MPL | Mini-MPL | optimal | 67 | 2.3 | 23.9 | 21.2 | 12.8 | 49
Toulouse-MPL_MPL | Mini-MPL | matlab | 74 | 2.3 | 32.1 | 25.1 | 20.5 | 63
Toulouse-MPL_MPL | Mini-MPL | calipso | 85 | 2.9 | 35.1 | 29.7 | 24.8 | 67
Toulouse-MPL_MPL | Mini-MPL | earlinet | 73 | 1.4 | 29.5 | 30.8 | 19.7 | 55
Toulouse-MPL_MPL | Mini-MPL | bellini | 58 | 1.4 | 28.3 | 27.0 | 18.5 | 56

### 3.5 Method-comparison conclusion

- **`eprof_v2` and `earlinet` give the cleanest, most stable long-term calibration** — lowest robust
  CV (31–32 %), lowest temporal CV (0.13–0.15), highest R² for `eprof_v2` (0.97). They reject the
  outlier nights that inflate the other methods' std-CV. **This is the recommendation for production**,
  where per-night quality and stability matter more than raw count (the Kalman handles sparsity).
- **`eprof_v1.2` remains a sound, simpler default**; use `eprof_v2` (or `earlinet` where a
  single-profile method is preferred) for the best precision.
- **Avoid `calipso` (retired) and `eprof_v1.1`/`main`** for quantitative calibration.
- **`bellini` (ALICENET)** is a solid, citable alternative, strongest on high-SNR CHM15k; on the
  lower-SNR 910 nm ALCs its 3–7 km band yields fewer points.
- **CL31 / CL51 (910 nm, noisy)** rarely yield a clean molecular window under any gated method —
  confirming they should use the **liquid-cloud calibration** as primary, Rayleigh only opportunistically.
- **Feed the per-night constants to the Kalman**, which absorbs the residual outliers.

---

## 4. v2 optimization — why clear nights fail, and the tuned C8

*Generated 2026-06-20; L2 recomputed 2026-06-21 after the coarse-cadence cloud-screening fix.
`validation/run_v2_sweep.py` + `analyze_v2_sweep.py`. 24 instruments (10 CHM15k, 4 Mini-MPL, 10 CL61)
× both levels × every clear night of 2026.*

**Why clear nights fail:** the **scattering-ratio gate (old `max_scattering_ratio = 1.10`) is the
dominant binding constraint on BOTH levels** — leave-one-gate-out recovers 77 % of failed L1-CL61
nights, 43 % L1-CHM15k, and 66–77 % of the (now far fewer) L2 failures. For **Mini-MPL** the binding
gate is the **window-start height** (its 532 nm signal runs out of SNR before 2 km). R² and ratio-std
never bind.

![Leave-one-gate-out: % of baseline-failed clear nights recovered by relaxing only that gate. L1 left, L2 right. Scattering ratio dominates on both levels; Mini-MPL is window-start-limited.](figs_extracted/v2_optimization_report_01.png)

| level · type | failed nights | #1 binding gate | #2 |
|---|---|---|---|
| L1 CHM15k | 263 | **scattering 43 %** | temporal_cv 31 % |
| L1 Mini-MPL | 154 | **window start 21 %** | scattering 13 % |
| L1 CL61 | 171 | **scattering 77 %** | residual 5 % |
| L2 CHM15k | 83 | **scattering 66 %** | residual 36 % |
| L2 CL61 | 146 | **scattering 77 %** | residual 13 % |

The estimated scattering ratio sits just above 1.10 on many genuinely clear nights, so the old gate
threw away usable molecular windows — overwhelmingly for **CL61**.

**The configuration trade-off (nine configs; gates not listed sit at the C0 baseline;
`max_rel_error = 15` always):**

| config | scatter | temporal_cv | residual | R² | start (m) | ratio_std |
|---|---|---|---|---|---|---|
| **C0** baseline (old v2) | 1.10 | 0.50 | 12 | 0.50 | 2000 | 0.30 |
| C3 looser shape/ratio | **1.15** | 0.50 | 20 | 0.50 | 2000 | 0.40 |
| C6 balanced | 1.12 | 0.80 | 16 | 0.40 | 1500 | 0.40 |
| C7 aggressive | 1.25 | 1.50 | 25 | 0.30 | 1000 | 0.50 |
| **C8 RECOMMENDED** | **1.15** | 0.80 | 16 | 0.40 | 1500 | 0.40 |
| *(C1 tcv≤0.8, C2 tcv≤1.2, C4 R²≥0.35, C5 start≥1.2 km also swept)* | | | | | | |

**valid % / σ_SD % (median over instruments of each type):**

**L1** (native, noisy — relaxation matters most)

| config | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| C0 baseline | 38.2 / 10.9 | 65.0 / 8.5 | 45.4 / 7.6 |
| C3 looser shape/ratio | 46.6 / 10.3 | 66.3 / 8.8 | 64.9 / **6.5** |
| C6 balanced | 67.5 / 9.6 | 68.9 / 8.8 | 68.5 / 7.3 |
| C7 aggressive | **93.4** / 9.2 | 72.8 / 10.9 | **85.9** / 7.1 |
| **C8 recommended** | 71.6 / 10.0 | 69.7 / 8.9 | 75.1 / **7.0** |

**L2** (averaged, clean — baseline already high; C8 neutral)

| config | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| C0 baseline | 91.0 / 9.3 | 100 / 6.4 | 75.2 / 7.0 |
| C3 looser shape/ratio | 93.2 / 9.1 | 100 / 6.4 | 77.1 / 7.8 |
| C6 balanced | 90.4 / **8.2** | 100 / 6.4 | 71.2 / 7.7 |
| C7 aggressive | **97.0** / 10.1 | 100 / 6.4 | 85.7 / 8.2 |
| **C8 recommended** | 91.0 / 9.3 | 100 / 6.4 | 75.2 / **7.0** |

![Yield vs short-term variability for every config, per type. Best is bottom-right. Circles = L1, squares = L2.](figs_extracted/v2_optimization_report_02.png)

![Valid-calibration fraction on clear nights by config, grouped by type (L1 left, L2 right).](figs_extracted/v2_optimization_report_03.png)

On **L1**, C0 → C8 pushes every type up-and-right (more valid, σ_SD flat/lower) — the gate relaxation
unlocks the noisy native data. On **L2**, the baseline is already ~91 % (CHM15k) / 75 % (CL61) because
averaging produces clean profiles, so C8 sits on top of C0 (neutral); only the aggressive C7 buys more
yield, at a σ_SD cost. C8 captures the L1 gains without disturbing L2.

**Recommendation — C8, the `eprof_v2` default.** Gates set to
`min_window_start_m=1500, min_r2=0.40, max_residual_pct=16, max_scattering_ratio=1.15,
max_ratio_std=0.40, max_temporal_cv=0.8` (weights & time-cell flagging unchanged; verified live in
`calibration/rayleigh/molecular_methods.py`). Net valid % vs the old baseline:

| | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| L1 | 38 → **72** (+34) | 65 → **70** (+5) | 45 → **75** (+30) |
| L2 | 91 → 91 (≈) | 100 → 100 | 75 → 75 (≈) |

σ_SD held flat or improved everywhere (L1-CL61 7.6 → 7.0; L2-CL61 7.0 → 7.0). C8 is a clear win on L1
and neutral on L2 — strictly safe as the global default.

**Adversarial review & caveats.** An independent 4-analyst review verified the (then L1-correct)
numbers and endorsed C8. The one substantive caveat it raised — "C8 inflates L2-CL61 σ_SD (6.9 → 8.7)"
— was **resolved**: it was the coarse-cadence screening bug admitting marginal cloudy nights; corrected
C8 L2-CL61 σ_SD is 7.0. Remaining caveats: (1) **no out-of-year / out-of-station holdout** — thresholds
tuned on 2026; validate on a second year and a station holdout before treating them as permanent; (2)
the leave-one-out diagnostic is **univariate** (correlated gates), so the 43–77 % recoveries are upper
bounds, but the *ranking* (scattering #1) is robust; (3) **per-instrument-type tuning is the natural
v2.1** — Mini-MPL is window-start-limited, CL61/CHM15k are scattering-limited; a per-type policy would
be Pareto-better than one global config, but ship global C8 now.

---

## 5. Network diagnosis — problematic-station causes, L1 vs L2

*Generated 2026-06-21. `validation/run_rayleigh_diag_light.py` + `analyze_rayleigh_diag.py`, on the v2
constants from the network run. All 153 CHM15k + Mini-MPL streams, 2026. The cause diagnosis uses L1
(instrument housekeeping lives only in native L1); the "L1 vs L2" section compares the calibration
outcome.*

### 5.1 Metrics and causes

Per stream (optimized v2): **clear%** = clear (fit-reaching) nights ÷ archive days;
**valid% (of clear)** = valid calibrations ÷ clear nights (method yield); **valid% (of archive)** =
availability × yield (overall productivity); **σ_SD**, **outlier%** = short-term variability and
drift-aware outlier rate. Four diagnostics map to candidate failure causes:

| cause | diagnostic | "bad" signature (vs healthy CHM15k median) |
|---|---|---|
| not enough clear sky | clear% | **11 %** vs 33 % |
| lots of FT aerosol | median window scattering ratio | **1.18** vs 1.12 |
| low laser | near-range signal strength → SNR proxy (+ `laser_life_time`) | SNR **13** vs 33 |
| electronic background | far-range noise → SNR proxy | SNR **22** vs 33 (high noise) |

A station is **problematic** if (within its type) it falls in the worst quartile of overall yield, or
the worst decile of σ_SD or outlier rate. Each is assigned the **dominant** cause (largest
standardized exceedance of the per-type threshold); secondary causes are listed too.

### 5.2 Result — 59 of 153 streams problematic

cause | # streams | what it looks like
---|---|---
ok | 94 | clear 33 %, scat 1.12, SNR 33, σ_SD 8.9 % — healthy
electronic background | 19 | low SNR from high far-range noise (clear sky fine)
not enough clear sky | 16 | clear% ~11 % — Arctic/maritime, persistent cloud
low laser | 12 | low SNR from weak near-range signal; often old laser
FT aerosol | 6 | scattering ratio ≥1.18 — persistent free-tropospheric aerosol
other | 6 | flagged on σ_SD/outliers without a single dominant metric

![CHM15k network diagnosis: (top-left) overall yield sorted, coloured by cause; (top-right) clear-sky availability drives yield; (bottom-left) FT aerosol vs σ_SD; (bottom-right) SNR vs outlier rate.](figs_extracted/rayleigh_network_diagnosis_report_01.png)

**What is wrong with the problematic stations:**

- **Not enough clear sky (16).** clear% ≈ 11 % (a third of healthy). High-latitude/maritime sites with
  persistent cloud — **Jan Mayen, Bjørnøya, Hopen, Tórshavn** (Arctic), **Bonaire** (trade-cumulus),
  **Camborne**. Not an instrument fault — too few clear nights.
- **Low laser (12).** Low SNR from *weak near-range signal*, often with high `laser_life_time`:
  **Vásárosnamény** (laser ≈ 53 900 h), **Nieuwkoop**, **Gottfrieding**. Ageing/under-powered laser →
  low SNR in the 2–6 km band. Action: laser service / re-collimation. **Persists on L2.**
- **Electronic background (19).** Low SNR from *high far-range noise* (the fit band itself is fine):
  **Friesoythe, Klippeneck, Bern**. Elevated detector/electronic background raises the noise floor.
  Action: check detector baseline / dark-current. **Persists on L2.**
- **FT aerosol (6).** Persistently elevated window scattering ratio (≥1.18): **Payerne (1.24), Bern**
  — frequent free-tropospheric aerosol. v2's scattering gate rejects nights even when skies look clear.
  Action: relaxed-scattering v2 (C8) + higher fit windows / per-site gate.

> Payerne: clear% 43 % (plenty of clear nights) yet valid ≈ 1 %, scat 1.24 AND SNR 20 — limited by
> **both** FT aerosol and background, not clear-sky availability.

- **Mini-MPL (5):** Aléria (83 %) and Lille (55 %) healthy; Brest *background*, Trappes *FT-aerosol*,
  the mobile unit *low-laser*.

### 5.3 L1 vs L2 — same streams, same v2 method

> **Screening-bug fix (2026-06-21).** An earlier version reported L2 reaching the fit on ~98 % of days
> vs ~32 % on L1 — physically impossible for the same data. Root cause: `filter_cloudy_profiles`
> rounded `profiles_per_min` to an integer; at L2's 5-min cadence that is `round(0.2)=0`, which zeroed
> the minimum-clear-profile threshold and **switched cloud screening off**, so cloudy L2 nights passed
> as clear. (It is *not* missing data, and `vertical_visibility` *is* read for L2.) After the fix
> (`profiles_per_min` kept as a float), L1 and L2 agree to **6 pp** median clear-reach difference (was
> ~66 pp). The numbers below are post-fix.

Medians over streams (v2):

| | clear% L1→L2 | valid/clear% L1→L2 | σ_SD% L1→L2 | outlier% L1→L2 |
|---|---|---|---|---|
| CHM15k (n=148) | 32 → 24 | 73 → 88 | 9.2 → 9.5 | 4 → 3 |
| Mini-MPL (n=5) | 13 → 13 | 100 → 100 | 6.5 → 6.4 | 4 → 4 |

![Paired per-stream L1 vs L2 (v2): clear-sky reach, yield-on-clear, σ_SD, outlier rate. Colour = L1-diagnosed cause. Panel 1 now sits on the diagonal (L1≈L2); panel 2 above it (L2 higher yield-on-clear).](figs_extracted/rayleigh_network_diagnosis_report_02.png)

**L1 and L2 now agree on clear-sky availability** — the clear-reach points lie on the diagonal (CHM15k
32 % vs 24 %; Mini-MPL identical at 13 %, as it must be since its L1 and L2 share the same 5-min ×
30-m grid). The remaining, *real* difference is **yield-on-clear: L2 is higher** — its time/range
averaging produces a cleaner molecular profile, so more clear nights pass the rel-error QC — at
**comparable σ_SD**.

**Per cause, L1 → L2 (CHM15k median clear% / valid-on-clear% / σ_SD):**

cause | L1 | L2 | reading
---|---|---|---
ok | 33 / 74 / 8.9 | 27 / 90 / 9.1 | same clear-sky; L2 turns more clear nights into valid ones; σ_SD ~equal
not enough clear sky | 11 / 82 / 12.2 | 8 / 92 / 8.5 | clear-sky-limited on **both** levels (atmosphere, not instrument)
low laser | 34 / 43 / 10.2 | 28 / 67 / 11.3 | **persists** — L2 averaging recovers some nights but σ_SD rises; weak laser still limits
electronic background | 30 / 35 / 10.0 | 23 / 69 / 12.3 | **persists** — L2 recovers nights but σ_SD *worsens* (10.0 → 12.3): noisy detector
FT aerosol | 35 / 76 / 9.7 | 26 / 81 / 10.9 | aerosol limits both levels similarly

**Takeaway on L1 vs L2:** they are consistent (same measurement). L2's only genuine edge is **higher
valid-on-clear from averaging** (cleaner fits), at similar precision. The two instrument-health causes
(**low-laser, electronic-background**) still show low yield and the σ_SD signature on **both levels**,
confirming they are hardware issues, not artefacts (averaging cannot recover a weak laser or a noisy
detector; for background it even raises σ_SD). For maximum usable calibrations, **L2** is slightly
preferable; for the instrument-health diagnosis you need **L1** (it carries the housekeeping).

### 5.4 Time series — full network

Every CHM15k stream's nightly lidar constant (v2, L1), sorted worst-yield first; green = robust
median, red = outliers, pink panels = problematic. 6 pages cover all 148 CHM15k; Mini-MPL separately.

![CHM15k time series page 1/6 (lowest-yield streams)](figs_extracted/rayleigh_network_diagnosis_report_03.png)
![CHM15k time series page 2/6](figs_extracted/rayleigh_network_diagnosis_report_04.png)
![CHM15k time series page 3/6](figs_extracted/rayleigh_network_diagnosis_report_05.png)
![CHM15k time series page 4/6](figs_extracted/rayleigh_network_diagnosis_report_06.png)
![CHM15k time series page 5/6](figs_extracted/rayleigh_network_diagnosis_report_07.png)
![CHM15k time series page 6/6 (highest-yield streams)](figs_extracted/rayleigh_network_diagnosis_report_08.png)
![Mini-MPL time series](figs_extracted/rayleigh_network_diagnosis_report_09.png)

**Problematic-station list (worst by yield; full 59-stream table in
`figs_paper_validation/rayleigh_diag/problematic_stations.md`):**

site | type | valid% | clear% | σ_SD% | scat | SNR | laser(h) | cause
---|---|---|---|---|---|---|---|---
Tórshavn | CHM15k | 0 | 0 | – | – | 91 | 16901 | clear-sky
QUALAIR | CHM15k | 0 | 0 | – | – | 30 | 42629 | clear-sky
Friesoythe | CHM15k | 0 | 31 | – | 1.18 | 61 | 26136 | background
Bonaire | CHM15k | 0 | 4 | – | 1.15 | 76 | 16012 | background→clear-sky
Payerne | CHM15k | 1 | 43 | – | 1.24 | 20 | 8934 | FT-aerosol + background
Vásárosnamény | CHM15k | 1 | 59 | – | 1.19 | 12 | 53948 | low-laser
Jan Mayen | CHM15k | 2 | 2 | – | 1.08 | 25 | 49507 | clear-sky
Nieuwkoop | CHM15k | 2 | 34 | – | 1.13 | 13 | 9390 | low-laser

**Takeaways:** ~60 % of the CHM15k network calibrates healthily; the problematic 40 % split into four
physically distinct, separable causes. **Two are not instrument faults** — *clear-sky-limited*
(Arctic/maritime) and *FT-aerosol* sites are atmosphere-limited (need longer accumulation or per-site
gate relaxation, not repair). **Two are instrument-health flags worth acting on** — *low-laser* and
*electronic-background* — and **both persist on L2**, confirming they are hardware, not screening,
issues; the far-range-noise and laser-age metrics can drive an automated network health alert.

---

## 6. Night-to-night variability — L1 2026 (CHM15k / Mini-MPL / CL61)

*MeteoSwiss E-PROFILE ALC paper — M. Hervo. Generated 2026-06-19.* This section estimates the
**night-to-night variability of the Rayleigh calibration constant** for the selected CHM15k, the
Mini-MPL, and **all CL61**, processed directly from the **L1 2026** daily archive
(`D:/E-PROFILE_L1_2026`, Feb–May 2026, every night). Methods are named by E-PROF version; the
historical **E-PROF v1.0 (sign error)** is included, and the CL61 Rayleigh result is cross-checked
against an **independent liquid-cloud calibration**.

### 6.1 Data and method

E-PROFILE **L1** daily files, every available night; darkness-adaptive night selection (SZA > 100°),
cloud/fog screening, profile-outlier screening, then the night-mean range-normalised signal + per-
profile stack feed the molecular-window search. **910 nm units (CL61) get the WV correction** (bundled
910 nm LUT + CAMS 2026); CHM15k (1064 nm) and Mini-MPL (532 nm) are unaffected. Because the prepared
profile is built *before* window selection and is method-independent, all series (the live methods +
**E-PROF v1.0**, the legacy window run through the full pipeline with `sign_error_v10`) come from **one
load per night**. Instruments: CHM15k — Payerne, Lindenberg, Aosta, Palaiseau, Granada, Magurele,
Bergen, Oslo, Hamburg, Hohenpeissenberg; Mini-MPL — Brest, Mini-MPL-Mobile, Aleria, Lille; CL61 —
Camborne, Zeebrugge, Birkenes, Lauder, Aosta, Uccle, Lanzhot, Sion/EPFL, Temelín, Edmonton.

Variability metrics are identical to §3.3 (σ_SD headline; σ_detrend; σ_within-month; σ_night; rob_CV;
valid %).

### 6.2 Which molecular selection is best

*(Figure — molecular-window methods on L1 2026: yield (left) and drift-insensitive precision (right).
Source figure `l1_2026_variability/method_precision_l1_2026.png` not committed to this repo.)*

Mean over all 24 instruments:

| method | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_within-month % | rob_CV % | CV % |
|---|--:|--:|--:|--:|--:|--:|--:|
| **E-PROF v2** | 44 | 6.7 | **9.5** | 9.2 | 9.9 | 13.0 | 13 |
| E-PROF v1.0 (sign error) | 60 | 4.4 | **10.3** | 9.2 | 10.2 | 13.3 | 37 |
| **EARLINET/SCC** | 51 | 4.9 | **10.5** | 10.1 | 10.7 | 13.9 | 16 |
| E-PROF v0.25 (MATLAB) | 63 | 21.6 | **11.4** | 10.6 | 11.5 | 14.8 | 23 |
| E-PROF v1.1 (sign cor) | 68 | 23.9 | **11.6** | 10.9 | 12.1 | 15.5 | 31 |
| E-PROF v1.2 (improved) | 54 | 12.7 | **12.9** | 10.8 | 12.3 | 15.3 | 59 |
| Bellini/ALICENET | 34 | 12.6 | **14.1** | 15.4 | 14.3 | 17.4 | 26 |

- **E-PROF v2 is the most precise** (σ_SD 9.5 %); **EARLINET** is a close second (10.5 %) with the best
  yield among the quality methods (51 %).
- **E-PROF v1.1 / v0.25 look stable night-to-night (σ_SD 11–12 %) but carry a huge within-night spread
  (σ_night 22–24 %)**: lacking an above-aerosol gate they admit aerosol windows — a *systematic*
  accuracy risk σ_SD alone doesn't reveal.
- **E-PROF v1.2 (the current production default) is mid-pack here (σ_SD 12.9 %, CV 59 %)** — a few
  spurious-window nights inflate its CV; v2/EARLINET's stricter gates avoid them. This is the case for
  the repo switching its default to **E-PROF v2**.
- **Bellini/ALICENET is last (σ_SD 14.1 %, yield 34 %)** and **fails entirely on CL61** (3–7 km band,
  no SNR at 910 nm).

**Per instrument type:**

| group | most precise (σ_SD) | E-PROF v2 | note |
|---|---|---|---|
| **CHM15k** (n=10) | EARLINET 9.2 % | 10.1 % (valid 36 %) | EARLINET/v2 lead |
| **Mini-MPL** (n=4) | **E-PROF v2** 11.1 % | 11.1 % (valid 62 %) | v2 clearly best; others 16–18 % |
| **CL61** (n=10) | E-PROF v0.25 7.7 % | 8.1 % (valid 45 %) | v0.25/v1.1 σ_SD low **but σ_night 22–24 %** (aerosol); **E-PROF v2 gives 8.1 % at σ_night 9 %** — the best *clean* choice. Bellini yields 0. |

*(Figure — per-group method precision σ_SD, lower = better, bar label = valid%. Source figure
`l1_2026_variability/method_precision_by_group.png` not committed to this repo.)*

For CL61, v0.25/v1.1 reach the lowest σ_SD only because they take a window every night, but their
22–24 % within-night spread shows the windows include aerosol; **E-PROF v2 matches their night-to-
night precision while keeping the within-night spread at 9 %**, so it is the method to use.

### 6.3 Per-instrument variability

*(Figure — per-instrument σ_SD and yield, recommended method per instrument, colour = type. Source
figure `l1_2026_variability/per_instrument_variability.png` not committed to this repo.)*

Each instrument is characterised with its **recommended method** — the most precise quality method
(preference E-PROF v2 → EARLINET → E-PROF v1.2) with ≥ 6 valid nights.

| instrument | type | fit-nights | method | valid | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_within-month % | rob_CV % |
|---|---|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| Hohenpeiss | CHM15k | 40 | E-PROF v2 | 16 | 40 | 5.9 | **4.5** | 6.5 | 6.1 | 10.6 |
| Hamburg | CHM15k | 34 | E-PROF v2 | 19 | 56 | 6.7 | **7.9** | 6.4 | 5.3 | 7.6 |
| Magurele | CHM15k | 35 | E-PROF v2 | 25 | 71 | 10.5 | **8.6** | 8.0 | 7.6 | 10.7 |
| Oslo | CHM15k | 18 | E-PROF v2 | 14 | 78 | 4.3 | **10.9** | 7.5 | 8.1 | 9.9 |
| Granada | CHM15k | 54 | E-PROF v2 | 13 | 24 | 5.8 | **11.4** | 9.6 | 8.7 | 11.3 |
| Lindenberg | CHM15k | 43 | E-PROF v2 | 22 | 51 | 7.3 | **13.3** | 11.0 | 12.7 | 14.8 |
| Aosta | CHM15k | 44 | E-PROF v2 | 16 | 36 | 6.8 | **13.9** | 14.8 | 14.9 | 19.9 |
| Bergen | CHM15k | 16 | EARLINET | 7 | 44 | 7.0 | **4.5** | 5.7 | 6.4 | 6.1 |
| Palaiseau | CHM15k | 53 | EARLINET | 10 | 19 | 8.4 | **7.6** | 2.9 | 4.8 | 8.0 |
| Payerne | CHM15k | 51 | E-PROF v1.2 | 7 | 14 | 31.7 | **14.2** | 15.3 | 8.4 | 9.3 |
| Mini-MPL-Mobile | Mini-MPL | 119 | E-PROF v2 | 85 | 71 | 2.2 | **6.8** | 5.7 | 8.6 | 13.6 |
| Aleria | Mini-MPL | 117 | E-PROF v2 | 95 | 81 | 1.6 | **6.8** | 7.5 | 15.6 | 17.7 |
| Lille | Mini-MPL | 75 | E-PROF v2 | 44 | 59 | 1.9 | **10.2** | 12.1 | 11.8 | 12.3 |
| Brest | Mini-MPL | 109 | E-PROF v2 | 42 | 39 | 2.4 | **20.5** | 23.5 | 25.9 | 42.1 |
| Lanzhot | CL61 | 49 | E-PROF v2 | 35 | 71 | 8.8 | **4.7** | 7.2 | 7.5 | 8.1 |
| Lauder | CL61 | 55 | E-PROF v2 | 42 | 76 | 6.9 | **6.0** | 7.7 | 10.1 | 11.4 |
| Sion / EPFL | CL61 | 50 | E-PROF v2 | 22 | 44 | 10.1 | **6.3** | 6.6 | 5.8 | 7.7 |
| Aosta | CL61 | 56 | E-PROF v2 | 41 | 73 | 7.8 | **7.5** | 8.1 | 7.6 | 11.5 |
| Temelín | CL61 | 32 | E-PROF v2 | 15 | 47 | 10.8 | **7.7** | 5.8 | 5.8 | 5.9 |
| Uccle | CL61 | 32 | E-PROF v1.1 | 9 | 28 | 16.0 | **7.7** | 9.3 | 9.7 | 6.0 |
| Birkenes | CL61 | 17 | E-PROF v2 | 14 | 82 | 9.3 | **7.9** | 5.7 | 5.2 | 6.2 |
| Camborne | CL61 | 23 | E-PROF v2 | 8 | 35 | 10.6 | **10.4** | 4.6 | 5.1 | 5.1 |
| Zeebrugge | CL61 | 31 | E-PROF v2 | 6 | 19 | 7.5 | **14.5** | 15.6 | 15.5 | 20.0 |
| Edmonton | CL61 | 9 | — | 0 | 0 | — | **—** | — | — | — |

Per-instrument **median σ_SD by type**: **CHM15k 9.7 %** (4.5–14.2), **Mini-MPL 8.5 %** (6.8–20.5),
**CL61 7.7 %** (4.7–14.5).

- **An irreducible night-to-night floor of ≈ 5 %** is reached by the best units of every type
  (Hohenpeissenberg & Bergen CHM15k 4.5 %, Lanzhot CL61 4.7 %).
- **CL61 is the most precise group** (median 7.7 %) — the headline new result (§6.4).
- **Problem units stand out**: Brest Mini-MPL (σ_SD 20.5 %), Zeebrugge CL61 (14.5 %), and Payerne
  CHM15k (only E-PROF v1.2 on 7 nights, σ_night 31.7 % — its spring-2026 aerosol defeats the strict
  gates).

*(Figure — per-instrument calibration-constant time series, E-PROF v1.2 vs E-PROF v2, title colour =
type. Source figure `l1_2026_variability/timeseries_l1_2026.png` not committed to this repo.)*

The E-PROF v2 series (red) sit in a tighter band than E-PROF v1.2 (blue), whose spikes are the
spurious-window nights that inflate its CV.

### 6.4 CL61 cross-check — Rayleigh vs liquid-cloud calibration

The CL61 Rayleigh calibration is cross-checked against the **independent liquid-water-cloud
calibration** (O'Connor/Hopkin) run on the monthly L2 archive (the cloud method needs a month of
profiles to accumulate enough liquid-cloud returns). The two methods calibrate *different physical
constants* — the Rayleigh L1 lidar constant vs a multiplier on the L2 attenuated backscatter — so only
their **precision and consistency** are comparable.

*(Figure — CL61 cross-check: liquid-cloud vs Rayleigh calibration precision. Source figure
`l1_2026_variability/cloud_vs_rayleigh.png` not committed to this repo.)*

| instrument | cloud months | cloud σ_within-month % | cloud σ_month-to-month % | Rayleigh n | Rayleigh σ_night % | Rayleigh σ_SD % |
|---|--:|--:|--:|--:|--:|--:|
| Birkenes | 2 | 11.4 | 0.9 | 14 | 9.3 | 7.9 |
| Camborne | 2 | 45.7 | 26.5 | 8 | 10.6 | 10.4 |
| Lanzhot | 2 | 7.2 | 4.0 | 35 | 8.8 | 4.7 |
| Lauder | 2 | 10.7 | 11.8 | 42 | 6.9 | 6.0 |
| Aosta | 2 | 16.8 | 0.7 | 41 | 7.8 | 7.5 |
| Sion / EPFL | 2 | 8.4 | 5.6 | 22 | 10.1 | 6.3 |
| Temelín | 1 | 16.5 | – | 15 | 10.8 | 7.7 |
| Uccle | 4 | 27.3 | 21.4 | 0 | – | – |
| Zeebrugge | 4 | 36.0 | **105.3** | 6 | 7.5 | **14.5** |
| Edmonton | 1 | 3.0 | – | 0 | – | – |

The cloud calibration runs with the **above-cloud aerosol two-way transmission correction applied**;
its effect is a physically sensible −4 to −6 % on the coefficient.

- **On the well-exposed CL61 the two independent methods agree at the ~10 % level**: Lanzhot (cloud
  σ_within-month 7.2 %, Rayleigh σ_night 8.8 %), Sion (8.4 / 10.1), Birkenes (11.4 / 9.3), Lauder
  (10.7 / 6.9). Strong corroboration — a method that does *not* use a molecular window reproduces the
  Rayleigh-derived calibratability.
- **Both methods independently flag Zeebrugge as anomalous** (cloud month-to-month 105 %, cloud
  coefficient ~10× the others; Rayleigh σ_SD 14.5 %, the worst CL61) — an instrument/site issue, not a
  method artefact.
- **Camborne and Uccle are noisier in the cloud method** (46 %, 27 %) on only 2–4 months with few clean
  liquid clouds; their Rayleigh values are the more reliable there.
- The mean cloud within-month scatter (18 %, inflated by Camborne/Uccle/Zeebrugge) is larger than the
  Rayleigh σ_SD (8 %), as expected for a sparse monthly cloud sample — but where both have data they
  **agree on which CL61 are stable and which are not**, which is the point of a cross-check.

> *Historical fix note.* An earlier `apply_transmission_correction` step produced a zero-median
> coefficient on L2 input — it integrated the **uncalibrated** stored backscatter (~1e6× too large for
> L2's units), so `T² = exp(−2·LR·B_aerosol)` underflowed to 0. The integral is now converted to a
> physical optical depth with the per-profile coefficient `AOD = LR·(S/S_THEORETICAL)·∫β dr` before
> the Beer–Lambert factor, giving the −4 to −6 % correction. (The legacy MATLAB reference carried the
> same latent bug.) *For the current authoritative liquid-cloud method and its fixed 100–2400 m gate,
> multiple-scattering η tables and mandatory WV correction, see the liquid-cloud calibration report.*

### 6.5 The CL61 fleet, and the sign error

**CL61 (new).** The long-run study (§3) had no CL61 (its L2-monthly archive and the old
`instruments.json` predate the fleet). Running the **10 CL61 from L1 2026 with the WV correction**
gives **CL61 median per-instrument σ_SD 7.7 %**, the *most precise* of the three groups, with high
yield on the well-exposed units (Lauder 76 %, Birkenes 82 %, Aosta 73 %, Lanzhot 71 %) — corroborated
by the liquid-cloud cross-check (§6.4).

**E-PROF v1.0 (sign error).** Including the historical sign-error baseline shows that **the Klett sign
error does not change the calibration-constant stability**: E-PROF v1.0 σ_SD 10.3 % vs E-PROF v1.1
(sign-corrected) 11.6 % — essentially the same (both use the legacy `main` window). The sign error
corrupts the *downstream attenuated-backscatter / AOD product*, **not** the Rayleigh calibration
constant `C_L`. This isolates the sign error's impact to β_att and confirms `C_L` is unaffected — a
useful, citable result given the operational network still runs v1.0.

**`instruments.json` rebuilt.** The shipped manifest predated the CL61 fleet (zero CL61) and carried
unused fields. `scripts/data/rebuild_instruments_json.py` rebuilds it from each stream's NetCDF
metadata (one entry per `(WMO, identifier)`), keeping only fields **in the files** or **used by the
code** (`WMO, Identifier, Type, SiteName, Latitude, Longitude, Altitude, Serial, Calibrated`) and
dropping unused ones (`Reference, FLength, NWS, Status`). Rebuilt from the **May 2026** fleet: **415
instruments** incl. **10 CL61** + 145 CHM15k + 4 Mini-MPL (+ 208 CL31, 48 CL51), 29 multi-instrument
stations. (May 2025 predates the CL61 rollout and yields 0 CL61, so the latest comprehensive month was
used.)

---

## Overall recommendations

1. **Repo production: use `eprof_v2` (config C8)** for CHM15k / CL61 / Mini-MPL Rayleigh calibration —
   most precise (σ_SD 9.5 % on L1 2026; 13.6 % long-run), robust to spurious windows. Use **`earlinet`**
   where a single-profile method is preferred (close second, higher yield) and as the fallback when v2
   yields too few nights. Both beat the older `eprof_v1.2` default on precision.
2. **Deployed operational network still runs E-PROF v1.0** (with the historical sign error) — distinct
   from the repo default; the sign error affects β_att, not `C_L`.
3. **Keep `calipso` retired** — physically inappropriate for ground-up ALCs (no stratospheric reference).
4. **Apply the WV correction on the 910 nm units (CL61)** — it brings them to CHM15k-class stability;
   it is mandatory for 910 nm (a no-WV mode is rejected). (See the liquid-cloud / WV report.)
5. **CL31 / CL51** rarely find a clean molecular window → calibrate them with the **liquid-cloud method**
   as primary, Rayleigh only opportunistically.
6. **The liquid-cloud calibration is a valid independent cross-check for CL61** — it agrees with
   Rayleigh at ~10 % on well-exposed units and flags the same problem instruments (Zeebrugge).
7. **Flag the outlier units**: Brest Mini-MPL, Zeebrugge & Camborne CL61, Payerne CHM15k (aerosol).
8. **Feed the per-night constants to the Kalman**, which absorbs the sparsity and residual outliers.
9. **Keep rejecting bad nights** — a "no calibration" (flag −2) the Kalman skips beats a spurious
   constant from an aerosol/noise window.
10. **Two network-health metrics are actionable** — far-range noise (electronic background) and
    near-range signal + laser age (low laser) both persist on L2 and can drive an automated alert.
