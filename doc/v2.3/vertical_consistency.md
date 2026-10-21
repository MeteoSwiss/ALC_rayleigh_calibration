# E-PROFILE ALC calibration v2.3 — a vertical-consistency test for the Rayleigh fallback tier

*Topic 2 of the v2.3 development (branch `v2.3-integration`). Method `molecular_method = "eprof_v2.3"`
(calibration version code 230), not the default. Code: `calibration/rayleigh/molecular_methods.py`
(`_binned_mean_noise`, `_vertical_consistency`, preset `eprof_v2.3`) and
`calibration/rayleigh/calibration.py` (`_window_search_setup`). Study: report 19, question q19;
scripts `validation/chm15k_stability/q19_*.py`, outputs `C:/DATA/Projects/202610_CHM15k_stability/q19/`.*

## Summary

* **Problem.** eprof_v2.2 added a noise-aware fallback tier ("tier 2") that rescues nights the strict v2
  gates reject: 14 % of the CHM15k nights, 32 % in June. These nights are aerosol-loaded. Their constants
  sit **+6.9 % (median; mean +10.8 %)** above the unit's tier-1 level, and they include junk values such as
  Payerne 2026-09-07 (C_L twice the unit median).
* **A bug disables the existing tier-2 gates.** The propagated noise `sigma_signal` is
  √(native profiles per 300-s bin) too large: ×4.2 on the 15-s CHM15k (median of 6633 nights). The
  tier-2 χ² test, both excess floors and the de-biased reference are therefore inert (χ²_red of the
  selected windows 0.09 instead of 1.28).
* **eprof_v2.3** keeps tier 1 exactly and repairs tier 2:
  1. the noise of the night mean is **measured** on the binned stack instead of propagated;
  2. a **vertical-consistency test**: a rescued window must not be brighter, in R(z) = signal/p_mol, than
     any 750-m block of the column above it or of the 1 km below it (by more than 10 % + 2σ), nor sit at
     the end of an aerosol tail (a block of the 1 km below brighter than the window by more than 10 % + 2σ).
     Candidates are tested best-first, so a failed first choice can be replaced by a consistent window.
* **Validation** (26 CHM15k units, all nights 2025-01 … 2026-10, exact recalibration; tier 1 bit-identical
  on 5848/5848 nights): the tier-2 offset falls from **+6.9 % to −0.7 % (mean +10.8 → +0.5 %)**. The
  rejected constants were the biased ones (+12.8 %, mean +18.3 %). Availability goes from 88.7 % to
  **83.5 %** of clear nights (tier 1 alone: 76.1 %): 59 % of v2.2's gain is kept. On verified-clean nights
  pushed into tier 2 by injected noise, 5.9 % are lost.
* **Payerne 2026-09-07** is rejected (flag −2): its window sits at the end of a decaying aerosol tail and no
  consistent alternative exists.
* **Scope.** CHM15k: measured noise and test. CL61: measured noise only (the test would read the CL61's
  instrumental slope of R(z) as aerosol and lower its tier-2 constants by 14 %). CL31, CL51 and Mini-MPL:
  the eprof_v2.2 fallback tier, unchanged.
* **Verification of the release code:** recalibrating 40 tier-2 nights with the integrated `eprof_v2.3`
  reproduces the study's exact runs on **40/40 nights** (same flag; constants equal to 2e-16).
* **Recommendation:** release it as `eprof_v2.3`; flag tier-2 nights of units with a strong characteristic
  gradient (> 4 %/km) in a follow-up; keep Payerne A (abnormal receiver) on tier 1 only.

## 1. What the fallback tier rescues

The strict tier of eprof_v2 (tier 1) rejects a window when its R², relative residual, ratio scatter or
scattering ratio exceed fixed thresholds. These are signal-to-noise proxies, so an ageing or noisy instrument
fails them even on clean nights. eprof_v2.2 therefore added tier 2: when no window passes tier 1, the three
SNR-driven gates are replaced by their noise-relative forms (χ²_red ≤ 2.5, ratio-scatter excess, de-biased
scattering-ratio reference). Tier 2 lifted the CHM15k availability from 76 % to 89 % of clear nights.

The rescued nights, however, are not clean nights that noise had rejected. They are aerosol-loaded:

* **Aerosol only adds backscatter.** At 1064 nm the extinction of a residual layer inside the 2–9 km band is
  negligible (R − 1 = 0.3 at 4 km gives α ≈ 1e-3 km⁻¹), so the cleanest part of the column has the lowest
  R = signal/p_mol. A window brighter than the air above it is contaminated, and its C_L is biased high.
* **Where tier 2 puts its windows.** Typically just above a decaying aerosol tail (≈ +3 … +5 % of residual
  aerosol); on dust or smoke nights inside a layer reaching 7–8 km, where clean air appears only at 8–9 km
  with R 0.35–0.7 of the window.

![Example profiles: a clean night, an accepted tier-2 window, a window in an aerosol tail, and a window brighter than the cleaner air above (published constant +81 %). Black: R(z)/R(window) in 300-m steps; squares: the 750-m test blocks (green passes, red is darker than the window = aerosol, purple is brighter just below = tail); shaded: ±10 %](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_profiles.png)

## 2. The noise bug

`calibration._sigma_on_fit_grid` divides the native first-difference noise by √(binned profiles) only. The
averaging of the 20–25 native profiles inside each 300-s bin is never divided out. The ratio of the
propagated `sigma_signal` to the noise of the night mean measured on the binned stack is:

| type | sampling | ratio (median) |
|---|---|---|
| CHM15k | 12–15 s | 4.2–4.8 |
| CHM15k | 30 s / 60 s | 3.0 / 2.1 (= √(300 s/dt) × 0.95) |
| CL61 | 5–10 s | 1.7 (1.2–2.9; its native noise is not white) |
| CL31 | 15 s | 2.1 |

With the measured noise the selected windows have χ²_red 1.58 (tier 1) and 1.28 (tier 2); with the
propagated one 0.12 and 0.09, so `max_chi2red = 2.5` never binds. This also explains the earlier finding
that the tier-2 χ² tolerance was "inert over 1.5–3.0".

![Propagated vs measured night-mean noise; χ²_red of the selected windows with both estimates; example SNR profile](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_sigma_check.png)

**The fix alone** (measured noise, thresholds unchanged) keeps 91.5 % of the tier-2 nights and lowers their
offset from +6.9 to +4.2 %. It is necessary for the test of §3: with the propagated noise the test detects
aerosol on 3 % of the first candidates instead of 14 %.

eprof_v2.3 measures the noise of the night mean on the masked binned stack (robust first differences of
consecutive 5-min profiles / √2 / √n, 11-gate running median). Where it cannot be measured (fewer than 8
binned profiles, < 40 min), it falls back to the propagation corrected by √(native profiles per bin). The
shared `sigma_signal` is left as it was: the v2.3 uncertainty budget applies its own time-binning correction,
so correcting it at the source would count that factor twice.

## 3. The vertical-consistency test

For a tier-2 candidate window, R(z) is averaged (noise-weighted) over the window and over 750-m blocks:

* above the window, up to the height where the expected 300-m block SNR falls below 3 (at most 9 km);
* below it, down to max(1.5 km, window bottom − 1 km).

The window is **rejected** when

* **aerosol:** some block is darker than the window by more than 10 % + 2σ (the window contains aerosol, or
  the profile is distorted in height); or
* **tail:** a block of the 1 km below is brighter than the window by more than 10 % + 2σ (the window sits at
  the end of an aerosol tail still decaying into it).

A block **above** that is brighter (cirrus, an elevated layer) never rejects. Candidates are taken
best-first by the operational objective; the first consistent one is kept, otherwise the night fails with
"tier-2 windows failed the vertical-consistency test".

| design (exact runs, CHM15k) | tier-2 kept | offset median / mean [%] | false rejection of clean windows (first candidate) |
|---|---|---|---|
| linear d ln R/dz over [window − 1 km, top], 5 %/km + 2σ (first proposal) | 61.3 % | +1.8 / +4.7 | 48 % (boundary-layer top below the window) |
| cleanest block only (no tail test) | 85.4 % | +3.3 / +4.7 | 6 % |
| **final: cleanest block + 1 km below, k = 2** | **58.6 %** | **−0.7 / +0.5** | re-selection, see §4b |
| final, k = 3 | 63.4 % | +0.1 / +1.6 | |

* **Why not a linear gradient:** it is dominated by the high-SNR gates next to the window and is fooled by
  U-shaped profiles (a decaying tail below plus an elevated layer above).
* **Why the 1 km below:** a window just above a still-decaying tail carries ~3–5 % residual aerosol. The same
  holds on tier-1 nights: tier-1 windows failing the tail test sit 2.5 % above those that pass; tier-1 windows
  failing the aerosol test (8.7 % of tier-1 nights) sit +7.7 % above (mean +12.6 %). Tier 1 is not changed
  here, but this is the size of the bias it still carries.
* **Thresholds:** ±10 % is the night-to-night scatter; 2σ applies to every block; 750-m blocks reach a block
  SNR of 10–20 at 7 km.
* **Windows with no measurable block above** (noisy nights whose wide windows reach the SNR limit) are judged
  on the 1 km below. A first design rejected them as "unverifiable" and lost 31.7 % of clean noisy nights in
  the injection test.

## 4. Validation

### 4a Offsets and availability (26 CHM15k units, Payerne A excluded)

Offsets are ln C minus the unit's tier-1 level (a local line through its tier-1 nights within ±45 d).

| series | tier-2 nights | offset median | mean | robust sd | p90 abs | availability (pooled) |
|---|---|---|---|---|---|---|
| tier 1 only (leave-one-out) | – | −0.2 | +1.2 | 10.6 | 19.6 | 76.1 % |
| v2.2 published | 858 | **+6.9** | **+10.8** | 17.9 | 40.2 | 88.7 % |
| measured noise only | 785 | +4.2 | +7.6 | 16.4 | 33.8 | 87.7 % |
| **eprof_v2.3 (k = 2)** | **503** | **−0.7** | **+0.5** | 14.6 | 27.0 | **83.5 %** |
| test with k = 3 | 544 | +0.1 | +1.6 | 14.6 | 27.5 | 84.1 % |

Where the bias goes (the 858 published tier-2 constants):

| outcome | nights | published offset median / mean | v2.3 median / mean |
|---|---|---|---|
| rejected | 423 | **+12.8 / +18.3** | – |
| kept, window re-selected | 393 | +3.8 / +4.3 | −0.2 / +0.7 |
| kept, same window | 42 | −2.0 / −1.5 | −2.0 / −1.5 |

68 nights are newly rescued. By first-candidate verdict the published constants sit at +23.0 % ("aerosol",
68 nights), +4.6 % ("tail", 622) and +1.1 % ("accept", 64).

* **Units with a strong characteristic gradient** (Oslo, Budapest B, Extremadura: +5.6, +5.8, +7.4 %/km,
  §6) keep +7.3 % after the test (75 nights); all other units −1.0 % (403 nights).
* **Scatter and annual cycle:** night-to-night scatter 9.5 → 9.4 %, annual amplitude 19.9 → 19.7 %. The
  tier-2 nights are not the cause of the CHM15k annual cycle (report 19).

![Tier-2 offsets by series and by month, and tier-2 nights kept per month](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_offsets.png)

### 4b Injection self-test (143 verified-clean CHM15k nights, 26 units)

Clean winter nights (fully clear, window passing the test at native noise) are modified on the native L1
profiles, before cloud screening, binning, noise estimation and window search. Aerosol is added as
Δrcs = C·β_a·T_mol² with its two-way transmission (LR 50 sr).

| test | result |
|---|---|
| forced into tier 2 at native noise | 100 % kept, bias 0.0 % (p90 2.9 %) |
| white noise added until tier 2 engages (101 cases): false rejection | v2.2 0 %, measured noise 3.0 %, **v2.3 5.9 %**, k = 3: 5.0 %, first design 31.7 % |
| decaying layer, scale height 1 km, R − 1 = 5–30 % at 4 km: detection | 91–98 % at native noise, 62–69 % at tier-2 noise |
| same, scale height 2 km | 37–93 % at native noise, 15–50 % at tier-2 noise |
| additive baseline B·r², B = −5 % | detected on 47 % |
| additive baseline B = +2 … +20 %, uniform aerosol | not detected (by design / by construction) |

An injected exponential tail never ends inside the search range, so a re-selected window still carries up to
~10 % + 2σ of it: the kept constants lose 20–40 % of their bias for steep tails, ≤ 15 % for slow ones. On real
nights the aerosol usually has a top above which the re-selected window is clean, consistent with the ~0 offset
of §4a. The test is a tolerance test, not an aerosol-free guarantee.

![Injection results](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_injection.png)

### 4c Co-located twins

Same-night twins (Schiphol A–D, Budapest A/B, Măgurele A/B) are not an independent truth here: on the nights
when one unit needs tier 2, its twin's own tier-1 constant is itself +6.2 % above its level (its window passes
the test on only 10 % of those nights). The twins cancel the common annual gain, not the common aerosol.

### 4d Verification of the release code

The study ran its rule through a patched copy of the package. The integrated `eprof_v2.3` (separate method,
v2.2 untouched, noise correction confined to the window search) was checked against those exact runs on 40
CHM15k tier-2 nights (15 kept on the same window, 15 re-selected, 10 rejected): **flags equal 40/40,
constants equal 40/40** (largest relative difference 2.2e-16). Unit tests: `tests/test_vertical_consistency_v23.py`.

## 5. Payerne CHM15k 2026-09-07

The published value is reproduced exactly with the classification screen (flag 0.5, C_L = 1.02e12, window
5.57–6.05 km above the instrument, tier 2). The window sits at the bottom end of a steeply decaying tail (R falls
by ~3.6× from 4 to 5.5 km): the block just below is 1.66× brighter ("tail"), and no block ≥ 1 km above is
measurable. No other tier-2 candidate passes, so eprof_v2.3 gives **flag −2**. The low junk value of 09-06
(−47 %) is rejected too; 09-05 is kept with a re-selected window.

![Payerne 2026-09-07 with and without the classification screen, and 2026-09-05 (range on the Y axis)](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_payerne_20260907.png)

## 6. Instrument scope

| type | measured noise | vertical-consistency test | reason |
|---|---|---|---|
| CHM15k (154 streams) | yes | yes | validated on 26 units (§4) |
| CL61 (13) | yes | **no** | noise alone: offset −2.5 → −2.6 %, scatter 16.1 → 14.1 %, availability 71.6 → 68.0 %. The test: tier-2 constants −14.4 %, because the CL61's R(z) falls 3–4.5 %/km on clean tier-1 nights (Lindenberg 3→8 km: 1.00 → 0.78, while the co-located CHM15k is flat) |
| CL31 (210), CL51 (52), Mini-MPL (5) | no | no | not validated: on the CL31 the measured noise opens tier 2 on ~19 % of the nights (v2.2: 0.7 %), with no tier-1 night to check them against |

`calibrate_rayleigh._window_search_setup` applies this scope; an explicit `molecular_params` entry
(`measured_noise`, `use_vc`) overrides it. The CL61 height dependence of R(z) is itself a finding for the CL61
Rayleigh calibration (water-vapour correction or range normalisation to be checked).

**Payerne A** (CHM15k TUB140016, abnormal receiver): tier 1 alone calibrates 16 % of its clear nights. v2.2
adds another 51 percentage points through tier 2, at −23.5 % below its tier-1 level; v2.3 still adds 13 points,
at −18 %. The test removes the junk rescues but cannot fix this unit: keep it on tier 1 only.

## 7. A unit's characteristic gradient

A unit's height dependence can be measured from the profiles of its clean winter tier-1 nights (d ln R/dz over
2.5–8 km, windows passing the aerosol test; 16–51 nights per unit, standard error ~1 %/km). It matches the
fixed-window recalibrations of report 19 with r = 0.92 (20 units) and recovers Oslo +5.6, Budapest B +5.8,
Extremadura +7.4 %/km. Payerne A comes out at −1.4 %/km on clean nights: its −14.7 %/km "station gradient" in
the availability study was mostly the aerosol bias of high rescued windows.

![Per-unit profile gradient vs the fixed-window recalibration](C:/DATA/Projects/202610_CHM15k_stability/q19/q19_unit_gradient.png)

Use (follow-up, not in v2.3): a two-pass unit flag — tier-2 nights of units with |g| > 4 %/km carry g × Δz.

## 8. Limits

* Blind to vertically uniform aerosol, to layers extending above the SNR limit (~7–9 km) and to positive
  additive baselines (a rising R(z) is deliberately allowed). Slow tails at tier-2 noise are detected 15–50 %.
* The reference is the tier-1 series, itself biased by ~+2 % (8.7 % of tier-1 windows fail the aerosol test at
  +7.7 %). Against verified tier-1 nights the v2.3 tier-2 nights sit at +1.6 %.
* Thresholds (10 %, 2σ, 750 m, 1 km) were chosen on this data set (26 units, 22 months); the trade-off is
  smooth and the injection test does not depend on the choice.
* Nights shorter than 40 min fall back to the corrected propagated noise.

## 9. How to run it

* `options.json`: `"molecular_method": "eprof_v2.3"`; or for a network run `ALC_MOLECULAR_METHOD=eprof_v2.3`.
* The rows carry `calibration_version = 230`; the dashboard labels them "v2.3".
* Thresholds can be overridden through `molecular_params` / `ALC_MOLECULAR_PARAMS`
  (`vc_k`, `vc_eps_tol`, `vc_tol_lo`, `vc_block_m`, `measured_noise`, `use_vc`).
* Every night that tier 1 calibrates is bit-identical to eprof_v2.2.

## 10. Recommendation

1. Release the repaired fallback tier as `eprof_v2.3` (measured noise + test for the CHM15k, measured noise for
   the CL61, v2.2 elsewhere).
2. Follow-ups: the unit-gradient flag of §7; Payerne A on tier 1 only; the CL61 R(z) slope; revisit the CL31
   Rayleigh tier 2 only with an independent reference (cloud calibration).
