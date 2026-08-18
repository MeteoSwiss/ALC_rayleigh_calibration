# Multiple-scattering correction of the liquid-cloud calibration

*Consolidated 2026-07-10 (M. Hervo, MeteoSwiss E-PROFILE ALC). Merges: multiple_scattering_check.md.*

Companion code: [`validation/multiple_scattering_eta.py`](../../validation/multiple_scattering_eta.py)
(full PVC implementation, validation, figure and table generator). Applied changes live in
`calibration/cloud/_filters.py` (`_ETA_*` tables + `apply_multiple_scattering_correction`) and
`calibration/cloud/calibration.py` (saturation warning) — the multiple-scattering code was split
out of `calibration.py` into `_filters.py` in the 2026-07 read-once refactor and is re-exported for
back-compat.

**Current state (2026-07-10, authoritative).** The multiple-scattering (MS) factor η(cloud-base) is
computed from the **reference-exact Hogan (2006) Photon Variance–Covariance (PVC) model at the
measured calibration-scene droplet radius a_G = 5.5 µm (≈ 11 µm diameter; α = 10 km⁻¹)**. **Each
instrument type has its OWN table** — including **CL61, which now has a dedicated table and no
longer borrows the CL51's**. This replaces the pre-2026-07 legacy Hopkin ladder (a_G ≈ 8 µm, flat
range dependence, CL61 borrowing CL51). At low cloud base η ≈ 0.95 (vs the legacy ≈ 0.83), i.e.
~9 % lower coefficient C_L for low clouds — the bulk of calibration scenes. The document below
records the deep check and the four independent validation campaigns that led to this state; the
legacy tables are retained only as clearly-labelled history.

## Contents

- [1. What was checked, and the current tables](#1-what-was-checked-and-the-current-tables)
- [2. Method — Hogan (2006) PVC model](#2-method--the-hogan-2006-pvc-model)
- [3. Validation vs the legacy CL31/CL51 tables + reference `multiscatter`](#3-validation-vs-the-legacy-cl31cl51-tables--reference-multiscatter)
- [4. η for all instruments (current PVC tables)](#4-η-for-all-instruments-current-pvc-tables)
- [5. Saturation warning for photon-counting instruments](#5-saturation-warning-for-photon-counting-instruments)
- [6. Impact on real CL31/CL51/CL61 data](#6-impact-on-real-cl31cl51cl61-data)
- [7. 2025–2026 confirmation on CL61 — slope, seasonality, Rayleigh agreement](#7-20252026-confirmation-on-cl61--slope-seasonality-rayleigh-agreement)
- [8. Co-located multi-FOV + Cloudnet-microphysics validation](#8-co-located-multi-fov--cloudnet-microphysics-validation)
- [9. CHM15k β intercomparison — absolute cross-instrument test](#9-chm15k-β-intercomparison--absolute-cross-instrument-test)
- [10. Droplet-size assumption vs Cloudnet climatology + literature](#10-droplet-size-assumption-vs-cloudnet-climatology--literature)
- [11. Assumptions of the liquid-cloud calibration — and how well they hold](#11-assumptions-of-the-liquid-cloud-calibration--and-how-well-they-hold)
- [12. Changes applied](#12-changes-applied)

## 1. What was checked, and the current tables

The O'Connor liquid-cloud calibration relies on the saturated integrated attenuated backscatter of a
fully-attenuating liquid cloud, **B = 1/(2 η S)**, where η ≤ 1 absorbs the multiple-scattering
enhancement of the apparent backscatter. The code corrects the signal gate-wise with η(range)
tables (`beta *= eta(range)`). Three questions were investigated:

1. Do the operational **CL31/CL51** tables fit O'Connor et al. (2004) and Hopkin et al. (2019)?
2. What η should the **CHM15k** (and CL61) use — and what do Le & O'Connor (2026) and
   Looschelders et al. (2025) imply?
3. What η applies to the **Mini-MPL / MPL**, and should the cloud method run on the
   photon-counting instruments at all?

**A code bug was found and fixed on the way** (review finding, 2026-07-05): the `otherwise`
branch of `apply_multiple_scattering_correction` interpolated the CL51 table endpoints for every
non-Vaisala instrument, while the MATLAB reference applies **no correction** (`ones()`) for them.
CHM15k / Mini-MPL cloud calibrations were therefore multiplied by a spurious 0.76–0.83
(≈ 15–20 % bias vs the MATLAB behaviour).

**The tables the code applies today** (verbatim from `calibration/cloud/_filters.py`,
`_ETA_*`; PVC at a_G = 5.5 µm, α = 10 km⁻¹; η at three representative cloud-base ranges):

| instrument | η(0.25 km) | η(1.125 km) | η(2.375 km) | table used |
|---|---|---|---|---|
| CL31 | 0.934 | 0.823 | 0.738 | `_ETA_CL31` (own PVC, 0.83 mrad FOV → strongest) |
| CL51 | 0.953 | 0.862 | 0.783 | `_ETA_CL51` (own PVC) |
| CL61 | 0.954 | 0.865 | 0.787 | `_ETA_CL61` (own PVC — **no longer borrows CL51**) |
| CHM15k | 0.981 | 0.939 | 0.893 | `_ETA_CHM15K` (own PVC) |
| Mini-MPL | 0.982 | 0.942 | 0.897 | `_ETA_MINIMPL` (own PVC) |
| MPL | 0.988 | 0.971 | 0.945 | `_ETA_MPL` (own PVC) |

The narrower the receiver FOV, the smaller the footprint at the cloud, the less forward-scattered
light is retained, and the closer η is to 1 — hence CL31 (widest, 0.83 mrad) gets the strongest
correction and MPL (narrowest) the weakest. CL61's marginally larger divergence (0.28 vs 0.21 mrad)
makes its correction fractionally weaker than the identically-FOV'd CL51 (η differs by < 0.01
everywhere). The CHM15k/Mini-MPL/MPL tables extend to 3.75 km; the Vaisala tables to 2.375 km.
Unknown instrument types get **no correction** (restores the MATLAB `otherwise` behaviour).

**Superseded 2026-07-06:** the pre-2026-07 code applied the legacy Hopkin single-ladder (a_G ≈ 8 µm,
CL31 ≈ CL51, CL61 := CL51 "as approximation") which is flat in range and physically implausible; the
values and the reasoning that retired it are kept below (§3–§4, §6, §10) as history.

## 2. Method — the Hogan (2006) PVC model

η was recomputed from first principles with the Photon Variance–Covariance model of
Hogan (2006, *Appl. Opt.* 45, 5984) — the same fast MS model Hopkin et al. (2019) used:

- half of the droplet extinction is diffracted into a forward Gaussian lobe of 1/e half-width
  **Θ = λ/(π a_G)** (a_G = equivalent-area droplet radius; Θ ≈ 26 mrad for a_G = 11 µm at 910 nm
  — much wider than any ceilometer FOV, so what matters is the FOV **footprint** at the cloud);
- the two-way problem folds into one-way transport in an equivalent medium (forward-lobe
  scattering rate = full extinction σ; return journey in vacuum — reciprocity);
- three photon populations per gate (unscattered / singly / multiply forward-scattered), each
  tracked by relative energy, lateral second moment, angular second moment and covariance
  (Hogan Eqs. 14–26); a telescope of half-angle FOV ρ_t accepts backscattered photons with
  lateral offset s ≤ ρ_t R, i.e. an accepted fraction 1 − exp(−(ρ_t R)²/s²);
- η(cbh) = [2 ∫ σ e^(−2τ) E dz]⁻¹ for a deep homogeneous cloud based at range cbh.

The implementation is a **line-for-line port of the reference** `multiscatter` 1.2.11
`small_angle.c` ('original' algorithm), including Eloranta's exact per-origin double-scattering
integration and within-gate multiple scattering, and **matches the reference binary to
< 1e-4 in η at every tested point** (see §3b).

**Instrument optics used** (half-angles; sources in the script header):

| instrument | receiver FOV ρ_t | divergence ρ_l | λ | source |
|---|---|---|---|---|
| CL31 | 0.83 mrad | ±0.4×±0.7 → 0.57 mrad | 910 nm | Wiegner et al. 2014 Tab. 1 / Vaisala |
| CL51 | 0.56 mrad | ±0.15×±0.25 → 0.21 mrad | 910 nm | idem |
| CL61 | **0.56 mrad** (identical to CL51) | ±0.2×±0.35 → 0.28 mrad | 910.55 nm | CL61 User Guide M212475EN (FOV); Le & O'Connor 2026 Tab. 1 (divergence) |
| CHM15k | 0.23 mrad | 0.15 mrad | 1064 nm | Wiegner et al. 2014 Tab. 1 |
| Mini-MPL | 0.11 mrad (FOV 220 µrad) | ~0.055 mrad | 532 nm | MiniMPL-532 ops manual |
| MPL | ~0.05 mrad | ~0.025 mrad | 532 nm | Campbell et al. 2002 |

## 3. Validation vs the legacy CL31/CL51 tables + reference `multiscatter`

### 3a. Fit to the legacy operational tables

The free microphysics parameters (droplet size, in-cloud extinction) were fitted jointly to the
legacy operational tables: best agreement at **a_G = 8 µm (16 µm diameter — inside the 8–20 µm band
Le & O'Connor use) and α = 10 km⁻¹**, giving rms deviations of **5.1 % (CL31 and CL51)** in η.

![PVC model vs legacy operational CL31/CL51 tables](figs_multiple_scattering/fig_eta_validation_cl31_cl51.png)

Verdict on question 1 — **the legacy tables are consistent with the literature envelope but not with
a per-instrument small-angle calculation**:

- their 0.76–0.83 range sits squarely in Hopkin's stated envelope ("η is usually between 0.7 and
  0.85 for 905–1064 nm in liquid water clouds"), and brackets O'Connor's constant η = 0.7;
- Hopkin et al. (2019) never published numeric tables or the microphysical configuration of
  their multiscatter runs, so an exact reproduction is not possible. The PVC recomputation with
  a single droplet size brackets the legacy values but has a **steeper range dependence**
  (model 0.95→0.67 vs table 0.83→0.76 over 0.25→2.4 km). Plausible causes: a droplet-size
  *distribution* (mixes lobe widths and flattens the slope), a different B-integration depth, or
  a different in-cloud extinction profile in the original runs.

### 3b. Cross-check against Hogan's reference `multiscatter` code (2026-07-06)

To rule out an implementation error, Hogan's own `multiscatter` 1.2.11 (compiled in WSL,
`-algorithms original none` = the exact Hogan 2006 PVC algorithm, no wide-angle) was run on the
same homogeneous-cloud profiles:

![Reference multiscatter vs this implementation vs legacy tables](figs_multiple_scattering/fig_eta_multiscatter_reference.png)

Findings:

1. **This repo's PVC implementation matches the reference exactly** (max |Δη| < 1e-4 over
   both instruments and all cloud-base ranges): after the first comparison showed a 0.03–0.06
   offset, Eloranta's exact per-origin double-scattering integration and the within-gate
   multiple scattering were ported line-for-line from `small_angle.c`, closing the gap
   completely. The Gaussian-moment treatment is retained only for triple-and-higher orders,
   exactly as in the reference.
2. **The reference code confirms the two features the legacy tables lack**: a strong range
   dependence (CL31 0.88 → 0.65 over 0.25→2.4 km) and a clear CL31–CL51 separation (0.03–0.04,
   from the 0.83 vs 0.56 mrad FOV). The legacy tables are nearly flat (0.83 → 0.76) and
   nearly identical for the two instruments (< 0.005 apart) — physically implausible for a
   per-instrument small-angle MS calculation.
3. **No single configuration reproduces the flat tables**: scanning droplet radius 5–25 µm and
   extinction 10–60 km⁻¹ in the reference code always yields a 2–4× steeper range dependence than
   the tables. A config matching the table at 0.25 km (a=5 µm, α=60 km⁻¹ → 0.822) drops to 0.58
   at 2.4 km instead of 0.76.
4. **Provenance of the legacy ladder — RESOLVED (2026-07-06) with E. Hopkin's original code**
   (`Emma.zip`, `VaisCeil_DayCalibration_EProf_EH.py::scatter_correct_Vais`). Her E-PROFILE Python
   applies **one single η ladder to ALL Vaisala ceilometers** — and its eight values (0.82881,
   0.82445, 0.81752, 0.81021, 0.80241, 0.79356, 0.78595, 0.77877 for 0.25-km bins up to 1.875 km)
   are **digit-for-digit the legacy "CL51" table**. Her docstring cites the source: *"Apply multiple
   scattering correction — source: http://www.met.reading.ac.uk/~swr99ejo/lidar_calibration/"*
   (Ewan O'Connor's Cloudnet calibration page), with the comment *"Note these values are
   instrument dependent"* — acknowledging the limitation. So the chain is:
   **O'Connor's page (one Hogan-code run, instrument/config not stated — the page's "beam
   divergence of 1 to 1.5 mrad" wording suggests a CT25K/CT75K-era configuration, matching
   neither the CL31 nor CL51 optics) → Hopkin's single all-Vaisala ladder (step-wise per 0.25-km
   bin, NOTHING above 2 km — clouds higher than 1.875 km got no correction in her code) →
   the E-PROFILE MATLAB (ladder relabeled "CL51"; a near-identical "CL31" variant ≤0.4 % apart;
   two rows APPENDED at 2.125/2.375 km; CL61 := CL51 "as approximation"; linear interpolation
   instead of steps) → this port.** The flat range-dependence and the CL31≈CL51 identity are
   inherited from that single legacy ladder — exactly what the physics comparison predicted.
   The thesis confirms the method (per-gate Hogan-code η, B·η target 0.025 sr⁻¹, S=18.8 sr) and
   the η→0.5 full-retention limit (citing Rogers et al. 1997).
5. **Where it matters**: at the dominant stratocumulus cloud bases (0.8–1.6 km) legacy table and
   physics agree within ±0.05 (≤ 6 % on C_L). The disagreement is at the extremes: below ~0.5 km the
   table likely over-corrects (η_true ≈ 0.88–0.95 vs 0.83), above ~2 km it under-corrects
   (η_true ≈ 0.65–0.70 vs 0.76) — an opposite-sign, cbh-dependent bias of up to ~10–15 % on the
   per-profile coefficient at the tails. The data-driven checks in §6–§9 confirm the PVC
   range-dependence is closer to reality.

### 3c. Literature review — is Hogan (2006) still the right model? (2026-07-06)

The natural worry is that the reference model itself (Hogan 2006, now 20 years old) might be
superseded by something more accurate. A survey of the multiple-scattering (MS) literature says
**no — for this specific regime (narrow-FOV, ground-based, near-IR, liquid water cloud) Hogan's
photon variance-covariance (PVC) method remains the accurate, standard forward model**, and the
newer work either targets a different regime or re-confirms it.

**Hogan's own validation covers exactly our regime.** Hogan (2006, *Appl. Opt.* 45, 5984, §5)
validates the PVC algorithm against a high-order Eloranta calculation — Eloranta's model having
itself been validated against Monte Carlo (Eloranta 1998; refs 5, 9 therein). The agreement is
**"within 4 %"**, and critically: *"Calculations have been performed using the same cloud profile
but with the laser divergence and telescope field of view varied between 0.005 and 50 mrad, and
the agreement between the new algorithm and high-order Eloranta calculations is equally good."*
That 0.005–50 mrad sweep **brackets every ceilometer FOV** (0.23–0.83 mrad), so the ±4 % accuracy
is established precisely where we use it — and our port reproduces this reference to < 1e-4 (§3b).

**Recent MS models target other regimes, not narrow-FOV ground-based near-IR liquid cloud:**

- **Space lidar (CALIOP 532 nm, ATLID 355 nm), wide beam.** [Shcherbakov et al. 2022, AMT 15,
  1729](https://amt.copernicus.org/articles/15/1729/2022/) and [2024, AMT 17,
  3011](https://amt.copernicus.org/articles/17/3011/2024/) build an empirical 3-parameter MS model
  and an Eloranta/McRALI Monte Carlo study — but only for 532/355 nm and satellite/wide-FOV
  geometry, where wide-angle MS dominates. Not applicable to a 0.5-mrad 910-nm ceilometer.
- **Oceanic lidar** (e.g. [Zhai et al. 2019, RS 11, 1870](https://doi.org/10.3390/rs11161870);
  MDPI RS 13, 3677 2021): semi-analytic Monte Carlo for water-penetrating green lidar — different
  medium and geometry.
- **Monte Carlo reference codes** (McRALI, semi-analytic MC, e.g. [Appl. Opt. 64, 6097
  2025](https://opg.optica.org/ao/abstract.cfm?uri=ao-64-21-6097)) are *reference* forward models,
  not faster operational replacements; where they overlap Hogan they agree.

**The recent work that IS in our regime re-confirms the physics we rely on.** The 2024 AMT
multi-layer study finds, for a narrow-FOV (≈1 mrad) ground-based lidar in liquid water cloud, that
the MS contribution is only **~5–30 % within 1 km of the cloud base and < 5 % beyond**, governed by
the *"escape effect"* — forward-scattered photons leave the narrow receiver cone — and that model
accuracy **increases** as the FOV/footprint shrinks (the small-angle-dominated limit, i.e. exactly
the PVC method's design regime). This is the same monotonic FOV dependence (narrower FOV → η closer
to 1, growing correction with penetration depth) that our recomputation shows and that the flat
legacy operational ladder lacks. [Le & O'Connor 2026](https://doi.org/10.5194/egusphere-2025-6331) —
the most recent ceilometer calibration paper — still uses **the same Hogan code** to compute their
CL61 ηS band, confirming it remains the community-standard tool today.

**Conclusion.** No newer model changes the numbers for a 910-nm sub-mrad ceilometer in liquid
water cloud; if anything the recent literature strengthens the case, because (a) Hogan's validated
±4 % accuracy is explicitly demonstrated across the ceilometer FOV range, and (b) the independent
2024 study reproduces the strong range- and FOV-dependence that the legacy operational ladder is
missing. The reference-exact PVC curves in this report are therefore the trustworthy physics; the
flat CL31≈CL51 legacy table is the outlier, for the provenance reason established in §3b.

## 4. η for all instruments (current PVC tables)

![Recomputed eta for all instruments](figs_multiple_scattering/fig_eta_all_instruments.png)

The figure and the §1 table give the reference-exact PVC η now applied in code (band in the figure
= 4–10 µm droplet radius). The physics is monotone in the receiver FOV: the narrower the FOV, the
smaller the footprint at the cloud, the less forward-scattered light is retained, the closer η is to
1. For the photon-counting instruments the MS correction is a **few percent**, not the 17–24 % the
CL51 endpoints wrongly imposed before the fix.

**CL61 has its own table.** The CL61 receiver FOV is **±0.56 mrad** (CL61 User Guide M212475EN,
receiving specifications) — **identical to the CL51**. The only optical difference is the slightly
larger CL61 divergence (0.28 vs 0.21 mrad), which changes η by **< 0.01** everywhere in the PVC
model (e.g. 0.865 vs 0.862 at 1.1 km). The code therefore carries a dedicated `_ETA_CL61` (marginally
weaker than CL51). Le & O'Connor (2026) compute a CL61-specific theoretical ηS band (droplet
diameters 8–20 µm, CL61 optics) — the same approach as here; Looschelders et al. (2025) do not
derive η (they cite Filioglou et al. 2023: WV and MS factors are candidate explanations for
inter-unit calibration differences, and quote Hopkin's < ±5 % coefficient stability for CL31/CHM15k).

**Superseded 2026-07-06 (history).** Before the 2026-07 migration, the code applied the a_G = 8 µm
PVC values with **CL61 := CL51 table** as an approximation of convenience. Those intermediate
8-µm values were:

| instrument | η(0.25 km) | η(1.125 km) | η(2.375 km) | table then applied |
|---|---|---|---|---|
| CL31 | 0.910 | 0.779 | 0.691 | legacy operational table |
| CL51 | 0.934 | 0.823 | 0.736 | legacy operational table |
| CL61 | 0.936 | 0.826 | 0.741 | CL51 table (borrowed) |
| CHM15k | 0.974 | 0.918 | 0.860 | 8-µm PVC table |
| Mini-MPL | 0.976 | 0.921 | 0.865 | 8-µm PVC table |
| MPL | 0.986 | 0.960 | 0.926 | 8-µm PVC table |

The migration to a_G = 5.5 µm (§10) and to a dedicated CL61 table is the current state in §1.

## 5. Saturation warning for photon-counting instruments

Hopkin et al. (2019) did cloud-calibrate Lufft CHM15k operationally (< ±5 % coefficients), but
our own experience is that the photon-counting detectors (CHM15k/CHM8k, Mini-MPL, MPL)
**saturate in the strong liquid-cloud return**, biasing the integrated backscatter low and making
the derived coefficient unreliable. The cloud calibration emits an explicit `UserWarning`
("photon-counting detector saturates in liquid clouds…") whenever it runs on these types; their
reference method remains the Rayleigh calibration. (Enforced by
`tests/test_all_instruments_run.py::test_cloud_runs_for_every_instrument`.)

## 6. Impact on real CL31/CL51/CL61 data

To quantify what actually changes if the legacy tables were swapped for the reference-exact PVC
tables, the operational cloud calibration was run **twice per station-day** (identical inputs, only
the η table swapped) over **Mar–May 2026** on a geographically diverse selection — 3 CL31
(Tenerife 28 °N / Lake Constance 48 °N / N Sweden 69 °N), 3 CL51 (Cyprus / Vienna / Iceland),
4 CL61 (Payerne / Camborne / Lindenberg / N Italy) — giving **343 station-days** with a valid
calibration in both modes (`validation/_eta_impact_experiment.py`). *(This experiment used the
intermediate a_G = 8 µm PVC table vs the legacy ladder; the final a_G = 5.5 µm table further
reduces the CL31 magnitude — see §9–§10.)*

![Per-profile ratio vs cloud-base height](figs_multiple_scattering/fig_eta_impact_profiles.png)

**The change is modest and entirely explained by the η-table difference.** Per-profile the
coefficient ratio C_pvc / C_legacy falls exactly on the analytic expectation η_legacy(z)/η_pvc(z)
(red line above) — a clean validation that nothing else moved. Median ratio **0.99**, P05–P95
**[0.91, 1.09]**; the two tables cross near **1.0–1.5 km**, so the sign of the change depends on
cloud-base height:

| cloud-base band | CL31 | CL51 | CL61 |
|---|---|---|---|
| 0.5–1.0 km | 0.974 | 0.938 | 0.930 |
| 1.0–1.5 km | 1.032 | 0.978 | 0.977 |
| 1.5–2.0 km | 1.067 | 1.000 | 1.003 |
| 2.0–3.0 km | 1.088 | 1.022 | 1.016 |

Day-level medians per type: **CL31 +2.6 %, CL51 −2.6 %, CL61 −3.1 %** (IQRs ±3–5 %). Because the
operational CBH gate is 500–2400 m, the largest tail disagreement (below 0.5 km, where the tables
diverge most) never enters an operational calibration — the low-cloud extreme is gated out.

**The physically decisive test — does the PVC range-dependence remove a spurious dependence of
the calibration constant on cloud-base height?** A correctly-calibrated instrument's constant must
not depend on where the cloud happens to sit, so a *smaller* |correlation(coefficient, CBH)| is
better:

| type | n profiles | legacy \|r\| | PVC \|r\| | verdict |
|---|---|---|---|---|
| CL51 | 84 | 0.314 | 0.205 | **PVC flatter (better)** |
| CL61 | 121 | 0.185 | 0.072 | **PVC flatter (better)** |
| CL31 | 137 | 0.169 | 0.346 | legacy flatter |

So for **CL51 and CL61 the reference-exact PVC correction measurably reduces the constant's
spurious CBH dependence** — direct evidence its range-dependence is closer to reality than the
flat legacy ladder. For **CL31 it goes the other way** with the *fixed 8 µm* table (PVC's steeper
wide-FOV correction over-rotates the trend); the CL31 sample is dominated by one Arctic station
(69 °N, 69/137 profiles) whose cloud regime may confound the test. §8 and §10 show this CL31
over-correction is an artefact of the 8 µm size choice, resolved by the a_G = 5.5 µm table.

![Per-station day-level ratio distributions](figs_multiple_scattering/fig_eta_impact_stations.png)

## 7. 2025–2026 confirmation on CL61 — slope, seasonality, Rayleigh agreement

Extending to **17 months (2025-01..2026-05)** on the three CL61 sites with an independent Rayleigh
reference and a co-located CHM15k (Payerne, Lindenberg, Aosta;
`validation/_ms_seasonal_experiment.py` + `_ms_seasonal_analyze.py`, ~822 cloud-days, ~55 000
in-cloud profiles). Two independent tests, both favouring PVC at **every** site:

**(i) Calibration-independent slope test** — for a fully-attenuating cloud the observed apparent
ratio is ηS_obs(z) = S·η_true(z)/f with f a *constant* calibration error, so
theoretical_ηS(z)/observed must be flat vs height for the correct model. |corr with height|:

| site | n profiles | legacy | PVC |
|---|---|---|---|
| Payerne | 4 996 | 0.242 | **0.057** |
| Lindenberg | 37 850 | 0.169 | **0.107** |
| Aosta | 12 773 | 0.552 | **0.388** |

![Slope discriminator](figs_multiple_scattering/fig_ms_slope_discriminator.png)

The legacy line is nearly *vertical* (it predicts almost no height dependence); the data clearly
decrease with height and follow the PVC curve. This is model-independent (the calibration level is
divided out) — the legacy flat ladder is falsified by the data's real height slope.

**(ii) Agreement with the independent CL61 Rayleigh constant** — PVC brings the cloud C_L closer to
the Rayleigh reference at all three sites (cloud C_L / Rayleigh C_L, closer to 1 = better):

| site | legacy/ray | PVC/ray |
|---|---|---|
| Payerne | 0.90 | **0.93** |
| Lindenberg | 0.96 | **0.98** |
| Aosta | 0.84 | **0.87** |

![Seasonal + agreement](figs_multiple_scattering/fig_ms_seasonal.png)

**Seasonality:** the legacy and PVC monthly-C_L curves are near-identical in shape at each site — so
**the MS-table choice does not create or remove a seasonal cycle**. The residual month-to-month
excursions (e.g. an April spike at Lindenberg/Aosta) appear in both methods and are sampling /
cloud-regime artefacts, not MS-driven; a clean seasonal assessment needs Cloudnet drizzle/high-cloud
screening (§8). PVC (at 8 µm here) still under-shoots the Rayleigh level by 7–16 %, so a residual
(droplet size, or a genuine cloud-vs-Rayleigh offset) remained — the Cloudnet-microphysics test (§8)
and the a_G = 5.5 µm migration (§9–§10) target it.

## 8. Co-located multi-FOV + Cloudnet-microphysics validation

Using co-located different-FOV ceilometers at Payerne (CL61+CL31+CHM15k), Lindenberg
(CL61+CHM15k) and Palaiseau/SIRTA (CL31+CHM15k), each with a Cloudnet `der` (droplet effective
radius + number concentration) product, 2025-2026, ~220 000 in-cloud profiles
(`validation/_ms_fov_experiment.py`, `_ms_fov_analyze.py`, `_ms_cloudnet_tests.py`). Lindenberg
CL61 excluded before 2025-05 (a ~2× hardware step-change). CHM15k is 1064 nm / narrowest FOV.

**Test 1 — co-located FOV.** Confounded by an instrument effect: the CHM15k (photon-counting)
**saturates in liquid clouds**, giving the *largest* apparent-ηS swing with height (2.5×→0.75×) —
backwards for MS (narrowest FOV should be flattest) — so it is not a clean MS probe. The clean
same-detector, same-wavelength pair (Payerne CL31/CL61) has too small a FOV contrast (0.83 vs
0.56 mrad) to strongly separate the models (both predict a near-flat ratio; data consistent with
PVC). Individual slopes are FOV-ordered (CL31 wider → steeper than CL61), as PVC predicts.

**C_L vs altitude (operational check).** Does the *calibrated* constant still depend on cloud
altitude? |corr(C_L, altitude)|, legacy → PVC: **Payerne CL61 0.20→0.01, Lindenberg CL61
0.53→0.33** (PVC removes the altitude dependence). CL31 with the *fixed* a_G = 8 µm table slightly
over-corrected (Payerne 0.19→0.26, Palaiseau 0.16→0.24). CHM15k stays ~0.90 under both (saturation
— unfixable by any MS model).

**Test 2 — first-principles closure with Cloudnet-MEASURED droplet radius + extinction (no
fitting).** Residual |corr(theoretical/observed ηS, height)|, legacy → PVC(measured):

| stream | n | legacy | PVC (measured) |
|---|---|---|---|
| Payerne CL61 | 4 234 | 0.198 | **0.065** |
| Lindenberg CL61 | 24 998 | 0.533 | **0.234** |
| Palaiseau CL31 | 21 720 | 0.218 | **0.025** |
| Payerne CL31 | 31 680 | 0.175 | 0.182 (tie) |
| CHM15k (all) | — | 0.86–0.91 | 0.70–0.82 (saturation) |

**With measured microphysics PVC wins decisively for CL61 and Palaiseau CL31** — the earlier CL31
"over-correction" was an artefact of the *fixed* droplet-size assumption, not the physics
(Palaiseau CL31 → 0.025 is near-perfect closure). Payerne CL31 is a tie — a unit-specific issue
(that CL31 has its own near-range/digitiser artefacts; cf. `network-offset-models`). CHM15k cannot
be cloud-calibrated (saturation).

**Test 3 — droplet-size dependence.** Apparent ηS **decreases with the Cloudnet-measured droplet
effective radius** at every site (slopes −0.001…−0.04 /µm), as PVC/physics predicts (larger
droplets → narrower forward lobe → more retained MS); the legacy table has **no** droplet
dependence. Present but modest/noisy.

![Test 2 first-principles closure](figs_multiple_scattering/fig_test2_cloudnet.png)

**Verdict.** The Vaisala-APD ceilometers (CL61, and CL51 by identical optics) are decisively
better served by PVC — first-principles closure with measured microphysics reduces the residual
altitude trend 3–4× and flattens the calibrated C_L. CL31 also follows PVC once the droplet size
is right; the fixed-table over-correction is a parameter choice, not a model error. CHM15k must
use Rayleigh (cloud saturation). This is the multi-instrument, measured-microphysics,
first-principles confirmation of the reference-exact PVC model over the flat legacy ladder.

Additional figures for this section: `fig_test1_fov_slopes.png`, `fig_test1_fov_ratio.png`,
`fig_cl_vs_altitude.png`, `fig_test3_dropletsize.png`.

## 9. CHM15k β intercomparison — absolute cross-instrument test

CHM15k is not cloud-calibrated (it saturates); it is the independent Rayleigh-calibrated 1064 nm
reference. For each site the co-located Vaisala is cloud-calibrated (legacy vs PVC C_L), its β
gets the 910 nm WV correction and the 910→1064 nm normalisation, is screened to clear air, and
compared to the CHM15k β in the 0.8–3.5 km aerosol/molecular band (reuses the paper pipeline in
`validation/paper/intercompare.py`; driver `validation/_ms_beta_intercompare.py`). Median relative
bias vs CHM15k (lower |·| = better absolute calibration):

| stream | legacy | PVC 8 µm | better | n |
|---|---|---|---|---|
| Payerne CL61 | +21.1% | **+18.9%** | PVC | 108k |
| Lindenberg CL61 | +28.5% | **+23.9%** | PVC | 773k |
| Payerne CL31 | **+74.5%** | +80.5% | legacy | 130k |
| Palaiseau CL31 | **+9.5%** | +15.0% | legacy | 142k |

![β vs CHM15k](figs_multiple_scattering/fig_beta_vs_chm15k.png)

**Two findings.** (1) **PVC (8 µm) improves the absolute CHM15k agreement for CL61 at both sites**,
but **over-corrects CL31 at both sites** — consistent split with the C_L-vs-altitude test.
Reconciling with §8 Test 2 (Palaiseau CL31 altitude-flatness 0.025 with *measured* microphysics):
PVC's η(z) *shape* is right for CL31, but the *fixed*-table PVC magnitude (a_G = 8 µm) is too strong
for CL31's wide FOV and over-corrects the absolute level. (2) The **absolute biases are large** (CL61
20–28%, CL31 +9…+80%) and dominated by *instrumental* offsets — the CL61 near-range electronic
undershoot ([[network-offset-models]]), CL31 unit-specific scale (Payerne CL31 +74% is a bad unit;
Palaiseau CL31 +9.5% is fine), wavelength/WV residuals, and the intrinsic cloud-vs-Rayleigh offset.
The MS-table choice moves them only ~2–6 %: **multiple scattering is a real second-order term, not
the dominant one, in cross-instrument β agreement.**

**The calibration-scene droplet size a_G = 5.5 µm wins at EVERY stream (this is the size now in
code).** Re-running the β-vs-CHM15k with a third table at the calibration-scene droplet radius (§10):

| stream | legacy | PVC 8 µm | **PVC 5.5 µm (current)** |
|---|---|---|---|
| Payerne CL61 | +21.1% | +18.9% | **+12.6%** |
| Lindenberg CL61 | +28.5% | +23.9% | **+16.8%** |
| Payerne CL31 | +74.5% | +80.5% | **+69.8%** |
| Palaiseau CL31 | +9.5% | +15.0% | **+7.6%** |

a_G = 5.5 µm is the **best of the three at all four streams** — it improves CL61 further AND
**fixes the CL31 over-correction** (below legacy at both CL31 sites). (The residual CL31 absolute
bias is instrumental — Payerne CL31's +70% is the bad unit / Kotthaus near-range artefacts, not MS.)

**Consolidated verdict (all tests §6–§10) — implemented in code:**
- **All 910 nm Vaisala (CL61, CL51, CL31) use the PVC table at the calibration-scene droplet size
  a_G ≈ 5.5 µm** (the reference-exact Hogan physics + the measured droplet climatology). It is
  validated on every independent test and best against the CHM15k reference for both instrument
  types. The earlier CL31 caveat is resolved: it was the *fitted* a_G = 8 µm, not the physics.
- **CHM15k → Rayleigh only** (cloud saturation); ideal β reference.

## 10. Droplet-size assumption vs Cloudnet climatology + literature

**What a_G is.** a_G is the *equivalent-area* droplet **radius** (√⟨r²⟩) that sets the forward
diffraction-lobe width θ = λ/(π a_G) in Hogan's model. It is a *radius*: a_G = 8 µm = **16 µm
diameter**; a_G = 5.5 µm = **11 µm diameter**. The original a_G = 8 µm was **fitted** (§3) to
reproduce the legacy tables — not measured.

**Literature.** All published ceilometer cloud-calibration methods use a **fixed** climatological
droplet size (or band) fed to Hogan's code with the instrument optics: O'Connor et al. 2004
(η ≈ 0.7–1, single representative re), Hopkin et al. 2019 (per-gate Hogan η ≈ 0.7–0.85, re not
published), Le & O'Connor 2026 (a **band** of 8–20 µm effective *diameter* = 4–10 µm effective
radius). **None use per-cloud measured microphysics** — so our §8 Test 2 (per-profile Cloudnet re)
is more rigorous than the literature. (Kotthaus et al. 2016 is about CL31 firmware/near-range
*signal artefacts*, not MS — but it explains the CL31 instrumental noise behind the +74 % Payerne
CL31 offset in §9.)

**Cloudnet climatology (8.7 M pixels/site, 2025-2026), THREE samples** — the distinction is the
whole point:

| sample | Payerne | Palaiseau | Lindenberg |
|---|---|---|---|
| all cloud pixels | 10.5 | 8.9 | 7.4 µm |
| cloud base (all clouds) | 9.8 | 8.1 | 6.6 µm |
| **valid-calibration scenes** | **5.8** | **5.6** | **5.4 µm** |

![Droplet climatology](figs_multiple_scattering/fig_droplet_climatology.png)

**The selection effect is decisive.** The clouds the O'Connor method *actually calibrates on*
(fully-attenuating, drizzle-free, homogeneous warm Sc — the S = 18.8 regime) have a droplet
effective radius of **~5.4–5.8 µm (≈ 11 µm diameter), remarkably uniform across all three sites**
and a *narrow* distribution (gold curve), because drizzle (large droplets) and heterogeneous
clouds are explicitly rejected. This is far smaller than the general cloud population (7–10 µm,
broad, site-dependent). **And it is the sample the MS correction is applied to** — so it is the
correct droplet size for the table.

**Conclusion (implemented).** The right MS-correction droplet size is the **calibration-scene value,
a_G ≈ 5–5.5 µm (~11 µm diameter)** — at the lower part of Le's 8–20 µm-diameter band and consistent
with continental-Sc literature (Miles 2000; Frisch 2002, ~5–9 µm). The earlier fitted **a_G = 8 µm
was too large for the calibration clouds** (it was tuned to the suspect legacy tables, not the
microphysics), which over-estimated MS and over-corrected — most for the wide-FOV CL31 (FOV²·a²).
An intermediate draft that argued 8 µm was fine had been misled by the *full* climatology, which
includes the many clouds we never calibrate on; extracting the calibration-scene subset (this
analysis) resolves it. The PVC tables were re-derived at **a_G = 5.5 µm** and are the ones in code
(§1); this removes the CL31 over-correction (the CL31 residual is *also* partly instrumental —
Kotthaus — so it does not vanish entirely).

## 11. Assumptions of the liquid-cloud calibration — and how well they hold

The O'Connor/Hopkin method rests on one identity: for a **fully-attenuating liquid water cloud**,
the integrated attenuated backscatter **B = ∫β dz = 1/(2·η·S)**. Every symbol must be true and
known. Assessment (● well-founded · ◐ real, quantified, partially fixable · ○ violated for some
instruments):

**A. Target physics**
- ● **S = 18.8 sr constant for liquid droplets** (±0.8 sr ≈ 4% floor, O'Connor 2004). Breaks for
  *drizzle* (larger drops, lower S) → depends entirely on the drizzle screen firing; wavelength
  independence 905/910/1064 nm assumed ("very similar", not zero).
- ◐ **η(z) known** — the weakest physics link (this whole report). 0.5–1, depends on droplet size /
  FOV / λ / height; the legacy table is flat and wrong, PVC is right in shape, and a single fixed
  droplet size (now a_G = 5.5 µm) is an approximation for a distribution. 15–25% at cloud-height
  extremes.
- ◐ **Cloud fully attenuating** — enforced only via the ×20 peak-sharpness test; leaky thin/broken
  Sc that look opaque bias C low, hard to catch per-profile.
- ● **Warm liquid only** (T > −20 °C, drizzle-free) — proxies; supercooled mixed-phase can slip through.

**B. Instrument (is ∫β unbiased?)**
- ○ **Detector linearity** — FALSE for CHM15k (photon-counting *saturates* in cloud → uncalibratable);
  Vaisala APDs far more linear but not perfect.
- ○ **No electronic/dark-signal bias in 100–2400 m** — FALSE: CL61 near-range undershoot, CL31
  firmware/near-range artefacts (Kotthaus 2016) bias ∫β altitude-dependently (worst low), seen as
  large offsets (+70% on a bad CL31 unit) that swamp the MS term in absolute β.
- ◐ **Constant absorbs the true window/blower attenuation** — only if window transmission is stable
  between calibration and application (cleaning cycles break it).

**C. Transfer & stability**
- ○ **Constant is slowly-drifting (Kalman)** — FALSE across hardware/firmware steps (Lindenberg CL61
  halved at 2025-05); window cleanings/laser aging → jumps not drifts.
- ◐ **Calibration-scene clouds represent where C is applied** — the method selects a *biased subset*
  (small-droplet drizzle-free warm Sc, ~5.5 µm vs 7–10 µm general). C should be a pure instrument
  gain that transfers, but any condition-dependent term (η, saturation, near-range bias) makes a
  clean-Sc C off elsewhere — the 7–16% CL61 cloud-vs-Rayleigh gap is this.
- ◐ **Horizontal homogeneity** (temporal-consistency filter) — rejects variable scenes at a yield cost.

**Bottom line.** The method is a *relative* technique only as good as its weakest assumption for a
given instrument. For clean Vaisala APDs on warm Sc it is ~10% accurate and the MS refinement
matters (PVC at a_G ≈ 5.5 µm); for a saturating or electronically-biased instrument no η tuning
rescues it — an independent reference (Rayleigh, or a clean co-located instrument) is what exposes
that. Ranked by residual impact: **instrument linearity/near-range bias (largest, uncorrectable) >
constant stability > multiple scattering (correctable, done here) > S/drizzle floor (~4%)**.

## 12. Changes applied

- `calibration/cloud/_filters.py` (multiple-scattering code, split out of `calibration.py` in the
  2026-07 read-once refactor; re-exported from `calibration/cloud/calibration.py` for back-compat):
  - `_ETA_CL31`, `_ETA_CL51`, `_ETA_CL61`, `_ETA_CHM15K`, `_ETA_MINIMPL`, `_ETA_MPL` — reference-exact
    PVC tables at **a_G = 5.5 µm** (each instrument type its own table; CL61 dedicated);
  - `apply_multiple_scattering_correction`: explicit per-instrument table map, unknown types →
    **no correction** (restores the MATLAB `otherwise` behaviour).
- `calibration/cloud/calibration.py`: saturation `UserWarning` for CHM15k/CHM8k/Mini-MPL/MPL cloud
  calibrations.
- `tests/test_all_instruments_run.py`: asserts the saturation warning fires exactly for the
  photon-counting types.
- `validation/multiple_scattering_eta.py`: full PVC implementation, validation, figures,
  table generator (rerun to regenerate everything under `figs_multiple_scattering/`).

**Impact on products.** Moving the CL31/CL51/CL61 tables from the flat legacy ladder to the
reference-exact PVC tables at a_G = 5.5 µm is a documented step-change: ≈ +3 %/−3 %/−3 % day-level
median for CL31/CL51/CL61 at the 8 µm stage (§6), with the 5.5 µm size further reducing the CL31
magnitude and improving the absolute CL61 agreement (§9). Because the CBH gate is 500–2400 m, the
divergent low-cloud tail is gated away. CHM15k and Mini-MPL cloud coefficients had also carried a
spurious 0.76–0.83 CL51 factor (now removed, replaced by their own few-percent PVC tables) — these
are non-primary products, flagged unreliable via the saturation warning (§5); their reference is the
Rayleigh calibration.
