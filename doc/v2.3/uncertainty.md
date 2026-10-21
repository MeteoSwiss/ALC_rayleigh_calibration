# E-PROFILE ALC calibration v2.3 — a realistic uncertainty for the Rayleigh lidar constant

*Topic 1 of the v2.3 development (branch `v2.3-dev`). Option `uncertainty_model = "v2.3"`, off by
default. Code: `calibration/rayleigh/uncertainty.py`; validation: `validation/v23_uncertainty/`.*

## Summary

* **Problem.** The published nightly uncertainty is about three times too small: reported 1 σ 2.7 %
  against an observed night-to-night 8.9–11 % on the CHM15k network; 63 % of the nights fall outside
  2 σ, and even same-night twins disagree beyond 2 σ on 32 % of the nights. It is blind to the
  dependence of the constant on the window height, to the molecular model and to most of the noise.
* **v2.3 budget** (option `uncertainty_model = "v2.3"`, the constant itself unchanged): the ensemble
  term plus a **vertical-consistency term** — C_L recomputed in 500 m sub-windows around the window,
  |d ln C_L/dz| × 1 km plus the scatter the line does not explain, **saturated at 15 %** — a
  molecular-model climatology, a residual-aerosol allowance and the photon noise **measured as the
  noise of the window mean**, in quadrature, reported at 2 σ. The three bold items come from an
  independent review (§8), which also re-validated everything on 37 fresh units.
* **Validation** (exact code values). Design units (36 units, 7 736 nights, 2025-01 … 2026-07, 99.0 %
  bit-identical to the published constants): CHM15k time-series pull 4.05 → **1.25** (63 % → 17 %
  outside 2 σ), same-night twins 1.69 → **0.60** (32 % → 3.5 %), CL61 3.13 → 1.08 (budget as first
  shipped). **Fresh units** (27 CHM15k, 8 CL61, 2 Mini-MPL, eprof_v2.3): CHM15k 3.91 → 1.17 as first
  shipped, **1.30 reviewed** (18 % outside 2 σ); CL61 2.85 → 1.08 / 1.20; twins 1.66 → 0.63 / **0.69**.
  The two truths bracket the answer, but loosely: a perfectly calibrated budget would give 1.00–1.46 on
  the time series (which contains real gain changes) and 0.40–0.77 on the twins (where the atmospheric
  error cancels) — about ±30 % on the level.
* **Conditional calibration** (review): as first shipped, the budget over-covered its steepest-gradient
  decile 2.2 × (pull 0.46) and under-covered the quiet bulk; the error does not grow like |g| × D. With
  the saturation the top decile is calibrated out of sample (pull 1.02) and the spread of the decile
  pulls falls by 60 %. The quiet nights stay above 1 on the time-series truth only — on the twin truth
  they are calibrated: that floor is real instrument gain change.
* **Large where it must be:** Messina, Payerne CL61, Vásárosnamény, Payerne CHM15k, Gottfrieding and
  Montsec rank 1–6 of 34 (observed ranks 1–7); the unit ranking improves (Spearman 0.56 → 0.72 design,
  0.40 → 0.58 fresh) but not significantly (95 % CI of the gain includes 0). The Payerne 2026-09-07
  night (window 6.06–6.54 km ASL, C_L twice the unit median) goes from 42 % to 162 % (2 σ; reviewed:
  71 %) and is rejected by the QC; eprof_v2.3 rejects it in the window search anyway. Schleiz (level
  steps ×2) and the Mini-MPL (gain drift within the night) are instrument problems no per-night
  uncertainty covers.
* **QC:** `u > C` fired on 0.08 % of the nights; the relative limit u/C_L > 60 % (2 σ) rejected 7.7 %
  with the budget as first shipped, 3.4 × enriched in nights more than 25 % off but half of them good
  nights; with the Kalman weighting on, rejecting them made the best estimates less consistent. With
  the reviewed budget the limit fires on < 0.5 % of the nights (extreme noise only).
* **Kalman** weighting by the uncertainty (eprof_v2.3): co-located best estimates 13.6 % (rms) and
  13.1 % (95th percentile) more consistent than with equal weights (eprof_v2.2: 5–8 %, marginal).
  **Cost:** 2.3 ms per night (0.2 %); the reviewed noise estimator adds ~0.4 ms.
* **Side finding:** v2.2's propagated noise (`sigma_signal`) is √k too large on the CHM15k (time
  binning missing; k = native profiles per 300-s bin); on the CL61 the native noise is not white and
  √k over-corrects by 1.9 ×. v2.2 is left untouched; the reviewed budget measures the noise instead.
* **Recommendation:** adopt the reviewed budget with `eprof_v2.3` and the Kalman weighting on, keep the
  relative QC only as a safety net, and carry the Kalman standard deviation into L2 instead of the fixed
  25 % (§5, §8.4).

## 1. The problem

The nightly Rayleigh constant is published with `lidar_constant_uncertainty` = 2 × the robust
standard deviation of a 125-member ensemble: aerosol lidar ratio 52 ± 10/± 20 sr × window shift
0/± 100/± 200 m × time sample (all profiles + four random 70 % subsets). Report 19 (q18) compared it
with what the constant actually does from night to night on 147 CHM15k units:

| | reported 1 σ (median) | observed night-to-night 1 σ | pull robust std | nights outside 2 σ |
|---|---|---|---|---|
| CHM15k network (147 units) | 2.7 % | 8.9 % | 3.3 | 51 % (Gaussian: 4.6 %) |
| Payerne A | 6.1 % | 22 % | — | — |

The ensemble is about three times too small, for three reasons:

1. **The window-shift dimension only sees the local gradient.** With a C_L-versus-height gradient
   g (%/km) the ±200 m shifts give a robust σ of only ≈ 0.15 g (2σ ≈ 0.3 g), whereas the molecular
   window moves by ~1 km between nights (median window-centre σ 0.86 km, summer vs winter 1.2 km;
   §3.6). A true lidar constant does not depend on the height of the window. When it does — an
   electronic profile distortion (dark/baseline offset), a faint aerosol layer, a window that is not
   purely molecular — nothing in a single profile can reliably tell these causes apart, and the
   honest answer is a large uncertainty.
2. **The molecular model is not in the budget.** US Standard 1976 instead of the real temperature and
   pressure biases C_L by ~2 % RMS at 2–6 km (−3 % in a southern summer; report 19, q15d).
3. **Photon noise is only seen through 70 % subsets of the same profiles** (they share 70 % of the
   data, so their spread is ≈ 0.65 of the noise), mixed into a median with the other dimensions.

Downstream, the L2 product does not carry the calibration uncertainty at all
(`uncertainties_att_backscatter_0` is a fixed 25 % of |β_att|), and the operational Kalman best
estimate gives every night the same weight.

## 2. Method

### 2.1 The budget

With `uncertainty_model = "v2.3"` the reported constant is **unchanged** (same ensemble median, same
QC except the gate in §2.7) and its uncertainty becomes

  σ² = σ_ens² + σ_vert² + σ_mol² + σ_aer² + σ_noise²,  `lidar_constant_uncertainty` = 2 σ · C_L (k = 2)

with all σ relative (1 σ). The five components are written to the calibration NetCDF as absolute
2-σ values (units of C_L) — `lidar_constant_uncertainty_ensemble`, `_vertical`, `_molecular`,
`_aerosol`, `_noise` — together with the height gradient `lidar_constant_height_gradient`
(d ln C_L/dz, km⁻¹).
The variable attributes state the 2-σ convention. None of this is written when the option is off,
so v2.2 files keep their exact structure.

### 2.2 (a) Ensemble term σ_ens

The existing 125-member ensemble, unchanged (it is the v2.2 `lidar_constant_uncertainty` / 2).

### 2.3 (b) Vertical-consistency term σ_vert

C_L is recomputed with the operational machinery (`_compute_cl_for_perturbation`: Klett from the
ground with the reference in the window, aerosol lidar ratio 52 sr, all profiles) in contiguous
500 m sub-windows that tile the usable molecular range:

* aligned on the operational window, from 1 km under its bottom — **but never under 4 km AGL** — to
  1 km above its top. Under a low window the layer below is the boundary-layer aerosol the window
  search deliberately avoided, and probing it measures that aerosol, not the window (with a 1.5 km
  floor the median gradient is 13.7 %/km instead of 3.4 %/km and the twins are over-covered 3×,
  §3.4). Under a high window it is free troposphere, where a strong gradient is a genuine warning
  (Payerne 2026-09-07, §3.8). 14 % of the CHM15k windows start above 4 km;
* above the window only while the mean per-gate SNR of the sub-window is ≥ 3 (30 m gates);
* each sub-window k gives C_k and its photon-noise σ_k (§2.5).

A straight line is fitted through ln C_k against height, weighted by 1/(σ_k² + 1 %²) (the 1 % floor
keeps the high-SNR low sub-windows from dominating). Its slope g = d ln C_L/dz is the vertical
inconsistency of the night, and

  σ_vert = √[ (|g| · D)² + s_res² ]

where s_res = √(max(χ² − ν, 0)/Σw) is the scatter of the C_k about the line that their noise (and
the 1 % floor) does not explain — a layer or a kink inside the probed range — and
**D = 1 km** is the representative displacement of the molecular window between nights
(§3.6: window-centre σ 0.86 km across nights, 0.78 km within a month, 1.2 km summer vs winter,
2.2 km between the 10th and 90th percentile). If fewer than three sub-windows are usable the width is
halved; with two the gradient is their difference; with none the whole ensemble spread is attributed
to the ±200 m shifts (g ≈ 2σ_ens/0.3 per km), a conservative fallback (on the validation nights the
fallback never fired and the two-sub-window branch 1–5 times in 2 200–2 400 nights).

**Saturation (review).** The term is then saturated, σ_vert → σ_vert / √(1 + (σ_vert / c)²) with
c = 15 % (`vertical_cap_rel = 0.15`). The displacement picture behind D is not what the data show:
the anomaly of a night is uncorrelated with g × (its window centre minus the unit's 31-day median
window centre) (slope +0.01, Spearman −0.001, 1 979 CHM15k nights). |g| is a good indicator of a
worse night, but the error does not grow like |g| × D: across octiles of |g| the observed time-series
scatter rises only from 9 to 17 % (eprof_v2.2; 9 to 12 % with eprof_v2.3) while |g| × D rises from
3 to 30 %. Since the time-series truth is an upper bound of the calibration error, the steepest decile
was over-covered by at least 2.3 ×. The window search keeps the cleanest part of the column, so the
window constant's error does not grow without bound with the gradient measured across the probed
range; the cap (chosen on the design units, validated on fresh units, §8) restores calibration in the
top decile and barely touches the typical night (a 7 % vertical term becomes 6.3 %, a 10 % one 8.3 %).

What the term measures — and what it does not:

* an **electronic distortion** (an additive baseline B·r² on the range-corrected signal) makes C_k
  grow (B > 0) or fall (B < 0) with height roughly as r²/β_mol; a synthetic molecular profile with
  B equal to 5 % of the molecular signal at 5 km gives g = +1.6 %/km and a vertical term 12 × that of
  the clean profile (§3.9);
* **aerosol inside a sub-window** raises its C_k (the reference assumes it is molecular); aerosol
  between two sub-windows is corrected by the Klett step with the assumed lidar ratio;
* **noise** alone makes g scatter by its standard error; `|g|` therefore carries the noise of the
  gradient estimate on average (E[g²] = g_true² + σ_g²) — measured with the odd and even profiles of
  each night, that noise is 1.3–1.5 %/km on the CHM15k, 1 % of E[g²]: |g| is signal;
* the term does **not** see a height-independent bias (a uniform free-tropospheric aerosol, a
  calibration-pulse or gain error) — those are invisible to any single-profile method.

### 2.4 (c) Molecular-model term σ_mol

The ln-bias of C_L from using US-1976 instead of the real T/p was computed per night with CAMS
(0.4°, L137) in each night's own window for 22 618 nights of 127 CHM15k units at 37–55 °N (report 19,
q15d). It is well described by a climatology (fit residual 1.35 %):

  b [%] = −0.421 + 1.516 cos ωd + 0.799 sin ωd + 0.080 (φ−47) − 0.017 (φ−47) cos ωd
          − 0.493 (z_c−4.5) − 0.210 (z_c−4.5) cos ωd + 0.008 h_stn

(d = day of year, ω = 2π/365.25, φ = latitude in °, z_c = window centre in **km above sea level**,
h_stn = station altitude in km; southern stations use the mirrored season; the latitude is clipped
to 35–62 °N). The fit is reproduced to the third decimal by `validation/v23_uncertainty/
u23r_molclim.py`, which also shows that the coefficients were fitted with z_c ASL: until the review
the module evaluated them with z_c above ground, which shifts b by +0.49 % per km of station altitude
(≤ 1 % at the 2-km stations; option `molecular_height_reference`, now `"asl"`). σ_mol = √(b² + 1.35 %²):
the expected bias of the night plus its synoptic scatter. Reading CAMS T/p per night would cost ~1 s
(daily file) to ~3 s (8 GB monthly file) — as much as the whole calibration — while the climatology
already reproduces the realised RMS (1.98 %). With `molecular_source = "cams"` the reference is
already CAMS and only a 0.5 % floor remains.

As an ABSOLUTE uncertainty the term is calibrated: |b| ≤ 2 σ_mol on 97 % of the nights of the fitted
units and on 96–99.7 % of the nights of units outside the fit up to 62 °N (`u23r_molecular_oos.py`).
North of the clip it under-covers: 79 % of the nights of five arctic units (Andøya, Jan Mayen,
Bjørnøya, Hopen, Tórshavn) are within 2 σ, with a mean realised bias 1.0 % above the climatology
(2.7 % at 76.5 °N). The time-series truth cannot test the slowly varying part of b (the 31-day median
removes it): replacing σ_mol by its synoptic part alone (1.35 %) changes the CHM15k pull by 2 %.

### 2.5 (d) Photon-noise term σ_noise

The constant is a window MEAN of RCS/β_mol (the gate-to-gate noise cancels in the Klett
normalisation), i.e. the time mean of the per-profile window means of signal/p_mol. Since the review
(`noise_estimator = "window_mean"`, default) σ_noise is MEASURED as exactly that: the window mean of
signal/p_mol is formed for every binned 300-s profile, its robust first difference between consecutive
profiles / √2 is the noise of one profile, / √n that of the night mean (× 1.2533 for a median
aggregation). The same estimator gives the noise of every sub-window constant of §2.3. Consecutive
differences also carry the atmospheric change over 5 min: an upper bound, like the eprof_v2.3
tier-2 noise (`molecular_methods._binned_mean_noise`). It needs ≥ 8 profiles (40 min); shorter nights
keep the per-gate estimate below.

The pre-review estimate (`noise_estimator = "per_gate"`) combines the per-gate noise that eprof_v2.2
measures on the native profiles, referred to the expected molecular signal ĉ·p_mol, as INDEPENDENT
gates: σ_noise = √(Σ e_i²)/N, e_i = σ_i/(ĉ p_mol,i). Split-night constants show that this is too small:
the scatter of (C_odd − C_even)/2, from the odd and the even 300-s profiles of the same night (the
atmosphere is common to both halves), is 1.6 × σ_noise on the CHM15k, 3.4–4.2 × on the CL61 and
2.8–3.1 × on the Mini-MPL (design and fresh units). Noise common to the gates of a profile (a
per-profile background-subtraction error, range-correlated detector noise) does not average over the
window; the window-mean estimator keeps it and passes the same test: 0.96 (CHM15k), 0.88 (CL61), 0.99
(Mini-MPL) on the fresh units.

**A propagation error in v2.2.** `calibrate_rayleigh` divides the noise of ONE NATIVE profile (12–15 s
on the CHM15k) by √(number of 300-s BINNED profiles), so its `sigma_signal` is too large by √k,
k = native profiles per bin. On the CHM15k (k = 5–30 for 60–10 s sampling) the per-gate σ/√k matches
the per-gate noise measured on the binned stack within 2–8 % (ratio 0.94–0.98, 4 300 nights): the √k
correction is right there. On the CL61 the ratio is 0.25–0.53 (median 0.52): its native noise is not
white and dividing by √k over-corrects by ~1.9 × — one more reason for the window-mean estimator,
which needs no assumption on k. The v2.2 gates are left untouched (bit-identical output), but their χ²_red ≤ 2.5
tolerance is a true χ²_red of ~50 on the CHM15k, and the per-gate SNR quoted in report 19 (q16a/q16c,
e.g. "SNR crosses 1 at ~6 km") is ~4.5 × too pessimistic.

### 2.6 (e) Residual uniform aerosol σ_aer

A vertically uniform aerosol in the window leaves R(z) flat and biases C_L high; no single-profile
test can see it. A fixed allowance covers it: 1 σ = 2 % at 1064 nm, the night-to-night effect of the
CAMS in-window aerosol scaled by its observed transfer to C_L (2.2 %, report 19 q15d), scaled to other
wavelengths as (λ/1064 nm)³ (β_mol ∝ λ⁻⁴, β_aer ∝ ~λ⁻¹): 1.25 % at 910 nm, 0.25 % at 532 nm.

### 2.7 QC: a relative limit replaces `u > C`

Flag −6 ("uncertainty exceeds value") required u > C_L, i.e. a 1-σ uncertainty above 50 %; it fired
on 2 of 2 436 nights (0.08 %). In v2.3 mode it becomes a relative limit on the full budget,
`u / C_L > max_rel_uncertainty` (2 σ), default **0.6** (1 σ = 30 %, about 4 × the median night). The
value is tied to the Kalman weighting: `kalman_best_estimate` clips the night-to-night weight ratio at
4 in σ, so a night beyond ~4 × the median uncertainty is no longer down-weighted — it is rejected
instead of being kept at a capped weight (§3.7).

**Review (2026-10): rejection does not pay.** On the design units the gate (pre-review budget)
rejected 7.6 % of the eprof_v2.2 nights: 33 % of them were more than 25 % off (3.7 × the kept nights),
but 40 % were good nights (< 10 % off), and as a QC score the budget separates bad from good nights
with an AUC of only 0.76 (0.70 with eprof_v2.3). With eprof_v2.3 the gate rejected 5.1 % of the
nights, half of them good, and 79 % of its tier-2 rejections were already removed by the
vertical-consistency test. With the Kalman weighting on, rejecting these nights made the co-located
best estimates LESS consistent than keeping them at their low weight (eprof_v2.3: twin rms −6 %
instead of −17 % against equal weights; eprof_v2.2: +4 % instead of −7 %). With the saturated vertical
term (§2.3) the 0.6 limit fires only on extreme noise or ensemble spread (< 0.5 % of the nights): the
weighting, not the rejection, handles the uncertain nights (§8).

### 2.8 Kalman weighting

`monitoring/kalman.kalman_best_estimate` already accepts `uncertainties=`: it keeps the empirical
measurement-variance scale (day-to-day scatter) and takes only the night-to-night RATIO of the
relative uncertainties (clipped to [0.5, 4] in σ). The network runner now passes them when
`kalman_uncertainty_weighting = 1` (options.json) or `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`; the
default stays the operational equal weighting.

## 3. Validation on real data

*This section is the validation of the budget as first delivered (eprof_v2.2 constants, per-gate
noise, no saturation; the vertical term emulated offline from the recorded sub-windows). The
independent review (§8) re-ran it with the exact code values, with eprof_v2.3, against
leakage-controlled truths and on fresh units, and changed three things; its numbers supersede the
ones below where they differ.*

### 3.1 Data and harness

* **Harness**: the q16a recalibration harness (`validation/chm15k_stability/q16a_recalibration.py`),
  which reproduces the published nightly constants bit-identically: `calibrate_rayleigh` unmodified,
  `options.json` + `eprof_v2.2`, the NetCDF write patched out. Each night's L1 is read once; the
  operational result and every budget ingredient come from that one in-memory night. The run records
  the constants of every 500 m sub-window over a WIDE range (from 1 km above the ground and 1.5 km
  under the window to a per-gate SNR of 1), so each design choice of §3.4 is re-evaluated offline
  without recalibrating; the shipped defaults were then re-run through the code (`verify`, §3.11).
  The CL61 run with the operational water-vapour line (910.55 nm, FWHM 0.188 nm, as the v2.2
  operational configuration; the code default 910.74/1.0 nm is the retired one).
* **Units** (36): 30 CHM15k — the twins Schiphol A–D, Budapest A/B, Măgurele A/B; Lindenberg,
  Granada, Évora, Potenza, Sofia, Aosta, Hohenpeißenberg, Palaiseau, Kleine Scheidegg, Wattisham,
  Schleswig, Oslo, Aldergrove, Zagreb, Tórshavn, Vásárosnamény, Gottfrieding, Madrid, Montsec,
  **Payerne A** (abnormal receiver, reported separately), **Schleiz**, **Messina** — plus 3 CL61
  (Lindenberg, Payerne, Aosta) and 3 Mini-MPL (Lille, Brest, Aleria).
* **Nights**: 2025-01-01 … 2026-07-31 in blocks of 7 consecutive nights (a night reads the L1 files of
  D−1 and D, so consecutive nights share files), one block every 2 weeks for the twins and key units,
  every 3 weeks otherwise; every unit uses the same calendar blocks, so twins share their nights.
  7 736 unit-nights, 2 436 successful calibrations, 2 090 of them with a usable time-series truth
  (Brest and Tórshavn have too few). Reproduction of the published constants: **99.0 %
  bit-identical** on 2 181 common nights (the remainder are night-key mismatches of the published
  file, as in q16a).

### 3.2 Two truths

1. **Time series** — the anomaly of ln C_L from the 31-day running median of the unit's published
   (all-night) series, after removing κ_u·(T_int − its 31-day running median), the sub-daily
   instrument-temperature term of report 19 (κ_u from the intra-night analysis, fallback 1.5 %/K,
   0 for Payerne, the CL61 and the Mini-MPL). Pull = anomaly / reported σ.
2. **Same-night twins** — the anomaly of ln(C_A/C_B) from its own 31-day running median, for the
   8 twin pairs. Atmosphere, molecular model and any gain change common to the site cancel; the
   pull is anomaly / √(σ_A² + σ_B²) with σ without the molecular term (it cancels too).

A realistic uncertainty gives a robust pull std of ~1 and ~5 % of the nights outside 2 σ.

### 3.3 The two truths disagree — by construction

Co-located units share most of their night-to-night scatter (published series, all nights):

| pair | type | nights | correlation of the 31-day anomalies | shared σ |
|---|---|---|---|---|
| Schiphol (6 pairs) | CHM15k / CHM15k | 122–153 | 0.55–0.80 | 8.7–11.9 % |
| Budapest A/B | CHM15k / CHM15k | 250 | 0.58 | 8.6 % |
| Măgurele A/B | CHM15k / CHM15k | 236 | 0.79 | 8.8 % |
| Lindenberg | CHM15k / CL61 | 142 | 0.36 | 5.4 % |
| Aosta | CHM15k / CL61 | 155 | 0.63 | 5.6 % |
| Payerne | CHM15k / CL61 | 39 | 0.62 | 25 % |

The co-located CL61 (another type, another wavelength, the same atmosphere) shares less with the
CHM15k than a CHM15k twin does. Roughly, a CHM15k's night-to-night variance splits into thirds:
**atmospheric** (shared with the CL61 — aerosol, molecular model: a calibration error the budget must
cover, but one that cancels in a twin ratio), **a CHM15k-type gain change common to the site** (shared
with the twin, not with the CL61: a real change of the instrument, not a calibration error), and
**unit-specific**. A budget covering exactly the calibration errors (atmospheric + unit-specific,
~0.7–0.8 of the variance) should therefore give time-series pulls of ~1.1–1.2 and twin pulls of
~0.6–0.7, not 1 and 1.

### 3.4 Where to probe: the design choice that matters

CHM15k without Payerne (29 units, 1 894 nights) and the 8 twin pairs (543 evaluable pair-nights);
vertical term |g|·D with D = 1 km plus the unexplained scatter; all other components identical.

| sub-window range of the vertical term | time-series pull | > 2 σ | twin pull | > 2 σ | median 1 σ |
|---|---|---|---|---|---|
| none (v2.2 ensemble only) | 4.05 | 63 % | 1.71 | 32 % | 2.7 % |
| the window only | 1.58 | 22 % | 0.79 | 6.8 % | 6.8 % |
| window + ≤ 1 km above (SNR ≥ 3) | 1.38 | 19 % | 0.67 | 4.4 % | 7.4 % |
| **+ 1 km below, floor 4 km AGL (shipped)** | **1.26** | **17 %** | **0.59** | **3.9 %** | 7.7 % |
| + 1 km below, floor 1.5 km (proposal in the brief) | 0.55 | 7 % | 0.25 | 0.7 % | 20.7 % |

Probing 1 km under every window (floor 1.5 km) triples the budget: the median gradient becomes
13.7 %/km instead of 3.4 %/km and the residual scatter 11.6 % instead of 3.2 %, because the layer under a low window is
the boundary-layer aerosol that the window search deliberately avoided. The twins are then
over-covered 3× (pull 0.25) — that range measures the aerosol, not the reliability of the window.
Limiting the downward extension to the free troposphere (above 4 km AGL) keeps the benefit for high
windows: lowering the floor to 3.5/3.0/2.5 km moves the twin pull to 0.56/0.53/0.48 and the time-series
pull to 1.18/1.05/0.88 (`validation/v23_uncertainty/u23_floor_test.py`). Using |g| rather than
√(g² − σ_g²) or √(g² + σ_g²) changes the pulls by < 0.1; the unexplained-scatter term is what brings
the twin tails from 11 % to 4 % outside 2 σ. On the twins the gradients of the two units correlate at
r = 0.76 and their window centres at 0.73: much of the vertical term is atmospheric and cancels in the
ratio, as §3.3 anticipates.

### 3.5 Results with the shipped defaults

| | nights | reported σ (median) v2.2 → v2.3 | observed σ | pull v2.2 → v2.3 | > 2 σ v2.2 → v2.3 |
|---|---|---|---|---|---|
| CHM15k, 29 units (Payerne excluded) | 1 894 | 2.7 → 7.7 % | 11.0 % | 4.05 → **1.26** | 63 → **17 %** |
| Payerne A | 46 | 6.1 → 15.3 % | 26.1 % | 2.95 → 1.47 | 67 → 15 % |
| CL61, 3 units | 111 | 2.7 → 7.5 % | 9.2 % | 3.10 → **1.02** | 55 → 11 % |
| Mini-MPL, 2 units | 39 | 1.8 → 3.1 % | 6.2 % | 3.19 → 1.86 | 49 → 33 % |
| Twins (Schiphol, Budapest, Măgurele) | 543 | — | 6.8 % (pair) | 1.71 → **0.59** | 32 → **3.9 %** |

Per twin site: Schiphol 1.66 → 0.62 (31 → 3.6 % outside 2 σ), Budapest 2.19 → 0.59 (36 → 4.4 %),
Măgurele 1.33 → 0.47 (29 → 4.3 %).

![Pull distributions (time series, twins) and reported vs observed σ per unit](C:/DATA/Projects/202610_v23/uncertainty/u23_pulls.png)

* **Coverage.** On the twins — the cleanest truth for the calibration error proper — the v2.3 budget
  over-covers slightly (pull 0.59, 3.9 % outside 2 σ, close to the Gaussian 4.6 %), as expected
  since its atmospheric part cancels in the ratio. On the time series it under-covers slightly (pull
  1.26) with heavy tails (17 % outside 2 σ): 45 % of the > 2 σ nights of the twin units have their
  twin also > 2 σ, same sign, the same night — common events (real site-wide gain changes or
  height-independent atmospheric effects) that no single-profile test can see.
* **Ranking.** Spearman correlation between the reported and the observed σ across units: 0.58 →
  0.67 (34 units; CHM15k only 0.57 → 0.61). Messina 1st (observed 1st), Payerne CL61 2nd (3rd),
  Vásárosnamény 3rd (7th), Payerne CHM15k 4th (2nd), Gottfrieding 5th (4th), Montsec 6th (5th).
  Messina carries the strongest gradient of the CHM15k network (median −10 %/km) and the Payerne CL61
  (−15 %/km) the known range-growing baseline of that unit.
* **Schleiz is not caught** (rank 15, observed 6th, pull 2.6): its within-window gradient is small
  (median −0.7 %/km); its scatter comes from level steps of about ×2 (2025-02, 2025-05, 2026-05, see
  §3.10) and multi-night episodes (2025-05-07/08 −112…−124 %, 2025-08-14…16 +73…+89 % with T_int
  4–6 K above its 31-day median). That is instrument instability, for the swap/step detection and the
  Kalman restart, not a per-night calibration error.
* **Mini-MPL** remains under-covered (pull 1.86): its photon noise and vertical term are tiny, and its
  scatter (6 %) comes from something the budget does not model; two units only.

Per unit (sorted by the v2.3 reported σ):

| unit | type | nights | observed σ | reported σ v2.2 | reported σ v2.3 | rank v2.3 / observed | median g [%/km] | pull v2.2 → v2.3 | > 2σ v2.2 → v2.3 |
|---|---|---|---|---|---|---|---|---|---|
| Messina A | CHM15k | 9 | 26.7 % | 14.6 % | 36.5 % | 1 / 1 | -10.0 | 1.7 → 0.51 | 33 → 11 % |
| Payerne C | CL61 | 14 | 25.7 % | 5.6 % | 17.9 % | 2 / 3 | -15.1 | 3.3 → 0.84 | 64 → 14 % |
| Vásárosnamény A | CHM15k | 56 | 16.0 % | 6.3 % | 15.6 % | 3 / 7 | -3.3 | 2.7 → 1.29 | 43 → 18 % |
| Payerne A | CHM15k | 46 | 26.1 % | 6.1 % | 15.3 % | 4 / 2 | -5.0 | 2.9 → 1.47 | 67 → 15 % |
| Gottfrieding 0 | CHM15k | 28 | 20.8 % | 6.3 % | 14.5 % | 5 / 4 | -5.0 | 2.9 → 1.11 | 57 → 18 % |
| Montsec A | CHM15k | 60 | 17.8 % | 5.5 % | 12.6 % | 6 / 5 | -7.0 | 2.2 → 0.98 | 55 → 15 % |
| Palaiseau B | CHM15k | 74 | 12.0 % | 3.8 % | 9.3 % | 7 / 11 | -3.6 | 3.2 → 1.27 | 55 → 15 % |
| Lindenberg C | CL61 | 58 | 8.5 % | 3.1 % | 8.4 % | 8 / 25 | -1.2 | 2.9 → 0.91 | 48 → 10 % |
| Oslo A | CHM15k | 24 | 9.9 % | 1.9 % | 8.2 % | 9 / 22 | -3.9 | 4.3 → 1.09 | 71 → 8 % |
| Schiphol C | CHM15k | 62 | 10.5 % | 2.9 % | 8.0 % | 10 / 15 | -1.7 | 4.2 → 1.20 | 63 → 6 % |
| Schiphol B | CHM15k | 70 | 11.2 % | 3.3 % | 8.0 % | 11 / 12 | -2.0 | 4.0 → 1.15 | 61 → 16 % |
| Budapest A | CHM15k | 99 | 15.3 % | 2.9 % | 7.9 % | 12 / 8 | -3.3 | 4.3 → 1.42 | 70 → 24 % |
| Măgurele RADO B | CHM15k | 99 | 10.1 % | 3.0 % | 7.8 % | 13 / 21 | +0.1 | 4.3 → 1.23 | 67 → 13 % |
| Schiphol D | CHM15k | 67 | 11.2 % | 3.2 % | 7.7 % | 14 / 13 | +0.2 | 3.5 → 1.30 | 57 → 13 % |
| Schleiz 0 | CHM15k | 95 | 17.6 % | 2.5 % | 7.7 % | 15 / 6 | -0.7 | 8.5 → 2.56 | 84 → 43 % |
| Măgurele A | CHM15k | 99 | 9.6 % | 3.0 % | 7.7 % | 16 / 23 | -0.1 | 3.7 → 1.09 | 63 → 10 % |
| Aldergrove A | CHM15k | 16 | 6.1 % | 3.0 % | 7.4 % | 17 / 33 | +1.1 | 1.8 → 0.55 | 25 → 0 % |
| Budapest B | CHM15k | 106 | 13.9 % | 2.1 % | 7.3 % | 18 / 9 | +1.6 | 6.9 → 1.50 | 77 → 25 % |
| Schiphol A | CHM15k | 72 | 10.2 % | 3.0 % | 7.3 % | 19 / 19 | -2.2 | 3.5 → 1.14 | 54 → 22 % |
| Potenza A | CHM15k | 53 | 8.5 % | 2.5 % | 7.2 % | 20 / 27 | -2.9 | 3.7 → 1.02 | 57 → 13 % |
| Aosta 0 | CHM15k | 84 | 10.4 % | 2.1 % | 7.2 % | 21 / 18 | -0.5 | 4.9 → 1.55 | 71 → 18 % |
| Sofia A | CHM15k | 65 | 8.4 % | 3.1 % | 7.1 % | 22 / 28 | -3.4 | 3.0 → 1.29 | 55 → 17 % |
| Kleine Scheidegg A | CHM15k | 77 | 7.6 % | 1.3 % | 7.0 % | 23 / 30 | -1.1 | 5.3 → 0.91 | 68 → 10 % |
| Zagreb A | CHM15k | 51 | 10.5 % | 2.6 % | 6.9 % | 24 / 16 | -1.6 | 3.6 → 0.97 | 57 → 8 % |
| Schleswig 0 | CHM15k | 47 | 7.1 % | 2.1 % | 6.8 % | 25 / 32 | -0.6 | 4.1 → 0.91 | 64 → 13 % |
| Lindenberg 0 | CHM15k | 94 | 10.2 % | 2.0 % | 6.8 % | 26 / 20 | -1.8 | 4.6 → 1.34 | 66 → 14 % |
| Évora B | CHM15k | 90 | 8.5 % | 3.0 % | 6.8 % | 27 / 26 | -1.9 | 2.8 → 1.08 | 52 → 13 % |
| Madrid A | CHM15k | 96 | 10.7 % | 2.4 % | 6.8 % | 28 / 14 | -3.4 | 4.4 → 1.23 | 67 → 21 % |
| Granada A | CHM15k | 73 | 12.0 % | 2.5 % | 6.5 % | 29 / 10 | -3.3 | 4.3 → 1.41 | 67 → 26 % |
| Wattisham A | CHM15k | 64 | 10.4 % | 1.9 % | 6.4 % | 30 / 17 | -0.7 | 4.1 → 1.47 | 69 → 19 % |
| Hohenpeißenberg 0 | CHM15k | 64 | 8.2 % | 1.6 % | 6.1 % | 31 / 29 | -2.4 | 5.4 → 1.36 | 67 → 16 % |
| Aosta B | CL61 | 39 | 7.3 % | 2.2 % | 6.0 % | 32 / 31 | +1.6 | 4.1 → 1.22 | 62 → 10 % |
| Lille A | Mini-MPL | 6 | 9.3 % | 4.6 % | 4.9 % | 33 / 24 | +0.0 | 2.0 → 1.92 | 17 → 17 % |
| Aleria A | Mini-MPL | 33 | 5.3 % | 1.7 % | 2.8 % | 34 / 34 | +0.1 | 3.0 → 1.81 | 55 → 36 % |

![Budget per unit: median nightly components against the observed variance](C:/DATA/Projects/202610_v23/uncertainty/u23_components.png)

### 3.6 How far does the window move? (the choice of D)

![Window-centre displacement statistics per unit](C:/DATA/Projects/202610_v23/uncertainty/u23_window_displacement.png)

From the published 2025–26 windows of the 33 CHM15k/CL61 units (network medians): the window centre
has a σ of 0.86 km from night to night (0.78 km about its 31-day median), consecutive nights differ
by 0.69 km, the summer window sits 1.2 km above the winter one, and the 10th–90th percentile span is
2.2 km. D = 1 km is the round value between the within-month σ and the seasonal shift; the
eprof_v2.2 tier-2 nights move further (1.5–2 km), which the vertical term covers only partly — those
nights are the ones with the largest gradients anyway. With the twin truth as the reference, D = 1 km
gives a twin pull of 0.59; D = 0.75 or 1.25 km would give ~0.65 or ~0.53 (window + 1 km above range:
0.72 / 0.67 / 0.60 for D = 0.75 / 1 / 1.25 km).

### 3.7 The relative QC limit

All units, 2 090 nights with a time-series truth; "off" = |anomaly| > 25 %:

| limit on u/C_L (2 σ) | nights rejected | of which off | kept nights off | Payerne A rejected |
|---|---|---|---|---|
| 40 % | 14.1 % | 28 % | 8.1 % | 37 % |
| 50 % | 10.1 % | 30 % | 8.7 % | 30 % |
| **60 % (default)** | **7.7 %** | **31 %** | **9.2 %** | **20 %** |
| 80 % | 5.0 % | 34 % | 9.7 % | 13 % |
| 100 % | 3.5 % | 34 % | 10.1 % | 11 % |
| v2.2 `u > C` (ensemble) | 0.08 % | — | — | — |

![Anomaly vs reported uncertainty, and the QC trade-off](C:/DATA/Projects/202610_v23/uncertainty/u23_qc.png)

The rejected nights are 3.4× more often off by more than 25 % than the kept ones at every limit; the
separation does not sharpen with a stricter limit, so the limit is set by the Kalman argument of §2.7
(60 % = 4 × the median 2 σ of ~15 %). An operator who values availability more can take 80–100 %.

### 3.8 Payerne 2026-09-07

The published constant of that night is 1.02 × 10¹² against a 2026 unit median of 5.26 × 10¹¹ (×1.94),
from a window 5.57–6.05 km AGL (6.06–6.54 km ASL) at the top of the usable range of a noisy receiver;
it calibrates only with the classification screen (flag −2 without it), which the validation
reproduces with the runner's own classification. The v2.2 ensemble gave 42 % (2 σ).

| night | C_L / unit median | v2.2 2 σ | v2.3 2 σ | of which vertical | gradient | v2.3 QC |
|---|---|---|---|---|---|---|
| 2026-09-06 | 0.55 | 28 % | 63 % | 27 % | −14 %/km | rejected |
| **2026-09-07** | **1.94** | **42 %** | **162 %** | **155 %** | **−77 %/km** | **rejected** |
| 2026-09-12 | 0.99 | 8 % | 25 % | 21 % | −2 %/km | kept |
| 2026-09-23 | 1.00 | 6 % | 26 % | 23 % | −8 %/km | kept |

![Sub-window constants against height: a typical Schiphol night, a typical Payerne night, the steepest Schleiz night and Payerne 2026-09-07 (range on the Y axis)](C:/DATA/Projects/202610_v23/uncertainty/u23_profiles.png)

On 2026-09-07 the constant is 1.23 × the window value 500 m under the window and 2.0 × 1 km under
it: the free-troposphere extension of the vertical test sees it (−77 %/km), the window-only range
would not (two half-windows, +23 %/km, 68 %). Both odd nights are rejected, both ordinary nights
kept; the noisy 09-06 night is caught mainly by its photon noise (49 % 2 σ).

### 3.9 Synthetic test

A pure molecular CHM15k profile (US-1976, 30 m gates) with a background-limited noise of 5 % at 5 km,
window 3–4 km: the clean profile gives g = −0.1 %/km and σ_vert = 0.13 %. Adding a constant
electronic offset B·r² on the range-corrected signal equal to ±5 % of the molecular power at 5 km
biases the window constant by ±2.0–2.5 % and gives g = +1.6 / −1.9 %/km, σ_vert = 1.6 / 1.9 %
(12–14 × the clean value; 2σ covers the bias); ±10 % gives g ≈ ±3.2 %/km. Tests:
`tests/test_uncertainty_v23.py::test_synthetic_distortion_inflates_vertical_term`.

### 3.10 Kalman weighting

`kalman_best_estimate` on each unit's recalibrated series (30 units with ≥ 30 nights), equal weights
(operational) against weights from the v2.3 uncertainty:

| metric | weighted vs equal (median over units) | better on |
|---|---|---|
| twin best estimates: rms deviation of ln(K_A/K_B) from its 31-day median | −4.8 % | 7 of 8 pairs |
| same, 95th percentile | −8.2 % | 6 of 8 pairs |
| day-to-day jitter of the best estimate | −3.8 % | 22 of 30 units |
| one-step-ahead innovation, the more certain half of the nights | +0.2 % | 15 of 30 units |
| one-step-ahead innovation, all nights | +1.9 % | 12 of 30 units |

![Kalman best estimate with equal weights and weighted by the v2.3 uncertainty (Payerne, Schleiz, Schiphol B)](C:/DATA/Projects/202610_v23/uncertainty/u23_kalman.png)

Weighting makes co-located best estimates more consistent and the estimate smoother, without
degrading the prediction of reliable nights (it follows the uncertain ones less, by design). The
effect is moderate because the filter keeps its empirical variance scale and caps the weight ratio at
4; it grows with the share of uncertain nights (eprof_v2.2 tier-2 recoveries, Payerne-type units).

### 3.11 Code check and cost

* The shipped defaults were re-run through the code on 41 nights of 5 units (Schiphol B, Schleiz,
  Payerne A, Granada, Lille): same constant on 100 % of the nights; σ_vert identical to the offline
  evaluation within 0.1 %-point on 85 % of them (correlation 0.998), the rest differing by 2–5 points
  where the code clips a partial sub-window at the 4 km floor, which the offline evaluation
  approximates.
* Runtime: **2.3 ms per night** (95th percentile 3.2 ms), about ten extra Klett runs, 0.2 % of the
  1.2–1.5 s a night takes (L1 read included). No extra input is read.

## 4. Limits

* **The vertical term does not identify the cause** — by design. It is large for an electronic
  distortion, for a faint aerosol layer near the window and for a window that is not purely
  molecular alike, and it is also large on aerosol-structured nights whose chosen window happens to
  be clean (a conservative error). Within-night gradients are mostly aerosol on ordinary units (they
  steepen in summer; altitude-independence audit), a distortion on Payerne-type units.
* **Height-independent errors are invisible** to it: a vertically uniform free-tropospheric aerosol,
  a gain or calibration-pulse error, the common-mode CHM15k gain cycle of report 19.
* **The two truths bracket the answer — loosely.** Part of a CHM15k's night-to-night variance is a
  gain change common to co-located twins but not to the co-located CL61 (§3.3), i.e. a real change of
  the instrument, not a calibration error, so the time-series pull over-states the error; the twin
  ratio cancels the atmospheric part of the error, so the twin pull under-states it. The direction is
  sound, but the review (§8.2) propagated the measured variance shares with their confidence
  intervals: a perfectly calibrated budget would give a time-series pull anywhere in 1.00–1.46 and a
  twin pull in 0.40–0.77. The bracket pins the level of the budget to about ±30 %, not to the
  "~1.2 and ~0.65" a strict reading of the thirds suggests.
* **Tails.** Even with the budget, 17 % of the CHM15k nights fall outside 2 σ of the time series
  (Gaussian 4.6 %); nearly half of them are site-wide events seen by both twins. A per-night budget
  cannot cover those; the Kalman process noise and the step detection must.
* **Instrument instability is not calibration uncertainty.** Units with level steps or episodic gain
  jumps (Schleiz) keep pulls well above 1; that belongs to the module-swap/step detection and the
  Kalman restart, not to the nightly uncertainty.
* **D is one network constant** (1 km), but the displacement picture behind it is not supported: a
  night's anomaly does not follow g × (its window displacement from the unit's usual height) (§2.3).
  The vertical term works as an empirical indicator of a worse night; its scale (D, the saturation at
  15 %) is set empirically, and the truths pin it to about ±30 % (§8.2). The 4 km floor of the
  downward extension is a physical choice (free troposphere), checked but not tuned (3–4 km give
  similar results).
* **Conditional calibration.** Even reviewed, the quietest decile of the CHM15k nights is
  under-covered on the time-series truth (pull ~1.5–2). On the twin truth the same nights are
  calibrated (pull 1.0): the floor is a real gain change of the instrument (Kalman process noise), not
  a calibration error.
* **The molecular climatology** was fitted on CHM15k nights at 37–55 °N; the latitude is clipped to
  35–62 °N (southern stations use the mirrored season). North of 62 °N it under-covers (79 % of the
  nights within 2 σ, mean bias +1 %, §2.4). Per-night CAMS is a drop-in replacement if its read cost
  is accepted.
* **The noise term** is measured on 300-s profiles; it includes the atmospheric change over 5 min (an
  upper bound) and not slower noise. Nights shorter than 8 binned profiles fall back to the per-gate
  estimate, which is 1.6–4 × too small.
* **Mini-MPL**: the constant of the Mini-MPL drifts within a night (first vs second half) by ~45 × its
  photon noise, an instrument gain drift that no calibration uncertainty covers; with an every-night
  truth the two fresh Mini-MPL units are under-covered (pull ~2.4). This belongs to the instrument
  monitoring and the Kalman process noise. A causal unit-level repeatability term (the unit's own
  recent scatter beyond its budget) would bring every time-series pull to ~0.9 but over-covers the
  twins (0.59 → 0.41): it absorbs real gain changes and is therefore not added to the nightly budget.
* **The validation truths leak and are noisy.** The u23 time-series truth excludes only the night
  itself; neighbouring nights share weather and aerosol regimes and pull the running median towards
  the night: excluding ±3 days raises the CHM15k pull by 15 % and the twin pull by 18 %. The running
  median also has its own error (median 3.6 % for the CHM15k), which inflates the raw pull by ~10 %.
  The two effects nearly cancel; §8 quotes both truths.

## 5. Proposal for L2

Today the L2 hub writes `uncertainties_att_backscatter_0` as a fixed 25 % of |β_att|. That is too large
for a well-calibrated instrument at high SNR, too small where the signal is noise-dominated, and blind
to the calibration (a typical Payerne CHM15k night, 31 % at 2 σ in v2.3, gets the same 25 % as a
typical Schiphol night at 16 %). Proposal, in two steps:

1. **Carry the calibration uncertainty.** Add `calibration_constant_0_uncertainty` (same units as
   `calibration_constant_0`, **1 σ standard uncertainty**, attribute `coverage_factor = 1`; the
   calibration NetCDF keeps its historical 2 σ) to L2. Its value is the uncertainty of the constant
   that L2 actually applies: the Kalman best estimate's standard deviation (`kalman_std` of the
   runner, which already combines the nightly uncertainties — when weighting is on — with the process
   noise accumulated since the last calibration), or u/2 of the nightly v2.3 constant when the
   nightly value is applied directly.
2. **Propagate it instead of the fixed 25 %.**
   u(β_att)/β_att = √[ (u_C/C)² + (σ_P/P)² ],
   with σ_P/P the per-gate relative noise of the profile (temporal first differences on the native
   profiles — the estimator v2.2 already runs, computed on the L2 time grid — or the instrument's own
   noise estimate). For the CHM15k a third, daytime term κ_u·ΔT_int (report 19: +5 % median, +7.5 %
   in summer, >10 % on a third of the units) should be added until the temperature correction is
   applied.

## 6. How to enable

| what | where | value |
|---|---|---|
| the v2.3 budget (C_L unchanged, u = full budget, components in the NetCDF, relative QC) | `options.json` | `"uncertainty_model": "v2.3"` |
| tune it (keys of `calibration.rayleigh.uncertainty.DEFAULT_PARAMS`) | `options.json` | `"uncertainty_params": {"window_displacement_m": 1000, "vertical_cap_rel": 0.15, "max_rel_uncertainty": 0.6, ...}` |
| the pre-review v2.3 budget (per-gate noise, climatology at the window centre AGL, no saturation) | `options.json` | `"uncertainty_params": {"noise_estimator": "per_gate", "molecular_height_reference": "agl", "vertical_cap_rel": null}` |
| same, for a network run without editing options.json | environment | `ALC_UNCERTAINTY_MODEL=v2.3` |
| weight the Kalman best estimate by the nightly uncertainty | `options.json` / environment | `"kalman_uncertainty_weighting": 1` / `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1` |

Defaults keep v2.2 exactly (`"ensemble"`, weighting off); the output of the default configuration is
bit-identical to v2.2 (tested: constants, uncertainties, flags, messages and every NetCDF variable
and attribute on the CHM15k, Mini-MPL and CL61 fixtures; `tests/test_uncertainty_v23.py`). The
integration defines molecular method `eprof_v2.3` with this option on.

## 7. Reproduction

```
# recalibration of the unit list with the v2.3 recording configuration (q16a harness, L1 read once
# per night, 8 workers): C:/DATA/Projects/202610_v23/uncertainty/runs/
python validation/v23_uncertainty/u23_validation.py run --workers 8
python validation/v23_uncertainty/u23_analysis.py windows     # window displacement (D)
python validation/v23_uncertainty/u23_common_mode.py          # twins vs co-located CL61
python validation/v23_uncertainty/u23_noise_check.py          # sigma_signal propagation check
python validation/v23_uncertainty/u23_analysis.py sweep       # design choices, pulls, twins
python validation/v23_uncertainty/u23_floor_test.py           # floor of the downward extension
python validation/v23_uncertainty/u23_analysis.py report      # final tables + figures
python validation/v23_uncertainty/u23_analysis.py payerne     # Payerne 2026-09-07 (+ classification)
python validation/v23_uncertainty/u23_analysis.py profiles    # sub-window profile figure
python validation/v23_uncertainty/u23_analysis.py kalman      # Kalman weighting
python validation/v23_uncertainty/u23_analysis.py verify      # shipped defaults in the code vs offline
python validation/v23_uncertainty/u23_analysis.py cost        # runtime

# independent review (section 8): a second harness reads each night's L1 once and computes both window
# methods, the exact shipped and revised budgets, split-night constants and the ensemble dimensions
# (12 workers); outputs in C:/DATA/Projects/202610_v23/uncertainty_review/
python validation/v23_uncertainty/u23r_harness.py run --set orig  --workers 12   # the 36 design units
python validation/v23_uncertainty/u23r_harness.py run --set fresh --workers 12   # 37 fresh units
python validation/v23_uncertainty/u23r_harness.py run --set origx --workers 12 --units <twins, Payerne A>
python validation/v23_uncertainty/u23r_truth.py                       # reproduction, truth leakage
python validation/v23_uncertainty/u23r_reliability.py orig eprof_v2.3  # conditional calibration
python validation/v23_uncertainty/u23r_reliability2.py orig eprof_v2.3 # twins, |g| octiles, displacement
python validation/v23_uncertainty/u23r_components.py fresh eprof_v2.2  # overlaps, split-night noise test
python validation/v23_uncertainty/u23r_bracket.py                      # variance shares, perfect-budget box
python validation/v23_uncertainty/u23r_saturation.py orig eprof_v2.3 shipped   # choice of the cap
python validation/v23_uncertainty/u23r_oos.py fresh eprof_v2.3         # before / after out of sample
python validation/v23_uncertainty/u23r_qc.py orig                      # QC cost, ROC, eprof_v2.3 interplay
python validation/v23_uncertainty/u23r_kalman.py origx eprof_v2.3 200  # Kalman: permutation test, variants
python validation/v23_uncertainty/u23r_ranking.py orig eprof_v2.2      # unit ranking with bootstrap CIs
python validation/v23_uncertainty/u23r_unitterm.py orig eprof_v2.2     # causal unit-level term
python validation/v23_uncertainty/u23r_molclim.py                      # climatology refit (ASL convention)
python validation/v23_uncertainty/u23r_molecular_oos.py                # climatology out of sample
python validation/v23_uncertainty/u23r_payerne.py                      # Payerne 2026-09-07, both methods
python validation/v23_uncertainty/u23r_figures.py                      # review figures
```

## 8. Independent review (2026-10)

An adversarial review attacked every claim of this page with a second harness
(`validation/v23_uncertainty/u23r_harness.py`: each night's L1 read once; both window methods,
eprof_v2.2 and eprof_v2.3; the EXACT shipped budget instead of an offline emulation; split-night
constants; the three ensemble dimensions separately), with leakage-controlled truths, and on
**37 fresh units never used in the design** (27 CHM15k in 14 countries, 12–69 °N, two mountain
sites; 8 CL61; 2 Mini-MPL with an every-night series). Full report, tables and figures:
`C:/DATA/Projects/202610_v23/uncertainty_review/u23r_review.md`.

### 8.1 What survived

* **The pooled pulls**, with the exact code values: CHM15k time series 4.05 → 1.25 (17 % outside
  2 σ), twins 1.69 → 0.60 (3.5 %), CL61 3.13 → 1.08 (design units, eprof_v2.2; the doc's 1.26 / 0.59 /
  1.02 came from the offline emulation). Out of sample (fresh units, eprof_v2.3): CHM15k 3.91 → 1.17
  (16 % outside 2 σ), CL61 2.85 → 1.08. Time hold-out: 1.27 in 2025, 1.33 in 2026 (eprof_v2.3). The
  ranking of the design candidates of §3.4 is the same on the fresh units and in both years, and the
  weight floor, the aerosol allowance and the molecular term move the pull by ≤ 0.15: no sign of
  in-sample over-tuning of the pooled number.
* **The √k correction** of the per-gate noise on the CHM15k (propagated / measured 0.965).
* **Payerne 2026-09-07** with eprof_v2.2: 42 % → 162 % (2 σ), rejected — reproduced exactly. With
  eprof_v2.3 the window search itself rejects that night (flag −2) — and also the ordinary 09-23 night
  at the unit median.
* **The molecular term** as an absolute uncertainty (|b| ≤ 2 σ_mol on 97 % of the fitted nights,
  96–99.7 % out of sample up to 62 °N).
* **No significant double counting.** The ensemble's shift and time dimensions overlap the vertical and
  noise terms by 3 % of the variance at the median night (total σ 1.7 % too large); the noise inside
  |g| is 1 % of E[g²] (odd/even halves); the climatological part of σ_mol moves the pull by 2 %.
* **The Kalman weighting** — and it is worth more than claimed with eprof_v2.3 (§8.3).

### 8.2 What did not survive (ranked by impact)

1. **Conditional calibration.** In deciles of the reported σ the observed CHM15k scatter rose only from
   7.6 to 12.8 % while the reported σ rose from 4.0 to 29.8 % (eprof_v2.3): pull 1.98 in the quietest
   decile, 0.44 in the noisiest. The pooled 1.25 averaged an under-covered bulk with an over-covered top.
   Per-unit pulls were over-dispersed (permutation p = 0.0015). The twin truth splits the causes: the
   quiet nights are calibrated there (pull 1.0 in the bottom quintile: the time-series floor is a real
   gain change), the top quintile is over-covered 3–4 ×. The top is the vertical term growing like
   |g| × D (§2.3), which the time-series truth — an upper bound of the calibration error — refutes.
   The same pattern holds out of sample (fresh units: 1.65 … 0.46).
2. **The photon-noise term** was 1.6 × too small on the CHM15k, 3.4–4.2 × on the CL61, 2.8–3.1 × on the
   Mini-MPL (split-night test): gates are not independent, and on the CL61 the √k white-noise
   correction itself over-corrects by ~1.9 × (§2.5).
3. **The relative QC** rejected good nights (40–50 % of its rejections less than 10 % off), was mostly
   redundant with eprof_v2.3's vertical-consistency test, and made the weighted Kalman best estimates
   less consistent (§2.7).
4. **The validation truths.** Excluding only the night itself lets neighbouring nights of the same
   weather regime into the running median: excluding ±3 days raises the CHM15k pull 1.25 → 1.43 and the
   twin pull 0.60 → 0.71. The median's own error (3.6 %) inflates the raw pull by ~10 %. The bracket
   argument holds in direction only: a perfect budget would give 1.00–1.46 / 0.40–0.77 (§4).
5. **The molecular climatology** was evaluated at the window centre AGL while fitted ASL (+0.49 % per km
   of station altitude); its fit script was not in the repository; it under-covers north of 62 °N.
6. **The unit ranking gain** 0.58 → 0.67 is not significant: 0.56 → 0.72 with the exact code, difference
   +0.16 (95 % CI −0.12…+0.31); fresh units 0.40 → 0.58 (−0.14…+0.40); the unit's median |g| alone
   ranks the design units as well (0.75). The budget does beat a unit's own past scatter as a predictor
   of next year's (0.70 vs 0.50 design, 0.65 vs 0.35 fresh).
7. **The Mini-MPL** under-coverage is confirmed with an every-night truth (pull 2.3 on two fresh
   units): its constant drifts within a night by ~45 × its photon noise — instrument drift, not a
   calibration error the budget could model.

### 8.3 What changed

The reported constant is unchanged, `uncertainty_model = "ensemble"` stays bit-identical to v2.2, and
the pre-review v2.3 budget remains reproducible through `uncertainty_params` (§6). Three defaults of the
v2.3 model changed, none of them a new fitted parameter except the cap:

| change | parameter (default) | why |
|---|---|---|
| photon noise measured as the noise of the window mean on the binned stack | `noise_estimator = "window_mean"` | the independent-gate estimate is 1.6–4 × too small (split-night test) |
| vertical term saturates | `vertical_cap_rel = 0.15` | the error of the steepest nights is ≥ 2.3 × below |g| × D (time-series truth = upper bound) |
| molecular climatology at the window centre ASL | `molecular_height_reference = "asl"` | the convention the coefficients were fitted with |

Before / after on the SAME held-out nights (fresh units, eprof_v2.3, exact code values; "±3 d" = the
leakage-controlled truth):

| fresh units, eprof_v2.3 | nights | v2.2 ensemble | v2.3 as first shipped | v2.3 reviewed |
|---|---|---|---|---|
| CHM15k: pooled pull (> 2 σ) | 1 753 | 3.91 (59 %) | 1.17 (16 %) | 1.30 (18 %) |
| CHM15k: pull in the bottom / top decile of the reported σ | | 6.52 / 1.63 | 1.65 / **0.46** | 1.59 / **1.02** |
| CHM15k: spread of the decile pulls (sd of log) | | 0.44 | 0.38 | **0.15** |
| CHM15k: median reported 1 σ | | 2.4 % | 7.5 % | 7.4 % |
| CHM15k, ±3 d truth: pooled pull (> 2 σ) | 1 586 | 4.15 (62 %) | 1.29 (19 %) | 1.44 (21 %) |
| CHM15k, ±3 d truth: top decile / spread | | 1.69 / 0.44 | 0.46 / 0.42 | 1.00 / 0.20 |
| CL61: pooled pull (> 2 σ) | 299 | 2.85 (50 %) | 1.08 (12 %) | 1.20 (16 %) |
| Mini-MPL (every-night own series) | 299 | 3.06 (52 %) | 2.30 (39 %) | 2.32 (39 %) |
| nights over the 0.6 QC limit (CHM15k / CL61) | | — | 5.6 % / 7.4 % | 0 % / 0 % |
| unit ranking (Spearman, 32 units) | | 0.40 | 0.58 | 0.59 |
| split-night noise test, CHM15k / CL61 / Mini-MPL (1 = right) | | — | 1.58 / 3.38 / 2.84 | 0.96 / 0.88 / 0.99 |
| *twins (design units, 520 pair-nights)* | | 1.66 (31 %) | 0.63 (3.8 %) | 0.69 (4.2 %) |

![Conditional calibration: time-series truth per decile (left), twin truth per quintile (middle), the vertical term against |g| (right)](C:/DATA/Projects/202610_v23/uncertainty_review/u23r_fig_reliability.png)

![Split-night test of the photon-noise term per instrument type](C:/DATA/Projects/202610_v23/uncertainty_review/u23r_fig_noise.png)

The pooled time-series pull rises from 1.17 to 1.30: the shipped number was lowered by the
over-covered top decile. Both, and the twin pull of 0.69, lie inside the range a perfectly calibrated
budget would give (time series 1.00–1.46, twins 0.40–0.77, §4). What the revision buys is conditional
calibration — the steepest decile goes from 2.2 × over-covered to calibrated and the spread of the
decile pulls falls by 60 %, out of sample and with both truths — and a noise term that passes the
split-night test for every instrument type. The bottom decile (pull ~1.6 on the time-series truth) is
a real gain change: on the twin truth it is calibrated. With the Kalman weighting on (eprof_v2.3, twin
units), the reviewed budget makes the co-located best estimates 13.6 % (rms) and 13.1 % (95th
percentile) more consistent than equal weights; rejecting on u/C would cut that to ~6 %. Payerne
2026-09-07 with eprof_v2.2 windows: 162 % → 71 % (still rejected, now carried by the noise of that
receiver, 24 → 48 %).

### 8.4 Recommendation

Release the reviewed budget (`uncertainty_model = "v2.3"`, the new defaults) together with
`eprof_v2.3` and the Kalman weighting (`kalman_uncertainty_weighting = 1`). Treat the weighting, not
the u/C rejection, as the way uncertain nights are handled: keep `max_rel_uncertainty = 0.6` only as a
safety net for extreme noise (with the saturated vertical term it fires on < 0.5 % of the nights) and
do not lower it. For L2 (§5) carry the Kalman standard deviation, which adds the real gain changes the
nightly budget deliberately leaves out (they are why the time-series pull stays above 1 and why the
Mini-MPL and Schleiz are under-covered). Open items: the vertical term's scale rests on one network D
and one cap, pinned only to about ±30 % by the available truths; north of 62 °N the molecular term
should use per-night CAMS T/p; the southern-hemisphere mirroring is untested.

## References

* Wiegner, M. and Geiß, A.: Aerosol profiling with the Jenoptik ceilometer CHM15kx, Atmos. Meas. Tech.,
  5, 1953–1964, 2012.
* Klett, J. D.: Lidar inversion with variable backscatter/extinction ratios, Appl. Opt., 24,
  1638–1643, 1985; sign correction: Speidel, J. and Vogelmann, H., Appl. Opt., 62, 861–868, 2023.
* JCGM 100:2008, Evaluation of measurement data — Guide to the expression of uncertainty in
  measurement (GUM).
* Mattis, I., D'Amico, G., Baars, H., Amodeo, A., Madonna, F., and Iarlori, M.: EARLINET Single
  Calculus Chain — technical — Part 2: Calculation of optical products, Atmos. Meas. Tech., 9,
  3009–3029, 2016.
