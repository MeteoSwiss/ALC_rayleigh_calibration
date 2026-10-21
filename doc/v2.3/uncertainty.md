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
  |d ln C_L/dz| × 1 km plus the scatter the line does not explain — a molecular-model climatology,
  a residual-aerosol allowance and the photon noise, in quadrature, reported at 2 σ.
* **Validation** on 36 units (30 CHM15k, 3 CL61, 3 Mini-MPL; 7 736 recalibrated nights 2025-01 …
  2026-07, 99.0 % bit-identical to the published constants): CHM15k time-series pulls 4.05 → **1.26**
  (63 % → 17 % outside 2 σ), same-night twins 1.71 → **0.59** (32 % → **3.9 %**), CL61 3.10 → 1.02.
  The two truths bracket the answer by construction: about a third of a CHM15k's scatter is a real
  gain change common to co-located twins (not a calibration error) and another third is atmospheric
  (an error that cancels in a twin ratio), so ~1.2 and ~0.65 are the expected values.
* **Large where it must be:** Messina, Payerne CL61, Vásárosnamény, Payerne CHM15k, Gottfrieding and
  Montsec rank 1–6 of 34 (observed ranks 1–7). The Payerne 2026-09-07 night (window 6.06–6.54 km ASL,
  C_L twice the unit median) goes from 42 % to **162 %** (2 σ) and is rejected by the new QC.
  Schleiz does not come out large (rank 15): its scatter is ×2 level steps, an instrument problem
  no per-night uncertainty can see.
* **QC:** `u > C` fired on 0.08 % of the nights; the relative limit u/C_L > 60 % (2 σ) rejects 7.7 %,
  3.4 × enriched in nights more than 25 % off.
* **Kalman** weighting by the new uncertainty: twin best estimates 5–8 % more consistent, the best
  estimate 4 % smoother, one-step-ahead prediction unchanged. **Cost:** 2.3 ms per night (0.2 %).
* **Side finding:** v2.2's propagated noise (`sigma_signal`) is 4.1–4.9 × too large (time binning
  missing); v2.2 is left untouched, the budget corrects it, the tier-2 χ² gates should be revisited.
* **Recommendation:** adopt the budget with `eprof_v2.3`, switch the Kalman weighting on with it, and
  carry the calibration uncertainty into L2 instead of the fixed 25 % (§5).

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
to the ±200 m shifts (g ≈ 2σ_ens/0.3 per km), a conservative fallback.

What the term measures — and what it does not:

* an **electronic distortion** (an additive baseline B·r² on the range-corrected signal) makes C_k
  grow (B > 0) or fall (B < 0) with height roughly as r²/β_mol; a synthetic molecular profile with
  B equal to 5 % of the molecular signal at 5 km gives g = +1.6 %/km and a vertical term 12 × that of
  the clean profile (§3.9);
* **aerosol inside a sub-window** raises its C_k (the reference assumes it is molecular); aerosol
  between two sub-windows is corrected by the Klett step with the assumed lidar ratio;
* **noise** alone makes g scatter by its standard error; `|g|` therefore carries the noise of the
  gradient estimate on average (E[g²] = g_true² + σ_g²);
* the term does **not** see a height-independent bias (a uniform free-tropospheric aerosol, a
  calibration-pulse or gain error) — those are invisible to any single-profile method.

### 2.4 (c) Molecular-model term σ_mol

The ln-bias of C_L from using US-1976 instead of the real T/p was computed per night with CAMS
(0.4°, L137) in each night's own window for 22 618 nights of 127 CHM15k units (report 19, q15d). It
is well described by a climatology (fit residual 1.35 %):

  b [%] = −0.421 + 1.516 cos ωd + 0.799 sin ωd + 0.080 (φ−47) − 0.017 (φ−47) cos ωd
          − 0.493 (z_c−4.5) − 0.210 (z_c−4.5) cos ωd + 0.008 h_stn

(d = day of year, ω = 2π/365.25, φ = latitude in °, z_c = window centre in km AGL, h_stn = station
altitude in km; southern stations use the mirrored season). σ_mol = √(b² + 1.35 %²): the expected
bias of the night plus its synoptic scatter. Reading CAMS T/p per night would cost ~1 s (daily file)
to ~3 s (8 GB monthly file) — as much as the whole calibration — while the climatology already
reproduces the realised RMS (1.98 %). With `molecular_source = "cams"` the reference is already CAMS
and only a 0.5 % floor remains.

### 2.5 (d) Photon-noise term σ_noise

From the per-gate noise profile that eprof_v2.2 measures on the native profiles (temporal first
differences / √2), referred to the expected molecular signal ĉ·p_mol (ĉ = median signal/p_mol in the
window: the measured signal crosses zero at a per-gate SNR below ~1). The constant is a window MEAN
of RCS/β_mol (the gate-to-gate noise cancels in the Klett normalisation), so
σ_noise = √(Σ e_i²)/N over the N window gates, e_i = σ_i/(ĉ p_mol,i).

**A propagation error in v2.2.** `calibrate_rayleigh` divides the noise of ONE NATIVE profile (12–15 s
on the CHM15k) by √(number of 300-s BINNED profiles), so its `sigma_signal` is too large by √k,
k = native profiles per bin = 20–25 on the CHM15k. Measured against the scatter of consecutive binned
profiles: ratio 4.1–4.9 on five CHM15k nights (√k = 4.5–5.0; `validation/v23_uncertainty/
u23_noise_check.py`). The v2.3 budget divides by √k. The v2.2 gates are left untouched (bit-identical
output), but note that their χ²_red ≤ 2.5 tolerance is therefore a true χ²_red of ~50, and that the
per-gate SNR quoted in report 19 (q16a/q16c, e.g. "SNR crosses 1 at ~6 km") is ~4.5 × too pessimistic.

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

### 2.8 Kalman weighting

`monitoring/kalman.kalman_best_estimate` already accepts `uncertainties=`: it keeps the empirical
measurement-variance scale (day-to-day scatter) and takes only the night-to-night RATIO of the
relative uncertainties (clipped to a factor 4 in σ). The network runner now passes them when
`kalman_uncertainty_weighting = 1` (options.json) or `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`; the
default stays the operational equal weighting.

## 3. Validation on real data

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
* **The two truths bracket the answer.** About a third of a CHM15k's night-to-night variance is a
  gain change common to co-located twins but not to the co-located CL61 (§3.3), i.e. a real change
  of the instrument, not a calibration error; the time-series pull therefore over-states the error
  and is expected slightly above 1. Conversely the twin ratio cancels the atmospheric part of the
  error (another third), so the twin pull is expected below 1.
* **Tails.** Even with the budget, 17 % of the CHM15k nights fall outside 2 σ of the time series
  (Gaussian 4.6 %); nearly half of them are site-wide events seen by both twins. A per-night budget
  cannot cover those; the Kalman process noise and the step detection must.
* **Instrument instability is not calibration uncertainty.** Units with level steps or episodic gain
  jumps (Schleiz) keep pulls well above 1; that belongs to the module-swap/step detection and the
  Kalman restart, not to the nightly uncertainty.
* **D is one network constant** (1 km). A per-unit D (that unit's window-height scatter) or a term
  growing with the distance of the night's window from the unit's usual height would sharpen the
  budget; neither is implemented. The 4 km floor of the downward extension is a physical choice
  (free troposphere), checked but not tuned (3–4 km give similar results).
* **The molecular climatology** was fitted on CHM15k nights at 36–60 °N; elsewhere it is clipped to
  that latitude range (southern stations use the mirrored season). Per-night CAMS is a drop-in
  replacement if its read cost is accepted.
* **The noise term assumes white noise** after the √k correction; the scatter of consecutive 300-s
  profiles is 5–10 % larger than that (a small correlated component), absorbed by the 1 % weight floor.
* **CL61 and Mini-MPL** are validated on three units each (two Mini-MPL with a usable truth, on their
  own sampled series); the Mini-MPL scatter (6 %) is not explained by the budget (3 %).
* **The validation sample** is 2025-01 … 2026-07 with one 7-night block every 2–3 weeks; the time-series
  truth uses the published all-night series, so its 31-day median is not degraded by the sampling.

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
| tune it (keys of `calibration.rayleigh.uncertainty.DEFAULT_PARAMS`) | `options.json` | `"uncertainty_params": {"window_displacement_m": 1000, "max_rel_uncertainty": 0.6, ...}` |
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
```

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
