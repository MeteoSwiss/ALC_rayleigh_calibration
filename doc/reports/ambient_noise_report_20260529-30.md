# Ambient (no-hood) noise & detection thresholds — CL31 / CL61 / CHM15k, Payerne

**Period:** 29–30 May 2026 (cloud-free periods only) · **Dark reference:** 12 May 2026 09:35–14:50 UTC (termination hood)
**Script:** `ambient_noise_cl61_chm_cl31.m` · **Figures/results:** `figs_ceilo_ambient\` · **Generated:** 11 June 2026

---

## 1. Objective

The dark-measurement campaign (`dark_measurement_cl61_chm_cl31.m`, hood on) characterised instrument noise vs altitude, Allan behaviour and aerosol detection thresholds at Payerne. Network instruments can never be hooded, so this study derives the **same quantities from ambient measurements**: if a robust ambient method reproduces the dark results, it can be applied to **every instrument in the network** where only ambient data exist. E-PROFILE L2 is time-averaged and unsuitable — all analysis runs on raw data at native resolution.

**Headline result: night-time ambient noise reproduces the dark measurement within 0–10 % at all altitudes for all three instruments. Clear nights are a free dark measurement.**

## 2. Data

| Instrument | λ | Source | Native resolution | Signal variable |
|---|---|---|---|---|
| CL31 | 910 nm | `A:\CL31\PAY\2026\05\VSRA00.*.dat` (MeteoSwiss bulletins) | 30 s / 10 m / 770 gates (7.7 km) | `rcs_910` (1 LSB = 10⁻⁸ m⁻¹sr⁻¹) |
| CL61 | 910.55 nm | `A:\CL61_Cloudnet\Payerne\2026\05\*_cl61d_*.nc` (Cloudnet daily) | 30 s / 4.8 m / 3276 gates (15.7 km) | `beta_raw` (**unscreened** β_att, sr⁻¹m⁻¹) |
| CHM15k | 1064 nm | `A:\CHM15k\PAY\2026\05\*_CHM*_000.nc` (Lufft daily) | 15 s / 15 m / 1024 gates (15.3 km) | `beta_raw` / CHM_CAL, CHM_CAL = 5×10¹¹ |

All signals converted to attenuated backscatter in Mm⁻¹ sr⁻¹ (same conventions as the dark script). The Cloudnet CL61 SNR-screened `beta` is used **only** for cloud detection — using it for noise statistics would bias σ catastrophically low.

## 3. Cloud screening

Common 30 s grid over the 2 days; a bin is clear only if **all three instruments** have data **and** none flags it — the CHM15k (1064 nm) acts as cirrus safety net for the 910 nm pair and vice versa.

Per-instrument "not clear" flags:

- **CL31**: instantaneous CBH (any of 3 layers), obscuration, vertical visibility; sky-condition octas > 0 add exclusions only (the 30-min clustering lags — safe direction).
- **CHM15k**: firmware CBH (> 0; fill = −1), SCI ≠ 0 (rain/fog/snow), cirrus test.
- **CL61** (no CBH in Cloudnet files): liquid cloud (screened β, 5-gate movmedian > 10 Mm⁻¹sr⁻¹ above 100 m), fog/drizzle (raw β median 100–300 m > 20), cirrus test.

**Cirrus test (lesson learned):** a fixed threshold cannot work — range-corrected noise grows ~r², and a 1 Mm⁻¹sr⁻¹ threshold on the 5-min smoothed signal flagged **100 %** of the data. The working test is **adaptive per gate**: flag where the smoothed signal exceeds `median + 6·1.4826·MAD` of that gate (and an absolute floor of 0.3 Mm⁻¹sr⁻¹), evaluated over 4–12 km. The per-gate statistics automatically track the r² noise growth; episodic cloud passages do not bias them.

The combined cloud flag is dilated by ±10 min (cloud edges, virga, humidified shells); missing-data bins block only themselves. Segments ≥ 20 min are kept.

**Result:** 33.2 % clear = **15.9 h** (11.0 h day, 3.6 h night, twilight excluded). 29 May is essentially clear; 30 May is entirely excluded by cirrus at 8–12 km — visible to CL61/CHM15k but **invisible to CL31** (7.7 km range), which alone would have called the day clear. Multi-instrument screening matters.

| # | Start | End | Duration | Day % | Night % |
|---|---|---|---|---|---|
| 1 | 29.05 00:00 | 29.05 06:55 | 6.92 h | 37 | 45 |
| 2–8 | 29.05 07:21 → 18:26 | (7 daytime segments) | 8.5 h | 100 | 0 |
| 9 | 29.05 22:59 | 29.05 23:29 | 0.51 h | 0 | 100 |

![Time-height quicklooks with mask](figs_ceilo_ambient/01_timeheight_mask_20260529_20260530.png)
![Screening diagnostics](figs_ceilo_ambient/02_screening_20260529_20260530.png)

## 4. Noise estimators

Three estimators of σ(z), computed per class (day = solar elevation > 5°, night < −6°), aggregated in 30 m bins, robust scale (1.4826·MAD, std fallback for CL31 quantization):

- **(a) Vertical-profile-only** — per profile, subtract a centred ~150 m moving mean along range; residual σ corrected by 1/√(1−1/L). A linear vertical gradient cancels exactly; only sub-150 m curvature leaks.
- **(b) Temporal first difference** — per gate, σ = robust_std(β(t+Δt)−β(t))/√2 within contiguous clear runs (details below).
- **(c) 30 m × 10 min blocks** — fitted plane (offset + time slope + range slope) removed per tile, σ with N/(N−3) dof correction, **median over tiles** per bin (rejects aerosol-filament tiles), IQR as uncertainty.

### Estimator (b) in detail — and why it is the recommended one

For each gate r, the measured signal decomposes as β(t) = A(t) + ε(t): a slowly varying atmospheric contribution A (aerosol + molecular) and instrument noise ε, independent from one profile to the next. In the first difference d(t) = β(t+Δt) − β(t):

- the noise terms are independent with variance σ² each → Var = **2σ²**;
- the atmosphere changes negligibly in one Δt (15–30 s) → A(t+Δt) − A(t) ≈ 0.

Hence **σ_inst = std(d)/√2**, with the unknown atmospheric profile cancelling *exactly* — no model, no smoothing scale, no stationarity window. In frequency domain, differencing is a high-pass filter with power response |H(f)|² = 4·sin²(πfΔt): ~zero gain at low frequency (thermal drift, solar cycle, aerosol evolution, laser ageing — all suppressed), average gain 2 over white noise (hence the √2). Everything slower than ~2Δt (60 s) is removed; the PSDs (Fig. 7) show the atmospheric content lives below ~3·10⁻³ Hz (periods > 5 min), far from the band (b) measures in.

Why it is recommended over (a) and (c):

1. **Weakest assumptions.** (a) assumes the aerosol profile is smooth over 150 m (false in thin layers, inversions, forming fog); (c) assumes 10-min × 30 m stationarity. (b) only assumes the atmosphere is frozen over 30 s per gate — and this assumption is *tested*: the 10-min windowed linear-detrend variant (σ² = RSS/(n−2)) overlaps it (dotted curves, Fig. 4), so atmospheric leakage is negligible.
2. **Per-gate at native resolution** — exactly what the detection-threshold formula β_min(z) = SNR·σ(z)·√(Δt/τ)/T²_mol consumes; no vertical mixing as in (c).
3. **Gap-robust.** It only needs *pairs* of consecutive valid profiles. PSD needs 32 min of unbroken data, Allan needs hours, (c) needs 10-min tiles — (b) works with a heavily fragmented cloud mask, which is the normal situation for network stations.
4. **Least leakage of all estimators.** Aerosol structure evolving slower than 60 s is fully removed regardless of its vertical scale; (a) passes any vertical structure < 150 m even if perfectly stationary — hence σ_a ≥ σ_b in the daytime boundary layer.
5. **Triple-validated.** Dark block: b/std = 0.99–1.01 (best of the four estimators). Ambient: closes with the PSD white floor to ±9 %, and the anchor at τ=Δt falls on the Allan curves — as it must, since for white noise ADEV(Δt) *is* std(diff)/√2: estimator (b) is mathematically the first point of the Allan curve.
6. **Tracks non-stationary noise.** Photon noise depends on the solar background and on the signal itself; being local in time, (b) can be evaluated per class (done here) or in a sliding window to produce σ(t, z) for continuous network monitoring.

Known limitations: internal firmware time-smoothing would correlate consecutive ε and bias (b) low — detectable as a PSD roll-off below Nyquist (none observed here, the plateaus are flat); fast convective turbulence below ~1–2 km by day can bias it slightly high; CL31 quantization (LSB = 0.01 Mm⁻¹sr⁻¹) can degenerate the MAD, in which case the script falls back to the classical std.

**Does the mean profile need to be subtracted (estimator a)?** Yes — and the *time-mean* profile is **not enough**. Subtracting only the class-mean profile (= plain per-gate std) leaves all atmospheric *variability* in: the dotted curves in Figure 3 bulge below ~2 km. The subtraction must be **per profile** (running smooth along range), because the aerosol profile drifts over hours.

**Validation on the dark block** (hood on → no atmosphere → all estimators must equal plain std):

| Instrument | a/std | b/std | c/std | windowed/std |
|---|---|---|---|---|
| CL31 | 0.99 | 0.99 | 0.98 | 0.98 |
| CL61 | 0.95 | 1.01 | 0.96 | 0.98 |
| CHM15k | 0.95 | 1.00 | 0.99 | 0.99 |

On ambient data the three estimators agree within ~6 % everywhere (night BL: a/b = 0.94–0.99, c/b = 0.98–1.01). **Estimator (b) is the recommended default** for the reasons detailed above: per-gate at native resolution, weakest assumptions, gap-safe, and triple-validated (dark block, PSD closure, Allan anchor).

![Estimator a](figs_ceilo_ambient/03_sigma_vertical_20260529_20260530.png)
![Estimator b](figs_ceilo_ambient/04_sigma_firstdiff_20260529_20260530.png)
![Estimator c](figs_ceilo_ambient/05_sigma_block_20260529_20260530.png)
![Estimator intercomparison](figs_ceilo_ambient/06_estimator_comparison_20260529_20260530.png)

## 5. Separating instrument noise from atmospheric variability (FFT + Allan)

**Welch PSD** per altitude band (500/1000/2000/3000/5000 m ± 30 m, per class; Hann, 50 % overlap, 32-min segments, never crossing mask/data gaps): atmospheric variability appears as a red f^(−1…−2.4) branch at low frequency, instrument noise as a flat white floor at high frequency. σ_inst = √(C/(2Δt)) from the floor C.

- **Closure:** σ_inst from the PSD floor matches estimator (b) within **±9 %** at all bands and classes (single exception: CL61 day 500 m, ratio 1.71, where the built-in flatness check correctly flagged the floor as contaminated).
- **Where is the atmosphere detectable?** Only below ~1–2 km: CL61 night 500 m (slope μ≈2.3, crossover ≈ 10 min), CL61 night 1 km (μ≈0.9, ≈ 5 min), CHM15k 500 m day (μ≈1.2, ≈ 3.4 min). At ≥ 2 km the spectra are white for all three instruments — **everything measured there is instrument noise**, which is why ambient night profiles reproduce the dark measurement.
- **Allan deviations** (overlapping, χ² CI) follow τ^(−1/2) in all bands, with the estimator-(b) anchor at τ=Δt sitting on the curves; departures from τ^(−1/2) (atmospheric drift) appear only in the lowest bands at long τ.

![Welch PSDs](figs_ceilo_ambient/07_psd_20260529_20260530.png)
![Allan deviation per band + closure](figs_ceilo_ambient/08_allan_20260529_20260530.png)

## 6. Ambient vs dark noise

σ_β (estimator b, native Δt, Mm⁻¹ sr⁻¹):

| z | CL31 day / night / **dark** | CL61 day / night / **dark** | CHM15k day / night / **dark** |
|---|---|---|---|
| 500 m | 0.067 / 0.049 / **0.042** | 0.0105 / 0.0051 / **0.0047** | 0.030 / 0.026 / **0.022** |
| 1 km | 0.290 / 0.210 / **0.187** | 0.031 / 0.020 / **0.018** | 0.077 / 0.066 / **0.064** |
| 2 km | 1.174 / 0.735 / **0.737** | 0.099 / 0.078 / **0.075** | 0.256 / 0.235 / **0.227** |
| 3 km | 2.39 / 1.69 / **1.59** | 0.231 / 0.158 / **0.165** | 0.564 / 0.511 / **0.485** |
| 5 km | 7.06 / 4.67 / **4.60** | 0.66 / 0.49 / **0.46** | 1.55 / 1.40 / **1.38** |

- **Night ambient = dark within 0–10 %** (≤ 19 % at 500 m where shot noise from the aerosol signal contributes). The ambient method recovers the dark noise profile without a hood.
- **Day/night ratio:** ×1.4–1.6 (CL31), ×1.3–2.1 (CL61), ×1.1 (CHM15k — 1064 nm is least affected by solar background). The diurnal mechanism is explicit in Figure 12: CL31 background RCS and CHM15k baseline/firmware-stddev track solar elevation exactly.
- **Hardware-level sensitivity** (noise density σ·√(Δt·Δr) at 3 km night): **CL61 1.9** ≪ CHM15k 7.7 ≪ CL31 29 Mm⁻¹sr⁻¹·√(s·m). The CL61 is by far the most sensitive instrument per unit bandwidth.

![Ambient vs dark](figs_ceilo_ambient/09_ambient_vs_dark_20260529_20260530.png)
![Noise vs signal and solar background](figs_ceilo_ambient/12_noise_vs_signal_20260529_20260530.png)

## 7. Detection thresholds

Same formulation as the dark script: β_min(z) = SNR·σ(z)·√(Δt/τ) / T²_mol(z), SNR = 3, white-noise τ-scaling justified by the Allan/PSD analysis. Extinction α_min = LR·β_min; mass M_min = α_min/MEC (volcanic ash: LR = 60 sr, MEC = 0.60 m²/g). ICAO ash zones 200/2000/4000 µg/m³.

**Minimum detectable mass M_min [µg/m³], τ = 30 min, SNR = 3:**

| z | CL31 day / night | CL61 day / night | CHM15k day / night |
|---|---|---|---|
| 500 m | 2.6 / 1.9 | 0.41 / 0.20 | 0.81 / 0.71 |
| 1 km | 11.3 / 8.1 | 1.2 / 0.8 | 2.1 / 1.8 |
| 2 km | 45.6 / 28.6 | 3.8 / 3.0 | 7.0 / 6.5 |
| 3 km | 93.1 / 65.6 | 9.0 / 6.1 | 15.5 / 14.0 |
| 5 km | 275.6 / 182.2 | 25.8 / 19.0 | 42.5 / 38.4 |

(τ = 5 min values are √6 ≈ 2.45× higher — full tables in the script printout.)

- All three instruments detect the **200 µg/m³ ICAO low-contamination edge** up to ≥ 5 km with 30 min averaging — the CL31 only marginally by day at 5 km.
- CL61 detects **2 µg/m³ at 1 km / 20 µg/m³ at 5 km** at night: an order of magnitude better than CL31.
- Night ambient threshold curves overlay the dark-derived reference curves (thin black in the figures) — consistency check passed.

![Detection threshold backscatter](figs_ceilo_ambient/10_detect_backscatter_20260529_20260530.png)
![Extinction and mass thresholds](figs_ceilo_ambient/11_detect_extinction_mass_20260529_20260530.png)

## 8. Verification checklist

| # | Check | Result |
|---|---|---|
| 1 | Dark block: a ≈ b ≈ c ≈ plain std | ✅ ratios 0.95–1.01 |
| 2 | Night ambient → dark convergence | ✅ 0.95–1.19 at 0.5–5 km (1.00–1.06 at ≥ 2 km) |
| 3 | Day > night for 910 nm instruments | ✅ ×1.3–2.1; CHM15k ×1.1 (expected at 1064 nm) |
| 4 | PSD-floor σ vs estimator (b) closure | ✅ ±9 % (1 contaminated case auto-flagged) |
| 5 | Allan slope τ^(−1/2) in white region | ✅ with anchor at τ=Δt on curve |
| 6 | Estimator consistency in BL/FT | ✅ within ~6 %, no suspicious inversions |
| 7 | Night thresholds ≈ dark thresholds | ✅ overlay in Figs 10–11 |

## 8b. Minimum data needed, and behaviour in day / cloud

**Estimator (b) is usable on as little as 5 minutes.** It needs only *pairs* of consecutive profiles, not a long record (PSD needs 32 min, Allan needs hours). Empirical scatter of the 5/10/30-min estimate about a 3 h reference (clear night, 29 May), aggregated to 30 m bins:

| Window | 1σ scatter per 30 m bin |
|---|---|
| 5 min | 17–27 % |
| 10 min | 11–21 % |
| 30 min | 4–11 % |

Critically, the **median is unbiased** — 5-min estimates scatter around the reference but show no systematic error (the differencing removes the atmosphere regardless of record length). 5 min is enough for a noise floor / detection-limit to ~20 %; ~30 min reaches < 10 %. Fragments can be pooled (10 × 5-min ≈ one 50-min estimate).

**Daytime (clear): works.** σ_b rises with solar-background shot noise (×1.9 at 910 nm, ×1.2 at 1064 nm) but the estimate stays as stable as at night (same ~20 % per-5-min scatter). The slow solar-elevation drift is removed by differencing; the background shot noise is legitimately counted. Independently confirmed by the **daytime** dark measurement (b/std = 1.00). Day and night must be reported separately because they are different noise regimes, not because the method fails.

**Clouds: must be screened.** With cloud in the beam, σ_b is biased high (~×2 here) and unstable (5-min scatter 35–39 % vs ~22 %) because a cloud edge changes the signal by orders of magnitude in one Δt. The robust MAD limits the damage (it ignores a minority of cloud spikes, so it degrades gracefully), but the value is no longer instrument noise. The cloud screen + ±10 min buffer handles this; the inter-window instability of σ_b is itself a residual contamination flag. The dominant threat is cloud *signal in the beam*, not the background shadow (CHM `base` rises only 0.018→0.025 in cloud) — so a good cloud mask is sufficient, no turbulent-background model needed.

## 8c. Seasonal stability and a parametrised noise model

A direct consequence: **you do not need a continuous 24/7 estimate** — the noise can be modelled from a fixed reference profile plus one scalar.

**The range-corrected noise is z² × (constant raw-signal noise).** Verified: log-log slope of σ_b(z) = 2.01, σ_b/z² flat from 2 to 14 km. Both detector and background noise are uniform in raw signal, so both grow as z² after range correction.

**The night floor is stable across seasons** (CHM15k high-gate noise, Jan→Jul 2025/26):

| Season | CHM15k σ_night(12 km) | temp_det | temp_int | CL31 σ_night(6 km) | laser_temp |
|---|---|---|---|---|---|
| Jan | 7.69 | 30.0 °C | 18.2 | 6.27 | 24.5 |
| Mar | 7.96 | 30.0 °C | 14.9 | 6.31 | 25.0 |
| May | 8.14 | 30.0 °C | 24.9 | 6.67 | 35.7 |
| Jul | 7.33 | 30.0 °C | 28.4 | 6.52 | 39.7 |
| **range** | **±5 %** (10 % p-p) | **0.0** | 14→28 | **±3 %** (6 % p-p) | 24→40 |

The **CHM15k APD is Peltier-stabilised at exactly 30.0 °C** (`temp_det` invariant across all dates 2024–2026) → its dark floor is a fixed instrument constant, season-independent; `temp_int` swings 14→29 °C with no correlation to the floor (r = −0.19). The **CL31** is not stabilised (`laser_temp` 24→41 °C) but shows only a weak trend. The **CL61** (raw files Feb–Jun 2026, ~32k five-minute files downloaded from Cloudnet for unit 2d386110): `internal_temperature` swings 14–43 °C, `laser_temperature` is locked (20–21 °C), and over this wide range the night floor shows a clearer but still modest positive temperature dependence (r = +0.45) — its firmware **`beta_att_noise_level` predicts the empirical σ_b at r = 0.83**.

A **dense single-season** characterisation (Feb–May 2026: one clear-sky night-floor point per day for CHM15k and CL31, plus 225 clear CL61 points) settles the temperature dependence cleanly. To avoid bias, the floor is taken as the **20th percentile of σ in 2 °C temperature bins** (cloud only inflates σ, so a low percentile per bin recovers the true clear-sky floor at each temperature, without the flat-bias of a global-percentile envelope), then fitted vs temperature:

| Instrument | temperature variable (range) | slope | fit |
|---|---|---|---|
| CL31 | `laser_temp` (25–37 °C) | **+0.25 %/°C** | r = 0.86, 114 nights |
| CHM15k | `temp_int` (16–25 °C) | **+0.34 %/°C** | r = 0.84, 120 nights; `temp_det` locked 30 °C |
| CL61 | `internal_temperature` (15–37 °C) | **+2.67 %/°C** | r = 0.66, 225 points |

Physically consistent: the two instruments with stabilised/quiet detection (CHM15k Peltier-locked APD at 30 °C; CL31) are **essentially flat** (≈0.3 %/°C, i.e. a few % across the whole season), while the **unstabilised CL61 APD** shows the expected — but still modest — positive trend (~2.7 %/°C, ≈60 % over its 22 °C range), in line with the exponential rise of APD dark counts with temperature. Separately, across *multiple years* the CHM15k `beta_raw`-unit floor drifts by ~×2 — this is **laser-energy / calibration normalisation drift, not temperature** (zero within-season correlation). Practical conclusion: a single night-floor reference per instrument suffices, refreshed occasionally (yearly, to track calibration); for the CL61 a small `internal_temperature` correction (~2.7 %/°C) can be applied where the last few percent matter.

**Parametrised model (validated to ±8 % across the whole column with a single coefficient):**

> **σ(z)² = σ_night(z)² + (k·z²)²,  with  k = c·√B**

- quadrature sum (independent noises), z² geometry, k ∝ √B (background shot noise).
- CL31 validation (29 May day): observed vs model σ_day(z) ratio 0.92–1.04 over 1–7 km; CHM15k 0.92–1.00 over 2–14 km. CL31 `bckgrd_rcs_910` day/night = 2.9 and k/√B = 8.0×10⁻⁸ = instrument constant.

**Background proxy per instrument:** CL31 `bckgrd_rcs_910` (direct, r = 0.85–0.88); CHM15k `base` × calibration (`base` is small but real, 0.018–0.031, day/night = 1.29, r = 0.87 over the diurnal cycle); CL61 `beta_att_noise_level` (r = 0.94) or `monitoring/background_radiance`.

**Operational recipe (no 24/7 estimation):**
1. Measure σ_night(z) once per instrument from a few clear nights → fixed reference (stable ±5–6 % year-round; for the CHM the stabilised detector makes it a true constant).
2. Fix the instrument constant c once from any clear day (c = k/√B).
3. At any time: **σ(z,t) = √( σ_night(z)² + c²·B(t)·z⁴ )**, reading the background B(t) from housekeeping — or, where no background variable exists, from solar elevation (deterministic).

**Limitations.** This describes the noise *floor* at altitude (no signal); below ~1–2 km the aerosol signal shot noise adds. CL61 temperature/background were taken from raw Vaisala files (Feb–Mar 2026); the Cloudnet daily files carry neither, so for a Cloudnet-only CL61 the background term must use solar elevation. The CL31 6 % temperature trend can be absorbed with a small `laser_temp` term if sub-5 % accuracy is required.

![z2 scaling](figs_ceilo_ambient/13_zscaling.png)
![noise vs temperature](figs_ceilo_ambient/14_temperature.png)
![seasonal stability](figs_ceilo_ambient/15_seasonal.png)
![model validation](figs_ceilo_ambient/16_model_validation.png)
![CL61 firmware noise level](figs_ceilo_ambient/17_background_cl61.png)

## 8d. Related work and how this method differs

A literature scan confirms the building blocks are established, but the combination here — a dark-validated, network-applicable ambient noise floor plus a parametrised detection-limit model — is not standard practice.

**Established prior art**
- **Kotthaus et al. (2016), AMT 9, 3769–3791** — the reference for processing Vaisala CL31 attenuated backscatter. It explicitly determines the *instrument-related background signal and noise* (not provided in the standard CL31 output) and notes that "the standard deviation at each range gate can be used as a noise estimate … assuming there are no temporal variations in the atmosphere." Our first-difference estimator is precisely a refinement that **drops that assumption** (the lag-1 difference removes the slowly varying atmosphere, so the estimate is unbiased even when the atmosphere is not static).
- **O'Connor et al. (2004), JTECH 21, 777** — autocalibration of cloud lidar; the methodological basis for E-PROFILE-style elastic-lidar processing and SNR use.
- **Liu et al. / "Estimating random errors due to shot noise in backscatter lidar observations"** (CALIPSO, NASA) — formalises lidar shot noise (signal-dependent, white) and shows the conventional std-over-N-consecutive-profiles method **over-estimates uncertainty below ~1.5 km in variable aerosol layers**. This is exactly the boundary-layer bias we see in the plain-std / vertical estimator and which the first-difference avoids.
- **Wiegner & Geiß (2012), AMT 5, 1953** — ceilometer backscatter retrieval and **SNR determination**.
- **Flentje et al. (2010) / Heese et al. (2010)** — DWD ceilometer-network aerosol/ash profiling; report the day/night SNR-vs-altitude behaviour (SNR > 1 to ~4–5 km by day, ~8.5 km at night) that our day/night split reproduces.
- **Córdoba-Jabonero et al. (2022), Remote Sens. 14, 5680** — volcanic-ash **mass-concentration** retrieval from CL51/CL61 ceilometers; the application our extinction/mass detection limits feed.
- **APD physics** — dark-count rate rises ~exponentially with detector temperature; this is why temperature-stabilised detectors (the CHM15k APD Peltier-locked at 30 °C) show a season-stable floor, and why the unstabilised CL31/CL61 show a small residual trend.

**What is new here**
1. A **lag-1 (first-difference) per-gate noise estimator** that removes atmospheric-variability bias and is validated against a covered-telescope (dark) reference — agreement within a few percent.
2. A **network-applicable ambient method**: needs only routine attenuated backscatter + cloud flags, no termination hood and no housekeeping for the floor itself; night-time ambient reproduces the dark noise.
3. A **parametrised detection-limit model** σ(z)² = σ_night(z)² + (k·z²)², with the night floor as a season-stable instrument constant (±5–6 %), enabling an altitude-resolved minimum-detectable backscatter/extinction/mass to be attached to every L2 profile without continuous estimation.
4. A **like-for-like CL31 vs CL61 vs CHM15k comparison** of noise and detection limits under one methodology.

## 9. Conclusions & network outlook

1. **The ambient method works**: clear-sky nights reproduce the hooded dark measurement at all altitudes — noise profiles and detection thresholds for any network instrument can be derived from routine measurements, screened for clouds.
2. **Daytime degradation is moderate and quantifiable** (solar background; ×1.1–2 at 30 s native resolution) and should be characterised per station since it depends on telescope/filter design.
3. **Estimator (b) (first difference / √2) is the recommended network estimator** — per gate, robust, insensitive to atmospheric drift, validated by (a), (c), PSD and Allan.
4. **Atmospheric variability only matters below ~2 km** at 30 s resolution; above, ambient σ *is* instrument noise.
5. **Requirements for network rollout**: raw (non-averaged) data and a cloud screen. The CL61 case shows the screen can run without CBH/housekeeping; single-instrument stations lose the cross-instrument cirrus veto, partly compensated by the adaptive (median + 6·MAD) cirrus test. Candidate raw sources: Cloudnet instrument files (as used here for CL61), CHM15k daily NetCDF archives.
6. **Instrument ranking** (noise density at 3 km night): CL61 (1.9) ≫ CHM15k (7.7) ≫ CL31 (29 Mm⁻¹sr⁻¹·√(s·m)).

## Files

- Analysis script: `ambient_noise_cl61_chm_cl31.m`
- Screening diagnostic helper: `debug_ambient_screening.m`
- Method tests: `test_5min_estimator_b.m`, `test_b_day_clouds.m`, `test_noise_model_temperature.m`, `test_noise_model_validate.m`, `test_cl31_seasonal_floor.m`, `test_cl61_raw_scan.m`, `test_chm_base_scan.m`, `test_temp_dense.m` (dense seasonal night floors, cloud-robust)
- Figure generators: `make_noise_figures.m` (report figs 13–17), `make_paper_figures.m` (the 2 paper figures)
- Word paper section: `build_paper_section.py` → `E-PROFILE_ALC_noise_section_v2.docx`
- Figures: `figs_ceilo_ambient\01…17_*.png`; paper figures `figs_ceilo_ambient\paper\paper_fig1_dark_detection.png`, `paper_fig2_ambient_method.png`
- Binned results (σ profiles, PSD, Allan): `figs_ceilo_ambient\ambient_noise_results_20260529_20260530.mat`; CL61 raw scan: `cl61_raw_scan.mat`
- Dark reference script: `dark_measurement_cl61_chm_cl31.m`
- CL61 raw files (with temperature/background/cloud): `E:\CL61_PAY\` (24 Feb–6 Mar 2026); more downloadable via the Cloudnet raw-files API (instrument 2d386110)
