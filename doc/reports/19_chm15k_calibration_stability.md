# 19 — CHM15k calibration stability: temperature dependence and diurnal behaviour

Date: 2026-10-20. Branch `remote-dark`. Plan: [19_chm15k_calibration_stability_plan.md](19_chm15k_calibration_stability_plan.md).
Code: `validation/chm15k_stability/`. Outputs: `C:/DATA/Projects/202610_CHM15k_stability/`.

## Summary

> **Notation.** **κ** is the temperature coefficient of the lidar constant, κ = d ln C_L / dT_int
> (%/K), and **κ_u** is its value for one unit (optical module). β is reserved for the backscatter
> coefficient (β_att = attenuated backscatter).

**Q1. Does the CHM15k Rayleigh constant depend on temperature? Yes.** It rises with the instrument's
**internal temperature** `temp_int`. The effect survives every control for season and atmosphere:

| Estimator | Units | Slope d ln C_L / dT_int | Significant (>3σ) |
|---|---|---|---|
| Network 2025–26, 31-day anomaly model (M2) | 158 unit-eras | **+2.06 ± 0.04 %/K** (median +2.08) | 100 positive / 0 negative |
| Same + aerosol & window-height controls (M2x) | 151 | **+2.00 ± 0.05 %/K** | 83 / 0 |
| 11-year archive 2013–2024, M2 | 22 unit-eras | **+1.43 ± 0.09 %/K** | 12 / 0 |
| Within a night, elapsed time held fixed | 150 units | **+1.46 ± 0.02 %/K** (median +1.31) | 112 / 0 |
| Same-wavelength pairs, night, vs ΔT_int between units | 8 pairs | +0.4 … +2.1 %/K | atmosphere cancels |

The coefficient is **unit-specific**, from about 0 to +3.5 %/K. It **grows with the optical-module
generation**: about +1.0 %/K for TUB06–TUB15 modules and +2.0 to +2.7 %/K for TUB20–TUB24.

**Q2. Is the constant free of diurnal / short-term patterns? No.** The internal temperature rises by
day: median +3.7 K at noon, +5.7 K in summer, up to ~10 K. The effective constant rises with it, by
**+4.8 % at noon (median), +7.5 % in summer**. A third of the units exceed +10 % on summer noons.
The Rayleigh calibration is night-only, so it never sees this, and daytime L2 attenuated backscatter
of the CHM15k is biased high by that amount.

The night-time temperature coefficient alone predicts the observed daytime behaviour:
* co-located CHM15k pairs: r = 0.98, observed = 1.16 × predicted;
* CHM15k/CL61 at three of four sites (Aosta +26 % predicted vs +29 % observed; Lindenberg +8 vs
  +12 %; Camborne +5 vs +6 %).
* Payerne is the exception, now explained (§7.1): its receiver is abnormally noisy (night background
  5–25× the network), so the night κ is unreliable. Its optical module also loses regulation on 37 % of
  daytime hours, with an opposite-sign effect.

Within a night, once temperature is accounted for, the constant is stable. There is no residual drift
(−0.08 %/h), and the hourly scatter is close to photon noise (2.5 % excess).

**Mechanism.** The effect is a **multiplicative, height-independent gain change in the raw photon
counts**. Excluded: the detector, the internal reference pulse, an additive electronic baseline,
optics/overlap, and **aerosol** (AERONET control: +2.36 → +2.01 %/K with AOD + Ångström, +2.13 %/K on the
cleanest nights; §7.4). The responsible component follows the **housing** temperature while the
regulated optical module stays at 33.85 °C.
* **Laser pulse energy is unlikely.** In theory it rises with crystal temperature (+0.2…+0.7 %/K), but the
  crystal is regulated, and the pulse rate moves only −0.06 %/K, uncorrelated with κ.
* **Leading hypothesis: the laser line vs the 1 nm receiver filter**, with a thermal alignment drift
  as a co-equal rival. The manual puts the filter inside the heated laser module, so only Lufft
  information can decide (§7.8, twelve questions).
* `temp_lom` adds a unit-specific term only when the module leaves regulation on hot afternoons (§7.3).

**The daytime drop of the internal calibration pulse** (−15 % at noon, −24 % in summer, tracking the
solar background at r = −0.96) is a real detector sensitivity loss. The firmware's division of the
signal by `p_calc` compensates it: the CHM15k/CL61 ratio does not follow it at Camborne. It is not
the driver of the temperature effect.

**Recommendation.** Add a per-unit temperature term to the CHM15k calibration chain:

C_L(t) = C_L,ref · exp(κ_u · (T_int(t) − T_ref))

* κ_u is estimated per unit from the intra-night method (≥ ~50 clear nights; generation median as a
  fall-back).
* Apply it (i) to normalise the nightly constants before the Kalman filter, and (ii) per profile in
  the L1→L2 step.

**Is it better than the housekeeping ML?** Yes, on the 2024 test year of the 11-year archive (§7.5):
* housekeeping only: exp model 20.2 % vs bagged-tree ML 22.3 %;
* operationally, as the 30-day running median of temperature-normalised constants × exp(κ ΔT):
  **13.0 %**, against 13.6 % for the plain running median and 14.2 % for the ML. Best at 9 of 11
  stations.

**Does it improve agreement with independent instruments? Yes, at 6 of 7 comparisons** (§7.6):
* the day−night bias against CL61 and EARLINET collapses from 7–23 pts to ≤ 4.5 pts;
* the spread drops at Aosta (27 → 18), Lindenberg (27.5 → 24) and Măgurele (41 → 37, 35 → 31.5);
* Payerne is over-corrected, as expected from §7.1.

Section 6 gives the gains and the open points.

## 1. Data and method

* **Housekeeping + signal cache** (`hk_cache.py`). The L1 files of all 154 CHM15k streams of the
  census, 2025-01 → 2026-06 (`D:/E-PROFILE_L1_2026`), are reduced to 10-min blocks:
  * medians of the housekeeping variables: `temp_int`, `temp_ext`, `temperature_detector`,
    `temperature_optical_module`, `calibration_pulse`, `window_transmission`, `status_laser`,
    `status_detector`, `laser_life_time`, `laser_pulses`, `bckgrd_rcs_0`, `stddev`;
  * clear-profile counts (no cloud base, `sci` = 0);
  * mean clear-sky `rcs_0` in range bins (100 m to 1 km, 250 m to 8 km);
  * `optical_module_id` per day.

  The cache holds 0.88 GB. The four CL61 co-located with a CHM15k are cached the same way.
* **Constants.** The operational nightly Rayleigh constants (dashboard NetCDF, E-PROFILE v2.2,
  `calibration_method` 0) are matched to the cache over each night's calibration window.
  Each unit is split into eras at optical-module changes.
* **11-year archive.** `A:\CHM15k\<stn>\chm_<stn>_1d_2013-2024.mat` (13 stations), the dataset behind
  `analyze_chm15k_stability.m` / `estimate_calib_from_housekeeping_v4.m`, exported to CSV
  (`export_11yr_archive.m`) and split at `serlom` changes.
* **Season control.** `temp_int` is ~70 % annual cycle (median seasonal fraction 0.69), so a naive
  regression mixes temperature with season (aerosol in the window, night length, SNR). Three models:
  * **M0** naive;
  * **M1** with an annual harmonic and a trend;
  * **M2** on anomalies against a 31-day running median. M2 keeps synoptic warm/cold spells only.
* **Within-night test** (`q2a`). The clear-sky signal in each night's own Rayleigh window is
  followed hourly. Over one night the molecular term varies by < 1 %, so S(h)/S_night = C(h)/C_night.
  The fit is joint in ΔT_int and elapsed time, because the instrument cools steadily through the
  night (r(T, h) ≈ −0.4 to −0.6).
* **Same-wavelength pairs** (`q2b`): Schiphol 06240 A/B/C/D (4 units, a few km apart; the L1 and L2
  files carry one common station position), Budapest 12843 A/B, Măgurele INO A/B.
  * The atmosphere cancels in the clear-sky signal ratio at 1.0–1.5 / 1.5–2.5 / 2.5–4.0 km.
  * Anomalies are taken against the night-time ratio within ±15 days, within the same module pair.
* **CHM15k vs CL61** (`q2d`) at Payerne, Lindenberg, Camborne and Aosta. The CL61's analog receiver
  gives an absolute daytime reference. It is affected by the 1064/910 nm Ångström difference and by
  water vapour at 910 nm, so it is sensitive only to effects of ≳ 10 %.

## 2. Q1 — the constant depends on the internal temperature

The raw scatter (top row) mixes the season with temperature. The anomaly model (bottom row) keeps only
warm/cold spells against the local seasonal state, and the positive slope remains:

![Raw and anomaly view of C_L vs T_int at five units](figs_chm15k_stability/q1_examples_raw_vs_anomaly.png)

![Distribution of the slopes per unit-era, network and 11-year archive](figs_chm15k_stability/q1_slope_histograms.png)

![Per-unit slopes, models M0 / M1 / M2](figs_chm15k_stability/q1_slopes_per_unit.png)

* **Removing the season shrinks the naive slope by about a third**:
  * network: +2.78 → +1.94 (M1) / +2.06 (M2) %/K;
  * archive: +2.33 → +1.17 / +1.43 %/K.

  So report 09's "r = +0.38 with internal temperature" (`analyze_chm15k_stability.m`, raw correlation)
  was part season, but most of the dependence is real.
* **Not haze.** Warm spells are hazier: r(T anomaly, low-level aerosol contrast) = +0.22. Adding the
  anomalies of the aerosol contrast ln(S_1–2.5 km / S_window) and of the window bottom height leaves the
  slope unchanged (+2.06 → +2.00 %/K).
* External temperature gives about half the slope (+1.04 %/K, M2). The relevant variable is the
  internal one.
* **Detector and optical-module temperatures** never vary by 5 K (both regulated; detector fixed at
  24.85 °C), so no slope can be fitted on them. The optical module never leaves regulation at night:
  `lom_exc_frac` = 0 for every unit. The "regulation loss" of Roşu et al. (2026) is an afternoon
  phenomenon on hot days.
* Payerne (TUB140016, 136 nights) is too noisy to resolve a slope on its own:
  * operational constant: M2 −2.2 ± 1.5 %/K;
  * constant re-run with the measured hood dark subtracted: +0.5 ± 0.9 %/K.

## 3. Q2a — within a night the constant follows the temperature, nothing else

![Intra-night slopes, drift and excess scatter per unit](figs_chm15k_stability/q2a_intranight_summary.png)

![Intra-night C(h)/C_night against T_int change, examples](figs_chm15k_stability/q2a_intranight_examples.png)

* Per unit, the within-night slope is median +1.41 %/K. With elapsed time held fixed it is
  **+1.31 %/K** (pooled +1.46), and **112 of 150 units are significant, none negative**.
* The apparent drift through the night (−0.25 %/h) is the instrument cooling. With temperature in the
  model it drops to −0.08 %/h, so aerosol settling out of the window is not the explanation.
* The hourly scatter of C(h)/C_night (median 5.2 %) is mostly photon noise (4.3 %). The excess is
  2.5 %. Apart from temperature, the constant has no short-term variability worth modelling.
* The within-night coefficient is smaller than the multi-day one (median +1.35 vs +2.08 %/K;
  per-unit r = +0.55). That fits a slow thermal response: the `temp_int` sensor reacts faster than the
  responsible component, so hour-scale swings are only partly transmitted.

## 4. Q2b — the daytime cycle

### 4.1 The internal calibration pulse drops by day (but is compensated)

![Diurnal cycle of the housekeeping, network](figs_chm15k_stability/q2c_hk_diurnal.png)

* At every unit, `calibration_pulse` falls by day: median −15 % at noon (−24 % in summer; p10–p90
  −34 … −12 %). It correlates with the solar background (r = −0.96 median) and not with `temp_int`
  (r = +0.08).
* The CHM15k firmware computes `beta_raw = ((P/lp) − b) / (cs · o(r) · p_calc) · r²`
  ([07_actris_chm15k_beta_att.md](paper_v08_review/07_actris_chm15k_beta_att.md)).
* If the drop were an artefact of the reference, the division would inflate the daytime signal, and
  ln(signal) would follow ln(p_cal) with slope −1. Against the CL61 at Camborne the slope is
  +0.18 ± 0.23, and the CHM15k/CL61 ratio stays within a few percent while p_cal falls 17 %.
* So the drop is a real detector sensitivity loss under daylight (photon-counting dead time), and the
  firmware compensates it.
* At the other three sites the regression on p_cal is not interpretable: p_cal and the temperature
  both follow the sun, so the two are collinear.

### 4.2 Same-wavelength pairs: the atmosphere cancels, a diurnal cycle remains

![Diurnal cycle of the signal ratio of co-located CHM15k pairs](figs_chm15k_stability/q2b_pairs_diurnal.png)

| Pair | Day − night (10–16 h solar, 1.5–2.5 km) | 95 % CI | Night slope vs ΔT_int |
|---|---|---|---|
| Schiphol A/B | −14.6 % | [−15.3, −13.0] | +0.5 ± 0.4 %/K |
| Schiphol A/C | +6.7 % | [+6.0, +7.6] | +1.2 ± 0.4 |
| Schiphol A/D | +2.6 % | [+1.4, +3.5] | +0.8 ± 0.3 |
| Schiphol B/C | +21.2 % | [+19.5, +22.6] | +2.1 ± 0.5 |
| Schiphol B/D | +17.2 % | [+15.9, +19.2] | +1.7 ± 0.4 |
| Schiphol C/D | −3.4 % | [−4.2, −2.1] | +0.4 ± 0.4 |
| Budapest A/B | −5.1 % | [−5.3, −4.8] | +2.0 ± 0.2 |
| Măgurele A/B | +0.3 % | [+0.1, +0.5] | +0.4 ± 0.1 |

Four points show the cycle is instrumental:
* The cycles are smooth and centred on solar noon.
* They are **the same in all three height bands**. Above the boundary layer (2.5–4 km) the air is
  horizontally homogeneous, so units a few km apart cannot differ through the atmosphere.
* They are absent in winter.
* The Schiphol pairs close exactly: A/B and A/C predict B/C at +21.3 %, observed +21.2 %.

So unit B reads ~15 % higher by day relative to its night than A, and C ~7 % lower.

![Pair ratio against the internal-temperature difference, night and day](figs_chm15k_stability/q2b_pairs_vs_dT.png)

### 4.3 CHM15k against CL61

![CHM15k / CL61 signal ratio and the CHM15k calibration pulse](figs_chm15k_stability/q2d_cross_cl61.png)

| Site | Day − night CHM/CL61, 1.0–1.5 / 1.5–2.5 / 2.5–4.0 km |
|---|---|
| Payerne | +4.2 / +2.7 / +1.1 % |
| Lindenberg | +11.4 / +11.5 / +9.0 % |
| Camborne | +6.7 / +5.5 / −0.2 % |
| Aosta | +29.1 / +29.1 / +25.6 % (summer +46 %) |

At Lindenberg and Aosta the excess peaks in the afternoon, not at noon. That follows the instrument
temperature, not the calibration pulse.

### 4.4 Synthesis — the night-time coefficient predicts the day

![Synthesis: estimator consistency and daytime prediction](figs_chm15k_stability/q3_synthesis.png)

* For each unit, predicted day−night = κ_u (within-night, elapsed time held fixed) × the rise of
  `temp_int` from night to day over the same clear hours.
* **Pairs:** the predictions fit at r = 0.98, with observed = 1.16 × predicted. For example:
  * Schiphol B/C: +17 predicted vs +21 % observed;
  * Schiphol A/D: +2.6 vs +2.6 %;
  * Budapest: −9.7 vs −5.1 %;
  * Măgurele: −0.7 vs +0.3 %.
* The 16 % excess matches the multi-day coefficient being larger than the within-night one (§3).
  The day-long warming is a longer, slower forcing than night cooling.
* **CL61 sites:** Aosta +26 vs +29 %, Lindenberg +8 vs +12 %, Camborne +5 vs +6 %.
* **Payerne** is the exception: +20 ± 8 % predicted, +2.7 % observed. It is explained in §7.1: a very
  noisy receiver makes its night-time coefficient unreliable, and its optical module loses regulation
  on 37 % of daytime hours with an opposite-sign effect.

![Coefficient by optical-module generation, and the resulting daytime excess](figs_chm15k_stability/q3_generation_daytime.png)

* Across the network, the noon excess of the effective constant over the night calibration is:
  * **median +4.8 %** (p10–p90 +1.1 … +11.2 %);
  * **summer +7.5 %** (+2.1 … +16.5 %);
  * above +10 % for **33 % of units** in summer.

## 5. Mechanism

![Height dependence and housekeeping response of the temperature coefficient](figs_chm15k_stability/q4_mechanism.png)

| Candidate | Fingerprint tested | Result | Verdict |
|---|---|---|---|
| Internal reference (signal ÷ p_cal) | d ln p_cal/dT should be ≈ −κ | **+0.21 %/K**. The un-divided signal keeps **+1.53 %/K**. Per-unit corr(κ, κ_pcal) = −0.03. | excluded |
| Detector sensitivity (photon counting, APD) | Regulated (24.85 °C, flat). Background enters through the same receiver, so a receiver gain change should move the pair **background** ratio with ΔT. | Budapest: signal +2.22 ± 0.23, background +0.44 ± 0.14 %/K. Măgurele: +1.55 ± 0.33 vs −1.05 ± 0.27 %/K. Schiphol too noisy (±0.8). | not the main driver |
| Additive electronic baseline / dark | Relative effect should shrink under strong signal | Clean nights **+1.40**, aerosol-loaded **+1.22 %/K**. Night background moves only +0.37 %/K. | excluded network-wide (Rome 0-380-6-1 is an exception: loaded +4.1 vs clean +1.4) |
| Optics / overlap / alignment (Hervo 2016 type) | Effect should concentrate below ~600 m | Flat **+1.2–1.5 %/K from 0.5 to 8 km**, small bump at ~0.9 km | excluded for the constant |
| Window | `window_transmission` vs T | no intra-night change | excluded |
| Laser pulse energy (driver / pump electronics at housing T) | Multiplicative, height-independent, absent from the background. For a passively Q-switched laser the pulse rate is the energy's mirror. | Pulse rate moves only **−0.06 %/K** (p10–p90 −0.15 … +0.03), uncorrelated with κ (r = −0.08). Literature: E ∝ 1/σ gives ≈ +0.2 %/K of *crystal* temperature, and the crystal is regulated (§7.2). | **unlikely** as the main driver |
| **Laser line vs narrowband filter** | Same as above. The background integrates over the whole band, the laser line sits at one point of it. Predicts curvature over wide temperature ranges. | Quadratic term mostly negative (§7.2). But the manual puts the 1 nm filter inside the heated laser-optical module, so it needs a heat leak (§7.8). | **leading, unproven** |
| Far-field alignment (boresight) drift | Multiplicative, concave, invisible in the background; FOV 0.45 mrad ≈ divergence | Should distort the incomplete-overlap zone, where κ is flat | **plausible rival** (§7.8) |

Other evidence:
* The optical module stays at 33.85 °C all night while the signal follows `temp_int` at 1–2 %/K.
  The responsible part therefore sits outside the regulated module, or is not represented by its
  sensor.
* The trend with module generation (≈ +1.0 %/K up to TUB15, +2.0–2.7 %/K for TUB20–24) points to a
  design change.

**Questions for Lufft:**
1. Where is the interference filter, and is it temperature-controlled? What is its bandwidth?
2. Can the internal laser-energy photodiode be read raw (the 100 % status is clipped)?
3. Are the laser driver / pump diodes on the main board at housing temperature?
4. What changed in the TUB16+ and TUB20+ modules?

## 6. Consequences and recommendations

1. **The "one constant per night, valid all day" assumption is wrong for the CHM15k.** Daytime L2 β_att
   is biased high by +5 % (median), +7.5 % in summer, and +15–30 % for the most sensitive units
   (Aosta, Schiphol B).
2. **Part of the seasonal cycle of the constant is instrument temperature.** Nights are ~11 K warmer
   (p5–p95 range of `temp_int`) in summer than in winter. With ~2 %/K, that is ~20 % of the 34 %
   seasonal amplitude (report 09 §2.3). The Kalman annual-cycle mode absorbs it only on average.
3. **Correction:** C_L(t) = C_L,ref · exp(κ_u (T_int(t) − T_ref)).
   * **κ_u:** estimated per unit by the intra-night method on the operational archive. It needs ~50
     clear nights, refreshed at module swaps. Use the generation median when there are too few.
     Use the intra-night κ (it is conservative). The multi-day κ is ~1.5× larger but more exposed to
     residual atmospheric confounds.
   * **Nightly constants:** normalising them to T_ref before the Kalman filter reduces the whole-record
     scatter from 18.8 to 15.2 % (130/148 units improve). Night-to-night scatter goes from 8.9 to 8.4 %.
   * **L2:** apply per profile in the L1→L2 step. This removes the daytime bias the night-only
     calibration cannot see.
   * **Validation:** the Schiphol/Budapest/Măgurele pairs give a direct check. The day−night pair
     ratio should collapse to ~0 after correction.
4. **Payerne CHM15k** (TUB140016): check the receiver (night background 5–25× the network). Do not use
   its night-derived κ: it over-corrects (§7.1, §7.6). Use a coefficient from the daytime CL61 comparison
   or leave it uncorrected.
5. **Report 09 §2.3** should be read with §2 here: its correlation with internal temperature was a
   raw (season-mixed) correlation. The dependence is real but ~⅓ smaller than the raw slope.

## 7. Follow-up investigations (2026-10-20, afternoon)

### 7.1 Why Payerne does not fit

![Payerne: instrument triangle and housekeeping cycles](figs_chm15k_stability/q6_payerne.png)

* **The Payerne CL31 cannot arbitrate.** Its day−night signal is 50–70 % off relative to *both* the
  CHM15k and the CL61 (CL61/CL31 −61 %, CHM15k/CL31 −50 %), its own daytime artefact.
* **Direct daytime coefficient.** Day-to-day variations of the CHM15k/CL61 ratio at a fixed solar hour
  (mean diurnal cycle removed) against the CHM15k `temp_int` anomaly:

  | Site | Slope |
  |---|---|
  | **Payerne** | **+0.7 ± 0.4 %/K** |
  | Lindenberg | +3.0 ± 0.2 %/K |
  | Camborne | +3.5 ± 0.5 %/K |
  | Aosta | +4.9 ± 0.2 %/K |

  Payerne's net daytime response really is small; the +20 % prediction was wrong, not the data.
* **Why Payerne's night coefficient is wrong.** The Payerne CHM15k receiver is in an exceptional state:
  * night background `bckgrd_rcs_0` 0.021, 5–25× the other units;
  * noise `stddev` 2–5× higher;
  * the lowest `calibration_pulse` (0.036);
  * hourly noise in the Rayleigh window of 21.5 %, against a network median of 4.3 %.

  Its estimators disagree accordingly: intra-night +2.0 ± 0.8 %/K, multi-day −2.3 ± 1.4 %/K. This
  matches the strong measured dark of this unit (report 16).
* **Regulation loss with the opposite sign.** Payerne's optical module leaves regulation on **37 % of
  daytime hours**, more than any other site tested (anomaly sd 1.6 K). With `temp_lom` in the daytime
  fit (§7.3), the two effects separate: `temp_int` +1.8 ± 0.3 %/K (in line with the network) and
  `temp_lom` −2.0 ± 0.5 %/K. The two partly cancel.
* **Module history.** TUB140016 was at Payerne 2015-11 → 2021-02 (archive M2 +4.45 ± 0.45 %/K).
  TUB200009 followed, 2021-03 → 2025-04 (+3.67 ± 0.48 %/K), and TUB140016 has been back since 2025-04,
  presumably after servicing; its current coefficient is poorly determined.
* **CL61 side.** Le et al. (2026) report temperature-dependent CL61 background/bias but no temperature
  dependence of the calibration. The Payerne CL61's laser energy rises only +2.5 % by day (Lindenberg
  +1.9 %). There is no sign that the CL61 compensates.
* **Consequence:** at Payerne the correction with the night κ over-corrects (§7.6). Payerne needs:
  * a receiver check (background/dark);
  * a coefficient from the daytime CL61 comparison rather than from its noisy nights;
  * or exclusion until then.

### 7.2 Laser energy or filter? Literature and two new tests

* **Laser pulse energy (literature).**
  * For a passively Q-switched Nd:YAG/Cr:YAG laser the pulse energy scales as E ∝ 1/σ (Degnan, 1995).
  * The 1064 nm emission cross-section falls by 0.20 %/K (Sato and Taira, 2012). Pulse energy therefore
    **rises** with *crystal* temperature: ≈ +0.2 %/K from σ; +0.4…+0.7 %/K measured over 100 K
    (Kimmelma et al., 2008; Bass et al., 2003).
  * In the CHM15k the laser sits in the regulated optical module (33.85 °C), so this mechanism is frozen.
* **Laser pulse energy (data).** In a passively Q-switched laser the pulse energy and the pulse rate
  trade off at fixed pump power. The CHM15k pulse rate moves only **−0.06 %/K** with `temp_int`
  (p10–p90 −0.15 … +0.03), uncorrelated with κ (r = −0.08). That bounds a pulse-energy change to
  ≈ +0.06 %/K, twenty times too small.
* **Filter mismatch (literature).**
  * The Lufft datasheet gives a **1 nm** filter ("1 nm or less" in the manual). A 0.1 nm figure on the
    Uni Köln instrument page is not supported by the manufacturer (§7.8).
  * Dielectric filters drift +0.002…+0.02 nm/K.
  * With a 1 nm passband and the laser line on its flank, the transmission at the line changes by
    about 1–4 %/K, unit-specific in size and sign.
  * This also predicts **curvature**: the slope must fall as the passband reaches the line.
* **Curvature (data).**

  ![Curvature of the temperature response](figs_chm15k_stability/q10_curvature.png)

  The quadratic term of ln C vs `temp_int` (season and trend controlled) is negative:
  * network −0.074 ± 0.006 %/K² (9 units significantly negative, 1 positive);
  * archive −0.164 ± 0.012 %/K² (9 negative, 0 positive).

  The response flattens at the warm end, as a passband approaching the line predicts.
* **Verdict (revised in §7.8 after the adversarial review).** The laser-line vs filter mismatch is the
  leading hypothesis, with a thermal alignment drift as a co-equal rival; pulse energy is unlikely. Madonna et al. (2018) reported the CHM15k pulse count varying ~10 %
  seasonally with enclosure temperature, but `beta_raw` is normalised by the pulse count. The PROBE COST
  report of Van Hove and Diémoz (2024) links the network-wide seasonal cycle to pulse numbers, detector
  sensitivity and background. Buxmann (2024) proposes aerosol in the window; §7.4 bounds that to ~15 %.

### 7.3 `temp_lom` (laser optical module temperature)

![temp_lom: night dither and daytime regulation loss](figs_chm15k_stability/q7_tlom.png)

* **Night.** `temp_lom` only dithers by 1 K steps (115 units). Its partial slope is
  +0.30 ± 0.05 %/K (significant at 4 units), against `temp_int` +1.57 %/K in the same fit.
* **Day, regulation lost** (hot afternoons, 3–37 % of daytime hours). Extra slope on `temp_lom` with
  `temp_int` held fixed:

  | Site | Extra slope on `temp_lom` |
  |---|---|
  | Aosta | +4.2 ± 0.5 %/K |
  | Budapest A−B | +1.4 ± 0.2 %/K |
  | Măgurele A−B | +0.4 ± 0.2 %/K |
  | Payerne | −2.0 ± 0.5 %/K |
* **Answer.** `temp_int` remains the better single predictor. `temp_lom` carries a real but
  unit-specific second term, of either sign, only when the module leaves regulation. It is worth adding
  as a second regressor at hot sites, not as a replacement.

### 7.4 AERONET: is it aerosol?

![AERONET control of the temperature coefficient](figs_chm15k_stability/q5_aeronet.png)

The data: AERONET Level 1.5 daily AOD (v3) at sites within 15 km of a CHM15k, 2025–2026. 27 units have
enough matched nights. Night AOD is the mean of the evening and morning days. AOD at 1064 nm comes from
AOD_1020 and the 440–870 Ångström exponent.
* Warm spells are slightly hazier: r(T anomaly, AOD anomaly) = +0.24 (median). The constant does respond
  to aerosol (+5.2 % per +0.1 AOD at 1064 nm, median).
* The temperature coefficient (31-day anomalies, same nights):

  | Model | Pooled slope |
  |---|---|
  | No aerosol term | +2.36 ± 0.11 %/K |
  | + AOD | +2.07 ± 0.11 %/K |
  | + AOD + Ångström | **+2.01 ± 0.11 %/K** |
  | Cleanest third of nights only | **+2.13 ± 0.15 %/K** |
* About 15 % of the multi-day slope is aerosol; ~2 %/K remains. Aerosol does not enter the intra-night
  and same-wavelength-pair estimates at all. **Aerosol is ruled out as the main cause.**

### 7.5 Housekeeping ML (`estimate_calib_from_housekeeping_v4.m`) vs the exponential model

![ML vs exponential model on the 2024 test year](figs_chm15k_stability/q8_ml_vs_exp.png)

* **Data:** the 11-year archive (13 stations), cleaned as in v4 but with the outlier screen applied
  **per optical module**. A station-wide screen flags whole module eras across a swap: Payerne 2024
  went from 21 to 43 valid nights.
* **Split:** train 2013–2023, test 2024 (the split v4 defines). The ML is the v4 bagged tree (30 trees,
  min leaf 8; predictors `p_calc`, `temp_int`, `base`, `laser_pulses`, `state_optics`, `cs`), retrained
  on the training years and averaged over 10 seeds.
* **The exponential is per module (`serlom`).** Nothing crosses a TUB swap. 104 test nights on modules
  swapped in during 2024 cannot be predicted by it and are excluded from the common comparison.

| Scenario | Model | Robust spread of ln(C_obs/C_pred), 1040 common 2024 nights |
|---|---|---|
| A — housekeeping only | ML bagged tree (retrained) | 22.3 % |
| A — housekeeping only | **C_ref · exp(κ T_int)** | **20.2 %** |
| A — housekeeping only | ML saved model 202502 (saw 2024: in-sample) | 22.9 % |
| B — operational (last 30 d of calibrations known) | running median 30 d | 13.6 % |
| B — operational | **running median of T-normalised C × exp(κ ΔT)** | **13.0 %** (best at 9 of 11 stations) |
| B — operational | ML × its 30-day bias | 14.2 % |

* **The simple exponential wins in both scenarios.** The ML adds variance: in scenario B it is worse than
  the plain running median.
* Housekeeping-only prediction (A) fails mainly on level changes it cannot see (ageing, servicing).
  Payerne's constant ramps from 4.5 to 6.5 ×10¹¹ during 2024. That is why a recent calibration is
  indispensable and the temperature term is a correction on top of it.
* The nightly Rayleigh noise (~9–10 %) is the floor for every model.

### 7.6 Agreement with CL61 and EARLINET after correction

![Paper-v0.9 validation with and without the temperature correction](figs_chm15k_stability/q9_validation_tcorr.png)

**Setup**
* The paper-v0.9 validation (`run_paper_v09`) was run unchanged, with its caches written to this
  study's folder; the paper outputs are untouched.
* Three CL61 sites were added: Lindenberg, Camborne and Aosta, with the CL61 cloud constant as
  reference.
* The correction exp(−κ_u(T_int(t) − T_ref(t))) is applied per hour to the CHM15k:
  * κ_u is the intra-night coefficient;
  * T_ref is a ±15-day running median of the calibration-night `temp_int`.
* Bands: stations 500–3000 m, EARLINET 500–5000 m. Night = 22–02 UTC, day = 10–14 UTC.

**Results**

| Comparison | κ [%/K] | Day − night bias, operational → corrected [pts] | IQR, operational → corrected [pts] | r_log |
|---|---|---|---|---|
| Aosta, CHM15k vs CL61 | +3.25 | +22.5 → **+4.5** | 27.0 → **18.4** | 0.941 → 0.958 |
| Lindenberg, CHM15k vs CL61 | +1.58 | +10.1 → **+2.2** | 27.5 → **23.8** | 0.865 → 0.869 |
| Camborne, CHM15k vs CL61 | +1.07 | +1.9 → −4.0 | 25.7 → 26.2 | 0.973 → 0.973 |
| Leipzig, CHM15k vs EARLINET PollyXT | +0.99 | +7.1 → **+1.9** | 28.6 → 28.2 | 0.920 → 0.918 |
| Măgurele B, CHM15k vs EARLINET RALI | +1.47 | +12.6 → **+1.1** | 41.0 → **37.0** | 0.857 → 0.857 |
| Măgurele A, CHM15k vs EARLINET RALI | +1.60 | +13.2 → **+2.1** | 35.0 → **31.5** | 0.899 → 0.899 |
| **Payerne, CHM15k vs CL61** | +2.02 | +6.1 → **−18.9** | 33.9 → 36.8 | 0.929 → 0.922 |

* **6 of 7 comparisons improve.** The day−night bias collapses from 7–23 pts to ≤ 4.5 pts, and the
  spread drops by up to 9 pts (Aosta).
* **Payerne is over-corrected**, as §7.1 predicts: night κ unreliable, opposite-sign regulation loss.
* **Some overall median biases move** (e.g. Aosta −7.8 → −15.3 %). The operational near-zero overall
  biases partly came from night-negative and day-positive hours cancelling. After correction the bias
  is the plain night-calibration level difference to the reference, the same by day and by night.

### 7.7 Schiphol unit B: the "daytime artefact" is the temperature coefficient

Schiphol unit B (module TUB200002) is a known problem unit:
* **report 05 §3.6:** +17 % against the other units, +25 % by day vs +11 % at night, concentrated near
  range;
* **memory note:** its static applied overlap is ~25 % too low at 500–1000 m, back to 1 only by
  ~1.6–2 km;
* **paper v0.9, Sect. 6.3:** "an inaccurate factory overlap function combined with a daytime artefact
  whose origin remains to be investigated".

The day-minus-night ratio separates the two. A static overlap error is there by night and by day, so it
cancels. Day−night of the pair signal ratio (10–16 h solar, %):

| Pair | 0.5–1 km | 1–1.5 km | 1.5–2.5 km | 2.5–4 km | 4–6 km |
|---|---|---|---|---|---|
| **B / C** | +21.0 | +23.4 | +21.2 | +20.1 | +19.2 |
| A / B | −15.0 | −16.0 | −14.6 | −13.1 | −14.2 |
| A / C (healthy pair) | +5.1 | +6.5 | +6.7 | +6.5 | +6.3 |

* **B's daytime excess is flat from 0.5 to 6 km**, so it is a pure gain, not a near-range effect.
* It matches B's temperature coefficient, the highest at Schiphol (κ_B = +2.8 %/K, intra-night;
  C ≈ 0, A +0.8, D +0.5): +17 % predicted, +21 % observed for B/C (§4.4).
* The near-range weighting reported earlier came from comparing absolute profiles, which mix the static
  overlap deficit (night and day) with this daytime gain.
* **Consequences.**
  * Schiphol B has two independent defects: a static overlap deficit (near range) and the largest
    temperature coefficient of the site (all heights, by day).
  * The temperature correction removes the second (B/C day−night 21 → 4 %-pt predicted residual).
    The first needs a static overlap correction.
  * The paper's "daytime artefact whose origin remains to be investigated" (Sect. 6.3) can now be
    attributed to the temperature dependence of the constant.
* The D pairs fall at 4–6 km (A/D −9.7, C/D −13.5) because of daytime noise in unit D's upper range.
  That band is not used for the conclusions.

### 7.8 Filter physics, an adversarial review and a filter model

**What the manufacturer says (Lufft CHM15k manual R19 2022, datasheet 2024, FAQ 2022)**
* The laser-optical module (LOM) is the "transmitting **and receiving** unit".
* The receiver has "a 1 nm (or less) bandwidth filter in front of the detector". The datasheet says
  **1 nm**. The 0.1 nm quoted on the Uni Köln instrument page is not supported by the manufacturer.
* The laser is a diode-pumped Nd:YAG "microchip" laser (passive Q-switching is likely but not stated),
  divergence < 0.5 mrad, receiver field of view 0.45 mrad.
* The LOM is **heater-only** (250 W, valid range 25–49 °C, no cooling). That is why it leaves
  regulation on hot afternoons (§7.3).
* The avalanche photodiode has its own Peltier stage and a bias compensated at 2.4 V/K.
* The `p_calc` reference light pulse is injected into the receiver chain (visible at 16 km equivalent
  range). Its wavelength, and whether it passes through the filter, are not documented.

So the premise "filter at housing temperature" is **contradicted** by the documentation. The filter
sits in the heated module, possibly on the Peltier detector head. The filter hypothesis survives only
through a heat leak from the housing into the module.

**What shifts a dielectric filter's centre wavelength** (literature: Newport/MKS, Semrock, Andover,
Fraunhofer)
* **Temperature:** layer expansion and dn/dT give a red shift on warming. Typical rates:
  * hard ion-sputtered coatings: +2…+5 pm/K;
  * soft evaporated coatings: up to +20 pm/K;
  * laser-line filters: < 20 pm/K.
  * measured on commercial narrow-band filters (Blifford 1966, NCAR): +8…+29 pm/K at normal incidence
    from −60 to +60 °C, linear and reversible, with the passband shape unchanged — the behaviour the
    model of §7.8 assumes. Those were visible-range evaporated filters. The shift grows roughly with
    wavelength, so the fitted ~20 pm/K at 1064 nm is in the normal range for such coatings.
  * modern hard-coated filters are 5–10× more stable: < 3 pm/K for reactive ion-assisted sputtered
    filters (Tsai et al. 2001); 2–5 pm/K quoted for current lidar filters (Zenteno-Hernández et al. 2024).
  * The data constrain only α/W ≈ 0.02 K⁻¹, not α and W separately. Two filters fit equally well: a
    1 nm soft-coated filter at ~20 pm/K, or a ~0.2 nm hard-coated filter at ~4 pm/K. The manual's
    '1 nm or less' allows both, so the filter specification (Lufft question 2) is needed to decide.
  * Sensitivity of a Gaussian passband: d ln T/dT = 8 ln2 · (Δ/W) · (α/W). A catalogue 10 nm hard-coated
    1064 nm filter (2–5 pm/K) gives only 0.05–0.14 %/K, even with the laser at half maximum, so it is
    15–40× too stable. A filter that wide is incompatible with the observed κ. To give ~2 %/K, a hard
    coating needs W ≈ 0.1–0.25 nm.
  * **Working assumption: W = 1 nm** (Lufft datasheet; operator decision 2026-10-20). The fit then fixes
    the filter drift at **α ≈ 20 pm/K (IQR 15–24)**. That is exactly the coefficient Jenoptik publishes for its
    JENOFILT narrow-band filters: ΔCWL = (T − 22) × 0.02 nm. Jenoptik made the CHM15k before Lufft took it
    over. JENOFILT 602 is an 11 nm, 400–1000 nm filter of cemented glass slides, not necessarily the
    CHM15k part, but it shows that this filter technology drifts at the rate the data require. Hard-coated
    filters (2–5 pm/K) would give only 0.2–0.5 %/K at 1 nm. Under this assumption the filter is
    conventional (soft) coated, and the rise of κ with module generation means a larger detuning: the
    laser line sits Δ0 ≈ κ × 0.09 nm/(%/K) off the filter centre, i.e. ≈ 0.09 nm (≤ 2015 modules), 0.15 nm
    (2016–19) and 0.18 nm (≥ 2020), if α is the same in all generations.
  * **Where is the laser line?** The manual says the source is a diode-pumped Nd:YAG microchip laser, without
    an exact wavelength.
    * The Nd:YAG 1064 nm line is known to metrological precision: the iodine-stabilised, frequency-doubled
      Nd:YAG standard is 532.245036 nm in vacuum (CIPM 2007). The fundamental is therefore 1064.490 nm in
      vacuum, **≈ 1064.20 nm in air**. Microchip Nd:YAG lasers are reported at 1064.17 nm.
    * The gain line is ~5 cm⁻¹ ≈ 0.57 nm wide and its centre shifts by +3…+5.5 pm/K. A microchip cavity mode
      tunes faster, ~13 pm/K (dn/dT plus expansion), with mode hops. Inside the heated module the line
      should be fixed.
    * Consistency check, with no fit involved:
      * catalogue 1064 nm filters are specified at a CWL of **1064.0 nm** (+0.2/−0 nm) at room temperature
        (22–23 °C);
      * the laser at 1064.2 nm therefore sits **~0.2 nm on the red side** at 22 °C;
      * at 0.02 nm/K the filter reaches the laser at **≈ 32 °C**.
      * The fit finds Δ0/W = +0.19 (laser red of the filter centre) and a best-match temperature T* ≈ 31 °C
        (IQR 26–43). The unit-to-unit spread of Δ0 (0.11–0.32 nm) is of the order of a 1 nm filter's CWL
        tolerance.
    * This makes the filter picture self-consistent, but it is not a proof: the CWL tolerance is as large
      as the offset.
  * **No Jenoptik 1064 nm narrow-band filter is in the public catalogue.**
    * JENOFILT 605 covers 1060–1800 nm, but only as 28 nm wide blocking filters, without a temperature
      coefficient.
    * The 0.02 nm/K figure comes from JENOFILT 602 (400–1000 nm, 11 nm), which is made of cemented glass
      slides, as is JENOFILT 601.
    * A 1 nm 1064 nm filter in the CHM15k would have been a custom part. Lufft question 2 should ask for its
      supplier and datasheet.
  * Commercial 1064 nm filters (supplier specifications, 2026-10). κ is predicted at the fitted
    detuning Δ/W = 0.19, i.e. κ ≈ 1.05 · α/W:

    | Filter | FWHM W | Coating | α (dλ/dT) | α/W [K⁻¹] | Predicted κ |
    |---|---|---|---|---|---|
    | **Jenoptik JENOFILT 602** (Jena; datasheet 1996) | 11 nm | conventional, cemented glass slides | **20 pm/K** (spec: ΔCWL = (T − 22) × 0.02 nm) | — | — |
    | Andover standard bandpass (e.g. 064FS02-25) | 1 nm | conventional (type not stated) | ≈ 30 pm/K (read off the generic coefficient-vs-wavelength chart, Andover selection guide p. 4; not a part spec) | 0.030 | ≈ 3 %/K |
    | Alluxa ultra-narrow (1064-0.9 … 1064-0.4) | 0.4–1 nm | hard, sputtered | 2–5 pm/K (Alluxa) | 0.002–0.013 | 0.2–1.3 %/K |
    | Semrock MaxLine laser-line 1064 | a few nm | hard | < 5 ppm/K ≈ 5 pm/K | < 0.002 | < 0.2 %/K |
    | Newport 05LF10-1064 / Alluxa 1064-10 | 10 nm | — / hard | < 20 / 2–5 pm/K | < 0.002 | < 0.2 %/K |
    | **Observed (CHM15k fit)** | — | — | — | **0.020 (IQR 0.015–0.024)** | **+1.5…+2 %/K** |

    A conventional 1 nm filter of the Andover type matches the observed α/W within a factor of 1.5.
    A hard-coated filter matches only if it is ≤ 0.4 nm wide. The rise of κ in newer modules then
    points to a narrower or more detuned filter, not to a harder coating, which would lower κ.
* **Angle of incidence:** λ_c(θ) = λ_c0·√(1 − sin²θ/n_eff²), a blue shift. It is second-order at
  normal incidence: about −29 pm/mrad for a filter angle-tuned to 5°, a few pm/K for a plausible
  thermal tilt. Minor.
* **Humidity in porous (evaporated) coatings:** a red shift of up to several nm between dry and humid
  air, with a time constant of hours to days. At fixed absolute humidity, housing RH falls ~6 %/K as
  `temp_int` rises, so a soft coating could mimic a large temperature coefficient. It does not apply
  to sealed or sputtered filters.
* **Mechanical stress and ageing:** soft coatings drift red over years. An ageing component could
  masquerade as the module-generation trend.

**Magnitude check.** With a Lorentzian passband the maximum relative slope is about 2/FWHM, i.e.
200 %/nm for 1 nm. A 1 nm filter drifting 5–20 pm/K with the laser on its flank therefore gives
1–4 %/K; a hard-coated filter at normal temperature coupling gives ≤ 1 %/K.

**Adversarial review — the competing explanations**

| Hypothesis | For | Against | Verdict |
|---|---|---|---|
| Filter vs laser line | multiplicative, height-flat, invisible in the solar background, unit-specific, concave | filter inside the heated LOM (needs a leak); 1 nm datasheet; one-signed κ needs a systematic offset; Hopkin et al. (2019) saw no annual cycle at Aberporth | plausible, unproven |
| **Far-field alignment (boresight) drift** | FOV 0.45 mrad ≈ divergence < 0.5 mrad; multiplicative and concave; background-blind; Hervo et al. (2016) prove thermo-sensitive optomechanics | should distort the 0.5–1.5 km incomplete-overlap zone | plausible, co-equal rival |
| Laser pulse energy | no direct energy monitor; Madonna et al. (2018) saw a seasonal pulse count | pulse rate flat (−0.06 %/K) unless servoed; crystal in the heated LOM | unlikely |
| Detector / electronics / window | — | Peltier detector; background and `p_calc` do not follow | refuted / unlikely |

Our data offer two discriminators against the alignment rival:
* κ is the same from 0.5 to 8 km, including the incomplete-overlap zone (§5);
* the Schiphol B/C day−night ratio is flat from 0.5 to 6 km (§7.7).

A boresight drift should reshape the near range. This favours a spectral (filter) mechanism, but it
is not conclusive for a tilt perpendicular to the transmitter–receiver baseline.

**A physical filter model** (`q11_filter_model.py`)

![Filter model: peak temperature, detuning, thermal lag](figs_chm15k_stability/q11_filter_model.png)

* Assumptions:
  * fixed laser line;
  * Gaussian passband of FWHM W;
  * centre drifting by α per K of filter temperature;
  * filter temperature = `temp_int` low-pass filtered with lag τ.
* ln C is then exactly quadratic in temperature. Each unit's slope b and curvature c (season-
  controlled, §7.2) give:
  * α/W = √(−c / 4 ln2);
  * detuning Δ0/W = b / (8 ln2 α/W);
  * the temperature T* at which the filter would be centred on the line.

**Results**
* **Curvature:** negative in 110 of 141 unit-eras, but significant in only 9. The inversion is
  therefore weakly constrained.
* **α/W ≈ 0.020 K⁻¹** (all c < 0), 0.027 K⁻¹ for the significant set. With the datasheet W = 1 nm this
  needs α ≈ 20–27 pm/K at the filter: a soft-coated filter with near-full temperature coupling.
  A hard-coated filter at 5 pm/K would need W ≈ 0.2 nm.
* **Detuning:** almost all units have the laser line on the red side of the passband centre
  (Δ0/W ≈ +0.1…+0.3), with a median transmission loss of ~10 % at the mean operating temperature.
  The near-uniform sign requires a systematic design offset, e.g. a filter specified for a warmer
  temperature.
* **Peak temperature:** T* ≈ 26 °C (significant set) to 31 °C (all), against a mean operating
  temperature of 20.5 °C and the LOM set-point of 33.85 °C. Suggestive of a filter specified near the
  module set-point but operating colder, but noisy.
* **Thermal lag:** the intra-night fit is best with **no lag** (τ = 0–0.5 h) and degrades steadily for
  τ ≥ 1 h. The responsive element follows the housing sensor within ~30 min. That is hard to
  reconcile with a component deep in the heated, thermally massive LOM, and it is a point for a
  housing-coupled element (a mount, a tube, or the alignment).
* Consequently the multi-day vs intra-night difference in κ (§3) is **not** explained by a simple
  lag. It is more likely residual atmospheric or seasonal content in the multi-day estimate.

**Verdict.**
* The cause is a multiplicative optical-throughput change in the transmitter–receiver path, tied to
  the housing temperature.
* Filter-vs-line and boresight drift both fit; the height-flatness slightly favours the filter, the
  no-lag result slightly favours a housing-coupled element.
* Pulse energy, detector, electronics and aerosol are excluded or unlikely.
* Only Lufft information (or a climate-chamber test) can decide. The correction itself (§8) does not
  depend on which mechanism it is.

**Open questions for Lufft**
1. Where exactly is the receiver interference filter mounted: on the Peltier detector head, in the
   heated LOM body, or in a part that follows the housing temperature?
2. Filter specification: CWL and tolerance, FWHM (1 nm or less?), passband shape (cavities), coating
   type (sputtered or evaporated; sealed?), dλ/dT, and the temperature at which the CWL is specified.
3. The factory-measured laser emission wavelength at the LOM set-point, per optical module.
4. Is the laser crystal temperature controlled beyond the LOM heater? Is the PRF closed-loop
   regulated through the pump current?
5. Does the `p_calc` reference light pass through the filter? Is it at 1064 nm (an LED or the laser)?
6. Can the transmitted pulse energy be read raw? `status_laser` is clipped at 100 %, and what does
   `LaserPower` (LAP) measure?
7. What changed in the optical module (filter, coating, supplier, optics, alignment) around 2015–2016
   and 2020?
8. Is the transmitter–receiver alignment athermal? What tilt (µrad/K) is expected between the beam and
   the 0.45 mrad field of view?
9. Payerne TUB140016: service history since its return in 2025 (abnormally high night background).
10. Can a climate-chamber test be run: housing temperature varied at constant LOM temperature, with an
    external energy meter at the window and a tunable 1064 nm source mapping the receiver passband?
12. Is there any compensation in the firmware before `beta_raw` (beyond p_calc, cs, overlap) that could depend on
    the season or the sun? Is the transmitter/receiver alignment sensitive to direct solar heating of the head
    (housing orientation)? The CHM15k shows a ~±17 % solstice-symmetric gain cycle every year (§7.14).
11. Why did the default detector set-point (`ApdTemp`) move from 25 to 30 °C around 2019–2020, and
    does Lufft recommend a value for a stable calibration?

### 7.9 The paper's Kalman figure with the temperature correction

![Kleine Scheidegg Kalman, operational vs temperature-corrected](figs_chm15k_stability/q12_kalman_tcorr_kse.png)

* The paper's operational Kalman filter (`kalman_improve_alc`) is run on the Kleine Scheidegg nightly
  constants as published (top) and after normalisation to 18 °C internal temperature (bottom):
  C_norm = C·exp(−κ_u (T_cal − T_ref)), with κ_u = +1.04 %/K (intra-night) and T_cal the
  calibration-window `temp_int`.
* The comparison uses the same 231 nights; the HK cache ends in July 2026.

| Version | Annual-cycle amplitude | Scatter around the model |
|---|---|---|
| Operational | ±17.6 % | 11.1 % |
| Corrected, κ = 1.04 %/K | ±14.0 % | 10.4 % |
| Corrected, multi-day κ = 1.45 %/K | ±12.7 % | 9.5 % |

* Temperature explains about a quarter of this unit's seasonal cycle, a weak case (§7.10). The rest
  has other causes (aerosol in the window, ageing, …).
* Script: `presentation/chm15k_temperature/make_kalman_tcorr.py`.

### 7.10 The seasonal cycle across the network — why Kleine Scheidegg was a weak example

![Seasonal amplitude before and after normalisation, all units](figs_chm15k_stability/q13_kalman_network.png)

Same Kalman annual-cycle model as §7.9, all CHM15k (main module, ≥ 100 nights; `q13_kalman_network.py`):

| | Operational | Within-night κ | Multi-day κ |
|---|---|---|---|
| Median seasonal amplitude (127 units) | ±16.3 % | ±11.6 % | ±10.0 % |
| Units with a smaller amplitude | — | 90 % | 87 % |
| Median scatter around the model | 12.9 % | 11.7 % | 11.5 % |

* About **35 %** (median) of the seasonal amplitude is removed. The removal grows with κ × the unit's
  seasonal swing of the calibration-night internal temperature (median 4.4 K).
* The multi-day κ comes from 31-day anomalies with the season removed, so applying it to the seasonal
  amplitude is not circular.
* Best cases:
  * Budapest B: ±19.8 → ±7.0 % (within-night κ +2.33 %/K) → ±2.8 % (multi-day κ);
  * Plovdiv: ±28.4 → ±12.2 → ±6.7 %;
  * Vásárosnamény: ±19.1 → ±5.8 → ±1.8 %.
* **Kleine Scheidegg** (§7.9) is a weak case: an old module (κ +1.04 %/K) and a small swing (high
  alpine site, heater-held winter temperature). Its remaining seasonal cycle has other causes.

![Budapest B, operational vs temperature-corrected](figs_chm15k_stability/q12_kalman_tcorr_bud.png)

**Can the pulse energy be estimated from the number of pulses?** In principle yes. For a passively
Q-switched laser, average power = E × PRF, so at constant pump power E ∝ 1/PRF. In practice:
* the pulse count per record moves only ±0.5 % (31-day anomalies) and −0.07 %/K with temperature, so
  the energy it would imply changes by about +0.07 %/K, about 1/25 of κ;
* the per-unit regression of ln C on ln(pulses) is unstable (median −1.5, IQR −4.6…+0.2) for lack of
  pulse-count variance;
* if the pulse rate is servo-controlled through the pump current (Lufft question 4), the proxy fails.

The firmware already divides the signal by the pulse count.

### 7.11 What differs between optical-module generations

![Optical-module generations](figs_chm15k_stability/q14_generation.png)

Production year from the TUB serial (127 streams, 2025–26; `q14_generation.py`). Medians for ≤ 2015 /
2016–2019 / ≥ 2020 modules:

| Feature | ≤ 2015 | 2016–19 | ≥ 2020 | Spearman vs year |
|---|---|---|---|---|
| κ intra-night [%/K] | +0.95 | +1.65 | +2.02 | **+0.42** (p 1e-6) |
| κ multi-day [%/K] | +1.67 | +2.38 | +2.42 | +0.28 |
| Night background (dark) | 0.0038 | 0.0020 | 0.0017 | **−0.53** (p 2e-10) |
| Night noise (`stddev`) | 2.1e-4 | 1.5e-4 | 1.4e-4 | **−0.53** |
| Detector set-point [°C] | 24.85 | 24.85 | 29.85 in most | **+0.52** |
| Optical-module set-point [°C] | 33.85 | 33.85 | 34.35–34.85 | **+0.50** |
| Curvature, noon `p_calc` drop, laser age*, `temp_int` | — | — | — | not significant |

* **Break 1, 2015.** The night background (detector dark counts) falls about fourfold, from yearly
  medians ~0.01 (2011–14) to ~0.003. This matches the hardware change of June 2015 in the Lufft FAQ
  (firmware `_552` for units built since then). κ rises from 2016 (+1.5 %/K).
* **Break 2, ~2019–2020.** The detector Peltier set-point moves from 24.85 to 29.85 °C in most new
  modules, and the optical-module heater set-point from 33.85 to 34.35–34.85 °C. Modules with the
  29.85 °C detector have κ +2.02 %/K, against +1.34 for those at 24.85 °C.
* \*Laser age here ignored the 16-bit wrap of `laser_life_time` at 65 536 h (§7.13); true ages are 28–51 kh.
* **Firmware is not the cause; κ follows the module hardware.** Firmware comes from the L1 attribute
  `instrument_firmware_version`.
  * Spearman κ vs firmware = 0.00, against +0.42 vs production year.
  * Eleven ≤ 2015 modules running the newest firmware 1.15 keep κ +0.95 %/K, the old-module value.
  * The low-κ firmware 0.743 group (8 modules ≥ 2020, κ +1.13) is the UK Met Office network. It is
    confounded with the maritime UK sites, not a firmware effect.
* **Which settings can the user change?** From the Lufft manual R19 (Table 5, service mode):
  * `ApdTemp` (APT): the detector Peltier set-point, default 30 °C. The manual says to change it
    "only after contacting the Lufft service".
  * `PeltierMode` (0–1024) and `ApdControlMode` are also service settings.
  * No parameter sets the optical-module temperature. The module is heated only, with a 25–49 °C valid
    range, so its set-point is factory-defined.
  * The 24.85 → 29.85 °C detector change is therefore a factory default, not an operator choice. Among
    2019–24 modules alone, the set-point changes κ only slightly: +2.02 %/K at 30 °C (n 21) against
    +1.87 at 25 °C (n 32). The module generation matters more than the detector temperature.
* Pulses per record and the constant level are not compared: they depend on each site's averaging
  time and on the factory scaling.
* Reading: the receiver (detector, filter or their thermal design) changed in 2015 and again around
  2020, and each change increased the temperature sensitivity. Lufft question 7.

### 7.12 Proposed test: step the detector set-point

**Question.** Does the detector temperature change the calibrated signal? If the filter sits on the
Peltier detector head (Lufft question 1), the same test also scans the filter passband.

**Why not a monotonic 1 °C-per-hour staircase.** Over a night, `temp_int` falls by 1–2 K. With
κ ≈ 2 %/K, that alone moves the signal by a few per cent, and so do aerosol and boundary-layer
changes. A one-way ramp would alias all of this into the detector response. Short blocks that return
to the baseline cancel these slow drifts.

**Set-up**
* **Instrument.** Prefer a healthy ≥ 2020 module next to a second CHM15k that stays unchanged (the
  Schiphol, Budapest or Măgurele pairs): the ratio cancels the atmosphere. At Payerne, use the CL61 as
  the atmospheric reference. Note that the Payerne CHM15k (TUB140016) is the network outlier (§7.1).
* **Setting.** `ApdTemp` (service mode). Agree the range with Lufft first, as the manual requires.
  Record the original value and restore it afterwards.
* **Sequence.** Clear, calm night, 22–04 UTC, 20-min blocks: 5 min to settle (check
  `temperature_detector` is within 0.1 K of the set-point), then 15 min of data. Alternate the baseline
  B = 30 °C with randomised levels across 26–34 °C, for example: B 26 B 34 B 28 B 32 B 27 B 33 B 29 B 31 B.
  That is 18 blocks and 8 levels. Repeat on a second night in reverse order.
* **Record** per block:
  * raw `rcs_0` and `beta_raw` at 2–4 km (above the boundary layer);
  * `calibration_pulse` (`p_calc`), `bckgrd_rcs_0` and `stddev`;
  * `temperature_detector`, `temperature_optical_module` and `temp_int`;
  * the reference instrument over the same band.
* **Analysis.**
  * Compare each level with the mean of its two neighbouring baseline blocks.
  * Regress the log ratio on the detector temperature. Fit a quadratic as well, in case the passband
    is being scanned.
  * Exclude the test nights from the operational calibration and flag them to the E-PROFILE hub.

**Reading the outcome**

| Raw signal | `beta_raw` (÷ `p_calc`) | Interpretation |
|---|---|---|
| changes | flat | Detector gain changes, but `p_calc` compensates it. The detector is excluded as the source of κ. |
| changes | changes, ~1–2 %/K, curved | Something on the detector head is temperature-sensitive and not seen by `p_calc`: most likely the filter. Map the passband from the curvature. |
| flat | flat | The detector head is insensitive. The filter is elsewhere, so look at a heat leak into the module or at alignment (§7.8). |

**Detection limit.** For a 15-min block ratio to a co-located CHM15k, assume ~1 % noise. Eight levels,
each bracketed by baselines, over two nights then resolve about ±0.1–0.2 %/K: well below the expected
1–2 %/K if the filter sits on the head.

### 7.13 Why the temperature correction leaves half of the annual cycle — attribution of all factors

Six parallel analyses (q15a–f) quantified every candidate factor on 126–136 CHM15k units (main
optical module, 2025–2026, Payerne excluded). Each factor is given as the effect on C_L, as network
medians. Reference values: annual amplitude ±16.3 % (peak DOY ~177); night-to-night scatter 12.5 %.

![Variance and amplitude budget](figs_chm15k_stability/q15f_budget.png)

| Rank | Factor | Sign | Sensitivity | Night-to-night | Day vs night (L2) | Annual amplitude, in phase | Evidence |
|---|---|---|---|---|---|---|---|
| 1 | **Instrument temperature, κ × T_int** (filter / line) | + | κ 1.3–2.1 %/K | ~4 % of variance | **+5 to +8 %** (summer +10 %) | **+5.7 to +7.4 %** (25–45 %) | 135 units, pairs, 4 estimators |
| 2 | **Unexplained common-mode gain, solstice-phased** | + | — | most of the shared ~58 % | ? | **~+8 to +9 %** (peak DOY ~150–177) | pairs + CL61 + 7–8 km band |
| 3 | Molecular model (US-1976 vs real T/p) | − in summer | 1:1 | 1.3 % | 0 | **−1.6 %** (anti-phase) | 126 units, CAMS T/p |
| 4 | Aerosol (AOD regression / inside the window) | + | +7 %/0.1 AOD | 2.2–2.5 % | 0 | +0.5 to +1.2 % | AERONET 22 units, CAMS |
| 5 | Window transmission | + | d ln C/d ln w ≈ 0.34 | ~1 % | −0.4 % | −0.7 to −0.9 % (peaks in winter) | 125 units, 1083 steps |
| 6 | Pulses-per-record normalisation / pulse energy | − / + | −0.07 %/K PRF | 0.1 % | ≤ +0.5 % | +0.3 to +0.4 % | 126 units |
| 7 | Window height, short nights, selection | ≈ 0 | +0.1 to +0.6 %/km | 0.5–1.2 % | 0 | ≤ +0.3 % | 134 units |
| 8 | Noise, background, afterpulse, overlap | ≈ 0 | — | < 0.5 % | residual solar baseline **+6 %-pt at 4 km** (daytime L2 aloft, not C_L) | < 0.3 % | 134 units, 27 in native L1 to 15 km |
| 9 | Aerosol extinction below the window | − | Klett corrects 0.93 of it | ~1 % | 0 | −0.3 % | 21 units |
| 10 | `p_calc` reference: seasonal / night-to-night | ≈ 0 | d ln C/d ln p ≈ +0.1 (division works) | ≤ 0.7 % | daytime −18 % drop is compensated | ~0.2 % | 125 units, 3 pairs |
| — | `p_calc` slow drift (ageing reference) | imported, − | cross-unit slope ≈ −0.7 | — | — | trend only: 5 units drift 12–16 %/yr | 5 units, tentative |
| 11 | H₂O absorption at 1064 nm (HITRAN + CAMS) | − | −0.012 %/mm PWV | 0.04 % | ≈ 0 | 0.1 % (anti-phase) | 29 units |
| 12 | Laser ageing (trend), detector regulation, λ⁻⁴ | — | −0.4 %/1000 h | 0 | ≤ 0.05 % | ≈ 0 (no aliasing) | archive 11 yr |
| — | Retrieval noise (u/C ≈ 5 %) | random | — | ~24 % of variance | — | 0 | — |

**What the analyses show**
* **The residual is a gain, not a distortion.** Clear-sky molecular-normalised signal at 7–8 km, where there
  is no aerosol, carries 14.8 of the 17.1 % annual cycle, in phase. κ is flat with altitude from 2.5 to 8 km.
* **It is common to co-located CHM15k but not instrument-specific.** After the κ normalisation, differences
  between twin units vanish: the Schiphol B/C ratio amplitude falls from 12.9 to 1.1 %. A common 8–10 %
  remains at Schiphol, Măgurele and Budapest. The co-located CL61 does not share it: the CHM15k/CL61 ratio
  keeps an 8–15 % amplitude.
* **It is not atmospheric.** Aerosol in the window, aerosol below it, the molecular model, water vapour and
  window/selection effects sum to about 0 % in phase. The lidar's own profiles show no summer aerosol excess.
  Only a vertically uniform free-tropospheric aerosol cycle (ΔAOD ~0.002) cannot be excluded with these
  data.
* **It is not temperature-phased.** C_L peaks at DOY ~177, but the calibration-night T_int peaks at DOY ~206.
  Any driver tracking a late-July temperature can remove at most ~half of the cycle, whatever its
  coefficient. Lags (best τ 12–24 h: +6 points), temp_ext, temp_lom, the detector temperature and a
  Gaussian passband do not do better.
* **What remains.** The remaining ~±9 % is a CHM15k-type instrumental gain that peaks near the **summer
  solstice**. It is shared by co-located units and identical at all heights. That points to a driver
  tracking **insolation or day length**, for example solar heating of parts the T_int sensor does not see,
  thermal gradients, or the short twilight-bounded summer nights. This has not been tested yet.
* Open: synoptic T_int anomalies give κ ≈ 0 at 7–8 km but +2 %/K in the window, while within nights the
  7–8 km band gives +1.44 %/K. Part of the multi-day κ may be atmospheric. The intra-night κ is the clean
  estimate.

![Annual cycle by altitude band and the flat κ(z)](figs_chm15k_stability/q15e_profiles.png)

![Driver temperatures, lags and phase](figs_chm15k_stability/q15c_factor1_drivers.png)

Side findings:
* `laser_life_time` is a 16-bit counter that wraps at 65 536 h. 18 units wrapped in 2025–26, and §7.11's
  "laser age" ignored the wrap.
* Every real laser replacement coincided with an optical-module swap. A swap that is not split adds ~3 %
  of fake annual amplitude.
* Slow `p_calc` drifts above ~5 %/yr are worth a QC flag.

Code: `validation/chm15k_stability/q15{a..f}_*.py`; results and per-family summaries:
`C:/DATA/Projects/202610_CHM15k_stability/q15/q15{a..f}_summary.md`.

### 7.14 The annual cycle is a real CHM15k gain: follow-up tests (q16a–g)

Seven further analyses tested the remaining hypotheses with the real calibration code. The operator asked about
window altitude, twilight, sunlight, the molecular detection limit, the lidar ratio, aerosol correction, code
errors, processing, humidity and a literature review.

**Answer: the annual cycle is a real, CHM15k-specific instrumental gain.** It is about ±17 %, symmetric about the
summer solstice, the same at every height and latitude, and present in every year since 2014. The instrument
temperature (κ) explains only its within-day part. No atmospheric, retrieval, code or processing factor explains
it.

![Same cycle in every year of the independent 11-year archive](figs_chm15k_stability/q16g_archive.png)

| Test (agent) | Result | Effect on the annual cycle |
|---|---|---|
| Real recalibration with fixed windows 3–4, 4–5, 5–6, 6–8 km (a; 20 units, reproduces published constants bit-identically) | amplitude 18–22 %, peak DOY 171–181 in every window; 6–8 km keeps 85–90 % | window altitude/boundary layer: ≤ 3 % |
| Twilight: dusk + dawn vs deep night, threshold 108° vs 100° (a) | sunlight in background ≤ 0.2 % of dark; dusk−dawn difference = κ×ΔT | < 0.5 %; do not adopt 108° (removes all summer nights north of ~51.5°N) |
| Sunlight drivers: day length, insolation 1–21 d, heating, LOM excursions (b) | timing solar, cause not: latitude test, synoptic sign or lag fail | ≤ 20 % of the cycle at best |
| Mirror test: spring vs mirrored autumn (b, g) | ln C −0.3…−1 % while T_int differs 4.5–5.2 K | excludes every temperature/humidity-phased driver |
| Molecular detectability (c) | SNR < 3 at 6.8 km (30 m gate); window 2.5–5.6 km; 7–8 km band is molecular (X = 1.00–1.03 to 10 km) | noise bias ≤ 0.2 %, annual 0.03 % |
| Lidar ratio / Klett correction (d) | whole correction 2τ ≈ 4.5 %; dC/dLR 0.8 %/10 sr; effective LR ~46–60 sr | ≤ 1.3 % in phase; a realistic seasonal LR error has the opposite sign |
| Code review (e) | published constants reproduced within 2 % (768 nights); raw molecular-fit slope already 21 % / DOY 176 | code adds +1…+2 % in phase |
| Processing: raw2l1, firmware, independent 11-yr chain (g) | raw2l1 copies `beta_raw` bit-identically; MATLAB + NCEP chain 2013–24: 17.1 %, peak DOY 171, 58/58 station-years in phase; persists on unchanged firmware for 8–9 yr | not processing, not firmware version |
| Absolute humidity (g) | peaks DOY 209, fails mirror, synoptic and climate tests | rejected (a window-bottom shape effect only) |
| CL61 control (g; operational constants, cloud method) | CL61 cycle 1–6 %; CHM15k/CL61 ratio 17–33 %; Aosta 7–8 km band 26.7 % vs 4.5 % on the same nights | the cycle is CHM15k-specific |
| Literature review, ~50 mechanisms (f) | no published mechanism fits all signatures | — |

**Code issues found** (q16e):
* **major:** the v2.2 fallback tier admits nights +9 % above the strict-tier constants. These are 33 % of June
  nights against 3–8 % in winter. Effect: +4 % on the June level, +0.2…0.8 pts on the amplitude.
* **minor:**
  * beam tilt is ignored: +1 % static on the 12 units at 15°;
  * the molecular depolarisation is 0.0301 instead of 0.0273 at 1064 nm (−0.35 %);
  * masked `rcs_0` cells are lost in a concatenate (no current effect);
  * nights with no usable profiles after 00 UTC are stamped with the previous day's date.

**Consequences**
* **Apply κ only to sub-daily T_int deviations.** Model the annual gain cycle separately, as the Kalman annual
  term does. Applying κ×T_int seasonally over-corrects: it creates a spring/autumn asymmetry and a June dip.
* **Remaining candidates, both instrumental and untested:**
  1. **Sun-driven alignment of the biaxial head.** The noon solar elevation swings ±23.4° at every latitude.
     The single-lens CL61 would be immune. Possible tests: housing orientation, the 15° vs 3–5° tilted units, a
     shading experiment.
  2. **An undocumented firmware compensation** applied before `beta_raw`. Ask Lufft.

  Both go to Lufft (new question 12).

Results and per-test summaries: `C:/DATA/Projects/202610_CHM15k_stability/q16/q16{a..g}_summary.md`; literature
review `q16f_literature_review.md`; code `validation/chm15k_stability/q16*_*.py`.

![Fixed-window recalibration: the cycle is the same in every window](figs_chm15k_stability/q16a_windows.png)

![CL61 control](figs_chm15k_stability/q16g_cl61.png)

## 8. Recommendation: a temperature-dependent C_L for the CHM15k — SWOT

**Recommendation: yes, introduce C_L(t) = C_L,ref · exp(κ_u (T_int(t) − T_ref)), but as a
calibration-side change first and an L1→L2 change second, behind guards.** It is not a change to
apply network-wide in one step.

### 8.1 Proposed design

1. **Calibration side** (this repo, low risk, reversible).
   * Write the calibration-window `temp_int` of every nightly constant into the calibration NetCDF
     (today only `temperature_optical_module` is stored).
   * Estimate κ_u per optical module with the intra-night method, with elapsed time held fixed (§3).
     It needs ≥ 50 clear nights and se ≤ 0.4 %/K. It is refreshed at every module swap.
   * Normalise the nightly constants to a fixed network reference T_ref = 25 °C before the Kalman
     filter.
   * Publish C_ref, κ_u, T_ref and a κ quality flag next to `lidar_constant`.
   * Units without a valid κ get κ_u = 0 (today's behaviour). Units flagged as outliers, like Payerne
     TUB140016, are forced to 0.
2. **Shadow L2** (3–6 months). Produce corrected L2 in parallel at the evaluation sites:
   * Schiphol A–D, Budapest, Măgurele (the pairs);
   * Payerne, Lindenberg, Camborne, Aosta (CL61);
   * Leipzig, Măgurele (EARLINET).

   Track the day−night pair ratio and the day−night bias vs CL61/EARLINET. Go/no-go criterion: both
   ≤ 5 pts at ≥ 80 % of these comparisons, and no comparison degraded by more than 5 pts.
3. **L1→L2 hub.** Apply exp(−κ_u (T_int(t) − T_ref)) per profile. The hub already reads `temp_int` and
   keys a CHM15k temperature correction on `optical_module_id` (the overlap model
   `1 + (a(z) T + b(z))/100`). The gain term is the same machinery with a height-independent
   coefficient.
   * Smooth `temp_int` (1 K quantisation → 1.5–2.5 % steps) with a ~1 h running mean.
   * Cap the factor (e.g. ±30 %).
   * Record κ_u and T_ref in the L2 metadata.

### 8.2 SWOT

| | |
|---|---|
| **Strengths** | **Weaknesses** |
| • Robust evidence. Four independent estimators agree (+1.3…+2.1 %/K): synoptic anomalies on 158 unit-eras, intra-night on 150 units, the 11-year archive, and atmosphere-free same-wavelength pairs. Aerosol is controlled with AERONET (−15 %). | • The physical mechanism is not confirmed. The laser line vs filter mismatch is favoured, but Lufft has not confirmed it. |
| • Simple: one parameter per optical module, interpretable, cheap to estimate from data already collected. | • `temp_int` is a proxy (sensor ≠ component). The coefficient depends on the time scale (hours 1.3, day ~1.5, weeks 2.1 %/K), so a static κ under-corrects slow changes. |
| • Beats the housekeeping ML: 13.0 vs 14.2 % operationally, 20.2 vs 22.3 % HK-only. | • The response is curved at the warm end, so the linear-in-T model over-corrects hot afternoons slightly. |
| • Validated against independent instruments. The day−night bias vs CL61/EARLINET goes from 7–23 to ≤ 4.5 pts at 6/7 comparisons, and the spread drops by up to 9 pts. | • Per-unit κ is noisy (±0.2–0.4 %/K, 1 K quantisation). New modules stay uncorrected until ~50 clear nights accumulate. |
| • It removes a systematic daytime bias the night-only calibration cannot see: median +4.8 %, summer +7.5 %, up to +30 %. | • Failure mode demonstrated: an abnormal unit (Payerne) is over-corrected (−19 pts). Per-unit QC is needed. |
| • The hub hook already exists (CHM15k T correction keyed by module on `temp_int`). | • Regulation loss on hot afternoons (`temp_lom`) is a second, unit-specific term that is not captured. |
| | • Validation covers 7 comparisons over 2025–26 only. |
| **Opportunities** | **Threats** |
| • Normalising before the Kalman removes part of the seasonal cycle of the constant (~20 % of the 34 % amplitude) and cuts record scatter from 18.8 to 15.2 %. That means fewer false Kalman alarms and cleaner swap detection. | • A wrong κ injects a bias that did not exist. That is worse than doing nothing for that unit; see Payerne. |
| • Day/night-consistent L2 benefits every daytime user: boundary-layer height, aerosol layers, model assimilation (OmB). | • Governance: the hub belongs to the E-PROFILE Data Hub. The change needs agreement, versioning, and possibly reprocessing; L2 users (Cloudnet, NWP) get a time-varying `calibration_constant_0`. |
| • κ as a health indicator: a sudden κ change flags receiver problems (Payerne-type) or an undocumented swap. | • A firmware, filter or design change by Lufft (cf. the generation trend) can silently change κ. Continuous monitoring is mandatory. |
| • Feedback to Lufft (filter temperature control) could fix it at the source for new units. | • Interaction with the existing CHM15k overlap temperature correction must be checked to avoid double-counting in the near range. The overlap model is anchored above ~800 m, so it should carry no gain term, but it needs verifying. |
| • A publishable result (first day/night test of ceilometer calibration with atmosphere-free pairs). It should go into the paper as a limitation / outlook. | • Daytime-heavy users may see apparent "changes" at go-live; communication is needed. |

### 8.3 Decision

* **Do now** (no external dependency, no product change):
  * calibration-side computation and publication of `temp_int`, κ_u and C_ref in the calibration
    NetCDF;
  * the κ QC;
  * the dashboard display;
  * the shadow L2 at the evaluation sites.
* **Ask Lufft** the four questions of §5 in parallel.
* **Switch the L1→L2 correction on** only after the shadow period meets the go/no-go criterion. Do it
  unit by unit (valid-κ units only), keeping κ = 0 as the default.
* **Do not** switch it on network-wide with a single generation-median κ. The unit spread (0…+3.5 %/K)
  is too large: a common κ would over-correct a third of the units.

## 9. Caveats

* `temp_int` is reported in whole degrees. Within-night ranges are only ~1–2 K, so per-unit κ carries
  ±0.2–0.4 %/K uncertainty.
* The temperature sensor is not the component. The coefficient depends on the time scale
  (hours < day < weeks); a lag model would fit better than a static κ.
* Schiphol unit positions are unknown (common station coordinates in L1 and L2). The 2.5–4 km band
  and the exact closure argue against an atmospheric cause, but a coastal effect on 1–1.5 km cannot be
  excluded.
* The CHM15k/CL61 test assumes a stable CL61 and a day-night-invariant 910/1064 nm Ångström ratio.
  It is only sensitive to ≳ 10 % effects. The Payerne mismatch is explained in §7.1.
* The network arm covers 2025–2026 only; the archive arm covers 13 stations, 2013–2024.

## 10. Reproduction

```
python validation/chm15k_stability/hk_cache.py --workers 16      # ~70 min, I/O-bound (D: SATA SSD)
python validation/chm15k_stability/q2d_cross_cl61.py --build      # CL61 caches + CHM/CL61 test
matlab: validation/chm15k_stability/export_11yr_archive.m         # 11-yr archive -> CSV
python validation/chm15k_stability/q1_temperature.py
python validation/chm15k_stability/q2a_intranight.py
python validation/chm15k_stability/q2b_pairs.py
python validation/chm15k_stability/q2c_hk_diurnal.py
python validation/chm15k_stability/q4_mechanism.py
python validation/chm15k_stability/q3_synthesis.py                 # needs q1 + q2a + q2c outputs
python validation/chm15k_stability/q6_payerne.py                   # Payerne investigation
python validation/chm15k_stability/q7_tlom.py                      # temp_lom
python validation/chm15k_stability/q5_aeronet.py --aeronet-dir <dir>   # AERONET L1.5 daily AOD (aeronet_match.json)
python validation/chm15k_stability/q10_curvature.py                # curvature (laser vs filter)
matlab: validation/chm15k_stability/compare_ml_vs_exp.m            # ML (v4 bagged tree) vs exp model
python validation/chm15k_stability/q9_validation_tcorr.py          # paper-v0.9 validation with/without T correction
```

Tables: `q1_slopes.csv`, `q2a_intranight.csv`, `q2b_pairs.csv`, `q2c_hk_diurnal.csv`,
`q2d_cross_cl61.csv`, `q3_synthesis.csv`, `q3_daytime_excess.csv`, `q3_normalisation_gain.csv`,
`q4_mechanism.csv`, `q4_pairs_background.csv`, `q5_aeronet.csv`, `q6_payerne_triangle.csv`,
`q6_daytime_kappa.csv`, `q7_tlom_night.csv`, `q7_tlom_day.csv`, `q9_validation_tcorr.csv`, `q10_curvature.csv`,
`ml_vs_exp_scores.csv`, `ml_vs_exp_per_station.csv`, `ml_vs_exp_predictions.csv` in
`C:/DATA/Projects/202610_CHM15k_stability/`.

### References

* CIPM (2007): Recommended values of standard frequencies, 127I2 a10 R(56) 32-0 at 532 nm
  (frequency-doubled Nd:YAG), https://www.bipm.org/documents/20126/41549514/M-e-P_I2_532.pdf
* Kushida, T. (1969): Linewidths and thermal shifts of spectral lines in neodymium-doped yttrium aluminum
  garnet and calcium fluorophosphate, Phys. Rev., 185, 500–508.
* Jenoptik Laser, Optik, Systeme GmbH (1996): JENOFILT 602 bandpass filter data sheet, doc.
  fil601-99-14-0196-en, https://www.jenoptik.com/-/media/websitedocuments/optics/coatings/jenofilt/filt_602.pdf
* Tsai, R.-Y., et al. (2001): Thermally stable narrow-bandpass filter prepared by reactive ion-assisted
  sputtering, Appl. Opt., 40, 1593–1598, https://doi.org/10.1364/AO.40.001593
* Zenteno-Hernández, J. A., et al. (2024): On the temperature stability requirements of free-running
  Nd:YAG lasers for atmospheric temperature profiling through the rotational Raman technique, Atmos.
  Meas. Tech., 17, 4687–4694, https://doi.org/10.5194/amt-17-4687-2024
* Blifford, I. H., Jr. (1966): Factors affecting the performance of commercial interference filters,
  Appl. Opt., 5, 105–111, https://doi.org/10.1364/AO.5.000105

* Hervo, M., Poltera, Y., and Haefele, A.: An empirical method to correct for temperature-dependent
  variations in the overlap function of CHM15k ceilometers, Atmos. Meas. Tech., 9, 2947–2959, 2016.
* Roşu, B., Roşu, A., Voiculescu, M., Constantin, D.-E., Belegante, L., and Nemuc, A.: Thermal behavior
  of CHM15k ceilometers and its statistical association with attenuated backscatter: a multi-site
  ACTRIS study, Preprints, doi:10.20944/preprints202610.0850.v1, 2026 (not peer-reviewed).
* Madonna, F., Rosoldi, M., Lolli, S., Amato, F., Vande Hey, J., Dhillon, R., Zheng, Y., Brettle, M., and
  Pappalardo, G.: Intercomparison of aerosol measurements performed with multi-wavelength Raman lidars,
  automatic lidars and ceilometers in the framework of INTERACT-II, Atmos. Meas. Tech., 11, 2459–2475, 2018.
* Madonna, F., Amato, F., Vande Hey, J., and Pappalardo, G.: Ceilometer aerosol profiling versus Raman lidar
  in the frame of the INTERACT campaign of ACTRIS, Atmos. Meas. Tech., 8, 2207–2223, 2015.
* Bass, M., Weichman, L. S., Vigil, S., and Brickeen, B. K.: The temperature dependence of Nd3+ doped
  solid-state lasers, IEEE J. Quantum Electron., 39, 741–748, 2003.
* Sato, Y. and Taira, T.: Temperature dependencies of stimulated emission cross section for Nd-doped
  solid-state laser materials, Opt. Mater. Express, 2, 1076–1087, 2012.
* Kimmelma, O., Tittonen, I., and Buchter, S. C.: Thermal tuning of laser pulse parameters in passively
  Q-switched Nd:YAG lasers, Appl. Opt., 47, 4262–4266, 2008.
* Degnan, J. J.: Optimization of passively Q-switched lasers, IEEE J. Quantum Electron., 31, 1890–1901, 1995.
* Hopkin, E., Illingworth, A. J., Charlton-Perez, C., Westbrook, C. D., and Ballard, S.: A robust automated
  technique for operational calibration of ceilometers using the integrated backscatter from totally
  attenuating liquid clouds, Atmos. Meas. Tech., 12, 4131–4147, 2019.
* Le, V., O'Connor, E. J., Filioglou, M., and Vakkari, V.: Operational performance of the Vaisala CL61
  ceilometer for atmospheric profiling, Atmos. Meas. Tech., 19, 4923–4941, 2026.
* Van Hove, M. and Diémoz, H.: Seasonal variation in the Rayleigh calibration factor of Automatic
  Lidar-Ceilometers: amplitude across Europe and possible explanations, PROBE COST report, Zenodo,
  doi:10.5281/zenodo.11074353, 2024.
* Buxmann, J.: Investigating the seasonal fluctuations of the CHM15K ceilometer calibration constant,
  PROBE report, Zenodo, doi:10.5281/zenodo.11108621, 2024.
