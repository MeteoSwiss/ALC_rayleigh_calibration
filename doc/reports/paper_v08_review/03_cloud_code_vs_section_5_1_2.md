# Adversarial review — paper v0.8 Section 5.1.2 "Cloud", Appendices 3–4 and the WV part of Appendix 1 vs the liquid-cloud calibration code

Branch `remote-dark` (HEAD 18876dd), reviewed 2026-10-10. Nothing in the repository was modified. Line numbers refer to the files as they are on this branch.

Files read: `calibration/cloud/calibration.py`, `calibration/cloud/_filters.py`, `calibration/water_vapor_correction/water_vapor.py`, `calibration/water_vapor_correction/cloud_water_vapor.py`, `calibration/config.py`, `calibration/flags.py`, `calibration/main.py` (Rayleigh CLI only — the cloud path lives in `scripts/run_network_calibration.py::_do_cloud`), `calibration/io/data_loader.py`, `calibration/io/instrument_day.py`, `monitoring/kalman.py`, `validation/multiple_scattering_eta.py`, `validation/_showcase_cloud_figure.py` (Figure 4), `validation/_cloudnet_droplet_climatology.py` + `validation/_ms_cloudnet_tests.py` (Figure A4-1), `validation/paper/run_paper_v08.py`, `validation/paper/fig_method_comparison.py`, `options.json`, `ops/config.sh`, `doc/reports/02_cloud_calibration.md`, `04_multiple_scattering.md`, `doc/WATER_VAPOR_AUDIT.md`, O'Connor et al. (2004) pp. 777–781, Hopkin et al. (2019) pp. 4131–4140, Wiegner et al. (2014) Table 1; figures image7.png (Fig. 4), image28.png (Fig. A3-1), image29.png (Fig. A4-1).

Two read-only computations were run with the project venv to settle numerical questions (scripts in the scratchpad, nothing written to the repo):

* `eta_check.py` — interpolates the operational η tables at the Table A4-1 heights and re-derives them with `validation.multiple_scattering_eta.eta_table`.
* `run_palaiseau_check.py` / `aer_split.py` — runs the **operational** cloud pipeline of this branch on the Figure-4 day (Palaiseau CL31, L1 file `D:/E-PROFILE_L1_2026/0-250-1001-07151/2026/06/L1_0-250-1001-07151_A20260612.nc`, CAMS `D:/CAMS/CAMS_Beta_20260612.nc`) and splits the sub-cloud aerosol integral into its parts.

---

## 0. Headline findings (ranked by impact on the paper)

1. **Sign error in the multiple-scattering sentence (Section 5.1.2 and Appendix 4).** The paper says uncorrected multiple scattering "biases the integrated backscatter low and the calibration coefficient high". From O'Connor's own Eq. 6, B = 1/(2ηS) with η < 1, the uncorrected integral is biased **high** and the O'Connor coefficient C **low** (the lidar constant C_L = 1/C high). The code confirms it: `beta * factor_profile` with η < 1 (`_filters.py:86`) *reduces* B and *raises* C by 1/η (≈ +12 % for a CL31 at 0.6 km). The sentence is self-contradictory as written ("enhances the apparent backscatter … biases the integrated backscatter low").

2. **The CL61 water-vapour spectrum in the code is still 910.74 nm / 1.0 nm, not the 910.55 / 0.188 nm the paper declares operational** (`water_vapor.py:35-41`: `"CL61": (910.74, 1.0)`). The manufacturer spectrum is only reachable through the `ALC_WV_SPECTRUM` environment override (`water_vapor.py:49-79`), which is **commented out** in `ops/config.sh:34` (the whole v2.2 cut-over block: "Do NOT uncomment before scripts/check_v22_archive.py passes") and set only in the CSCS batch scripts (`ops/cscs/alc_v22_release.sbatch:47`, `_payloads`, `_sens`, `_addons`). Whether the daily server run uses the paper's spectrum therefore depends on server-side configuration that the repository does not show as enabled. (`validation/paper/run_paper_v08.py:5-10` states that the dashboard "v2.2" constants used in Section 6 were computed with 910.55/0.188 — consistent with the CSCS reprocessing, not with the committed default.)

3. **The aggregation described in 5.1.2 ("aggregated monthly … median of all individual monthly medians … inter-monthly standard deviation") does not exist in the code.** The product is a **daily median** of the per-profile C (`calibration.py:794`, `run_network_calibration.py:422`), turned into C_L = C_applied / C (`calibration.py:825`), with uncertainty C_L·σ/median (`run_network_calibration.py:447`); the dashboard then runs a **Kalman filter on the daily cloud constants** (`run_network_calibration.py:497-520 _kalman_rows`, for `method in ("rayleigh", "cloud")`, via `monitoring/kalman.py::kalman_best_estimate`). Section 6.3/Appendix 6 use that daily Kalman cloud C_L (`run_paper_v08.py:5-7`), Section 6.1 uses ±15-day medians of the daily values (`fig_method_comparison.py:8-9,32`). The flow-chart A3-1 itself says "Daily median C → Kalman-smoothed calibration constant", contradicting the 5.1.2 text. Section 5.1.3 describes the Kalman for Rayleigh only.

4. **Figure 4's headline value is not what the operational pipeline produces for that day.** Running the operational code on 2026-06-12 gives **n = 248, C = 1.29** (C_L = 7.73e7), whereas the figure (reproduced exactly by `_showcase_cloud_figure.py`) gives **n = 210, C = 1.39**. The figure (i) omits the sub-cloud aerosol two-way transmission correction that 5.1.2 presents as part of the method (−7 % on this day: per-profile T²_aer median 0.932, range 0.871–0.958), and (ii) runs the cloud tests on β' *without* the η correction, whereas the code applies η first (`calibration.py:761-770`) — hence the different n. The figure's "C_L = 0.72" is 1/C, i.e. C_L/1e8, not the Wiegner C_L the dashboard stores (7.7e7).

5. **"C = 1/C_L" holds only in normalised units.** The code defines C_L = calibration_constant_applied / C (`calibration.py:819-825`); the applied constant is the E-PROFILE default 1e8 for CL31/CL51 L1 (`calibration.py:294-300, 414`) and 1.0 for CL61. So C = 1/C_L is exact for the CL61 only; for CL31/CL51 it is C = 1e8/C_L. The 5.1.2 statement and Fig. 4(f) caption should say so (or define C_L as dimensionless there).

6. **The L2 quality-flag filter described in 5.1.2 is dead code.** `build_cloud_input` sets `quality_flag=None` (`calibration.py:327`) on the only reader path, so `apply_instrument_filters` never enters the quality-flag branch (`_filters.py:120`). The operational input is L1 anyway (`options.json: "data_level": "L1"`; runner `_l1_file`).

7. **The aerosol transmission correction is not "residual": it is a 5–10 % term with an unpropagated ±3–4 % assumption.** On the Figure-4 day (a clean calibration day) it lowers C by 7 % (median), up to 13 % per profile. The 20/70 sr bounds are computed (`_filters.py:361-367`) but discarded (`calibration.py:801`: `C_corrected, _C_low, _C_high = …`) and never enter the reported uncertainty. Its integral starts at gate index 4 (40 m, inside the <100 m zone the calibration otherwise excludes) and ends 5 gates below the **peak**, i.e. it includes the cloud base: on that day 8 % of B_aer comes from 40–100 m and 14 % from the cloud-base gates (T² 0.932 instead of 0.946 if only the genuine sub-cloud column were used, ≈ 1.5 % of C).

8. **Appendix 4 omits the in-cloud extinction assumption that dominates η.** The tables are reproduced bit-exactly by `eta_table(k, 5.5, alpha=10e-3)` (max deviation 0.00000), i.e. PVC at a_G = 5.5 µm **and α = 10 km⁻¹**, τ_max = 6, dr = 2 m. Sensitivity (CL31, 1.125 km): a_G 5.0→6.0 µm moves η by only ±1.2 % (0.833→0.813), α 10→30 km⁻¹ moves it by −16 % (0.823→0.695). O'Connor (2004) used 15–20 km⁻¹; report 04 §3a shows α = 10 km⁻¹ was *fitted* to the legacy Hopkin tables jointly with a_G = 8 µm, then a_G was revised to 5.5 µm from Cloudnet while α was kept. The paper lists FOV, divergence, wavelength and r_e as the inputs and never mentions α.

9. **Upper CBH limit equals the top of the integration window (2400 m = 2400 m).** `cbh_maxheight = cal_maxheight = 2400` (`calibration.py:127,129`). A cloud with base at 2.1–2.4 km has part (up to all) of its return above the integration top; the 5 % ratio test only catches peaks above ≈ 2.35 km. The paper's rationale "to avoid … insufficient integration at high altitudes" is not realised by equal limits; it biases C high for the highest accepted clouds (a few %, up to ~15 % for slowly decaying returns that just pass the 1/20 test).

10. **Minor table/number issues.** Table A4-1's CL61 row is the CL51 row (code: 0.924/0.875/0.807 → 0.92/0.88/0.81, marginally weaker than CL51 because of its larger divergence, `_filters.py:65-66`); CHM15k at 1 km is 0.9450 (0.94/0.95 borderline); CHM15k and Mini-MPL are **not** identical in the code (Δη = 0.001–0.005), only at two decimals — physically because λ and FOV both halve so ρ_t/θ is 0.0037 vs 0.0036. Appendix 4 uses the Cloudnet *effective* radius r_e directly as Hogan's *equivalent-area* radius a_G (`multiple_scattering_eta.py:13, 222-223`; report 04 l.475); for a stratocumulus gamma distribution r_e/a_G ≈ 1.1–1.2, worth a sentence.

---

## 1. Statement-by-statement table

Verdict key: MATCHES / DIFFERS / NOT IN CODE / NOT IN PAPER. "Server config" marks items that depend on the environment of the machine running the daily cron.

| # | Paper statement (5.1.2 / App. 1 / App. 3 / App. 4) | Code location | Verdict | Actual value / behaviour |
|---|---|---|---|---|
| 1 | S_liq = 18.8 sr (O'Connor 2004) | `_filters.py:17` `S_THEORETICAL = 18.8` | MATCHES | Hopkin uses the same value for 910 and 1064 nm. |
| 2 | B = ∫₀^{z_cb} β' dz = 1/(2 S_liq) | `_filters.py:151-184` (`calculate_lidar_ratio`), limits `calibration.py:126-127` | DIFFERS (limits) | Trapezoid over the gates in **[100, 2400] m** (`cal_minheight`/`cal_maxheight`), neither from 0 nor *to the cloud base* — the integral must run *through* the cloud (O'Connor Eq. 4: 0→∞). Flow-chart A3-1 repeats "∫_{0.1 km}^{z_cb}". Profile skipped if <3 valid gates or >10 % NaN; S = 1/(2B) only if B > 0. |
| 3 | C = 1/(2 S_liq ∫β'_L2 dz), β_true = C β'_L2 | `calibration.py:776` `calibration_coefficients = S_consistent / S_THEORETICAL`; `_filters.py:183` `S = 1/(2B)` | MATCHES (formula) | C = S_app/18.8 = 1/(2·18.8·B). Operational input is **L1** (CL31/51: rcs_0/1e8, `calibration.py:294-300,414`; CL61: physical β_att), not L2 — "β'_L2" is only notation. |
| 4 | "C is the reciprocal of C_L"; C > 1 ⇒ instrument under-estimates β | `calibration.py:812-825`: `calibration_factor = 1/cal_median`; `lidar_constant = cc_applied / cal_median`; runner `run_network_calibration.py:442-452` | DIFFERS (partly) | Stored product = Wiegner **C_L = calibration_constant_applied / C** (`cal_value`), uncertainty = C_L·σ_C/median. C = 1/C_L only when the applied constant is 1 (CL61); CL31/CL51: C = 1e8/C_L. The sign statement (C > 1 ⇒ file's β too small) is correct for the O'Connor multiplier. Fig. 4(f) "C_L = 0.72" is 1/C (= C_L/1e8); the operational C_L for that day is 7.7e7. |
| 5 | MS: β'_corr = η(z) β'; PVC (Hogan 2006) from FOV, divergence, λ, droplet r_e | `_filters.py:59-86` (`beta * factor_profile[:, None]`, table interpolated on `range_km` with linear **extrapolation**); tables `_filters.py:20-56`; derivation `multiple_scattering_eta.py:66-233` | MATCHES (mechanics) | η(range) multiplies the whole profile, aerosol included, *before* the filters (`calibration.py:761`). Extrapolated below 0.25 km: η(0.1 km) = 0.959 (CL31) — a MS "correction" of the aerosol layer, ≤0.3 % effect. Gate-wise η(z) vs η(CBH) differ by 0.05 % on the Figure-4 day (checked). |
| 6 | Uncorrected MS "biases the integrated backscatter low and the calibration coefficient high" | — | DIFFERS (sign) | See §2a: B_meas = 1/(2ηS) > 1/(2S) → B **high**, C **low** (C_L high). Code raises C by 1/η. Appendix 4 repeats the error. |
| 7 | Table A4-1 η: CL31 0.89/0.84/0.76; CL51 0.92/0.87/0.80; CL61 0.92/0.87/0.80; CHM15k 0.97/0.95/0.91; Mini-MPL 0.97/0.95/0.91 | `_filters.py:20-49`, interpolated as the code does (`_interp1_linear_extrap`) | MATCHES (CL31, CL51, Mini-MPL); DIFFERS slightly (CL61, CHM15k) | Code at 0.5/1/2 km: CL31 0.894/0.835/0.759; CL51 0.922/0.873/0.803; **CL61 0.924/0.875/0.807** (→ 0.92/0.88/0.81; the paper's CL61 row is a copy of the CL51 row); **CHM15k 0.969/0.9450/0.906**; Mini-MPL 0.971/0.948/0.909. CHM15k − Mini-MPL = −0.001…−0.005 at the nodes: not identical, only at two decimals (ρ_t/θ = 0.0037 vs 0.0036). |
| 8 | Table A4-1 optics (half-angles): FOV 0.83/0.56/0.56/0.23/0.11 mrad; div ±0.40×±0.70, ±0.15×±0.25, ±0.20×±0.35, 0.15, 0.055; λ 910/910/910.55/1064/532 | `multiple_scattering_eta.py:205-215` `OPTICS` (ρ_t, ρ_l half-angles, 1/e), `WAVELENGTH`; comments 36-47 | MATCHES | ρ_t are half-angles (Wiegner 2014 Table 1 "RFOV (half angle)": CL31 0.83, CL51 0.56, CHM15k 0.23 ✓; CL61 = CL51 per user guide; Mini-MPL 220 µrad full → 0.11). The two divergence half-angles are combined as their RMS (0.57, 0.21, 0.28 mrad) and used as a 1/e Gaussian half-angle — not stated in the paper. λ = 910.55 only in the η table (WV uses `LASER_SPECTRUM`). |
| 9 | η evaluated "for the droplet effective radius of the calibration clouds" (≈5.5 µm) | tables == `eta_table(k, 5.5, alpha=10e-3)` bit-exactly | MATCHES (5.5 µm) / NOT IN PAPER (α = 10 km⁻¹, dr = 2 m, τ_max = 6; a_G is the equivalent-area radius, fed with the Cloudnet effective radius) | α is the dominant sensitivity (§2e). The derivation script's "validation" step (`multiple_scattering_eta.py:60-62, 253-275`) now fits the model to its own output, since `_ETA_CL31/_ETA_CL51` imported "to validate against" *are* the PVC tables since 2026-07. |
| 10 | WV spectra: CL31 909.7 nm / ~6 nm; CL51 910.0 / 3.4; CL61 910.55 / 0.188 (App. 1: manufacturer spectrum "used operationally") | `water_vapor.py:35-41` `LASER_SPECTRUM`; `laser_spectrum_for` 69-79; `set_defaults` `calibration.py:222-238`; override `ops/config.sh:34` (commented), `ops/cscs/alc_v22_*.sbatch` | MATCHES (CL31, CL51); **DIFFERS for CL61 — server config** | Committed default **CL61 = (910.74, 1.0)**; 910.55/0.188 only via `ALC_WV_SPECTRUM`, which is commented out in `ops/config.sh`. Spectrum is one source of truth for Rayleigh and cloud (`laser_spectrum_for`). `ALC_WV_DISABLE=1` switches the correction off (runner l.390). |
| 11 | T²_WV from CAMS forecast T/q collocated to the station + HITRAN LUT "spectrally averaged absorption cross-section … convolved with the laser spectrum" | `cloud_water_vapor.py:134-282`; CAMS `ensure_cams_file` 183-192 (0.4° monthly/daily, 1° fallback), `cams_point_too_far` 194-201, `cams_levels_all_times` 202 (L137, ideal-gas n_w); LUT `water_vapor.py:85-123` (903–918 nm, 0.0083 nm; 1000 heights at 10 m) | DIFFERS (description) | The code does **not** average the cross-section: it computes τ(λ,z) per LUT wavelength (262-267) and takes the emission-weighted mean of the monochromatic two-way transmissions, T²_eff = Σ g(λ) e^{−2τ(λ)} / Σ g (268-269, Gaussian σ = FWHM/2.355, ±3σ band) — the Wiegner & Gasteiger (2015) effective transmission, which is the correct quantity (Jensen: T²_eff ≥ e^{−2τ̄}, matters for the 6 nm CL31 line). CAMS: "forecast" product (`scripts/download_cams_oob.py:98`), nearest grid point, **nearest 3-h step** (271-278). |
| 12 | β'_corr = β'_corr / T²_WV | `calibration.py:755` `data.beta = data.beta / trans2` | MATCHES | Direction correct (§2b). Gates where T² is NaN/≤0 are silently set to 1 (`cloud_water_vapor.py:281`). |
| 13 | "a night without a matching CAMS profile is left un-calibrated" | `calibration.py:739-756` (raise, no fallback); runner 480-493 → flags −4 / −10 / −99; `flags.py:33,39` | MATCHES | It is the whole **UTC day** (`slice_to_date`, runner l.35 "full daytime of D"), not a night. −10 when the station is outside the CAMS domain, −99 for any other WV failure (missing LUT, bug). |
| 14 | Health: discard if the **L2 quality flag** is invalid in the integration column | `calibration.py:327` `quality_flag=None`; `_filters.py:119-129` | NOT IN CODE (dead) | No reader path ever fills `quality_flag`; the branch cannot execute. Operational input is L1. |
| 15 | L1 laser energy < 90 % of nominal → discard | `_filters.py:140-146` (`le < 90`), `energy_threshold` `calibration.py:143`; loader `data_loader.py:330` (`laser_energy`/`laser_pulse_energy`, NaN for CHM15k/Mini-MPL) | MATCHES | NaN → kept (gate inactive on photon counters and on L2). Block-mean over the 30 s working grid. |
| 16 | Window transmission: no correction; discard below 50 % | `_filters.py:131-138`, `window_correction_threshold = 50` `calibration.py:142` | MATCHES | Reject-only gate as described; NaN kept. (Hopkin: 90 %.) Legacy `window_threshold = 90` and `temp_threshold = −20` fields (`calibration.py:133,136`) are unused. |
| 17 | Test 1 "totally attenuating": signal 300 m above the peak ≤ 1/20 of the peak | `_filters.py:196-198` (`gate_300m = round(300/dr)` = 30 gates), 219-241 | MATCHES | Peak = `nanargmax` of the η-corrected profile with gates ≤ 100 m zeroed (l.222, verified equal to MATLAB `beta_roi(1:lower_gate)=0`, no off-by-one), **no upper bound** (η extrapolation suppresses high gates). The test samples one gate exactly 300 m above the peak. Reject if `beta_above * 20 > max_beta`. |
| 18 | Test 2 "isolated base": 300 m below ≤ 1/20 | `_filters.py:243-249` | MATCHES | Same single-gate test. |
| 19 | Test 3 "low aerosol": sub-cloud integral < 5 % of the 100–2400 m integral | `_filters.py:252-267`, `ratio_filter = 0.05` | MATCHES (numbers) / DIFFERS (what "beneath the cloud" is) | Plain sums × dr. "Beneath the cloud" = from 100 m to **5 gates (50 m) below the β peak**, not below the CBH, so the lowest cloud gates count as "aerosol"; the denominator includes the cloud. It is the dominant rejection reason (Figure-4 day: 963 profiles). Median value on accepted profiles 2.4 %. |
| 20 | Test 4: cloud base 500–2400 m a.g.l. | `_filters.py:269-277`; `cbh_minheight/maxheight` `calibration.py:128-129`; CBH source `calibration.py:303-308` (`cbh[:, 0]` lowest layer), `data_loader.py:704-728` (block **minimum** over 30 s) | MATCHES | Instrument-reported lowest layer (Vaisala = backscatter maximum, median +10 m above the peak on the Figure-4 day); falls back to the peak height when no CBH. Upper limit = integration top (§2f-ii). |
| 21 | Test 5: 5 consecutive profiles within ±10 % of the group mean | `_filters.py:282-306`; `n_consecutive = 5`, `consistency_range = 10` | MATCHES | Applied to S (≡ C/18.8 before the aerosol correction). "Consecutive" = adjacent columns of the **30 s** working grid (`instrument_day.py:287-288`: 30 s / 10 m) → a 2.5-min window; gaps are not detected. Overlapping windows do not double-count (`S_consistent[i] = S[i]` once). Hopkin: three profiles either side and ≥ 10 profiles/day. |
| 22 | Aerosol two-way transmission: B_aer = ∫_{below cloud} β' dz, C_corr = exp(−2 S_aer B_aer), S_aer = 50 sr, bounds 20/70 sr "to characterise the uncertainty" | `_filters.py:309-369`; `aerosol_lidar_ratio = 50` runner l.391; `calibration.py:800-808` | MATCHES (formula, 50 sr) / NOT IN CODE (bounds in the output) | C_corrected = C_base · exp(−2·50·B_aer) with B_aer = C_base·Σβ[**gate 4 … peak−5**]·dr (40 m → 50 m below the peak; peak = `nanargmax` of the *full* profile, no 100 m floor; skipped if peak index < 9). C_low/C_high computed and **discarded** (`_C_low, _C_high`); the reported uncertainty is only C_L·σ/median. Figure-4 day: factor 0.932 median (0.871–0.958). |
| 23 | Integration limits | `calibration.py:126-127` | fixed 100–2400 m | Not from the ground, not from the CBH, not to the cloud top (Hopkin: below 200 m and above 2.4 km not used). |
| 24 | "computed for each accepted profile and aggregated monthly … median of all individual monthly medians … inter-monthly standard deviation" | `calibration.py:793-808` (mean/median/std ddof=1/mode per **day**); runner 420-452; `_kalman_rows` 497-520; `monitoring/kalman.py` | NOT IN CODE | Daily **median** of per-profile C (Hopkin: mode) → C_L = applied/median → per-stream daily CSV + NetCDF (method 1, version 1000). Dashboard Kalman per (stream, method) on the daily cloud C_L (constant/random-walk predict model, rolling-IQR outliers, relative process noise 4 %/d, 15 %/yr). No monthly statistic anywhere in the operational chain. |
| 25 | Section 6 uses "daily Kalman best-estimate constants" for the cloud method | `_kalman_rows` for `method == "cloud"` when ≥ 5 successful days; `run_paper_v08.py:5-7`; `fig_method_comparison.py:32` (±15 d medians) | MATCHES (Section 6) / NOT IN PAPER (methods) | A Kalman on cloud constants exists and is what 6.3/App. 6 consume; 5.1.3 says "Daily Rayleigh calibration constants are combined…" only; 5.1.2 says monthly. The flow-chart A3-1 says "Daily median C → Kalman-smoothed". |
| 26 | Fig. 4(e): "theoretical ηS curve computed from the droplet effective radius measured by Cloudnet that day" | `_showcase_cloud_figure.py:44-66, 106-126` (Cloudnet `der` median over CBH−30…+300 m at the valid scenes; `eta_for_cloud_base(…, alpha=10e-3, dr=5.0)`) | MATCHES (figure) / NOT ACKNOWLEDGED | The operational code uses the fixed 5.5 µm table; for this day the Cloudnet r_e happens to be 5.5 µm, so the curves coincide, but the caption suggests a per-day r_e. The figure's per-profile C also uses η at the CBH (l.111-113) rather than gate-wise (0.05 % difference, checked). |
| 27 | Fig. 4: n = 210 profiles, C = 1.39, "C_L = 0.72" | `_showcase_cloud_figure.py:91-118` | DIFFERS from the operational pipeline | Operational code on the same file: **n = 248, C = 1.294, C_L = 7.73e7**; without the transmission correction C = 1.393. The figure skips the aerosol transmission correction and runs the filters on un-η-corrected β'. "C_L = 0.72" = 1/C = C_L/1e8. Panels (d)/(e) plot S_app/C_med ("S_phys", l.118), i.e. the apparent S *after* dividing by the day's median C. |
| 28 | App. 4: calibration clouds have r_e "close to 5.5 µm … consistent across the three sites" | `_cloudnet_droplet_climatology.py:49-67` (`match_cache.npz`), built by `_ms_cloudnet_tests.py:85-95` | MATCHES | Medians of Cloudnet `der` over [CBH−30 m, CBH+300 m] at the nearest Cloudnet time (≤ 0.1 h) for profiles that passed all cloud tests (Payerne CL61+CL31, Palaiseau CL31, Lindenberg CL61): **5.8 / 5.6 / 5.4 µm** (Fig. A4-1), "5.5" = median of medians; "cloud base" 7.7/7.3/5.7; "all" 9.1/8.5/6.9. Caveats: r_e (effective) ≠ a_G (equivalent-area); the subset is conditioned on calibration success. |
| 29 | (Hopkin) minimum 10 profiles per day | `flags.py:116-130` `min_profiles_success = 3` | NOT IN PAPER | One 5-profile window (2.5 min) already yields a flag-1 daily constant; n is 0 or ≥ 5, so the 0.5 "partial" flag is unreachable for the cloud method. |
| 30 | Rejection flags −20…−26, −6, CBH-based −1 override | `flags.py:43-49, 55-90`; runner 424-441 | NOT IN PAPER | Dominant-filter attribution when no profile survives; −6 when σ_C > median C; day forced to −1 if the instrument saw no cloud. |

---

## 2. Physics / sign checks

### (a) Multiple scattering — the paper has the sign wrong

O'Connor (2004) Eqs. 3–6: dτ = η S β_true dz, B = ∫β_obs dz = 1/(2ηS), η ≤ 1 ("values between 0.5 and 1"). Multiple scattering *reduces the apparent extinction* of the cloud, so the observed signal is extinguished more slowly and the integrated backscatter is **larger** than the single-scattering value: B_meas = B_ss/η > B_ss. Hence, uncorrected,

* B is biased **high** by 1/η,
* C = 1/(2 S_liq B) is biased **low** by η (CL31 at 0.6 km: η ≈ 0.88 → C 12 % too low),
* the lidar constant C_L = C_applied/C is biased **high**.

The code does exactly this: `return beta * factor_profile[:, None]` with η < 1 (`_filters.py:86`) shrinks the integral and raises C. The paper sentence (5.1.2 and the first paragraph of Appendix 4) — "this multiply-scattered contribution enhances the apparent backscatter and, if left uncorrected, biases the integrated backscatter low and the calibration coefficient high" — is wrong on both counts (and internally inconsistent: an enhanced backscatter cannot integrate to a smaller B). It would be right for C_L, not for C. Suggested wording: "…biases the integrated backscatter high and the calibration coefficient C low (the lidar constant C_L high)".

### (b) Water vapour — direction right, altitude bookkeeping almost right

`data.beta = data.beta / trans2` with T² < 1 raises β', raises B, lowers S_app and C (lowers C_L): correct — WV absorption makes the file's β' too small. Height handling in `compute_wv_transmission`:

* n_w profile: CAMS geopotential height is converted to AGL by `z_agl = cams_z[:, i] - data.station_altitude` (l.245) and interpolated onto the instrument range grid (AGL) — **correct**, with the lowest-level value filled below the model surface (248-254) and 0 above the top.
* τ is accumulated from the first range gate (tau = 0 at `range_col[0]`, 262-267) — correct for a ground-based lidar; the ≤10 m offset is irrelevant.
* **The LUT is indexed by range AGL, not altitude ASL** (`height_indices = argmin |range_col − abs_cs_height|`, 229-232), while the LUT heights (0.005–9.995 km) parametrise the pressure/temperature broadening of the cross-section. Measured on the bundled LUT: the spectrum-mean cross-section changes with altitude by +0.7 % (CL31) and +0.4 % (CL51) at 2 km, but by −3.7 % at 0.5 km, −6.8 % at 1 km and −12.7 % at 2 km for the narrow CL61 line (910.55/0.188). So for a CL61 at ~500 m a.s.l. (Payerne, Aosta, Sion) τ_WV is overestimated by ~4 %; since the CL61 two-way transmission at 2 km is ~0.9, the effect on C is < 0.5 % — a real inconsistency but numerically negligible. (The 910.74/1.0 spectrum would be insensitive, −0.7 %.)
* Tilt: no zenith-angle handling anywhere in the cloud code; range is used as height. This is fine for B (O'Connor's integral is along the line of sight) and for the η tables (range-based); it makes the WV path ~1/cos θ longer than assumed (2 % for 12°), i.e. ≤ 0.5 % in T² — negligible.
* Time: the nearest 3-hourly CAMS step is used (`_interp1_nearest_extrap_cols`), no interpolation. Fine, not stated in the paper.

### (c) Integral: dz, trapz vs sum, 1e6

* The calibration integral is a trapezoid on the working grid: **10 m** (`instrument_day.py:287 target_range_m = 10.0`; the runner passes `average_range_m=0` because the shared read is already coarsened, `run_network_calibration.py:396-397`). For the L2 path it would be the 30 m L2 grid (dz from `np.diff(x)`, so any grid is handled). The sub-cloud and aerosol integrals use plain sums × dr (`_filters.py:258-261, 344`) — immaterial.
* Units: for L2 the loader folds the 1e-6 (`rcs = attbsc_0·cc·1e-6`, then `beta = rcs/cc`, `calibration.py:288`) so β is physical m⁻¹ sr⁻¹ and S_app ≈ 18 sr; for L1 CL31/CL51 `beta = rcs_0 / 1e8` (`calibration.py:294-300`, `INSTRUMENT_CAL_DEFAULT`, l.414); CL61 L1 is already β_att (factor 1). Consistent: B_total on the Figure-4 day = 0.0266 sr⁻¹ = 1/(2·18.8) after calibration.
* The aerosol transmission correction uses the calibrated integral (`B_aerosol = C_base * B_aerosol_raw`, `_filters.py:352-353`), so AOD is physical; the comment at 345-351 about L2 units is stale but harmless.
* Negative β' (Vaisala cosmetic offset, Kotthaus 2016) is integrated as-is over the ~1.5 km of noise above the cloud inside the 100–2400 m window — typically < 0.5 % of B.

### (d) Figure 4 vs the operational code

Reproduced exactly by `validation/_showcase_cloud_figure.py` (its n = 210, C = 1.3865 with a 5.5 µm η at the CBH). The operational pipeline of this branch on the same L1 file gives n = 248, C = 1.2945 (transmission correction on) / 1.3934 (off). Differences:

1. The figure never calls `apply_transmission_correction` (−7 % median on this day; the paper presents the correction as part of the method).
2. The figure applies the cloud tests to β'/T²_WV without η (`S_app_all` from `beta_wv`, l.93-96), the code applies η first — different accepted set (210 vs 248).
3. η at the CBH vs gate-wise: 0.05 %, irrelevant.
4. The theory curve uses the Cloudnet r_e of the day (5.5 µm) with α = 10 km⁻¹ and dr = 5 m; the operational table has the same a_G and α (dr = 2 m) — the caption does not say the operational correction is a fixed table.
5. "C_L = 0.72" in panel (f) is C_L/1e8; the dashboard reports 7.7e7 for this day.
6. The figure uses the 1° `D:/CAMS` daily files (flat site, negligible), fine, but note CLAUDE.md's warning about 1° monthlies.

### (e) Appendix 4 droplet radius and the η tables

* "Calibration profiles" in Fig. A4-1 = profiles that passed all five tests in the MS experiment runs (`results.jsonl` → `_ms_cloudnet_tests.build_jobs`), matched to Cloudnet `der` (median over CBH−30…+300 m, time within 0.1 h). The site medians are 5.8/5.6/5.4 µm; "5.5" is a median, not a mode, and "consistent across sites" is fair (spread 0.4 µm).
* a_G in Hogan's model is the equivalent-area radius √⟨r²⟩ (`multiple_scattering_eta.py:13`; report 04 l.475); Cloudnet `der` is r_e = ⟨r³⟩/⟨r²⟩ ≥ a_G. For a stratocumulus gamma distribution r_e/a_G ≈ 1.1–1.2 → a_G ≈ 4.6–5.0 µm; η(1.125 km, CL31) at 5.0 µm is 0.833 vs 0.823 at 5.5 µm (+1.2 %, C −1.2 %). Small but should be stated as an approximation (report 04 l.532 does).
* **α = 10 km⁻¹ is the hidden parameter.** η depends on the lateral spread accumulated over the penetration depth ~τ/α: α = 30 km⁻¹ (O'Connor used 15–20) gives η(1.125 km, CL31) = 0.695 instead of 0.823, i.e. C +18 %. Report 04 §3a: α = 10 km⁻¹ came out of a joint fit (with a_G = 8 µm) to the legacy Hopkin ladder; a_G was then replaced by the Cloudnet 5.5 µm but α was kept. The paper's statement that η is computed from FOV, divergence, λ and r_e is incomplete, and the uncertainty discussion ("η is most sensitive to droplet size for wide-FOV instruments") misses the dominant sensitivity. The consistency of the CL61 Rayleigh/cloud ratio in Section 6.1 (0.94) is partly tunable through α and should be reported with the α assumption.
* Table A4-1 "CL61 = CL51": the code deliberately gives the CL61 its own table (larger divergence → η +0.001…+0.004); the printed row should be 0.92/0.88/0.81.
* CHM15k and Mini-MPL are not identical (Δη up to 0.005); the two-decimal coincidence is physical (λ/FOV ratio) and worth one sentence rather than a silent identity.

### (f) Bug / weakness checklist

| Item | Finding |
|---|---|
| Off-by-one in the 300 m tests | None. `gate_300m = round(300/dr)` = 30 on the 10 m grid, `idx_above = max_idx + 30`; ROI zeroing `beta_roi[:lower_gate+1] = 0` equals MATLAB `beta_roi(1:lower_gate)=0`. |
| Peak detection grid | Working grid 30 s / 10 m for L1 (`instrument_day.py:287-288`); the 30 m L2 grid gives 10-gate tests. The peak search has no upper bound (gates > 2.4 km only damped by the η extrapolation). |
| CBH source | Instrument lowest layer (`cbh[:, 0]`), **minimum** over each 30 s block (`_block_reduce_cloud_base`), fallback to the peak height. OK. The CBH gate is only an acceptance window (report 02 §1), fine. |
| **CBH upper limit = integration top** | `cbh_maxheight = cal_maxheight = 2400 m`. Return of clouds based at 2.1–2.4 km is partly outside the window; the ratio test only rejects peaks above ≈ 2.35 km. Biases C high in the top ~300 m of the CBH window (few %, up to ~15 % for a return that just meets the 1/20-at-300 m criterion). Recommend `cbh_maxheight ≤ cal_maxheight − 300 m` or integrating to peak + 300 m. |
| **Quality-flag filter dead** | `quality_flag=None` on the only reader path (`calibration.py:327`). |
| Tilt | Not handled; negligible (see 2b). |
| Negative β' | Integrated as-is; ≤ 0.5 %. |
| **5 % aerosol test when the integral is cloud-dominated** | Denominator includes the cloud, so the criterion is "sub-cloud < 5 % of (sub-cloud + cloud)". The numerator runs to 50 m below the *peak* and therefore contains the lowest cloud gates (14 % of it on the Figure-4 day, where the Vaisala CBH sits 10 m above the peak). It is the dominant rejection (963 of the cloudy profiles that day) and is partly a cloud-base-sharpness test, not an aerosol test. |
| **Transmission correction window** | Starts at gate index 4 (= 40 m, inside the overlap zone excluded from everything else) and uses `nanargmax` over the whole profile (no 100 m floor), inconsistent with the cloud-filter peak; 8 % of B_aer from 40–100 m that day. Hard-coded gate indices (`[4 : max_idx-5]`, skip if `max_idx < 9`) make the window grid-dependent (120 m on a 30 m L2 grid). |
| Transmission correction magnitude | −7 % median, −13 % worst profile on a clean day; up to −12.5 % at the 5 % admission limit with 50 sr (−5 % … −17 % over the 20–70 sr bounds). The bounds are not propagated. The corrected β' is η-weighted (η < 1 below the cloud) → ~0.5 % under-correction. |
| ±10 % group test on S vs C | Equivalent (C ∝ S) before the transmission correction; the correction is applied afterwards, so final per-profile C within a window can spread beyond ±10 %. Harmless. |
| Overlapping windows / duplicates | No double counting: `S_consistent[i:i+5] = S[i:i+5]` is idempotent; daily statistics are over unique profiles (`valid_idx`). |
| Minimum sample | One window (5 profiles, 2.5 min) → flag 1; partial flag 0.5 unreachable; Hopkin's ≥ 10 not enforced. |
| Daily statistic | Median (Hopkin: mode; both computed, `calculate_mode` reproduces MATLAB `histcounts`). |
| Legacy/unused config | `temp_threshold`, `window_threshold` unused; module docstring (`calibration.py:48-72`) still describes the retired RH round-trip WV chain (commit b3b838e moved to the ideal-gas route). |
| Mini-MPL / CHM15k | WV gated by wavelength (`in_water_vapor_band`, 900–920 nm), saturation warning only (`calibration.py:722-728`); η tables applied. Consistent with "included for reference only". |

---

## 3. What the code does that the paper omits, and vice-versa

### In the code, not in the paper
* Input is **L1** (`rcs_0` divided by the default E-PROFILE constant 1e8 for CL31/CL51; β_att for CL61) on a **30 s / 10 m** working grid shared with the other passes; the L2 path exists but is not operational. The stored product is the Wiegner **C_L = 1e8/C** (CL31/51) or 1/C (CL61), with uncertainty C_L·σ_C/median_C.
* Full **UTC day** (day and night) per calibration; daily **median** of per-profile C; a single 5-profile window suffices.
* **Kalman filter on the daily cloud constants** (dashboard: constant predict model, rolling-IQR outliers, relative process noise), consumed by Section 6.3 and Appendix 6; ±15-day medians in Section 6.1.
* Flag taxonomy: 1 / 0.5 (unreachable) / 0 no data / −1 no liquid cloud / −4 CAMS missing / −6 σ > median / −10 outside CAMS domain / −20…−26 dominant rejection (window, energy, above, below, ratio, CBH, consistency) / −99 exception; a day is forced to −1 if the instrument reported no cloud base.
* η: PVC with **α = 10 km⁻¹**, τ_max = 6, a_G = 5.5 µm used as equivalent-area radius, applied gate-wise to the whole profile with linear extrapolation outside 0.25–2.375 km; CL61 has its own table.
* WV: emission-weighted mean of monochromatic transmissions (not a mean cross-section); CAMS **forecast** model levels at the nearest 0.4° point with a 1° fallback, nearest 3-h step; LUT 903–918 nm at 0.0083 nm, indexed by AGL range; below-surface n_w fill; NaN transmissions set to 1; env switches `ALC_WV_SPECTRUM`, `ALC_WV_DISABLE`.
* Health gates are NaN-tolerant (no energy/window gate on photon counters or L2).
* Aerosol correction integrates from 40 m to 50 m below the peak; bounds computed and dropped.

### In the paper, not in the code
* The L2 quality-flag filter (dead code).
* Monthly aggregation / median of monthly medians / inter-monthly standard deviation.
* The 20–70 sr uncertainty characterisation (computed, never reported).
* The integral limits 0 → z_cb (code: 100–2400 m) — also on the flow-chart.
* "night" (code: day).
* The statement C = 1/C_L (true only for the CL61 convention).
* Operational CL61 spectrum 910.55/0.188 (code default 910.74/1.0 unless the environment override is set — server configuration).

### Flow-chart A3-1 specifically
Correct: the three corrections, the five tests with their numbers, S_aer = 50 (20–70) sr, C = S_app/S_liq, "daily median C → Kalman". Wrong/misleading: "∫_{0.1 km}^{z_cb}" (should be 0.1–2.4 km); "quality flag" in the health box (dead); "20–70 sr bounds" implies they are used.

---

## 4. Items that depend on server-side configuration

* **CL61 WV spectrum**: `ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'` — commented out in `ops/config.sh:34` (v2.2 cut-over block), exported in `ops/cscs/alc_v22_*.sbatch`. The committed default is 910.74/1.0. The paper's Section 6 constants come from the dashboard "v2.2" record, which `run_paper_v08.py` says was computed with 910.55/0.188 (CSCS reprocessing). Whether the live daily cron on zueub434 has the override enabled cannot be determined from the repository (`doc/OPERATIONS_v22_cutover.md:116` lists the export as a cut-over step).
* CAMS archive: `ALC_CAMS_DIR` (0.4° monthly) with `ALC_CAMS_DIR_FALLBACK` (1°); the Figure-4 script used `D:/CAMS` daily files.
* `ALC_WV_DISABLE`, `ALC_CLOUD_TYPES_EXTRA`, `aerosol_lidar_ratio` (hard-coded 50 in the runner, not in `options.json`).

---

## 5. Suggested text corrections (minimal)

1. 5.1.2 and App. 4: fix the MS sign ("biases B high and C low; C_L high").
2. 5.1.2 eq. for B: write ∫ over the whole profile (0.1–2.4 km in practice, through the cloud), not 0 → z_cb; same on Fig. A3-1.
3. 5.1.2: "C = C₀/C_L with C₀ the constant applied to the input (1e8 for CL31/CL51 L1, 1 for CL61)"; same in the Fig. 4(f) caption, and report C_L in absolute units or state the normalisation.
4. 5.1.2 health filtering: drop the L2 quality-flag sentence (or implement it); say the gates are NaN-tolerant.
5. 5.1.2 aerosol test: "integrated from 100 m to 50 m below the backscatter peak".
6. 5.1.2 aerosol transmission: state its typical magnitude (5–10 %) and either propagate the 20–70 sr bounds into the uncertainty or remove the claim.
7. 5.1.2 aggregation: replace the monthly text with "daily median of the accepted per-profile coefficients (≥ 5 profiles), Kalman-filtered as in 5.1.3 (applied to both methods)"; extend 5.1.3 accordingly.
8. 5.1.2 WV: "the two-way transmission is the emission-spectrum-weighted average of the monochromatic transmissions computed from a HITRAN look-up table (Wiegner and Gasteiger, 2015)"; "day" instead of "night"; mention the nearest 3-h CAMS step and 0.4° point.
9. App. 1 / 5.1.2: state that 910.55/0.188 is operational only where `ALC_WV_SPECTRUM` is set (or make it the code default).
10. App. 4: add the in-cloud extinction α = 10 km⁻¹ (with its origin and sensitivity), the a_G vs r_e distinction, the CL61 row (0.92/0.88/0.81), and that the operational correction is a fixed table (Fig. 4(e) caption).
11. Fig. 4: regenerate from the operational pipeline (n = 248, C = 1.29 with the aerosol correction) or state that the aerosol correction is not applied in the figure.
12. Consider `cbh_maxheight = 2100 m` (or integrate to peak + 300 m) and start the aerosol-correction integral at 100 m / stop it at the CBH; both are one-line changes in `_filters.py`.
