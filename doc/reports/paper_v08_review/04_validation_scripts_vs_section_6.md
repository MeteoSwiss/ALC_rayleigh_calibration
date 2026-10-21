# Adversarial review of Section 6 "Validation" (paper v0.8) against the repository code

Reviewer role: hostile referee + code auditor. Repository `C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration`
(branch `remote-dark`, nothing modified). Paper text: `scratchpad/paper/paper_v08.md` lines 313–385. Figures 10–15 =
`media/image15.png` … `image20.png`; their panel titles are byte-for-byte the titles produced by `run_paper_v08.py`
and the saved `fig_v08_*.png`, so the paper figures ARE the saved outputs.

Reproducibility status of the inputs used for this review (all on disk):

| Item | Where | Status |
|---|---|---|
| Section 6.2/6.3 statistics | `C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/summary_v08.json` (written 2026-10-09 16:03) | present; every quoted number cross-checked |
| Section 6.1 statistics | `.../paper_v08/method_comparison_summary.json` + per-stream `method_test_2026/<key>/<key>_cal.csv` | present; flags re-tabulated |
| Hourly matched arrays | `.../paper_v08/cache_<site>.npz` | present; used for the recomputations below (quarterly medians, day/night, estimator sensitivity, Angstrom sensitivity) |
| Dashboard constants | `.../paper_v08/calib_dashboard/<WMO>_<ident>_<rayleigh|cloud>.csv`, `SOURCE.txt` = "fetched 2026-10-09 12:54 UTC" | present |
| EARLINET SCC files | `A:/EARLINET/ari` (1834 b1064 files), `A:/EARLINET/ino` (2900) | present; metadata inspected |
| Dashboard HTML, bucket NetCDFs | remote | not re-fetched (not reproducible offline; the CSVs are taken as given) |
| Appendix 5/6 inputs (dark npz, `diag_v22_dark_v2`) | referenced by `run_paper_v08.py:108–118` | not audited (out of scope) |

Run logs show the full pipeline ran twice (`run_full.log` crashed in a plotting call at an earlier code version;
`run_pay_sir.log`/`run_amsterdam.log` completed with `PAPER_V08_DONE`). The JSON on disk is from the completed runs.

---

## 1. Processing-chain table (Sections 6.2 and 6.3)

Legend: OK = paper statement is what the code does; PARTIAL = true in substance but the description omits a
material detail; WRONG = description does not match the implementation.

### 1.1 Shared chain (Section 6.3 stations; same functions feed the CHM15k side of 6.2)

| # | Paper statement | Code (file:line) | Verdict | Actual behaviour |
|---|---|---|---|---|
| 1 | "Starting from the Level-1 range-corrected signal" | `intercompare.read_l1` 262–307, `_l1_day` 228–259, `retime_hourly` 618–650 | PARTIAL | `rcs_0` is **hourly-median aggregated at read time** (one `retime_hourly` per daily file, coverage `n·dt ≥ 1800 s` at 643). Every later step (overlap, calibration, WV, conversion, screening) runs on hourly medians, not on native profiles. `cloud_base_height`, `vertical_visibility` and `temp_int` are hourly-medianed the same way (257). |
| 2 | "temperature-dependent overlap correction is applied to the CHM15k (Section 5.2)" before division | `run_paper_validation.channel_beta` 117–122 → `overlap.correct_rcs` 87–107; models `D:/TEMP_MODELS/202606` (39) | OK | Multiplicative `rcs·(1+(a(z)·T+b(z))/100)`, T = hourly-median internal temperature (K→°C at `intercompare` 245), a,b → 0 above the model top (`right=0.0`, 99–100), factor clamped 0.1–10 (106). Applied before the division by C_L (136). A missing model only prints a line (122); no such line in the run logs, so all CHM15k were corrected. Consistent with the constants: the operational Rayleigh fit starts ≥ 1.5 km (`molecular_methods.py:131`), where the correction is a no-op. EARLINET side identical (`earlinet.compare` 154–157). |
| 3 | "divided by the operational lidar constant … daily Kalman best-estimate … as published by the dashboard" | `channel_beta` 132–137; `intercompare.load_calib_series` 198–214 (falls through `<key>_L1.csv` → `<key>_L2.csv` → `<key>.csv`, i.e. the dashboard CSV); `interp_calib` 217–224; `fetch_dashboard_calib.kalman_series` 55–65, `write_stream` 93–112; `run_paper_v08._configure` 166–171 | OK, with caveats | **Which constant**: the Plotly traces "Kalman estimate (v2.2)" (Rayleigh) and "Kalman estimate (O'Connor)" (cloud) of each station page, fetched **2026-10-09 12:54 UTC** (`SOURCE.txt`). β = rcs_0/C_L·1e6 with C_L **linearly interpolated** between daily Kalman values at the hourly times and **clamped** outside the series (`np.interp`). **Kalman type**: causal predict/update (`monitoring/kalman.py` 31–41; the dashboard plots the precomputed `<key>_kalman.csv`, `charts.py` 478–488): the value holds between calibrations and jumps after each new one (CSV evidence: Payerne CL31 cloud 3.871e7 held 3–14 Jul, 6.50e7 on 15 Jul, 7.64e7 on 17 Jul; Payerne CL61 cloud 1.087 held 1–9 Mar on 2 cloud days). **Days without a value**: inside the record 1–3 days per stream are missing and are linearly interpolated; before the first / after the last Kalman day the end value is held (Payerne CL61 cloud ends 2026-10-02 → 3–8 Oct use the 2 Oct value; Amsterdam C starts 2025-01-13 → 11–12 Jan use the 13 Jan value). No rule excludes days far from any calibration. |
| 4 | "corrected for water-vapour absorption with CAMS humidity profiles and the same emission spectra as their calibration" | `channel_beta` 140–147 → `intercompare.apply_wv` 366–394 → `water_vapor.two_way_wv_transmission` 478–554 (Gaussian line, ±3σ); spectra `intercompare.WV_PARAMS` 91–92 = `water_vapor.LASER_SPECTRUM` 35–39 for CL31 (909.7/6.0) and CL51 (910.0/3.4); CL61 forced to (910.55, 0.188) by `run_paper_v08.py:42` (env) and `:170` | OK | Same spectra as the dashboard constants: verified here by comparing the local v2.2 method-test rerun (same env) with the dashboard nightly values of the Payerne CL61 — median ratio 1.0000 on 60/60 cloud days and on 26 common Rayleigh nights. **Omitted detail**: the humidity is the **monthly-mean** CAMS profile for every month that has a 0.4° monthly file (`cams_periods` 338–352; `D:/CAMS_Monthly_04` covers 2025-01…2026-06) and the daily 0.4° file from 2026-07 on (354–363); the calibration constants themselves used the humidity of each calibration day. The 1° fallback (`find_cams_month` 435–441, CLAUDE.md orography trap) was NOT triggered for any month of the three station periods (0.4° monthlies exist for all of them). Days without usable CAMS are NaN-masked for the 910 nm channels only (2026-08-31, 2026-10-09 per logs) → the matched sample of the Vaisala channels is 1–2 days shorter than the CHM15k's. |
| 5 | "converted to 1064 nm by separating the molecular contribution, computed from CAMS temperature and pressure, from the aerosol contribution, to which an Ångström exponent of 1 is applied" | `channel_beta` 154–156 → `intercompare.wavelength_correct_molecular` 495–516 (formula at 515); `_mol_att_cams` 468–492; `_molecular_beta` 398–417 (US-std fallback); `ALPHA = 1.0` at `run_paper_validation.py:69` | PARTIAL | Exactly computed: `β1064 = βm,att(1064) + [β910 − βm,att(910)]·(910.55/1064)^1`, with `βm,att = βm·T²m` from the **period-mean** CAMS T/p (monthly mean or daily), hydrostatic extrapolation below the lowest CAMS level, US-std where CAMS is missing. Only the "residual" `β − βm T²m` is scaled (factor 0.856). **Not converted**: (i) the two-way *particle* transmission, kept at its 910 nm value inside the residual; (ii) the cross term "molecular return attenuated by particles", `−βm T²m (1−T²p)`, which is scaled by 0.856 although at 1064 nm it is ≈0.46 of its 910 nm value. Net: in an aerosol-free gate above a layer of particle optical depth τp(910) the converted β_att is **low by ≈1.5·τp** (−3 % for τp = 0.02, −8 % for 0.05, −12 % for 0.08); the error scales with the molecular fraction of the gate (≈1 % in the boundary layer, largest at the top of the 500–3000 m band). See §4.5. α is a *backscatter* Ångström exponent applied to an *attenuated* residual; it is fixed, no uncertainty is propagated. |
| 6 | "Profiles affected by clouds, fog or precipitation (any cloud base below 20 km or a reported vertical visibility, ±15 min around them) are removed" | `channel_beta` 162–164 → `intercompare.screen` 520–538 | WRONG in the details | `screen` runs on the **hourly** arrays. `has_cloud` (524) = any finite CBH in (0, 20 km) in the hourly `nanmedian` of the L1 CBH layers → one native profile with a cloud base flags the whole hour. The "±15 min" window is `win = 2·round(900/3600)+1 = 1` (529–533) → **no expansion** into adjacent hours. Vertical visibility (525–526): hourly median of the native values *including* the −1/−99 "no fog" sentinels, sentinel removed afterwards → fog is detected only when reported for more than half of the hour. No precipitation variable is used (precipitation is removed only through the instrument's own CBH/VV). Each channel is screened with its **own** cloud detection; a pair therefore uses the union of two instruments' detections. `qf` is all-zero for L1 (`read_l1` 305) → the quality-flag line is inert. Effective rule to print: "any hour in which the instrument reported a cloud base below 20 km, or vertical visibility for most of the hour, is removed". |
| 7 | "averaged to hourly medians (≥ 30 min of data per hour) on a common altitude grid" | `_l1_day` 257 / `retime_hourly` 625–648; `grid_and_stats` 192 (second retime, no-op), 205–209; `build_common_grid` 653–659; `regrid` 662–670 | OK | Hourly bins labelled by their **start**; coverage `n_samples·dt_native ≥ 1800 s`. Common grid: step = coarsest channel (15 m at Payerne/Palaiseau, CL61 4.8 m and CL31 10 m block-averaged within ±7.5 m), from the highest first gate to the lowest last gate. |
| 8 | "No noise filtering and no dark-signal correction" | `intercompare.noise_filter_mode` 67–72 ("none" unless env); `grid_and_stats` 188 passes `snr_idx=(0,)` but `snr_mask` 555–556 returns all-True for < 5 samples per bin (hourly input → 1 sample); `channel_beta` 126–129 (`bdark` only for appendix channels, `run_paper_v08.py` 134–141) | OK | True, but the SNR≥3 mask is inert only because the data are already hourly — a fragile way of "not filtering". Worth one sentence in the code, not in the paper. |
| 9 | "median of the matched profiles, computed per range gate over the hours where all instruments provide data" | `run_paper_v08.fig_station` 395–398 (`common = np.all(isfinite)`), 396 (`n_common_h`), `_median_iqr` 290–298 with `minn = max(30, 2 %·hours)` (393) | OK for panel (a) only | The relative-difference profile (panel b), the bin-scatter panels and **all statistics** are **pairwise** (405–410, 430), not restricted to the common hours: Payerne CHM15k pair 1704 h, CL31 pair 1583 h, common 1623; Amsterdam pairs 3067/3272/3234 h, common 2768. Recomputed on the common hours the medians move ≤ 0.6 % (Payerne CHM15k +30.4 → +31.0; Amsterdam A +4.2 → +4.3, B +20.3 → +20.2, C +4.9 → +5.0), so the conclusions hold, but the N printed in the panel titles is not the sample of panel (a). |
| 10 | "bin-scatter diagrams of the hourly values of each instrument against the reference between 500 and 3000 m" | `panel_binscatter` 344–375 | OK | 2-D histogram of pairs with **x > 0 and y > 0** (346), 66 log bins from 10^−2.3 to 10 (287); binned median/IQR in 22 bins with ≥ 50 pairs (357). The N in the title is `pair_stats["n"]` (434) = pairs with ref > 0 *including negative test values*, so the title N exceeds the number of plotted points (Palaiseau 648 489 vs 447 327 plotted; Payerne CL31 264 314 vs 222 351). |
| 11 | "the median relative difference 100·(X−ref)/ref over 500–3000 m (stations) / 500–5000 m (EARLINET) is the headline metric" | `rel_diff` 248–253; `pair_stats` 256–275 | OK | Population = **all hour×gate pairs pooled** (N ≈ 1.6–8·10⁵), both finite **and ref > 0**; test values ≤ 0 are kept (d ≤ −100 %); pairs with ref ≤ 0 are dropped (0.0–2.7 % of pairs, 2.5–2.7 % at Payerne where the CL61 reference is at the molecular level at 2–3 km). IQR = pooled percentiles. Layer medians (270–274): same rule per layer (`z >= lo & z < hi`). No weighting by hour → long clear episodes dominate; night/day subsets at 402–408 (22–02 UTC / 10–14 UTC). |
| 12 | "the correlation is computed on log-transformed values (log r)" | `pair_stats` 262, 267 | OK | Pearson r of log10 over the pairs with **both** values > 0 (a smaller population than N: Palaiseau 447 327 of 648 489; Payerne CL31 222 351 of 264 314), pooled over all altitudes of the band → the value is dominated by the common vertical gradient (β spans > 1 decade from 500 m to 3–5 km), not by instrument agreement at a given height. |

### 1.2 EARLINET-specific chain (Section 6.2)

| # | Paper statement | Code (file:line) | Verdict | Actual behaviour |
|---|---|---|---|---|
| 13 | "particle backscatter coefficients at 1064 nm were converted forward to β_att by adding the molecular backscatter and applying the two-way particle and molecular transmittance, with the lidar ratio assumed by the SCC retrieval" | `earlinet.read_earlinet` 76–80 (molecular), 108–114 (LR), 121–137 (forward model) | OK with caveats | **S**: `assumed_particle_lidar_ratio[0]` of each file; 55 sr in all 1834 Leipzig files, 50 sr in 2602/2900 Măgurele files, the other 298 fall back to `LR_1064 = 50` (29) — coincidentally the same value, so no inconsistency today, but the fallback is silent. **T/p for βm and αm**: the US Standard Atmosphere 1976 (`load_standard_atmosphere`, 78) evaluated at `grid + station_alt`, NOT CAMS; the SCC itself used the US-std at Leipzig (`atmospheric_molecular_calculation_source = 0` in all files) but ECMWF at Măgurele (2517/2900 files). **Transmission**: trapezoidal OD of `ext = βp·S` plus αm from the ground (134–136), `T²` applied to `(βp + βm)` (137). Below the overlap cut the extinction of the lowest trusted gate is **extended to the ground** (134); zero above the top gate. Interior NaN gaps of βp are bridged linearly by `np.interp` (121). |
| 14 | "EARLINET range bins below the full-overlap height of each reference system were discarded" | `SITES` 36–43 (`ari` 800 m, `ino` 1100 m hard-coded), cut at 125 | PARTIAL | The 800 m at Leipzig coincides with the SCC product's own first gate (822 m a.g.l. w.r.t. the lidar in every file; first grid point 810 m w.r.t. the CHM15k altitude) → the cut is a no-op and the "full-overlap height" is whatever ELDA was configured with. At Măgurele the SCC product starts at ≈ 800 m a.g.l. (p50 802 m) but the code discards everything below **1100 m**; this number has no documented source (operator choice) and the discarded 800–1100 m gates are not used for the OD integral either (they are NaN before 126). |
| 15 | a.g.l. / a.s.l. alignment, station altitude | 117 (`rng = alt_asl − station_alt`), `compare` 168 passes the **CHM15k L1** station altitude | PARTIAL | The lidar's own `station_altitude` variable is ignored: Leipzig lidar 125 m vs CHM15k 138 m (13 m shift of the EARLINET profile), Măgurele lidar 93 m (77/83 m in some files) vs CHM15k L1 70 m (23 m shift). Negligible for the statistics, but **horizontal co-location is not what the text says**: the Leipzig CHM15k (`0-20000-0-10471`, 51.320 N 12.450 E, 138 m) is **3.8 km** from the PollyXT (51.3526 N 12.4345 E, "OCEANET … PollyXT_next_generation (arielle) as mobile facility"); Măgurele 0.1–1.4 km (two lidar positions in the files). |
| 16 | "profiles affected by clouds or fog (±15 min around them) were removed" (CHM15k) | `compare` 167 → `intercompare.screen` on hourly data | WRONG in the details | Same as row 6 (whole-hour exclusion on the CHM15k's own CBH, no ±15 min expansion). EARLINET side: the SCC automatic cloud mask (`cloud_mask_type = 2`); the SCC `cirrus_contamination` flag is not used (Leipzig: 830 of 1834 files flagged "cirrus_detected", Măgurele 720 of 2900). |
| 17 | "For each EARLINET profile (averaging time ≥ 30 min), the CHM15k profiles within ±30 min were median-averaged and interpolated onto the EARLINET range grid" | 99–107 (`time_bounds` midpoint; profiles with `tb[1]−tb[0] < 1800 s` skipped: 296 at Leipzig, 567 at Măgurele), `compare` 175–197 | WRONG in the details | The CHM15k series is **hourly** (`read_l1`). The ±30 min window selects the 1–2 hourly bins whose **start label** lies in [t_e−30 min, t_e+30 min] (181–182), i.e. bins covering [t_e, t_e+90 min]: the CHM15k sample is centred ≈30 min **after** the EARLINET window on average and can lie entirely outside it (bin start = t_e+30 min → zero overlap with a 60-min EARLINET average). "Median-averaged" = median of 1–2 hourly medians. The "≥ 30 min of CHM data" test at 185 is always satisfied (1 bin × 3600 s). Interpolation: linear `np.interp` of the hourly profile onto the 15 m EARLINET grid (195). |
| 18 | "Statistics are computed over 500–5000 m a.g.l." | `compare` 201; `ZBAND_EARLINET` `run_paper_v08.py:60`; `pair_stats` 450 | PARTIAL | Effective band = **810–5000 m** (Leipzig) and **1110–5000 m** (Măgurele) because the reference is NaN below; the "500–1000 m" layer is 810–1000 m at Leipzig and empty (NaN in the JSON) at Măgurele. |
| 19 | Daytime profiles | no time-of-day selection anywhere in `earlinet.py` / `run_paper_v08.fig_earlinet` | PARTIAL (undisclosed) | Daytime EARLINET profiles are included: 247 of the 1151 Leipzig matches and 281 of the 636 Măgurele matches have their midpoint at 08–16 UTC (623 / 226 at 20–04 UTC). The CHM15k constant is nightly (Kalman-interpolated to the day). Day and night medians differ by 5–11 % (§4.3). |

---

## 2. Number cross-check

Source for "found": `summary_v08.json` unless marked *recomputed* (from `cache_<site>.npz` / the dashboard CSVs with
the paper's own definitions) or *figure* (panel title). ✓ = reproduced; ≈ = reproduced after rounding; ✗ = not supported.

### 2.1 Section 6.2 and Figures 11–12

| Paper | Found | Verdict |
|---|---|---|
| Leipzig period 30 May 2025 – 2 May 2026 | `period` "2025-05-30 to 2026-05-02" | ✓ (set by the EARLINET data availability, the configured window is 2025-01-01…2026-10-08, `run_paper_v08.py:157`) |
| 1151 temporally matched profile pairs | `n_profiles` 1151, figure "1151 matched profiles" | ✓ |
| N = 320 254 range bins, log r = 0.92 | 320254, 0.9216 | ✓ |
| median +4.3 % (IQR −9 to +20 %) | +4.28, q25 −8.9, q75 +19.8 | ✓ |
| Leipzig layers +7.7 / +5.2 / +0.7 % (1–2, 2–3, 3–5 km) | +7.75 / +5.17 / +0.71 | ✓ |
| +13 % between 800 and 1000 m | layer "500–1000" = +13.28 (effectively 810–1000 m); 810–1000 m *recomputed* +13.3 | ✓ |
| EARLINET truncated below ~810 m | `lowest_earlinet_m` 810 | ✓ |
| Măgurele period 6 Jan 2025 – 10 Jun 2026 | "2025-01-06 to 2026-06-10" | ✓ |
| 636 matched pairs; N = 159 962; +0.3 %; log r 0.86 | 636; 159962; +0.35; 0.857 | ✓ |
| IQR −20 to +21 % | −19.6, +21.3 | ✓ |
| Măgurele layers +1.1 / +2.0 / −1.9 % | +1.07 / +1.98 / −1.88 | ✓ |
| truncated below ~1110 m | 1110 (first grid point ≥ the hard-coded 1100 m) | ✓ |
| unit A +2.5 % (log r 0.90) | `magurele_a` +2.54, 0.899 | ✓ |
| "median difference stays within about ±5 % between 1 and 5 km at both stations" | Leipzig 1–2 km = **+7.7 %** (quoted in the same sentence) | ✗ self-contradictory; write "within ±8 %" or "within ±5 % above 2 km" |
| CHM15k median within ~0.01 Mm⁻¹sr⁻¹ of EARLINET above 2 km (Leipzig) | *recomputed* max |ΔMedian| 0.006, median 0.002 | ✓ |
| quarterly medians +1 to +6 % (Leipzig) | *recomputed* 2025Q2 +1.9, Q3 +5.7, Q4 +1.2, 2026Q1 +4.8, Q2 +5.4 | ✓ |
| quarterly −5 to +6 % (Măgurele) | *recomputed* +5.5, −4.5, −1.2, +0.5, −1.5, +5.2 | ✓ |
| "Kalman lidar constants of these CHM15k vary seasonally by up to ±25 %" | *recomputed* Kalman C_L over the matched period, relative to its median: Leipzig −13 % / +36 % (max/min = 1.56); Măgurele B −30 % / +40 % (max/min = 2.0); Măgurele A −25 % / +43 % | ✗ understated and asymmetric; "seasonally" and the temperature attribution (Bellini et al.) are asserted, not shown |
| "about 11 months at Leipzig and 17 months at Măgurele" | 11.1 and 17.1 months | ✓ |
| Fig. 11 caption: "(b) … median and IQR of the CHM15k values in bins of the EARLINET value", "logarithmic colour scale", "(c) dotted lines ±10 %" | code 344–375, 332–333 | ✓ (bins need ≥ 50 pairs) |
| Fig. 12 caption "As Figure 10, but for the RALI lidar…" | Figure 10 is the method comparison; the mirrored layout is **Figure 11** | ✗ cross-reference off by one |

### 2.2 Section 6.3 and Figures 13–15

| Paper | Found | Verdict |
|---|---|---|
| Payerne 4 Mar – 8 Oct 2026, 1623 common hours | period "2026-03-04 to 2026-10-08", `n_common_hours` 1623 | ✓ (4 Mar = first CL61 *Rayleigh* Kalman day, code comment `run_paper_v08.py:70`; the CL61 *cloud* record — the reference — starts 2026-02-28, so "the start of the CL61 calibration record" is imprecise) |
| CHM15k +13 / +33 / +40 % (500–1000, 1–2, 2–3 km) | +13.09 / +32.64 / +39.55 | ✓ |
| +30.4 % over 500–3000 m, log r 0.93 | +30.37, 0.927; figure "+30.4 %, log r = 0.93, N = 284350" | ✓ |
| night +28 %, day +33 % | 27.5 / 33.1 | ✓ |
| CHM15k −3 %, CL31 +3 % in 100–500 m | −2.46 / +2.85 | ✓ |
| CL31 median −0.6 %, log r 0.32 | −0.61, 0.315; figure "−0.6 %, log r = 0.32, N = 264314" | ✓ |
| CL31 −24 % (500–1000 m) to about +100 % (2–3 km) | −23.6 / +97.4 | ✓ |
| CL31 exclusion "7–17 July" = "days between the replacement and the first post-replacement cloud calibration" | code excludes [2026-07-07 00:00, 2026-07-18 00:00) (`run_paper_v08.py:79–80`) = 7–17 July inclusive; dashboard CSV: last pre-swap cloud day 3 Jul, **first post-swap cloud calibration 15 Jul** (Kalman 3.87e7 → 6.50e7), next 17 Jul (→ 7.64e7) | ✓ consistent; strictly the first post-swap calibration is 15 Jul, so 15–17 Jul are excluded in addition (conservative, fine — the Kalman was still 20 % below the new level on 17 Jul) |
| "the dark lowers the operational Rayleigh constant of the Payerne CHM15k by about 20 % (the constant obtained after dark subtraction is 25 % higher)" | `chm15k_over_operational.dark` = 1.248 (1/1.248 = 0.80) | ✓ |
| "CL61 and CHM15k agree within 7 % after dark subtraction of both" | appendix pair "CL61 cloud + dark vs CHM15k + dark" = −6.96 % | ✓ |
| dark ~10⁻² Mm⁻¹sr⁻¹ = 10–20 % of the molecular signal | molecular 1064 nm medians 0.065 (3 km), 0.053 (5 km) → 15–19 % | ✓ (order of magnitude) |
| Amsterdam 11 Jan 2025 – 8 Oct 2026, 2768 common hours | "2025-01-11 to 2026-10-08", 2768 | ✓ |
| A +4.2 %, C +4.9 % (log r 0.95, 0.94) | +4.24 (0.951), +4.88 (0.938) | ✓ |
| B +20.3 % | +20.33; figure log r 0.95 (0.945) | ✓ |
| B +29 % (500–1000 m), +12 % (2–3 km) | +29.0 / +12.3 | ✓ |
| night, 3–5 km: A, B, C within ±3 % | +0.6 / −1.3 / −2.7 | ✓ |
| B +30 % day / +15 % night (500–3000 m) | 30.1 / 14.6 | ✓ |
| A and C "within about ±10 % between 300 m and 4 km" | *recomputed* median rel-diff profile 300–4000 m: A −7.2…+10.3, C **−12.1**…+9.8 | ≈ (C reaches −12 % near 300 m) |
| Fig. 14 caption "Panels as in Figure 12" | the two-row layout is that of **Figure 13** | ✗ cross-reference off by one |
| Fig. 14 caption flattening below ~0.07 Mm⁻¹sr⁻¹ = reference noise | consistent with panels c–e (binned medians flatten at ≈0.06–0.07) | ✓ (interpretation plausible; it is the "ref ≤ 0 / ref noisy" selection effect of §4.4) |
| Palaiseau 3 Jan 2025 – 8 Oct 2026, 3890 common hours | "2025-01-03 to 2026-10-08", 3890 | ✓ |
| CL31 within −1 % (+3 % at night) in the lowest 500 m | 100–500 m −0.74; night +2.88 | ✓ |
| −17 % (500–1000), −34 % (1–2 km) | −17.2 / −33.7 | ✓ |
| −28.5 % over 500–3000 m, log r 0.65 | −28.52, 0.650; figure "−28.5 %, log r = 0.65, N = 648489" | ✓ |
| "median approaches zero above ~2 km" | *recomputed* CL31 median 0.024 at 2.5 km, 0.005 at 3 km (molecular ≈ 0.07) | ✓ |
| Fig. 15 caption "Layout as in Figure 10" | the one-row (a, b, c) layout is that of **Figure 11** | ✗ cross-reference off by one |
| Figure titles vs captions (panel letters, bands, N) | all identical to the JSON | ✓ |

### 2.3 Section 6.1 and Figure 10 (from `method_comparison_summary.json` and the `_cal.csv` flag census)

| Paper | Found | Verdict |
|---|---|---|
| 33 instruments, 16 countries, 7 CL31 / 5 CL51 / 10 CL61 / 7 CHM15k / 4 Mini-MPL | 33 keys in `method_test_streams.json`; countries CH, GB, FR, ES, CZ, SE, GR, BE, IS, AT, CY, DE, NL, RO, NO, IT = 16; type counts 7/5/10/7/4 | ✓ (note: `run_method_test.ps1:11` ran a 34th stream, Lauder NZ `0-20008-0-LAU_A`, which failed on every day with flag −4 "No CAMS" and was silently dropped from the figure) |
| 1 Feb – 29 Jun 2026 | ps1 `--start 20260201 --end 20260629`; `N_DAYS = 149` | ✓ |
| CL31: 0–3 valid nights, 0–2 % of days | n_ray 0/1/0/3/1/0/0; 0–2.0 % | ✓ |
| CL31 "mostly because no candidate window passes the validity gates" | flag −2 ("No molecular window passed the validity gates") vs −1 ("Not a clear night", `calibration.py:605–621` = cloud/fog in or below the window): Athens 113/36, Getafe 88/45, Palaiseau 83/64, Payerne 79/69, Temelín 68/**79**, Gävle 48/**58** (+43 nights without a night-time window), Camborne 23/**123** | ≈ true for 4 of 7 units; at Camborne the clear-night screen dominates |
| CL31 Rayleigh/cloud factors 0.3–12 | 0.33 (Palaiseau, 1 night), 7.0 (Getafe, 3), 12.5 (Camborne, 1) | ✓ |
| CL51 Rayleigh 3–44 %, median 10 % | 2.7 / 4.7 / 10.1 / 34.2 / 44.3 → median 10.1 | ✓ (2.7 rounded to 3) |
| CL51 ratio median 1.38; 0.99 Keflavík; 1.2–2.5 others | 1.384; 0.991; 1.20 / 1.38 / 1.65 / 2.49 | ✓ |
| CL51 scatter 20 % vs 10 % | type medians 20.5 (4 units; Praha NaN, n = 4 < 5) vs 10.4 | ✓ |
| cloud 19–70 % at the 12 CL31/CL51, scatter 6–17 % | 19.5–70.5 %; 5.9–16.9 % | ✓ |
| CHM15k cloud 15–54 % of days | 15.4–54.4 | ✓ |
| CHM15k ratios 1.47–1.72, 1.62 overall ("32–42 % below") | 1.472–1.720, median of medians 1.624; 1/1.472 = 0.68, 1/1.720 = 0.58 | ✓ |
| CHM15k scatter 18–41 % (cloud) vs 15–25 % (Rayleigh) | 17.5–40.8; 14.7–24.9 | ✓ |
| Mini-MPL: 0 cloud days; 68–91 % of days fail the factor-20 test; Rayleigh 6–41 %; scatter 5–15 % | n_cld = 0 (all four); flag −22 "peak not sharp above (300 m)": Lille 90, Brest 104, Mobile 102, Aléria 125 days = **60.4 / 69.8 / 68.5 / 83.9 % of calendar days** but **79.6 / 70.7 / 68.5 / 90.6 % of days with data**; Rayleigh 6.0–40.9 %; scatter 5.4–14.5 % | ≈ the "68–91 %" uses *days with data* as denominator whereas every other percentage of the section uses the 149 calendar days (`fig_method_comparison.py:28,80`) — inconsistent within one paragraph |
| CL61 Rayleigh 11–49 %, cloud 3–48 % | 10.7–49.0; 2.7–48.3 | ✓ (rounded) |
| CL61 median ratio 0.94; eight units within ±10 % (0.91–1.10) | 0.944; excluding Payerne/Uccle: 0.910 … **1.102** (Aosta) | ≈ Aosta is marginally outside ±10 % |
| CL61 scatter 4–10 % cloud, 6–17 % Rayleigh (26 % Payerne) | 4.4–9.9 (Lanžhot NaN, n = 4); 5.6–16.9; Payerne 25.8 | ✓ |
| Payerne 0.83, Uccle 0.79; "17–21 % below" | 0.830 (n = 29), 0.786 (n = 17) | ✓ |
| "cloud method … limited by the frequency of suitable liquid clouds, especially at Lanžhot and at the Alpine sites of Aosta and Sion" | Lanžhot: 121 of 133 days flagged **−21 "laser energy too low"** (1 day −22); Aosta: −21 32, −22 51, −24 33; Sion: −22 61, −24 38; Camborne CL61 (13 %): −24 "aerosol below cloud > limit" 94 days | ✗ for Lanžhot (an instrument-energy gate, not cloud frequency); loose for Aosta/Sion (cloud-quality gates −22/−24, not absence of clouds) |
| Fig. 10 caption: pairing convention, 10–90 % whiskers, ±10 % band, "number of pairs at the bottom", 1.48·MAD(ln C_L) | `fig_method_comparison.py` 50–60, 63–67, 75–77, 130–136 | ✓ (1.4826) |
| "0.83 … rises from 0.78 to 0.94 (Appendix 6)" | see §3.4 | explained: period effect |

---

## 3. Section 6.1 — settings, pairing, denominators, and the 0.83 vs 0.78 question

**3.1 Which settings were used.** `run_method_test.ps1` runs `scripts/run_streams_parallel.py` with
`ALC_MOLECULAR_METHOD=eprof_v2.2`, `ALC_L1_ROOT=D:/E-PROFILE_L1_2026`, `ALC_WV_SPECTRUM={"CL61":[910.55,0.188]}`,
`ALC_RAYLEIGH_TYPES_EXTRA=CL31,CL51`, `ALC_CLOUD_TYPES_EXTRA=CHM15k,Mini-MPL` (lines 6–10), `--methods rayleigh,cloud`,
CAMS `D:/CAMS_Monthly_04;D:/CAMS_daily`, date range chunked in 5 (lines 11–13). The type routing is widened by
`run_network_calibration.py:131–137`. Everything else is `options.json` (Rayleigh: window centre searched 2–6 km,
`range_start_m/range_end_m` lines 22–23; v2.2 window start ≥ 1.5 km, `molecular_methods.py:131`; molecular reference =
**US Standard Atmosphere**, `options.json:32,34`; `LRaer: 52`, line 30). "Operational settings" is justified: for the
Payerne CL61 the rerun reproduces the dashboard's nightly constants exactly (median local/dashboard = 1.0000 on 60/60
cloud days and 26 common Rayleigh nights; 3 nights exist only locally, 2 only on the dashboard — CAMS daily/monthly
differences at the margin). This was checked for one stream only.

**3.2 CHM15k / Mini-MPL cloud run.** The cloud code applies each type's own PVC multiple-scattering table
(`_filters.py:38–80`, a_G = 5.5 µm for all types, comment `calibration.py:455–470`), no water-vapour correction for
1064/532 nm (`calibration.py:736–737`), and only **warns** about detector saturation (`calibration.py:718–727`) — it
never rejects, so the paper's "the generic cloud screening … contains no such test" is correct. The "factor of 20 within
300 m" is `attenuation_factor = 20` (`calibration.py:144`) applied ±300 m around the peak (`_filters.py:197, 234–247`),
flags −22 (above) / −23 (below). The paper should say that the CHM15k cloud constants *are* multiple-scattering
corrected with the CHM15k's own η table (a_G = 5.5 µm), so the 32–42 % deficit is net of that correction.

**3.3 Pairing and its asymmetry.** `paired_ratio` (`fig_method_comparison.py:50–60`): each value of the *tested*
method is divided by the median of the *reference* method within ±15 days, **or by the whole-period median when fewer
than 3 reference values fall in the window** (silent fallback, line 59). The Vaisala/photon-counting asymmetry
(lines 74–77) is a design choice — the non-operational method is always the numerator-side sample — so that the box
width in Fig. 10b reflects the scatter of the method under test. Consequences the text should state: (i) box widths
are not comparable across the two groups (CL61 boxes = Rayleigh nightly scatter, CHM15k boxes = cloud daily scatter);
(ii) with the fallback, Lanžhot (4 cloud days), Camborne CL31 (1 night), Palaiseau CL31 (1), Getafe (3) and Praha CL51
(4 Rayleigh nights) are ratios to a whole-period median or single values, not ±15-day pairs.

**3.4 "Fraction of days with a valid calibration".** Denominator = **149 calendar days** (`N_DAYS`, lines 28, 80),
not days with data. Streams with fewer data days: Temelín CL61 104, Lille 113, Keflavík 123, Sion 124, Payerne CL61
126, Lanžhot 133, Aléria 138, Getafe/Praha/Uccle-CL61 139. "Valid" = flag 1 **or 0.5** (line 44); flag 0.5 is a
*partially clear night* — the fit was made on the clear subset of the night after cloudy profiles were removed
(`calibration/rayleigh/calibration.py:1334`, `calibration/io/data_loader.py:957–1003`). The share of 0.5 among
"valid" Rayleigh nights is large: Aléria 50/50, Mobile 51/61, Granada 45/75, Payerne CHM15k 36/54, Aosta 48/73,
Lanžhot 40/67, Wien 33/51. Neither point is stated. (The Kalman series on the dashboard also ingest 0.5 nights —
`fetch_dashboard_calib.nightly_values` 68–90 keeps every finite `lidar_constant` regardless of flag.)

**3.5 Scatter metric.** `scatter()` lines 63–67: 100·1.4826·MAD(ln C) over the valid values, NaN below 5 values —
matches "1.48 times the MAD of ln C_L". Type-level medians therefore skip Praha CL51 Rayleigh (4 nights) and Lanžhot
cloud (4 days).

**3.6 0.83 (Section 6.1) versus 0.78 (Appendix 6).** Recomputed from the dashboard CSVs (`calib_dashboard/
0-20000-0-06610_C_{rayleigh,cloud}.csv`), nightly Rayleigh C_daily over the cloud reference:

| period | nightly / Kalman-cloud (appendix definition, `fig_appendix_dark` 563–564) | nightly / median cloud within ±15 d (Section 6.1 definition) |
|---|---|---|
| 2026-03-04 … 2026-10-08 (appendix) | **0.783** (n = 62) | 0.785 |
| 2026-02-01 … 2026-06-29 (Section 6.1) | 0.825 (n = 28) | 0.821 |
| local method-test rerun, Feb–Jun | — | **0.830** (n = 29) |

The difference is a **period effect** (the ratio drops to ≈0.75 in Jul–Oct), not a pairing (±0.004) or source
(local = dashboard, 1.0000) effect. The paper should say "0.83 over February–June; 0.78 over the March–October period of
Appendix 6, the summer nights lowering the ratio", and may note that a lower summer ratio is what a
temperature-dependent dark (Le et al.) predicts.

---

## 4. Referee-level methodological issues (with code evidence and what to change)

### 4.1 Circularity of the EARLINET free-troposphere agreement
*Evidence.* The CHM15k Rayleigh constant is the ratio of the clear-night signal to the US-std molecular profile in a
window whose centre is searched between 2 and 6 km and that starts ≥ 1.5 km (`options.json:22–23`, `molecular_methods.py:131`,
`options.json:32,34` = US-std). The EARLINET β_att above the aerosol is, by construction of the forward model
(`earlinet.py:76–80, 137`), `βm(US-std)·T²` where the SCC retrieval set βp = 0 ("assumed backscatter-ratio value in
calibration range" = 1 in every file). The SCC normalisation height is 10–12 km at Leipzig (p5 6.1–8.0 km) but
**5.5–5.9 km at Măgurele (p5 3.5–4.0 km)**, i.e. inside the CHM15k window. So at Măgurele the 3–5 km layer compares two
quantities both pinned to the same molecular model on the same clean scenes; at Leipzig the Klett solution propagates
down from 10–12 km, so the 3–5 km βp is a genuine retrieval but the absolute scale of both instruments is still the same
US-std βm. The "+0.7 %" and "−1.9 %" at 3–5 km therefore cannot detect a common error of the molecular reference
(density ±3–5 % seasonally versus the real atmosphere; the Măgurele SCC itself used ECMWF T/p in 2517/2900 files while the
forward model re-adds US-std βm) and are partly tautological. The genuine test is 1–3 km (aerosol-dominated), where the
numbers are +7.7/+5.2 % (Leipzig) and +1.1/+2.0 % (Măgurele).
*Change.* State the SCC normalisation ranges and the CHM15k window; reword "agreement over more than four kilometres …
obtained without any tuning to the reference instrument" — both are normalised to the same molecular model; present
the 1–3 km layer as the validation of the transfer and the 3–5 km layer as a consistency check.

### 4.2 Boundary-layer transmission assumption of the forward model
*Evidence.* `earlinet.py:134` extends the extinction of the lowest trusted gate (810/1110 m) to the ground; this
two-way transmission scales the **entire** EARLINET profile above. For the Leipzig median scene (βp ≈ 0.2 Mm⁻¹sr⁻¹ at 810 m,
S = 55 sr) T² = 0.98 (2 %); for a loaded boundary layer (βp = 1) T² = 0.92 (8 %); at Măgurele the extension runs from
1110 m (3 % median, > 10 % loaded) although the SCC product contains data down to ≈ 800 m that could feed the OD integral.
Since the real lowest layer is usually more loaded than its top, the assumption biases the EARLINET β_att high and the
CHM15k relative difference negative by a few per cent — of the same size as the headline numbers.
*Change.* Quantify (e.g. use the CHM15k's own BL backscatter × S for the OD below the cut, as a sensitivity), use the
800–1100 m SCC gates at Măgurele for the OD, document the origin of the 800/1100 m cut-offs.

### 4.3 Day and night are mixed, and they differ
*Evidence (recomputed from the caches, midpoint 08–16 UTC vs 20–04 UTC, 500–5000 m).* Leipzig: day **+8.0 %** (247
profiles) vs night **+3.0 %** (623); Măgurele B: **+5.9 %** (281) vs **−4.9 %** (226); Măgurele A: +8.1 % vs −3.2 %.
Payerne CHM15k: +33 vs +28 (quoted); Amsterdam B: +30 vs +15 (quoted). The CHM15k reads systematically higher by day
relative to every reference, by 5–15 %. Nothing in 6.2 mentions that daytime SCC 1064 nm products (low SNR) are
included, nor the day/night split. The "+0.3 %" at Măgurele is the average of −4.9 % and +5.9 %.
*Change.* Add the day/night split to 6.2 (one sentence + numbers), discuss the daytime CHM15k excess as a common
feature (solar background/afterpulse?), and state that the nightly Rayleigh constant is applied to daytime data.

### 4.4 The headline metric: pooled hour×gate population, asymmetric exclusion, estimator sensitivity
*Evidence.* `rel_diff` keeps test values ≤ 0 and drops ref ≤ 0 (`run_paper_v08.py:252`); `log r` uses both > 0 (267);
`ratio_of_medians` is also computed (268) but never reported. Recomputed on the same pairs:

| pair (band) | fraction of test values ≤ 0 | ref ≤ 0 excluded | median rel. diff. (paper) | median log-ratio (both > 0) | ratio of pooled medians |
|---|---|---|---|---|---|
| Payerne CHM15k vs CL61 | 0.4 % | 2.7 % | +30.4 | +30.5 | +29.0 |
| Payerne CL31 vs CL61 | **15.9 %** | 2.5 % | −0.6 | +14.0 | +32.1 |
| Amsterdam A / B / C vs D | 0.0 % | 0.2–0.3 % | +4.2 / +20.3 / +4.9 | +4.2 / +20.3 / +4.9 | +2.9 / +16.7 / +4.2 |
| Palaiseau CL31 vs CHM15k | **31.0 %** | 0.2 % | −28.5 | −10.6 | +3.7 |
| Leipzig | 0.1 % | 0.0 % | +4.3 | +4.3 | +5.8 |
| Măgurele B / A | 2.0 / 0.7 % | 0.0 % | +0.3 / +2.5 | +1.1 / +2.8 | **+12.3 / +9.4** |

Per layer, Palaiseau 1–2 km is −33.7 % (paper) but −14.0 % on positive pairs; 2–3 km is −77 % vs +37 %: in the
noise-dominated regime no estimator is meaningful, and the quoted layer values are estimator-dependent by tens of per
cent. The direction of the ref > 0 conditioning is subtle: dropping the negative *reference* noise realisations
retains reference values biased high, so the pairwise median is biased **negative** where the reference is noisy
(Amsterdam 2–3 km, Payerne CL61 at the molecular level), not positive — the Fig. 14 caption describes exactly this
effect. Keeping negative *test* values, by contrast, is the unbiased choice when the reference is clean (Palaiseau).
Finally, the pooled log r over 500–3000/5000 m mostly measures the shared vertical gradient.
*Change.* State the population (pooled hour×gate pairs, N), the fraction of non-positive values, that negative test values
are kept, that log r is on positive pairs only (give that N) and pooled over altitude; add the ratio of medians (or the
median of differences normalised by the median reference) as a second estimator in the tables; quote Palaiseau/Payerne
CL31 layer numbers above ~1.5 km as "noise-dominated" rather than as percentages.

### 4.5 The 910 → 1064 nm conversion: Ångström sensitivity and an uncorrected cross term
*Ångström (recomputed on the caches with the US-std molecular profile, factor (910.55/1064)^(α−1) on the residual).*

| α | Payerne CHM15k vs CL61, 500–3000 m | Payerne 100–500 m (CHM15k) | Palaiseau CL31 vs CHM15k, 500–3000 m | Palaiseau 100–500 m |
|---|---|---|---|---|
| 0.5 | +23.9 % | −8.4 % | −25.0 % | +6.0 % |
| 1.0 (paper) | +30.4 % | −2.5 % | −28.5 % | −0.7 % |
| 1.5 | +36.8 % | +3.8 % | −31.8 % | −6.9 % |

The aerosol fraction of the converted channel's median profile is 0.8 below 1 km and ≈ 0 at 2.5–3 km (Payerne CL61:
0.456 vs βm 0.087 at 300 m; 0.074 vs 0.069 at 2.5 km), so a ±0.5 uncertainty of the *backscatter* Ångström exponent
moves the boundary-layer agreement by **±6–7 %** and the 500–3000 m headline by ±3–6 %. The statements "CHM15k and CL31
lie within −3 % and +3 % of the CL61 [100–500 m] despite their different wavelengths" and "the CL31 agrees with the
CHM15k to within −1 %" are therefore agreement *within the conversion uncertainty*, not evidence of per-cent closure;
Section 6.3 says so for Palaiseau only.
*Cross term (from `intercompare.py:515`).* For an aerosol-free gate above a layer of particle OD τp(910),
Output − True ≈ −0.9·βm(1064)T²m·(1−T²p(1064)), i.e. a **relative error ≈ −1.5·τp(910)**: −3 % (τp = 0.02), −8 % (0.05),
−12 % (0.08). It vanishes where the aerosol fraction is high and is largest at the top of the band; at Payerne it
makes the CL61 reference several per cent too low at 2–3 km (part of the "+40 %" and of the +30.4 % headline; comparable
to the −7 % residual quoted after dark correction), at Palaiseau it pushes the CL31 further below the CHM15k.
*Change.* Give an uncertainty budget for the conversion (α range from the literature for 910–1064 nm backscatter;
include the τp term, or correct it with the CAMS/AERONET AOD or the residual itself × an assumed S as the EARLINET
forward model does); reword the "few per cent" boundary-layer claims accordingly. Also state that βm for the conversion
comes from CAMS T/p while the Rayleigh constants are normalised to the US-std atmosphere (`options.json:32,34`), and
that the humidity for the WV correction is the monthly mean before 2026-07.

### 4.6 The Payerne reference is calibrated by a different method and starts in a Kalman transient
*Evidence.* The reference CL61 carries the liquid-cloud method's systematics (η table a_G = 5.5 µm, assumed droplet
lidar ratio, the 100–2400 m integration, the WV correction inside the cloud integral): any common-mode error of that
method shifts every Payerne number; no budget is given. The CL61 cloud Kalman is based on **2 cloud days** until
10 March (CSV: 1.159 and 1.055 on 28 Feb/1 Mar, held at 1.087 through 9 Mar, then 1.002 → 0.955 → 1.045 → 1.094 within
five days), i.e. a ±8 % start-up transient over the first ~2 weeks of a 7-month record. The paper's reason for 4 March
("start of the CL61 calibration record") is in fact the first *Rayleigh* Kalman day (`run_paper_v08.py:70`).
*Change.* Quote the cloud-method systematic budget for the reference; either start mid-March or state the transient's
size; correct the sentence on 4 March.

### 4.7 CL31 optical-module swap: the exclusion is consistent, but the Kalman lags for weeks
*Evidence.* The exclusion [7, 17 July] matches the CSV (first post-swap calibration 15 July). But the causal Kalman
reaches only 94 % of the post-swap nightly median (8.12e7 vs pre-swap 3.28e7, step ×2.47) on 18–24 July and is < 95 %
on 10 of the 81 days up to 8 October; the CL31 β_att is correspondingly 5–6 % too high on those days. Small compared
with the CL31's range-dependent artefacts, but the text claims the excluded days are the only ones affected.
*Change.* One sentence: "the Kalman estimate converges to the post-replacement level within ~3 weeks (≤ 6 % low on
18–24 July)", or extend the exclusion to the convergence date; cite the memory finding that an operational filter
should be restarted at a documented optical-module swap.

### 4.8 Cloud screening: selection effects and own-detector screening
*Evidence.* Rows 6/16 of §1: whole-hour exclusion on each instrument's own CBH, no ±15 min expansion, fog only if
reported > 30 min of the hour, no precipitation variable, SCC cirrus flag unused. The matched sample is thus
conditioned on clear hours as seen by each instrument (the CL31 misses high clouds that the CL61 detects; a CHM15k-vs-CL61
pair uses the union). Clear-sky conditioning selects dry, high-pressure, often hazy-BL situations — the same regime
the Rayleigh calibration uses — so the validation never exercises humid or post-frontal conditions where the WV
correction matters most, and the 500–3000 m medians are dominated by long clear episodes (no per-day weighting).
*Change.* Describe the screening as implemented; add the number of hours removed per channel; consider a per-day
weighted median as a robustness check.

### 4.9 "Operational constants as published": fetch date, reprocessed history, distributed L2
*Evidence.* Constants fetched 2026-10-09 (`SOURCE.txt`) for L1 up to 2026-10-08; the dashboard series are a
*reprocessed* history under the current configuration (e.g. the CL61 spectrum switch of 2026-08 changed earlier CL61
nights retroactively; the local rerun matches it night-for-night), not what was displayed on each date. Project context
(CLAUDE.md) indicates that the distributed E-PROFILE Level-2 files at the time of writing still carried the v1.0 Rayleigh
constants (CHM15k, CL61) and a default constant for the CL31; a reader cannot reproduce Section 6 from the L2 product.
*Change.* State "constants as published on 9 October 2026 (reprocessed with the current configuration)", that β_att was
rebuilt from L1 `rcs_0` with those constants, and whether/when the distributed L2 product carries them.

### 4.10 Kalman is causal and interpolated/held
*Evidence.* §1 row 3. A day before a jump carries the old constant; the last/first days of a record use a held value;
no maximum distance to the nearest calibration is enforced (Payerne CL31 cloud: 11 days without calibration in July;
Lanžhot-type gaps would be worse).
*Change.* Say "causal daily Kalman estimate, linearly interpolated to the hour, held at the record ends"; report the
median distance to the nearest calibration per stream.

### 4.11 EARLINET temporal matching is shifted by the hourly labels
*Evidence.* §1 row 17. The CHM15k sample is centred ~30 min after the EARLINET midpoint on average and can lie outside
the EARLINET averaging window; "profiles within ±30 min, median-averaged" is not implemented as described.
*Change.* Either fix the selection to bin *centres* (`t_b + 30 min ∈ [t_e−30, t_e+30]`, or re-read native profiles for
the EARLINET path) and re-run, or describe the matching as implemented ("the hourly CHM15k median whose hour starts
within ±30 min of the EARLINET mid-time"). Expected impact: scatter, not bias, but it must be stated.

### 4.12 "Co-located" at Leipzig
*Evidence.* 3.8 km between the E-PROFILE CHM15k (DWD site 10471, 51.320 N 12.450 E) and the PollyXT (TROPOS, 51.3526 N
12.4345 E, "mobile facility" per the SCC files); different station altitudes (138 vs 125 m) not reconciled in the code.
*Change.* Write "3.8 km apart" and treat BL heterogeneity as part of the Leipzig 1–2 km excess (+7.7 %, +13 % at 810–1000 m).

### 4.13 The "±25 % seasonal" Kalman statement and the Bellini attribution
*Evidence.* §2.1: the Kalman C_L spans −13/+36 % (Leipzig) and −30/+40 % (Măgurele) around its median; seasonality and
its temperature origin are not shown, and a factor-2 change is also what a drifting detector, a window-contamination
episode or a clear-night-selection bias would produce.
*Change.* Show the C_L time series with internal temperature (or cite the dashboard figure), quote the actual range,
and keep the conclusion to what the quarterly medians support: the smoothed constant tracks whatever the instrument did,
within ±5 %.

### 4.14 Minor / editorial
* Three caption cross-references are off by one (Fig. 12 → "As Figure 11", Fig. 14 → "Panels as in Figure 13",
  Fig. 15 → "Layout as in Figure 11"): a figure was inserted after the captions were written.
* Panel-title N (pairs with ref > 0) ≠ plotted points (both > 0); "N matched range bins" should say which.
* "Statistics over 500–5000 m" → "from the reference's lowest usable gate (810/1110 m) to 5000 m".
* Fig. 10 x-labels are 16-character truncations ("Athens Observato", "Keflavikurflugvo", "Amsterdam Ap Sch", "Sion
  Epfl Alpole", "Mini Mpl Mobile"; `fig_method_comparison.py:78`).
* Fig. 10b draws IQR boxes for single values (Camborne, Palaiseau CL31).
* Section 6 intro: "Except in Section 6.1 … only the operational calibration is used" — 6.1 itself uses locally rerun
  nightly values (no Kalman), not dashboard values; say so.
* Section 6.2 first paragraph: a verb is missing after the Leipzig/Măgurele sentences ("… at this urban continental
  site" refers to a site not yet named).

---

## 5. What must change before submission (priority order)

1. **Describe the processing as implemented** (hourly medians before every step; whole-hour cloud exclusion on the
   instrument's own CBH, no ±15 min expansion; EARLINET matching on hourly bins; effective EARLINET band 810/1110–5000 m;
   negative test values kept, ref ≤ 0 dropped; log r on positive pairs; pairwise, not common-hour, statistics) — or
   change the code to match the text (items 4.11, 4.8) and re-run.
2. **Add the day/night split to 6.2** and discuss the 5–11 % daytime excess of the CHM15k (4.3).
3. **Give an uncertainty budget for the 910 → 1064 nm conversion** (α ±0.5 → ±6–7 % in the BL; τp cross term ≈ −1.5·τp
   above the BL) and reword the "few per cent" boundary-layer agreements at Payerne and Palaiseau (4.5).
4. **Qualify the EARLINET free-troposphere agreement** (same molecular model on both sides; Măgurele normalisation
   inside the CHM15k window; BL transmission extension; 3.8 km separation at Leipzig) (4.1, 4.2, 4.12).
5. **Report a second estimator** (ratio of medians or median log-ratio) and the fraction of non-positive values; label
   CL31 layers above ~1.5 km as noise-dominated (4.4).
6. **Fix the factual slips**: Leipzig 1–2 km +7.7 % is not "within ±5 %"; Lanžhot's cloud failures are laser-energy
   rejections, not cloud frequency; the Mini-MPL 68–91 % uses a different denominator; Aosta 1.102 is outside ±10 %;
   "0.83 vs 0.78" is a period effect (Feb–Jun vs Mar–Oct); 4 March is the first *Rayleigh* Kalman day; the Kalman range
   at the EARLINET sites is −30/+40 %, not ±25 %; three caption cross-references.
7. **State the calibration provenance precisely**: causal Kalman, interpolated/held, dashboard values fetched 2026-10-09
   as a reprocessed history, flag-0.5 (partially clear) nights included, 149-calendar-day denominators in 6.1, the
   Payerne CL61 cloud start-up transient and the CL31 post-swap Kalman lag (4.6, 4.7, 4.9, 4.10, 3.4).
