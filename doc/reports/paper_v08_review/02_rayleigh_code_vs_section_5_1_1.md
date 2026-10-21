# Adversarial review — paper v0.8 Section 5.1.1 "Rayleigh" vs. the operational Rayleigh code

Repository: `C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration`, branch `remote-dark` (HEAD 18876dd; the
working tree carries uncommitted edits to `options.json` (only `cams_folder`), `scripts/run_network_calibration.py`
(research env override for instrument routing) and `calibration/incremental.py` (Windows file-lock retry) — none
of them change the Rayleigh algorithm). Paper text: `scratchpad/paper/paper_v08.md` lines 151–187 (Section 5.1.1
incl. the Figure 3 caption on line 157). Nothing in the repository was modified for this review.

---

## 0. What is actually "operational" — the configuration the review is made against

The daily cron (`ops/run_daily.sh` → `ops/ops_daily.py` → `scripts/run_network_calibration.py --sens --omb`) builds the
Rayleigh options like this (`scripts/run_network_calibration.py:313-337`):

```python
o = CalibrationOptions.from_json(REPO / "options.json")      # :313  <- options.json is the operational source
o.folder_root = L1_ROOT; o.data_level = DataLevel.L1; o.cams_folder = CAMS
if MOLECULAR_METHOD: o.molecular_method = MOLECULAR_METHOD   # :324  ALC_MOLECULAR_METHOD env (optional)
if os.environ.get("ALC_DARK_PROFILE"): o.dark_profile_file = ...   # :329  unset operationally
if os.environ.get("ALC_WV_DISABLE") == "1": o.apply_wv_correction = False  # :336  unset operationally
```

Instrument routing: `RAYLEIGH_TYPES = {"CL61", "CHM15k", "Mini-MPL"}` (`:131`), `CLOUD_TYPES = {"CL31","CL51","CL61"}` (`:132`).

Operational values, from `options.json` (committed) + `calibration/config.py` defaults for keys options.json does not set:

| parameter | value | source |
|---|---|---|
| `molecular_method` | **`eprof_v2`** | `options.json:29` (no `ALC_MOLECULAR_METHOD` in the committed `ops/config.sh`) |
| `molecular_source` / `use_std_atm` | **`standard`** (US Std 1976) / 1 | `options.json:32-34` |
| `use_sza_night` / `sza_night_threshold` | 1 / **100°** | `options.json:10-11` |
| `min_time_range` | 3 h | `options.json:14` |
| `time_aggregation` / `screen_profile_outliers` / `profile_outlier_nmad` | mean / 1 / 4.0 | `options.json:16-18` |
| `z_low_cloud` | **6000 m** (dataclass default would be 4000) | `options.json:20` |
| `max_ratio_cloudy` | 0.5 — **defined but never read anywhere** | `config.py:224`; grep finds no consumer |
| `range_start_m` / `range_end_m` (window *centres*) | 2000 / 6000 m | `options.json:22-23` |
| `min_window_start_m / min_window_r2 / max_window_rel_error` | 2000 / 0.5 / 50 — **ignored by eprof_v2** (only v1.2 reads them, `rayleigh_fit.py:209-222`) | `options.json:25-27` |
| eprof_v2 gates (registered defaults, "C8") | start ≥ 1500 m, R² ≥ 0.40, residual ≤ 16 %, scattering ratio ≤ 1.15, ratio_std ≤ 0.40, temporal CV ≤ 0.8, rel_error ≤ 15 %, weights 0.25/0.20/0.10/0.10/0.35/0.20, cell flag 4 MAD + 25 % | `molecular_methods.py:103-107` |
| `LRaer` / `threshold_quality` | 52 sr / 15 % | `options.json:30-31` |
| `apply_wv_correction` | 1 | `options.json:38` |
| CL61 laser spectrum for WV | **(910.74 nm, FWHM 1.0 nm)** unless `ALC_WV_SPECTRUM` is exported | `water_vapor.py:35-41`, override `:49-79` |
| `dark_profile_file` / `subtract_background` / `sign_error_v10` | "" (OFF) / False / False | `config.py:240,230,300`; `ops/config.sh:37-39` says "Deliberately NOT set: ALC_DARK_PROFILE" |
| `half_length_options_m` / `fit_range_increment_bins` | 250,490,…,1930 m / 8 bins | `config.py:274,277` |
| LR deltas / window shifts / time subsets | (−20,−10,0,+10,+20) sr / (−200,−100,0,+100,+200) m / 4 subsets × 70 %, seed 42 | `calibration.py:67-70, 912` |
| `aerosol_scattering_threshold` (flag −9) | 2.0 | `config.py:209` |

**The v2.2 question.** `ops/config.sh:19-40` contains a *commented-out* "v2.2 CUTOVER BLOCK" (`ALC_FULLCAL_DIR=…v2.2`,
`ALC_MOLECULAR_METHOD="eprof_v2.2"`, `ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'`). Commit d8788a5 (2026-08-17) says
"deploying the code must not flip anything; the flip is uncommenting the block", and `doc/OPERATIONS_v22_cutover.md`
states that at the time of writing "le cron du serveur tourne encore `eprof_v2`". `CLAUDE.md` asserts (2026-08-31)
that the NetCDFs E-PROFILE publishes carry `calibration_version` 220 (= eprof_v2.2, `output.py:185-190`) and
reproduce the 910.55/0.188 reruns night-for-night. **Whether the server's `ops/config.sh` has the block uncommented
is server-side state that cannot be determined from the repository.** Consequences for the paper:

* The paper's **"version 2.0"** + **"910.55 nm / 0.188 nm"** is not a configuration that exists in the code: with
  the committed configuration (eprof_v2) the CL61 is WV-corrected with **910.74 nm / FWHM 1.0 nm**; the 910.55/0.188
  spectrum is only ever set *together with* `eprof_v2.2` by the cutover block. Either the paper describes v2.2 (and
  must then also describe the v2.2 noise-aware fallback tier, see §2 B6) or it describes v2.0 with the wrong
  spectrum. Appendix 1 ("the manufacturer spectrum … is therefore used operationally") agrees with v2.2.
* `doc/reports/01_rayleigh_calibration.md:123,1144` still says the deployed network runs E-PROF v1.0 — stale.
* Branch `main` (which the memory notes call "operational code only") has `options.json` → `eprof_v1.2` and
  `/data/pay/...` paths: it is frozen and far behind; it is **not** what the paper describes.

Everything below is checked against **eprof_v2 as committed**, with v2.2 differences called out explicitly.

---

## 1. Statement-by-statement table (Section 5.1.1 incl. Figure 3 caption)

Verdicts: **MATCHES** / **DIFFERS** / **NOT IN CODE** (claimed by the paper, absent or inactive in code) /
**NOT IN PAPER** (present in code, undocumented). Line numbers refer to the files named.

| # | Paper statement | Code location | Verdict | Actual behaviour in code |
|---|---|---|---|---|
| 1 | Rayleigh calibration for CHM15k, CL61 and Mini-MPL; CL31/CL51 excluded | `run_network_calibration.py:131-132`; `config.py:70-79` | MATCHES | Runner routes exactly {CL61, CHM15k, Mini-MPL} to Rayleigh. `supports_calibration` would also allow CHM8k/MPL (none in census). `calibrate_rayleigh` still *runs* on CL31/51 if called directly (warning only, `calibration.py:510-514`). |
| 2 | "operational E-Profile Rayleigh calibration version 2.0 … common Python package" | `options.json:29`, `output.py:185-190`, `ops/config.sh:32-35` | DIFFERS (version label) | Committed config = `eprof_v2` (code 200). Published files reportedly carry 220 (eprof_v2.2). See §0. |
| 3 | Night window = real darkness, SZA, "Sun more than about 10° below the horizon", latitude/season adaptive | `data_loader.py:898-920`; `options.json:11` | MATCHES | `keep = in_night & (sza > 100°)`, i.e. Sun ≥ 10° below the horizon (nautical). SZA from a simple declination `−23.44°·cos(2π(doy+10)/365.25)` and hour angle from solar time `UTC + lon/15` (no equation of time, no refraction: boundary error ≲ 15 min). Bounded to [12:00 solar D−1, 12:00 solar D]. Falls back to a 20–04 solar-clock window if lat/lon are non-finite. At high latitudes in summer the window is empty → flag 0 "No profiles in nighttime window" (`calibration.py:572-579`). |
| 4 | "A minimum of 3 hours of valid data is required" | `data_loader.py:997-1003, 1016, 1027, 1042, 1096` | **DIFFERS** | The 3-h test (`min_profiles = 3·60·profiles_per_min`) is applied **only on nights that contain cloud/fog**. A completely clear night returns at `:1000-1003` *before* any duration check: the only minima then are ≥ 10 binned profiles for the MAD screen (`calibration.py:797`, else the screen is skipped, not the night) and ≥ 5 for cell flagging. A clear 1-h night (e.g. short high-latitude summer nights, data gaps) is calibrated. Also, the test counts *profiles*, not hours of valid data, using `profiles_per_min` from the first two time stamps only (`:1010-1014`). |
| 5 | Profiles "converted to geometric height using the station altitude" | `data_loader.py:377`; `average_ceilometer_data:828`; `calibration.py:658, 707` | MATCHES | `altitude_grid = range_alc + station_altitude` (m ASL) is what the US-Std atmosphere is interpolated onto; `range_alc` (AGL) is used for 1/r² and all optical-depth integrals. No double counting (see §2 check A). Tilt is ignored — L1 files carry `tilt_angle` (Payerne CHM15k 3°, CL61 3.8°, Hohenpeissenberg 0°) which nothing reads: ≤ 0.2 % at 4°. |
| 6 | "the signal is first binned onto the common grid (30 m, 5 min). Cloud screening is then performed" | `calibration.py:554-558, 605-606, 633-647`; `data_loader.py:744-833` | DIFFERS (order) | Order is deliberately **cloud screen at native resolution first, then bin** (comment `:554-558`), then the MAD screen on the binned profiles. Binning = block means of consecutive *kept* profiles (`round(300 s / native dt)` and `round(30 m / native dr)` native samples per block, `:792, :825`): CHM15k → 29.96 m × 300 s; **CL61 (4.8 m) → factor 6 = 28.8 m**, not 30 m; Mini-MPL 15 m → 30 m, 30 m → no-op. Blocks are index-based, not clock-aligned, and may straddle gaps left by the screen. |
| 7 | Cloud screening with CBH: contaminated profiles "(and ±15 min around them) are removed" | `data_loader.py:1019-1046` | MATCHES | Profiles with `cbh[:,0] <= z_low_cloud` (or Vaisala fog = finite `vertical_visibility`, CL31/51/61 only, `:987-992`) are dilated by `maximum_filter1d(size=2w−1)`, w = round(15·profiles_per_min) → ±(15 min − 1 profile). |
| 8 | "the signal above the cloud base minus 500 m is masked" | `data_loader.py:1081-1094` | MATCHES | `rcs[i, range >= cbh−500] = NaN` for every remaining profile with a cloud; profile counted "usable" only if `cbh − 500 > range_end_m` (6000 m). |
| 9 | "with a maximum cloudy fraction of 50 % allowed below 4000 m" | `options.json:20`; `config.py:223-224`; `data_loader.py:1016-1028, 1040-1043, 1096-1097` | **NOT IN CODE** (both numbers) | (a) Operational altitude threshold is **6000 m** (`z_low_cloud: 6000`; the 4000 m is only the dataclass default). (b) There is **no cloudy-fraction rule**: `max_ratio_cloudy` exists in `config.py` but is never read. The actual rule is "≥ 3 h-worth of profiles must survive" the low-cloud removal, the ±15-min dilation and the usability test — i.e. the tolerated cloudy fraction is 1 − 3 h/night length (≈ 57 % on a 7-h summer night, ≈ 79 % on a 14-h winter night). |
| 10 | "residual outlier profiles … rejected with a MAD test before the retained profiles are averaged" | `calibration.py:796-810, 825-826` | MATCHES (details undocumented) | On the *binned* profiles: metric = median of `rcs/r²` over 2–6 km per profile; reject `|metric − median| > 4 × 1.4826 × MAD`; applied only if ≥ 10 profiles and if ≥ max(3, 50 %) survive. Collapse = `np.nanmean` (no clipping of negative values — L1 `rcs_0` is ~40 % negative at range and is kept, so no positive noise bias). |
| 11 | Figure 3: Hohenpeissenberg CHM15k, 974 m a.s.l. | L1 file `L1_0-20000-0-10962_0*.nc` | DIFFERS (trivial) | `station_altitude = 977.0 m` in the L1 metadata (the value the code uses). |
| 12 | Figure 3: molecular profile "computed from co-located temperature and pressure data" | `calibration.py:668-711`; `options.json:34` | **DIFFERS** | Operationally `molecular_source = "standard"` → `load_standard_atmosphere(US Std 1976)` (`atmosphere.py:177-221`). The CAMS T/p path (`load_cams_atmosphere`, `:292-363`) exists but is not enabled. The caption contradicts the paper's own text on line 159. The diagnostic plot labels the curve "molecular" (`plotting.py:889, 1144`). |
| 13 | Figure 3: "grey, profiles excluded as signal outliers in the absence of low cloud"; "343 of 761 profiles retained" | `calibration.py:1275-1288`; `plotting.py:978-1000, 1043-1044` | DIFFERS | In the compact diagnostic (the operational PNG) the hatching is drawn over the **native** profiles using the *cloud-screen* keep mask: red = `cbh < z_low_cloud`, grey = "unused but not low-cloud-flagged" = the **±15-min contamination margin and fog profiles**, *not* the MAD outliers (the MAD screen acts on 5-min bins and is not displayed). The native display is strided to ≤ 1000 columns (`_st = ceil(n/1000)`, `:1277`): a 9.5-h CHM15k night at 15 s has ≈ 2300 profiles, so "761" (≈ 2300/3) looks like the strided display count, not the true number. Check the source of these numbers. |
| 14 | Molecular reference "derived from the US Standard Atmosphere 1976", T/p "interpolated to the lidar height grid" | `atmosphere.py:177-221`; `calibration/data/standard_atmosphere_US_1976_50km.csv` | MATCHES | Linear interpolation of T (K) and p (Pa) from a 7.5-m-step table onto `altitude_grid` (ASL). |
| 15 | Bucholtz (1995) with "wavelength-dependent refractive index and depolarization corrections" | `atmosphere.py:29-34, 55-97, 128-148` | DIFFERS (depolarisation) | Refractive index: Edlén/Peck–Reeder dispersion at λ — wavelength-dependent ✓. **Depolarisation `DEPOLARIZATION_RATIO = 0.0301` is a fixed constant** (≈ Bucholtz's ~350 nm value; at 910–1064 nm Bucholtz gives ≈ 0.0279): King factor (6+3ρ)/(6−7ρ) and phase function (`:93-96`) are *not* wavelength-dependent. Effect: β_mol +0.26 % → C_L −0.26 %. Also the molecular lidar ratio used for α_mol is 8π/3 (`:34, :151`) while the depolarisation-corrected phase function implies 8π/3·(1+2γ)/(1+γ) = 8.50 sr → α_mol 1.5 % low → < 0.1 % on C_L. Wavelengths: CL61 reference at **910.0 nm** (`config.py:63`; +0.24 % β_mol vs 910.55), CHM15k 1064.0 (Appendix 1 measures 1064.47: +0.18 %). All negligible but the paper's wording overstates the code. |
| 16 | "expected attenuated molecular signal … two-way molecular transmission and geometric range correction" | `atmosphere.py:153-166` | MATCHES | `p_mol = β_mol · exp(−2∫α_mol dz) / r²`; left-Riemann cumsum on the binned grid from the first gate. |
| 17 | WV correction for 910 nm from CAMS T/q collocated to the station | `calibration.py:732-783`; `water_vapor.py:358-412` | MATCHES (window undocumented) | Nearest CAMS grid point (0.4° Europe box or a regional box, `cams.py:95-212`), L137 model levels integrated hydrostatically per time step and averaged on height. **The averaging window is fixed 20:00 UTC D−1 → 04:00 UTC D** (`calibration.py:664-666`, from `hour_min/hour_max`), *not* the SZA night and not solar time — wrong by ~8–12 h for the out-of-domain stations (Canada, New Zealand). |
| 18 | HITRAN LUT convolved with a Gaussian; CL61 at 910.55 nm / FWHM 0.188 nm | `water_vapor.py:35-41, 69-79, 522-556`; `calibration/data/abs_cross_wv_910nm.nc` | **DIFFERS unless `ALC_WV_SPECTRUM` is exported** | Code default `LASER_SPECTRUM["CL61"] = (910.74, 1.0)`; 910.55/0.188 only through the env override that ships commented-out with the v2.2 block (§0). Gaussian truncated at ±3σ (`:548`), weights normalised over the full LUT grid; T² (not T) is spectrally averaged (`:554`) ✓. LUT: 903–918 nm at 0.0083 nm (58 points within ±3σ of the 0.188-nm line — adequate), 1000 height levels. Cross-section height index uses **AGL range** (`:505-508`) while the humidity profile is interpolated on **ASL** altitude (`:512`) — minor inconsistency for high stations. |
| 19 | "This transmission is divided out of the range-corrected signal before the molecular fit" | `calibration.py:782` | MATCHES | `data.rcs = data.rcs / t2_wv[None, :]` after binning, before the MAD screen and the fit. |
| 20 | 910 nm night without collocated CAMS is not calibrated | `calibration.py:740-778` | MATCHES | Flag **−4** (no file / unusable profile / grid mismatch) or **−10** (nearest CAMS point > max(1°, 1.5 cells) away = outside the regional domain). |
| 21 | "scanning candidate windows of different widths between approximately 1.5 and 6 km" | `molecular_methods.py:240-245`; `config.py:274-277` | DIFFERS (upper bound) | Centres 1980–5820 m (step 240 m, 17), half-widths 240–1920 m (8) → 136 windows, widths 0.48–3.84 km. Window **start** ≥ 1500 m (gate), window **top up to ≈ 7.7 km**. |
| 22 | "time–altitude cells contaminated … part of the night are flagged and excluded, and the fit is performed on the time-cleaned mean profile" | `molecular_methods.py:450-487, 845-868` | MATCHES for the window/slope, **DIFFERS for the constant** | Cells with `ratio > max(median + 4·MAD, 1.25·median)` per altitude are flagged (upper tail only, so a layer present most of the night defines the median and is kept). The *window search, slope and all gates* use the cleaned mean (`clean_signal`, `:854-856`), but `cell_flag` is never returned to the caller: the **Klett/ensemble constant is computed from the uncleaned `rcs_mean`** (`calibration.py:826, 918-925, 1053-1072`). See §2 B1. |
| 23 | Composite score: R², near-unity scattering ratio, small shape residual, low temporal variability | `molecular_methods.py:679-690, 705` | MATCHES (incomplete) | `quality = R² − 0.25·|SR−1| − 0.20·resid% /100 − 0.10·ratio_std − 0.35·tCV − 0.20·rel_err/100 + 0.10·n_pts/max`; the SNR (ratio_std), slope-vs-median (rel_err) and window-size terms are not mentioned. |
| 24 | "If no window qualifies, the night is flagged as non-calibrating" | `calibration.py:967-980` | MATCHES | Flag −2. |
| 25 | Gate: slope > 0 | `molecular_methods.py:589`; `rayleigh_fit.py:57-59`; `calibration.py:984-997` | MATCHES | Eligibility gate; the post-hoc flag −7 can only fire if the method returned a non-positive slope (it cannot for v2). |
| 26 | Gate: "intercept must remain smaller than the slope" | `molecular_methods.py:590`; `rayleigh_fit.py:59`; `calibration.py:998-1011` | MATCHES literally, **vacuous** | `abs(intercept) < slope` compares b (in signal units, ~10⁻³ for CHM15k) with a (= C_L ≈ 10¹¹): dimensionally meaningless and always true → flag −8 is dead code. The physically meaningful offset test is the rel_error gate (slope vs. median ratio ≤ 15 %). |
| 27 | Gate: slope within 15 % of the robust median-ratio estimate | `molecular_methods.py:364-366, 591`; `calibration.py:1013-1028` | MATCHES | `rel_error = |a − median(y/x)|/median(y/x)` ≤ 15 % as a gate (v2) and again post hoc with `threshold_quality` → flag −2. |
| 28 | Window start above ~1.5 km; R² ≥ 0.4 | `molecular_methods.py:103-104, 588, 604` | MATCHES | 1500 m, 0.40 (eprof_v2 and v2.2). |
| 29 | Rayleigh-shape residual, scattering ratio relative to the cleanest window, temporal variability "must be small" | `molecular_methods.py:604-608` | MATCHES (thresholds undocumented) | residual_pct ≤ **16 %** (RMSE of the forced-zero-intercept fit / mean), ratio_std ≤ **0.40**, scattering ratio ≤ **1.15** relative to the *minimum* median ratio over clean windows (R² ≥ 0.5, `:380-386`), temporal CV of the window-mean ratio ≤ **0.8**. v2.2 adds a fallback tier (chi²_red ≤ 2.5, noise-de-biased reference, `:612-634`) — not in the paper. |
| 30 | Klett inversion, LR_aer 52 sr, LR_mol 8π/3, reference in the molecular window | `atmosphere.py:366-478`; `calibration.py:386-409` | MATCHES | Corrected Fernald/Klett signs (verified, §2 check B). Reference = window **mean** of beta_att/β_mol (not a point), reference index = window centre. |
| 31 | "Extinction is integrated starting at 100 m altitude" | `calibration.py:395` (`i_start_ext = 0  # from ground`); `rayleigh_fit.py:447-456`; `config.py:228` | **NOT IN CODE** | `z_start_ext = 100` exists in the dataclass but is read nowhere (grep). Klett and the transmission integral both start at the **first range gate** (≈ 15 m CHM15k; CL61 first native gate is 0 m → 12 m after binning). |
| 32 | "Points where the retrieved total backscatter falls below the molecular backscatter are excluded" | `calibration.py:415-419`; `config.py:231`; `rayleigh_fit.py:458` | **NOT IN CODE** (default) | `consider_points_lower_than_molecular = True` → the reset block is skipped; such gates are *kept* (negative β_aer; only its extinction is clipped to 0, `atmosphere.py:470`). `calculate_lidar_constant` excludes only gates with β_tot ≤ 0. |
| 33 | Per-gate C_L = signal / retrieved total backscatter, corrected for two-way extinction | `rayleigh_fit.py:443-472` | MATCHES | `cl = rcs / β_tot · exp(2·OD_tot(0→z))`, OD by cumulative trapezoid from gate 0; C_L = **median** over the window gates. |
| 34 | Ensemble: LR ±20 sr, window ±200 m, "several random subsets" | `calibration.py:67-70, 912-925, 1042-1072` | MATCHES (details undocumented) | LR ∈ {32,42,52,62,72} sr × shifts {−200,−100,0,+100,+200} m × {all, 4 random 70 % subsets (seed 42)} = 125 members, unweighted; members that fail (`_compute_cl_for_perturbation` → None: shifted window outside the range, < 3 gates, C_L ≤ 0) drop out silently. |
| 35 | Daily constant = median over the ensemble | `calibration.py:1093-1096` | MATCHES | `np.median(cl_values)`. |
| 36 | Uncertainty = robust 2σ from IQR and MAD | `calibration.py:1097-1106` | MATCHES | `u = 2 · max(IQR/1.349, MAD/0.6745)`. The per-member statistical term (`uncertainty_single`, `rayleigh_fit.py:476-477`) is computed and discarded; the NetCDF attribute `lidar_constant_uncertainty.comments = "(np.std(CL_all) + error_fit*np.median(CL_all))*2"` (`output.py:118`) describes that *discarded* formula — stale metadata. |
| 37 | "independent cross-check" with a slope-based estimate; reject if > 15 % | `calibration.py:1149-1199` | MATCHES mechanically; "independent" is overstated | `cl_slope = slope · exp(2·OD_aer(0 → window start))` (aerosol OD only, from the nominal Klett) vs `cl_median`. Both scale with the window signal; the ensemble value reduces to the window-mean ratio of the *uncleaned* profile, the slope value to the fit on the *cleaned* profile (§2 B1), so the test mostly measures intercept + time-cleaning effects. Flag −3. |
| 38 | Reject if uncertainty > constant | `calibration.py:1200-1202` | MATCHES | Flag −6 (checked after −3). |
| 39 | "A dedicated scattering-ratio test additionally rejects nights in which an aerosol layer is present within or just below the selected molecular window" | `calibration.py:1203-1225`; `config.py:208-209` | DIFFERS | Flag −9: ratio of the chosen window's slope to the **10th percentile of slopes of all clean candidate windows anywhere in 2–6 km** must be ≤ **2.0**. Nothing specific to "just below"; threshold is very loose (a window carrying up to 2× the cleanest window's signal passes; the real guard is the 1.15 gate inside the selection). |
| 40 | Flag 1 clear night / 0.5 partially clear | `calibration.py:1334`; `data_loader.py:1100` | DIFFERS (semantics) | `is_partial = n_cloud_masked > 0`: 0.5 means "at least one *retained* profile had a (high, > 6 km) cloud masked above cbh−500". A night where low-cloud periods were *removed* but no high cloud remains in the retained profiles gets flag **1**. |
| 41 | Negative codes: unsuitable conditions; not proportional; disagreement; missing CAMS; noisy/all-missing; uncertainty > value; invalid fit; aerosol below window | `flags.py:26-51`; `calibration.py` | MATCHES, incomplete | −1, −2 (two causes: no eligible window / rel_error > 15 %), −3, −4 (three causes), −5 (RCS all-NaN *or* all 125 perturbations failed), −6, −7, −8 (dead), −9. **Not in paper:** 0 "no data / no profiles in window", −10 "CAMS point too far (outside domain)", −11 classification veto (needs `--classify`; `ALC_CLASSIFY=0` operationally), −99 exception. |
| 42 | Constant, uncertainty, window limits and **quality flag** written to a NetCDF file per instrument and day | `output.py:106-180, 198-306`; `calibration.py:1326-1329, 1361-1373` | **DIFFERS** | The yearly NetCDF (`ALC_calibration_<wmo>_<id><year>.nc`) has **no flag variable** and **only successful nights are written** (a rejection returns before `write_calibration_result`). It carries lidar_constant, lidar_constant_uncertainty, bottom/top height (**m ASL**, `calibration.py:902-903`), calibration_method (0), calibration_version (200/220), laser_wavelength, housekeeping. Flags (incl. 0.5 vs 1) live only in the runner's per-stream `<key>_cal.csv`. |

---

## 2. Bugs and methodological weaknesses that can bias the paper's results

Ordered by expected impact on the published numbers. "Sign" = direction of the error on C_L.

### B1. The constant is computed from the *uncleaned* night mean, while the window/slope/gates use the time-cleaned mean — **sign +, up to +15 %**

* `molecular_methods.py:847-856` flags contaminated cells and fits `clean_signal = nanmean(masked)`; the result
  travels back as `MethodWindow.cell_flag` but `rayleigh_fit.py:72-103` (`_result_from_method_window`) drops it.
* `calibration.py:826` `rcs_mean = _agg(rcs_use, axis=0)` and `:918-925` (time subsets) never see the flag; the
  ensemble (`:1053-1072` → `_compute_cl_for_perturbation(rcs_mean=rcs_mean_t, …)`) and `calculate_lidar_constant`
  therefore use the mean *with* the flagged cells.
* Algebra: with `reference_value = ⟨beta_att/β_mol⟩_win` (`:388-391`) the Klett output in the window satisfies
  `cl(z) ≈ ⟨rcs_mean/(β_mol r²)⟩_win · exp(2·OD_aer(0→z))` — the per-gate constants are essentially the
  window-mean ratio of the **uncleaned** profile. `cl_slope` (`:1190`) uses the **cleaned** fit slope. The −3 test
  therefore passes whenever the cleaning changed the window signal by < 15 %, and the *reported* constant keeps the
  contamination. Scenario: cirrus/aerosol cells at 2× molecular during 10 % of the kept 5-min bins → cleaned fit
  unaffected, `cl_median` +10 %, error_pct ≈ 10 % < 15 % → flag 1, C_L **+10 %**. The paper's sentence "the fit is
  performed on the time-cleaned mean profile" is true only for the window selection.
* Fix: apply `cell_flag` to `rcs_use` before `rcs_mean`/the subsets (and recompute `signal`), or at least document.

### B2. The 3-h minimum is not enforced on clear nights — sign ±, precision

`data_loader.py:1000-1003` returns before the duration test. Short dark periods (SZA > 100° lasts 1–2 h at
55–60°N in June, polar-adjacent stations in May/August), or nights where the L1 file is partial, are calibrated
from < 1 h of data. No bias, but noisier constants concentrated in summer, which feeds the seasonal-amplitude and
availability statistics. The test is also profile-count based with `profiles_per_min` from the first two samples.

### B3. Cloud screen parameters differ from the paper — mostly availability, some sign +

Operational `z_low_cloud = 6000 m` (paper 4000 m) and no 50 % rule (§1 #9). With 6000 m, profiles under
mid-level cloud (4–6 km) are removed rather than masked — fewer calibrations on such nights (conservative). The
"3 h must survive" rule lets up to ~80 % of a winter night be cloudy. Since the `has_low_cloud` test looks only at
`cbh[:,0]`, a profile whose *lowest* layer is above 6 km but with a second layer is handled correctly (lowest layer
is layer 0), fine. The paper must simply state the real rule.

### B4. A uniform aerosol background inside the window is undetectable by construction — sign +, 5–15 %

All aerosol gates are *relative*: scattering ratio vs. the cleanest candidate window (`molecular_methods.py:380-386`),
flag −9 vs. the 10th-percentile slope (`calibration.py:1220`), temporal CV (fluctuations only), rel_error/residual
(curvature only). A free troposphere with a spatially uniform, steady scattering ratio R passes every gate with
R² high and `SR = 1.0`. Then the Klett normalisation (`reference_value`) absorbs the excess, β_aer ≈ 0 is retrieved
in the window, and **C_L = R · C_L,true**: a 5 % scattering ratio gives **+5 %**, independent of the lidar ratio.
The ±20 sr ensemble does *not* probe this: the LR only rescales the aerosol OD *below* the window (for AOD_below =
0.05 at 52 sr, 32–72 sr spans exp(±2·0.019) ≈ ±4 %, which is what the "sensitivity 2σ" reports). At 1064 nm
(β_mol(4 km) ≈ 1.1·10⁻⁷ m⁻¹sr⁻¹) a "clean" free troposphere typically has R ≈ 1.05–1.2, at 910 nm ≈ 1.03–1.1,
so the Rayleigh constants of the CHM15k/CL61 carry a plausible **+5–15 % high bias with a seasonal/latitudinal
modulation** (summer FT aerosol) that is absent from the stated uncertainty. The gate ceiling 1.15 means even a
*detected* 15 % excess relative to the cleanest window is accepted. This should be stated as the method's dominant
systematic (Wiegner & Geiß 2012 discuss it).

### B5. Fixed US Standard Atmosphere → season-locked bias — sign: C_L low in summer, ≈ 2–3 % peak-to-peak

β_mol ∝ p/T evaluated from the 1976 standard profile. With representative mid-latitude climatology:

| altitude | US Std p/T (hPa/K) | January (≈259 K/606 hPa at 4 km) | July (≈273 K/626 hPa at 4 km) |
|---|---|---|---|
| 3 km | 2.610 | 2.597 (**−0.5 %**) | 2.537 (**−2.8 %**) |
| 4 km | 2.351 | 2.340 (−0.5 %) | 2.293 (−2.5 %) |
| 5 km | 2.113 | 2.103 (−0.5 %) | 2.067 (−2.2 %) |

The standard atmosphere is ~1 % denser than the mid-latitude annual mean at 3–5 km, ~0.5 % denser than winter and
~2.5 % denser than summer → C_L **−0.5 % (winter) to −2.5 % (summer)**, i.e. ≈ 2–3 % peak-to-peak with a
summer minimum (≈ 3–5 % at Nordic stations with very cold winters). This cannot explain the 30–40 % cycle of
Section 5.1.3, but it is a network-wide, season-locked error with a defined phase that should be removed or
quoted before the residual is attributed to the instrument; the code already has the fix (`molecular_source:
"cams"`, `calibration.py:669-705`, one line of `options.json`). The `options.json:33` comment "<1 % difference" is
right for the annual mean only. The molecular-transmission part is < 0.1 %.

### B6. (eprof_v2.2 only) Noise propagation to the fit grid is wrong by the time-binning factor — tier-2 gates far more permissive than documented

* `calibration.py:646` measures `sigma_native` per **native** profile (15 s / 30 s) before binning (`noise.py:67-99`
  returns the per-profile σ). `_sigma_on_fit_grid` (`:167-203`) reduces by √(native gates per fit gate) and by
  √`n_profiles` where `n_profiles = rcs_use.shape[0]` is the number of **5-min bins** — the √(300 s/dt_native)
  (√20 CHM15k, √10 CL61) from the time binning is missing → `sigma_signal` too large by 3–4.5×, **chi²_red
  (`molecular_methods.py:318-323`) too small by 10–20×**. The gate `max_chi2red = 2.5` therefore admits windows
  whose true reduced χ² is ≈ 25–50, and the `ratio_std_excess`/`temporal_cv_excess` floors are inflated likewise.
* For the CL61 the WV correction (`:782`) divides the signal by T²_wv ≈ 0.6–0.8 *after* σ was measured → σ
  additionally too small by T² (opposite direction, ×1.5–2.8 on χ²). Net for CL61 ≈ 5× too permissive.
* Consequence: the "noise-aware recovery" tier admits nights v2 rejected, with aerosol-curved windows passing a
  physically meaningless χ² threshold → recovered nights biased high (B4 mechanism). Only matters if the server runs
  v2.2 (the CLAUDE.md claim); then the paper must describe this tier and the C8-vs-v2.2 validation should be
  re-read with this in mind (the 2.5 threshold was presumably *tuned* on the mis-scaled statistic, which is why it
  "works" empirically but is not interpretable as "2.5× the photon noise").

### B7. Klett/transmission integrals start at the first gate, extinction clipped one-sided — sign ±1–2 %

`calibration.py:395` (`i_start_ext = 0`), `rayleigh_fit.py:447-456`. The paper's "from 100 m" is not what runs. In
the incomplete-overlap/afterpulse region the manufacturer-overlap-corrected `rcs_0` is unreliable: under-corrected
overlap → β_aer < 0 → `ext_aer = max(0, β_aer)·S_a` (`atmosphere.py:470`) = 0 → aerosol OD of the lowest
≈100–300 m missing → C_L low (α = 0.1 km⁻¹ over 200 m → −4 % … typically −1 to −2 %); over-corrected overlap →
the opposite. The one-sided clipping also biases ext_aer positive wherever noise dominates (negligible inside the
window: ~0.1 %). Starting at 100 m would *exclude* real BL extinction (−1.7 % for α = 0.1 km⁻¹ over 85 m) — neither
choice is "right"; the point is that the paper and code disagree and the near-range OD is a ±1–2 % systematic
that the ensemble does not capture.

### B8. The CAMS humidity averaging window is 20–04 UTC, not the night actually used — CL61 only, out-of-domain stations

`calibration.py:664-666`. For European stations the solar offset ≤ 2 h is harmless (3-hourly CAMS). For Edmonton,
Ontario, Bonaire, Lauder/Auckland the WV column is taken from local afternoon/day; with a ~15–20 % WV effect and a
~10 % diurnal PWV cycle this is a ≈ 1–2 % bias on those CL61 constants. Also an inconsistency with the paper's
"computed for each night".

### B9. Uncertainty semantics — what the "2σ sensitivity 3.3 %" in Figure 3 does and does not contain

The 125-member spread (`calibration.py:1093-1106`) is dominated by the ±200 m window shift (vertical gradient of
the signal/molecular ratio, i.e. residual aerosol/baseline) and by the correlated 70 % subsets; the LR axis
contributes ≈ 0 in a clean window; photon noise of the slope (`std_err`) and the discarded `uncertainty_single`
are not included; B4/B5/B7 systematics are not included. It is a sensitivity index, not a 95 % confidence
interval; the paper calls it "uncertainty" and the NetCDF metadata describes a different formula (§1 #36).

### B10. Minor inconsistencies / dead gates (each < 0.5 % or availability-only)

* Fixed depolarisation 0.0301 (UV value) instead of ≈ 0.0279 at 910–1064 nm; CL61 molecular reference at 910.0 nm;
  CHM15k at 1064.0 nm; α_mol with 8π/3 vs the 8.50 sr implied by the phase function — combined ≲ 0.5 % (§1 #15).
* `|b| < a` gate dimensionally vacuous (§1 #26) → the paper lists a check that cannot fail.
* Flag −9 threshold 2.0 (§1 #39) — effectively redundant with the 1.15 selection gate.
* Tilt ignored (`tilt_angle` in L1 = 3–4° at Payerne): height overestimated by 0.2 %, slant OD +0.2 %.
* CL61 binned to 28.8 m, not 30 m; blocks index-based (no effect on a mean of a stationary signal).
* `profiles_per_min` from the first two samples; irregular cadence → wrong ±15 min / 3 h counts.
* Ensemble members dropping out (shifted window beyond the range, C_L ≤ 0) would skew the median; cannot happen with
  the 15-km E-PROFILE L1 range, but there is no guard.
* The flag-0.5 semantics (§1 #40) — "partially clear" is under-reported.
* LUT cross-sections indexed by AGL height while humidity is ASL (`water_vapor.py:505-512`): pressure-broadening
  evaluated ~1 km too low at a 1000-m station — second-order on a 15–20 % correction (≲ 1 %).
* `klett_inversion` uses rectangle sums (`atmosphere.py:420-422`), `calculate_lidar_constant` trapezoids — ~0 at 30 m.

### Checks requested that came out clean

* **A. a.g.l./a.s.l. consistency.** Station altitude enters exactly once: `altitude_grid = range + alt` for the
  standard atmosphere (`data_loader.py:377`, rebinned at `:828`) and `wv_alt_grid = alt + range` for the humidity
  interpolation (`calibration.py:766`, `water_vapor.py:512`). `range_alc` (AGL) is used for 1/r², molecular OD,
  Klett and the LUT index. No double counting, no omission.
* **B. Klett/Fernald signs and integration direction.** `atmosphere.py:444-470`: `qt = cum[ref] − cum[R]` =
  ∫_R^{Rc} β_mol(S_a−S_m)dz; numerator `X·exp(+2qt)`; denominator `X(Rc)/β(Rc) + 2·S_a·∫_R^{Rc} X·exp(+2qt)dr`
  with the signed cumulative so gates above the reference get the forward solution. This is the correct Fernald
  backward/forward form (Speidel & Vogelmann 2023 Eq. 10); the historical −2 exponent survives only under
  `sign_error_v10=True` (default False, not set operationally). ✓
* **C. WV applied as T² to the range-corrected signal.** `trans**2` averaged over the Gaussian (`water_vapor.py:554`)
  then `rcs / t2` (`calibration.py:782`). ✓ Grid aligned ASL ✓ (LUT index AGL, minor, above).
* **D. Night-mean construction.** nanmean of MAD-screened 5-min bins, no clipping of negative values (L1 `rcs_0`
  contains ~40 % negative samples at range and has no masked entries in the files inspected). ✓
* **E. Time-cleaning direction.** `flag_contaminated_cells` flags only `ratio > max(med+4·MAD, 1.25·med)`
  (`molecular_methods.py:484-487`): it cannot remove aerosol-free cells (no lower-tail flagging); it *keeps* a layer
  present > 50 % of the night (median = aerosol level) — acknowledged in the code (`calibration.py:211-215`). Not a
  bias source by itself; B1 is the problem.

---

## 3. Code behaviour the paper should mention, and paper content with no code counterpart

### 3.1 In the code, not in the paper

1. **Which version and which CL61 spectrum are live** (§0): the 910.55/0.188 spectrum exists only with the
   eprof_v2.2 cutover; v2.2's noise-aware fallback tier (`max_chi2red 2.5`, de-biased scattering reference, noise
   floors; `molecular_methods.py:131-137, 612-634`) is undocumented and currently mis-scaled (B6).
2. **No background/baseline subtraction**: the free intercept is used only as a gate/diagnostic
   (`subtract_background=False`, `calibration.py:379-382`); the constant is a ratio estimator that absorbs any
   additive baseline — the mechanism behind the known CHM15k near-top-of-window bias and the dark-profile option,
   which is **OFF operationally** (`dark_profile_file=""`, `ops/config.sh:37-39`, verified; `_dark_on_grid` returns
   None at `calibration.py:310-312`). The paper should say the Rayleigh constants are dark-free/baseline-uncorrected.
3. **CL61 units**: the E-PROFILE L1 `rcs_0` of the CL61 is Vaisala's calibrated β_att (`units = m^-1.sr^-1`,
   first gate 0 m, 4.8 m, 30 s — verified on `L1_0-20000-0-06610_C20260810.nc`), so the CL61 "lidar constant" is a
   dimensionless correction factor to the manufacturer calibration (≈ 1–2), while the output NetCDF labels it
   "V*m^3/sr" (`config.py:102-104`). The molecular reference for the CL61 is evaluated at 910.0 nm (two
   wavelengths coexist on purpose: 910.0 Rayleigh / 910.55 WV).
4. Ensemble details: 5 LR × 5 shifts × (1 + 4 random 70 % subsets, fixed seed) = 125 unweighted members; failed
   members dropped; nominal member = LR 52, shift 0, all profiles.
5. Flag −9 definition (10th-percentile reference, factor 2.0) and the extra codes 0 / −10 / −11 / −99; the −5 code
   covers "all-NaN" *and* "all perturbations failed".
6. Minimum sample sizes: ≥ 10 bins for the MAD screen, ≥ 5 for cell flagging, ≥ 5 finite gates per window.
7. CAMS handling: monthly file preferred, daily fallback; regional boxes for out-of-domain stations; nearest grid
   point; 20–04 UTC averaging; a station outside every box → −10.
8. Screening order (native-resolution cloud screen → binning → MAD screen) and the CL61 28.8-m grid.
9. Output: only successes are written to the yearly NetCDF; same-day/same-method rows are overwritten; every row
   carries `calibration_version`; window limits are m ASL; flags live in the per-stream CSV. Operational plumbing
   worth one sentence in Section 5.1.3/5.2: per-day caches + regression guard (`calibration/incremental.py`) and
   the daily merge with a 5-day self-healing backfill (`ops/config.sh:71-72`).
10. Tilt is ignored; polar-summer nights produce no calibration (flag 0).
11. Section 5.1.3 should state that the fixed standard atmosphere contributes a ≈ 2–3 % summer-low cycle (B5) and
    that an undetected uniform FT aerosol contributes a positive, possibly seasonal, bias (B4), before attributing
    the 30–40 % to the instrument temperature.

### 3.2 In the paper, not in the code

* Extinction from 100 m (gate 0 in code); exclusion of gates with β_tot < β_mol (kept); the 50 % cloudy-fraction
  rule and the 4000 m threshold (none / 6000 m); the 3-h minimum for clear nights (none); the quality flag in the
  NetCDF (absent; failures not written); the molecular profile "from co-located T and p" in the Figure 3 caption
  (US Std 1976); grey hatching = MAD outliers in Figure 3 (= ±15-min margin/fog); "windows between 1.5 and 6 km"
  (tops to 7.7 km); "aerosol just below the window" test (not implemented as such); an *independent* cross-check
  (not independent); "wavelength-dependent depolarization correction" (constant 0.0301).

### 3.3 Suggested minimal edits to Section 5.1.1 (if the code is kept as is)

* Replace "version 2.0" by the version actually deployed and, if v2.2, add one paragraph on the noise-aware fallback;
  state the CL61 spectrum consistently with it.
* Cloud screen: "profiles with a cloud base below 6 km (or Vaisala vertical visibility) and ±15 min around them are
  removed; the signal above CBH−500 m is masked in the remaining profiles; the night is rejected if fewer than 3 h of
  usable profiles remain (fully clear nights are not subject to the 3-h test)".
* Replace "Extinction is integrated starting at 100 m" by "from the first range gate" and drop the β_tot < β_mol
  exclusion sentence (or enable the option).
* Figure 3 caption: molecular profile from the US Standard Atmosphere 1976; grey = profiles within ±15 min of a
  low cloud; verify the profile counts (native vs. displayed/strided vs. 5-min bins) and the 977 m station altitude.
* Flags: add 0 (no data), −10 (station outside the CAMS domain), −99 (processing error); say the NetCDF contains
  successful calibrations only and that flags are kept in the per-instrument time-series file.
* Add the two dominant unquantified systematics (uniform FT aerosol, +5–15 %; standard-atmosphere season cycle,
  ≈ 2–3 %) to the uncertainty discussion, and note that the reported 2σ is a perturbation spread.
