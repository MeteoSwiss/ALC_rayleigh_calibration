# Adversarial review of "The E-Profile Automatic Lidar and Ceilometer network", draft v0.8

Reviewed on 2026-10-10 against the code of this repository (branch `remote-dark`, HEAD 18876dd), the saved
outputs of the paper scripts (`C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08/`),
the local E-PROFILE archives and the cited literature. The review was done by one coordinating session and six
delegated reviewers (manuscript, Rayleigh code, cloud code, validation scripts, Kalman/overlap/dark/availability,
literature); every headline finding below was re-checked in the code by the coordinator. Nothing in the repository
was modified except the files of this folder.

Detailed reports (statement-by-statement tables with file:line evidence):

| File | Scope |
|---|---|
| [01_manuscript_level.md](01_manuscript_level.md) | 48 manuscript-level findings: status/transparency, notation, internal contradictions, overclaims, structure |
| [02_rayleigh_code_vs_section_5_1_1.md](02_rayleigh_code_vs_section_5_1_1.md) | Section 5.1.1 vs `calibration/rayleigh` (42 statements, 10 ranked bugs/weaknesses) |
| [03_cloud_code_vs_section_5_1_2.md](03_cloud_code_vs_section_5_1_2.md) | Section 5.1.2, App. 3–4 vs `calibration/cloud` (30 statements, physics/sign checks, Figure 4 rerun) |
| [04_validation_scripts_vs_section_6.md](04_validation_scripts_vs_section_6.md) | Section 6 vs `validation/paper/*` (chain table, every number cross-checked, 14 referee issues) |
| [05_kalman_overlap_dark_availability.md](05_kalman_overlap_dark_availability.md) | Section 5.1.3, Figure 7, App. 2, 5, 6 vs `monitoring/kalman.py`, `improve_alc_calib`, dark scripts |
| [06_literature.md](06_literature.md) | Reference-list integrity, claim-by-claim verification, missing literature, consistency with published results |
| [07_actris_chm15k_beta_att.md](07_actris_chm15k_beta_att.md) | ACTRIS Cloudnet raw CHM15k files: who writes `beta_att`, what `c_cal` contains, comparison with the E-Profile Rayleigh constant |

Convention: "C_L" is the lidar constant (signal = C_L · β_att · ...), "C" the O'Connor cloud coefficient
(C = C_applied / C_L). Percentages are relative differences unless stated.

---

## 1. Verdict

The science is largely sound and the validation numbers in Section 6 reproduce exactly from the saved outputs,
but the manuscript is not yet an accurate description of what runs. A referee who reads the code (the paper
promises reproducibility of L1 → L2) will find that several methods are described as they were planned rather
than as implemented, that three different Kalman filters and three different cloud aggregations are described or
used, that one appendix is built on a bug, and that a few physical statements have the wrong sign or an
unstated dominant assumption. None of the headline conclusions collapses, but the paper needs one more full pass
in which every method paragraph is rewritten from the code, the operational status is stated with dates, and
the known systematics (uniform free-tropospheric aerosol, dark signals, in-cloud extinction, Ångström exponent,
standard atmosphere) are given numbers.

## 2. Must-fix findings (ranked)

1. **Appendix 2 is built on a bug.** `fig_appendix2_availability.py:92-93` accepts only letters as instrument
   identifier (`[A-Za-z]+`); the 60 DWD CHM15k streams with identifier `0` (5653 station-month files in the hourly
   archive) are silently dropped. Consequences: CHM15k 113 ever / 89 active instead of ≈173 / ≈147, network
   458 / 350 instead of ≈518 / ≈408 (which then agrees with Fig. 1 and Section 3), every per-type statistic,
   the dips and the Finnish-drop sentence. Fix the two regexes, `--rescan`, regenerate text and figure.
   (05, §C.3)
2. **What is operational is not what the paper says, and the paper must say it with dates.** The committed
   configuration runs `eprof_v2` with the CL61 water-vapour spectrum 910.74 nm / FWHM 1.0 nm; the paper's
   "version 2.0" with "910.55 nm / 0.188 nm" is not a configuration that exists: the manufacturer spectrum only
   ships with the v2.2 cut-over block, which is commented out in `ops/config.sh:19-40`. The distributed Level-2
   product (CEDA/GTS) carries the v1.0 Rayleigh constants for CHM15k/CL61 and the default 1e8 for CL31/CL51
   (`doc/OPERATIONS.md`, CLAUDE.md), while every constant used in Section 6 is the dashboard's v2.2 Kalman
   series fetched on 2026-10-09 (a reprocessed history). Section 3's "the CL31 is calibrated by the manufacturer
   using liquid cloud signal" is wrong on both counts. State: which Rayleigh version and which spectrum are
   applied in the distributed L2 at the time of writing, whether cloud constants are applied or only published,
   and the switch date. (02 §0, 03 §4, 04 §4.9, 01 #1)
3. **Three Kalman filters, one description.** Section 5.1.3 describes the equations of the hub's NRT driver
   (`improve_alc_calib`: sinusoid if > 20 valid points, one global 1.5·IQR, per-type process noise with no CL61
   entry and a probable copy-paste error for MPL, older optical modules discarded as "outliers"); Figure 6 was
   produced by the reprocessing driver (constant model, rolling IQR, 87 outliers removed, non-causal, caption
   colours of predictions/updates swapped, "30–40 % peak-to-peak" is the half-amplitude); Section 6 uses the
   dashboard variant `monitoring/kalman.py` (constant model, rolling IQR, one relative noise set 4 %/step, no
   per-night weights) applied to BOTH Rayleigh and cloud constants (`run_network_calibration.py:498-520`). The
   per-night Rayleigh uncertainty of 5.1.1 is never used (author's own TODO), so "physically consistent
   uncertainty propagation" is not supported. Describe the filter actually used for the published constants,
   give its parameters, and say that it is applied to the cloud constants too. (05 §A)
4. **Cloud aggregation is described three ways, none of them implemented.** 5.1.2 "monthly ... median of monthly
   medians, inter-monthly SD" (no code does this); App. 3 "daily median C → Kalman"; Section 6 "daily Kalman".
   The product is the daily median of per-profile C (one 5-profile window, i.e. 2.5 min, suffices for a
   flag-1 day; Hopkin required ≥ 10 profiles) and the dashboard Kalman. (03 §0.3, #24, #29)
5. **Multiple-scattering sentence has the wrong sign** (5.1.2 and App. 4): B = 1/(2ηS) with η < 1 means the
   uncorrected integral is biased HIGH and the coefficient C LOW (C_L high); the code multiplies β′ by η < 1,
   shrinking B and raising C. The sentence is also self-contradictory ("enhances the apparent backscatter ...
   biases the integrated backscatter low"). (03 §2a)
6. **The η tables hide their dominant parameter.** The operational tables are reproduced bit-exactly by the PVC
   model at a_G = 5.5 µm AND in-cloud extinction α = 10 km⁻¹ (a value fitted to the legacy Hopkin tables, kept
   when a_G was revised); the paper lists FOV, divergence, λ and r_e only. Checked here: CL31 η at 1 km is 0.835
   (α = 10), 0.76 (α = 20), 0.71 (α = 30 km⁻¹, the function's own default and O'Connor's 15–20), i.e. C changes
   by +10 to +18 %. Droplet radius ±0.5 µm moves η by only ±1 %. The CL61 Rayleigh/cloud ratio 0.94 of 6.1 is
   partly tunable through α. Also: Table A4-1's CL61 row is a copy of the CL51 row (code gives 0.92/0.88/0.81),
   and r_e (effective) is used as Hogan's equivalent-area radius a_G (r_e/a_G ≈ 1.1–1.2). (03 §0.8, §2e)
7. **Rayleigh code bug: the constant is computed from the uncleaned night mean.** The time–altitude cell cleaning
   that 5.1.1 presents is applied to the window search, slope and gates only; `cell_flag` is set
   (`molecular_methods.py:866`) and never consumed, so the Klett/ensemble constant uses `rcs_mean` with the
   contaminated cells (`calibration.py:826, 918-925, 1053-1072`). The "method disagreement" test (flag −3,
   15 %) then measures the cleaning effect, and contamination below 15 % is published as flag 1 (sign +,
   up to +10–15 % on such nights). Fix the code or describe the cleaning as a window-selection step only.
   (02 §2 B1)
8. **Rayleigh method description vs code, several items.** Cloud screen: threshold 6000 m (paper 4000 m) and
   there is no 50 % cloudy-fraction rule (`max_ratio_cloudy` is never read; the real rule is "≥ 3 h of usable
   profiles survive"); the 3 h minimum is not applied on fully clear nights (`data_loader.py:1000-1003`);
   extinction is integrated from the first gate, not from 100 m (`z_start_ext` is dead); gates with
   β_tot < β_mol are kept (option default True); windows reach 7.7 km (not 6); the "intercept < slope" gate is
   dimensionally vacuous and cannot fail; the yearly NetCDF has no flag variable and holds only successful
   nights (flags live in the per-stream CSV); the molecular reference is the US Standard Atmosphere 1976
   (the Figure 3 caption says "co-located temperature and pressure"); Figure 3's "343 of 761 profiles" is a
   display-stride count (CHM15k 15 s → ≈2300 native profiles, ≤ 120 five-minute bins) and its grey hatching is
   the ±15 min cloud margin, not the MAD outliers. (02 §1 rows 4, 9, 12, 13, 21, 26, 31, 32, 42)
9. **Figure 4 is not the operational result.** The figure script omits the sub-cloud aerosol transmission
   correction that 5.1.2 presents as part of the method (−7 % that day) and applies the cloud tests before η;
   the operational pipeline on the same file gives n = 248, C = 1.29 (figure: 210, 1.39). "C_L = 0.72" in
   panel (f) is C_L/1e8. Panel (e)'s caption ("r_e measured by Cloudnet that day") describes a per-day η while
   the operational correction is a fixed 5.5 µm table. Regenerate from the operational code. (03 §2d)
10. **Section 6 processing is not as described.** Every step runs on hourly medians formed at read time; the
    cloud screen excludes whole hours on each instrument's own CBH with NO ±15 min expansion (window = 1 at
    hourly cadence, `intercompare.py:529-533`); the EARLINET match selects the hourly bins whose START lies
    within ±30 min of the EARLINET mid-time (sample centred ≈30 min late); the EARLINET band is effectively
    810/1110–5000 m; the statistics are pairwise, not common-hour; negative test values are kept but reference
    values ≤ 0 are dropped; log r is pooled over altitude and mostly measures the vertical gradient. Describe
    as implemented or fix and rerun. (04 §1, §4.11)
11. **The headline estimator is fragile where the CL31 is noise-dominated.** Palaiseau 500–3000 m: −28.5 %
    (paper, median relative difference) vs −10.6 % (median log-ratio on positive pairs) vs +3.7 % (ratio of
    medians); Payerne CL31: −0.6 % vs +14.0 % vs +32.1 %; 31 % of the Palaiseau CL31 test values are ≤ 0. Where
    the reference is noisy the pairwise median is biased negative (Amsterdam 2–3 km, the Fig. 14 caption
    describes this effect). Report a second estimator and the fraction of non-positive values; label the CL31
    layers above ≈1.5 km as noise-dominated rather than quoting percentages. (04 §4.4)
12. **EARLINET agreement is partly circular and mixes day and night.** Both the CHM15k constant and the SCC
    1064 nm elastic product are normalised to the same US-Std molecular model; at Măgurele the SCC
    normalisation height (5.5–5.9 km) lies inside the CHM15k window. Day/night split (recomputed): Leipzig
    +8.0 % day vs +3.0 % night; Măgurele B +5.9 % vs −4.9 % (the "+0.3 %" is their average); the CHM15k reads
    5–15 % higher by day against every reference. The sentence "an upper bound on the Rayleigh calibration
    error" is logically invalid (two biased quantities can cancel); "within about ±5 % between 1 and 5 km"
    contradicts the +7.7 % quoted in the same sentence; "±25 % seasonal" understates the Kalman range
    (−13/+36 % Leipzig, −30/+40 % Măgurele); the Leipzig instruments are 3.8 km apart. (04 §4.1–4.3, 4.12–4.13)
13. **910 → 1064 nm conversion has no uncertainty budget.** Ångström ±0.5 moves the boundary-layer agreement by
    ±6–7 % and the 500–3000 m headlines by ±3–6 % (Payerne CHM15k +23.9/+30.4/+36.8 % for α = 0.5/1/1.5); an
    uncorrected cross term (molecular return attenuated by particles, scaled by 0.856 instead of ≈0.46) gives a
    relative error ≈ −1.5·τ_p above the aerosol layer. The "within −3 % and +3 %" and "within −1 %"
    boundary-layer closures are agreement within the conversion uncertainty, not per-cent closure. (04 §4.5)
14. **Unstated dominant Rayleigh systematics.** (a) A spatially uniform, steady aerosol scattering ratio R in
    the window passes every gate (all gates are relative to the cleanest candidate window; ceiling 1.15) and
    gives C_L = R · C_true: +5–15 % plausible at 1064/910 nm, season- and latitude-dependent; the ±20 sr ensemble
    does not probe it (it only rescales the OD below the window, ±4 % for AOD 0.05). (b) The fixed US-Std
    atmosphere gives a season-locked bias of ≈2–3 % peak-to-peak (C_L ≈ 0.5 % low in winter, 2.5 % low in
    summer at 3–5 km); the code has `molecular_source: "cams"` ready. (c) The near-range OD from the first
    gate with one-sided extinction clipping is a ±1–2 % systematic. None is in the uncertainty discussion; the
    "2σ sensitivity" is a 125-member perturbation spread, not a confidence interval. (02 §2 B4, B5, B7, B9)
15. **Figure 7 shows a research method the text never describes.** Panel (a) is the clear-sky-noise "M2"
    retrieval of the manufacturer overlap (`validation/overlap_clearsky/`, not operational, made with an
    uncommitted script variant: committed chain gives 115/17 not 116/16); panel (b) is the operational
    Hervo-2016 temperature model evaluated at T_int = 25 °C, which Section 5.2 never mentions (the text stops at
    the daily overlap function, Section 4.3 says the hub applies only the a(z), b(z) model). Either describe M2
    with its validation and QC criteria, or replace panel (a). (05 §D)
16. **Payerne reference and transient.** The Payerne CL61 reference carries its own altitude-dependent dark
    (−7 to −28 % of molecular at 3–6 km) and the cloud method's systematics (η with α = 10 km⁻¹, S_liq, WV inside
    the cloud integral); its cloud Kalman starts on 2 cloud days (±8 % transient over the first two weeks);
    4 March is the first Rayleigh Kalman day, not "the start of the CL61 calibration record" (cloud record from
    28 Feb). The "+30.4 %" is not a CHM15k-only error, and the Conclusions' "discrepancies ... rather than by
    the calibration" contradicts App. 6, where the dark biases the operational Rayleigh CONSTANT by ≈20 %.
    (04 §4.6, 01 #27–28)
17. **Internal contradictions in numbers and cross-references.** Instrument counts: 423 (text) / 424 (Fig. 1
    legend) / CL61 6 (§3), 10 (Fig. 1, §6.1), 9 (App. 2), 13 (App. 6) / CL51 52 vs 50 / CHM15k 152 vs 151; Fig. 1b
    ≈420 stations in 2023 vs App. 2 368. Payerne CL61 ratio 0.83 (6.1, Feb–Jun) vs 0.78 (same paragraph and
    App. 6, Mar–Oct) is a period effect (the ratio falls to ≈0.75 in Jul–Oct) and must be said. Mini-MPL
    "68–91 %" uses days-with-data as denominator while every other percentage in 6.1 uses 149 calendar days.
    Lanžhot's missing cloud calibrations are laser-energy rejections (121 of 133 days flag −21), not "frequency
    of suitable liquid clouds"; Aosta's ratio 1.102 is outside the "±10 %" claimed for eight units. Captions:
    Fig. 12 "As Figure 10" → 11, Fig. 14 "Panels as in Figure 12" → 13, Fig. 15 "Layout as in Figure 10" → 11.
    App. 5: the CL61 session range omits the +2.4 % session of 29 Sep; "CL31 decreased steadily" is false
    (309 → 109 → 168 → 125 → 80 → 9 %); the CHM15k "−18 to −31 %" of App. 6 is not the gate-wise range
    (−14 to −54 %) and the statistic is undefined; A5-1 (MATLAB, hard-coded 5e11, nominal 1e8, hand-set window
    whose end time differs between script and rendered title) is inconsistent with A5-2's convention
    (CL31 dark 3.5× smaller). (01 #10–13, 04 §2–3, 05 §B–C)
18. **App. 6 statistics.** Spearman ρ = 0.68 on n = 10 falls to 0.58 (n.s.) without Payerne or Uccle, i.e. the
    "ranking agrees" rests on the two units it set out to explain; Edmonton's "20–25 % of molecular" is NaN in the
    output (no cloud constant); the ratio 0.94 is computed on 102 nights vs 0.78 on 62 (report the 62 common
    nights too); "~0.001 % measured" cloud immunity has no script; the far-field offset is "a third of the effect"
    in the text but "negligible (operator decision 2026-10-09)" in `run_paper_v08.py:116-117`; two hood products
    (v2 for A6-1, v3 for A6-2's 1.67) coexist. (05 §B.4)

## 3. Should-fix findings

- Notation: C is the lidar constant in 2.1 and its reciprocal in 5.1.2; for CL31/CL51 "C = 1/C_L" holds only
  in units of the 1e8 default (code: C_L = C_applied / C). Define C_L, C and the normalisation once. (01 #4–5)
- The L2 quality-flag health filter of 5.1.2 is dead code (`quality_flag=None`); the operational input is L1; the
  20/70 sr aerosol-transmission bounds are computed and discarded, so "to characterise the uncertainty" is not
  realised; the correction is a 5–13 % term, not "residual", and integrates from 40 m to 50 m below the PEAK
  (includes the lowest cloud gates; 8 % of B_aer from 40–100 m on the Figure-4 day); the upper CBH limit equals
  the integration top (2400 m), so the return of the highest accepted clouds is truncated (C biased high by a few
  %, up to ≈15 %); the 5 % "aerosol" test is partly a cloud-base-sharpness test and the dominant rejection;
  "a NIGHT without CAMS" should be "a day" (cloud calibration uses the full UTC day). (03 §1 #14, #19, #22, §2f)
- Rayleigh: v2.2's noise-aware fallback tier is undocumented and its σ is mis-scaled by the time-binning factor
  (√20 CHM15k, √10 CL61; CL61 σ also measured before the WV division), so `max_chi2red = 2.5` is 10–50× more
  permissive than its documented meaning — relevant if v2.2 is what runs. CAMS humidity for the Rayleigh WV
  correction is averaged over a fixed 20–04 UTC window (wrong by 8–12 h for Edmonton/Ontario/Lauder). Fixed
  depolarisation 0.0301 (UV value) instead of ≈0.0279. Tilt ignored (≤ 0.2 %). Flags 0, −10, −11, −99 missing
  from the paper; flag 0.5 means "a retained profile had a > 6 km cloud masked", not "partially clear night";
  0.5 nights are a large share of the "valid" nights in 6.1 (Aléria 50/50, Mobile 51/61, Granada 45/75) and the
  dashboard Kalman ingests them. (02 B6, B8, B10; 04 §3.4)
- Section 6.1: constants come from a local v2.2 rerun (verified identical to the dashboard night-for-night for
  Payerne), not from the dashboard — the intro sentence "except in 6.1 ... only the operational calibration" is
  right but should say so; the ±15-day pairing silently falls back to the whole-period median when < 3 reference
  values exist (Lanžhot, Camborne CL31, Palaiseau CL31, Getafe, Praha); a 34th stream (Lauder NZ) was run and
  dropped (flag −4 every day). The CHM15k cloud constants ARE η-corrected with the CHM15k's own table, so the
  32–42 % deficit is net of that correction — say it; the forced CHM15k/Mini-MPL cloud run is a demonstration of
  the operational configuration, not a method test (no saturation screening as in Hopkin). (04 §3, 01 #26)
- Section 5.2 must describe what is operational: the per-instrument temperature model a(z), b(z) keyed by optical
  module and driven by the internal temperature, how it is fitted from the daily overlap functions, and the
  near-range validity (Fig. 7b "negligible above 500 m" vs §3 "incomplete overlap up to 800 m"). (01 #16)
- Section 5.1.3 seasonal cycle: cite Buxmann (2024), Van Hove & Diémoz (2024), Bellini et al. (2024); remove the
  standard-atmosphere (≈2–3 %) and the uniform-aerosol contributions before attributing the cycle to the
  instrument temperature; show the housekeeping evidence or drop the claim. (01 #31, 02 B4–B5)
- Kalman as applied in Section 6 is causal, linearly interpolated to the hour and held at record ends; the CL31
  post-swap Kalman stays ≤ 6 % low for ≈3 weeks after the excluded 7–17 July (not only the excluded days are
  affected). (04 §4.7, 4.10)
- Appendix 4: give the MS uncertainty budget (α, a_G, r_e vs a_G, marine vs continental droplet sizes — the
  figure title itself says marine r_e ≈ 8–13 µm while 5.5 µm is applied network-wide); CHM15k and Mini-MPL η
  differ by 0.001–0.005 (identical only at two decimals). (03 §2e, 01 #34)
- Appendix 1: the ±0.10 nm spectrometer uncertainty and the 0.19 nm CL61 offset are not reconciled (resolution
  limits the width, not the centroid); the wavelength scale was checked only at 760 nm; the real evidence is the
  cloud CBH-slope test — make that the argument; quantify the CL31 ±1 nm peak wander in transmission. (01 #40–41)
- Section 7.1: the network climatology mixes CL31 (7.7 km range) with CHM15k (15 km) and Vaisala/Lufft CBH
  definitions that differ by up to 300 m (§2.2) — latitude and type are confounded; use medians of a skewed
  distribution; "late summer 3,500–3,800 m" does not match Fig. 16 (Sep ≈ 3150 m). (01 #20, #42)
- Section 4.3: the fixed 25 % L2 uncertainty needs a budget table (Rayleigh nightly spread, Kalman, cloud
  per-profile spread, S_liq ±0.8 sr, WV/CAMS humidity, η(α, a_G), S_aer 20–70 sr, overlap, dark). (01 #39)
- Missing sections: Code and data availability (package name/version/DOI, L12L2 not public, overlap_probe,
  raw2l1, constants DOI 10.5281/zenodo.18739645, dashboard URL, EARLINET/Cloudnet/CAMS/CEDA); a parameter table
  for 5.1.1/5.1.2/5.1.3. (01 #46)
- Completeness: author placeholders, affiliations, acknowledgements placeholder, Section 7 title missing,
  Figures 1, 2, 3, 5, 7, 9 never referenced in the text, Figures 8–9 undescribed, Figure 5 placeholder caption
  (the network ranking with median C = 1.434 is a key result that deserves a paragraph and a comparison with
  Hopkin 2019), broken sentence in 4.4 ("queryable through an, which"), duplicated sentences in 2.1, Section 4.3
  reads as a code audit, Section 7.1 reads auto-generated, "AI-generated content may be incorrect" alt-texts,
  61 unresolved co-author comments of Feb–Apr 2026 (several still unanswered: Kalman parameters/reproducibility,
  overlap details, Zenodo DOI, CL51 removal). (01 #21–23, #47–48)
- Facts: CL61 L1 gate spacing is 4.8 m (3276 gates to 15.7 km, 10 s at Lindenberg — verified on a local L1 file),
  not "15 m"; Hohenpeissenberg station altitude in L1 is 977 m (paper 974 m); CHM15k "6.5 KHz" → kHz;
  "E-Profile" (57×) vs the official "E-PROFILE" (1×); WIGOS ID "0-20000-0-06736A" (the A is the instrument
  identifier). (01 #44, 02 row 11)

## 4. What held up

- Every number of Sections 6.1–6.3 and of the captions of Figures 10–15 reproduces from the saved outputs
  (`summary_v08.json`, `method_comparison_summary.json`, per-stream CSVs), and the figure panel titles are the
  script's own strings. The Appendix 5/6 numbers reproduce with the exceptions listed above. (04 §2, 05 §B)
- The Klett/Fernald signs and integration direction, the a.g.l./a.s.l. bookkeeping, the two-way WV transmission
  applied to the range-corrected signal, the night-mean construction without clipping, the direction of the WV
  correction in the cloud method, the 300 m / factor-20 tests, the five cloud tests' numbers, S_liq, S_aer and
  the FOV/divergence values are all as stated. (02 "checks that came out clean", 03 §1)
- The η tables are bit-exact PVC output; the Cloudnet droplet climatology (5.8/5.6/5.4 µm medians) supports
  the 5.5 µm choice at the three continental sites. (03 §1 #9, #28)
- The 0.83 vs 0.78 Payerne discrepancy is a period effect, not an error in either number. (04 §3.6)
- The local v2.2 rerun reproduces the dashboard constants night-for-night (median ratio 1.0000), so the
  constants used in Section 6 are what the published pipeline produces. (04 §3.1)

## 5. Suggested order of work

1. Fix and rerun Appendix 2 (regex); reconcile every instrument count and period in the paper with one table.
2. Decide and state the operational status (version, spectrum, hub vs dashboard, cloud constants in L2).
3. Rewrite 5.1.1, 5.1.2, 5.1.3 and 5.2 from the code (use the statement tables of reports 02, 03 and 05 as
   checklists); add a parameter table and the flag table; fix the MS sign; add α to App. 4.
4. Decide on the two code fixes that change results (cell-flag cleaning in the Rayleigh constant; the ±15 min /
   bin-centre matching in the validation chain) and rerun, or describe the current behaviour honestly.
5. Add the missing systematics with numbers (uniform FT aerosol, standard atmosphere, η(α), Ångström, cross
   term, dark) and an uncertainty budget table; add the day/night split to 6.2.
6. Regenerate Figures 3 (caption), 4 (operational pipeline), 6 (caption/colours, state the driver), 7 (replace
   or describe), A5-1 (same convention as A5-2); fix the caption cross-references.
7. References: see [06_literature.md](06_literature.md) (missing entries, uncited entries, bibliographic errors,
   must-cite literature).

---

## 6. Actions taken on 2026-10-10 (manuscript v0.9)

Operator decisions: the paper describes the **dashboard version** (Rayleigh v2.2, CL61 water-vapour spectrum
910.55 nm / 0.188 nm, O'Connor cloud calibration of this repository); the only Kalman filter of the paper is
**improve_alc_calib**, to be run on top of both the cloud and the Rayleigh constants (the dashboard Kalman is for
display only); Figure 5 shows the lidar constant C_L; Figure 7 says "overlap estimation".

`doc/E-PROFILE ALC Paper_v0.8.docx` (the version shared for review) now carries all these changes as tracked
changes by "Claude (AI review)" plus Word comments for the points that need the authors' input; your own
tracked revisions are kept untouched. The untouched original is `doc/E-PROFILE ALC Paper_v0.8_before_review_2026-10-10.docx`
and `doc/E-PROFILE ALC Paper_v0.9.docx` is an identical copy of the new v0.8. Parts 1–4 of the edits carry the
timestamp 18:00 UTC, the reference fixes 20:00 UTC. Figure swaps are not tracked by Word and are listed here.

| Item | What changed | Source |
|---|---|---|
| Appendix 2 | Identifier pattern fixed (letters or one digit, `fig_appendix2_availability.py`); hourly archive rescanned, supplement re-read for the numeric identifiers (`appendix2_incremental_cache.py`); text, caption and figure regenerated with panel (a) only: 518 instruments at 457 sites (173 CHM15k), 408 active in March 2026 (was 458 / 113 / 350). The old MATLAB panels were already deleted by your own tracked change | `paper_v08/appendix2/` |
| Section 6 (Figs 12–16, A6-1, all numbers) | Rerun with the dashboard nightly constants smoothed by the improve_alc filter | `validation/paper/kalman_improve_alc.py`, `run_paper_v09.py`, `figs_paper_validation/paper_v09/` |
| New Figure 4 (Section 5.1.1) | Network ranking of the CHM15k Rayleigh C_L (154 instruments, median 2.96e11 vs default 3e11), with description paragraph; Figures 4–17 renumbered 5–18 (captions and cross-references as tracked changes) | `validation/paper/fig_v09_chm15k_network.py` |
| Figure 5 (was 4) | Rebuilt from the operational cloud pipeline (n = 248, C = 1.29, C_L = 7.7e7, equal to the dashboard) | `validation/paper/fig_v09_cloud_example.py` |
| Figure 6 (was 5) | Network ranking of the dashboard-version cloud C_L (272 instruments) | `validation/paper/fig_v09_cloud_network.py` |
| Figure 7 (was 6) | Kleine Scheidegg 2025–2026, dashboard nightly v2.2 + improve_alc filter | `paper_v09/fig_v09_kalman_kleine_scheidegg.png` |
| Figure 8 (was 7) | Relabelled "overlap estimation"; committed data give 115 kept / 17 rejected (v0.8 image: 116/16, uncommitted script) | `validation/overlap_clearsky/fig_paper_network.py` |
| Sections 5.1.1–5.1.3, 5.2, App. 4 | Rewritten from the code (see reports 02, 03, 05), MS sign fixed, α = 10 km⁻¹ stated, CL61 η row fixed, overlap temperature model and overlap estimation described | tracked changes |
| CL61 water vapour, Section 5.1.1 | "biases high by 15–20 %" replaced: with the operational line the uncorrected constant is 3–4 % LOW (Payerne −3.8 %, Lindenberg −2.7 %) | `diag_v22_l55s008` vs `diag_v22_nowv` |
| References and citations (report 06, Parts A and B) | List rebuilt in Copernicus style, alphabetical, with DOIs (78 entries): entries added for every work cited in the text (D'Amico, Engelmann, Freudenthaler [discussion paper], Frisch, Madonna 2015/2018, Miles, Nicolae 2013/2026, Wandinger, US Standard Atmosphere, CAMS, Cloudnet, WMO-No. 8 and 306, ISO/FDIS 28902-4, A-Profiles, Sauvageat 2026 dataset, Le 2026 reply to RC1, Hogan 2008, Eloranta, Klett 1981, Fernald, Speidel & Vogelmann, Heese, Campbell 2002, Flentje 2010, Emeis, Illingworth 2007/2015, Cimini, Roininen SOP, Muñiz-Rosado, Peuch); phantom Hogan (2006) and Klett (1983) replaced; Stevens (2019) corrected to its DOI; Hoffmann 2023, Pinto 2021 and UNFCCC 2022 removed (never cited, unrelated). In-text fixes: 2010 ash-network sentence (DWD ceilometer network), CBH differences 50–100 m (Wagner) and Martucci 160/569 m, Spänkuch thresholds, CL61 analogue APD and factory cloud calibration, CHM15k seasonal cycle (Madonna, Van Hove & Diémoz, Buxmann instead of Bellini 2024), Kotthaus reason, BUFR 0 02 121 once per profile, 0.188 nm attributed to the Le et al. reply, citation forms (Münkel, Wiegner & Gasteiger, An 2020, Bellini et al.) | `edits_part5_references.py` (session scratchpad), report 06 |
| Word comments (12) | Points needing the authors: CHM15k CBH algorithm change source, Görsdorf 2016 location and "300 m", ICAO quotation, WMO quotation wording, hub-guide version, 2017 WMO agreement, HITRAN edition, Nicolae 2026 identity, Freudenthaler 2018 status, EARLINET 10–20 % source, Vaisala DA10 guide year/URL, summary of the reference changes | comments by "Claude (AI review)" |

## 7. CHM15k manufacturer calibration (c_cal) — ACTRIS check

The E-Profile Level-1 files do not contain the CHM15k firmware calibration: raw2l1 converts `beta_att` back to
`beta_raw = beta_att / c_cal` and does not write `c_cal` (0 of 152 L1 files checked). On ACTRIS Cloudnet, 4 of the
24 CHM15k with data since 2023 write the `beta_att` format (Galati, Limassol, Cluj, Munich LMU). In all four,
`c_cal` is the generic Lufft default 3.2e-12 (1/c_cal = 3.125e11), not a factory calibration of the instrument,
and none of the four is an E-Profile instrument. For the 154 E-Profile CHM15k, C_L × 3.2e-12 has a median of 0.95
(p10–p90 0.57–1.28; 56 % within ±20 %): the default is right for the network mean but off by up to a factor of
about 2 for single instruments. A per-instrument comparison of manufacturer and E-Profile calibration is therefore
not possible with the available data, and no comparison was added to the paper. Details in report 07.

Open points for the authors are listed in the final message of the session and in reports 01–07.
