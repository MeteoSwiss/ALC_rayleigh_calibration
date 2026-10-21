# 19 — CHM15k calibration stability: temperature effect and short-term / diurnal behaviour (PLAN)

Status: **executed 2026-10-20 → results in [19_chm15k_calibration_stability.md](19_chm15k_calibration_stability.md)**. Branch `remote-dark`.

## 0. Questions

- **Q1 — temperature.** Does the CHM15k Rayleigh lidar constant C_L depend on instrument temperature, across
  *several* units? We want the effect separated from season, with a %/K slope per unit and for the network.
- **Q2 — short-term / diurnal.** Is the working assumption *"C_L has no diurnal or short-term pattern"* correct?

## 1. What we already know (memory, reports, code, literature)

| Source | Finding | Limitation |
|---|---|---|
| Hervo et al. 2016 (AMT 9, 2947) §5, Fig. 6–7 | Overlap error depends linearly on **internal temperature** (more correction when warm). It is up to −45 % at 250 m (35 °C) versus −20 % (<20 °C), and it vanishes by ~600 m. The manufacturer says the T-sensitive part is the **photodetector spatial sensitivity**. The model is explicitly meant to fix day/night variation. | **By construction it is blind to a range-independent (multiplicative) T effect**: the fit anchors at R_OK, so a far-range gain change is absorbed in the fit constant. The paper says nothing about C_L. |
| Report 09 §2.3 (A:\CHM15k, 13 stations, 2013–2024, MATLAB `analyze_chm15k_stability.m`) | Calibration CV 23 %, seasonal peak-to-peak 34 %. **corr(C, T_int) = +0.38** (median); housekeeping (T_int + laser pulses + optics) explains R² = 0.38. | **Season is not separated** from T_int. The figure and scripts are missing from the repo (only the digest survives). |
| Report 09 §3.2 Allan | CHM15k PAY: σ_A(1 d) 13 %, floor 11 % at ~21 d, ×1.2 only. LIN: 15 % → 8 % at ~1 yr. | The PAY floor ≈ 1-day scatter means real variability on day-to-month scales, *not* averageable noise. The source is unknown. |
| `rayleigh_availability/chm15k_constant_vs_temperature.py` (9 Oct, Payerne only, 2025–26) | Operational C: **+0.56 ± 0.59 %/K (r 0.07)**. Dark + background corrected: **+1.16 ± 0.42 %/K (r 0.21)**. | 1 unit. **TUB200009 → TUB140016 swap (Apr-2025) not segmented.** Season not removed. `temp_int` quantised to 1 K. The dark-corrected slope is larger, which hints at a **T-dependent dark** (below). |
| [remote-dark-campaign] M3 hood truth | Payerne A hood dark: D_slow 60–500 m = **21 rcs/°C**. The slow thermal dynamics lag by hours. | The dark is additive and moves the Rayleigh constant (+24 % at PAY), so a T-dependent dark gives a **T-dependent C_L**. This is a concrete mechanism. |
| [overlap-temperature-confound] | T_int is 74 % annual cycle at PAY. A naive overlap–T slope (+0.9 %/°C) drops to ~0 once day-of-year is controlled. Night-only methods only see *seasonal* T. `temperature_optical_module` is Peltier-locked (~33.9 °C): **use `temp_int`**. | The same trap applies to C_L. |
| [chm-cl61-800m-wv-extrapolation] / [amsterdam-b-overlap-deficit] | A real **diurnal** T effect exists in the near range: Amsterdam B bias 8 → 27 % over 18 → 35 °C. At PAY the CL61/CHM ratio is ~8 % diurnal at 300 m and anti-correlated with T_int. | Near range (<~800 m) only. Nothing yet tests the far range / C_L by day. |
| Roşu et al. 2026, preprints 202610.0850 (3 CHM15k: RADO-Galati, RADO-Bucharest/Măgurele, AGORA-Granada; Jul-2023 → Jul-2025, Cloudnet `beta_smooth`, 1-min HK) | T_int tracks T_ext (r 0.96–0.98), with a heater floor ~15–20 °C when T_ext < 0 °C. **Optical module is regulated at ~34.4 °C but its regulation SATURATES**: it climbs to 40–48 °C on hot afternoons once T_int > ~35 °C. The **detector** is stable (Galati ±0.02 °C) but rises by 1.8–5.7 °C at Bucharest and Granada when T_int > ~40 °C. 51 strong κ–T correlation regions, 47 of them in the warm season, mostly < 2 km. | **Not a calibration study.** Pixel-wise Pearson of κ against T over 10 days is confounded by the diurnal BL cycle, which the authors concede ("do not demonstrate direct temperature-induced modification"). No same-wavelength reference, no far-range test, no C_L. What we take from it: (i) the **regulation-loss regime** (hot afternoons) is the hypothesis to test, and (ii) Galati and Măgurele are in our data (A:\CHM15k\GAL, E-PROFILE `0-20008-0-INO` A/B). |
| Web literature scan | Rayleigh-based CHM15k system constant stable to ~7 % over multi-day (MDPI). Ambient T affects ceilometer calibration and "check frequently" (Raman comparison). LED calibration pulse unreliable below −8 °C. | Secondary sources; to be cited only after reading the originals. |

**Code facts:**
- The Rayleigh fit is **night-only** (`use_sza_night=1`, SZA > 100°, `calibration/rayleigh/calibration.py:561`), with one C per night. **So the pipeline cannot see a diurnal C_L cycle at all.** If one exists, it lands uncorrected in the daytime L2.
- The overlap-T correction (`validation/paper/overlap.py`, models in `ALC_OVERLAP_DIR`) is a no-op above ~700–840 m, so it cannot touch the Rayleigh window (≥ 2 km).
- CHM15k L1 carries useful housekeeping: `temp_int`, `temp_ext`, `temperature_detector`, `temperature_optical_module`, **`calibration_pulse`** (internal reference, photons/shot), `status_laser`, `status_detector`, `window_transmission`, `laser_life_time`, `laser_pulses`, `bckgrd_rcs_0`, `stddev`.

**Key physical hypothesis (from Hervo 2016 + Roşu 2026 + the hood):** C_L is stable **while the optical-module and detector regulation holds**. It can drift (i) through the T-dependent dark, at any T, and (ii) when regulation saturates on hot afternoons (T_int > ~35 °C).

**Regime (ii) only occurs by day. The night-only Rayleigh pipeline never samples it.**

### Provisional answer to Q2 (to be tested, not concluded)
The assumption is **probably wrong in part**.
- Near range (<~800 m): a diurnal, T-driven effect exists (Hervo 2016; Amsterdam B).
- Far range / C_L: there is **no evidence either way**. Night-only Rayleigh structurally cannot test it, and the PAY Allan floor (~11 % not averaging down) says something real moves on ≤ 3-week scales.
- Candidate mechanisms that would make C_L diurnal or short-term:
  - (a) photodetector sensitivity vs temperature (manufacturer statement), and optical-module/detector regulation loss above T_int ≈ 35–40 °C (Roşu 2026)
  - (b) T-dependent dark/baseline (hood: 21 rcs/°C)
  - (c) daytime solar background → photon-counting pile-up/dead-time loss
  - (d) window transmission (dew, rain, dirt)
  - (e) laser pulse energy vs temperature or ageing

## 2. Strategy: break the season ↔ temperature confound

Night-only Rayleigh constants alone cannot separate T from season, because aerosol in the window, night length and molecular SNR are all seasonal. Three independent levers are used:

1. **Synoptic T anomalies** (days): deseasonalise both C and T_int per unit, then regress anomaly on anomaly. A warm spell in March versus a cold one is not "season".
2. **Unit contrast**: units with different T_int amplitude in the same climate. Indoor/sheltered vs outdoor installs decouple T_int from season.
3. **Same-wavelength co-located CHM15k pairs** (aerosol cancels, works **day and night**): **Schiphol 06240 A/B/C/D** (4 units) and **Budapest 12843 A/B**. The hourly ratio κ_A/κ_B at 1–4 km (above every overlap) against ΔT_int(A−B) and against local hour/solar background. **This is the only clean test of diurnal C_L.**

## 3. Steps

### Step 1 — Literature & prior-work recovery (½ day)
- Digest the preprint once it is supplied.
- Re-read Hervo 2016 §5 for per-K numbers (Fig. 6: R² 0.45; Fig. 7).
- Recover `analyze_chm15k_stability.m` / `estimate_calib_from_housekeeping` and the missing `calib_stability_chm15k_drivers.png`, to see whether the 11-yr analysis removed season. If not, Step 3 redoes it.

### Step 2 — Housekeeping cache for all CHM15k (compute, home PC)
- Script `validation/chm15k_stability/hk_cache.py`. For 154 CHM15k streams, 2025–2026 (`D:/E-PROFILE_L1_2025|2026`), extract **10-min medians** of the housekeeping above plus `bckgrd_rcs_0`. Output: one parquet/npz per stream under `C:/DATA/Projects/202610_CHM15k_stability/hk/`. Size: about 154 × 2 yr × 105 k rows × 12 vars, roughly 1–2 GB.
- Read HK variables only, never `rcs_0`, so it stays cheap. Use a 16–24 worker pool.
- Era segmentation per stream:
  - `laser_life_time` resets (laser swap)
  - steps in `calibration_pulse` / `stddev` noise floor (optical-module swap; the Schiphol D trick from [remote-dark-campaign])
  - known events: PAY Apr-2025 TUB140016, Schiphol D 2026-07-11

### Step 3 — Q1: C_L vs temperature, network panel
Inputs: the dashboard nightly constants (`paper_v08/calib_dashboard/nc`, start/end_time), matched to the Step 2 HK over each night's window.

Per unit and era, three models of log C_n:
- **M0 (naive)**: κ·T_int + const. This reproduces the 9-Oct number and report 09.
- **M1 (season-controlled)**: κ·T_int + annual harmonic (sin/cos DOY) + linear trend (ageing). Report the **partial** κ.
- **M2 (anomaly)**: deseasonalised (30-d running median removed) C anomaly vs T_int anomaly. This captures the synoptic-scale T effect.

Additional analyses:
- Repeat with `temperature_detector`, `temp_ext`, `calibration_pulse`, `window_transmission`, `status_laser` as regressors, plus a multivariate version.
- Add a **regulation-loss** regressor: excursion of `temperature_optical_module` above its unit's set-point (night median), and of `temperature_detector` above its own. At night these are mostly 0, so Step 3 mainly verifies that, and Step 5 carries the test.
- **Selection-bias guard**: clear-night availability is seasonal, so weight by night and check the T-coverage per unit.
- **Dark arm (PAY only, the only hood)**: compare the operational C against the dark-corrected C to quantify how much of κ the T-dependent dark carries (hood D_slow 21 rcs/°C).
- Pooled meta-analysis: distribution of the per-unit partial κ, its median ± IQR, and how many units are significant (≥ 3σ). Check heterogeneity against module generation (TUB12/14/15/16/20…) and installation type.
- **Success criterion**: a partial κ distinguishable from 0 network-wide (e.g. median |κ| > 2·SE across ≥ 20 units). Otherwise conclude "no detectable multiplicative T effect beyond X %/K".

### Step 4 — Q2a: intra-night short-term stability (nights only)
- For high-SNR clear nights (PAY, LIN, Schiphol A–D, Budapest, plus ~10 geographically spread CHM15k [prefer-geographically-diverse-samples]), re-run the Rayleigh fit on **2-h sub-windows** with the **molecular window frozen** to the full-night window. Use a driver calling `calibration.rayleigh` with a per-sub-window time mask; core code unchanged.
- Metric: within-night C(t)/C_night against the within-night ΔT_int (night cooling, typically 2–6 K). The noise expectation comes from the per-gate σ propagation already in `calibration.py:168`; within-night scatter in excess of that is real.
- Known confound (M3 lesson): atmosphere evolution co-varies with nocturnal cooling. Mitigation: use only nights where the frozen window's R² stays ≥ the gate in every sub-window, and cross-check with Step 5 pairs, which do not share this confound.

### Step 5 — Q2b: diurnal stability from same-wavelength pairs (decisive)
- Sites: Schiphol 06240 A/B/C/D (6 pairs), Budapest 12843 A/B, plus Măgurele `0-20008-0-INO` A/B (the census gives identical lat/lon 44.3448/26.0123 for both).
- ⚠ Census coordinates are identical per WMO for Schiphol and Budapest too. The Schiphol units likely stand at different runway positions a few km apart, so the 1.5–4 km ratio is still the right layer (horizontally homogeneous aloft), but the actual positions need checking. Măgurele is the hot-summer site where Roşu 2026 saw regulation loss, so it is the most diagnostic pair if real.
- Key composite: R against **optical-module / detector excursion above set-point**. If R moves only when the regulation saturates, it is the regulation-loss signature. Hot-day vs normal-day split.
- Hourly mean `rcs_0` (the hourly `_streams_*` cache is fine for means, not for noise), cloud- and precipitation-screened (`cloud_base_height`, `sci`). Take each pair's ratio R(t) = median over 1.5–4 km (above all overlaps) of rcs_A/rcs_B, normalised by the ratio of their nightly Rayleigh constants.
- Analyses:
  - (i) **diurnal composite** of R by local solar hour, per season
  - (ii) R vs ΔT_int(A−B) and vs `bckgrd_rcs_0` (solar background → photon-counting dead-time)
  - (iii) day-minus-night ΔR, with a bootstrap CI over days
- Range-dependence check: R(z) at 1, 2, 3 km. Flat in z means a gain / C_L effect; height-varying means dark/background (additive) or overlap.
- Controls: units with identical module generation; the 2026-07-11 Schiphol D swap as a natural experiment (the T-response before and after).
- Cross-wavelength sanity only, not decisive: PAY CHM15k vs CL61 and LIN CHM15k vs CL61, day vs night at 1.5–3 km. The Ångström and WV confound is documented in [chm-cl61-800m-wv-extrapolation].

### Step 6 — Report
- `doc/reports/19_chm15k_calibration_stability.md` with inline landscape figures:
  - per-unit κ (M0/M1/M2) forest plot
  - network κ histogram
  - PAY/LIN/Schiphol C vs T anomaly
  - intra-night C(t) vs ΔT
  - Schiphol pair diurnal composites
  - R vs ΔT_int and vs background
- **Verdict** on Q1 (%/K, network ± spread) and Q2 (diurnal amplitude upper bound, short-term drivers).
- Recommendation: either keep one C per night, or add a T-dependent C_L term or a daytime correction to the L1→L2 hub.
- Update report 09 §2.3 (season confound), CLAUDE.md, and memory.

## 4. Practical

- **Machine**: home PC (32 cores). Step 2 is I/O-bound on D: (L1 2025+2026). Steps 3–5 are light. No CSCS needed. If CSCS is used for the 2025 L1, **CPU partitions only**.
- **New code** in `validation/chm15k_stability/` (R&D, branch `remote-dark`, not `main`). The operational pipeline is untouched.
- **Writes**: ~1–2 GB to `C:/DATA/Projects/202610_CHM15k_stability/` (HK cache + sub-window fits). Please confirm this location before the run.
- **Effort**: Step 2 ≈ 2–3 h wall-clock. Steps 3–5 ≈ 1 working day of analysis each. Report ≈ ½ day.

## 5. Open inputs needed from you
1. ~~Preprint~~ — received and digested (Roşu et al. 2026, §1).
2. The location of `analyze_chm15k_stability.m` / `estimate_calib_from_housekeeping` (A:\CHM15k\Calibration?), if you have it.
3. Any known module/laser swap dates for Schiphol A–C and Budapest A/B (L1 has no module id).
4. OK to write the ~1–2 GB cache under `C:/DATA/Projects/202610_CHM15k_stability/`?
