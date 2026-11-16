# E-PROFILE ALC calibration v2.3 — Rayleigh calibration in polar summer (low-sun fallback)

*Draft for Confluence, 2026-10-21. Topic 3 of the v2.3 release. Option `polar_night_fallback`, default
off. Code: `calibration/rayleigh/polar.py`. Evaluation scripts: `validation/polar_summer/p1_geometry.py`,
`p2_recalibrate.py`, `p3_analyse.py`. Outputs: `C:/DATA/Projects/202610_v23/polar/`.*

## Summary

* **Problem.** The Rayleigh calibration only uses profiles with the sun more than 10° below the horizon
  (SZA > 100°). North of 56.6° N there is no such profile around the June solstice. Eleven Rayleigh-calibrated
  E-PROFILE streams (ten CHM15k and the Birkenes CL61) lose 30 (South Uist) to 165 (Hopen) nights a year, and
  stay without any calibration for up to seven months. No southern station is affected.
* **Option.** `polar_night_fallback` (default **off**; v2.2 results stay bit-identical — checked on the test
  fixtures and on 6775 real nights). On a night **without any darkness** it calibrates on the **6 h centred on
  solar midnight**, never with the sun more than 12° above the horizon (SZA ≥ 78°). Everything else is the
  operational chain. Such nights get flag **0.25** ("Low-sun success (polar summer)"); a night whose low-sun
  window is too short gets **−14**.
* **Availability.** Of 1599 nights without darkness (2025-01 … 2026-08), the fallback calibrates **154**
  (operational: 0; best global threshold, 90°: 111; any variant: 181, the rest is cloud). Calibrations rise by
  +52 % at Hopen and Andøya, +81 % at Jan Mayen, +27 % at Oslo; the longest summer gap shrinks e.g. from 152 to
  28 days at Andøya and from 216 to 51 days at Jan Mayen.
* **Solar background.** CHM15k background in the low-sun window: 1.0–1.2 × dark south of 63° N, 2.1–2.7 × dark
  at the four Arctic stations. Calibrating the SAME night on SZA bands from 80° to 100° (CHM15k: 964 band-night
  pairs at the polar stations, 1784 at the controls) shows **no solar bias beyond ~2 %**; all twilight subsets
  sit 1.3–3.4 % above the dark core, including the background-free 98–100° band, so this offset is not solar.
  The solar residual is visible only above the calibration window (r²/β_mol growth at 7–11 km).
* **Continuity.** At the edges of the no-night season, low-sun constants match the neighbouring dark nights to
  +1.5 % (CI −3.1…+6.1) once their warmer housing is accounted for (κ); scatter 16 % vs 12 %; reported
  uncertainty 4.4 % vs 4.1 %.
* **QC.** The far-range residual has the same distribution on low-sun and dark nights (cirrus and noise
  dominate), so it is kept as a diagnostic; the background is controlled by the SZA floor.
* **Recommendation.** Enable the fallback in v2.3 with the defaults and let the 0.25 nights into the Kalman
  (already wired). Do **not** lower the global threshold instead: it changes every short summer night of the
  network (+0.3…+0.8 %) and still leaves the Arctic stations without nights.

## 1. The problem: no night, no Rayleigh calibration

The Rayleigh calibration only uses night-time profiles. With `use_sza_night = 1` it keeps the profiles whose
solar zenith angle (SZA) exceeds `sza_night_threshold = 100°` (sun 10° below the horizon), inside a 24-h
window centred on solar midnight (`calibration/io/data_loader.py::filter_time_range`). The lowest sun of a
night is reached at solar midnight, where the SZA is

  SZA_max = 180° − |latitude + declination|.

Around the June solstice (declination +23.44°) SZA_max stays below 100° north of **56.6° N**. No profile then
qualifies, the night ends with flag 0 ("No profiles in nighttime window"), and the station has no Rayleigh
calibration for weeks to months. The same holds south of 56.6° S in December; no E-PROFILE station is that far
south (Lauder, 45.0° S, keeps SZA_max ≥ 111.5°; Auckland ≥ 119.7°).

![Lowest sun per night and nights lost per year](C:/DATA/Projects/202610_v23/polar/p1_geometry.png)

*Left: SZA at solar midnight through 2025 for the eleven Rayleigh-calibrated streams north of 56.6° N.
Right: nights per year without a single qualifying profile, as a function of latitude, for thresholds 100°
to 90°.*

Affected streams of the E-PROFILE L1 census (433 streams) that are Rayleigh-calibrated (CHM15k, CL61,
Mini-MPL). Nights per year without a single profile above the threshold (2025 geometry, operational formula):

| Station | Type | Lat [°N] | Lowest sun at the solstice (SZA) | No profile > 100° (period) | > 98° | > 96° | > 94° | > 92° | > 90° | Nights with < 1 h / < 3 h of SZA > 100° |
|---|---|---|---|---|---|---|---|---|---|---|
| Hopen | CHM15k | 76.5 | 80.1 | **165** (01 Apr – 12 Sep) | 155 | 145 | 134 | 123 | 112 | 166 / 170 |
| Bjørnøya | CHM15k | 74.5 | 82.1 | **155** (06 Apr – 07 Sep) | 145 | 134 | 123 | 112 | 98 | 156 / 161 |
| Jan Mayen | CHM15k | 70.9 | 85.6 | **136** (15 Apr – 28 Aug) | 126 | 114 | 102 | 88 | 72 | 138 / 144 |
| Andøya | CHM15k | 69.3 | 87.2 | **128** (19 Apr – 24 Aug) | 116 | 104 | 90 | 76 | 57 | 128 / 136 |
| Tórshavn | CHM15k | 62.0 | 94.5 | **81** (13 May – 01 Aug) | 64 | 42 | 0 | 0 | 0 | 82 / 96 |
| Flesland (Bergen) | CHM15k | 60.3 | 96.3 | **66** (20 May – 24 Jul) | 45 | 0 | 0 | 0 | 0 | 68 / 84 |
| Lerwick | CHM15k | 60.1 | 96.4 | **65** (21 May – 24 Jul) | 43 | 0 | 0 | 0 | 0 | 67 / 82 |
| Oslo | CHM15k | 59.9 | 96.6 | **63** (22 May – 23 Jul) | 40 | 0 | 0 | 0 | 0 | 66 / 82 |
| Birkenes | CL61 | 58.4 | 98.2 | **46** (30 May – 14 Jul) | 0 | 0 | 0 | 0 | 0 | 49 / 69 |
| Stornoway | CHM15k | 58.2 | 98.3 | **44** (31 May – 13 Jul) | 0 | 0 | 0 | 0 | 0 | 47 / 68 |
| South Uist | CHM15k | 57.4 | 99.2 | **30** (07 Jun – 06 Jul) | 0 | 0 | 0 | 0 | 0 | 35 / 60 |

Another 36 CL31/CL51 streams north of 56.6° N lose their nights in the same way, but they are not
Rayleigh-calibrated operationally (liquid-cloud method), so they are not affected. Relaxing the threshold
to 94° would rescue the seven streams south of 63° N but not the four Arctic stations, which still lose
90–134 nights at 94° and 57–112 nights at 90°. Full table: `C:/DATA/Projects/202610_v23/polar/p1_stations.md`.

## 2. How much sun is left at the lowest point of the night

The CHM15k reports its solar background per profile (`bckgrd_rcs_0`). Over the 10-min housekeeping cache of
report 19 (clear blocks, 2025–2026), relative to each unit's dark median (SZA ≥ 105°):

| SZA [°] | ≥ 94 | 92–94 | 90–92 | 88–90 | 86–88 | 84–86 | 82–84 | 80–82 | 78–80 | 70–75 |
|---|---|---|---|---|---|---|---|---|---|---|
| background / dark (median of 10 polar CHM15k) | 0.98–1.01 | 1.12 | 1.37 | 1.68 | 2.02 | 2.19 | 2.40 | 2.66 | 2.94 | 3.32 |

The solar part of the background collapses by two orders of magnitude between SZA 90° and 94°; above ~95° it is
below the detection floor (< 1e-5 photons per shot), consistent with report 19 (q16a: ≤ 0.2 % of the dark level
at SZA ≥ 98°). For comparison, at noon it is ~14 × dark.

![CHM15k solar background versus SZA](C:/DATA/Projects/202610_v23/polar/p1_background.png)

*Left: background relative to the dark median, per polar CHM15k. Right: the solar part in photons per shot;
it is similar for all units (1–4e-3 at SZA 80–88°) while the dark level differs by ×5 between units.*

Consequences for the eleven affected streams at the solstice (the worst night of the year), with the
recommended low-sun window of section 3 (6 h centred on solar midnight, SZA ≥ 78°), time-averaged over the
window:

| | lowest sun (SZA) | window | background / dark, window mean (edge) |
|---|---|---|---|
| Hopen | 80.1° | SZA > 78.0°, 4.4 h (floor reached) | 2.63 (2.75) |
| Bjørnøya | 82.1° | SZA > 78.0°, 5.9 h | 2.72 (2.98) |
| Jan Mayen | 85.6° | SZA > 80.6°, 6 h | 2.29 (2.72) |
| Andøya | 87.2° | SZA > 81.8°, 6 h | 2.12 (2.30) |
| Tórshavn, Flesland, Lerwick, Oslo, Stornoway, South Uist | 94.5–99.2° | SZA > 87.3–90.9°, 6 h | 0.99–1.16 (1.2–1.9) |
| Birkenes (CL61, no background variable) | 98.2° | SZA > 90.1°, 6 h | – |

* **Seven streams south of 63° N**: the window centre is as dark as a normal night; only its last hour on
  each side sees a sun up to 3° above the horizon. Only the 100° convention excludes these nights.
* **Four Arctic stations**: the sun stays 3–12° above the horizon and the background is **2.1–2.7 × dark** on
  average (solar part 0.9–2.6e-3 photons per shot). This is the regime the bias test of section 4.3 has to
  cover.

## 3. Method: the low-sun fallback

The operational selection is left untouched. The fallback only acts on nights that have **no darkness**, and it
replaces the empty selection by the profiles with the **lowest sun of that night**.

1. **Trigger (geometry, not data).** For the night key date D the module computes, with the same declination
   and hour-angle formula as `filter_time_range`, the lowest sun SZA_max and the dark duration (hours with
   SZA > `sza_night_threshold`). The fallback replaces the operational selection only if the dark duration is
   at most `polar_trigger_dark_hours` (default 0 h: no profile can qualify). A data gap during a real night
   therefore never pulls daylight profiles in.
2. **Selection = the lowest sun.** Keep the profiles with SZA above a floor
   `max( min(SZA_max − polar_sza_delta, SZA(±polar_min_hours/2)), polar_min_sza )`, where SZA(±h/2) is the
   SZA h/2 hours either side of solar midnight. With the defaults (`polar_sza_delta` 0, `polar_min_hours` 6,
   `polar_min_sza` 78) this is **the 6 h centred on solar midnight, never with the sun more than 12° above the
   horizon**. The window is contiguous and centred on solar midnight. It is applied by calling the
   operational `filter_time_range` itself with the floor as threshold, so the SZA formula, the night window and
   the profile bookkeeping are exactly the operational ones.
3. **Background control = the SZA floor.** The solar background is a function of the SZA (section 2), so
   `polar_min_sza` caps it (≤ 3 × dark at 78°). Where the floor shortens the window (Hopen around the solstice:
   4.4 h) the window is still used if it lasts at least the operational minimum clear duration
   (`min_time_range`, 3 h); otherwise the night is rejected ("sun too high").
4. **Everything downstream is operational.** Cloud screen (`min_time_range` 3 h for partly cloudy nights),
   L2-grid binning, molecular-window search (eprof_v2 / v2.2 gates), Klett, the 125-member sensitivity
   ensemble, all QC gates.
5. **Far-range residual (diagnostic; optional QC).** By day the CHM15k under-subtracts the solar background:
   an additive residual, constant in raw signal, that biases the pointwise constant by δ/(slope·p_mol(z)) and
   grows as r²/β_mol (report 19, q15e). The module measures it in the far-range band `polar_residual_band_m`
   (8–12 km; median over the profiles of each profile's median residual), converts it into the implied bias at
   the window centre and writes it into the message. It rejects a night only if `polar_max_residual_pct` > 0;
   the default is 0 (off) because the statistic does not separate low-sun from dark nights (section 4.5).
6. **Flags.** A low-sun calibration that passes every gate gets flag **0.25** ("Low-sun success (polar
   summer)", `calibration.flags.LOW_SUN_FLAG`) and a message such as
   `Low-sun (polar summer): SZA > 81.8 deg over 6.0 h (night max 87.2, dark 0.0 h), far-range residual
   +1.4 %; Successful (...)`. A rejected one gets **−14** ("Low-sun night rejected (polar summer)"). Neither
   value had a meaning before. 0.25 is a usable calibration: it is part of `SUCCESS_FLAGS`
   (`calibration/flags.py`, the dashboard mirror `monitoring/config.py`, and `_SUCCESS` in
   `scripts/run_network_calibration.py`, i.e. the Kalman input), while staying distinguishable from 1 / 0.5 in
   the CSV and on the dashboard (own label and colour). −14 is shown in the dashboard class "No night at this
   latitude". With the option off neither value can occur, so these one-token additions change nothing.

When the option is off nothing of this runs (`calibrate_rayleigh` imports the module only behind the option).
Two hooks in `calibration/rayleigh/calibration.py` (after the time filter and in the QC block) and one line at
the result are the only changes to existing code paths.

## 4. Evaluation with real data

### 4.1 Data, harness and validation

Every night 2025-01-01 … 2026-08-31 of the eleven affected streams (CHM15k and the Birkenes CL61), plus ten
control streams for April–September (CHM15k Eskdalemuir, List, Schleswig, Lindenberg, Schiphol A, Palaiseau,
Budapest A, Granada; CL61 Lindenberg and Uccle), was recalibrated from the native L1 (D:/E-PROFILE_L1_2026)
with the operational code: `calibrate_rayleigh`, options.json, molecular method `eprof_v2.2` (= the published
`calibration_version` 220), as in the report-19 harness (q16a): NetCDF writing patched out, nothing else
changed. Each L1 file was read once; all variants of a night ran on the same in-memory copy (assembly checked
bit-identical to `load_l1_data`). Variants per night:

* `op` — operational (SZA > 100°);
* `final` — the fallback with the recommended defaults (6-h window, SZA ≥ 78°); `polar` (4-h window or
  3° below the lowest sun) and `pol_h2` (2 h) as alternatives;
* `t98` … `t90` — the operational selection with the threshold lowered to 98, 96, 94, 92, 90°;
* `b80_85` … `b98_100`, `b100` — the same night restricted to the SZA bands [80, 85), [85, 90), [90, 94),
  [94, 98), [98, 100) and > 100° (on nights with ≥ 1 h of darkness, `min_time_range` 1 h as in q16a);
* `pol_force` — the low-sun window forced on nights with < 3 h of darkness (trigger study).

**Validation.** The `op` run reproduces the published v2.2 constants (`load_constants`, dashboard NetCDF):
**98.7 % of the 1507 nights common to both are bit-identical** (median ratio 1.000000; 21 streams, 14 of
them 100 %). The few mismatches are the known night-key mapping of the published files (q16a). With
the option **on**, every one of the **6775** real nights that have a dark period gives exactly the same flag
and constant as with the option off. The unit tests (`tests/test_polar_fallback.py`) check the same on the
bundled fixtures against constants recorded on v2.3-dev before the change, and that the polar code is never
called when the option is off.

### 4.2 Availability: nights calibrated when there is no night

Nights without darkness in the record (2025-01 … 2026-08, nights with L1 data) and how many of them each
variant calibrates (flag 1 / 0.5 / 0.25). Operationally the answer is 0 everywhere.

| Station | nights without darkness | SZA > 98 | > 96 | > 94 | > 92 | > 90 | low-sun 4 h | **low-sun 6 h (default)** | any variant |
|---|---|---|---|---|---|---|---|---|---|
| Hopen | 268 | 1 | 2 | 3 | 4 | 5 | 12 | **15** | 16 |
| Bjørnøya | 253 | 0 | 2 | 1 | 3 | 3 | 5 | **5** | 8 |
| Jan Mayen | 225 | 0 | 4 | 5 | 8 | 10 | 13 | **17** | 18 |
| Andøya | 213 | 2 | 4 | 3 | 6 | 10 | 30 | **34** | 37 |
| Tórshavn | 142 | 2 | 5 | 7 | 10 | 12 | 9 | **13** | 14 |
| Flesland | 120 | 4 | 10 | 13 | 17 | 21 | 13 | **18** | 21 |
| Lerwick | 72 | 1 | 3 | 3 | 3 | 4 | 3 | **5** | 5 |
| Oslo | 115 | 9 | 21 | 28 | 35 | 32 | 28 | **35** | 37 |
| Birkenes (CL61, 2026 only) | 44 | 12 | 10 | 8 | 8 | 8 | 9 | **6** | 16 |
| Stornoway | 87 | 3 | 5 | 4 | 6 | 5 | 5 | **5** | 7 |
| South Uist | 60 | 0 | 0 | 0 | 0 | 1 | 0 | **1** | 2 |
| **total** | **1599** | 34 | 66 | 75 | 100 | 111 | 127 | **154** | 181 |

![Availability](C:/DATA/Projects/202610_v23/polar/p3_availability.png)

*Nights without darkness in the record (grey labels) and nights calibrated with a lowered threshold (orange
to red), with the 4-h and the recommended 6-h low-sun window (blue). The black bar is the number of nights
calibratable by any variant: the limit set by the clouds.*

* The default low-sun window calibrates **154 of the 1599 nights without darkness** (operational: 0), 85 % of
  what any variant can reach (181): the rest is cloud. Lowering the global threshold reaches at most 111
  (at 90°), and the 7 streams south of 63° N gain nothing more below 94°: their nights are short, not absent.
* Over the whole record the number of Rayleigh calibrations rises by **+52 % at Hopen and Andøya, +81 % at Jan
  Mayen, +68 % at Tórshavn, +32 % at Birkenes, +27 % at Oslo**, +16–21 % at Bjørnøya, Flesland and Lerwick,
  +8 % at Stornoway and +2 % at South Uist (whose no-night season is one month).
* The longest summer stretch without any calibration (March–October 2025) shrinks from 216 to 51 days at Jan
  Mayen, 152 → 28 d at Andøya, 220 → 92 d at Hopen, 207 → 85 d at Bjørnøya, 86 → 28 d at Flesland, 66 → 19 d at
  Oslo, 104 → 53 d at Stornoway; in 2026: Hopen 124 → 38 d, Jan Mayen 124 → 45 d, Andøya 91 → 26 d, Birkenes
  44 → 26 d, Oslo 67 → 18 d. Tórshavn (171 → 131 d) and Bjørnøya 2026 stay long: their summers were almost
  always cloudy.
* Only the CL61 at Birkenes did better with a short window (2 h: 13 nights, 6 h: 6). A longer window more
  often contains a cloud, which makes the night "partly cloudy" and then requires 3 h of clear profiles
  (`min_time_range`); a fully clear 2-h window has no such requirement. With one summer of one unit this is
  not conclusive; over the ten CHM15k the 6-h window gives 148 nights against 77 for 2 h.

### 4.3 Bias versus solar elevation (same night)

To measure what the sun does to the constant, the SAME night is calibrated twice: once on its dark core
(SZA > 100°) and once on the profiles of one SZA band only (dusk and dawn together). The ratio is corrected
for the internal temperature with each unit's intra-night κ. Polar stations contribute their spring and autumn
nights (long twilights), the controls their summer nights.

| SZA band | 80–85° | 85–90° | 90–94° | 94–98° | 98–100° |
|---|---|---|---|---|---|
| **CHM15k, polar stations** (10 units) | +2.7 % [1.0, 4.7] | +1.3 % [−0.1, 3.2] | +2.4 % [1.0, 4.3] | +1.8 % [1.1, 2.8] | +2.1 % [1.0, 3.3] |
| nights | 182 | 203 | 205 | 227 | 147 |
| background / dark | 2.73 | 2.00 | 1.30 | 1.01 | 1.00 |
| **CHM15k, controls 37–55° N** (8 units) | +3.4 % [2.5, 4.6] | +2.9 % [1.8, 3.9] | +3.0 % [1.7, 3.7] | +2.2 % [1.2, 3.1] | +1.6 % [0.8, 2.5] |
| nights | 362 | 409 | 348 | 377 | 288 |
| background / dark | 1.64 | 1.40 | 1.09 | 1.01 | 1.00 |
| CL61, Lindenberg + Uccle (no κ) | +8.9 % [−8.8, 30] | +6.6 % [−13.8, 16.4] | +2.6 % [−0.9, 10.7] | +1.5 % [−7.8, 13.7] | +2.2 % [0.1, 9.4] |
| CL61, Birkenes (no κ) | +0.2 % [−6.7, 5.3] | −0.5 % [−2.5, 0.7] | +3.1 % [−0.1, 3.6] | −3.1 % [−5.6, 3.8] | −1.3 % [−7.8, 8.5] |

ln(C_band / C_dark), median [95 % CI]; background = same-night `bckgrd_rcs_0` ratio (CHM15k).

* **Every twilight subset is 1.3–3.4 % above the dark core — including SZA 98–100°, where the background is
  exactly the dark level.** This common offset is therefore not solar. It matches the twilight/deep-night offset
  found in report 19 (q16a: +2.1 % for SZA 100–106° against > 110°, +0.5 % in its dusk/dawn-symmetric
  estimator), attributed there to intra-night drift that a single κ does not fully remove; the band subsets mix
  dusk and dawn in proportions set by the clouds.
* **The solar part is the change ACROSS the bands**, relative to 98–100°: −0.8…+0.6 % at the polar stations
  down to SZA 80–85° (background 2.7 × dark), +0.6…+1.8 % at the controls (1.6 × dark). It does not grow with
  the background (the polar stations, with the higher background, show less), so **the solar effect on the
  constant is ≤ ~2 % and not resolved** at the night-to-night scatter (8–11 % robust sd per band).
* The **mechanism is visible** where it should be: the molecular-normalised profile of the SZA 80–85° subsets
  carries an excess that grows with height (+1 %-pt at 4.5 km, +3 at 5.5 km, +10 at 7 km, +18 at 9 km, +42 at
  11 km; right panel; at the polar stations alone +5.5 at 5.5 km and +29 at 7 km) — the r²/β_mol signature of
  the under-subtracted solar background (report 19, q15e). Inside the calibration window (2.5–5.5 km) it stays
  at 0…+3 %-pt (polar stations: up to +5.5 at the window top). From SZA 90° upwards it is within the ±2 %-pt
  that the 98–100° subset (no sun) also shows, up to 7 km.
* The CL61 samples are small (41–57 control nights, 10–17 at Birkenes) and noisy; they show no SZA trend
  either. Birkenes' low-sun window (SZA ≥ 90–94°) lies where even the CHM15k shows no background increase.

![Bias versus SZA](C:/DATA/Projects/202610_v23/polar/p3_bias_sza.png)

*Left: ln(C of the SZA subset / C of the dark core, SZA > 100°) on the same night, T_int-corrected, median and
95 % CI; dots are single nights. Middle: same-night background ratio (CHM15k) and change of the far-range
residual. Right: change of the molecular-normalised profile X(z) relative to the dark core (CHM15k), range on
the vertical axis; the grey band is the typical calibration window.*

### 4.4 Window length and threshold relaxation

**Window length** (nights without darkness, same night, T_int-corrected, relative to the 6-h default):

| window | nights calibrated | C vs 6-h window, same night | relative uncertainty (median) | median of the lowest SZA used |
|---|---|---|---|---|
| 2 h | 90 | −0.45 % (CI −1.0 … +0.3), n = 74 | 4.7 % | 96.3° |
| 4 h or 3° below the lowest sun | 127 | −0.06 % (CI −0.75 … +0.20), n = 118 | 4.0 % | 93.3° |
| **6 h (default)** | **154** | — | 3.9 % | 89.1° |

The window length does not change the constant; the 6-h window calibrates 21 % more nights than the 4-h one
and 71 % more than the 2-h one, at the same uncertainty.

**Lowering the global threshold instead** changes nights that already calibrate. On the nights with a short
dark period (< 6 h, n ≈ 355), relative to the operational SZA > 100° and T_int-corrected:

| threshold | 98° | 96° | 94° | 92° | 90° |
|---|---|---|---|---|---|
| ln(C_threshold / C_op) | +0.25 % (0.05…0.48) | +0.43 % (0.21…0.70) | +0.50 % (0.32…0.71) | +0.75 % (0.43…1.08) | +0.64 % (0.41…1.07) |

Each lowering dilutes the +1.5…+3 % twilight offset of section 4.3 into the night mean. Small, but it is a
systematic change of the whole network's summer series, for no gain at the four Arctic stations.

### 4.5 The far-range residual is not a usable per-night background QC

A per-night background QC would compare the low-sun night with the station's dark nights. The CHM15k background
variable (`bckgrd_rcs_0`) is not carried by the calibration's data structure, the CL61 has none, and in polar
summer there is no dark profile in the night's own files to compare with. The rcs-based candidate is the
far-range residual (section 3, item 5): reference-free, since it is zero on a dark night by construction. On
real data:

| | dark nights (`op`) | low-sun nights (default window) |
|---|---|---|
| n | 1565 | 154 |
| median implied window bias | +0.7 % | +1.4 % |
| robust sd | 4.4 % | 3.4 % |
| rejected at \|bias\| > 10 / 20 / 30 / 50 % | 18 / 7.3 / 3.6 / 1.9 % | 12 / 8.4 / 5.2 / 0.6 % |

(time-median estimator; the night-mean median gives the same picture.) Upper-tropospheric cirrus and noise
dominate the statistic: any threshold rejects dark and low-sun nights at the same rate, i.e. for reasons that
have nothing to do with the sun. The +0.7 %-pt shift of the median is the solar residual, consistent with
section 4.3. Same-night, the change of the residual correlates only weakly with the change of the constant
(Spearman +0.11, n = 2838; the top fifth of residual changes goes with +4 % in C, the bottom fifth with
+0.3 %), too weakly to act on single nights. **The residual is therefore reported in the message (diagnostic) and the QC is
off by default** (`polar_max_residual_pct = 0`). The background is controlled by the SZA floor (section 2),
whose validity is the band test of section 4.3.

![Far-range residual on dark and low-sun nights](C:/DATA/Projects/202610_v23/polar/p3_qc.png)

### 4.6 Continuity with the dark-night series, scatter and uncertainty

The CHM15k constant has a real instrumental annual cycle (±17 %, peak near the solstice; report 19 §7.13–7.14),
so polar-summer constants are **expected** to sit above the spring and autumn level, and the dashboard Kalman
(random walk) simply carries the spring level through the gap. The test therefore removes the network seasonal
course (median annual harmonic of 134 mid-latitude CHM15k: 17.2 %, peak DOY 176) and compares each low-sun
constant with the dark nights of the same unit and optical module within ±45 days, i.e. at the edges of the
no-night season. The same statistic for every dark night against its dark neighbours (leave-one-out) is the
control.

| | n | median [95 % CI] | with κ·ΔT_int removed | robust sd | reported 2σ uncertainty (median) | ΔT_int vs neighbours |
|---|---|---|---|---|---|---|
| dark nights (control) | 543 | −0.1 % [−1.0, +1.4] | −0.6 % [−1.5, +0.4] | 11.8 % | 4.1 % | 0.0 K |
| low-sun nights | 74 | +7.3 % [+1.0, +14.3] | **+1.5 % [−3.1, +6.1]** | 16.3 % | 4.4 % | +1.9 K |

* The low-sun constants are continuous with the dark-night series once the warmer housing is accounted for
  (the low-sun window includes evening hours: +1.9 K, times κ 1.1–3.2 %/K per unit). The κ term partly double
  counts the seasonal temperature effect already contained in the network curve, so the truth lies between the
  two columns; neither is significantly different from the dark-night control at the night-to-night scatter.
* Per unit: Andøya −1.0 % (n = 10), Hopen +3.3 % (n = 6), Tórshavn −8.4 % (n = 6), Oslo +15 % and Flesland
  +14 % (n = 26 and 17; +1.3 and +2.8 % after κ: both have a housing temperature that swings 9–10 K over the
  season, κ 3.2 and 2.5 %/K). Single nights at Jan Mayen, Bjørnøya and South Uist are outliers (−98, −56, +74 %)
  of the kind the Kalman outlier screen removes.
* Night-to-night scatter of low-sun constants is 16 % against 12 % for dark nights (more noise in short,
  cloud-broken windows); the reported ensemble uncertainty barely changes (4.4 % vs 4.1 %) — like for every
  night, it underestimates the real scatter (topic 1 of v2.3).
* In the time series (figure) the low-sun constants continue the seasonal curve of the dark nights. Oslo 2025
  (monthly medians relative to the module's dark-night median): April 1.33 and May 1.39 (dark) → June 2.17
  (low-sun) → July 1.83 (mostly low-sun) → August 1.54 (dark); the last dark nights of July (2.25, 3.11) are as
  high as the June low-sun ones. The June nights that look high used low windows (1.7–2.7 km, summer aerosol),
  and a lowered threshold (SZA > 90°) gives the same constants on those nights. Andøya 2025: low-sun May–July
  1.18–1.36 against 1.18 for the dark nights of August–September.
* The Arctic island units (Hopen, Bjørnøya, Jan Mayen) have large non-seasonal level steps through the year
  (window transmission 70–98 %, icing) which limit any continuity test there.

![Continuity](C:/DATA/Projects/202610_v23/polar/p3_continuity.png)

*Four series (CHM15k, each point divided by the median dark-night constant of its optical module; the yellow
bands are the nights without darkness). Right: each low-sun constant against the dark nights of the same unit
within ±45 days, after removing the network seasonal course of the CHM15k constant (17 %, peak DOY 176).*

### 4.7 Short nights (trigger)

The fallback is triggered only on nights without any darkness (`polar_trigger_dark_hours = 0`). Forcing the
6-h low-sun window on nights that have a short dark period (`pol_force`, dark < 3 h; 414 nights at the polar
and 52–55° N stations) calibrates more of them — 16 % instead of 11 % with 1–2 h of darkness, 14 % instead of
9 % with 2–3 h — and gives +1.2 % (CI −0.2…+2.5, n = 37) relative to the dark-only constant on nights where both
succeed (the twilight offset of section 4.3). Raising the trigger would therefore help short nights, but it
would also replace nights that calibrate today, network-wide down to ~52° N in June. It is left at 0 h; the
parameter allows it to be raised later if the operators want more summer nights at 55–57° N.

## 5. Limits

* **Short record.** One full summer (2025) and January–August 2026; Birkenes (CL61) only from November 2025. The
  availability numbers are per record, not climatological.
* **Clouds dominate.** Polar summer is cloudy on the Arctic islands (low stratus, fog). The fallback can only
  calibrate the nights that are clear around solar midnight; most summer nights at Hopen, Bjørnøya and Jan
  Mayen stay without a calibration.
* **Bias test coverage.** The same-night SZA subsets cover SZA 80–100°. The Hopen and Bjørnøya windows at the
  solstice reach SZA 78–80° (background 2.6–2.7 × dark on average, up to 3.0 × at the window edge), at the edge
  of the lowest tested band (80–85°, 2.7 × dark same-night): a mild extrapolation.
* **Temperature.** Same-night ratios are corrected with each unit's intra-night κ (report 19, q2a; 1.5 %/K
  where unknown). Twilight subsets are split between dusk (warmer) and dawn (cooler), so the correction is
  small, but κ errors of ±1 %/K leave ±0.5 % on the band ratios. No κ is applied to the CL61.
* **No per-night background QC.** See section 4.5: the far-range residual is dominated by upper-tropospheric
  cirrus and noise, so it could only catch gross failures. The protection against the solar background is the
  SZA floor, validated statistically, not a per-night measurement.
* **Continuity is confounded.** The CHM15k has a real instrumental gain cycle of ±17 % peaking at the solstice
  (report 19 §7.13–7.14), and the Arctic island units show large non-seasonal level changes (window icing,
  window transmission 70–98 % at Hopen). Low-sun summer constants are therefore EXPECTED to sit above the
  spring/autumn level; the comparison has to remove the seasonal course, and remains noisy.
* **SZA formula.** The operational formula has no equation of time: the window can be centred up to ~16 min
  away from the true solar midnight, i.e. ≲ 1° of SZA at the window edges at these latitudes (the lowest sun
  itself is unaffected).
* **The published NetCDF carries no flag.** A low-sun constant is written as a normal row of
  `ALC_calibration_<key><yyyy>.nc`; only the runner CSV and the dashboard see the 0.25 flag. If downstream
  users of the NetCDF must tell them apart, a `calibration_flag` variable has to be added at integration.

## 6. Recommendation

**Enable `polar_night_fallback` in v2.3 with the defaults, and count the 0.25 nights in the Kalman.**

1. **Strategy: a fallback, not a lower threshold.** Keep `sza_night_threshold = 100°` for every night that has
   a dark period, and use the lowest-sun window only on nights that have none. Lowering the global threshold
   changes every short summer night of the network (twilight subsets carry a non-solar +1…+2 % offset relative
   to the dark core, section 4.3) and still leaves the four Arctic stations without calibration for 57–112
   nights a year even at 90°. The fallback touches only the nights that today have no calibration at all.
2. **Window: 6 h centred on solar midnight, sun never above 12° (SZA ≥ 78°).** The window length does not
   change the constant (2 h / 4 h / 6 h agree within 0.5 % on the same nights), and the 6-h window calibrates the
   most nights (154 of 1599, against 127 for 4 h, 90 for 2 h and 111 for a global 90° threshold). The 78°
   floor binds only at Hopen and Bjørnøya around the solstice, where it shortens the window to 4.4–6 h.
3. **Background QC: the SZA floor, not a per-night statistic.** The background is a function of the SZA: in
   the recommended windows it is 1.0–1.2 × dark at the seven stations south of 63° N and 2.1–2.7 × dark at the
   four Arctic stations. Over that range the same-night test finds no solar bias beyond ~2 % (≤ 0.6 % at the
   polar stations at 2.7 × dark). The far-range residual cannot single out bad nights (same distribution on
   dark and low-sun nights), so it stays a diagnostic in the message (`polar_max_residual_pct = 0`); it can be
   switched on (e.g. 50 %, which rejects ~2 % of dark nights) as a gross-error guard if wanted.
4. **Flag 0.25 = usable, counted in the Kalman and on the dashboard.** The low-sun constants connect with
   the dark-night series (+1.5 %, CI −3.1…+6.1, once the warmer housing is accounted for), at a larger
   night-to-night scatter (16 % vs 12 %). The distinct flag keeps them identifiable (own colour on the
   dashboard; the CSV carries flag and message). For a series that must exclude them, filter on flag 0.25.
5. **Watch.** (a) Hopen, Bjørnøya and Jan Mayen are calibrated in summer only on the few clear nights — the
   gain is real but the summer level there rests on 5–17 nights per season. (b) A future station north of
   ~78° N would hit the SZA floor for weeks around the solstice and stay uncalibrated (by design). (c) If the
   NetCDF product must distinguish low-sun constants, add a `calibration_flag` variable at integration.
   (d) The remaining +1…+3 % twilight offset is common to every subset of a night that is not its dark core;
   it is part of the CHM15k intra-night behaviour (report 19), not of the polar fallback.

## 7. How to enable

`options.json` now documents the keys (`_comment_polar`) with the option OFF. To enable, set
`polar_night_fallback` to 1; the other values are the recommended defaults (they are also the
`CalibrationOptions` defaults, so absent keys behave the same):

```json
"polar_night_fallback": 1,
"polar_trigger_dark_hours": 0,
"polar_sza_delta": 0,
"polar_min_hours": 6,
"polar_min_sza": 78,
"polar_max_residual_pct": 0,
"polar_residual_band_m": [8000, 12000]
```

or in Python: `options.polar_night_fallback = True` on a `CalibrationOptions`. The option only acts when the
SZA night selection is used (`use_sza_night: 1`) and the station coordinates are known. For the v2.3 method
(`eprof_v2.3`, defined at integration) it is one switch of the three topics; no other change is needed: the
0.25 flag is already counted as a success by the runner's Kalman and by the dashboard.

Parameters, if they have to be changed:

| key | meaning | default |
|---|---|---|
| `polar_trigger_dark_hours` | the fallback replaces the operational selection when the night has no more than this many hours with SZA > `sza_night_threshold` (0 = only nights without any darkness) | 0 |
| `polar_min_hours` | length of the low-sun window centred on solar midnight | 6 h |
| `polar_sza_delta` | alternative window definition: SZA ≥ lowest sun − delta (the wider of the two is used) | 0 |
| `polar_min_sza` | never use a profile with a higher sun (SZA below this); caps the solar background | 78° |
| `polar_max_residual_pct` | reject when the far-range residual implies more than this bias at the window centre (0 = off, diagnostic only) | 0 |
| `polar_residual_band_m` | far-range band of that residual (m above the instrument) | 8000–12000 |

## 8. Reproduction

| step | command | output |
|---|---|---|
| geometry, background vs SZA (CPU only) | `python validation/polar_summer/p1_geometry.py` | `p1_stations.csv/.md`, `p1_geometry.png`, `p1_background.csv/.png` |
| recalibration of every night (8 workers, each L1 file read once) | `python validation/polar_summer/p2_recalibrate.py run --workers 8` | `runs/<stream>/<first-night>.csv` |
| the recommended defaults on the nights where the SZA floor shortens the 6-h window (Hopen, Bjørnøya; elsewhere `final` = the `pol_h6` variant) | `python validation/polar_summer/p2_recalibrate.py final` | `runs_final/<stream>.csv` |
| check that the read-once night equals `load_l1_data` | `python validation/polar_summer/p2_recalibrate.py selftest` | console |
| analysis, tables and figures | `python validation/polar_summer/p3_analyse.py` | `p3_*.csv`, `p3_*.png`, `p3_digest.txt` |
| unit tests (bit-identity off, synthetic polar night) | `pytest tests/test_polar_fallback.py` | — |

All outputs are under `C:/DATA/Projects/202610_v23/polar/`. Interpreter: the main checkout's `.venv`, run from
the repository root with `PYTHONIOENCODING=utf-8`.
