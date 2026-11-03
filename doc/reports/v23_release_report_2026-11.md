# Calibration release v2.3: the network recompute and its differences to the v2.2 release

Date 2026-11-03. Status: recompute complete (Rayleigh + cloud, 434 streams, 2025-01-01 to 2026-11-01); the
OmB / sensitivity add-ons run behind it; nothing transferred to zueub434 yet. Companion documents: the plan
`doc/v2.3/network_rerun_plan_2026-11.md`, the eta decision `eta_switch_decision_2026-11-02.md`, the three v2.3
design notes `doc/v2.3/uncertainty.md`, `vertical_consistency.md`, `polar_summer.md`, and the cutover procedure
`doc/OPERATIONS_v22_cutover.md`. Every number below comes from `validation/v23_release/compare_v22_v23.py` run on the
local copies of the two release trees (`C:/DATA/Projects/202611_v23_release/calout_v22_rel`, `calout_v23_rel`);
its tables are in `C:/DATA/Projects/202611_v23_release/compare/` and its figures in `figs_v23_release/`.

## 1. Summary

The release bundles four changes decided in October and November 2026: the repaired Rayleigh fallback tier
(`eprof_v2.3`: measured noise for CHM15k and CL61, vertical-consistency test for the CHM15k, tier-2 window top capped
at 5 km AGL for the CL61, Payerne CHM15k on tier 1 only), the uncertainty budget (`uncertainty_model v2.3`, five
components written to the NetCDF, Kalman weighted by it), the polar-summer low-sun fallback (flag 0.25, version
231), and the second generation of multiple-scattering tables for the cloud calibration (PVC, a_G 3.3 um, alpha
20 /km for the three Vaisala types, S unchanged at 18.8 sr, version 1100). Compared with the v2.2 release over
the common window 2025-01-01 to 2026-08-13:

* **Rayleigh, strict-tier nights are bit-identical**, as designed: 22 828 CHM15k, 653 CL61 and 577 Mini-MPL nights
  (83 %, 78 % and 100 % of the v2.2 successes) carry the same constant in both releases. Every night that moved
  or disappeared came from the fallback tier, verified by reruns with the fallback switched off (appendix B).
* **CHM15k: 8 % fewer calibrated nights** (27 640 to 25 375). The vertical-consistency test removes 2 749 tier-2
  rescues that v2.2 accepted and adds 330 new ones; the 2 063 tier-2 nights kept in both releases moved by a
  median of -2.7 % (P10-P90 -11 to +2 %). The Kalman best estimate changes by -0.5 % (median over 153 streams).
* **CL61: 16 % fewer Rayleigh nights** (842 to 707) from the 5 km cap on the fallback windows; the common nights are
  unchanged (31 moved by a median of -0.1 %). Per stream the Kalman moves by +0.4 % (median), +11 % at Payerne,
  where the cap removes the summer high windows that are biased low without a dark correction.
* **Mini-MPL: identical** (577 nights, Kalman +0.1 %).
* **About 160 polar-summer low-sun nights** (flag 0.25) at 11 stations north of 57 deg N that had no calibration at
  all before. Their level is 1.2 to 1.5 times the dark-night level at the Arctic stations (summer against winter
  nights of the CHM15k annual cycle), 2.6 at Oslo and 3.0 at South Uist (one night): those two need a look before
  the constants are trusted.
* **Uncertainty: the budget is three times the ensemble spread.** On identical nights the 2-sigma relative
  uncertainty is 14 % (CHM15k, CL61) and 7 % (Mini-MPL) instead of 5-6 %; the vertical-consistency term dominates
  (10 % of 14.7 % at the CHM15k). The Kalman posterior standard deviation is unchanged (9 % at 1 sigma).
* **Cloud: -1.6 to -1.9 %** on every Vaisala day (medians 0.981 CL31, 0.983 CL51, 0.984 CL61), exactly the table
  change (-1.0 to -2.6 % between cloud base 0.5 and 2 km). The day count moves by -0.35 %: 288 days fall out
  (244 of them to the aerosol-below-cloud filter, which sees a 2 % higher calibrated aerosol load), 31 come in.
* **Three Czech streams changed instrument inside the window** (Plzen 2025-07-16, Tusimice 2026-06-25, Libus
  2026-06-26: CL31 replaced by CL51 under the same stream identifier). One census type per stream cannot describe
  them: v2.3 calibrates their whole record as CL51 (the newer unit), so the CL31 part carries the CL51 eta table
  (+3 to +4 % on C_L against the CL31 treatment). This is the one real defect found; it needs a per-day instrument
  type in the runner and a Kalman restart at the swap (section 10).
* **The archive now runs to 2026-11-01 for the whole network** (the v2.2 release ended 2026-08-13), so the cutover
  gap fill on zueub434 is one or two days.

Verdict: the recompute does what the plan says and nothing else; the three Czech streams and the two low-sun
outliers are documented limits, not blockers. Ready for the archive check and the transfer (section 11).

## 2. What was run

| item | value |
|---|---|
| code | commit `46bb407` of `v2.3-dev`, staged read-only at `/scratch/mch/mhrvo/alc_v23_code` (`PYTHONPATH`), symlinked `~/alc_v23_code` on Linard |
| jobs | Linard (CSCS, CPU partition `postproc`), per-stream arrays: **337836** (434 streams, 2026-11-02 23:00Z), **339769** (the 275 Vaisala streams again, 2026-11-03 08:23Z), **339770** (OmB + sensitivity, `--no-cal --sens --omb`, chained `afterany`); 337837 cancelled |
| scripts | `ops/cscs/alc_v23_release_array.sbatch`, `alc_v23_addons_array.sbatch`; keylists `/scratch/mch/mhrvo/keys_v23_release.txt`, `keys_v23_vaisala.txt` |
| runner call | `scripts/run_network_calibration.py --stream <key> --start 20250101 --end 20261101 --methods rayleigh,cloud --force --ignore-coverage --workers 1` |
| switches | `ALC_MOLECULAR_METHOD=eprof_v2.3`, `ALC_UNCERTAINTY_MODEL=v2.3`, `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`, `ALC_POLAR_NIGHT_FALLBACK=1`, `ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'`; `ALC_ETA_TABLES` unset (= generation v2); `ALC_DARK_PROFILE` unset (network dark-free, operator decision 2026-08-17); `station_molecular_params` in `options.json`: Payerne CHM15k A `use_noise: false` |
| data | L1 `/scratch/mch/mhrvo/E-PROFILE_L1` (12 000 files per month, 380-390 stations, complete to 2026-10-31); CAMS `CAMS_Monthly_04` (0.4 deg monthlies) then `CAMS` (dailies) |
| output | `/scratch/mch/mhrvo/E_PROFILE_calout_v23_rel/<key>/` (`<key>_cal.csv`, `_kalman.csv`, `_hk*.csv`, `_status.csv`, `<year>/ALC_calibration_<key><year>.nc`), `RELEASE_PROVENANCE.txt` |
| v2.2 reference | `E_PROFILE_calout_v22_rel` (`alc_v22_release.sbatch`, `eprof_v2.2`, same WV spectrum and CAMS, 2025-01-01 to 2026-08-13) |
| comparison | `validation/v23_release/compare_v22_v23.py` (this report), `scan_type_changes.py` (section 8.3) |

**Incident.** The code was staged with an rsync that excluded `*.nc`, which dropped the water-vapour absorption
table `calibration/data/abs_cross_wv_910nm.nc`. In the first array every 910 nm calibration (CL61 Rayleigh and
cloud, CL31 and CL51 cloud: 275 streams) failed with flag -99 (`FileNotFoundError`), while the CHM15k and Mini-MPL
Rayleigh results were valid. The file was restored in both code copies and the 275 streams were recomputed in place
by array 339769 (`KEYS_FILE` subset, `--force`); the add-ons were re-chained behind both arrays. The sbatch header,
the handoff and the project memory now carry the rule: never exclude `calibration/data/*.nc` when staging the code.

## 3. Completeness and coverage

* Streams: v2.3 434 (the whole 2026 census: 154 CHM15k, 210 CL31, 52 CL51, 13 CL61, 5 Mini-MPL); v2.2 433. The
  Zeebrugge CL51 (`0-20000-0-06418_A`) is new since the v2.2 census. Nothing is missing.
* Coverage: 392 of 434 streams have data to 2026-11-01, 28 more to October; 10 streams end earlier because their
  L1 does (the Lindenberg CL61 `0-20000-0-10393_C` on 2026-06-15, a few on 2026-02 to 2026-06). The v2.2 release
  ends on 2026-08-13 for 402 streams, so the comparison window is 2025-01-01 to 2026-08-13 and the v2.3 tree adds
  2.5 months: in September and October 2026 it holds about 1 650 Rayleigh and 4 200 cloud successes per month.
* Every Vaisala stream of the second pass finished with return code 0 and no traceback (275 of 275; the two
  Belgian CL51 streams, Uccle and Diepenbeek, took 3.5 h each over the 22-month window).

![Rayleigh calibrations per month and type, both releases](figs_v23_release/fig_availability_month.png)

*Figure 1. Calibrated Rayleigh nights per month over the common window (grey v2.2, colour v2.3). CHM15k: the
repaired fallback tier accepts fewer nights every month. CL61: the summer months lose the most (the 5 km cap).
Mini-MPL: identical.*

## 4. Rayleigh calibration

### 4.1 What changed in the method

`eprof_v2.3` keeps tier 1 (the strict molecular-window search) bit-identical to `eprof_v2.2` and repairs tier 2,
the noise-aware fallback that v2.2 used on nights without a strict window: the CHM15k and CL61 use the noise measured
in the profile instead of the noise model, the CHM15k fallback must pass the vertical-consistency test (the
constant retrieved in two sub-windows has to agree within the budget), the CL61 fallback window may not end above
5 km AGL (`TIER2_MAX_END_M`), and the Payerne CHM15k A runs tier 1 only (`station_molecular_params`), because its
noisy receiver made its tier-2 nights -18 % low. The flags do not encode the tier: flag 0.5 marks a night that was
only partially clear (`filter_cloudy_profiles`), on either tier, and flag 1 a fully clear one, on either tier;
low-sun nights get flag 0.25; the NetCDF `calibration_version` is 230 (231 for low-sun nights). Neither the CSV nor
the NetCDF records which tier produced a night, so the comparison recovers it from the constants (identical =
strict tier, bit-identical by construction) and from reruns with the fallback switched off (appendix B).

### 4.2 Counts and the moved nights

| type | v2.2 successes | v2.3 successes | identical | tier-2 nights in both, ratio v2.3/v2.2 | dropped | added | low-sun |
|---|---|---|---|---|---|---|---|
| CHM15k | 27 640 | 25 375 (-8.2 %) | 22 828 | 2 063, median 0.973 (P10-P90 0.890-1.021) | 2 749 | 330 | 154 |
| CL61 | 842 | 707 (-16.0 %) | 653 | 31, median 0.999 (0.958-1.049) | 158 | 17 | 6 |
| Mini-MPL | 577 | 577 | 577 | 0 | 0 | 0 | 0 |

The fallback tier produced 17 % of the v2.2 CHM15k successes (4 812 of 27 640) and 22 % of the CL61 ones (189 of
842); v2.3 keeps 2 063 and 31 of them, moved, and adds 330 and 17. The CHM15k nights that the vertical-consistency
test removes were the tier-2 rescues the q19 study identified as the low-biased ones; the nights kept in both releases move down by 2.7 % on median, the direction the twin and the CL61
cloud references predicted. The CL61 common nights do not move (the measured noise reproduces the model where the
window is unchanged); what changes is which nights exist.

Payerne CHM15k A: 162 successes in v2.2, 44 in v2.3 (tier 1 only). The unit is the reference of the dashboards and
of the paper; its summer record is thin (about four nights per summer month) and its constants are now the clean
ones.

![tier-2 ratio, CL61 window tops, uncertainty](figs_v23_release/fig_rayleigh_changes.png)

*Figure 2. Left: the ratio v2.3 / v2.2 on the tier-2 nights present in both releases. Middle: CL61 molecular-window
tops (ASL); the cap acts on the fallback windows only, so tops above 5.5 km remain from tier 1. Right: relative
2-sigma uncertainty of the nightly constant on identical nights, v2.2 ensemble spread (dashed) against the v2.3
budget (solid), as densities.*

### 4.3 CL61 per stream (common window)

| stream | site | identical | changed | dropped | added | low-sun | Kalman v2.3/v2.2 (median, P5-P95) |
|---|---|---|---|---|---|---|---|
| 0-20000-0-03808_C | Camborne | 47 | 3 | 7 | 0 | 0 | 1.014 (0.980-1.215) |
| 0-20000-0-06418_B | Zeebrugge | 17 | 2 | 16 | 3 | 0 | 1.047 (0.911-1.202) |
| 0-20000-0-06447_B | Uccle | 8 | 5 | 17 | 0 | 0 | 1.005 (0.965-1.175) |
| 0-20000-0-06610_C | Payerne | 15 | 2 | 26 | 7 | 0 | 1.114 (0.985-1.295) |
| 0-20000-0-10393_C | Lindenberg | 163 | 1 | 23 | 1 | 0 | 0.998 (0.945-1.090) |
| 0-20000-0-11538_B | Temelin | 43 | 3 | 8 | 3 | 0 | 0.996 (0.926-1.037) |
| 0-20000-0-14015_B | Ljubljana | 18 | 1 | 13 | 1 | 0 | 1.033 (0.922-1.200) |
| 0-20008-0-BIR_A | Birkenes | 16 | 1 | 5 | 0 | 6 | 0.998 (0.850-1.003) |
| 0-203-10-LNG_A | Lanzhot | 134 | 7 | 19 | 1 | 0 | 0.997 (0.901-1.032) |
| 0-380-5-1_B | Aosta | 134 | 3 | 8 | 0 | 0 | 0.995 (0.934-1.014) |
| 0-756-4-EERLCL61_A | Sion | 58 | 3 | 16 | 1 | 0 | 1.004 (0.967-1.046) |

The Payerne CL61 loses 26 of its 43 common nights, most of them summer nights with 5-7 km windows, whose dark-free
constants were 15-30 % low (the hood-dark leverage, `eta_payerne_gap_2026-11-02.md`); its Kalman rises by 11 %
towards the cloud constant. Zeebrugge and Ljubljana move by +3 to +5 % for the same reason with fewer nights. The
CL61 Rayleigh constant stays a diagnostic in the network (dark-free); the cloud constant is the CL61 product
(policy item 5 of the plan).

### 4.4 Flags over the common window (Rayleigh rows)

| flag | meaning | v2.2 | v2.3 |
|---|---|---|---|
| 1 | successful (tier 1) | 11 691 | 10 551 |
| 0.5 | partial success (tier 2) | 17 368 | 15 948 |
| 0.25 | low-sun success (polar summer) | 0 | 160 |
| 0 | no data | 2 211 | 439 |
| -1 | unsuitable conditions | 53 632 | 55 209 |
| -2 | signal not proportional to molecular | 2 234 | 5 364 |
| -3 | method disagreement | 2 480 | 1 939 |
| -4 | missing CAMS | 209 | 209 |
| -6 | uncertainty exceeds value | 6 | 17 |
| -9 | another layer with lower signal | 24 | 20 |
| -99 | exception | 2 | 2 |

The fallback nights that v2.3 rejects land in -2 (the vertical-consistency failure is reported as a non-molecular
signal) and -1; -3 falls because the measured noise removes method disagreements that the model noise created. The
0 rows drop from 2 211 to 439 because the scratch L1 archive was completed between the two runs. Flags 1 and 0.5
both fall because fallback nights carry either flag (a partially clear night is 0.5 on both tiers).

## 5. Polar-summer low-sun nights

| station | type | lat | low-sun nights | dark nights | months | level, low-sun / dark |
|---|---|---|---|---|---|---|
| Hopen | CHM15k | 76.5 | 19 | 29 | Apr-Jul 2025, Apr-Aug 2026 | 1.29 |
| Bjornoya | CHM15k | 74.5 | 6 | 26 | Jul-Aug 2025, Jul-Aug 2026 | 1.25 |
| Jan Mayen | CHM15k | 70.9 | 15 | 25 | Apr-Aug 2025, Apr-Jun 2026 | 1.40 |
| Andoya | CHM15k | 69.3 | 36 | 76 | Apr-Jul 2025, Apr-Aug 2026 | 1.19 |
| Torshavn | CHM15k | 62.0 | 15 | 25 | May-Jun 2025, Jun-Jul 2026 | 1.18 |
| Flesland | CHM15k | 60.3 | 17 | 108 | May-Jul 2025, May-Jul 2026 | 1.26 |
| Lerwick | CHM15k | 60.1 | 3 | 40 | Jul 2025 | 1.52 |
| Oslo | CHM15k | 59.9 | 39 | 151 | May-Jul 2025, May-Jul 2026 | 2.63 |
| Birkenes | CL61 | 58.4 | 6 | 26 | Jun-Jul 2026 | 1.00 |
| Stornoway | CHM15k | 58.2 | 5 | 68 | Jun-Jul 2025 | 1.29 |
| South Uist | CHM15k | 57.4 | 1 | 67 | Jun 2025 | 3.03 |

The fallback does what it was built for: the Arctic CHM15k now have a summer record. The level ratios of 1.2 to 1.5
are expected for a CHM15k, whose constant peaks near the solstice (+17 % annual gain, plus the temperature
coefficient) while its dark nights at these latitudes are all winter nights; the CL61 at Birkenes, which has no such
cycle, gives 1.00. Oslo (2.63 on 39 nights) and South Uist (3.03 on one night) are outside that range and must be
inspected (residual solar background, window choice) before their low-sun constants are used; the flag 0.25 and
version 231 make them easy to exclude downstream (`polar_summer.md` section 6).

## 6. Uncertainty budget

On the identical tier-1 nights the 2-sigma relative uncertainty of the nightly constant is:

| type | nights | v2.2 ensemble spread, median (P10-P90) | v2.3 budget, median (P10-P90) |
|---|---|---|---|
| CHM15k | 22 828 | 4.8 % (2.3-10.0) | 14.2 % (8.9-26.6) |
| CL61 | 653 | 5.8 % (2.4-11.6) | 13.8 % (8.2-28.9) |
| Mini-MPL | 577 | 4.8 % (1.4-11.2) | 6.7 % (3.9-13.5) |

Components written to the NetCDF (`lidar_constant_uncertainty_<component>`, 2 sigma; all Rayleigh nights of the
v2.3 tree, median relative value):

| type | nights | total | ensemble | vertical | molecular | aerosol | noise |
|---|---|---|---|---|---|---|---|
| CHM15k | 28 652 | 14.7 % | 5.2 | 10.5 | 3.7 | 4.0 | 3.8 |
| CL61 | 949 | 13.7 % | 5.7 | 9.7 | 3.6 | 2.5 | 4.0 |
| Mini-MPL | 665 | 6.8 % | 4.9 | 1.2 | 3.6 | 0.5 | 0.2 |

The vertical-consistency term is the budget: it carries the height dependence of the constant inside the window
(overlap, residual aerosol, dark), which the ensemble spread never saw. The relative QC (flag -6, uncertainty above
the value) rejects 17 nights instead of 6. The budget is what the Kalman weights (section 7) and what the L2
proposal carries (`uncertainty.md` section 5).

## 7. Kalman best estimate

Per stream, median ratio v2.3 / v2.2 on common dates, and the posterior relative standard deviation (1 sigma):

| type, method | streams | ratio, median of stream medians (P10-P90) | Kalman std v2.2 | v2.3 |
|---|---|---|---|---|
| CHM15k Rayleigh | 153 | 0.995 (0.989-1.002) | 9.0 % | 9.2 % |
| CL61 Rayleigh | 11 | 1.004 (0.996-1.047) | 9.1 % | 8.9 % |
| Mini-MPL Rayleigh | 5 | 1.001 (1.000-1.007) | 11.4 % | 11.8 % |
| CL31 cloud | 210 | 0.982 (0.980-0.985) | 5.6 % | 5.7 % |
| CL51 cloud | 51 | 0.985 (0.982-0.988) | 5.2 % | 5.5 % |
| CL61 cloud | 11 | 0.984 (0.980-0.990) | 8.4 % | 8.8 % |

The weighting by the budget and the cleaner tier-2 record move the CHM15k best estimate by half a percent; the
cloud series move by the eta-table factor and nothing else. The posterior spread is unchanged: the budget
re-weights nights against each other but the series are long enough that the process noise sets the width.

## 8. Cloud calibration: eta tables, generation v2

### 8.1 The constant

| type | common days (v2.2 / v2.3) | ratio v2.3 / v2.2, median (P5-P95) | expected from the tables at CBH 0.5 / 1 / 2 km |
|---|---|---|---|
| CL31 | 58 325 (58 553 / 58 350) | 0.9808 (0.974-0.989) | 0.990 / 0.986 / 0.977 |
| CL51 | 14 725 (14 782 / 14 731) | 0.9836 (0.976-1.001); without the two retyped streams 0.9832 (0.976-0.991) | 0.991 / 0.984 / 0.976 |
| CL61 | 792 (795 / 792) | 0.9837 (0.977-0.992) | 0.991 / 0.985 / 0.977 |

The shift is the table ratio at the cloud bases that dominate the calibration scenes (1 to 1.5 km), uniform across
types, with no residual structure: the calibration code around the tables is unchanged.

![cloud ratio, Kalman ratios, low-sun nights](figs_v23_release/fig_cloud_kalman_lowsun.png)

*Figure 3. Left: cloud constant per day, v2.3 / v2.2, by type (densities); the CL51 tail at 1.03 is the two
retyped streams of section 8.3. Middle: per-stream median of the Kalman ratio. Right: polar-summer low-sun nights
per station with their level against the dark nights.*

### 8.2 Days that change state

Of the days calibrated in v2.2, 288 are rejected in v2.3 and 31 new ones are accepted (net -0.35 % of 74 000 days,
spread evenly over the months). The lost days go to the aerosol-below-cloud filter (-24: 244 days), cloud base out
of range (-25: 20), peak shape (-22/-23: 20), neighbours (-26: 2), uncertainty (-6: 2). The mechanism is the
constant itself: the filter compares the calibrated aerosol backscatter below the cloud with a fixed limit, and the
1.6-1.9 % lower constant lifts every profile's calibrated aerosol load by the same amount, so the marginal days cross
the limit. No stream loses more than 7 days.

### 8.3 Streams whose instrument changed inside the window

The comparison flagged two CL51 streams (Tusimice `0-20000-0-11438_A`, Libus `0-20000-0-11520_A`) whose cloud ratio is
1.03-1.04 instead of 0.98 on every day: v2.2 calibrated them as CL31, v2.3 as CL51 (the ratio is exactly
eta_v2(CL51) / eta_2026-07(CL31)). The census refresh of 2026-10 had retyped them from the newest L1 file. A scan of
the first L1 file of every month of every stream (`validation/v23_release/scan_type_changes.py`, 10 571
stream-months) and of the daily files around the changes shows a real replacement of the instrument under the same
stream identifier at three Czech stations: the range extent jumps from 770 gates (7.7 km, CL31) to 1 540 gates
(15.4 km, CL51), and the L1 `instrument_type` attribute follows one to twelve weeks later:

| stream | site | new range extent from | L1 type attribute CL51 from | v2.2 typed | v2.3 typed |
|---|---|---|---|---|---|
| 0-20000-0-11450_A | Plzen | 2025-07-16 | 2025-10-09 | CL51 | CL51 |
| 0-20000-0-11438_A | Tusimice | 2026-06-25 | 2026-07-03 | CL31 | CL51 |
| 0-20000-0-11520_A | Libus | 2026-06-26 | 2026-07-03 | CL31 | CL51 |

No other stream changes type or range extent inside the window. Consequences: in v2.3 the CL31 part of these records
(Plzen 2025-01 to 2025-07, the other two 2025-01 to 2026-06) is corrected with the CL51 eta table, 3-4 % too high;
in v2.2 the CL51 part of Tusimice and Libus was 3-4 % too low. The step in C_L at the swap (a different unit) is
smoothed by the Kalman in both releases. The three events are now in `monitoring/instrument_events.json` (dashboard
markers). The fix is a per-day instrument type in the runner, taken from the L1 range extent rather than from the
census or the lagging attribute, with a Kalman restart at the swap (section 10, item 1); the three streams are then
recomputed.

## 9. NetCDF and CSV content

| | v2.2 tree | v2.3 tree |
|---|---|---|
| files | 843 | 846 |
| `calibration_version` | 220 (28 531 rows), 1000 (74 130) | 230 (30 105), 231 (161), 1100 (85 095) |
| uncertainty components | none | five, in the 335 files that hold Rayleigh rows (CHM15k, CL61, Mini-MPL); the 511 cloud-only files (CL31, CL51) have none |
| `lidar_constant_height_gradient` | no | yes (Rayleigh rows) |
| Kalman CSV | `kalman`, `kalman_std` | `kalman`, `kalman_std` (weighted) |

One inconsistency of the CSV side: the `version` column of `<key>_cal.csv` is filled for Rayleigh rows (230 / 231)
but empty for cloud rows, although the NetCDF carries 1100; the runner now writes it (change made after the release
run, so the release CSVs keep the empty column; the NetCDF is the product).

## 10. Known limits and follow-ups

1. **Instrument swaps inside a stream** (section 8.3): per-day type from the L1 range extent in the runner, Kalman
   restart at unit changes, recompute of the three Czech streams. Code change, one day; until then the three
   streams carry a 3-4 % bias on their CL31 part.
2. **CL61 Rayleigh** stays dark-free and capped; it is a diagnostic in the network, the cloud constant is the
   product. The hub-side policy (which of our constants the operational coefficients carry for the CL61) is to be
   written into `extract_l2_opcoeff.py` at the cutover.
3. **Low-sun outliers**: Oslo (2.6 times the dark level, 39 nights) and South Uist (one night at 3.0) to be inspected.
4. **Autumn decline** of both Payerne Rayleigh constants (Aug to Oct 2026, -22 %, cloud flat): still open, watched
   on the dashboard.
5. **L2 uncertainty**: the hub applies a fixed 25 %; the budget and the Kalman standard deviation are now available
   to replace it (`uncertainty.md` section 5).
6. **Lindenberg CL61** ends 2026-06-15 in the scratch L1; the paper and the inter-comparison use the local copy.
7. **The incident rule**: the code staging rsync must include `calibration/data/*.nc` (`alc_v23_release_array.sbatch`
   header). A startup self-check of the data files in the runner would make this impossible to repeat.
8. **CSV cloud version** (section 9): fixed in the runner for the next runs.
9. **The tier is not recorded**: nothing in the CSV message or the NetCDF says whether a night came from the strict
   or the fallback tier (flag 0.5 means a partially clear night). A tier token in the message, or a
   `calibration_tier` variable, would make the next comparison and the dashboard filters direct.

## 11. Cutover checklist (phases 2-4 of the plan)

- [ ] add-ons array 339770 finished (`_omb_cache.npz`, `_sens_cache.npz` per stream), logs without traceback
- [ ] archive check: adapt `scripts/check_v22_archive.py` to a v2.3 version (expected codes 230 / 231 / 1100, flag
      0.25, the five components, no -99 on the Vaisala streams, the three Czech streams listed as known) and run it on
      the scratch tree
- [ ] copy the complete tree down (scratch is purged): CSV, NC, TXT, JSON and the two cache files per stream
- [ ] transfer to zueub434 as `/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.3` (rsync through WSL, the
      transfer section of `OPERATIONS_v22_cutover.md`)
- [ ] `ops/config.sh`: output tree, `ALC_MOLECULAR_METHOD=eprof_v2.3`, `ALC_UNCERTAINTY_MODEL=v2.3`,
      `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`, `ALC_POLAR_NIGHT_FALLBACK=1`, `ALC_OPCOEFF_CSV`; the deployed
      `options.json` with the polar switch and the Payerne override; gap fill from 2026-11-01 to yesterday
- [ ] payloads, a FULL dashboard build (no `--changed-only`), publish, the J+1 and J+7 checks; rollback = re-comment
      the block and rebuild
- [ ] `main` gets the operational snapshot (never a wholesale merge)
- [ ] `extract_l2_opcoeff.py` with the CL61 policy; the uncertainty proposal to the hub
- [ ] inter-comparison dashboard (cloud streams of Payerne and Lindenberg and the CL31 ladder on the v2 tables), the
      paper v1.0 (section 6 from these constants), the Confluence pages, and a note to the E-PROFILE users: Vaisala
      beta_att up by 1.6-1.9 % against the 2026-07 tables, Rayleigh constants unchanged on tier-1 nights, fewer
      fallback nights (-8 % CHM15k, -16 % CL61), low-sun nights in the Arctic, five new uncertainty variables, new
      version codes 230 / 231 / 1100

## Appendix A: files

* Report figures: `doc/reports/figs_v23_release/fig_availability_month.png`, `fig_rayleigh_changes.png`,
  `fig_cloud_kalman_lowsun.png` (regenerated by the comparison script).
* Tables: `C:/DATA/Projects/202611_v23_release/compare/` (`summary.txt`, `rayleigh_nights.csv`,
  `rayleigh_nights_per_stream.csv`, `lowsun_nights.csv`, `kalman_ratio_per_stream.csv`, `cloud_days.csv`),
  `type_changes.csv`.
* Scripts: `validation/v23_release/compare_v22_v23.py`, `scan_type_changes.py`; `ops/cscs/alc_v23_*.sbatch`.
* Release trees: `/scratch/mch/mhrvo/E_PROFILE_calout_v23_rel` (CSCS), local copy
  `C:/DATA/Projects/202611_v23_release/calout_v23_rel`.

## Appendix B: verification that only the fallback tier moved (local reruns, 2026-11-03)

The comparison cannot read the tier from the archives (section 4.1), so the invariance of the strict tier was
checked by rerunning nights locally (same L1 files as on CSCS, md5-identical; CAMS `D:/CAMS_Monthly_04`):

* the code of the v2.2 release (git worktree at `c2cbc6c`, 2026-08-17, `eprof_v2.2`) reproduces the v2.2 release
  constants bit for bit on Payerne CHM15k A (2025-04-01, 2025-06-12, 2026-03-02) and Deuselbach CHM15k
  (`0-20000-0-10615_0`: 2025-02-05, 2025-03-06, 2025-06-12);
* the v2.3 release code (`46bb407`) under `eprof_v2.2` with the fallback on reproduces them bit for bit as well
  (every commit between the two releases was bisected on the Payerne night: none changes the strict tier; the
  only change on Payerne is the intended tier-1-only override of `0025ced`);
* with the fallback switched off (`ALC_MOLECULAR_PARAMS='{"use_noise": false}'`, `eprof_v2.2`) the same code
  rejects all six nights with -2 "No molecular window passed the validity gates": they were fallback nights, with
  flag 1 because the nights were fully clear;
* under the release configuration the Deuselbach nights become 2025-02-05 moved (window 3 657-7 493 m to
  5 096-6 055 m, -7.9 %), 2025-03-06 moved (-5.0 %) and 2025-06-12 rejected, the v2.3 behaviour of a fallback
  night, and the Payerne nights are rejected by the tier-1-only override.

Scratch outputs of these runs are not kept; the commands are the runner calls of section 2 with
`--start <night> --end <night>` and the environment switches named above.
