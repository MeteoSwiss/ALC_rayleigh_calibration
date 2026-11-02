# Multiple-scattering correction of the Vaisala cloud calibration, v2 (draft)

Draft of 2026-11-02. Companion to `04_multiple_scattering.md` (v1), which is unchanged. Nothing here is
applied to the operational code (`calibration/cloud/_filters.py` is untouched). Every number below is
read from a CSV in `C:/DATA/Projects/202611_eta_study/` (path given with each table); if a run is
repeated, take the numbers from the CSV, not from this text. The CL61 water-vapour (WV) line is the
operational 910.55 nm / FWHM 0.188 nm throughout; the few inputs that came from older 910.74 nm runs are
flagged where they appear.

## 1. Purpose and what v1 assumed

The O'Connor cloud calibration needs the multiple-scattering factor eta(range). v1 replaced the legacy
Hopkin ladders by PVC tables (Hogan 2006) computed for each Vaisala type with two fixed cloud parameters:

| parameter | v1 value | where it came from |
|---|---|---|
| droplet equivalent-area radius a_G | 5.5 um | "calibration-scene size" from an earlier Cloudnet comparison |
| in-cloud extinction alpha | 10 /km | lower end of a fit scan against the legacy tables |
| receiver FOV / divergence | Wiegner et al. 2014 Table 1, Vaisala manuals | see section 3.3 |
| liquid lidar ratio S | 18.8 sr | O'Connor et al. 2004 |

v2 asks whether these inputs are right, measures them where possible, and proposes new tables.

## 2. Principle

**Every input is fixed from measurement, physics or manufacturer specification. The CBH slope of the
cloud lidar constant C_L and the agreement with the Rayleigh constant are independent checks only; they
are never used to choose or tune an input.** (Project decision, 2026-11-02: a flat C_L is a diagnostic,
not a target for finding the inputs.) Section 6 therefore gives sensitivities, not fits.

## 3. Inputs

### 3.1 Droplet size a_G: Cloudnet, lowest 100 m above the liquid base

Cloudnet `der` product (Payerne, Palaiseau, Lindenberg; 2025-01..2026-06), matched to every accepted
calibration profile (nearest 30-s Cloudnet profile within 6 min, liquid base = lowest liquid gate within
+-150 m of the ceilometer CBH, followed down through contiguous liquid). Only pixels with
`der_retrieval_status == 1` (reliable) are used; status 4 ("surrounding ice") never carries a radius.
r_e from `der_scaled` (LWP-scaled), a_G = r_e exp(-1.5 sigma^2) = 0.832 r_e (lognormal, sigma 0.35).

Source: `step1d/summary_msfov.csv` (scene list of the July ms_fov run, see flag below).

| site / instrument | profiles with status-1 data in lowest 100 m (days) | a_G median [P10-P90] (um) | alpha median [P10-P90] (/km) |
|---|---|---|---|
| Payerne CL61 | 2142 (41) | 3.48 [1.7-5.0] | 32.8 [12.1-241] |
| Payerne CL31 (distorted unit) | 17470 (210) | 3.60 [2.1-5.5] | 26.6 [7.3-163] |
| Palaiseau CL31 | 10469 (227) | 3.10 [2.1-4.9] | 35.4 [11.4-97.6] |
| Lindenberg CL61 (>= 2025-05) | 12447 (192) | 3.04 [1.9-5.5] | 12.8 [2.2-41.2] |
| all Vaisala | 42528 (386) | 3.31 [2.0-5.3] | 24.0 [5.2-108] |

- **a_G = 3.3 um** is adopted (all-Vaisala median, lowest 100 m). It grows with depth: 3.0 / 3.3 / 3.6 um
  in the lowest 50 / 100 / 150 m (all Vaisala). Uncertainty used below: 2.5-5 um.
- Only 30-46 % of matched profiles have status-1 pixels in the lowest 100 m.
- Cloudnet alpha = 2 pi N_scaled r_e^2 exp(-3 sigma^2) is reported for comparison only (section 3.2): its
  P90 is driven by unphysical N_scaled (> 5000 cm-3 at Payerne).
- *910.74 flag:* the scene list comes from the July ms_fov run (910.74 nm, 1 deg CAMS); only the CL61
  profile selection can depend on it. Re-matched on the operational 910.55 nm dump
  (`step1d/summary_dump.md`) every median agrees within 2 %.

![Cloudnet microphysics at the calibration scenes](figs_04v2/fig_microphysics_msfov.png)

### 3.2 Extinction alpha: the ceilometer's own tail decay (alpha_k)

**Method.** For every accepted calibration profile, today's operational calibration was re-run and the
attenuated backscatter stored on the gates from 60 m below to 300 m above the peak
(`validation/eta_study/extract_incloud_profiles.py`; reproduces the operational per-profile C exactly).
alpha_k is the log-slope of beta / beta_peak on +20..+150 m above the peak (gates above 0.5 % of the peak),
inverted with the PVC forward model computed for the same gates and estimator (k_app ~ 1.87-1.92 alpha).

**Why the tail and not the full shape.** A full-shape fit of the PVC model does not identify alpha: the
observed peak is broader than PVC plus the nominal pulse (110 ns CL31/CL51, 160 ns CL61) plus gate
integration. The fit prefers an effective response of about 200-250 ns, and the fitted alpha follows the
assumed response (`step2/response_width_test.csv`). A Gaussian response preserves an exponential slope,
so alpha_k stays within +-1.5 % for any response from 110 to 400 ns.

**Robustness** (`step2b/`):

| test | result | source |
|---|---|---|
| fit window +20..+80 / +20..+150 / +40..+150 m | alpha_k changes by +-5-10 % | `window_sensitivity_by_unit.csv` |
| noise floor | SNR at +150 m is 0.3-1.1 for CL31/CL51 (6 for CL61): the slope effectively uses +20..+100 m | idem |
| causal-receiver test | across alpha_k quintiles the 10-90 % rise below the peak stays at ~29-32 m while the tail k triples -> tail = atmosphere, rise = instrument | `quintile_rise_tail.csv` |
| pulse x causal receiver x PVC fit | receiver time constant and alpha are not separable from the profile shape | `receiver_fit.csv` |
| CL61 native 4.8 m vs 9.6 m working gates | alpha ratio 1.00 (P10-P90 0.86-1.11) | `step2/cl61_native_vs_working.csv` |
| 20 m-gate CL31 (41 units, real L1 configuration) | consistent within ~2 /km | `network/summary.txt` |

**Network values** (270 units, 2025-01..2026-06; `network/unit_summary_all.csv`):

| type | units | alpha_k, median of unit medians | P10-P90 of unit medians |
|---|---|---|---|
| CL31 | 211 | 27.8 /km | 24.0-31.3 |
| CL51 | 48 | 30.2 /km | 26.1-31.7 |
| CL61 | 11 | 31.0 /km | 27.2-32.2 |

**Regional pattern** (`network/summary.txt`): with the type median removed, subarctic maritime units sit at
-5.6 /km, northern Europe -3.4, NW Atlantic -1.2, central and western Europe +0.7..+0.9. Region explains
37 % of the between-unit variance, firmware 15 % (19 levels), gate spacing 3 %. Within one firmware (CL31
fw 205, 62 units) the regional pattern persists; within central + western Europe (126 CL31) firmware
explains 2 %. alpha_k is higher in summer than in winter (+7..+13 % in August). This behaves like cloud
microphysics (alpha ~ N^(1/3) LWC^(2/3)), not like an instrument setting.

**Cloudnet comparison** (`step2/closure_cloudnet.csv`, +-15 min, status 1, lowest 100 m): median ratio
alpha_k / Cloudnet alpha = 0.92 (Payerne CL31), 0.86 (Payerne CL61), 0.85 (Palaiseau CL31) and 2.51
(Lindenberg CL61, the only 35 GHz radar site). Profile and daily correlations are weak (Spearman <= 0.2).
The two agree on the typical level at the two 94 GHz sites and share no variability.

![Mean normalised profile vs PVC](figs_04v2/fig_a_mean_profile_vs_lut.png)
![alpha_k robustness](figs_04v2/fig_alpha_k_robustness.png)
![Network alpha_k](figs_04v2/fig_network_alpha_k.png)
![Closure with Cloudnet](figs_04v2/fig_c_closure_cloudnet.png)

### 3.3 Receiver FOV and laser divergence

| type | receiver FOV rho_t | divergence rho_l | source and convention |
|---|---|---|---|
| CL31 | 0.83 mrad | +-0.4 x +-0.7 -> 0.57 mrad (RMS) | Wiegner et al. 2014 Table 1: "RFOV the receiver field of view (half angle) ... from user manuals"; divergence from Vaisala |
| CL51 | 0.56 mrad | +-0.15 x +-0.25 -> 0.21 mrad | idem |
| CL61 | 0.56 mrad | +-0.2 x +-0.35 -> 0.28 mrad | CL61 User Guide M212475EN; Le et al. 2026 Table 1 |

The code treats rho_t as a top-hat half-angle and rho_l as a 1/e Gaussian half-angle. No source says
whether the manufacturer "+-" value is an edge, a 1/e or a 1/e^2 width (or a full angle), no answer from
Vaisala is available, and no temperature dependence of the optics is documented. **Central value:** the
documented Wiegner et al. (2014) half-angle, used as a top-hat edge (rho_t) and the Vaisala divergence as a
1/e half-angle (rho_l). **Uncertainty:** the convention is carried explicitly in the eta band (section 4):
rho_t and rho_l are each drawn between / sqrt(2) and x sqrt(2) of the documented value. It is the largest
single term (section 6).

### 3.4 Liquid lidar ratio S: Mie at the measured droplet size

Mie (miepython; water n = 1.3275 - 4.86e-7 i at 910.55 nm, Hale & Querry 1973) for a lognormal
population (sigma 0.35), S = 4 pi <Q_ext r^2> / <Q_back r^2> (`step4/mie_lidar_ratio.csv`):

| r_e (um) | 3.5 | 3.97 (= a_G 3.3 / 0.832) | 4.0 | 4.5 | 5.5 | 8.0 | 10.0 |
|---|---|---|---|---|---|---|---|
| S (sr) | 17.46 | **18.35** | 18.41 | 19.21 | 19.87 | 19.18 | 18.86 |

The 10 um value (18.86 sr) reproduces O'Connor et al. (2004, 18.8 sr). v2 uses S = 18.35 sr; uncertainty
+-0.5 sr (the spread over r_e 3.5-4.5 um).

## 4. Proposed eta tables

PVC at a_G 3.3 um and the network alpha_k of each type, on the operational CBH nodes, with the documented
optics (Wiegner 2014 half-angle as top-hat edge). +-1 sigma band (16th / 84th percentiles) from 300
draws: a_G uniform 2.5-5 um; alpha_k uniform between the P10 and P90 of the type's unit medians; receiver
FOV rho_t and divergence rho_l each log-uniform between / sqrt(2) and x sqrt(2) of the documented value
(unknown edge / 1/e / 1/e^2 / full-angle convention); PVC model error 2 % on eta. Tables in
the `_ETA_<type>` format: `step3/eta_tables_proposed.py.txt` (`_PROPOSED`, `_LO`, `_HI`); values in
`step3/eta_tables_proposed.csv`; script `validation/eta_study/eta_tables_proposed.py`.

C_L scales with eta(CBH) x S. Change against the operational configuration (PVC 5.5 um / 10 /km,
S 18.8), from `step3/eta_change_vs_operational.csv` combined with S = 18.35 / 18.8:

| type | CBH | eta operational | eta proposed [-1 sigma, +1 sigma] | C_L change, eta only | C_L change, eta and S [band] |
|---|---|---|---|---|---|
| CL31 | 0.5 km | 0.894 | 0.856 [0.808, 0.875] | -4.2 % | -6.5 % [-11.8, -4.5] |
| CL31 | 1 km | 0.835 | 0.782 [0.727, 0.804] | -6.3 % | -8.5 % [-15.0, -6.0] |
| CL31 | 2 km | 0.759 | 0.699 [0.645, 0.719] | -7.9 % | -10.1 % [-17.0, -7.5] |
| CL51 | 0.5 km | 0.922 | 0.884 [0.846, 0.907] | -4.1 % | -6.4 % [-10.5, -4.0] |
| CL51 | 1 km | 0.873 | 0.816 [0.771, 0.844] | -6.5 % | -8.7 % [-13.8, -5.6] |
| CL51 | 2 km | 0.803 | 0.733 [0.686, 0.765] | -8.7 % | -10.9 % [-16.6, -7.1] |
| CL61 | 0.5 km | 0.924 | 0.884 [0.842, 0.908] | -4.3 % | -6.6 % [-11.0, -4.1] |
| CL61 | 1 km | 0.875 | 0.817 [0.766, 0.845] | -6.7 % | -8.9 % [-14.6, -5.8] |
| CL61 | 2 km | 0.807 | 0.735 [0.682, 0.764] | -9.0 % | -11.2 % [-17.5, -7.6] |

The band is asymmetric because the a_G range (2.5-5 um) is not centred on 3.3 um; the optics convention
widens it mostly towards lower eta (a wider FOV collects more multiply scattered light).

![Proposed eta tables](figs_04v2/fig_eta_tables_proposed.png)

## 5. Checks (not fits)

### 5.1 CBH slope, network-wide

Slope of ln C_L against cloud-base height. Per accepted profile, the C_L of each eta variant is obtained
from the stored in-cloud windows on the same profiles, C_L,v = C_L,oper x Int(beta/T2 eta_v) / Int(beta/T2
eta_oper); this matches exact re-runs of the calibration to a median 0.998 (P5-P95 +-1 %,
`slope/network/validation_windowed_vs_exact.txt`). Unit x month fixed effects remove the seasonal cycle;
CIs are clustered by unit-day. Source: `slope/network/cbh_slope_network.csv`. S does not change slopes.

| type (units, profiles) | operational | legacy Hopkin | proposed (a_G 3.3, alpha_k) |
|---|---|---|---|
| CL31 (211, 8.8 M) | +0.7 [0.6, 0.9] %/km | +7.3 | **-1.7 [-1.8, -1.5]** -> fails |
| CL51 (48, 1.8 M) | +2.8 [2.5, 3.2] (daily data +4.8) | +7.9 | **-0.3 [-0.6, +0.05]** -> passes |
| CL61 (11, 58 k) | +2.3 [0.9, 3.6] | +7.2 | **-0.9 [-2.3, +0.4]** -> passes |

Removing the distorted Payerne CL31 changes the CL31 rows by less than 0.1 %/km. The result does not
depend on the statistic: daily Theil-Sen, per-profile unit fixed effects and unit x month fixed effects
give the same ranking (`slope/network/cbh_slope_network.csv`).

### 5.2 Level against the CL61's own Rayleigh constant

Exact re-run of the operational cloud calibration per unit-day and eta variant
(`validation/eta_study/step4_cloud_eta_variants.py`; the operational variant reproduces the operational
910.55 nm cloud constants to 2e-16). Ratio of the median daily cloud C_L to the median daily Rayleigh C_L
(v2.2, 910.55 nm) over the common period, 7-day block bootstrap 95 % CI (`step4/ratio_cloud_rayleigh.csv`).
"Proposed" = eta (a_G 3.3, alpha_k) and Mie S.

| unit, reference | days cloud / Rayleigh | operational | legacy Hopkin | proposed eta | proposed eta and S |
|---|---|---|---|---|---|
| Payerne CL61, Rayleigh + measured dark | 61 / 49 | 1.053 [0.97, 1.13] | 0.986 [0.91, 1.06] | 0.967 [0.89, 1.05] | 0.944 [0.87, 1.02] |
| Payerne CL61, Rayleigh without dark | 59 / 29 | 1.230 [1.09, 1.35] | 1.155 | 1.129 | 1.102 [0.98, 1.22] |
| Lindenberg CL61, Rayleigh without dark | 225 / 142 | 1.037 [1.01, 1.07] | 0.951 [0.93, 0.99] | 0.969 [0.94, 1.01] | 0.946 [0.92, 0.98] |
| 9 other CL61 (no dark), median [unit range] | | 1.048 [0.95, 1.19] | 0.959 | 0.986 | 0.963 [0.86, 1.10] |

**Not decisive.** At Payerne the measured dark moves the Rayleigh reference by 17 %, more than any eta
choice. Lindenberg has no dark measurement. **The Lindenberg CL61 was swapped on 2025-04-28**: the daily
quick-look shows a sudden drop of the background at ~12 UTC that day (the ~x2 step previously described as
"at 2025-05"). The Lindenberg window used here starts on 2025-05-01 and therefore covers only the new unit;
the cloud C_L shows no further step later in 2025 (`step4/daily_cloud_variants.csv`). Le et al. (2026,
AMT 19, 4923) report hood measurements only for their Finnish units, with a dark residual of ~1-2 % of the
molecular signal (figure as quoted to us, not re-checked here); the laser-power decline they report for
Lindenberg (40 % to ~10 %, cloud factor down ~x3; report 09 section 2.7) concerns the **old** unit. The
unit-to-unit spread of the other CL61 (~0.86-1.19) is larger than the 7-9 % separating the variants.

![Step 4 time series](figs_04v2/fig_step4_timeseries.png)

## 6. Sensitivities (not fits)

One input perturbed at a time around the central values (a_G 3.3 um, network alpha_k, type optics,
operational WV line), all 270 units, same accepted profiles; the S row is analytic (C_L scales with S;
+-0.5 sr is +-2.7 % at S ~ 18.4-18.8 sr). Source:
`sensitivity/sensitivity_table.csv`; script `validation/eta_study/sensitivity_table.py`.

| input change | CL31 slope (%/km) / level at 1 km | CL51 | CL61 |
|---|---|---|---|
| a_G -0.5 / +0.5 um | +0.67 / -0.51; +2.5 / -2.2 % | +0.81 / -0.66; +2.2 / -2.0 % | +0.80 / -0.66; +2.2 / -1.9 % |
| alpha -5 / +5 /km | +0.84 / -0.61; +2.9 / -2.5 % | +0.92 / -0.74; +2.4 / -2.1 % | +0.88 / -0.71; +2.3 / -2.0 % |
| rho_t / sqrt(2), x sqrt(2) | +1.00 / -1.33; +3.4 / -5.4 % | +1.99 / -1.77; +4.9 / -5.5 % | +1.73 / -1.87; +4.2 / -5.5 % |
| rho_l / sqrt(2), x sqrt(2) | -0.15 / -0.34; +0.1 / -1.9 % | -0.20 / +0.30; -0.4 / +0.5 % | -0.31 / +0.05; -0.5 / -0.2 % |
| CL31 WV line centre -1 / +1 nm | -0.32 / +0.19; -0.6 / +0.4 % | | |
| CL31 WV line FWHM 3 / 9 nm | +0.47 / -0.30; +0.9 / -0.6 % | | |
| S -0.5 / +0.5 sr | 0; -2.7 / +2.7 % | 0; -2.7 / +2.7 % | 0; -2.7 / +2.7 % |

The receiver-FOV convention is the largest term, followed by alpha and a_G; the CL31 WV line is small.
WV caveat: humidity is not stored with the profiles, so a line change enters as T2 -> T2^R, with R taken
from the HITRAN table on a standard humidity profile at 1 km (R spans 0.92-1.14 over 908.7-910.7 nm and
2-8 nm FWHM, `slope_null/cl31_wv_line_to_r.csv`).

## 7. Open issues

| issue | why it matters | what would settle it |
|---|---|---|
| Vaisala FOV / divergence convention (edge, 1/e, 1/e^2; half or full angle) | largest sensitivity: +-1-2 %/km and up to 5.5 % in level; no answer available | central = documented half-angle as top-hat edge; carried in the eta band (/ or x sqrt(2)); a written answer from Vaisala would remove it |
| CL31 residual slope (-1.7 %/km with the proposed eta) | the CL31 check fails while CL51 and CL61 pass | the FOV answer above; a CL31-specific range term (an earlier report found the flat CL31 median to be a cancellation, `cbh_slope_heterogeneity.md`) |
| Lindenberg CL61 dark (new unit since the 2025-04-28 swap) | no dark reference for the unit in the window; Le et al.'s laser decline concerns the old unit | hood measurement on the new unit |
| why the Payerne CL61 dark is worth 17 % of the Rayleigh constant | it dominates the level check | remote-dark reports 14-17, `dark_clearsky_method.md`; repeat hood |
| Payerne CL61 slope negative for every variant (-3 to -12 %/km, 61 days) | site-specific, not eta (`step4/cbh_slope_by_variant.csv`) | more days; near-range CL61 effects at Payerne |
| full in-cloud shape unusable (effective response 200-250 ns) | only the tail slope carries alpha | the Vaisala impulse response, or a hard-target test |

## 8. Bugs in the July Cloudnet scripts (Test 2 is void)

`validation/_ms_cloudnet_tests.py::match_one` (July 2026) computed the Cloudnet extinction with three
errors (sized in `step1d/summary_msfov.csv`):

| error | effect |
|---|---|
| factor 1.4 instead of the lognormal moment exp(-3 sigma^2) = 0.69 | alpha x2.02 |
| `der` (fixed-N retrieval) paired with `N_scaled` (LWP-scaled retrieval) | a further x1.0-1.4; inconsistent pair |
| `N_scaled` read with `np.asarray`, keeping the 9.97e36 fill value | in 32-59 % of scenes the median N was the fill value; Test 2 then silently replaced those alpha by the 10 /km default |

The new matcher (`validation/eta_study/cloudnet_microphysics_match.py`) uses `der_scaled` with
`N_scaled`, the correct moment, masked fill values and status-1 pixels only. Its replica of the old
formula reproduces the July match cache (same r_e and alpha medians, profile count within 0.3 %). **The July Test 2 closure result should not be used.**
It also used the 910.74 nm configuration.

## 9. Compute and reproducibility

Scripts: `validation/eta_study/` (not committed at the time of writing). Data and figures:
`C:/DATA/Projects/202611_eta_study/` (`REPORT.md` holds the working log). L1 `D:/E-PROFILE_L1_2026`,
CAMS 0.4 deg (`D:/CAMS_Monthly_04;D:/CAMS_daily`), Cloudnet `D:/Cloudnet_Droplet_radius`.

| step | script | where / time |
|---|---|---|
| Cloudnet microphysics | `cloudnet_microphysics_match.py` | PC, 29 s (12 processes) |
| PVC forward table | `forward_lut.py` | PC, 81-92 s (6 processes) |
| in-cloud extraction, 43 study units | `extract_incloud_profiles.py` | PC, 22 min (24 processes), 1.07 M profiles |
| in-cloud extraction, network | same + `run_network_pc.py`, `cscs/extract_fit_array.sbatch` | CSCS Linard postproc (CPU only) arrays 335747 / 335908 / 336141 (8 CPU, 12 GB per task, throttled to keep <= 6 of 14 nodes) for 219 units; PC for 51; split in `network/split.csv` |
| shape fit, alpha_k | `fit_alpha_shape.py`, `response_width_test.py`, `alpha_k_robustness.py` | PC, ~2 min each |
| per-unit summary | `collect_unit_summary.py`, `network_alpha_k_analysis.py` | Linard pp-short 336043 / 336154; PC |
| CBH-slope test | `cbh_slope_windowed.py` | Linard pp-short 336378 (55 s) + PC |
| level check | `step4_cloud_eta_variants.py`, `step4_cbh_slope.py` | PC, ~20 min (16 processes) |
| proposed tables | `eta_tables_proposed.py` | PC, 9 s |
| sensitivities | `sensitivity_table.py` | Linard pp-short 336463 + PC |

The Linard dry run (job 335708, pp-short) is bit-identical to the PC on every per-profile column. The
CL61 WV line is checked at start-up by every script (`(910.55, 0.188)`). On Linard the code was staged
under `/scratch/mch/mhrvo/eta_study/code` (the Linard `$HOME` is over quota).

### 5.3 Payerne gap and per-unit check (2026-11-02)

See `eta_payerne_gap_2026-11-02.md`: the +7.8 pct of section 5.2 context was the v1 hood dark; with dark v3 the CL61 Rayleigh+dark and the CHM15k agree +1 +-3 pct Jun-Oct 2026, and the proposed eta does not improve the cloud/Rayleigh agreement per unit.

### 5.5 CBH binscatters per unit and the switch decision (2026-11-02)

See `eta_switch_decision_2026-11-02.md`: per type the proposed table flattens the CL51 (+2.8 -> -0.3 pct/km) and the CL61 (+2.3 -> -0.9) and over-corrects the CL31 (+0.7 -> -1.7); per unit the residual shapes (+-3-4 pct/km) exceed the difference between tables (units above 2 pct/km: 174 -> 145 of 270). Recommendation: adopt the proposed inputs as a versioned table with a full rerun (section 3 there lists the reruns, dashboards and paper work).
