# The Payerne "8 %" gap between the CL61 Rayleigh + dark and the CHM15k, and the per-unit eta check (2026-11-02)

Companion to `04_multiple_scattering_v2.md` (sections 5.2 and 7, open item 1 of `eta_study_HANDOFF.md`).
Scripts: `validation/eta_study/payerne_extension.py` (the inter-comparison dashboard arithmetic replayed on an
L1-only hourly cache, 2026-06-01..10-31, 1136 paired hours; `--dark v1|v3`) and `payerne_gap_figures.py`.
Data: `C:/DATA/Projects/202611_eta_study/ext/` (`summary.txt`, `summary_darkv1.txt`, `monthly_gap*.csv`,
`nights_*.csv`, runs `c_nodark`, `a_v22`, `b_cloud`, cache `_streams_payerne_ext.npz`).

## 1. What the +7.8 % is

The +7.8 % of the handoff is the dashboard's default state: CL61 Rayleigh + measured dark (`l55s008dark`, hood
profile **v1**, 6 logbook sessions) against the CHM15k on OUR `v2.2dark` constants (same v1 dark), 2026-06-11..08-14,
678 h, manufacturer WV line on both sides. It is not a comparison with the distributed L2: against the L2 v1.0 the
same channel is −14.8 %, and the four CHM15k references span 23 points (L2 −14.8, v2.0 +2.4, v2.2 −0.4,
v2.2dark +7.8). The CHM15k TUB140016 is not an absolute anchor at Payerne in summer (per-night scatter ±40 %, the
strict v2.0 gates keep 2 nights in the window, its own hood dark is worth +11..+20 % of the constant).

It is a mixture of months: June +18.6 %, July +6.2 %, August −2.8 % (night +21.1 / +7.2 / −0.9). It is flat across
altitude bands apart from the known TUB140016 near-range overlap deficit (+10..+17 points in 0.5–1 km only) and
flat across the molecular fraction of the CL61 signal, i.e. a scale difference, not an Angstrom or WV signature.
Sensitivities: Angstrom 0.5–1.5 moves it ±5 points; the crossed WV state (constants at 910.55 nm, profile divided
by the 910.74/1.0 transmission) adds +23 points, the trap behind several older numbers.

## 2. Cause: the 6-session dark profile

Replaying the identical arithmetic with the 10-session, window-detected, P-smoothed hood profile **v3**
(`dark_profiles_payerne_v3.npz`, reruns `diag_v22_dark_v3` for BOTH instruments):

| CL61 Rayleigh + dark vs CHM15k v2.2 + dark, 500–3000 m, medrel % | Jun | Jul | Aug | Sep | Oct |
|---|---|---|---|---|---|
| hood dark v1 (dashboard state; Sep–Oct constants clamped, ignore) | +18.2 | +6.1 | +0.2 | — | — |
| **hood dark v3, all hours** | **+2.2** | **−0.3** | **−3.1** | **+2.5** | **+3.9** |
| hood dark v3, night only | +0.5 | +0.7 | −0.5 | +4.3 | +4.6 |

The v1 → v3 change raises the June CL61 constant (per-night medians 0.944 → 0.990, Kalman 0.88 → 0.97) and
lowers the CHM15k one (Kalman 7.0 → 6.5e11), both closing the gap. The dependence of the dark-corrected CL61 constant
on the molecular-window centre falls from −5.7 %/km (v1, n = 73) to −2.3 %/km [−5.2, +1.2] (v3, n = 100), against
−16.9 %/km [−21.6, −11.6] without any dark. **The two independent Rayleigh calibrations agree to +1 ± 3 % in every
month from June to October; there is no 8 % offset to explain.** The dashboard still runs on the v1 dark; rebuilding
it on v3 is a separate task.

![monthly gap](figs_04v2/fig_payerne_gap_monthly.png)
*Left/middle: monthly medrel of the CL61 channels against the CHM15k (v2.2 + dark v3, v2.2 no dark); right:
uncalibrated signal ratios (not visually checked before the session ended).*

![CL61 nights](figs_04v2/fig_payerne_cl61_nights.png)
*Per-night CL61 Rayleigh constants against the window centre (dark v3 / no dark), monthly medians of the CL61 and
CHM15k series.*

## 3. What the eta level check becomes

CL61 cloud (operational eta) vs CHM15k v2.2 + dark v3: −5.1 / −5.0 / −9.5 % (Jun / Jul / Aug); proposed eta
(× 0.90): +5.5 / +5.6 / +0.6 %. Cloud / Rayleigh+dark(v3) per month: Mar 1.03, Apr 0.96, May 0.99, Jun 1.06,
Jul 1.08, Aug 1.09 (operational) → 0.93–0.98 (proposed). Symmetric: the check does not discriminate at Payerne
(still "not decisive"), but it no longer favours the proposal either.

## 4. Per CL61 unit: the proposal does not improve the cloud / Rayleigh agreement

`step4/ratio_cloud_rayleigh.csv`, 12 unit-references. Operational → proposed: every unit moves by −9 %; 4 units
improve (Uccle, Temelin, Ljubljana, Zeebrugge), 4 worsen (Birkenes, Lanzhot, Aosta, Sion), Lindenberg, Camborne and
Payerne flip sign (+4..+5 % → −4..−6 %). Median |ln ratio| over the ten dark-free units 0.052 → 0.075, mean
0.074 → 0.086, units within 5 %: 4 → 4. The unit spread (0.95–1.19) exceeds the shift, and part of it is the dark of
the dark-free Rayleigh reference, unit-specific (opaque-low-cloud remote estimate at 3–6 km, `cl61_cloud_dark.log`):
the two strongly negative units are the two highest ratios (Payerne −25 %: 1.23 dark-free → 1.05 hood-corrected;
Uccle −20 %: 1.19), the positive-dark units are the lowest (Aosta +17, Lanzhot +23, Birkenes +60 %: 0.95). With the
dark accounted for, the operational table gives ratios ≈ 1.0–1.05 at the dark-characterised units and the proposal
0.90–0.95. The Payerne hood measurement removes the electronic-offset ambiguity only there; elsewhere the Rayleigh
reference carries an unknown ±5..20 % dark.

![units](figs_04v2/fig_cl61_cloud_rayleigh_units.png)

## 5. New open point: the autumn decline of both Rayleigh constants

From late August both Rayleigh constants fall together (CL61 1.064 → 0.948 → 0.822, CHM15k 7.8 → 6.5 → 6.1e11,
Aug / Sep / Oct), so their ratio holds, while BOTH cloud constants are flat (CL61 1.16 → 1.15 → 1.01 on 7 / 5 / 2
days; CL31 0.81 → 0.79 → 0.87e8). Housekeeping: CL61 window transmission 84 → 94 % and reported laser energy
96.5 → 90.5 % on 2026-07-28 (an intervention), laser 87–88 % and temp_int 33 → 21 °C by October; CHM15k temp_int
32 → 23 °C (its kappa of +1.3..+2.1 %/K predicts −12..−19 %). Common far-range causes (temperature-dependent darks on
cooler electronics with a fixed May–September profile; an attenuating layer the Klett term misses) are not
separated. The Sep–Oct cloud / Rayleigh ratio (1.21) must not be used for eta until this is understood; an autumn
hood session would settle the dark part.

## 6. Data added

Payerne L1 2026-10-09..31 synced from CSCS (`/scratch/mch/mhrvo/E-PROFILE_L1`) into
`D:/E-PROFILE_L1_2026/0-20000-0-06610/2026/10/` (69 files, 4.8 GB). CAMS daily exists only to 2026-10-08 (locally
and on CSCS), so CL61 nights after that are flagged −4 (no WV); the CHM15k calibrates.
