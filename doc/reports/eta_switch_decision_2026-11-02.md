# Switching the Vaisala cloud-calibration eta inputs: CBH check per unit, recommendation, and what a switch implies (2026-11-02)

Companion to `04_multiple_scattering_v2.md` (proposed tables, section 4; checks, section 5) and
`eta_payerne_gap_2026-11-02.md` (the Payerne level check, closed). Script
`validation/eta_study/cbh_binscatter_units.py`; per-unit table
`C:/DATA/Projects/202611_eta_study/slope/network/cbh_slope_per_unit.csv`.

The question: should the operational eta tables (`calibration/cloud/_filters.py`, `_ETA_CL31/_CL51/_CL61`:
PVC, a_G 5.5 um, alpha 10 /km, S 18.8 sr) be replaced by the study's inputs (PVC, a_G 3.3 um, alpha = the type's
in-cloud decay alpha_k 27.8 / 30.2 / 31.0 /km, Mie S 18.35 sr), and does that remove the cloud-base-height (CBH)
dependence of C_L unit by unit? The CBH flatness is a diagnostic, not a fitting target (operator rule): the inputs
were fixed from Cloudnet droplet sizes and the ceilometers' own in-cloud decay, never from these checks.

## 1. CBH dependence of C_L, per type and per unit

Daily medians of ln C_L per eta variant from the stored in-cloud windows of every accepted profile (270 Vaisala
units, 69 776 unit-days, `cbh_slope_windowed.py --stats`), centred on each unit's own median and binned in 0.25 km
of CBH. Variants: a = operational, b = legacy Hopkin, c = proposed.

![per type](figs_04v2/fig_cbh_binscatter_types.png)
*Pooled per type. The operational residual is a rise of C_L from 0.5 to about 1.5 km then a plateau; the proposed
table flattens it for the CL51 and the CL61 and bends the CL31 down above 1.5 km.*

| slope of ln C_L vs CBH | CL31 (211 units) | CL51 (48) | CL61 (11) |
|---|---|---|---|
| per profile, unit x month fixed effects: operational | +0.7 [0.6, 0.9] %/km | +2.8 [2.5, 3.1] | +2.3 [0.9, 3.6] |
| same, proposed | **−1.7 [−1.8, −1.5]** | **−0.3 [−0.6, 0.0]** | **−0.9 [−2.3, +0.4]** |
| same, legacy Hopkin | +7.3 | +7.9 | +7.2 |
| per unit, daily Theil-Sen, median: operational → proposed | +1.6 → −0.7 | +5.1 → +2.2 | +2.4 → −1.5 |
| per unit, IQR of the slope (proposed) | −3.2 .. +1.1 | +0.1 .. +3.4 | −2.4 .. +2.0 |
| units within ±2 %/km: operational → proposed | 86 → 103 | 6 → 18 | 4 → 4 |
| median |slope| per unit: operational → proposed | 2.3 → 2.1 | 5.1 → 2.6 | 2.8 → 2.4 |

(The daily statistic keeps the seasonal co-variation of CBH and C_L, hence the larger CL51 values; the fixed-effects
statistic removes it.)

![CL61 units](figs_04v2/fig_cbh_binscatter_CL61_units.png)
*One panel per CL61 unit. Lindenberg (225 days) and Aosta keep a residual +2 %/km with the proposed table, Temelin
+7 %/km on 34 days, Payerne is negative with every variant (−7 → −10 %/km, a site effect, not eta). The distorted
Payerne CL31 is shown for reference only.*

![slope histograms](figs_04v2/fig_cbh_slope_per_unit_hist.png)
*Per-unit slope distributions. The proposed table recentres the CL51 and CL61 distributions on zero and shifts the
CL31 one slightly below it; the width of each distribution (unit-specific shapes) exceeds the shift.*

Per-unit galleries: `figs_04v2/fig_cbh_binscatter_CL51_units.png` (48 units),
`figs_04v2/fig_cbh_binscatter_CL31_units_p1..p6.png` (211 units, 42 per page).

**Answer to "does it remove the altitude dependence?"** For the type averages, yes for the CL51 and the CL61
(slopes within the CI of zero), and the CL31 goes from +0.7 to −1.7 %/km, i.e. slightly over-corrected. Unit by unit,
no table makes every unit flat: the unit-specific residuals (±3–4 %/km, Lindenberg +2, Temelin +7, Payerne −10)
are larger than the difference between the tables, and they are not an eta property (section 7 of report v2:
FOV convention, near-range instrument terms, site effects). The proposed table reduces the number of units with a
slope above 2 %/km from 174 to 145 (of 270).

## 2. Recommendation

**Adopt the proposed inputs (a_G 3.3 um, type alpha_k, Mie S) as a new, versioned eta table, with a full rerun.**
The reasons are the inputs, not the checks:

- a_G = 5.5 um in the operational table was a Cloudnet effective radius (`der`, −30..+300 m above base, all
  statuses) used as if it were the equivalent-area radius; the chain 5.5 → 4.4 (layer) → 4.0 (`der_scaled`) → 3.3 um
  (r_e → a_G) corrects a known mis-assignment (report v2, section 3.1).
- alpha = 10 /km was a generic value; the ceilometers' own in-cloud tail decay gives 28–31 /km on 12.6 M profiles,
  consistent with Cloudnet at the 94 GHz sites (section 3.2). It is the input that moves C_L (−6 / −9 / −11 % at
  CBH 0.5 / 1 / 2 km); smaller drops alone would raise C_L.
- S: the Mie value at the measured droplet size (18.35 sr) replaces O'Connor's 18.8 sr (−2.4 %); the same code
  reproduces 18.9 sr at 10 um.
- The checks do not contradict the change: CBH flatness improves for CL51 and CL61 and is marginally worse for the
  CL31 (−1.7 %/km); the level checks at Payerne are symmetric (cloud vs Rayleigh+dark 1.03–1.09 operational vs
  0.93–0.98 proposed; vs CHM15k −5..−9.5 % vs +0.6..+5.6 %) and per unit the proposal is a −9 % shift with no net
  improvement (`eta_payerne_gap_2026-11-02.md`, section 4). None of them has the power to decide ±9 %.
- What remains uncertain either way: the Vaisala FOV / divergence convention (±5 % in level, ±1–2 %/km), the
  CL31 residual −1.7 %/km, and the autumn Rayleigh decline at Payerne (not an eta matter). A written answer from
  Vaisala on the FOV convention is worth asking for before or alongside the switch; it would not reverse the a_G and
  alpha corrections.

Against switching now: the change moves every Vaisala L2 attenuated backscatter up by 7–12 % (CBH-dependent) for
some 270 units and users; the CL31 overshoot; the FOV ±5 %. If the operator prefers to wait for Vaisala, the
operational table stays as is and this note documents the pending change.

## 3. What a switch implies, concretely

1. **Code:** replace `_ETA_CL31/_CL51/_CL61` in `calibration/cloud/_filters.py` by the tables of
   `step4/eta_tables.csv` (variant c; `eta_tables_proposed.py` writes them) and `S` by 18.35 sr in the cloud
   configuration; bump the cloud `calibration_version` written to the NetCDFs (currently 1000) so L2 users can tell
   the two generations apart; update `options.json` comments, report 02 and CLAUDE.md.
2. **Recalculate every Vaisala cloud constant from L1.** The eta enters the per-profile integral and production
   keeps no per-profile cache (`calibration/incremental.py` caches OmB and sensitivity only), so the exact
   per-profile rescaling of the study (`cbh_slope_windowed.py`, median ratio 0.998) is a diagnostic, not a product.
   Scope: all CL31 / CL51 / CL61 streams (about 270) over the whole archive (2025-01 to date), Rayleigh untouched
   (CHM15k, CL61 Rayleigh). Compute: a network month is 3–5 min as a SLURM array on the CSCS CPU partitions
   (`pp-short` / `postproc`, never GPU), so the full history is a few hours plus the merge; on the PC
   `scripts/run_streams_parallel.py` at 1–4 s per instrument-day is about half a day. Then the Kalman series,
   `extract_l2_opcoeff.py` → `operational_coefficients.csv`, and the per-night NetCDFs republished to the EWC bucket.
3. **Dashboards.** The operational dashboard: full `build_dashboard.py` (not `--changed-only`) + `ops/publish.sh`.
   The inter-comparison dashboard v3: rerun the Payerne / Lindenberg cloud streams of `diag_v22_l55s008` (and the
   CL31 ladder `diag_v22_cl31l910` / `_cl31wieg`), rebuild the payloads (`build_v3.py`, `render_v3.py`, `check_v3.py`);
   the Rayleigh variants and the dark do not change. The Hopkin heatmap dumps (`cloud_profile_dump_*`) must be
   regenerated for the new table.
4. **Paper.** Section 6 uses the dashboard constants (v2.2, 910.55/0.188, operational eta): the CL61 cloud rows at
   Payerne / Lindenberg / Camborne / Aosta and the CL51 / CL31 rows move by −6..−11 % in C_L (+7..+12 % in beta_att),
   so `run_paper_v09.py` must be rerun and the methods text must describe the new table (a_G, alpha_k, S) and the
   version. The Rayleigh rows and the CHM15k references are unchanged. The paper validation at 910.55 nm
   (`figs_paper_validation/paper_python_l55`) is the starting point.
5. **Order of work:** FOV question to Vaisala (parallel) → code + version bump → network rerun on CSCS → dashboards →
   paper → ops switch (`main` snapshot per the branch doctrine). Nothing of this has been started.

## 7. Operator decisions 2026-11-02 (after the note): one alpha, S unchanged

- **One alpha = 30 /km for the three types.** alpha is a property of the cloud, not of the instrument; the type medians (27.8 / 30.2 / 31.0 /km) differ by less than the between-unit scatter (SD 2.2-2.6 /km) and follow a regional pattern. Against the type-specific tables the common value moves C_L by -0.8..-1.3 pct (CL31), +0.1 pct (CL51), +0.3..+0.6 pct (CL61), inside the +-2-3 pct that the +-5 /km uncertainty of alpha already carries; the CL31 CBH slope moves by about -0.3 pct/km (from -1.7 to about -2.0 pct/km).
- **S stays at 18.8 sr** (O Connor et al. 2004, Hopkin et al. 2019, Cloudnet). The Mie value at the measured droplet size (18.35 sr at r_e 4.0 um) differs by 2.4 pct, inside the +-0.8 sr of the standard value, and the standard keeps the method comparable with the community. The C_L effect of the proposal is therefore the eta part only.
- **Effect of the proposal as decided** (a_G 3.3 um, alpha 30 /km, S 18.8 sr) on C_L against the operational table, at CBH 0.5 / 1 / 1.5 / 2 km: CL31 -5.0 / -7.4 / -8.5 / -9.1 pct; CL51 -4.1 / -6.4 / -7.8 / -8.6 pct; CL61 -4.0 / -6.3 / -7.6 / -8.5 pct. Sensitivity alpha 25 / 35 /km: -4.8 / -9.6 pct at 1 km (CL31), -4.2 / -8.4 (CL51), -4.1 / -8.2 (CL61).
- Tables (eta at the operational CBH nodes 0.25 .. 2.375 km; `validation.multiple_scattering_eta.eta_table(itype, 3.3, alpha=0.030)`):
  CL31 0.9046 0.8735 0.8252 0.7883 0.7588 0.7346 0.7142 0.6968 0.6818 0.6686;
  CL51 0.9312 0.9054 0.8637 0.8306 0.8031 0.7796 0.7594 0.7417 0.7260 0.7120;
  CL61 0.9318 0.9071 0.8664 0.8339 0.8069 0.7838 0.7639 0.7464 0.7310 0.7171.
- The per-unit cloud/Rayleigh ratios of section 4 for this configuration are the step-4 variant c (eta only): median over the ten dark-free units 0.978, Payerne hood-corrected 0.967; the CBH slopes of section 1 are unchanged (S does not enter them).

## 8. Alpha ladder 10 / 20 / 30 /km (2026-11-02, `validation/eta_study/alpha_ladder.py`, 43 study units, 1.0 M profiles, exact per-profile rescaling)

![alpha ladder](figs_04v2/fig_alpha_ladder.png)

C_L change against the 2026-07 tables at CBH 0.5 / 1 / 2 km (median per-profile ratio, median over units):

| table | CL31 | CL51 | CL61 |
|---|---|---|---|
| a_G 3.3 um, alpha 10 | +4.5 / +5.8 / +7.9 pct | +3.4 / +4.6 / +6.8 | +3.3 / +4.4 / +6.6 |
| a_G 3.3 um, alpha 20 | -1.4 / -1.9 / -2.6 | -1.1 / -1.6 / -2.4 | -1.0 / -1.5 / -2.3 |
| a_G 3.3 um, alpha 30 | -5.8 / -7.3 / -9.1 | -4.8 / -6.3 / -8.6 | -4.6 / -6.0 / -8.4 |
| a_G 3.7 um, alpha 20 | -2.6 / -3.4 / -4.5 | -2.1 / -2.9 / -4.1 | -2.0 / -2.7 / -4.0 |
| a_G 3.7 um, alpha 30 | -7.3 / -9.0 / -10.9 | -6.0 / -7.8 / -10.4 | -5.7 / -7.4 / -10.1 |

CBH slope of ln C_L (unit x month fixed effects, these 43 units): operational -0.2 / +3.5 / +2.3 pct/km (CL31 / CL51 / CL61); a_G 3.3: alpha 10 +2.2 / +5.9 / +4.6, alpha 20 -1.1 / +2.6 / +1.3, alpha 30 -2.9 / +0.5 / -0.8; a_G 3.7: alpha 20 -1.6 / +2.0 / +0.7, alpha 30 -3.2 / -0.1 / -1.3. CL61 cloud / Rayleigh (ratio of medians): Payerne vs Rayleigh + hood dark v3: operational 1.032, a_G 3.3 alpha 10 / 20 / 30 = 1.096 / 1.012 / 0.952, a_G 3.7: 1.085 / 0.996 / 0.933; median of the 10 dark-free units: 1.042, a_G 3.3: 1.091 / 1.026 / 0.981, a_G 3.7: 1.078 / 1.014 / 0.968. Reading: the whole alpha range 10-30 spans about +-6 pct in C_L at 1 km; the CBH slope is flattest at 30 /km for CL51 and CL61 and at 10-20 for the CL31; the cloud/Rayleigh ratio is centred closest to 1 at 20 /km and sits 2-5 pct below 1 at 30 (within the +-5 pct the dark-free references and the Payerne window dependence carry). None of these checks selects alpha; the inputs are the measured tail decay (28-31 /km, a lower bound) and the literature assumption (15-20 /km, O Connor 2004).

Measuring alpha during the calibration: the per-unit median tail decay alpha_k (CL61 26.6-37.8 /km, p10-p90 per unit about 17-44) was already used as step-4 variant d (per-unit alpha): it differs from the type value by at most 1.5 pct in C_L, because the between-unit spread (SD 2.2-2.6 /km) is small against the +-5 /km that move C_L by +-2.5 pct. A per-scene alpha would add the per-profile scatter (p10-p90 17-44 /km, +-7 pct in eta) without a measurable gain, and the tail decay is smoothed by the receiver response (effective 200-250 ns), so it under-estimates the true extinction: the measured value is a lower bound, not a per-night truth.

## 9. Extinction of the calibration clouds measured two ways, and the final choice alpha = 20 /km (2026-11-02)

`validation/eta_study/extinction_histograms.py`; outputs `C:/DATA/Projects/202611_eta_study/extinction/`.

![extinction histograms](figs_04v2/fig_extinction_histograms.png)

1. Forward inversion of the calibrated profiles (S = 18.8 sr, per-gate eta of the calibration, started below the cloud onset at beta_att = 20 Mm-1 sr-1): mean extinction over the first 50 m above the onset 6-7 /km (P10-P90 4-10), over the first 100 m 9 /km (6.5-11), at the signal peak 12-14 /km (9-17), optical depth 0.18-0.19 at the peak; optical depth 1 is reached 110-125 m above the onset (mean 8-9 /km to there); the backscatter-weighted extinction over the first two optical depths is 10-11.5 /km. Same numbers at Payerne CL61 / CL31, Palaiseau CL31 and Lindenberg CL61 (n 7 430 / 43 224 / 31 852 / 24 532 profiles).
2. Cloudnet 2 pi N a_G^2 in the lowest 50 / 100 / 150 m above the Cloudnet liquid base (the Vaisala cloud base, 30-60 m ABOVE the signal peak): Payerne CL61 26 / 33 / 36 /km, Payerne CL31 22 / 27 / 29, Palaiseau 31 / 35 / 40, Lindenberg 11 / 13 / 17 (medians; P10-P90 spans a factor 10).
3. Tail decay above the peak (step 2, +20..+150 m): 27-32 /km, a lower bound of the local extinction there.

Reading: the extinction of these clouds rises steeply with depth, about 5 /km at the onset, 13 /km at the signal peak 50 m higher, 25-35 /km another 100 m up. The three measurements describe one profile, not three answers. The PVC table assumes a homogeneous cloud; the layer that forms the integrated backscatter (optical depth 0-2) spans 10-30 /km, so one homogeneous value is a compromise. **Best estimate adopted for the release: alpha = 20 /km** (a_G 3.3 um, S 18.8 sr) - the middle of the measured range, the backscatter-weighted value of the first optical depth (10) and the deeper-layer values (30) bracketing it, and the value O Connor et al. (2004) assumed (15-20). It is not chosen on the checks: against the 2026-07 tables C_L moves by -1.0..-1.4 / -1.5..-1.9 / -2.3..-2.6 pct at CBH 0.5 / 1 / 2 km (section 8); the residual CBH slopes become -1.1 / +2.6 / +1.3 pct/km (43 units) and the Payerne cloud / Rayleigh+dark ratio 1.01. The remaining uncertainty of alpha (10-30 /km) is worth about +-5 pct in C_L, the FOV convention another +-5 pct.
