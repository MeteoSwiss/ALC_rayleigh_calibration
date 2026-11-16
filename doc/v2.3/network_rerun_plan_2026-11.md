# Plan: one combined network recompute — Rayleigh v2.3 (uncertainty budget, repaired tier 2, polar-summer fallback) and the cloud eta switch (2026-11-02, plan only, nothing started)

Scope: the whole 2026 L1 census (434 streams: 154 CHM15k, 210 CL31, 52 CL51, 13 CL61, 5 Mini-MPL) over
2025-01-01 to date, recomputed once with every pending change, validated on two months first, then flipped on
zueub434 with the v2.2 cutover procedure (`doc/OPERATIONS_v22_cutover.md`) as the template.

## 1. What goes into the release

| change | switch / place | state today |
|---|---|---|
| Rayleigh tier 2 repaired: measured noise (CHM15k, CL61) + vertical-consistency test (CHM15k) | `molecular_method: eprof_v2.3` (`calibration_version` 230) | code on `v2.3-dev` = `v2.3-integration`, off by default |
| Uncertainty budget (ensemble + vertical + molecular + aerosol + noise), reviewed defaults, relative QC, components written to the NetCDF | `uncertainty_model: v2.3`, `kalman_uncertainty_weighting: 1` | code merged (46993af), off by default |
| Polar-summer low-sun fallback (6 h around solar midnight, SZA >= 78 deg, flag 0.25 / -14) | `polar_night_fallback: 1` in `options.json` | code merged (0101912), off by default |
| CL61 water-vapour line as the code default (910.55 / 0.188) instead of the `ALC_WV_SPECTRUM` export | `water_vapor.LASER_SPECTRUM` | changed today, UNCOMMITTED (`water_vapor.py`, `tests/test_water_vapor.py`) |
| Cloud eta tables: PVC a_G 3.3 um, one alpha = 20 /km for the three types (best estimate from the measured extinction profile of the calibration clouds, `eta_switch_decision_2026-11-02.md` section 9), S unchanged at 18.8 sr; C_L -1.0..-1.4 / -1.5..-1.9 / -2.3..-2.6 % at CBH 0.5 / 1 / 2 km | `calibration/cloud/_filters.py` `_ETA_CL31/_CL51/_CL61` (13 nodes), cloud `calibration_version` 1100 | implemented 2026-11-02 (tests `tests/test_release_v23_switches.py`) |
| Dashboard: optical-unit / firmware change marks, HK columns | `monitoring/unit_changes.py`, `instrument_events.json` | merged; the full recompute regenerates the HK, so no separate `--hk-only` backfill |
| Operational switches kept as they are | no `ALC_DARK_PROFILE` (network homogeneity, operator decision 2026-08-17); CAMS 0.4 deg monthlies + daily, 1 deg per-month fallback | unchanged |

Versioning: output tree `ALC_calibration_v2.3`; NetCDF `calibration_version` 230 (Rayleigh) and 1100 (cloud) so
that the two eta generations are distinguishable downstream; the Kalman best estimate carries its standard
deviation for the L2 proposal.

## 2. Open items on the Rayleigh v2.3 list — what blocks the network recompute and what does not

Sources: `doc/v2.3/uncertainty.md` sections 4, 5, 8.4; `doc/v2.3/vertical_consistency.md` sections 6, 7, 10;
`doc/v2.3/polar_summer.md` sections 5, 6; the ref23 agreement study; the to-do of 2026-10-21.

**Must be closed before the recompute (they change what the NetCDFs contain):**

1. **Commit the working tree.** The CL61 spectrum default and its test, `options.json`, the eta switch, the
   version bumps. Run the test suite: the v2.2 bit-identity tests (`tests/test_uncertainty_v23.py`,
   `tests/test_default_method.py`, the sample-data tests) must stay green with the defaults, and a v2.3 fixture
   run must write the uncertainty components, the 0.25 flag and the new version numbers.
2. **Eta decision and implementation** (section 1): the cloud rows of every Vaisala stream change by -6..-11 %, so
   it has to be in the same recompute, not a second one.
3. **Polar-summer nights identifiable in the NetCDF.** The writer (`calibration/io/output.py`) already has two
   flag-type variables, `calibration_method` (0-3) and `calibration_version` (method codes, `VERSION_CODES`,
   flag_values / flag_meanings), but no per-night QC flag: a low-sun constant is written as a normal row and only
   the CSV and the dashboard see 0.25 (`polar_summer.md` section 5, last item). Lightest fix: a dedicated
   `calibration_version` code for the low-sun nights (231 next to 230), which reuses the existing flag machinery
   and needs no new variable; do it before the recompute, because the NetCDFs are what the hub consumes.
4. **Payerne CHM15k A on tier 1 only.** Recommended in `vertical_consistency.md` (section 6) but no per-station
   override exists in the runner (`grep` finds none). Either add a small per-key `molecular_params` override
   (`station_molecular_params` in `options.json`, read by the runner) or accept v2.3's behaviour for this unit
   (still +13 points of tier-2 nights at -18 %). It is the reference unit of every dashboard and of the paper, so
   decide it now.
5. **Which of OUR CL61 constants the hub applies.** The hub's own Rayleigh v1.0 on the CL61 is known to be wrong
   (Payerne L2 applies 2.126 since 2026-07-29 and halves its beta_att; replacing it is one motivation of this
   work), so this is not a question about the hub: the release has to state which constant the operational
   coefficients carry for the CL61 — the cloud constant (dark-immune, agrees with the CHM15k and Vaisala's factory
   scale within a few %) or the Rayleigh one (dark-free in the network, see item 6) — and the hub takes that.
6. **CL61 Rayleigh at high windows without a dark.** In v2.3 the CL61 keeps tier 2 with the measured noise and no
   vertical-consistency test. Its dark-free constants from 5-7 km windows are biased low by 15-30 % (the hood
   leverage, `eta_payerne_gap_2026-11-02.md`; ref23: Payerne fallback-tier nights +28.8 % against the cloud).
   Operationally the network stays dark-free (item 7 below), so either cap the CL61 window top for tier 2
   (`molecular_params` for the type, e.g. 5 km) or document that the CL61 Rayleigh constant is a diagnostic and the
   cloud constant is the product. Decide before the recompute; the cap is a one-line parameter.

**Should be done, can follow in the validation months or after:**

7. Molecular term of the budget north of 62 deg N: per-night CAMS T/p instead of the climatology (under-covers
   there: 79 % of nights within 2 sigma). The polar fallback adds nights exactly there. Cheap to evaluate on the
   two validation months (the Arctic streams), costs one CAMS read per CHM15k night if adopted.
8. Southern-hemisphere mirroring of the climatology is untested: check on Lauder and Auckland in the validation
   months.
9. Unit-gradient flag (|g| > 4 %/km, two-pass) — follow-up of `vertical_consistency.md` section 7, not in v2.3.
10. L2: carry `calibration_constant_0_uncertainty` and propagate it instead of the fixed 25 % (`uncertainty.md`
    section 5). Hub-side change; the recompute only has to write the Kalman standard deviation it needs.
11. The CL61 R(z) slope on clean nights (3-4.5 %/km) and the autumn decline of both Payerne Rayleigh constants
    (Aug -> Oct -22 %, cloud constants flat; `eta_payerne_gap_2026-11-02.md` section 5): watch them in the
    validation months (an autumn month), not blocking.
12. Report 19 section 7.15 (q17) and the attribution deck: documentation only.
13. Promotion to `main` was declined on 2026-10-21 ("not yet"); local `main` is at 497ab95, `origin/main` carries the
    v2.2 operational line and 388 commits of `v2.3-dev` are not on `main`. Promotion = snapshot of the operational
    subset (never a wholesale merge), after the validation months pass.

## 3. The recompute, phase by phase

**Phase 0 — freeze and smoke test (PC, 1 day).** Items 1-6 above; `pytest`; a three-stream run of Payerne A/B/C
over 2026-06-01..2026-10-31 with the release configuration (`ALC_MOLECULAR_METHOD=eprof_v2.3`,
`ALC_UNCERTAINTY_MODEL=v2.3`, `ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`, `polar_night_fallback: 1`) into a fresh
`ALC_FULLCAL_DIR`; inspect the NetCDFs (components, versions, flag) and the CSV flags. Sync the snapshot to
balfrin (`$HOME/alc_v23_cscs`, read-only via `PYTHONPATH`; `git fetch && git reset --hard` as usual).

**Phase 1 — validation months on CSCS (CPU partitions only; 1 day).** Two months over the full census, as a
per-stream array (`ops/cscs/alc_v22_*.sbatch` as templates, `partition=postproc` / `pp-short`, never a GPU
partition): a winter month (2026-01) and a polar-summer month (2026-06), each run twice — the current v2.2
configuration and the release configuration — into separate trees. Diff report (one script,
`validation/v23_reference/` harness as the base): yield per type and tier; distribution of the constant ratios
v2.3/v2.2 (tier-1 nights must be bit-identical, tier-2 CHM15k offsets near the q19 numbers, CL61 as decided in
item 6); uncertainty pulls against the twin and the CL61-cloud references (ref23); polar stations: nights added,
flag 0.25 share, level continuity; cloud: the -6 / -9 / -11 % shift at CBH 0.5 / 1 / 2 km and the per-type CBH
slopes (`cbh_binscatter_units.py` on the stored windows). Gate: no change outside the documented ones; otherwise
fix and repeat the month.

**Phase 2 — release recompute on CSCS (1 day wall).** 2025-01-01 to the last complete day, 434 streams,
`--methods rayleigh,cloud --force --ignore-coverage`, fresh tree `/scratch/mch/mhrvo/E_PROFILE_calout_v23`
(never two jobs in one tree). One 256-core postproc node took the v2.2 release within its 12-h limit (peak
147 GiB); a per-stream array is the alternative if the node is busy. Then OmB + sensitivity (`--no-cal --omb --sens`,
`--mem=24G`, per month) chained with `--dependency=afterany`, then the Kalman with weighting, the archive check
(`scripts/check_v22_archive.py` adapted to versions 230 / 1100, flag 0.25, the uncertainty components), a dashboard
build on balfrin viewed through the SSH tunnel for a last inspection, and the transfer to zueub434 as
`/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.3` (rsync through WSL, the v2.2 transfer section of the
cutover doc).

**Phase 3 — cutover on zueub434 (half a day).** Exactly the v2.2 procedure: one edit of `ops/config.sh` (the
output tree, `ALC_MOLECULAR_METHOD=eprof_v2.3`, `ALC_UNCERTAINTY_MODEL=v2.3`,
`ALC_KALMAN_UNCERTAINTY_WEIGHTING=1`, `ALC_OPCOEFF_CSV`; `ALC_WV_SPECTRUM` becomes redundant once the default is
committed, keep it one more cycle), the polar switch in the deployed `options.json`, gap fill from the archive's
last day to yesterday, payloads, a FULL dashboard build (no `--changed-only`), publish, the J+1 and J+7 checks,
rollback = re-comment the block and rebuild. `main` gets the operational snapshot at this point.

**Phase 4 — downstream (1-2 days).** `extract_l2_opcoeff.py` -> operational coefficients -> hub (with the CL61
decision of item 5 and the uncertainty proposal of item 10); the inter-comparison dashboard (rerun its
Payerne / Lindenberg cloud streams and the CL31 ladder with the new tables, rebuild the payloads; the dark-v3
rebuild running now is independent); the paper (section 6 from the new constants, methods text for the eta
table and the v2.3 Rayleigh, `run_paper_v09.py`); the Confluence pages of the three topics (ids in
`v23-release` memory) and a note to the E-PROFILE users: Vaisala beta_att up by 7-12 %, Rayleigh constants
unchanged on tier-1 nights, new uncertainty variables, new version numbers.

## 4. Compute and data

Timings measured on balfrin (single core): Rayleigh 2.2-4.1 s per instrument-day, cloud ~1 s. Full history:
172 Rayleigh streams x ~660 days x 3 s ~ 95 CPU-h; 275 cloud streams x 660 x 1 s ~ 50 CPU-h; both fit one
256-core node in a few hours, OmB + sensitivity are the larger part (native L1 month loads, 1-2.5 min per
stream-month, ~300 CPU-h). Data on scratch: L1 to date (Payerne complete to 2026-10-31; check the other streams'
front edge with `refresh_census.py`), CAMS 0.4 deg monthlies (`CAMS_Monthly_04`, verify the 2026-07..10 months
exist — the daily `CAMS/` stops at 2026-10-08), 1 deg monthlies as fallback. Scratch is purged: copy the release
tree down as soon as it is checked.

## 5. Decisions needed from the operator before phase 0

- eta switch: go / wait for Vaisala's FOV answer (inputs decided 2026-11-02: a_G 3.3 um, alpha 30 /km for all
  three types, S 18.8 sr).
- CL61 Rayleigh in the release: window cap for tier 2, or diagnostic-only with the cloud constant as the product
  (item 6), and the hub-side CL61 constant policy (item 5).
- Payerne A: per-station override or v2.3 as is (item 4).
- Version numbers and tree name (230 / 1100, `ALC_calibration_v2.3`).
