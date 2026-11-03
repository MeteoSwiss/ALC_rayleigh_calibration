# η study: handoff to a local Claude Code session (2026-11-02)

Read this file and `CLAUDE.md` before doing anything. You run on Max's desktop PC (32 cores, 128 GB, data local), in the repo `C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration`. Everything below was done today from a cloud project thread. This file replaces that thread's memory.

## 1. What the study is

The study reviews the multiple-scattering correction η (PVC, Hogan 2006) used in the O'Connor liquid-cloud calibration of the Vaisala CL31, CL51 and CL61. It also reviews the lidar ratio S. Conventions in the code:
- `beta *= eta`, B = ∫β′dz = 1/(2ηS), so **C_L ∝ η·S**.
- Report 04 and CLAUDE.md speak of the O'Connor multiplier C ∝ 1/η. "Higher η → lower C" means the multiplier, not C_L.

## 2. Where things are (all UNCOMMITTED; Max commits and pushes himself, never `git push`)

- Scripts: `validation/eta_study/`, for example `cloudnet_microphysics_match.py`. The operational code is untouched.
- Report v2: `doc/reports/04_multiple_scattering_v2.md`, figures in `doc/reports/figs_04v2/`.
- Report v1 `doc/reports/04_multiple_scattering.md`: two correction notes were added today on the 5.5 µm (lines ~22 and ~636). The 910.55 nm audit flags were also added today.
- WV audit log: `doc/reports/wv_cl61_91074_audit.md`.
- Data and outputs: `C:\DATA\Projects\202611_eta_study\` (`REPORT.md` with inline figures; `fr_figs\` holds 8 PNGs and `captions.txt`).
- French explainer page (Max's claude.ai, private): https://claude.ai/artifact/QEsNDNhJm4E5JnRwVWESux. It has 9 figures, the CHM15k check and the 5.5 → 3.3 µm chain.

## 3. Results (current recommendation, revised 15:10)

**Inputs, each fixed from physics, measurements or specs (never from the checks):**

| Input | Value | Source |
|---|---|---|
| a_G (equivalent-area RADIUS) | 3.3 µm (r_e ≈ 4.0 µm) | Cloudnet `der_scaled`, `der_retrieval_status==1`, lowest 100 m above liquid base, 3 sites. a_G = 0.832 r_e (lognormal σ 0.35; CloudnetPy der default N = 200 cm⁻³, σ 0.35). 50/100/150 m give 3.0/3.3/3.6 µm. |
| α | CL31 27.8, CL51 30.2, CL61 31.0 km⁻¹ | In-cloud tail decay (+20..+150 m above peak), 270 units, 12.6 M profiles. It is a lower bound with respect to the receiver response. Its pattern is regional (subarctic −5.6 km⁻¹), not firmware. It matches Cloudnet at 94 GHz sites (0.85–0.92); Lindenberg 35 GHz is off (2.5). |
| FOV / divergence | Wiegner 2014 Table 1 half-angles (as in code) | The convention (edge, 1/e, 1/e²) is UNKNOWN (Max doesn't know). This gives ±1–2 %/km on the slope and about ±5 % on the level. |
| S | 18.35 sr (Mie, 910.55 nm, r_e 4.0) | The Mie code gives 18.9 sr at 10 µm, matching O'Connor's 18.8. |

The 5.5 µm of July was a Cloudnet r_e (`der`, CBH −30..+300 m, all statuses), used as if it were a_G. The chain is 5.5 → 4.4 (layer) → 4.0 (der_scaled) → 3.3 µm (r_e→a_G). Smaller drops RAISE C_L. The C_L drop comes from α (10 → ~30 km⁻¹).

**Effect:** C_L −6 / −9 / −11 % vs operational at CBH 0.5 / 1 / 2 km, the same for all three types (η part −4 / −6.5 / −9 %, S part −2.4 %). ±1σ is about ±3–5 %.

**Checks (diagnostics only):**
- De-seasoned CBH slope, proposed: CL31 −1.7 %/km (fails), CL51 −0.3 (passes), CL61 −0.4 to −0.9 (passes). Operational: +0.7 / +2.8 / +2.3. Hopkin: +7 to +8.
- CL61 cloud/Rayleigh C_L, operational → proposed: Payerne + dark 1.053 → 0.944, Lindenberg (no dark) 1.037 → 0.946, 9 other CL61 1.048 → 0.963. Not decisive.
- Payerne vs CHM15k L2 v1.0 (v3 inter-comparison dashboard, June–Aug 2026, 678 h): CL61 cloud −2.3 % → about +9 % (scaled). CL61 Rayleigh + dark is +7.8 %. CHM15k v2.2 is −13 % from its own v1.0. This is the most unfavourable check: the proposal moves the CL61 cloud from the CHM15k onto its own Rayleigh + dark.

## 4. Hard rules (from Max)

- **A flat C_L vs CBH and the Rayleigh/CHM15k level are DIAGNOSTICS.** Never choose an input (a_G, α, FOV, S, WV line), and never keep a type's table, because a check passes. Scans may appear only as sensitivity tables, never as a "best" or "null" value.
- CAMS on the PC: `D:/CAMS_Monthly_04;D:/CAMS_daily`. NEVER `D:/CAMS` or `D:/CAMS_run_v20` (1° orography, which truncates the WV column). The code defaults to `A:` and silently drops days.
- CL61 WV line: 910.55 nm / FWHM 0.188 nm. Results with 910.74 / 1.0 are obsolete.
- Rerunning July scripts: re-apply `tests/_ceilo_reader.average_ceilo_data` (averaging moved into the loader on 2026-07-10). The old η tables are in git `83f451f^`.
- Lindenberg CL61: the unit was swapped on 2025-04-28. Use data from 2025-05-01 only. It has no dark measurement.
- Payerne CL31 is a distorted unit: flag it and show results with and without it, never as the only evidence.
- Compute: prefer this PC. If several jobs share it, the total is 32 cores. CSCS: CPU partitions only (Linard: postproc, pp-serial, pp-short, pp-long, xfer; never GPU, never `--gres/--gpus`). Keep ≤ 6 of 14 Linard nodes; `~/squ` double-counts nodes. Don't edit `~/bin/squ_gentle.sh` or clean the over-quota Linard `$HOME`.
- Ask before moving data. Use `/tmp` for scratch on zueub434.

## 5. Open items, in recommended order

1. **Explain the 8 % gap between the CL61 Rayleigh + dark and the CHM15k at Payerne.** It does not depend on η. Check which CHM15k constant (v1.0 L2, v2.2 or v2.2 + dark, which spread 13–21 %) and the 910 → 1064 nm conversion (report 07 §8). Do this before any operational η switch.
2. FOV convention of the CL31, CL51 and CL61 (Vaisala).
3. The CL31 −1.7 %/km residual, explained from independent information: FOV convention, CL31 laser spectrum 909.7 / 6 nm and its temperature drift, a range-dependent instrument term. No tuning.
4. Lindenberg CL61 dark (none for the new unit) and Payerne dark (17 %, far above Le et al. 2026).
5. Max decides whether and when to change `_ETA_CL31/_CL51/_CL61` and S in operations. Nothing is applied.

## 6. Session 2026-11-02 (local PC): open item 1 closed, per-unit eta question answered

See doc/reports/eta_payerne_gap_2026-11-02.md. The +7.8 pct was the 6-session v1 hood dark; with dark v3 on both instruments the CL61 Rayleigh+dark and the CHM15k agree +1 +-3 pct in every month Jun-Oct 2026. The proposed eta shifts every CL61 unit by -9 pct and does not improve cloud/Rayleigh agreement (median |ln ratio| 0.052 -> 0.075). New open item: both Rayleigh constants fall ~22 pct Aug->Oct while cloud constants stay flat (CL61 intervention 2026-07-28, cooling 9-12 K); do not use Sep-Oct ratios for eta. Payerne L1 synced to 2026-10-31; CAMS daily ends 2026-10-08. Dashboard still on dark v1 (rebuild = separate task). Remaining open items: 2, 3, 4, 5 of section 5.

- 2026-11-02, later: per-unit CBH binscatters and the switch recommendation in `doc/reports/eta_switch_decision_2026-11-02.md` (`validation/eta_study/cbh_binscatter_units.py`). Recommendation: adopt a_G 3.3 um / alpha_k / Mie S as a versioned eta table with a full cloud rerun of all Vaisala streams, dashboard rebuilds and paper section 6 redo; Max decides.

- 2026-11-02, later: combined network recompute PLAN (Rayleigh v2.3 + uncertainty + polar fallback + eta switch) in `doc/v2.3/network_rerun_plan_2026-11.md`, with the open v2.3 items split into blocking (commit, eta decision, polar flag in NetCDF, Payerne A tier-1 override, CL61 hub constant / 2.126, CL61 tier-2 window cap) and follow-ups. Nothing started.

- 2026-11-02 evening: phase 0 done (switch tests, Payerne + polar smoke); final eta inputs a_G 3.3 um, **alpha 20 /km** (extinction histograms: 5 -> 13 -> 30 /km with depth), S 18.8; release launched on Linard as a per-stream array (see the session report / RELEASE_PROVENANCE.txt in /scratch/mch/mhrvo/E_PROFILE_calout_v23_rel).

- 2026-11-02 23:00 LAUNCHED on Linard (balfrin CPU partitions down): release array **337836** (0-433 pct 96, `ops/cscs/alc_v23_release_array.sbatch`, code 46bb407 staged at `/scratch/mch/mhrvo/alc_v23_code`, keylist `/scratch/mch/mhrvo/keys_v23_release.txt`, output `/scratch/mch/mhrvo/E_PROFILE_calout_v23_rel`, logs `/scratch/mch/mhrvo/logs/alc_v23_arr_337836_<task>.log`), add-ons array **337837** chained afterany (OmB + sensitivity). Tomorrow: `squeue --me`, count `<key>/<key>_cal.csv` and RELEASE_PROVENANCE.txt, grep the logs for Traceback, then adapt `scripts/check_v22_archive.py` to versions 230/231/1100 and run it; the network L1 on CSCS ends 2026-08-13 (Payerne to 11-01), the tail is filled on zueub434 at the cutover.

- 2026-11-03 09:00 INCIDENT: the first array (337836) ran without `calibration/data/abs_cross_wv_910nm.nc` (my rsync excluded *.nc), so every 910 nm calibration (CL61 Rayleigh+cloud, CL31/CL51 cloud, 275 streams) failed with flag -99 while the CHM15k / Mini-MPL Rayleigh results (159 streams) are valid. File restored in both code copies; the Vaisala streams were resubmitted as array **339769** (keylist `/scratch/mch/mhrvo/keys_v23_vaisala.txt`, `KEYS_FILE` env of the sbatch), the add-ons re-chained as **339770** (afterany 337836 and 339769); the old add-ons 337837 cancelled. The three Dutch high-resolution CHM15k tasks of 337836 were still running (7 h). Comparison framework: `validation/v23_release/compare_v22_v23.py` on the local copies `C:/DATA/Projects/202611_v23_release/calout_v2{2,3}_rel`.

- 2026-11-03 midday: Vaisala rerun 339769 complete except the two CL51 streams Uccle/Diepenbeek (3.5 h+); comparison script extended (retyped streams, cloud-day transitions, budget components, Kalman std) and the release report written: `doc/reports/v23_release_report_2026-11.md` (figures `figs_v23_release/`). Findings: tier-1 nights identical, CHM15k -8 % nights (tier-2 rescues removed), CL61 -16 % (cap), cloud -1.6..-1.9 %, 160 low-sun nights (Oslo 2.6x and South Uist 3.0x the dark level = outliers to inspect), budget 14 % vs 5 % ensemble. DEFECT FOUND: three Czech streams swapped CL31 -> CL51 under the same ident (Plzen 2025-07-16, Tusimice 2026-06-25, Libus 2026-06-26; `validation/v23_release/scan_type_changes.py`): the census carries one type per stream (newest file), so v2.3 corrected their CL31 part with the CL51 eta table (+3-4 %); events recorded in `monitoring/instrument_events.json`; fix = per-day type from the L1 range extent + Kalman restart at the swap. Also: the scratch L1 covers the whole network to 2026-10-31 (the v2.3 tree runs to 2026-11-01; the v2.2 one ended 2026-08-13). Flag 0.5 means a PARTIALLY CLOUDY night (filter_cloudy_profiles), not the fallback tier: tier-2 nights carry flag 1 or 0.5 like tier-1 nights; the CSV/NetCDF do not record the tier.
- 2026-11-03 13:30: rerun 339769 finished (275/275, rc 0, no traceback), local v2.3 copy final, report numbers final (`compare/summary_final.txt`). Add-ons array 339770 (OmB + sensitivity, 434 tasks) running on Linard. Next: archive check adapted to 230/231/1100, copy of the complete tree, transfer to zueub434, cutover (report section 11).
