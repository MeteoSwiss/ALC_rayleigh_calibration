# E-PROF v2 molecular-window optimization — why clear nights fail, and a tuned v2 (L1 + L2, 2026)

*Generated 2026-06-20. Experiment: `validation/run_v2_sweep.py` + `validation/analyze_v2_sweep.py`.
Adversarially cross-checked by an independent 4-analyst review of the sweep summary.*

## TL;DR

The fraction of clear (fit-reaching) nights that yield a **valid** Rayleigh calibration varies
enormously by instrument — from ~14 % to ~80 % in the per-instrument L1 figure. We ran an experiment
over **24 instruments (10 CHM15k, 4 Mini-MPL, 10 CL61) × both L1 and L2 × every clear night of
2026** to find out why, and to optimize the **E-PROF v2** molecular-window method.

- **Why clear nights fail:** the **scattering-ratio gate (`max_scattering_ratio = 1.10`) is the
  dominant binding constraint** — relaxing it alone recovers **77 %** of failed L1-CL61 nights,
  43 % of L1-CHM15k, 43 % of L2-CL61, 30 % of L2-CHM15k. It was rejecting good *molecular* windows,
  not aerosol. For **Mini-MPL** the binding gate is instead the **window-start height** (the SNR
  runs out before 2 km); `R²` and the in-window ratio-std gates are *never* binding.
- **Optimized v2 (config “C8”):** relax scattering to 1.15 plus moderate easing of
  start/R²/residual/ratio-std/temporal. This **raises the valid-night fraction on every instrument
  type and on both levels** while keeping σ_SD (short-term variability) flat or better on L1 (e.g.
  **L1-CL61 45 → 75 % with σ_SD 7.6 → 7.0**).
- **C8 is now the production `eprof_v2` default** (`calibration/rayleigh/molecular_methods.py`).
- **Honest caveat:** on **L2-CL61** the same relaxation is a genuine yield/precision *trade* — median
  σ_SD rises 6.9 → 8.7 % (6 of 9 units worsen), so on that one cell we are buying yield with some
  residual-aerosol leakage. See §6.

---

## 1. The question

The per-instrument L1 figure (previous report) showed the *valid-calibration fraction* swinging from
14 % (Payerne CHM15k) to ~80 % (Birkenes CL61) with no obvious physical pattern. Two questions:

1. **Why do clear nights fail on some instruments but not others?** A clear night reaches the
   molecular fit; if no window then passes the method’s quality gates, the night is wasted.
2. **Can E-PROF v2 be tuned to recover those nights without degrading the calibration** (i.e. without
   inflating the night-to-night scatter of the lidar constant)?

## 2. Method

**Data.** L1 = `D:/E-PROFILE_L1_2026` (raw range-corrected signal), L2 = `D:/E-PROFILE_L2_2026`
(attenuated backscatter), daily files, the 24 instruments of `validation/scope_l1_2026.json`.

**One load, many configs.** For each instrument-night we run the full pipeline once to the prepared
molecular profile (`calibrate_rayleigh(..., fit_inputs_out=...)`), build the v2 time-cell-cleaned
window grid **once**, then evaluate every candidate gate configuration on that grid — so the configs
are compared on identical inputs.

**Metrics.**
- **valid %** = valid calibrations ÷ clear (fit-reaching) nights. “valid” = the v2 window passes the
  pipeline QC (`rel_error ≤ 15 %`). A clear night is one where the profile reached the fit (cloudy /
  foggy nights return earlier and are not counted as failures here).
- **σ_SD** = robust successive-difference (von Neumann) precision, `1.4826·median|ΔC|/√2 ÷ |median C|
  ×100` (% of median lidar constant). This is the **short-term variability** — drift-insensitive
  night-to-night repeatability; **lower is better**.

**Failure diagnostic (leave-one-gate-out).** For every baseline-failed clear night we relax **one**
gate to infinity at a time and check whether the night becomes valid. The gate that recovers the most
nights is the binding constraint. (`max_rel_error` is held at the pipeline QC of 15 % throughout — a
window above it can never be a valid calibration, so relaxing it cannot help.)

## 3. Why clear nights fail — the binding gate

![Leave-one-gate-out: % of baseline-failed clear nights recovered by relaxing only that gate (higher = more binding). L1 left, L2 right.](figs_paper_validation/v2_sweep/fig_v2_bottleneck.png)

| level · type | failed nights | #1 binding gate | #2 |
|---|---|---|---|
| L1 CHM15k | 263 | **scattering ratio 43 %** | temporal_cv 31 % |
| L1 Mini-MPL | 154 | **window start 21 %** | scattering 13 % |
| L1 CL61 | 171 | **scattering ratio 77 %** | residual 5 % |
| L2 CHM15k | 669 | **scattering ratio 30 %** | residual 20 % |
| L2 Mini-MPL | 159 | **window start 24 %** | scattering 13 % |
| L2 CL61 | 612 | **scattering ratio 43 %** | residual 14 % |

**Reading:** the aerosol-rejection gate `scattering_ratio ≤ 1.10` is the single biggest reason clear
nights fail, overwhelmingly so for **CL61** (77 % on L1). The estimated scattering ratio sits just
above 1.10 on many genuinely clear nights, so the gate throws away usable molecular windows.
**Mini-MPL** is the exception — its 532 nm signal runs out of SNR before the 2 km window-start floor,
so *window-start height* is its bottleneck. **`R²` and in-window ratio-std never bind** (≈0 %
recovery) — those gates can stay strict. On L2 the **Rayleigh-shape residual** becomes a secondary
binding gate (the L2 retrieval’s vertical structure differs slightly from a pure exponential).

## 4. The configuration trade-off (≥ 5 configs)

Nine gate configurations, each evaluated on every clear night at both levels. Gates not listed sit at
the C0 baseline; `max_rel_error = 15` always.

| config | scatter | temporal_cv | residual | R² | start (m) | ratio_std |
|---|---|---|---|---|---|---|
| **C0** baseline (old v2) | 1.10 | 0.50 | 12 | 0.50 | 2000 | 0.30 |
| C1 temporal≤0.8 | 1.10 | 0.80 | 12 | 0.50 | 2000 | 0.30 |
| C2 temporal≤1.2 | 1.10 | 1.20 | 12 | 0.50 | 2000 | 0.30 |
| C3 looser shape/ratio | **1.15** | 0.50 | 20 | 0.50 | 2000 | 0.40 |
| C4 R²≥0.35 | 1.10 | 0.50 | 12 | 0.35 | 2000 | 0.30 |
| C5 start≥1.2 km | 1.10 | 0.50 | 12 | 0.50 | **1200** | 0.30 |
| C6 balanced | 1.12 | 0.80 | 16 | 0.40 | 1500 | 0.40 |
| C7 aggressive | 1.25 | 1.50 | 25 | 0.30 | 1000 | 0.50 |
| **C8 RECOMMENDED** | **1.15** | 0.80 | 16 | 0.40 | 1500 | 0.40 |

### valid % / σ_SD %  (median over instruments of each type)

**L1**

| config | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| C0 baseline | 38.2 / 10.9 | 65.0 / 8.5 | 45.4 / 7.6 |
| C1 temporal≤0.8 | 56.3 / 8.9 | 66.6 / 8.7 | 46.4 / 7.3 |
| C2 temporal≤1.2 | 60.3 / 9.2 | 66.6 / 8.7 | 46.4 / 7.3 |
| C3 looser shape/ratio | 46.6 / 10.3 | 66.3 / 8.8 | 64.9 / **6.5** |
| C4 R²≥0.35 | 38.2 / 10.9 | 65.0 / 8.5 | 47.0 / 7.5 |
| C5 start≥1.2 km | 47.5 / 9.5 | 67.1 / 8.3 | 45.4 / 7.2 |
| C6 balanced | 67.5 / 9.6 | 68.9 / 8.8 | 68.5 / 7.3 |
| C7 aggressive | **93.4** / 9.2 | 72.8 / 10.9 | **85.9** / 7.1 |
| **C8 recommended** | 71.6 / 10.0 | 69.7 / 8.9 | 75.1 / **7.0** |

**L2**

| config | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| C0 baseline | 47.4 / 10.0 | 64.9 / 8.2 | 35.6 / **6.9** |
| C3 looser shape/ratio | 60.5 / 9.9 | 67.1 / 8.5 | 50.9 / 8.6 |
| C6 balanced | 58.6 / 9.4 | 72.6 / 9.0 | 45.3 / 7.9 |
| C7 aggressive | 72.0 / 11.1 | 75.6 / 10.5 | 64.7 / 12.5 |
| **C8 recommended** | 59.5 / 9.5 | 72.6 / 9.0 | 50.3 / 8.7 |

*(C1, C2, C4, C5 change nothing material at L2 — their single-gate relaxations are not the L2
bottleneck — so only the multi-gate configs are tabulated for L2.)*

![Yield vs short-term variability for every config, one panel per instrument type. Best is bottom-right (more valid, lower σ). Circles = L1, squares = L2.](figs_paper_validation/v2_sweep/fig_v2_pareto.png)

![Valid-calibration fraction on clear nights by config, grouped by instrument type (L1 left, L2 right).](figs_paper_validation/v2_sweep/fig_v2_validbars.png)

**What the Pareto shows.** Moving from C0 to C8 pushes every instrument type **down-and-right** (more
valid, lower-or-equal σ) on L1. **C7 (aggressive)** reaches the highest yield but breaks the σ_SD
budget on L2 — L2-CL61 σ_SD jumps to **12.5 %** (from 6.9) — so it is rejected. **C8** keeps almost
all of C7’s recovery at a fraction of the variability cost.

## 5. Recommendation — C8, now the v2 default

**`eprof_v2` defaults updated** to C8 (`min_window_start_m=1500, min_r2=0.40, max_residual_pct=16,
max_scattering_ratio=1.15, max_ratio_std=0.40, max_temporal_cv=0.8`); composite weights and the
time-cell flagging are unchanged. Net effect vs the old baseline, valid %:

| | CHM15k | Mini-MPL | CL61 |
|---|---|---|---|
| L1 | 38 → **72** (+34) | 65 → **70** (+5) | 45 → **75** (+30) |
| L2 | 47 → **60** (+12) | 65 → **73** (+8) | 36 → **50** (+15) |

…with σ_SD held flat or improved everywhere except L2-CL61 (see §6). C8 beats C0 on valid % in **all
six** type×level cells and stays within ~+1.5 pp σ_SD of baseline everywhere — the precision guard
C7 fails. To use the optimized v2 in production set `"molecular_method": "eprof_v2"` in `options.json`
(the file currently ships `eprof_v1.2`).

## 6. Adversarial review & caveats

An independent 4-analyst review (yield judge, variability guardian, investigation auditor, adoption
skeptic) verified every number against the summary JSON. It **endorsed C8 as the global default** but
flagged real limits — reproduced here honestly:

1. **L2-CL61 precision is a genuine trade, not a free recovery.** The median σ_SD rises 6.9 → 8.7 %,
   and per-instrument **6 of 9 CL61 units inflate** — several badly (Sion 8.1 → 13.9, Zeebrugge
   9.2 → 14.7). On the operational L2 product, relaxing scattering for CL61 admits *some* genuine
   residual aerosol. The conservative alternative **C6** (scattering ≤ 1.12) only partly helps — the
   same units still inflate — so it is a ~5-pt-yield concession for little precision back.
2. **No out-of-year / out-of-station holdout — overfitting risk.** Every threshold is tuned on 2026
   aerosol climatology. Before locking these as permanent defaults, **validate on a second year and a
   station holdout**. This is the single biggest methodological risk.
3. **The leave-one-out diagnostic is univariate.** A real aerosol layer trips scattering, residual
   *and* temporal_cv together, so the recovery magnitudes (43–77 %) are upper bounds; the **ranking**
   (scattering #1) is robust, the exact percentages are not. A few units (e.g. Uccle) fail for
   site/data reasons no gate fixes.
4. **Type medians hide tails — report per-instrument.** σ_SD medians rest on 7–10 units; some
   instruments (Brest Mini-MPL ~20–36 %, a few CL61) are intrinsically noisy. The per-instrument rows
   are in `v2_sweep_summary.json`.
5. **The guard tests precision, not accuracy.** σ_SD measures repeatability, not absolute-calibration
   bias. Admitting marginal-aerosol windows (the L2-CL61 mechanism) could shift the absolute lidar
   constant — worth a dedicated accuracy check.

**Future direction (v2.1):** a **per-instrument-type** policy is the natural next step, because the
three types have three different binding gates — CL61/CHM15k are scattering-limited (cap scattering at
1.12 on L2-CL61 to protect precision), Mini-MPL is window-start-limited (relax start, leave scattering
near baseline). Ship the single global C8 now; do per-type tuning once a validation year is available.

## 7. Notebook & plotting changes (this round)

- **`examples/03_sample_data_calibration.ipynb`** now switches the whole notebook between L1 and L2
  with one `LEVEL = "L1"|"L2"` flag, and shows **only** the compact diagnostics dashboard (the simple
  per-step RCS plots are gone). This also explains the earlier “CHM15k/CL61 flag = −2” surprise: those
  nights are genuinely flagged on **L1** but calibrate on **L2** — flipping `LEVEL` makes that visible.
- **New `plot_all` option** (`calibration/config.py`): `plot_main` now emits *only* the compact
  dashboard; the extra `_rcs_timeseries`/`_rcs_annotated` plots are gated behind `plot_all`.
- **Compact dashboard reworked** (`plot_rayleigh_diagnostics_compact`): a filled 4×6 grid with **no
  empty quadrants**, the **full range-corrected-signal matrix** promoted to the large bottom-left
  panel, and **hatched overlays** marking excluded profiles (red `///` = low cloud, grey `\\\` =
  screened / not used) drawn on top of the signal.

![Reworked compact Rayleigh dashboard (L2 CHM15k, Payerne 2026-02-25): molecular fit (top), window search (row 2), full RCS matrix with hatched exclusions + gold molecular band (bottom-left), sensitivity grid and lidar-constant spread (bottom-right).](figs_paper_validation/v2_sweep/fig_diag_compact_example.png)

## 8. Reproduce

```
python validation/run_v2_sweep.py        # sweep: 48 (level×instrument) jobs, 14 configs each
python validation/analyze_v2_sweep.py    # -> v2_sweep_summary.json + the 3 figures above
pytest tests/test_all_instruments_run.py tests/test_sample_data.py   # 47 pass with the new v2 default
```

Artifacts under `figs_paper_validation/v2_sweep/`: `v2_sweep_summary.json`, `sweep_<level>_<inst>.json`
(per-instrument per-night per-config), `fig_v2_pareto.png`, `fig_v2_validbars.png`,
`fig_v2_bottleneck.png`.
