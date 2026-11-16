# Remote dark M1+M3 — implementation state after the first session (2026-08-20)

Branch `remote-dark`, code `remote_dark/`, outputs
`C:/DATA/Projects/202606_E-PROFILE_calibration/remote_dark/`. Everything below ran on real data;
every number is from a run, not a projection. Plan: `14_remote_dark_methodologies_plan.md`.

## What stands

**Model families (M1 step 1) — the truth ceilings are established.** Fit of each type's
parametric family to the Payerne hood truth (`remote_dark/models.py`):

| type | family | ceiling |
|---|---|---|
| CL31 | flat + two damped resonances (9 par) | **R²(P) = 0.998 at 60–1500 m** — the near-range ring is fully representable |
| CHM15k | flat + negative exponential (3 par) | **R²(b, structured) = 0.942**, Rayleigh-band misfit 6 % |
| CL61 | hood template × amplitude | exact by construction |

The CL31's 96 % Rayleigh-band "misfit" is the gate-locked ripple — it belongs to the gate-fold
channel (anchor 3), not to the smooth family, and is already retrievable network-wide.

**Hood-dynamics probe (`probe_hood_dynamics.py`) — three design facts, all new:**

* the background variable **sweeps under the hood** (0→18, same order as open sky): `bckgrd_rcs_0`
  is partly detector-driven, not a pure solar readout;
* the dark's modulated response is **near-range concentrated**: |corr(dark, B or T)| up to
  0.68–0.72 at 60–500 m, ~0.1 above 2 km — the modulated component lives exactly where the static
  CL31 dark lives;
* it is **unit-dependent**: the post-swap CL31 optic block responds ~2× more than the pre-swap
  one. D_bg/D_T are per-unit products; a type constant would be wrong.

**M1 anchored fit (`m1_fit.py`) — shape recovery works at first attempt.** Family projected onto
the v4 per-gate evidence, nightly amplitudes re-fit with a soft C_ref anchor, three alternations,
on the existing v4 night caches:

| stream | nights | vs hood (2–8 km) |
|---|---|---|
| Payerne A (CHM15k, post-swap era) | 216 | **corr +0.94** — the cuvette SHAPE is recovered |
| Payerne A (era 20250601-) | 190 | corr +0.89 |
| Payerne C (CL61) | 105 | corr +0.86 |
| Payerne A (pre-swap, 27 nights) | 27 | corr +0.34 — too few nights, consistent with v4 experience |

This is the thing the free-b v4 could not do on the CHM15k (the collinear cuvette). The family
constraint is doing its job on shape.

## What is broken or refuted, named plainly

* **M1 amplitude is wrong** — gain +2.4e7 on the CHM15k, ~0 on the CL61 template step. A
  units/scale bug in the anchor prior (`C_ref` units vs the cache's A_n) and/or the template WLS,
  to be instrumented and fixed next session. Until the amplitude closes, M1 delivers shape only.
* **M3 twilight with the RAW background regressor is refuted** — sky-vs-hood corr ≈ −0.2. The
  terminator's background is locked to the boundary-layer evolution itself (dawn convection rises
  with the sun), so the regression attributes real atmospheric change to B. Kept in the module
  docstring as a negative result.
* **High-passed regressor is insufficient as validated so far** (corr −0.09/−0.19), and the deeper
  problem is on the TRUTH side: under the hood there are no cloud shadows, so the hood sessions
  contain almost no fast-B variation — the fast-B coupling truth is weakly defined. M3 validation
  must be redesigned as a JOINT model (D_bg·B + D_T·T against the hood's slow drift, where the
  25 h sessions have full leverage) rather than component-wise.

## Next session, in order

1. Fix the M1 amplitude: print/trace `C_ref`, median A_n before/after anchoring, and the template
   normalisation; then re-run the three Payerne streams. Success = amplitude gain within
   ~±30 % with the shape correlations held.
2. M1 injection self-test (the v4 machinery generalised): inject the hood family at known
   amplitude into a control stream, require unbiased recovery.
3. M3 joint fit on the two 25 h/21 h hood sessions: dark(t, z) = D_bg(z)·B(t) + D_T(z)·T(t) with
   the collinearity handled by the two sessions' different B/T phase relationships; then the same
   joint model on sky twilights.
4. CL31 sky test for M1: build the missing Payerne B night cache (v4 `--start 20250101`), fit,
   compare near range — the band where the family ceiling is 0.998.

---

## Session 2 (same day, continued) — the amplitude closed where it can close

Chronology of the amplitude chase, each step forced by a measurement:

1. The "units bug" was three bugs. (a) The hood npz stores `_b` (P-view) and `_b_rcs` (rcs view);
   the loader had them semantically swapped, so every session-1 comparison was against the wrong
   view — the "+2.4e7 amplitude" was z² at 5 km. (b) The bounds of the family step walked with the
   iterate instead of anchoring on the hood prior (τ contracted ×0.6 per iteration to the floor).
   (c) The CHM15k family itself was wrong: the hood P-view PEAKS at ~3 km; a monotone exponential
   cannot, so its prior τ was garbage. Replaced by −A·(z/1km)^k·e^(−z/τ) (peak at kτ = 2.7 km,
   matching the hood; fit band from 1 km — below that the hood shows an overlap-region +P spike
   that is not subtractable dark).
2. The anchor went through three designs: era-scoped Kalman (rejected — the Kalman smooths across
   the unit swap, the documented CL31 failure mode); era-scoped raw-cal median (rejected — the v4
   cache amplitudes and the archive constants differ by a CONVENTION factor ~2.25, so an absolute
   import drags everything); final: **per-night relative anchor** A_n ~ N(k̂·C_n, σ) with the
   units factor k̂ estimated once — the archive contributes its independent night-to-night
   structure, the absolute scale stays with the data.
3. The anchored alternation was dropped for the template and gexp families (their far-range shapes
   brush the molecular/z² directions; every night step drained them — measured, not assumed).
   Single family step on the v4-initialised state.
4. **The injection self-test is now part of every run**: the same estimator is re-run on the same
   nights + 1× the hood dark; the recovered increment is the estimator's gain, and
   θ_corr = θ_raw / gain.

**Verdicts:**

| stream | shape corr (2–8 km) | injection gain | θ_corr | verdict |
|---|---|---|---|---|
| Payerne C (CL61, 105 nights) | **+1.00** | 0.76 | **+1.13** | **closed** — remote dark within ~13 % of the hood |
| Payerne A (CHM15k, 216 nights) | +0.24 (near-range +0.98) | 0.48 | +5.2 ≠ 1 | **not closable from one station** — the raw θ carries a contamination component NOT proportional to any true dark (v4's "aerosol/misfit ~2–3×", now measured per-stream); the far-gate anchor cannot help (the gexp hump is e-suppressed at 10 km). This is the case the plan reserved for M2. |

The CHM15k negative is a *result*, not a failure: the self-test now quantifies, per stream, the
exact thing v4 could only bound globally — and it says the cuvette amplitude needs the network.

**Still open for the next session:** the M3 joint B+T fit on the two long hood sessions; the CL31
sky test (its night cache must be built with the co-located CHM15k's clear-night list — the CL31
has no Rayleigh nights of its own in the archive, which is why v4 called it untestable).

