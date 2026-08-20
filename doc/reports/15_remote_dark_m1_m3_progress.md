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
