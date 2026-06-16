# Water-vapor correction — code audit & validation

**Date:** 2026-06-15
**Scope:** cross-validation and bug audit of the water-vapor (WV) two-way transmission
used for the Rayleigh calibration of 910 nm ALCs.

**Codes reviewed**

| | File | Role |
|---|---|---|
| Python (new) | `rayleigh_calibration/water_vapor.py` | the deliverable under test |
| MATLAB | `liquid_cloud_calibration/wv_t2eff.m` | Gaussian-weighted T² core |
| MATLAB | `compute_wv_transmission.m` | full pipeline (CAMS → T²(z,t)) |
| MATLAB | `Water_absorption_correction/get_water_vapor_number_concentration_from_RH.m` | RH → n_wv |
| MATLAB (dep.) | `MDA/monitoring_alc_monthly/get_Beta_CAMS_oper_monthly.m` | CAMS read + L137 geopotential |
| MATLAB (dep.) | `Basics functions/convert_humidity/convert_humidity.m` + `calculate_saturation_vapor_pressure_liquid.m` | q → RH |

**External reference:** ACTRIS-Cloudnet `atmoslib` (thermodynamics). See §6.

---

## 1. Method

Three layers of validation, all in `tests/test_water_vapor.py` (12 tests, all passing):

1. **vs MATLAB, hardcoded references** — the MATLAB routines were run on identical
   synthetic and real inputs (driver `tmp_wv_ref.m`) and their outputs frozen into the
   tests. This pins the Python to the validated MATLAB to ≤2·10⁻⁴ on the T² core and
   ≤1.2 % on the full CAMS→T² chain.
2. **vs an independent external library** — `atmoslib` (ACTRIS-Cloudnet) for the
   thermodynamic primitives (vapour pressure, saturation pressure, absolute humidity).
3. **internal invariants** — bounds (0<T²≤1), monotonicity, zero-absorption→unity, and a
   **regression test for the one latent bug found** (§5).

---

## 2. Verdict

The Python port is **correct and, at the CAMS→n_wv step, slightly more accurate than the
MATLAB** it replaces. The two T²-core implementations are line-for-line equivalent; the
humidity path differs by design (Python is the cleaner route). One genuine latent bug,
shared by both codes, was found and fixed in the Python (flagged for MATLAB).

| ID | Where | Severity | Status |
|---|---|---|---|
| **F3** | both (`cams_water_vapor_profile` / `get_Beta_CAMS_oper_monthly`) | latent bug | **fixed in Python**, flagged for MATLAB |
| F1 | MATLAB only | ≤0.3 % bias (cold), ~0.01 % in warm BL | documented (Python avoids it) |
| F2 | MATLAB only | +0.097 % bias | documented (Python avoids it) |
| F4 | Python only | benign safety net | keep |
| F5a/b | both | shared modelling approximation | documented, acceptable |

Net systematic Python−MATLAB difference in the boundary layer where WV matters: **~0.1 %
on n_wv, <0.05 % on T²**. The ~1.2 % spread seen in the full-chain test is dominated by
grid/interpolation/time-averaging choices, not by these constants — i.e. the codes
genuinely agree.

---

## 3. Confirmed correct (faithful ports)

- **`wv_t2eff_core` ≡ `wv_t2eff.m`.** Same σ=FWHM/(2√(2ln2)), same Gaussian, same ±3σ
  band mask, same `ext = abscs·n_wv/1e4` (cm²→m²), same `cumtrapz` over range, same
  `T² = Σ T²·g / Σ g`. Verified to atol 2·10⁻⁴ on two synthetic cases
  (`test_core_vs_matlab_constant`, `test_core_vs_matlab_varying`).
- **L137 geopotential ≡ `get_Beta_CAMS_oper_monthly.m`.** Same Rd=287.06, same a/b
  half-level coefficients, same `T_moist=T·(1+0.609133·q)`, same
  `z_f=z_h+T·Rd·α`, `z_h+=T·Rd·dlogP` recursion, same /9.80665.
- **q→e identical.** `convert_humidity` uses `c = M_wet/M_dry = 18.0152/28.9644 =
  0.621981`, which is **exactly** the Python `EPS`; the e=(q·P)/(c+(1−c)q) formula is the
  same. So the two codes diverge only *after* e (see F1).
- **Full Payerne CL61 chain** (CAMS q → n_wv → T²) matches `compute_wv_transmission.m` to
  ≤1.2 % at 1/3/6 km (`test_full_profile_vs_matlab_payerne`).

---

## 4. Findings

### F3 — latent bug (both codes): top-of-atmosphere guard assumes ascending level order
**Location:** `cams_water_vapor_profile` (Python), `get_Beta_CAMS_oper_monthly.m:259`
(MATLAB, `if i == 1`).

The hydrostatic integration walks from the surface (last array index) upward, and the
top-of-atmosphere singularity (upper half-level pressure = 0 ⇒ `log(0)`) is handled by a
special case. **Both codes keyed that special case off the loop position** (`i==0` Python,
`i==1` MATLAB) instead of the physical top half-level. This is correct *only* if the
`level` axis is stored ascending (1=top … 137=surface), which today's CAMS files happen to
be. On a reordered or subset axis the guard fires on the wrong level — dividing by a
nonzero `Ph_lev` at the true top (→ `log(p/0)=∞`) or applying the 0.1 Pa replacement at a
spurious level — corrupting the whole geopotential profile (the mutation test produced a
surface "height" of 4955 m instead of 501 m).

**Fix (Python):** sort all level-indexed arrays ascending on entry (no-op for today's
files) **and** key the special case on `idx == 0` (the half-level whose pressure is
identically zero), not on loop position. Locked by `test_geopotential_invariant_to_level_order`,
which builds ascending and descending synthetic CAMS files and asserts identical output —
and which **fails** if either half of the fix is removed (verified by mutation).

**MATLAB:** the same fragility exists in `get_Beta_CAMS_oper_monthly.m`. It is not on the
critical path for current files, so it is *flagged, not changed* (that routine is shared by
other pipelines; change only with the user's go-ahead). Suggested one-line hardening: sort
by `level` ascending before the `t` loop, or guard on `Ph_lev==0` instead of `i==1`.

### F1 — MATLAB only: saturation-formula round-trip (q→RH→Pw uses two different formulas)
The MATLAB derives n_wv as q →(`convert_humidity`)→ RH →(`get_water_vapor…RH`)→ Pw. But
`convert_humidity` uses **Murphy & Koop 2005** for saturation while
`get_water_vapor_number_concentration_from_RH` uses **Wagner–Pruß IAPWS-95**. The two do
not perfectly cancel, leaving a residual factor `Pw_final/e_true = Pws_WP/es_MK`:

| T [K] | 240 | 250 | 273.15 | 290 | 300 |
|---|---|---|---|---|---|
| bias | +0.32 % | +0.11 % | 0.00 % | −0.011 % | −0.001 % |

So **~0.01 % in the warm boundary layer** (where WV absorption dominates) and ≤0.3 % only
at cold upper levels where WV→0. Real but negligible. **The Python avoids it entirely** by
going straight q → e = (q·P)/(ε+(1−ε)q) → n_wv = e/(k_B·T); no saturation formula and no
phase ambiguity (the MATLAB also silently uses over-*liquid* saturation at all
temperatures — its over-ice branch is dead code because T is in Kelvin).

### F2 — MATLAB only: number-density coefficient 7.25e22 vs 1/k_B
`get_water_vapor_number_concentration_from_RH` ends with `nw = 7.25e22·Pw[Pa]/T` (the Rw
factors cancel). The exact value is `1/k_B = 7.24297e22`, so the MATLAB n_wv is **+0.097 %
high**, flat. The Python uses `KB = 1.380649e-23` exactly. Captured by
`test_n_wv_constant_matches_matlab` (rel 2·10⁻³).

### F4 — Python only (benign): T² clipped to 1 on degenerate input
`wv_t2eff_core` sets `T²≤0` or non-finite to 1.0 (= no correction). This mirrors the final
clip in `compute_wv_transmission.m:141` (absent from `wv_t2eff.m`, which clips in its
caller). It only fires on degenerate input and is a safe fallback. Minor asymmetry: it does
not clip a hypothetical `T²>1` — but a positive WV profile can never produce that, and both
codes share the gap. Keep.

### F5 — shared modelling approximations (not bugs; note in the paper)
- **(a) Out-of-band weight.** Wavelengths with |λ−λ₀|>3σ but still inside the LUT get
  T²=1 (no absorption) yet keep their Gaussian weight in the normalised average. The mass
  beyond ±3σ is <0.27 %, so T² is biased high by <0.3 %. Identical in both codes.
- **(b) Cross-section height mapping.** LUT cross-sections are mapped to range gates by
  nearest **AGL** height, and their pressure/temperature broadening is the LUT's built-in
  standard atmosphere, not the actual CAMS p,T. Weak effect; identical in both codes.
- **Design difference (not a discrepancy):** the Python returns one **time-averaged** T²
  profile per calibration window (correct for a night-mean Rayleigh fit), whereas
  `compute_wv_transmission.m` returns T²(z,t) resolved per profile (for L2 products). Same
  physics, different temporal granularity.

---

## 5. Fix applied

Only **F3** warranted a code change (it is the sole genuine bug). Applied to
`rayleigh_calibration/water_vapor.py`:

1. sort `level`, `T`, `q` ascending on entry to `cams_water_vapor_profile`;
2. key the top-of-atmosphere special case on `idx == 0` (zero-pressure half-level).

Regression test `test_geopotential_invariant_to_level_order` builds ascending + descending
synthetic CAMS files and asserts byte-identical profiles; mutation-tested to confirm it
goes red without the fix. F1/F2 are MATLAB-side and need no Python change; they are recorded
here and the Python already does the right thing.

---

## 6. External reference (the "find a repository online" request)

- **ACTRIS-Cloudnet `atmoslib`** (the cloud-radar/lidar processing org's thermodynamics
  library, v2.4.1) is used as an *independent* check of the humidity primitives:
  `vapor_pressure(P,q)`, `saturation_vapor_pressure(T)`, `absolute_humidity(T,e)`. Three
  tests cross-check our Pw and n_wv against it (`test_vapor_pressure_matches_atmoslib` to
  ~2 ppm — the only gap is atmoslib's slightly different M_w/M_d constant;
  `test_n_wv_matches_atmoslib`; `test_n_wv_from_RH_matches_matlab` ties atmoslib's
  saturation pressure to the MATLAB IAPWS reference values to 1 %).
- **CloudnetPy itself does *not* implement a Wiegner-style spectral WV transmission
  correction.** Cloudnet calibrates ceilometers against liquid-cloud returns
  (O'Connor 2004), and its 910 nm instruments are not WV-corrected the way we do here. So
  CloudnetPy is *not* a drop-in reference for `water_vapor.py`; `atmoslib` is the right,
  authoritative external anchor for the thermodynamic core, and the MATLAB
  (Wiegner & Gasteiger 2015) remains the reference for the spectral transmission itself.

---

## 7. Test inventory (`tests/test_water_vapor.py`, 12 tests)

| Test | Validates against | Checks |
|---|---|---|
| `test_core_vs_matlab_constant` | MATLAB `wv_t2eff.m` | T² core, constant inputs (atol 2e-4) |
| `test_core_vs_matlab_varying` | MATLAB `wv_t2eff.m` | T² core, λ-varying abscs + decaying n_wv |
| `test_core_bounds_and_monotonic` | invariant | 0<T²≤1, non-increasing |
| `test_no_absorption_gives_unity` | invariant | abscs=0 ⇒ T²=1 |
| `test_n_wv_constant_matches_matlab` | MATLAB convention | 1/k_B ≈ 7.25e22 (F2) |
| `test_vapor_pressure_matches_atmoslib` | atmoslib | Pw(q,P) to ~2 ppm |
| `test_n_wv_matches_atmoslib` | atmoslib | n_wv via absolute humidity (0.2 %) |
| `test_n_wv_from_RH_matches_matlab` | atmoslib + MATLAB | RH→n_wv (1 %) |
| `test_band_and_laser_spectrum` | spec | band membership, λ₀/FWHM table |
| `test_full_profile_vs_matlab_payerne` | MATLAB `compute_wv_transmission.m` | full CAMS→T² chain (≤1.2 %) |
| `test_lut_sanity` | data | HITRAN LUT shape/coverage |
| `test_geopotential_invariant_to_level_order` | invariant | **F3 regression** (level-axis order) |

Run: `python -m pytest tests/test_water_vapor.py -v -o addopts=""`
(`-o addopts=""` disables the repo's default `--cov`, which needs pytest-cov.)
Data-dependent tests (CAMS/LUT) skip cleanly when those files are absent; the atmoslib
tests skip if the library is not installed.
