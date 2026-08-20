# Remote retrieval of the ceilometer electronic baseline — three candidate methodologies

**Plan only** (2026-08-20). Aim: estimate the mean electronic dark structure `b(z)` of every network
instrument **without a site visit**, verified against the Payerne termination-hood measurements.
Nothing below is implemented.

---

## 0. Problem statement, and what "verified" means

The target is the **mean additive baseline** `b(z)` (rcs_0 units) — the thing a terminal-hood
session measures — not the noise *variance* σ(z), which is already solved (the per-gate
`σ² = β + αP + γP²` intercept method, validated 1–4 % against the hood; that machinery is an asset
here, not a goal).

Per type, the shapes are known from the Payerne hood library and they are *very* different — which
is why one method will not fit all:

| type | hood-measured dark structure | difficulty class |
|---|---|---|
| CL31 | **two under-damped resonances** (AC-coupling ring ~142 kHz + transmitter ripple ~30 kHz, R²=0.98), large near range (−14…−60 % of signal at 450–600 m), decaying with range; plus a gate-locked far ripple | mixed: high-frequency parts easy, damped near-range part hard |
| CL51 | CL31-like ripple family (period 40–80 m, unit-specific) | as CL31 |
| CHM15k | **smooth negative-exponential "cuvette"** — the worst case: nearly collinear with the mean molecular shape (the proven identifiability wall) | hardest |
| CL61 | small and stable: above ~350 m ±1.5·10⁻¹⁴ m⁻¹sr⁻¹ around zero (Looschelders 2025, 6 units × 16 hood sessions, type-common); structure concentrated <350 m (120 m protrusion, afterpulsing rise <1.2 km — Le & O'Connor); slight T-dependence (15–31 °C); our v4 retrieves the smooth part (θ_corr = 0.96 vs hood) | easiest above 350 m; near range is the open part |

**Verification standard** (the user's requirement): the Payerne hood library —
`cl31_b_dark.npz`, `chm15k_b_dark.npz`, the CL61 npz; four CL31 sessions (12 May, **26–27 May
~25 h**, 09 Jun, 23 Jun 2026). Protocol in §5: methods are frozen on synthetic + non-Payerne data
first, then compared to the hood **once**, blind. The ~25 h session is special: it spans a full
diurnal cycle *under the hood*, so it contains the dark's own temperature and daylight-background
modulation — the ground truth for Methodology 3.

## 1. What we already have (assets and proven dead ends)

Assets, all reusable:

* **v4 clear-sky estimator** (`rayleigh_availability/dark_from_clearsky.py`): per-gate robust line
  across nights `S_n(z) = A_n·m_n(z)·g(z) + B_n·a_n(z) + d_n·z² + b(z)` with CAMS molecular + CAMS
  aerosol regressors; injection **self-test** that measures, per stream, how much of a given dark
  shape is recoverable — the honesty gauge every methodology below must pass through.
  Validated: Payerne CL61 θ_corr = 0.96; CHM15k detected at ~3σ but +50 % amplitude; free-profile
  amplitude unusable (aerosol/misfit contamination 2–3×).
* **The identifiability wall, proven**: any dark component collinear with the mean molecular shape
  is invisible *by principle* to single-station night-to-night regression — a uniform shift of the
  nightly amplitudes slides along the regression lines. CHM15k's cuvette is the worst case.
* **Gate-fold ripple retrieval** (`validation/paper/_offset_lib.py`): split-half reproducibility +
  gate-locked folding recovers the periodic component network-wide (CL51 4/10, CL31 6/10 units) —
  but *only* the periodic component.
* **Noise-variance intercept method** (overlap work): per-gate electronic σ, validated vs hood.
* **Cloud calibration is dark-immune** (measured: +0.08 % CL31, ~0.001 % CL61) — an *absolute*
  amplitude anchor that exists precisely for the instruments whose Rayleigh path is compromised.
* **Natural experiments**: unit swaps with L1/L2 serial attributes (Payerne CL31 Feb-2025, CHM15k
  Apr-2025 already characterised); station elevations spanning ~0–2000 m ASL.

Dead ends, proven internally — do not revisit:

* **Fourier in range** for the smooth components: 99–100 % of the Rayleigh-window bias integral
  lives in the smooth band, inseparable from molecular. FT remains the right tool for the
  gate-locked ripple only.
* **Plain temporal PCA/EOF** on night profiles: free loadings absorb the static dark
  (λ′ 0.71 → 0.33). PCA *as such* cannot break a static-vs-static degeneracy.
* **Amplitude screens** on the nightly A_n: they remove the bright nights that anchor the
  regression and degrade the detector.
* Overlap-probe-style window accumulation: yield-dead (~2 usable nights per multi-year archive).

### 1bis. Direct answers to the questions asked

**Is Fourier the right tool?** Only for the gate-locked ripple (40–80 m periods — already solved).
For everything smooth, no — and the proposed fallback "aerosol vs molecular+noise" does not hold
either: aerosol is as smooth as molecular in range-frequency, so FT separates neither pair. What
separates aerosol from molecular is not frequency but **forward knowledge and dynamics**: molecular
is computable per night from CAMS T/p (shape known, amplitude tied to C_L), aerosol varies
night-to-night with its own CAMS regressor and is non-negative. That is regression structure, not
spectral structure.

**Is PCA the right tool?** Not by itself — proven above; the degeneracy is *in the data*, and no
rotation of a degenerate basis removes it. PCA becomes powerful the moment an external structure
breaks the symmetry first, and each methodology below is one way of doing that: across **units**
(M2 — the dark lives in gate space, the atmosphere in altitude space), across **modulators**
(M3 — CompCor-style PCA restricted to signal-free gates, then projected back), or after an
**absolute anchor** removes the amplitude freedom (M1). The astronomy literature reached the same
conclusion from the same failure: KLIP/ADI subtracts a quasi-static speckle pattern only because
the *scene rotates* while the pattern stays fixed; without that motion it self-subtracts, exactly
our wall.

## 2. The physical levers available (inventory the methodologies draw from)

1. **Absolute amplitude anchor**: cloud-calibrated C_L (dark-immune) × CAMS molecular predicts the
   absolute clear-night signal → the residual is dark + aerosol + response error. Kills the
   "uniform A_n shift" freedom for CL31/CL51/CL61.
2. **Gate space vs altitude space**: electronics are locked to time-after-trigger (gate index);
   the atmosphere is locked to altitude. Station elevation (0–2000 m), tilt, and range-resolution
   settings shift the atmosphere in gate coordinates while the dark stays put.
3. **Signal-free gates**: beyond the SNR horizon (CL31 ≳ 6 km at night; CHM15k ≳ 9–10 km at
   1064 nm) the mean clear-night profile *is* the dark plus firmware-baseline residual — a direct
   measurement of `b(z)` at far range, needing only extrapolation inward.
4. **Type-common parametric form**: the hood library gives the *model family* per type (two
   resonances; negative exponential; growing baseline) — a handful of parameters per unit instead
   of a free profile, and a family that is largely *not* collinear with molecular.
5. **Instrumental modulators the atmosphere does not share**: solar background (sweeps orders of
   magnitude at twilight while the atmosphere is frozen; drives the firmware baseline error),
   internal temperature (within-night, guarding the seasonal confound), laser energy (shared with
   the signal — usable only on band-passed ripple components).
6. **Unit swaps**: before/after residual difference at a fixed site = difference of two units'
   darks, seasonally adjusted — absolute-free *pairwise* constraints.

---

## 3. The three methodologies

### M1 — Physics-anchored parametric retrieval (« model the electronics, not the profile »)

**Idea.** Stop estimating a free `b(z)` (2 000 unknowns, degenerate) and fit the **type-specific
electronic model** (3–8 parameters per unit) under multiple weak but *absolute* constraints:

* model family per type from the hood physics: CL31/CL51 `b(z) = Σ_{k=1,2} A_k e^{−z/τ_k}
  sin(2πf_k z/c + φ_k)` (priors on f_k from Payerne: ~142 kHz and ~30 kHz), CHM15k
  `b(z) = −A e^{−z/τ}` (+ the cuvette curvature term), CL61 smooth monotone spline with a
  growth prior;
* **anchor 1 — cloud C_L** (CL31/CL51/CL61): on clear nights, predicted absolute molecular signal
  = C_L,cloud × m_n(z); the mid-range residual constrains the dark where the signal still
  dominates. Converts the fatal "dark vs amplitude" degeneracy into a benign "dark vs
  multiplicative response g(z)" one — benign because g is smooth-multiplicative while the CL31
  family is oscillatory, and because the parametric form removes the per-gate freedom;
* **anchor 2 — signal-free far gates**: direct `b(z)` beyond the SNR horizon (lever 3), the
  strongest constraint on the exponential tails. Not a fixed range cut: adopt Le & O'Connor's
  stationary-wavelet noise-gate identification (SWT bior1.1, minimax threshold) to find the
  aerosol/hydrometeor-free gates per profile — published, tested on CL61, directly portable;
* **anchor 3 — gate-fold ripple** (existing): fixes the periodic component exactly;
* fit per unit by robust weighted least squares across the clear-night ensemble, the v4 machinery
  providing the per-night CAMS regressors and caching.

**Why it can break the wall.** The wall applies to a *free* per-gate baseline. A two-resonance
CL31 model is nearly orthogonal to molecular by construction; the CHM15k exponential is the honest
test — its collinearity with molecular is *quantified per unit* by the existing injection
self-test, and the far-gate anchor pins the tail that the regression cannot see.

**Per type.** CL31/CL51: complete solution plausible (all three anchors live). CL61: already
solved by v4; M1 adds the parametric smoothing. CHM15k: partial — no cloud anchor; relies on far
gates + form + self-test-corrected θ; if the self-test says the cuvette stays >50 % invisible,
M1 alone will not close CHM15k and M2 takes over.

**Verification.** Blind Payerne A/B/C; gate-by-gate comparison to the hood npz over 0.3–10 km;
success = CL31 near-range structure recovered within ~20 % of the −14…−60 % hood effect, CHM15k
window-integrated Rayleigh bias predicted within a factor ~1.5.

**Risks.** g(z)-vs-dark leakage at 910 nm if the WV correction is imperfect (mitigate: fit 1064
CHM15k with no WV term as the control); firmware baseline subtraction making "signal-free" gates
not exactly dark (characterise on the hood data first — the hood profiles contain the same
firmware behaviour).

**Effort.** ~1–2 weeks; reuses v4, the cloud constants, `_offset_lib`, the hood library.

---

### M2 — « The network is the hood »: joint cross-unit decomposition

**Idea.** The degeneracy is unbreakable per station but *not across stations*. For one type,
model all units jointly:

```
S_{u,n}(z) = A_{u,n}·m_{u,n}(z)·g_u(z) + B_{u,n}·a_{u,n}(z) + d_{u,n}·z² + s_u·D_type(z) + r_u(z)
```

with `D_type(z)` a **shared** dark basis (1–3 components), `s_u` unit amplitudes, `r_u(z)` a small
unit-specific remainder (the ripple, from anchor 3). Identification comes from lever 2: across
units, `m_{u,n}` shifts and rescales in gate coordinates (station elevation 0–2000 m, climate,
pressure) while `D_type` does not move — the network plays the role the hood plays at one site.
Unit-swap events (lever 6) add absolute pairwise differences; the Payerne units inject the anchor
that converts the shared basis from relative to absolute.

Two implementations, in order:

1. *exploratory* — pooled robust PCA on the per-unit v4 residuals (after CAMS forward removal),
   [unit-night × gate]: the leading components that are common in gate space and do **not**
   regress on any atmospheric covariate are the candidate dark modes; directly comparable to the
   hood templates. This is the defensible version of "use PCA".
2. *confirmatory* — the hierarchical fit above (alternating robust LS), with a **network-scale
   injection self-test**: inject hood-template darks of known amplitude into a random subset of
   units and require unbiased recovery before believing anything.

**Why it can break the wall.** A CHM15k at Payerne (491 m), at sea level, and at 1 600 m see
*different* molecular shapes in gate space; their common gate-locked residual cannot be atmosphere.
This is the one methodology whose information content grows with the network (~140 CHM15k, ~90
CL31, ~50 CL51, ~60 CL61) rather than with one site's archive.

**Per type.** The CHM15k answer — the type M1 cannot fully close, and the type with the largest
fleet. Also the natural frame for CL51/CL31 unit-amplitude screening (which units need attention).
The shared-basis assumption now has direct published support for CL61: Looschelders et al. (2025)
find six co-located units × 16 hood sessions give *similar* dark profiles across instruments above
40 m — one shared `D_type(z)` with unit amplitudes is what their data shows. For CL31 the Kotthaus
"sensor-specific frequency" means the ripple parameters stay per-unit; the shared basis carries
only the smooth components.

**Verification.** Hold Payerne out of the fit entirely; predict its units' darks from the network
solution; compare to the hood. Sanity anchors from the v4 scan: Messina must come out ~5× Payerne,
Montsec/Twenthe ~0.

**Risks.** Site systematics that are *also* gate-locked (firmware version families, range-gate
settings) masquerading as dark — mitigate by fitting firmware/hardware cohorts separately and by
the injection test; elevation leverage concentrated in few mountain stations (check the design
matrix conditioning before trusting the result).

**Effort.** ~2–3 weeks; needs one balfrin sweep to build the per-unit residual cache network-wide
(the v4 cache machinery exists), then the joint fit is cheap.

---

### M3 — Modulation and covariance: catch the dark *moving*, not sitting

**Idea.** The static mean is the hard part; the dark's **modulated** parts are separable because
their drivers are instrumental, not atmospheric (lever 5). Estimate each driven component from its
driver, then reassemble:

* **Twilight background sweep** (primary): at dawn/dusk the solar background sweeps orders of
  magnitude in <1 h while the atmosphere is essentially frozen. Per gate, regress the signal on
  the measured background around the terminator → `D_bg(z)`, the background-driven baseline error
  (the CHM15k firmware over-subtraction is exactly this). Twice a day, every station, no clear
  night needed.
* **CompCor-style noise-ROI PCA** (the second defensible PCA): PCA of the *time series* of the
  signal-free far-gate block → temporal loadings of the electronic fluctuation modes; regress
  every gate on these loadings → full-range spatial patterns of everything that co-fluctuates
  with the electronics. Recovers the shapes of all modulated components at once, including ones
  without a housekeeping proxy.
* **Band-passed laser-energy regression**: the transmitter ripple shares its driver with the
  atmospheric signal (both ∝ pulse energy), so regress only the 40–80 m band-passed signal on the
  energy proxy (laser status/lifetime housekeeping) — in that band the atmosphere is quiet.
* **Within-night temperature regression** for the T-dependent part, with day-of-year partialled
  out (the seasonal-confound guard from the overlap work). Published motivation now exists: FMI
  correct the CL61 with a hood-measured `P_instrument(r, T)` **lookup table by internal
  temperature** (Le & O'Connor §3.1.3) — M3's ambition is that table *without* the hood.
* **Laser-power epochs** (new lever, from Le & O'Connor §4.1): background noise anti-tracks the
  laser-power decline — stepwise, jumping at transmitter swaps, compensation ceasing below ~40 %
  power. Laser power is in the housekeeping; segment every retrieval by power epoch and regress
  the noise-normalised components on it. An aging laser is a slow modulator sweeping the dark's
  amplitude over months — usable like the twilight sweep, on another timescale.

One caveat both papers force on the twilight sweep: CL31 and CL61 firmware applies a **zero-bias
compensation** for the solar pedestal (Kotthaus 2016; confirmed for CL61 by Le & O'Connor — the
mean above 5 km follows β_mol). The sweep therefore retrieves the *residual* background error
`B_bk` after that compensation — precisely the quantity that pollutes the data, but not to be
read as the raw solar pedestal.

**What it cannot do, honestly.** The fully static, never-modulated remainder is invisible to M3 by
construction. M3 is not a stand-alone solution; it is (i) the only route to the
*background-dependent* dark — which the clear-night methods never see, since they select for
darkness, yet which matters for every daytime product — and (ii) a set of independent components
that reduce what M1/M2 must explain.

**Per type.** CHM15k: `D_bg` is expected to be the dominant driven term (baseline subtraction
error). CL31/CL51: the energy-band ripple + background pedestal. CL61: small everything (its dark
is small and smooth) — mostly a null check.

**Verification — the 25 h hood session is the goldmine.** That session spans a full diurnal cycle
*under the hood*: the dark's own background and temperature response is in the truth data. Freeze
the M3 regressions on sky data, then compare `D_bg(z)` and the T-slope directly to what the hood
session shows. No other methodology has ground truth for its *dynamics*; M3 does.

**Risks.** Twilight is not perfectly frozen atmosphere (BL transition at dawn) — mitigate by using
both terminators, high gates, and the falling/rising asymmetry; housekeeping proxies coarse on
some firmware.

**Effort.** ~1 week for the twilight + CompCor core; the 25 h hood validation is a day.

---

## 4. Applicability matrix

| component | CL31/CL51 | CHM15k | CL61 |
|---|---|---|---|
| gate-locked ripple | solved (gate-fold) | n/a | n/a |
| near-range damped structures | **M1** (anchored resonances) + M3 energy band | — | — |
| smooth collinear pedestal (cuvette) | — | **M2** (network) + M1 far-gate/form | v4 already (θ=0.96) |
| background-driven baseline | M3 twilight | **M3 twilight** (dominant) | M3 (null check) |
| absolute amplitude | cloud C_L anchor | network + swaps + far gates | v4 + cloud C_L |

## 5. Common verification protocol (agreed before any implementation)

1. **Freeze first**: each method is developed on synthetic injections + non-Payerne streams; the
   Payerne hood comparison is run once, blind, per method.
2. **Injection self-test everywhere**: no retrieved amplitude is reported without its per-stream
   recoverability factor (the v4 affine correction generalised to M1/M2/M3).
3. **Metrics**: gate-wise agreement with the hood npz (0.3–10 km, per band); window-integrated
   predicted Rayleigh bias (what the calibration actually feels); for M3, dynamic agreement
   against the 25 h session. For CL61 the discriminating band is **below 350 m** (the 120 m
   protrusion, the afterpulsing rise) — above it the dark is ±1.5·10⁻¹⁴ m⁻¹sr⁻¹ and any method
   "succeeds" trivially. Hood-protocol detail when comparing: discard the first ~20 min of a
   session (near-range transient — Kotthaus 2016, applied by Looschelders 2025).
4. **Network sanity**: Messina ≈ 5× Payerne, Montsec/Twenthe ≈ 0, and the doctrine stands — only
   subtract what is *measured* (hood) *or* detected strongly (θ > 0.5); weak-θ subtraction was
   shown to improve nothing.

## 6. Recommendation and sequencing

Run **M1 and M3 in parallel first** (they are cheap, independent, and each provides anchors the
other lacks: M1 the static amplitude, M3 the driven components + the only daytime-relevant term),
with the blind Payerne comparison at the end of both. **M2 second**, seeded by whatever M1 leaves
unexplained on CHM15k — it is the only route through the cuvette wall and the only one that scales
to "all instruments" by construction, but it is also the most expensive to validate honestly.

The end product is a fused per-unit `b̂(z)` with a per-component provenance tag (measured ripple /
anchored parametric / network-shared / driven) and a recoverability factor — the same honesty
structure the v4 θ already has, extended to the full profile.

---

### Literature notes (checked 2026-08-20)

**The two closest works, read in full (both in `doc/`):**

* **Le & O'Connor et al., egusphere-2025-6331** (CL61 operational performance, FMI): the state of
  the art for *hood-based* correction — `P_instrument(r, T)` measured under the hood at a range of
  internal temperatures over two years, stored as a lookup table, subtracted by current T
  (their Eqs. 12–15), solar variance separated at 10–12 km (Eq. 11). Also: the SWT noise-gate
  identification (§3.1, adopted by our M1/M3); background noise anti-tracks laser power with
  steps at transmitter swaps (§4.1 — our M3 laser-power lever); afterpulsing rise below 1.2 km;
  firmware zero-bias compensation confirmed on CL61. **Our aim restated in their terms: obtain
  their lookup table without their hood.**
* **Looschelders et al. 2025, Met. Apps** (6 co-located CL61, SIRTA): the CL61 dark is small
  above 350 m (±1.5·10⁻¹⁴ m⁻¹sr⁻¹), structured below (120 m protrusion, ppol > xpol), *stable in
  time and similar across all six units* (16 hood sessions), slight T-dependence — direct
  evidence for M2's shared-basis assumption on CL61, and the definition of where the CL61
  problem actually lives (<350 m). Their SNR filter explicitly *assumes noise constant with
  height and neglects height-dependent electronic noise* — the gap our variance-intercept
  method already fills.

No published *remote* dark retrieval exists; the field measures it (termination hood) and
corrects what it can reach: [Kotthaus et al. 2016, AMT](https://amt.copernicus.org/articles/9/3769/2016/)
(CL31 background procedure and the "sensor-specific frequency" ripple — our gate-fold is its
network generalisation), the [ARM ceilometer handbook](https://www.arm.gov/publications/tech_reports/handbooks/ceil_handbook.pdf)
(hood as an optional accessory; CL61 median dark profile from hood sessions). The transferable
methods are extramural: exoplanet high-contrast imaging
([KLIP/ADI PSF subtraction](https://iopscience.iop.org/article/10.3847/2515-5172/aa9d18),
[robust local-weighting variants](https://www.aanda.org/articles/aa/full_html/2020/02/aa35859-19/aa35859-19.html))
— a quasi-static pattern subtracted only because the scene *moves* while the pattern stays, with
self-subtraction as the failure mode, i.e. our identifiability wall under another name — and
[camera fixed-pattern-noise estimation](https://piscat.readthedocs.io/Tutorial3/Tutorial3.html)
(M2 is the "many detectors, one pattern" variant; M3's CompCor is the fMRI noise-ROI variant).
ICA/NMF blind-source-separation surveys were checked and offer nothing beyond what the anchored
regressions already encode — with only one quasi-static mixture per station, blind separation has
no diversity to work with; every workable route above supplies that diversity physically
(anchors, the network, modulators).
