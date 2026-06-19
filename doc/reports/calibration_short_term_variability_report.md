# Reducing the short-term variability of ALC calibration — diagnosis & recommendations

**Author:** M. Hervo · **Date:** 2026-06-17
**Scope:** Rayleigh + liquid-cloud calibration; CL31, CL51, CL61, CHM15k.

## Question & short answer

The daily/per-night calibration coefficients scatter a lot. **Why, and how do we reduce it?**

**Most of the short-term scatter is averageable measurement noise**, not real calibration
change. It is largest where the calibration signal is weakest (the per-night **Rayleigh**
molecular fit: ~**39 %** night-to-night), and smallest for the strong-signal **liquid-cloud**
target (~**8 %** day-to-day). Underneath the noise sits a **real, non-averageable floor
(~8–20 %)** that is the *instrumental* drift identified separately (laser ageing, internal
temperature, window) — see `calibration_stability_report.md`.

**We do NOT want to average over months — that would hide the instrument changes (window
cleaning, laser swaps, drift) the Kalman is meant to detect.** The right target is the
**shortest integration that yields one usable calibration point** at a useful precision; those
points are then fed to the **Kalman, which adaptively tracks drift and detects steps**. That
single-calibration integration time is **short and instrument-dependent**: ~**1 day** for the
liquid-cloud method (strong signal, 8–9 %), but ~**2–3 weeks** for the per-night Rayleigh
(weak molecular signal, 39 %/night). Vertical averaging (wider molecular fit window) gives a
secondary ×2 on the per-night Rayleigh noise; quality/SNR filtering mainly throws nights away
without improving the per-estimate noise. **§2 below is the operative section.**

## 1. Why is the short-term scatter so large? — it is mostly noise

The Allan deviation σ_A(τ) of each instrument's daily calibration separates **random noise**
(σ_A falls as τ^(−1/2): averaging helps) from **real drift** (σ_A rises with τ: averaging
hurts). The minimum marks the optimal averaging time.

![Allan deviation of the calibration vs integration time](figs_paper_validation/calib_integration_allan.png)

*(a) Allan deviation vs integration time τ (log-log); the dotted line is the pure-white-noise
τ^(−1/2) slope; the pentagram marks each series' optimum. (b) the short-term scatter a W-day
running median removes (rises to the total short-term variability).*

| Series | λ / method | σ_A at 1 day | optimal τ | noise floor at optimum | reduction |
|---|---|---|---|---|---|
| **CL61 Rayleigh** | 910, molecular | **39 %** | ~21 d | 9 % | **×4.6** |
| CL31 cloud | 910, cloud | 31 % | ~60 d | 19 % | ×1.6 |
| CHM15k LIN | 1064, molecular | 15 % | ~1 yr | 8 % | ×1.9 |
| CHM15k PAY | 1064, molecular | 13 % | ~21 d | 11 % | ×1.2 |
| CL51 cloud | 910, cloud | 9 % | ~90 d | 2 % | ×4.5 |
| CL61 cloud | 910, cloud | 8 % | ~30 d | 3 % | ×2.3 |

Two clear patterns:
- **Per-night Rayleigh is photon/SNR-limited** (39 % at 1 day): the molecular return at 3–6 km
  is ~100× weaker than a cloud target, so a single night is very noisy — but it follows the
  τ^(−1/2) line, i.e. it is **almost pure averageable noise** (×4.6 down by ~3 weeks).
- **The liquid-cloud target is intrinsically precise** (8–9 % at 1 day) because the in-cloud
  signal is strong; it still benefits from ~1-month averaging.
- **CHM15k PAY barely averages down (×1.2)**: its 1-day scatter is *already* dominated by the
  real instrumental drift floor (~11 %), not noise — consistent with it being a stable, mature
  unit whose variability is seasonal/temperature, not random.

## 2. Single-calibration integration time + Kalman change detection (the operative answer)

We want the **shortest integration that gives one usable calibration point**, then let the
**Kalman detect changes** — not a months-long average that would smear them.

![single-calibration integration time and change detection](figs_paper_validation/calib_single_time_changedetect.png)

*(a) precision of ONE calibration vs its integration time, with 10 % / 5 % targets. (b,c) a
+20 % step (e.g. window cleaned) on a slow drift: the responsive Kalman recovers it in days,
while a 90-day average smears it over months.*

**Shortest integration for one usable calibration point:**

| Instrument / method | precision of 1 day | τ for ≤10 % | τ for ≤5 % | noise floor |
|---|---|---|---|---|
| **CL61 cloud** | **8 %** | **1 day** | ~3–4 weeks | 3 % |
| **CL51 cloud** | **9 %** | **1 day** | ~10 days | 2 % |
| CHM15k (1064, Ray) | 13–15 % | 1–2 days / ~25 d | floor-limited | 9–11 % |
| **CL61 Rayleigh** | 39 % | **~18 days** | unreachable | 9 % |
| CL31 cloud | 31 % | floor-limited (>19 %) | — | 19 % |

So **the cloud method gives a usable single calibration every single day** (8–9 %) — that is
the integration time to use for change monitoring, *not* a monthly average. The per-night
Rayleigh is too noisy for a daily point (39 %) and effectively needs ~2–3 weeks; it is the
slow **absolute anchor**, not the fast change detector.

**The Kalman is the change detector, not an averager.** Feed it the single-calibration points
(daily for cloud, per-night for Rayleigh) with their uncertainties; tune its **measurement
noise = the single-calibration σ** and its **process noise to the expected real change rate**
(large enough to follow steps/drift). In the demo it recovers a 20 % step in **~4 days (cloud)**
/ ~12 days (Rayleigh); a 90-day average needs ~3 months and lags the drift. On the short CL61
records the production Kalman **diverged** for 3 units — guard it (reject non-physical jumps)
or fall back to a short robust running median, but keep the window short.

### 2b. Rayleigh night-averaging window — now darkness-adaptive (SZA), was fixed clock

Rayleigh works only at night, so the **within-night averaging window** is set by how much dark
time is used. The pipeline used a **fixed solar-clock window (solar 20:00→04:00, 8 h)** — it is
*not* darkness-adaptive (despite a memory that commit 6d684d8 added that; 6d684d8 only fixed the
solar-*time* clock, no sun-angle). At Payerne this **wastes ~2–5 h of true dark time in winter**
(13 h dark vs 8 h used) and **includes ~2 h of twilight in summer** (June dark = 5.7 h). Naively
widening the *clock* window backfires (CV 25 %→95 % at 12 h) because it pulls in twilight.

**Implemented** a solar-zenith-angle night selection in `data_loader.filter_time_range`
(`use_sza_night`, `sza_night_threshold`, default 100° = sun ~10° below horizon; clock fallback
when coordinates are missing; `use_sza_night=False` reproduces the old result exactly). On
Payerne CL61 2026 it **recovers more usable nights (25 vs 20)** and correctly tracks the dark
period across seasons/latitude. It does **not lower the per-night CV** (≈25–33 %): that scatter
is night-to-night atmospheric/fit variability, not within-night averaging — so the window is the
*correct* and *adaptive* choice (essential for high-latitude sites and summer, where the fixed
clock fails), but the per-night noise is still beaten only by the Kalman over nights (§2) and by
the wider fit window (§3). Recommended threshold 100° (more nights) to 108° (astronomical,
cleanest fit) — a fit-quality vs night-count trade-off.

## 3. Lever — vertical integration (molecular fit window)

Re-running the CL61 Payerne Rayleigh with a **wider molecular fit window** (2–9 km vs 2–6 km,
plus longer fit half-lengths) cut the per-night CV from **24.9 % → 11.1 % (×2.2)** — more
range gates averaged ⇒ less photon noise. The trade-off: it **also cut the usable nights
(20 → 6)** because a clean molecular signal to 9 km is rarer, and shifted the median ~−4 %
(sampling a higher, more aerosol-free column). Recommendation: a moderate widening (e.g. to
7–7.5 km) is worth testing operationally for the noisy 910 nm Rayleigh; the 1064 nm CHM15k
already has enough SNR that this matters less.

## 4. Lever — noise filtering (SNR / quality / outliers)

| filter | effect on CL61 Rayleigh |
|---|---|
| tighten method-agreement gate (quality 15→8) | nights 20→10, **CV unchanged (24.8 %)** |
| require ≥5 h night (vs 3 h) | no change |

Tighter gates **remove nights without lowering the per-estimate noise** — they improve
*reliability* (reject bad fits) but are not the route to lower variability. The effective
"filter" for variability is **temporal averaging** (§2). The
**housekeeping filters from the literature** (window transmission < 90 % reject, laser power
< 40 % flag — Hopkin 2019, Le 2026) remove the episodic excursions that no amount of averaging
will fix.

### 4a. Outlier *profiles* are the cause of the random-subset spread — screen, don't medianize

The Rayleigh sensitivity test (random 70 % time-subsets) showed big subset-to-subset spread.
Cause: the night's profiles were collapsed with a **mean** (`np.nanmean`), which is sensitive to
a few residual-aerosol / cloud-edge / noisy profiles — so each random subset that includes or
excludes them shifts. The per-altitude CL and the perturbation aggregate already use the median.

Tested on Payerne CL61 2026 (night-to-night CV):

| time collapse | night-to-night CV | note |
|---|---|---|
| mean, no screen (legacy) | 33.3 % | outlier-sensitive |
| **mean + outlier screen (now default)** | **23.8 %** | screen profiles by MAD on the molecular-band signal, then mean |
| median profile | 81.6 % | **worse** — median is ~1.5× noisier at the weak 3–6 km signal |

So the fix is **remove outlier profiles, then take the mean** — *not* a median profile (which
is robust to outliers but far noisier where the molecular signal is weak). **Implemented:**
`screen_profile_outliers` (default on, MAD threshold `profile_outlier_nmad`=4) +
`time_aggregation` ("mean" default; "median" available). This cut the per-night CV ×1.4.

### 4b. Molecular-window detection was selecting aerosol / low-R² windows — diagnosed & **fixed**

**The bug.** The Python window search chose the center (2–6 km) with the smallest **Σ|intercept|**,
then the largest R² at that center, with a **free** intercept and **no R² floor**. This is
degenerate: in the high-altitude noise region the signal → 0, so the regression fits it with
intercept **b ≈ 0** (trivially small) *and* **R² ≈ 0** (no real correlation). Minimising Σ|b|
therefore systematically **selects a high, noise-dominated window with near-zero R²**, and
nothing rejects it. This is exactly what the diagnostic figure for Payerne CL61 **2026-03-12**
showed: the "optimum" sat at center **4.62 km in a near-zero-R² region**, while genuine R²≈1 was
only at low center (where the boundary-layer aerosol makes signal and molecular both decay with
height → spurious correlation).

**The literature names this failure mode.** Mattis, D'Amico, Baars et al. (2016, *AMT* 9, 3009;
EARLINET Single Calculus Chain) state that a pure minimum-signal/minimum-background search
"does **not guarantee that there are no particles** … would find a minimum also in the case that
there are fewer particles than in other altitude regions only … may cause **large errors**." The
robust remedies they and others use are (i) an **SNR / standard-deviation gate**, (ii) a
**Rayleigh-shape ("molecular") test** — require the window's shape to match the computed Rayleigh
profile (Freudenthaler et al. 2018 residual ≤ 1 %; Baars et al. 2016), (iii) a **scattering-ratio
bound** R ≤ 1.1 (Wiegner & Geiß 2012; CALIOP 1.01 ± 0.01), and (iv) **prefer the lowest qualifying
window** (best SNR). R² of the signal-vs-molecular fit is a valid operationalisation of the shape
test — and crucially **R² is the right metric because it *collapses* in the noise region** (so a
noise window is rejected, not selected, the opposite of Σ|b|).

**The original MATLAB had all three guards; the Python port lost them.**
`Auto_Calib_25/Rayleigh/rayleigh_fit.m` (Hervo & Poltera 2014) **forces the intercept to zero**
(`opts.Lower/Upper=[-Inf 0]/[Inf 0]`), chooses the center by **minimum Σ RMSE** (not Σ|b|), and
the caller (`auto_calib_v23/24.m`) **rejects** the night when `best_r2 < min_r2_rfit` (=0.5) or
`best_rmse > max_rmse_rfit`. It also searches centers only to **5000 m** (Python: 6000 m).

| | Auto_Calib_25 (MATLAB, reference) | Python port (before) | Python (now, fixed) |
|---|---|---|---|
| Intercept | forced **b = 0** | free | free, gated `|b| < a` |
| Center criterion | min Σ **RMSE** | min Σ **\|intercept\|** | **max R²** among valid |
| R² floor | reject `R² < 0.5` | **none** | reject if no window `R² ≥ min_window_r2` |
| Above aerosol | center 2–5 km | center 2–6 km, start unconstrained | **window start ≥ `min_window_start_m`** |

**The fix (implemented in `rayleigh_fit.py:find_optimal_molecular_window`).** A window is now
**eligible** only if it (rec #1) **starts above** `min_window_start_m` (above the BL aerosol),
(R² fix) has `R² ≥ min_window_r2`, a positive slope, and `|b| < a`, and (rec #2) has the fit
slope consistent with the pointwise median ratio (`relative_error ≤ max_window_rel_error`,
rejecting aerosol curvature). Among eligible windows the **highest-R²** one is selected; if **none**
qualifies the night is **flagged non-calibration** (flag −2) instead of emitting a spurious
constant. New `options.json` knobs (production defaults): `min_window_start_m`=2000,
`min_window_r2`=0.5, `max_window_rel_error`=50. The window-search diagnostic plot now outlines the
eligible region (green) and shows the per-center max-R² (the new score) alongside the old,
degenerate Σ|b|.

**Verified (Payerne CL61, March 2026).** The diagnosed night **2026-03-12 is now correctly
rejected** (no molecular window above the aerosol passes); a clean night like **2026-03-28 selects
center 3.42 km** (window ≈ 2.21–4.63 km AGL, R²≈1) where the old Σ|b| rule would have pushed to
high center. Calibrated-night count is stricter but honest:

| gate | calibrated nights / 31 |
|---|---|
| start ≥ 1500 m, R² ≥ 0.5 | 5 |
| **start ≥ 2000 m, R² ≥ 0.5 (default)** | **3** |
| start ≥ 2500 m, R² ≥ 0.5 | 3 |
| start ≥ 2000 m, R² ≥ 0.4 | 3 |

Lowering the R² floor (0.5→0.4) adds **zero** nights — the rejected nights genuinely lack a
molecular signal above the BL, it is **not** a marginal-R² artefact. The **start altitude** is the
real lever (1500 m→5, 2000 m→3): the classic SNR-vs-aerosol-purity trade-off. 3–5 clean Rayleigh
nights/month is consistent with the literature (ceilometer Rayleigh needs an aerosol-free,
high-SNR column — Wiegner 2014) and with Rayleigh's role here as the **slow absolute anchor** (§2,
needs ~2–3 weeks of points), not the fast detector. `min_window_start_m` can be lowered to ~1500 m
where more points are needed and the BL is low.

Before/after (regenerated diagnostics, `figs_paper_validation/rayleigh_diag/`):

![2026-03-12 — now rejected: no eligible (green) region; R²≈1 only below 2 km in the aerosol, R²→0 above](figs_paper_validation/rayleigh_diag/plots/0-20000-0-06610/2026/20260312_0-20000-0-06610_window_search.png)

![2026-03-28 — clean night: optimum (red ✗) at center 3.42 km inside the eligible high-R² region; grey dotted = old degenerate Σ|b|](figs_paper_validation/rayleigh_diag/plots/0-20000-0-06610/2026/20260328_0-20000-0-06610_window_search.png)

Remaining options (not needed given the above, but available): set the lower bound dynamically
from the L1 `layer`/`layer_aerosol` aerosol-top, or widen the window (§3) to trade nights for
lower per-night noise.

## 5. Per-instrument recommendations

| Instrument | single-calibration integration | change-detection capability | notes / other levers |
|---|---|---|---|
| **CL61 cloud** (910) | **1 day (≈8 %)** | detects a ~20 % step in **~3–4 days** | best fast monitor; but carries the +21 % offset → use Rayleigh/CHM for the absolute scale |
| **CL51 cloud** (910) | **1 day (≈9 %)** | ~few days for a ~20 % step | clean, responsive |
| **CHM15k** (1064, Rayleigh) | 1–2 days (≈13 %) | floor-limited (~9–11 %) → only ≳15 % steps fast | variability is **real instrumental** → **HK correction** (internal-T overlap, Hervo 2016) + recalibrate on laser-age steps; averaging won't help below the floor |
| **CL61 Rayleigh** (910) | **~2–3 weeks** for ≈10 % (per-night 39 %) | slow (~2 weeks) | the **absolute anchor**, not the fast detector; **wider fit window** halves the per-night noise |
| **CL31** (910, cloud) | floor-limited (~19–31 %) | only large (>~20 %) changes | needs attention: tighten cloud-target selection + HK/window monitoring |

**Cross-cutting recommendations**
1. **Keep the single-calibration integration SHORT** — daily for the cloud method, per-night
   for Rayleigh — and feed every point to the **Kalman**. Do **not** pre-average over months.
2. **Tune the Kalman as a change tracker, not an averager:** measurement noise = the
   single-calibration σ (8 % cloud / 39 % night Rayleigh), process noise set to follow real
   drift/steps. **Guard it** against divergence on short records (it failed for 3 CL61 units).
3. **Combine the two methods:** cloud for day-to-day *responsiveness* (catches window/laser
   events fast), molecular/CHM for the *absolute* scale (validated anchor, updated more slowly).
4. **Widen the 910 nm Rayleigh fit window** moderately (×2 less per-night noise) so it too can
   give a point on a useful timescale.
5. **Drive/monitor with housekeeping** (laser power, internal temperature, window transmission)
   — this catches the episodic step events directly and removes the non-averageable
   instrumental floor, which no integration can.
6. Keep robust outlier rejection and the strict water-vapour correction (910 nm).

![Cross-instrument scatter (precision)](figs_paper_validation/cross_instrument_stability.png)
![Instrumental driver of the residual floor (CHM15k)](figs_paper_validation/calib_stability_chm15k_drivers.png)

---
*Scripts: `analyze_calib_integration.m` (Allan deviation), `run_cl61_rayleigh_params.py`
(fit-window / quality / night-length test), `analyze_chm15k_stability.m` +
`cross_instrument_stability.py` (drivers & cross-instrument). Companion:
`calibration_stability_report.md` (driver attribution + literature).*
