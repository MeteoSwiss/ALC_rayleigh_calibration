# Molecular-window methods — full-archive multi-site comparison

*MeteoSwiss E-PROFILE ALC paper — M. Hervo. 2026-06-18.*

The 7 molecular-window detection methods (`main`, `improved`, `matlab`, `calipso`, `earlinet`,
`optimal`, `bellini`) run over the **entire E-PROFILE L2-monthly archive** for **14 instruments** —
10 CHM15k (Payerne, Lindenberg, Aosta, Palaiseau, Granada, Magurele, Bergen, Oslo, Hamburg,
Hohenpeissenberg) + 4 Mini-MPL (Brest, Toulouse, Corsica, SIRTA) — sampling **5 nights/month** over
each site's full period (~80–113 months; ~400–520 fit-nights per site). This stresses every method on
years of real, varied atmospheres. Driver: `longrun_methods.py` (per-instrument JSON checkpoints).

## Headline

Over hundreds of nights the **classic std/mean CV explodes (100–1500 %) for the higher-yield methods**
(`main`, `calipso`, `improved`, `matlab`, `optimal`) — they each admit a few nights with a spurious
window and an extreme constant. The **robust scatter (MAD-based) collapses to 31–42 % for all methods**
(the inflation was rare outliers, which the Kalman would smooth). The two views together give the real
ranking:

| method | calibrated-frac | robust CV % | std CV % (outlier-sensitive) | med R² | med temporal_cv |
|---|---|---|---|---|---|
| **optimal** | 0.45 | **31** | 128 | **0.97** | 0.15 |
| **earlinet** | 0.45 | **32** | 41 | 0.90 | **0.13** |
| bellini | 0.42 | 33 | 43 | 0.87 | 0.21 |
| matlab | 0.52 | 35 | 380 | 0.94 | 0.22 |
| improved | 0.46 | 38 | 293 | 0.97 | 0.16 |
| main | 0.57 | 37 | 445 | 0.91 | 0.25 |
| calipso | 0.79 | 42 | 385 | 0.65 | 0.27 |

(Full per-instrument numbers: [`ranking_robust_longrun.md`](figs_paper_validation/molecular_methods_longrun/ranking_robust_longrun.md),
[`method_comparison_multisite.md`](figs_paper_validation/molecular_methods_longrun/method_comparison_multisite.md).)

![Usable nights + robust night-to-night CV per method, full archive (14 sites)](figs_paper_validation/molecular_methods_longrun/summary_robust_longrun.png)

![Calibration-constant time series per method, full archive, per site](figs_paper_validation/molecular_methods_longrun/timeseries_longrun.png)

## Precision — the right way to evaluate (drift-insensitive)

**CV is the wrong metric**: it mixes measurement *precision* with the real *seasonal + laser-ageing
drift* a calibration is supposed to track, so it penalises a precise instrument for having a seasonal
cycle. We instead use metrics that remove slow drift (full numbers:
[`precision_longrun.md`](figs_paper_validation/molecular_methods_longrun/precision_longrun.md)):

- **σ_SD** — *successive-difference* precision (von Neumann): robust scatter of |CLᵢ₊₁−CLᵢ| between
  time-ordered consecutive calibrations ÷√2. Slow drift cancels in the difference ⇒ short-term noise only.
- **σ_detrend** — robust scatter of CL minus its ~2-month rolling median (seasonal + trend removed).
- **σ_within-month** — robust scatter of CL pooled within each calendar month (season ≈ const).
- **σ_night** — average single-night spread (in-window std of signal/molecular ÷ CL).
- **valid %** — fraction of sampled nights yielding a valid calibration.

| method | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_month % | CV % (ref) |
|---|---|---|---|---|---|---|
| **optimal** | 45 | 8.6 | **13.6** | 12.5 | 8.4 | 128 |
| **earlinet** | 45 | **4.2** | 15.9 | 13.6 | 10.0 | 41 |
| bellini | 42 | 10.7 | 18.3 | 15.4 | 11.0 | 43 |
| matlab | 52 | 16.7 | 18.5 | 15.3 | 12.2 | 380 |
| improved | 46 | 9.4 | 19.2 | 16.1 | 12.3 | 293 |
| main | 57 | 20.0 | 20.4 | 16.0 | 12.1 | 444 |
| calipso | **79** | 14.9 | 24.9 | 19.9 | 16.4 | 385 |

![Yield and drift-insensitive precision metrics per method](figs_paper_validation/molecular_methods_longrun/precision_longrun.png)

**This resolves the earlier metric ambiguity.** The decisive facts:
1. **CV (128–444 %) ≫ σ_SD (14–25 %) for every method** — i.e. *most of the CV was real drift, not
   noise*. The genuine per-night precision is ~14–25 %, far better than CV implied. (The two methods with
   "low" CV — earlinet 41 %, bellini 43 % — only looked good because they reject the drift-revealing
   nights; their σ_SD is mid-pack.)
2. **`optimal` is the most precise method** by all three drift-insensitive measures (σ_SD 13.6 %,
   σ_detrend 12.5 %, σ_within-month 8.4 %) and has a low within-night spread (8.6 %).
3. **`earlinet` is second on σ_SD and has the tightest within-night spread (4.2 %)** — its narrow,
   low-altitude windows are internally very consistent.
4. **`calipso` is the least precise (σ_SD 24.9 %)** despite the highest yield (79 %) — quantity at the
   cost of precision (its noisy high windows).
5. **The drift CV reveals is physical and useful** — the seasonal/laser-ageing cycle (visible in the
   time series) is exactly what the Kalman should track; a precision metric must not conflate it with noise.

## Findings

1. **`optimal` and `earlinet` give the cleanest, most stable long-term calibration** — the lowest
   robust CV (31–32 %), the lowest temporal CV (0.13–0.15), and (for `optimal`) the highest R² (0.97).
   They reject the outlier nights that inflate the other methods' std-CV. **This is the recommendation
   for a production calibration**, where per-night quality and stability matter more than raw count
   (the Kalman handles sparsity).
2. **`calipso` maximises yield (0.79) but at the cost of quality** — lowest R² (0.65, noisy high
   windows) and the highest robust scatter (42 %). Quantity over quality; not advised for ceilometers
   (consistent with the EarthCARE/ATLID point: "highest clean layer" ≠ "best" without a stratosphere).
3. **`main` is the highest-yield/highest-scatter legacy method** — it calibrates the most raw nights
   (0.57) but with the worst std-CV (445 %): the degenerate behaviour, confirmed at archive scale.
4. **The "best by automated score" flips with the metric and period** (improved on a short clean test,
   earlinet by std-CV, calipso by yield-weighted robust score) — i.e. *there is no single winner*; the
   choice is a yield-vs-quality trade-off. For E-PROFILE production we weight **quality + stability →
   `optimal` (or `earlinet`)**.
5. **CHM15k seasonal cycle is visible** in the time series (the known internal-temperature / laser
   cycle) for the stable methods — exactly what a calibration monitor should track; the unstable
   methods bury it under scatter.
6. **Mini-MPL (532 nm) is the easy case** — all methods reach R²≈0.99 and low temporal CV; the clean,
   high-SNR molecular column means even `main`/`calipso` behave. `optimal` still has the lowest CV.

## Recommendation (confirmed at archive scale)

Use **`optimal`** (or **`earlinet`** where a single-profile method is preferred) for CHM15k/CL61/
Mini-MPL Rayleigh calibration. `improved` remains a sound, simpler default. Avoid `calipso`/`main` for
quantitative calibration. `bellini` (ALICENET) is a solid, citable alternative, strongest on high-SNR
CHM15k. Feed the per-night constants to the Kalman, which absorbs the residual outliers.
