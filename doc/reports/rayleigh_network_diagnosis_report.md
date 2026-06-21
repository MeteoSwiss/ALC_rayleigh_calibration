# Rayleigh calibration across the whole network — time series & problematic-station diagnosis (2026)

*Generated 2026-06-21. `validation/run_rayleigh_diag_light.py` + `analyze_rayleigh_diag.py`, on the v2
constants from the network run. All 153 CHM15k + Mini-MPL streams, L1, 2026.*

## Metrics

Per stream (optimized v2, L1):

- **clear%** = clear (fit-reaching) nights ÷ archive days — *clear-sky availability* (atmosphere).
- **valid% (of clear)** = valid calibrations ÷ clear nights — *method yield on clear sky*.
- **valid% (of archive)** = valid ÷ archive days = availability × yield — *overall productivity* (the
  sorted bars below; low here can be either too few clear nights OR low yield).
- **σ_SD**, **outlier%** = short-term variability and drift-aware outlier rate of the nightly constant.

Four diagnostic metrics, sampled per station, map to the candidate failure causes:

| cause | diagnostic | "bad" signature (vs healthy CHM15k median) |
|---|---|---|
| not enough clear sky | clear% | **11%** vs 33% |
| lots of FT aerosol | median window scattering ratio | **1.18** vs 1.12 |
| low laser | near-range signal strength → SNR proxy (+ `laser_life_time`) | SNR **13** vs 33 |
| electronic background | far-range noise → SNR proxy | SNR **22** vs 33 (high noise) |

A station is **problematic** if (within its type) it falls in the worst quartile of overall yield, or
the worst decile of σ_SD or outlier rate. Each problematic station is assigned the **dominant** cause
(largest standardized exceedance of the per-type threshold); secondary causes are listed too.

## Result — 59 of 153 streams problematic

cause | # streams | what it looks like
---|---|---
✅ ok | 94 | clear 33%, scat 1.12, SNR 33, σ_SD 8.9% — healthy
🟫 electronic background | 19 | low SNR from high far-range noise (clear sky fine)
🟪 not enough clear sky | 16 | clear% ~11% — Arctic/maritime, persistent cloud
🟧 low laser | 12 | low SNR from weak near-range signal; often old laser
🟥 FT aerosol | 6 | scattering ratio ≥1.18 — persistent free-tropospheric aerosol
⬜ other | 6 | flagged on σ_SD/outliers without a single dominant metric

![CHM15k network diagnosis: (top-left) overall yield sorted, coloured by cause; (top-right) clear-sky availability drives yield; (bottom-left) FT aerosol vs σ_SD; (bottom-right) SNR vs outlier rate.](figs_paper_validation/rayleigh_diag/fig_cause_overview.png)

## What is wrong with the problematic stations

**Not enough clear sky (16).** Defining feature: clear% ≈ 11% (a third of healthy). These are
high-latitude / maritime sites with persistent cloud — **Jan Mayen, Bjørnøya, Hopen, Tórshavn**
(Arctic), **Bonaire** (trade-cumulus), **Camborne**. Nothing is wrong with the instrument; there
simply are too few clear nights to calibrate. Not fixable by tuning — only by accumulating over a
longer period.

**Low laser (12).** Defining feature: low SNR driven by *weak near-range signal*, frequently with a
high `laser_life_time`. Clearest cases: **Vásárosnamény** (laser ≈ 53 900 h — far beyond a typical
~20 k h service life), **Nieuwkoop**, **Gottfrieding**. Aging/under-powered laser → low SNR in the
2–6 km fit band → noisy or failed molecular fits. Action: flag for laser service / re-collimation.

**Electronic background (19).** Defining feature: low SNR driven by *high far-range noise* (signal in
the fit band fine). E.g. **Friesoythe, Klippeneck, Bern**. Elevated detector/electronic background
raises the noise floor, inflating the window scattering ratio and σ_SD. Action: check detector
baseline / dark-current and overlap.

**FT aerosol (6).** Defining feature: persistently elevated window scattering ratio (≥1.18). E.g.
**Payerne (1.24), Bern (1.19)** — sites with frequent free-tropospheric aerosol (Saharan dust,
boundary-layer venting). The fit window rarely finds aerosol-free molecular air, so v2's
scattering-ratio gate rejects nights even when skies are "clear". Action: these benefit most from the
relaxed-scattering v2 (C8) and from higher fit windows; consider per-site scattering gate.

> Note: Payerne is interesting — clear% 43% (plenty of clear nights) but valid ≈ 1%, with scat 1.24
> AND low SNR (20): it is limited by **both** FT aerosol and background, not clear-sky availability.

**Mini-MPL (5):** Aléria (83%) and Lille (55%) healthy; Brest flagged *background*, Trappes
*FT-aerosol*, the mobile unit *low-laser* (high σ_SD despite many nights).

## Time series — full network

Every CHM15k stream's nightly lidar constant (v2, L1), sorted worst-yield first; green = robust
median, red = outliers, pink panels = problematic. 6 pages cover all 148 CHM15k; Mini-MPL separately.

![CHM15k time series page 1/6 (lowest-yield streams)](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p1.png)
![CHM15k time series page 2/6](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p2.png)
![CHM15k time series page 3/6](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p3.png)
![CHM15k time series page 4/6](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p4.png)
![CHM15k time series page 5/6](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p5.png)
![CHM15k time series page 6/6 (highest-yield streams)](figs_paper_validation/rayleigh_diag/fig_ts_CHM15k_p6.png)
![Mini-MPL time series](figs_paper_validation/rayleigh_diag/fig_ts_MiniMPL_p1.png)

## Problematic-station list

The full table (59 streams, with valid%, clear%, σ_SD, outlier%, scattering, SNR, laser age, and the
assigned cause) is in `figs_paper_validation/rayleigh_diag/problematic_stations.md`. Worst by yield:

site | type | valid% | clear% | σ_SD% | scat | SNR | laser(h) | cause
---|---|---|---|---|---|---|---|---
Tórshavn | CHM15k | 0 | 0 | – | – | 91 | 16901 | clear-sky
QUALAIR | CHM15k | 0 | 0 | – | – | 30 | 42629 | clear-sky
Friesoythe | CHM15k | 0 | 31 | – | 1.18 | 61 | 26136 | background
Bonaire | CHM15k | 0 | 4 | – | 1.15 | 76 | 16012 | background→clear-sky
Payerne | CHM15k | 1 | 43 | – | 1.24 | 20 | 8934 | FT-aerosol + background
Vásárosnamény | CHM15k | 1 | 59 | – | 1.19 | 12 | 53948 | low-laser
Jan Mayen | CHM15k | 2 | 2 | – | 1.08 | 25 | 49507 | clear-sky
Nieuwkoop | CHM15k | 2 | 34 | – | 1.13 | 13 | 9390 | low-laser

## Takeaways

- ~60% of the CHM15k network calibrates healthily; the problematic 40% split into four physically
  distinct, separable causes.
- **Two are not instrument faults** — *clear-sky-limited* (Arctic/maritime) and *FT-aerosol* sites are
  atmosphere-limited; they need longer accumulation or per-site gate relaxation, not repair.
- **Two are instrument health flags worth acting on** — *low-laser* (esp. > ~40 k h lasers like
  Vásárosnamény, QUALAIR, Jan Mayen) and *electronic-background* stations should be flagged for
  service. The far-range-noise and laser-age metrics here can drive an automated network health alert.

## Reproduce
```
python validation/run_rayleigh_diag_light.py   # sampled HK + scattering per stream
python validation/analyze_rayleigh_diag.py     # summary + figures + problematic_stations.md
```
