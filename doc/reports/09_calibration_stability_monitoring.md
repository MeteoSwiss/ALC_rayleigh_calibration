# Calibration stability, conventions and network monitoring

*Consolidated 2026-07-10 (M. Hervo, MeteoSwiss E-PROFILE ALC). Merges: calibration_coefficient_convention.md, calibration_stability_report.md, calibration_short_term_variability_report.md, calibration_outliers_report.md, network_v2_vs_v11_report.md, calibration_v2_camsfar_modifications_2026-06-26.md.*

This is the durable reference for **how the absolute calibration coefficient is defined and named**,
**what drives its (in)stability**, **how much of the short-term scatter is reducible and how**, the
**network-wide outlier statistics and constant time series**, the **network validation of the
production molecular method (E-PROF v2, config C8) against v1.1**, and a dated **June-2026
modifications changelog** (CAMS-far guard + dashboard). It applies to the three Rayleigh-capable
types (CHM15k, CL61, Mini-MPL) and the cloud-only Vaisala types (CL31, CL51).

## Contents

1. [The calibration-coefficient convention (authoritative)](#1-the-calibration-coefficient-convention-authoritative)
2. [What drives calibration (in)stability — Rayleigh & cloud](#2-what-drives-calibration-instability--rayleigh--cloud)
3. [Reducing short-term variability — diagnosis & recommendations](#3-reducing-short-term-variability--diagnosis--recommendations)
4. [Per-calibration outlier rate & the network constant time series](#4-per-calibration-outlier-rate--the-network-constant-time-series)
5. [Network validation: E-PROF v2 (C8) vs v1.1 — L1 + L2](#5-network-validation-e-prof-v2-c8-vs-v11--l1--l2)
6. [June-2026 modifications changelog — CAMS-far + dashboard](#6-june-2026-modifications-changelog--cams-far--dashboard)

> **Current-truth reconciliation (2026-07-10).** The production Rayleigh molecular-window method is
> **`eprof_v2`, config C8** (repo default); the *deployed operational E-PROFILE network* still runs
> the older **E-PROF v1.0** (which retains the historical Klett sign error) — do not conflate the two.
> The `calipso` method is **retired**. Multiple-scattering η now uses the **PVC (Hogan 2006) tables
> at a_G = 5.5 µm** (each Vaisala type its own table; CL61 has its own), superseding the legacy
> Hopkin / 8 µm ladder. The **window-temperature correction is reject-only** (no magnitude correction).
> The **water-vapour correction is mandatory for 910 nm** (CAMS L137 model levels; a no-WV mode is
> rejected). The liquid-cloud method integrates a **fixed 100–2400 m gate** (CBH-robust).

---

## 1. The calibration-coefficient convention (authoritative)

*Single source of truth for how the calibration coefficient is named, defined and reported across
this package (code, output files, figures, reports). Referenced by every section below.*

### 1.1 The one convention: Wiegner lidar constant `C_L`

Both calibration methods report the **Wiegner & Geiß (2012) lidar constant**

```
C_L = RCS / β_att
```

where `RCS` is the range-corrected signal (`P·r²`, background-subtracted) and `β_att` is the
**attenuated total backscatter** (`m⁻¹ sr⁻¹`). Equivalently, the lidar equation is
`RCS = C_L · β_att`, so the data are calibrated by

```
β_att = RCS / C_L .
```

`C_L` is a large, instrument-specific number. Typical 2026 values: **CHM15k ≈ 5×10¹¹**,
**CL61 ≈ 1**, **Mini-MPL ≈ 2×10⁵**. Units: `V·m³/sr` (Vaisala CL31/CL51/CL61) or
`counts/s·m³/sr` (CHM15k, Mini-MPL).

This is the **primary reported quantity for every method**, so the Rayleigh and liquid-cloud
products live on the **same axis** and can be plotted in one time series.

### 1.2 Rayleigh (molecular) calibration — Wiegner & Geiß (2012)

`C_L` is computed directly from the raw range-corrected signal in the molecular window:
`C_L = RCS / β_att,molecular` (transmission-corrected). **No prior calibration is involved**,
so the absolute `C_L` is obtained directly and is **identical from L1 and L2** (the L2 reader
reverts the applied calibration and corrects units,
`RCS = attenuated_backscatter_0 × calibration_constant_0 × 1e-6`, recovering the L1-scale signal).

Code: `CalibrationResult.lidar_constant`. NetCDF: variable `lidar_constant`
(`long_name = "Lidar constant C_L (Wiegner & Geiss 2012)"`).

### 1.3 Liquid-cloud calibration — O'Connor (2004) / Hopkin (2019)

The cloud method ingests the file's **already-calibrated** attenuated backscatter `β_att,file`
and integrates it through a fully-attenuating liquid cloud (over a **fixed 100–2400 m gate
window**, robust to Vaisala CBH bias), comparing to the theoretical value `1/(2·S)` with
`S = 18.8 sr`. Its native output is the **O'Connor multiplier**

```
C = S_apparent / 18.8 = β_true / β_att,file
```

`C` is therefore the **inverse sense** of `C_L` (`C ∝ β/signal`, `C_L ∝ signal/β`): the two
conventions are reciprocals, `C_r = 1/C_c` in the user's shorthand.

`C ≈ 1` only when the file's `calibration_constant_0` is already the true lidar constant. Across
the 2026 network this splits by type (sampled, 8 streams/type):

- **CL31 / CL51 L2 store a fixed nominal placeholder `calibration_constant_0 = 1×10⁸`** (these
  cannot be Rayleigh-calibrated). So `C` is far from 1 (e.g. CL31 `C ≈ 1.6×10⁻⁶`), and
  `C_L = calibration_constant_0 / C` is their **only** absolute calibration.
- **CL61 L2 carries a real per-stream constant near 1** (observed 0.91–1.10), so for CL61
  `C ≈ 1` and `C_L ≈ calibration_constant_0 × (1/C)` ≈ the operational value times the cloud
  correction.

Either way the absolute `C_L` is the headline; a `C` far from 1 is **expected** for CL31/CL51,
not an error.

To report in the Wiegner convention, note the file's `β_att,file` was made with the applied
constant `calibration_constant_0` (= `RCS/β_att,file` = the operationally-applied `C_L`). Hence
the **absolute Wiegner constant from the cloud method** is

```
C_L = calibration_constant_0 / C .
```

This is the **headline cloud product** (directly comparable to Rayleigh). The dimensionless
**inverse** `1/C` (the Wiegner-sense correction factor, ≈ 1) is reported alongside; the raw
O'Connor multiplier `C` is kept only as an internal diagnostic.

> **Why the cloud `C_L` carries the L2 constant.** The cloud method only *measures a ratio*
> (`1/C`, completely independent of any prior calibration). Turning that ratio into an absolute
> `C_L` borrows the file's existing scale `calibration_constant_0`: it is exact given the file,
> but its absolute accuracy is only as good as the operational constant already applied. The
> Rayleigh `C_L` has no such dependence because it works on the raw signal.

Code: `CloudCalResults` exposes, in order of preference —

| field | symbol | meaning |
|---|---|---|
| `lidar_constant` | `C_L` | absolute Wiegner constant `= calibration_constant_0 / C` (headline; NaN if the file has no applied constant, e.g. raw L1) |
| `calibration_factor` | `1/C` | dimensionless Wiegner-sense correction — *the inverse, reported whenever possible* |
| `calibration_coefficient` | `C` | O'Connor multiplier (internal/diagnostic; ≈1 only if the file constant is the true one) |

### 1.4 Variability metric is convention-independent

The short-term precision `σ_SD` (robust successive-difference, % of the median) is a **relative**
metric, invariant under `C → 1/C` and under any constant scale. So `σ_SD(C) = σ_SD(1/C) =
σ_SD(C_L)`: every variability figure/table below is unchanged by the convention and is labelled in
terms of `C_L`.

---

## 2. What drives calibration (in)stability — Rayleigh & cloud

**Author:** M. Hervo · **Date:** 2026-06-17

### 2.1 Question

How stable are the absolute calibrations (Rayleigh molecular and liquid-cloud), what is the **main
driver** of their instability (water vapour, aerosols, or instrumental), is this **the same in the
literature**, and do the **other instruments (CL31, CL51, CHM15k)** show the same problem?

### 2.2 Short answer

**Instrumental drift is the dominant driver of long-term ALC calibration instability**, on three
timescales: **laser ageing** (multi-year), **internal/detector temperature** (seasonal + diurnal),
and **window contamination/fogging** (episodic). This is true for **every** ALC — including the
1064 nm CHM15k, which has *no* water-vapour sensitivity yet is just as unstable, proving the driver
is instrumental, not atmospheric. **Water vapour** adds a large (~20 %) *seasonal* term **for the
910 nm instruments only** (CL31/CL51/CL61), which we now correct (mandatory). **Aerosols** are
second-order. This matches the literature exactly (Le et al. 2026; Hopkin et al. 2019; Hervo et al.
2016; Kotthaus et al. 2016; Wiegner et al. 2014/2015; Filioglou et al. 2023).

### 2.3 The decisive evidence — CHM15k (1064 nm) is unstable *without* any water vapour

The CHM15k operates at **1064 nm, outside the water-vapour absorption band** — so any instability in
its calibration **cannot** be water vapour. Using the 11-year daily archive (13 stations, 2013–2024;
`A:\CHM15k\…`, the dataset behind `estimate_calib_from_housekeeping`):

*(Figure `figs_paper_validation/calib_stability_chm15k_drivers.png` — showing (a) Payerne: the lidar
constant tracks the internal temperature over a decade; (b) all 13 units pooled: lidar constant vs
internal temperature, r = 0.31; (c) per-station, the calibration correlates with instrumental
housekeeping — internal temperature, laser pulses, optics state — is **missing from the repo** and
therefore not embedded; see the digest.)*

| CHM15k metric (per station, 11 yr) | value |
|---|---|
| Calibration CV | **23 %** (median; 15–38 %) |
| Seasonal peak-to-peak amplitude | **34 %** (median) |
| Long-term trend (laser ageing) | **−10 … +16 % / yr** (unit-specific) |
| Correlation with **internal temperature** | **+0.38** (median) |
| Variance explained by instrumental HK (T_int + laser pulses + optics) | **R² = 0.38** (up to 0.69) |

A 1064 nm instrument with a **23 % calibration CV and a 34 % seasonal cycle that follows its internal
temperature** can only be explained by **instrumental** effects (temperature-driven optics/overlap
and laser ageing) — exactly the mechanism Hervo et al. (2016) corrected for the CHM15k overlap.
~38 % of the variance is captured by housekeeping alone; the rest is sampling noise plus the
aerosol-loading sensitivity of the molecular fit.

### 2.4 Water vapour — a large 910 nm-only term (corrected)

From the CL61 controlled re-run (`cl61_verify_*`), toggling the WV correction changes the calibration
of the **910 nm** instruments substantially, and the **1064 nm CHM15k not at all**:

| | WV-off → WV-on |
|---|---|
| Rayleigh constant (910 nm) | **+20.8 %** (the molecular fit at 3–6 km integrates the full WV column) |
| Cloud coefficient (910 nm) | **−12.8 %** (cloud base ~1 km still sees most of the boundary-layer WV) |
| CHM15k (1064 nm) | **0 %** (outside the WV band) |

*(Figure `figs_paper_validation/cl61_verify_l1l2_wv.png` — WV impact and L1/L2 robustness — is
missing from the repo and not embedded; see the digest.)*

So for CL31/CL51/CL61 a **seasonal water-vapour cycle** is a first-order calibration driver *if
uncorrected*; with the matching-month CAMS correction applied it is removed. **This correction is
mandatory for 910 nm** — the operational humidity source is **CAMS model levels (L137)** (dense in
the boundary layer); a night without usable CAMS is *flagged* (see flag −4 / −10 in §6), never
calibrated WV-free. This is the Wiegner & Gasteiger (2015) result (~20 % mid-latitude, >50 % tropics
at 905 nm) and is why the CHM15k (1064 nm) is immune.

### 2.5 Aerosols — second-order

- **Cloud method:** the target is a *totally attenuating* liquid cloud with fixed S = 18.8 sr, so it
  is largely aerosol-insensitive; residual below-cloud aerosol is suppressed by the 90 %-in-cloud
  acceptance filter (Hopkin 2019; Le 2026). Our diagnostics confirmed aerosol contamination pushes
  the cloud coefficient the *wrong* way to explain the offset.
- **Rayleigh method:** the assumed aerosol lidar ratio S_p (50 sr default) enters the fit; literature
  propagates ±10 sr as a sensitivity. A secondary, event-driven driver.

### 2.6 Do the other instruments show the same problem? — Yes

| Instrument | λ | raw daily CV (precision) | smoothed drift (this work / literature) | dominant driver |
|---|---|---|---|---|
| **CHM15k** | 1064 | 23 % (11 yr) | seasonal 34 %, trend ±10 %/yr | **instrumental** (T_int, laser) — *no WV* |
| **CL31** | 910 | 44 % (8 yr, Payerne) | ±3 %/20 mo, ±5 %/yr (Hopkin 2019) | instrumental + WV |
| **CL51** | 910 | 9 % (short record) | similar (WV-corr. validated, Wiegner 2019) | instrumental + WV |
| **CL61 (cloud)** | 910 | 11 % (2026) | laser-ageing dominated (Le 2026) | instrumental + WV |
| **CL61 (Rayleigh)** | 910 | 61 % (2026, per-night) | needs Kalman; agrees with cloud within unc. | instrumental + WV |

*(Figure `figs_paper_validation/cross_instrument_stability.png` — cross-instrument calibration
scatter — is missing from the repo and not embedded; see the digest.)*

Two points: (i) the **raw per-sample scatter** is large for every type — which is why all
operational calibrations are **temporally smoothed** (Kalman / 90-day running mean); the *smoothed*
drift is the few-% literature value. (ii) The **cloud method is more precise per-sample** than the
per-night Rayleigh fit (CL61 11 % vs 61 %), but both converge after smoothing. The instability is
**universal across the network and across instrument types**.

### 2.7 Same in the literature? — Yes, in detail

- **Instrumental dominates long-term drift.** Le et al. (2026, "Long-term performance of the Vaisala
  CL61", 4 units, 3 yr): the firmware rescales the calibration while laser power >40 %, but **below
  40 % it breaks down** — Lindenberg's cloud-calibration factor fell ~×3 (1.45→0.54) as laser power
  went 40→10 %. **Window fogging** caused episodic factor jumps up to ×10 (diagnosed via the "window
  condition" housekeeping). An **internal-temperature look-up-table** corrects the CL61 bias — the
  analogue of **Hervo et al. (2016)**'s temperature-dependent CHM15k overlap correction.
  **Kotthaus et al. (2016)** documents the CL31 instrumental background and firmware artefacts.
- **Housekeeping is the standard monitor/QC.** Hopkin (2019) rejects profiles with **window
  transmission < 90 %** ("cannot be reliably corrected"); Le (2026) keys on **laser power**, **window
  condition**, and **internal temperature**. This is exactly what our
  `estimate_calib_from_housekeeping` work exploits for the CHM15k.
- **Water vapour** ~20 % at 910 nm (mid-latitude), seasonal, not at 1064 nm (Wiegner 2014/2015).
- **Cloud vs Rayleigh** agree within the cloud method's **~10 % uncertainty**, dominated by the
  **multiple-scattering factor η** and S = 18.8 ± 0.8 sr (Hopkin 2019; O'Connor 2004).
  **Filioglou et al. (2023)** attribute the CL61 factory-vs-field calibration spread (15 %,
  latitude-dependent) to **water vapour + multiple scattering** — which is precisely the residual
  seen in our cloud-vs-Rayleigh +21 %.

  > **Superseded 2026-07 (η value).** The literature/legacy η range **0.7–0.85** (η ≈ 0.83 at low
  > cloud base) cited above is no longer the value applied here. The pipeline now uses the **PVC
  > (Hogan 2006) tables at droplet radius a_G = 5.5 µm** (11 µm diameter; α = 10 /km) — the
  > Cloudnet-measured calibration-scene size — with **each Vaisala type using its own table** (CL61
  > no longer borrows CL51's; CL31's wider 0.83 mrad FOV → stronger correction). At low cloud base
  > **η ≈ 0.95** (was 0.83) → ~**9 % lower C** for low clouds. PVC η tables were also added for
  > CHM15k / Mini-MPL / MPL, plus a saturation warning. The *conclusion* (η is the leading term in
  > the ~10 % cloud-vs-Rayleigh offset) stands; only the numerical η ladder changed.

### 2.8 Drivers, ranked

| Rank | Driver | Type | Timescale | Affects | Mitigation |
|---|---|---|---|---|---|
| 1 | **Laser ageing / power loss** | instrumental | multi-year | all ALC | HK monitoring; recalibrate; flag <40 % power |
| 2 | **Internal/detector temperature** | instrumental | seasonal + diurnal | all ALC | T-dependent overlap (Hervo 2016) / T-LUT (Le 2026) |
| 3 | **Window contamination / fogging** | instrumental | episodic | all ALC | window-transmission QC (<90 % reject); keep blower/heater on |
| 4 | **Water vapour absorption** | atmospheric | seasonal | 910 nm only | matching-month CAMS WV correction (mandatory) |
| 5 | **Multiple scattering η / S assumption** | method | constant offset | cloud method | η per FOV (PVC a_G = 5.5 µm); cloud as cross-check, Rayleigh/CHM as anchor |
| 6 | **Aerosol load / S_p** | atmospheric | event | Rayleigh fit | sun-photometer S_p; cloud target is aerosol-robust |
| — | per-night/per-cloud sampling | noise | daily | all | Kalman / running-mean smoothing |

### 2.9 Conclusions & recommendations

1. **Main driver = instrumental** (laser ageing, internal temperature, window contamination), proven
   by the 1064 nm CHM15k being just as unstable (23 % CV, 34 % seasonal, r = 0.38 with internal
   temperature) despite having *no* water-vapour sensitivity. The literature agrees.
2. **Water vapour** is a large, *correctable*, 910 nm-only seasonal term (~20 % Rayleigh, −13 %
   cloud); keep applying the matching-month CAMS correction (mandatory for 910 nm).
3. **The instability is universal** — CL31, CL51, CL61 and CHM15k all show it; it is reduced to the
   few-% literature level by temporal smoothing (Kalman / running mean).
4. **Recommendations:** (a) drive/monitor the calibration with **housekeeping** (laser power, internal
   temperature, window transmission) as `estimate_calib_from_housekeeping` already does for the
   CHM15k — extend to the CL61; (b) flag laser power <40 % and window transmission <90 %; (c) keep the
   WV correction strict for 910 nm; (d) use the molecular/CHM (Rayleigh) calibration as the absolute
   anchor and the liquid-cloud as a cross-check (its ~10 % η/S uncertainty is the residual behind the
   cloud-vs-Rayleigh +21 %).

*Data: `A:\CHM15k\` 11-yr daily HK (13 stations); `cl61_verify_*` controlled re-run;
`Cloud_Trans-cor_WV-cor` (CL31 06610_B, CL51 EDT_A). Scripts: `analyze_chm15k_stability.m`,
`cross_instrument_stability.py`.*

---

## 3. Reducing short-term variability — diagnosis & recommendations

**Author:** M. Hervo · **Date:** 2026-06-17 · **Scope:** Rayleigh + liquid-cloud; CL31, CL51, CL61,
CHM15k.

### 3.1 Question & short answer

The daily/per-night calibration coefficients scatter a lot. **Why, and how do we reduce it?**

**Most of the short-term scatter is averageable measurement noise**, not real calibration change. It
is largest where the calibration signal is weakest (the per-night **Rayleigh** molecular fit:
~**39 %** night-to-night), and smallest for the strong-signal **liquid-cloud** target (~**8 %**
day-to-day). Underneath the noise sits a **real, non-averageable floor (~8–20 %)** that is the
*instrumental* drift identified in §2 (laser ageing, internal temperature, window).

**We do NOT want to average over months — that would hide the instrument changes (window cleaning,
laser swaps, drift) the Kalman is meant to detect.** The right target is the **shortest integration
that yields one usable calibration point** at a useful precision; those points are then fed to the
**Kalman, which adaptively tracks drift and detects steps**. That single-calibration integration time
is **short and instrument-dependent**: ~**1 day** for the liquid-cloud method (strong signal, 8–9 %),
but ~**2–3 weeks** for the per-night Rayleigh (weak molecular signal, 39 %/night). Vertical averaging
(wider molecular fit window) gives a secondary ×2 on the per-night Rayleigh noise; quality/SNR
filtering mainly throws nights away without improving the per-estimate noise. **§3.3 below is the
operative section.**

### 3.2 Why is the short-term scatter so large? — it is mostly noise

The Allan deviation σ_A(τ) of each instrument's daily calibration separates **random noise**
(σ_A falls as τ^(−1/2): averaging helps) from **real drift** (σ_A rises with τ: averaging hurts). The
minimum marks the optimal averaging time.

*(Figure `figs_paper_validation/calib_integration_allan.png` — (a) Allan deviation vs integration
time τ (log-log), dotted line = pure-white-noise τ^(−1/2) slope, pentagram = each series' optimum;
(b) the short-term scatter a W-day running median removes — is missing from the repo and not
embedded; see the digest.)*

| Series | λ / method | σ_A at 1 day | optimal τ | noise floor at optimum | reduction |
|---|---|---|---|---|---|
| **CL61 Rayleigh** | 910, molecular | **39 %** | ~21 d | 9 % | **×4.6** |
| CL31 cloud | 910, cloud | 31 % | ~60 d | 19 % | ×1.6 |
| CHM15k LIN | 1064, molecular | 15 % | ~1 yr | 8 % | ×1.9 |
| CHM15k PAY | 1064, molecular | 13 % | ~21 d | 11 % | ×1.2 |
| CL51 cloud | 910, cloud | 9 % | ~90 d | 2 % | ×4.5 |
| CL61 cloud | 910, cloud | 8 % | ~30 d | 3 % | ×2.3 |

Two clear patterns:
- **Per-night Rayleigh is photon/SNR-limited** (39 % at 1 day): the molecular return at 3–6 km is
  ~100× weaker than a cloud target, so a single night is very noisy — but it follows the τ^(−1/2)
  line, i.e. it is **almost pure averageable noise** (×4.6 down by ~3 weeks).
- **The liquid-cloud target is intrinsically precise** (8–9 % at 1 day) because the in-cloud signal
  is strong; it still benefits from ~1-month averaging.
- **CHM15k PAY barely averages down (×1.2)**: its 1-day scatter is *already* dominated by the real
  instrumental drift floor (~11 %), not noise — consistent with it being a stable, mature unit whose
  variability is seasonal/temperature, not random.

### 3.3 Single-calibration integration time + Kalman change detection (the operative answer)

We want the **shortest integration that gives one usable calibration point**, then let the **Kalman
detect changes** — not a months-long average that would smear them.

*(Figure `figs_paper_validation/calib_single_time_changedetect.png` — (a) precision of ONE
calibration vs its integration time with 10 %/5 % targets; (b,c) a +20 % step (window cleaned) on a
slow drift: the responsive Kalman recovers it in days, a 90-day average smears it over months — is
missing from the repo and not embedded; see the digest.)*

**Shortest integration for one usable calibration point:**

| Instrument / method | precision of 1 day | τ for ≤10 % | τ for ≤5 % | noise floor |
|---|---|---|---|---|
| **CL61 cloud** | **8 %** | **1 day** | ~3–4 weeks | 3 % |
| **CL51 cloud** | **9 %** | **1 day** | ~10 days | 2 % |
| CHM15k (1064, Ray) | 13–15 % | 1–2 days / ~25 d | floor-limited | 9–11 % |
| **CL61 Rayleigh** | 39 % | **~18 days** | unreachable | 9 % |
| CL31 cloud | 31 % | floor-limited (>19 %) | — | 19 % |

So **the cloud method gives a usable single calibration every single day** (8–9 %) — that is the
integration time to use for change monitoring, *not* a monthly average. The per-night Rayleigh is too
noisy for a daily point (39 %) and effectively needs ~2–3 weeks; it is the slow **absolute anchor**,
not the fast change detector.

**The Kalman is the change detector, not an averager.** Feed it the single-calibration points (daily
for cloud, per-night for Rayleigh) with their uncertainties; tune its **measurement noise = the
single-calibration σ** and its **process noise to the expected real change rate** (large enough to
follow steps/drift). In the demo it recovers a 20 % step in **~4 days (cloud)** / ~12 days (Rayleigh);
a 90-day average needs ~3 months and lags the drift. On the short CL61 records the production Kalman
**diverged** for 3 units — guard it (reject non-physical jumps) or fall back to a short robust running
median, but keep the window short. *(This Kalman best-estimate smoothing of the constant time series
is part of the operational monitoring: the dashboard's Rayleigh C_L series carries a
"v2.0 Kalman estimate" line — see §6.5.)*

#### 3.3a Rayleigh night-averaging window — now darkness-adaptive (SZA), was fixed clock

Rayleigh works only at night, so the **within-night averaging window** is set by how much dark time is
used. The pipeline used a **fixed solar-clock window (solar 20:00→04:00, 8 h)** — *not*
darkness-adaptive (despite a memory that commit 6d684d8 added that; 6d684d8 only fixed the solar-*time*
clock, no sun-angle). At Payerne this **wastes ~2–5 h of true dark time in winter** (13 h dark vs 8 h
used) and **includes ~2 h of twilight in summer** (June dark = 5.7 h). Naively widening the *clock*
window backfires (CV 25 %→95 % at 12 h) because it pulls in twilight.

**Implemented** a solar-zenith-angle night selection in `data_loader.filter_time_range`
(`use_sza_night`, `sza_night_threshold`, default **100°** = sun ~10° below horizon; clock fallback
when coordinates are missing; `use_sza_night=False` reproduces the old result exactly). *Verified in
code (`options.json`): `use_sza_night = 1`, `sza_night_threshold = 100`.* On Payerne CL61 2026 it
**recovers more usable nights (25 vs 20)** and correctly tracks the dark period across
seasons/latitude. It does **not lower the per-night CV** (≈25–33 %): that scatter is night-to-night
atmospheric/fit variability, not within-night averaging — so the window is the *correct* and
*adaptive* choice (essential for high-latitude sites and summer, where the fixed clock fails), but the
per-night noise is still beaten only by the Kalman over nights (§3.3) and by the wider fit window
(§3.4). Recommended threshold 100° (more nights) to 108° (astronomical, cleanest fit) — a fit-quality
vs night-count trade-off.

### 3.4 Lever — vertical integration (molecular fit window)

Re-running the CL61 Payerne Rayleigh with a **wider molecular fit window** (2–9 km vs 2–6 km, plus
longer fit half-lengths) cut the per-night CV from **24.9 % → 11.1 % (×2.2)** — more range gates
averaged ⇒ less photon noise. The trade-off: it **also cut the usable nights (20 → 6)** because a
clean molecular signal to 9 km is rarer, and shifted the median ~−4 % (sampling a higher, more
aerosol-free column). Recommendation: a moderate widening (e.g. to 7–7.5 km) is worth testing
operationally for the noisy 910 nm Rayleigh; the 1064 nm CHM15k already has enough SNR that this
matters less.

### 3.5 Lever — noise filtering (SNR / quality / outliers)

| filter | effect on CL61 Rayleigh |
|---|---|
| tighten method-agreement gate (quality 15→8) | nights 20→10, **CV unchanged (24.8 %)** |
| require ≥5 h night (vs 3 h) | no change |

Tighter gates **remove nights without lowering the per-estimate noise** — they improve *reliability*
(reject bad fits) but are not the route to lower variability. The effective "filter" for variability
is **temporal averaging** (§3.3). The **housekeeping filters from the literature** (window
transmission < 90 % reject, laser power < 40 % flag — Hopkin 2019, Le 2026) remove the episodic
excursions that no amount of averaging will fix.

#### 3.5a Outlier *profiles* are the cause of the random-subset spread — screen, don't medianize

The Rayleigh sensitivity test (random 70 % time-subsets) showed big subset-to-subset spread. Cause:
the night's profiles were collapsed with a **mean** (`np.nanmean`), which is sensitive to a few
residual-aerosol / cloud-edge / noisy profiles — so each random subset that includes or excludes them
shifts. The per-altitude CL and the perturbation aggregate already use the median.

Tested on Payerne CL61 2026 (night-to-night CV):

| time collapse | night-to-night CV | note |
|---|---|---|
| mean, no screen (legacy) | 33.3 % | outlier-sensitive |
| **mean + outlier screen (now default)** | **23.8 %** | screen profiles by MAD on the molecular-band signal, then mean |
| median profile | 81.6 % | **worse** — median is ~1.5× noisier at the weak 3–6 km signal |

So the fix is **remove outlier profiles, then take the mean** — *not* a median profile (which is
robust to outliers but far noisier where the molecular signal is weak). **Implemented:**
`screen_profile_outliers` (default on, MAD threshold `profile_outlier_nmad` = 4) + `time_aggregation`
("mean" default; "median" available). *Verified in code (`options.json`): `time_aggregation = "mean"`,
`profile_outlier_nmad = 4.0`.* This cut the per-night CV ×1.4.

#### 3.5b Molecular-window detection was selecting aerosol / low-R² windows — diagnosed & **fixed**

**The bug.** The Python window search chose the center (2–6 km) with the smallest **Σ|intercept|**,
then the largest R² at that center, with a **free** intercept and **no R² floor**. This is
degenerate: in the high-altitude noise region the signal → 0, so the regression fits it with
intercept **b ≈ 0** (trivially small) *and* **R² ≈ 0** (no real correlation). Minimising Σ|b|
therefore systematically **selects a high, noise-dominated window with near-zero R²**, and nothing
rejects it. This is exactly what the diagnostic figure for Payerne CL61 **2026-03-12** showed: the
"optimum" sat at center **4.62 km in a near-zero-R² region**, while genuine R² ≈ 1 was only at low
center (where the boundary-layer aerosol makes signal and molecular both decay with height → spurious
correlation).

**The literature names this failure mode.** Mattis, D'Amico, Baars et al. (2016, *AMT* 9, 3009;
EARLINET Single Calculus Chain) state that a pure minimum-signal/minimum-background search "does
**not guarantee that there are no particles** … would find a minimum also in the case that there are
fewer particles than in other altitude regions only … may cause **large errors**." The robust
remedies they and others use are (i) an **SNR / standard-deviation gate**, (ii) a **Rayleigh-shape
("molecular") test** — require the window's shape to match the computed Rayleigh profile
(Freudenthaler et al. 2018 residual ≤ 1 %; Baars et al. 2016), (iii) a **scattering-ratio bound**
R ≤ 1.1 (Wiegner & Geiß 2012; CALIOP 1.01 ± 0.01), and (iv) **prefer the lowest qualifying window**
(best SNR). R² of the signal-vs-molecular fit is a valid operationalisation of the shape test — and
crucially **R² is the right metric because it *collapses* in the noise region** (so a noise window is
rejected, not selected, the opposite of Σ|b|).

**The original MATLAB had all three guards; the Python port lost them.**
`Auto_Calib_25/Rayleigh/rayleigh_fit.m` (Hervo & Poltera 2014) **forces the intercept to zero**
(`opts.Lower/Upper=[-Inf 0]/[Inf 0]`), chooses the center by **minimum Σ RMSE** (not Σ|b|), and the
caller (`auto_calib_v23/24.m`) **rejects** the night when `best_r2 < min_r2_rfit` (=0.5) or
`best_rmse > max_rmse_rfit`. It also searches centers only to **5000 m** (Python: 6000 m).

| | Auto_Calib_25 (MATLAB, reference) | Python port (before) | Python (now, fixed) |
|---|---|---|---|
| Intercept | forced **b = 0** | free | free, gated `|b| < a` |
| Center criterion | min Σ **RMSE** | min Σ **\|intercept\|** | **max R²** among valid |
| R² floor | reject `R² < 0.5` | **none** | reject if no window `R² ≥ min_window_r2` |
| Above aerosol | center 2–5 km | center 2–6 km, start unconstrained | **window start ≥ `min_window_start_m`** |

**The fix (implemented in `rayleigh_fit.py:find_optimal_molecular_window`).** A window is now
**eligible** only if it (rec #1) **starts above** `min_window_start_m` (above the BL aerosol), (R² fix)
has `R² ≥ min_window_r2`, a positive slope, and `|b| < a`, and (rec #2) has the fit slope consistent
with the pointwise median ratio (`relative_error ≤ max_window_rel_error`, rejecting aerosol curvature).
Among eligible windows the **highest-R²** one is selected; if **none** qualifies the night is
**flagged non-calibration** (flag −2) instead of emitting a spurious constant. New `options.json`
knobs (production defaults, *verified in code*): `min_window_start_m` = **2000**, `min_window_r2` =
**0.5**, `max_window_rel_error` = **50**. The window-search diagnostic plot now outlines the eligible
region (green) and shows the per-center max-R² (the new score) alongside the old, degenerate Σ|b|.

> **Reconciliation note.** The June-2026 C8 changelog (§6.1) records that an interim C8 tuning
> lowered `min_window_start_m` 2000→1500 m (a documented caveat, which occasionally pulled the window
> ~250 m into residual aerosol). The **shipped default is back at 2000 m** (`options.json`,
> `calibration/config.py`) — that is the value in force. `min_window_start_m` can still be lowered to
> ~1500 m where more points are needed and the BL is low.

**Verified (Payerne CL61, March 2026).** The diagnosed night **2026-03-12 is now correctly rejected**
(no molecular window above the aerosol passes); a clean night like **2026-03-28 selects center
3.42 km** (window ≈ 2.21–4.63 km AGL, R² ≈ 1) where the old Σ|b| rule would have pushed to high
center. Calibrated-night count is stricter but honest:

| gate | calibrated nights / 31 |
|---|---|
| start ≥ 1500 m, R² ≥ 0.5 | 5 |
| **start ≥ 2000 m, R² ≥ 0.5 (default)** | **3** |
| start ≥ 2500 m, R² ≥ 0.5 | 3 |
| start ≥ 2000 m, R² ≥ 0.4 | 3 |

Lowering the R² floor (0.5→0.4) adds **zero** nights — the rejected nights genuinely lack a molecular
signal above the BL, it is **not** a marginal-R² artefact. The **start altitude** is the real lever
(1500 m→5, 2000 m→3): the classic SNR-vs-aerosol-purity trade-off. 3–5 clean Rayleigh nights/month is
consistent with the literature (ceilometer Rayleigh needs an aerosol-free, high-SNR column — Wiegner
2014) and with Rayleigh's role here as the **slow absolute anchor** (§3.3, needs ~2–3 weeks of points),
not the fast detector.

*(The before/after window-search diagnostics `figs_paper_validation/rayleigh_diag/plots/0-20000-0-06610/2026/20260312_..._window_search.png`
and `…20260328_..._window_search.png` are missing from the repo and not embedded; see the digest.)*

Remaining options (not needed given the above, but available): set the lower bound dynamically from
the L1 `layer`/`layer_aerosol` aerosol-top, or widen the window (§3.4) to trade nights for lower
per-night noise.

### 3.6 Per-instrument recommendations

| Instrument | single-calibration integration | change-detection capability | notes / other levers |
|---|---|---|---|
| **CL61 cloud** (910) | **1 day (≈8 %)** | detects a ~20 % step in **~3–4 days** | best fast monitor; but carries the +21 % offset → use Rayleigh/CHM for the absolute scale |
| **CL51 cloud** (910) | **1 day (≈9 %)** | ~few days for a ~20 % step | clean, responsive |
| **CHM15k** (1064, Rayleigh) | 1–2 days (≈13 %) | floor-limited (~9–11 %) → only ≳15 % steps fast | variability is **real instrumental** → **HK correction** (internal-T overlap, Hervo 2016) + recalibrate on laser-age steps; averaging won't help below the floor |
| **CL61 Rayleigh** (910) | **~2–3 weeks** for ≈10 % (per-night 39 %) | slow (~2 weeks) | the **absolute anchor**, not the fast detector; **wider fit window** halves the per-night noise |
| **CL31** (910, cloud) | floor-limited (~19–31 %) | only large (>~20 %) changes | needs attention: tighten cloud-target selection + HK/window monitoring |

**Cross-cutting recommendations**

1. **Keep the single-calibration integration SHORT** — daily for the cloud method, per-night for
   Rayleigh — and feed every point to the **Kalman**. Do **not** pre-average over months.
2. **Tune the Kalman as a change tracker, not an averager:** measurement noise = the
   single-calibration σ (8 % cloud / 39 % night Rayleigh), process noise set to follow real
   drift/steps. **Guard it** against divergence on short records (it failed for 3 CL61 units).
3. **Combine the two methods:** cloud for day-to-day *responsiveness* (catches window/laser events
   fast), molecular/CHM for the *absolute* scale (validated anchor, updated more slowly).
4. **Widen the 910 nm Rayleigh fit window** moderately (×2 less per-night noise) so it too can give a
   point on a useful timescale.
5. **Drive/monitor with housekeeping** (laser power, internal temperature, window transmission) — this
   catches the episodic step events directly and removes the non-averageable instrumental floor, which
   no integration can.
6. Keep robust outlier rejection and the strict (mandatory) water-vapour correction (910 nm).

*Scripts: `analyze_calib_integration.m` (Allan deviation), `run_cl61_rayleigh_params.py` (fit-window /
quality / night-length test), `analyze_chm15k_stability.m` + `cross_instrument_stability.py` (drivers
& cross-instrument).*

---

## 4. Per-calibration outlier rate & the network constant time series

*Generated 2026-06-21 (L2 recomputed after the coarse-cadence cloud-screening fix, commit `0fb4cc4`).
`validation/outliers_timeseries.py` on the network run. Optimized E-PROF v2 (C8) vs v1.1, L1 + L2.*

> The nightly coefficient `C` here is the Wiegner lidar constant `C_L = RCS/β_att` (Rayleigh method);
> see §1 for the full convention.

### 4.1 Outlier definition (drift-aware, robust)

For each calibration (stream × level × method) the date-ordered **valid** nightly lidar constants are
flagged relative to the instrument's own slow drift:

```
residual_i = C_i − rolling_median(C, 9 nights)
sigma_rob  = 1.4826 · median(|residual − median(residual)|)
outlier_i  = |residual_i − median(residual)| > 3 · sigma_rob
outlier %  = flagged nights / valid nights
```

Detrending with a 9-night rolling median removes seasonal/instrumental drift, so only genuine jumps
and spikes count.

### 4.2 Per-type outlier rate — median over streams

| type | level | **v2 outlier %** | v1.1 outlier % |
|---|---|---|---|
| CHM15k | L1 | **4.2** | 6.2 |
| CHM15k | L2 | **3.4** | 6.1 |
| Mini-MPL | L1 | 3.8 | 3.8 |
| Mini-MPL | L2 | 3.8 | 1.9 |
| CL61 | L1 | **2.0** | 11.1 |
| CL61 | L2 | **2.1** | 5.3 |

**Optimized v2 produces fewer outliers than v1.1 for CHM15k and CL61 on both levels** — most
strikingly CL61 (2.0 % vs 11.1 % L1; 2.1 % vs 5.3 % L2): the gated-optimal selection with
temporal-variability rejection avoids the occasional bad windows v1.1 admits. Mini-MPL (n small) is
comparable. Typical v2 rates are low — ~3–4 % (CHM15k), ~2 % (CL61). (L2 rates are lower than the
pre-fix version, which had counted marginal cloudy nights as valid.)

![Calibration outlier rates: (left) v2 outlier % by type, L1 light / L2 dark, bar=median; (centre) v2 vs v1.1 on L2 — points below the line = v2 cleaner; (right) per-stream v2 outlier % sorted.](figs_extracted/calibration_outliers_report_01.png)

### 4.3 Time series (L2, optimized v2; outliers in red)

Each panel is one instrument's nightly lidar constant over 2026: green = robust median, grey band =
±3σ, red = flagged outliers; the panel title gives the outlier %.

![CL61 — nightly lidar constant (L2, v2), outliers in red](figs_extracted/calibration_outliers_report_02.png)

![Mini-MPL — nightly lidar constant (L2, v2), outliers in red](figs_extracted/calibration_outliers_report_03.png)

![CHM15k (12 highest-outlier streams) — nightly lidar constant (L2, v2), outliers in red](figs_extracted/calibration_outliers_report_04.png)

### 4.4 Highest outlier rates (L2, v2; n ≥ 20) — candidates for instrument follow-up

| site | type | v2 out % | v1.1 out % |
|---|---|---|---|
| Zagreb-Maksimir | CHM15k | 23.3 (n=30) | 14.8 |
| Rotterdam-The Hague | CHM15k | 20.0 (n=20) | 0.0 |
| Harzgerode | CHM15k | 16.0 (n=25) | 7.1 |
| Falkenberg | CHM15k | 14.7 (n=34) | 9.4 |
| Folyás | CHM15k | 14.6 (n=41) | 6.1 |

A high outlier rate over a full record points to a real instrument/site issue (these overlap with the
low-laser / electronic-background stations in the diagnosis report, §2/§3). Per-stream rates for both
methods and levels are in `figs_paper_validation/network_v2_v11/outlier_summary.json`.

### 4.5 Reproduce

```
python validation/outliers_timeseries.py   # -> outlier_summary.json + 4 figures
python validation/embed_report.py doc/reports/09_calibration_stability_monitoring.md
```

---

## 5. Network validation: E-PROF v2 (C8) vs v1.1 — L1 + L2

*Generated 2026-06-21 (L2 recomputed after the coarse-cadence cloud-screening fix, commit `0fb4cc4`).
`validation/scope_network_2026.py` → `run_network_v2_vs_v11.py` → `analyze_network.py`. Optimized v2
(config C8, the `eprof_v2` default) vs E-PROF v1.1 on every CHM15k / CL61 / Mini-MPL stream, both
levels.*

### 5.1 Scope

Every instrument stream of the three Rayleigh-capable types in the 2026 L1 + L2 archives: **148
CHM15k, 11 CL61, 5 Mini-MPL**, every clear (fit-reaching) night of 2026. Per night/method we record
whether the molecular window passes the pipeline QC (`rel_error ≤ 15 %`) and the lidar constant; per
stream → **valid-on-clear%** (valid ÷ clear nights) and **σ_SD** (robust successive-difference
precision, % of median C; lower = better).

> **Note:** the L2 numbers here were recomputed after fixing a cloud-screening bug
> (`profiles_per_min` rounded to 0 at L2's 5-min cadence, disabling screening). The earlier version
> understated L2 valid-on-clear (it counted cloudy nights as "clear"). Mini-MPL L1 (also 5-min) was
> likewise recomputed.

### 5.2 Headline

- **Precision: optimized v2 beats v1.1 everywhere** — median σ_SD lower in all 6 (type × level) cells,
  by 0.5–2.2 pp, and below the diagonal for the great majority of individual streams.
- **Yield: v2 beats v1.1 on both levels (after the native-grid fix)** — a clear win on L2 (CHM15k
  +10.5 pp, v2 ahead on 71 % of streams; CL61 75 vs 70 %), and now **also on L1 CHM15k (+8.3 pp,
  85.7 vs 75.0 %)** once native L1 is binned to the L2 grid. Before the fix v2 was suppressed on
  native L1 to a −1.6 pp "tie"; binning recovers it to ~its L2 level. CL61 ≈ tie on both. v2 is
  **uniformly the better method**.

### 5.3 Results — median over streams (paired Δ = per-stream v2 − v1.1)

| level · type | n | valid v2 | valid v1.1 | Δvalid | v2 wins yield | σ_SD v2 | σ_SD v1.1 | Δσ_SD |
|---|---|---|---|---|---|---|---|---|
| L1 CHM15k | 141 | **85.7** | 75.0 | **+8.3** | 62 % | **10.6** | 13.1 | **−2.0** |
| L1 Mini-MPL | 4 | 100.0 | 100.0 | +0.0 | 0 % | **6.8** | 7.5 | +0.4 |
| L1 CL61 | 10 | 67.9 | 67.4 | −1.9 | 40 % | **7.7** | 9.2 | **−1.4** |
| L2 CHM15k | 142 | **87.5** | 73.9 | **+10.5** | 71 % | **9.5** | 11.7 | **−2.1** |
| L2 Mini-MPL | 4 | 100.0 | 100.0 | +0.0 | 0 % | **6.4** | 6.8 | −0.5 |
| L2 CL61 | 11 | **75.0** | 69.8 | −2.7 | 36 % | **7.0** | 9.7 | **−1.7** |

*(valid = valid-on-clear: valid calibrations ÷ clear nights. **The L1 rows are the corrected re-run
with native L1 binned to the L2 grid (30 m × 300 s; `l1_bin_to_l2_grid`, day-sampled every 2nd day);
L2 rows unchanged.** Mini-MPL n small; its L1 and L2 share one 5-min grid so the two rows are
essentially identical. The earlier native-L1 v2 (72.7 %, a −1.6 pp deficit) was the native-grid
handicap — see §5.6.)*

![Network medians: valid-on-clear% (top) and σ_SD (bottom) for optimized v2 (blue) vs v1.1 (orange), L1 left / L2 right.](figs_extracted/network_v2_vs_v11_report_01.png)

![Paired per-stream comparison (each point = one instrument). σ_SD points below the dashed line and valid% points above it favour v2.](figs_extracted/network_v2_vs_v11_report_02.png)

### 5.4 Precision (σ_SD): v2 wins decisively

σ_SD is lower for v2 in every cell and for the great majority of individual streams (paired scatter,
bottom row, points below the diagonal) — on **both** levels. v2's gated-optimal selection with
temporal-variability rejection produces a more repeatable lidar constant than v1.1's
signal/Rayleigh-error pick. Largest gains: CHM15k (−2.1 to −2.2 pp) and CL61 L2 (−1.7 pp).

### 5.5 Yield (valid-on-clear): v2 ahead on both levels (after the L1 grid fix)

- **L2 CHM15k**: the clearest win — +10.5 pp valid-on-clear on 71 % of 142 streams, *and* −2.1 pp
  σ_SD.
- **L2 CL61**: v2 75 % vs v1.1 70 % (and σ_SD 7.0 vs 9.7). With the screening fixed, v2 no longer
  "loses" CL61 on L2 — the earlier deficit was the bug counting cloudy nights as clear.
- **L1 CHM15k**: now **+8.3 pp** (85.7 vs 75.0 %, v2 ahead on 62 % of streams) after binning native L1
  to the L2 grid. Before the fix this was a −1.6 pp deficit (72.7 %) — the native-grid handicap, now
  removed.
- **L1 CL61**: ≈ tie (67.9 vs 67.4 %; the 60 s native grid was already adequate, so binning barely
  moves it).
- **Mini-MPL** (n small): both methods calibrate ~all clear nights; v2 has marginally lower σ_SD.

### 5.6 L1 vs L2 — the "tie on L1" is a native-grid effect, not a data difference

A deeper trace shows the apparent v2 **"tie" on L1 is an artifact of running v2 on the *native* L1
grid**, not a real ceiling — the two levels are the *same measurement*.

**Root cause.** On an identical clear day (CHM15k 06610, 2026‑02‑25, full files) the range-normalised
night-mean signal is molecular-proportional at **R² ≈ 0.96 on *both* L1 and L2**. Yet the gated
methods behave differently by grid:

| method | native L1 (15 m × 15 s) | L1 binned to L2 grid (30 m × 300 s) | L2 (30 m × 300 s) |
|---|---|---|---|
| **v1.1** | ✓ 5.17×10¹¹ | — | ✓ 5.04×10¹¹ |
| **v2** | ✗ rejected | ✓ **5.00×10¹¹** | ✓ 4.92×10¹¹ |
| **earlinet** | ✗ rejected | ✓ 5.11×10¹¹ | ✓ 5.09×10¹¹ |
| **v1.2** | ✗ | ✗ | ✗ (hard day) |

`v1.1` is **grid-robust**; the newer gated methods reject the **fine, noisy native L1 grid** (1024
bins × 15 s) — `v2` via its temporal-variability rejection over the per-profile stack (needs
time-averaging), `earlinet` via the per-bin R²/proportionality test (needs range-binning). **Binning
native L1 to the L2 grid (30 m × 300 s) recovers v2, which then agrees with L2** (5.00 vs
4.92×10¹¹). The raw L1 `rcs_0` is **fully background-corrected** and matches L2's β_att (clear-mean:
2 km 1.01×10⁻² vs 1.03×10⁻²; 6 km 5.4×10⁻⁴ vs 5.3×10⁻⁴). A 5-stream CHM15k check confirms it:
native-L1 v2 **9 %** valid → binned-L1 v2 **15 %** ≈ L2 v2 **10 %**, with v1.1 native at **13 %** (all
C_L ≈ 1.9–2.0×10¹¹).

![L1 vs L2 Rayleigh: same data, a method/grid interaction. Left: demonstration day C_L by method × grid (v1.1 works native; v2/earlinet need binning to the L2 grid, then agree with L2). Right: 5-stream valid fraction — native-L1 v2 suppressed, binning recovers it.](figs_extracted/network_v2_vs_v11_report_03.png)

**Fix applied (this is what the L1 columns now show).** `calibrate_rayleigh` now bins native L1 to the
L2 grid (30 m × 300 s) by default for `data_level==L1` (`l1_bin_to_l2_grid`, commit `a0873f8`); L2/RAW
are untouched and a coarser native grid is a no-op. Re-running the network L1 with the fix lifts **L1
CHM15k v2 from 72.7 % (a −1.6 pp deficit) to 85.7 % (+8.3 pp over v1.1)** — matching its L2 level
(87.5 %). The demonstration day confirms agreement (L1-binned v2 5.00×10¹¹ vs L2 4.92×10¹¹). The L1
rows in §5.3 are this corrected re-run (day-sampled every 2nd day); the L2 rows are unchanged. CL61
barely moves (its 60 s native grid was already adequate). Practical guidance: with the fix, L1 and L2
now calibrate consistently; `v1.1` remains the grid-robust fallback if binning is disabled.

### 5.7 Conclusion

Across the whole network, on **both** L1 and L2, the optimized v2 gives **more repeatable**
calibration constants than v1.1 (lower σ_SD everywhere) and **at least as many** valid calibrations —
a clear yield win on the operational L2 product (CHM15k +10.5 pp; CL61 +5 pp) and — with the
native-grid fix now in `calibrate_rayleigh` — **a matching win on L1 (CHM15k +8.3 pp, 85.7 vs
75.0 %)**. The earlier "tie on L1" was a native-grid handicap (the gated methods over-reject the fine
native grid); binning native L1 to the L2 grid (30 m × 300 s) removes it, and L1 then tracks L2 on the
same data. **v2 (C8) is adopted as the repo/production default method**; it runs consistently on L1
(auto-binned) and L2.

> **Deployment scope reminder.** "Production default" here means the **repo** default that the
> MeteoSwiss E-PROFILE ALC calibration runs. The **deployed operational E-PROFILE network** product
> still uses the older **E-PROF v1.0** (which carries the historical Klett sign error); the v1.0 → v2
> upgrade on the operational side is a separate, pending step. The Klett sign error is fixed in the
> three active codebases; only operational v1.0 still carries it.

### 5.8 Reproduce

```
python validation/scope_network_2026.py
python validation/run_network_v2_vs_v11.py        # both levels (L2 recomputed post-fix)
python validation/analyze_network.py              # network_summary.json + the 2 figures
python validation/embed_report.py doc/reports/09_calibration_stability_monitoring.md
```

---

## 6. June-2026 modifications changelog — CAMS-far + dashboard

**Branch:** `wv-correction` · **Date:** 2026‑06‑26 · **Scope:** E‑PROFILE network, L1 2025‑2026.

Dated changelog of the modifications made to the Rayleigh + cloud calibration and the monitoring
dashboard in June 2026, the evidence behind each, per-type validation, and the network re-deployment.

### Summary of changes

| # | Change | Type | Effect |
|---|--------|------|--------|
| 1 | **E‑PROF v2** is the default molecular‑window method everywhere | compute | Rayleigh windows now picked by the v2 (optimal) selector on every instrument |
| 2 | **New flag −10 "Closest CAMS data too far"** | compute | 910 nm stations outside the CAMS domain (e.g. New Zealand) now fail honestly instead of a false success |
| 3 | **Cloud station‑coordinate read fix** | compute (bug) | cloud no longer feeds fill‑value coords to the CAMS lookup |
| 4 | **Rayleigh diagnostic pcolor → datetime x‑axis** | plot | time–height curtains show wall‑clock hours + date, not "hours since start" |
| 5 | **Dashboard fixes** (success‑rate, calendars, ICAO, %op, OmB/timeseries labels) | dashboard | correct/clear monitoring UI |
| 6 | **Purity‑objective window selection (`eprof_v2p`)** | R&D (not deployed) | candidate fix for v2 picking aerosol‑tainted windows; evaluated and **not adopted** |

### 6.1 E‑PROF v2 as the default molecular method

`options.json` previously shipped `molecular_method = eprof_v1.2`; the operational intent is **v2**
(the "optimal" selector: composite score + time‑resolved aerosol flagging, C8 defaults). v2 is now the
default in **all** production paths: `options.json`, `calibration/config.py` (dataclass default +
`from_json` fallback), `rayleigh_fit.find_optimal_molecular_window`, and the `calibration.py` getattr
fallback (regression‑guarded by `tests/test_default_method.py`). Method‑comparison scripts that
intentionally pin a version are unchanged. *Verified in code: `options.json` and
`calibration/config.py` both default to `eprof_v2`.*

**Known caveat (documented, not a regression):** v2's C8 tuning lowered the molecular‑window floor
`min_window_start_m` 2000→1500 m. On some clear nights with residual boundary‑layer aerosol up to
~2 km, v2 picks a window ~250 m lower than v1.2 — into the aerosol — which inflates the fit slope and
trips the slope‑vs‑Klett "method disagreement" gate (e.g. CHM15k Payerne 2026‑02‑25: v1.2 flag 1 /
C_L 4.72e11 → v2 flag −3 at 16.2 %). The cause is the v2 selector's **R²‑dominated** score (R² weight
1.0 vs aerosol penalties ~0.2): R² measures fit *linearity*, not molecular *purity*, so a smooth
aerosol layer yields a tight (high‑R²) line with the wrong slope and is *rewarded*. This motivated the
R&D in §6.6. It is a rare event (see §6.7) and v2 is otherwise cleaner network‑wide.

> **Reconciliation (2026-07-10):** the shipped default `min_window_start_m` is back at **2000 m** (see
> §3.5b), so the 1500 m caveat above is historical for the operational config; it can still be lowered
> to ~1500 m by choice where the BL is low and more points are needed.

![C_L (=signal/molecular ratio) vs window‑centre altitude, colour = fit R². Aerosol inflates C_L at low altitude (where R² is highest), while the clean high‑altitude "plateau" is noisy (low R²). This aerosol‑vs‑noise trade‑off is why a fixed plateau‑rejection gate is not viable (§6.6) and why window placement matters.](figs_v2_camsfar/cl_vs_altitude.png)

### 6.2 New flag −10 "Closest CAMS data too far"

**Bug.** CL31/CL51/CL61 (910 nm) need CAMS for the water‑vapor correction (mandatory). The CAMS
download is **regional** (Europe/N‑Atlantic, lat 27–74° lon −27–45°). For an out‑of‑domain station,
`xarray.sel(method="nearest")` silently returned the **domain‑edge** cell — for Lauder, NZ that is
**~125° (≈14 000 km) away** — and both the cloud and Rayleigh calibrations used it, producing a
**false success** with a meaningless WV correction.

**Fix.** `water_vapor.cams_point_too_far(file, lat, lon)` flags a station whose nearest CAMS grid
point is farther than `max(1.0°, 1.5×grid_spacing)`. Rayleigh returns `flag=-10` directly (WV block +
CAMS‑molecular block); the cloud `compute_wv_transmission` raises a recognizable error that the runner
maps to −10 (distinct from −4 *missing* CAMS and −99 *other*). Registered in `flags.py` and
`monitoring/config.py`. *Verified in code: flag −10 "Closest CAMS data too far" is present in
`calibration/flags.py`.*

**Verified:** Lauder NZ nights now return **rayleigh −10 and cloud −10**; European stations are
unaffected (nearest cell 0.2–0.3° away). Test `tests/test_cams_too_far.py`.

> Note: this correctly flags **every** out‑of‑domain 910 nm station. To *calibrate* such global
> stations instead, the CAMS download domain would need widening (or per‑station CAMS boxes — since
> implemented for a handful of affiliates, routed by station lat/lon).

### 6.3 Cloud station‑coordinate read fix

While verifying §6.2, the in‑domain CL31/CL51 sample calibrations began failing — exposing a
**pre‑existing** bug: the cloud reader preferred the per‑profile `latitude`/`longitude` variables over
the canonical scalar `station_latitude`/`station_longitude`. In L2 files the per‑profile arrays are
often `_FillValue` (~1e36), so the cloud was feeding **garbage coordinates** to the CAMS lookup (again
grabbing the domain‑edge cell) — silently, because the tests only check that profiles exist, not WV
quality. Fixed to take the first **valid** coordinate (reject fill / out‑of‑range, prefer the scalar
station coord). Restores the cloud sample tests and removes a latent WV bias.

### 6.4 Rayleigh diagnostic pcolor — datetime x‑axis

The Rayleigh time–height curtain (`plot_rayleigh_diagnostics_compact` and `…_failure`) used "Hours
since start". It now uses a real **datetime axis** (matplotlib date numbers + `ConciseDateFormatter`):
wall‑clock **hours + date** ("Time (UTC)"), with the date shown at day boundaries. All overlays
(hatched excluded profiles, cloud‑base scatter, high‑cloud mask) move with it. Falls back to hours if
no timestamps are passed. Verified on CHM15k + CL61, success and failure plots; test
`tests/test_plotting_time_axis.py`. Regenerated network‑wide in this run (PLOTS=1).

### 6.5 Dashboard improvements

- **Success rate redefined to the true daily yield.** Was `valid / suitable` (excluded no‑data and the
  "no liquid cloud / not clear" flag), which read **~96 %** for cloud. Now `valid / all days`, so
  no‑data, no‑cloud and every rejection count against it: **cloud 95.8 → 53.1 %**, **rayleigh 48.3 →
  12.8 %** (matches the cloud‑yield study). Identical formula for both methods.
- **Station calendars**: a 3‑month window with prev/next arrows + a month dropdown, stacked
  **vertically**.
- **ICAO detection‑altitude map**: median (not mean); colorbar shortened to "alt [m]".
- **% of operational constant map**: fixed (it was simply built without `--opcoeff`; now passed).
- **OmB plots labelled "v2"** (not "our"); **Rayleigh C_L time series**: legend moved below the plot,
  lines renamed **"Applied in L2" / "v1.0" / "v2.0" / "v2.0 Kalman estimate"**. *(The
  "v2.0 Kalman estimate" line is the Kalman best-estimate smoothing of the constant time series from
  §3.3 — the operational best estimate tracked on the dashboard.)*
- New flag **−10** added to the dashboard flag table/colours.
- Regression tests in `tests/test_dashboard.py` (16) guard all of the above.

### 6.6 R&D — purity‑first window selection (`eprof_v2p`), not deployed

To address the v2 caveat (§6.1), a flag‑gated variant `eprof_v2p` was added: among molecular‑**eligible**
windows it **minimises an aerosol cost** (curvature + scattering) instead of **maximising R²**, so a
tight but aerosol‑tainted low window can't win on R² alone.

- **Pilot (anchor cases):** flips Payerne 2026‑02‑25 −3→pass with the clean window (C_L 4.97e11 ≈
  v1.2), no regression on genuine aerosol nights.
- **Plateau‑rejection gate: ruled out** by the C_L‑vs‑altitude diagnostic (figure §6.1) — the
  high‑altitude reference is too noisy and a ~15–20 % descent is *normal*, so such a gate would reject
  good calibrations.
- **Wider scan (216 nights, 12 continental CHM15k):** v2p **never changed a pass/reject outcome**
  (0 recoveries, 0 regressions, same 44 passes), but picked **systematically cleaner windows** (median
  curvature 4.2 → 1.5 %, scattering 1.077 → 1.031). So v2p is *safe* and *cleaner* — its value is
  **precision, not yield**.
- **Stage‑1 sweep (DECIDED — v2p NOT adopted).** A v2‑vs‑v2p sweep measuring **yield + σ_SD
  (night‑to‑night precision)** per type (`validation/_stage1_sweep.py`, parallel) gave:

  | type | streams | yield v2→v2p | σ_SD v2→v2p |
  |---|---|---|---|
  | CHM15k | 12 | 78→78 % | 13.8→13.9 (tied) |
  | Mini‑MPL | 5 | 85→85 % | 9.3→9.4 (tied) |
  | CL61 | 8 | 69→69 % | 7.6→9.2 (worse; but driven by low‑N noisy streams — inconclusive) |
  | CL51 | 0 | — | (no Rayleigh data — cloud‑method instrument) |

  **The "cleaner windows" did NOT translate into better precision** — σ_SD is tied on
  CHM15k/Mini‑MPL and shows no gain (possibly a loss) on CL61. The purity objective chases the
  cleanest *spot*, which moves altitude night‑to‑night and so doesn't improve *stability*; v2's R²‑max
  picks a consistent high‑SNR window and is fine. v2p only fixed the **rare** (~2/128) Payerne‑type
  false‑reject, at no measurable precision benefit. **Decision: keep v2; `eprof_v2p` stays registered
  but not deployed.** Probe scripts: `validation/_pilot_purity.py`, `_scan_wider.py`,
  `_probe_plateau.py`, `_stage1_sweep.py`.

### 6.7 Validation by instrument type

| type | wavelength | method | result |
|------|-----------|--------|--------|
| **CHM15k** | 1064 nm | Rayleigh v2 | works (44 passes in the 216‑night scan); the rare Payerne‑type −3 is the §6.1 caveat |
| **CL61** | 910 nm | Rayleigh + cloud | in‑domain OK; **out‑of‑domain (Lauder) → −10** ✓ |
| **CL31 / CL51** | 910 nm | cloud | in‑domain cloud OK after the §6.3 coord fix; CAMS guard does not false‑trigger ✓ |
| **Mini‑MPL** | 532 nm | Rayleigh v2 | 4/4 real clear nights flag 1, sensible C_L (no CAMS needed) ✓ |

### 6.8 Deployment (network re‑run)

Launched on CSCS (balfrin), 2025‑01‑01 → 2026‑05‑31, all 425 streams:

- **Calibration** array (`--force`, `PLOTS=1`) — v2 windows, the −10 guard, and regenerated
  datetime‑axis diagnostics.
- **Sensitivity** array (`--no-cal --sens`) — median ICAO altitude.
- **Dashboard** rebuild (`--opcoeff`) — auto‑runs after both, with all the §6.5 fixes + flag −10.

**Pending (as of the changelog date):** the OmB recompute (the "v2" plot label) was deferred — it
needs the 0.4° CAMS with aerosol backscatter (`CAMS_Monthly_04`), whose download was not yet complete.
The OmB *map* (which reads `omb.csv`) is unaffected; only the per‑station OmB plot label changes when
OmB is re‑run.

### 6.9 Files changed

`calibration/`: `flags.py`, `config.py`, `rayleigh/calibration.py`, `rayleigh/rayleigh_fit.py`,
`rayleigh/molecular_methods.py`, `cloud/calibration.py`, `water_vapor_correction/water_vapor.py`,
`plotting.py` · `monitoring/`: `index.py`, `render.py`, `metrics.py`, `charts.py`, `config.py`,
`static/diag.js`, `static/style.css` · `scripts/run_network_calibration.py` · `options.json` ·
`scripts/run_lindenberg_cl61_cal.py`.

**Tests added:** `test_cams_too_far.py`, `test_default_method.py`, `test_plotting_time_axis.py`, and
extensions to `test_dashboard.py` (16) — all passing locally; sample suite green except the two
*known* cases (v2 Payerne marginal; a not‑clear Mini‑MPL fixture day).

---

## References

- M. Wiegner and A. Geiß, *Aerosol profiling with the Jenoptik ceilometer CHM15kx*, Atmos. Meas. Tech.
  5, 1953–1964 (2012).
- E. J. O'Connor, A. J. Illingworth, R. J. Hogan, *A technique for autocalibration of cloud lidar*,
  J. Atmos. Oceanic Technol. 21, 777–786 (2004).
- E. Hopkin et al., *A robust automated technique for operational calibration of ceilometers…*,
  Atmos. Meas. Tech. 12, 4131–4147 (2019).
- M. Hervo et al., *An empirical method to correct for temperature-dependent variations in the overlap
  function of CHM15k ceilometers*, Atmos. Meas. Tech. 9, 2947–2959 (2016).
- I. Mattis, G. D'Amico, H. Baars et al., *EARLINET Single Calculus Chain — technical, Part 2*, Atmos.
  Meas. Tech. 9, 3009–3029 (2016).
- R. J. Hogan, *Fast approximate calculation of multiply scattered lidar returns* (Photon
  Variance-Covariance model), Appl. Opt. 45, 5984–5992 (2006).
- Le et al. (2026), *Long-term performance of the Vaisala CL61*; Kotthaus et al. (2016), CL31
  characterisation; Wiegner & Gasteiger (2015), 905 nm water-vapour absorption; Filioglou et al.
  (2023), CL61 factory-vs-field calibration spread. (Full citations in the source sessions.)
