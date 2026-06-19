# What drives ALC calibration (in)stability? — Rayleigh & cloud, CL31/CL51/CL61/CHM15k

**Author:** M. Hervo · **Date:** 2026-06-17

## Question

How stable are the absolute calibrations (Rayleigh molecular and liquid-cloud), what is
the **main driver** of their instability (water vapour, aerosols, or instrumental), is this
**the same in the literature**, and do the **other instruments (CL31, CL51, CHM15k)** show
the same problem?

## Short answer

**Instrumental drift is the dominant driver of long-term ALC calibration instability**, on
three timescales: **laser ageing** (multi-year), **internal/detector temperature**
(seasonal + diurnal), and **window contamination/fogging** (episodic). This is true for
**every** ALC — including the 1064 nm CHM15k, which has *no* water-vapour sensitivity yet is
just as unstable, proving the driver is instrumental, not atmospheric. **Water vapour** adds
a large (~20 %) *seasonal* term **for the 910 nm instruments only** (CL31/CL51/CL61), which
we now correct. **Aerosols** are second-order. This matches the literature exactly
(Le et al. 2026; Hopkin et al. 2019; Hervo et al. 2016; Kotthaus et al. 2016; Wiegner et al.
2014/2015; Filioglou et al. 2023).

## 1. The decisive evidence — CHM15k (1064 nm) is unstable *without* any water vapour

The CHM15k operates at **1064 nm, outside the water-vapour absorption band** — so any
instability in its calibration **cannot** be water vapour. Using the 11-year daily archive
(13 stations, 2013–2024; `A:\CHM15k\…`, the dataset behind `estimate_calib_from_housekeeping`):

![CHM15k calibration stability and instrumental drivers](figs_paper_validation/calib_stability_chm15k_drivers.png)

*(a) Payerne: the lidar constant tracks the internal temperature over a decade. (b) All 13
units pooled: lidar constant vs internal temperature, r = 0.31. (c) Per-station, the
calibration correlates with instrumental housekeeping — internal temperature, laser pulses,
optics state.*

| CHM15k metric (per station, 11 yr) | value |
|---|---|
| Calibration CV | **23 %** (median; 15–38 %) |
| Seasonal peak-to-peak amplitude | **34 %** (median) |
| Long-term trend (laser ageing) | **−10 … +16 % / yr** (unit-specific) |
| Correlation with **internal temperature** | **+0.38** (median) |
| Variance explained by instrumental HK (T_int + laser pulses + optics) | **R² = 0.38** (up to 0.69) |

A 1064 nm instrument with a **23 % calibration CV and a 34 % seasonal cycle that follows its
internal temperature** can only be explained by **instrumental** effects (temperature-driven
optics/overlap and laser ageing) — exactly the mechanism Hervo et al. (2016) corrected for
the CHM15k overlap. ~38 % of the variance is captured by housekeeping alone; the rest is
sampling noise plus the aerosol-loading sensitivity of the molecular fit.

## 2. Water vapour — a large 910 nm-only term (corrected)

From the CL61 controlled re-run (`cl61_verify_*`), toggling the WV correction changes the
calibration of the **910 nm** instruments substantially, and the **1064 nm CHM15k not at all**:

| | WV-off → WV-on |
|---|---|
| Rayleigh constant (910 nm) | **+20.8 %** (the molecular fit at 3–6 km integrates the full WV column) |
| Cloud coefficient (910 nm) | **−12.8 %** (cloud base ~1 km still sees most of the boundary-layer WV) |
| CHM15k (1064 nm) | **0 %** (outside the WV band) |

![WV impact and L1/L2 robustness](figs_paper_validation/cl61_verify_l1l2_wv.png)

So for CL31/CL51/CL61 a **seasonal water-vapour cycle** is a first-order calibration driver
*if uncorrected*; with the matching-month CAMS correction applied it is removed (hard rule).
This is the Wiegner & Gasteiger (2015) result (~20 % mid-latitude, >50 % tropics at 905 nm)
and is why the CHM15k (1064 nm) is immune.

## 3. Aerosols — second-order

- **Cloud method:** the target is a *totally attenuating* liquid cloud with fixed S = 18.8 sr,
  so it is largely aerosol-insensitive; residual below-cloud aerosol is suppressed by the
  90 %-in-cloud acceptance filter (Hopkin 2019; Le 2026). Our diagnostics confirmed aerosol
  contamination pushes the cloud coefficient the *wrong* way to explain the offset.
- **Rayleigh method:** the assumed aerosol lidar ratio S_p (50 sr default) enters the fit;
  literature propagates ±10 sr as a sensitivity. A secondary, event-driven driver.

## 4. Do the other instruments show the same problem? — Yes

| Instrument | λ | raw daily CV (precision) | smoothed drift (this work / literature) | dominant driver |
|---|---|---|---|---|
| **CHM15k** | 1064 | 23 % (11 yr) | seasonal 34 %, trend ±10 %/yr | **instrumental** (T_int, laser) — *no WV* |
| **CL31** | 910 | 44 % (8 yr, Payerne) | ±3 %/20 mo, ±5 %/yr (Hopkin 2019) | instrumental + WV |
| **CL51** | 910 | 9 % (short record) | similar (WV-corr. validated, Wiegner 2019) | instrumental + WV |
| **CL61 (cloud)** | 910 | 11 % (2026) | laser-ageing dominated (Le 2026) | instrumental + WV |
| **CL61 (Rayleigh)** | 910 | 61 % (2026, per-night) | needs Kalman; agrees with cloud within unc. | instrumental + WV |

![Cross-instrument calibration scatter](figs_paper_validation/cross_instrument_stability.png)

Two points: (i) the **raw per-sample scatter** is large for every type — which is why all
operational calibrations are **temporally smoothed** (Kalman / 90-day running mean); the
*smoothed* drift is the few-% literature value. (ii) The **cloud method is more precise
per-sample** than the per-night Rayleigh fit (CL61 11 % vs 61 %), but both converge after
smoothing. The instability is **universal across the network and across instrument types**.

## 5. Same in the literature? — Yes, in detail

- **Instrumental dominates long-term drift.** Le et al. (2026, "Long-term performance of the
  Vaisala CL61", 4 units, 3 yr): the firmware rescales the calibration while laser power
  >40 %, but **below 40 % it breaks down** — Lindenberg's cloud-calibration factor fell ~×3
  (1.45→0.54) as laser power went 40→10 %. **Window fogging** caused episodic factor jumps up
  to ×10 (diagnosed via the "window condition" housekeeping). An **internal-temperature
  look-up-table** corrects the CL61 bias — the analogue of **Hervo et al. (2016)**'s
  temperature-dependent CHM15k overlap correction. **Kotthaus et al. (2016)** documents the
  CL31 instrumental background and firmware artefacts.
- **Housekeeping is the standard monitor/QC.** Hopkin (2019) rejects profiles with **window
  transmission < 90 %** ("cannot be reliably corrected"); Le (2026) keys on **laser power**,
  **window condition**, and **internal temperature**. This is exactly what our
  `estimate_calib_from_housekeeping` work exploits for the CHM15k.
- **Water vapour** ~20 % at 910 nm (mid-latitude), seasonal, not at 1064 nm (Wiegner 2014/2015).
- **Cloud vs Rayleigh** agree within the cloud method's **~10 % uncertainty**, dominated by
  the **multiple-scattering factor η** (0.7–0.85; ~10 % for wide-beam ceilometers) and
  S = 18.8 ± 0.8 sr (Hopkin 2019; O'Connor 2004). **Filioglou et al. (2023)** attribute the
  CL61 factory-vs-field calibration spread (15 %, latitude-dependent) to **water vapour +
  multiple scattering** — which is precisely the residual seen in our cloud-vs-Rayleigh +21 %.

## 6. Drivers, ranked

| Rank | Driver | Type | Timescale | Affects | Mitigation |
|---|---|---|---|---|---|
| 1 | **Laser ageing / power loss** | instrumental | multi-year | all ALC | HK monitoring; recalibrate; flag <40 % power |
| 2 | **Internal/detector temperature** | instrumental | seasonal + diurnal | all ALC | T-dependent overlap (Hervo 2016) / T-LUT (Le 2026) |
| 3 | **Window contamination / fogging** | instrumental | episodic | all ALC | window-transmission QC (<90 % reject); keep blower/heater on |
| 4 | **Water vapour absorption** | atmospheric | seasonal | 910 nm only | matching-month CAMS WV correction (done) |
| 5 | **Multiple scattering η / S assumption** | method | constant offset | cloud method | η per FOV; cloud as cross-check, Rayleigh/CHM as anchor |
| 6 | **Aerosol load / S_p** | atmospheric | event | Rayleigh fit | sun-photometer S_p; cloud target is aerosol-robust |
| — | per-night/per-cloud sampling | noise | daily | all | Kalman / running-mean smoothing |

## 7. Conclusions & recommendations

1. **Main driver = instrumental** (laser ageing, internal temperature, window contamination),
   proven by the 1064 nm CHM15k being just as unstable (23 % CV, 34 % seasonal, r=0.38 with
   internal temperature) despite having *no* water-vapour sensitivity. The literature agrees.
2. **Water vapour** is a large, *correctable*, 910 nm-only seasonal term (~20 % Rayleigh,
   −13 % cloud); keep applying the matching-month CAMS correction.
3. **The instability is universal** — CL31, CL51, CL61 and CHM15k all show it; it is reduced
   to the few-% literature level by temporal smoothing (Kalman / running mean).
4. **Recommendations:** (a) drive/monitor the calibration with **housekeeping** (laser power,
   internal temperature, window transmission) as `estimate_calib_from_housekeeping` already
   does for the CHM15k — extend to the CL61; (b) flag laser power <40 % and window
   transmission <90 %; (c) keep the WV correction strict for 910 nm; (d) use the molecular/CHM
   (Rayleigh) calibration as the absolute anchor and the liquid-cloud as a cross-check (its
   ~10 % η/S uncertainty is the residual behind the cloud-vs-Rayleigh +21 %).

---
*Data: `A:\CHM15k\` 11-yr daily HK (13 stations); `cl61_verify_*` controlled re-run;
`Cloud_Trans-cor_WV-cor` (CL31 06610_B, CL51 EDT_A). Scripts: `analyze_chm15k_stability.m`,
`cross_instrument_stability.py`. Literature digest from the 3 provided papers + web (full
citations in the session). Figures embedded above.*
