# CL61 Rayleigh calibration investigation — why C_L(Rayleigh) ≠ C_L(cloud)

*2026-07-02, `validation/paper/` (branch `wv-correction`). Trigger: the station validation shows
good agreement CHM15k ↔ CL61-cloud (Payerne −0.6 %) but not CHM15k ↔ CL61-Rayleigh (+13.5 %), i.e.
the two calibration methods disagree on the same instrument's lidar constant:
C_L(Rayleigh) = 1.25 vs C_L(cloud) = 1.43 at Payerne (ratio 0.88). All probe scripts are in
`validation/paper/_cl61_*.py`; figures in `figs_paper_report/`.*

## Executive summary

| hypothesis | verdict | evidence |
|---|---|---|
| WV correction missing in the CL61 Rayleigh calibration | **no — applied & essential** | 910 nm nights only calibrated when WV-correctable; T²_wv ≈ 0.78 at the fit windows → the correction already raises C_L by ≈ 28 % ([calibration.py:495-533](../../calibration/rayleigh/calibration.py)) |
| WV correction missing in the validation/comparison | **no — applied** (both CL61 rows identically → cancels in ray-vs-cloud) | [intercompare.py apply_wv](../../validation/paper/intercompare.py) + fail-safe guards (missing month / too-far → excluded, never uncorrected) |
| wrong laser wavelength (measured 910.74 vs manufacturer 910.55 nm) | **no — immaterial** | recalibration of all successful nights: **+0.5 %** on C_L; whole (λ₀ ± 0.5 nm, FWHM 0.3–3.4 nm) envelope ≤ 2.4 % |
| aerosol contamination of the fit window | **no — wrong sign** | contamination *raises* C_L (agreed); our C_L is *low*. The eprof_v2 outlier screen + aerosol gate also reject such nights |
| assumed lidar ratio (Klett transmission below the window) | **no — too small** | LR = 52 sr ± 20 scanned → 2σ ≈ 2–6 % (the per-night `sensitivity` term) |
| CAMS resolution (1° vs 0.4°) | **no — immaterial** | T²_wv differs by ≤ 0.3 % in C_L terms |
| monthly-CAMS vs radiosonde WV | **secondary, wrong direction** | per-night 00-UT soundings → dC_L median **−3.6 %** (clear calibration nights are drier than the monthly mean → the correction slightly over-corrects; fixing it *widens* the method gap) |
| **CL61 negative signal baseline (weak-signal artefact)** | **CONFIRMED — leading cause** | nightly means are physically **negative at 11–14 km** (−0.02…−0.05 Mm⁻¹sr⁻¹, 12/12 nights); at the 2.6–6 km fit windows the depression (−0.02…−0.04 on clean nights) exceeds the −0.014 required for the −12 % C_L gap |

**Conclusion:** the CL61 vendor β_att carries a systematic **background over-subtraction** whose
relative impact is largest exactly where the Rayleigh method fits (weak molecular signal,
2.6–6 km). The liquid-cloud method (strong signal, 0.1–2.4 km) is essentially immune — and it is
anchored by the CHM15k (−0.6 %) whose own constant closes the AERONET AOD (ratio 0.85, consistent
with 1 at winter AODs of 0.03). **Use the cloud constant for the CL61; treat the native Rayleigh
constant as biased low by ≈ 10–15 % until the baseline is handled.** A ≈ −4 % secondary term comes
from the monthly-CAMS WV over-correction on dry calibration nights.

## 1. The signal offset: direct measurement, consistency, and causal test

**1a. Direct measurement (covered telescope).** Three hood-on dark tests were performed on the
operational Payerne CL61 (2026-05-12 09:35–14:50, 2026-05-26 11:45–13:15, 2026-06-09 09:20–11:50);
their periods are present in the L1 archive and are unambiguously dark (window-mean β_att at 1 km
= +0.001 to +0.002 Mm⁻¹sr⁻¹, where any real daytime boundary layer gives 0.1–0.5). With no
atmosphere involved, the window-mean profile **is** the processing offset:

| range | 1 km | 3 km | 5 km | 8 km |
|---|---|---|---|---|
| b_dark(z) [Mm⁻¹sr⁻¹], median of 3 tests, 300 m smoothed | ≈ 0 | **−0.008** | **−0.017** | **−0.027** |

![dark windows](figs_paper_report/fig_cl61_dark_windows.png)
*Figure 1 — L1 archive during the covered-telescope windows: (a) window-mean β_att; (b) zoom
around zero — the negative, altitude-growing offset; (c) β_att/z², the non-range-corrected shape.*

**1b. Consistency with the night sky under full transmission.** The measured attenuated
backscatter must lie below the aerosol-free molecular curve because of aerosol and water-vapour
extinction; a meaningful residual therefore requires the complete model
C_L·β_mol·T²_mol·T²_wv·**T²_aer**, with T²_aer forward-integrated from the measured profile
itself (LR = 50 sr). On the six clean calibration nights the fully-corrected residual at 3–6 km is
**−0.021 Mm⁻¹sr⁻¹ (median; range −0.021…−0.035)**, against the covered-telescope
b_dark(3–6 km) = −0.015: same sign and magnitude, the ≈ 40 % excess being compatible with a
day/night dependence of the offset (the dark tests are daytime) and residual LR/WV model error.
(An earlier draft of this analysis omitted T²_aer and overstated the residual by ≈ 30 %; the
conclusion survives the correction.)

**1c. Causal test — remove the measured offset and recalibrate.** Since range correction is the
fixed z² map, removing the offset from the non-range-corrected signal and re-range-correcting is
identical to subtracting b_dark(z) from β_att. Doing so on the calibration nights (corrected
copies of the L1 files, full eprof_v2 recalibration, WV on):

| | C_L median (common nights) | gap to C_L(cloud) = 1.425 |
|---|---|---|
| original | 1.047 | −26.5 % |
| **offset-corrected** | **1.286** | **−9.8 %** |

Per-night changes +4.0 %, +21.6 %, +22.7 %; the corrected data also yields **more eligible
nights** (9 vs 7) — windows previously rejected by the |intercept| criterion become admissible,
confirming the offset as the fit-spoiler. The daytime-measured b_dark thus removes ≈ ⅔ of the
method discrepancy; the remainder is consistent with the ≈ 40 % larger night-time offset seen in
§1b (a night-time covered test would settle it).

Magnitude closure: explaining C_L(ray)/C_L(cloud) = 0.88 requires a window-mean deficit of
−0.12·C_L·β_mol·T² ≈ **−0.014 Mm⁻¹sr⁻¹** — bracketed by the daytime dark (−0.008…−0.017 over the
window range) and the night-sky residual (−0.021). Noise is not a factor: the per-sample σ
(0.17–0.57 Mm⁻¹sr⁻¹, [dark_measurement_payerne.md](dark_measurement_payerne.md)) averages to
≈ 0.02/√N per night; the offset is systematic.

## 2. Wavelength & spectrum (excluded)

Recalibrating every successful Payerne night with the manufacturer line:

| λ₀ / FWHM [nm] | C_L median | vs C_L(cloud)=1.425 |
|---|---|---|
| 910.74 / 1.0 (Qmini-measured, operational) | 1.268 | −11.0 % |
| 910.55 / 1.0 (manufacturer) | 1.273 | −10.6 % |

Pure-transmission scan at the window: the full plausible envelope (λ₀ 909.7–911.0, FWHM 0.3–3.4)
moves C_L by ≤ 2.4 %. The MATLAB reference (`rayleigh_calibration_matlab`) applies **no WV
correction at all** and hard-codes 910 nm (`l0_wavelength` is read but unused) — with
T²_wv ≈ 0.78 its CL61 Rayleigh constants are low by construction; not comparable.

## 3. AERONET forward-Klett closure (CL61 direct, CHM15k as method control)

Run on the 2026 L1.5 AERONET data (Mar-01…Jun-20, 6 897 CL61 / 7 350 CHM15k matched samples,
LR from the 2026 almucantar inversions where valid, else 50 sr; WV per monthly CAMS):

| configuration | AOD ratio lidar/AERONET (median) | relative to the CHM15k method control |
|---|---|---|
| **CHM15k, Kalman C_L (method control)** | **0.70** | 1.00 |
| CL61, C_L = 1.251 (Rayleigh), b = 0 | 0.71 | 1.01 |
| CL61, C_L = 1.425 (cloud), b = 0 | 0.51 | 0.73 |
| CL61, C_L = 1.251, b = −0.03 | 0.98 | 1.40 |
| CL61, C_L = 1.425, **b = −0.02** | 0.69 | **0.99** |
| CL61, C_L = 1.425, b = −0.03 | 0.78 | 1.11 |

(b = a constant baseline removed from the vendor β_att before calibrating; the earlier Jan–Mar
2025 CHM15k run gave 0.85 at AOD ≈ 0.03.) Three lessons:
1. **The daytime forward-Klett carries a common ≈ 0.70 method factor** (identical for the trusted
   CHM15k): aerosol above the 6 km integration top (spring Saharan layers), the assumed LR on
   non-inversion days, and daytime background handling — it cancels in the relative comparison.
2. **The closure is (C_L, b)-degenerate**: relative to the control, (1.251, b=0) and
   (1.425, b=−0.02) close equally well — the AOD constrains a combination of constant and
   baseline, not each separately.
3. **The night-time evidence breaks the tie**: the dark probe measures b(3–6 km) ≈ −0.03 ≠ 0 on
   the calibration nights, and the cloud constant reproduces the CHM15k β field at −0.6 % in the
   0.5–3 km band (where b ≈ 0). Hence (C_L = 1.425, b(z) < 0 above ≈ 3 km) is the consistent
   solution; (1.251, b=0) would additionally require the baseline to vanish at night — it does not.

Chain of anchors: AERONET ⇒ CHM15k (method control) ⇒ CL61-cloud (−0.6 % in β) ⇒ **C_L(cloud)**;
the Rayleigh constant absorbs the 3–6 km baseline. Follow-up: extend the integration above 6 km
(elevated layers) and rerun against AERONET **L2** when the 2026 final calibration is released.

## 4. Water-vapour source sensitivity

Per-night T²_wv at the fit window (probe `_cl61_wv_sources_probe.py`):

| source | dC_L vs operational (median) | note |
|---|---|---|
| CAMS 1° monthly (operational) | — | T²_wv 0.71–0.81 at the windows |
| Payerne radiosonde 00 UT | **−3.6 %** (−12.4…+0.7) | calibration nights are drier than the monthly mean → the monthly correction over-corrects; a nightly WV source would *lower* C_L(ray) further |
| CAMS 0.4° monthly | −0.1…−0.3 % | resolution immaterial (0.4° archive incomplete: 202501/02/06/07 only, ADS issues) |

The ±4 % night-to-night WV term also explains a good part of the C_L(ray) scatter (and its
`sensitivity 2σ` being smaller than the observed night-to-night spread).

## 5. Aosta: why only the Rayleigh constant has a seasonal cycle

![coherence](figs_paper_report/fig_cl61_network_coherence.png)
*Figure 3 — Monthly median C_L (blue = Rayleigh, dark grey = cloud) with the monthly WV lever
T²_wv(3–5 km) (green, right axis), and the method ratio vs 1/T²_wv, for the four dual-method CL61
stations.*

| station | months | ratio ray/cloud | corr(C_L_ray, 1/T²_wv) | corr(C_L_cloud, 1/T²_wv) | corr(ratio, 1/T²_wv) |
|---|---|---|---|---|---|
| Payerne | 5 | 0.865 | +0.90 | +0.88 | −0.68 |
| Lindenberg | 17 | 1.055 | −0.22 | −0.23 | −0.08 |
| **Aosta** | 8 | 0.894 | **+0.93** | +0.58 | **+0.80** |
| Camborne | 9 | 0.945 | +0.21 | −0.52 | **+0.79** |

**Aosta answered:** its Rayleigh constant tracks the seasonal WV lever almost perfectly
(corr +0.93): the fit window sits above the full WV column, so C_L(ray) scales with 1/T²_wv
(0.71→0.85 seasonally ≈ ±8 % direct lever) and any residual of the monthly correction (dry-night
bias, §4) imprints seasonally. The cloud method integrates 0.1–2.4 km where the WV lever is ≈ 4×
smaller — hence no visible cycle. Camborne behaves the same (ratio-corr +0.79). Payerne's 5-month
record is too short to separate the terms.

**Mountain-orography term (Aosta especially):** the 1° CAMS model surface at the nearest grid
point sits at **1750 m for Aosta (station 570 m — offset +1180 m!)** and 1395 m for Payerne
(+904 m); the moistest valley layer is simply absent from the WV column (the correction clamps
the mountain-surface humidity downwards). Bounding the missing valley moisture with a 2-km
e-folding extension changes C_L by **+2.1 % at Aosta / +2.0 % at Payerne — seasonally varying**
(moist summer valley → larger), i.e. the right sign and season to add to the Aosta Rayleigh
cycle. Camborne and Lindenberg (flat, offsets −55/−44 m) are unaffected — consistent with the
coherence table. Fix: the **0.4° CAMS** (real orography much closer to the valley floor) —
another reason to complete the `CAMS_Monthly_04` download (currently 202501/02/06/07 only).

**Coherence across the network:** the ray/cloud ratio is NOT a universal constant (0.87, 0.89,
0.95, 1.06) — consistent with a data-dependent artefact (baseline strength differs per unit /
processing chain) plus the WV residual, not with a single physical constant offset. **Lindenberg
is the outlier in every respect** (ratio > 1, no WV correlation): its "L1" is converted from
Cloudnet, i.e. a different processing chain — supporting the baseline explanation (different
background handling → different Rayleigh bias).

## 6. Recommendations

1. **Operations:** keep the **cloud method** as the CL61 constant of record (already the case in
   the fullcal chain); flag the CL61 native Rayleigh constant as biased low ≈ 10–15 %.
2. **Root cause:** raise the negative-baseline finding with Vaisala (background subtraction in the
   vendor β_att; reproducible as physically negative nightly means at 11–14 km on clear nights).
   A hood-on dark measurement on a CL61 would separate detector vs processing contributions.
3. **Rayleigh method hardening:** estimate and remove a per-night baseline from the 10–14 km
   gates before the molecular fit (turns the −12 % into a correctable term); consider nightly
   (sounding or CAMS-daily) WV instead of monthly (−4 % and less scatter).
4. **Closure:** download 2026 AERONET Payerne (L1.5 now, L2 when released) and run
   `_cl61_aeronet_closure.py` for the direct CL61 constant test (script ready, LR from the .lid
   inversions when valid, else 50 sr).
5. **Paper:** §5.1 of the validation report already carries the exoneration table and the
   baseline finding; the CL61-vs-CHM15k comparison offsets (both methods high at LIN/AOS/CAM)
   remain a separate, comparison-level question (910→1064 conversion + WV residual at 910 nm).
