# E-PROFILE ALC attenuated-backscatter validation — methodology and results

*Consolidated 2026-07-10 (M. Hervo, MeteoSwiss E-PROFILE ALC). Merges: attbsc_validation_technical.md, paper_validation_L1_report.md, paper_validation_report.md, paper_python_validation.md, l1_vs_l2_validation.md, l1_validation_cscs.md, payerne_noise_filter_comparison.md.*

This report is the durable reference for the validation of the E-PROFILE ALC calibrated
attenuated-backscatter product β_att. The **authoritative results** are those of the
**uniform Level-1 rewrite (`run_paper_validation.py`, generated 2026-07-08, extended to
30 June 2026)**, which folds in the 2026-07 cloud recalibration (per-Vaisala-type multiple-
scattering tables), the molecular-aware 910→1064 nm wavelength conversion, and the removal of
the SNR selection filter. Earlier result tables (the 2026-07-02 L2-based paper validation and
the operational-Python-vs-MATLAB benchmark) are retained only as one-line historical notes where
their numbers were superseded — they are **not** re-tabulated here.

Two companion investigations are cross-linked but not reproduced in depth:
[report 07](07_cl61_calibration.md) (CL61 Rayleigh-vs-cloud and the
hood-dark campaign; consolidated as report 07) and the near-range tilt / electronic-offset analysis
(consolidated as report 08). The **Payerne CL61–CHM15k near-range tilt (200–900 m) is RESOLVED**
(2026-07-08/09): it is a unit-specific near-range (overlap-normalization) difference, dominated by
the CHM15k TUB140016 static-overlap module — the water-vapour, molecular T/p, Ångström-centre and
electronics corrections are each quantified and **exonerated**. Any statement below that a
comparison is "confounded by WV / Ångström / T-p / electronics" reflects an earlier reading and is
corrected accordingly.

## Table of contents

1. [Methodology and conventions](#1-methodology-and-conventions)
2. [Paper validation — station intercomparison and EARLINET (uniform L1)](#2-paper-validation--station-intercomparison-and-earlinet-uniform-l1)
3. [Where the significant differences come from](#3-where-the-significant-differences-come-from)
4. [Water-vapour, electronic-offset and dark corrections in the validation](#4-water-vapour-electronic-offset-and-dark-corrections-in-the-validation)
5. [L1 vs L2 calibration constant](#5-l1-vs-l2-calibration-constant)
6. [L1 β_att validation (CSCS calibration applied to the raw signal)](#6-l1-att-validation-cscs-calibration-applied-to-the-raw-signal)
7. [Noise-filter sensitivity — keep `none`](#7-noise-filter-sensitivity--keep-none)
8. [Limitations, caveats and reproducibility](#8-limitations-caveats-and-reproducibility)

---

## 1. Methodology and conventions

### 1.1 Objective and scope

We validate the **calibrated attenuated-backscatter coefficient** β_att produced by the E-PROFILE
automatic lidars and ceilometers (ALC). All comparisons are of the *calibrated* product, so they
jointly test (i) the absolute calibration, (ii) the water-vapour and wavelength harmonisation, and
(iii) instrument hardware differences. Two independent validation families are run through **one
single methodology for every instrument at every site**:

1. **Multi-ceilometer station intercomparison** — co-located instruments at 7 stations, each
   calibrated independently (nightly Rayleigh `eprof_v2` or daily liquid-cloud O'Connor, both
   Kalman-smoothed), compared pairwise against the station reference over 500–3000 m AGL.
2. **Ceilometer vs EARLINET** — the CHM15k at Leipzig, Palaiseau and Magurele (both co-located
   units A and B at Magurele) against the EARLINET research-lidar reference over 500–5000 m AGL;
   the Trappes Mini-MPL against the SIRTA 532 nm channel at its native wavelength.

Everything is computed from the **native Level-1 `rcs_0`** — no provider L2 β_att, no provider
calibration constant, no MATLAB legacy reference. This removes the per-type unit conventions (the
CL61 `calibration_constant_0` exception of earlier drafts disappears by construction) and makes
the three physical corrections explicit and identical everywhere: overlap (CHM15k), water vapour
(910 nm family), and the molecular-aware 910→1064 nm conversion.

**Sites and channels.** Payerne (CHM15k reference; CL31 native + hood-dark corrected, CL61 ×2
methods), Amsterdam/Schiphol (4 × CHM15k), Uccle (CL51 reference native + offset-corrected; CL61),
Palaiseau/SIRTA (CHM15k reference; CL31, Mini-MPL 532 nm), Lindenberg, Aosta, Camborne (CHM15k
reference; CL61 ×2 methods each). EARLINET systems: `ari`/`lei` (Leipzig), `sir`/`sir_532`
(Palaiseau), `cbw` (Cabauw), `ino`/`ino_a` (Magurele).

### 1.2 Data

| Source | Content | Period | Notes |
|---|---|---|---|
| E-PROFILE L1 daily archive (`D:/E-PROFILE_L1_2026`) | native `rcs_0`, CBH, vertical visibility, internal temperature `temp_int` | per site config | the **only** β input |
| Calibration series (`calib/<key>_L1.csv`) | nightly Rayleigh + daily cloud lidar constants C_L, Kalman-smoothed | 2025-01 → 2026-06 | all L1-derived; same series feeding the monitoring dashboard |
| Overlap models (`ALC_OVERLAP_DIR` = `D:/TEMP_MODELS/202606`) | per-instrument `a(z)`, `b(z)`, `overlap_ref(z)` (`overlap_probe_eprofile`) | one model per unit, fitted on 53–778 days spanning 2020–2026 | **127 units — every CHM15k in this study covered** (Amsterdam ×4 and Magurele ×2 built 2026-07-03) |
| Electronic-offset patterns (`network_offset/<key>.npz`) | fixed digitizer-ripple `b_phys(z)` from clear nights | Mar–Jun 2026 | applied at Uccle (CL51), §4.2 |
| Terminal-hood dark (`cl31_b_dark.npz`) | measured covered-telescope background `P_dark(z)·z²`, two-resonance model | 4 hood sessions May–Jun 2026 | applied at Payerne (CL31), §4.3 |
| CAMS reanalysis (`D:/CAMS`, `D:/CAMS_Monthly_04`) | water-vapour + T/p profiles, monthly **1°** (June 2026: daily **0.4°**) | 2018-01 → 2026-06 | fail-safe: uncorrectable month ⇒ NaN, never fail-open |
| EARLINET SCC L2 (`A:/EARLINET`) | 1064 nm (532 nm for `sir_532`) particle backscatter + per-scene lidar ratio | 2025-01 → 2026-06 | cloud-screened upstream by the SCC |
| US standard atmosphere 1976 | molecular backscatter/extinction | — | fallback molecular where CAMS absent |

### 1.3 Calibration application and the Wiegner C_L convention

All calibration series are expressed as the **absolute Wiegner lidar constant `C_L = RCS/β_att`** —
the single physical constant of the instrument, for both the Rayleigh and the liquid-cloud methods,
and the same quantity the monitoring dashboard displays. This is the single naming/definition
convention across code, outputs, figures and reports; the cloud O'Connor coefficient `C`
(`β_true = C·β_L2`) maps to it as `C_L = calibration_constant_0 / C`.

In the uniform-L1 pipeline one formula covers both methods, applied directly on the `rcs_0` scale
(no provider constant to divide out):

```
β_att [Mm⁻¹ sr⁻¹] = rcs_0 / C_L · 1e6,   C_L the daily-Kalman constant interpolated to profile times.
```

The daily C_L is produced per night (Rayleigh, `eprof_v2`) or per day (liquid-cloud, O'Connor/Hopkin)
from the native L1 archive and smoothed with the operational E-PROFILE Kalman filter (predict/update,
process-noise floor, daily grid). Outliers are rejected before smoothing
(|log C_L − median(log C_L)| > 4·MAD — the raw series occasionally contains physically impossible
C_L that would wreck the filter); a short/gappy record the Kalman cannot handle falls back to the
constant median C_L. Because the cloud method is offset-immune, the offset-corrected twin keeps its
native C_L.

> **Superseded (earlier L2-based drafts):** the L2 comparison applied the constant as
> `β = β_L2·(calibration_constant_0 / C_L)` (Rayleigh) or `β = β_L2·(INSTRUMENT_CAL_DEFAULT / C_L)`
> (cloud), and carried a **CL61 exception** — the CL61's L2 `calibration_constant_0` is Vaisala's
> internal factor in a different unit system and had to be left un-divided (applied default = 1).
> The uniform-L1 path derives β from `rcs_0` directly, so this per-type convention exception
> **disappears by construction** and is no longer a caveat.

### 1.4 The uniform per-channel pipeline (eight steps)

Every channel at every site goes through the same steps
(`run_paper_validation.channel_beta` + `grid_and_stats`):

1. **Read L1** `rcs_0` from the daily archive (`intercompare.read_l1`), hourly-median retimed, with
   CBH / vertical visibility / internal temperature (K→°C) carried along.
2. **Overlap correction (CHM15k)** — `rcs_0 · (1 + (a(z)·T + b(z))/100)`, T = internal temperature
   in °C (`overlap.correct_rcs`). The Level-1 already carries the *static* reference overlap; the
   model corrects the *temperature-dependent residual* (Hervo et al. 2016, `overlap_probe_eprofile`).
   The convention was verified against the package source: the model fits
   `Dif = 100·(O_ref − O_daily)/O_daily` against T in °C, so the signal factor is exactly
   `1 + Dif/100` (the reference overlap cancels). The correction acts below ~720 m (full overlap) —
   in-band (≥ 500 m) it moves the statistics by < 0.1 % (verified by an on/off rerun), so it matters
   for the *near-range product*, not the numbers below. Each station figure shows the check directly:
   the CHM15k median profile is drawn corrected (solid) and uncorrected (dashed), separating only
   below ~700 m.
2b. **Electronic-offset correction (flagged units)** — `rcs_0 − b_phys(z)`, the clear-night
   digitizer-ripple pattern. Applied at Uccle (CL51, the network's strongest case), shown as an
   *additional channel* so native and corrected are both in the results (§4.2).
3. **Calibration** — the C_L formula above (both methods).
4. **Water vapour (CL31/CL51/CL61)** — divide by the two-way WV transmission from monthly CAMS at
   the instrument's laser line (`intercompare.apply_wv`). A month with no usable CAMS is NaN-masked,
   never passed through uncorrected. Impact quantified per instrument in §4.1.
5. **Molecular-aware 910→1064 nm conversion** (CL31/CL51/CL61 → 1064 nm): the analytic molecular
   (Rayleigh) attenuated backscatter is removed and re-added at the target line, and the Ångström
   law (α = 1) scales the **aerosol residual only**:
   `β₁₀₆₄ = β_mol·T²_mol|₁₀₆₄ + [β₉₁₀ − β_mol·T²_mol|₉₁₀]·(910/1064)^α`. The molecular part uses
   `β_mol ∝ p/T` from **CAMS T/p** (`intercompare.wavelength_correct_molecular`; 0.4° monthly with a
   1° fallback, same archive as the WV step; hydrostatic extrapolation below the lowest CAMS level
   for elevated/valley sites such as Aosta; US-standard fallback if a month has no CAMS). This
   replaces a single-Ångström-on-the-total step that mis-scaled the λ⁻⁴ molecular part (§3.4). The
   532 nm Mini-MPL keeps its two-component `molaer` form (US-std); the Uccle CL51↔CL61 pair is
   910-vs-910 (no conversion). *(Verified in code: `run_paper_validation.py` calls
   `wavelength_correct_molecular` for the 910 nm family and `wavelength_correct(..., 'molaer')` for
   the Mini-MPL.)*
6. **Screening (science stream)** — quality flag > 0; any cloud base 0–20 km; fog / finite vertical
   visibility; ±15 min expansion. The display stream keeps clouds visible (quality-flag masking only).
   EARLINET profiles are cloud-screened upstream by the SCC.
7. **Gridding** — native profiles (~30 s) median-aggregated onto a 60-min grid, bins kept only with
   ≥ 30 min coverage; channels aligned on the union time grid and a common altitude grid. **No SNR
   filter** on the compared samples (§7): the per-gate window SNR (σ_rob = 1.4826·MAD) is still
   computed but used *diagnostically* (CL31 detection-limit profile, §3.1). `ALC_VAL_L1_SNR=1`
   restores the legacy SNR filter for cross-checks.
8. **Statistics vs the site reference** over 500–3000 m AGL (EARLINET 500–5000 m): headline pair
   **median relative bias** and **log-space Pearson r**; linear relbias and r kept for continuity.

The EARLINET path (`earlinet.compare`) uses steps 1–3 (+ overlap for CHM15k), no WV/wavelength
conversion (native-wavelength comparisons), the same screening, ±30-min matching with ≥ 30-min
averages on both sides, and a physical EARLINET β_att reconstruction (per-scene lidar ratio,
below-overlap OD extension, gates below the instrument overlap excluded).

### 1.5 Metric definitions

Over all valid (time, altitude) pairs in the band (stations 500–3000 m AGL; EARLINET 500–5000 m):

| metric | definition | character |
|---|---|---|
| relbias | `100·mean(x−r)/mean(r)` | dominated by rare large values (clouds/plumes) |
| **med relbias** | `100·median((x−r)/r)`, r > 0 | robust central agreement |
| r | Pearson on **linear** β | dominated by large values |
| **log r** | Pearson on **log₁₀ β**, positive pairs | robust across the dynamic range |
| N | number of valid pairs | — |

The robust pair **(med relbias, log r)** is the headline; the linear pair is kept for continuity
and for sensitivity to the aerosol-event scaling. Panel (a) of every station figure shows medians
restricted to the **common hours** where every channel reports, so the profiles describe the same
atmospheric sample. On the linear axis the median is drawn into negative values — a channel whose
median crosses β = 0 has hit its noise/offset floor (that is the information); only a > 20 %-coverage
rule ever cuts a curve, so it is not an SNR filter.

### 1.6 Harmonisation corrections — the physics

**Water vapour (910 nm only).** CL31/CL51/CL61 emit in the 905–911 nm H₂O absorption band, so their
β_L2 is attenuated by water vapour that does **not** affect the 1064 nm CHM15k. The spectrally-
averaged two-way transmission T²_wv(z) is computed from CAMS humidity and a HITRAN/MT-CKD LUT,
weighted by each instrument's measured laser spectrum, and β is divided by it. The WV correction is
**mandatory** at 910 nm (a no-WV degraded mode is rejected — it worsens results); a 910 nm night
without usable CAMS is flagged, never calibrated WV-free. Per-instrument laser parameters
(2026-06-02 Qmini campaign):

| Instrument | λ₀ [nm] | FWHM [nm] | WV correction |
|---|---|---|---|
| CL31 | 909.7 | 6.0 | yes |
| CL51 | 910.0 | 3.4 | yes |
| CL61 | 910.74 | 1.0 | yes |
| CHM15k | 1064.47 | 0.5 | no (outside band) |

Typical effect: T²_wv ≈ 0.80–0.87, i.e. β boosted ~14–23 % (quantified per site in §4.1).

**Fail-safe, never fail-open.** A month with no CAMS file, or whose nearest grid point is farther
than `max(1°, 1.5 grid cells)` from the station (`cams_point_too_far`, the same guard that gives
flag −10 in the calibration), is **excluded (NaN) and reported** — data is never passed through
uncorrected or corrected with a domain-edge point. The 2025–2026 runs excluded **no** months (all
11 stations well inside the Europe box, nearest grid point ≤ 0.7°).

**Wavelength.** The molecular-aware component-separated conversion (step 5) supersedes the single-
Ångström-on-the-total step; the aerosol Ångström exponent is α = 1 (continental backscatter
Ångström 0.9–1.2; DeLiAn / Floutsi 2023, Haarig 2025) and the result is weakly sensitive to it.

---

## 2. Paper validation — station intercomparison and EARLINET (uniform L1)

*Authoritative results: `run_paper_validation.py`, generated 2026-07-08, six recent stations
extended to 30 June 2026, uniform Level-1 methodology, 2026-07 cloud recalibration folded in, no
SNR filter. Machine-readable numbers in `figs_l1_validation/summary_stats.csv` (28 rows) and
`discrepancy_analysis.json`.*

**Headline.** The CHM15k network is mutually consistent to a few % (Amsterdam C −4.8 % / D +1.4 %
vs A, both Magurele units vs EARLINET +1/+3 %, Leipzig +4 %, Palaiseau −8 %). After the 2026-07
cloud recalibration the **two independent CL61 calibrations converge** — liquid-cloud and Rayleigh
give the same lidar constant (Payerne C_L 1.207 vs 1.216, < 1 %). With the molecular-aware
910→1064 nm conversion in the pipeline, the 910 nm **CL61 agrees with the 1064 nm CHM15k to within
±3 %** at Payerne (+0.1 %), Lindenberg (−1.7 %) and Camborne (−2.9 %) — down from the
+21 / +32 / +28 % the retired single-Ångström step produced (§3.4). Two residuals remain, both
localised and *not* methodological: **Aosta** CL61 +12.0 % (a degraded 82 % window, §3.4) and
**Uccle** CL61 −16.5 % (a 910-vs-910 CL61-vs-CL51 calibration difference, §3.5).

### 2.1 Analysis flow

![flow intercompare](figs_paper_report/fig_flow_intercompare.png)
*Figure 1 — Station-intercomparison flow. Blue: inputs. Orange: calibration series. Purple:
auxiliary atmospheric data. Red: the screening / averaging / detection gates (science stream).
Green: outputs. The median-profile panel uses only hours where every channel reports; the
statistics are pairwise vs the reference channel.*

![flow earlinet](figs_paper_report/fig_flow_earlinet.png)
*Figure 2 — Ceilometer-vs-EARLINET flow. The EARLINET branch (left) converts SCC particle
backscatter to attenuated backscatter with the per-scene lidar ratio, a molecular profile from the
standard atmosphere, and a two-way transmission whose below-overlap extinction is extended from the
lowest trusted gate (backscatter itself stays NaN there). The CHM branch (right) is calibrated and
screened exactly like the station intercomparison, then matched within ±30 min.*

### 2.2 Calibration series

![calibration time series](figs_l1_validation/fig_calib_timeseries.png)
*Figure 3 — Lidar constant C_L for every instrument (raw nightly/daily ×, Kalman line ± 1σ;
blue = Rayleigh, dark grey = cloud; CL61 panels overlay both methods — after the 2026-07 cloud
recalibration the two now overlap instead of splitting, §3.3). Input to everything that follows.*

Stability of the 910 nm cloud series (`discrepancy_analysis.json`):

| series | raw days | daily scatter | dominant period | seasonal amplitude |
|---|---|---|---|---|
| Uccle CL51 (cloud) | 263 | 13 % | **272 d** | **34 %** |
| Payerne CL31 (cloud) | 259 | 52 % | 136 d | 87 % |
| Payerne CL61 (cloud) | 62 | 9 % | 61 d | 12 % |

The single-photodiode 910 nm Vaisalas oscillate seasonally even with the WV correction ON —
consistent with a laser centre-wavelength drift across the steep 910 nm absorption band that the
fixed-line WV model cannot cancel. Any result for a 910 nm cloud-calibrated channel therefore
depends on the season sampled (§3.3, §3.5).

### 2.3 Station intercomparison (500–3000 m AGL, screened, no SNR filter)

Reference channels (0 by construction) omitted; **med relbias / log r** is the headline pair.
WV = in-band β increase from the water-vapour correction (§4.1). The 910 nm channels use the
molecular-aware CAMS conversion (§1.4 step 5); the CHM15k references and the 910-vs-910 Uccle pair
are unaffected by it.

| station | channel | med relbias | log r | relbias | r | N | WV |
|---|---|---|---|---|---|---|---|
| Payerne | CL31 (cloud) | −7.9 % | 0.18 | +21.0 % | 0.38 | 199 231 | +23 % |
| Payerne | CL31 (cloud, hood-dark corr) | −19.0 % | **0.37** | −39.7 % | 0.52 | 199 231 | +23 % |
| Payerne | CL61 (cloud) | **+0.1 %** | **0.94** | −2.0 % | 0.98 | 157 815 | +30 % |
| Payerne | CL61 (Rayleigh) | **−0.6 %** | 0.94 | −0.9 % | 0.99 | 157 815 | +30 % |
| Amsterdam | CHM15k B | +17.3 % | 0.97 | +25.8 % | 0.95 | 195 250 | — |
| Amsterdam | CHM15k C | **−4.8 %** | 0.95 | −5.0 % | 0.92 | 186 250 | — |
| Amsterdam | CHM15k D | **+1.4 %** | 0.95 | +3.3 % | 0.95 | 196 500 | — |
| Uccle | CL61 (cloud) | −16.5 % | 0.69 | −20.3 % | 0.86 | 135 791 | +21 % |
| Uccle | **CL51 (cloud, offset-corr)** | **−0.5 %** | 0.76 | −0.2 % | 0.92 | 171 182 | +20 % |
| Palaiseau | CL31 (cloud) | −40.1 % | 0.69 | −41.2 % | 0.85 | 153 716 | +27 % |
| Palaiseau | Mini-MPL (532→1064) | −52.0 % | 0.80 | −48.1 % | 0.82 | 68 807 | — |
| Lindenberg | CL61 (cloud) | **−1.7 %** | 0.93 | −2.3 % | 0.98 | 862 000 | +29 % |
| Lindenberg | CL61 (Rayleigh) | **+2.5 %** | 0.93 | +0.2 % | 0.98 | 862 000 | +29 % |
| Aosta | CL61 (cloud) | +12.0 % | 0.95 | +8.2 % | 0.96 | 266 699 | +21 % |
| Aosta | CL61 (Rayleigh) | −5.0 % | 0.94 | −10.9 % | 0.95 | 266 699 | +23 % |
| Camborne | CL61 (cloud) | **−2.9 %** | **0.97** | −1.1 % | 0.98 | 82 999 | +35 % |
| Camborne | CL61 (Rayleigh) | **+1.5 %** | 0.97 | +1.1 % | 0.99 | 82 999 | +34 % |

**Median without the SNR filter.** The table is the **unfiltered** median over all ≥ 30-min-covered
gates (earlier drafts kept only per-gate SNR ≥ 3 samples). Two things follow, and both are the point
of the paper (see §7 for the full sensitivity):

- **The CL61 is unchanged by removing the SNR filter** (Payerne −0.5 → +0.1 %, Lindenberg
  −1.3 → −1.7 %, Aosta +11.6 → +12.0 %, Camborne −3.0 → −2.9 %; log r drops only ~0.03). Its ±3 %
  agreement is **not** an artefact of the SNR selection — it holds on the raw signal.
- **The CL31 median moves a lot** (Payerne −4.6 → −7.9 %, Palaiseau −13.4 → **−40.1 %**), because
  without the SNR filter its gates above the detection ceiling — where the CL31 reads noise but the
  CHM15k still reads real molecular signal — pull the ratio toward −100 %. That is the detection
  ceiling showing up in a column statistic, which is why the CL31 is reported as a **detection-limit
  profile** (§3.1), not a single median.

Uccle CL61 (Rayleigh) has **no usable L1 calibration series** — the native-signal molecular fit
converged on only **4 nights** over the whole 2025-12 → 2026-06 span (versus **200 days** for the
cloud method on the same unit), too few and too scattered to Kalman-smooth, so it is skipped
(reported, not hidden). The cause is **data delivery, not the sky**: this Uccle stream arrives in a
~50 % duty cycle (30 s bursts separated by regular ~5.5 min transmission gaps), so the Rayleigh
gate's ≥ 3 h-of-clear-profiles requirement is almost never met on short high-latitude summer nights
(a network-side collection issue, not correctable here). The Uccle CL61 (cloud) reads −16.5 %
against the CL51 reference (both 910 nm, so the wavelength conversion is a no-op) — a genuine
CL61-vs-CL51 calibration difference (§3.5), not a spectral effect.

> **Superseded — 2026-07-02 L2 paper validation** (`figs_paper_report/`; not re-tabulated): the same
> stations before the 2026-07 cloud recalibration and molecular-aware conversion. It read the CL61
> as a ~15 % method discrepancy (Payerne cloud +4.3 % vs Rayleigh +20.9 %) and the network CL61 far
> high (Lindenberg +35/+43 %, Aosta +11/+35 %, Camborne +39/+53 % cloud/Rayleigh) under a single
> Ångström conversion. Both readings are corrected above.
>
> **Superseded — operational-Python-vs-MATLAB benchmark** (`paper_python_validation.md`): the same
> Mar–May 2026 tables carrying a MATLAB legacy column (e.g. Payerne CL61 cloud +16.0 %, Rayleigh
> +1.7 % under the older WV/coordinate fixes). Kept only for provenance; the numbers are those of
> its own earlier calibration snapshot, not the current pipeline.

### 2.4 Ceilometer / lidar vs EARLINET (500–5000 m AGL)

Unfiltered, the same as the station intercomparison: the per-match SNR ≥ 3 removal is off (the
EARLINET profiles are already cloud-screened by the SCC, and we compare 30-min averages on both
sides).

| site | compared instrument | med relbias | log r | relbias | r | matched |
|---|---|---|---|---|---|---|
| Palaiseau `sir` | CHM15k (Rayleigh) | **−8.4 %** | 0.75 | −6.2 % | 0.91 | 196 |
| Magurele `ino` | CHM15k **B** (Rayleigh) | **+1.1 %** | **0.86** | −1.8 % | 0.42 | 636 |
| Magurele `ino_a` | CHM15k **A** (Rayleigh) | **+2.8 %** | **0.90** | −0.6 % | 0.41 | 636 |
| Leipzig `ari` | CHM15k (Rayleigh) | **+4.4 %** | **0.92** | +10.0 % | 0.91 | 1151 |
| Palaiseau `sir_532` | Mini-MPL (native 532 nm) | **−3.0 %** | 0.86 | −4.7 % | 0.87 | 37 |

Leipzig `lei` and Cabauw `cbw` have no EARLINET 1064 nm files in 2025–2026. The two co-located
**Magurele units agree to 1.7 %** (+1.1 / +2.8 %) against the same independent reference, each with
its own nightly calibration and its own overlap model — a direct unit-consistency validation of the
whole chain. At native 532 nm the Mini-MPL agrees with EARLINET to **−3.0 %** — the instrument and
its Rayleigh calibration are fine; its −52 % *station* entry (§2.3) is the 532→1064 conversion (§3.2).
Profile-level linear r is low by construction (different technique, ±30 min sampling, per-scene
EARLINET lidar ratio) — the median agreement and log r are the meaningful metrics; Magurele's
linear r 0.42 coexisting with log r 0.86 is the clearest argument for the log-space metric.

> **Superseded** (2026-07-02 L2 path): Palaiseau `sir` −4.4 % / log r 0.91 (177 matched), Magurele
> `ino` +2.1 % / 0.92 (625), Leipzig `ari` +5.8 % / 0.94 (933), `sir_532` −1.9 % / 0.88 (54). These
> carried the SNR ≥ 3 gate; the unfiltered uniform-L1 numbers above are the current ones (removing
> the filter shifts Magurele/Leipzig toward zero and moves the noisiest `sir` set, which reaches to
> 5 km, the other way — see §2.5).

### 2.5 Station and EARLINET figures

Layout per station: (a) median ± IQR profiles over common hours; (b) scatter vs the reference
(500–3000 m, log-log); (c) histogram of differences; (d–g) curtains, all data in greyscale with only
the kept (science-stream) gates in colour; black dots = cloud base. Channel colours: CHM15k red,
CL31 orange, CL51 purple, Mini-MPL green, CL61 blue (Rayleigh) / dark grey (cloud); Amsterdam's four
CHM15k units use a distinct palette.

![payerne](figs_l1_validation/fig_payerne.png)
*Figure 4 — Payerne, Mar–Jun 2026. CHM15k reference; in panel (a) the CHM15k median is drawn
overlap-corrected (solid red) and uncorrected (dashed red) — the curves separate only below ~700 m,
the direct check of the temperature-dependent overlap correction. The two CL61 entries are the same
instrument calibrated two ways: cloud **+0.1 %**, Rayleigh **−0.6 %** — agreeing with each other
(§3.3) and with the CHM15k under the molecular-aware 910→1064 nm conversion (§3.4), unchanged when
the SNR filter is dropped. The CL31 unfiltered median is −7.9 % with the network's lowest log r
(0.18); on the linear axis its median plunges through zero to ≈ −4 Mm⁻¹sr⁻¹ above ~4 km — the
detection ceiling made visible (§3.1). Its terminal-hood dark-corrected twin (light blue) subtracts
the measured covered-telescope background (§4.3): in 0.5–2 km it drops onto the CHM15k, log r rises
0.18 → 0.37, but the median moves the other way (−7.9 → −19.0 %) — that native "agreement" was
partly the positive dark pedestal inflating β.*

![amsterdam](figs_l1_validation/fig_amsterdam.png)
*Figure 5 — Amsterdam, four co-located CHM15k, unit A reference, all overlap-corrected with
unit-specific models. C −4.8 %, D +1.4 % vs A (unfiltered); unit B +17.3 % — a real unit effect,
§3.6.*

![amsterdam overlap impact](figs_l1_validation/fig_amsterdam_overlap_impact.png)
*Figure 5b — Impact of the temperature-dependent overlap correction on the four Amsterdam units (the
only site with four unit-specific models side by side). (a) The four model factors at each unit's
observed median internal temperature: genuinely unit-specific in shape AND sign below ~400 m.
(b) Median β_att corrected (solid) vs uncorrected (dashed). (c) Realised impact on the science
stream: ±2–3 % at 350 m, < 1 % at 500 m, tens of % below ~250 m. (d) Four-way agreement in the
300–700 m near-range band, with vs without the correction: the healthy units tighten slightly while
unit B's near-range excess is untouched (+34.8 → +35.1 %) — independent confirmation that B's problem
is not overlap (§3.6).*

**Day-only / night-only views** (solar elevation > 5°) of the identical synchronized matrices:

![amsterdam day](figs_l1_validation/fig_amsterdam_day.png)
*Figure 5c — Amsterdam, DAY only. Unit B's median profile separates from A/C/D through the lowest
~1.5 km; the B difference histogram is visibly displaced positive.*

![amsterdam night](figs_l1_validation/fig_amsterdam_night.png)
*Figure 5d — Amsterdam, NIGHT only. The four median profiles nearly collapse; B retains a reduced
+11 % floor.*

| channel | hours | med relbias | log r | relbias | r | N |
|---|---|---|---|---|---|---|
| CHM15k B | all | +17.3 % | 0.97 | +25.8 % | 0.95 | 195 250 |
| CHM15k B | **day** | **+25.2 %** | 0.96 | +33.0 % | 0.94 | 102 750 |
| CHM15k B | **night** | **+10.8 %** | 0.98 | +16.6 % | 0.96 | 92 500 |
| CHM15k C | all | −4.8 % | 0.95 | −5.0 % | 0.92 | 190 000 |
| CHM15k C | day | −7.4 % | 0.94 | −8.5 % | 0.89 | 101 250 |
| CHM15k C | night | −2.0 % | 0.97 | −0.5 % | 0.96 | 88 750 |
| CHM15k D | all | +1.4 % | 0.95 | +3.3 % | 0.95 | 196 500 |
| CHM15k D | day | +1.3 % | 0.94 | +2.4 % | 0.95 | 103 500 |
| CHM15k D | night | +1.5 % | 0.97 | +4.5 % | 0.96 | 93 000 |

Reading (balanced sampling, ≈ 103 k day / 92 k night pairs): unit B swings by **+14 points** between
night (+10.8 %) and day (+25.2 %) while its log r stays 0.96–0.98 in both — the daytime excess is a
*bias*, not added scatter. C shifts negative by day (−7.4 vs −2.0 % at night) while D stays flat
(+1.3 / +1.5 %). Practical network-QC consequence: a **day-minus-night split of the median relative
bias** is a cheap, calibration-independent detector of this failure mode — B's Δ(day−night) = +14
points stands out against ≤ 5 points for the healthy units.

![amsterdam without B](figs_l1_validation/fig_amsterdam_noB.png)
*Figure 5e — Amsterdam, healthy trio only (A reference, C, D). The three median profiles collapse
over the full column (601 common hours); C −4.8 % / log r 0.95, D +1.4 % / log r 0.95 — the CHM15k
unit-to-unit consistency in its cleanest form.*

![uccle](figs_l1_validation/fig_uccle.png)
*Figure 6 — Uccle, Mar–Jun 2026. CL51 (cloud) reference (blue), shown with its offset-corrected twin
(green, −0.5 % median, §4.2). The native CL51 median runs to 6 km with a huge IQR — its 40 m ripple,
a fixed additive error, swamps the weak signal aloft (driving the median through zero near ~2.9 km);
the offset-corrected twin rises cleanly. The CL61 (cloud, black) reads −16.5 % (log r 0.69) — a
910-vs-910 calibration difference (the CL61's own new multiple-scattering table + the removed window
correction, §3.5), not a wavelength effect; its Rayleigh twin is absent (only 4 usable molecular
nights).*

![sirta](figs_l1_validation/fig_sirta.png)
*Figure 7 — Palaiseau/SIRTA, Mar 2025–Feb 2026. CL31 (cloud) unfiltered median −40.1 % — its
detection ceiling: above ~0.7–1.5 km the CL31 reads noise while the CHM15k reference still sees
molecular signal (§3.1); Mini-MPL −52 % — the ill-conditioned 532→1064 conversion (§3.2): at native
532 nm the same instrument reads −3.0 % vs EARLINET (Figure 13).*

![lindenberg](figs_l1_validation/fig_lindenberg.png)
*Figure 8 — Lindenberg, 2025–2026 (18 months, N ≈ 862 k pairs). Both CL61 calibrations agree with
the CHM15k under the molecular-aware conversion: **cloud −1.7 %, Rayleigh +2.5 %**, log r ≈ 0.93 —
the flat-Ångström step it replaced gave +32/+35 % growing with altitude (§3.4).*

![aosta](figs_l1_validation/fig_aosta.png)
*Figure 9 — Aosta, 2025–2026. CL61 cloud +12.0 %, Rayleigh −5.0 % — the two disagree, a
station-specific degraded-window / cloud-calibration residual (82 % window, §3.4) that the wavelength
conversion cannot absorb.*

![camborne](figs_l1_validation/fig_camborne.png)
*Figure 10 — Camborne, 2025–2026. CL61 cloud **−2.9 %**, Rayleigh **+1.5 %**, log r 0.97 — both agree
with the CHM15k under the molecular-aware conversion (the flat step gave +28/+32 %, §3.4).*

![earlinet ari](figs_l1_validation/fig_earlinet_ari.png)
*Figure 11 — Leipzig vs CHM15k: 1151 matched ≥ 30-min profiles, med +4.4 %, log r 0.92 (unfiltered).*

![earlinet sir](figs_l1_validation/fig_earlinet_sir.png)
*Figure 12 — Palaiseau vs CHM15k: 196 matched, med −8.4 %, log r 0.75 (unfiltered; this system's
2000 m overlap limits the comparison to the free troposphere, and its 196 matches reaching to 5 km
are the noisiest set — the SNR-filtered value was −0.0 %).*

![earlinet ino](figs_l1_validation/fig_earlinet_ino.png)
*Figure 13 — Magurele unit B: 636 matched, med +1.1 %, log r 0.86. The linear r (0.42) is
event-driven — the clearest argument for the log-space metric.*

![earlinet ino_a](figs_l1_validation/fig_earlinet_ino_a.png)
*Figure 14 — Magurele unit A: 636 matched, med +2.8 %, log r 0.90 — 1.7 % from its co-located twin
against the same reference.*

![earlinet sir 532](figs_l1_validation/fig_earlinet_sir_532.png)
*Figure 15 — Native 532 nm: Mini-MPL Trappes vs SIRTA 532 (no wavelength conversion): 37 matched,
med −3.0 %, log r 0.86 — vindicates the instrument; the station −52 % is the conversion (§3.2).*

---

## 3. Where the significant differences come from

Three post-hoc decompositions on the exact paper matrices localize every significant difference in
time-of-day, season and altitude (`fig_report_*`; numbers in `discrepancy_analysis.json`).

![splits](figs_l1_validation/fig_report_daynight_seasonal.png)
*Figure 16 — Day/night + seasonal splits of med relbias (top) and log r (bottom) per channel.*

![altitude bands](figs_l1_validation/fig_report_altitude_bands.png)
*Figure 17 — Median relative bias per altitude band (0.5–1, 1–2, 2–3 km). A near-range mechanism
fades aloft; a conversion/noise mechanism grows where the aerosol fraction shrinks; a calibration
scale error is flat.*

![monthly](figs_l1_validation/fig_report_monthly_bias.png)
*Figure 18 — Monthly median relative bias. Flat = scale-like; seasonal cycle = WV/laser-drift
residual; step = instrument change.*

### 3.1 CL31: report a detection-limit profile, not a single median

The CL31 median is not a calibration number — it is dominated by *where the instrument stops
measuring*. Removing the SNR filter makes this explicit (§2.3): the Palaiseau CL31 unfiltered median
falls to −40.1 % because, above its detection ceiling, the CL31 reads noise while the CHM15k reference
still reads real molecular signal, so the column ratio is dragged toward −100 %. A single median
hides that; a **detection-limit profile** shows it directly, built like the operational sensitivity
product (`calibration.sensitivity.detection`). The minimum detectable attenuated backscatter at
averaging time τ is

```
β_att,min(r) = SNR · σ₀(r) · √(dt/τ),
```

with σ₀(r) the native-sampling noise (robust lag-1 profile-to-profile scatter, so the slow
atmospheric signal cancels), SNR = 3, τ = 30 min. The **maximum usable altitude** is the highest
range where β_att,min still sits below the clear-air molecular floor β_mol(r) (from CAMS T/p at the
instrument's own wavelength).

![CL31 detection-limit profile](figs_l1_validation/fig_sensitivity_profile_payerne.png)
*Figure 19 — Detection-limit profile at Payerne (June 2026, SNR 3 @ 30 min). Solid = each
instrument's noise floor β_att,min; black/grey = the clear-air molecular floor at 1064 nm (CHM) and
910 nm (CL31/CL61); dashed = the max usable altitude where the two cross. The CL31 and CL61 are both
910 nm, so they share the same molecular target: the CL61 clears it to 3.60 km, the CL31 only to
0.67 km — the same-wavelength, same-site comparison shows the CL61 is usable ~5× higher. The CHM15k
(1064 nm) reaches 2.07 km against its own intrinsically dimmer λ⁻⁴ floor.*

| instrument | λ | max usable altitude (SNR 3, 30 min) |
|---|---|---|
| **CL31** | 910 nm | **0.67 km** |
| CHM15k | 1064 nm | 2.07 km |
| **CL61** | 910 nm | **3.60 km** |

Caveats: this is the strict *molecular-detection* limit — a brighter aerosol layer (≳ 10× the
molecular floor) lets all three reach higher; and the absolute altitudes scale with the averaging
(at τ = 3 h the CL31 ceiling rises by ≈ √6). The ordering CL31 ≪ CHM15k < CL61 is the robust result,
and is why the CL31 is validated *only below ~0.7–1.5 km* (where its median is consistent with the
network) and reported as a profile above that.

### 3.2 Palaiseau Mini-MPL (−52 %): the 532→1064 conversion, not the instrument

Three convergent proofs. (i) **Native-wavelength closure**: −3.0 % vs EARLINET 532 nm (Figure 15).
(ii) **Altitude-band decomposition**: −35 / −47 / −74 % in the 0.5–1 / 1–2 / 2–3 km bands (Figure 17)
— growing exactly where the aerosol fraction shrinks, opposite of any instrument/overlap defect.
(iii) **Ångström insensitivity**: the median in-band aerosol fraction of the converted signal is
essentially zero (the extracted β_aer = β − β_mol·T² is a difference of two nearly equal numbers), so
no plausible α rescues the conversion.

![minimpl alpha](figs_l1_validation/fig_report_minimpl_alpha.png)
*Figure 20 — Mini-MPL bias vs assumed Ångström exponent: flat — the conversion, not α, is the
problem.*

A conversion physics error was also found and fixed in an earlier iteration (the molaer model
subtracted the *unattenuated* molecular from the *attenuated* 532 nm signal; using the attenuated
molecular recovered ≈ 5 points). The residual is the ill-conditioning: above the SNR gate the band
signal is ≈ 98 % molecular, so any residual few-% scale/model error is amplified ≈ 13× into the
converted 1064 nm value. **Consequence:** validate the Mini-MPL at its native wavelength; an elastic
532→1064 conversion by molecular subtraction is structurally unreliable in clean air.

### 3.3 Payerne CL61: the two calibration methods now agree

![cl61 methods](figs_l1_validation/fig_report_cl61ray_payerne.png)
*Figure 21 — The SAME physical C_L from the two methods. After the 2026-07 cloud recalibration (the
CL61 gets its own multiple-scattering table) the two series overlap instead of splitting.*

With the earlier cloud tables the two CL61 constants disagreed by ~12 % (C_L cloud 1.425 vs Rayleigh
1.251, ratio 0.885), which earlier drafts read as a genuine method discrepancy. **The 2026-07
recalibration removes it:** C_L(cloud) ≈ 1.18–1.21 vs C_L(Rayleigh) ≈ 1.22–1.25 — agreement to a few
percent (ratio ≈ 0.99–1.03 over the averaging window), and the two β rows now differ by only **0.7
points (+0.1 % cloud vs −0.6 % Rayleigh)**, both agreeing with the CHM15k once the molecular-aware
conversion is applied (§3.4). Two independent calibrations — strong-signal liquid-cloud O'Connor and
weak-signal Rayleigh molecular — land on the same lidar constant. The old "the cloud method is the
one to trust" reading was an artefact of the previous cloud table being ~14 % high, which happened to
cancel the WV correction into apparent agreement.

> **Superseded — CL61 Rayleigh-vs-cloud "~15 % method discrepancy"** (2026-07-02 drafts and the
> hood-dark AC-coupling offset hypothesis): the 2026-07 recalibration brings the two methods into
> agreement, so the entry is superseded — the (small) common CL61-vs-CHM15k offset is the wavelength
> conversion (§3.4), not either calibration. The full history — the termination-hood dark campaign,
> the AC-coupling pulse-response model, the multi-instrument hood significance test — is in the
> dedicated investigation [report 07](07_cl61_calibration.md)
> (report 07) and is **not** reproduced here.

### 3.4 CL61 vs CHM15k: the 910→1064 nm wavelength conversion — the fix, now applied

This is the paper's central methodological result; the full experiment is the companion report
[report 07](07_cl61_calibration.md). With the earlier
single-Ångström conversion the CL61 sat **+21 / +32 / +52 / +28 %** above the CHM15k (Payerne /
Lindenberg / Aosta / Camborne), growing with altitude and larger at night and in winter — exactly
the conditions of high molecular share, the fingerprint of a wavelength-conversion error rather than
of calibration or noise (Lindenberg the clean case: the bias ran +16.5 → +47 % across the 0.5–1 /
1–2 / 2–3 km bands, Figure 17).

The mechanism: a **single Ångström exponent on the total signal** mis-scales the molecular part,
which follows λ⁻⁴ (ratio β_mol(910)/β_mol(1064) = 1.877, King-corrected, Bucholtz 1995), not λ⁻ᵅ
with α ≈ 1. The flat exponent leaves the molecular part ≈ 1.6× over-scaled (0.855 applied where 0.534
is required), and the positive bias grows precisely where the molecular fraction is large.

**The component-separated molecular-aware conversion is now in the pipeline** (§1.4 step 5): the
analytic molecular part (β_mol ∝ p/T from CAMS) is removed and re-added at 1064 nm, and the Ångström
law scales the aerosol residual only. It brings the CL61 (cloud) to **+0.1 % (Payerne), −1.7 %
(Lindenberg), −2.9 % (Camborne)** — agreement to within ±3 %, with cloud and Rayleigh converging
(§3.3) — and it collapses the altitude tilt (the controlled treatment matrix confirms Lindenberg's
across-band spread falls 12.4 → 0.2 points; no multiplicative calibration-scale fix can do that). The
α ≈ 1 used is physical, not tuned, and the result is weakly sensitive to it.

![method summary](figs_l1_validation/fig_method_summary.png)
*Figure 22 — Wavelength-conversion treatment matrix across the CL61 sites: the single-Ångström step
(large, altitude-growing positive bias) vs the molecular-aware component-separated conversion
(near-zero, altitude-flat). Per-site panels: `fig_method_payerne.png`, `fig_method_lindenberg.png`,
`fig_method_aosta.png`, `fig_method_camborne.png`, `fig_method_uccle.png`.*

**Aosta** is the exception the conversion does *not* absorb: a residual **+12.0 %** with the two CL61
calibrations disagreeing (cloud +12.0 vs Rayleigh −5.0). The cause is a **window-transmission**
effect, not the wavelength conversion. Aosta's CL61 has a severely degraded window (median
transmission 82 %, p10 80 %) — the worst of the benchmark units. The *previous* cloud calibration
corrected β for the reported window transmission (β /= (T/100)², a +47 % inflation at Aosta) before
deriving the constant; that is physically wrong (the reported value is an arbitrary manufacturer-
scaled diagnostic, and the constant already absorbs the real window attenuation), so the 2026-07
recalibration reverted it to a **reject-only gate** (drop below 50 %, no β correction). That fix is
the bulk of Aosta's cloud shift (+16 % → +52 %) and of the cloud-vs-Rayleigh split (the Rayleigh
constant comes from a different source and never carried the window correction). The residual is
therefore a **degraded-window / cloud-calibration issue specific to this unit**, not the wavelength
methodology — which did its job. Recommendation: flag Aosta's degraded-window data.

> Consistent with the current-truth: the **window-temperature/transmission magnitude correction was
> removed 2026-07-06 — it is now reject-only**. Any older text proposing a window-T β correction is
> superseded.

### 3.5 Uccle CL61 (−16.5 %, log r 0.69): a 910-vs-910 calibration difference

Uccle is the one CL61 site with a **910 nm reference** (the CL51), so the wavelength conversion is a
no-op and the WV correction nearly cancels (WV on/off moves the median by only ~1 %, §4.1). The
−16.5 % is therefore a genuine **CL61-vs-CL51 calibration difference**, not a spectral or WV effect.
It sign-flipped from the +19 % of earlier drafts through two 2026-07 changes: the CL61 now carries
its **own** multiple-scattering table (a_G = 5.5 µm) instead of borrowing the CL51's (the dominant
part), and the removed window-transmission correction (§3.4) treated the CL51 reference (91 % window)
and the CL61 (94.6 %) differently — together a ~39-point relative shift between the two cloud
constants. With no usable CL61 Rayleigh series at Uccle (only 4 nights in six months, §2.3) there is
no independent tie-breaker, and the CL51 reference itself carries the network's strongest electronic
ripple (§4.2) and a 272-day / 34 % calibration oscillation (§2.2) sampled over only three months
(MAM). The Uccle CL61 number is the least transferable in the study: re-evaluate against the
offset-corrected CL51 over a full year, or against a 1064 nm reference through the §3.4 conversion.

### 3.6 Amsterdam CHM15k B (+17.3 %): a unit-specific daytime + near-range effect

Split day/night (§2.5): unit B reads +25 % by day against unit A but only +11 % at night. Split by
altitude band (Figure 17): +36 % at 0.5–1 km, +16 % at 1–2 km, +9 % at 2–3 km — concentrated in
daytime and in the lowest kilometre, fading aloft. A dedicated six-probe investigation pins it down:

![amsterdam B investigation](figs_l1_validation/fig_amsterdam_b_investigation.png)
*Figure 23 — Amsterdam unit B, six probes vs units A/C/D (unfiltered). (a) The bias follows a smooth
diurnal cycle peaking at 15 UT (+32 %), while C and D stay flat. (b) The daytime bias profile peaks
at ≈ +48 % near 500 m and fades aloft. (c) The daily bias is stable over the whole quarter — no
drift, no step. (d) The bias rises monotonically with solar elevation. (e) It rises equally with
internal temperature (18→35 °C: +8→+27 %). (f) The additive-vs-gain test: the median absolute
difference (B−A) by day does NOT follow a pure-gain shape — the excess is concentrated in the lowest
1.5 km beyond what a gain error produces.*

Reading: unit B carries (i) a ~+11 % night-time floor — notable because a pure gain difference would
be absorbed by its own nightly Rayleigh calibration, so even this floor is a profile-shape effect —
plus (ii) a daytime, near-range-weighted excess up to +20 points that tracks solar load and is
stationary over the quarter. Together they indict the unit's daytime signal handling (solar-
background / afterpulse subtraction), *not* the calibration (flat daily series, no step), *not* the
overlap model (< 1 % in-band, and the effect peaks near 500 m above the overlap region — Figure 5b(d)
shows B's 300–700 m excess is identical with and without the correction).

**Housekeeping cross-examination.** The full L1 housekeeping suite was compared four-ways:

![amsterdam B housekeeping](figs_l1_validation/fig_amsterdam_b_housekeeping.png)
*Figure 24 — Amsterdam A–D housekeeping. (a–d) diurnal composites of solar background, calibration
pulse, detector temperature and raw noise; (e) window transmission; (f) laser / detector quality;
(g) service-code frequency; (h) Spearman rank of unit-B daytime bias vs its own housekeeping.*

| marker | A | **B** | C | D |
|---|---|---|---|---|
| laser lifetime [h] | 15 207 | **49 904** | 21 154 | 50 371 |
| window transmission [%] | 95 | **86** | 95 | 86 |
| calibration pulse [ph/shot] | 0.068 | **0.051** | 0.074 | 0.056 |
| solar background, noon−night [ph/shot] | +0.038 | **+0.050** | +0.022 | +0.060 |
| records with service bits [%] | 2.0 | **0.0** | 0.01 | 5.5 |

The housekeeping sharpens the diagnosis by elimination. B and D form the **aged pair** — both lasers
at ≈ 50 000 h (3× A/C), both windows at 86 %, both with a weakened calibration pulse — yet **D agrees
with A to +1.4 % while B reads +17.3 %**: laser age, window contamination and calibration-pulse loss
are each acquitted as sufficient causes by the D control. Even the daytime solar-background load is
largest on D (+0.060 vs B's +0.050) with no bias. What is unique to B is that daytime load **converts
into a near-range positive β residual** — imperfect background/afterpulse subtraction inside this
unit's chain. Two practical conclusions: (i) B **self-reports clean** (0 % service bits, while the
healthy D flags 5.5 %), so status-based network QC cannot catch this defect — only co-location can;
(ii) B and D are both due for service on the aging markers, but only B's data is biased. C and D
remain the healthy pair, bounding the healthy CHM15k unit-to-unit spread near ±5 % (unfiltered).

### 3.7 Synthesis of mechanisms

| mechanism | fingerprint | affected results |
|---|---|---|
| **910→1064 nm conversion — molecular-aware (RESOLVED, now in pipeline)** | flat-Ångström over-scaled the λ⁻⁴ molecular part ≈ 1.6× → offset grew with altitude/night/winter; the analytic-Rayleigh (CAMS) form with α ≈ 1 removes it | CL61 vs CHM15k +21/+32/+28 % → +0.1/−1.7/−2.9 % (Payerne/Lindenberg/Camborne); §3.4 |
| **Window-transmission mishandling** (previous cal corrected β /= (T/100)², now reject-only) | flat in altitude; cloud vs Rayleigh disagree; worst at degraded windows | Aosta CL61 (82 % window, +47 % inflation removed) cloud +12.0 % vs Rayleigh −5.0 % (§3.4); part of Uccle |
| CL61 own multiple-scattering table (2026-07) | ~39-pt shift vs CL51; no wavelength/WV component (910-vs-910) | Uccle CL61 −16.5 % vs CL51 (§3.5) |
| 910 nm laser-line drift → WV residual + C_L oscillation | seasonal bias cycle (272 d Uccle, 136 d CL31); season-dependent medians | all CL31/CL51 results; Uccle reference (§3.5) |
| Digitizer fixed-pattern ripple | −11 % at 2–3 km on the Uccle CL51 itself; log r loss | Uccle above ~2 km (corrected twin provided, §4.2) |
| 532→1064 conversion ill-conditioning | grows with altitude; α-insensitive; native-λ closure fine | Mini-MPL station entry only (§3.2) |
| Unit-specific daytime/near-range behaviour | day ≫ night, fades with altitude | Amsterdam B (§3.6) |
| SNR floor / detection ceiling (single-diode) | unfiltered median → −40 % as CL31 reads noise above ~0.7 km while the reference still sees molecular signal; low log r | CL31 everywhere; detection-limit profile (§3.1) |
| CL31 detector-background (dark) | measured terminal-hood dark (P_dark·z², two-resonance R² = 0.98); removal raises log r 0.18 → 0.37 | Payerne CL31 hood-dark twin (§4.3) |

---

## 4. Water-vapour, electronic-offset and dark corrections in the validation

### 4.1 Water-vapour impact (910 nm family)

Each CL31/CL51/CL61 channel was run twice — WV correction ON vs OFF — through the otherwise
identical pipeline. Two numbers per instrument: the **in-band β increase** (median β_on/β_off − 1 ≈
median 1/T²_wv − 1) and the **change in median relative bias** vs the site reference.

![wv impact](figs_l1_validation/fig_wv_impact.png)
*Figure 25 — The WV correction raises in-band β_att by 14–23 % depending on site humidity and laser
line. Identical for the CL61 cloud and Rayleigh rows (an instrument/site property, not a method
property) — a built-in consistency check.*

| site | instrument | β increase | Δ med relbias |
|---|---|---|---|
| Payerne | CL31 | +13.7 % | +11.7 % |
| Payerne | CL61 (both methods) | +20.4 % | +15.8 / +18.0 % |
| Uccle | CL51 (reference) | +19.2 % | 0 (by construction) |
| Uccle | CL61 | +20.6 % | **+1.1 %** |
| Palaiseau | CL31 | +18.0 % | +15.4 % |
| Lindenberg | CL61 | +19.1 % | +24.1 / +23.1 % |
| Aosta | CL61 | +14.8 % | +14.0 / +16.6 % |
| Camborne | CL61 | +23.1 % | +23.2 / +23.1 % |

Two readings. (i) Against a **1064 nm reference** the correction shifts the comparison by its full
size (up to +24 points at Lindenberg) — without it every 910 nm instrument would sit ~15–25 % low:
the correction is not optional at 910 nm. (ii) Against a **910 nm reference** (Uccle: CL61 vs CL51)
the corrections nearly cancel (Δ +1.1 %) — the Uccle CL61 offset is **not** a WV artefact.

### 4.2 Uccle CL51 — clear-night ripple (no hood available)

The Uccle CL51 carries the strongest correctable electronic offset of the network scan
([report 08](08_overlap_nearrange_offset.md)): a fixed 40 m (4-gate) digitizer
ripple at 9.4 % of the mid-range signal, split-half reproducibility 1.00. The corrected twin
(`rcs_0 − b_phys`, same cloud C_L — the O'Connor method is offset-immune) vs the native unit:

| band | med relbias (corr vs native) |
|---|---|
| 500–3000 m (headline) | **−0.9 %** |
| 500–1000 m | −0.4 % |
| 1000–2000 m | −0.7 % |
| 2000–3000 m | **−11 %** |

The ripple is a *fixed additive* error: negligible against the strong low-level signal, an increasing
fraction of the weakening signal aloft (−11 % at 2–3 km, more above). Consequence: **Uccle results
that reach above ~2 km are reference-limited** — the corrected CL51 should be the reference for any
free-troposphere use (§3.5). The network-wide offset scan and the other flagged units (Diepenbeek,
Chilbolton, Kuopio, 6 × CL31) are covered in the offset report (consolidated as report 08).

### 4.3 Payerne CL31 — measured terminal-hood dark

Payerne is the **only** benchmark site with covered-telescope (terminal-hood) sessions, so its CL31
carries the *measured* dark rather than a clear-night proxy. Four hood sessions (2026-05-12, the
~25 h of 26–27 May, 09 June, 23 June; 4206 covered profiles) give the raw background in P = rcs_0/z²
space, and a two-resonance physical model — a fast AC-coupling amplifier ring (Λ ≈ 1053 m,
≈ 142 kHz, matching Kotthaus et al. 2016's 159 kHz high-pass corner) plus a slow transmitter ripple
(Λ ≈ 5080 m, ≈ 30 kHz) — fits it to **R² = 0.98** (a single damped sinusoid gives 0.74).

![CL31 hood dark model](figs_l1_validation/fig_cl31_offset_physical_model.png)
*Figure 26 — Payerne CL31 covered-telescope background. (a) the raw offset P = rcs_0/z² is the
superposition of two under-damped resonances; (b) the correction subtracted from L1 is that model
× z² (done in the raw-signal space, before calibration: `b_dark(z) = P_model(z)·z²`, the physically
correct place for a detector-background offset); (c) the two modes separated.*

**Applied correctly** — `rcs_0 − b_dark(z)` before the C_L scaling — the effect is diagnostic, not
cosmetic: **log r rises 0.18 → 0.37** and linear r 0.38 → 0.52 (the coherent dark that decorrelated
the CL31 from the CHM15k is gone), while the median moves *away* from zero (−7.9 → −19.0 %). The
reading is important: the native CL31's apparent "+21 %" high bias was largely its own **positive
dark pedestal** inflating β; removing the measured dark exposes that the CL31 genuinely does **not**
detect the weak molecular signal in the free troposphere (§3.1). The clear-night proxy captured only
~1/30 of it (|b| 0.5 vs 15.7 rcs_0 units) — the hood is what a single-diode CL31 needs. It is the
**correlation**, not the median-of-ratio, that shows the correction working, and the CL31 is still
best summarised by the detection-limit profile (§3.1), not by any single median.

---

## 5. L1 vs L2 calibration constant

*Native L1 (`D:/E-PROFILE_L1_2026`, binned to the 30 m L2 grid) vs the E-PROFILE L2 product, both
with the operational `eprof_v2` molecular-window method. The Rayleigh calibration spans the full L1
record (2025-01 … 2026-06); the cloud coefficient keeps the inter-comparison window.*

The historical L1↔L2 Rayleigh difference was a **method/grid interaction, not a data difference**:
the gated molecular-window methods over-reject the fine *native* grid. The fix is to bin native
L1/RAW to the L2 grid (30 m × 300 s, `l1_bin_to_l2_grid`) before calibrating. Switching from
`eprof_v1.2` to `eprof_v2` plus the longer window also raised the Payerne CHM series from 4 raw
nights (old figure) to 34 (L1) / 35 (L2).

![calibration L1 vs L2, 2025-2026](figs_extracted/l1_vs_l2_validation_01.png)
*Figure 27 — Calibration coefficient time series, L1 (native, binned to the L2 grid) vs L2, 2025–2026,
`eprof_v2`.*

| channel | calib | L1 median (n) | L2 median (n) | L1/L2 |
|---|---|---|---|---|
| Payerne CHM15k (Rayleigh) | Rayleigh | 6.070e+11 (479) | 6.240e+11 (475) | 0.973 |
| Payerne CL61 (Rayleigh) | Rayleigh | — (0) | 5.877e−01 (42) | — |
| Amsterdam CHM15k A | Rayleigh | 3.457e+11 (501) | 3.520e+11 (500) | 0.982 |
| Amsterdam CHM15k B | Rayleigh | 2.907e+11 (501) | 2.954e+11 (501) | 0.984 |
| Amsterdam CHM15k C | Rayleigh | 2.879e+11 (498) | 2.918e+11 (479) | 0.987 |
| Amsterdam CHM15k D | Rayleigh | 3.174e+11 (500) | 3.096e+11 (500) | 1.025 |
| Uccle CL61 (Rayleigh) | Rayleigh | — (0) | 1.120e+00 (155) | — |
| Palaiseau CHM15k (Rayleigh) | Rayleigh | — (0) | 4.765e+11 (360) | — |
| Palaiseau Mini-MPL (Rayleigh) | Rayleigh | 1.270e+05 (305) | 1.269e+05 (345) | 1.001 |
| Payerne CL31 (cloud) | cloud | 3.468e+00 (70) | 3.829e+00 (81) | 0.906 |
| Payerne CL61 (cloud) | cloud | — (0) | 9.804e−01 (71) | — |
| Uccle CL51 (cloud) | cloud | — (0) | 2.171e+00 (75) | — |
| Uccle CL61 (cloud) | cloud | 9.032e−01 (70) | 8.559e−01 (58) | 1.055 |
| Palaiseau CL31 (cloud) | cloud | 1.333e+00 (23) | 1.530e+00 (353) | 0.871 |

**Rayleigh** L1/L2 median = **0.985** (range 0.973–1.025) over 2025–2026 — binning the native L1 to
the L2 grid makes L1 and L2 calibrate identically. **Cloud** coefficients are physical O(1) (apparent
lidar ratio ~18–22 sr after the units fix): CL61 ~0.9, SIRTA CL31 ~1.5, the old Payerne CL31 ~3.5
(gain-confirmed); L1/L2 ~0.9–1.06. (Entries with L1 n = 0 are streams where only the L2 product
yielded a series in this run — reported, not hidden.)

---

## 6. L1 β_att validation (CSCS calibration applied to the raw signal)

*Generated by `run_l1_validation.py`. For every channel the native L1 range-corrected signal `rcs_0`
(2025–2026) is converted to β_att by applying the **CSCS calibration** —
`β_att = rcs_0 / C_L · 1e6`, with C_L the Kalman lidar constant from the CSCS calout
`<key>_kalman.csv` — then WV (910 nm), wavelength normalisation, cloud/fog screening, hourly +
common-altitude gridding, and statistics over 500–3000 m AGL vs the site reference (CHM15k Rayleigh).
No L2 product and no raw vendor files are used.*

This section demonstrates the calibration is **portable from the CSCS archive straight onto the raw
signal**, and confirms the median-C_L values that feed the operational dashboard. The relative biases
here predate the 2026-07 cloud recalibration and the molecular-aware wavelength conversion, so where
they differ from §2 the **§2 numbers are current**; this table is retained for the median-C_L column
and the linear-r / N sample sizes.

| site | channel | calib | rel. bias | r | N | median C_L |
|---|---|---|---|---|---|---|
| Payerne | CHM15k (Rayleigh) *(ref)* | rayleigh | +0.0 % | 1.000 | 623 411 | 6.41e+11 |
| Payerne | CL31 (cloud) | cloud | +51.8 % | 0.369 | 475 616 | 3.37e+07 |
| Payerne | CL61 (cloud) | cloud | −7.6 % | 0.986 | 152 638 | 1.46 |
| Payerne | CL61 (Rayleigh) | rayleigh | +6.9 % | 0.989 | 152 638 | 1.25 |
| Amsterdam | CHM15k A *(ref)* | rayleigh | +0.0 % | 1.000 | 627 500 | 3.44e+11 |
| Amsterdam | CHM15k B | rayleigh | +22.8 % | 0.958 | 582 250 | 2.82e+11 |
| Amsterdam | CHM15k C | rayleigh | −1.4 % | 0.931 | 556 750 | 2.87e+11 |
| Amsterdam | CHM15k D | rayleigh | −3.4 % | 0.966 | 584 000 | 3.16e+11 |
| Uccle | CL51 (cloud) *(ref)* | cloud | +0.0 % | 1.000 | 572 531 | 7.85e+07 |
| Uccle | CL61 (cloud) | cloud | +16.6 % | 0.854 | 149 094 | 1.17 |
| Uccle | CL61 (Rayleigh) | rayleigh | +28.1 % | 0.856 | 149 094 | 1.1 |
| Palaiseau | CHM15k (Rayleigh) *(ref)* | rayleigh | +0.0 % | 1.000 | 298 219 | 4.27e+11 |
| Palaiseau | CL31 (cloud) | cloud | −19.8 % | 0.867 | 228 997 | 7.91e+07 |
| Palaiseau | Mini-MPL (Rayleigh) | rayleigh | −61.0 % | 0.818 | 80 178 | 1.34e+05 |
| Lindenberg | CHM15k (Rayleigh) *(ref)* | rayleigh | +0.0 % | 1.000 | 850 750 | 3.93e+11 |
| Lindenberg | CL61 (cloud) | cloud | +24.4 % | 0.973 | 829 500 | 1.35 |
| Lindenberg | CL61 (Rayleigh) | rayleigh | +17.6 % | 0.974 | 829 500 | 1.45 |
| Aosta | CHM15k (Rayleigh) *(ref)* | rayleigh | +0.0 % | 1.000 | 780 558 | 2.47e+11 |
| Aosta | CL61 (cloud) | cloud | +3.8 % | 0.956 | 236 639 | 1.59 |
| Aosta | CL61 (Rayleigh) | rayleigh | +16.6 % | 0.949 | 236 639 | 1.3 |
| Camborne | CHM15k (Rayleigh) *(ref)* | rayleigh | +0.0 % | 1.000 | 385 436 | 2.2e+11 |
| Camborne | CL61 (cloud) | cloud | −1.2 % | 0.131 | 75 818 | 1.26 |
| Camborne | CL61 (Rayleigh) | rayleigh | +5.9 % | 0.131 | 75 818 | 1.16 |

The `figs_l1/` panels (`fig_l1_payerne.png`, `fig_l1_amsterdam.png`, `fig_l1_uccle.png`,
`fig_l1_sirta.png`, `fig_l1_lindenberg.png`, `fig_l1_aosta.png`, `fig_l1_camborne.png`) reproduce the
station layout for this CSCS-applied-to-L1 variant.

![payerne L1 validation](figs_l1/fig_l1_payerne.png)
*Figure 28 — Payerne, CSCS calibration applied to the raw L1 signal (CHM15k Rayleigh reference).*

---

## 7. Noise-filter sensitivity — keep `none`

*Window 2026-03-01 → 2026-06-30, Payerne (0-20000-0-06610). Reference: CHM15k (Rayleigh). Band
500–3000 m AGL. Driver: `_run_payerne_noisefilter.py`.*

Three selectable noise-filter modes (`ALC_VAL_NOISE_FILTER` / `intercompare.set_noise_filter`),
applied to the **native sub-hourly profiles** before the 60-min median. *(Verified in code:
`intercompare.py` exposes exactly `("none", "snr3", "cloudnet")` with the default resolving to
`none` unless `ALC_VAL_L1_SNR != 0`.)*

- **none** — no detection gate; the ungated signal (symmetric noise averages out in the median).
- **snr3** — per-gate window SNR ≥ 3 (|median| / (robust σ / √n)); the operational sensitivity gate.
- **cloudnet** — the Cloudnet / ceilopyter per-profile filter (`ceilopyter.noise.remove_noise`): each
  sample below 5× the per-profile far-range noise (std of the top 10 % of range gates, on the
  non-range-corrected signal `rcs₀/r²`) is masked, then the survivors are median-averaged.

**Median relative bias vs reference [%] (500–3000 m):**

| Channel | No filter | SNR ≥ 3 | Cloudnet filter |
|---|---|---|---|
| CL31 (cloud) | −7.9 | −4.6 | +85.6 |
| CL61 (cloud) | +0.1 | −0.5 | −10.5 |
| CL61 (Rayleigh) | −0.6 | −1.1 | −8.2 |
| CL31 (cloud, hood-dark corr) | −19.0 | +6.6 | +107.8 |

**Retained hourly-gate samples N (500–3000 m):**

| Channel | No filter | SNR ≥ 3 | Cloudnet filter |
|---|---|---|---|
| CHM15k (Rayleigh) *(ref)* | 233 132 | 210 800 | 121 580 |
| CL31 (cloud) | 199 231 | 78 240 | 7 469 |
| CL61 (cloud) | 157 815 | 140 544 | 77 221 |
| CL61 (Rayleigh) | 157 815 | 140 544 | 77 221 |
| CL31 (cloud, hood-dark corr) | 199 231 | 78 240 | 7 469 |

**Interpretation.**

- **Data retained collapses as the gate tightens.** The Cloudnet filter keeps ~31 % of the samples
  `none` keeps (SNR ≥ 3 ~68 %), and it is *instrument-dependent*: the noisiest unit (CL31 cloud)
  drops from N = 199 231 to N = 7 469.
- **Cloudnet is a per-profile filter, so it conditions on positive noise excursions before
  averaging.** For a low-SNR instrument this pushes the retained median UP with altitude: CL31 (cloud)
  median bias goes −8 % → +86 %. The tell-tale of censoring bias is the divergence between the median
  and the linear-mean relbias, and the drop in log r for the clean CL61 channel (its molecular points,
  which anchored the correlation, are removed).
- **SNR ≥ 3 is milder** — a window-median gate (drops a gate only if the hourly median is not 3σ-
  significant), not a per-profile cherry-pick — but it still cuts N and shifts the noisiest channels.
- **`none` preserves the symmetric-noise-averages-out property** that an unbiased instrument-vs-
  instrument bias comparison needs, which is why it is the validation default.

**Recommendation (the finding of this study): keep `none` for the calibration/bias validation.** The
Cloudnet filter is destructive for the bias comparison — it degrades the good CL61 (Payerne cloud
+0.1 → −10.5 %) and inflates the noisy CL31 (−7.9 → +85.6 %); it is faithful to Cloudnet's
*product-generation* pipeline (per-profile screening + range sparsity are acceptable there) but is
the wrong tool for estimating a mean bias from the survivors. Use `snr3` only when a detection-limited
view is explicitly wanted.

![Noise-filter comparison](figs_payerne_noise_filter/fig_payerne_noisefilter_comparison.png)
*Figure 29 — Median relative bias per channel under the three noise-filter modes.*

![Payerne none](figs_payerne_noise_filter/fig_payerne_none.png)
*Figure 30 — Payerne, no filter (the validation default).*

![Payerne snr3](figs_payerne_noise_filter/fig_payerne_snr3.png)
*Figure 31 — Payerne, SNR ≥ 3 (legacy operational sensitivity gate).*

![Payerne cloudnet](figs_payerne_noise_filter/fig_payerne_cloudnet.png)
*Figure 32 — Payerne, Cloudnet per-profile filter — visibly destructive for the bias comparison.*

---

## 8. Limitations, caveats and reproducibility

### 8.1 Limitations and caveats

- **Uccle CL61 (Rayleigh):** no usable L1 calibration series — the molecular fit converged on only
  4 nights in six months (the stream's ~50 % duty-cycle data delivery keeps the ≥ 3 h-of-clear
  Rayleigh gate from ever being met on short summer nights). Skipped, reported (§2.3).
- **Leipzig `lei` / Cabauw `cbw`:** no EARLINET 1064 nm files in the 2025–2026 window.
- **`sir_532`:** only 37 matched profiles (the L1 screening is stricter than the earlier L2 path) —
  the −3.0 % closure is robust in sign but has ~±3 % sampling uncertainty.
- **Overlap correction** is validated as *harmless in-band* here (< 0.1 %); its actual benefit is
  below 500 m, outside these statistics by design.
- **Electronic-offset patterns** are Mar–Jun 2026 extractions; a firmware/board change requires
  re-extraction. Only Uccle is corrected in this study (the strongest case); the other flagged units
  (Diepenbeek, Chilbolton, Kuopio, 6 × CL31) follow the same one-line recipe (report 08).
- **Three-month stations** (Payerne, Amsterdam, Uccle: Mar–Jun 2026) sample one season; every 910 nm
  number there inherits the §2.2 oscillation phase. The 18-month sites (Lindenberg, Aosta, Camborne,
  SIRTA) average it.
- **Absolute 532 nm ↔ 1064 nm intercomparison (Mini-MPL) is intrinsically limited** — molecular
  backscatter scales as ≈ λ⁻⁴, so the 532 nm signal is molecular-dominated and the absolute bias is
  bracketed by the wavelength model. The Mini-MPL is an honest 532 nm *tracker* (native-wavelength
  closure −3.0 % vs EARLINET, temporal correlation, time-height consistency); its converted station
  entry measures the conversion's conditioning, not the instrument (§3.2).
- **910 nm cloud-calibrated channels** depend on the season sampled because of the laser-line-drift
  C_L oscillation (§2.2, §3.5) that the fixed-line WV model cannot fully cancel.

### 8.2 Reproducibility

From the repo root (branch `wv-correction`), with the `calibration` package importable. Inputs are
env-overridable; the defaults are the paths used for this report.

```bash
# Inputs (defaults shown)
export ALC_VAL_CALOUT="C:/DATA/Projects/202606_E-PROFILE_calibration/E_PROFILE_calout_2025_2026"  # CSCS cloud+Rayleigh calout
export ALC_VAL_L1_ROOT="D:/E-PROFILE_L1_2026"          # native L1 rcs_0 archive (the ONLY beta input)
export ALC_OVERLAP_DIR="D:/TEMP_MODELS/202606"         # CHM15k temperature-dependent overlap models (127 units)
export ALC_VAL_CAMS_04="D:/CAMS_Monthly_04"            # 0.4 deg monthly CAMS for the molecular-aware conversion + WV
# 1 deg CAMS at D:/CAMS is the per-month fallback; L2_monthly at A:/E-PROFILE_L2_monthly (CL61-Rayleigh series only)
# SNR filter is OFF by default for stations (Section 7); set ALC_VAL_L1_SNR=1 to restore the legacy SNR filter.

# 1. Feed the operational calout into the report's calibration series (writes calib/<key>{,_L1,_L2}.csv)
python -m validation.paper.dashboard_to_calib          # CL61-Rayleigh (absent in calout) preserved from calib_benchmark

# 2. Regenerate the validation, everything UNFILTERED (per-site + EARLINET figures, fig_calib_timeseries,
#    fig_wv_impact, summary_stats.csv). Stations and EARLINET both run with no SNR filter.
python -m validation.paper.run_paper_validation

# 2b. (optional) legacy SNR-filtered comparison for a cross-check
ALC_VAL_L1_SNR=1 python -m validation.paper.run_paper_validation

# 3. Regenerate the discrepancy figures + discrepancy_analysis.json (day/night, seasonal, altitude bands, monthly)
python -m validation.paper.discrepancy_analysis full

# 4. The wavelength-conversion experiment (companion report): treatment matrix + fig_method_*.png
python -m validation.paper._cl61_methodology_experiment

# 5. The CL31 detection-limit profile (Figure 19): fig_sensitivity_profile_payerne.png
python -m validation.paper._sensitivity_profile_test

# 5b. The Payerne CL31 terminal-hood dark model (writes cl31_b_dark.npz; re-run only if the hood sessions change)
python -m validation.paper._cl31_offset_physical_model

# 6. The noise-filter sensitivity (Section 7)
python -m validation.paper._run_payerne_noisefilter
```

`dashboard_to_calib` reads the SAME `$ALC_VAL_CALOUT` as `intercompare.CALOUT`, so the report's
Kalman constants and the L1 pipeline's constants cannot drift apart. To rebuild the calibration
constants from scratch (no CSCS calout) instead of step 1, run
`python -m validation.paper.calib_benchmark`.

**Machine-readable outputs:** `figs_l1_validation/summary_stats.csv` (28 rows),
`discrepancy_analysis.json` (all split/band/monthly numbers).

**Fog exclusion** (CL31/CL51/CL61 `vertical_visibility` > 0 ⇒ fog ⇒ excluded) is applied in the
screening step (`intercompare.screen`), verified on Edmonton 2026-02-01 where 66/155 night profiles
were fog.
