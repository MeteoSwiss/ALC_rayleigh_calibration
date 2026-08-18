# Water-vapour absorption correction for ~910 nm ceilometers — literature, wavelength configuration and sensitivity

*Consolidated 2026-07-10 (M. Hervo, MeteoSwiss E-PROFILE ALC). Merges: attbsc_wv_literature_review.md, wv_fwhm_literature_review.md, wv_wavelength_sensitivity.md, payerne_cl61_wv_sensitivity.md. Incorporates the CAMS-resolution figure set (`figs_wv_resolution/`, §5), which had no separate markdown source.*

## Contents

- [0. Operational facts (cheat-sheet)](#0-operational-facts-authoritative)
- [1. Literature basis for the WV correction](#1-literature-basis-for-the-wv-correction)
- [2. Laser emission spectrum per instrument (λ₀ and FWHM)](#2-laser-emission-spectrum-per-instrument-λ-and-fwhm)
- [3. Sensitivity of T²_wv to the wavelength configuration](#3-sensitivity-of-twv-to-the-wavelength-configuration)
- [4. CL61 Rayleigh-calibration sensitivity to the WV correction (Payerne)](#4-cl61-rayleigh-calibration-sensitivity-to-the-wv-correction-payerne)
- [5. Resolution sensitivity of the WV correction (CAMS grid)](#5-resolution-sensitivity-of-the-wv-correction-cams-grid)
- [6. Implementation notes and reproduction](#6-implementation-notes-and-reproduction)
- [Sources](#sources)

---

## 0. Operational facts (authoritative)

- **The WV correction is MANDATORY for all ~910 nm instruments** (Vaisala CL31 / CL51 / CL61).
  A no-WV "degraded" mode is **rejected** — omitting the correction worsens the result (§4), so
  the pipeline refuses to emit a WV-biased constant.
- **Humidity source = CAMS model levels (L137)** (`wv_source='cams'`), dense in the boundary
  layer. An **ERA5 path exists but is research-only** (Earth Data Hub ERA5 is a 19-level pressure
  subset with a coarser BL); it stays out of operations.
- **CAMS resolution:** operational monthly **0.4°** with a per-month **1° fallback**; both are
  L137 model levels, so vertical resolution is preserved. The same CAMS field feeds the cloud-WV
  correction and OmB alike.
- **A 910 nm night without usable CAMS is flagged, never calibrated WV-free.** Verified in code
  (`calibration/rayleigh/calibration.py`, `calibration/cloud/calibration.py`): no CAMS file →
  `flag = -4`; station outside the CAMS domain → `flag = -10`; unusable WV profile or a
  transmission grid mismatch → `flag = -4`; in the cloud path a failed / all-ones / empty
  transmission raises a hard error ("no fallback: this period is NOT calibrated"). The correction
  is gated on the laser **wavelength** (900–920 nm), not the instrument name, so 1064 nm (CHM15k)
  and 532 nm (Mini-MPL) are never WV-corrected by physics rather than by a hard-coded list.
- **CHM15k (1064 nm) is the WV-free reference:** outside the H₂O absorption band, it needs no
  correction and is used as the molecular-vs-molecular anchor in the sensitivity tests below.

---

## 1. Literature basis for the WV correction

### 1.1 Bottom line

The **Wiegner & Gasteiger (2015)** "WAPL" scheme is the **singular, anchor method** for
water-vapour (WV) absorption correction of 905–910 nm ceilometer attenuated backscatter. Every
later work either **validates it** (Wiegner et al. 2019), **substitutes a cheaper empirical
parameterisation** for operational use (Hopkin et al. 2019, the E-PROFILE operational reference),
or **acknowledges but neglects it** (Kotthaus et al. 2016). It is the standard, and it is
essentially unique — there is no competing rigorous formulation.

**No open, runnable implementation of a Wiegner-style *spectral* (HITRAN/MT-CKD line-by-line or
LUT) WV correction exists in any public repository.** CloudnetPy does **not** implement it
(confirmed by direct code inspection — see §1.4); the E-PROFILE/EUMETNET operational lineage uses
the simpler empirical `Twv = 1 − 0.17·IWV^0.52`. **Our MATLAB code
(`compute_wv_transmission.m`, `wv_t2eff.m`, `apply_wv_correction_to_L2.m`) and its Python port
(now `calibration/water_vapor_correction/water_vapor.py`) appear to be the only runnable open
implementation of the spectral WAPL-style correction** — a genuine novelty worth stating in the
paper.

### 1.2 Summary table

| Paper | Instrument / λ | Humidity source | Spectroscopy | WV correction method & magnitude | Code |
|---|---|---|---|---|---|
| **Wiegner & Gasteiger 2015**, AMT 8, 3971 ([doi](https://doi.org/10.5194/amt-8-3971-2015)) | Vaisala CT25k/CL31/CL51, 905–910 nm (not 1064 nm) | (method paper; any T,q profile) | **HITRAN 2005 + MT-CKD continuum**, ARTS line-by-line 895–930 nm @0.01 cm⁻¹, stored netCDF LUT (0.1 cm⁻¹, 10 m) | Spectral LUT; **averages over Gaussian laser spectrum** (λ₀+FWHM≈3.4 nm) → effective T²_w,eff. **Ignoring it biases β_p by ~20 % mid-lat, >50 % tropics (worst case; ~35 % typical tropical)** | WAPL netCDF archive = **private/on-request**, not open |
| **Wiegner et al. 2019** (CeiLinEx2015), AMT 12, 471 ([doi](https://doi.org/10.5194/amt-12-471-2019)) | Vaisala CL51 etc. vs RALPH ref. lidar | radiosonde / model | uses W&G 2015 LUT | **Validation** of W&G 2015 (multiplies signal by T_w,eff⁻²); near-range transmission agreement ~1–5 % | none |
| **Hopkin et al. 2019**, AMT 12, 4131 ([doi](https://doi.org/10.5194/amt-12-4131-2019)) — *E-PROFILE operational ref.* | CL31/CL51/CT25k/CS135 @910 nm + CHM15k @1064 nm | **NWP** (Met Office UKV; ECMWF profiles *provided by M. Hervo, MeteoSwiss*) | **none** (empirical) | O'Connor liquid-cloud calib (S=18.8±0.8 sr) **+ empirical** `Twv = 1 − 0.17·IWV^0.52` (Markowicz 2008). Agrees with WAPL **to within 2 %**; ~12 % annual cycle if ignored | multiple-scattering code (Hogan 2006) public; **no WV repo** |
| **Kotthaus et al. 2016**, AMT 9, 3769 ([doi](https://doi.org/10.5194/amt-9-3769-2016)) | Vaisala CL31, 905±10 nm, FWHM~4 nm, 0.3 nm/K | — | — | **Acknowledges** WV sensitivity (cites W&G 2015, Markowicz 2008) but **neglects it** (c_absolute=1) | none |
| Markowicz et al. 2008 | CT25k, ~905 nm | — | — | Case-specific predecessor; basis for Hopkin's empirical fit | none |
| **CloudnetPy** (ACTRIS-Cloudnet) | CL31/CL51, CHM15k | — | **none** | **No WV correction.** Scalar site `calibration_factor` (Vaisala 1.0, Lufft 3e-12, "probably incorrect") + O'Connor/Hogan liquid-cloud classification | [github.com/actris-cloudnet/cloudnetpy](https://github.com/actris-cloudnet/cloudnetpy) — open, but **no WAPL** |

### 1.3 Key points

1. **W&G 2015 is the method to cite** as the rigorous benchmark. Its three defining ingredients —
   HITRAN+MT-CKD cross-sections, line-by-line LUT, and **averaging over the laser emission
   spectrum (λ₀ + FWHM)** — are exactly what our implementation reproduces.
2. **Magnitude**: ~20 % on retrieved backscatter at mid-latitudes (>50 % tropics, worst case).
   Quote the >50 % as an *upper bound* (W&G stress there is "no generally applicable value" — it
   depends on height, aerosol load, algorithm). The 2025 CL51 study (Jin et al., MDPI rs17122013)
   independently states ">20 % if water vapor correction is ignored." *Caveat:* that paper
   reports the magnitude but does **not** adopt the full WAPL pipeline — do not cite it as a WAPL
   application.
3. **Operational vs rigorous split**: Hopkin et al. (2019) deliberately chose the cheap empirical
   `Twv = 1 − 0.17·IWV^0.52` over WAPL "because it requires a radiative-transfer model or access to
   their WAPL database." Notably their WAPL-vs-empirical comparison (their Fig. 5) **used ECMWF
   water-vapor profiles provided by M. Hervo** — a direct link between this author and the
   operational reference. The two agree to 2 %.
4. **Our novelty**: a full **spectral WAPL-style correction, in both MATLAB and open Python**
   (validated against MATLAB to ≤0.4 % and against ACTRIS-Cloudnet `atmoslib` — see
   `tests/WATER_VAPOR_AUDIT.md`), driven by CAMS humidity + the ECMWF L137 model levels. This is
   the **only open runnable spectral implementation** found.
5. **CloudnetPy confirmation**: direct inspection of CloudnetPy (HEAD 0db960b8) found **zero**
   references to HITRAN / MT-CKD / WV transmission / spectral averaging in the ceilometer code; it
   uses a scalar calibration factor + O'Connor liquid-cloud classification. So CloudnetPy is *not*
   a WV-correction reference — `atmoslib` (its thermodynamics library) is the right external anchor
   for the humidity primitives only.

### 1.4 Open follow-ups (not resolved by the literature search)

- The **Vande Hey (2015) thesis** and **Madonna et al. (2018)** were named but did not surface as
  verified claims — worth a manual check for any independent spectral WV treatment.
- Whether the original **W&G WAPL netCDF archive** was ever published openly (Zenodo / LMU
  institutional) so it could be cited/reused rather than regenerated.
- State explicitly in the paper which **HITRAN edition + MT-CKD version** our
  `abs_cross_647_full_levels_1000.nc` LUT uses, vs W&G's HITRAN-2004/MT-CKD basis.

---

## 2. Laser emission spectrum per instrument (λ₀ and FWHM)

### 2.1 Why the FWHM matters

The WV correction weights the (sharply structured) H₂O absorption cross-section by the laser
**emission spectrum**, assumed Gaussian with central wavelength λ₀ and full width at half maximum
Δλ (FWHM). Wiegner & Gasteiger (2015) introduced this and rank the required inputs **by decreasing
relevance: (1) the water-vapour profile, (2) the central wavelength λ₀, (3) the spectral width
Δλ.** The spectral width is the *least* sensitive of the three — but it is still needed, and it is
the worst-documented, especially for the CL61.

Three different "widths" appear in the datasheets and must not be confused:
- **emission FWHM (Δλ, spectral)** — the quantity used by the WV correction (this section);
- **receiver optical-filter bandwidth** — 36 nm for CL31/CL51 (Wiegner et al. 2019); irrelevant to the weighting;
- **pulse duration (temporal FWHM)** — e.g. 160 ns for the CL61; not a spectral width.

### 2.2 Summary — literature/manufacturer values vs our operational config

| Instrument | central λ₀ | **emission FWHM Δλ** | source(s) | our config | Payerne Qmini (measured, 2026-06-02) |
|---|---|---|---|---|---|
| **CL31** | 905 ± 10 nm (older) → 910 ± 10 nm (newer), 25 °C | **≈ 4 nm** (typical); 2021 datasheet omits it | Vaisala CL31 datasheet/ARM VCEIL handbook; Kotthaus et al. 2016; Wiegner et al. 2019 | (909.7, **6.0**) | 909.7, **5–7 nm**, peak wandering 909.0–910.1 nm |
| **CL51** | 910 ± 10 nm @25 °C, drift 0.27 nm K⁻¹ | **3.4 nm** (Vaisala); 3.5 nm used; 1–4 nm explored | **Wiegner & Gasteiger 2015**; Wiegner et al. 2019; Cordoba-Jabonero / MDPI 2025 | (910.0, **3.4**) | — (not measured) |
| **CL61** | **910.55 nm** | **not documented** — only λ₀ is published | Vaisala CL61 User Guide M212475EN-E; Le et al. 2026; Looschelders et al. 2025; Laffineur et al. 2026 | (910.74, **1.0**) | **910.74 ± 0.10**, single narrow line (true ≤ 0.03 nm; apparent ≈1.5 nm spectrometer-limited) |
| CS135 (Campbell, context) | 912 nm (stable) | ±3.5 nm | Wiegner et al. 2019 | — | — |
| CHM15k (context) | 1064 nm | n/a — outside the H₂O band, no WV | — | (1064.47, 0.5) | 1064.47 ± 0.10 |

The operational values are set in `calibration/water_vapor_correction/water_vapor.py`
(`LASER_SPECTRUM`): CL31 (909.7, 6.0), CL51 (910.0, 3.4), CL61 (910.74, 1.0), CHM15k (1064.47,
0.5), CHM8k (1064.0, 0.5) — matching this table.

### 2.3 CL51 — the only well-documented case (Δλ ≈ 3.4 nm)

[Wiegner & Gasteiger (2015)](https://doi.org/10.5194/amt-8-3971-2015) is the reference for the WV
correction and is explicit (p. 3975):

> *"For the CL51 ceilometer, e.g., λ₀ = 910 ± 10 nm at 25 °C with a drift of 0.27 nm K⁻¹ is
> specified by Vaisala. We assume a Gaussian shape of the spectrum with λ₀ between 901 and 919 nm,
> and a full width at half maximum (FWHM … Δλ) between 1.0 and 4.0 nm. **According to Vaisala, Δλ is
> of the order of 3.4 nm.**"*

For their analysis they adopt a "realistic value of **Δλ = 3.5 nm**" and explore Δλ = 2.5, 3.0,
3.5, 4.0 nm (their Fig. 2). [Wiegner et al. (2019, CeiLinEx2015)](https://doi.org/10.5194/amt-12-471-2019)
re-state λ = 910 ± 10 nm for the CL31/CL51 and use the same framework; the recent CL51 WV study
[(Remote Sens. 17, 2013, 2025)](https://doi.org/10.3390/rs17122013) follows Wiegner & Gasteiger.
**Our operational value `CL51 = (910.0, 3.4)` is exactly the Vaisala/W&G figure** — no change
needed.

### 2.4 CL31 — broad multimode diode, ≈ 4 nm, poorly pinned

The CL31 uses a pulsed **multimode** InGaAs diode. Older specs and the ARM VCEIL handbook give
λ₀ = 905 ± 10 nm at 25 °C with a typical FWHM ≈ 4 nm; newer units are quoted at 910 ± 10 nm
(Wiegner et al. 2019). Notably the **current Vaisala CL31 datasheet (B210415EN, 2021) lists no
laser wavelength or spectral width at all** — only "pulsed diode laser, Class 1M". Kotthaus et al.
(2016), the standard CL31 processing reference, discusses the instrument but not a precise emission
linewidth.

Our Qmini measurement gives a **broader 5–7 nm FWHM with the peak wandering 909.0–910.1 nm** between
acquisitions — consistent with the multimode diode redistributing power among longitudinal modes,
and broader than the ≈ 4 nm in the older literature. Our operational `CL31 = (909.7, 6.0)` sits in
the measured range (and above the ~4 nm literature value); the breadth is real and is why a single
nominal wavelength is only approximate for the CL31.

### 2.5 CL61 — only the centre (910.55 nm) is published; the FWHM is undocumented

This is the key gap. **Every** source we found gives only the central wavelength, never an emission
linewidth:

- **Vaisala CL61 User Guide (M212475EN-E)**: laser wavelength **910.55 nm**, InGaAs diode; the only
  "FWHM" in the spec table is the **pulse duration, 160 ns (temporal)** — not spectral
  (as tabulated by Le et al. 2026, their Table 1).
- [Le et al. (2026)](https://doi.org/10.5194/egusphere-2025-6331) — the most thorough CL61
  performance study (4 ACTRIS sites) and its supplement: 910.55 nm only; no emission spectrum.
  (Their abstract loosely writes "905 nm" for the molecular discussion, but the spec is 910.55 nm.)
- [Looschelders et al. (2025)](https://doi.org/10.1002/met.70088) — "pulsed laser diode at a
  wavelength of 910.55 nm"; no FWHM.
- [Laffineur et al. (2026, CONIOPOL, Uccle)](https://doi.org/10.5194/egusphere-2026-948) —
  "910.55 nm InGaAs diode laser"; no FWHM.

So **the CL61 emission linewidth has not been reported in the literature or specified by Vaisala.**
To our knowledge our Qmini campaign (910.74 ± 0.10 nm, a single **narrow** line, true width
≤ 0.03 nm, apparent ≈ 1.5 nm limited by the spectrometer) is the **first measured characterisation
of the CL61 emission spectrum** — a genuine contribution. Two consequences for the WV correction:

1. The CL61 line is much narrower than the CL31/CL51 multimode diodes (3–7 nm). It is **not** safe
   to inherit their ~3–4 nm widths, but it is **also not** safe to treat the CL61 as monochromatic
   (see §3.3): the measured ±0.10 nm centre uncertainty would then dominate. Our operational
   `CL61 = (910.74, 1.0)` keeps a deliberate ~1 nm averaging for robustness.
2. We measure **910.74 nm**, 0.19 nm above the Vaisala spec **910.55 nm**. The literature uniformly
   uses 910.55; at a ~1 nm bandwidth the difference is < 0.5 % on the WV correction (§3.1), so
   either is adequate, but the measured value is the more defensible one to report.

Relatedly, Le et al. (2026, Fig. S2) show the CL61 **laser temperature cycling** (≈ 19–22 °C,
~127 s period). At the CL51 drift of 0.27 nm K⁻¹ such swings imply sub-nm centre drift, reinforcing
that a finite (~1 nm) spectral averaging — rather than a fixed monochromatic line — is the robust
choice for the CL61 WV correction.

### 2.6 Context — other instruments

- **Campbell CS135**: λ = 912 nm (stable), spectral width ±3.5 nm (Wiegner et al. 2019).
- **Lufft CHM15k/CHM8k**: 1064 nm — outside the H₂O band, **no WV correction**; our measured
  1064.47 nm matches the manufacturer.

### 2.7 How the literature treats the FWHM sensitivity

- **Relevance order (Wiegner & Gasteiger 2015):** water-vapour profile ≫ central wavelength λ₀ >
  spectral width Δλ. *"The variability with Δλ depends on λ₀ but is in most cases comparably
  small."* Absorption rises by ~×2 from λ₀ = 905→908 nm but varies < 10 % from 908→918 nm.
- **Central-wavelength sensitivity (Wiegner et al. 2019):** *"the transmission is much more
  sensitive to errors of the assumed wavelength λ_on than to errors of the water-vapour content …
  dT_w,eff/dλ > 0.02 nm⁻¹"* (> 2 % per nm).
- **Recurring caveat:** several papers note that *"the unknown emission spectrum of the diode
  laser … can introduce significant errors"* — i.e. the spectrum is treated as an unmeasured
  manufacturer input.

This literature consensus **matches our own sensitivity study** (§3, figure there) exactly: at fixed FWHM the
CL61 910.55-vs-910.74 nm choice is ≈ 0.4–0.5 %; the FWHM itself is a ≈ 1–4 % lever; and the CL31
broad band self-averages to < 1 % despite its wandering. The one regime the older literature does
not emphasise — because it predates the narrow-line CL61 — is that a **near-monochromatic**
treatment makes T²_wv swing ≈ 27 % within the ±0.10 nm centre uncertainty, which is why the ~1 nm
averaging must be kept.

---

## 3. Sensitivity of T²_wv to the wavelength configuration

**Question.** The WV correction weights the H₂O absorption cross-section by the laser emission
spectrum (λ₀ and FWHM) to get the two-way transmission T²_wv. Across 905–915 nm the H₂O
cross-section varies by orders of magnitude over fractions of a nanometre (Fig. a), so the assumed
λ₀ and FWHM matter. How much does the *wavelength configuration* change T²_wv — the manufacturer
value vs the measured one, the spectral width, and (for the CL31) the spectral breadth and
inter-acquisition wandering?

**Emission wavelengths** (Qmini fibre spectrometer, Payerne, 2 June 2026; scale verified on the
O₂ A-band to 0.1 nm): CHM15k 1064.47 ± 0.10 nm; CL61 910.74 ± 0.10 nm (single narrow line,
apparent ≈1.5 nm FWHM spectrometer-limited, true scatter ≤0.03 nm); CL31 ≈5–7 nm FWHM centred
≈909.7 nm, peak wandering 909.0–910.1 nm between acquisitions. Manufacturer CL61 value:
**910.55 nm**. The operational WV code uses the **measured** values (`CL61 = (910.74, 1.0)`,
`CL31 = (909.7, 6.0)`).

**Metric:** median T²_wv over the **500–3000 m** comparison band (drives the validation bias, since
β_corr = β / T²_wv) and over the **2–6 km** Rayleigh window (drives the lidar constant). "β impact"
is the resulting change in corrected backscatter (β ∝ 1/T²_wv). Profiles from CAMS at Payerne.

*(Figure — WV correction vs laser-wavelength configuration at Payerne — not archived in this
report tree; the source `figs_paper_validation/wv_wavelength_sensitivity.png` was not preserved.
It showed: (a) the H₂O absorption cross-section across 905–915 nm at ~3 km (log scale) with the
normalised laser spectra overlaid — CL61 manufacturer 910.55 nm, CL61 measured 910.74 nm, CL31
909.7 nm / 6 nm FWHM — the CL61 narrow lines sampling individual absorption features while the CL31
broad band averages over many; (b) two-way WV transmission T²_wv vs range for the three
configurations with the CL61 FWHM 0.1–1.5 nm envelope shaded and the 500–3000 m comparison band
greyed; (c) median T²_wv over 500–3000 m for every configuration. Reproduce with
`wv_wavelength_sensitivity.py`, §6.)*

### 3.1 CL61 — manufacturer 910.55 nm vs measured 910.74 nm

At the operational FWHM (1.0 nm) the 0.19 nm difference between the manufacturer and measured
central wavelengths changes T²_wv by **< 0.5 %**, stable across the season:

| month | T²_wv (910.55) | T²_wv (910.74) | β difference |
|---|---|---|---|
| 2026-02 | 0.848 | 0.852 | 0.43 % |
| 2026-03 | 0.861 | 0.864 | 0.41 % |
| 2026-04 | 0.810 | 0.814 | 0.49 % |
| 2026-05 | 0.836 | 0.840 | 0.45 % |

→ **Negligible (≈ 0.4–0.5 %).** With a ~1 nm laser bandwidth the WV correction averages over the
H₂O line structure, so the exact centre within ±0.2 nm hardly matters. The measurement *confirms*
the manufacturer wavelength is adequate for the CL61; either value gives essentially the same WV
correction.

### 3.2 CL61 — spectral width (FWHM) is the larger lever

| config (λ₀ = 910.74) | T²_wv (500–3000 m) | T²_wv (2–6 km) | β impact |
|---|---|---|---|
| FWHM 0.1 nm (≈monochromatic) | 0.831 | 0.795 | **+4.0 %** |
| FWHM 0.5 nm | 0.857 | 0.818 | +0.9 % |
| **FWHM 1.0 nm (operational)** | **0.864** | **0.826** | **0** |
| FWHM 1.5 nm | 0.873 | 0.836 | −1.0 % |

The assumed **width** moves T²_wv more than the 0.19 nm centre shift: a near-monochromatic
treatment raises the corrected β by ~4 % relative to FWHM 1.0 nm. A narrower line sees less
*average* absorption only if it sits in a micro-window — which leads to the key caveat below.

### 3.3 Key caveat — do **not** treat the CL61 as monochromatic

The measured CL61 line is intrinsically narrow (true width ≤ 0.03 nm), so one might use a
near-monochromatic spectrum. But then T²_wv becomes extremely sensitive to the **±0.10 nm
wavelength-calibration uncertainty** (and to any laser drift), because the line can fall in a
clear micro-window or directly on a strong H₂O line:

| FWHM | T²_wv at λ₀−0.1 | at λ₀ | at λ₀+0.1 | β swing over ±0.1 nm |
|---|---|---|---|---|
| 0.1 nm | 0.982 | 0.831 | 0.771 | **27 %** |
| 1.0 nm | 0.862 | 0.864 | 0.868 | **0.7 %** |

→ A monochromatic CL61 model would make the WV correction swing **≈ 27 %** within the measurement
uncertainty alone. The operational **FWHM = 1.0 nm** is therefore a deliberate, robust choice:
averaging over ~1 nm smooths the line structure so the correction is insensitive to the exact
(uncertain) centre. Using the true narrow line would require pinning λ₀ to < 0.01 nm and tracking
its drift — not warranted for a ~0.5 % gain.

### 3.4 CL31 — broad band self-averages despite the wandering

The CL31 has a large *nominal* spectral uncertainty (≈5–7 nm FWHM, peak wandering 909.0–910.1 nm),
yet the WV correction is **robust**: over the full λ₀ × FWHM range the median T²_wv stays within
**0.873–0.879** (a **< 1 %** β spread):

| config | T²_wv (500–3000 m) |
|---|---|
| 909.7 nm, FWHM 6 (operational) | 0.876 |
| 909.0 nm, FWHM 6 | 0.879 |
| 910.1 nm, FWHM 6 | 0.874 |
| 909.7 nm, FWHM 5 | 0.873 |
| 909.7 nm, FWHM 7 | 0.878 |

→ The broad band integrates over many H₂O lines, so the spectral breadth and the inter-acquisition
wandering **average out**: the CL31's spectral "messiness" does *not* propagate into WV-correction
uncertainty (< 1 %). The breadth that makes a single nominal wavelength ill-defined is exactly what
makes the WV correction insensitive to it.

### 3.5 CHM15k

At 1064.47 nm the CHM15k is **outside** the H₂O absorption band, so no WV correction is applied and
the wavelength configuration is irrelevant to it (it is the WV-free 1064 nm reference).

### 3.6 Summary of wavelength-configuration sensitivity

| Factor | β impact on WV correction | note |
|---|---|---|
| CL61 manufacturer (910.55) vs measured (910.74), FWHM 1.0 | **≈ 0.4–0.5 %** | negligible; manufacturer value adequate |
| CL61 FWHM 0.5 → 1.5 nm | **≈ 1–4 %** | width is the larger lever |
| CL61 treated as monochromatic, ±0.10 nm centre | **≈ 27 %** | **avoid**; FWHM ≈ 1 nm averaging is essential |
| CL31 λ₀ 909.0–910.1 nm × FWHM 5–7 nm | **< 1 %** | broad band self-averages; robust |
| CHM15k (1064 nm) | — | outside band, no WV |

**Bottom line.** With the operational configuration the wavelength *configuration* contributes
**< 1 %** to the WV correction for both 910 nm instruments — far below the WV correction itself
(≈ 16–18 %, see §4) and the other calibration sensitivities (method ≈ 14 pts, Ångström ≈ 14
pts/unit). The Qmini measurements confirm the operational values, and the manufacturer 910.55 nm
would have been equally fine for the CL61. The one thing to avoid is a **monochromatic** CL61
model: at the true narrow linewidth the ±0.10 nm wavelength uncertainty would swing the WV
correction by ~27 %, so the ~1 nm spectral averaging must be kept. For the CL31, the broad emission
band makes the WV correction insensitive to its (poorly defined, wandering) central wavelength.

---

## 4. CL61 Rayleigh-calibration sensitivity to the WV correction (Payerne)

**Question.** Does the WV absorption correction **in the CL61 Rayleigh calibration** matter, and by
how much? The CL61 emits at 910.74 nm, inside the H₂O absorption band, so the high molecular-
reference region used by the Rayleigh calibration is itself attenuated by WV. If that attenuation
is not removed, the retrieved lidar constant is biased low (Wiegner & Gasteiger, 2015).

### 4.1 Method — calibration run twice

The Rayleigh lidar constant `C_L = RCS/β_att` (Wiegner convention) for the Payerne CL61
(`0-20000-0-06610_C`, 46.8137° N 6.9425° E, Mar–May 2026, 30 m grid) was computed **twice**:

- **WV calib** — RCS divided by the two-way WV transmission T²_wv before the lidar-constant fit
  (`apply_wv_correction = 1`, the production setting);
- **no-WV calib** — same fit without that division (`apply_wv_correction = 0`).

**Everything else is identical**: same 30 m CL61 L2, same Kalman filter, same display-side WV
correction and Ångström scaling, same CHM15k (Rayleigh, 1064 nm) reference, same screening and
30 m comparison grid. The two CL61 channels therefore differ **only in the calibration constant**,
isolating the effect of the WV correction *in the calibration*.

### 4.2 Results

**Calibration constant `C_L`** (Kalman median over the period):

| | lidar constant C_L | ratio |
|---|---|---|
| WV calib | 0.595 | — |
| no-WV calib | 0.514 | **0.86** |

The ratio **0.86 ≈ T²_wv(ref)**: the WV-uncorrected constant is low by the two-way WV transmission
at the molecular reference, so the WV correction **raises the constant by ≈ 16 %** (1 / 0.86).

**Validation against the CHM15k (Rayleigh), 500–3000 m AGL:**

| CL61 Rayleigh | rel. bias vs CHM15k | RMSE [Mm⁻¹sr⁻¹] | r | N |
|---|---|---|---|---|
| **WV calib** (production) | **+1.7 %** | 0.060 | 0.987 | 77 522 |
| **no-WV calib** | **+19.3 %** | 0.083 | 0.988 | 77 522 |

*(Figure — Payerne CL61 Rayleigh WV-in-calibration sensitivity, Mar–May 2026 — not archived in
this report tree; the source `figs_paper_validation/sensitivity_payerne_wv.png` was not preserved.
It showed: (a) an example profile, (b) scatter vs CHM15k, (c–e) time-height attenuated backscatter
for CHM15k, CL61 (WV calib) and CL61 (no-WV calib), with the no-WV-calibrated CL61 (e)
systematically brighter/higher than the WV-calibrated one (d); black dots = cloud-base detections.
Reproduce with `sensitivity_payerne_wv.m`, §6.)*

### 4.3 Interpretation

- **The WV correction in the calibration is essential: a ≈ 16–18 % effect** (+1.7 % → +19.3 %, a
  17.6-percentage-point validation swing; 1/T²_wv ≈ 16 % on the constant). Without it the CL61
  looks ≈ +19 % *too high* relative to the 1064 nm reference, whereas the correctly WV-calibrated
  CL61 reads **+1.7 %** — essentially unbiased, the value expected from the clean molecular-vs-
  molecular (CHM Rayleigh ↔ CL61 Rayleigh) comparison.
- **The correlation is unchanged (r ≈ 0.987).** The calibration constant only rescales the profile;
  the range-dependent WV shape is removed by the display-side correction in *both* cases, so r is
  insensitive to the calibration-WV toggle. The WV-in-calibration effect is purely a ≈ 16 %
  **scale** (bias) effect.
- The magnitude (≈ 16–18 %) is consistent with Wiegner & Gasteiger (2015), who report a ~20 %
  backscatter bias at mid-latitudes when WV absorption is ignored at ~910 nm. *(An earlier run with
  the CL61 coordinates left at 0,0 over-estimated this at ≈ 25 %: the Gulf-of-Guinea CAMS column is
  far moister than Payerne — hence the mandatory station-coordinate patch, see §6.)*

**Conclusion.** The water-vapour correction must be applied in the Rayleigh calibration of 910 nm
ALCs; omitting it biases the CL61 lidar constant by ≈ 16 % and the validated attenuated backscatter
accordingly (+1.7 % → +19.3 % vs the 1064 nm reference). The production dataset uses it
(`apply_wv_correction = 1`), and — per §0 — a no-WV degraded mode is rejected outright.

---

## 5. Resolution sensitivity of the WV correction (CAMS grid)

Beyond the *spectral* configuration (§3), the WV correction depends on the **humidity profile
source** — its horizontal resolution, its vertical resolution in the boundary layer, and its
temporal sampling. This section quantifies those, comparing the operational **CAMS 0.4°** (L137
model levels) against the **CAMS 1° fallback**, against **ERA5 0.25°** and hourly IFS, and against
**radiosonde truth** at Payerne. The error metric is ΔC, the resulting error on the cloud-
calibration constant at a 2 km cloud base (β ∝ 1/T²_wv, so a T²_wv error maps directly to a C
error). Ten stations spanning flat lowland to deep Alpine valley are used (2025).

**Headline findings.**
- **Flat/lowland sites are ~1 % or better.** Where the model orography matches the station
  altitude (Payerne 490 m, Montpellier 1 m, Uccle 100 m, Falsterbo 2 m, Köln 50 m), CAMS 0.4°
  reproduces the radiosonde WV correction to within ≈ 1 % on C; the 0.4°→1° coarsening is a small
  additional change.
- **Valley/mountain orography error dominates.** The error scales with the **model–station
  orography mismatch**, not with resolution per se: the coarse grid cell sits at the wrong altitude
  and samples the wrong humidity column. **Aosta** (560 m station, but the CAMS 0.4° cell is
  **+1627 m** too high) is the extreme — a deep Alpine valley the model fills with mountain — and
  shows the largest ΔC of the set.
- **Below-surface handling.** Where the model surface is above the station (valley case), the
  humidity profile is extrapolated downward with a **constant number-density fill** below the model
  surface; this is a deliberate, bounded choice but is the root of the valley error.
- **Temporal sampling (1-hourly vs 3-hourly) is negligible.** Sampling the WV correction 3-hourly
  and interpolating (vs hourly truth) costs **RMS ≈ 0.40 % (winter) / 0.38 % (summer)** on C, with
  P95 < 0.85 % — well below the orographic and spectral terms. The 3-hourly CAMS cadence is
  therefore adequate.
- **Seasonality / IWV.** The correction magnitude grows with integrated water vapour (≈ 6–10 % on
  C at IWV ≈ 5 mm rising to ≈ 25–30 % at IWV ≈ 35–40 mm); the *source error* (CAMS − sonde) stays
  within ≈ ±2 % across the IWV range, i.e. the WV correction is applied robustly even in the moist
  summer boundary layer.

![Coarse-CAMS horizontal-resolution impact on the WV correction — 10 stations (2025). (a) grid-cell altitude (orography) mismatch, CAMS 1° (grey) vs CAMS 0.4° operational (green), cell − station [m]; (b) resulting cloud-calibration error ΔC at 2 km vs independent ERA5, summer; (c) |ΔC| grows with |CAMS 0.4° orography error| — flat sites (blue) cluster < 1 %, mountain/valley sites (red, Aosta/Samedan/Radstadt) rise with the mismatch.](figs_wv_resolution/fig_wv_10stations_impact.png)

![WV transmission by model resolution — flat plateau vs deep Alpine valley (2025-07-15 12 UT). Left: Payerne (490 m, CL61), CAMS 0.4° orography error +289 m — CAMS 1° (dotted), CAMS 0.4° operational (red dashed) and ERA5 0.25° (blue) T²_wv(range) essentially coincide. Right: Aosta (560 m, CL61), CAMS 0.4° orography error +1627 m — the operational CAMS profiles depart markedly from ERA5, the signature of the valley being filled with mountain in the coarse grid.](figs_wv_resolution/fig_wv_10stations_profiles.png)

![Payerne — seasonal WV correction error (CAMS 0.4° operational vs radiosonde truth). (a) ΔC at 2 km cloud base by source (CAMS, hourly IFS) and season (winter/summer) — box spread within ±1–1.5 %; (b) correction magnitude 100·(1−T²) at 2 km vs integrated water vapour IWV, winter (blue) and summer (orange) sondes — rising from ≈ 6 % at 5 mm to ≈ 30 % at 40 mm; (c) source error ΔC (CAMS − sonde) vs IWV — bounded within ≈ ±2 % across the IWV range.](figs_wv_resolution/fig_wv_payerne_seasonal.png)

![Payerne — water-vapour correction vs profile source and season (radiosonde = truth). Top row winter (2025-01-15, IWV = 10 mm), bottom row summer (2025-07-15, IWV = 29 mm). (a,d) humidity profile n_H₂O — radiosonde (black), IFS at point hourly (blue), CAMS 0.4° 3-hourly (red dashed); (b,e) resulting T²_wv(height); (c,f) cloud-calibration error ΔC vs cloud-base height for CAMS − sonde (red) and IFS − sonde (blue) — errors within ≈ ±1–2 % up to a 4 km cloud base at this flat site.](figs_wv_resolution/fig_wv_payerne_profiles.png)

![Payerne — benefit of 1-hourly vs 3-hourly WV sampling (hourly IFS as truth, error at 2 km cloud base). (a) diurnal T²_wv (summer 07-15): hourly IFS truth vs 3-hourly sampled+interpolated track closely; (b) winter ΔC histogram from 3-hourly sampling: RMS = 0.40 %, P95 = 0.79 %; (c) summer ΔC histogram: RMS = 0.38 %, P95 = 0.84 %.](figs_wv_resolution/fig_wv_payerne_temporal.png)

**Operational consequence.** The operational monthly **0.4°** CAMS (with the per-month **1°
fallback**) is adequate for flat and lowland stations. For **deep-valley Alpine sites** the coarse-
grid orography error is the leading WV-correction uncertainty (metres of altitude mismatch, not
resolution as such); this is a known limitation flagged for those stations, not a reason to reject
the correction. The 3-hourly cadence and the manufacturer-vs-measured wavelength choice are both
negligible by comparison.

---

## 6. Implementation notes and reproduction

- **Operational code path (current layout).** The spectral WAPL-style correction lives in
  `calibration/water_vapor_correction/` — `water_vapor.py` (the Rayleigh path:
  `two_way_wv_transmission`, `in_water_vapor_band`, `laser_spectrum_for`, `LASER_SPECTRUM`) and
  `cloud_water_vapor.py` (the liquid-cloud path). Both are driven by CAMS specific humidity q on
  the ECMWF L137 model levels and the HITRAN/MT-CKD LUT `abs_cross_647_full_levels_1000.nc` (a
  bundled 910 nm sub-band, `data/abs_cross_wv_910nm.nc`, is used out-of-the-box when no external
  LUT path is configured). The correction is a **faithful Python port** of the validated MATLAB
  routines `wv_t2eff.m`, `compute_wv_transmission.m`,
  `get_water_vapor_number_concentration_from_RH.m` — validated to ≤ 0.4 % (see
  `tests/WATER_VAPOR_AUDIT.md`).
  *Superseded note:* older reports referenced the module as `rayleigh_calibration/water_vapor.py`;
  it now lives under `calibration/water_vapor_correction/`.
- **Read-once integration.** Each instrument-day loads L1 once + CAMS once, coarsens to the shared
  30 s × 10 m working grid, and computes T²_wv once for reuse across the Rayleigh and cloud
  branches (`calibration/io/instrument_day.py`; `wv_source='cams'` default).
- **Enforcement (verified in code).** The correction is gated on the laser wavelength (900–920 nm),
  not the instrument name. A 910 nm night with no CAMS → `flag = -4`; station outside the CAMS
  domain → `flag = -10`; unusable WV profile or transmission grid mismatch → `flag = -4`. In the
  cloud path a failed / all-ones / empty transmission raises a hard `RuntimeError` ("no fallback:
  this period is NOT calibrated"). 1064 nm and 532 nm never enter the WV block.
- **Reproduce the wavelength-configuration study (§3):** `wv_wavelength_sensitivity.py` — builds
  T²_wv with the operational machinery (`water_vapor.two_way_wv_transmission`, the HITRAN LUT, CAMS
  humidity at Payerne) for each (λ₀, FWHM) and renders the figure.
- **Reproduce the CL61 Rayleigh sensitivity (§4):** `sensitivity_payerne_wv.m`. The WV-calibrated
  channel uses the correct-coordinates CSV (`…\fullcal_stdmol_check`); the WV-off channel uses
  `…\fullcal_all_noWV` (no CAMS → lat/lon-independent). Both use the Payerne CL61 L2 **patched to
  the real coordinates** (the raw files report 0,0). The WV-off CSV is a **diagnostic only**; the
  operational chain keeps `apply_wv_correction = 1`.
- **Resolution study (§5):** figures in `figs_wv_resolution/` compare CAMS 0.4° / CAMS 1° / ERA5
  0.25° / hourly IFS against radiosonde truth at Payerne and across 10 stations (2025); the
  metric ΔC is the cloud-calibration error at a 2 km cloud base. The below-surface fill (constant
  n_H₂O below the model surface) and the 3-hourly→interpolated sampling are the two source
  choices exercised.

---

## Sources

- [Wiegner & Gasteiger (2015), AMT 8, 3971–3984](https://doi.org/10.5194/amt-8-3971-2015) — WV correction (WAPL); HITRAN+MT-CKD LUT; averaging over the laser spectrum; CL51 Δλ ≈ 3.4 nm; input relevance order.
- [Wiegner et al. (2019, CeiLinEx2015), AMT 12, 471–490](https://doi.org/10.5194/amt-12-471-2019) — validation of W&G 2015; CL31/CL51 910±10 nm; CS135 912±3.5 nm; dT_w,eff/dλ > 0.02 nm⁻¹.
- [Hopkin et al. (2019), AMT 12, 4131–4147](https://doi.org/10.5194/amt-12-4131-2019) — E-PROFILE operational reference; O'Connor liquid-cloud calib + empirical `Twv = 1 − 0.17·IWV^0.52`; agrees with WAPL to 2 %.
- [Kotthaus et al. (2016), AMT 9, 3769–3791](https://doi.org/10.5194/amt-9-3769-2016) — CL31 processing; acknowledges but neglects WV.
- [O'Connor et al. (2004), JAOT 21, 777](https://journals.ametsoc.org/view/journals/atot/21/5/1520-0426_2004_021_0777_atfaoc_2_0_co_2.xml) — liquid-cloud calibration reference.
- [Water vapour correction with a 910 nm CL51 (Remote Sens. 17, 2013, 2025)](https://doi.org/10.3390/rs17122013) — independent ">20 % if WV ignored" (magnitude only; not full WAPL).
- [Le et al. (2026), Vaisala CL61 performance, EGUsphere](https://doi.org/10.5194/egusphere-2025-6331) (+ supplement) — CL61 910.55 nm; 160 ns pulse; laser-temp cycling; cloud C.
- [Looschelders et al. (2025), Meteorol. Appl. 32, e70088](https://doi.org/10.1002/met.70088) — CL61 inter-instrument; 910.55 nm; ±5 % calibration.
- [Laffineur et al. (2026, CONIOPOL), EGUsphere](https://doi.org/10.5194/egusphere-2026-948) — CL61 910.55 nm, Uccle.
- [CloudnetPy](https://github.com/actris-cloudnet/cloudnetpy) — open, but no WAPL/WV correction (scalar calibration factor only).
- [HITRAN / MT-CKD](https://hitran.org/mtckd/) — spectroscopy basis.
- Vaisala datasheets: CL31 B210415EN, CL51 B210861EN, CL61 User Guide M212475EN-E; ARM VCEIL handbook.
- Verification: 17 sources fetched, 79 claims extracted, 25 verified, 23 confirmed / 2 refuted (literature review, 2026-06-16).
