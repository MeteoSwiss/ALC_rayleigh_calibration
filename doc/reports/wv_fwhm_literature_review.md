# Laser emission FWHM for the water-vapour correction — literature review (CL31, CL51, CL61)

*M. Hervo (MeteoSwiss), 2026-06-16. Supports the WV-correction wavelength-configuration study
([wv_wavelength_sensitivity.md](wv_wavelength_sensitivity.md)) and the β_att validation
([attbsc_validation_technical.md](attbsc_validation_technical.md)).*

## Why the FWHM matters

The water-vapour (WV) correction weights the (sharply structured) H₂O absorption cross-section by
the laser **emission spectrum**, assumed Gaussian with central wavelength λ₀ and full width at half
maximum Δλ (FWHM). Wiegner & Gasteiger (2015) introduced this and rank the required inputs **by
decreasing relevance: (1) the water-vapour profile, (2) the central wavelength λ₀, (3) the spectral
width Δλ.** The spectral width is the *least* sensitive of the three — but it is still needed, and
it is the worst-documented, especially for the CL61.

Three different "widths" appear in the datasheets and must not be confused:
- **emission FWHM (Δλ, spectral)** — the quantity used by the WV correction (this review);
- **receiver optical-filter bandwidth** — 36 nm for CL31/CL51 (Wiegner et al. 2019); irrelevant to the weighting;
- **pulse duration (temporal FWHM)** — e.g. 160 ns for the CL61; not a spectral width.

## Summary — what the literature/manufacturer give

| Instrument | central λ₀ | **emission FWHM Δλ** | source(s) | our config | Payerne Qmini (measured, 2026-06-02) |
|---|---|---|---|---|---|
| **CL31** | 905 ± 10 nm (older) → 910 ± 10 nm (newer), 25 °C | **≈ 4 nm** (typical); 2021 datasheet omits it | Vaisala CL31 datasheet/ARM VCEIL handbook; Kotthaus et al. 2016; Wiegner et al. 2019 | (909.7, **6.0**) | 909.7, **5–7 nm**, peak wandering 909.0–910.1 nm |
| **CL51** | 910 ± 10 nm @25 °C, drift 0.27 nm K⁻¹ | **3.4 nm** (Vaisala); 3.5 nm used; 1–4 nm explored | **Wiegner & Gasteiger 2015**; Wiegner et al. 2019; Cordoba-Jabonero / MDPI 2025 | (910.0, **3.4**) | — (not measured) |
| **CL61** | **910.55 nm** | **not documented** — only λ₀ is published | Vaisala CL61 User Guide M212475EN-E; Le et al. 2026; Looschelders et al. 2025; Laffineur et al. 2026 | (910.74, **1.0**) | **910.74 ± 0.10**, single narrow line (true ≤ 0.03 nm; apparent ≈1.5 nm spectrometer-limited) |
| CS135 (Campbell, context) | 912 nm (stable) | ±3.5 nm | Wiegner et al. 2019 | — | — |
| CHM15k (context) | 1064 nm | n/a — outside the H₂O band, no WV | — | (1064.47, 0.5) | 1064.47 ± 0.10 |

## 1. CL51 — the only well-documented case (Δλ ≈ 3.4 nm)

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
**Our operational value `CL51 = (910.0, 3.4)` is exactly the Vaisala/W&G figure** — no change needed.

## 2. CL31 — broad multimode diode, ≈ 4 nm, poorly pinned

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

## 3. CL61 — only the centre (910.55 nm) is published; the FWHM is undocumented

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

1. The CL61 line is much narrower than the CL31/CL51 multimode diodes (3–7 nm). It is **not**
   safe to inherit their ~3–4 nm widths, but it is **also not** safe to treat the CL61 as
   monochromatic (see §5): the measured ±0.10 nm centre uncertainty would then dominate. Our
   operational `CL61 = (910.74, 1.0)` keeps a deliberate ~1 nm averaging for robustness.
2. We measure **910.74 nm**, 0.19 nm above the Vaisala spec **910.55 nm**. The literature uniformly
   uses 910.55; at a ~1 nm bandwidth the difference is < 0.5 % on the WV correction
   ([wv_wavelength_sensitivity.md](wv_wavelength_sensitivity.md)), so either is adequate, but the
   measured value is the more defensible one to report.

Relatedly, Le et al. (2026, Fig. S2) show the CL61 **laser temperature cycling** (≈ 19–22 °C,
~127 s period). At the CL51 drift of 0.27 nm K⁻¹ such swings imply sub-nm centre drift, reinforcing
that a finite (~1 nm) spectral averaging — rather than a fixed monochromatic line — is the robust
choice for the CL61 WV correction.

## 4. Context — other instruments

- **Campbell CS135**: λ = 912 nm (stable), spectral width ±3.5 nm (Wiegner et al. 2019).
- **Lufft CHM15k/CHM8k**: 1064 nm — outside the H₂O band, **no WV correction**; our measured
  1064.47 nm matches the manufacturer.

## 5. How the literature treats the FWHM sensitivity (and how it matches our study)

- **Relevance order (Wiegner & Gasteiger 2015):** water-vapour profile ≫ central wavelength λ₀ >
  spectral width Δλ. *"The variability with Δλ depends on λ₀ but is in most cases comparably
  small."* Absorption rises by ~×2 from λ₀ = 905→908 nm but varies < 10 % from 908→918 nm.
- **Central-wavelength sensitivity (Wiegner et al. 2019):** *"the transmission is much more
  sensitive to errors of the assumed wavelength λ_on than to errors of the water-vapour content …
  dT_w,eff/dλ > 0.02 nm⁻¹"* (> 2 % per nm).
- **Recurring caveat:** several papers note that *"the unknown emission spectrum of the diode
  laser … can introduce significant errors"* — i.e. the spectrum is treated as an unmeasured
  manufacturer input.

This literature consensus **matches our own sensitivity study**
([wv_wavelength_sensitivity.md](wv_wavelength_sensitivity.md)) exactly: at fixed FWHM the CL61
910.55-vs-910.74 nm choice is ≈ 0.4–0.5 %; the FWHM itself is a ≈ 1–4 % lever; and the CL31 broad
band self-averages to < 1 % despite its wandering. The one regime the older literature does not
emphasise — because it predates the narrow-line CL61 — is that a **near-monochromatic** treatment
makes T²_wv swing ≈ 27 % within the ±0.10 nm centre uncertainty, which is why the ~1 nm averaging
must be kept.

![WV correction vs laser-wavelength configuration at Payerne (companion sensitivity study): (a) H₂O absorption band with the CL31/CL61 laser spectra overlaid, (b) T²_wv profiles, (c) median T²_wv by configuration.](figs_paper_validation/wv_wavelength_sensitivity.png)

## 6. What is useful for our paper

- **CL51 = 3.4 nm** is the only manufacturer-stated Vaisala ALC emission FWHM (Wiegner & Gasteiger
  2015) — cite it; our config matches.
- **The CL61 emission linewidth is undocumented** in Vaisala specs and in all CL61 papers (Le 2026,
  Looschelders 2025, Laffineur 2026) — they give only 910.55 nm, and the "FWHM" in the spec is the
  160 ns pulse duration. **Our Qmini measurement appears to be the first reported CL61 emission
  spectrum** — worth stating as a contribution.
- The **measured CL61 centre (910.74 nm) differs from the spec (910.55 nm) by 0.19 nm**; negligible
  for the correction but the measured value is preferable.
- The literature's **relevance order (WV ≫ λ₀ > Δλ)** and the **dT/dλ > 0.02 nm⁻¹** sensitivity
  support our framing; our monochromatic-instability result extends it to the narrow-line CL61.
- The **CL61 laser-temperature cycling** (Le 2026, Fig. S2) is a citable physical argument for using
  a finite spectral-averaging width rather than a fixed monochromatic line.
- **Cross-checks for the cloud calibration:** CL61 O'Connor/Hopkin cloud factors C ≈ 1.0–1.4 across
  the FMI sites (Le 2026, Figs S4–S7) and ≈ ±5 % inter-instrument spread (Looschelders 2025)
  corroborate our Payerne CL61 cloud C ≈ 0.99.

---
### Sources
- [Wiegner & Gasteiger (2015), AMT 8, 3971–3984](https://doi.org/10.5194/amt-8-3971-2015) — WV correction; CL51 Δλ ≈ 3.4 nm; relevance order.
- [Wiegner et al. (2019, CeiLinEx2015), AMT 12, 471–490](https://doi.org/10.5194/amt-12-471-2019) — CL31/CL51 910±10 nm, CS135 912±3.5 nm, dT/dλ.
- [Kotthaus et al. (2016), AMT 9, 3769–3791](https://doi.org/10.5194/amt-9-3769-2016) — CL31 processing.
- [Water vapour correction with a 910 nm CL51 (Remote Sens. 17, 2013, 2025)](https://doi.org/10.3390/rs17122013).
- [Le et al. (2026), Vaisala CL61 performance, EGUsphere](https://doi.org/10.5194/egusphere-2025-6331) (+ supplement) — CL61 910.55 nm; 160 ns pulse; laser-temp cycling; cloud C.
- [Looschelders et al. (2025), Meteorol. Appl. 32, e70088](https://doi.org/10.1002/met.70088) — CL61 inter-instrument; 910.55 nm; ±5 % calibration.
- [Laffineur et al. (2026, CONIOPOL), EGUsphere](https://doi.org/10.5194/egusphere-2026-948) — CL61 910.55 nm, Uccle.
- Vaisala datasheets: CL31 B210415EN, CL51 B210861EN, CL61 User Guide M212475EN-E; ARM VCEIL handbook.
