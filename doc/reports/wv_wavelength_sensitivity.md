# Sensitivity of the water-vapour correction to the laser-wavelength configuration

*M. Hervo (MeteoSwiss), 2026-06-16. Companion to the β_att validation
([attbsc_validation_technical.md](attbsc_validation_technical.md)) and the calibration
sensitivity study ([payerne_cl61_calibration_sensitivity.md](payerne_cl61_calibration_sensitivity.md)).*

## Question

The water-vapour (WV) correction weights the H₂O absorption cross-section by the **laser
emission spectrum** (central wavelength λ₀ and spectral FWHM) to get the two-way transmission
T²_wv. Across 905–915 nm the H₂O cross-section varies by orders of magnitude over fractions of a
nanometre (Fig. a), so the assumed λ₀ and FWHM matter. How much does the *wavelength
configuration* change T²_wv — i.e. the manufacturer value vs the measured one, the spectral
width, and (for the CL31) the spectral breadth and inter-acquisition wandering?

**Emission wavelengths** (Qmini fibre spectrometer, Payerne, 2 June 2026; scale verified on the
O₂ A-band to 0.1 nm): CHM15k 1064.47 ± 0.10 nm; CL61 910.74 ± 0.10 nm (single narrow line,
apparent ≈1.5 nm FWHM spectrometer-limited, true scatter ≤0.03 nm); CL31 ≈5–7 nm FWHM centred
≈909.7 nm, peak wandering 909.0–910.1 nm between acquisitions. Manufacturer CL61 value:
**910.55 nm**. The operational WV code uses the **measured** values
(`CL61 = (910.74, 1.0)`, `CL31 = (909.7, 6.0)`).

Metric: median T²_wv over the **500–3000 m** comparison band (drives the validation bias, since
β_corr = β / T²_wv) and over the **2–6 km** Rayleigh window (drives the lidar constant). "β impact"
is the resulting change in corrected backscatter (β ∝ 1/T²_wv). Profiles from CAMS at Payerne.

![Water-vapour correction vs laser-wavelength configuration at Payerne. (a) H₂O absorption cross-section across 905–915 nm at ~3 km (grey, log scale) with the normalised laser spectra overlaid — CL61 manufacturer 910.55 nm (red), CL61 measured 910.74 nm (blue), CL31 909.7 nm/6 nm FWHM (green); the CL61 narrow lines sample individual absorption features while the CL31 broad band averages over many. (b) Two-way WV transmission T²_wv vs range for the three configurations, with the CL61 FWHM 0.1–1.5 nm envelope shaded (blue); the 500–3000 m comparison band is greyed. (c) Median T²_wv over 500–3000 m for every configuration.](figs_paper_validation/wv_wavelength_sensitivity.png)

## 1. CL61 — manufacturer 910.55 nm vs measured 910.74 nm

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

## 2. CL61 — spectral width (FWHM) is the larger lever

| config (λ₀ = 910.74) | T²_wv (500–3000 m) | T²_wv (2–6 km) | β impact |
|---|---|---|---|
| FWHM 0.1 nm (≈monochromatic) | 0.831 | 0.795 | **+4.0 %** |
| FWHM 0.5 nm | 0.857 | 0.818 | +0.9 % |
| **FWHM 1.0 nm (operational)** | **0.864** | **0.826** | **0** |
| FWHM 1.5 nm | 0.873 | 0.836 | −1.0 % |

The assumed **width** moves T²_wv more than the 0.19 nm centre shift: a near-monochromatic
treatment raises the corrected β by ~4 % relative to FWHM 1.0 nm. A narrower line sees less
*average* absorption only if it sits in a micro-window — which leads to the key caveat below.

## 3. Key caveat — do **not** treat the CL61 as monochromatic

The measured CL61 line is intrinsically narrow (true width ≤ 0.03 nm), so one might use a
near-monochromatic spectrum. But then T²_wv becomes extremely sensitive to the **±0.10 nm
wavelength-calibration uncertainty** (and to any laser drift), because the line can fall in a
clear micro-window or directly on a strong H₂O line:

| FWHM | T²_wv at λ₀−0.1 | at λ₀ | at λ₀+0.1 | β swing over ±0.1 nm |
|---|---|---|---|---|
| 0.1 nm | 0.982 | 0.831 | 0.771 | **27 %** |
| 1.0 nm | 0.862 | 0.864 | 0.868 | **0.7 %** |

→ A monochromatic CL61 model would make the WV correction swing **≈ 27 %** within the
measurement uncertainty alone. The operational **FWHM = 1.0 nm** is therefore a deliberate,
robust choice: averaging over ~1 nm smooths the line structure so the correction is insensitive
to the exact (uncertain) centre. Using the true narrow line would require pinning λ₀ to
< 0.01 nm and tracking its drift — not warranted for a ~0.5 % gain.

## 4. CL31 — broad band self-averages despite the wandering

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

→ The broad band integrates over many H₂O lines, so the spectral breadth and the
inter-acquisition wandering **average out**: the CL31's spectral "messiness" does *not* propagate
into WV-correction uncertainty (< 1 %). The breadth that makes a single nominal wavelength
ill-defined is exactly what makes the WV correction insensitive to it.

## 5. CHM15k

At 1064.47 nm the CHM15k is **outside** the H₂O absorption band, so no WV correction is applied
and the wavelength configuration is irrelevant to it (it is the WV-free 1064 nm reference).

## Summary

| Factor | β impact on WV correction | note |
|---|---|---|
| CL61 manufacturer (910.55) vs measured (910.74), FWHM 1.0 | **≈ 0.4–0.5 %** | negligible; manufacturer value adequate |
| CL61 FWHM 0.5 → 1.5 nm | **≈ 1–4 %** | width is the larger lever |
| CL61 treated as monochromatic, ±0.10 nm centre | **≈ 27 %** | **avoid**; FWHM ≈ 1 nm averaging is essential |
| CL31 λ₀ 909.0–910.1 nm × FWHM 5–7 nm | **< 1 %** | broad band self-averages; robust |
| CHM15k (1064 nm) | — | outside band, no WV |

**Bottom line.** With the operational configuration the wavelength *configuration* contributes
**< 1 %** to the WV correction for both 910 nm instruments — far below the WV correction itself
(≈ 16–18 %) and the other calibration sensitivities (method ≈ 14 pts, Ångström ≈ 14 pts/unit). The
Qmini measurements confirm the operational values, and the manufacturer 910.55 nm would have been
equally fine for the CL61. The one thing to avoid is a **monochromatic** CL61 model: at the true
narrow linewidth the ±0.10 nm wavelength uncertainty would swing the WV correction by ~27 %, so the
~1 nm spectral averaging must be kept. For the CL31, the broad emission band makes the WV
correction insensitive to its (poorly defined, wandering) central wavelength.

---
*Reproduce:* `ALC_rayleigh_calibration/wv_wavelength_sensitivity.py` — builds T²_wv with the
operational machinery (`water_vapor.two_way_wv_transmission`, HITRAN LUT
`abs_cross_647_full_levels_1000.nc`, CAMS humidity at Payerne) for each (λ₀, FWHM) and renders the
figure. Operational laser spectrum (`water_vapor.py LASER_SPECTRUM`): CL31 (909.7, 6.0),
CL51 (910.0, 3.4), CL61 (910.74, 1.0), CHM15k (1064.47, 0.5).
