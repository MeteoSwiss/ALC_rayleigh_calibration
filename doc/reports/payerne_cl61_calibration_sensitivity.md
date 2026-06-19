# Payerne CL61 calibration sensitivity study

*M. Hervo (MeteoSwiss), 2026-06-16. Companion to the β_att validation
([attbsc_validation_technical.md](attbsc_validation_technical.md)). Period Mar–May 2026,
CL61 `0-20000-0-06610_C` at Payerne (**46.8137° N, 6.9425° E, 491 m ASL**) on the E-PROFILE 30 m
grid, reference = colocated CHM15k (Rayleigh, 1064 nm), comparison band 500–3000 m AGL, hourly
screened. The station coordinates matter: they set the CAMS grid point for the water-vapour
correction — see the CL61 lat/lon fix at the end.*

How robust is the CL61 attenuated-backscatter calibration to the main processing choices?
Four base calibrations are computed, then four sensitivity axes are quantified.

## 1. Four CL61 calibrations vs CHM15k (Rayleigh)

Each calibration is applied identically downstream (same 30 m L2, Kalman filter, display-side
WV correction, Ångström scaling, screening); they differ **only in the calibration constant**.

| Calibration | corr. factor | rel. bias vs CHM15k | r |
|---|---|---|---|
| **cloud + WV** | 0.987 | **+16.0 %** | 0.987 |
| cloud − WV | 1.174 | +34.6 % | 0.986 |
| **Rayleigh + WV** | 0.833 | **+1.7 %** | 0.987 |
| Rayleigh − WV | 0.971 | +19.3 % | 0.988 |

![Payerne CL61 — four calibrations vs CHM15k (Mar–May 2026): (a) median profile, (b) scatter vs CHM15k, (c–g) time-height attenuated backscatter for CHM15k and the four CL61 calibrations. Black dots = cloud base.](figs_paper_validation/sensitivity_payerne_4calib.png)

- **Calibration method (cloud vs Rayleigh), WV applied:** +16.0 % vs +1.7 % → a **≈ 14-point
  spread**. The liquid-cloud method (O'Connor S = 18.8 sr) gives a higher β_att than the
  molecular-Rayleigh method; the CHM15k(Rayleigh)–CL61(Rayleigh) pair (**+1.7 %**) is the cleanest
  absolute check (both molecular, reference WV-free) — near-zero, as expected for two
  molecular-calibrated instruments. This is the dominant sensitivity.
- All four correlate with the CHM15k at **r ≈ 0.987** — the choices shift the *scale*, not the shape.
- *(All values use the corrected Payerne coordinates. With the earlier lat/lon = 0 the WV read
  sampled the Gulf of Guinea; the Rayleigh + WV bias was then −8.5 % and the cloud + WV +15.0 % —
  see the lat/lon fix note at the end.)*

## 2. Water-vapour correction (in the calibration)

From the table above, toggling WV **in the calibration**:

| Method | WV-off → WV-on | shift |
|---|---|---|
| Cloud | +34.6 % → +16.0 % | **−18.6 pts** |
| Rayleigh | +19.3 % → +1.7 % | **−17.6 pts** |

WV is **essential for both** methods (≈ 18 pts; consistent with Wiegner & Gasteiger 2015,
~20 %). Note the opposite sign of the constant change: WV **lowers** the cloud constant (it raises
the below-cloud integrated backscatter → lower apparent S) but **raises** the Rayleigh constant
(it raises the RCS at the molecular reference). The net bias moves the same way (down) for both.

## 3. Ångström exponent (wavelength correction 910 → 1064 nm)

The 910 nm β is scaled to 1064 nm by `β /= (λ/λ_target)^(−Å)`, i.e. **× (910.74/1064)^Å =
× 0.856^Å** — a uniform factor, so the bias scales **exactly**: bias(Å) = (1+bias₁)·0.856^(Å−1) − 1.
Effect ≈ **−14 % per unit Å**.

| Å | factor vs Å=1 | CL61 cloud (WV) | CL61 Rayleigh (WV) |
|---|---|---|---|
| 0.0 | ×1.168 | +35.5 % | +18.8 % |
| 0.5 | ×1.081 | +25.4 % | +9.9 % |
| **1.0** (used) | ×1.000 | **+16.0 %** | **+1.7 %** |
| 1.5 | ×0.925 | +7.3 % | −5.9 % |
| 2.0 | ×0.856 | −0.7 % | −12.9 % |

This is a **large** sensitivity: over a plausible Å = 0.5–1.5 the bias moves by ~±9 pts. Å = 1
(Haarig et al. 2025, high-RH boundary layer) is the adopted value; the result is sensitive to it,
so Å should be stated explicitly in the paper.

## 4. Lidar ratio

- **Cloud — droplet lidar ratio S** (O'Connor): the constant is `S_apparent / S_theoretical`, so
  **C ∝ 1/S** (exact). bias(S) = (1+bias₁)·18.8/S − 1.

  | S (sr) | CL61 cloud (WV) bias |
  |---|---|
  | 15 | +45.4 % |
  | 18.0 | +21.2 % |
  | **18.8** (used, ±0.8) | **+16.0 %** (±~5 pts) |
  | 20 | +9.0 % |
  | 22 | −0.9 % |

  Within the Hopkin et al. (2019) uncertainty (18.8 ± 0.8 sr) the cloud bias is ± ~5 pts; over a
  wider plausible S it is large. The cloud calibration's absolute level is tied to the assumed S.

- **Rayleigh — aerosol lidar ratio** (Klett): the calibration already perturbs the aerosol LR over
  52 ± 20 sr (plus altitude-shift); the resulting **constant uncertainty is only ± 5.1 %**, because
  the molecular reference region is nearly aerosol-free. The Rayleigh calibration is **robust** to
  the lidar ratio.

## 5. Cloud vs Rayleigh, multiple scattering and the field of view

The dominant sensitivity (§1) is the calibration *method*: the **same** CL61 reads **+14 %** higher
under the liquid-cloud (O'Connor) calibration than under the molecular-Rayleigh one (cloud +16.0 %
vs Rayleigh +1.7 % against the CHM15k). The molecular-Rayleigh CL61 agrees with the colocated CHM15k
(1064 nm, Rayleigh) to **+1.7 %** — two *independent* molecular calibrations — so the molecular scale
is the trustworthy anchor and the **cloud method is the outlier**. This section tests whether the
multiple-scattering (MS) correction explains it, characterises the CL61 MS **two independent ways**,
and keeps the CHM as the final arbiter.

### 5.1 A pure calibration-constant offset
The cloud/Rayleigh β ratio is **range-independent** — 1.143 at *every* level 500–3000 m (Mar–May
2026) — so the difference is entirely in the calibration constant, not a range-dependent artefact
(overlap, a range-dependent MS profile, WV and Ångström all cancel between two channels of the same
instrument).

### 5.2 How η is defined and how the CL31/CL51 were characterised (literature)
The liquid-cloud calibration (O'Connor et al. 2004; Hopkin et al. 2019) forces the integrated
attenuated backscatter through a fully attenuating liquid cloud to its theoretical value
**B = ∫β dz = 1/(2·η·S)**, with the droplet lidar ratio **S = 18.8 ± 0.8 sr** (essentially constant
905–1064 nm; target B = 0.0266 sr⁻¹). **η is the multiple-scattering factor** — the share of
forward-scattered photons recaptured by the receiver — which **depends on the laser beam divergence,
the receiver field of view (FOV) and altitude**, computed per range gate with the fast model of
**Hogan (2006, Appl. Opt. 45, 5984–5992)**. Hopkin et al. (2019) characterised it this way for the
Vaisala CL31/CL51 and report **η ≈ 0.7–0.85** for 905–1064 nm liquid clouds (smaller higher in the
cloud). Our `liquid_cloud_calibration.m` uses exactly these per-gate values (η ≈ 0.76–0.83).
Wiegner et al. (2014) neglect MS (aerosol focus); O'Connor's original instrument was a CT75K
(divergence 0.66 mrad, FOV 0.75 mrad).

### 5.3 The receiver FOV — the decisive number (this corrects an earlier hypothesis)
Multiple scattering is set by the **receiver FOV**, not the wavelength. The half-angle receiver FOVs:

| instrument | receiver FOV (half-angle) | source |
|---|---|---|
| CHM15k | 0.23 mrad | Wiegner et al. 2014, Table 1 |
| **CL61** | **0.56 mrad** | Vaisala spec M212475EN-E |
| CL51 | 0.56 mrad | Wiegner et al. 2014, Table 1 |
| CL31 | 0.83 mrad | Wiegner et al. 2014, Table 1 |

**The CL61 FOV (±0.56 mrad) is identical to the CL51's and wider than the CHM15k — it is *not* a
narrow lidar-class FOV.** (Its laser *beam divergence*, ±0.2 × 0.35 mrad, is narrow, but the
*receiver FOV*, which controls the MS capture, is CL51-class.) Its η is therefore expected in the
**same 0.7–0.85 regime** as the CL51, and **borrowing the CL51 η table for the CL61 is appropriate,
not an over-correction.** Le et al. (2026), who compute the CL61 MS from the CL61's *own* divergence
and FOV for 8–20 µm droplets, land at the same operating point (their CL61 cloud C ≈ 1.0–1.4).
*(This revises an earlier hypothesis in this note that the CL61 FOV was much narrower; the Vaisala
specification shows it equals the CL51's.)*

### 5.4 Independent characterisation of the CL61 multiple scattering
Two independent lines, **neither using the CHM**:

1. **FOV / Hogan model.** With FOV 0.56 mrad (= CL51) the Hogan-2006 per-gate η is the CL51 profile,
   **η ≈ 0.79** near 1–2 km cloud base — i.e. ≈21 % MS, *not* ≈10 %.
2. **In-cloud depolarisation (CL61-unique, data-driven).** Liquid droplets are spherical, so the
   single-scattering linear depolarisation δ = 0; any in-cloud δ is purely multiple scattering. Over
   **2273 fully-attenuating liquid clouds** at Payerne the CL61 δ rises from ≈ 0.05 at cloud base to
   a **peak ≈ 0.11 at ~120 m** depth, then falls as the signal attenuates — a clear but **moderate**
   MS signature, consistent with CL51-class scattering (η ≈ 0.8) and **inconsistent with η ≈ 1**
   (which would leave δ ≈ 0).

![Payerne CL61 in-cloud linear depolarisation vs depth above cloud base (March 2026, 2273 fully-attenuating liquid clouds). δ grows from ≈0.05 at cloud base to ≈0.11 at ~120 m — the multiple-scattering signature. Its moderate magnitude is consistent with CL51-class MS (η≈0.8), not with η≈1 (which would give δ≈0). Independent of the CHM and of the calibration constant.](figs_paper_validation/cl61_incloud_depol.png)

Both independent estimates give **η_CL61 ≈ 0.8 (CL51-like)** — so the code's MS factor is appropriate
and **multiple scattering does *not* explain the +14 %.**

### 5.5 So what is the +14 %? — and the CHM as final validation
With η ≈ 0.8 correct and S = 18.8 sr established, the cloud method finds the CL61 internally
self-consistent (apparent S ≈ 18.5 sr → C ≈ 0.99), yet the molecular scale says the CL61 should be
~14 % lower (C ≈ 0.87). Closing the gap would require **either** η ≈ 0.90 (but the FOV *and* the
depolarisation both say ≈ 0.79) **or** an effective droplet S ≈ 21.5 sr (but 18.8 ± 0.8 sr is well
established). **Neither independent line supports the value needed**, so the **+14 % is a genuine
cloud-vs-Rayleigh inter-method discrepancy, not a fixable CL61 η error.** It sits at the upper end of
the cloud-calibration uncertainty flagged by O'Connor (2004) (~5–10 % from MS/S) once real-cloud S
variability and the residual MS/S degeneracy are folded in; a small shared molecular-calibration
offset cannot be excluded.

**Final validation = CHM15k.** The CL61 Rayleigh matches the colocated CHM15k molecular calibration
to **+1.7 %** (two wavelengths, two instruments), so the molecular scale is the reliable absolute
reference. **Recommendation:** adopt the molecular-Rayleigh calibration for the CL61 absolute scale
and treat the liquid-cloud constant as carrying the larger (~10–15 %) systematic. The +14 % is an
open inter-method difference (candidate causes: droplet-S variability, residual MS/S degeneracy),
**not** a CL61 FOV/η problem.

*Sources: O'Connor et al. (2004), J. Atmos. Oceanic Technol. 21, 777; Hopkin et al. (2019), AMT 12,
4131; Hogan (2006), Appl. Opt. 45, 5984; Wiegner et al. (2014), AMT 7, 1979 (Table 1 FOVs);
Vaisala CL61 spec M212475EN-E (FOV); Le et al. (2026), EGUsphere, doi:10.5194/egusphere-2025-6331.
Reproduce: `cl61_incloud_depol.py` (depolarisation), `cl61_cloud_vs_rayleigh.m` (offset).*

### 5.6 Is it saturation or aerosol contamination? — diagnostics (no)

Two further artefact candidates for the cloud-high offset were tested with dedicated
diagnostics (à la the *Long-term CL61* supplement, Le et al. 2026, and Hopkin et al. 2019).

![Payerne CL61 cloud calibration diagnostics](figs_paper_validation/cl61_cloud_diagnostic.png)

*(a) Apparent droplet lidar ratio S = C·18.8 vs cloud-base height (2-D histogram) against the
theoretical 18.8 sr; (b) its distribution — tightly peaked at ≈ 18.5 sr, so the cloud method is
internally self-consistent; (c) raw per-profile peak β_att — a smooth, bimodal (clear-air ~1e-5 /
in-cloud ~3–5e-4) distribution whose maximum (6.8e-4) sits well above p99.9 (5.4e-4): **no ceiling
or pile-up**.*

- **Detector saturation — ruled out.** Panel (c) shows no hard clipping ceiling, and the apparent S
  correlates **negatively** with the in-cloud peak β (corr ≈ −0.48); saturation would require a
  strong **positive** correlation (high signal → suppressed apparent S).
- **Aerosol contamination — ruled out (wrong sign).** Below-cloud aerosol adds to the integral and
  would push the cloud coefficient **down**, the opposite of the observed high bias; the
  90 %-in-cloud acceptance filter removes contaminated profiles. (corr(S, cloud-base) ≈ −0.48 is a
  mild geometric/aerosol signature, but in the wrong direction to explain a cloud-*high* offset.)

The molecular fit is correspondingly clean — it tracks the reference through an aerosol-free
1.7–5.6 km window (see `figs_paper_validation/rayleigh_diag/.../*_molecular_fit.png`).

**Verdict:** neither saturation nor aerosol explains the offset; it is a genuine inter-method
difference, consistent with §5.5. *Reproduce: `cl61_cloud_diagnostic.py`,
`run_rayleigh_diag_payerne.py`.*

> **Network-wide follow-up (2026-06-17):** this Payerne offset was confirmed across **all nine
> network CL61** and stress-tested (Kalman on/off, WV on/off, **L1 vs L2**). Cloud runs higher than
> Rayleigh at **every** site: network median **+26 %**, Payerne **+21 %** (consistent with the
> CHM-anchored +14–19 % here). Both calibrations are strongly water-vapor dependent (Rayleigh
> **+20 %**, cloud **−13 %**), and the Rayleigh is data-level-independent (L1 = L2 to ~1 %). A
> `parfor` bug had initially produced a spurious "0 % cloud WV" and an apparent ~17 %
> irreproducibility; run serially with WV correctly applied, the cloud calibration **reproduces**
> this study's value (C = 0.968). Full details: **`cl61_calibration_verification_report.md`**.

## 6. Standard atmosphere vs CAMS T/p (molecular profile) — now implemented

The Rayleigh molecular reference (β_mol ∝ P/T) needs a temperature/pressure profile. The
**original** MATLAB (`Auto_Calib_25\Rayleigh\auto_calib_v24.m`) built it from **ECMWF/MACC
reanalysis** at the site and date (`get_TPRH_MACC`, T/p interpolated to the range grid). The Python
port dropped that and defaulted to the **US Standard 1976** atmosphere (the back-port
`calibrateRayleighWithData.m` even states *"ECMWF loading not implemented. Using standard
atmosphere"*); CAMS was kept only for the WV correction. The Python calibration now has a
**`molecular_source` option** (`'standard'` | `'cams'`): `'cams'` takes the molecular T/p from the
actual CAMS profile at the site (the same monthly CAMS file as the WV correction, via the shared
ECMWF-L137 read, interpolated to the 30 m grid) — i.e. it **restores the original reanalysis-based
molecular profile**, CAMS being the modern successor to MACC.

**Validation against the original MATLAB.** `rayleigh_calibration_matlab` is a back-port *from* the
Python, so it is not an independent check. Feeding identical CAMS (and std) T/p to the Python
`calculate_molecular_properties` and to the **original** `get_rayleigh_v3.m` (`Auto_Calib_25`) gives
β_mol equal to **floating-point round-off** (max relative difference ≈ 1 × 10⁻¹⁵, original/Python
ratio 1.0000000000 over the 2–6 km fit window). The Bucholtz-1995 formula and constants
(T₀ = 288.15 K, N = 2.547 × 10²⁵ m⁻³, ρ = 0.0301, Edlén index, linear T/p interpolation) are
identical across the original, the Python and the back-port; the new CAMS path is
molecular-formula-exact.

**Measured impact.** Re-running the Payerne CL61 Rayleigh calibration over Feb–May 2026 **both
ways** (WV on, correct site coordinates, identical in every respect except the molecular source),
on the 18 clear nights:

| molecular source | median lidar constant C |
|---|---|
| standard (US 1976) | 0.6049 |
| CAMS T/p | 0.6030 |

The per-night paired ratio C(CAMS)/C(std) is **1.000 ± 0.007** (mean ± SD; median 1.002; ratio of
medians 0.997). The molecular-source choice moves the lidar constant by **< 0.3 % on average —
smaller than the ±0.7 % night-to-night calibration scatter**, i.e. effectively zero. It propagates
to a uniform ≈ 0.3 % rescale of β_att (≈ 0.3 percentage points on the validation bias), far below
every other sensitivity here. This confirms the earlier analytic estimate (< 1 %; the 3–6 km
US-1976-vs-CAMS air-density ratio was 0.9995).

→ **Negligible** (< 0.5 %): at a mid-latitude site the US Standard 1976 density is an excellent
proxy for the molecular reference, so the standard-vs-CAMS choice for the *molecular* profile does
not matter. Production therefore keeps `molecular_source = 'standard'`: it is always available
(CAMS covers only selected months), whereas `'cams'` would skip any night without a matching CAMS
month for a < 0.3 % gain. `'cams'` is offered for reanalysis-fidelity / sensitivity work — and is
what the original MATLAB effectively did. **For water vapour the choice is not optional**: a
standard atmosphere is dry (no humidity), so WV *requires* CAMS (or radiosonde) humidity — running
"with a standard atmosphere" for WV is the WV-off case (§2), i.e. a 16–25-point error.

## Summary — sensitivities ranked by impact

| Factor | Impact on CL61 bias | Note |
|---|---|---|
| Calibration **method** (cloud vs Rayleigh) | **≈ 14 pts** | genuine inter-method gap; **not** the MS factor (CL61 FOV ≈ CL51 → η appropriate; §5). Rayleigh (CHM-validated) is the anchor |
| **WV correction** in calibration | **≈ 18 pts** | essential; needs CAMS humidity |
| **Ångström** exponent | ≈ 14 pts per unit Å | state Å explicitly (Å=1 used) |
| **Cloud** lidar ratio S | ±5 pts (±0.8 sr); large over wider S | sets cloud absolute level |
| **Rayleigh** aerosol lidar ratio | ± 5 % | robust (clean reference) |
| Std-atm vs CAMS **molecular** profile | **< 0.3 %** (measured) | negligible; now selectable (`molecular_source`) |
| WV **laser-wavelength config** | **< 1 %** | manuf≈measured; CL31 broad band self-averages — see [wv_wavelength_sensitivity.md](wv_wavelength_sensitivity.md) |

**Bottom line.** The calibration *method* (≈14 pts) and the *WV correction* (≈18 pts) dominate, the
*Ångström* exponent and *cloud S* are secondary-but-significant scale factors that must be stated,
the *Rayleigh aerosol lidar ratio* is robust (±5 %), and the *standard-atmosphere vs CAMS molecular
profile* is negligible (<0.3 %, measured by re-running both ways) — while WV itself genuinely
requires CAMS humidity.

---
*Reproduce:* `sensitivity_payerne_4calib.m` (4 calibrations) and `sensitivity_payerne_wv.m`
(WV detail). WV-off calibrations are diagnostics in `…\fullcal_all_noWV` and `…\Cloud_noWV`; the
Ångström and cloud-S axes are exact analytic scalings of the calibration constant. The molecular
source is the new Python `molecular_source` option: `run_camsmol_payerne.py` re-runs the CL61
Rayleigh both ways (writing `…\fullcal_stdmol_check` and `…\fullcal_camsmol`); `validate_cams_molecular.py`
saves the β_mol + T/p, and `validate_cams_molecular_original.m` cross-checks it against the
**original** `Auto_Calib_25\Rayleigh\get_rayleigh_v3.m` (≈1e-15 match), `validate_cams_molecular_matlab.m`
against the `rayleigh_calibration_matlab` back-port. Both re-runs use the corrected 30 m L2 and
the real Payerne coordinates.

**CL61 lat/lon fix.** The raw CL61 files report `latitude = longitude = 0`, and the L2 files +
station manifest inherited it (same bug class as the `station_altitude = 0` already fixed in
RAW2L2). At 0,0 every water-vapour correction sampled CAMS in the Gulf of Guinea. Fixed:
`process_CL61_month.m` / `raw2L2_CL61.m` now fall back to the configured coordinates; the existing
L2 monthly files + manifest were patched in place (`patch_cl61_coords.py`); and both WV-on CL61
calibrations were recomputed with the correct coordinates (cloud: `recompute_cl61_cloud_payerne.m`
→ `Cloud_…_fixlatlon`; Rayleigh: the correct-coords `run_camsmol_payerne.py` → `fullcal_stdmol_check`;
the WV-off diagnostics use no CAMS, so they are lat/lon-independent). The §1–§4 numbers and the
Payerne validation figure were regenerated with these — this is what moved the Rayleigh + WV bias
from −8.5 % to +1.7 %. Production uses `apply_wv_correction = 1`, `molecular_source = 'standard'`,
Å = 1, S = 18.8 sr.
