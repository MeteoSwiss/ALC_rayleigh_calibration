# Water-vapor absorption correction for ~910 nm ceilometers — literature & code review

*Compiled 2026-06-16 for the E-PROFILE ALC β_att validation paper. Multi-source search with
adversarial verification (17 sources fetched, 79 claims extracted, 25 verified, 23 confirmed /
2 refuted).*

## Bottom line

The **Wiegner & Gasteiger (2015)** "WAPL" scheme is the **singular, anchor method** for
water-vapor (WV) absorption correction of 905–910 nm ceilometer attenuated backscatter. Every
later work either **validates it** (Wiegner et al. 2019), **substitutes a cheaper empirical
parameterisation** for operational use (Hopkin et al. 2019, the E-PROFILE operational
reference), or **acknowledges but neglects it** (Kotthaus et al. 2016). It is the standard, and
it is essentially unique — there is no competing rigorous formulation.

**No open, runnable implementation of a Wiegner-style *spectral* (HITRAN/MT-CKD line-by-line or
LUT) WV correction exists in any public repository.** CloudnetPy does **not** implement it
(confirmed by direct code inspection — see below); the E-PROFILE/EUMETNET operational lineage
uses the simpler empirical `Twv = 1 − 0.17·IWV^0.52`. **Our MATLAB code
(`compute_wv_transmission.m`, `wv_t2eff.m`, `apply_wv_correction_to_L2.m`) and its new Python
port (`rayleigh_calibration/water_vapor.py`) appear to be the only runnable open implementation
of the spectral WAPL-style correction** — a genuine novelty worth stating in the paper.

## Summary table

| Paper | Instrument / λ | Humidity source | Spectroscopy | WV correction method & magnitude | Code |
|---|---|---|---|---|---|
| **Wiegner & Gasteiger 2015**, AMT 8, 3971 ([doi](https://doi.org/10.5194/amt-8-3971-2015)) | Vaisala CT25k/CL31/CL51, 905–910 nm (not 1064 nm) | (method paper; any T,q profile) | **HITRAN 2005 + MT-CKD continuum**, ARTS line-by-line 895–930 nm @0.01 cm⁻¹, stored netCDF LUT (0.1 cm⁻¹, 10 m) | Spectral LUT; **averages over Gaussian laser spectrum** (λ₀+FWHM≈3.4 nm) → effective T²_w,eff. **Ignoring it biases β_p by ~20 % mid-lat, >50 % tropics (worst case; ~35 % typical tropical)** | WAPL netCDF archive = **private/on-request**, not open |
| **Wiegner et al. 2019** (CeiLinEx2015), AMT 12, 471 ([doi](https://doi.org/10.5194/amt-12-471-2019)) | Vaisala CL51 etc. vs RALPH ref. lidar | radiosonde / model | uses W&G 2015 LUT | **Validation** of W&G 2015 (multiplies signal by T_w,eff⁻²); near-range transmission agreement ~1–5 % | none |
| **Hopkin et al. 2019**, AMT 12, 4131 ([doi](https://doi.org/10.5194/amt-12-4131-2019)) — *E-PROFILE operational ref.* | CL31/CL51/CT25k/CS135 @910 nm + CHM15k @1064 nm | **NWP** (Met Office UKV; ECMWF profiles *provided by M. Hervo, MeteoSwiss*) | **none** (empirical) | O'Connor liquid-cloud calib (S=18.8±0.8 sr) **+ empirical** `Twv = 1 − 0.17·IWV^0.52` (Markowicz 2008). Agrees with WAPL **to within 2 %**; ~12 % annual cycle if ignored | multiple-scattering code (Hogan 2006) public; **no WV repo** |
| **Kotthaus et al. 2016**, AMT 9, 3769 ([doi](https://doi.org/10.5194/amt-9-3769-2016)) | Vaisala CL31, 905±10 nm, FWHM~4 nm, 0.3 nm/K | — | — | **Acknowledges** WV sensitivity (cites W&G 2015, Markowicz 2008) but **neglects it** (c_absolute=1) | none |
| Markowicz et al. 2008 | CT25k, ~905 nm | — | — | Case-specific predecessor; basis for Hopkin's empirical fit | none |
| **CloudnetPy** (ACTRIS-Cloudnet) | CL31/CL51, CHM15k | — | **none** | **No WV correction.** Scalar site `calibration_factor` (Vaisala 1.0, Lufft 3e-12, "probably incorrect") + O'Connor/Hogan liquid-cloud classification | [github.com/actris-cloudnet/cloudnetpy](https://github.com/actris-cloudnet/cloudnetpy) — open, but **no WAPL** |

## Key points for the paper

1. **W&G 2015 is the method to cite** as the rigorous benchmark. Its three defining
   ingredients — HITRAN+MT-CKD cross-sections, line-by-line LUT, and **averaging over the laser
   emission spectrum (λ₀ + FWHM)** — are exactly what our implementation reproduces.
2. **Magnitude**: ~20 % on retrieved backscatter at mid-latitudes (>50 % tropics, worst case).
   Quote the >50 % as an *upper bound* (W&G stress there is "no generally applicable value" —
   it depends on height, aerosol load, algorithm). The 2025 CL51 study (Jin et al., MDPI
   rs17122013) independently states ">20 % if water vapor correction is ignored." *Caveat:* that
   paper reports the magnitude but does **not** adopt the full WAPL pipeline — do not cite it as
   a WAPL application.
3. **Operational vs rigorous split**: Hopkin et al. (2019) deliberately chose the cheap
   empirical `Twv = 1 − 0.17·IWV^0.52` over WAPL "because it requires a radiative-transfer model
   or access to their WAPL database." Notably their WAPL-vs-empirical comparison (their Fig. 5)
   **used ECMWF water-vapor profiles you (M. Hervo) provided** — a direct link between this
   paper's author and the operational reference. The two agree to 2 %.
4. **Our novelty**: a full **spectral WAPL-style correction, now in both MATLAB and open Python**
   (validated against MATLAB to ≤0.4 % and against ACTRIS-Cloudnet `atmoslib` — see
   `tests/WATER_VAPOR_AUDIT.md`), driven by CAMS humidity + the ECMWF L137 model levels. This is
   the **only open runnable spectral implementation** found.
5. **CloudnetPy confirmation** (your prior finding is correct): direct inspection of CloudnetPy
   (HEAD 0db960b8) found **zero** references to HITRAN / MT-CKD / WV transmission / spectral
   averaging in the ceilometer code; it uses a scalar calibration factor + O'Connor liquid-cloud
   classification. So CloudnetPy is *not* a WV-correction reference — `atmoslib` (its
   thermodynamics library) is the right external anchor for the humidity primitives only.

## Open follow-ups (not resolved by this search)
- The **Vande Hey (2015) thesis** and **Madonna et al. (2018)** were named but did not surface
  as verified claims — worth a manual check for any independent spectral WV treatment.
- Whether the original **W&G WAPL netCDF archive** was ever published openly (Zenodo / LMU
  institutional) so it could be cited/reused rather than regenerated.
- State explicitly in the paper which **HITRAN edition + MT-CKD version** our
  `abs_cross_647_full_levels_1000.nc` LUT uses, vs W&G's HITRAN-2004/MT-CKD basis.

## Verified sources
- Wiegner & Gasteiger 2015 — https://amt.copernicus.org/articles/8/3971/2015/
- Wiegner et al. 2019 (CeiLinEx2015) — https://amt.copernicus.org/articles/12/471/2019/
- Hopkin et al. 2019 — https://amt.copernicus.org/articles/12/4131/2019/
- Hopkin PhD thesis — https://centaur.reading.ac.uk/85509/1/15013194_Hopkin_thesis.pdf
- Kotthaus et al. 2016 — https://amt.copernicus.org/articles/9/3769/2016/
- O'Connor et al. 2004 — https://journals.ametsoc.org/view/journals/atot/21/5/1520-0426_2004_021_0777_atfaoc_2_0_co_2.xml
- CloudnetPy — https://github.com/actris-cloudnet/cloudnetpy
- HITRAN / MT-CKD — https://hitran.org/mtckd/
- Jin et al. 2025 (CL51, magnitude only) — https://www.mdpi.com/2072-4292/17/12/2013
