# Validation of E-PROFILE attenuated backscatter — technical documentation

*Maxime Hervo (MeteoSwiss), June 2026. Working notes for the E-PROFILE ALC paper.*

> **2026-06-21 update (Python calibration campaign).** Coefficients now follow the **Wiegner lidar
> constant `C_L = RCS/β_att`** throughout ([calibration_coefficient_convention.md](calibration_coefficient_convention.md));
> the cloud O'Connor coefficient `C` (`β_true = C·β_L2`) maps to it as `C_L = calibration_constant_0 / C`.
> Two findings from the same campaign are folded in below: **(i)** the L1↔L2 Rayleigh difference is a
> **method/grid interaction, not a data difference** — the gated methods over-reject the fine *native*
> grid; fixed by binning native L1/RAW to the L2 grid (30 m × 300 s, `l1_bin_to_l2_grid`; directly
> relevant to §7.6) — see [network_v2_vs_v11_report.md](network_v2_vs_v11_report.md). **(ii)** the cloud
> calibration prefers the **native cadence** (the temporal-consistency gate) — see
> [cloud_optimization_report.md](cloud_optimization_report.md). The Python package was renamed
> `rayleigh_calibration` → `calibration`; the molecular `calipso` method has been retired.

## 1. Objective and scope

We validate the **calibrated attenuated backscatter coefficient** β_att produced by the
E-PROFILE automatic lidars and ceilometers (ALC). The calibration under test is the
**corrected Python Rayleigh calibration** (a re-implementation that fixed a sign error in the
earlier code) combined with a **Kalman filter** that turns per-night lidar constants into a
smooth daily time series. Validation uses six independent strategies:

1. **EARLINET** — ceilometer (CHM15k) β_att against the colocated EARLINET research-lidar
   reference, at four sites.
2. **Payerne** — three different instruments at one site (CHM15k, CL31, CL61).
3. **Amsterdam/Schiphol** — four nominally identical CHM15k (instrument-to-instrument spread).
4. **A CL51 + CL61 site (Uccle)** — CL51 against CL61, with the CL61 shown under **two
   independent calibrations** (liquid-cloud and Rayleigh), i.e. treated as two instruments.
5. **SIRTA/Palaiseau** — three colocated E-PROFILE instruments at **three wavelengths**
   (CHM15k 1064 nm, CL31 910 nm, Mini-MPL 532 nm). The 532 nm Mini-MPL is an independent
   research-grade reference at a third wavelength, making SIRTA the strongest cross-wavelength
   check in the set (and SIRTA also hosts the EARLINET `sir` system of strategy 1).
6. **Lindenberg** — the E-PROFILE CHM15k (1064 nm) against an independent **ACTRIS-Cloudnet
   CL61** (910 nm) obtained from the Cloudnet archive, i.e. a cross-network / cross-data-source
   check of the same calibration applied to two different processing chains.

All comparisons are of the *calibrated* product, so they jointly test (i) the absolute
calibration, (ii) the water-vapor and wavelength harmonisation, and (iii) instrument
hardware differences.

## 2. Data

### 2.1 Ceilometer / lidar L2
- **E-PROFILE L2** attenuated backscatter, read as monthly-concatenated files
  `A:\E-PROFILE_L2_monthly\<WIGOS>\<yyyy>\L2_<WIGOS>_<id><yyyymm>.nc`
  (`read_L2_monthly.m`). Recent daily files (`D:\E-PROFILE_L2_2026`) are concatenated to the
  monthly archive with `Eprofile_alc_monthly_concatenation_L2.m`.
- **EARLINET** Level-2 1064 nm particle backscatter, converted to *attenuated* backscatter
  in `read_earlinet_att_backscatter.m` (assumed lidar ratio 50 sr, molecular two-way
  transmission from the US-1976 standard atmosphere, overlap hold below a site-specific
  minimum range). Stations: `sir` (Palaiseau), `lei`+`ari` (Leipzig), `cbw` (Cabauw),
  `ino` (Magurele).

### 2.2 Atmospheric ancillary
- **CAMS** monthly T/RH profiles (`D:\CAMS\CAMS_Beta_YYYYMM.nc`) for the water-vapor
  correction. Coverage currently ends 2026-02; later months fall back to the nearest
  available month.
- **HITRAN/MT-CKD** water-vapor absorption cross-section LUT
  `...\MDA\monitoring_alc_monthly\abs_cross_647_full_levels_1000.nc`.

## 3. Calibration

### 3.1 Rayleigh (Python, sign-corrected) + Kalman
Per-night lidar constants C_L are produced by the Python package in
`C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration` and written as
`D:\E-PROFILE_calibration_rayleigh\fullcal_all\<WIGOS>_<id>\<WIGOS>_<id>_cl.csv`
(flag ∈ {1, 0.5} = clear / partially-clear night). `load_rayleigh_python_kalman.m`:

1. keeps successful nights, **rejects outliers** with |log C_L − median(log C_L)| > 4·MAD
   (the raw series occasionally contains physically impossible C_L up to ~10¹⁷ that would
   wreck the filter),
2. **normalises** to a dimensionless coefficient C_L / median(C_L) because the Kalman bridge
   (`run_kalman_from_matlab.py`) is tuned for O(1) values,
3. runs the operational **Kalman filter** (daily grid, predict/update, process-noise floor),
4. rescales back to C_L units and returns the legacy `rayleigh_<WIGOS>_<id>.mat` schema
   (`daily_dates, daily_C_kalman, daily_C_std_kalman`). A short/gappy record that the Kalman
   cannot handle falls back to the constant median C_L.

The correction applied to the L2 product is
**β_true = β_L2 · (C_op / C_Rayleigh)**, where C_op is the operational
`calibration_constant_0`. Outside the calibrated period the daily series is held (clamped),
so recent profiles use the latest stable calibration.

### 3.2 Liquid-cloud calibration
For Vaisala instruments (CL31/CL51/CL61) the alternative is the O'Connor/Hopkin liquid-cloud
calibration (`A:\E-PROFILE_L2_Calibration\Cloud_Trans-cor_WV-cor\calibration_<WIGOS>_<id>.mat`,
Kalman already applied). Convention: **β_true = β_L2 · C** (dimensionless multiplier).

The calibration time series for every instrument (raw daily values and the Kalman estimate)
is shown below; it reveals seasonal calibration cycles tracked by the filter, a step change
in the Payerne CL31 around 2024, and the widening Kalman uncertainty across data gaps
(e.g. Magurele 2022–2024).

![Calibration coefficient time series for all instruments: raw daily values (× crosses) and the Kalman estimate (red line, ±1σ band). Rayleigh channels show the Wiegner lidar constant C_L = RCS/β_att (a.u.); cloud channels the O'Connor multiplier C (β_true = C·β_L2), which maps to the Wiegner constant as C_L = calibration_constant_0 / C (see calibration_coefficient_convention.md).](figs_paper_validation/calibration_timeseries.png)

## 4. Harmonisation corrections

### 4.1 Water vapor (910 nm only)
CL31/CL51/CL61 emit in the 905–911 nm water-vapor absorption band, so their β_L2 is
attenuated by H₂O that does **not** affect the 1064 nm CHM15k. `apply_wv_correction_to_L2.m`
(+ `compute_wv_transmission.m`) computes the spectrally-averaged two-way transmission
T²_wv(z) from CAMS humidity and the HITRAN LUT, weighted by each instrument's measured laser
spectrum, and divides β by T²_wv. Per-instrument laser parameters (2026-06-02 Qmini campaign):

| Instrument | λ₀ [nm] | FWHM [nm] | WV correction |
|---|---|---|---|
| CL31 | 909.7 | 6.0 | yes |
| CL51 | 910.0 | 3.4 | yes |
| CL61 | 910.74 | 1.0 | yes |
| CHM15k | 1064.47 | 0.5 | no (outside band) |

Typical effect: median T²_wv ≈ 0.80–0.87 (i.e. β boosted ~15–25 %), increasing with the
water-vapor column above the gate.

### 4.2 Wavelength
β is normalised to a common wavelength with an Ångström scaling
`(λ/λ_target)^(−Å)`, Å = 1 (Haarig 2025, 532/1064 nm at high RH). λ_target = 1064 nm for the
CHM/EARLINET comparisons; 910 nm for the all-Vaisala Uccle comparison (then a no-op).

## 5. Matching, gridding, statistics

- **Colocated continuous instruments** (`paper_val_process.m`): each channel is calibrated,
  WV- and wavelength-corrected, quality-/cloud-/rain-screened, retimed to 60-min medians on a
  common regular grid (union of channels), and averaged onto a common altitude grid. A
  *display* stream (quality-only) feeds the time-height pcolors; a *screened* stream
  (cloud/rain/fog excluded ±15 min) feeds the profile, scatter and statistics. One physical
  instrument may appear as several channels (e.g. CL61 cloud + Rayleigh).
- **EARLINET** (`paper_val_earlinet.m`): EARLINET provides sparse nighttime profiles and the
  CHM record is gappy, so the CHM is read only for months containing EARLINET profiles and
  matched profile-by-profile (CHM profiles within ±30 min of each EARLINET time, median),
  interpolated to the EARLINET altitude grid.
- **Statistics** over 500–3000 m AGL (EARLINET 500–5000 m): N, mean/median bias, RMSE,
  relative bias and Pearson r versus the reference channel.

## 6. Figures
- Colocated sections (`paper_val_figure.m`): one example profile (all channels overlaid),
  one time-height pcolor per channel, and a scatter of every channel vs the reference.
- EARLINET (`paper_val_earlinet_figure.m`): median matched profile ± IQR, scatter
  (CHM vs EARLINET), and two matched-profile "curtains" (EARLINET and CHM by date).
- Output: `figs_paper_validation/validation_*.png` (300 dpi).

## 7. Results

Comparison band 500–3000 m AGL (EARLINET 500–5000 m). "rel. bias" is the mean of
(comparison − reference) over reference; r is the Pearson correlation of the paired
hourly (colocated) or matched (EARLINET) attenuated-backscatter values.

> **Status (2026-06-16, updated).** All four sections are now final. CAMS humidity is
> available through 2026-05, so the **910 nm** sections (Payerne §7.1, Uccle §7.3) use a real
> per-month water-vapor correction (no fallback); months without CAMS (e.g. June 2026) are
> excluded, not fallback-corrected. The Python Rayleigh calibration applies the WV correction
> itself (calibration/water_vapor_correction/water_vapor.py, validated against the MATLAB reference to 0.4 %
> and against ACTRIS-Cloudnet atmoslib — see tests/WATER_VAPOR_AUDIT.md). The CL61 Rayleigh
> constants were regenerated WV-corrected: e.g. CL_L2/CL_Rayleigh = 0.502/0.638 ≈ 0.79 ≈
> T²_wv(ref) at Payerne, confirming the operational L2 constant was the WV-*uncorrected* one and
> the Rayleigh now recovers the true constant. **This materially changed the 910 nm numbers (see
> §7.1, §7.3).** A subsequent fix corrected the Payerne CL61 **station coordinates** (the raw CL61
> files report lat/lon = 0, which had placed the WV CAMS read at 0° N 0° E — the Gulf of Guinea):
> with the real Payerne coordinates the CL61 Rayleigh bias moves from −8.5 % to **+1.7 %**, the
> near-zero agreement expected between two molecular-calibrated instruments. Both CL61 calibrations
> and the Payerne validation were regenerated with the correct coordinates.

### 6.1 Dataset homogeneity (calibration applied per instrument)

The dataset is homogeneous: every channel is Kalman-filtered, every 910 nm channel is
water-vapor-corrected (910 nm sits in the H₂O absorption band; 1064 nm does not), and all
channels are compared on a common 30 m grid.

| Instrument | λ (nm) | Calibration | Kalman | WV correction | Grid |
|---|---|---|---|---|---|
| CHM15k (CHM15k ref, Amsterdam, EARLINET) | 1064 | corrected Python Rayleigh | ✓ | — (outside band) | 30 m |
| Mini-MPL (SIRTA) | 532 | corrected Python Rayleigh | ✓ | — (outside band) | 30 m |
| CL31 | 909.7 | liquid-cloud (O'Connor, WV-corrected) | ✓ | ✓ | 30 m |
| CL51 | 910.0 | liquid-cloud (O'Connor, WV-corrected) | ✓ | ✓ | 30 m |
| CL61 — instrument #1 | 910.74 | liquid-cloud (O'Connor, WV-corrected) | ✓ | ✓ | 30 m |
| CL61 — instrument #2 | 910.74 | corrected Python Rayleigh | ✓ | ✓ | 30 m |

For the 910 nm channels the WV correction is applied **once** physically: the calibration step
removes the WV bias from the constant (cloud C from WV-corrected integrated backscatter; Rayleigh
constant from WV-divided RCS), and the comparison step removes the range-dependent two-way WV
transmission T²_wv(z) from the profile so it is comparable to the WV-free 1064 nm reference.
Kalman: `load_rayleigh_python_kalman` for the Rayleigh channels, `daily_C_kalman` for the cloud
channels. (Confirmed in the run log: CHM channels carry no `[WV]` tag; every CL31/CL51/CL61
channel carries `[WV T²≈0.78–0.82]`; every channel loads/applies a Kalman series.)

### 7.1 Payerne — 4 channels, CL61 as two instruments (Mar–May 2026, ref = CHM15k Rayleigh)
| Comparison | calib | rel. bias | RMSE [Mm⁻¹sr⁻¹] | r | N |
|---|---|---|---|---|---|
| CL31 vs CHM15k | cloud + WV | +31.8 % | 0.46 | 0.48 | 79 348 |
| CL61 vs CHM15k | **cloud + WV** | **+16.0 %** | 0.07 | **0.987** | 77 522 |
| CL61 vs CHM15k | **Rayleigh + WV** | **+1.7 %** | 0.06 | **0.987** | 77 522 |

**The CL61 is shown as two instruments (cloud + Rayleigh), both WV-corrected, both correlating
with the CHM15k at r ≈ 0.987 on the matched 30 m grid.** Against the 1064 nm Rayleigh reference the
molecular-Rayleigh CL61 is essentially unbiased (**+1.7 %**) while the liquid-cloud CL61 reads
**+16.0 %** — a **≈14-point spread between the liquid-cloud (O'Connor S = 18.8 sr) and
molecular-Rayleigh calibrations** of the same instrument. The CHM15k(Rayleigh)–CL61(Rayleigh) pair
is the cleanest absolute check (both molecular, reference needs no WV): the two agree to **+1.7 %**,
the near-zero result expected for two molecular-calibrated instruments. **The +14 % cloud–Rayleigh
offset is a pure, range-independent calibration-constant difference; it is *not* explained by
multiple scattering (the CL61 receiver FOV ≈ the CL51's, so its η is appropriate — confirmed by the
CL61 in-cloud depolarisation) but is a genuine cloud-vs-Rayleigh inter-method discrepancy; see
[payerne_cl61_calibration_sensitivity.md](payerne_cl61_calibration_sensitivity.md) §5.**

Four CL61 fixes were needed to reach this:
1. **`station_altitude` was 0** in raw2L2 (raw CL61 `elevation`=0) → `altitude` was AGL not ASL →
   ~491 m vertical misalignment vs CHM. `process_CL61_month.m` now falls back to 491 m.
2. **Native 4.8 m grid** broke the O'Connor cloud calibration (the aerosol-ratio filter
   over-rejects on fine gates → 0 valid profiles). raw2L2 now **reproduces the E-PROFILE L2 grid
   (30 m vertical averaging, 5-min temporal)**, on which the cloud calibration succeeds (671
   profiles, C ≈ 0.99). Verified unbiased (30 m/native β ratio = 1.000) and consistent with the
   official E-PROFILE 30 m files (cloud C 1.05 vs 1.08 on the same days).
3. **`build_common_grid` used an LCM rule** that exploded to a 420 m comparison grid when one
   channel was slightly off 30 m (lcm(28,30)); replaced by the coarsest-channel resolution → a
   stable 30 m grid (N rose from ~6 k to ~83 k; r is now the honest 30 m value, not the inflated
   0.99 of the coarse grid).
4. **`station_latitude/longitude` were 0** in the raw CL61 (inherited by the L2 + manifest), so
   every water-vapor correction — the cloud constant, the Rayleigh constant, and the display-side
   T²_wv — sampled CAMS at 0° N 0° E (Gulf of Guinea) instead of Payerne. `process_CL61_month.m`/
   `raw2L2_CL61.m` now fall back to the real coordinates; the existing L2 files + manifest were
   patched (`patch_cl61_coords.py`) and both CL61 calibrations recomputed. **This moved the CL61
   Rayleigh bias from strongly negative (≈ −8.5 %, the wrong WV location) to near-zero (+1.7 %)** and
   raised the cloud accordingly — the earlier negative bias was largely an artefact of the wrong WV
   location.

**CL31 is degraded over the spring (Mar–May)** (+31.8 %, r 0.48): a known CL31 behaviour
(window/optics in the April–May high-load period); a clean late-winter window gives a much smaller
CL31 bias, so the spring figure here is the conservative case.

![Payerne four-channel validation (Mar–May 2026) in a 3×3 layout: profile median ± IQR (left column), scatter vs CHM15k (b), histogram of channel−CHM15k differences in Mm⁻¹ sr⁻¹ (c), and time-height attenuated backscatter for CHM15k (Rayleigh), CL31 (cloud), CL61 (cloud) and CL61 (Rayleigh) (d–g) — all WV-corrected, CL61 on the E-PROFILE-style 30 m grid. Black dots = cloud-base-height detections.](figs_paper_validation/validation_payerne.png)

### 7.2 Amsterdam/Schiphol — 4× CHM15k (Mar–May 2026, ref = CHM15k A)
| Comparison | rel. bias | RMSE | r | N |
|---|---|---|---|---|
| CHM15k B vs A | +24.0 % | 0.16 | 0.96 | 28 884 |
| CHM15k C vs A | −14.9 % | 0.13 | 0.93 | 32 868 |
| CHM15k D vs A | −15.6 % | 0.12 | 0.94 | 34 196 |

Four nominally identical CHM15k are highly correlated (r 0.92–0.96) but spread ≈±20 %
in absolute calibration — a direct measure of the instrument-to-instrument calibration
uncertainty after Rayleigh calibration.

![Four colocated CHM15k at Amsterdam/Schiphol (Mar–May 2026): profile median ± IQR, scatter vs instrument A, and time-height attenuated backscatter for A/B/C/D.](figs_paper_validation/validation_amsterdam.png)

### 7.3 Uccle — CL51 vs CL61, CL61 as two instruments (Mar–May 2026, ref = CL51 cloud)
| Comparison | calib | rel. bias | RMSE | r | N |
|---|---|---|---|---|---|
| CL61 vs CL51 | cloud + WV | −2.9 % | 0.29 | 0.94 | 54 365 |
| CL61 vs CL51 | **Rayleigh + WV** | **−34.7 %** | 0.31 | 0.94 | 54 365 |

*(Mar–May 2026, the common period with §7.1–7.2; reference is the CL51 liquid-cloud calibration.)*

The **cloud and Rayleigh calibrations of the same CL61 disagree** (−2.9 % vs −34.7 % relative to the
CL51). Crucially the reference here is a *cloud-calibrated* CL51, which carries the **same**
liquid-cloud over-correction (the CL51 multiple-scattering factor; §7.1 and
[payerne_cl61_calibration_sensitivity.md](payerne_cl61_calibration_sensitivity.md) §5), so this
comparison is **confounded** — both the CL51-cloud reference and the CL61-cloud channel are inflated
together, which is why the CL61-Rayleigh appears far below them. The **clean** CL61 cloud-vs-Rayleigh
test is therefore Payerne against the molecular 1064 nm CHM15k (§7.1: cloud +16.0 % vs Rayleigh
+1.7 %, a 14-point gap); §5 shows this is a genuine inter-method discrepancy, **not** a
multiple-scattering/FOV error (the CL61 receiver FOV equals the CL51's, and its in-cloud
depolarisation confirms CL51-class multiple scattering).
This warrants follow-up (cloud-calibration S-value / multiple-scattering for CL61) and is a
result in its own right.

![Uccle CL51 vs CL61, with the CL61 shown under both cloud and Rayleigh calibration (Mar–May 2026): profile median ± IQR, scatter vs CL51, and time-height attenuated backscatter. Black dots = cloud-base-height detections.](figs_paper_validation/validation_cl51_06447.png)

### 7.4 EARLINET — CHM15k vs research-lidar reference (vs CHM15k Rayleigh)
| Site | rel. bias | r | matched N |
|---|---|---|---|
| Palaiseau (sir) | −4.5 % | 0.37 | 362 |
| Leipzig (lei) | +1.7 % | 0.21 | 5 096 |
| Magurele (ino) | −7.3 % | 0.25 | 1 049 |
| Cabauw (cbw) | +23.2 % | 0.30 | 29 |
| Leipzig (ari) | +36.7 % | 0.16 | 1 667 |

The well-sampled sites (Palaiseau, Leipzig-lei, Magurele) agree with EARLINET to within
≈±7 % in the median; Cabauw (29 profiles) and the ari system are less reliable. Profile-level
r is low by construction (different technique, ±30 min sampling, assumed EARLINET lidar ratio
of 50 sr); the median-profile agreement and bias are the meaningful metrics.

![EARLINET validation at Palaiseau: median matched profile (±IQR), scatter, and matched EARLINET/CHM15k curtains.](figs_paper_validation/validation_earlinet_sir.png)

![EARLINET validation at Leipzig (lei).](figs_paper_validation/validation_earlinet_lei.png)

![EARLINET validation at Magurele.](figs_paper_validation/validation_earlinet_ino.png)

### 7.5 SIRTA/Palaiseau — three instruments at three wavelengths (Mar 2025–Feb 2026, ref = CHM15k Rayleigh)

SIRTA is the strongest **cross-wavelength** check in the set: three colocated E-PROFILE
instruments spanning **532, 910 and 1064 nm**. The CHM15k (1064 nm, Python Rayleigh) is the
reference; the CL31 (910 nm) is liquid-cloud + WV-corrected; the **Mini-MPL (532 nm)** is an
independent, research-grade micropulse lidar, Rayleigh-calibrated by the same Python pipeline.
A full **12-month, all-season** window is used (the Mini-MPL L2 ends 2026-02; CAMS covers it,
so the CL31 WV correction is real throughout).

| Comparison | calib | λ (nm) | wavelength model | rel. bias | RMSE [Mm⁻¹sr⁻¹] | r | N |
|---|---|---|---|---|---|---|---|
| CL31 vs CHM15k | cloud + WV | 910 | Ångström (Å=1) | **−16.4 %** | 0.22 | **0.81** | 191 232 |
| Mini-MPL vs CHM15k | Rayleigh | 532 | molecular + aerosol | −45.5 % | 0.19 | **0.79** | 88 063 |

- **CL31 (910 nm) reads −16.4 %** below the 1064 nm CHM15k over the full year (r 0.81) — within
  the network CL31 spread (Payerne CL31 ran +6 % in a clean window and +32 % in the high-load
  spring; SIRTA's −16 % is the all-season figure for a different unit), and consistent with the
  CL31 being the noisiest 910 nm class.
- **The Mini-MPL (532 nm) is an honest cross-wavelength tracker but its *absolute* bias is
  intrinsically uncertain.** Comparing 532 nm β_att to 1064 nm is hard because molecular
  backscatter scales as ≈λ⁻⁴, so β_mol(532) ≈ 17 × β_mol(1064) and the molecular term dominates
  the 532 nm signal. A single *aerosol* Ångström exponent (Å=1) under-normalises the molecular
  part and gives **+206 %**; the physically-correct two-component conversion (molecular β_att from
  the standard atmosphere via `get_rayleigh_v3`, aerosol remainder scaled by Å, recombined at
  1064 nm — `paper_val_process.m/convert_wavelength_molaer`) gives **−45.5 %**, because the
  aerosol remainder β(532)−β_mol(532) is a small difference of large numbers and the conversion
  neglects the (larger) aerosol two-way transmission at 532 nm. The true value is bracketed
  between these; the robust metric is therefore the **temporal correlation (r ≈ 0.79 over a full
  year)** and the time-height consistency, not the absolute bias. This is a genuine, citable
  limitation of absolute cross-wavelength ALC intercomparison, and is why the 910/1064 pairs
  (close wavelengths, molecular term small) give clean absolute biases while 532/1064 does not.
- The figure shows all three instruments resolving the **same** time-height structure over the
  year, with the profile overlay and the scatter/difference panels against the CHM15k.

![SIRTA/Palaiseau three-wavelength validation (Mar 2025–Feb 2026): profile median ± IQR for CHM15k 1064 nm (Rayleigh), CL31 910 nm (cloud) and Mini-MPL 532 nm (Rayleigh) (a); scatter vs CHM15k (b); histogram of channel−CHM15k differences (c); and time-height attenuated backscatter for CHM15k (d), CL31 (e) and the Mini-MPL (f) — the Mini-MPL converted to 1064 nm with the two-component molecular+aerosol model. Black dots = cloud-base-height detections.](figs_paper_validation/validation_sirta.png)

### 7.6 Lindenberg — E-PROFILE CHM15k vs ACTRIS-Cloudnet CL61 (cross-source)

Lindenberg carries an E-PROFILE **CHM15k** (`0-20000-0-10393`, 1064 nm, Python Rayleigh;
2025-01→2026-03) and **no** colocated E-PROFILE CL31/CL61. Its CL61 (910 nm) is instead
available from the **ACTRIS-Cloudnet** archive, so Lindenberg is a **cross-network /
cross-data-source** check: the *same* Python calibration applied to a different processing
chain. The Cloudnet CL61 raw was pulled through the open Cloudnet API
(`download_cloudnet_longrun.py`, PID `21.12132/3.695573e5981845d9`, 200+ day-folders into the
RAW per-day layout) and calibrated with the RAW reader (WV-corrected) and the seven
molecular-window methods (`calibrate_cloudnet_cl61.py`):

| method | nights calibrated / 49 | median C_L | robust CV | median rel. err |
|---|---|---|---|---|
| eprof_v1.2 (improved) | 18 | 2.77 | 46 % | 16 % |
| eprof_v1.1 (main) | 18 | 2.11 | 45 % | 10 % |
| eprof_v0.25 (matlab) | 16 | 1.86 | 49 % | 12 % |
| earlinet | 2 | 2.21 | 30 % | 8 % |
| eprof_v2 / bellini | 0 | — | — | — |

*(Original native-grid run, relabelled to the current method names; the `calipso` row is removed —
that method is retired. These are **native RAW 4.8 m × 60 s** results.)*

The Cloudnet CL61 **calibrates cleanly from raw via the same pipeline** (the high-yield methods
on 16–18 of 49 sampled nights, median relative error 9–16 %), demonstrating that the calibration
is portable across data sources. The strict gated methods (**`eprof_v2`, `bellini`**) calibrate **0**
nights here — but this is now understood to be the **native-grid handicap**, not a property of the
CL61 data: the gated v2 over-rejects the fine native grid (its temporal-variability gate on the noisy
per-profile stack), exactly as on native L1 CHM15k. The fix (`l1_bin_to_l2_grid`: bin native L1/RAW to
the L2 30 m × 300 s grid before the fit, commit `a0873f8`) recovers v2 — see
[network_v2_vs_v11_report.md](network_v2_vs_v11_report.md). **Re-running `calibrate_cloudnet_cl61.py`
with the grid fix (RAW → L2 grid) is the immediate next step and is expected to lift the v2 yield to
the v1.1/v1.2 level.** The **profile-level β_att intercomparison** against the
colocated E-PROFILE CHM15k (overlap 2025-01→2025-04) is the cross-source validation enabled by
this calibration; the calibration constants and their night-to-night stability are reported here,
and the gridded β_att comparison is the immediate next step.

## 8. Data limitations / caveats
- **Payerne CL61**: E-PROFILE L2 regenerated from raw (E:\CL61_PAY via RAW2L2, Feb–Jun 2026,
  ~4 months). Rayleigh now succeeds (23 clear nights) and is WV-corrected, but only February
  has CAMS, so the WV-valid Rayleigh is presently limited to Feb (the rest is held for CAMS).
  The liquid-cloud calibration detects clouds but the O'Connor opacity/aerosol-ratio filters
  reject every Payerne CL61 profile (CL61 multiple scattering + 5-min averaging), so no cloud C
  is produced — the robust CL61 cloud-vs-Rayleigh comparison is therefore at Uccle.
- **EARLINET vs ceilometer** profile-level correlation is intrinsically low (different
  technique, ±30 min sampling, assumed EARLINET lidar ratio); the median-profile agreement
  and bias are the meaningful metrics.
- **Water-vapor rule**: 910 nm instruments are only calibrated/compared when the matching
  month's CAMS is available — no nearest-month fallback (apply_wv_correction_to_L2 excludes
  months without CAMS; the Python Rayleigh skips such nights). CAMS currently ends 2026-02, so
  spring 910 nm work is held until Mar–Jun CAMS is downloaded.
- **calibration_constant_0** conventions and WIGOS folder/filename can be inconsistent across
  the network (e.g. Palaiseau files are `0-250-1001-07151` inside a `0-20000-0-07151` folder).
- **Site selection for the added comparisons.** A survey of the E-PROFILE L2 archive for further
  multi-instrument sites retained **SIRTA** (CHM15k + CL31 + Mini-MPL, three wavelengths) and
  **Lindenberg** (E-PROFILE CHM15k + ACTRIS-Cloudnet CL61). **Ljubljana** (`0-20000-0-14015`) was
  examined but carries a single ALC (one CL51) in the archive — no CL61 and no second colocated
  instrument — so it does not meet the "CL61 + a second instrument over several months" criterion
  and is not used. At SIRTA the **Mini-MPL L2 ends 2026-02**, so its tri-wavelength window is
  Mar 2025–Feb 2026 (12 months; CAMS covers it, so the CL31 WV correction is real, no fallback).
  An all-season window is also *required* for the Mini-MPL: in clean winter air the 532→1064
  molecular+aerosol conversion is ill-posed (the aerosol residual is a tiny difference of large
  numbers), so a window including the aerosol-rich summer is needed (§7.5).
- **Absolute 532 nm ↔ 1064 nm intercomparison (Mini-MPL) is intrinsically limited** — molecular
  backscatter scales as ≈λ⁻⁴, so the 532 nm signal is molecular-dominated and the absolute bias is
  bracketed by the wavelength model (+206 % aerosol-Ångström vs −45 % two-component, §7.5). The
  Mini-MPL is used as an independent 532 nm *tracker* (temporal correlation, time-height
  consistency); its absolute calibration is not a clean validation reference the way the close
  910/1064 pairs are.

## 9. Reproducibility
Run `make_validation_figures.m` (`sectionsToRun = {payerne, amsterdam, cl51, sirta, earlinet}`).
SIRTA alone: `run_sirta_validation.m` (avoids re-running the four heavy sections). Key functions:
`paper_val_process.m` (now supports per-channel `WIGOS_ID` for a colocated instrument on its own
WIGOS id, and a per-channel `wavelengthModel = 'molaer'` two-component conversion for far-apart
wavelengths such as the 532 nm Mini-MPL), `paper_val_figure.m`, `paper_val_figure_payerne.m`,
`paper_val_earlinet.m`, `paper_val_earlinet_figure.m`, `load_rayleigh_python_kalman.m`,
`compute_wv_transmission.m`, `apply_wv_correction_to_L2.m`, `interp_calib.m`, `get_rayleigh_v3.m`,
`get_std_atm.m`, `read_earlinet_att_backscatter.m`, `read_L2_monthly.m`.

The Lindenberg cross-source CL61 (§7.6) is in `C:\Users\hervo\OneDrive\Documents\ALC_rayleigh_calibration`:
`download_cloudnet_longrun.py` (Cloudnet open API → RAW per-day layout) then
`calibrate_cloudnet_cl61.py` (RAW reader, WV-corrected, seven molecular-window methods).
**Fog exclusion** (CL31/CL51/CL61 `vertical_visibility` > 0 ⇒ fog ⇒ excluded) is applied in both
pipelines: MATLAB `paper_val_process.m/screen_profiles`, and Python
`calibration/io/data_loader.py/filter_cloudy_profiles` (verified on Edmonton 2026-02-01,
where 66/155 night profiles were fog).
