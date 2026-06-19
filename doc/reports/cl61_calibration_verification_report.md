# CL61 calibration verification — Rayleigh vs liquid-cloud, with Kalman / water-vapor / L1–L2 sensitivity

**Author:** M. Hervo · **Date:** 2026-06-17
**Scope:** every CL61 in the E-PROFILE network (9 instruments), year 2026.

## 1. What was done and why

The paper validates calibrated attenuated backscatter (β_att). For the CL61 the two
independent absolute calibrations — the **molecular (Rayleigh)** fit and the
**liquid-cloud** (O'Connor/Hopkin) method — disagree at Payerne, with cloud running
higher. To make sure this is a real, properly-computed result and not a processing
artefact, **all nine network CL61 were re-calibrated from scratch** and the result
was stress-tested along three axes:

| Axis | Variants | How |
|---|---|---|
| Calibration method | Rayleigh (molecular) · liquid-cloud | Python pipeline · MATLAB `liquid_cloud_calibration` |
| Water-vapor correction | **on** (production) · **off** | `apply_wv_correction` 1/0 in both pipelines |
| Kalman smoothing | **with** · **without** | both pipelines store the raw daily series *and* the Kalman series |
| Data level | **L2** monthly · **L1** daily (`D:\E-PROFILE_L1_2026`) | `DataLevel.L2_MONTHLY` / `DataLevel.L1` |

All runs are restricted to **2026** (the period the L1 archive covers); the
cloud-vs-Rayleigh agreement is reported for the **Mar–May 2026 validation window**
(where the direct CHM-15k ground-truth exists), with full-2026 as a period check.
Instruments are identified by the **`instrument_type` global attribute** ("CL61"),
never by wavelength (a co-located CL51 is also ~910 nm).

The nine CL61:

| WIGOS_id | site | L2 id | L1 id (CL61) | network L1? |
|---|---|---|---|---|
| 0-20000-0-06418 | Zeebrugge (BE) | B | B | yes |
| 0-20000-0-06447 | Uccle (BE) | B | B | yes (co-located CL51 = "A") |
| 0-20000-0-06610 | **Payerne (CH)** | C | — | no — research CL61, only the operational CL31 feeds network L1 (raw at `E:\CL61_PAY`) |
| 0-20008-0-BIR | Birkenes (NO) | A | A | yes |
| 0-20008-0-EDT | Edmonton (CA) | B | B | yes (sparse) |
| 0-20008-0-LAU | Lauder (NZ) | A | A | yes |
| 0-203-10-LNG | Langenlebarn (AT) | A | A | yes |
| 0-380-5-1 | Aosta St-Christophe (IT) | B | B | yes |
| 0-756-4-EERLCL61 | (CH) | A | A | yes (2-month record) |

## 2. Bugs found and fixed during verification

Re-running from scratch exposed several real bugs (this *was* the point of the exercise):

1. **L1 reader — housekeeping name** (`data_loader.py`): read `temperature_optical_module`
   unconditionally; Vaisala CL61 L1 uses `temperature_laser`. → every L1 night crashed.
   Fixed with a manufacturer-aware fallback (HK only feeds diagnostics, never the fit).
2. **L1 reader — no-cloud sentinel** (`data_loader.py`): L1 `cloud_base_height` uses
   **−999.9 / −1000** for "no cloud"; the clear-night test expected −9.0. → every L1
   profile read as cloudy (0 valid nights). Fixed by normalizing the L1 sentinel.
3. **Cloud instrument-type detection** (`liquid_cloud_calibration.m`): the file-attribute
   detection handled `CL31`/`CL51` but **not `CL61`**. Added CL61 (and CHM15k); WV
   applicability now keys on `instrument_type`, not wavelength.
4. **Cloud WV error handling**: made **strict** — a 910 nm WV failure now raises an error
   and the period is **excluded, never silently calibrated uncorrected** (hard rule). A
   `dispo_error` typo in the read catch was also fixed.
5. **Cloud parfor WV bug** (the big one — see §5): under `parfor` the WV correction
   silently produced no effect, inflating the cloud coefficient by ~15 %. Fixed by
   running the cloud driver serially.

## 3. Headline result — cloud is higher than Rayleigh at every CL61

![CL61 network cloud vs Rayleigh, controlled re-run](figs_paper_validation/cl61_verify_cloud_vs_rayleigh.png)

*(a) Per-instrument cloud-vs-Rayleigh relative difference 100·(β_cloud/β_Rayleigh − 1)
(L2, WV-on, robust median, Mar–May 2026). (b) WV-on vs WV-off.*

**The liquid-cloud calibration is systematically higher than the molecular calibration
at every CL61 — 8/8 sites with Mar–May data** (EDT has no clear Mar–May night). The
direction is unambiguous and universal.

**Network-median offset = +25.8 %** (range +3 % Uccle … +59 % Aosta; Payerne **+21 %**).
The Payerne value is **fully consistent with the direct CHM-15k validation**
(cloud +16 % vs CHM, Rayleigh +1.7 % vs CHM → +14–19 %). With the parfor WV bug fixed
(§5), the network cloud calibration **reproduces the validation exactly** (Payerne
cloud C = 0.968, Kalman 0.987 — identical to the `fixlatlon` validation value). The
per-site spread is driven mainly by Rayleigh-night sampling (some sites have only
8–10 clear nights) and real inter-site differences, not by method irreproducibility.

> **cloud-vs-Rayleigh for CL61 = +21 % at Payerne (CHM-anchored), +26 % network median;
> the molecular/CHM calibration is the reliable absolute reference.**

## 4. Data level — L1 and L2 give the same Rayleigh calibration (~1 %)

![L1 vs L2 consistency and WV impact](figs_paper_validation/cl61_verify_l1l2_wv.png)

*(a) Rayleigh constant from L1 vs L2 (WV-on); points on the 1:1 line, per-night ratio
**1.01 ± 1 %**. (b) WV impact per method.*

Calibrating from raw **L1** (4.8 m native, `rcs_0` in V·m²) and from **L2** monthly
(reconstructed `rcs = β_att·calConst·1e-6`, ~30 m) yields the **same physical lidar
constant to ≈ 1 %** (per-night ratio 1.01, CV ≈ 1 %; median over sites 1.008). This
confirms the L2 reconstruction and the stored `calibration_constant_0` are sound and
that **the Rayleigh calibration is independent of data level**. (L1 and L2 do not always
flag the same nights as clear, so only the overlap is compared; Payerne's CL61 is not in
the network L1 archive and uses its from-raw L2.)

## 5. Water-vapor correction — ~20 % on Rayleigh AND ~−13 % on the cloud (both matter)

| Quantity | WV-off → WV-on change (network median) |
|---|---|
| **Rayleigh** lidar constant | **+20.8 %** (range +15 … +33 %) |
| **Cloud** calibration coefficient | **−12.8 %** (range −6 … −14 %) |

**Both calibrations are strongly water-vapor dependent.** The molecular fit at 3–6 km
integrates the full WV column (~+20 %); the cloud integral at a ~1 km cloud base still
sees most of the boundary-layer WV column below it (~−13 %, with T²_wv(1 km) ≈ 0.87,
T²_wv(3 km) ≈ 0.78). So the **hard rule — never calibrate/compare 910 nm without a
matching-month WV correction — is essential for both methods.**

> **This corrects an earlier (wrong) finding of "~0 % cloud WV sensitivity".** That 0 %
> was a **parfor concurrency bug**: under `parfor` the cloud WV step silently produced no
> effect (WV-on `.mat` == WV-off `.mat`, byte-identical), and additionally threw spurious
> "all-ones" failures for some Jan/Feb months. Run **serially**, the WV correction works
> correctly (verified single-month test: Payerne April C 1.19→1.03, −13 %, with the
> expected T²_wv profile), and the cloud calibration reproduces the CHM-validated value.
> The strict error-handling (bug #4 above) now guarantees a 910 nm period is never
> calibrated without a valid WV correction — only June 2026 (genuinely no CAMS file) is
> excluded.

### 5b. The cloud calibration does NOT normalize to 1064 nm internally (verified)

Checked because the (buggy) cloud coefficient had run ~17 % high, ≈ the 910→1064 Ångström
factor (1.169). **No internal 1064-normalization exists:**
- **Code:** every β-handling step in `liquid_cloud_calibration.m` is at the native
  wavelength — read (`×1e-6` unit only), WV (`÷T²_wv`), multiple-scattering (`×η`,
  0.76–0.83), `C = S_apparent/18.8`. No Ångström/1064 factor; `raw2L2` sets λ=1064 only
  for the CHM-15k.
- **Empirical:** the WV-corrected cloud gives Payerne C = 0.968 → apparent droplet
  **S = C·18.8 = 18.2 sr ≈ 18.8** (self-consistent at 910 nm). A 1064-normalization would
  force C ≈ 1.13 / S ≈ 21 sr.

The ~17 % "drift" was **the parfor WV bug, not a normalization and not irreproducibility**
— with WV applied serially the cloud calibration is reproducible (Payerne 0.968 = the
validation's 0.968). The cloud-high offset is therefore a **genuine inter-method
difference**, most plausibly in the cloud method's effective S·η assumptions for the CL61.

## 6. Kalman smoothing

**Kalman is numerically unstable for short/irregular CL61 records** — it diverged on 3 of 9:

| instrument | raw daily-median cloud C | Kalman cloud C |
|---|---|---|
| Zeebrugge 06418 | 1.13 (robust) | **13.9** (diverged — March `calibration_constant_0 = 45.22` glitch) |
| Lauder LAU | 1.09 (robust) | **2.37** (diverged) |
| EERLCL61 | 1.04 (robust) | **NaN** (2-month record) |
| other 6 | — | track the raw median (stable; e.g. Payerne raw 0.968 / Kalman 0.987) |

→ **the raw daily median is the robust estimator** (used throughout this report); the
operational Kalman should be guarded against non-physical excursions on short records.
On stable instruments Kalman and raw agree to ≲2 %, so Kalman is not a magnitude driver.

## 7. Is the cloud-high offset an artefact? (saturation / aerosol)

No — it survives every controlled check (WV, Kalman, L1/L2, 1064) and is not a detector
or aerosol artefact:

![Payerne CL61 cloud calibration diagnostics](figs_paper_validation/cl61_cloud_diagnostic.png)

*(a) Apparent droplet S = C·18.8 vs cloud-base height — peaked near 18.8 sr (internally
self-consistent). (b) Its distribution. (c) Raw peak β_att: smooth, bimodal, no ceiling/
pile-up → no saturation.*

- **Saturation: ruled out.** No β ceiling; apparent S correlates **negatively** with peak
  β (−0.48) — saturation needs a strong *positive* correlation.
- **Aerosol: ruled out (wrong sign).** Below-cloud aerosol would push the cloud coefficient
  **down**, opposite to the observed high bias; the 90 %-in-cloud filter removes it anyway.

The Rayleigh side is clean too — the molecular fit tracks the reference through an
aerosol-free 1.7–5.6 km window:

![Rayleigh molecular fit, Payerne CL61, 13 Mar 2026](figs_paper_validation/rayleigh_diag/plots/0-20000-0-06610/2026/20260313_0-20000-0-06610_molecular_fit.png)

## 8. Conclusions

1. **Cloud > Rayleigh at every CL61** (8/8 with Mar–May data). Robust and universal.
2. **Magnitude +21 % (Payerne, CHM-anchored) / +26 % network median** — consistent across
   methods once the parfor WV bug is fixed; the cloud calibration **reproduces** the
   validation (no real irreproducibility).
3. **Rayleigh is reproducible and data-level-independent** (L1 = L2 to ~1 %; agrees with
   CHM-15k +1.7 %) → the molecular/CHM calibration is the reference.
4. **Water vapor matters for BOTH methods** — Rayleigh +20 %, cloud −13 %. Never compare
   910 nm without a matching-month WV correction (now enforced strictly).
5. The cloud-high offset is **not** a saturation, aerosol, 1064-normalization, Kalman, or
   data-level artefact — it is a genuine inter-method difference (likely the cloud
   method's effective S·η for the CL61).
6. Five real bugs were found and fixed (§2), the most consequential being the **parfor WV
   bug**, which had inflated the cloud coefficient and produced a spurious "0 % cloud WV"
   and an apparent "+17 % irreproducibility".

### Recommendations / follow-up
- Use the molecular/CHM calibration as the CL61 absolute reference; treat liquid-cloud as
  a cross-check carrying the larger systematic.
- Guard the operational Kalman against non-physical excursions on short CL61 records.
- Root-cause the parfor WV statefulness if parallel cloud processing is wanted (serial is
  correct meanwhile).
- Commit the L1-reader, instrument-type, strict-WV and serial fixes.

---

### Artefacts produced
- Rayleigh re-runs: `D:\E-PROFILE_calibration_rayleigh\cl61_verify\{L2,L1}_{WVon,WVoff}\`
- Cloud re-runs (serial, WV-correct): `A:\E-PROFILE_L2_Calibration\cl61_verify_cloud\{WVon,WVoff}\`
- Master table: `figs_paper_validation\cl61_verify_summary.csv`
- Figures: `cl61_verify_cloud_vs_rayleigh.png`, `cl61_verify_l1l2_wv.png`,
  `cl61_cloud_diagnostic.png`, `rayleigh_diag\...\*_molecular_fit.png`
- Scripts: `run_cl61_variants.py`, `rerun_cl61_cloud_variants.m`, `cl61_verify_analysis.py`,
  `test_cloud_wv.m`
- Code fixes: `rayleigh_calibration/data_loader.py` (L1 HK + cbh sentinel);
  `liquid_cloud_calibration.m` (CL61 type detection, type-based + strict WV, `disp_error`);
  `run_cl61_variants.py` (instrument_type L1 id); `rerun_cl61_cloud_variants.m` (serial)
