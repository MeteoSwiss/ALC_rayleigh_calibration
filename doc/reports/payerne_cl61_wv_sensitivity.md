# Sensitivity of the Payerne CL61 Rayleigh calibration to the water-vapour correction

*M. Hervo (MeteoSwiss), 2026-06-16. Companion to the β_att validation
([attbsc_validation_technical.md](attbsc_validation_technical.md)).*

## Question

Does the water-vapour (WV) absorption correction **in the CL61 Rayleigh calibration**
matter, and by how much? The CL61 emits at 910.74 nm, inside the H₂O absorption band, so the
high molecular-reference region used by the Rayleigh calibration is itself attenuated by WV. If
that attenuation is not removed, the retrieved lidar constant is biased low (Wiegner &
Gasteiger, 2015).

## Method — calibration run twice

The corrected-Python Rayleigh lidar constant for the Payerne CL61 (`0-20000-0-06610_C`,
46.8137° N 6.9425° E, Mar–May 2026, 30 m grid) was computed **twice**:

- **WV calib** — RCS divided by the two-way WV transmission T²_wv before the lidar-constant fit
  (`apply_wv_correction = 1`, the production setting);
- **no-WV calib** — same fit without that division (`apply_wv_correction = 0`).

**Everything else is identical**: same 30 m CL61 L2, same Kalman filter, same display-side WV
correction and Ångström scaling, same CHM15k (Rayleigh, 1064 nm) reference, same screening and
30 m comparison grid. The two CL61 channels therefore differ **only in the calibration
constant**, isolating the effect of the WV correction *in the calibration*.

## Results

**Calibration constant** (Kalman median over the period):

| | lidar constant C | ratio |
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

![Payerne CL61 Rayleigh WV-in-calibration sensitivity (Mar–May 2026): (a) example profile, (b) scatter vs CHM15k, (c–e) time-height attenuated backscatter for CHM15k, CL61 (WV calib) and CL61 (no-WV calib). The no-WV-calibrated CL61 (e) is systematically brighter/higher than the WV-calibrated one (d). Black dots = cloud-base detections.](figs_paper_validation/sensitivity_payerne_wv.png)

## Interpretation

- **The WV correction in the calibration is essential: a ≈ 16–18 % effect** (+1.7 % → +19.3 %, a
  17.6-percentage-point validation swing; 1/T²_wv ≈ 16 % on the constant). Without it the CL61
  looks ≈ +19 % *too high* relative to the 1064 nm reference, whereas the correctly WV-calibrated
  CL61 reads **+1.7 %** — essentially unbiased, the value expected from the clean
  molecular-vs-molecular (CHM Rayleigh ↔ CL61 Rayleigh) comparison.
- **The correlation is unchanged (r ≈ 0.987).** The calibration constant only rescales the
  profile; the range-dependent WV shape is removed by the display-side correction in *both* cases,
  so r is insensitive to the calibration-WV toggle. The WV-in-calibration effect is purely a
  ≈ 16 % **scale** (bias) effect.
- The magnitude (≈ 16–18 %) is consistent with Wiegner & Gasteiger (2015), who report a ~20 %
  backscatter bias at mid-latitudes when WV absorption is ignored at ~910 nm. *(An earlier run with
  the CL61 coordinates left at 0,0 over-estimated this at ≈ 25 %: the Gulf-of-Guinea CAMS column is
  far moister than Payerne — see the lat/lon note below.)*

**Conclusion.** The water-vapour correction must be applied in the Rayleigh calibration of
910 nm ALCs; omitting it biases the CL61 lidar constant by ≈ 16 % and the validated attenuated
backscatter accordingly (+1.7 % → +19.3 % vs the 1064 nm reference). The production dataset uses it
(`apply_wv_correction = 1`).

---
*Reproduce:* `sensitivity_payerne_wv.m`. The WV-calibrated channel uses the correct-coordinates CSV
`…\fullcal_stdmol_check`; the WV-off channel uses `…\fullcal_all_noWV` (no CAMS →
lat/lon-independent). Both use the Payerne CL61 L2 patched to the real coordinates (the raw files
report 0,0 — see [payerne_cl61_calibration_sensitivity.md](payerne_cl61_calibration_sensitivity.md)
for the lat/lon fix). The WV-off CSV is a **diagnostic only**; the operational chain keeps
`apply_wv_correction = 1`.
