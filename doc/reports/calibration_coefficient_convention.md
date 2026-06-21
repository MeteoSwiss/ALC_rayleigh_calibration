# Calibration coefficient — definition and convention

*Single source of truth for how the calibration coefficient is named, defined and reported
across this package (code, output files, figures, reports). Referenced by the other reports.*

## The one convention: Wiegner lidar constant `C_L`

Both calibration methods report the **Wiegner & Geiß (2012) lidar constant**

```
C_L = RCS / β_att
```

where `RCS` is the range-corrected signal (`P·r²`, background-subtracted) and `β_att` is the
**attenuated total backscatter** (`m⁻¹ sr⁻¹`). Equivalently, the lidar equation is
`RCS = C_L · β_att`, so the data are calibrated by

```
β_att = RCS / C_L .
```

`C_L` is a large, instrument-specific number. Typical 2026 values: **CHM15k ≈ 5×10¹¹**,
**CL61 ≈ 1**, **Mini-MPL ≈ 2×10⁵**. Units: `V·m³/sr` (Vaisala CL31/CL51/CL61) or
`counts/s·m³/sr` (CHM15k, Mini-MPL).

This is the **primary reported quantity for every method**, so the Rayleigh and liquid-cloud
products live on the **same axis** and can be plotted in one time series.

## Rayleigh (molecular) calibration — Wiegner & Geiß (2012)

`C_L` is computed directly from the raw range-corrected signal in the molecular window:
`C_L = RCS / β_att,molecular` (transmission-corrected). **No prior calibration is involved**,
so the absolute `C_L` is obtained directly and is **identical from L1 and L2** (the L2 reader
reverts the applied calibration and corrects units, `RCS = attenuated_backscatter_0 ×
calibration_constant_0 × 1e-6`, recovering the L1-scale signal).

Code: `CalibrationResult.lidar_constant`. NetCDF: variable `lidar_constant`
(`long_name = "Lidar constant C_L (Wiegner & Geiss 2012)"`).

## Liquid-cloud calibration — O'Connor (2004) / Hopkin (2019)

The cloud method ingests the file's **already-calibrated** attenuated backscatter `β_att,file`
and integrates it through a fully-attenuating liquid cloud, comparing to the theoretical value
`1/(2·S)` with `S = 18.8 sr`. Its native output is the **O'Connor multiplier**

```
C = S_apparent / 18.8 = β_true / β_att,file
```

`C` is therefore the **inverse sense** of `C_L` (`C ∝ β/signal`, `C_L ∝ signal/β`): the two
conventions are reciprocals, `C_r = 1/C_c` in the user's shorthand.

`C ≈ 1` only when the file's `calibration_constant_0` is already the true lidar constant. Across
the 2026 network this splits by type (sampled, 8 streams/type):

- **CL31 / CL51 L2 store a fixed nominal placeholder `calibration_constant_0 = 1×10⁸`** (these
  cannot be Rayleigh-calibrated). So `C` is far from 1 (e.g. CL31 `C ≈ 1.6×10⁻⁶`), and
  `C_L = calibration_constant_0 / C` is their **only** absolute calibration.
- **CL61 L2 carries a real per-stream constant near 1** (observed 0.91–1.10), so for CL61
  `C ≈ 1` and `C_L ≈ calibration_constant_0 × (1/C)` ≈ the operational value times the cloud
  correction.

Either way the absolute `C_L` is the headline; a `C` far from 1 is **expected** for CL31/CL51,
not an error.

To report in the Wiegner convention, note the file's `β_att,file` was made with the applied
constant `calibration_constant_0` (= `RCS/β_att,file` = the operationally-applied `C_L`). Hence
the **absolute Wiegner constant from the cloud method** is

```
C_L = calibration_constant_0 / C .
```

This is the **headline cloud product** (directly comparable to Rayleigh). The dimensionless
**inverse** `1/C` (the Wiegner-sense correction factor, ≈ 1) is reported alongside; the raw
O'Connor multiplier `C` is kept only as an internal diagnostic.

> **Why the cloud `C_L` carries the L2 constant.** The cloud method only *measures a ratio*
> (`1/C`, completely independent of any prior calibration). Turning that ratio into an absolute
> `C_L` borrows the file's existing scale `calibration_constant_0`: it is exact given the file,
> but its absolute accuracy is only as good as the operational constant already applied. The
> Rayleigh `C_L` has no such dependence because it works on the raw signal.

Code: `CloudCalResults` exposes, in order of preference —

| field | symbol | meaning |
|---|---|---|
| `lidar_constant` | `C_L` | absolute Wiegner constant `= calibration_constant_0 / C` (headline; NaN if the file has no applied constant, e.g. raw L1) |
| `calibration_factor` | `1/C` | dimensionless Wiegner-sense correction — *the inverse, reported whenever possible* |
| `calibration_coefficient` | `C` | O'Connor multiplier (internal/diagnostic; ≈1 only if the file constant is the true one) |

## Variability metric is convention-independent

The short-term precision `σ_SD` (robust successive-difference, % of the median) is a **relative**
metric, invariant under `C → 1/C` and under any constant scale. So `σ_SD(C) = σ_SD(1/C) =
σ_SD(C_L)`: every variability figure/table is unchanged by the convention and is labelled in
terms of `C_L`.

## References

- M. Wiegner and A. Geiß, *Aerosol profiling with the Jenoptik ceilometer CHM15kx*,
  Atmos. Meas. Tech. 5, 1953–1964 (2012).
- E. J. O'Connor, A. J. Illingworth, R. J. Hogan, *A technique for autocalibration of cloud
  lidar*, J. Atmos. Oceanic Technol. 21, 777–786 (2004).
- E. Hopkin et al., *A robust automated technique for operational calibration of ceilometers…*,
  Atmos. Meas. Tech. 12, 4131–4147 (2019).
