# Dark-measurement reproduction - Payerne L1

Python port of `dark_measurement_cl61_chm_cl31.m` run on the operational L1 files for **Payerne** (0-20000-0-06610), window **2026-05-12T09:35 .. 2026-05-12T14:50 UTC** (the primary date from the MATLAB code).

> **Data caveat.** The MATLAB experiment used hood-on (covered) raw files; here we reuse the operational L1 network archive for the same dates, so the telescopes are uncovered and the low-altitude field contains real atmosphere. The per-gate noise below ~1 km is therefore an *upper bound*. Estimator (b) (temporal first difference) is used because it cancels the static atmosphere and slow variability and isolates the white detector noise; the 4-7 km floor and the whole detection-threshold methodology reproduce the MATLAB experiment.

> **Absolute-scale caveat.** CL61 L1 is already physical attenuated backscatter (m^-1 sr^-1), so its noise floor and thresholds are absolute. CL31 and CHM15k are operationally **uncalibrated**: their raw rcs_0 is divided by the default constants (CL31 1e8, CHM15k 3e11), so their beta, beta_min, M_min and ICAO altitudes scale inversely with those defaults - treat the CL31/CHM15k absolute numbers as default-constant-dependent (the MATLAB used CHM_CAL=5e11, a 5/3 shift). The relative behaviour and the methodology are unaffected.

## Instruments

| Instrument | Stream | lambda [nm] | Profiles | Gates | dt [s] | Top [m] |
| --- | --- | --- | --- | --- | --- | --- |
| CL31 | B | 910.0 | 630 | 770 | 30 | 7700 |
| CL61 | C | 910.5 | 630 | 3276 | 30 | 15720 |
| CHM15k | A | 1064.0 | 1260 | 1024 | 15 | 15345 |

## 1. Dark-window attenuated backscatter

![time-height](dark_measurement/figs/01_timeheight.png)

## 2. Noise profiles

![profiles](dark_measurement/figs/02_profiles.png)

![estimator-b](dark_measurement/figs/03_estimator_b_noise.png)

### Estimator (b) noise sigma(beta_att) at probe altitudes [Mm^-1 sr^-1]

*In this uncovered-L1 reproduction the 500-2000 m rows remain **upper bounds** - the first difference does not fully cancel fast atmospheric variability between consecutive profiles. The **4-7 km mean** is the representative detector-noise floor.*

| Altitude [m] | CL31 | CL61 | CHM15k |
| --- | --- | --- | --- |
| 500 | 0.04193 | 0.004721 | 0.03579 |
| 1000 | 0.1887 | 0.01803 | 0.1079 |
| 2000 | 0.7863 | 0.07637 | 0.3926 |
| 3000 | 1.531 | 0.1745 | 0.8256 |
| 5000 | 4.393 | 0.483 | 2.256 |
| mean 4-7 km | 5.47 | 0.5731 | 2.924 |

## 3. Allan deviation (tau scaling)

![allan](dark_measurement/figs/04_allan.png)

## 4. Non-range-corrected signal

![rawsignal](dark_measurement/figs/05_rawsignal.png)

## 5. Aerosol detection thresholds

![detect-backscatter](dark_measurement/figs/06_detect_backscatter.png)

![detect-extinction-mass](dark_measurement/figs/07_detect_extinction_mass.png)

### Detection thresholds at SNR=3, Ash: LR=60 sr, MEC=0.60 m^2/g

**tau = 30min** (beta in Mm^-1 sr^-1, alpha in Mm^-1, M in ug/m^3)

| Instrument | z [m] | beta_min | alpha_min | M_min |
| --- | --- | --- | --- | --- |
| CL31 | 500 | 1.626e-02 | 9.754e-01 | 1.6 |
| CL31 | 1000 | 7.322e-02 | 4.393e+00 | 7.3 |
| CL31 | 2000 | 3.056e-01 | 1.834e+01 | 30.6 |
| CL31 | 3000 | 5.958e-01 | 3.575e+01 | 59.6 |
| CL31 | 5000 | 1.714e+00 | 1.028e+02 | 171.4 |
| CL61 | 500 | 1.829e-03 | 1.097e-01 | 0.2 |
| CL61 | 1000 | 6.994e-03 | 4.196e-01 | 0.7 |
| CL61 | 2000 | 2.967e-02 | 1.780e+00 | 3.0 |
| CL61 | 3000 | 6.787e-02 | 4.072e+00 | 6.8 |
| CL61 | 5000 | 1.884e-01 | 1.130e+01 | 18.8 |
| CHM15k | 500 | 9.806e-03 | 5.884e-01 | 1.0 |
| CHM15k | 1000 | 2.959e-02 | 1.775e+00 | 3.0 |
| CHM15k | 2000 | 1.077e-01 | 6.462e+00 | 10.8 |
| CHM15k | 3000 | 2.267e-01 | 1.360e+01 | 22.7 |
| CHM15k | 5000 | 6.200e-01 | 3.720e+01 | 62.0 |

### ICAO-threshold detection altitude (night, tau=30min, SNR=3, Ash: LR=60 sr, MEC=0.60 m^2/g)

| Instrument | 200 ug/m3 | 2000 ug/m3 | 4000 ug/m3 |
| --- | --- | --- | --- |
| CL31 | 5500 m | 7700 m | 7700 m |
| CL61 | 15710 m | 15720 m | 15720 m |
| CHM15k | 8841 m | 15345 m | 15345 m |

*This per-station ICAO detection altitude is the headline scalar that will populate the network map on the dashboard.*
