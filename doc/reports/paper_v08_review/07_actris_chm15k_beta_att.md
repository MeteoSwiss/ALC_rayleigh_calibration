# ACTRIS Cloudnet raw CHM15k files: the beta_att format, c_cal, and the E-PROFILE Rayleigh constant

> Session note: paths under `scratchpad/` were in the temporary session folder and are not kept; the raw and product files can be re-downloaded from the URLs in Section 5.

Check date: 2026-10-10. Sources: the Cloudnet data-portal API, one raw file per instrument
(read with netCDF4 only), the CloudnetPy and ceilopyter sources on GitHub, five Cloudnet lidar
product files, E-PROFILE L1 headers (`D:/E-PROFILE_L1_2026`) and the E-PROFILE v2.2 nightly Rayleigh
constants (`calout_v22_04`, method `rayleigh`, flag 1 or 0.5, 2025-01-01 to 2026-08-13).

## Summary

1. **Cloudnet does hold beta_att-format raw files.** 4 of the 24 CHM15k/CHM15k-x instruments
   with data since 2023 write it (global attribute `netcdf_mode = 1`). All four are still active:
   **Galati CHM210108** (firmware 1.100), **Limassol CHM230107** (1.110), **Cluj CHM230128** (1.120)
   and **Munich (LMU) CHM230134** (1.140). The other 20 write `beta_raw` (`netcdf_mode = 2`, or no
   attribute for firmware 1.040 and older). Many `beta_raw` units run firmware 1.07 to 1.15, so the
   format is a setting the operator chooses, not something the firmware forces.
2. **c_cal is the generic Lufft default in all four files:** a scalar float32 equal to 3.2e-12
   (`long_name = "calibration constant (c)"`, no units), so 1/c_cal = **3.125e11**. It is not a
   per-instrument factory calibration. The part that differs per instrument is `scaling` (cs, the
   "device-dependent system constant"). It is present in both formats and is already divided out
   inside `beta_raw` and `beta_att`.
3. **None of the four beta_att instruments is in E-PROFILE.** I read the serial of all 154 E-PROFILE
   CHM15k from their L1 headers and none matches. No E-PROFILE CHM15k lies within 30 km of Galati,
   Cluj or Limassol. The closest one to Munich (DWD München-Stadt, 2.8 km away) is a different unit
   (TUB180007, not TUB240017). **So c_cal cannot be compared with an E-PROFILE constant of the same
   instrument.**
4. **Substitute comparison.** 11 Cloudnet units are also E-PROFILE instruments, all in `beta_raw`
   format. For them, C_L,E-PROFILE / (1/3.2e-12) ranges from **0.57 to 1.67 (median 1.08)**. Over
   the last 3 months it ranges from 0.59 to 1.91 (median 1.17). Across all 154 E-PROFILE CHM15k the
   ratio has a median of **0.95** (p10 to p90: 0.57 to 1.28), and only 56 % of instruments are within
   ±20 %. The generic c_cal is therefore about right for the network mean, but for a single
   instrument it is off by a factor of up to about 2. It is not a usable instrument calibration.
5. **CloudnetPy never keeps c_cal as the calibration.** Since 2026-10-01 (CloudnetPy ≥ 1.96.7, reader
   delegated to ceilopyter), `netcdf_mode == 1` leads to `beta_raw = beta_att / c_cal`, which is then
   multiplied by Cloudnet's own calibration factor (default 3e-12). I checked this pixel by pixel:
   product / raw = **0.9375 exactly** for Cluj and Limassol. Before that, CloudnetPy read `beta_att`
   as if it were `beta_raw`, and these instruments were given a calibration factor of about 1
   (verified: 1.000 for Munich and 1.001 for Galati, products of 2026-09-15), so in practice the
   firmware c_cal was used. As a result, Cloudnet β for the four units **dropped by 6.25 %** at the
   2026-10-01 switch.

## 1. Cloudnet CHM15k inventory (API)

`/api/instrument-pids` lists 30 CHM15k-family PIDs. 24 of them have a `lidar` product since 2023;
6 have none (CHM060012 DWD, CHM080085 CNR, CHM130102 UGR, CHM210109 IAC, CHM230126 Met Éireann,
CHMstddrd UGR).

Raw uploads (`/api/raw-files?instrument=chm15k`) numbered 1440 on 2026-09-15, 1056 on 2026-10-01 and
1330 on 2026-10-08. Most of these are 5-minute files from Cabauw, Evora, Granada, Melpitz and
Payerne. 17 CHM15k uploaded on all three days, except Melpitz, which is missing on 2026-10-01. The
two CHM15k-x (`instrument=chm15kx`: UFS Schneefernerhaus, and TROPOS LACROS now at Invercargill)
uploaded on 2026-10-08. Five more instruments are recently inactive: CHM150101 (last 2026-09-08),
CHM200114 (last 2026-09-01), CHM170119 (replaced by CHM250105 on 2026-02-12), CHX090103 (last
2025-08-01) and CHM170158 (RV Meteor, 2024).

| Serial (Cloudnet) | Instrument | PID (hdl.handle.net/…) | Site(s) | First lidar day ≥2023 | Latest |
|---|---|---|---|---|---|
| CHM100110 | DWD CHM 15k (2010) | 21.12132/3.cdf99c536bd04146 | lindenberg | 2023-01-01 | 2026-10-10 |
| CHM120109 | UoC CHM 15k (2012) | 21.12132/3.324507c6ce8e489c | juelich | 2023-01-01 | 2026-10-10 |
| CHM150101 | IPSL CHM 15k (2015) | 21.12132/3.b14e15ec9c9b43fb | palaiseau | 2023-01-01 | 2026-09-08 |
| CHM150118 | KNMI CHM 15k (2015) | 21.12132/3.c07f533cd86041f1 | cabauw | 2023-09-01 | 2026-10-10 |
| CHM160146 | IFUSP CHM 15k (2016) | 21.12132/3.f3c16fc50e8a4964 | campina (ATTO) | 2023-01-01 | 2026-10-09 |
| CHM170103 | TROPOS CHM 15k (Melpitz) | 21.12132/3.75d6aaafbf1a461d | melpitz | 2026-01-29 | 2026-10-10 |
| CHM170119 | UGR CHM 15k (2017) | 21.12132/3.77a75f3b32294855 | granada | 2023-01-01 | 2026-02-12 |
| CHM170137 | INOE CHM 15k (2017) | 21.12132/3.c60c931fac9d43f0 | bucharest | 2023-01-01 | 2026-10-10 |
| CHM170158 | MPI-M CHM 15k (2017) | 21.12132/3.d20b2f985a124a7e | rv-meteor | 2024-08-17 | 2024-09-22 |
| CHM190157 | ENEA CHM 15k (2019) | 21.12132/3.a5d479cd32a14a88 | lampedusa | 2023-01-16 | 2026-10-09 |
| CHM190159 | NUIG CHM 15k | 21.12132/3.5c5d0a39d2cd403d | mace-head | 2023-01-01 | 2026-10-10 |
| CHM200110 | MeteoSwiss CHM 15k (2020) | 21.12132/3.5ffe1a5df37145cc | payerne | 2023-01-01 | 2026-10-10 |
| CHM200114 | LIM CHM 15k (CEILIM) | 21.12132/3.7d8e02584f8541a9 | falkenberg, leipzig-lim | 2023-01-01 | 2026-09-01 |
| CHM200127 | IPSL CHM 15k (2020) | 21.12132/3.9b914f8c94c4438a | palaiseau | 2026-02-05 | 2026-10-10 |
| **CHM210108** | DJUG CHM 15k (2021) | 21.12132/3.079cdfced7a7438e | **galati** | 2023-01-01 | 2026-10-10 |
| **CHM230107** | ECoE CHM 15k (2023) | 21.12132/3.d7b65c4b72b6432d | **limassol** | 2024-06-26 | 2026-10-10 |
| CHM230109 | PULS CHM 15k (2023) | 21.12132/3.eb35ee6de17b49b4 | rzecin | 2023-07-11 | 2026-10-08 |
| **CHM230128** | UBB CHM 15k (2023) | 21.12132/3.d6bf209bc48b48a4 | **cluj** | 2024-04-11 | 2026-10-09 |
| **CHM230134** | LMU CHM 15k (2023) | 21.12132/3.8584159a7b3a4545 | **munich** | 2024-12-15 | 2026-10-10 |
| CHM230136 | UÉ CHM 15k (2023) | 21.12132/3.982ea3f82c82491b | evora | 2024-07-25 | 2026-10-10 |
| CHM250105 | UGR CHM 15k (2026) | 21.12132/3.25a1c14a153b4706 | granada | 2026-02-12 | 2026-10-10 |
| CHX090103 | LMU CHM 15k (2009) | 21.12132/3.60eef0a44a7542a0 | munich | 2023-01-01 | 2025-08-01 |
| CHX090105 | UFS CHM 15k-x (2009) | 21.12132/3.c31cb1ea90be473d | schneefernerhaus | 2023-01-01 | 2026-10-09 |
| CHX090111 | TROPOS CHM 15k-x (LACROS) | 21.12132/3.cd3578e9f68b42c0 | eriswil, invercargill, leipzig | 2023-01-01 | 2026-10-10 |

## 2. Native format of each instrument (one raw file each)

`device_name` is the instrument serial (CHM…). `serlom` is the laser-optical-module serial (TUB…),
which raw2l1 writes into the E-PROFILE `instrument_serial_number`. Firmware is the third token of
`software_version`. The "Cloudnet cf" column is the calibration factor returned by
`/api/calibration` ("–" means none is configured, so the default 3e-12 is used). All files are
NETCDF3_CLASSIC.

| Device | serlom | Site | File date | Firmware | netcdf_mode | Signal variable | c_cal | Pulse norm. | scaling (cs) | Cloudnet cf | In E-PROFILE (serlom match) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **CHM210108** | TUB210031 | galati | 2026-09-15 | 1.100 | **1** | **beta_att** | **3.2e-12** | p_cal | 0.5578 | – | no |
| **CHM230107** | TUB230017 | limassol | 2026-10-01 | 1.110 | **1** | **beta_att** | **3.2e-12** | p_cal | 0.2696 | – | no |
| **CHM230128** | TUB230022 | cluj | 2026-10-08 | 1.120 | **1** | **beta_att** | **3.2e-12** | p_cal | 0.3578 | – | no |
| **CHM230134** | TUB240017 | munich | 2026-09-15 | 1.140 | **1** | **beta_att** | **3.2e-12** | p_cal | 0.1930 | – | no |
| CHM100110 | TUB190005 | lindenberg | 2026-09-15 | 1.150 | 2 | beta_raw | – | p_calc | 0.3050 | 4e-12 | 0-20000-0-10393_0 |
| CHM120109 | TUB160062 | juelich | 2026-09-15 | 1.150 | 2 | beta_raw | – | p_calc | 0.4396 | 2.3e-12 | 0-276-0-10508_A (from 2026-06-18) |
| CHM150101 | TUB140013 | palaiseau | 2026-09-08 | 1.070 | 2 | beta_raw | – | p_calc | 0.3007 | 2.3e-12 | 0-250-1001-07151_B |
| CHM200127 | TUB200032 | palaiseau | 2026-09-15 | 1.080 | 2 | beta_raw | – | p_calc | 0.3611 | – | 0-250-1001-75105006_A (QUALAIR, Paris, until 2026-02-04) |
| CHM150118 | TUB150008 | cabauw | 2026-10-08 | 1.040 | absent | beta_raw | – | p_calc | 0.2356 | – | 0-20000-0-06348_A |
| CHM170103 | TUB160061 | melpitz | 2026-10-08 | 1.040 | absent | beta_raw | – | p_calc | 0.3073 | – | 0-20008-0-MEL_A |
| CHM170119 | TUB170006 | granada | 2026-02-12 | 1.080 | 2 | beta_raw | – | p_calc | 0.4905 | 5.2e-12 | 0-20008-0-UGR_A (until 2026-02-11) |
| CHM250105 | TUB250006 | granada | 2026-10-01 | 1.150 | 2 | beta_raw | – | p_calc | 0.1829 | – | 0-20008-0-UGR_A (from 2026-02-12) |
| CHM170137 | TUB170068 | bucharest | 2026-10-01 | 1.130 | 2 | beta_raw | – | p_calc | 0.4210 | 5e-12 | 0-20008-0-INO_A |
| CHM200110 | TUB140016 | payerne | 2026-10-08 | 1.130 | 2 | beta_raw | – | p_calc | 0.3336 | – | 0-20000-0-06610_A (from 2025-04-02) |
| CHM230136 | TUB240008 | evora | 2026-10-08 | 1.140 | 2 | beta_raw | – | p_calc | 0.1877 | – | 0-620-3704-5480_B |
| CHM160146 | TUB170019 | campina | 2026-10-09 | 1.040 | absent | beta_raw | – | p_calc | 0.3233 | – | no |
| CHM170158 | TUB170065 | rv-meteor | 2024-09-22 | 1.080 | 2 | beta_raw | – | p_calc | 0.4966 | – | no |
| CHM190157 | TUB190057 | lampedusa | 2026-10-01 | 1.040 | absent | beta_raw | – | p_calc | 0.4403 | – | no |
| CHM190159 | TUB190058 | mace-head | 2026-09-15 | 0.754 | absent | beta_raw | – | p_calc | 0.6966 | 4.5e-12 | no |
| CHM200114 | TUB200014 | falkenberg | 2026-09-01 | 1.150 | 2 | beta_raw | – | p_calc | 0.3696 | – | no (E-PROFILE Falkenberg is DWD unit TUB170036) |
| CHM230109 | TUB230015 | rzecin | 2026-09-15 | 1.110 | 2 | beta_raw | – | p_calc | 0.2920 | – | no |
| CHX090103 | TUB230039 | munich | 2025-08-01 | 1.110 | 2 | beta_raw | – | p_calc | 0.2826 | – | no |
| CHX090105 | TUB080029 | schneefernerhaus | 2026-10-08 | 1.150 | 2 | beta_raw | – | p_calc | 0.1524 | 3e-12 | no |
| CHX090111 | TUB090013 | invercargill | 2026-10-08 | 1.140 | 2 | beta_raw | – | p_calc | 0.1158 | 2.897e-12 | no |

How the beta_att files are laid out (identical in all four):

- `beta_att(time, range)` has units `1/m 1/sr` and
  `long_name = "attenuated backscatter ((P_raw / lp) - b) / (cs * o(r) * p_cal) * r * r * c"`.
  The `beta_raw` long name is the same formula without `* c` and with `p_calc`. This means
  `beta_att = c_cal × beta_raw` exactly, as raw2l1 and ceilopyter both assume.
- `c_cal` is a 0-d float32 equal to 3.1999999e-12, with no units attribute.
- The variable set depends on the mode, not on the firmware. Mode 1 has `beta_att`, `beta_att_hr`,
  `c_cal`, `p_cal`, `fan` and `nn4`. Mode 2 has `beta_raw`, `beta_raw_hr`, `p_calc`, `nn2` and
  `nn3`. Mode-2 files at firmware 1.150 still have the mode-2 set.
- `wmo_id = 0` and `wigos_id = ""` in all four, so nothing links them to a WIGOS station.

## 3. c_cal versus the E-PROFILE Rayleigh constant

Definitions: C_L,manufacturer = 1/c_cal and ratio = C_L,E-PROFILE / (1/c_cal) = C_L,E-PROFILE × c_cal.
A ratio above 1 means the firmware β_att would be higher than the E-PROFILE Rayleigh-calibrated β_att
by that factor.

### 3a. The beta_att instruments (the requested comparison)

| Instrument | Device / serlom | Site | Firmware | c_cal | 1/c_cal | E-PROFILE median C_L | Ratio | n nights |
|---|---|---|---|---|---|---|---|---|
| DJUG CHM 15k | CHM210108 / TUB210031 | Galati (RO) | 1.100 | 3.2e-12 | 3.125e11 | not in E-PROFILE | – | – |
| ECoE CHM 15k | CHM230107 / TUB230017 | Limassol (CY) | 1.110 | 3.2e-12 | 3.125e11 | not in E-PROFILE | – | – |
| UBB CHM 15k | CHM230128 / TUB230022 | Cluj (RO) | 1.120 | 3.2e-12 | 3.125e11 | not in E-PROFILE | – | – |
| LMU CHM 15k | CHM230134 / TUB240017 | Munich (DE) | 1.140 | 3.2e-12 | 3.125e11 | not in E-PROFILE (DWD München-Stadt 0-20000-0-10865_0, 2.8 km, is TUB180007) | – | – |

**c_cal is equal to the Lufft default (3.2e-12) for every instrument.** It contains no
instrument-specific factory information.

### 3b. Cloudnet units that are also E-PROFILE instruments (beta_raw format): what c_cal = 3.2e-12 would give

The E-PROFILE statistics are restricted to the period in which the station carried the same serlom.
Serial-change dates come from a bisection over the daily L1 files.

| Device / serlom | Cloudnet site | E-PROFILE key | Firmware | Same-serial period | n nights | Median C_L [p25–p75] | **Ratio** | Last-3-month median C_L (n, window) | Ratio (3 mo) | Cloudnet cf | C_L × Cloudnet cf |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CHM100110 / TUB190005 | lindenberg | 0-20000-0-10393_0 | 1.150 | 2025-01-01–2026-08-13 | 212 | 4.27e11 [3.56e11–4.81e11] | **1.37** | 4.93e11 (33, 05-14–08-13) | 1.58 | 4e-12 | 1.71 |
| CHM120109 / TUB160062 | juelich | 0-276-0-10508_A | 1.150 | 2026-06-18–2026-08-13 | 27 | 2.09e11 [1.87e11–2.29e11] | **0.67** | 2.09e11 (27, same) | 0.67 | 2.3e-12 | 0.48 |
| CHM150101 / TUB140013 | palaiseau | 0-250-1001-07151_B | 1.070 | 2025-01-01–2026-08-13 | 215 | 4.69e11 [3.99e11–5.37e11] | **1.50** | 5.53e11 (44, 05-14–08-13) | 1.77 | 2.3e-12 | 1.08 |
| CHM200127 / TUB200032 | palaiseau (was QUALAIR, Paris) | 0-250-1001-75105006_A | 1.080 | 2025-01-01–2026-01-29 | 140 | 3.82e11 [3.20e11–4.46e11] | **1.22** | 2.99e11 (24, 2025-10-30–2026-01-29) | 0.96 | – (3e-12) | 1.15 |
| CHM150118 / TUB150008 | cabauw | 0-20000-0-06348_A | 1.040 | 2025-01-01–2026-08-13 | 198 | 3.52e11 [3.18e11–3.97e11] | **1.13** | 4.20e11 (32, 05-14–08-13) | 1.34 | – (3e-12) | 1.06 |
| CHM170103 / TUB160061 | melpitz | 0-20008-0-MEL_A | 1.040 | 2025-01-01–2026-08-13 | 193 | 1.77e11 [1.52e11–1.97e11] | **0.57** | 1.86e11 (37, 05-14–08-13) | 0.59 | – (3e-12) | 0.53 |
| CHM170119 / TUB170006 | granada (old unit) | 0-20008-0-UGR_A | 1.080 | 2025-01-01–2026-02-11 | 203 | 2.49e11 [2.29e11–2.85e11] | **0.80** | 2.30e11 (48, 2025-10-22–2026-01-21) | 0.73 | 5.2e-12 | 1.30 |
| CHM250105 / TUB250006 | granada (new unit) | 0-20008-0-UGR_A | 1.150 | 2026-02-12–2026-08-13 | 88 | 3.37e11 [3.10e11–3.79e11] | **1.08** | 3.64e11 (51, 05-13–08-12) | 1.17 | – (3e-12) | 1.01 |
| CHM170137 / TUB170068 | bucharest | 0-20008-0-INO_A | 1.130 | 2025-01-01–2026-08-13 | 234 | 3.38e11 [2.88e11–3.91e11] | **1.08** | 4.16e11 (43, 05-13–08-12) | 1.33 | 5e-12 | 1.69 |
| CHM200110 / TUB140016 | payerne | 0-20000-0-06610_A | 1.130 | 2025-04-02–2026-08-13 | 138 | 5.21e11 [4.42e11–6.45e11] | **1.67** | 5.97e11 (34, 05-13–08-12) | 1.91 | – (3e-12) | 1.56 |
| CHM230136 / TUB240008 | evora | 0-620-3704-5480_B | 1.140 | 2025-01-01–2026-08-13 | 333 | 3.19e11 [2.71e11–3.58e11] | **1.02** | 3.62e11 (67, 05-14–08-13) | 1.16 | – (3e-12) | 0.96 |

Network context: for all 154 E-PROFILE CHM15k with at least 10 Rayleigh nights (27,640 nights,
without serial restriction), C_L × 3.2e-12 has a median of 0.95, p10/p25/p75/p90 of
0.57/0.77/1.11/1.28 and an extreme range of 0.10 to 2.93. The share of instruments within ±10 %,
±20 %, ×1.5 and ×2 of the default is 32 %, 56 %, 79 % and 90 % respectively.

Notes:

- Serial changes found by bisection over the L1 files:
  - Payerne went from TUB200009 to TUB140016 on 2025-04-02.
  - Granada went from TUB170006 to TUB250006 on 2026-02-12, the same day the Cloudnet PID changed.
  - Jülich was TUB120017 (firmware 0.743, later 1.13) until 2025-11-18. There are no L1 files from
    then until 2026-06-18, when TUB160062 appears.
  - The SIRTA unit TUB200032 sent E-PROFILE data as QUALAIR (Paris) until 2026-02-04. It has
    appeared on Cloudnet at Palaiseau since 2026-02-05.
- The last-3-month C_L is 13–23 % above the full-period median at most sites (Lindenberg,
  Palaiseau, Cabauw, Bucharest, Payerne, Evora). That window is a summer-only sample. I did not
  investigate the cause further.
- The Cloudnet calibration factors set per instrument (last column) are not systematically closer
  to E-PROFILE than the Lufft default. Palaiseau improves from 1.50 to 1.08, but Lindenberg goes
  from 1.37 to 1.71, Bucharest from 1.08 to 1.69 and Jülich from 0.67 to 0.48. For Jülich, however,
  the factor may date from the earlier optics module.

## 4. What CloudnetPy does with beta_att / c_cal

From the source code:

- **Since 2026-10-01** the commit "Process CHM 15k using Ceilopyter" (CloudnetPy ≥ 1.96.7; 1.97.1 is
  deployed) makes `cloudnetpy/instruments/lufft.py` call `ceilopyter.read_chm15k`. In
  `ceilopyter/instruments/chm15k.py::_get_beta` (version 0.2.2; the logic has been there since
  ceilopyter 0.2.0 in 2025-06), `if nc.netcdf_mode == 1: return nc["beta_att"][:] / nc["c_cal"][:]`.
  The code comment reads "netcdf_mode is 1 (introduced in firmware version 1.050) … beta_att =
  beta_raw * c_cal", and mode 2 is described as "recommended in ACTRIS CCRES ALC SOP". The result is
  then multiplied by `calibration_factor`, which is the per-PID value from the Cloudnet calibration
  database or, if none, a default of 3e-12 with a log warning. Note that this default (3e-12)
  differs from the Lufft default (3.2e-12). c_cal is used only to undo the firmware calibration.
- **Before 2026-10-01**, from 2022-03-04 (commit "Use fallback variable name") up to 1.96.x, the
  native reader did `beta_raw = self._getvar("beta_raw", "beta_att")`. In other words, it read
  `beta_att` as if it were `beta_raw`, ignored `c_cal`, and multiplied by `calibration_factor`. A
  correct result therefore depended on the processing configuration assigning a factor of about 1
  to these instruments.

Checked on real products (Cloudnet product `beta_raw` divided by the raw `beta_att` of the same day,
with profiles matched in time, gates 50 and above):

| Product | CloudnetPy | calibration_factor in product | product / raw beta_att |
|---|---|---|---|
| 20260915_galati_chm15k_079cdfce.nc | 1.96.1 | 1.001 | 1.001 (constant) |
| 20260915_munich_chm15k_8584159a.nc | 1.96.1 | 1.0 | 1.000 (constant) |
| 20261001_limassol_chm15k_d7b65c4b.nc | 1.96.7 | 3e-12 | 0.9375 = 3e-12 / 3.2e-12 |
| 20261008_cluj_chm15k_d6bf209b.nc | 1.97.1 | 3e-12 | 0.9375 = 3e-12 / 3.2e-12 |
| 20261008_galati_chm15k_079cdfce.nc | 1.97.1 | 3e-12 | (no raw file of that day; same factor) |

Before the switch, Cloudnet β for these units equalled the firmware β_att, i.e. it used c_cal.
After the switch it is 6.25 % lower (the Cloudnet default replaces the Lufft default). Today
`/api/calibration` returns no factor for the four PIDs, even for dates before the switch. This
suggests (inferred, not confirmed) that the earlier factors of about 1 were removed when the reader
changed. raw2l1, the E-PROFILE reader, does the equivalent of ceilopyter: `beta_raw = beta_att / c_cal`
(3.2e-12 if `c_cal` is missing), and `instrument_serial_number = serlom`. If one of the four
instruments joined E-PROFILE, its `rcs_0` would therefore be on the same scale as the `beta_raw`
instruments.

## 5. Downloaded files

The raw files are in `scratchpad/actris_chm15k_files/` and the products in its `cloudnet_products/`
subfolder. All 24 raw downloads were verified against the MD5 checksum from the API. In total I
downloaded 431.3 MB of raw files and 121.1 MB of products, plus small API JSON responses
(`api/`) and source files (`cloudnetpy_src/`).

| File | Site | Device | Size (B) | URL |
|---|---|---|---|---|
| 20260915_lindenberg_CHM100110.nc | lindenberg | CHM100110 | 50,150,932 | https://cloudnet.fmi.fi/api/download/raw/6a7cc2fa-cd00-415d-a367-ff9a52b6f660/20260915_lindenberg_CHM100110.nc |
| 20260915_JOYCE_CHM120109_000.nc | juelich | CHM120109 | 25,169,832 | https://cloudnet.fmi.fi/api/download/raw/8f111976-9d65-4960-903a-714480f46f43/20260915_JOYCE_CHM120109_000.nc |
| chm15k_0a_Lz1Lb87Ppr2R15mF15s_v01_20260908_000001.nc | palaiseau | CHM150101 | 8,580,148 | https://cloudnet.fmi.fi/api/download/raw/5c2de683-494b-4f07-8cc4-81006c414ae3/chm15k_0a_Lz1Lb87Ppr2R15mF15s_v01_20260908_000001.nc |
| 06348_A202610081110_CHM15k.nc | cabauw | CHM150118 | 172,540 | https://cloudnet.fmi.fi/api/download/raw/bbf44320-78d2-4cb0-a3a6-7035f322bae2/06348_A202610081110_CHM15k.nc |
| 20261009_ATTO_CHM160146.nc | campina | CHM160146 | 9,591,864 | https://cloudnet.fmi.fi/api/download/raw/62334fcf-fe93-4851-a8bf-1d3b128d84fc/20261009_ATTO_CHM160146.nc |
| 20261008_Melpitz_CHM170103_0205_000.nc | melpitz | CHM170103 | 53,788 | https://cloudnet.fmi.fi/api/download/raw/1fcb445b-0a2c-41dc-a303-96a164de5b33/20261008_Melpitz_CHM170103_0205_000.nc |
| 20260212_Granada_CHM170119_0745_000.nc | granada | CHM170119 | 97,500 | https://cloudnet.fmi.fi/api/download/raw/3a260067-aad9-4246-a392-48402dbdd284/20260212_Granada_CHM170119_0745_000.nc |
| 20261001_Magurele_CHM170137_000.nc | bucharest | CHM170137 | 12,585,620 | https://cloudnet.fmi.fi/api/download/raw/848f60c4-9bf9-4dc9-aa3e-0e7e20820485/20261001_Magurele_CHM170137_000.nc |
| 20240922_FSMETEOR_CHM170158.nc | rv-meteor | CHM170158 | 37,749,652 | https://cloudnet.fmi.fi/api/download/raw/7f458e5f-eec3-439a-84ba-13963beb4375/20240922_FSMETEOR_CHM170158.nc |
| 20261001_LMP_CHM190157_000.nc | lampedusa | CHM190157 | 24,981,936 | https://cloudnet.fmi.fi/api/download/raw/e8b4ae8a-b4bf-42c0-92bf-c32a0c7d4648/20261001_LMP_CHM190157_000.nc |
| 20260915_MaceHead_CHM190159.nc | mace-head | CHM190159 | 25,169,780 | https://cloudnet.fmi.fi/api/download/raw/717cd8b0-41e6-401a-b83d-ed25340d6106/20260915_MaceHead_CHM190159.nc |
| 20261008_pay_CHM200110_1145_000.nc | payerne | CHM200110 | 97,512 | https://cloudnet.fmi.fi/api/download/raw/4509cb88-57e6-447a-a6ee-bf13a8cead27/20261008_pay_CHM200110_1145_000.nc |
| 20260901_Falkenberg_CHM200114_000.nc | falkenberg | CHM200114 | 25,059,108 | https://cloudnet.fmi.fi/api/download/raw/a2fb1c1a-591a-4d79-9cf8-6c815fdff680/20260901_Falkenberg_CHM200114_000.nc |
| chm15k_0a_Lz1Lb87Ppr2R15mF15sI200127_v01_20260915_000006.nc | palaiseau | CHM200127 | 25,169,828 | https://cloudnet.fmi.fi/api/download/raw/548c33ea-f71b-4d28-a79c-f8a5f49689d0/chm15k_0a_Lz1Lb87Ppr2R15mF15sI200127_v01_20260915_000006.nc |
| 20260915_UGAL-REXDAN_CHM210108_000.nc | galati | CHM210108 | 25,146,968 | https://cloudnet.fmi.fi/api/download/raw/6b46780d-b796-4795-9c89-d0be244e0e9e/20260915_UGAL-REXDAN_CHM210108_000.nc |
| 20261001_Limassol_CHM230107_000.nc | limassol | CHM230107 | 25,142,576 | https://cloudnet.fmi.fi/api/download/raw/7441be84-2e2f-4142-9d46-38202544a3e4/20261001_Limassol_CHM230107_000.nc |
| 20260915_Rzecin_CHM230109_000.nc | rzecin | CHM230109 | 37,177,444 | https://cloudnet.fmi.fi/api/download/raw/09165189-e689-4b0e-a87a-5a09c0c2b47f/20260915_Rzecin_CHM230109_000.nc |
| 20261008_Cluj_CHM230128_000.nc | cluj | CHM230128 | 22,589,620 | https://cloudnet.fmi.fi/api/download/raw/dca8eebd-7641-4f33-825a-753bedd6a756/20261008_Cluj_CHM230128_000.nc |
| CHM15kLMU_20260915.nc | munich | CHM230134 | 38,237,040 | https://cloudnet.fmi.fi/api/download/raw/a2589a8d-9160-48a5-b33a-29fb3203e971/CHM15kLMU_20260915.nc |
| 20261008_EVO_CHM230136_0710_000.nc | evora | CHM230136 | 75,652 | https://cloudnet.fmi.fi/api/download/raw/86dbf2e8-4a97-46d4-adb6-1107924ad94a/20261008_EVO_CHM230136_0710_000.nc |
| 20261001_Granada_CHM250105_2355_000.nc | granada | CHM250105 | 88,764 | https://cloudnet.fmi.fi/api/download/raw/6fc23cb6-6bab-47df-8ae3-99dd966dfaf0/20261001_Granada_CHM250105_2355_000.nc |
| CHM15kxLMU_20250801.nc | munich | CHX090103 | 13,027,960 | https://cloudnet.fmi.fi/api/download/raw/c22c5493-65ba-45f9-b513-2d75fd251958/CHM15kxLMU_20250801.nc |
| CHM15kUFS_20261008.nc | schneefernerhaus | CHX090105 | 25,170,984 | https://cloudnet.fmi.fi/api/download/raw/829afca9-0fce-468c-bc3d-ac6d46eb98ad/CHM15kUFS_20261008.nc |
| 20261008_Invercargill_CHX090111_1135_000.nc | invercargill | CHX090111 | 53,820 | https://cloudnet.fmi.fi/api/download/raw/9ccfdcce-14bb-481d-90d5-96cf7c28bb2c/20261008_Invercargill_CHX090111_1135_000.nc |
| cloudnet_products/20260915_galati_chm15k_079cdfce.nc | galati | CHM210108 | 27,632,994 | https://cloudnet.fmi.fi/api/download/product/3909d606-08cb-40d2-bd9b-6aeb83cc63f3/20260915_galati_chm15k_079cdfce.nc |
| cloudnet_products/20261008_galati_chm15k_079cdfce.nc | galati | CHM210108 | 23,496,387 | https://cloudnet.fmi.fi/api/download/product/1195223e-ae90-4014-a349-84d10b92ca47/20261008_galati_chm15k_079cdfce.nc |
| cloudnet_products/20260915_munich_chm15k_8584159a.nc | munich | CHM230134 | 25,534,443 | https://cloudnet.fmi.fi/api/download/product/8f1f7893-892d-4211-be89-64420231e009/20260915_munich_chm15k_8584159a.nc |
| cloudnet_products/20261001_limassol_chm15k_d7b65c4b.nc | limassol | CHM230107 | 23,533,902 | https://cloudnet.fmi.fi/api/download/product/d29b7c7e-467a-4169-a323-e0757d3ac458/20261001_limassol_chm15k_d7b65c4b.nc |
| cloudnet_products/20261008_cluj_chm15k_d6bf209b.nc | cluj | CHM230128 | 20,912,518 | https://cloudnet.fmi.fi/api/download/product/25a11066-09fc-4914-8303-7db5e5275f83/20261008_cluj_chm15k_d6bf209b.nc |

## 6. Method and reproducibility

The scripts are in `scratchpad/actris_check_scripts/`. All were run with
`.venv/Scripts/python.exe -I`, and none was run from the download folder.

| Script | Purpose |
|---|---|
| `query_raw_chm15k.py` / `summarise_raw_chm15k.py` | Raw-upload listings for the sample days |
| `list_chm15k_pids.py` / `latest_per_pid.py` | PID inventory and the latest day per PID |
| `select_and_download.py` | Picks one file per PID and downloads it with an MD5 check; writes `download_manifest.json` |
| `inspect_raw_files.py` | Format markers, `c_cal`, `scaling`, firmware, serials and all global attributes; writes `inspect_results.json` |
| `match_census.py`, `eprofile_side.py`, `eprofile_all_chm_serials.py`, `find_swap_dates.py` | E-PROFILE serials from L1 headers opened by exact path (no directory scan), serial-change dates and Rayleigh statistics |
| `build_comparison.py` | Builds `comparison.json` |
| `cloudnet_calibration.py`, `check_cloudnet_product.py`, `check_cloudnet_pairs.py` | Cloudnet calibration factors and the product-versus-raw check |

Nothing in the repository was modified.

Caveats:

- I inspected one raw file per instrument. c_cal and the mode could in principle have changed at
  other dates (for example after an operator reconfiguration).
- The E-PROFILE station list is the 2026 L1 census (154 CHM15k). An instrument that reported to
  E-PROFILE only before 2026 would not be seen.
- The "ratio" for the beta_raw units is hypothetical: it is what the firmware default c_cal would give
  if those units were switched to beta_att mode.
