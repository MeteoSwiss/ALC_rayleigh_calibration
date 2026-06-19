# Rayleigh Calibration for Automated Lidars and Ceilometers

A Python package for Rayleigh (molecular) calibration of ceilometers and lidars used in the MeteoSwiss / E-PROFILE network. Includes water-vapour correction for 910 nm instruments, six selectable molecular-window detection methods, liquid-cloud calibration, and a suite of analysis / validation scripts for the E-PROFILE ALC paper.

## Features

- **Multi-instrument support**: CHM15k, CHM8k, CL51, CL61, Mini-MPL, and more
- **Water-vapour correction** (mandatory for 910 nm): spectral two-way WV absorption from CAMS monthly means + HITRAN cross-section LUT; nights without matching CAMS are excluded (flag −4), never silently skipped
- **Six molecular-window detection methods** (selectable via `molecular_method` in options.json):
  - `main` = E-PROF v1.1 — original E-PROFILE method
  - `improved` = E-PROF v1.2 — aerosol-robust grid search (production default)
  - `optimal` = E-PROF v2 — temporal aerosol rejection + layer flagging
  - `calipso`, `earlinet`, `bellini` — alternative reference methods
- **Liquid-cloud calibration** (`cloud_calibration/`): Python port of MATLAB O'Connor method, bit-for-bit validated; reads E-PROFILE L2 or Cloudnet CL61 raw; mandatory WV correction
- **Kalman smoothing bridge**: shells out to `run_kalman_from_matlab.py` to produce a smoothed daily lidar constant
- **Three input data levels**: L1 (`rcs_0`), L2 daily, L2 monthly (see *Input data level*)
- **CF-compliant output**: CSV (`<WMO>_<id>_cl.csv`) + optional NetCDF4

## Critical operational constraints

- `apply_wv_correction` **must stay = 1** in `options.json` (production). Never calibrate 910 nm without a valid matching-month CAMS file — no fallback.
- `molecular_method` defaults to `"improved"` (E-PROF v1.2) in production. Do not change without re-running the long-run validation.
- `molecular_source` stays `"standard"` (US Standard Atmosphere 1976) for production.

## Installation

```bash
pip install -e .
```

Diagnostic plots need matplotlib (optional): `pip install matplotlib`.

**Required:** Python ≥ 3.9, NumPy, SciPy, netCDF4, pandas  
**Optional:** matplotlib (plots), scipy.ndimage (smoothing)

## Quick Start

```python
from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from rayleigh_calibration.config import InstrumentType

options = CalibrationOptions.from_json("options.json")
info = InstrumentInfo(site_name="Payerne", wmo_id="0-20000-0-06610",
                      identifier="A", instrument_type=InstrumentType.CHM15k,
                      latitude=46.82, longitude=6.95, altitude=491)
calibrate_rayleigh("20240115", info, options)
```

## Configuration

### Key options.json fields

| Field | Default | Description |
|---|---|---|
| `molecular_method` | `"improved"` | Molecular-window detection method (see Methods section) |
| `apply_wv_correction` | `1` | **Must stay 1** — WV correction for 910 nm |
| `molecular_source` | `"standard"` | Atmospheric profile source (`"standard"` or `"cams"`) |
| `cams_folder` | `"D:/CAMS/"` | Path to CAMS monthly NetCDF files |
| `abs_cs_lookup_table` | (path) | HITRAN WV cross-section LUT |
| `data_level` | `"L1"` | Input data level (`L1`, `L2_daily`, `L2_monthly`) |
| `folder_root` | — | Root of input data archive |
| `hour_min` / `hour_max` | 20 / 4 | Nighttime window (UTC) |
| `z_low_cloud` | 4000 | Max cloud height for clear-night selection (m) |
| `LRaer` | 52 | Aerosol lidar ratio assumed in molecular window (sr) |
| `threshold_quality` | 15 | Max rel. error (%) to accept a calibration |

### Input data level

`data_level` selects which product is read; `folder_root` must point to the matching archive:

| `data_level` | Variable | Layout | Example root |
|---|---|---|---|
| `L1` (default) | `rcs_0` | `<WMO>/YYYY/MM/L1_<WMO>_<id><YYYYMMDD>.nc` | `D:/E-PROFILE_L1` |
| `L2_daily` | `attenuated_backscatter_0` + `calibration_constant_0` | `<WMO>/YYYY/MM/L2_<WMO>_<id><YYYYMMDD>.nc` | `D:/E-PROFILE_L2_2021-2025` |
| `L2_monthly` | `attenuated_backscatter_0` + `calibration_constant_0` | `<WMO>/YYYY/L2_<WMO>_<id><YYYYMM>.nc` | `A:/E-PROFILE_L2_monthly` |
| `RAW` | instrument-native (CL61 `beta_att`, CHM `beta_raw`) | `<WMO>/YYYYMMDD/*.nc` or `<WMO>/YYYYMMDD.nc` | `R:/CL61/RAW_cloudnet_dl` |

For L2 levels: `rcs = attenuated_backscatter_0 × calibration_constant_0 × 1e-6`. The `1e-6` matches the MATLAB reference; stored β_att is always micro-scaled (`1E-6/(m·sr)`) regardless of the `units` attribute. Height variable can be `range`, `height`, or `altitude` (ASL — subtracted from `station_altitude` to get AGL).

## Molecular window methods (E-PROF versioning)

| Key | Label | Description |
|---|---|---|
| `main` | E-PROF v1.1 | Original E-PROFILE grid-search method |
| `improved` | E-PROF v1.2 | Aerosol-robust grid search; production default |
| `optimal` | E-PROF v2 | Temporal aerosol rejection + per-layer flagging; highest precision, lower yield |
| `calipso` | CALIPSO-type | Fixed high-altitude molecular window |
| `earlinet` | EARLINET/SCC-type | SCC-style window search |
| `bellini` | Bellini/ALICENET | ALICENET gradient-based window |
| `eprof_v10` | E-PROF v1.0 | Pre-a4e7140 baseline (sign-error in Rayleigh slope); for historical comparison only — **not selectable in production** |

`METHODS_DISPLAY` in `validation/compare_molecular_methods.py` is the ordered tuple used for display/aggregation (includes `eprof_v10`). `METHODS` in `rayleigh_calibration/rayleigh/molecular_methods.py` lists only the live-selectable methods.

**Long-run results (14 sites, full archive):** `optimal` is most precise (robust CV 10–14 %); `improved` best balances precision and yield; `main` has the highest yield (57 %) but noisiest nights; `calipso` highest yield (79 %) but weakest signal-selectivity. See the reports in `doc/reports/` (figures under `C:\DATA\Projects\202606_E-PROFILE_calibration\figs_paper_validation\molecular_methods_longrun\`).

## Algorithm overview

1. **Data loading**: L1 `rcs_0`, or L2 reconstructed, or RAW instrument signal
2. **Fog exclusion**: profiles with `vertical_visibility` set → excluded (ceilometer fog flag)
3. **Nighttime selection**: 20:00–04:00 UTC (configurable)
4. **Clear-night filtering**: remove profiles with low clouds or precipitation
5. **Water-vapour correction** (910 nm only): compute two-way T²_wv(r) from CAMS monthly q/T/lnsp + HITRAN LUT (`water_vapor_correction/water_vapor.py`); divide RCS; skip night if no CAMS (flag −4)
6. **Atmospheric model**: US Standard Atmosphere 1976 (or CAMS for molecular density)
7. **Molecular calculation**: Bucholtz (1995) β_mol(λ, z)
8. **Molecular window detection**: selected method (see Methods section)
9. **Rayleigh fit**: linear regression signal vs β_mol in the window
10. **Validation**: slope method vs Klett cross-check
11. **Output**: append to `<WMO>_<id>_cl.csv` (date, flag, lidar_constant, uncertainty, …)

## Flag meanings

| Flag | Meaning |
|------|---------|
| 1 | Successful calibration |
| 0.5 | Partially clear night (some clouds removed) |
| 0 | No data available |
| −1 | Not a clear night |
| −2 | Signal not proportional to molecular scattering (fog night) |
| −3 | Poor agreement between calibration methods |
| −4 | Missing WV data (CAMS); 910 nm night skipped — **no fallback** |
| −5 | RCS contains only NaN values |
| −6 | Uncertainty exceeds calibration value |
| −7 | Negative Rayleigh fit slope |
| −8 | Rayleigh fit intercept exceeds slope |

## Repository layout

```
rayleigh_calibration/         installable package (public API unchanged)
├── config.py  plotting.py  main.py    shared config, plotting, CLI entry point
├── rayleigh/                molecular calibration: calibration.py, rayleigh_fit.py,
│                            atmosphere.py, molecular_methods.py
├── io/                      data_loader.py (L1/L2/RAW readers) + output.py (CSV/NetCDF)
├── water_vapor_correction/  water_vapor.py (spectral two-way WV transmission)
├── cloud_calibration/       cloud_calibration.py (liquid-cloud O'Connor method)
└── data/                    standard_atmosphere_US_1976_50km.csv (package data)
scripts/        runnable calibration runners + data/ (acquisition) + doc_tools/
validation/     comparison / sensitivity / reproduction studies  (validation/README.md)
examples/       Jupyter quick-start notebooks  (examples/README.md)
tests/          pytest suite
doc/            reference papers, method descriptions, reports/  (doc/README.md)
lost_and_found/ one-off R&D scratch + old logs (unmaintained)
```

`options.json` and `instruments.json` stay at the repo root (read from the working
directory — run scripts from the repo root). Figures and result archives are written
**outside** the repo to `C:\DATA\Projects\202606_E-PROFILE_calibration`.

## Key source files

### Core package (`rayleigh_calibration/`)

| Module | Purpose |
|---|---|
| `rayleigh/calibration.py` | Main calibration pipeline; calls WV correction, molecular calc, window detection |
| `rayleigh/molecular_methods.py` | Six molecular-window detection methods |
| `rayleigh/atmosphere.py` | Standard-atmosphere loader + Bucholtz molecular properties |
| `rayleigh/rayleigh_fit.py` | Molecular-window grid search + lidar-constant fit |
| `water_vapor_correction/water_vapor.py` | WV two-way transmission from CAMS; port of MATLAB `wv_t2eff` + `compute_wv_transmission` |
| `cloud_calibration/cloud_calibration.py` | Liquid-cloud calibration (O'Connor method); bit-for-bit vs MATLAB |
| `io/data_loader.py` | Multi-level reader (L1/L2/RAW/Cloudnet CL61); fog-flag threading |
| `io/output.py` | CF-compliant CSV + NetCDF writers |

The public API is re-exported from the top level — `from rayleigh_calibration import
calibrate_rayleigh, CalibrationOptions, InstrumentInfo, ...` is unchanged.

### Scripts & studies

| Location | Contents |
|---|---|
| `scripts/` | Calibration runners: `run_all_l2monthly.py`, `run_calibration.py`, `run_targeted_rayleigh.py`, `calibrate_cloudnet_cl61.py`, `run_lindenberg_cl61_cal.py` |
| `scripts/data/` | Data acquisition & prep: Cloudnet downloads, `build_l2_manifest.py`, `build_cop_lookup.py`, converters |
| `validation/` | Method comparison, MATLAB↔Python parity, long-run & sensitivity studies (see `validation/README.md`) |
| `examples/` | Quick-start notebooks (see `examples/README.md`) |

### Test suite (`tests/`)

| File | Coverage |
|---|---|
| `test_water_vapor.py` | WV vs MATLAB refs + atmoslib; validates geopotential, level-order invariance |
| `test_cloud_calibration_vs_matlab.py` | Cloud calib bit-for-bit vs MATLAB (max rel diff 1.4e-5) |
| `test_calibration.py`, `test_methods_smoke.py` | Config / atmosphere unit tests; molecular-methods smoke test |

## Output

Primary output: CSV at `C:\DATA\Projects\202606_E-PROFILE_calibration\E-PROFILE_calibration_rayleigh\fullcal_all\<WMO>_<id>\<WMO>_<id>_cl.csv` (the output root is configured per script and via `options.folder_output`)

| Column | Description |
|---|---|
| `date` | Night date (YYYYMMDD) |
| `flag` | Success flag (1 = ok, 0.5 = partial, negative = failure) |
| `lidar_constant` | Daily lidar constant C_L |
| `uncertainty` | Calibration uncertainty (same units as C_L) |
| `bottom_height` | Bottom of molecular window (m AGL) |
| `top_height` | Top of molecular window (m AGL) |
| `message` | Human-readable status message |

## References

- Bucholtz, A. (1995). Rayleigh-scattering calculations for the terrestrial atmosphere. *Applied Optics*, 34(15), 2765–2773.
- Wiegner, M., & Geiß, A. (2012). Aerosol profiling with the Jenoptik ceilometer CHM15kx. *AMT*, 5(8), 1953–1964.
- Wiegner, M., & Gasteiger, J. (2015). Correction of water vapor absorption for aerosol remote sensing with ceilometers. *AMT*, 8(9), 3971–3984. [Spectral WV correction basis]
- Hopkin, E., et al. (2019). A robust automated technique for operational calibration of ceilometers using the integrated backscatter from overcast stratocumulus. *AMT*, 12(7), 4131–4147. [Cloud calibration O'Connor method]
- E-PROFILE Programme: https://e-profile.eu

## Changelog

### 2026 — repository reorganization

- Package split into themed subpackages: `rayleigh/` (calibration, fit, atmosphere, molecular methods), `io/` (readers + writers), `water_vapor_correction/`, `cloud_calibration/`; shared `config.py` / `plotting.py` / `main.py` stay at the package top level. The public API (`from rayleigh_calibration import …`) is unchanged.
- US Standard Atmosphere 1976 table shipped as package data (`rayleigh_calibration/data/`); `load_standard_atmosphere(None, grid)` resolves it regardless of the working directory.
- Loose scripts sorted into `scripts/` (runners + `data/`), `validation/` (studies), `examples/` (notebooks) and `lost_and_found/` (R&D scratch + logs); root `test_*.py` consolidated into `tests/`.
- Figures and result archives moved out of the repo to `C:\DATA\Projects\202606_E-PROFILE_calibration`; hardcoded output paths in the scripts updated accordingly.

### 2026 (this session)

- **WV correction** added (`water_vapor.py`): spectral two-way T²_wv from CAMS + HITRAN LUT; mandatory for 910 nm (flag −4 on missing CAMS, no fallback). Effect: Payerne CL61 Feb-28 C_L +20 %.
- **Six molecular methods** (`molecular_methods.py`): `main`/`improved`/`optimal`/`calipso`/`earlinet`/`bellini`; selectable via `molecular_method` in options.json.
- **E-PROF versioning**: `main`→v1.1, `improved`→v1.2, `optimal`→v2; historical baseline `eprof_v10` (pre-sign-error) tracked in `METHODS_DISPLAY` for paper comparison.
- **Cloud calibration port** (`cloud_calibration.py`): bit-for-bit vs MATLAB O'Connor; reads E-PROFILE L2 and Cloudnet CL61 RAW; full WV correction included.
- **Fog exclusion**: `vertical_visibility` threaded through all readers; flagged profiles excluded from Rayleigh + cloud calibration.
- **RAW data level**: native instrument reader for CHM15k β_raw and Cloudnet CL61 β_att day-folders/daily files.
- **Cloudnet CL61 pipeline**: download → homogenize → calibrate scripts for Lindenberg and Hyytiala (ACTRIS-Cloudnet API).
- **Long-run validation**: 14 sites (10 CHM15k + 4 Mini-MPL), full archive, 7 methods; results in `figs_paper_validation/molecular_methods_longrun/`.
- **Sign-error fix** (commit a4e7140): corrected sign convention in Rayleigh slope; all post-fix CSVs in `fullcal_all/`.

### v2.0.0 (2024)

- Complete rewrite with modern Python practices; type hints, dataclasses, vectorised calculations.

### v1.0.0 (2015–2024)

- Original E-PROFILE implementation by hem.
