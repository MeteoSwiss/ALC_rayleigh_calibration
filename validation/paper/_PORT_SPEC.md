I'll now compile a comprehensive, precise port specification document based on all the MATLAB code I've read. This will be detailed enough to reproduce the numbers to within a few percent.

---

## MATLAB CEILOMETER/LIDAR VALIDATION PIPELINE — PRECISE PORT SPECIFICATION

This document specifies **every formula, threshold, unit, NaN-handling choice, file path, and variable** needed to port the validation pipeline to Python.

---

### 1. CORE PIPELINE: `paper_val_process.m`

#### 1.1 Purpose
Multi-channel inter-comparison of attenuated-backscatter lidar/ceilometer profiles. Applies calibration corrections, water-vapor absorption, wavelength normalization, quality screening, and computes bias/RMSE/correlation statistics versus a reference channel.

#### 1.2 Configuration (Defaults in cfg struct)

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `WIGOS_ID` | string | required | E-PROFILE station ID (e.g., '0-20000-0-06610') |
| `startDate`, `endDate` | datetime | required | Date range |
| `channels` | struct array | required | See §1.3 |
| `referenceChannel` | int | 1 | Channel index (1-based) for stats baseline |
| `folderL2` | string | 'A:\E-PROFILE_L2_monthly\' | E-PROFILE monthly L2 root |
| `folderCloudKalman` | string | 'A:\E-PROFILE_L2_Calibration\Cloud_Trans-cor_WV-cor\' | Cloud calibration files |
| `lambda_target` | float | 1064 | Target wavelength [nm] for normalization |
| `angstromExponent` | float | 1 | Ångström exponent α (532/1064 @ high RH = 1) |
| `zMin`, `zMax` | float | 500, 3000 | Altitude band for statistics [m AGL] |
| `dtBin` | duration | minutes(60) | Temporal averaging bin size |
| `nightOnly` | bool | false | Filter to night profiles only |
| `szaThreshold` | float | 95 | Solar zenith angle threshold [deg] for night |
| `snrThreshold` | float | -1 | SNR threshold [1] (≤0 disables); see §1.4.6 |
| `wvcfg` | struct | empty | Fields: `.abs_cs_lookup_table`, `.cams_folder` |
| `rayleighOpts` | cell | {} | Name-value pairs for `load_rayleigh_python_kalman` |

#### 1.3 Channel Configuration (channels(k) struct)

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `source` | string | required | 'eprofile' \| 'earlinet' |
| `identifier` | string | '' | Instrument ID (e.g., 'A', 'B', 'C') for E-PROFILE |
| `earlinetStation` | string | '' | Station code ('sir', 'lei', 'cbw', 'ino', 'ari') for EARLINET |
| `calib` | string | required | 'rayleigh_python' \| 'cloud' \| 'none' |
| `label` | string | required | Display name (e.g., 'CL61 (Rayleigh)') |
| `applyWV` | bool | auto | Auto: `true` if 900 ≤ wavelength ≤ 920 nm |
| `bgCorrection` | bool | false | Apply instrument background correction |
| `bgMinRangeAGL` | float | 2400 | [m] minimum range for background estimation |
| `wvWavelength` | float | auto | Override WV laser wavelength [nm] |
| `wvFwhm` | float | auto | Override WV laser FWHM [nm] |
| `wavelengthModel` | string | 'angstrom' | 'angstrom' \| 'molaer' (two-component) |

#### 1.4 Processing Pipeline

**Step 1: Load L2 Data (lines 93–103)**

For each channel:
- **E-PROFILE**: `read_L2_monthly(WIGOS_ID, identifier, startDate, endDate, folderL2)`
  - Files: `{folderL2}/{WIGOS_ID}/{yyyy}/L2_{WIGOS_ID}_{identifier}{yyyymm}.nc`
  - Fields: `time`, `altitude` [m ASL], `attenuated_backscatter_0` [Mm⁻¹ sr⁻¹], `calibration_constant_0`, `quality_flag`, `cloud_base_height`, `vertical_visibility`, `l0_wavelength` [nm], `station_altitude` [m], `station_latitude/longitude` [deg], `instrument_type` (string), `site_location` (string)
- **EARLINET**: `read_earlinet_att_backscatter(station, startDate, endDate)` (see §3)
  - Cached per physical source key: `{source}|{WIGOS_ID}|{identifier}` (EARLINET: `earlinet|{station}`)

**Step 2: Solar Zenith Angle (line 106)**

`sza = solar_zenith_angle(L2.time, L2.station_latitude, L2.station_longitude)` (see §5)
- Returns [deg], shape (Nprofiles, 1)

**Step 3: Cloud-Base Height (lines 109–114)**

```matlab
if has field cloud_base_height:
    cbh = min(L2.cloud_base_height, [], 2, 'omitnan')  % lowest of 3 layers
    cbh(cbh <= 0 | cbh > 20000) = NaN                   % guard fill values
else:
    cbh = NaN(numel(L2.time), 1)
```
Units: [m] (absolute, not AGL); shape (Nprofiles, 1).

**Step 4: Calibration Correction (lines 116–118, detail in §1.4.1)**

`[correction, calibInfo] = get_calibration(cfg, ch, L2, WIGOS_ID)`

Returns: `correction` (Nprofiles, 1) multiplicative factor; `calibInfo` struct.

```matlab
beta = beta_L2 .* correction   % element-wise multiplication
```

**Step 4.1: get_calibration Details (lines 262–306)**

**Case 'none':**
- `correction = ones(Nprofiles, 1)`

**Case 'rayleigh_python':**
- `kal = load_rayleigh_python_kalman(wig, identifier, rayleighOpts...)` (see §2)
- Returns struct with `daily_dates`, `daily_C_kalman` [Rayleigh lidar constant, m sr / Mm⁻¹]
- `CL_kal = interp_calib(kal.daily_dates, kal.daily_C_kalman, L2.time)` (see §1.6)
  - Interpolates to profile times, clamping (holding first/last value beyond range)
- **Formula:**
  ```
  correction = L2.calibration_constant_0 / CL_kal   [element-wise division]
  ```
  - Interpretation: `beta_corrected = beta_L2 * (CL_L2 / CL_Rayleigh)` — converts from L2's calibration to Rayleigh

**Case 'cloud':**
- Load: `{folderCloudKalman}/calibration_{wig}_{identifier}.mat`
- Extract: `kal.daily_dates`, `kal.daily_C_kalman` [cloud lidar constant]
- `C_kal = interp_calib(kal.daily_dates, kal.daily_C_kalman, L2.time)` (clamping)
- **Formula:**
  ```
  correction = C_kal   [direct application, already in Mm⁻¹ sr⁻¹ units]
  ```

**Step 5: Optional Background Correction (lines 121–126)**

If `ch.bgCorrection == true`:
- `[bgProfile, ~] = estimate_background_profile(L2, sza, bgMinRangeAGL)` (see §4)
- `[beta_bg, ~] = apply_background_correction(beta_L2, L2, sza, bgProfile)` (see §4.2)
- `beta = beta_bg .* correction` (apply correction to background-corrected signal)

**Step 6: Water-Vapor Absorption Correction (lines 128–136)**

If `900 ≤ wavelength ≤ 920` and not manually disabled:
- `[lam0, fwhm] = wv_params_for_type(L2.instrument_type)` (see §1.5)
- Allow per-channel override: `ch.wvWavelength`, `ch.wvFwhm`
- `[beta, wvInfo] = apply_wv_correction_to_L2(L2, beta, lam0, fwhm, wvcfg)` (see §1.7)
  - **Formula:** `beta_corr = beta / T2_wv` where T²_wv is two-way water-vapor transmission [0,1]
  - Units: returns corrected backscatter [Mm⁻¹ sr⁻¹]
  - **Key rule:** Months without CAMS data → profiles set to NaN (not approximated)
  - Returns `wvInfo` struct with: `.applied` (bool), `.wavelength`, `.laser_fwhm`, `.months`, `.months_excluded`, `.median_t2`, `.n_excluded`

**Step 7: Wavelength Normalization to Target (lines 138–151)**

1. Read channel wavelength: `lambda = L2.l0_wavelength` [nm]
2. Choose model: default 'angstrom' or per-channel override
3. **If 'molaer' AND |λ − λ_target| > 1 nm:**
   - `beta = convert_wavelength_molaer(beta, L2, lambda, cfg)` (see §1.8)
4. **Else (standard Ångström):**
   ```matlab
   wl_corr = (lambda / lambda_target)^(-angstromExponent)
   beta = beta / wl_corr
   ```
   - **Formula:**
     ```
     β_target = β_λ × (λ / λ_target)^α
     ```
   - Example: λ=910 nm, λ_target=1064 nm, α=1 → `β_1064 = β_910 × (910/1064)^1 = β_910 × 0.856`
   - **Sign convention:** Exponent α > 0 means shorter wavelength → larger backscatter (Ångström scaling).

**Step 8: Two Data Streams (lines 153–156)**

**Display stream (quality flagged only):**
```matlab
beta_disp = beta
beta_disp(L2.quality_flag > 0) = NaN   % mask quality_flag != 0
```

**Screened stream (full filtering):**
```matlab
beta_scr = screen_profiles(beta, L2, sza, cfg)
```
(See §1.4.2 for screening logic)

**Step 1.4.2: screen_profiles Detail (lines 309–344)**

Applied to `beta_scr` only (used for profile overlay, scatter, and statistics).

**1. Quality flag masking:**
```matlab
beta(quality_flag > 0) = NaN
```

**2. Cloud/rain/fog exclusion (±15 min):**
```matlab
excl_margin = minutes(15)
has_cloud = any(~isnan(cloud_base_height), 2)  % any layer with CBH
vv = vertical_visibility
vv(vv == -1) = NaN                               % convert sentinel to NaN
has_excl = has_cloud | ~isnan(vv)              % cloud OR fog/precip indicator

dt = median(diff(L2.time))
if dt > 0:
    win = 2 * round(excl_margin / dt) + 1       % symmetric window
    expanded = movmax(double(has_excl), win) > 0 % dilate by ±15 min
else:
    expanded = has_excl
beta(expanded, :) = NaN
```

**3. Night filtering (optional):**
```matlab
if cfg.nightOnly:
    beta(sza < cfg.szaThreshold, :) = NaN      % keep sza >= 95°
```

**4. SNR filtering (optional):**
```matlab
if cfg.snrThreshold > 0:
    edges = min(L2.time) : minutes(60) : max(L2.time) + minutes(60)
    [~, ~, binIdx] = histcounts(L2.time, edges)
    for b = 1:max(binIdx):
        idx = find(binIdx == b)
        if numel(idx) < 3: continue
        blk = beta(idx, :)
        med = median(blk, 1, 'omitnan')
        mad = 1.4826 * median(abs(blk - med), 1, 'omitnan')  % MAD = 1.4826 * median(|x - med|)
        low = abs(med) ./ max(mad, eps) < cfg.snrThreshold
        beta(idx, low) = NaN
```
- **Interpretation:** Within each hourly bin, compute robust noise level (MAD at each altitude) and mask gates where the signal-to-noise ratio (median / MAD) is below threshold.

#### 1.5 wv_params_for_type (lines 378–397)

Measured emission wavelength and FWHM per instrument type (Payerne Qmini campaign 2026-06-02):

| Type | λ₀ [nm] | FWHM [nm] |
|------|---------|-----------|
| CL31 | 909.7 | 6.0 |
| CL51 | 910.0 | 3.4 |
| CL61 | 910.74 | 1.0 |
| CHM15K, CHM15KX, CHM8K | 1064.47 | 0.5 |
| MINI-MPL, MINIMPL, MPL | 532.0 | 0.1 |
| EARLINET | 1064.0 | 0.5 |
| Default (unrecognized) | 910.0 | 3.4 |

#### 1.6 interp_calib (§2 reference)

Interpolates daily calibration series onto profile times with clamping:
```matlab
function v = interp_calib(dd, val, t)
dd, val cleaned: remove NaN/NaT
if numel(val) == 1:
    v = repmat(val, size(t))
else:
    v = interp1(dd, val, t, 'linear')   % linear interpolation
v(t < dd(1)) = val(1)                   % clamp before first date
v(t > dd(end)) = val(end)                % clamp after last date
```

#### 1.7 apply_wv_correction_to_L2 (§1.8 reference; full details in §1.7 of this spec)

Returns: `[beta_corr, info]` where
- `beta_corr` [Mm⁻¹ sr⁻¹]: corrected backscatter (or NaN for months without CAMS)
- `info`: struct with `.applied`, `.wavelength`, `.laser_fwhm`, `.months` (cell of 'yyyymm'), `.months_excluded`, `.median_t2`, `.n_excluded`

**Key logic:**
- Only works for 900 ≤ λ ≤ 920 nm
- Requires exact-month CAMS file: `{cams_folder}/CAMS_Beta_{yyyymm}.nc`
- Missing CAMS → **set entire month to NaN** (no fallback)
- Calls `compute_wv_transmission(data, cfg)` per month (see §1.8)

#### 1.8 convert_wavelength_molaer (lines 400–436)

Two-component (molecular + aerosol) wavelength conversion for large wavelength separations.

**Inputs:**
- `beta` [Mm⁻¹ sr⁻¹]: attenuated backscatter at wavelength λ
- `L2`: struct (uses `.altitude`, `.station_altitude`)
- `lambda`: actual wavelength [nm]
- `cfg`: struct (uses `.lambda_target`, `.angstromExponent`)

**Algorithm:**
1. **Load US-1976 standard atmosphere once (cached):**
   - `[ZS, TS, PS] = get_std_atm()` [m, K, Pa]
   
2. **Create uniform AGL grid for molecular calculations:**
   - `zc = 0:30:15000` [m AGL]
   - Interpolate T, P to this grid
   
3. **Compute molecular backscatter at both wavelengths:**
   - `[~, ~, ~, bml] = get_rayleigh_v3(zc, Pc, Tc, lambda)` [m⁻¹ sr⁻¹]
   - `[~, ~, ~, bmt] = get_rayleigh_v3(zc, Pc, Tc, lambda_target)` [m⁻¹ sr⁻¹]
   - Convert: `bml *= 1e6`, `bmt *= 1e6` → [Mm⁻¹ sr⁻¹]
   
4. **Interpolate molecular profiles to native altitude grid (AGL):**
   - `z_agl = L2.altitude - L2.station_altitude` [m]
   - `bml_z = interp1(zc, bml, z_agl, 'linear', NaN)` [1 × Nalt]
   - `bmt_z = interp1(zc, bmt, z_agl, 'linear', NaN)` [1 × Nalt]
   
5. **Split aerosol component and rescale:**
   - `wl_corr = (lambda / lambda_target)^(-angstromExponent)`
   - `beta_aer = beta - bml_z` [Mm⁻¹ sr⁻¹, broadcasts over time]
   - `beta_t = bmt_z + beta_aer / wl_corr`
   
**Formula:**
```
β_target = β_mol(λ_target, z) + [β(λ, z) − β_mol(λ, z)] × (λ_target/λ)^α
         = β_mol(λ_target, z) + [β(λ, z) − β_mol(λ, z)] / wl_corr
```

#### 1.9 Temporal Synchronization (lines 175–190)

Build one regular time grid (union of all channels):

```matlab
for k = 1:nCh:
    vScr = sprintf('scr_%d', k)
    vDisp = sprintf('disp_%d', k)
    vCbh = sprintf('cbh_%d', k)
    TT = timetable(chData{k}.time, chData{k}.beta_scr, chData{k}.beta_disp, chData{k}.cbh, ...)
    TTs{k} = retime(TT, 'regular', 'median', 'TimeStep', cfg.dtBin)
end
TTsync = TTs{1}
for k = 2:nCh:
    TTsync = synchronize(TTsync, TTs{k}, 'union')
```

- Each channel: resample to regular `dtBin` grid using **median aggregation**
- Then take the **union** of all time grids across channels
- Result: `time_sync` (Nt, 1) datetime, common to all channels

#### 1.10 Common Altitude Grid & Averaging (lines 192–205)

```matlab
altGrid = build_common_grid({alt₁, alt₂, ..., alt_nCh})
```

**Algorithm (lines 347–363):**
```matlab
z_start = -inf; z_end = inf; dz_common = 0
for each alt list:
    dz = median(diff(alt))
    dz_common = max(dz_common, dz)          % coarsest native resolution
    z_start = max(z_start, min(alt))
    z_end = min(z_end, max(alt))
dz_common = max(round(dz_common), 1)
altGrid = z_start : dz_common : z_end      % uniform grid [m ASL]
```

**Averaging onto common grid (lines 366–375):**
```matlab
half_dz = median(diff(altGrid)) / 2
G = nan(nT, nZ)
for iz = 1:nZ:
    sel = abs(alt_src - altGrid(iz)) < half_dz
    if any(sel):
        G(:, iz) = mean(B(:, sel), 2, 'omitnan')
```
- For each altitude in common grid, average all source-grid gates within ±half_dz

#### 1.11 Statistics vs Reference Channel (lines 209–238)

```matlab
zMin_asl = cfg.zMin + R.station.altitude      % convert AGL to ASL
zMax_asl = cfg.zMax + R.station.altitude
zMask = (altGrid >= zMin_asl) & (altGrid <= zMax_asl)
iref = cfg.referenceChannel
ref = R.beta{iref}(:, zMask)                  % reference channel, masked

for k = 1:nCh:
    cur = R.beta{k}(:, zMask)                 % current channel
    m = isfinite(cur) & isfinite(ref)
    a = cur(m); b = ref(m)
    n = numel(a)
    
    if n > 2:
        d = a - b                              % differences
        bias = mean(d, 'omitnan')
        medbias = median(d, 'omitnan')
        rmse = sqrt(mean(d.^2, 'omitnan'))
        std = std(d, 'omitnan')
        relbias_pct = 100 * mean(d, 'omitnan') / mean(b, 'omitnan')
        cc = corrcoef(a, b)                    % Pearson correlation matrix (2x2)
        r = cc(1, 2)                           % off-diagonal element
```

**Definitions:**
- **Bias:** `mean(current − reference)`
- **Median bias:** `median(current − reference)`
- **RMSE:** `sqrt(mean((current − reference)²))`
- **Std dev:** sample standard deviation of differences
- **Relative bias:** `100 × mean(current − reference) / mean(reference)` [%]
- **Correlation r:** Pearson coefficient (linear correlation, range [−1, 1])

#### 1.12 Output Structure R

```matlab
R.channels(k)      % metadata: label, calib, source, wavelength, instrument_type, 
                   % identifier, median_correction, wv
R.altGrid          % (1 x Nz) m ASL, common altitude grid
R.time_sync        % (Nt x 1) datetime, synchronized time grid
R.beta{k}          % (Nt x Nz) Mm⁻¹ sr⁻¹, screened corrected backscatter
R.beta_disp{k}     % (Nt x Nz) Mm⁻¹ sr⁻¹, display (quality-flagged only)
R.cbh_disp{k}      % (Nt x 1) m, cloud-base height [m] on time_sync
R.station          % struct: site_location, latitude, longitude, altitude
R.stats(k)         % struct: n, bias, medbias, rmse, std, relbias_pct, r
R.cfg              % resolved configuration
```

---

### 2. RAYLEIGH CALIBRATION: `load_rayleigh_python_kalman.m`

#### 2.1 Purpose
Load corrected Python Rayleigh calibration CSV, apply Kalman filter, return daily lidar constant with uncertainty.

#### 2.2 Input / Output

**Inputs:**
- `WIGOS_ID`: E-PROFILE station WIGOS id (e.g., '0-20000-0-06610')
- `identifier`: instrument identifier (e.g., 'A', 'C')
- **Options (name-value pairs):**
  - `'CsvRoot'`: root of Python output tree (default: `'D:\E-PROFILE_calibration_rayleigh\fullcal_all'`)
  - `'CacheFolder'`: Kalman result cache dir (default: `'A:\E-PROFILE_L2_Calibration\rayleigh_python_kalman'`)
  - `'PyScript'`: path to `run_kalman_from_matlab.py` (default: `'C:\Users\hervo\OneDrive\Documents\Python\improve_alc_calib\run_kalman_from_matlab.py'`)
  - `'PyExe'`: Python executable (default: `'python'`)
  - `'StartDate'`, `'EndDate'`: datetime, filter nights within range (default: NaT = no limit)
  - `'Recompute'`: ignore cache, recompute (default: false)

**Output:** `kal` struct or `[]` (empty if no usable data)
```matlab
kal.daily_dates        % (N x 1) datetime, daily Kalman grid
kal.daily_C            % (N x 1) raw daily-median CL [lidar constant units]
kal.daily_C_std        % (N x 1) raw daily spread
kal.daily_C_kalman     % (N x 1) Kalman-smoothed CL [same units]
kal.daily_C_std_kalman % (N x 1) Kalman uncertainty
kal.meta               % struct: wigos_id, identifier, median_CL, csv_file, 
                       % n_nights, source='python_rayleigh_sign_corrected', 
                       % kalman_normalised=true, kalman_fallback
```

#### 2.3 Processing Steps

**1. Cache lookup (lines 63–68):**
- Cache file: `{CacheFolder}/rayleigh_{WIGOS_ID}_{identifier}.mat`
- If exists and `Recompute=false`, load and return cached struct

**2. Locate CSV (lines 71–77):**
- CSV path: `{CsvRoot}/{WIGOS_ID}_{identifier}/{WIGOS_ID}_{identifier}_cl.csv`
- Fail if not found

**3. Read & filter CSV (lines 79–126):**
- **Columns expected:** `date` (yyyyMMdd), `flag`, `lidar_constant`, `uncertainty` (optional)
- **Parse night date:** convert 'date' column to datetime (yyyyMMdd format)
- **Filter successful nights:** `flag == 1` (clear) OR `flag == 0.5` (partially clear)
- **Keep valid nights:** `isfinite(lidar_constant) & (lidar_constant > 0) & ~isnat(nightDate)`
- **Apply date range:** if `StartDate` or `EndDate` given, keep `nightDate >= StartDate` and `nightDate <= EndDate`

**4. Robust outlier rejection (lines 108–126):**
```matlab
if numel(clVal) >= 5:
    logCL = log(clVal)
    medlog = median(logCL)
    madlog = 1.4826 * median(abs(logCL - medlog))
    if madlog > 0:
        keep = abs(logCL - medlog) <= 4 * madlog   % reject > 4 MAD in log space
        if any(~keep):
            clVal = clVal(keep); uncVal = uncVal(keep); etc.
```
- **Interpretation:** Remove nights with physically impossible lidar constants (e.g., ~1e17 on poor nights)

**5. Min 3 nights required (lines 128–132):**
```matlab
if numel(clVal) < 3:
    return []
```

**6. Normalise before Kalman (lines 139–155):**
```matlab
medCL = median(clVal, 'omitnan')
cNorm = clVal / medCL              % dimensionless O(1) coefficient
cNormStd = uncVal / medCL
inTime = nightDate
inC = cNorm
inCstd = cNormStd
% append trailing NaN row to extend grid to EndDate (if given)
if ~isnat(EndDate):
    inTime = [inTime; EndDate + days(1)]
    inC = [inC; NaN]
    inCstd = [inCstd; NaN]
```
- **Reason:** Kalman filter tuned for O(1) coefficients; rescale back after filtering

**7. Run Kalman filter via Python (lines 157–176):**
- Write temp CSV: `{tempname}_rayleigh_kal_in.csv` with columns `time` (yyyy-mm-dd), `C`, `C_std`
- Call: `python {PyScript} {tmpIn} {tmpOut}`
- Read result: `{tmpOut}` with columns `time`, `C_daily`, `C_daily_std`, `C_kalman`, `C_kalman_std`
- Fail if script not found or execution error

**8. Rescale Kalman output (lines 181–187):**
```matlab
kal.daily_dates = datetime(Tout.time, 'InputFormat', 'yyyy-MM-dd')
kal.daily_C = Tout.C_daily * medCL
kal.daily_C_std = Tout.C_daily_std * medCL
kal.daily_C_kalman = Tout.C_kalman * medCL       % back to CL units
kal.daily_C_std_kalman = Tout.C_kalman_std * medCL
```

**9. Fallback (lines 194–203):**
```matlab
if nnz(isfinite(kal.daily_C_kalman)) < 2:    % Kalman output all NaN
    kal.daily_C_kalman = repmat(medCL, size(kal.daily_dates))
    kal.daily_C_std_kalman = repmat(std(clVal, 'omitnan'), size(kal.daily_dates))
    kal.meta.kalman_fallback = 'constant_median'
```

**10. Cache result (lines 205–211):**
- Save to `{CacheFolder}/rayleigh_{WIGOS_ID}_{identifier}.mat`

---

### 3. EARLINET DATA: `read_earlinet_att_backscatter.m`

#### 3.1 Purpose
Read EARLINET Level-2 1064 nm backscatter profiles, compute attenuated backscatter using Rayleigh formula, return structure compatible with E-PROFILE L2.

#### 3.2 Station Metadata

| Code | Folder | Overlap Min [m AGL] | Name |
|------|--------|---------------------|------|
| 'sir' | A:\EARLINET\sir\ | 2000 | Palaiseau (SIRTA) |
| 'lei' | A:\EARLINET\lei\ | 800 | Leipzig |
| 'cbw' | A:\EARLINET\cbw\ | 1000 | Cabauw |
| 'ino' | A:\EARLINET\ino\ | 1100 | Magurele |
| 'ari' | A:\EARLINET\ari\ | 800 | Leipzig (ARI) |

#### 3.3 File Discovery & Deduplication (lines 70–111)

**Pattern:** `{folder}/*{station}*b1064*20*.nc`

**Filename:** `EARLINET_AerRemSen_{stn}_Lev02_b1064_{yyyyMMddHHmm}_{yyyyMMddHHmm}_{vXX}_{qcYY}.nc`

**Deduplication logic:**
```matlab
for each file:
    split on '_'; extract parts[6] (start), parts[7] (end), parts[8] (version)
    key = parts[6] '_' parts[7]                         % "yyyyMMddHHmm_yyyyMMddHHmm"
    version = sscanf(parts[8], 'v%d')                  % e.g., v02 -> 2
[unique_keys, ~, key_idx] = unique(keys)
for each unique key:
    keep only highest version
```

#### 3.4 Molecular Rayleigh Calculation (lines 121–125)

```matlab
[z_atm, T_std, P_std] = get_std_atm()      % US-1976 standard atmosphere
range_ref = (0:15:15000)'                  % uniform AGL grid [m], column vector
P_interp = interp1(z_atm, P_std, range_ref + station_altitude_m)  % interpolate
T_interp = interp1(z_atm, T_std, range_ref + station_altitude_m)
[alpha_mol, beta_mol] = get_rayleigh_v3(range_ref, P_interp, T_interp, 1064)
% alpha_mol, beta_mol: [m⁻¹], [m⁻¹ sr⁻¹] (molecular extinction & backscatter)
```

#### 3.5 Per-File Processing (lines 134–181)

For each file (parallelized `parfor`):

**1. Read data:**
```matlab
altitude_raw_asl = ncread(filepath, 'altitude')'       % [1 x N_alt] m ASL
time_raw = ncread(filepath, 'time')'                   % [1 x N_time] seconds since 1970-01-01
backscatter_raw = ncread(filepath, 'backscatter')'     % [1 x N_alt] m⁻¹ sr⁻¹ (L2 particle)
```

**2. Filter by time (line 143):**
```matlab
file_time = datetime(filename_parts{6}, 'InputFormat', 'yyyyMMddHHmm')
if file_time < startDate || file_time > endDate + days(1):
    continue
```

**3. Interpolate to uniform range grid (line 155):**
```matlab
valid_mask = ~isnan(altitude_raw_asl)
range_raw = altitude_raw_asl(valid_mask) - station_altitude_m  % [m AGL]
backscatter_interp = interp1(range_raw, backscatter_raw(valid_mask), range_ref')
% interp1 returns NaN outside interpolation domain
```

**4. Overlap correction (lines 158–161):**
```matlab
idx_ov = find(range_ref' < overlap_min)
if ~empty(idx_ov):
    backscatter_interp(idx_ov) = backscatter_interp(idx_ov(end))  % fill with last valid
```

**5. Compute attenuated backscatter (lines 163–170):**
```matlab
LR_1064 = 50                                % assumed lidar ratio [sr]
extinction = backscatter_interp * LR_1064  % aerosol extinction [m⁻¹]
extinction(isnan(extinction)) = 0           % treat NaN as zero
sum_ext = cumtrapz(range_ref', extinction + alpha_mol', 2)  % cumulative optical depth [unitless]
trans = exp(-sum_ext)                       % transmission [0,1]
att_bsc = (backscatter_interp + beta_mol') .* trans .* trans  % [m⁻¹ sr⁻¹]
% att_bsc = (particle_backscatter + rayleigh) × T²
```

**6. Store result:**
```matlab
att_backscatter(f, :) = att_bsc
time_datenum(f) = time_raw / 86400 + datenum(1970, 1, 1)  % convert seconds to MATLAB datenum
valid(f) = true
```

#### 3.6 Post-Processing (lines 183–226)

```matlab
keep = valid                                   % filter to valid files
att_backscatter = att_backscatter(keep, :)
time_datenum = time_datenum(keep)

if isempty(time_datenum):
    return []                                  % no valid data

time_dt = datetime(time_datenum, 'ConvertFrom', 'datenum')
[time_dt, sortIdx] = sort(time_dt)            % sort by time
att_backscatter = att_backscatter(sortIdx, :)
[time_dt, uniqIdx] = unique(time_dt, 'stable')  % remove duplicates
att_backscatter = att_backscatter(uniqIdx, :)

att_backscatter = att_backscatter * 1e6      % m⁻¹ sr⁻¹ -> Mm⁻¹ sr⁻¹
```

#### 3.7 Output Structure

```matlab
L2.time                     % (Np x 1) datetime
L2.altitude                 % (Nalt x 1) m ASL = range_ref + station_altitude_m
L2.attenuated_backscatter_0 % (Np x Nalt) Mm⁻¹ sr⁻¹
L2.calibration_constant_0   % (Np x 1) all ones (already calibrated)
L2.quality_flag             % (Np x Nalt) int32, all zeros
L2.cloud_base_height        % (Np x 3) NaN (not available)
L2.vertical_visibility      % (Np x 1) NaN (not available)
L2.l0_wavelength            % scalar: 1064 [nm]
L2.station_altitude         % station altitude [m ASL]
L2.station_latitude         % latitude [deg]
L2.station_longitude        % longitude [deg]
L2.instrument_type          % 'EARLINET'
L2.site_location            % station name string
```

---

### 4. BACKGROUND CORRECTION

#### 4.1 estimate_background_profile.m

**Purpose:** Estimate instrument background noise from nocturnal clear-sky observations.

**Inputs:**
- `L2`: struct with `.time`, `.altitude` [m ASL], `.attenuated_backscatter_0`, `.cloud_base_height`, `.station_altitude`
- `sza`: solar zenith angle [deg]
- `minRangeAGL`: [m] minimum range for estimation (default: 2400)

**Parameters (hard-coded):**
- `szaNight = 108`: SZA threshold for deep night classification [deg]
- `maxCloudFrac = 0.10`: max fraction of cloudy profiles per hourly bin
- `smoothWindow = 21`: rolling mean window [gates]

**Algorithm:**

**1. Select nighttime profiles (sza > 108°)**

**2. Restrict to range > minRangeAGL**

**3. Group into hourly bins, filter by cloud fraction**
- If cloud_base_height non-NaN for any layer → has cloud
- Keep hourly bins with cloudFrac < 0.10

**4. Convert to pseudo-signal:**
```matlab
PSGN = beta ./ range_agl²   % range_agl in [m]
```

**5. Compute global statistics:**
```matlab
medProfile = median(PSGN, 1, 'omitnan')
pct25 = prctile(PSGN, 25, 1)
pct75 = prctile(PSGN, 75, 1)
bg_global = -medProfile      % negate (background is negative noise)
bg_global_smooth = movmean(bg_global, smoothWindow, 'omitnan')
```

**6. Temporal analysis (if ≥3 months and ≥100 valid profiles):**
- Call `analyze_background_temporal(PSGN, validNightTimes, range_high)` → returns segments and monthly profiles per segment
- Build time-varying `bgProfile` (Nprofiles × Nalt)
- Else: use static `bg_global_smooth` for all times

**7. Extend to full altitude range:**
```matlab
for each altitude:
    if range > minRangeAGL: use smoothed profile
    if range <= minRangeAGL: constant = profile(first gate)
    if range > highest gate: constant = profile(last gate)
```

**Output:**
- `bgProfile`: (Nprofiles × Nalt) OR (1 × Nalt) static
- `diagInfo`: struct with diagnostic data

#### 4.2 apply_background_correction.m

**Purpose:** Apply background correction in pseudo-signal space, optionally with diurnal modulation.

**Inputs:**
- `beta`: (Nprofiles × Nalt) attenuated backscatter [Mm⁻¹ sr⁻¹]
- `L2`: struct with `.time`, `.altitude`, `.station_altitude`
- `sza`: solar zenith angle [deg]
- `bgProfile`: (1 × Nalt) static OR (Nprofiles × Nalt) time-varying

**Parameters (hard-coded):**
- `topGateDepth_m = 300`: depth of top gates for noise estimation [m]
- `cirrusSmoothWin = 7`: smooth window for cirrus detection
- `smoothMinutes = 25`: temporal smoothing for noise [min]
- `szaNight = 108`: night threshold [deg]

**Algorithm:**

**1. Convert to pseudo-signal:**
```matlab
R² = range_agl²
PSGN = beta / R²     % element-wise
```

**2. Expand bgProfile if static:**
```matlab
if size(bgProfile, 1) == 1:
    ADD = repmat(bgProfile, nT, 1)
else:
    ADD = bgProfile
```

**3. Check if dynamic modulation needed:**
```matlab
meanProfile = mean(ADD, 1, 'omitnan')
biasSignificant = abs(meanProfile(1)) > 3 * std(meanProfile, 'omitnan')
```

**4. If biasSignificant, apply dynamic modulation:**

**a) Extract top gates:**
```matlab
nTop = round(topGateDepth_m / dR)
if even: nTop += 1
nTop = min(nTop, nAlt)
topCols = (nAlt - nTop + 1) : nAlt
TOP = PSGN(:, topCols)   % (nT × nTop)
```

**b) 2D smooth + compute relative variance:**
```matlab
TOPm = smooth2d_nanmean(TOP, cirrusSmoothWin, cirrusSmoothWin)  % smooth
TOPs = smooth2d_nanstd(TOP, cirrusSmoothWin, cirrusSmoothWin)
RV = (TOPs ./ TOPm)²     % relative variance
```

**c) Mask cirrus (high variance → cloud signal):**
```matlab
cirrus = isnan(RV) | (RV < 1)      % RV < 1 = noise-dominated, keep
TOP(cirrus) = NaN
badTime = sum(cirrus, 2) > 10       % timesteps with too many masked gates
```

**d) Temporal smoothing:**
```matlab
dt_sec = median(seconds(diff(L2.time)))
wt = round((60 / dt_sec) * smoothMinutes) + 1
if even: wt += 1
topSmooth = smooth2d_nanmean(TOP, wt, nTop)
centralCol = floor((nTop - 1) / 2) + 1
noiseMean = topSmooth(:, centralCol)
noiseMean(badTime) = NaN
noiseMean = fillmissing(noiseMean, 'linear')
noiseMean = fillmissing(noiseMean, 'nearest')
Fback = -noiseMean
```

**e) Nocturnal baseline:**
```matlab
nightMask = sza(:) > szaNight
if any(nightMask & ~isnan(Fback)):
    base = mean(Fback(nightMask), 'omitnan')
else:
    base = mean(Fback, 'omitnan')
ADD = ADD + (Fback - base)      % modulate climatology with dynamic noise
```

**5. Apply correction:**
```matlab
BSGN = PSGN + ADD              % add background in pseudo-signal space
beta_corrected = BSGN .* R²    % convert back to backscatter units
```

**Output:**
- `beta_corrected`: (Nprofiles × Nalt) [Mm⁻¹ sr⁻¹]
- `diagInfo`: struct with diagnostic info

---

### 5. WATER-VAPOR CORRECTION

#### 5.1 apply_wv_correction_to_L2.m

**Purpose:** Remove 910 nm water-vapor absorption to match 1064 nm instruments.

**Inputs:**
- `L2`: struct (uses `.time`, `.altitude`, `.station_altitude/latitude/longitude`)
- `beta`: (nTime × nAlt) attenuated backscatter [Mm⁻¹ sr⁻¹]
- `wavelength`: laser wavelength [nm]
- `laser_fwhm`: laser spectral width [nm]
- `wvcfg`: struct with `.abs_cs_lookup_table` (path), `.cams_folder` (path)

**Algorithm (lines 36–99):**

**1. Check wavelength range:**
```matlab
if wavelength < 900 or wavelength > 920:
    return beta (unchanged), info.applied = false
```

**2. Parse time into (YYYY)MM groups:**
```matlab
range = L2.altitude - L2.station_altitude    % [m AGL]
t = L2.time
ym_all = cellstr(datestr(t, 'yyyymm'))        % 'yyyymm' per profile
ym_unique = unique(ym_all, 'stable')
```

**3. For each month:**
```matlab
for each ym:
    idx = strcmp(ym_all, ym)  % indices of profiles in this month
    
    % MUST have exact-month CAMS file
    camsFile = {cams_folder}/CAMS_Beta_{ym}.nc
    if ~isfile(camsFile):
        beta_corr(idx, :) = NaN         % EXCLUDE entire month
        n_excluded += sum(idx)
        months_excluded{end+1} = ym
        continue
    
    % Build month's data struct
    data.range = range
    data.time = t(idx)
    data.station_altitude/latitude/longitude = ...
    data.beta = beta(idx, :)'                % [nAlt × nMonthTimes]
    
    % Compute transmission
    trans2 = compute_wv_transmission(data, cfg)   % [nAlt × nMonthTimes]
    beta_corr(idx, :) = beta(idx, :) / trans2'   % divide by T²
    t2_accum = [t2_accum; trans2(:)]
    months{end+1} = ym
```

**4. Summary:**
```matlab
info.applied = ~isempty(months)
info.median_t2 = median(t2_accum, 'omitnan')
info.months = months
info.months_excluded = months_excluded
info.n_excluded = n_excluded
```

**Output:**
- `beta_corr`: (nTime × nAlt) [Mm⁻¹ sr⁻¹] with WV correction applied (or NaN for excluded months)
- `info`: struct with applied flag, wavelength, FWHM, months processed, median T², exclusion count

#### 5.2 compute_wv_transmission.m

**Purpose:** Compute two-way water-vapor transmission from CAMS T/RH and HITRAN absorption cross-sections.

**Inputs:**
- `data`: struct with `.range` [m AGL], `.time` [datetime], `.station_altitude/latitude/longitude`, `.beta` (nAlt × nTime)
- `config`: struct with `.abs_cs_lookup_table` (path to LUT NetCDF), `.cams_folder` (path), `.wavelength` [nm], `.laser_fwhm` [nm], `.date_str` ('yyyymm')

**Algorithm (lines 48–143):**

**1. Load absorption cross-section LUT:**
```matlab
lut = ncread(abs_cs_lookup_table, ...)
abs_cs_height = lut.height_in_km * 1000       % km -> m
abs_cs_wl = lut.lambda                        % wavelength grid [nm]
abs_cs_all_values = lut.abscs_ave             % [n_wl × n_height]
```

**2. Load CAMS T/RH profiles:**
```matlab
[time_CAMS, CAMS_z, CAMS_data] = get_Beta_CAMS_oper_monthly(
    date_str, station_lat, station_lon, info_cams, [], cams_folder, 'CAMS_Beta')
% CAMS_data.T [K], CAMS_data.RH [0-1], shape (n_height × n_cams_time)
% CAMS_z [m ASL], shape (n_height × n_cams_time)
```

**3. Precompute laser spectrum (Gaussian):**
```matlab
sigma = laser_fwhm / (2 * sqrt(2 * ln(2)))
gauss = normpdf(abs_cs_wl(:), wavelength, sigma)
index_sum = (abs_cs_wl(:) > wavelength - 3*sigma) & (abs_cs_wl(:) < wavelength + 3*sigma)
sum_gauss = sum(gauss)
```

**4. Map ceilometer range to LUT heights:**
```matlab
[~, height_indices] = min(abs(range_col' - abs_cs_height(:)), [], 1)
abs_cs = abs_cs_all_values(:, height_indices)  % [n_wl × n_range]
```

**5. Compute water-vapor density for each CAMS time step:**
```matlab
nw_all = get_water_vapor_number_concentration_from_RH(CAMS_data.T, CAMS_data.RH, CAMS_z)
% [n_height × n_cams] number concentration [m⁻³]

wv_density_all = zeros(n_range, n_cams)
for i = 1:n_cams:
    z_agl = CAMS_z(:,i) - station_altitude
    wv_density = interp1(z_agl, nw_all(:,i), range_col)
    % fill NaN below lowest model level with nearest valid value
    nan_mask = isnan(wv_density)
    if any(nan_mask):
        first_valid = find(~nan_mask, 1, 'first')
        if ~empty(first_valid):
            wv_density(nan_mask & range_col < max(z_agl)) = wv_density(first_valid)
    wv_density(isnan(wv_density)) = 0
    wv_density_all(:, i) = wv_density
```

**6. Compute transmission for each CAMS time step (vectorized):**
```matlab
for i = 1:n_cams:
    ext = abs_cs .* (wv_density_all(:, i)' / 1e4)  % [n_wl × n_range]
    sum_ext = zeros(size(ext))
    sum_ext(index_sum, :) = cumtrapz(range_col, ext(index_sum, :), 2)
    trans = exp(-sum_ext)                          % [n_wl × n_range]
    trans2_cams(:, i) = sum(trans.^2 .* gauss, 1)' / sum_gauss
    % Gaussian-weighted average of trans² over wavelength
```

**7. Interpolate from CAMS to ceilometer time grid:**
```matlab
if n_cams == 1:
    trans2 = repmat(trans2_cams, 1, n_profiles)
else:
    ceilo_dn = datenum(data.time)
    cams_dn = datenum(time_CAMS)
    trans2 = interp1(cams_dn, trans2_cams', ceilo_dn, 'nearest', 'extrap')'
```

**8. Safety check:**
```matlab
trans2(trans2 <= 0 | isnan(trans2)) = 1   % clamp: no correction if invalid
```

**Output:**
- `trans2`: (n_range × n_profiles) two-way transmission [0,1]

---

### 6. RAYLEIGH MOLECULAR PROPERTIES: `get_rayleigh_v3.m`

#### 6.1 Purpose
Compute molecular (Rayleigh) backscatter, extinction, and attenuated backscatter from P, T, wavelength.

#### 6.2 Inputs

| Param | Type | Units | Notes |
|-------|------|-------|-------|
| `range` | (N,1) | m | Range from instrument (AGL-like, NOT altitude) |
| `P` | (N,1) | Pa | Pressure; must all be > 15000 Pa (check) |
| `T` | (N,1) | K | Temperature; must all be > 50 K (check) |
| `lambda_rec` | scalar | nm | Reception wavelength (usually same as emission) |
| `lambda_em` | scalar (optional) | nm | Emission wavelength (for Raman); if omitted, λ_rec = λ_em |

#### 6.3 Outputs

| Var | Shape | Units | Definition |
|-----|-------|-------|------------|
| `alpha_mol` | (N,1) | m⁻¹ | Molecular extinction coefficient |
| `beta_mol` | (N,1) | m⁻¹ sr⁻¹ | Molecular (Rayleigh) backscatter |
| `lidar_signal` | (N,1) | m⁻¹ sr⁻¹ m² | Simulated range-corrected signal = β_mol × T_em × T_rec / r² |
| `beta_att` | (N,1) | m⁻¹ sr⁻¹ | Molecular attenuated backscatter = β_mol × T_em × T_rec |
| `density` | (N,1) | m⁻³ | Molecular number density |
| `trans` | (N,1) | — | Two-way transmission (receive) at λ_rec |
| `trans_em` | (N,1) | — | Two-way transmission (emit) at λ_em |

#### 6.4 Constants (lines 61–68)

```matlab
T0 = 273.15 + 15 = 288.15 K          % reference temperature
P0 = 101325 Pa                       % reference pressure
rho = 0.0301                         % depolarisation ratio
N = 2.547e25 m⁻³                    % reference molecular density
lambda_rec *= 1e-9                   % nm -> m
lambda_em *= 1e-9
```

#### 6.5 Refractive Index (lines 69–70)

**Cauchy formula (Ciddor, 2007):**
```matlab
m = ( 5791817 / (238.0185 - 1/(lambda_rec*1e6)²)
    + 167909 / (57.362 - 1/(lambda_rec*1e6)²) ) * 1e-8 + 1
```
- Input: λ_rec in m; computes using λ in μm

#### 6.6 Molecular Backscatter (lines 72–98)

**Number density:**
```matlab
density = N_A * P / R / T   where N_A = 6.022005e23 mol⁻¹, R = 8.314418 J mol⁻¹ K⁻¹
```

**Molecular backscatter (Bucholtz 1995):**
```matlab
for i = 1:length(range):
    beta_mol_tot(i) = 24π³ × (m² − 1)² × (6 + 3ρ) / [λ_rec⁴ × N² × (m² + 2)² × (6 − 7ρ)] × N × T₀ × P(i) / P₀ / T(i)
    beta_mol_tot_em(i) = 24π³ × (m² − 1)² × (6 + 3ρ) / [λ_em⁴ × N² × (m² + 2)² × (6 − 7ρ)] × N × T₀ × P(i) / P₀ / T(i)

gamma = rho / (2 - rho)
Pray = 3/4 / (1 + 2γ) × (1 + 3γ + (1 − γ) cos²(π))   % cos(π) = -1
beta_mol = beta_mol_tot / (4π) × Pray                 % [m⁻¹ sr⁻¹]
alpha_mol = beta_mol × 8π / 3                         % [m⁻¹]
beta_mol_em = beta_mol_tot_em / (4π) × Pray
alpha_mol_em = beta_mol_em × 8π / 3
```

#### 6.7 Transmission (lines 103–112)

```matlab
step_z = range(2) - range(1)          % uniform step size [m]
for i = 1:length(range):
    trans_em(i) = exp(-sum(alpha_mol_em(1:i) × step_z))       % cumulative extinction (emit)
    trans(i) = exp(-sum(alpha_mol(1:i) × step_z))             % cumulative extinction (receive)
    lidar_signal(i) = beta_mol(i) × trans_em(i) × trans(i) / range(i)²
    beta_att(i) = beta_mol(i) × trans_em(i) × trans(i)        % attenuated backscatter
```

---

### 7. STANDARD ATMOSPHERE: `get_std_atm.m`

**Purpose:** Load US-1976 standard atmosphere.

**File:** `standard_atmosphere_US_1976_50km.csv` (in working directory)

**Output:**
```matlab
z   % (n,1) [m] altitude MSL
T   % (n,1) [K] temperature
P   % (n,1) [Pa] pressure
```

**Source:** http://www.digitaldutch.com/atmoscalc/

---

### 8. SOLAR ZENITH ANGLE: `solar_zenith_angle.m`

#### 8.1 Purpose
Compute solar zenith angle for datetime array (Spencer 1971 / NOAA method).

#### 8.2 Inputs
- `t`: datetime array (any shape)
- `lat`: station latitude [deg, +N]
- `lon`: station longitude [deg, +E]

#### 8.3 Output
- `sza`: [deg], same shape as `t`; range [0, 180]

#### 8.4 Algorithm (lines 18–45)

**1. Fractional day of year:**
```matlab
doy = day(t, 'dayofyear') + (hour(t) + minute(t)/60 + second(t)/3600) / 24
```

**2. Fractional year [rad]:**
```matlab
yr = year(t)
nDays = 365 + (1 if leap year else 0)
gamma = 2π × (doy − 1) / nDays
```

**3. Equation of time [min] and solar declination [rad] (Spencer 1971):**
```matlab
eqtime = 229.18 × ( 0.000075
                  + 0.001868 cos(γ) − 0.032077 sin(γ)
                  − 0.014615 cos(2γ) − 0.04089 sin(2γ) )
decl = 0.006918
       − 0.399912 cos(γ) + 0.070257 sin(γ)
       − 0.006758 cos(2γ) + 0.000907 sin(2γ)
       − 0.002697 cos(3γ) + 0.00148 sin(3γ)
```

**4. True solar time:**
```matlab
utc_min = hour(t) × 60 + minute(t) + second(t) / 60
tst = utc_min + eqtime + 4 × lon    % true solar time [min]
ha = deg2rad(tst / 4 − 180)         % hour angle [rad]
```

**5. Solar zenith angle:**
```matlab
lat_r = deg2rad(lat)
cos(sza) = sin(lat_r) sin(decl) + cos(lat_r) cos(decl) cos(ha)
cos(sza) = clamp(cos(sza), -1, 1)    % for acos stability
sza = rad2deg(acos(cos(sza)))        % [deg]
```

---

### 9. READ E-PROFILE L2: `read_L2_monthly.m`

#### 9.1 Purpose
Read monthly E-PROFILE L2 NetCDF files and concatenate.

#### 9.2 File Path Pattern
```
{folderRoot}/{WIGOS_ID}/{yyyy}/L2_{WIGOS_ID}_{identifier}{yyyymm}.nc
```

#### 9.3 WIGOS ID Conversion
```matlab
if length(WIGOS_ID) == 5 & all(isstrprop(WIGOS_ID, 'digit')):
    WIGOS_ID = ['0-20000-0-' WIGOS_ID]     % convert WMO to WIGOS
```

#### 9.4 NetCDF Fields Read

| Field | Shape | Dtype | Notes |
|-------|-------|-------|-------|
| `time` | (Np,) | double | days since or seconds since epoch (parsed from units) |
| `altitude` | (Nz,) | double | m ASL |
| `attenuated_backscatter_0` | (Np, Nz) | (after transpose check) | Mm⁻¹ sr⁻¹, can be (Nz, Np) or (Np, Nz) |
| `calibration_constant_0` | (Np,) | single | varies by calib method |
| `cloud_base_height` | (Np, Nlayers) or (Np,) | single | m AGL, NaN if none; up to 3 layers |
| `vertical_visibility` | (Np,) | single | m, NaN if no fog/precip |
| `quality_flag` | (Np, Nz) or (Np,) | int32 | 0 = good; > 0 = reject |
| `l0_wavelength` | scalar or (Np,) | double | nm |
| `station_altitude` | scalar or (Np,) | double | m ASL |
| `station_latitude` | scalar or (Np,) | double | deg |
| `station_longitude` | scalar or (Np,) | double | deg |

#### 9.5 Processing

**1. Month range:** All months from startDate to endDate inclusive

**2. Transpose check (lines 164–167):**
```matlab
if size(attBsc, 1) == nAlt && size(attBsc, 2) ~= nAlt:
    attBsc = attBsc'    % ensure (Np × Nz)
```

**3. Time conversion:**
```matlab
refDateStr = extractAfter(units, 'days since ') or 'seconds since '
refDate = datetime(refDateStr(1:10), 'InputFormat', 'yyyy-MM-dd')
timeDT = refDate + days(timeRaw) or seconds(timeRaw)
```

**4. Fill value handling (lines 260–264):**
```matlab
fillVal = single(9.96921e+36)   % NetCDF default fill
attBsc(attBsc > fillVal * 0.9) = NaN
calConst(calConst > fillVal * 0.9) = NaN
cbh(cbh > fillVal * 0.9) = NaN
vv(vv > fillVal * 0.9) = NaN
```

**5. Date range trim, sort, remove duplicates**

#### 9.6 Output Structure
Same as L2 in §1.4 (per-channel processing), but from monthly files.

---

### 10. FIGURE: `paper_val_figure_payerne.m`

#### 10.1 Layout (3 × 3 tiles)

```
┌─────────────┬──────────┬──────────┐
│             │ scatter  │ histogram│  (tiles 1-3)
│   profile   ├──────────┼──────────┤
│   median    │ pcolor 1 │ pcolor 2 │  (tiles 4-6)
│   ±IQR      ├──────────┼──────────┤
│   (3 rows)  │ pcolor 3 │ pcolor 4 │  (tiles 7-9)
└─────────────┴──────────┴──────────┘
```

#### 10.2 Tile 1 (left, 3 rows): Median Profile ± IQR

**Algorithm (lines 34–75):**

For each channel:
```matlab
B = R.beta{k}(:, zMask)           % (Nt × Nz) screened data
med = median(B, 1, 'omitnan')     % median per altitude
q1 = prctile(B, 25, 1)            % 25th percentile (IQR low)
q3 = prctile(B, 75, 1)            % 75th percentile (IQR high)
nz = sum(isfinite(B), 1)          % count of valid profiles per alt
good = isfinite(med) & (med > 0) & (nz >= max(10, minProfFrac*Nt))
firstBad = find(~good & z > zMin, 1, 'first')
keep = [true ... true, false ... false]  % truncate above firstBad
```
- **Plot:**
  - Solid line: median profile (kept region, positive only)
  - Dotted lines: 25th / 75th percentile (IQR, kept region, positive only)

#### 10.3 Tile 2 (top middle): Scatter vs Reference

**Algorithm (lines 77–104):**
```matlab
for k = 1:nCh:
    if k == iref: continue
    cur = R.beta{k}(:, band)
    m = isfinite(cur) & isfinite(ref) & cur > 0 & ref > 0
    a = cur(m); b = ref(m)
    if numel(a) > 6000: subsample to 6000
    scatter(a, b, label=sprintf('%s (r=%.2f, %+.0f%%)', ...
        R.channels(k).label, R.stats(k).r, R.stats(k).relbias_pct))
plot 1:1 line (black dashed)
xlim, ylim = [max(min(allv), 1e-2), prctile(allv, 99.8)]
```
- **Axes:** log-log
- **Data:** all channels (except reference) vs reference, altitudes zMin to zMax

#### 10.4 Tile 3 (top right): Histogram of Differences

**Algorithm (lines 106–131):**
```matlab
for k = 1:nCh:
    if k == iref: continue
    d = R.beta{k}(:, band) - ref
    d = d(isfinite(d))
    histogram(d, edges, Normalization='pdf', DisplayStyle='stairs', ...
        label=sprintf('%s (med %+.2f)', R.channels(k).label, median(d)))
xline(0, 'k--')
xlim = [-xmax, xmax] where xmax = prctile(abs(dall), 99)
```
- **Units:** [Mm⁻¹ sr⁻¹] (absolute difference)
- **Overlay:** vertical line at x = 0 (perfect agreement)

#### 10.5 Tiles 5, 6, 8, 9 (2 × 2 pcolor block): Time-Altitude Contours

**Algorithm (lines 133–158):**
```matlab
for k = 1:min(nCh, 4):
    B = R.beta_disp{k}(:, zMask)'         % (Nz × Nt)
    B(B < 1e-3) = 1e-3                    % floor for log
    pcolor(time_sync, z_agl(zMask), log10(abs(B)))
    shading('flat')
    clim(opts.clim)                        % default [-2 1]
    if CBH available:
        plot(time_sync, cbh_disp{k}, '.', color='k', markersize=3)
```
- **Colormap:** log₁₀(β) with limits [clim(1), clim(2)]
- **Overlay:** cloud-base height as black dots

---

### 11. EARLINET VALIDATION: `paper_val_earlinet.m`

#### 11.1 Purpose
Validate CHM15k (Rayleigh-calibrated) against colocated EARLINET using profile-by-profile matching (not continuous sync).

#### 11.2 Inputs

```matlab
site                         % EARLINET station code
chmLabel                     % display label for CHM channel
chmWigos, chmId              % E-PROFILE WIGOS & id of CHM
startDate, endDate           % date range
opts.zMin, opts.zMax         % [m AGL] comparison band (default: 500, 5000)
opts.matchMinutes            % [min] matching window (default: 30)
opts.folderL2                % [default: 'A:\E-PROFILE_L2_monthly\']
```

#### 11.3 Processing (lines 34–79)

**1. Read EARLINET reference:**
```matlab
E = read_earlinet_att_backscatter(site, startDate, endDate)
zE = E.altitude(:)'                         % m ASL
betaE_all = double(E.attenuated_backscatter_0)  % (Np × Nz)
```

**2. Load CHM Rayleigh calibration:**
```matlab
kal = load_rayleigh_python_kalman(chmWigos, chmId)
```

**3. Loop over months with EARLINET profiles:**
```matlab
months = unique(dateshift(E.time, 'start', 'month'))
for mi = 1:numel(months):
    m0 = months(mi)
    m1 = m0 + calmonths(1) - seconds(1)
    sel = E.time >= m0 & E.time <= m1
    if ~any(sel): continue
    
    L2 = read_L2_monthly(chmWigos, chmId, m0, m1, folderL2)
    if isempty(L2): continue
    
    % Apply Rayleigh calibration to CHM
    CL = interp_calib(kal.daily_dates, kal.daily_C_kalman, L2.time)
    betaC = double(L2.attenuated_backscatter_0) .* (double(L2.calibration_constant_0) ./ CL)
    betaC(L2.quality_flag > 0) = NaN
    
    % Match profiles
    for j = 1:numel(etimes):
        etime = E.time(sel)(j)
        eprof = betaE_all(sel(j), :)
        
        w = abs(L2.time - etime) <= minutes(matchMinutes)
        if ~any(w): continue
        
        cprof = median(betaC(w, :), 1, 'omitnan')            % CHM on CHM grid
        cprofE = interp1(double(L2.altitude), cprof, zE, 'linear', NaN)  % to EARLINET grid
        
        mE = [mE; eprof]
        mC = [mC; cprofE]
        mT = [mT; etime]
```

#### 11.4 Statistics (lines 86–104)

```matlab
band = altAGL >= zMin & altAGL <= zMax
a = betaC(:, band)
b = betaE(:, band)
ok = isfinite(a) & isfinite(b) & a > 0 & b > 0
av = a(ok); bv = b(ok)
d = av - bv

stats.n = numel(av)
stats.bias = mean(d)
stats.medbias = median(d)
stats.rmse = sqrt(mean(d.^2))
stats.relbias_pct = 100 * mean(d) / mean(bv)
cc = corrcoef(av, bv)
stats.r = cc(1, 2)
```

#### 11.5 Output Structure R

```matlab
R.altitude          % (1 × Nz) m ASL (EARLINET grid)
R.altAGL            % (1 × Nz) m AGL
R.time              % (Np × 1) matched profile times
R.betaE             % (Np × Nz) EARLINET attenuated backscatter
R.betaC             % (Np × Nz) CHM on EARLINET grid
R.station           % struct: site_location, altitude, latitude, longitude
R.labels            % cell: {'EARLINET ...', chmLabel}
R.stats             % struct: n, bias, medbias, rmse, relbias_pct, r
R.cfg               % struct: site, chmWigos, chmId, startDate, endDate, zMin, zMax, n_earlinet, n_matched
```

---

### 12. EARLINET FIGURE: `paper_val_earlinet_figure.m`

#### 12.1 Layout (2 × 2 tiles)

```
┌──────────────┬──────────────┐
│  (a) med±IQR │ (b) scatter  │
├──────────────┼──────────────┤
│ (c) EARLINET │ (d) CHM      │
│  curtain     │  curtain     │
└──────────────┴──────────────┘
```

#### 12.2 Tile 1 (a): Median Profile ± IQR

**Algorithm (lines 40–51):**
```matlab
plot_med_iqr(ax1, bE(:, zMask), z(zMask), colE, 'EARLINET')
plot_med_iqr(ax1, bC(:, zMask), z(zMask), colC, chmLabel)
```

Where:
```matlab
med = median(B, 1, 'omitnan')
q = prctile(B, [25 75], 1)
fill([q(1,:) fliplr(q(2,:))], [z fliplr(z)], color, alpha=0.15)
plot(med, z, '-', linewidth=2)
```

#### 12.3 Tile 2 (b): Binned Density Scatter

**Algorithm (lines 54–77):**
```matlab
band = z >= zMin & z <= zMax
a = bC(:, band); b = bE(:, band)    % a = CHM (y), b = EARLINET (x)
m = isfinite(a) & isfinite(b) & a > 0 & b > 0
a = a(m); b = b(m)
binscatter(log10(b), log10(a), 60)   % bin in log space, 60 bins
plot 1:1 line (black dashed)
xlim = ylim = [log10(max(min([a;b]), 1e-2)), log10(prctile([a;b], 99.8))]
xticks = yticks = ceil(lim(1)):floor(lim(2))
```
- **Axes:** log₁₀ scale (both)
- **Data:** matched profiles within comparison band

#### 12.4 Tiles 3 & 4 (c, d): Curtains

**Algorithm (lines 79–88):**
```matlab
[ts, si] = sort(R.time)
bE = R.betaE(si, :)
bC = R.betaC(si, :)
np = numel(ts)

for ax in [ax3, ax4]:
    B = bE or bC (selected);  B = B';  B(B < 1e-3) = 1e-3
    pcolor(1:np, z, log10(abs(B)))
    shading('flat')
    clim(clim_)
    xticks = round(linspace(1, np, min(5, np)))
    xticklabels = datestr(ts(xticks), 'yy-mm-dd')
```
- **X-axis:** profile index (sorted by time), with date labels
- **Y-axis:** altitude [m AGL]
- **Color:** log₁₀(β)

---

## SUMMARY TABLE: Key Parameters & Thresholds

| Parameter | Value | Unit | Where |
|-----------|-------|------|-------|
| Default comparison band | 500–3000 | m AGL | `zMin`, `zMax` |
| EARLINET comparison band (CHM) | 500–5000 | m AGL | `paper_val_earlinet.m` |
| Cloud/rain/fog exclusion margin | ±15 | min | `screen_profiles` |
| SNR window (hourly bins) | 60 | min | `screen_profiles` |
| Night threshold (SZA) | 95 | deg | `szaThreshold` |
| Deep-night threshold (background est.) | 108 | deg | `estimate_background_profile` |
| Background smoothing window | 21 | gates | `estimate_background_profile` |
| Background min range | 2400 | m AGL | `bgMinRangeAGL` |
| Ångström exponent (default) | 1.0 | — | `angstromExponent` |
| Temporal bin (default) | 60 | min | `dtBin` |
| Rayleigh outlier rejection threshold | 4 | MAD (log space) | `load_rayleigh_python_kalman` |
| Rayleigh min nights | 3 | — | `load_rayleigh_python_kalman` |
| EARLINET assumed lidar ratio | 50 | sr | `read_earlinet_att_backscatter` |
| EARLINET uniform range grid | 0:15:15000 | m AGL | `read_earlinet_att_backscatter` |
| US-1976 atm. uniform grid (molaer) | 0:30:15000 | m AGL | `convert_wavelength_molaer` |
| WV correction band | 900–920 | nm | `apply_wv_correction_to_L2` |
| WV top gates for background | 300 | m | `apply_background_correction` |
| WV relative variance threshold | 1.0 | — | `apply_background_correction` (cirrus mask) |
| Scatter plot subsample | 6000 | — | `paper_val_figure_payerne` |
| EARLINET matching window | 30 | min | `paper_val_earlinet` |
| Histogram percentile | 99.0 | % | `paper_val_figure_payerne` |
| Pcolor floor value | 1e-3 | Mm⁻¹ sr⁻¹ | `paper_val_figure_payerne` |

---

This specification captures all formulas, thresholds, units, file layouts, and processing steps needed to faithfully reproduce the MATLAB pipeline in Python.