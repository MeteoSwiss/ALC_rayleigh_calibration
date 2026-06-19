# Rayleigh Calibration for MATLAB

A function-based MATLAB package for performing Rayleigh (molecular) calibration of ceilometers and lidars with parallel processing support.

## Features

- **Multi-instrument support**: CHM15k, CHM8k, CL51, CL61, Mini-MPL
- **Parallel processing**: Uses `parfor` for batch calibration
- **Function-based**: Simple functions, no classes required
- **JSON configuration**: Compatible with E-PROFILE configuration files
- **Vectorized calculations**: Optimized for MATLAB

## Requirements

- MATLAB R2016b or later
- Parallel Computing Toolbox (optional, for parallel processing)

## Installation

1. Extract the package to your desired location
2. Add the folder to MATLAB path:
```matlab
addpath('/path/to/rayleigh_calibration_matlab');
```

## Quick Start

### Run calibration for all stations (yesterday):
```matlab
main_rayleigh_calibration
```

### Edit configuration in the script:
```matlab
% In main_rayleigh_calibration.m:
INSTRUMENTS_FILE = 'E-PROFILE_instruments_L2_auto.json';
OPTIONS_FILE = 'E-PROFILE_options_calibration_rayleigh.json';
DATE_TO_PROCESS = '20240115';  % or '' for yesterday
USE_PARALLEL = true;
NUM_WORKERS = 8;  % 0 = auto-detect
```

### Programmatic usage:
```matlab
% Load configurations
instruments = loadInstruments('E-PROFILE_instruments_L2_auto.json');
options = loadCalibrationOptions('E-PROFILE_options_calibration_rayleigh.json');

% Enable parallel processing
options.UseParallel = true;
options.NumWorkers = 0;  % auto-detect

% Run batch calibration with PARFOR
results = runCalibrationBatch('20240115', instruments, options);

% Check results
for i = 1:length(results)
    fprintf('%s: CL = %.3e (Flag: %.1f)\n', ...
        results(i).SiteName, results(i).LidarConstant, results(i).Flag);
end
```

### Single instrument calibration:
```matlab
% Get one instrument
instruments = loadInstruments('instruments.json');
info = instruments(1);

% Load options
options = loadCalibrationOptions('options.json');

% Run calibration
result = calibrateRayleigh('20240115', info, options);

if result.Flag >= 0.5
    fprintf('Success! CL = %.4e +/- %.4e\n', ...
        result.LidarConstant, result.Uncertainty);
else
    fprintf('Failed: %s\n', result.Message);
end
```

## File Structure

```
rayleigh_calibration_matlab/
├── main_rayleigh_calibration.m   % Main entry script
├── loadInstruments.m             % Load instrument JSON
├── loadCalibrationOptions.m      % Load options JSON
├── loadL1Data.m                  % Load NetCDF L1 files
├── buildFilePaths.m              % Build file paths
├── filterTimeRange.m             % Filter to nighttime
├── filterCloudyProfiles.m        % Remove cloudy data
├── loadStandardAtmosphere.m      % Load/approximate std atm
├── calcMolecularProperties.m     % Molecular backscatter
├── findOptimalMolecularWindow.m  % Rayleigh fit grid search
├── klettInversion.m              % Klett aerosol retrieval
├── calcLidarConstant.m           % Calculate lidar constant
├── calibrateRayleigh.m           % Single station calibration
├── runCalibrationBatch.m         % Parallel batch processing
└── README.md
```

## Flag Meanings

| Flag | Meaning |
|------|---------|
| 1 | Successful calibration |
| 0.5 | Partially clear night |
| 0 | No data available |
| -1 | Not a clear night |
| -2 | Signal not proportional to molecular |
| -3 | Poor method agreement |
| -4 | Missing model data |
| -5 | RCS contains only NaN |
| -6 | Uncertainty exceeds value |
| -7 | Negative Rayleigh fit |
| -8 | Fit intercept exceeds slope |

## Parallel Processing

The package automatically uses `parfor` when:
- Parallel Computing Toolbox is available
- `options.UseParallel = true`
- More than one instrument needs calibration

```matlab
% Enable with specific workers
options.UseParallel = true;
options.NumWorkers = 8;

% Or auto-detect
options.NumWorkers = 0;
```

## Configuration Files

Uses the same JSON format as the Python version - fully compatible with existing E-PROFILE configurations.

### Options JSON example:
```json
{
    "folder_root": "/data/L1_FILES/",
    "folder_output": "/data/calibration/",
    "hour_min": 20,
    "hour_max": 4,
    "min_time_range": 3,
    "z_low_cloud": 4000,
    "LRaer": 52,
    "threshold_quality": 15,
    "use_std_atm": 1
}
```

## License

MIT License
