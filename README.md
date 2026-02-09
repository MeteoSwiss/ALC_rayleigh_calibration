# Rayleigh Calibration for Automated Lidars and Ceilometers

A modern, maintainable Python package for performing Rayleigh (molecular) calibration of ceilometers and lidars. This is a complete rewrite of the original E-PROFILE calibration scripts with significant improvements in code organization, performance, and maintainability.

## Features

- **Multi-instrument support**: CHM15k, CHM8k, CL51, CL61, Mini-MPL, and more
- **Modern Python**: Type hints, dataclasses, and clean architecture
- **Optimized performance**: Vectorized NumPy operations replace explicit loops
- **Comprehensive logging**: Track calibration progress and diagnose issues
- **Flexible configuration**: JSON-based configuration with Python dataclass validation
- **CF-compliant output**: NetCDF4 output following Climate and Forecast conventions
- **Optional visualization**: Matplotlib-based diagnostic plots

## Installation

### From PyPI (when published)

```bash
pip install rayleigh-calibration
```

### From source

```bash
git clone https://github.com/e-profile/rayleigh-calibration.git
cd rayleigh-calibration
pip install -e ".[all]"
```

### Dependencies

**Required:**
- Python >= 3.9
- NumPy >= 1.20
- SciPy >= 1.7  
- netCDF4 >= 1.5

**Optional:**
- matplotlib >= 3.4 (for plotting)

## Quick Start

### Command Line Usage

Process yesterday's data for all configured stations:

```bash
rayleigh-calibration --verbose
```

Process a specific date:

```bash
rayleigh-calibration --date 20240115 --verbose
```

Process a single station:

```bash
rayleigh-calibration --station 0-20000-0-06610 --date 20240115
```

### Python API

```python
from rayleigh_calibration import (
    calibrate_rayleigh,
    load_instruments,
    CalibrationOptions,
)

# Load configuration
instruments = load_instruments("instruments.json")
options = CalibrationOptions.from_json("options.json")

# Run calibration for a single instrument
result = calibrate_rayleigh("20240115", instruments[0], options)

if result.is_successful:
    print(f"Lidar constant: {result.lidar_constant:.4e}")
    print(f"Uncertainty: {result.uncertainty:.4e}")
else:
    print(f"Calibration failed: {result.flag_meaning}")
```

## Configuration

### Calibration Options (JSON)

```json
{
    "folder_root": "/data/L1_FILES/",
    "folder_output": "/data/calibration/rayleigh/",
    "folder_ECMWF": "/data/ECMWF/",
    "hour_min": 20,
    "hour_max": 4,
    "min_time_range": 3,
    "z_low_cloud": 4000,
    "LRaer": 52,
    "threshold_quality": 15,
    "use_std_atm": 1,
    "plot_main": 0,
    "plot_all": 0
}
```

### Instrument Configuration (JSON)

```json
[
    {
        "SiteName": "PAYERNE, Switzerland",
        "WMO": "0-20000-0-06610",
        "Identifier": "A",
        "Type": "CHM15k",
        "Latitude": 46.82,
        "Longitude": 6.95,
        "Altitude": 491,
        "Calibrated": "1"
    }
]
```

## Output

Calibration results are written to CF-compliant NetCDF4 files with the following structure:

- **time**: Central time of calibration period
- **lidar_constant**: Calculated lidar constant
- **lidar_constant_uncertainty**: Calibration uncertainty
- **calibration_bottom_height**: Bottom of molecular window (m ASL)
- **calibration_top_height**: Top of molecular window (m ASL)
- **calibration_method**: 0 = Rayleigh calibration

## Algorithm Overview

1. **Data Loading**: Read L1 NetCDF files from current and previous day
2. **Time Filtering**: Select nighttime profiles (default: 20:00-04:00 UTC)
3. **Cloud Filtering**: Remove profiles contaminated by low clouds
4. **Atmospheric Model**: Load temperature/pressure profiles (standard atmosphere or ECMWF)
5. **Molecular Calculation**: Calculate molecular backscatter using Bucholtz (1995)
6. **Rayleigh Fit**: Find optimal molecular window via grid search
7. **Klett Inversion**: Retrieve aerosol backscatter and extinction
8. **Lidar Constant**: Calculate using Wiegner & Geiss (2012) method
9. **Validation**: Cross-check slope method vs. Klett method
10. **Output**: Write results to NetCDF

## Flag Meanings

| Flag | Meaning |
|------|---------|
| 1 | Successful calibration |
| 0.5 | Partially clear night (some clouds removed) |
| 0 | No data available |
| -1 | Not a clear night |
| -2 | Signal not proportional to molecular scattering |
| -3 | Poor agreement between calibration methods |
| -4 | Missing atmospheric model data |
| -5 | RCS contains only NaN values |
| -6 | Uncertainty exceeds calibration value |
| -7 | Negative Rayleigh fit slope |
| -8 | Rayleigh fit intercept exceeds slope |

## Migration from Original Code

The new package is a complete rewrite with improved architecture. Key changes:

### Code Organization

| Original | New |
|----------|-----|
| `ALC_CAL_rayleigh_auto.py` | `rayleigh_calibration/main.py` |
| `ALC_cal_rayleigh_tools.py` | Split into: `calibration.py`, `atmosphere.py`, `data_loader.py`, `rayleigh_fit.py`, `output.py` |
| `utils.py` | Removed (replaced by NumPy builtins) |

### Performance Improvements

- **Vectorized calculations**: Molecular backscatter calculation is now ~10x faster
- **Efficient I/O**: NetCDF files are read with proper chunking
- **Memory optimization**: Data arrays are pre-allocated

### Maintainability Improvements

- **Type hints**: All functions have type annotations
- **Dataclasses**: Configuration objects with validation
- **Docstrings**: Comprehensive NumPy-style documentation
- **Logging**: Structured logging instead of print statements
- **Testing**: Unit tests with pytest

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Run the test suite: `pytest`
5. Submit a pull request

## References

- Bucholtz, A. (1995). Rayleigh-scattering calculations for the terrestrial atmosphere. Applied Optics, 34(15), 2765-2773.
- Wiegner, M., & Geiß, A. (2012). Aerosol profiling with the Jenoptik ceilometer CHM15kx. Atmospheric Measurement Techniques, 5(8), 1953-1964.
- E-PROFILE Programme: https://e-profile.eu

## License

MIT License - see LICENSE file for details.

## Changelog

### v2.0.0 (2024)

- Complete rewrite with modern Python practices
- Added type hints and dataclasses
- Vectorized molecular calculations
- Improved error handling and logging
- Added command-line interface
- Split monolithic file into modular architecture
- Added comprehensive documentation

### v1.0.0 (2015-2024)

- Original E-PROFILE implementation by hem
