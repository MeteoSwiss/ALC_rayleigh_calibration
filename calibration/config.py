"""
Configuration dataclasses for Rayleigh calibration.

This module defines typed configuration objects that replace the previous
JSON-based configuration with proper type hints and validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional
import json


class DataLevel(str, Enum):
    """Input data product the calibration reads from.

    - ``L1``        : native L1 daily files with the range-corrected signal
                      ``rcs_0`` (e.g. ``D:/E-PROFILE_L1``, ``F:/E_PROFILE_ALC/L1_FILES``).
                      Layout: ``<root>/<WMO>/YYYY/MM/L1_<WMO>_<id><YYYYMMDD>.nc``.
    - ``L2_daily``  : daily L2 files with ``attenuated_backscatter_0`` +
                      ``calibration_constant_0`` (e.g. ``D:/E-PROFILE_L2_2021-2025``).
                      Layout: ``<root>/<WMO>/YYYY/MM/L2_<WMO>_<id><YYYYMMDD>.nc``.
    - ``L2_monthly``: monthly L2 files, same variables as ``L2_daily`` but one file
                      per month (e.g. ``A:/E-PROFILE_L2_monthly``).
                      Layout: ``<root>/<WMO>/YYYY/L2_<WMO>_<id><YYYYMM>.nc``.

    For both L2 levels the range-corrected signal is reconstructed as
    ``rcs = attenuated_backscatter_0 * calibration_constant_0 * 1e-6`` (fixed factor,
    matching the MATLAB reference ``loadL2Data.m``; the stored attenuated backscatter is
    always micro-scaled even when the ``units`` attribute mislabels it), recovering the
    L1-equivalent signal so the rest of the pipeline is unchanged.
    """
    L1 = "L1"
    L2_DAILY = "L2_daily"
    L2_MONTHLY = "L2_monthly"
    # Native manufacturer/Cloudnet RAW NetCDF (CHM15k 'beta_raw' daily files; CL61 'beta_att'
    # Cloudnet 5-min files). Layout: <root>/<WMO>/YYYYMMDD/*.nc (one daily file or many 5-min
    # files per day-folder; the reader concatenates them). See load_raw_data.
    RAW = "RAW"


class InstrumentType(str, Enum):
    """Supported ceilometer/lidar instrument types."""
    CHM15k = "CHM15k"
    CHM8k = "CHM8k"
    CL31 = "CL31"
    CL51 = "CL51"
    CL61 = "CL61"
    MINI_MPL = "Mini-MPL"
    MPL = "MPL"

    @property
    def wavelength_nm(self) -> float:
        """Return the laser wavelength in nanometers for this instrument type."""
        wavelengths = {
            InstrumentType.CHM15k: 1064.0,
            InstrumentType.CHM8k: 1064.0,
            InstrumentType.CL31: 910.0,
            InstrumentType.CL51: 910.0,
            InstrumentType.CL61: 910.0,
            InstrumentType.MINI_MPL: 532.0,
            InstrumentType.MPL: 532.0,
        }
        return wavelengths[self]

    @property
    def supports_calibration(self) -> bool:
        """Check if this instrument type supports Rayleigh calibration."""
        # CL31 and CL51 have signal distortion issues
        return self in {
            InstrumentType.CHM15k,
            InstrumentType.CHM8k,
            InstrumentType.CL61,
            InstrumentType.MINI_MPL,
            InstrumentType.MPL,
        }

    @property
    def no_cloud_value(self) -> float:
        """Return the no-cloud flag value for this instrument type."""
        if self == InstrumentType.MINI_MPL:
            return -999.9
        return -9.0

    @property
    def lidar_constant_units(self) -> str:
        """Return the units for lidar constant for this instrument type."""
        if self in {InstrumentType.CL31, InstrumentType.CL51, InstrumentType.CL61}:
            return "V*m^3/sr"
        return "counts/s*m^3/sr"


@dataclass
class InstrumentInfo:
    """Information about a single ceilometer/lidar instrument."""
    site_name: str
    wmo_id: str
    identifier: str
    instrument_type: InstrumentType
    latitude: float
    longitude: float
    altitude: float
    calibrated: bool = True       # used by the pipeline (skip uncalibrated instruments)
    serial: str = ""              # carried from the file (instrument_serial_number)

    @classmethod
    def from_dict(cls, data: dict) -> InstrumentInfo:
        """Create InstrumentInfo from a dictionary (e.g., from instruments.json).

        Only the fields the calibration code uses are read. The legacy Reference/FLength/
        NWS/Status fields were dropped from instruments.json (not in the source files and
        not used); ``.get`` defaults keep older manifests loading unchanged.
        """
        return cls(
            site_name=data["SiteName"],
            wmo_id=data["WMO"],
            identifier=data.get("Identifier", "A"),
            instrument_type=InstrumentType(data["Type"]),
            latitude=float(data["Latitude"]),
            longitude=float(data["Longitude"]),
            altitude=float(data.get("Altitude", 0)),
            calibrated=str(data.get("Calibrated", "1")) == "1",
            serial=data.get("Serial", ""),
        )

    def to_legacy_dict(self) -> dict:
        """Convert back to the instruments.json dictionary format."""
        return {
            "WMO": self.wmo_id,
            "Identifier": self.identifier,
            "Type": self.instrument_type.value,
            "SiteName": self.site_name,
            "Latitude": self.latitude,
            "Longitude": self.longitude,
            "Altitude": self.altitude,
            "Serial": self.serial,
            "Calibrated": "1" if self.calibrated else "0",
        }


@dataclass
class CalibrationOptions:
    """Configuration options for Rayleigh calibration."""
    # Paths
    folder_root: Path = field(default_factory=lambda: Path("/data/zue/E_PROFILE/ALC/L1_FILES/"))
    folder_output: Path = field(default_factory=lambda: Path("/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh/"))

    # Input data product to read (see DataLevel). folder_root must point at the
    # matching archive. L2 levels reconstruct rcs from attenuated backscatter.
    data_level: DataLevel = DataLevel.L1

    # Time selection (solar time)
    hour_min: int = 20  # Start hour on previous day (clock fallback)
    hour_max: int = 4   # End hour on current day (clock fallback)
    min_time_range: int = 3  # Minimum hours of clear data required
    # Darkness-adaptive night selection (solar zenith angle). When enabled, the
    # nighttime averaging window is the real dark period (SZA > threshold) rather
    # than the fixed solar-clock window -> more dark hours in winter, no twilight
    # in summer, latitude-adaptive. Falls back to hour_min/hour_max if disabled
    # or station coordinates are missing.
    use_sza_night: bool = True
    sza_night_threshold: float = 100.0  # deg; 100 = sun ~10 deg below horizon (nautical)
    # Time collapse before the molecular fit: "mean" (efficient for the weak-signal
    # photon noise at 3-6 km) or "median" (robust but ~1.5x noisier there).
    time_aggregation: str = "mean"
    # Optional pre-averaging before the Rayleigh fit. When set, the calibration reads
    # the full nightly data, then block-averages it in time and range before any
    # filtering/fitting. Disabled by default; the Lindenberg runner enables 60 s / 30 m.
    average_time_s: Optional[float] = None
    average_range_m: Optional[float] = None
    # Drop residual-aerosol/cloud/noise outlier PROFILES before the mean (MAD-based on
    # the molecular-band signal) -> robust without the median's noise penalty.
    screen_profile_outliers: bool = True
    profile_outlier_nmad: float = 4.0

    # Plotting flags. plot_main -> the single compact Rayleigh diagnostics dashboard only.
    # plot_all -> additionally emit the simple per-step RCS plots (time-series + annotated).
    plot_main: bool = False
    plot_all: bool = False

    # Cloud detection
    z_low_cloud: float = 4000.0  # Maximum altitude for low cloud detection (m)
    max_ratio_cloudy: float = 0.5

    # Extinction calculation
    lidar_ratio_aerosol: float = 52.0  # Aerosol lidar ratio (sr)
    z_start_ext: float = 100.0  # Start altitude for extinction calc (m)
    calc_ext_above_molecular: bool = False
    subtract_background: bool = False
    consider_points_lower_than_molecular: bool = True

    # Quality thresholds
    threshold_quality: float = 15.0

    # Atmosphere model
    use_std_atm: bool = True

    # Source of the molecular reference profile (beta_mol proportional to P/T):
    #   'standard' (default) — US Standard 1976 atmosphere (matches the MATLAB
    #                          reference, which always uses the standard atmosphere);
    #   'cams'              — actual CAMS T/p at the site (needs the monthly CAMS
    #                          file, same archive as the WV correction).
    # At mid-latitudes the two differ by <1 % in molecular density. When set to
    # 'cams' this takes precedence over use_std_atm.
    molecular_source: str = "standard"

    # Water-vapor correction (910 nm instruments only). Requires monthly CAMS +
    # the HITRAN cross-section LUT. A 910 nm night without a matching CAMS month is
    # skipped (never calibrated without a valid WV correction).
    apply_wv_correction: bool = False
    cams_folder: Path = field(default_factory=lambda: Path("D:/CAMS/"))
    abs_cs_lookup_table: Path = field(default_factory=lambda: Path(""))
    # Auto-download a missing CAMS file from the ADS instead of skipping the night.
    # Off by default (calibration stays offline and reproducible). When on, a missing
    # CAMS_Beta_*.nc is fetched via calibration.io.download_cams_beta (needs cdsapi +
    # cfgrib + ADS credentials in ~/.cdsapirc; the ADS endpoint is selected automatically).
    auto_download_cams: bool = False
    # What to fetch on a miss: 'day' (CAMS_Beta_<YYYYMMDD>.nc, just the night — light,
    # and the only choice valid for a recent date whose month is not yet complete) or
    # 'month' (CAMS_Beta_<YYYYMM>.nc, the whole completed month — heavy but reused all month).
    cams_download_scope: str = "day"

    # Rayleigh fit window search parameters (not typically changed)
    half_length_options_m: tuple = field(default_factory=lambda: tuple(range(250, 2000, 240)))
    range_start_m: float = 2000.0
    range_end_m: float = 6000.0
    fit_range_increment_bins: int = 8
    # Molecular-window validity gates (fix the high-altitude low-R² selection; see
    # rayleigh_fit.find_optimal_molecular_window). A window is eligible only if it
    # starts above the boundary-layer aerosol (min_window_start_m), has a genuine
    # linear signal~molecular relation (R² >= min_window_r2, cf. MATLAB min_r2_rfit),
    # and its slope is consistent with the pointwise median ratio
    # (<= max_window_rel_error %, rejecting aerosol curvature). The best window is the
    # highest-R² eligible one; if none qualifies the night is flagged non-calibration.
    min_window_start_m: float = 2000.0
    min_window_r2: float = 0.5
    max_window_rel_error: float = 50.0
    # Molecular-window detection strategy, keyed by E-PROF version (see molecular_methods.py):
    # 'eprof_v1.2' (default, production = improved), 'eprof_v1.1' (legacy main, sign-corrected),
    # 'eprof_v0.25' (MATLAB Auto_Calib_25), 'earlinet', 'eprof_v2' (optimal), 'bellini'.
    # Legacy aliases (improved/main/matlab/optimal) are still accepted. The min_window_*/
    # max_window_* gates above apply to 'eprof_v1.2'; other methods use their own defaults.
    molecular_method: str = "eprof_v1.2"

    # E-PROF v1.0 baseline ONLY: reproduce the pre-a4e7140 calibration (the Klett
    # sign error + total-OD reference). Default False = the production, corrected
    # calibration. Used solely to regenerate the historical v1.0 series for the
    # method comparison (see atmosphere.klett_inversion / calibrate_rayleigh).
    sign_error_v10: bool = False

    @classmethod
    def from_json(cls, filepath: Path) -> CalibrationOptions:
        """Load options from a JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)

        return cls(
            folder_root=Path(data.get("folder_root", "/data/zue/E_PROFILE/ALC/L1_FILES/")),
            folder_output=Path(data.get("folder_output", "/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh/")),
            data_level=DataLevel(data.get("data_level", "L1")),
            hour_min=data.get("hour_min", 20),
            hour_max=data.get("hour_max", 4),
            min_time_range=data.get("min_time_range", 3),
            use_sza_night=bool(data.get("use_sza_night", 1)),
            sza_night_threshold=float(data.get("sza_night_threshold", 100.0)),
            time_aggregation=str(data.get("time_aggregation", "mean")),
            average_time_s=(None if data.get("average_time_s", None) is None else float(data.get("average_time_s"))),
            average_range_m=(None if data.get("average_range_m", None) is None else float(data.get("average_range_m"))),
            screen_profile_outliers=bool(data.get("screen_profile_outliers", 1)),
            profile_outlier_nmad=float(data.get("profile_outlier_nmad", 4.0)),
            plot_main=bool(data.get("plot_main", 0)),
            plot_all=bool(data.get("plot_all", 0)),
            z_low_cloud=float(data.get("z_low_cloud", 4000)),
            max_ratio_cloudy=float(data.get("max_ratio_cloudy", 0.5)),
            lidar_ratio_aerosol=float(data.get("LRaer", 52)),
            z_start_ext=float(data.get("z_start_ext", 100)),
            calc_ext_above_molecular=bool(data.get("calc_ext_above_molecular", 0)),
            subtract_background=bool(data.get("subtract_background", 0)),
            consider_points_lower_than_molecular=bool(data.get("consider_points_lower_than_molecular", 1)),
            threshold_quality=float(data.get("threshold_quality", 15)),
            use_std_atm=bool(data.get("use_std_atm", 1)),
            molecular_source=str(data.get("molecular_source", "standard")),
            range_start_m=float(data.get("range_start_m", 2000)),
            range_end_m=float(data.get("range_end_m", 6000)),
            min_window_start_m=float(data.get("min_window_start_m", 2000)),
            min_window_r2=float(data.get("min_window_r2", 0.5)),
            max_window_rel_error=float(data.get("max_window_rel_error", 50.0)),
            molecular_method=str(data.get("molecular_method", "eprof_v1.2")),
            sign_error_v10=bool(data.get("sign_error_v10", 0)),
            apply_wv_correction=bool(data.get("apply_wv_correction", 0)),
            cams_folder=Path(data.get("cams_folder", "D:/CAMS/")),
            abs_cs_lookup_table=Path(data.get("abs_cs_lookup_table", "")),
            auto_download_cams=bool(data.get("auto_download_cams", 0)),
            cams_download_scope=str(data.get("cams_download_scope", "day")),
        )


@dataclass
class CalibrationResult:
    """Result of a Rayleigh calibration attempt."""
    # Wiegner & Geiss (2012) lidar constant C_L = RCS / beta_att (calibrate via
    # beta_att = RCS / C_L). Same definition/units as the cloud product's lidar_constant.
    lidar_constant: float
    flag: int
    uncertainty: float
    calibration_bottom_height: Optional[float] = None
    calibration_top_height: Optional[float] = None
    message: str = ""

    @property
    def is_successful(self) -> bool:
        """Check if calibration was successful."""
        return self.flag == 1 or self.flag == 0.5

    @property
    def flag_meaning(self) -> str:
        """Get human-readable meaning of the flag value."""
        meanings = {
            1: "Successful",
            0.5: "Partially clear night",
            0: "No data",
            -1: "Not a clear night",
            -2: "Signal not proportional to molecular scattering",
            -3: "Bad quality ratio between methods",
            -4: "Missing model data",
            -5: "RCS contains only NaN",
            -6: "Uncertainties higher than CL values",
            -7: "Negative Rayleigh fit",
            -8: "Rayleigh fit issue: |b| > a",
        }
        return meanings.get(int(self.flag), "Unknown flag")


def load_instruments(filepath: Path) -> list[InstrumentInfo]:
    """Load instrument configuration from JSON file."""
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        data = json.load(f)
    return [InstrumentInfo.from_dict(item) for item in data]
