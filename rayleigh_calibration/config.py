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
    calibrated: bool = True
    reference: bool = False
    serial: str = ""
    focal_length: float = 5.0
    network: str = ""
    status: bool = True

    @classmethod
    def from_dict(cls, data: dict) -> InstrumentInfo:
        """Create InstrumentInfo from a dictionary (e.g., from JSON)."""
        return cls(
            site_name=data["SiteName"],
            wmo_id=data["WMO"],
            identifier=data.get("Identifier", "A"),
            instrument_type=InstrumentType(data["Type"]),
            latitude=float(data["Latitude"]),
            longitude=float(data["Longitude"]),
            altitude=float(data.get("Altitude", 0)),
            calibrated=data.get("Calibrated", "1") == "1",
            reference=data.get("Reference", "0") == "1",
            serial=data.get("Serial", ""),
            focal_length=float(data.get("FLength", 5)),
            network=data.get("NWS", ""),
            status=data.get("Status", "1") == "1",
        )

    def to_legacy_dict(self) -> dict:
        """Convert back to legacy dictionary format for compatibility."""
        return {
            "SiteName": self.site_name,
            "WMO": self.wmo_id,
            "Identifier": self.identifier,
            "Type": self.instrument_type.value,
            "Latitude": self.latitude,
            "Longitude": self.longitude,
            "Altitude": self.altitude,
            "Calibrated": "1" if self.calibrated else "0",
            "Reference": "1" if self.reference else "0",
            "Serial": self.serial,
            "FLength": str(self.focal_length),
            "NWS": self.network,
            "Status": "1" if self.status else "0",
        }


@dataclass
class CalibrationOptions:
    """Configuration options for Rayleigh calibration."""
    # Paths
    folder_ecmwf: Path = field(default_factory=lambda: Path("E:/ECMWF/"))
    folder_root: Path = field(default_factory=lambda: Path("/data/zue/E_PROFILE/ALC/L1_FILES/"))
    folder_output: Path = field(default_factory=lambda: Path("/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh/"))

    # Input data product to read (see DataLevel). folder_root must point at the
    # matching archive. L2 levels reconstruct rcs from attenuated backscatter.
    data_level: DataLevel = DataLevel.L1

    # Time selection (hours UTC)
    hour_min: int = 20  # Start hour on previous day
    hour_max: int = 4   # End hour on current day
    min_time_range: int = 3  # Minimum hours of clear data required

    # Plotting flags
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

    # Water-vapor correction (910 nm instruments only). Requires monthly CAMS +
    # the HITRAN cross-section LUT. A 910 nm night without a matching CAMS month is
    # skipped (never calibrated without a valid WV correction).
    apply_wv_correction: bool = False
    cams_folder: Path = field(default_factory=lambda: Path("D:/CAMS/"))
    abs_cs_lookup_table: Path = field(default_factory=lambda: Path(""))

    # Rayleigh fit window search parameters (not typically changed)
    half_length_options_m: tuple = field(default_factory=lambda: tuple(range(250, 2000, 240)))
    range_start_m: float = 2000.0
    range_end_m: float = 6000.0
    fit_range_increment_bins: int = 8

    @classmethod
    def from_json(cls, filepath: Path) -> CalibrationOptions:
        """Load options from a JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)

        return cls(
            folder_ecmwf=Path(data.get("folder_ECMWF", "E:/ECMWF/")),
            folder_root=Path(data.get("folder_root", "/data/zue/E_PROFILE/ALC/L1_FILES/")),
            folder_output=Path(data.get("folder_output", "/data/pay/REM/ACQ/E_PROFILE_ALC/Calibration/rayleigh/")),
            data_level=DataLevel(data.get("data_level", "L1")),
            hour_min=data.get("hour_min", 20),
            hour_max=data.get("hour_max", 4),
            min_time_range=data.get("min_time_range", 3),
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
            range_start_m=float(data.get("range_start_m", 2000)),
            range_end_m=float(data.get("range_end_m", 6000)),
            apply_wv_correction=bool(data.get("apply_wv_correction", 0)),
            cams_folder=Path(data.get("cams_folder", "D:/CAMS/")),
            abs_cs_lookup_table=Path(data.get("abs_cs_lookup_table", "")),
        )


@dataclass
class CalibrationResult:
    """Result of a Rayleigh calibration attempt."""
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
