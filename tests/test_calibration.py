"""
Basic tests for the rayleigh_calibration package.
"""

import pytest
import numpy as np
from pathlib import Path

from rayleigh_calibration.config import (
    InstrumentType,
    InstrumentInfo,
    CalibrationOptions,
    CalibrationResult,
)
from rayleigh_calibration.atmosphere import (
    calculate_refractive_index,
    calculate_rayleigh_phase_function,
    calculate_molecular_properties,
)


class TestInstrumentType:
    """Tests for InstrumentType enum."""

    def test_wavelengths(self):
        """Test that wavelengths are correctly defined."""
        assert InstrumentType.CHM15k.wavelength_nm == 1064.0
        assert InstrumentType.CL51.wavelength_nm == 910.0
        assert InstrumentType.MINI_MPL.wavelength_nm == 532.0

    def test_supports_calibration(self):
        """Test calibration support flags."""
        assert InstrumentType.CHM15k.supports_calibration is True
        assert InstrumentType.CL31.supports_calibration is False
        assert InstrumentType.CL51.supports_calibration is False
        assert InstrumentType.CL61.supports_calibration is True

    def test_no_cloud_value(self):
        """Test no-cloud flag values."""
        assert InstrumentType.CHM15k.no_cloud_value == -9.0
        assert InstrumentType.MINI_MPL.no_cloud_value == -999.9


class TestInstrumentInfo:
    """Tests for InstrumentInfo dataclass."""

    def test_from_dict(self):
        """Test creation from dictionary."""
        data = {
            "SiteName": "Test Station",
            "WMO": "0-20000-0-12345",
            "Identifier": "A",
            "Type": "CHM15k",
            "Latitude": "46.82",
            "Longitude": "6.95",
            "Altitude": "491",
            "Calibrated": "1",
        }
        info = InstrumentInfo.from_dict(data)
        
        assert info.site_name == "Test Station"
        assert info.wmo_id == "0-20000-0-12345"
        assert info.instrument_type == InstrumentType.CHM15k
        assert info.latitude == 46.82
        assert info.calibrated is True

    def test_to_legacy_dict(self):
        """Test conversion back to legacy format."""
        info = InstrumentInfo(
            site_name="Test",
            wmo_id="12345",
            identifier="A",
            instrument_type=InstrumentType.CHM15k,
            latitude=46.0,
            longitude=7.0,
            altitude=500.0,
        )
        legacy = info.to_legacy_dict()
        
        assert legacy["SiteName"] == "Test"
        assert legacy["Type"] == "CHM15k"


class TestCalibrationResult:
    """Tests for CalibrationResult dataclass."""

    def test_is_successful(self):
        """Test success detection."""
        success = CalibrationResult(lidar_constant=1e12, flag=1, uncertainty=1e10)
        assert success.is_successful is True
        
        partial = CalibrationResult(lidar_constant=1e12, flag=0.5, uncertainty=1e10)
        assert partial.is_successful is True
        
        failed = CalibrationResult(lidar_constant=-1, flag=-1, uncertainty=0)
        assert failed.is_successful is False

    def test_flag_meaning(self):
        """Test flag meaning strings."""
        result = CalibrationResult(lidar_constant=1e12, flag=1, uncertainty=1e10)
        assert result.flag_meaning == "Successful"
        
        result = CalibrationResult(lidar_constant=-1, flag=-1, uncertainty=0)
        assert result.flag_meaning == "Not a clear night"


class TestAtmosphere:
    """Tests for atmospheric physics functions."""

    def test_refractive_index(self):
        """Test refractive index calculation."""
        # At 1064 nm, refractive index should be close to 1
        m = calculate_refractive_index(1064e-9)
        assert 1.0 < m < 1.001
        
        # Shorter wavelength should have higher refractive index
        m_532 = calculate_refractive_index(532e-9)
        m_1064 = calculate_refractive_index(1064e-9)
        assert m_532 > m_1064

    def test_rayleigh_phase_function(self):
        """Test Rayleigh phase function calculation."""
        p = calculate_rayleigh_phase_function()
        # Phase function at 180° should be positive and around 1.5
        assert 1.0 < p < 2.0

    def test_molecular_properties_shape(self):
        """Test that molecular properties have correct shapes."""
        temperature = np.full(100, 288.0)  # K
        pressure = np.full(100, 101325.0)  # Pa
        range_alc = np.arange(0, 10000, 100)  # 0-10 km
        wavelength = 1064e-9
        
        props = calculate_molecular_properties(
            temperature, pressure, range_alc, wavelength
        )
        
        assert props.beta_mol.shape == (100,)
        assert props.alpha_mol.shape == (100,)
        assert props.transmission.shape == (100,)
        
        # Check physical constraints
        assert np.all(props.beta_mol > 0)
        assert np.all(props.alpha_mol > 0)
        assert np.all(props.transmission <= 1)
        assert np.all(props.transmission > 0)

    def test_molecular_backscatter_decreases_with_altitude(self):
        """Test that molecular backscatter decreases with altitude."""
        # Standard atmosphere: pressure decreases with altitude
        altitude = np.arange(0, 10000, 100)
        temperature = 288.0 - 0.0065 * altitude  # Lapse rate
        pressure = 101325 * (temperature / 288.0) ** 5.256
        
        props = calculate_molecular_properties(
            temperature, pressure, altitude, 1064e-9
        )
        
        # Beta_mol should decrease with altitude (lower pressure)
        assert props.beta_mol[0] > props.beta_mol[-1]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
