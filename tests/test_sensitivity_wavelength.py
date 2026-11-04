# -*- coding: utf-8 -*-
"""The sensitivity product must not divide by a zero wavelength.

Two L1 days of the v2.3 release add-ons (Oslo CHM15k, Libus CL51) carried l0_wavelength = 0, and
``alpha_mol0_for_wavelength`` raised ZeroDivisionError from its lambda^-4 scaling, which cost both streams
their sensitivity and OmB products. Now the per-day reader falls back to the type's nominal wavelength when the
file value is not a positive number, and the function itself refuses such a value with a clear message.
"""
import numpy as np
import pytest

from calibration.sensitivity.detection import ALPHA_MOL0, alpha_mol0_for_wavelength


def test_tabulated_and_nearby_wavelengths():
    assert alpha_mol0_for_wavelength(910.0) == ALPHA_MOL0[910]
    assert alpha_mol0_for_wavelength(910.55) == ALPHA_MOL0[910]
    assert alpha_mol0_for_wavelength(1064.5) == ALPHA_MOL0[1064]


@pytest.mark.parametrize("bad", [0.0, -5.0, float("nan")])
def test_non_positive_wavelength_is_refused(bad):
    with pytest.raises(ValueError, match="positive number"):
        alpha_mol0_for_wavelength(bad)


def test_reader_falls_back_to_the_nominal_wavelength():
    from calibration.io.instrument_day import usable_wavelength
    for bad in (0.0, -1.0, float("nan"), None, "x"):
        assert usable_wavelength(bad, 1064.0) == 1064.0
    assert usable_wavelength(910.55, 1064.0) == 910.55
    assert usable_wavelength(np.float32(1064.3), 1064.0) == pytest.approx(1064.3)
