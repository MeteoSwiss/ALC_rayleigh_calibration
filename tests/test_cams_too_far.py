"""Regression: a station OUTSIDE the regional CAMS domain (e.g. New Zealand) must be detected as
'closest CAMS data too far' (-> flag -10), instead of silently using the domain-edge cell and
producing a false 910 nm calibration. Uses a synthetic Europe/N-Atlantic CAMS grid (no big files)."""
import numpy as np
import pytest

xr = pytest.importorskip("xarray")

from calibration.water_vapor_correction.water_vapor import cams_point_too_far, cams_nearest_offset_deg


def _make_cams(tmp_path, lat0=27.0, lat1=74.0, lon0=-27.0, lon1=45.0, step=1.0):
    """Minimal CAMS-like file: just latitude/longitude coords (what the guard reads)."""
    lats = np.arange(lat0, lat1 + step, step)
    lons = np.arange(lon0, lon1 + step, step)
    ds = xr.Dataset({"t": (("latitude", "longitude"), np.zeros((lats.size, lons.size)))},
                    coords={"latitude": lats, "longitude": lons})
    p = tmp_path / "CAMS_Beta_synthetic.nc"
    ds.to_netcdf(p)
    return str(p)


def test_out_of_domain_station_is_too_far(tmp_path):
    f = _make_cams(tmp_path)
    assert cams_point_too_far(f, -45.04, 169.68) is True      # Lauder, New Zealand
    assert cams_point_too_far(f, -33.9, 18.4) is True         # Cape Town
    off, _ = cams_nearest_offset_deg(f, -45.04, 169.68)
    assert off > 100                                          # nearest cell is ~125 deg away


def test_in_domain_stations_are_not_too_far(tmp_path):
    f = _make_cams(tmp_path)
    for lat, lon in [(46.81, 6.94),   # Payerne
                     (52.21, 14.12),  # Lindenberg
                     (50.22, -5.33),  # Camborne
                     (60.0, 10.0)]:   # southern Norway
        assert cams_point_too_far(f, lat, lon) is False, (lat, lon)


def test_nonfinite_location_not_flagged(tmp_path):
    # Unknown coords must not be hard-failed by the guard (let the normal path decide).
    f = _make_cams(tmp_path)
    assert cams_point_too_far(f, float("nan"), 6.9) is False
