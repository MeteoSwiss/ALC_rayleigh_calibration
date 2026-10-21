# -*- coding: utf-8 -*-
"""Tests of the polar-summer (low-sun) fallback of the Rayleigh night selection (v2.3 option
``polar_night_fallback``, calibration/rayleigh/polar.py).

  * Option OFF (the default) -> the operational eprof_v2 / v2.2 results are bit-identical: the
    bundled clear-night fixtures give exactly the constants recorded on v2.3-dev before the change,
    and the polar module is never even called.
  * Option ON but a real night exists -> still bit-identical (the fallback only triggers on nights
    without darkness).
  * A synthetic high-latitude solstice night (no profile reaches SZA 100) -> operational flag 0
    ("No profiles in nighttime window"); with the option: the fallback triggers, the night is
    calibrated on the lowest-sun profiles, flagged LOW_SUN_FLAG (0.25), and the true constant is
    recovered.
  * Sun too high (no window long enough above polar_min_sza) and a large residual solar background
    are rejected with LOW_SUN_REJECT_FLAG (-14).
  * Geometry: the closed-form dark duration and lowest sun agree with the operational per-profile
    selection of filter_time_range.
"""
from __future__ import annotations

import copy
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType, DataLevel
from calibration.flags import (FLAG_MEANINGS, LOW_SUN_FLAG, LOW_SUN_REJECT_FLAG, SUCCESS_FLAGS,
                               flag_label, is_success)
from calibration.io.data_loader import CeilometerData, filter_time_range
from calibration.rayleigh import polar
from calibration.rayleigh.atmosphere import (DEFAULT_STANDARD_ATMOSPHERE, load_standard_atmosphere,
                                             calculate_molecular_properties)

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "examples" / "data"
OPTIONS = REPO / "options.json"

# Constants produced by v2.3-dev (35e7c23) BEFORE the polar code was added, on the bundled L2
# clear-night fixtures with options.json (identical for eprof_v2 and eprof_v2.2):
# (wmo, date, type) -> (flag, lidar_constant, uncertainty, bottom_height, top_height)
REFERENCE = {
    ("0-20000-0-06610", "20260225", InstrumentType.CHM15k):
        (1, 492307124267.0247, 31669922985.400066, 2706.9575130939443, 5583.010503053656),
    ("0-20000-0-07014", "20260424", InstrumentType.MINI_MPL):
        (1, 262415.5271364225, 16687.475435226414, 3464.150001525879, 7304.150001525879),
}


# --------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------
def _fixture_info(wmo, itype, date):
    p = DATA / "L2" / wmo / "2026" / date[4:6] / f"L2_{wmo}_A{date}.nc"
    with Dataset(p) as ds:
        lat = float(ds.variables["station_latitude"][...])
        lon = float(ds.variables["station_longitude"][...])
        alt = float(ds.variables["station_altitude"][...])
    return InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier="A", instrument_type=itype,
                          latitude=lat, longitude=lon, altitude=alt)


def _fixture_options(tmp_path, method):
    o = CalibrationOptions.from_json(OPTIONS)
    o.data_level = DataLevel.L2_DAILY
    o.folder_root = DATA / "L2"
    o.folder_output = tmp_path
    o.plot_main = False
    o.plot_all = False
    o.molecular_method = method
    return o


def _synthetic_night(lat, lon, date_str, c_true=3.0e11, delta_frac=0.0, seed=0, dt_s=60.0):
    """A clear, purely molecular CHM15k-like record over the night [D-1, D] (48 h, 60-s profiles).

    rcs = C * beta_mol * T_mol^2 (+ photon-like noise constant in raw signal, i.e. ~r^2 in rcs),
    plus an optional additive residual baseline ``delta`` constant in raw signal (the daylight
    under-subtraction model), expressed as a fraction of the molecular signal at 4 km.
    """
    rng = np.random.default_rng(seed)
    z = np.arange(1, 501) * 30.0                                 # 30 m .. 15 km
    alt = 10.0
    atm = load_standard_atmosphere(DEFAULT_STANDARD_ATMOSPHERE, z + alt)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, z, 1064e-9)
    rcs0 = c_true * mol.beta_att_mol                             # range-corrected, noise-free
    sig4 = float(np.interp(4000.0, z, rcs0 / z ** 2))            # raw molecular signal at 4 km
    d1 = datetime.strptime(date_str, "%Y%m%d") - timedelta(days=1)
    n = int(48 * 3600 / dt_s)
    times = [d1 + timedelta(seconds=dt_s * (k + 0.5)) for k in range(n)]
    sigma_raw = 0.5 * sig4                                       # per-profile SNR 2 at 4 km
    noise = rng.normal(0.0, sigma_raw, (n, z.size)) * z[None, :] ** 2
    rcs = rcs0[None, :] + noise + delta_frac * sig4 * z[None, :] ** 2
    t_days = np.array([(t - datetime(1970, 1, 1)).total_seconds() / 86400.0 for t in times])
    nan = np.full(n, np.nan)
    return CeilometerData(
        time=t_days, time_datetime=times, hours_since_start=(t_days - t_days.min()) * 24.0,
        range_alc=z, altitude=alt, altitude_grid=z + alt, latitude=lat, longitude=lon,
        rcs=rcs, cbh=np.full((n, 3), -9.0),
        temperature_optical_module=nan.copy(), window_transmission=nan.copy(),
        status_laser=nan.copy(), status_detector=nan.copy(), laser_life_time=nan.copy(),
        calibration_pulse=nan.copy(), vertical_visibility=nan.copy())


def _synthetic_options(tmp_path, **polar_kw):
    o = CalibrationOptions.from_json(OPTIONS)
    o.data_level = DataLevel.L1
    o.folder_output = tmp_path
    o.plot_main = False
    o.plot_all = False
    o.apply_wv_correction = False
    for k, v in polar_kw.items():
        setattr(o, k, v)
    return o


def _synthetic_info(lat, lon):
    return InstrumentInfo(site_name="SYNTH", wmo_id="0-99999-0-POLAR", identifier="A",
                          instrument_type=InstrumentType.CHM15k, latitude=lat, longitude=lon,
                          altitude=10.0)


# --------------------------------------------------------------------------------------------
# option OFF: bit-identical
# --------------------------------------------------------------------------------------------
def test_option_defaults_off():
    o = CalibrationOptions()
    assert o.polar_night_fallback is False
    oj = CalibrationOptions.from_json(OPTIONS)
    assert oj.polar_night_fallback is False
    assert oj.sza_night_threshold == 100.0


@pytest.mark.parametrize("method", ["eprof_v2", "eprof_v2.2"])
@pytest.mark.parametrize("key", list(REFERENCE), ids=lambda k: k[2].value)
def test_off_bit_identical_to_reference(tmp_path, monkeypatch, key, method):
    """Option off: exactly the pre-change constants, and the polar module is never called."""
    def _boom(*a, **k):
        raise AssertionError("polar code must not run when the option is off")
    monkeypatch.setattr(polar, "apply_low_sun_fallback", _boom)
    monkeypatch.setattr(polar, "low_sun_residual_qc", _boom)
    wmo, date, itype = key
    r = calibrate_rayleigh(date, _fixture_info(wmo, itype, date), _fixture_options(tmp_path, method))
    flag, c, u, zb, zt = REFERENCE[key]
    assert r.flag == flag
    assert r.lidar_constant == c
    assert r.uncertainty == u
    assert r.calibration_bottom_height == zb
    assert r.calibration_top_height == zt
    assert r.message.startswith("Successful")


@pytest.mark.parametrize("key", list(REFERENCE), ids=lambda k: k[2].value)
def test_on_with_real_night_is_bit_identical(tmp_path, key):
    """Option ON on a night WITH darkness: the fallback does not trigger, nothing changes."""
    wmo, date, itype = key
    o = _fixture_options(tmp_path, "eprof_v2")
    o.polar_night_fallback = True
    fi = {}
    r = calibrate_rayleigh(date, _fixture_info(wmo, itype, date), o, fit_inputs_out=fi)
    flag, c, u, zb, zt = REFERENCE[key]
    assert (r.flag, r.lidar_constant, r.uncertainty) == (flag, c, u)
    assert (r.calibration_bottom_height, r.calibration_top_height) == (zb, zt)
    assert "low_sun" not in fi


# --------------------------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize("lat,date", [(70.0, "20250621"), (60.1, "20250621"), (52.0, "20250621"),
                                      (57.4, "20250601"), (46.8, "20251221"), (-45.0, "20251221"),
                                      (76.5, "20250315")])
def test_geometry_matches_operational_selection(lat, date):
    """Closed-form lowest sun and dark duration == what filter_time_range keeps profile by profile."""
    lon = 10.0
    d = _synthetic_night(lat, lon, date, dt_s=60.0)
    o = CalibrationOptions()
    for thr in (100.0, 96.0, 90.0, 80.0):
        o.sza_night_threshold = thr
        kept = filter_time_range(d, date, o)
        geo = polar.night_geometry(date, lat, thr)
        assert abs(len(kept.time) / 60.0 - geo.dark_hours) <= 2.0 / 60.0 + 1e-9, (lat, date, thr)
    # the lowest sun of the night is reached at solar midnight: 180 - |lat + decl|
    geo = polar.night_geometry(date, lat, 100.0)
    decl = np.degrees(polar.declination_rad(date))
    assert geo.sza_max == pytest.approx(180.0 - abs(lat + decl), abs=1e-9)


def test_floor_window_rules():
    lat, date = 69.3, "20250621"
    geo = polar.night_geometry(date, lat, 100.0)
    assert geo.dark_hours == 0.0 and geo.sza_max < 90.0
    # delta only: the floor is SZA_max - delta
    f = polar.low_sun_floor(geo, lat, 3.0, 0.0, 0.0)
    assert f == pytest.approx(geo.sza_max - 3.0)
    # min_hours only: the window SZA > floor lasts exactly min_hours
    f = polar.low_sun_floor(geo, lat, 0.0, 4.0, 0.0)
    assert polar.hours_above_sza(lat, geo.decl_rad, f) == pytest.approx(4.0, abs=1e-6)
    # both: the larger window wins; min_sza caps it
    f_both = polar.low_sun_floor(geo, lat, 3.0, 4.0, 0.0)
    assert f_both == pytest.approx(min(geo.sza_max - 3.0, f))
    assert polar.low_sun_floor(geo, lat, 3.0, 4.0, 86.0) == 86.0


# --------------------------------------------------------------------------------------------
# synthetic polar-summer night: the fallback triggers, is flagged, and recovers C
# --------------------------------------------------------------------------------------------
LAT, LON, DATE, C_TRUE = 70.0, 15.0, "20250621", 3.0e11


def test_synthetic_polar_night_off_has_no_profiles(tmp_path):
    d = _synthetic_night(LAT, LON, DATE, C_TRUE)
    r = calibrate_rayleigh(DATE, _synthetic_info(LAT, LON), _synthetic_options(tmp_path),
                           preloaded_data=d)
    assert r.flag == 0
    assert "nighttime window" in r.message


def test_synthetic_polar_night_fallback_flagged_and_unbiased(tmp_path):
    d = _synthetic_night(LAT, LON, DATE, C_TRUE)
    o = _synthetic_options(tmp_path, polar_night_fallback=True)
    fi = {}
    r = calibrate_rayleigh(DATE, _synthetic_info(LAT, LON), o, preloaded_data=d,
                           fit_inputs_out=fi)
    assert r.flag == LOW_SUN_FLAG, r.message
    assert r.message.startswith("Low-sun (polar summer)")
    assert abs(r.lidar_constant / C_TRUE - 1.0) < 0.03, r.lidar_constant / C_TRUE
    ls = fi["low_sun"]
    assert ls.dark_hours == 0.0 and not ls.rejected
    assert ls.sza_floor < ls.sza_max < 100.0
    assert ls.window_hours >= o.polar_min_hours - 1e-9
    assert ls.n_profiles == pytest.approx(60.0 * ls.window_hours, abs=2)
    assert abs(ls.residual_pct) < 3.0          # no solar residual was simulated
    assert not is_success(r.flag)              # kept apart from 1/0.5 until promoted


def test_synthetic_sun_too_high_rejected(tmp_path):
    d = _synthetic_night(LAT, LON, DATE, C_TRUE)
    o = _synthetic_options(tmp_path, polar_night_fallback=True, polar_min_sza=88.0)
    r = calibrate_rayleigh(DATE, _synthetic_info(LAT, LON), o, preloaded_data=d)
    assert r.flag == LOW_SUN_REJECT_FLAG
    assert "sun too high" in r.message


def test_synthetic_solar_residual_rejected(tmp_path):
    """An additive residual baseline (daylight under-subtraction) trips the far-range QC."""
    d = _synthetic_night(LAT, LON, DATE, C_TRUE, delta_frac=0.05)
    o = _synthetic_options(tmp_path, polar_night_fallback=True, polar_max_residual_pct=1.0)
    fi = {}
    r = calibrate_rayleigh(DATE, _synthetic_info(LAT, LON), o, preloaded_data=d, fit_inputs_out=fi)
    assert r.flag == LOW_SUN_REJECT_FLAG, r.message
    assert "far-range residual" in r.message
    assert fi["low_sun"].residual_pct > 1.0


def test_far_range_residual_recovers_offset():
    z = np.arange(1, 501) * 30.0
    p_mol = np.exp(-z / 8000.0) / z ** 2
    a, zc = 2.0e11, 4000.0
    delta = 0.05 * a * float(np.interp(zc, z, p_mol))
    bias, x_far = polar.far_range_residual(a * p_mol + delta, p_mol, z, a, zc)
    assert bias == pytest.approx(5.0, rel=1e-9)
    assert x_far > 1.0
    assert np.isnan(polar.far_range_residual(a * p_mol, p_mol, z, np.nan, zc)[0])


# --------------------------------------------------------------------------------------------
# flags
# --------------------------------------------------------------------------------------------
def test_flags_distinct_and_labelled():
    assert LOW_SUN_FLAG not in SUCCESS_FLAGS and LOW_SUN_REJECT_FLAG not in SUCCESS_FLAGS
    assert flag_label(LOW_SUN_FLAG) == FLAG_MEANINGS[0.25] != FLAG_MEANINGS[0]
    assert flag_label(LOW_SUN_REJECT_FLAG) == FLAG_MEANINGS[-14]
    # the pre-existing labels are unchanged
    for f in (1, 1.0, 0.5, 0, -1, -3, -11, -20, -99):
        assert flag_label(f) == FLAG_MEANINGS[0.5 if f == 0.5 else int(f)]
    # the dashboard mirror labels the new flags identically
    from monitoring import config as mcfg
    assert mcfg.flag_label(0.25, "rayleigh") == flag_label(0.25)
    assert mcfg.flag_label(-14, "rayleigh") == flag_label(-14)
    assert mcfg.cal_class(0.25) == "ok" and mcfg.cal_class(-14) == "nonight"


def test_options_json_roundtrip(tmp_path):
    import json
    cfg = json.loads(OPTIONS.read_text(encoding="utf-8"))
    cfg.update(polar_night_fallback=1, polar_min_sza=80, polar_residual_band_m=[7000, 11000])
    p = tmp_path / "o.json"
    p.write_text(json.dumps(cfg), encoding="utf-8")
    o = CalibrationOptions.from_json(p)
    assert o.polar_night_fallback is True and o.polar_min_sza == 80.0
    assert o.polar_residual_band_m == (7000.0, 11000.0)
    o2 = copy.copy(o)
    o2.polar_night_fallback = False
    assert o.polar_night_fallback is True
