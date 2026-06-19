"""
Robustness tests: the Rayleigh and the liquid-cloud calibration must RUN for every instrument
type at both L1 and L2, producing a reasonable result object (never an exception), and emitting a
suitability *warning* for the combinations the method is not validated for:

  * Rayleigh  — CL31/CL51 (signal distortion) warn; CHM15k/CL61/Mini-MPL are suitable.
  * Cloud     — only the 910 nm CL61 is suitable; CL31/CL51 (distortion) and the non-910 nm
                CHM15k/Mini-MPL warn but still run.

Uses the CL31/CL51/CHM15k/Mini-MPL/CL61 sample files bundled under examples/data (L1 and L2). The
Rayleigh tests run fully offline; the cloud tests need CAMS for the mandatory 910 nm water-vapour
correction (skipped per-file when that month's CAMS_Beta is absent, e.g. in CI). The non-910 nm
instruments do not use the WV correction, so they run without CAMS.
"""
from __future__ import annotations

import glob
import os
import warnings
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
from calibration.cloud import liquid_cloud_calibration, CloudCalConfig

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "examples" / "data"
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
CAMS_DIR = Path("D:/CAMS")
WV_BAND = {"CL31", "CL51", "CL61"}           # 910 nm -> mandatory WV correction (needs CAMS)


def _discover():
    """[(level, instrument_type, file, wmo, ident, lat, lon, alt, yyyymm), ...] from examples/data."""
    out = []
    for level in ("L1", "L2"):
        for f in sorted(glob.glob(str(DATA / level / "*" / "2026" / "*" / "*.nc"))):
            with Dataset(f, "r") as d:
                it = getattr(d, "instrument_type", "?")
                g = lambda n: (float(np.asarray(d.variables[n][:]).ravel()[0])
                               if n in d.variables else 0.0)
                rec = (level, it, f, getattr(d, "wigos_station_id", ""),
                       getattr(d, "instrument_id", "A"),
                       g("station_latitude"), g("station_longitude"), g("station_altitude"))
            out.append(rec + (os.path.basename(f)[-11:-5],))   # YYYYMM
    return out


CASES = _discover()
IDS = [f"{lvl}-{it}" for (lvl, it, *_rest) in CASES]


def _suitability_warned(records):
    return any("indicative" in str(w.message) or "not well-suited" in str(w.message)
               for w in records)


@pytest.mark.skipif(not CASES, reason="no bundled sample data under examples/data")
@pytest.mark.parametrize("case", CASES, ids=IDS)
def test_rayleigh_runs_for_every_instrument(case):
    level, it, f, wmo, ident, lat, lon, alt, _ym = case
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = DATA / level
    o.data_level = DataLevel.L1 if level == "L1" else DataLevel.L2_DAILY
    o.plot_main = False
    o.apply_wv_correction = False                 # keep the "does it run" test offline (no CAMS)
    info = InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier=ident,
                          instrument_type=InstrumentType(it),
                          latitude=lat, longitude=lon, altitude=alt)
    ds = os.path.basename(f)[-11:-3]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        result = calibrate_rayleigh(ds, info, o)   # must not raise
    # a reasonable result object with a numeric flag (calibrated or flagged, both fine)
    assert result is not None
    assert np.isreal(result.flag)
    assert np.isreal(result.lidar_constant)
    # suitability warning iff the instrument is not well-suited to Rayleigh (CL31/CL51)
    expect = not InstrumentType(it).supports_calibration
    assert _suitability_warned(rec) is expect, (it, level, [str(w.message) for w in rec])


@pytest.mark.skipif(not (CASES and WV_LUT.is_file()),
                    reason="no bundled sample data / WV LUT")
@pytest.mark.parametrize("case", CASES, ids=IDS)
def test_cloud_runs_for_every_instrument(case):
    level, it, f, wmo, ident, lat, lon, alt, ym = case
    if it in WV_BAND and not (CAMS_DIR / f"CAMS_Beta_{ym}.nc").is_file():
        pytest.skip(f"CAMS_Beta_{ym}.nc needed for the 910 nm WV correction not present")
    cfg = CloudCalConfig(
        nc_file=str(f), instrument=it, apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=str(CAMS_DIR), abs_cs_lookup_table=str(WV_LUT),
        station_latitude=lat, station_longitude=lon)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = liquid_cloud_calibration(cfg)        # must not raise for any instrument x level
    # a reasonable result: ran to completion; clouds may or may not be present that day
    assert res is not None
    assert int(res.n_profiles) >= 0
    if res.n_profiles > 0:
        assert np.isfinite(res.cal_median) and res.cal_median > 0
    # only the 910 nm CL61 is "suitable"; everything else should warn
    expect = it != "CL61"
    assert _suitability_warned(rec) is expect, (it, level, [str(w.message) for w in rec])
