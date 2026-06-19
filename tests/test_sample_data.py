"""Calibration tests on the bundled sample data (``examples/data/``).

These exercise the readers and both calibrations on small real fixtures committed to
the repo, so the package can be smoke-tested without the full E-PROFILE archive.

- 1064/532 nm (CHM15k, Mini-MPL) Rayleigh runs fully self-contained.
- 910 nm cases (CL61, CL31, CL51) need a monthly CAMS file for the mandatory
  water-vapour correction, so they SKIP when CAMS is absent. The HITRAN WV LUT is
  bundled (``calibration/data/abs_cross_wv_910nm.nc``), so no external LUT is needed.

Selected days: clear nights for the molecular instruments, overcast days for the cloud
instruments (see ``examples/data/README.md``).
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo
from calibration.config import InstrumentType, DataLevel
from calibration.cloud import liquid_cloud_calibration, CloudCalConfig

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "examples" / "data"
OPTIONS = REPO / "options.json"

# CAMS lives outside the repo (large, machine-specific). 910 nm tests skip if absent.
CAMS = Path(CalibrationOptions.from_json(OPTIONS).cams_folder) if OPTIONS.exists() else Path("D:/CAMS")
needs_cams = pytest.mark.skipif(
    not CAMS.exists(), reason=f"CAMS folder not available ({CAMS}); needed for 910 nm WV"
)

# (type, wmo, date, InstrumentType, needs_cams) — Rayleigh fixtures
RAYLEIGH = [
    ("CHM15k", "0-20000-0-06610", "20260225", InstrumentType.CHM15k, False),
    ("Mini-MPL", "0-20000-0-07014", "20260423", InstrumentType.MINI_MPL, False),
    ("CL61", "0-756-4-EERLCL61", "20260304", InstrumentType.CL61, True),
]
# (type, wmo, date, expect_valid_constant) — cloud fixtures
CLOUD = [
    ("CL31", "0-20000-0-06602", "20260220", True),
    ("CL51", "0-20000-0-02998", "20260116", True),
    # Sion CL61: cloud profiles are detected but the site rarely yields the warm-liquid
    # stratocumulus O'Connor needs, so we only assert that the cloud path runs on CL61.
    ("CL61", "0-756-4-EERLCL61", "20260414", False),
]


def _l2(wmo: str, date: str) -> Path:
    return DATA / "L2" / wmo / "2026" / date[4:6] / f"L2_{wmo}_A{date}.nc"


def _info(wmo: str, itype: InstrumentType, ncpath: Path) -> InstrumentInfo:
    with Dataset(ncpath) as ds:
        lat = float(ds.variables["station_latitude"][...])
        lon = float(ds.variables["station_longitude"][...])
        alt = float(ds.variables["station_altitude"][...])
    return InstrumentInfo(site_name=wmo, wmo_id=wmo, identifier="A",
                          instrument_type=itype, latitude=lat, longitude=lon, altitude=alt)


# --------------------------------------------------------------------------- #
# Readers: every committed fixture loads and carries a signal variable
# --------------------------------------------------------------------------- #
def test_sample_data_present():
    ncs = list(DATA.rglob("*.nc"))
    assert len(ncs) >= 10, f"expected the bundled fixtures under {DATA}, found {len(ncs)}"


@pytest.mark.parametrize("ncfile", sorted(DATA.rglob("*.nc")), ids=lambda p: p.name)
def test_fixture_readable(ncfile):
    with Dataset(ncfile) as ds:
        assert "time" in ds.dimensions and len(ds.dimensions["time"]) > 0
        assert any(v in ds.variables for v in ("rcs_0", "attenuated_backscatter_0"))
        assert getattr(ds, "instrument_type", "") != ""


# --------------------------------------------------------------------------- #
# Rayleigh (molecular) calibration on the bundled L2 fixtures
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("typ,wmo,date,itype,need_cams", RAYLEIGH, ids=[r[0] for r in RAYLEIGH])
def test_rayleigh_sample(typ, wmo, date, itype, need_cams):
    if need_cams and not CAMS.exists():
        pytest.skip(f"CAMS not available ({CAMS}); needed for {typ} 910 nm WV")
    o = CalibrationOptions.from_json(OPTIONS)
    o.data_level = DataLevel.L2_DAILY
    o.folder_root = DATA / "L2"
    o.cams_folder = CAMS
    o.abs_cs_lookup_table = Path("")          # use the bundled WV LUT
    o.folder_output = Path(tempfile.mkdtemp())
    for k in ("plot_main", "plot_all"):
        if hasattr(o, k):
            setattr(o, k, 0)
    res = calibrate_rayleigh(date, _info(wmo, itype, _l2(wmo, date)), o)
    assert res.flag in (1, 1.0, 0.5), f"{typ} flag={res.flag} ({res.flag_meaning})"
    assert res.lidar_constant > 0


# --------------------------------------------------------------------------- #
# Liquid-cloud calibration on the bundled L2 fixtures (910 nm -> needs CAMS)
# --------------------------------------------------------------------------- #
@needs_cams
@pytest.mark.parametrize("typ,wmo,date,expect_valid", CLOUD, ids=[c[0] for c in CLOUD])
def test_cloud_sample(typ, wmo, date, expect_valid):
    cfg = CloudCalConfig(
        nc_file=str(_l2(wmo, date)),
        cams_folder=str(CAMS),
        apply_wv_correction=1,            # abs_cs_lookup_table omitted -> bundled WV LUT
        aerosol_lidar_ratio=50,
    )
    res = liquid_cloud_calibration(cfg)
    assert res.n_profiles > 0, f"{typ}: no in-cloud profiles detected"
    if expect_valid:
        assert np.isfinite(res.cal_median) and res.cal_median > 0, \
            f"{typ}: cal_median={res.cal_median}"
