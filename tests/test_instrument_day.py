"""Tests for the read-once shared loader (``calibration.io.instrument_day``).

Phase B foundation (B1 + B2): one L1 read + one CAMS closest-cell read + one WV
transmission per instrument-day, with a dynamically bin-averaged working grid.

Two layers:
  * fast UNIT tests of the per-file coarsening rule + the read_cams=False path on the
    bundled CL61 L1 — always run;
  * an INTEGRATION test that the once-computed WV transmission is bit-identical to the
    production ``compute_wv_transmission`` — skipped when the CAMS file is absent (CI).
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from calibration.io.instrument_day import _applied_factors, load_instrument_day

REPO = Path(__file__).resolve().parents[1]
BUNDLED_L1 = REPO / "examples/data/L1/0-756-4-EERLCL61/2026/03/L1_0-756-4-EERLCL61_A20260304.nc"
CAMS_DIR = Path("D:/CAMS_Monthly_04")
CAMS_FALLBACK = Path("D:/CAMS")


def _cams_available(yyyymm: str) -> bool:
    return (CAMS_DIR / f"CAMS_Beta_{yyyymm}.nc").is_file() or (
        CAMS_FALLBACK / f"CAMS_Beta_{yyyymm}.nc"
    ).is_file()


# --------------------------------------------------------------------------- unit


def test_coarsen_factor_is_dynamic_per_file():
    """The block factor is round(target / median native step), clamped >=1 — derived from
    the file's own grid, never hard-coded per instrument type."""
    # CL61-like: 4.8 m / 30 s -> range averaged x2 (9.6 m), time untouched.
    fine = SimpleNamespace(
        time_datetime=[dt.datetime(2026, 3, 6) + dt.timedelta(seconds=30 * i) for i in range(50)],
        range_alc=np.arange(3276) * 4.8,
    )
    assert _applied_factors(fine, target_range_m=10.0, target_time_s=15.0) == (1, 2)

    # CHM15k-like: 15 m / 15 s -> already at/above target, nothing averaged.
    coarse = SimpleNamespace(
        time_datetime=[dt.datetime(2026, 3, 6) + dt.timedelta(seconds=15 * i) for i in range(50)],
        range_alc=np.arange(1024) * 15.0,
    )
    assert _applied_factors(coarse, target_range_m=10.0, target_time_s=15.0) == (1, 1)

    # A high-rate stream: 5 s -> time averaged x3.
    fast = SimpleNamespace(
        time_datetime=[dt.datetime(2026, 3, 6) + dt.timedelta(seconds=5 * i) for i in range(50)],
        range_alc=np.arange(770) * 10.0,
    )
    assert _applied_factors(fast, target_range_m=10.0, target_time_s=15.0) == (3, 1)


@pytest.mark.skipif(not BUNDLED_L1.is_file(), reason="bundled CL61 L1 not present")
def test_load_once_and_coarsen_no_cams():
    """read_cams=False: L1 is read once and the working grid is a consistent bin-average."""
    idd = load_instrument_day([BUNDLED_L1], "CL61", CAMS_DIR, read_cams=False)
    assert idd is not None
    tf, rf = idd.coarsen
    assert tf >= 1 and rf >= 1
    # Working grid is coarser-or-equal, and the shapes match the reported factors.
    assert idd.working.rcs.shape[0] == idd.native.rcs.shape[0] // tf
    assert idd.working.rcs.shape[1] == idd.native.rcs.shape[1] // rf
    assert idd.working.range_alc.size <= idd.native.range_alc.size
    # No CAMS requested -> no cell, no WV.
    assert idd.cams_temperature is None
    assert idd.wv_transmission is None


# -------------------------------------------------------------------- integration


@pytest.mark.skipif(
    not (BUNDLED_L1.is_file() and _cams_available("202603")),
    reason="bundled L1 or CAMS 202603 absent (e.g. CI)",
)
def test_wv_transmission_once_matches_production():
    """The WV transmission computed once in the loader is bit-identical to a direct
    production ``compute_wv_transmission`` call on the same working grid."""
    from calibration.cloud.calibration import (
        CloudCalConfig,
        compute_wv_transmission,
        set_defaults,
    )

    idd = load_instrument_day(
        [BUNDLED_L1], "CL61", CAMS_DIR, cams_folder_fallback=CAMS_FALLBACK
    )
    assert idd is not None and idd.wv_transmission is not None
    wv = idd.wv_transmission
    n_time, n_range = idd.working.rcs.shape
    assert wv.shape == (n_range, n_time)
    assert np.isfinite(wv).all()
    assert (wv > 0).all() and (wv <= 1.0 + 1e-9).all()

    w = idd.working
    view = SimpleNamespace(
        range=np.asarray(w.range_alc, float),
        time=np.asarray(w.time_datetime, dtype="datetime64[ns]"),
        time_num=np.asarray(w.time, float) + 719529.0,
        beta=np.asarray(w.rcs, float).T,
        station_latitude=w.latitude,
        station_longitude=w.longitude,
        station_altitude=w.altitude,
    )
    cfg = set_defaults(
        CloudCalConfig(
            instrument="CL61",
            apply_wv_correction=True,
            cams_folder=str(CAMS_DIR),
            cams_folder_fallback=str(CAMS_FALLBACK),
            station_latitude=w.latitude,
            station_longitude=w.longitude,
        )
    )
    ref = compute_wv_transmission(view, cfg)
    assert np.array_equal(ref, wv)


# ------------------------------------------------ cloud: _ceilo_from_shared


@pytest.mark.skipif(not BUNDLED_L1.is_file(), reason="bundled CL61 L1 not present")
def test_ceilo_from_shared_beta_matches_reader():
    """_ceilo_from_shared rebuilds cloud's beta from the shared read, bit-identical to
    read_ceilometer_data (bundled CL61 -> physical units, factor 1, no /C)."""
    from calibration.cloud.calibration import (
        CloudCalConfig,
        _ceilo_from_shared,
        read_ceilometer_data,
        set_defaults,
    )

    ref, status = read_ceilometer_data(
        str(BUNDLED_L1), set_defaults(CloudCalConfig(instrument="CL61"))
    )
    assert status == 0
    idd = load_instrument_day([BUNDLED_L1], "CL61", CAMS_DIR, read_cams=False)
    mine = _ceilo_from_shared(idd, set_defaults(CloudCalConfig(instrument="CL61")))

    assert mine.beta.shape == ref.beta.shape  # (range, time)
    fa, fb = np.isfinite(ref.beta), np.isfinite(mine.beta)
    assert np.array_equal(fa, fb)
    assert np.array_equal(ref.beta[fa], mine.beta[fb])
    assert np.array_equal(np.asarray(ref.range), np.asarray(mine.range))


@pytest.mark.skipif(not BUNDLED_L1.is_file(), reason="bundled CL61 L1 not present")
def test_slice_to_date_matches_full_day():
    """slice_to_date keeps a single UTC day's profiles, and the cloud beta rebuilt from the
    slice still matches read_ceilometer_data -- so a night [D-1,D] read can feed the day-D
    cloud pass."""
    from calibration.cloud.calibration import (
        CloudCalConfig,
        _ceilo_from_shared,
        read_ceilometer_data,
        set_defaults,
    )

    idd = load_instrument_day([BUNDLED_L1], "CL61", CAMS_DIR, read_cams=False)
    # UTC date from time (days since 1970); time_datetime may be cftime, not datetime.
    day0 = int(np.floor(float(np.asarray(idd.native.time)[0])))
    date = dt.date(1970, 1, 1) + dt.timedelta(days=day0)
    sliced = idd.slice_to_date(date)
    assert sliced.native.rcs.shape == idd.native.rcs.shape  # all profiles kept

    ref, _ = read_ceilometer_data(
        str(BUNDLED_L1), set_defaults(CloudCalConfig(instrument="CL61"))
    )
    mine = _ceilo_from_shared(sliced, set_defaults(CloudCalConfig(instrument="CL61")))
    fa, fb = np.isfinite(ref.beta), np.isfinite(mine.beta)
    assert np.array_equal(fa, fb)
    assert np.array_equal(ref.beta[fa], mine.beta[fb])


_PAYERNE_MAR = Path(r"D:\E-PROFILE_L1_2026\0-20000-0-06610\2026\03")


@pytest.mark.skipif(not _PAYERNE_MAR.is_dir(), reason="L1 archive not present (e.g. CI)")
@pytest.mark.parametrize("label,letter", [("CHM15k", "A"), ("CL31", "B")])
def test_ceilo_from_shared_beta_raw_path(label, letter):
    """The raw-signal path: CHM15k (counts) / CL31 (V*m^2) divide by the calibration
    constant -- still bit-identical to read_ceilometer_data."""
    from calibration.cloud.calibration import (
        CloudCalConfig,
        _ceilo_from_shared,
        read_ceilometer_data,
        set_defaults,
    )

    f = _PAYERNE_MAR / f"L1_0-20000-0-06610_{letter}20260306.nc"
    if not f.is_file():
        pytest.skip(f"{f.name} absent")
    ref, status = read_ceilometer_data(str(f), set_defaults(CloudCalConfig(instrument=label)))
    assert status == 0
    idd = load_instrument_day([f], label, CAMS_DIR, read_cams=False)
    mine = _ceilo_from_shared(idd, set_defaults(CloudCalConfig(instrument=label)))
    fa, fb = np.isfinite(ref.beta), np.isfinite(mine.beta)
    assert np.array_equal(fa, fb)
    assert np.array_equal(ref.beta[fa], mine.beta[fb])
