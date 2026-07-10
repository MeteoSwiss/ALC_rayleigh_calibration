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
