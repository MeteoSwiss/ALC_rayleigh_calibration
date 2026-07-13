"""Chunk-outer read-once (doc/reports/11_batch_readonce_design.md): unit + parity tests.

The parity test needs the local Payerne L1 archive (D:/E-PROFILE_L1_2026) and is skipped
where that data is absent (CI): it asserts that a night served as ``slice_night(D)`` of a
multi-day batched read is identical to today's fresh per-night ``[D-1, D]`` read, for the
native grid, the re-coarsened day slice and the OmB/sens day dict.
"""
import datetime as dt
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from calibration.io.instrument_day import load_instrument_day  # noqa: E402
import run_network_calibration as R  # noqa: E402

_L1_DIR = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610/2026/06")
_HAVE_L1 = all((_L1_DIR / f"L1_0-20000-0-06610_A202606{d:02d}.nc").exists()
               for d in (13, 14, 15, 16))


def test_chunks_split_covers_window_once():
    s, e = datetime(2026, 5, 1), datetime(2026, 6, 30)
    ch = list(R._chunks(s, e, 31))
    assert ch[0] == (datetime(2026, 5, 1), datetime(2026, 5, 31))
    assert ch[1] == (datetime(2026, 6, 1), datetime(2026, 6, 30))
    days = [d for cs, ce in ch for d in R._days(cs, ce)]
    assert len(days) == 61 and len(set(days)) == 61          # every day exactly once
    assert list(R._chunks(s, s, 7)) == [(s, s)]               # 1-day window = 1 chunk (daily run)
    # chunk sizes never exceed ndays and concatenate seamlessly
    for (a0, a1), (b0, b1) in zip(ch, ch[1:]):
        assert (a1 - a0).days < 31 or (a1 - a0).days == 30
        assert b0 == a1 + dt.timedelta(days=1)


@pytest.mark.skipif(not _HAVE_L1, reason="local Payerne L1 archive not available")
def test_slice_night_matches_per_night_read():
    f = {d: str(_L1_DIR / f"L1_0-20000-0-06610_A202606{d:02d}.nc") for d in (13, 14, 15, 16)}
    kw = dict(cams_folder="D:/CAMS_Monthly_04", read_cams=False, build_working=False)
    chunk = load_instrument_day([f[13], f[14], f[15], f[16]], "CHM15k", **kw)
    night = load_instrument_day([f[14], f[15]], "CHM15k", **kw)     # night of the 15th = [14, 15]
    assert chunk is not None and night is not None

    view = chunk.slice_night(dt.date(2026, 6, 15))
    # native parity: same profiles, same signal, same cloud bases
    assert np.array_equal(np.asarray(view.native.time), np.asarray(night.native.time))
    assert np.array_equal(np.asarray(view.native.rcs, dtype="float64"),
                          np.asarray(night.native.rcs, dtype="float64"), equal_nan=True)
    assert np.array_equal(np.asarray(view.native.cbh, dtype="float64"),
                          np.asarray(night.native.cbh, dtype="float64"), equal_nan=True)
    # `working` mirrors the per-night read (build_working=False -> alias of native)
    assert view.working is view.native and night.working is night.native

    # day slice built from the night view == day slice from the real night read (coarse grid)
    v_day = view.slice_to_date(dt.date(2026, 6, 15))
    n_day = night.slice_to_date(dt.date(2026, 6, 15))
    assert np.array_equal(np.asarray(v_day.working.time), np.asarray(n_day.working.time))
    assert np.array_equal(np.asarray(v_day.working.rcs, dtype="float64"),
                          np.asarray(n_day.working.rcs, dtype="float64"), equal_nan=True)

    # OmB/sens day dict parity
    a, b = v_day.to_omb_dict(), n_day.to_omb_dict()
    assert np.array_equal(a["time"], b["time"])
    assert np.array_equal(a["rcs"], b["rcs"], equal_nan=True)
    assert np.array_equal(a["cbh"], b["cbh"], equal_nan=True)
    assert a["wl"] == b["wl"] and a["alt"] == b["alt"]


@pytest.mark.skipif(not _HAVE_L1, reason="local Payerne L1 archive not available")
def test_chunk_reader_serves_nights_and_boundaries():
    s = dict(wmo="0-20000-0-06610", ident="A", type="CHM15k",
             lat=46.81, lon=6.94, alt=490.0, site="PAYERNE")
    start, end = datetime(2026, 6, 14), datetime(2026, 6, 16)
    reader = R._make_chunk_reader(s, start, end, chunk_days=2)    # boundary between 15 and 16
    for ds, d in (("20260615", dt.date(2026, 6, 15)), ("20260616", dt.date(2026, 6, 16))):
        got = reader(ds)
        assert got is not None
        ref = load_instrument_day(
            [str(_L1_DIR / f"L1_0-20000-0-06610_A{(d - dt.timedelta(days=1)):%Y%m%d}.nc"),
             str(_L1_DIR / f"L1_0-20000-0-06610_A{d:%Y%m%d}.nc")],
            "CHM15k", cams_folder="D:/CAMS_Monthly_04", read_cams=False, build_working=False)
        # 20260616 is the FIRST day of the second chunk: its night must still include the 15th
        assert np.array_equal(np.asarray(got.native.time), np.asarray(ref.native.time))
        assert np.array_equal(np.asarray(got.native.rcs, dtype="float64"),
                              np.asarray(ref.native.rcs, dtype="float64"), equal_nan=True)
