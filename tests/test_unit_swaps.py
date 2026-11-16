"""Unit swaps inside a stream and per-stream Kalman options (v2.3 release, known limit 1).

* the Vaisala type of one L1 file comes from its range extent (CL31 7.7 km, CL51 15.4 km), every other
  type passes through, an unreadable file falls back to the census type;
* the Kalman starts afresh on every instrument_events.json entry with kalman_restart, so no value of
  the old unit leaks into the new unit's best estimate;
* station_kalman_params.exclude_flags keeps those nights out of one stream's Kalman only.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from calibration.config import CalibrationOptions
from calibration.io import vaisala_type as vt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_network_calibration as rnc  # noqa: E402


def _l1(path: Path, n_gates: int, step_m: float) -> Path:
    with Dataset(path, "w") as nc:
        nc.createDimension("range", n_gates)
        v = nc.createVariable("range", "f4", ("range",))
        v[:] = (np.arange(n_gates) + 1) * step_m
    return path


@pytest.mark.parametrize("n,step,census,expected", [
    (770, 10.0, "CL31", "CL31"),      # CL31, 10 m gates
    (385, 20.0, "CL31", "CL31"),      # CL31, 20 m gates
    (1540, 10.0, "CL31", "CL51"),     # a CL51 still labelled CL31 (Usti, Tusimice, Libus, Plzen)
    (1540, 10.0, "CL51", "CL51"),
    (770, 10.0, "CL51", "CL31"),      # the CL31 part of a stream the census already calls CL51
    (1540, 10.0, "CHM15k", "CHM15k"),  # non-Vaisala types pass through untouched
    (770, 10.0, "CL61", "CL61"),
])
def test_type_from_range_extent(tmp_path, n, step, census, expected):
    f = _l1(tmp_path / f"L1_{n}_{step}.nc", n, step)
    assert vt.file_instrument_type(f, census) == expected


def test_unreadable_file_keeps_census_type(tmp_path):
    assert vt.file_instrument_type(tmp_path / "missing.nc", "CL31") == "CL31"
    bad = tmp_path / "bad.nc"
    bad.write_bytes(b"not a netcdf")
    assert vt.file_instrument_type(bad, "CL51") == "CL51"


def _rows(start, values, flag=1.0, method="cloud"):
    d0 = datetime.strptime(start, "%Y%m%d")
    return [dict(date=(d0 + timedelta(days=i)).strftime("%Y%m%d"), method=method, flag=flag,
                 cal_value=v, uncertainty=0.05 * v) for i, v in enumerate(values)]


def _events(tmp_path, monkeypatch, key, date):
    f = tmp_path / "events.json"
    f.write_text(json.dumps({key: [{"date": date, "label": "swap", "kalman_restart": True}]}))
    monkeypatch.setattr(rnc, "EVENTS_FILE", f)


def test_kalman_restarts_at_a_unit_change(tmp_path, monkeypatch):
    key = "0-99999-0-SWAP_A"
    rows = _rows("20260601", [1.0e8] * 30) + _rows("20260701", [2.8e8] * 30)
    _events(tmp_path, monkeypatch, key, "2026-07-01")
    plain = {r["date"]: float(r["kalman"]) for r in rnc._kalman_rows(rows, weighted=False)}
    split = {r["date"]: float(r["kalman"]) for r in rnc._kalman_rows(rows, weighted=False, key=key)}
    # Without the restart the filter drags the old unit into July; with it, July starts on 2.8e8.
    assert plain["20260702"] < 2.0e8
    assert split["20260702"] == pytest.approx(2.8e8, rel=1e-6)
    assert split["20260629"] == pytest.approx(1.0e8, rel=1e-6)
    # No grid day between the last old-unit night and the first new-unit night is invented.
    assert all(not ("20260630" < d < "20260701") for d in split)


def test_restart_needs_kalman_restart_true(tmp_path, monkeypatch):
    key = "0-99999-0-MARK_A"
    f = tmp_path / "events.json"
    f.write_text(json.dumps({key: [{"date": "2026-07-01", "label": "window cleaned"}]}))
    monkeypatch.setattr(rnc, "EVENTS_FILE", f)
    assert rnc._kalman_restarts(key) == []
    rows = _rows("20260601", [1.0e8] * 30) + _rows("20260701", [2.8e8] * 30)
    assert rnc._kalman_rows(rows, weighted=False, key=key) == rnc._kalman_rows(rows, weighted=False)


def test_repo_events_restart_the_four_vaisala_swaps_and_payerne_b():
    keys = ["0-20000-0-11450_A", "0-20000-0-11438_A", "0-20000-0-11520_A", "0-20000-0-11502_A",
            "0-20000-0-06610_B"]
    for k in keys:
        assert len(rnc._kalman_restarts(k)) == 1, k


def test_excluded_flags_only_for_the_listed_stream(tmp_path, monkeypatch):
    rows = (_rows("20260501", [1.0e11] * 20, flag=1.0, method="rayleigh")
            + _rows("20260521", [2.6e11] * 20, flag=0.25, method="rayleigh"))
    opts = CalibrationOptions.from_json(REPO / "options.json")
    opts.station_kalman_params = {"0-99999-0-OSLO_A": {"exclude_flags": [0.25]}}
    monkeypatch.setattr(rnc.CalibrationOptions, "from_json", classmethod(lambda cls, p: opts))
    other = {r["date"]: float(r["kalman"]) for r in rnc._kalman_rows(rows, weighted=False,
                                                                       key="0-99999-0-OTHER_A")}
    oslo = {r["date"]: float(r["kalman"]) for r in rnc._kalman_rows(rows, weighted=False,
                                                                      key="0-99999-0-OSLO_A")}
    assert max(other.values()) > 1.5e11           # low-sun nights pull an unlisted stream up
    assert max(oslo.values()) == pytest.approx(1.0e11, rel=1e-6)
    assert max(oslo) == "20260520"                 # the excluded nights are not in its grid


def test_repo_options_exclude_oslo_low_sun_nights():
    o = CalibrationOptions.from_json(REPO / "options.json")
    assert o.station_kalman_params["0-20000-0-01492_A"]["exclude_flags"] == [0.25]
