"""Tests for CAMS file resolution and the download script's calibration compatibility.

The pure-Python parts run fully offline (no network, no eccodes, no archive):
  * ``download_cams_beta.parse_args`` (month / day / range),
  * ``build_output`` writes ``z``/``lnsp`` as SURFACE fields (the layout both
    calibration readers require — the historical bug overwrote ``z`` with the full
    model-level geopotential, which moved the surface ~24 km),
  * the monthly-or-daily resolver (``find_cams_file`` / ``ensure_cams_file``).

One end-to-end compatibility test reconstructs a download-style file from the real
monthly archive and confirms the calibration reader returns an identical T/p/H
profile; it skips when ``D:/CAMS`` is unavailable (e.g. CI).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from calibration.io import download_cams_beta as dl
from calibration.io.cams import (
    candidate_cams_files,
    ensure_cams_file,
    find_cams_file,
    _month_dates,
    _night_dates,
)

# The CAMS_Beta model-level set: model level 1 (surface-field slot) + 38..137.
LEVELS = np.array([1] + list(range(38, 138)), dtype="int32")


# --------------------------------------------------------------------------- #
# parse_args
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "argv,label,ndays",
    [
        (["202606"], "202606", 30),
        (["20260615"], "20260615", 1),
        (["20260614", "20260615"], "20260614_20260615", 2),
    ],
)
def test_parse_args(argv, label, ndays):
    dates, lab = dl.parse_args(argv)
    assert lab == label
    assert len(dates) == ndays
    assert dates[0].count("-") == 2          # YYYY-MM-DD form


def test_parse_args_swaps_reversed_range():
    dates, _ = dl.parse_args(["20260615", "20260614"])
    assert dates[0] == "2026-06-14" and dates[-1] == "2026-06-15"


# --------------------------------------------------------------------------- #
# build_output: z / lnsp must be SURFACE fields (finite only at level index 0)
# --------------------------------------------------------------------------- #
def _synthetic_raw(nt: int = 1, ny: int = 2, nx: int = 2):
    """A minimal grib_to_netcdf-style raw dataset: t/q full profiles, z/lnsp
    broadcast surface fields (as the ADS delivers them)."""
    import xarray as xr

    nl = LEVELS.size
    shape = (nt, nl, ny, nx)
    dims = ("time", "level", "latitude", "longitude")
    rng = np.random.default_rng(0)

    t = np.broadcast_to(np.linspace(280.0, 220.0, nl)[None, :, None, None], shape).copy()
    q = np.full(shape, 1e-3)
    z_surf = 5000.0          # m^2 s^-2  (~510 m ASL)
    lnsp = float(np.log(9.0e4))   # ~900 hPa surface pressure

    data = {
        "t": (dims, t),
        "q": (dims, q),
        "z": (dims, np.full(shape, z_surf)),       # broadcast across all levels
        "lnsp": (dims, np.full(shape, lnsp)),      # broadcast across all levels
    }
    for v in ("aerext355", "aerext532", "aerext1064",
              "aerbackscatgnd355", "aerbackscatgnd532", "aerbackscatgnd1064"):
        data[v] = (dims, rng.random(shape) * 1e-6)

    times = np.datetime64("2026-01-01T03:00") + np.arange(nt) * np.timedelta64(3, "h")
    coords = {
        "time": ("time", times),       # datetime64, like real grib_to_netcdf output
        "level": ("level", LEVELS),
        "latitude": ("latitude", np.linspace(46.5, 46.0, ny)),
        "longitude": ("longitude", np.linspace(7.0, 7.5, nx)),
    }
    return xr.Dataset(data, coords=coords), z_surf, lnsp


def test_build_output_writes_surface_z_lnsp(tmp_path):
    import xarray as xr

    raw, z_surf, lnsp = _synthetic_raw()
    out = tmp_path / "CAMS_Beta_20260101.nc"
    dl.build_output(raw, str(out))

    with xr.open_dataset(out) as o:
        for v in ("t", "q", "z", "lnsp", "altitude", "altitude_agl"):
            assert v in o, f"{v} missing from download output"
        z = o["z"].isel(time=0, latitude=0, longitude=0).values
        s = o["lnsp"].isel(time=0, latitude=0, longitude=0).values

    # The whole point of the fix: z and lnsp are finite ONLY at level index 0
    # (the surface slot), NaN on every other level — exactly how the archive stores
    # them and how _cams_levels / _cams_levels_all_times read them.
    assert np.isfinite(z[0]) and np.isfinite(s[0])
    assert not np.any(np.isfinite(z[1:])), "z must be NaN above the surface slot"
    assert not np.any(np.isfinite(s[1:])), "lnsp must be NaN above the surface slot"
    assert abs(float(z[0]) - z_surf) < 1.0      # surface geopotential preserved
    assert abs(float(s[0]) - lnsp) < 1e-3


# --------------------------------------------------------------------------- #
# Resolver: monthly OR daily
# --------------------------------------------------------------------------- #
def test_candidate_order_and_names():
    cands = candidate_cams_files("/cams", "20260615")
    assert [p.name for p in cands] == ["CAMS_Beta_202606.nc", "CAMS_Beta_20260615.nc"]


def test_find_prefers_monthly_then_daily(tmp_path):
    assert find_cams_file(tmp_path, "20260615") is None
    daily = tmp_path / "CAMS_Beta_20260615.nc"
    daily.write_bytes(b"x")
    assert find_cams_file(tmp_path, "20260615") == daily     # daily accepted
    monthly = tmp_path / "CAMS_Beta_202606.nc"
    monthly.write_bytes(b"x")
    assert find_cams_file(tmp_path, "20260615") == monthly   # monthly preferred


def test_ensure_no_autodownload_returns_none(tmp_path):
    assert ensure_cams_file(tmp_path, "20991231", auto_download=False) is None


def test_grib_to_netcdf_prefers_cli_then_cfgrib(monkeypatch):
    """Conversion uses the eccodes CLI when present, else the pure-Python cfgrib path
    (no conda / no system eccodes needed)."""
    calls = []

    # CLI on PATH -> shell out to grib_to_netcdf
    monkeypatch.setattr(dl.shutil, "which", lambda name: "grib_to_netcdf")
    monkeypatch.setattr(dl.subprocess, "run", lambda *a, **k: calls.append("cli"))
    dl.grib_to_netcdf("x.grib", "y.nc")
    assert calls == ["cli"]

    # CLI absent -> fall back to cfgrib (_cfgrib_to_dataset)
    class _FakeDS:
        def to_netcdf(self, p):
            calls.append("cfgrib")

        def close(self):
            pass

    monkeypatch.setattr(dl.shutil, "which", lambda name: None)
    monkeypatch.setattr(dl, "_cfgrib_to_dataset", lambda p: _FakeDS())
    dl.grib_to_netcdf("x.grib", "y.nc")
    assert calls == ["cli", "cfgrib"]


def test_ensure_autodownload_failure_is_graceful(tmp_path, monkeypatch):
    """A failed download (no eccodes, ADS error, ...) must not crash the calibration:
    ensure_cams_file returns None, removes any partial file, and swallows the error so
    the night is simply skipped (flag -4)."""
    def boom(dates, out_path, **kw):
        Path(out_path).write_bytes(b"partial")      # simulate a partial write
        raise RuntimeError("simulated ADS/eccodes failure")

    monkeypatch.setattr(dl, "download_to_netcdf", boom)
    res = ensure_cams_file(tmp_path, "20991231", auto_download=True, scope="day")
    assert res is None
    assert list(tmp_path.iterdir()) == []           # partial file cleaned up


def test_night_and_month_dates():
    assert _night_dates("20260615") == ["2026-06-14", "2026-06-15"]
    md = _month_dates("202602")
    assert md[0] == "2026-02-01" and md[-1] == "2026-02-28" and len(md) == 28


# --------------------------------------------------------------------------- #
# End-to-end: a download-style file is read IDENTICALLY to the source archive
# --------------------------------------------------------------------------- #
_CAMS = Path("D:/CAMS")
_SRC = _CAMS / "CAMS_Beta_202605.nc"


@pytest.mark.skipif(not _SRC.exists(), reason=f"real CAMS archive not available ({_SRC})")
def test_download_style_file_reads_identically(tmp_path):
    import xarray as xr

    from calibration.water_vapor_correction.water_vapor import cams_temperature_pressure_profile

    lat, lon = 46.8, 6.9
    t0 = np.datetime64("2026-05-14T20:00")
    t1 = np.datetime64("2026-05-15T04:00")
    h_ref, t_ref, p_ref = cams_temperature_pressure_profile(_SRC, lat, lon, t0, t1)

    # Reconstruct a download-style file (build_output on a monthly slice with the
    # surface z/lnsp broadcast across levels, as the raw ADS GRIB delivers them).
    ds = xr.open_dataset(_SRC)
    sub = ds.sel(latitude=slice(47.5, 46.0), longitude=slice(6.0, 7.5))
    m = (ds.time.values >= np.datetime64("2026-05-14T00")) & (ds.time.values <= np.datetime64("2026-05-15T12"))
    sub = sub.isel(time=np.where(m)[0])
    z0 = sub["z"].isel(level=0).values
    s0 = sub["lnsp"].isel(level=0).values
    raw = sub[["aerext355", "aerext532", "aerext1064",
               "aerbackscatgnd355", "aerbackscatgnd532", "aerbackscatgnd1064", "t", "q"]].copy()
    raw["z"] = sub["t"].dims, np.broadcast_to(z0[:, None, :, :], sub["t"].shape).copy()
    raw["lnsp"] = sub["t"].dims, np.broadcast_to(s0[:, None, :, :], sub["t"].shape).copy()
    ds.close()

    out = tmp_path / "CAMS_Beta_20260515.nc"
    dl.build_output(raw, str(out))

    h_d, t_d, p_d = cams_temperature_pressure_profile(out, lat, lon, t0, t1)
    assert np.nanmax(np.abs(h_d - h_ref)) < 0.5                       # sub-meter (float32)
    assert np.nanmax(np.abs(t_d - t_ref)) < 0.01                      # K
    assert np.nanmax(np.abs(p_d - p_ref)) / np.nanmax(p_ref) < 1e-4   # relative
