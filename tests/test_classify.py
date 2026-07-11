"""ceiloclass integration: the E-PROFILE adapter + the read-once shared-loader path classify
consistently, and the shared loader carries CL61 depolarization.

Skipped when ``ceiloclass`` / ``ceilopyter`` are not installed (they are optional deps, not needed
for calibration) or when the Payerne L1/CAMS sample is absent (CI). Validates that the in-repo
adapter (`calibration/classify/`) reproduces the standalone port and that classifying off the
shared ``InstrumentDayData`` (coarse working grid + depol) matches a direct native file read —
including the depol-ice contamination the Rayleigh scattering-ratio gate misses (the Q5 finding).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ceiloclass")
pytest.importorskip("ceilopyter")

from ceilopyter import average_time  # noqa: E402
from ceiloclass.classification import Target, classify  # noqa: E402

from calibration.classify import (  # noqa: E402
    cams_to_model,
    ceilo_from_shared,
    read_eprofile_l1,
    station_altitude,
)
from calibration.io.instrument_day import load_instrument_day  # noqa: E402

_L1_DIR = Path(r"D:\E-PROFILE_L1_2026\0-20000-0-06610\2026\03")
_CAMS = Path(r"D:\CAMS_Monthly_04\CAMS_Beta_202603.nc")
_CAMS_FB = Path(r"D:\CAMS\CAMS_Beta_202603.nc")
LAT, LON = 46.8137, 6.9425


def _cams() -> Path:
    return _CAMS if _CAMS.is_file() else _CAMS_FB


def _l1(letter: str) -> Path:
    return _L1_DIR / f"L1_0-20000-0-06610_{letter}20260306.nc"


def _fractions(res) -> dict[str, float]:
    tot = res.target.size
    return {t.name: float((res.target == t.value).sum()) / tot * 100 for t in Target}


_HAVE_CAMS = _cams().is_file()
_REPO = Path(__file__).resolve().parents[1]
_L1_CHM = _L1_DIR / "L1_0-20000-0-06610_A20260304.nc"


@pytest.mark.skipif(not (_L1_CHM.is_file() and _HAVE_CAMS), reason="Payerne L1/CAMS sample absent (CI)")
def test_rayleigh_classification_screen_flags_contaminated_window():
    """Q5 Rayleigh screen (calibrate_rayleigh ``contam_profile``): a successful night is rejected with
    flag -11 when the contamination profile covers its selected molecular window, and is UNCHANGED when
    the contamination is elsewhere or absent (so it never touches a clean night)."""
    import tempfile

    from calibration import CalibrationOptions, InstrumentInfo, calibrate_rayleigh
    from calibration.config import DataLevel, InstrumentType

    with tempfile.TemporaryDirectory() as td:
        def opts():
            o = CalibrationOptions.from_json(str(_REPO / "options.json"))
            o.folder_root = _L1_DIR.parents[2]     # D:/E-PROFILE_L1_2026
            o.data_level = DataLevel.L1
            o.cams_folder = _cams().parent
            o.plot_main = o.plot_all = False
            o.folder_output = Path(td)
            return o

        info = InstrumentInfo(site_name="PAY", wmo_id="0-20000-0-06610", identifier="A",
                              instrument_type=InstrumentType("CHM15k"),
                              latitude=LAT, longitude=LON, altitude=490.0)
        base = calibrate_rayleigh("20260304", info, opts())
        if base.flag not in (1, 1.0, 0.5) or not np.isfinite(base.calibration_bottom_height):
            pytest.skip("baseline night is not a Rayleigh success on this archive")

        a, b = base.calibration_bottom_height, base.calibration_top_height
        h = np.linspace(0.0, 8000.0, 200)
        over = np.column_stack([h, ((h >= a - 50) & (h <= b + 50)).astype(float)])       # 100% in window
        outside = np.column_stack([h, (h < max(a - 800.0, 500.0)).astype(float)])        # 100% below it

        assert calibrate_rayleigh("20260304", info, opts(), contam_profile=over).flag == -11
        assert calibrate_rayleigh("20260304", info, opts(), contam_profile=outside).flag == base.flag
        assert calibrate_rayleigh("20260304", info, opts(), contam_profile=None).flag == base.flag


@pytest.mark.skipif(not (_l1("C").is_file() and _HAVE_CAMS), reason="Payerne CL61 L1 / CAMS sample absent (CI)")
def test_cl61_shared_loader_matches_file_and_detects_contamination():
    """CL61 Payerne 2026-03-06: the read-once shared loader (coarse working grid + depol) gives the
    same classification as a direct native file read, and both flag the depol-ice contamination."""
    cams = str(_cams())
    f = str(_l1("C"))
    alt = station_altitude(f)

    ceilo_f = average_time(read_eprofile_l1(f), 30.0)
    res_f = classify(
        ceilo_f, cams_to_model(cams, LAT, LON, ceilo_f.time, ceilo_f.range, alt),
        altitude=alt, use_wet_bulb=False,
    )

    idd = load_instrument_day([f], "CL61", cams, read_cams=False, target_time_s=30.0, target_range_m=10.0)
    ceilo_s = ceilo_from_shared(idd)
    assert ceilo_s.depol is not None, "CL61 depolarization must survive the shared loader"
    res_s = classify(
        ceilo_s, cams_to_model(cams, LAT, LON, ceilo_s.time, ceilo_s.range, alt),
        altitude=alt, use_wet_bulb=False,
    )

    ff, fs = _fractions(res_f), _fractions(res_s)
    # The depol-ice layer (2.5-3.7 km, 0.2 depol) is detected on both paths (Q5).
    assert ff["ICE"] > 0.5 and fs["ICE"] > 0.5
    # File (native range) vs shared (coarse 30 s/10 m) agree closely on every dominant class.
    for t in ("CLEAR", "AEROSOL", "ICE"):
        assert abs(ff[t] - fs[t]) < 1.0, f"{t}: file {ff[t]:.2f}% vs shared {fs[t]:.2f}%"


@pytest.mark.skipif(not (_l1("A").is_file() and _HAVE_CAMS), reason="Payerne CHM15k L1 / CAMS sample absent (CI)")
def test_single_channel_adapter_classifies_via_shared_loader():
    """A single-channel instrument (CHM15k) has no depol yet still classifies through the shared
    loader path -- guards that the adapter handles the depol-absent case."""
    cams = str(_cams())
    f = str(_l1("A"))
    alt = station_altitude(f)
    idd = load_instrument_day([f], "CHM15k", cams, read_cams=False, target_time_s=30.0, target_range_m=10.0)
    ceilo = ceilo_from_shared(idd)
    assert ceilo.depol is None
    res = classify(ceilo, cams_to_model(cams, LAT, LON, ceilo.time, ceilo.range, alt),
                   altitude=alt, use_wet_bulb=False)
    assert res.target.size > 0 and np.isfinite(res.strong_beta)
