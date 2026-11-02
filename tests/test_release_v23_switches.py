# -*- coding: utf-8 -*-
"""Release v2.3 switches (2026-11): eta tables v2 and their generations, version codes, low-sun version,
per-station molecular_params override, CL61 tier-2 window cap, polar env override."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]


def test_eta_tables_v2_match_the_generator():
    import calibration.cloud._filters as CF
    assert (CF.ETA_V2_AG_UM, CF.ETA_V2_ALPHA_KM) == (3.3, 20.0)
    from validation.multiple_scattering_eta import CBH_GRID_KM, eta_table
    for it in ("CL31", "CL51", "CL61"):
        tab = getattr(CF, f"_ETA_{it}")
        ref = eta_table(it, CF.ETA_V2_AG_UM, alpha=CF.ETA_V2_ALPHA_KM * 1e-3)[:13]
        assert np.allclose(tab[:, 0], CBH_GRID_KM[:13])
        assert np.allclose(tab[:, 1], ref, atol=6e-6)
        old = getattr(CF, f"_ETA_PVC55_{it}")
        assert np.all(tab[:10, 1] < old[:, 1])       # stronger correction everywhere (common nodes)


def test_eta_generation_switch(monkeypatch):
    import calibration.cloud._filters as CF
    from calibration.cloud import CloudCalConfig
    rng = np.array([250.0, 875.0, 2375.0])
    beta = np.ones((3, 1))
    cfg = CloudCalConfig(instrument="CL31")
    monkeypatch.delenv("ALC_ETA_TABLES", raising=False)
    assert CF.eta_table_generation() == "v2" and CF.cloud_version_name() == "cloud_oconnor_v2"
    v2 = CF.apply_multiple_scattering_correction(beta, rng, cfg)[:, 0]
    monkeypatch.setenv("ALC_ETA_TABLES", "pvc55")
    assert CF.eta_table_generation() == "pvc55" and CF.cloud_version_name() == "cloud_oconnor"
    old = CF.apply_multiple_scattering_correction(beta, rng, cfg)[:, 0]
    assert old[0] == pytest.approx(0.93354, abs=1e-5) and old[2] == pytest.approx(0.73795, abs=1e-5)
    assert np.all(v2 < old)
    monkeypatch.setenv("ALC_ETA_TABLES", "hopkin_legacy")
    leg = CF.apply_multiple_scattering_correction(beta, rng, cfg)[:, 0]
    assert leg[0] == pytest.approx(0.82854, abs=1e-5)


def test_version_codes_and_lowsun():
    from calibration.io.output import VERSION_CODES, version_code, version_code_for
    assert VERSION_CODES["eprof_v2.3"] == 230 and VERSION_CODES["eprof_v2.3_lowsun"] == 231
    assert VERSION_CODES["eprof_v2.2_lowsun"] == 221
    assert VERSION_CODES["cloud_oconnor"] == 1000 and VERSION_CODES["cloud_oconnor_v2"] == 1100
    assert len(set(VERSION_CODES.values())) == len(VERSION_CODES)
    assert version_code_for("eprof_v2.3", 1.0) == 230
    assert version_code_for("eprof_v2.3", 0.5) == 230
    assert version_code_for("eprof_v2.3", 0.25) == 231
    assert version_code_for("eprof_v2.2", 0.25) == 221
    assert version_code_for("eprof_v2", 0.25) == 200          # no low-sun entry -> the method's code
    assert version_code_for("eprof_v2.3", None) == version_code("eprof_v2.3")


def test_station_override_in_options_json():
    from calibration.config import CalibrationOptions
    o = CalibrationOptions.from_json(REPO / "options.json")
    assert o.station_molecular_params["0-20000-0-06610_A"] == {"use_noise": False}
    assert CalibrationOptions().station_molecular_params == {}
    merged = dict(o.molecular_params or {}, **o.station_molecular_params["0-20000-0-06610_A"])
    assert merged["use_noise"] is False


@pytest.mark.parametrize("method,itype,expect", [
    ("eprof_v2.3", "CL61", 5000.0), ("eprof_v2.3", "CHM15k", None), ("eprof_v2.2", "CL61", None)])
def test_cl61_tier2_cap_scope(method, itype, expect):
    from calibration.config import CalibrationOptions, InstrumentType
    from calibration.rayleigh.calibration import _window_search_setup
    o = CalibrationOptions()
    o.molecular_method = method
    params, _ = _window_search_setup(o, InstrumentType(itype), None, 1.0)
    got = (params or {}).get("tier2_max_end_m")
    assert got == expect


def test_tier2_cap_removes_high_windows_only_in_tier2():
    from calibration.rayleigh import molecular_methods as MM
    import inspect
    sig = inspect.signature(MM._select_optimal)
    assert sig.parameters["tier2_max_end_m"].default == np.inf
    src = inspect.getsource(MM._select_optimal)
    assert "tier2_max_end_m" in src.split("relaxed = (base")[1].split("elig = _tvar(relaxed")[0]
    assert "tier2_max_end_m" not in src.split("strict = (base")[1].split("strict = _tvar(strict")[0]
