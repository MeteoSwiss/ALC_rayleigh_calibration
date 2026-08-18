"""Regression: the production molecular-window method is eprof_v2 EVERYWHERE — the operational
options.json, the CalibrationOptions dataclass default, and the from_json fallback. Guards against
silently reverting to eprof_v1.2."""
import json
from pathlib import Path

from calibration import CalibrationOptions
from calibration.rayleigh.rayleigh_fit import find_optimal_molecular_window

REPO = Path(__file__).resolve().parents[1]


def test_dataclass_default_is_v2():
    assert CalibrationOptions().molecular_method == "eprof_v2"


def test_operational_options_json_is_v2():
    oj = json.loads((REPO / "options.json").read_text(encoding="utf-8"))
    assert oj["molecular_method"] == "eprof_v2"


def test_from_json_loads_v2():
    assert CalibrationOptions.from_json(REPO / "options.json").molecular_method == "eprof_v2"


def test_find_optimal_window_default_is_v2():
    import inspect
    assert inspect.signature(find_optimal_molecular_window).parameters["method"].default == "eprof_v2"
