"""Lock the instrument-type -> calibration-method mapping.

Scientific invariant: the weak 910 nm CL31/CL51 must NEVER be Rayleigh-calibrated
(the molecular return is unusable), and CHM15k / Mini-MPL must never be
cloud-calibrated. Rayleigh -> {Mini-MPL, CHM15k, CL61}; Cloud -> {CL31, CL51, CL61}.
These tests bind to the real gate (`_gate_methods`) so the mapping can't silently drift.
"""
import importlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
rnc = importlib.import_module("scripts.run_network_calibration")

ALL_TYPES = ["CL31", "CL51", "CL61", "CHM15k", "Mini-MPL"]
BOTH = {"rayleigh", "cloud"}

# (type, expect_do_rayleigh, expect_do_cloud) when BOTH methods are requested
GATE_TABLE = [
    ("CL31", False, True),      # 910 nm -> cloud only, never rayleigh
    ("CL51", False, True),      # 910 nm -> cloud only, never rayleigh
    ("CHM15k", True, False),    # 1064 nm molecular -> rayleigh only
    ("Mini-MPL", True, False),  # molecular -> rayleigh only
    ("CL61", True, True),       # does both
]


def test_type_sets_are_exact():
    """Locking membership catches the most likely regression: editing the sets."""
    assert rnc.RAYLEIGH_TYPES == {"CL61", "CHM15k", "Mini-MPL"}
    assert rnc.CLOUD_TYPES == {"CL31", "CL51", "CL61"}


@pytest.mark.parametrize("itype", ["CL31", "CL51"])
def test_cl31_cl51_are_never_rayleigh(itype):
    """The headline guarantee the whole file exists for."""
    assert itype not in rnc.RAYLEIGH_TYPES
    do_ray, _ = rnc._gate_methods(itype, BOTH)
    assert do_ray is False


@pytest.mark.parametrize("itype", ["CHM15k", "Mini-MPL"])
def test_chm15k_minimpl_are_never_cloud(itype):
    assert itype not in rnc.CLOUD_TYPES
    _, do_cld = rnc._gate_methods(itype, BOTH)
    assert do_cld is False


@pytest.mark.parametrize("itype,exp_ray,exp_cld", GATE_TABLE)
def test_gate_table_with_both_methods(itype, exp_ray, exp_cld):
    assert rnc._gate_methods(itype, BOTH) == (exp_ray, exp_cld)


@pytest.mark.parametrize("itype", ALL_TYPES)
def test_no_cal_disables_everything(itype):
    assert rnc._gate_methods(itype, BOTH, no_cal=True) == (False, False)


@pytest.mark.parametrize("itype", ALL_TYPES)
def test_method_request_is_honoured(itype):
    # asking for only cloud never runs rayleigh (even on rayleigh types)...
    assert rnc._gate_methods(itype, {"cloud"})[0] is False
    # ...and asking for only rayleigh never runs cloud (even on cloud types).
    assert rnc._gate_methods(itype, {"rayleigh"})[1] is False


def test_sens_omb_method_matches_the_gate():
    """The Kalman method used for OmB/sensitivity must agree with the gate:
    rayleigh types -> 'rayleigh', 910 nm cloud-only types -> 'cloud'."""
    for itype, method in rnc.SENS_OMB_METHOD.items():
        if method == "rayleigh":
            assert itype in rnc.RAYLEIGH_TYPES
        elif method == "cloud":
            assert itype in rnc.CLOUD_TYPES and itype not in rnc.RAYLEIGH_TYPES
        else:
            pytest.fail(f"unexpected SENS_OMB_METHOD for {itype}: {method}")
