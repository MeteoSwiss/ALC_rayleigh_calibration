# -*- coding: utf-8 -*-
"""The OmB coarse vertical grid must reject a corrupt range axis instead of asking numpy for an
astronomical array.

Lerwick's L1 file of 2026-04-01 carries range values from -3e30 to 4e17 m; in the v2.3 release add-ons
``np.arange`` raised "Maximum allowed size exceeded" and the uncaught error cost the stream its whole
OmB pass. ``vertical_edges`` now checks the axis and raises a ValueError the per-day loop of the runner
logs and skips.
"""
import numpy as np
import pytest

from calibration.omb.omb import vertical_edges


def test_regular_axis_gives_the_anchored_grid():
    rng = np.arange(15.0, 15000.0, 15.0)
    edges = vertical_edges(rng, y_max=6000.0, resol_z=120.0)
    assert edges[0] == 15.0
    assert np.allclose(np.diff(edges), 120.0)
    assert edges[-1] >= 6000.0


def test_corrupt_axis_is_rejected_with_a_message():
    rng = np.arange(15.0, 15000.0, 15.0)
    rng[0] = -3.1e30
    rng[-1] = 4.2e17
    with pytest.raises(ValueError, match="degenerate range axis"):
        vertical_edges(rng, y_max=6000.0, resol_z=120.0)


def test_all_nan_axis_is_rejected():
    with pytest.raises(ValueError, match="fewer than two finite gates"):
        vertical_edges(np.full(100, np.nan), y_max=6000.0, resol_z=120.0)


def test_non_increasing_axis_is_rejected():
    rng = np.arange(15.0, 3000.0, 15.0)
    rng[10] = rng[9]
    with pytest.raises(ValueError, match="not strictly increasing"):
        vertical_edges(rng, y_max=6000.0, resol_z=120.0)
