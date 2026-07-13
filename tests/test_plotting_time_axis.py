"""Regression: the Rayleigh diagnostic pcolor uses a real datetime (hours + date) x-axis
instead of 'hours since start'. Tests the helpers + that both diagnostic functions accept
the timestamps (so the runner can keep passing them in future versions)."""
import inspect
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import dates as mdates  # noqa: E402

from calibration.plotting import (  # noqa: E402
    _time_x, _label_time_x,
    plot_rayleigh_diagnostics_compact, plot_rayleigh_diagnostics_failure,
)


def test_time_x_falls_back_to_hours_without_datetimes():
    hrs = np.array([0.0, 0.5, 1.0])
    x, is_dt = _time_x(hrs, None)
    assert is_dt is False and np.allclose(x, hrs)
    # empty list also falls back (a failure plot for a no-data night)
    x2, is_dt2 = _time_x(hrs, [])
    assert is_dt2 is False


def test_time_x_uses_datetimes_when_available():
    t0 = datetime(2026, 3, 15, 18, 0, 0)
    times = [t0 + timedelta(minutes=30 * i) for i in range(5)]   # crosses nothing, monotonic
    x, is_dt = _time_x(np.arange(5.0), times)
    assert is_dt is True
    assert np.allclose(x, mdates.date2num(times))               # matplotlib date numbers
    assert np.all(np.diff(x) > 0)


def test_label_time_x_datetime_axis():
    fig, ax = plt.subplots()
    _label_time_x(ax, True)
    assert ax.get_xlabel() == "Time (UTC)"
    assert isinstance(ax.xaxis.get_major_formatter(), mdates.ConciseDateFormatter)
    plt.close(fig)


def test_label_time_x_hours_fallback():
    fig, ax = plt.subplots()
    _label_time_x(ax, False)
    assert ax.get_xlabel() == "Hours since start"
    plt.close(fig)


def test_diag_functions_accept_time_datetime():
    assert "time_datetime" in inspect.signature(plot_rayleigh_diagnostics_compact).parameters
    assert "time_datetime" in inspect.signature(plot_rayleigh_diagnostics_failure).parameters
