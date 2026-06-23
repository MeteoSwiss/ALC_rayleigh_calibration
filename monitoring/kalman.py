"""Vendored E-PROFILE operational Kalman best-estimate, for the dashboard.

The primitives (kalman_predict / kalman_update / constant) are a faithful copy of the
operational routine in
    improve_alc_calib/src/improve_alc_calib/cal_best_estimate.py
used by run_lindenberg_cl61_cal.py. They are vendored here (rather than imported from
that sibling project) so the dashboard stays self-contained and portable to the server
where improve_alc_calib is not checked out.

Magnitude note
--------------
The operational noise defaults are tuned for one instrument's C_L magnitude. fullcal_all
mixes instruments spanning ~1e10..1e12, so kalman_best_estimate() normalizes each station
by its own median, runs the filter with *relative* process noise, then de-normalizes. The
relative parameters match the Lindenberg run (4% drift, 15%/yr seasonal accumulation).
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

# Relative process-noise (series normalized to ~1 before filtering); see module docstring.
_VAR_EPS_CONST = 0.04 ** 2            # day-to-day random-walk drift variance
_VAR_EPS_TEMP = (0.15 / 365.0) ** 2   # seasonal variance accumulated per gap-day


# --- Operational primitives (faithful copy) ---------------------------------

def constant(t, y):
    """Time-independent predict model: f(t)=mean(y), df/dt=0 (pure random walk)."""
    const = float(np.mean(y))
    return (lambda x: const), (lambda x: 0.0), const


def kalman_predict(t, t_0, x_0, var_0, predictfunc,
                   var_eps_const=_VAR_EPS_CONST, var_eps_temp=_VAR_EPS_TEMP):
    x_a = x_0 + (predictfunc(t) - predictfunc(t_0))
    var_eps = var_eps_const + var_eps_temp * (t - t_0).days
    var_a = var_0 + var_eps
    return x_a, var_a


def kalman_update(meas, x_a, var_meas, var_a):
    gain = var_a / (var_a + var_meas)
    x_t = x_a + gain * (meas - x_a)
    var_t = var_a - gain * var_a
    return x_t, var_t


# --- Best-estimate wrapper (normalized) -------------------------------------

def kalman_best_estimate(times, values):
    """Daily Kalman best estimate of a noisy (times, values) calibration series.

    Mirrors run_lindenberg_cl61_cal.kalman_best_estimate: daily-median aggregation,
    rolling-IQR outlier rejection, measurement noise from rolling-mean residuals, and a
    predict-only step on gap days. Operates on a median-normalized copy so one set of
    relative noise parameters fits every instrument magnitude.

    Returns (grid_dates: datetime64[ns], state, std) in the ORIGINAL C_L units, or three
    empty arrays if there are too few points.
    """
    empty = (np.array([], dtype="datetime64[ns]"), np.array([]), np.array([]))
    values = np.asarray(values, dtype=float)
    times = list(times)
    good = np.isfinite(values)
    values = values[good]
    times = [t for t, g in zip(times, good) if g]
    if values.size < 5:
        return empty

    # Normalize by the median so relative process noise applies to any instrument.
    scale = float(np.median(values))
    if not np.isfinite(scale) or scale <= 0:
        return empty
    norm = values / scale

    order = np.argsort([_as_dt(t) for t in times])
    times = [_as_dt(times[i]) for i in order]
    norm = norm[order]

    # --- daily-median aggregation -----------------------------------------
    day_keys = [t.date() for t in times]
    uniq_days = sorted(set(day_keys))
    daily_t = [datetime(d.year, d.month, d.day) for d in uniq_days]
    daily_v = np.array([np.median(norm[[k == d for k in day_keys]]) for d in uniq_days])
    n_days = daily_v.size
    if n_days < 5:
        return empty

    # --- rolling-IQR outlier flag -----------------------------------------
    win = 30 if n_days > 100 else 10
    half = win // 2
    keep = np.ones(n_days, dtype=bool)
    for k in range(n_days):
        lo, hi = max(0, k - half), min(n_days, k + half + 1)
        w = daily_v[lo:hi]
        med = np.median(w)
        q25, q75 = np.percentile(w, [25, 75])
        iqr = q75 - q25
        if daily_v[k] < med - 1.5 * iqr or daily_v[k] > med + 1.5 * iqr:
            keep[k] = False
    clean_t = [t for t, k in zip(daily_t, keep) if k]
    clean_v = daily_v[keep]
    if clean_v.size < 5:
        clean_t, clean_v = daily_t, daily_v

    # --- predict model + measurement-noise variance -----------------------
    predict_func, _, _ = constant(clean_t, clean_v)
    roll = np.copy(clean_v)
    for k in range(clean_v.size):
        lo, hi = max(0, k - 5), min(clean_v.size, k + 6)
        roll[k] = np.mean(clean_v[lo:hi])
    var_meas = float(np.mean((clean_v - roll) ** 2))
    if not np.isfinite(var_meas) or var_meas <= 0:
        var_meas = float(np.var(clean_v)) or 1.0

    # --- contiguous daily grid + Kalman loop ------------------------------
    grid = [clean_t[0] + timedelta(days=k)
            for k in range((clean_t[-1] - clean_t[0]).days + 1)]
    obs = {t.date(): v for t, v in zip(clean_t, clean_v)}
    t_prev = clean_t[0] - timedelta(days=10)
    x_est, var_est = float(roll[0]), float(var_meas)

    state, variance = [], []
    for t_cur in grid:
        y = obs.get(t_cur.date(), np.nan)
        x_a, var_a = kalman_predict(t_cur, t_prev, x_est, var_est, predict_func)
        if np.isfinite(y):
            x_est, var_est = kalman_update(float(y), x_a, var_meas, var_a)
            t_prev = t_cur
            if not np.isfinite(x_est):
                x_est, var_est = x_a, var_a
        else:
            x_est, var_est = x_a, var_a
        state.append(x_est)
        variance.append(var_est)

    # De-normalize back to original C_L units.
    grid_np = np.array([np.datetime64(t) for t in grid])
    return grid_np, np.array(state) * scale, np.sqrt(np.array(variance)) * scale


def _as_dt(t) -> datetime:
    """Coerce a pandas Timestamp / numpy datetime64 / datetime to a python datetime."""
    if isinstance(t, datetime):
        return t
    if hasattr(t, "to_pydatetime"):
        return t.to_pydatetime()
    return datetime.utcfromtimestamp(np.datetime64(t, "s").astype("int64"))
