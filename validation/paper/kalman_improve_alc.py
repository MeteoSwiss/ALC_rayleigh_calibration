"""kalman_improve_alc.py - the operational E-PROFILE Kalman best estimate (improve_alc_calib) applied to
a daily calibration series.

This is the near-real-time filter of improve_alc_calib/run_cal_best_estimate_evaluation_nrt.py,
reproduced step by step so that the paper figures and the Section 6 constants use the SAME filter that
the E-PROFILE service runs on top of the Rayleigh and liquid-cloud calibrations (the dashboard's own
smoothing is for display only):

  1. outlier flagging: one global inter-quartile-range criterion over the whole series
     (value < q25 - 1.5 IQR or > q75 + 1.5 IQR); flagged nights are kept but not assimilated;
  2. prediction model: annual sinusoid A sin(w doy + phi) + C0 (least squares on day-of-year) when
     more than THRESHOLD_ANNUAL_CYCLE valid nights exist, otherwise a constant;
  3. measurement variance: mean squared residual of the valid nights to the fitted model (one scalar);
  4. process noise: var_eps_const per step + var_eps_temp per day of gap, instrument-type specific;
  5. sequential filter on a daily grid from 30 days before the first valid night: predict every day,
     update on days with a valid (non-flagged) night; flag 2 = measured, 1 = predicted.

The primitives (kalman_predict, kalman_update, fitfunc_sin) are imported from improve_alc_calib when the
sibling project is on disk, otherwise an identical local copy is used.

The operational driver carries absolute process-noise values for the CHM15k, Mini-MPL and MPL only.
For the other instrument types (the liquid-cloud constants of the CL31, CL51 and CL61, whose constants
are of order 1e8 and 1) the CHM15k values are applied RELATIVE to the series median, i.e. the same
fractional day-to-day drift (5e9 / 3e11 = 1.7 % per step) and seasonal accumulation
((3e11/365) / 3e11 = 0.27 % per gap-day). This choice must be confirmed when the filter is deployed
on the cloud constants.
"""
from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

THRESHOLD_ANNUAL_CYCLE = 20           # more valid nights than this -> annual-cycle model
IQR_FACTOR = 1.5
START_OFFSET_DAYS = 30
REFERENCE_CL = 3e11                   # CHM15k magnitude the absolute parameters were tuned for
VAR_EPS_CONSTANT = {"CHM15k": 5e9 ** 2, "Mini-MPL": 5e3 ** 2, "MPL": 5e4 ** 2}
VAR_EPS_TEMPORAL = {"CHM15k": (3e11 / 365) ** 2, "Mini-MPL": (3e5 / 365) ** 2, "MPL": (3e11 / 365) ** 2}

_IMPROVE = Path("C:/Users/hervo/OneDrive/Documents/Python/improve_alc_calib/src/improve_alc_calib")
try:
    if _IMPROVE.is_dir() and str(_IMPROVE) not in sys.path:
        sys.path.insert(0, str(_IMPROVE))
    from cal_best_estimate import kalman_predict, kalman_update, fitfunc_sin   # type: ignore
    SOURCE = "improve_alc_calib"
except Exception:                                                              # pragma: no cover
    SOURCE = "vendored copy"
    DOYFREQ = 2 * np.pi / 365

    def fitfunc_sin(doy, amp, phase, offset):
        return amp * np.sin(DOYFREQ * doy + phase) + offset

    def kalman_predict(t, t_0, x_0, var_0, predictfunc, var_eps_const=5e10 ** 2,
                       var_eps_temp=(3e11 / 365) ** 2):
        x_a = x_0 + (predictfunc(t) - predictfunc(t_0))
        var_eps = var_eps_const + (var_eps_temp * (t - t_0).days)
        return x_a, var_0 + var_eps

    def kalman_update(meas, x_a, var_meas, var_a):
        gain = var_a / (var_a + var_meas)
        return x_a + gain * (meas - x_a), var_a - gain * var_a


def flag_outliers(values):
    """Global IQR criterion of flag_outliers_lom (optical-module handling not applied here)."""
    v = np.asarray(values, float)
    q25, q75 = np.nanpercentile(v, 25), np.nanpercentile(v, 75)
    cut = IQR_FACTOR * (q75 - q25)
    return (v < q25 - cut) | (v > q75 + cut)


def fit_model(times, values):
    """(predict function of datetime, parameters, model name) as in the NRT driver."""
    if len(values) > THRESHOLD_ANNUAL_CYCLE:
        doy = np.array([t.timetuple().tm_yday for t in times], float)
        scale = float(np.median(values))                    # a-priori scaled to the series magnitude
        p0 = [scale / 3.0, 3 * np.pi / 2, scale]
        params, _ = curve_fit(fitfunc_sin, doy, np.asarray(values, float), p0, maxfev=20000)
        return (lambda t: fitfunc_sin(t.timetuple().tm_yday, *params)), params, "annual cycle"
    const = float(np.mean(values))
    return (lambda t: const), np.array([const]), "constant"


def process_noise(itype, scale):
    """(var_eps_const, var_eps_temp) for the instrument type; relative scaling for other types."""
    if itype in VAR_EPS_CONSTANT:
        return VAR_EPS_CONSTANT[itype], VAR_EPS_TEMPORAL[itype]
    rel_c = np.sqrt(VAR_EPS_CONSTANT["CHM15k"]) / REFERENCE_CL
    rel_t = np.sqrt(VAR_EPS_TEMPORAL["CHM15k"]) / REFERENCE_CL
    return (rel_c * scale) ** 2, (rel_t * scale) ** 2


def kalman_improve_alc(times, values, itype="CHM15k"):
    """Operational best estimate of a daily calibration series.

    times  : sequence of datetime (one per calibration, any order)
    values : the daily constants (NaN allowed)
    itype  : instrument type for the process-noise parameters

    Returns dict(grid, state, std, flag, outlier_dates, model, params, var_meas): `grid` is the daily
    datetime grid from 30 days before the first valid night to the last one, `flag` is 2 on days with an
    assimilated measurement and 1 on predicted days.
    """
    t = np.array([datetime(x.year, x.month, x.day) if hasattr(x, "year") else x for x in times])
    v = np.asarray(values, float)
    good = np.isfinite(v)
    order = np.argsort(t[good])
    t, v = t[good][order], v[good][order]
    if v.size < 3:
        raise ValueError("fewer than 3 valid calibrations")
    out = flag_outliers(v)
    t_fit, y_fit = t[~out], v[~out]
    predict, params, model = fit_model(list(t_fit), y_fit)
    res = y_fit - np.array([predict(x) for x in t_fit])
    var_meas = float(np.mean(res ** 2))
    var_c, var_t = process_noise(itype, float(np.median(y_fit)))

    obs = {x.date(): (y, o) for x, y, o in zip(t, v, out)}
    t0 = t_fit[0] - timedelta(days=START_OFFSET_DAYS)
    x0, var0 = float(predict(t0)), var_meas
    grid = [t0 + timedelta(days=k) for k in range((t[-1] - t0).days + 1)]
    state, variance, flag = [], [], []
    for tc in grid:
        y, is_out = obs.get(tc.date(), (np.nan, False))
        x_a, var_a = kalman_predict(tc, t0, x0, var0, predict, var_eps_const=var_c, var_eps_temp=var_t)
        if np.isfinite(y) and not is_out:
            x_t, var_t_ = kalman_update(float(y), x_a, var_meas, var_a)
            flag.append(2)
        else:
            x_t, var_t_ = x_a, var_a
            flag.append(1)
        t0, x0, var0 = tc, x_t, var_t_
        state.append(x_t)
        variance.append(var_t_)
    return dict(grid=np.array(grid), state=np.array(state), std=np.sqrt(np.array(variance)),
                flag=np.array(flag), outlier_dates=[x.date() for x in t[out]], model=model,
                params=params, var_meas=var_meas, source=SOURCE)
