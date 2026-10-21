# -*- coding: utf-8 -*-
"""q20 - shared definitions for the per-unit correction budget of the CHM15k Rayleigh constant (report 19).

Twenty geographically spread CHM15k units (plus Payerne A, reported separately as the abnormal unit). For
each unit the operational nightly Rayleigh constant C_L (published NetCDF, `common.load_constants`) is
corrected for every effect quantified in report 19 (q2a, q15a-f, q16a-g, q17a-b), one at a time and
cumulatively, and the annual harmonic and the night-to-night scatter are compared before/after.

All corrections act on ln C_L (multiplicative effects). Every correction term is centred on its median
over the unit's nights, so the corrected series keeps the operational level; only the variability is
changed. The order of the cumulative chain follows the task numbering (1 temperature ... 8 other).

Coefficients and their sources
------------------------------
1 temperature   kappa_u (intra-night, elapsed time held fixed; q2a_intranight.csv b_T_ctrl) applied to
                the SYNOPTIC deviation of the calibration-window T_int from its 31-day running median only
                (q16b: kappa applied to the seasonal T swing over-corrects).
2 aerosol       CAMS AOD550 by species, synoptic multiple regression with unit fixed effects (q17a):
                dust 72.9, organic matter 40.3, black carbon 55.3, sulphate 56.7, sea salt 17.3 % per unit
                AOD. Sensitivities: CAMS ln total AOD550 5.34 %/ln AOD; AERONET ln AOD1064 5.15 %/ln AOD.
3 window        d ln C / d ln w = +0.34 (q15b, pooled anomalies, cleaning steps and pairs).
4 module        era offsets at optical_module_id changes, fitted jointly with the annual harmonic (q15a:
                every real laser replacement coincided with a module swap; laser_life_time wraps at
                65 536 h and is not an event).
5 p_calc drift  imported reference ageing: d ln C = -0.69 d ln p_calc on the slow (linear-trend) part only
                (q15b cross-unit Theil-Sen), applied where |trend| >= 5 %/yr (q15b QC threshold).
6 molecular     ln of the slope of the real (CAMS 0.4 deg L137 T/p) molecular signal against the US-1976
                model in each night's own window (q15d F3, same code path).
7 tier-2        nights rescued by the eprof_v2.2 noise-aware fallback tier (q16e: +9.4 % median vs
                tier-1 nights within +-45 d) identified night by night from the q19 record run (which
                reproduces the published constants exactly); removed from the series.
8 other         window geometry and sampling (window bottom/top, duration, clear-profile count) tested
                with a pooled synoptic regression on the 20 units; kept only if significant.
"""
from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
from common import OUT  # noqa: E402
from q1_temperature import running_median  # noqa: E402
from validation.paper.kalman_improve_alc import flag_outliers, kalman_improve_alc  # noqa: E402

Q20 = OUT / "q20"
EPOCH = datetime(1970, 1, 1)
T_START = float((datetime(2025, 1, 1) - EPOCH).days)        # first night start (inclusive)
T_END = float((datetime(2026, 8, 1) - EPOCH).days)          # last night start (exclusive)
W = 2 * np.pi / 365.25
ANOM_WINDOW = 31.0

# (stream key, short label, country, setting) - the 20 study units
UNITS = [
    ("0-20000-0-06240_A", "Schiphol A", "NL", "coastal"),
    ("0-20000-0-06240_B", "Schiphol B", "NL", "coastal"),
    ("0-20000-0-06240_C", "Schiphol C", "NL", "coastal"),
    ("0-20000-0-06240_D", "Schiphol D", "NL", "coastal"),
    ("0-20000-0-12843_A", "Budapest A", "HU", "continental"),
    ("0-20000-0-12843_B", "Budapest B", "HU", "continental"),
    ("0-20000-0-10393_0", "Lindenberg", "DE", "continental"),
    ("0-20000-0-10035_0", "Schleswig", "DE", "coastal"),
    ("0-20000-0-10962_0", "Hohenpeissenberg", "DE", "mountain (977 m)"),
    ("0-20000-0-06736_A", "Kleine Scheidegg", "CH", "alpine (2060 m)"),
    ("0-380-5-1_0", "Aosta", "IT", "alpine valley (560 m)"),
    ("0-20008-0-POT_A", "Potenza", "IT", "hill (760 m)"),
    ("0-20000-0-01492_A", "Oslo", "NO", "coastal (fjord)"),
    ("0-20000-0-03590_A", "Wattisham", "UK", "near coast"),
    ("0-250-1001-07151_B", "Palaiseau", "FR", "continental"),
    ("0-20008-0-UGR_A", "Granada", "ES", "continental (680 m)"),
    ("0-724-60000-0000000CIEMAT001_A", "Madrid", "ES", "continental (669 m)"),
    ("0-620-3704-5480_B", "Evora", "PT", "continental"),
    ("0-100-20000-0000001_A", "Sofia", "BG", "continental (590 m)"),
    ("0-20008-0-INO_A", "Magurele", "RO", "continental"),
]
PAYERNE = ("0-20000-0-06610_A", "Payerne", "CH", "continental (490 m)")
ALL_UNITS = UNITS + [PAYERNE]
LABEL = {k: lab for k, lab, _, _ in ALL_UNITS}
TWIN_PAIRS = [("0-20000-0-06240_A", "0-20000-0-06240_B"), ("0-20000-0-06240_A", "0-20000-0-06240_C"),
              ("0-20000-0-06240_A", "0-20000-0-06240_D"), ("0-20000-0-06240_B", "0-20000-0-06240_C"),
              ("0-20000-0-06240_B", "0-20000-0-06240_D"), ("0-20000-0-06240_C", "0-20000-0-06240_D"),
              ("0-20000-0-12843_A", "0-20000-0-12843_B")]

# Coefficients (see the module docstring for the sources)
K_WIN = 0.34
K_PCAL = 0.69
PCAL_MIN_TREND = 5.0                     # %/yr
AER_SPECIES = {"duaod550": 72.9, "omaod550": 40.3, "bcaod550": 55.3, "suaod550": 56.7, "ssaod550": 17.3}
AER_CAMS_LN = 5.34                       # % per ln AOD550
AER_AERONET_LN = 5.15                    # % per ln AOD1064
TIER2_BIAS = 9.4                         # % (q16e), used only for the 'scale instead of drop' sensitivity

STEPS = ["raw", "kappa", "aerosol", "window", "module", "pcal", "molecular", "tier2", "other"]
STEP_LABEL = {"raw": "operational", "kappa": "1 temperature (synoptic κ)", "aerosol": "2 aerosol (CAMS species)",
              "window": "3 window transmission", "module": "4 module eras", "pcal": "5 p_calc drift",
              "molecular": "6 molecular model", "tier2": "7 tier-2 nights removed", "other": "8 window geometry"}


def doy(t):
    """Day of year (1..366) of days-since-1970 values."""
    return np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in np.atleast_1d(t)], float)


def to_datetime(t):
    return [EPOCH + timedelta(days=float(x)) for x in np.atleast_1d(t)]


def harmonic(t, y, with_trend=False):
    """OLS fit of y = m + a cos(w doy) + b sin(w doy) [+ trend]. Returns (phasor a+ib, mean, fitted, ok)."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    ok = np.isfinite(t) & np.isfinite(y)
    if ok.sum() < 30:
        return np.nan + 1j * np.nan, np.nan, np.full(len(y), np.nan), ok
    d = doy(t)
    cols = [np.ones(len(t)), np.cos(W * d), np.sin(W * d)]
    if with_trend:
        cols.append((t - np.nanmean(t[ok])) / 365.25)
    X = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    return c[1] + 1j * c[2], c[0], X @ c, ok


def peak_doy(p):
    return float(np.mod(np.angle(p) / W, 365.25)) if np.isfinite(p) else np.nan


def inphase(p, pref):
    """Component of phasor p along the direction of the reference phasor pref (same units as p)."""
    if not (np.isfinite(p) and np.isfinite(pref)) or abs(pref) == 0:
        return np.nan
    return float(np.real(p * np.conj(pref)) / abs(pref))


def screen(c):
    """Kalman global-IQR outlier screen (True = kept), as the operational best estimate does."""
    c = np.asarray(c, float)
    good = np.isfinite(c) & (c > 0)
    keep = np.zeros(len(c), bool)
    if good.sum() >= 3:
        keep[good] = ~flag_outliers(c[good])
    return keep


def n2n_scatter(t, lnc):
    """Robust night-to-night scatter [%]: 1.4826 MAD of ln C about its 31-day running median."""
    t = np.asarray(t, float)
    y = np.asarray(lnc, float)
    ok = np.isfinite(t) & np.isfinite(y)
    if ok.sum() < 30:
        return np.nan
    order = np.argsort(t[ok])
    tt, yy = t[ok][order], y[ok][order]
    r = yy - running_median(tt, yy, ANOM_WINDOW)
    r = r[np.isfinite(r)]
    return 100 * 1.4826 * np.median(np.abs(r - np.median(r)))


def series_metrics(t, lnc, pref=None):
    """Annual amplitude [%], peak DOY, in-phase amplitude vs pref [%], n2n scatter [%], model scatter [%].

    The series is first screened with the Kalman global IQR criterion (on C), as the best-estimate filter
    does; all metrics are computed on the kept nights."""
    t = np.asarray(t, float)
    lnc = np.asarray(lnc, float)
    keep = screen(np.exp(lnc)) & np.isfinite(t)
    p, m, fit, ok = harmonic(t[keep], lnc[keep])
    res = lnc[keep] - fit
    out = dict(n=int(keep.sum()), amp=100 * abs(p) if np.isfinite(p) else np.nan, peak=peak_doy(p),
               n2n=n2n_scatter(t[keep], lnc[keep]),
               scat_model=100 * float(np.sqrt(np.nanmean(res ** 2))) if keep.sum() > 30 else np.nan,
               phasor=p)
    out["inphase"] = 100 * inphase(p, pref) if pref is not None else out["amp"]
    return out


def kalman_metrics(t, c):
    """Paper Kalman (kalman_improve_alc) on a nightly series: (grid days, state, model params, amp %, peak DOY)."""
    t = np.asarray(t, float)
    c = np.asarray(c, float)
    ok = np.isfinite(t) & np.isfinite(c) & (c > 0)
    if ok.sum() < 30:
        return None
    k = kalman_improve_alc(to_datetime(t[ok]), c[ok], itype="CHM15k")
    grid = np.array([(g - EPOCH).total_seconds() / 86400.0 for g in k["grid"]])
    amp = peak = np.nan
    if k["model"] == "annual cycle":
        a, ph, off = k["params"]
        amp = 100 * abs(a) / off
        peak = ((np.pi / 2 if a > 0 else 3 * np.pi / 2) - ph) / (2 * np.pi / 365) % 365
    return dict(grid=grid, state=k["state"], std=k["std"], amp=amp, peak=peak, var_meas=k["var_meas"])


def centred(x):
    """x minus its median over finite values (0 where x is not finite)."""
    x = np.asarray(x, float)
    out = np.zeros(len(x))
    ok = np.isfinite(x)
    if ok.any():
        out[ok] = x[ok] - np.median(x[ok])
    return out
