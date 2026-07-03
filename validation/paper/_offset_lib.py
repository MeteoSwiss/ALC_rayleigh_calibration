# -*- coding: utf-8 -*-
"""Shared library for network-scale ceilometer electronic-offset extraction from E-PROFILE L1.

The fixed additive electronic offset (Kotthaus et al. 2016 instrument background P^bgi) is estimated
from CLEAR-NIGHT per-gate medians of P = rcs_0 / z^2 (homoscedastic-noise space): the offset is a
FIXED range pattern, so it survives the per-gate median while the VARIABLE atmosphere averages toward
a smooth molecular/aerosol baseline. A gaussian high-pass removes that smooth baseline and isolates
the oscillatory ripple component (the part cleanly recoverable without a covered/hood measurement).

Screening: nighttime (low solar background), cloud-free (all CBH layers == fill), cleanest aerosol
percentile. Internal temperatures (temperature_laser; CL61 also temp_int) are pooled per profile so
the offset can be split into temperature bins (Kotthaus's laser heat-sink temperature classes)."""
import glob
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

L1ROOT = Path("D:/E-PROFILE_L1_2026")
C_LIGHT = 2.99792458e8


def _to_dt(tv):
    return np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in tv])


def read_file(f, temp_var="temperature_laser"):
    """Read one L1 daily file -> (datetimes, range, rcs_0[time,range], cbh[time,nlayer], temp[time])."""
    with Dataset(f) as nc:
        tt = _to_dt(np.asarray(nc.variables["time"][:], "f8"))
        r = np.asarray(nc.variables["range"][:], "f8")
        x = np.asarray(nc.variables["rcs_0"][:], "f8")
        if x.shape[0] != tt.size:
            x = x.T
        cbh = None
        if "cloud_base_height" in nc.variables:
            cbh = np.asarray(nc.variables["cloud_base_height"][:], "f8")
            if cbh.ndim == 1:
                cbh = cbh[:, None]
            if cbh.shape[0] != tt.size:
                cbh = cbh.T
        temp = None
        if temp_var and temp_var in nc.variables:
            temp = np.asarray(nc.variables[temp_var][:], "f8").ravel()
            if temp.size != tt.size:
                temp = None
    return tt, r, x, cbh, temp


def list_files(wmo, ident, months, year="2026"):
    fs = []
    for mm in months:
        fs += sorted(glob.glob(str(L1ROOT / wmo / year / mm / f"L1_{wmo}_{ident}*.nc")))
    return fs


def load_clearnight(wmo, ident, months=("03", "04", "05", "06"), year="2026",
                    night=(22, 3), aerosol_keep_pct=50, temp_var="temperature_laser",
                    aer_lo=1500.0, aer_hi=4000.0, max_profiles=200000):
    """Pool clear-night profiles. Returns dict(rng, X, temp, n_pool, n_kept, n_files) or None."""
    pool_x, pool_t, rng = [], [], None
    n_pool = 0
    for f in list_files(wmo, ident, months, year):
        try:
            tt, r, x, cbh, temp = read_file(f, temp_var)
        except Exception:
            continue
        if x.ndim != 2 or x.shape[1] < 50:
            continue
        hr = np.array([d.hour + d.minute / 60.0 for d in tt])
        lo, hi = night
        night_m = (hr >= lo) | (hr <= hi) if lo > hi else (hr >= lo) & (hr <= hi)
        nocloud = np.nanmax(cbh, axis=1) < 0 if cbh is not None else np.ones(tt.size, bool)
        sel = night_m & nocloud
        if not sel.any():
            continue
        pool_x.append(x[sel]); rng = r
        pool_t.append(temp[sel] if temp is not None else np.full(int(sel.sum()), np.nan))
        n_pool += int(sel.sum())
    if not pool_x:
        return None
    X = np.vstack(pool_x); T = np.concatenate(pool_t)
    if X.shape[0] > max_profiles:                # cap for memory on very long records
        idx = np.linspace(0, X.shape[0] - 1, max_profiles).astype(int)
        X, T = X[idx], T[idx]
    zkm = rng / 1000.0
    P = X / zkm[None, :] ** 2
    sm = uniform_filter1d(np.nan_to_num(P), 9, axis=1)
    aer = np.nanmedian(sm[:, (rng >= aer_lo) & (rng <= aer_hi)], axis=1)
    keep = aer < np.nanpercentile(aer, aerosol_keep_pct)
    return dict(rng=rng, X=X[keep], temp=T[keep], n_pool=n_pool, n_kept=int(keep.sum()),
                n_files=len(list_files(wmo, ident, months, year)))


def offset_stats(X, rng):
    """Per-gate robust median of P=rcs_0/z^2 and its MAD-based precision."""
    zkm = rng / 1000.0
    P = X / zkm[None, :] ** 2
    Pmed = np.nanmedian(P, axis=0)
    Pmad = 1.4826 * np.nanmedian(np.abs(P - Pmed[None, :]), axis=0)
    N = np.isfinite(P).sum(axis=0)
    se = Pmad / np.sqrt(np.maximum(N, 1))
    return Pmed, Pmad, se, N


def highpass(y, dr, sigma_m=350.0):
    """Remove the smooth molecular/aerosol baseline; keep the electronic ripple."""
    return y - gaussian_filter1d(np.nan_to_num(y), sigma_m / dr)


def autocorr(y, dr):
    yy = np.nan_to_num(y - np.nanmean(y))
    ac = np.correlate(yy, yy, "full")[len(yy) - 1:]
    return np.arange(len(ac)) * dr, ac / (ac[0] if ac[0] else 1.0)


def dominant_period(hp, rng, dr, lo, hi, pmin=20.0, pmax=1500.0):
    """First strong autocorrelation peak in [lo,hi] -> (Lambda_m, peak_height, ripple_rms)."""
    from scipy.signal import find_peaks
    m = (rng >= lo) & (rng <= hi)
    if m.sum() < 50:
        return np.nan, np.nan, np.nan
    lag, ac = autocorr(hp[m], dr)
    good = (lag >= pmin) & (lag <= pmax)
    if not good.any():
        return np.nan, np.nan, np.nan
    pk, props = find_peaks(ac, height=0.15)
    pk = [p for p in pk if pmin <= lag[p] <= pmax]
    rms = float(np.sqrt(np.nanmean(hp[m] ** 2)))
    if not pk:
        return np.nan, float(np.max(ac[good])), rms
    p0 = pk[int(np.argmax(ac[pk]))]
    return float(lag[p0]), float(ac[p0]), rms


def fit_ripple(hp, rng, Lam, lo, hi):
    """Least-squares period-Lam ripple (fundamental + Nyquist harmonic), undamped. Returns dict."""
    m = (rng >= lo) & (rng <= hi)
    r = rng[m]; y = np.nan_to_num(hp[m])
    A = np.c_[np.ones_like(r), np.cos(2 * np.pi * r / Lam), np.sin(2 * np.pi * r / Lam),
              np.cos(2 * np.pi * r / (Lam / 2)), np.sin(2 * np.pi * r / (Lam / 2))]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ coef
    R2 = 1 - np.sum((fit - y) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-30)
    return dict(coef=coef, amp_fund=float(np.hypot(coef[1], coef[2])),
                amp_harm=float(np.hypot(coef[3], coef[4])), R2=float(R2),
                f_MHz=C_LIGHT / (2 * Lam) / 1e6)
