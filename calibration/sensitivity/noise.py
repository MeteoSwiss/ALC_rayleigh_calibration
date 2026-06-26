"""Per-gate noise estimators and supporting statistics.

Faithful ports of the MATLAB helper functions used by
``dark_measurement_cl61_chm_cl31.m`` and ``ambient_noise_cl61_chm_cl31.m``.

Conventions
-----------
* ``beta``    : attenuated backscatter, time x range, any consistent unit
                (we use Mm^-1 sr^-1 everywhere downstream).
* ``dtime_s`` : profile time stamps as seconds (float) from an arbitrary
                origin; only differences matter.
* ``dt``      : nominal sampling step [s] (median of diff(dtime_s)).
* ``r``       : range / altitude AGL of each gate [m].

All estimators return *per-gate* sigmas (one value per range gate) plus the
number of samples that contributed, so they can be aggregated into altitude
bins with :func:`bin_rms`.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "robust_std",
    "first_difference_sigma",
    "plain_std",
    "bin_rms",
    "solar_elevation",
    "classify_day_night",
    "overlapping_adev_ci",
    "allan_time",
    "band_timeseries",
]


def robust_std(X: NDArray, axis: int = 0) -> NDArray:
    """Robust scale 1.4826 * MAD along ``axis``.

    Falls back to the plain standard deviation wherever the MAD is degenerate
    (e.g. CL31 quantisation can make it exactly 0). Port of
    ``local_robust_std_dim``.
    """
    X = np.asarray(X, dtype=float)
    med = np.nanmedian(X, axis=axis, keepdims=True)
    s = 1.4826 * np.nanmedian(np.abs(X - med), axis=axis)
    with np.errstate(invalid="ignore"):
        s_fallback = np.nanstd(X, axis=axis)
    s = np.asarray(s, dtype=float)
    bad = ~(s > 0)
    s[bad] = s_fallback[bad]
    return s


def plain_std(beta: NDArray) -> NDArray:
    """Plain per-gate temporal std about the time-mean profile (time x range).

    Removes the static atmosphere but keeps all atmospheric variability, so it
    upper-bounds the noise. Returns one value per range gate.
    """
    beta = np.asarray(beta, dtype=float)
    with np.errstate(invalid="ignore"):
        return np.nanstd(beta, axis=0)


def first_difference_sigma(
    beta: NDArray, dtime_s: NDArray, dt: float, sel: NDArray | None = None
) -> tuple[NDArray, NDArray]:
    """Estimator (b): temporal first difference / sqrt(2), per gate.

    Differences are formed only between consecutive profiles that are both
    selected and contiguous in time (gap < 1.5*dt), never across a data gap or
    a cloud-mask gap. White noise is uncorrelated between profiles, so the
    variance of the difference is twice the per-profile variance; dividing the
    robust sigma of the difference by sqrt(2) recovers the per-profile noise.

    Port of ``local_pairs`` + the first-difference block.

    Returns ``(sigma_gate, n_gate)`` (both length n_range).
    """
    beta = np.asarray(beta, dtype=float)
    dtime_s = np.asarray(dtime_s, dtype=float)
    n = beta.shape[0]
    if sel is None:
        sel = np.ones(n, dtype=bool)
    sel = np.asarray(sel, dtype=bool)

    both = sel[:-1] & sel[1:]
    gap = np.diff(dtime_s)
    pair_i = np.where(both & (gap < 1.5 * dt))[0]
    if pair_i.size == 0:
        nr = beta.shape[1]
        return np.full(nr, np.nan), np.zeros(nr, dtype=int)

    D = beta[pair_i + 1, :] - beta[pair_i, :]
    sigma_gate = robust_std(D, axis=0) / np.sqrt(2.0)
    n_gate = np.sum(np.isfinite(D), axis=0)
    return sigma_gate, n_gate


def bin_rms(
    sig_gate: NDArray,
    n_gate: NDArray,
    r: NDArray,
    z_edges: NDArray,
    min_n: int = 30,
) -> tuple[NDArray, NDArray]:
    """Aggregate per-gate sigmas into altitude bins as the RMS (variance mean).

    Port of ``local_bin_rms``. Bins with fewer than ``min_n`` contributing
    samples are set to NaN. Returns ``(sig_bin, n_bin)`` of length
    ``len(z_edges) - 1``.
    """
    sig_gate = np.asarray(sig_gate, dtype=float).ravel()
    n_gate = np.asarray(n_gate, dtype=float).ravel()
    r = np.asarray(r, dtype=float).ravel()
    nb = len(z_edges) - 1

    sig_bin = np.full(nb, np.nan)
    n_bin = np.zeros(nb, dtype=float)

    idx = np.digitize(r, z_edges) - 1  # 0..nb-1, -1 / nb fall outside
    # Close the top bin like MATLAB discretize: a gate exactly on the last edge
    # belongs to the last bin (digitize would push it to nb and drop it). Use an
    # exact == test so genuinely out-of-range gates (r > z_edges[-1]) stay out.
    idx[r == z_edges[-1]] = nb - 1
    ok = np.isfinite(sig_gate) & (idx >= 0) & (idx < nb)
    if not np.any(ok):
        return sig_bin, n_bin

    var = sig_gate**2
    for b in range(nb):
        m = ok & (idx == b)
        if np.any(m):
            sig_bin[b] = np.sqrt(np.nanmean(var[m]))
            n_bin[b] = np.nansum(n_gate[m])
    sig_bin[n_bin < min_n] = np.nan
    return sig_bin, n_bin


def solar_elevation(dtime, lat: float, lon: float) -> NDArray:
    """Solar elevation [deg] from the Spencer (1971) declination / equation of
    time, accurate to ~0.3 deg - plenty for a day/night classification.

    Port of ``local_solar_elevation``. ``dtime`` may be a numpy datetime64
    array or any array convertible to pandas datetimes (interpreted as UTC).
    """
    import pandas as pd

    t = pd.to_datetime(np.asarray(dtime))
    doy = t.dayofyear.to_numpy(dtype=float)
    # Include sub-second part (pandas .second is integer; MATLAB second() is
    # fractional) so frac_h matches MATLAB exactly on sub-second timestamps.
    sec = t.second.to_numpy(dtype=float) + t.microsecond.to_numpy(dtype=float) / 1e6
    frac_h = (
        t.hour.to_numpy(dtype=float)
        + t.minute.to_numpy(dtype=float) / 60.0
        + sec / 3600.0
    )
    g = 2.0 * np.pi / 365.0 * (doy - 1.0 + (frac_h - 12.0) / 24.0)
    decl = (
        0.006918
        - 0.399912 * np.cos(g)
        + 0.070257 * np.sin(g)
        - 0.006758 * np.cos(2 * g)
        + 0.000907 * np.sin(2 * g)
        - 0.002697 * np.cos(3 * g)
        + 0.001480 * np.sin(3 * g)
    )  # [rad]
    eqtime = 229.18 * (
        0.000075
        + 0.001868 * np.cos(g)
        - 0.032077 * np.sin(g)
        - 0.014615 * np.cos(2 * g)
        - 0.040849 * np.sin(2 * g)
    )  # [min]
    tst = frac_h * 60.0 + eqtime + 4.0 * lon  # true solar time [min]
    ha = (tst / 4.0) - 180.0  # hour angle [deg]
    lat_r = np.deg2rad(lat)
    sin_el = np.sin(lat_r) * np.sin(decl) + np.cos(lat_r) * np.cos(decl) * np.cos(
        np.deg2rad(ha)
    )
    return np.rad2deg(np.arcsin(np.clip(sin_el, -1.0, 1.0)))


def classify_day_night(
    dtime, lat: float, lon: float, elev_day: float = 5.0, elev_night: float = -6.0
) -> NDArray:
    """Per-profile class: 1 = day (elev > elev_day), 2 = night
    (elev < elev_night), 0 = twilight / unused. Matches the ambient script."""
    elev = solar_elevation(dtime, lat, lon)
    cls = np.zeros(elev.shape, dtype=int)
    cls[elev > elev_day] = 1
    cls[elev < elev_night] = 2
    return cls


def band_timeseries(beta: NDArray, r: NDArray, z0: float, halfwin: float) -> NDArray:
    """Mean over the gates of a +/- halfwin band centred on z0 (port of
    ``local_band_ts``). Returns a length-n_time series."""
    r = np.asarray(r, dtype=float)
    sel = (r >= z0 - halfwin) & (r <= z0 + halfwin)
    if not np.any(sel):
        return np.full(beta.shape[0], np.nan)
    return np.nanmean(beta[:, sel], axis=1)


def overlapping_adev_ci(
    y: NDArray, dt: float, tau: float, conf: float = 0.683, min_edf: float = 4.0
) -> tuple[float, float, float]:
    """Overlapping Allan deviation of a 1-D series with a chi-squared CI.

    Port of ``local_oadev_ci``. Refs: Allan & Barnes (1981); IEEE Std
    1139-2008; Riley, NIST SP 1065 (2008); white-FM edf from Howe, Allan &
    Barnes (1981). Returns ``(adev, lo, hi)`` (NaN where undefined).
    """
    from scipy.stats import chi2

    y = np.asarray(y, dtype=float).ravel()
    # linear gap fill, no extrapolation of end values
    y = _fill_linear_interior(y)
    ok = np.isfinite(y)
    if ok.mean() < 0.5:
        return np.nan, np.nan, np.nan
    y = y[ok]
    M = y.size  # frequency-type samples
    # round half away from zero, to match MATLAB round() exactly
    m = max(1, int(np.floor(tau / dt + 0.5)))  # averaging factor
    if M < 2 * m + 1:
        return np.nan, np.nan, np.nan
    x = np.concatenate([[0.0], np.cumsum(y)]) * dt  # integrated (phase), M+1
    idx = np.arange(0, x.size - 2 * m)
    d = x[idx + 2 * m] - 2 * x[idx + m] + x[idx]
    avar = np.sum(d**2) / (2.0 * tau**2 * idx.size)
    adev = np.sqrt(avar)
    edf = (3 * (M - 1) / (2 * m) - 2 * (M - 2) / M) * (4 * m**2) / (4 * m**2 + 5)
    if not np.isfinite(edf) or edf < min_edf:
        return adev, np.nan, np.nan
    alpha = 1.0 - conf
    hi = np.sqrt(edf * avar / chi2.ppf(alpha / 2.0, edf))
    lo = np.sqrt(edf * avar / chi2.ppf(1.0 - alpha / 2.0, edf))
    return adev, lo, hi


def allan_time(x: NDArray, dt: float, tau: float) -> tuple[float, int]:
    """Non-overlapping (classic) Allan deviation of a 1-D series. Port of
    ``local_allan_time``. Returns ``(adev, nseg)``."""
    x = np.asarray(x, dtype=float).ravel()
    m = max(1, int(np.floor(tau / dt + 0.5)))  # round half away from zero (MATLAB)
    nseg = x.size // m
    if nseg < 2:
        return np.nan, nseg
    x = x[: nseg * m]
    seg_means = np.nanmean(x.reshape(nseg, m), axis=1)
    d = np.diff(seg_means)
    return float(np.sqrt(0.5 * np.nanmean(d**2))), nseg


def _fill_linear_interior(y: NDArray) -> NDArray:
    """Linearly interpolate interior NaNs, leave leading/trailing NaNs in place
    (equivalent to MATLAB ``fillmissing(y,'linear','EndValues','none')``)."""
    y = np.asarray(y, dtype=float).copy()
    ok = np.isfinite(y)
    if ok.sum() < 2 or ok.all():
        return y
    n = y.size
    xi = np.arange(n)
    first, last = np.argmax(ok), n - 1 - np.argmax(ok[::-1])
    interior = np.zeros(n, dtype=bool)
    interior[first : last + 1] = True
    fillpos = interior & ~ok
    if np.any(fillpos):
        y[fillpos] = np.interp(xi[fillpos], xi[ok], y[ok])
    return y
