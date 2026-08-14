# -*- coding: utf-8 -*-
"""Indicators for the Rayleigh availability work — ONE implementation, used by every phase.

The user's minimum set is availability, Kalman-outlier count and change in calibration constants;
the rest exist because availability alone can always be bought with outliers, and a gate change is
only acceptable if precision, continuity and agreement with independent references all hold.

  1. availability      A = n_valid / n_clear     (a "clear" night = the pipeline reached the fit,
                       i.e. flag not in {0, -1}; that is the denominator the gates actually act on)
  2. kalman_outliers   count of nights the operational rolling-IQR test rejects (monitoring.kalman
                       rule: window 30 if >100 days else 10, median +/- 1.5*IQR)
  3. continuity        on nights valid under BOTH configs: median and p95 of |C_new/C_ref - 1|
  4. sigma_sd          1.4826*median(|diff C|)/sqrt(2)/|median C| -- the existing short-term
                       variability metric, so numbers stay comparable with the C8 study
  5. flag_mix          share of each rejection flag over clear nights (shows WHERE nights moved to;
                       a gate change that merely converts -2 into -9 is not a gain)
  6. noise_regression  yield gain vs the stream's measured night noise -- the mechanism check: gains
                       must concentrate at NOISY streams (admitting noisy-clean nights), not at
                       clean ones (which would mean we simply loosened the aerosol defence)

All functions take plain dicts/arrays so they work on baseline JSONs, sweep JSONs and archive CSVs
alike. Nothing here plots.
"""
from __future__ import annotations
import numpy as np

VALID_FLAGS = (1.0, 0.5)
NOT_CLEAR_FLAGS = (0.0, -1.0)      # no data / not a clear night: never reached the gates


# --------------------------------------------------------------------------- basic per-stream
def is_valid(flag):
    return flag in VALID_FLAGS


def is_clear(flag):
    """Did this night reach the molecular fit (i.e. is it in the denominator the gates act on)?"""
    return flag is not None and flag not in NOT_CLEAR_FLAGS and flag != -99.0


def availability(rec):
    """rec: {date: [flag, cl, ...]} -> dict(n_nights, n_clear, n_valid, availability_pct)."""
    flags = [v[0] for v in rec.values()]
    n_clear = sum(is_clear(f) for f in flags)
    n_valid = sum(is_valid(f) for f in flags)
    return dict(n_nights=len(flags), n_clear=n_clear, n_valid=n_valid,
                availability_pct=(100.0 * n_valid / n_clear) if n_clear else float("nan"))


def flag_mix(rec):
    """Share of every flag over CLEAR nights -> {flag: pct}. Shows where rejected nights went."""
    flags = [v[0] for v in rec.values()]
    clear = [f for f in flags if is_clear(f)]
    if not clear:
        return {}
    out = {}
    for f in set(clear):
        out[f] = 100.0 * sum(1 for x in clear if x == f) / len(clear)
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def constants(rec):
    """(dates, C) for the valid nights, date-sorted."""
    items = sorted((d, v[1]) for d, v in rec.items()
                   if is_valid(v[0]) and v[1] is not None and np.isfinite(v[1]) and v[1] > 0)
    return [d for d, _ in items], np.array([c for _, c in items], float)


def sigma_sd(rec):
    """Short-term variability [%]: 1.4826*median|diff C|/sqrt(2) / |median C| (the C8 definition).

    Night-to-night differences, so a slow seasonal drift does not count as noise. Needs >=3 nights.
    """
    _, c = constants(rec)
    if c.size < 3:
        return float("nan")
    med = np.median(c)
    if not np.isfinite(med) or med == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2.0) / abs(med) * 100.0)


def window_heights(rec):
    """(dates, window bottom height [m]) for the valid nights."""
    items = sorted((d, v[3]) for d, v in rec.items()
                   if is_valid(v[0]) and v[3] is not None and np.isfinite(v[3]))
    return [d for d, _ in items], np.array([h for _, h in items], float)


def altitude_gradient(rec, min_n=25):
    """How much the retrieved C_L depends on WHERE the molecular window was placed.

    A true lidar constant cannot depend on the fit altitude, so a non-zero gradient is an
    unmodelled profile defect (background over/under-subtraction, overlap residual, or aerosol the
    molecular model does not capture). Measured two ways over a stream's valid nights:
      rho        Spearman of C_L against window bottom height (robust, unit-free)
      slope_pct_per_km  log-linear slope, i.e. the % change in C_L per km of window height

    This matters here because a gate change that recovers nights by fitting HIGHER will shift the
    constant by (gradient x altitude shift) at any station whose gradient is non-zero -- an
    altitude confound, not a noise effect. Healthy co-located instruments show ~0 gradient.
    """
    dates, c = constants(rec)
    _, h = window_heights(rec)
    if c.size < min_n or h.size != c.size:
        return dict(n=int(c.size), rho=float("nan"), slope_pct_per_km=float("nan"))
    ok = np.isfinite(h) & np.isfinite(c) & (c > 0)
    if ok.sum() < min_n:
        return dict(n=int(ok.sum()), rho=float("nan"), slope_pct_per_km=float("nan"))
    h, c = h[ok], c[ok]
    rho = float(np.corrcoef(_rank(h), _rank(c))[0, 1])
    slope = float(np.polyfit(h / 1000.0, np.log(c), 1)[0] * 100.0)
    return dict(n=int(ok.sum()), rho=rho, slope_pct_per_km=slope)


def altitude_shift(rec_new, rec_ref):
    """Median window-height difference between newly admitted and retained nights [m].

    Combined with altitude_gradient this predicts how much of a constant offset is an altitude
    confound rather than a real change: expected_offset_pct ~ gradient * shift_km.
    """
    added = newly_admitted(rec_new, rec_ref)
    kept = [d for d in rec_new if is_valid(rec_new[d][0])
            and d in rec_ref and is_valid(rec_ref[d][0])]
    ha = np.array([rec_new[d][3] for d in added
                   if rec_new[d][3] is not None and np.isfinite(rec_new[d][3])], float)
    hk = np.array([rec_new[d][3] for d in kept
                   if rec_new[d][3] is not None and np.isfinite(rec_new[d][3])], float)
    if ha.size < 3 or hk.size < 3:
        return dict(n_added=int(ha.size), n_kept=int(hk.size), shift_m=float("nan"))
    return dict(n_added=int(ha.size), n_kept=int(hk.size),
                shift_m=float(np.median(ha) - np.median(hk)))


def offset_added_vs_kept(rec_new, rec_ref, month_matched=True):
    """Relative difference between the newly admitted and the retained constants [%].

    month_matched=True compares within each calendar month before aggregating, which removes the
    seasonal confound (a gate change that admits summer nights would otherwise appear biased
    simply because C_L has a seasonal cycle).
    """
    added = {d: rec_new[d][1] for d in newly_admitted(rec_new, rec_ref) if rec_new[d][1]}
    kept = {d: rec_new[d][1] for d in rec_new
            if is_valid(rec_new[d][0]) and d in rec_ref and is_valid(rec_ref[d][0]) and rec_new[d][1]}
    if len(added) < 3 or len(kept) < 3:
        return float("nan")
    if not month_matched:
        return 100.0 * (np.median(list(added.values())) / np.median(list(kept.values())) - 1.0)
    ratios = []
    for m in sorted({d[:6] for d in kept} & {d[:6] for d in added}):
        a = [v for d, v in added.items() if d[:6] == m]
        k = [v for d, v in kept.items() if d[:6] == m]
        if len(a) >= 3 and len(k) >= 3:
            ratios.append(np.median(a) / np.median(k))
    return 100.0 * (np.median(ratios) - 1.0) if ratios else float("nan")


def continuity(rec_new, rec_ref):
    """|C_new/C_ref - 1| on nights valid under BOTH -> dict(n, median_pct, p95_pct, max_pct).

    Guards against a gate change that silently moves the constant on nights that were already fine.
    """
    common = [d for d in rec_new if d in rec_ref
              and is_valid(rec_new[d][0]) and is_valid(rec_ref[d][0])
              and rec_new[d][1] and rec_ref[d][1]]
    if not common:
        return dict(n=0, median_pct=float("nan"), p95_pct=float("nan"), max_pct=float("nan"))
    r = np.array([abs(rec_new[d][1] / rec_ref[d][1] - 1.0) for d in common]) * 100.0
    return dict(n=len(common), median_pct=float(np.median(r)),
                p95_pct=float(np.percentile(r, 95)), max_pct=float(r.max()))


def newly_admitted(rec_new, rec_ref):
    """Dates valid under new but NOT under ref -- the nights a change actually buys."""
    return sorted(d for d, v in rec_new.items()
                  if is_valid(v[0]) and not (d in rec_ref and is_valid(rec_ref[d][0])))


def lost(rec_new, rec_ref):
    """Dates valid under ref but NOT under new -- the nights a change costs."""
    return sorted(d for d, v in rec_ref.items()
                  if is_valid(v[0]) and not (d in rec_new and is_valid(rec_new[d][0])))


# --------------------------------------------------------------------------- Kalman outliers
def kalman_outliers(rec):
    """Nights the operational rolling-IQR test rejects, via monitoring.kalman's own rule.

    Reimplemented here (rather than importing) so the indicator works on any {date: C} series
    without a DataFrame, but the rule is byte-for-byte the one in monitoring/kalman.py:92-102:
    window 30 if n_days > 100 else 10, centred, point included in its own window, fence
    median +/- 1.5*IQR, and the whole rejection discarded if fewer than 5 points survive.
    """
    dates, c = constants(rec)
    n = c.size
    if n < 5:
        return dict(n_valid=n, n_outliers=0, outlier_dates=[], rate_per_100=float("nan"))
    scale = np.median(c)
    v = c / scale if scale else c
    win = 30 if n > 100 else 10
    half = win // 2
    keep = np.ones(n, bool)
    for k in range(n):
        w = v[max(0, k - half):min(n, k + half + 1)]
        med = np.median(w)
        q25, q75 = np.percentile(w, [25, 75])
        iqr = q75 - q25
        if v[k] < med - 1.5 * iqr or v[k] > med + 1.5 * iqr:
            keep[k] = False
    if keep.sum() < 5:                                  # the module's safety valve
        keep = np.ones(n, bool)
    bad = [dates[i] for i in range(n) if not keep[i]]
    return dict(n_valid=n, n_outliers=len(bad), outlier_dates=bad,
                rate_per_100=100.0 * len(bad) / n)


# --------------------------------------------------------------------------- aggregation
def by_month(rec):
    """{YYYYMM: {...availability...}} -- the seasonal view (summer is the failure mode)."""
    months = {}
    for d, v in rec.items():
        months.setdefault(d[:6], {})[d] = v
    return {m: availability(r) for m, r in sorted(months.items())}


def season_of(yyyymm):
    """'summer' = Jun-Aug (the C8 blind spot), 'winter' = Nov-Feb, else 'shoulder'."""
    m = int(yyyymm[4:6])
    return "summer" if m in (6, 7, 8) else ("winter" if m in (11, 12, 1, 2) else "shoulder")


def by_season(rec):
    out = {}
    for d, v in rec.items():
        out.setdefault(season_of(d[:6]), {})[d] = v
    return {s: availability(r) for s, r in sorted(out.items())}


def summarise(rec, rec_ref=None):
    """One row of every indicator for one stream+config."""
    row = availability(rec)
    row["sigma_sd_pct"] = sigma_sd(rec)
    row.update({f"kal_{k}": v for k, v in kalman_outliers(rec).items() if k != "outlier_dates"})
    row["flag_mix"] = flag_mix(rec)
    row["by_season"] = {s: a["availability_pct"] for s, a in by_season(rec).items()}
    if rec_ref is not None:
        row["continuity"] = continuity(rec, rec_ref)
        row["n_newly_admitted"] = len(newly_admitted(rec, rec_ref))
        row["n_lost"] = len(lost(rec, rec_ref))
    return row


def noise_regression(gain_by_stream, noise_by_stream):
    """Spearman correlation of per-stream availability GAIN vs measured night noise.

    The mechanism check for the whole hypothesis: if the gates were failing because of photon
    noise, the streams that gain must be the NOISY ones. A gain uncorrelated with noise would mean
    we simply weakened the aerosol defence everywhere instead.
    """
    keys = [k for k in gain_by_stream if k in noise_by_stream
            and np.isfinite(gain_by_stream[k]) and np.isfinite(noise_by_stream[k])]
    if len(keys) < 5:
        return dict(n=len(keys), spearman=float("nan"))
    g = np.array([gain_by_stream[k] for k in keys], float)
    s = np.array([noise_by_stream[k] for k in keys], float)
    rg, rs = _rank(g), _rank(s)
    rho = float(np.corrcoef(rg, rs)[0, 1])
    return dict(n=len(keys), spearman=rho)


def _rank(a):
    order = np.argsort(a, kind="mergesort")
    r = np.empty(a.size, float)
    r[order] = np.arange(a.size, dtype=float)
    return r
