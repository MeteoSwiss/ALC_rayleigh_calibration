# -*- coding: utf-8 -*-
"""Payerne dark profiles v3: same 10 sessions and covered intervals as v2, but smoothed in the
NON-range-corrected signal P = rcs_0 / r^2 with a range-adaptive window, then multiplied back by r^2.

Why: v2 smooths the range-corrected dark with a FIXED 450 m window. The dark noise is roughly
stationary in P, so in rcs units it grows as r^2 and a fixed window leaves gate-to-gate wiggles of
the order of 10 % of the molecular signal in the Rayleigh fit range (3-6 km) and ~50 % at 10 km:
subtracting that profile injects noise into the calibration. In P the dark is a slowly varying
function of range, so it can be smoothed with a window that widens with range (w = max(W_MIN,
W_FRAC * r)) without blurring real structure, and the r^2 factor is applied only afterwards.

Pooling: ONE median per gate over ALL covered profiles of all sessions pooled together (operator
choice 2026-10-09: not a median of per-session medians). The per-session medians are kept only to
give the between-session standard error, i.e. the uncertainty from the evolution of the dark.

Output: <DATA>/dark_profiles_payerne_v3.npz, same keys as v2 (drop-in for ALC_DARK_PROFILE), plus
<ident>_P (smoothed P) and <ident>_P_raw (unsmoothed pooled P) for diagnostics.

Run:  python rayleigh_availability/build_dark_v3.py
"""
from __future__ import annotations
import json
import sys
from datetime import date, datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dark_profiles as DP   # noqa: E402

WINDOWS_JSON = DP.DATA / "dark_windows_detected.json"
V2 = DP.DATA / "dark_profiles_payerne_v2.npz"
OUT = DP.DATA / "dark_profiles_payerne_v3.npz"
MEDIAN_M = 150.0     # median filter first: removes isolated outlier gates
W_MIN = 300.0        # running-mean width at short range [m]
W_FRAC = 0.25        # ... widening as W_FRAC * r at long range (1.25 km at 5 km, 2.5 km at 10 km)


def median_filter(x, k):
    half = k // 2
    pad = np.pad(x, half, mode="edge")
    return np.array([np.median(pad[i:i + k]) for i in range(x.size)])


def smooth_P(rng, P):
    """Median filter, then a running mean whose width grows with range, in P units."""
    dz = float(np.median(np.diff(rng)))
    x = np.asarray(P, float).copy()
    good = np.isfinite(x)
    x[~good] = np.interp(rng[~good], rng[good], x[good])
    x = median_filter(x, max(3, int(round(MEDIAN_M / dz)) | 1))
    c = np.concatenate([[0.0], np.cumsum(x)])
    out = np.empty_like(x)
    for i, r in enumerate(rng):
        h = max(1, int(round(max(W_MIN, W_FRAC * r) / dz / 2)))
        a, b = max(0, i - h), min(x.size, i + h + 1)
        out[i] = (c[b] - c[a]) / (b - a)
    out[~good] = np.nan
    return out


def collect_all(ident):
    """All covered profiles in P units (rcs_0 / r^2), pooled over sessions, plus the session of each."""
    blocks, sess, rng_ref = [], [], None
    for k, (w0, w1) in enumerate(DP.WINDOWS[ident]):
        t0, t1 = DP._parse(w0), DP._parse(w1)
        if ident == "B" and t0 < DP.CL31_SWAP:
            continue                                     # CL31: post optical-module swap only
        day = t0.date()
        while day <= t1.date():
            got = DP._read_day(ident, datetime.combine(day, datetime.min.time()))
            if got is not None:
                times, rng, rcs = got
                if rng_ref is None:
                    rng_ref = rng
                m = (times >= t0) & (times <= t1)
                if m.any() and rng.size == rng_ref.size:
                    blocks.append(rcs[m] / np.where(rng > 0, rng ** 2, np.nan)[None, :])
                    sess.append(np.full(int(m.sum()), k))
            day = date.fromordinal(day.toordinal() + 1)
    return rng_ref, np.vstack(blocks), np.concatenate(sess)


def main():
    win = json.loads(WINDOWS_JSON.read_text())
    DP.MARGIN_S = 0.0
    DP._parse = lambda s: datetime.strptime(s, "%Y-%m-%d %H:%M:%S")
    v2 = dict(np.load(V2))
    out = {}
    for ident, name in (("A", "CHM15k"), ("B", "CL31"), ("C", "CL61")):
        DP.WINDOWS[ident] = [tuple(w) for w in win[ident]]
        rng, X, sess = collect_all(ident)
        r2 = np.where(rng > 0, rng ** 2, np.nan)
        P_raw = np.nanmedian(X, axis=0)                                  # one median over ALL profiles
        P = smooth_P(rng, P_raw)
        # per-session medians: only for the between-session standard error
        S = np.vstack([smooth_P(rng, np.nanmedian(X[sess == k], axis=0)) for k in np.unique(sess)])
        per_window = [dict(n=int((sess == k).sum())) for k in np.unique(sess)]
        b = P * r2
        sem = np.nanstd(S, axis=0) / np.sqrt(S.shape[0]) * r2
        out.update({f"{ident}_range": rng, f"{ident}_b_rcs": b, f"{ident}_b": P, f"{ident}_sem": sem,
                    f"{ident}_nwin": len(per_window), f"{ident}_nprof": int(sum(w["n"] for w in per_window)),
                    f"{ident}_P": P, f"{ident}_P_raw": P_raw, f"{ident}_P_sessions": S})
        band = (rng >= 3000) & (rng <= 6000)
        # gate-to-gate roughness: residual to a 1.5 km running mean, in rcs units
        def rough(x):
            k = max(3, int(round(1500 / np.median(np.diff(rng)))) | 1)
            sm = np.convolve(np.pad(x, k // 2, mode="edge"), np.ones(k) / k, mode="valid")
            return np.nanstd((x - sm)[band])
        print(f"{name}: {len(per_window)} sessions, {out[f'{ident}_nprof']} profiles | 3-6 km: "
              f"mean b v2 {np.nanmean(v2[f'{ident}_b_rcs'][band]):+.3e}  v3 {np.nanmean(b[band]):+.3e} | "
              f"roughness v2 {rough(v2[f'{ident}_b_rcs']):.2e}  v3 {rough(b):.2e} | "
              f"sem v2 {np.nanmean(v2[f'{ident}_sem'][band]):.2e}  v3 {np.nanmean(sem[band]):.2e}")
    np.savez(OUT, **out)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
