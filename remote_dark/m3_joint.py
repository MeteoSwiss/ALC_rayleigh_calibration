# -*- coding: utf-8 -*-
"""M3 — joint decomposition of the dark's dynamics, hood side and sky side.

The redesign the twilight failure forced. Within one diurnal cycle the background B and the
internal temperature T are collinear (both sun-driven), so "the B-coupling" and "the T-coupling"
are not separately identifiable from one session — and pretending otherwise is how the first
twilight estimator produced noise. What IS identifiable:

  * the SLOW coupling  D_slow(z)  [per °C]  — the dark's response to the diurnal thermal state,
    regressed on T with B's slow part deliberately absorbed into it (stated, not hidden: this is
    "the temperature-attributed diurnal response", the quantity Le & O'Connor tabulate as
    P_instrument(r, T));
  * the FAST coupling  D_fast(z)  [per B-unit] — the response to background changes too fast for
    the thermal mass (cloud shadows, ramps; high-passed B at 45 min).

Hood side (`run_hood`): the four multi-hour sessions (A, B pre-swap, B post-swap, C — all ≥20 h)
give the TRUTH for both couplings, with no atmosphere in front.

Sky side (`run_sky_T`): on clear nights, the same slow regression per gate on (signal vs T within
the night, linear drift removed) — the atmosphere is the contaminant, handled by night-median
aggregation and by only claiming gates where many nights agree. Compared against the hood D_slow:
that comparison is the M3 verdict.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta

import numpy as np

from remote_dark.common import PAYERNE, ensure_out, read_day
from remote_dark import hood
from remote_dark.m3_twilight import _highpass

MIN_PROF = 60


def _slowfast_fit(times, Y, Bv, Tv, rng):
    """Per gate: Y(t,z) = c0 + D_slow·T̃ + D_fast·B_hp + c1·τ  (2 IRLS passes, 4-MAD clip).

    T̃ is centred T in °C; B_hp the 45-min high-passed background; τ centred hours (absorbs
    monotone drifts that are neither: laser aging within the session, pressure).
    Returns dict with D_slow, D_fast, their formal SEs, and the T–B_hp collinearity.
    """
    hrs = np.array([(x - times[0]).total_seconds() / 3600.0 for x in times])
    if not np.isfinite(Bv).any():
        # CL61: no background variable in the L1 -> far-gate proxy (mean rcs/z^2 over the top 8 %
        # of gates, the CompCor noise-ROI idea). Under the hood this is a slice of the response
        # itself -- mild endogeneity, accepted and stated; those gates are excluded from the fit
        # output by construction of the band metrics.
        top = rng >= np.nanpercentile(rng, 92)
        with np.errstate(invalid="ignore"):
            Bv = np.nanmean(Y[:, top] / rng[top] ** 2, axis=1)
    good = np.isfinite(Bv) & np.isfinite(Tv)
    if good.sum() < MIN_PROF:
        return None
    Bh = _highpass(Bv, times)
    good &= np.isfinite(Bh)
    if good.sum() < MIN_PROF:
        return None
    t = hrs[good] - hrs[good].mean()
    Tt = Tv[good] - np.nanmean(Tv[good])
    Bn = Bh[good]
    X = np.column_stack([np.ones_like(t), Tt, Bn, t])
    Y = Y[good]
    W = np.isfinite(Y)
    Yf = np.where(W, Y, 0.0)
    beta = None
    for _ in range(2):
        XtX = np.einsum("ti,tj,tz->zij", X, X, W.astype(float))
        XtY = np.einsum("ti,tz->zi", X, Yf * W)
        beta = np.full((rng.size, 4), np.nan)
        var = np.full((rng.size, 4), np.nan)
        for z in range(rng.size):
            if W[:, z].sum() < MIN_PROF:
                continue
            try:
                Ai = np.linalg.inv(XtX[z] + 1e-12 * np.eye(4))
            except np.linalg.LinAlgError:
                continue
            beta[z] = Ai @ XtY[z]
            resid_z = Y[:, z] - X @ beta[z]
            s2 = np.nanvar(np.where(W[:, z], resid_z, np.nan))
            var[z] = s2 * np.diag(Ai)
        resid = Y - X @ beta.T
        mad = np.nanmedian(np.abs(resid - np.nanmedian(resid, axis=0)), axis=0) + 1e-30
        W = np.isfinite(Y) & (np.abs(resid) < 4.0 * mad)
        Yf = np.where(W, Y, 0.0)
    return {"D_slow": beta[:, 1], "D_fast": beta[:, 2],
            "se_slow": np.sqrt(var[:, 1]), "se_fast": np.sqrt(var[:, 2]),
            "corr_T_Bhp": float(np.corrcoef(Tt, Bn)[0, 1]), "n_prof": int(good.sum())}


def run_hood():
    """The truth: slow/fast couplings of every multi-hour hood session."""
    out = {}
    for ident in "ABC":
        for s in hood.session_frames(ident, min_hours=10):
            r = _slowfast_fit(s["times"], s["dark"], s["hk"]["bckgrd"], s["hk"]["t_int"],
                              s["rng"])
            if r is None:
                continue
            key = f"{ident}_{s['era']}_{s['t0']:%Y%m%d}"
            r["rng"] = s["rng"]
            out[key] = r
            m = (s["rng"] >= 60) & (s["rng"] <= 500)
            print(f"  hood {key}: {r['n_prof']} prof  corr(T,B_hp)={r['corr_T_Bhp']:+.2f}  "
                  f"|D_slow| 60-500m median={np.nanmedian(np.abs(r['D_slow'][m])):.3g}/degC  "
                  f"|D_fast| ={np.nanmedian(np.abs(r['D_fast'][m])):.3g}/B")
    return out


def run_sky_T(ident: str, nights: list[str], max_nights: int = 120):
    """The sky-side slow coupling: within-night T regression on clear nights, night-median pooled.

    The atmosphere is present, so each night's D_slow is contaminated by whatever atmospheric
    evolution correlates with T that night; the MEDIAN over ~100 nights keeps only what repeats
    with the instrument. No B_hp term at night (no sun); τ handles the monotone part.
    """
    wmo = PAYERNE["wmo"]
    per = []
    rng = None
    used = 0
    for ds in nights[:max_nights]:
        day = datetime.strptime(ds, "%Y%m%d")
        d = read_day(wmo, ident, day)
        if d is None:
            continue
        # the same morning-night segment convention as the v4 extraction: solar elevation < -6 deg
        from remote_dark.common import solar_elevation_deg
        el = solar_elevation_deg(PAYERNE["lat"], PAYERNE["lon"], d["times"])
        m = el < -6.0
        if m.sum() < MIN_PROF:
            continue
        times, Y = d["times"][m], d["rcs"][m]
        Tv = d["hk"]["t_int"][m]
        if not np.isfinite(Tv).any() or np.nanstd(Tv) < 0.3:
            continue                       # no thermal sweep this night -> no information
        rng = d["rng"]
        hrs = np.array([(x - times[0]).total_seconds() / 3600.0 for x in times])
        t = hrs - hrs.mean()
        Tt = Tv - np.nanmean(Tv)
        X = np.column_stack([np.ones_like(t), Tt, t])
        W = np.isfinite(Y)
        Yf = np.where(W, Y, 0.0)
        XtX = np.einsum("ti,tj,tz->zij", X, X, W.astype(float))
        XtY = np.einsum("ti,tz->zi", X, Yf * W)
        D = np.full(rng.size, np.nan)
        for z in range(rng.size):
            if W[:, z].sum() < MIN_PROF:
                continue
            try:
                D[z] = np.linalg.solve(XtX[z] + 1e-12 * np.eye(3), XtY[z])[1]
            except np.linalg.LinAlgError:
                continue
        per.append(D)
        used += 1
    if not per:
        return None
    P = np.vstack(per)
    return {"rng": rng, "D_slow_med": np.nanmedian(P, axis=0),
            "D_slow_mad": np.nanmedian(np.abs(P - np.nanmedian(P, axis=0)), axis=0),
            "n_nights": used}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sky", action="store_true", help="also run the sky-side T regression (slow)")
    ap.add_argument("--max-nights", type=int, default=120)
    a = ap.parse_args()

    print("== hood joint decomposition (truth)")
    hd = run_hood()
    out = ensure_out("m3")
    payload = {}
    for k, v in hd.items():
        for f in ("D_slow", "D_fast", "se_slow", "se_fast"):
            payload[f"{k}_{f}"] = v[f]
        payload[f"{k}_rng"] = v["rng"]
        payload[f"{k}_corr"] = np.array(v["corr_T_Bhp"])
    np.savez(out / "hood_joint.npz", **payload)
    print(f"   -> {out / 'hood_joint.npz'}")

    if not a.sky:
        return
    from rayleigh_availability import dark_from_clearsky as v4
    for ident in ("A", "C"):
        nights = v4.clear_nights_from_csv(v4.DEFAULT_CSV, PAYERNE["wmo"],
                                          "A" if ident == "B" else ident,
                                          "20250101", "20260813")
        print(f"== sky T regression {ident}: {len(nights)} clear nights available")
        r = run_sky_T(ident, nights, a.max_nights)
        if r is None:
            print("   nothing usable")
            continue
        np.savez(out / f"sky_T_{ident}.npz", **r)
        # the verdict: sky D_slow vs the hood D_slow of the same instrument
        truth_keys = [k for k in hd if k.startswith(ident)]
        for tk in truth_keys:
            hv = hd[tk]
            m = ((r["rng"] >= 60) & (r["rng"] <= 1500)
                 & np.isfinite(r["D_slow_med"]) & np.isfinite(np.interp(
                     r["rng"], hv["rng"], hv["D_slow"])))
            hint = np.interp(r["rng"], hv["rng"], hv["D_slow"])
            if m.sum() < 10:
                continue
            cc = float(np.corrcoef(r["D_slow_med"][m], hint[m])[0, 1])
            gain = float(np.nansum(r["D_slow_med"][m] * hint[m])
                         / max(np.nansum(hint[m] ** 2), 1e-30))
            print(f"   sky vs hood {tk} (60-1500 m): corr {cc:+.2f}  gain {gain:+.2f}  "
                  f"({r['n_nights']} nights)")


if __name__ == "__main__":
    main()
