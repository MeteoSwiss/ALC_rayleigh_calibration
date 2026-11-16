# -*- coding: utf-8 -*-
"""M3 — twilight background sweep: the background-driven dark component, from sky data alone.

Idea. At the terminator the optical background sweeps orders of magnitude in tens of minutes while
the atmosphere and the instrument's thermal state barely move. Per gate, regressing the signal on
the measured background B(t) over such a window isolates the background-COUPLED part of the
baseline, D_bg(z) — the firmware background-compensation residual that clear-night methods never
see (they select for darkness) and that pollutes every daytime product.

Identifiability, learned from the hood probe (probe_hood_dynamics):
  * B and internal temperature are collinear on DIURNAL timescales — but not at the terminator,
    where B moves ~10x faster than the thermal mass. The per-event model therefore carries the
    slow drift explicitly (τ, τ²) and claims only the fast coupling for B.
  * The dark's modulated response is near-range concentrated (|corr| up to 0.7 at 60–500 m under
    the hood) — D_bg(z) is expected to live exactly where the static CL31 dark lives.
  * The response is UNIT-dependent (the post-swap CL31 optic responds ~2x more), so D_bg is a
    per-unit product, never a type constant.

Per event e (dawn or dusk), all gates at once, two IRLS passes with 4-MAD clipping:

    rcs(t, z) = c0(z) + D_e(z)·B̃(t) + c1(z)·τ + c2(z)·τ²     B̃ = B centred/scaled, τ = hours

First-run lesson (kept, because it is the finding): with the RAW background as regressor the sky
estimate anticorrelates with the hood truth (corr ≈ −0.2). At the terminator the background is
locked to the boundary-layer evolution itself — dawn convection rises with the sun — so the
regression attributes real atmospheric change to B. Identification must come from SPEED, not from
the sweep: the default regressor is therefore the HIGH-PASSED background (B minus a ~45-min
running median): cloud-shadow flickers and fast ramps that the boundary layer cannot follow, but
the firmware's background compensation follows instantly.

and the campaign product is the across-event median D(z) with MAD, dawn/dusk separated (a real
instrumental coupling must agree in sign and magnitude between the two terminators; residual
atmospheric drift does not).

Validation (the point of doing this at Payerne first): the two multi-hour hood sessions contain
their own dawn/dusk, with no atmosphere in front — running the SAME estimator on them yields the
true D_bg(z). The sky estimate is compared to it blind.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime

import numpy as np

from remote_dark.common import (PAYERNE, days_between, ensure_out, read_day,
                                solar_elevation_deg)
from remote_dark import hood

#: Twilight = solar elevation inside this band (deg). Below −12° the sky is dark (no sweep);
#: above +25° the background saturates its diurnal plateau and clouds dominate the variance.
EL_BAND = (-12.0, 25.0)
MIN_PROF = 40          #: fewer profiles than this cannot constrain 4 coefficients robustly
MIN_BSPAN = 5.0        #: B must actually sweep (its MAD-normalised span) or the event is void


def _events_from_day(d, lat, lon, exclude=()):
    """Split one day into dawn/dusk twilight events; honour exclusion windows (hood periods)."""
    el = solar_elevation_deg(lat, lon, d["times"])
    rising = np.gradient(el) > 0
    out = []
    for name, m in (("dawn", (el >= EL_BAND[0]) & (el <= EL_BAND[1]) & rising),
                    ("dusk", (el >= EL_BAND[0]) & (el <= EL_BAND[1]) & ~rising)):
        if m.sum() < MIN_PROF:
            continue
        tt = d["times"][m]
        if any((tt[0] <= t1 and tt[-1] >= t0) for t0, t1 in exclude):
            continue
        out.append({"kind": name, "times": tt, "rcs": d["rcs"][m],
                    "B": d["hk"]["bckgrd"][m], "t_int": d["hk"]["t_int"][m]})
    return out


HP_MIN = 45.0          #: high-pass window (minutes) for the default regressor


def _highpass(B, times, minutes=HP_MIN):
    """B minus its running median over *minutes* — the fast component only."""
    t = np.array([(x - times[0]).total_seconds() / 60.0 for x in times])
    out = np.full_like(B, np.nan, dtype=float)
    half = minutes / 2.0
    for i in range(B.size):
        m = (t >= t[i] - half) & (t <= t[i] + half) & np.isfinite(B)
        if m.sum() >= 5:
            out[i] = B[i] - np.median(B[m])
    return out


def _fit_event(ev, rng, regressor="hp"):
    """One event -> D_e(z) + diagnostics, or None. Two IRLS passes, 4-MAD residual clipping."""
    B = ev["B"].astype(float)
    if not np.isfinite(B).any():
        # No background variable (CL61): far-gate proxy — the raw-signal mean over the top 8 % of
        # gates, which is background-dominated on every type.
        top = rng >= np.nanpercentile(rng, 92)
        with np.errstate(invalid="ignore"):
            B = np.nanmean(ev["rcs"][:, top] / rng[top] ** 2, axis=1)
    if regressor == "hp":
        B = _highpass(B, ev["times"])
    good = np.isfinite(B)
    if good.sum() < MIN_PROF:
        return None
    bmad = np.nanmedian(np.abs(B - np.nanmedian(B))) + 1e-30
    if (np.nanmax(B) - np.nanmin(B)) / bmad < MIN_BSPAN and np.nanstd(B) < 0.5:
        return None
    t0 = ev["times"][good][0]
    tau = np.array([(x - t0).total_seconds() / 3600.0 for x in ev["times"][good]])
    tau -= tau.mean()
    bs = float(np.nanstd(B[good])) + 1e-30
    Bn = (B[good] - np.nanmean(B[good])) / bs
    X = np.column_stack([np.ones_like(tau), Bn, tau, tau ** 2])
    Y = ev["rcs"][good]
    W = np.isfinite(Y)
    Yf = np.where(W, Y, 0.0)

    beta = None
    for _ in range(2):                                     # IRLS: plain, then clipped
        XtX = np.einsum("ti,tj,tz->zij", X, X, W.astype(float))
        XtY = np.einsum("ti,tz->zi", X, Yf * W)
        beta = np.full((rng.size, 4), np.nan)
        okz = np.linalg.matrix_rank(XtX.sum(axis=0)) == 4   # cheap global sanity
        for z in range(rng.size):
            if W[:, z].sum() < MIN_PROF:
                continue
            try:
                beta[z] = np.linalg.solve(XtX[z] + 1e-12 * np.eye(4), XtY[z])
            except np.linalg.LinAlgError:
                continue
        resid = Y - X @ beta.T
        mad = np.nanmedian(np.abs(resid - np.nanmedian(resid, axis=0)), axis=0) + 1e-30
        W = np.isfinite(Y) & (np.abs(resid) < 4.0 * mad)
        Yf = np.where(W, Y, 0.0)

    # collinearity of the fast regressor with the slow drift — the identifiability diagnostic
    c_bt = float(np.corrcoef(Bn, tau)[0, 1])
    return {"D": beta[:, 1] / bs, "kind": ev["kind"], "n_prof": int(good.sum()),
            "b_span": float(np.nanmax(B) - np.nanmin(B)), "corr_B_tau": c_bt,
            "t0": ev["times"][good][0]}


def run_sky(wmo: str, ident: str, d0: datetime, d1: datetime, lat: float, lon: float,
            exclude=(), regressor="hp"):
    """All twilight events of a period -> aggregated D_bg(z)."""
    evs, rng = [], None
    for day in days_between(d0, d1):
        d = read_day(wmo, ident, day)
        if d is None:
            continue
        rng = d["rng"]
        for ev in _events_from_day(d, lat, lon, exclude):
            r = _fit_event(ev, rng, regressor)
            if r is not None:
                evs.append(r)
    return _aggregate(evs, rng)


def run_hood(ident: str, regressor="hp"):
    """The SAME estimator on the multi-hour hood sessions -> the true D_bg(z), per era."""
    out = {}
    for s in hood.session_frames(ident, min_hours=10):
        lat, lon = PAYERNE["lat"], PAYERNE["lon"]
        d = {"times": s["times"], "rcs": s["dark"],
             "hk": {"bckgrd": s["hk"]["bckgrd"], "t_int": s["hk"]["t_int"]}}
        evs = []
        for ev in _events_from_day(d, lat, lon):
            r = _fit_event(ev, s["rng"], regressor)
            if r is not None:
                evs.append(r)
        agg = _aggregate(evs, s["rng"])
        if agg is not None:
            out[s["era"] + "_" + f"{s['t0']:%Y%m%d}"] = agg
    return out


def _aggregate(evs, rng):
    if not evs or rng is None:
        return None
    D = np.vstack([e["D"] for e in evs])
    dawn = np.vstack([e["D"] for e in evs if e["kind"] == "dawn"]) if any(
        e["kind"] == "dawn" for e in evs) else None
    dusk = np.vstack([e["D"] for e in evs if e["kind"] == "dusk"]) if any(
        e["kind"] == "dusk" for e in evs) else None
    return {
        "rng": rng,
        "D_med": np.nanmedian(D, axis=0),
        "D_mad": np.nanmedian(np.abs(D - np.nanmedian(D, axis=0)), axis=0),
        "n_events": len(evs),
        "D_dawn": (np.nanmedian(dawn, axis=0) if dawn is not None else np.full(rng.size, np.nan)),
        "D_dusk": (np.nanmedian(dusk, axis=0) if dusk is not None else np.full(rng.size, np.nan)),
        "n_dawn": 0 if dawn is None else dawn.shape[0],
        "n_dusk": 0 if dusk is None else dusk.shape[0],
        "corr_B_tau_med": float(np.nanmedian([e["corr_B_tau"] for e in evs])),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ident", default="B")
    ap.add_argument("--start", default="20260401")
    ap.add_argument("--end", default="20260813")
    ap.add_argument("--regressor", default="hp", choices=["hp", "raw"])
    a = ap.parse_args()
    wmo, lat, lon = PAYERNE["wmo"], PAYERNE["lat"], PAYERNE["lon"]
    ident = a.ident
    d0 = datetime.strptime(a.start, "%Y%m%d")
    d1 = datetime.strptime(a.end, "%Y%m%d")

    # hood periods must never count as sky
    excl = [(hood._parse(t0), hood._parse(t1)) for t0, t1 in hood.WINDOWS.get(ident, [])]

    print(f"== sky sweep {wmo}_{ident} {a.start}..{a.end}")
    sky = run_sky(wmo, ident, d0, d1, lat, lon, exclude=excl, regressor=a.regressor)
    if sky is None:
        print("   no usable events")
        return
    print(f"   events: {sky['n_events']} (dawn {sky['n_dawn']} / dusk {sky['n_dusk']}), "
          f"median corr(B, tau) = {sky['corr_B_tau_med']:.2f}")

    print(f"== hood truth {ident}")
    hd = run_hood(ident, regressor=a.regressor)
    for k, v in hd.items():
        print(f"   {k}: {v['n_events']} events (dawn {v['n_dawn']} / dusk {v['n_dusk']})")

    out = ensure_out("twilight")
    payload = {"rng": sky["rng"], "sky_D_med": sky["D_med"], "sky_D_mad": sky["D_mad"],
               "sky_D_dawn": sky["D_dawn"], "sky_D_dusk": sky["D_dusk"],
               "sky_n_events": sky["n_events"]}
    for k, v in hd.items():
        payload[f"hood_{k}_D_med"] = v["D_med"]
        payload[f"hood_{k}_n"] = v["n_events"]
    np.savez(out / f"{wmo}_{ident}.npz", **payload)

    # the comparison that decides M3: sky vs hood, near range, in P-view
    rng = sky["rng"]
    for k, v in hd.items():
        m = (rng >= 60) & (rng <= 1500) & np.isfinite(sky["D_med"]) & np.isfinite(v["D_med"])
        if m.sum() < 10:
            continue
        num = float(np.nansum(sky["D_med"][m] * v["D_med"][m]))
        den = float(np.nansum(v["D_med"][m] ** 2)) + 1e-30
        gain = num / den
        cc = float(np.corrcoef(sky["D_med"][m], v["D_med"][m])[0, 1])
        print(f"   sky vs {k} (60-1500 m): corr {cc:+.2f}   amplitude gain {gain:+.2f} "
              f"(1 = perfect)")
    print(f"   npz -> {out / f'{wmo}_{ident}.npz'}")


if __name__ == "__main__":
    main()
