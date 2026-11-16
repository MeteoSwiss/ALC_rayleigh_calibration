# -*- coding: utf-8 -*-
"""Is the measured Payerne "dark" a residual of the far-field background subtraction?

The L1 rcs_0 is background-subtracted and range-corrected. A small error c in the subtracted background
(a constant in the NON-range-corrected signal P = rcs_0 / r^2) appears in rcs_0 as c * r^2, i.e. a dark
that grows quadratically with range -- the shape seen under the hood above ~1.5 km.

Check 1 (hood sessions): P_dark(r) = rcs_0/r^2, time median per session. If P_dark is flat beyond a few
km, the dark is a background offset c; fit rcs_dark = c r^2 over 2-15 km and report how much of the
dark it explains.
Check 2 (normal operation): P(r) at the far range (12-15 km) on clear-sky night hours of the comparison
period. Without any atmospheric signal left there, P should average to zero; a systematic offset there
means the operational background estimate is biased, and that bias can be estimated from the far field
in every profile (no hood needed).

Run:  python rayleigh_availability/dark_background_check.py
"""
from __future__ import annotations
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
OUT_FIG = REPO / "doc" / "reports" / "figs_remote_dark" / "dark_background_check.png"
OUT_JSON = DATA / "dark_background_check.json"
NAMES = {"A": "CHM15k", "B": "CL31", "C": "CL61"}
FAR = {"A": (12000.0, 15000.0), "B": (6000.0, 7600.0), "C": (12000.0, 15000.0)}   # far-field windows


def read_day(ident, day, with_cbh=False):
    f = L1 / f"{day:%Y}" / f"{day:%m}" / f"L1_0-20000-0-06610_{ident}{day:%Y%m%d}.nc"
    if not f.exists():
        return None
    with netCDF4.Dataset(f) as ds:
        t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
        r = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
        x = np.ma.filled(ds["rcs_0"][:].astype("f8"), np.nan)
        cbh = np.ma.filled(ds["cloud_base_height"][:].astype("f8"), np.nan) if with_cbh else None
    if x.shape[0] != t.size:
        x = x.T
    tt = np.array([datetime(1970, 1, 1) + timedelta(days=float(v)) for v in t])
    return tt, r, x, cbh


def p_profile(r, x):
    """Non-range-corrected signal P = rcs_0 / r^2 (r in km), NaN at r <= 0."""
    rk = r / 1000.0
    with np.errstate(all="ignore"):
        return np.where(rk[None, :] > 0, x / rk[None, :] ** 2, np.nan)


def binned(r, y, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    v = np.array([np.nanmedian(y[(r >= a) & (r < b)]) if ((r >= a) & (r < b)).any() else np.nan
                  for a, b in zip(edges[:-1], edges[1:])])
    return c, v


def hood_sessions():
    win = json.loads((DATA / "dark_windows_detected.json").read_text())
    edges = np.arange(250.0, 15500.0, 250.0)
    out = {}
    for ident in ("A", "B", "C"):
        out[ident] = []
        for w0, w1 in win[ident]:
            t0 = datetime.strptime(w0, "%Y-%m-%d %H:%M:%S")
            t1 = datetime.strptime(w1, "%Y-%m-%d %H:%M:%S")
            blocks, r = [], None
            day = t0.date()
            while day <= t1.date():
                got = read_day(ident, datetime.combine(day, datetime.min.time()))
                if got:
                    tt, r, x, _ = got
                    m = (tt >= t0) & (tt <= t1)
                    if m.any():
                        blocks.append(x[m])
                day += timedelta(days=1)
            X = np.vstack(blocks)
            rcs_med = np.nanmedian(X, axis=0)
            P = np.nanmean(p_profile(r, X), axis=0)                 # time MEAN of P (offset estimate)
            rc, Pb = binned(r, P, edges)
            # fit rcs_dark = c * r^2 (r in km) over 2-15 km (CL31: 2-7.6 km)
            fit = (r >= 2000) & (r <= (7600 if ident == "B" else 15000)) & np.isfinite(rcs_med)
            rk = r[fit] / 1000.0
            c = float(np.sum(rcs_med[fit] * rk ** 2) / np.sum(rk ** 4))
            resid = rcs_med[fit] - c * rk ** 2
            frac = 1.0 - float(np.sum(resid ** 2) / np.sum(rcs_med[fit] ** 2))
            far = (r >= FAR[ident][0]) & (r <= FAR[ident][1])
            out[ident].append(dict(start=w0[:10], c=c, explained=frac, P_far=float(np.nanmean(P[far])),
                                   P_3_6=float(np.nanmean(P[(r >= 3000) & (r <= 6000)])),
                                   rc=rc, Pb=Pb))
    return out


def clear_nights(ident, d0, d1):
    """Far-field P on clear-sky night hours (21-03 UTC, no cloud base reported in the hour)."""
    rows = []
    day = d0
    while day <= d1:
        got = read_day(ident, day, with_cbh=True)
        if got:
            tt, r, x, cbh = got
            h = np.array([v.hour for v in tt])
            night = (h >= 21) | (h < 3)
            cb = cbh if cbh.ndim == 2 else cbh[:, None]
            clear = ~np.any(np.isfinite(cb) & (cb > 0), axis=1)
            m = night & clear
            if m.sum() > 60:
                P = p_profile(r, x[m])
                far = (r >= FAR[ident][0]) & (r <= FAR[ident][1])
                mid = (r >= 3000) & (r <= 6000)
                rows.append((day, float(np.nanmean(P[:, far])), float(np.nanmean(P[:, mid])), int(m.sum())))
        day += timedelta(days=1)
    return rows


def main():
    hood = hood_sessions()
    summary = {}
    for ident in ("A", "B", "C"):
        print(f"== {NAMES[ident]}: hood sessions  (P = rcs_0/r^2, r in km)")
        for s in hood[ident]:
            print(f"   {s['start']}  c (fit rcs=c r^2) = {s['c']:+.3e}   explained {100 * s['explained']:5.1f} %"
                  f"   <P> 3-6 km = {s['P_3_6']:+.3e}   <P> far = {s['P_far']:+.3e}")
        summary[ident] = dict(hood_c=[s["c"] for s in hood[ident]],
                              hood_explained=[s["explained"] for s in hood[ident]],
                              hood_P_far=[s["P_far"] for s in hood[ident]])
    night = {}
    for ident in ("A", "C"):
        rows = clear_nights(ident, datetime(2026, 3, 4), datetime(2026, 10, 8))
        night[ident] = rows
        pf = np.array([r[1] for r in rows])
        pm = np.array([r[2] for r in rows])
        print(f"== {NAMES[ident]}: clear-sky night hours, {len(rows)} nights: <P> far = "
              f"{np.nanmedian(pf):+.3e} (IQR {np.nanpercentile(pf, 25):+.3e} .. {np.nanpercentile(pf, 75):+.3e})")
        summary[ident].update(night_P_far_median=float(np.nanmedian(pf)),
                              night_P_far_q25=float(np.nanpercentile(pf, 25)),
                              night_P_far_q75=float(np.nanpercentile(pf, 75)), n_nights=len(rows),
                              hood_c_median=float(np.median(summary[ident]["hood_c"])))
    OUT_JSON.write_text(json.dumps(summary, indent=1, default=str))
    figure(hood, night)


def figure(hood, night):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    cmap = plt.get_cmap("viridis")
    for k, ident in enumerate(("A", "B", "C")):
        ax = axs[0, k]
        ss = hood[ident]
        for i, s in enumerate(ss):
            ax.plot(s["Pb"], s["rc"] / 1000, color=cmap(i / max(len(ss) - 1, 1)), lw=1.2, label=s["start"])
            ax.axvline(s["c"], color=cmap(i / max(len(ss) - 1, 1)), lw=0.6, ls=":")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_ylabel("Range [km]")
        ax.set_xlabel("P = rcs$_0$ / r$^2$  (r in km, L1 units)")
        ax.set_title(f"{NAMES[ident]}: hood dark, NON-range-corrected (dotted: fitted constant c)", fontsize=10)
        lo = np.nanpercentile(np.concatenate([s["Pb"] for s in ss]), 2)
        hi = np.nanpercentile(np.concatenate([s["Pb"] for s in ss]), 98)
        ax.set_xlim(lo - 0.3 * abs(hi - lo), hi + 0.3 * abs(hi - lo))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    for k, ident in enumerate(("A", "C")):
        ax = axs[1, 2 * k] if k == 0 else axs[1, 2]
        rows = night[ident]
        t = [r[0] for r in rows]
        ax.plot(t, [r[1] for r in rows], "o", ms=3, label=f"clear night, {FAR[ident][0]/1e3:.0f}-{FAR[ident][1]/1e3:.0f} km")
        hc = [s["c"] for s in hood[ident]]
        ht = [datetime.strptime(s["start"], "%Y-%m-%d") for s in hood[ident]]
        ax.plot(ht, hc, "s", color="r", ms=6, label="hood sessions (fitted c)")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"{NAMES[ident]}: far-field P in normal operation vs hood offset", fontsize=10)
        ax.set_ylabel("P [L1 units]")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axs[1, 1].axis("off")
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=110)
    print(f"-> {OUT_FIG}")


if __name__ == "__main__":
    main()
