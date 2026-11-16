# -*- coding: utf-8 -*-
"""Runs ON balfrin — P1 sweep v3: find SHORT covered segments (~30 min), the length of a real
operator termination test. The full-day sweep is blind to them (its p10 statistic bottoms out
at ~2.4 h of a day); Payerne's own logbook is mostly short sessions.

Detector, per day: per-profile band MEAN of P over 0.4-3 km (averaging ~170 gates kills the
noise that defeated the p95 statistic) and first-gates mean over 40-400 m; rolling 30-min
medians of both; a covered segment collapses BOTH to ~zero (the first-gates veto is built in:
fog/snow keep the first gates bright). Stored per day: the minimum rolling ratios vs the
day's own level, the below-deck duration, and the hour of the minimum (for verification).
Flagging and threshold calibration happen locally against the short Payerne hood sessions.
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

L1 = Path("/scratch/mch/mhrvo/E-PROFILE_L1")
OUT = Path("/scratch/mch/mhrvo/remote_dark_net")
D0, D1 = datetime(2015, 1, 1), datetime(2026, 8, 13)
DECIM = 3
WORKERS = 64
WIN_S = 1800.0


def l1_path(wmo, ident, day):
    a = L1 / wmo / f"{day:%Y}" / f"{day:%m}" / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
    if a.exists():
        return a
    b = L1 / wmo / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
    return b if b.exists() else None


def rolling_median(x, w):
    if x.size < w:
        return np.array([np.nanmedian(x)])
    out = np.empty(x.size - w + 1)
    for i in range(out.size):
        out[i] = np.nanmedian(x[i:i + w])
    return out


def day_stats(task):
    wmo, ident, ds = task
    import netCDF4
    day = datetime.strptime(ds, "%Y%m%d")
    f = l1_path(wmo, ident, day)
    if f is None:
        return None
    try:
        with netCDF4.Dataset(f) as nc:
            t = np.ma.filled(nc.variables["time"][::DECIM].astype("f8"), np.nan)
            rng = np.ma.filled(nc.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(nc.variables["rcs_0"][::DECIM, :].astype("f4"), np.nan)
    except Exception:                                                       # noqa: BLE001
        return None
    if rcs.ndim != 2 or rcs.shape[0] < 40 or rcs.shape[1] != rng.size:
        return None
    z = np.where(rng > 0, rng, np.nan)
    P = np.where(np.isfinite(rcs) & (np.abs(rcs) < 1e30), rcs, np.nan) / z[None, :] ** 2
    band = (rng >= 400) & (rng <= 3000)
    first = (rng >= 40) & (rng <= 400)
    with np.errstate(all="ignore"):
        bm = np.abs(np.nanmean(P[:, band], axis=1))
        fg = np.abs(np.nanmean(P[:, first], axis=1))
    dt = float(np.nanmedian(np.diff(t))) * 86400.0
    if not np.isfinite(dt) or dt <= 0:
        return None
    w = max(8, int(round(WIN_S / dt)))
    bm30, fg30 = rolling_median(bm, w), rolling_median(fg, w)
    med_bm = float(np.nanmedian(bm))
    med_fg = float(np.nanmedian(fg))
    if not (med_bm > 0 and med_fg > 0):
        return None
    n = min(bm30.size, fg30.size)
    bm30, fg30 = bm30[:n], fg30[:n]
    with np.errstate(all="ignore"):
        joint = np.fmax(bm30 / med_bm, fg30 / med_fg)      # BOTH must be low
    if not np.isfinite(joint).any():
        return None
    imin = int(np.nanargmin(joint))
    hour = float((datetime(1970, 1, 1) + timedelta(days=float(t[imin]))).hour
                 + (datetime(1970, 1, 1) + timedelta(days=float(t[imin]))).minute / 60.0)
    low = (bm / med_bm < 0.1) & (fg / med_fg < 0.15)
    dur_min = float(low.sum() * dt / 60.0)
    return (f"{wmo},{ident},{ds},{rcs.shape[0]},{np.nanmin(joint):.4f},"
            f"{bm30[imin] / med_bm:.4f},{fg30[imin] / med_fg:.4f},{hour:.2f},{dur_min:.0f}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    streams = json.load(open(Path.home() / "chm_streams.json"))
    keys = {(s["wmo"], s["ident"]) for s in streams}
    keys.add(("0-20000-0-06610", "A"))
    tasks = []
    for wmo, ident in sorted(keys):
        day = D0
        while day <= D1:
            tasks.append((wmo, ident, f"{day:%Y%m%d}"))
            day += timedelta(days=1)
    print(f"{len(tasks)} stream-days (short-segment sweep)", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for i, r in enumerate(ex.map(day_stats, tasks, chunksize=16)):
            if r is not None:
                rows.append(r)
            if (i + 1) % 20000 == 0:
                print(f"  {i + 1}/{len(tasks)}, {len(rows)} with data", flush=True)
    with open(OUT / "p1_short_sweep.csv", "w") as fh:
        fh.write("wmo,ident,date,n_prof,min_joint30,min_bm30,min_fg30,hour_min,dur_min\n")
        fh.write("\n".join(rows) + "\n")
    print(f"P1C_DONE {len(rows)} rows", flush=True)


if __name__ == "__main__":
    main()
