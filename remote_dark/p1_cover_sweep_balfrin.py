# -*- coding: utf-8 -*-
"""Runs ON balfrin — P1 archive sweep of the CHM15k dark plan (report 18): find ACCIDENTAL
covered-telescope periods (maintenance covers) in the 2025-2026 L1 archive.

A covered receiver is unmistakable: the 0.4-3 km band, which always carries boundary-layer
aerosol on an open instrument (clear, cloudy or precipitating alike — fog and precip make it
BRIGHTER, never dimmer), collapses to the dark/noise scale, and the firmware reports no cloud.

The sweep only computes cheap per-day summary statistics; the actual flagging happens locally
on each stream's own distribution (a covered day is a bottom outlier of s_near relative to the
stream's normal range), with the known Payerne hood sessions (2026-05-12, 2026-05-26, ...) as
the positive control.

Per (stream, day): daily median and 10th percentile of the per-profile near-band signal
s_near = p95 of |P| over 0.4-3 km (P = rcs/z^2), the far-band noise scale (12-15 km), the
fraction of profiles with no cloud base, and n_prof. Output: one CSV.
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
D0, D1 = datetime(2015, 1, 1), datetime(2026, 8, 13)   # archive: Payerne 2015-, most 2020-
DECIM = 3                      #: profile decimation for IO (statistics, not spectra)
WORKERS = 64


def l1_path(wmo, ident, day):
    a = L1 / wmo / f"{day:%Y}" / f"{day:%m}" / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
    if a.exists():
        return a
    b = L1 / wmo / f"L1_{wmo}_{ident}{day:%Y%m%d}.nc"
    return b if b.exists() else None


def day_stats(task):
    wmo, ident, ds = task
    import netCDF4
    day = datetime.strptime(ds, "%Y%m%d")
    f = l1_path(wmo, ident, day)
    if f is None:
        return None
    try:
        with netCDF4.Dataset(f) as nc:
            rng = np.ma.filled(nc.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(nc.variables["rcs_0"][::DECIM, :].astype("f4"), np.nan)
            cbh = (np.ma.filled(nc.variables["cloud_base_height"][::DECIM, :]
                                .astype("f4"), np.nan)
                   if "cloud_base_height" in nc.variables else None)
    except Exception:                                                       # noqa: BLE001
        return None
    if rcs.ndim != 2 or rcs.shape[0] < 10 or rcs.shape[1] != rng.size:
        return None
    z = np.where(rng > 0, rng, np.nan)
    P = np.abs(rcs / z[None, :] ** 2)
    near = (rng >= 400) & (rng <= 3000)
    far = (rng >= 12000) & (rng <= 15000)
    with np.errstate(all="ignore"):
        s_near = np.nanpercentile(P[:, near], 95, axis=1)
        s_far = float(np.nanmedian(P[:, far]))
        med = float(np.nanmedian(s_near))
        p10 = float(np.nanpercentile(s_near, 10))
    nocb = (float(np.mean(~np.isfinite(cbh[:, 0]) | (cbh[:, 0] <= 0)))
            if cbh is not None else np.nan)
    return f"{wmo},{ident},{ds},{rcs.shape[0]},{med:.4e},{p10:.4e},{s_far:.4e},{nocb:.3f}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    streams = json.load(open(Path.home() / "chm_streams.json"))
    keys = {(s["wmo"], s["ident"]) for s in streams}
    keys.add(("0-20000-0-06610", "A"))          # Payerne = positive control (hood sessions)
    tasks = []
    for wmo, ident in sorted(keys):
        day = D0
        while day <= D1:
            tasks.append((wmo, ident, f"{day:%Y%m%d}"))
            day += timedelta(days=1)
    print(f"{len(tasks)} stream-days to sweep", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for i, r in enumerate(ex.map(day_stats, tasks, chunksize=16)):
            if r is not None:
                rows.append(r)
            if (i + 1) % 5000 == 0:
                print(f"  {i + 1}/{len(tasks)} scanned, {len(rows)} with data", flush=True)
    with open(OUT / "p1_cover_sweep.csv", "w") as fh:
        fh.write("wmo,ident,date,n_prof,s_near_med,s_near_p10,s_far,frac_no_cbh\n")
        fh.write("\n".join(rows) + "\n")
    print(f"P1_DONE {len(rows)} rows", flush=True)


if __name__ == "__main__":
    main()
