# -*- coding: utf-8 -*-
"""q15d step 1: extract CAMS columns at every CHM15k station (report 19, attribution question q15, agent d).

For each CHM15k stream of the census, the nearest CAMS 0.4 deg grid cell is sampled for every time step of
the monthly files D:/CAMS_Monthly_04 (2025-01 .. 2026-06) and of the daily files D:/CAMS_daily
(2026-07 .. 2026-08): aerosol backscatter at 1064 nm (attenuated, from the ground), temperature, model-level
altitude (m ASL) and full-level pressure (L137 a/b coefficients x surface pressure). Only the lowest
61 model levels (surface .. ~11 km) are kept, which covers every Rayleigh window (zt < 8 km).

Reading the whole field per time block and indexing the station points is ~100x cheaper than slicing one
grid point at a time (chunks are 62 x 26 x 34 x 47).

Output: C:/DATA/Projects/202610_CHM15k_stability/q15/q15d_cams/cams_<tag>.npz
Run:    python validation/chm15k_stability/q15d_cams_extract.py [--workers 4]
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import OUT, chm15k_streams, stream_key  # noqa: E402
from calibration.water_vapor_correction.water_vapor import _A137, _B137  # noqa: E402

DST = OUT / "q15" / "q15d_cams"
LEV0 = 40            # level index 40 = model level 77 (~11 km); index 100 = level 137 (surface)
VARS = ("aerbackscatgnd1064", "t", "altitude")


def _points(lat_grid, lon_grid):
    keys, ii, jj = [], [], []
    for wmo, ident, site, lat, lon in chm15k_streams():
        keys.append(stream_key(wmo, ident))
        ii.append(int(np.argmin(np.abs(lat_grid - lat))))
        jj.append(int(np.argmin(np.abs(lon_grid - lon))))
    return keys, np.array(ii), np.array(jj)


def extract(path: str):
    import netCDF4
    path = Path(path)
    out = DST / f"cams_{path.stem.replace('CAMS_Beta_', '')}.npz"
    if out.exists():
        return f"skip {out.name}"
    with netCDF4.Dataset(path) as ds:
        lat, lon = ds["latitude"][:], ds["longitude"][:]
        keys, ii, jj = _points(np.asarray(lat), np.asarray(lon))
        level = np.asarray(ds["level"][:], int)[LEV0:]
        tv = ds["time"]
        time = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
        tdays = np.array([(np.datetime64(x) - np.datetime64("1970-01-01")) / np.timedelta64(1, "D") for x in time])
        nt = len(tdays)
        tchunk = ds["t"].chunking()[0]
        res = {v: np.full((len(keys), nt, len(level)), np.nan, np.float32) for v in VARS}
        lnsp = np.full((len(keys), nt), np.nan)
        # surface pressure is stored in one finite level slot
        probe = np.ma.filled(ds["lnsp"][0, :, ii[0], jj[0]], np.nan)
        slot = int(np.where(np.isfinite(probe))[0][0])
        for t0 in range(0, nt, tchunk):
            t1 = min(nt, t0 + tchunk)
            for v in VARS:
                blk = np.ma.filled(ds[v][t0:t1, LEV0:, :, :].astype(np.float32), np.nan)
                res[v][:, t0:t1, :] = np.moveaxis(blk[:, :, ii, jj], 2, 0)
            lb = np.ma.filled(ds["lnsp"][t0:t1, slot, :, :], np.nan)
            lnsp[:, t0:t1] = lb[:, ii, jj].T
    ps = np.exp(lnsp)                                               # (n_pt, n_t) Pa
    ph_lo = _A137[level][None, None, :] + _B137[level][None, None, :] * ps[:, :, None]
    ph_hi = _A137[level - 1][None, None, :] + _B137[level - 1][None, None, :] * ps[:, :, None]
    p = (0.5 * (ph_lo + ph_hi)).astype(np.float32)
    DST.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, keys=np.array(keys), time=tdays, level=level, aer1064=res["aerbackscatgnd1064"],
                        T=res["t"], alt=res["altitude"], p=p, lat=np.asarray(lat)[ii], lon=np.asarray(lon)[jj])
    return f"done {out.name} ({nt} steps)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    files = sorted(Path("D:/CAMS_Monthly_04").glob("CAMS_Beta_2*.nc"))
    daily = sorted(Path("D:/CAMS_daily").glob("CAMS_Beta_2026*.nc"))
    files += [f for f in daily if f.stem[10:16] in ("202607", "202608") and len(f.stem) == 18]
    with ProcessPoolExecutor(a.workers) as ex:
        for msg in ex.map(extract, [str(f) for f in files]):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
