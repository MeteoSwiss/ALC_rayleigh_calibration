# -*- coding: utf-8 -*-
"""Build the 10-min housekeeping + binned-signal cache for every CHM15k stream (report 19, step 2).

For each daily L1 file it keeps, per 10-min block:
  * the median of the housekeeping variables (temperatures in degC, calibration pulse, window
    transmission, laser/detector status, laser life time, background, noise),
  * the number of profiles and of CLEAR profiles (no cloud base reported, sky condition index 0),
  * the mean rcs_0 of the clear profiles in coarse range bins (100 m to 1 km, 250 m to 8 km).

Only clear profiles enter the signal mean, so the cache can serve both the intra-night stability
test (signal in the Rayleigh window) and the same-wavelength pair ratios by day and by night. The
optical_module_id global attribute is stored per day to segment the record at module swaps.

Run:  python validation/chm15k_stability/hk_cache.py --workers 20
      python validation/chm15k_stability/hk_cache.py --streams 0-20000-0-06240_A --force
"""
from __future__ import annotations

import argparse
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BLOCKS_PER_DAY, CACHE, HK_VARS, L1_ROOT, PAIR_SITES, RANGE_EDGES, TEMP_VARS, YEARS,  # noqa: E402
                    cache_file, chm15k_streams, stream_key)

N_BINS = len(RANGE_EDGES) - 1


def _read_var(ds, name, n):
    """Variable as float64 with masked -> NaN; an all-NaN column when the stream lacks it."""
    if name not in ds.variables:
        return np.full(n, np.nan)
    return np.ma.filled(ds[name][:].astype("f8"), np.nan)


def _block_median(values, block, valid_blocks):
    """Median of `values` per 10-min block (NaN-aware); NaN where a block has no data."""
    out = np.full(BLOCKS_PER_DAY, np.nan)
    for b in valid_blocks:
        v = values[block == b]
        v = v[np.isfinite(v)]
        if v.size:
            out[b] = np.median(v)
    return out


def process_day(path):
    """Reduce one L1 day to (hk[144, nvar], n_all[144], n_clear[144], rcs[144, nbin], module, tilt)."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        t = np.ma.filled(ds["time"][:].astype("f8"), np.nan)
        n = t.size
        rng = np.ma.filled(ds["range"][:].astype("f8"), np.nan)
        # scalar on CHM15k, per-profile on CL61
        tilt = float(np.nanmedian(np.ma.filled(ds["tilt_angle"][...].astype("f8"), np.nan))) \
            if "tilt_angle" in ds.variables else np.nan
        module = str(ds.getncattr("optical_module_id")) if "optical_module_id" in ds.ncattrs() else ""
        hk = {v: _read_var(ds, v, n) for v in HK_VARS}
        cbh = np.ma.filled(ds["cloud_base_height"][:, 0].astype("f8"), np.nan) if "cloud_base_height" in ds.variables \
            else np.full(n, np.nan)
        sci = _read_var(ds, "sci", n)
        igate = np.searchsorted(rng, RANGE_EDGES[-1])
        rcs = np.ma.filled(ds["rcs_0"][:, :igate].astype("f4"), np.nan)
    for v in TEMP_VARS:
        hk[v] = hk[v] - 273.15
    day0 = np.floor(np.nanmedian(t))
    block = np.clip(((t - day0) * BLOCKS_PER_DAY).astype(int), 0, BLOCKS_PER_DAY - 1)
    ok_t = np.isfinite(t)
    valid_blocks = np.unique(block[ok_t])

    hk_out = np.stack([_block_median(hk[v], block, valid_blocks) for v in HK_VARS], axis=1)
    # Clear = no cloud base reported (fill/negative/NaN) and no precipitation/fog flagged by sci.
    clear = ok_t & ~(np.isfinite(cbh) & (cbh > 0)) & ((sci == 0) | ~np.isfinite(sci))
    n_all = np.bincount(block[ok_t], minlength=BLOCKS_PER_DAY).astype("i2")
    n_clear = np.bincount(block[clear], minlength=BLOCKS_PER_DAY).astype("i2")

    # Range-bin means per profile (NaN-aware), then mean over the clear profiles of each block.
    gbin = np.digitize(rng[:igate], RANGE_EDGES) - 1
    rcs_out = np.full((BLOCKS_PER_DAY, N_BINS), np.nan, dtype="f4")
    if clear.any():
        x = rcs[clear]
        fin = np.isfinite(x)
        x0 = np.where(fin, x, 0.0).astype("f8")
        prof = np.full((x.shape[0], N_BINS), np.nan)
        for k in range(N_BINS):
            sel = gbin == k
            if not sel.any():
                continue
            cnt = fin[:, sel].sum(axis=1)
            s = x0[:, sel].sum(axis=1)
            prof[:, k] = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
        bclear = block[clear]
        for b in np.unique(bclear):
            rcs_out[b] = np.nanmean(prof[bclear == b], axis=0)
    return day0, hk_out.astype("f4"), n_all, n_clear, rcs_out, module, tilt


def day_files(wmo, ident):
    """All daily L1 files of the stream in YEARS, in date order."""
    out = []
    d = date(YEARS[0], 1, 1)
    end = date(YEARS[-1], 12, 31)
    while d <= end:
        f = L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
        if f.exists():
            out.append(f)
        d += timedelta(days=1)
    return out


def build_stream(wmo, ident, force=False):
    out = cache_file(wmo, ident)
    if out.exists() and not force:
        return stream_key(wmo, ident), "exists", 0
    files = day_files(wmo, ident)
    if not files:
        return stream_key(wmo, ident), "no L1", 0
    rows = []
    errors = 0
    for f in files:
        try:
            rows.append(process_day(f))
        except Exception:  # a corrupt day must not lose the stream; count it and move on
            errors += 1
    if not rows:
        return stream_key(wmo, ident), f"all {errors} days failed", 0
    tilt = np.nanmedian([r[6] for r in rows])
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.npz")
    np.savez_compressed(
        tmp,
        day=np.array([r[0] for r in rows]),
        hk=np.stack([r[1] for r in rows]),
        n_all=np.stack([r[2] for r in rows]),
        n_clear=np.stack([r[3] for r in rows]),
        rcs=np.stack([r[4] for r in rows]),
        module=np.array([r[5] for r in rows], dtype="U32"),
        tilt=np.float64(tilt),
        hk_vars=np.array(HK_VARS, dtype="U32"),
        range_edges=RANGE_EDGES,
    )
    tmp.replace(out)
    return stream_key(wmo, ident), f"ok ({errors} bad days)", len(rows)


def _worker(args):
    wmo, ident, force = args
    try:
        return build_stream(wmo, ident, force)
    except Exception:
        return stream_key(wmo, ident), "FAILED " + traceback.format_exc(limit=2).replace("\n", " | "), 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--streams", nargs="*", help="restrict to these <wmo>_<ident> keys")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    streams = chm15k_streams()
    if a.streams:
        streams = [s for s in streams if stream_key(s[0], s[1]) in a.streams]
    CACHE.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    # Key streams first: the same-wavelength pairs and the CHM15k units co-located with a CL61.
    priority = {f"{w}_{u}" for w, units in PAIR_SITES.items() for u in units}
    priority |= {"0-20000-0-06610_A", "0-20000-0-10393_0", "0-20000-0-03808_A", "0-380-5-1_0"}
    streams = sorted(streams, key=lambda s: stream_key(s[0], s[1]) not in priority)
    jobs = [(s[0], s[1], a.force) for s in streams]
    print(f"{len(jobs)} CHM15k streams, {a.workers} workers -> {CACHE}", flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            key, status, ndays = fu.result()
            print(f"[{i:3d}/{len(jobs)}] {time.time() - t0:7.0f}s  {key:28s} {ndays:4d} d  {status}", flush=True)
    print("HK_CACHE_DONE", flush=True)


if __name__ == "__main__":
    main()
