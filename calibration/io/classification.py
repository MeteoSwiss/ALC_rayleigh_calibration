# -*- coding: utf-8 -*-
"""Read Cloudnet target-classification curtains (the ceiloclass product) back in.

The classification is written per stream-DAY, while a Rayleigh night straddles two of them
(evening of d-1 + morning of d). :func:`read_classification_curtain` stitches the days a night can
touch into one (time, range, code) curtain; the consumer selects the night itself, on the profile
times the fit actually used -- which is what keeps a "fraction of the night" from being diluted by
the daylight hours in the same files.

Reading the product needs nothing but netCDF4: ``ceiloclass``/``ceilopyter`` are required to
PRODUCE a classification, never to consume one.
"""
from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np

__all__ = ["read_classification_curtain", "classification_files_for_night"]


def classification_files_for_night(root: Path, key: str, wmo: str, date_str: str) -> list[Path]:
    """The classification NetCDFs a Rayleigh night for ``date_str`` can draw on (d-1 and d)."""
    d = datetime.strptime(date_str, "%Y%m%d")
    out = []
    for dd in (d - timedelta(days=1), d):
        ds = dd.strftime("%Y%m%d")
        f = Path(root) / key / "classification" / wmo / ds[:4] / f"{key}_{ds}_classification.nc"
        if f.is_file():
            out.append(f)
    return out


def read_classification_curtain(
    files: Iterable[Path],
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Stitch classification files into ``(times, range_agl, target_classification)``.

    ``times`` are ``datetime64[s]``; ``range_agl`` is the range grid of the FIRST file read, with
    later files interpolated (nearest) onto it if their grid differs; codes follow the Cloudnet
    convention (0 clear, 1 droplet, 2 drizzle/rain, 3 ice, 4 supercooled, 5 aerosol, ...).

    Returns None when no file could be read -- the caller's screening is then a no-op.
    """
    import netCDF4

    times: list[np.ndarray] = []
    codes: list[np.ndarray] = []
    grid: Optional[np.ndarray] = None
    for f in files:
        try:
            with netCDF4.Dataset(str(f)) as nc:
                rng = np.asarray(nc.variables["range"][:], dtype=float)          # AGL
                tc = np.asarray(nc.variables["target_classification"][:])        # (time, range)
                tv = nc.variables["time"]
                t = _decode_time(np.asarray(tv[:], dtype=float), getattr(tv, "units", ""), f.name)
        except (OSError, KeyError, IndexError):
            continue
        if tc.ndim != 2 or tc.shape != (t.size, rng.size):
            continue
        if grid is None:
            grid = rng
        elif rng.shape != grid.shape or not np.allclose(rng, grid):
            j = np.clip(np.searchsorted(rng, grid), 0, rng.size - 1)
            tc = tc[:, j]
        times.append(t)
        codes.append(tc)
    if grid is None or not times:
        return None
    t_all = np.concatenate(times)
    c_all = np.concatenate(codes, axis=0)
    order = np.argsort(t_all)
    return t_all[order], grid, c_all[order]


def _decode_time(vals: np.ndarray, units: str, fname: str) -> np.ndarray:
    """Decode the time axis to datetime64[s].

    Cloudnet writes ``time`` as decimal hours since midnight UTC of the file's day, but the units
    attribute is what actually says so -- honour it when it carries an epoch, and fall back to the
    filename date otherwise (the ``<key>_<YYYYMMDD>_classification.nc`` convention).
    """
    u = str(units).lower()
    if "since" in u:
        base = u.split("since", 1)[1].strip().replace("t", " ").replace("z", "").strip()
        t0 = None
        for fmt, n in (("%Y-%m-%d %H:%M:%S", 19), ("%Y-%m-%d %H:%M", 16), ("%Y-%m-%d", 10)):
            try:
                t0 = np.datetime64(datetime.strptime(base[:n], fmt), "s")
                break
            except ValueError:
                continue
        if t0 is not None:
            per_s = {"second": 1.0, "minute": 60.0, "hour": 3600.0, "day": 86400.0}
            scale = next((v for k, v in per_s.items() if u.startswith(k)), 3600.0)
            return t0 + (vals * scale).astype("timedelta64[s]")
    day = np.datetime64(datetime.strptime(_date_from_name(fname), "%Y%m%d"), "s")
    return day + (vals * 3600.0).astype("timedelta64[s]")          # decimal hours (Cloudnet default)


def _date_from_name(fname: str) -> str:
    for part in Path(fname).stem.split("_"):
        if len(part) == 8 and part.isdigit():
            return part
    raise ValueError(f"no YYYYMMDD in classification filename {fname!r}")
