# -*- coding: utf-8 -*-
"""Vaisala CL31 / CL51 type of ONE L1 file, from its range extent.

The station census carries one type per stream, taken from the newest L1 ``instrument_type``
attribute, and that attribute lags a unit swap by one to twelve weeks. At four Czech stations a CL31
was replaced by a CL51 under the same stream identifier (Plzen 2025-07-16, Tusimice 2026-06-25,
Libus 2026-06-26, Usti nad Labem 2026-09-16), so one type per stream calibrates one side of the swap
with the wrong eta table (3-4 % on C_L).

The range extent is unambiguous across the network (scan of the newest file of all 258 Vaisala
streams, 2026-11-16): every CL31 ends at 7.7 km (385 or 770 gates), every CL51 at 15.4 km
(1540 gates). A file whose last gate lies beyond 10 km is therefore a CL51, otherwise a CL31.
Other instrument types are returned unchanged.
"""
from __future__ import annotations

from functools import lru_cache
from os import PathLike

import numpy as np
from netCDF4 import Dataset

VAISALA_TYPES = frozenset({"CL31", "CL51"})
#: Last-gate range (m) above which a Vaisala file is a CL51 (CL31: 7.7 km, CL51: 15.4 km).
CL51_MIN_RANGE_M = 10_000.0


@lru_cache(maxsize=8192)
def _max_range_m(path: str) -> float | None:
    try:
        with Dataset(path) as nc:
            r = np.asarray(nc.variables["range"][:], dtype=float)
    except (OSError, KeyError, RuntimeError, ValueError):
        return None
    r = r[np.isfinite(r)]
    return float(r.max()) if r.size else None


def file_instrument_type(path: str | PathLike, census_type: str) -> str:
    """Instrument type of one L1 file: CL31 or CL51 from the range extent for a Vaisala stream,
    ``census_type`` for every other type and for a file whose range cannot be read."""
    if census_type not in VAISALA_TYPES:
        return census_type
    top = _max_range_m(str(path))
    if top is None:
        return census_type
    return "CL51" if top > CL51_MIN_RANGE_M else "CL31"
