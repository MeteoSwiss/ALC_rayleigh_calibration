"""Lean multi-day L1 loader for the OmB and sensitivity passes.

``read_ceilometer_data`` is geared to a single calibration night (and filters to
the solar night for Rayleigh). OmB needs all cloud-free periods over a window and
sensitivity needs clear day *and* night, so both want the raw range-corrected
signal across many days without any night/cloud pre-filtering. This loader
concatenates the native ``rcs_0`` (kept as float32 to bound memory) plus the
cloud base height and station metadata for a list of daily L1 files.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import netCDF4
from numpy.typing import NDArray


def load_l1_window(paths: Iterable[str]) -> Optional[dict]:
    """Concatenate the daily L1 files in *paths* (chronological).

    Returns a dict with ``time`` (datetime64[ns]), ``rcs`` (time x range,
    float32), ``cbh`` (time, lowest cloud base AGL [m], NaN = clear), ``range``
    (m AGL), ``wl`` [nm], ``lat``, ``lon``, ``alt``; or ``None`` if no file
    yields a usable ``rcs_0``. Days whose range grid differs from the first are
    skipped (a grid change mid-window would misalign the matrix).
    """
    times, rcs, cbhs = [], [], []
    rng = wl = lat = lon = alt = None
    for fp in paths:
        try:
            ds = netCDF4.Dataset(str(fp))
        except OSError:
            continue
        try:
            if "rcs_0" not in ds.variables or "range" not in ds.variables:
                continue
            r = np.asarray(ds.variables["range"][:], dtype="float64")
            if rng is None:
                rng = r
                wl = float(ds.variables["l0_wavelength"][...])
                lat = float(ds.variables["station_latitude"][...])
                lon = float(ds.variables["station_longitude"][...])
                alt = float(ds.variables["station_altitude"][...])
            elif r.size != rng.size:
                continue  # range grid changed -> cannot concatenate this day
            days = np.asarray(ds.variables["time"][:], dtype="float64")
            t = np.datetime64("1970-01-01") + (days * 86400.0 * 1e9).astype("timedelta64[ns]")
            sig = np.asarray(ds.variables["rcs_0"][:], dtype="float32")
            cb = _lowest_cbh(ds)
        finally:
            ds.close()
        if sig.shape[0] != t.size or sig.shape[1] != rng.size:
            continue
        times.append(t)
        rcs.append(sig)
        cbhs.append(cb)
    if not times:
        return None
    time = np.concatenate(times)
    rcs_all = np.concatenate(rcs, axis=0)
    cbh_all = np.concatenate(cbhs)
    order = np.argsort(time)
    return dict(time=time[order], rcs=rcs_all[order], cbh=cbh_all[order],
                range=rng, wl=wl, lat=lat, lon=lon, alt=alt)


def _lowest_cbh(ds: netCDF4.Dataset) -> NDArray:
    """Lowest cloud base height per profile [m AGL], NaN where clear/none.

    For Vaisala ceilometers (CL31/CL51/CL61) a reported ``vertical_visibility``
    (fog or strong precipitation, when the sky is obscured and no cloud base is
    returned) must be treated as a low cloud base so those profiles are screened
    out of the clear-sky OmB. CHM15k/Mini-MPL do not report it (CHM15k's ``vor``
    is a different quantity), so the ``vertical_visibility`` variable name keys
    exactly the Vaisala instruments."""
    nt = ds.variables["rcs_0"].shape[0]
    if "cloud_base_height" in ds.variables:
        cb = np.asarray(ds.variables["cloud_base_height"][:], dtype="float64")
        cb = np.where(cb > 0, cb, np.nan)
        with np.errstate(invalid="ignore"):
            cbh = np.nanmin(cb, axis=1) if cb.ndim == 2 else cb
    else:
        cbh = np.full(nt, np.nan)
    if "vertical_visibility" in ds.variables:        # Vaisala fog/precip -> low cloud base
        vis = np.asarray(ds.variables["vertical_visibility"][:], dtype="float64").reshape(-1)
        vis = np.where(vis > 0, vis, np.nan)
        if vis.size == cbh.size:
            cbh = np.fmin(cbh, vis)                   # fmin ignores NaN: vis where clear, min where both
    return cbh
