"""Read one instrument-day's L1 + CAMS ONCE and share a coarsened working grid.

Phase B foundation (steps B1 + B2). This module only *orchestrates* existing, tested
functions so that, per instrument-day:

* the L1 file(s) are read **once**              -> :func:`load_l1_data`
* a shared **working grid** is derived          -> :func:`average_ceilometer_data`
    (dynamic, per-file bin-averaging to a config target -- ``round(target / median
     native step)``, simple block means, cloud-base by *min* -- i.e. the pipeline's own
     reduction, so it is identical to what the cloud step already produces)
* the **CAMS closest cell** is read once         -> :func:`_cams_levels_all_times`
    (lazy nearest-cell read of the 8.5 GB monthly file, LRU-cached by (file, lat, lon))
* the two-way **WV transmission** is computed once on the working grid (910 nm only)
                                                 -> :func:`compute_wv_transmission`

It is purely **additive**: no existing step is modified. Downstream steps (Rayleigh,
cloud, OmB, sensitivity, classification, housekeeping) can later consume the one
:class:`InstrumentDayData` instead of each re-reading L1/CAMS and recomputing WV
(Phase B5). Note that the sensitivity/OmB *noise* metrics should read ``native`` (not
``working``): block-averaging N gates drops per-gate noise by ~sqrt(N) and would
otherwise understate the instrument's true detection limit.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import List, Optional, Tuple, Union

import netCDF4
import numpy as np
from numpy.typing import NDArray

from ..config import InstrumentType
from .cams import ensure_cams_file
from .data_loader import CeilometerData, average_ceilometer_data, load_l1_data

# MATLAB datenum of 1970-01-01: `compute_wv_transmission` compares CAMS datenums against
# ``data.time_num`` on this convention (see ``_cams_levels_all_times``).
_MATLAB_DATENUM_1970 = 719529.0

# The water-vapour correction applies to 910 nm ceilometers only; 1064 nm (CHM15k) has
# negligible WV absorption and skips it.
_WV_WAVELENGTHS_NM = {910.0}

# Backscatter variable names in the order read_ceilometer_data picks them; for an L1 file
# this resolves to rcs_0 (the same variable load_l1_data reads), so its units drive cloud's
# beta reconstruction (see _ceilo_from_shared).
_BETA_VARS = (
    "attenuated_backscatter_0", "rcs_0", "beta", "beta_raw",
    "attenuated_backscatter", "beta_att",
)


def _read_meta(path: Path) -> Tuple[Optional[str], Optional[float]]:
    """rcs_0 units (for cloud's beta) and l0_wavelength (for OmB/sensitivity), attribute-only."""
    units: Optional[str] = None
    wl: Optional[float] = None
    try:
        with netCDF4.Dataset(str(path)) as nc:
            for nm in _BETA_VARS:
                if nm in nc.variables:
                    units = getattr(nc.variables[nm], "units", None)
                    break
            if "l0_wavelength" in nc.variables:
                wl = float(np.asarray(nc.variables["l0_wavelength"][...]).ravel()[0])
    except OSError:
        pass
    return units, wl


@dataclass
class InstrumentDayData:
    """One instrument-day, read once and shared across every processing step.

    ``native`` is full resolution; ``working`` is the dynamically bin-averaged grid every
    downstream step should consume -- except the sensitivity/OmB noise metrics, which read
    ``native`` so block-averaging does not understate the per-gate noise. The CAMS closest
    cell (all time steps) and the once-computed WV transmission ride along.
    """

    instrument_type: InstrumentType
    wavelength_nm: float
    rcs_units: Optional[str]              # units of the L1 rcs_0 (for cloud's beta reconstruction)
    wavelength_nm_file: Optional[float]   # the file's l0_wavelength (for OmB/sensitivity)
    native: CeilometerData
    working: CeilometerData
    coarsen: Tuple[int, int]              # (time_factor, range_factor) applied to `working`
    cams_file: Optional[str]
    cams_time_num: Optional[NDArray]      # (n_t,)        MATLAB datenum
    cams_z_asl: Optional[NDArray]         # (n_lev, n_t)  m ASL
    cams_temperature: Optional[NDArray]   # (n_lev, n_t)  K
    cams_nw: Optional[NDArray]            # (n_lev, n_t)  water-vapour number density m^-3
    wv_transmission: Optional[NDArray]    # (n_range_working, n_time_working) two-way T^2; None @1064 nm

    @property
    def latitude(self) -> float:
        return self.native.latitude

    @property
    def longitude(self) -> float:
        return self.native.longitude

    @property
    def altitude(self) -> float:
        return self.native.altitude

    def slice_to_date(self, date: dt.date) -> "InstrumentDayData":
        """A view whose ``native`` profiles are restricted to one UTC date.

        Day-scoped consumers (cloud, monitoring) take this slice of the night ``[D-1, D]``
        union read, so the files are opened once for the whole night yet each service sees
        only its day. Only ``native`` is sliced; the working grid / CAMS cell / WV transmission
        are dropped (rebuild if a sliced consumer needs them -- cloud does not).
        """
        n = self.native
        t_days = np.asarray(n.time, dtype="float64")  # UTC days since 1970-01-01
        day0 = int((np.datetime64(date) - np.datetime64("1970-01-01")).astype("timedelta64[D]").astype(int))
        idx = np.nonzero(np.floor(t_days).astype("int64") == day0)[0]

        def _sl(a):
            if a is None:
                return None
            arr = np.asarray(a)
            return arr[idx] if arr.ndim >= 1 and arr.shape[0] == t_days.size else a

        sliced = dataclasses.replace(
            n,
            time=_sl(n.time),
            time_datetime=[n.time_datetime[i] for i in idx],
            hours_since_start=_sl(n.hours_since_start),
            rcs=_sl(n.rcs),
            cbh=_sl(n.cbh),
            temperature_optical_module=_sl(n.temperature_optical_module),
            window_transmission=_sl(n.window_transmission),
            status_laser=_sl(n.status_laser),
            status_detector=_sl(n.status_detector),
            laser_life_time=_sl(n.laser_life_time),
            calibration_pulse=_sl(n.calibration_pulse),
            vertical_visibility=_sl(n.vertical_visibility),
        )
        return dataclasses.replace(
            self, native=sliced, working=sliced,
            cams_time_num=None, cams_z_asl=None, cams_temperature=None, cams_nw=None,
            wv_transmission=None,
        )

    def to_omb_dict(self) -> dict:
        """The ``load_l1_window`` dict (OmB / sensitivity input) built from this native read.

        Matches ``calibration.io.l1_window.load_l1_window`` field for field: float32 ``rcs``,
        lowest cloud base (with the Vaisala ``vertical_visibility`` floor), the file
        ``l0_wavelength``, time-sorted -- so OmB/sensitivity can consume the shared read
        instead of re-opening the file.
        """
        n = self.native
        days = np.asarray(n.time, dtype="float64")
        time = np.datetime64("1970-01-01") + (days * 86400.0 * 1e9).astype("timedelta64[ns]")
        rcs = np.asarray(n.rcs, dtype="float32")
        cb = np.asarray(n.cbh, dtype="float64")
        cb = np.where(cb > 0, cb, np.nan)
        with np.errstate(invalid="ignore"):
            cbh = np.nanmin(cb, axis=1) if cb.ndim == 2 else cb
        vv = n.vertical_visibility
        if vv is not None:
            vis = np.asarray(vv, dtype="float64").reshape(-1)
            vis = np.where(vis > 0, vis, np.nan)
            if vis.size == cbh.size:
                cbh = np.fmin(cbh, vis)  # fmin ignores NaN: vis where clear, min where both
        order = np.argsort(time)
        wl = self.wavelength_nm_file if self.wavelength_nm_file is not None else self.wavelength_nm
        return dict(
            time=time[order], rcs=rcs[order], cbh=cbh[order],
            range=np.asarray(n.range_alc, dtype="float64"), wl=wl,
            lat=float(n.latitude), lon=float(n.longitude), alt=float(n.altitude),
        )


def _applied_factors(
    native: CeilometerData, target_range_m: float, target_time_s: float
) -> Tuple[int, int]:
    """The (time, range) block factors ``average_ceilometer_data`` will actually apply.

    Same rule it uses internally -- ``round(target / median native step)``, clamped to >=1
    -- recomputed here only so the result is reportable. Dynamic per file, so it adapts to
    each stream's configured resolution and to changes over time (nothing is hard-coded).
    """
    time_factor = range_factor = 1
    dts = (
        np.diff(np.asarray(native.time_datetime, dtype="datetime64[ns]").astype("int64"))
        .astype("float64")
        / 1e9
    )
    dts = dts[np.isfinite(dts) & (dts > 0)]
    if target_time_s and dts.size:
        time_factor = max(1, int(np.round(float(target_time_s) / float(np.median(dts)))))
    dr = np.diff(np.asarray(native.range_alc, dtype="float64"))
    dr = dr[np.isfinite(dr) & (dr > 0)]
    if target_range_m and dr.size:
        range_factor = max(1, int(np.round(float(target_range_m) / float(np.median(dr)))))
    return time_factor, range_factor


def _wv_transmission_once(
    working: CeilometerData,
    itype: InstrumentType,
    cams_folder: str,
    cams_folder_fallback: str,
    wv_lut: str,
    wv_source: str,
    era5_cache: str,
) -> NDArray:
    """Two-way WV transmission on the WORKING grid, computed once.

    Reuses the production :func:`compute_wv_transmission` verbatim -- so it is bit-identical
    to what the cloud/OmB path produces on the same grid -- via a light data view. That
    function only reads ``.range``, ``.time``, ``.time_num``, ``.beta.shape[1]`` and the
    station coords off the object, and its CAMS read is the same cached
    ``_cams_levels_all_times``.
    """
    # Lazy import to avoid any import cycle (cloud.calibration -> io.cams).
    from ..cloud.calibration import CloudCalConfig, compute_wv_transmission, set_defaults

    time_dt64 = np.asarray(working.time_datetime, dtype="datetime64[ns]")
    time_num = np.asarray(working.time, dtype="float64") + _MATLAB_DATENUM_1970
    view = SimpleNamespace(
        range=np.asarray(working.range_alc, dtype="float64"),
        time=time_dt64,
        time_num=time_num,
        beta=np.asarray(working.rcs, dtype="float64").T,  # (range, time): only .shape[1] is read
        station_latitude=float(working.latitude),
        station_longitude=float(working.longitude),
        station_altitude=float(working.altitude),
    )
    cfg = set_defaults(
        CloudCalConfig(
            instrument=itype.value,
            apply_wv_correction=True,
            cams_folder=str(cams_folder),
            cams_folder_fallback=str(cams_folder_fallback),
            abs_cs_lookup_table=str(wv_lut),
            wv_source=wv_source,
            era5_cache=era5_cache,
            station_latitude=float(working.latitude),
            station_longitude=float(working.longitude),
        )
    )
    return compute_wv_transmission(view, cfg)


def load_instrument_day(
    l1_paths: List[Union[str, Path]],
    instrument_type: Union[InstrumentType, str],
    cams_folder: Union[str, Path],
    *,
    cams_folder_fallback: Union[str, Path] = "",
    wv_lut: Union[str, Path] = "",
    wv_source: str = "cams",
    era5_cache: str = "",
    target_range_m: float = 10.0,
    target_time_s: float = 30.0,
    read_cams: bool = True,
    auto_download_cams: bool = False,
    build_working: bool = True,
) -> Optional[InstrumentDayData]:
    """Read L1 (once) + CAMS (once, closest cell) and build the shared working grid.

    Parameters
    ----------
    l1_paths : list of str/Path
        L1 file(s) for the instrument-day (e.g. a night is ``[D-1, D]``).
    instrument_type : InstrumentType or str
        e.g. ``"CL61"``, ``"CHM15k"``, ``"CL31"``.
    cams_folder : str/Path
        Primary CAMS dir (0.4 deg); ``cams_folder_fallback`` is the 1 deg set.
    wv_lut : str/Path
        Absorption cross-section LUT; ``""`` uses the bundled 910 nm LUT.
    target_range_m, target_time_s : float
        Working-grid config (default 10 m / 30 s). The block factor is derived per file.
        30 s makes every instrument uniform (only CHM15k is finer, at 15 s) and shifts its
        classification by <1 pp; calibration/OmB/sensitivity run on `native`, so it does not
        affect them.
    read_cams : bool
        Read the CAMS cell + compute WV. Set ``False`` to only load/coarsen L1.
    auto_download_cams : bool
        Fetch a missing CAMS file from the ADS (default ``False`` -- offline/local).

    Returns
    -------
    InstrumentDayData or None
        ``None`` if the L1 file(s) could not be read.
    """
    itype = (
        instrument_type
        if isinstance(instrument_type, InstrumentType)
        else InstrumentType(str(instrument_type))
    )

    native = load_l1_data([Path(p) for p in l1_paths], itype)
    if native is None:
        return None
    rcs_units, wl_file = _read_meta(Path(l1_paths[0]))

    # The coarsened working grid is only needed by consumers that classify (or otherwise want
    # the reduced grid). Rayleigh/cloud/OmB/sensitivity all use `native`, so building it there
    # is pure overhead (block-averaging the full rcs) -- skip it unless asked (build_working).
    factors = _applied_factors(native, target_range_m, target_time_s)
    working = (
        average_ceilometer_data(native, average_time_s=target_time_s, average_range_m=target_range_m)
        if build_working
        else native
    )

    cams_file: Optional[str] = None
    cams_cell: Tuple[Optional[NDArray], ...] = (None, None, None, None)
    wv: Optional[NDArray] = None
    if read_cams:
        day8 = native.time_datetime[0].strftime("%Y%m%d")
        path = ensure_cams_file(
            str(cams_folder), day8, auto_download=auto_download_cams,
            latitude=native.latitude, longitude=native.longitude,
        )
        if path is None and str(cams_folder_fallback):
            path = ensure_cams_file(
                str(cams_folder_fallback), day8, auto_download=auto_download_cams,
                latitude=native.latitude, longitude=native.longitude,
            )
        if path is not None:
            cams_file = str(path)
            # Lazy import to avoid any import cycle.
            from ..cloud.calibration import _cams_levels_all_times

            cams_cell = _cams_levels_all_times(cams_file, native.latitude, native.longitude)
            if itype.wavelength_nm in _WV_WAVELENGTHS_NM:
                wv = _wv_transmission_once(
                    working, itype, str(cams_folder), str(cams_folder_fallback),
                    str(wv_lut), wv_source, era5_cache,
                )

    return InstrumentDayData(
        instrument_type=itype,
        wavelength_nm=itype.wavelength_nm,
        rcs_units=rcs_units,
        wavelength_nm_file=wl_file,
        native=native,
        working=working,
        coarsen=factors,
        cams_file=cams_file,
        cams_time_num=cams_cell[0],
        cams_z_asl=cams_cell[1],
        cams_temperature=cams_cell[2],
        cams_nw=cams_cell[3],
        wv_transmission=wv,
    )
