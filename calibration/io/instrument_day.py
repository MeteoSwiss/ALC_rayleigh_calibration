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

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import List, Optional, Tuple, Union

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
    target_time_s: float = 15.0,
    read_cams: bool = True,
    auto_download_cams: bool = False,
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
        Working-grid config (default 10 m / 15 s). The block factor is derived per file.
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

    working = average_ceilometer_data(
        native, average_time_s=target_time_s, average_range_m=target_range_m
    )
    factors = _applied_factors(native, target_range_m, target_time_s)

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
