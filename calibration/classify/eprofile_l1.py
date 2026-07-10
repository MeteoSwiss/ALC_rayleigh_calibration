"""Read E-PROFILE L1 ALC data into a ceilopyter ``Ceilo`` for the Cloudnet classifier.

``ceiloclass`` classifies a ceilopyter ``Ceilo`` (screened attenuated backscatter ``beta`` in
sr-1 m-1, plus ``depol`` for a CL61). Its native readers parse the *raw* instrument files;
E-PROFILE instead harmonises every ceilometer to one netCDF layout whose signal lives in
``rcs_0`` (the range-corrected signal). This module bridges the gap, two ways:

* :func:`read_eprofile_l1` — read L1 file(s) directly (native resolution), for standalone use;
* :func:`ceilo_from_shared` — build the same ``Ceilo`` from an already-loaded
  :class:`~calibration.io.instrument_day.InstrumentDayData` (the read-once working grid), so the
  classifier rides on the shared L1+CAMS load with no extra file read.

Physics of the conversion
-------------------------
The lidar identity ``RCS = P r^2 = C * beta_att`` gives ``beta_att = rcs_0 * calibration_factor``,
with the per-type factor chosen so ``beta`` lands in the same sr-1 m-1 scale the native ceilopyter
readers produce (so their noise floors and the classifier's adaptive threshold keep working):

======  ==================  =====================
type    E-PROFILE rcs_0     calibration_factor
======  ==================  =====================
CHM15k  range-corr. signal  3e-12  (== native)
CL31    beta * 1e8          1e-8   (rcs_0 / 1e8)
CL51    beta * 1e8          1e-8   (rcs_0 / 1e8)
CL61    beta_att (physical) 1.0    (rcs_0 already beta)
======  ==================  =====================
"""

from __future__ import annotations

from os import PathLike

import netCDF4
import numpy as np
from cftime import num2pydate
from numpy import ma

from ceilopyter.ceilo import Ceilo
from ceilopyter.ceilo_raw import CeiloRaw, concatenate_raw
from ceilopyter.noise import NOISE_FLOORS, screen_noise

# Per E-PROFILE ``instrument_type``: ceilopyter reader key (selects the noise floor), fallback
# wavelength (nm) if the file lacks ``l0_wavelength``, and the rcs_0 -> native-scale beta factor.
_SPECS: dict[str, tuple[str, float, float]] = {
    "CHM15k": ("chm15k", 1064.0, 3e-12),
    "CHM15kx": ("chm15k", 1064.0, 3e-12),
    "CL31": ("cl31", 910.0, 1e-8),
    "CL51": ("cl51", 910.0, 1e-8),
    "CL61": ("cl61", 910.55, 1.0),
}


def instrument_type(file: str | PathLike) -> str:
    """Return the E-PROFILE ``instrument_type`` global attribute of an L1 file."""
    with netCDF4.Dataset(file) as nc:
        return str(nc.instrument_type)


def station_altitude(file: str | PathLike) -> float | None:
    """Site altitude (m a.s.l.) from ``station_altitude`` (E-PROFILE's name for it).

    ceiloclass's own ``read_altitude`` looks for a variable literally called ``altitude`` and so
    returns ``None`` on E-PROFILE files; we read the right variable and hand the value to
    ``classify`` to align the model profile.
    """
    with netCDF4.Dataset(file) as nc:
        if "station_altitude" in nc.variables:
            return float(np.asarray(nc["station_altitude"][:]).ravel()[0])
    return None


def _spec_for(itype: str) -> tuple[str, float, float]:
    if itype not in _SPECS:
        raise ValueError(f"Unsupported E-PROFILE instrument_type {itype!r}; known: {sorted(_SPECS)}")
    return _SPECS[itype]


def _screen(raw: CeiloRaw, key: str, calibration_factor: float) -> Ceilo:
    """Common tail of both readers: scale to native beta, noise-screen, mask depol, wrap."""
    beta_raw = raw.beta * calibration_factor
    beta = screen_noise(beta_raw, raw.range, noise_floor=NOISE_FLOORS[key])
    if raw.depol is not None:
        # Raw depol is meaningless where there is no signal; keep only values under the
        # backscatter mask, exactly as ceilopyter's CL61 reader does.
        raw.depol = ma.masked_where(ma.getmaskarray(beta), raw.depol)
    return Ceilo(raw, beta_raw, beta, calibration_factor)


def _read_one(file: str | PathLike, wavelength_fallback: float) -> CeiloRaw:
    with netCDF4.Dataset(file) as nc:
        time = num2pydate(nc["time"][:], nc["time"].units)
        rng = np.asarray(nc["range"][:], dtype=float)
        beta = ma.asarray(nc["rcs_0"][:], dtype=float)  # un-calibrated raw signal (scaled later)
        if "l0_wavelength" in nc.variables:
            wavelength = float(np.asarray(nc["l0_wavelength"][:]).ravel()[0])
        else:
            wavelength = wavelength_fallback
        depol = (
            ma.asarray(nc["linear_depol_ratio"][:], dtype=float)
            if "linear_depol_ratio" in nc.variables
            else None
        )
        internal_temperature = (
            np.asarray(nc["temp_int"][:], dtype=float) if "temp_int" in nc.variables else None
        )
    # zenith_angle left None: E-PROFILE ceilometers are near vertical (tilt <~3 deg) and
    # classify()/plot use range as height directly, so it is never consumed downstream.
    return CeiloRaw(time, rng, beta, wavelength, zenith_angle=None, depol=depol,
                    internal_temperature=internal_temperature)


def read_eprofile_l1(
    files: str | PathLike | list[str | PathLike],
    calibration_factor: float | None = None,
) -> Ceilo:
    """Read E-PROFILE L1 file(s) for one instrument into a ceilopyter ``Ceilo`` (native grid).

    Concatenate the raw ``rcs_0``, scale to native-units beta with the per-type
    ``calibration_factor``, screen with the instrument's noise floor, and (CL61) reuse the
    backscatter mask for the otherwise-noisy depolarization.
    """
    if not isinstance(files, list):
        files = [files]
    key, wl_fallback, default_factor = _spec_for(instrument_type(files[0]))
    factor = default_factor if calibration_factor is None else calibration_factor
    concat = concatenate_raw([_read_one(f, wl_fallback) for f in files])
    return _screen(concat, key, factor)


def ceilo_from_shared(idd, calibration_factor: float | None = None) -> Ceilo:
    """Build the ceilopyter ``Ceilo`` from a shared ``InstrumentDayData`` (read-once path).

    The operational equivalent of :func:`read_eprofile_l1`: it takes the already-loaded, coarse
    ``working`` grid (rcs_0, range, time, and — for a CL61 — depolarization) instead of re-reading
    the file, so the classifier consumes the same L1 read as every calibration pass. ``idd`` is
    duck-typed (only ``.instrument_type`` / ``.working`` / ``.wavelength_nm_file``).
    """
    key, wl_fallback, default_factor = _spec_for(idd.instrument_type.value)
    factor = default_factor if calibration_factor is None else calibration_factor
    w = idd.working
    time = list(np.asarray(w.time_datetime))
    rng = np.asarray(w.range_alc, dtype=float)
    beta = ma.masked_invalid(np.asarray(w.rcs, dtype=float))          # (time, range), raw rcs_0
    wavelength = float(idd.wavelength_nm_file) if idd.wavelength_nm_file else wl_fallback
    depol = None if w.depol is None else ma.masked_invalid(np.asarray(w.depol, dtype=float))
    raw = CeiloRaw(time, rng, beta, wavelength, zenith_angle=None, depol=depol,
                   internal_temperature=None)
    return _screen(raw, key, factor)
