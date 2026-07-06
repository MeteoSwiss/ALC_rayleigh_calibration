"""
overlap.py — CHM15k temperature-dependent overlap correction (Hervo et al., 2016).

The CHM15k near-range signal is attenuated by an incomplete, temperature-dependent laser/telescope
overlap. The E-PROFILE Level-1 `rcs_0` already carries the instrument's STATIC reference overlap
`overlap_ref`; what remains is a residual that depends on the internal temperature. The
`overlap_probe_eprofile` package (github.com/martin-obs) derives, per instrument, a temperature model

    Dif(z) = a(z) * T + b(z)                          [percent]   with T = internal temperature [degC]

fitted as the relative difference between the reference and the daily (true) overlap,

    Dif/100 = (overlap_ref - overlap_daily) / overlap_daily      =>   overlap_daily = overlap_ref / (1 + Dif/100)

so the true (temperature-dependent) overlap is `overlap_ref / (1 + Dif/100)`. Because the L1 `rcs_0`
was formed by dividing the raw signal by `overlap_ref`, correcting it to the true overlap is simply

    rcs_corrected(z) = rcs_0(z) * overlap_ref(z) / overlap_daily(z) = rcs_0(z) * (1 + Dif(z)/100).

The reference overlap cancels, so only `a`, `b` and the internal temperature are needed. This
convention (Celsius, and the exact arithmetic) was verified against `build_temp_model.py` /
`overlap_utils.write_temp_model_to_netcdf` in the source package; the same formula is stored in each
model file's `description` attribute. The package only *builds* the model — there is no signal-apply
code there — so this small kernel is the application step.

The models live in `D:/TEMP_MODELS/202606` as `Overlap_correction_model_<SITE>_<serlom>_<id>.nc`
(one per instrument; each is fitted over 200+ days spanning years, so a single model plus the daily
temperature covers the whole record). A model carries `range` [m AGL], `a`, `b` and `overlap_ref`;
we use `range`, `a`, `b`. The correction is confined to the near range (a and b vanish once overlap
is complete, ~700 m for the CHM15k), so it is a no-op above the ~500 m analysis floor.
"""
from __future__ import annotations
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

MODEL_DIR = Path(os.environ.get("ALC_OVERLAP_DIR", "D:/TEMP_MODELS/202606"))

# Guard rails on the multiplicative factor (1 + Dif/100). Above full overlap Dif = 0 -> factor = 1
# (true no-op). The clamp only matters for the lowest, sub-overlap gates (well below the 500 m
# analysis floor) where the empirical Dif can be large; it keeps the corrected signal finite and
# same-signed rather than letting a Dif < -100 % flip the sign.
_FACTOR_MIN, _FACTOR_MAX = 0.1, 10.0

_index = None            # {(wigos_station_id, instrument_id): path} built once
_cache = {}              # {(wmo, ident): model dict or None}


def _build_index():
    """Map every overlap model file by its (wigos_station_id, instrument_id) global attributes."""
    global _index
    _index = {}
    if not MODEL_DIR.is_dir():
        return
    for fp in sorted(MODEL_DIR.glob("Overlap_correction_model_*.nc")):
        try:
            with Dataset(fp) as nc:
                # files from 2026-07 spell the attribute "insturment_id" (upstream typo) — accept both
                iid = getattr(nc, "instrument_id", getattr(nc, "insturment_id", ""))
                key = (getattr(nc, "wigos_station_id", ""), str(iid))
        except Exception:
            continue
        _index.setdefault(key, fp)


def load_overlap_model(wmo, ident):
    """Return the overlap model for one instrument as dict(range, a, b) in m AGL / percent, or None
    when no model exists for this (WIGOS id, instrument id). Results are cached."""
    if (wmo, ident) in _cache:
        return _cache[(wmo, ident)]
    if _index is None:
        _build_index()
    fp = _index.get((wmo, ident))
    model = None
    if fp is not None:
        with Dataset(fp) as nc:
            model = dict(range=np.asarray(nc["range"][:], "f8"),
                         a=np.asarray(nc["a"][:], "f8"),
                         b=np.asarray(nc["b"][:], "f8"),
                         file=fp.name)
    _cache[(wmo, ident)] = model
    return model


def correct_rcs(rcs, range_agl, temp_C, model):
    """Apply the temperature-dependent overlap correction to a range-corrected signal.

    rcs        : (time, range) range-corrected signal (L1 rcs_0), reference overlap already applied.
    range_agl  : (range,) instrument range above ground [m] (CHM15k L1 range is AGL).
    temp_C     : (time,) internal temperature [degC] per profile.
    model      : dict from load_overlap_model (range, a, b in percent).

    Returns the corrected signal rcs * (1 + (a(z)*T + b(z)) / 100), with a, b interpolated onto
    range_agl (0 above the model's fully-overlapped top gate, so a no-op there) and the factor
    clamped to a physical range. Input NaNs propagate unchanged.
    """
    a = np.interp(range_agl, model["range"], model["a"], left=model["a"][0], right=0.0)
    b = np.interp(range_agl, model["range"], model["b"], left=model["b"][0], right=0.0)
    T = np.asarray(temp_C, "f8")
    dif = np.outer(T, a) + b[None, :]                            # (time, range), percent
    # Profiles with no internal temperature keep factor 1 (no correction) instead of being blanked;
    # the near-range correction is negligible above the analysis floor anyway.
    dif[~np.isfinite(T), :] = 0.0
    factor = np.clip(1.0 + dif / 100.0, _FACTOR_MIN, _FACTOR_MAX)
    return rcs * factor
