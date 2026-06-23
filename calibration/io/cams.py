"""CAMS_Beta file resolution and optional on-demand download.

The Rayleigh (molecular + water-vapor) and liquid-cloud calibrations need a CAMS
model-level file (``CAMS_Beta_*.nc``) covering the processed night. Historically the
archive held one monthly file per month (``CAMS_Beta_YYYYMM.nc``). This module lets a
calibration:

* accept **either** a monthly file (``CAMS_Beta_YYYYMM.nc``) **or** a daily file
  (``CAMS_Beta_YYYYMMDD.nc``) for the night, and
* optionally **download the missing file on the fly** from the Atmosphere Data Store
  (see :mod:`calibration.io.download_cams_beta`).

A night ending on date ``D`` spans ``D-1`` evening .. ``D`` morning, so the daily file
is downloaded for **both** ``D-1`` and ``D`` and stored as ``CAMS_Beta_<D>.nc`` (the
calibration reader then averages whichever forecast steps fall inside the night window).

Resolution order is monthly first (it covers the whole night and is the archive's
native granularity), then daily. ``find_cams_file`` never downloads; ``ensure_cams_file``
downloads only when ``auto_download=True`` and nothing is found.
"""
from __future__ import annotations

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional, Union

_PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def monthly_name(date_str: str) -> str:
    """``CAMS_Beta_YYYYMM.nc`` for the month containing *date_str* (YYYYMM[DD])."""
    return f"CAMS_Beta_{date_str[:6]}.nc"


def daily_name(date_str: str) -> str:
    """``CAMS_Beta_YYYYMMDD.nc`` for the day *date_str* (YYYYMMDD)."""
    return f"CAMS_Beta_{date_str[:8]}.nc"


def candidate_cams_files(cams_folder: _PathLike, date_str: str) -> List[Path]:
    """Candidate CAMS paths for the night ending on *date_str* (YYYYMMDD), in
    preference order: monthly first, then daily. Existence is **not** checked here."""
    folder = Path(cams_folder)
    cands = [folder / monthly_name(date_str)]
    if len(date_str) >= 8:
        cands.append(folder / daily_name(date_str))
    return cands


def find_cams_file(cams_folder: _PathLike, date_str: str) -> Optional[Path]:
    """Return the first existing CAMS file (monthly preferred, then daily) for the
    night ending on *date_str*, or ``None`` if neither exists. Never downloads."""
    for path in candidate_cams_files(cams_folder, date_str):
        if path.exists():
            return path
    return None


def ensure_cams_file(
    cams_folder: _PathLike,
    date_str: str,
    *,
    auto_download: bool = False,
    scope: str = "day",
    log: Optional[logging.Logger] = None,
) -> Optional[Path]:
    """Resolve the CAMS file for the night ending on *date_str*, downloading it from
    the ADS if it is missing and *auto_download* is enabled.

    Parameters
    ----------
    cams_folder : path-like
        Directory holding (and receiving) ``CAMS_Beta_*.nc``.
    date_str : str
        Night-end date ``YYYYMMDD`` (``YYYYMM`` is accepted for a monthly-only lookup).
    auto_download : bool
        If ``True`` and no file exists, fetch it from the ADS. Requires ``cdsapi`` +
        ADS credentials, and the eccodes ``grib_to_netcdf`` tool on PATH.
    scope : {'day', 'month'}
        What to download when missing:

        * ``'day'`` (default) → ``CAMS_Beta_<D>.nc`` covering just the night (dates
          ``D-1`` and ``D``). Light, and the only option that works for a *recent*
          date whose month is not yet complete.
        * ``'month'`` → ``CAMS_Beta_<YYYYMM>.nc`` for the whole month. Heavy (~GBs) but
          reused by every night that month; only valid once the month is complete.

    Returns
    -------
    Path or None
        Path to a usable CAMS file, or ``None`` if it is missing and was not (or could
        not be) downloaded.
    """
    log = log or logger

    existing = find_cams_file(cams_folder, date_str)
    if existing is not None:
        return existing
    if not auto_download:
        return None

    folder = Path(cams_folder)
    folder.mkdir(parents=True, exist_ok=True)

    if scope == "month":
        out_path = folder / monthly_name(date_str)
        dates = _month_dates(date_str)
    elif scope == "day":
        out_path = folder / daily_name(date_str)
        dates = _night_dates(date_str)
    else:
        raise ValueError(f"Unknown CAMS download scope {scope!r} (expected 'day' or 'month')")

    log.warning(
        "CAMS file missing for %s; auto-downloading %s (%s..%s) from the ADS",
        date_str, out_path.name, dates[0], dates[-1],
    )
    # Import lazily: cdsapi / cfgrib are only needed for an actual download, so the
    # resolver stays importable (and find_cams_file usable) without them.
    from .download_cams_beta import download_to_netcdf, CALIBRATION_VARIABLES

    try:
        # Only the 4 model-level fields the calibration reads (t/q/z/lnsp) — keeps the
        # request under the ADS cost limit and the daily file small.
        download_to_netcdf(dates, out_path, variables=CALIBRATION_VARIABLES)
    except Exception as exc:  # noqa: BLE001 — surface any download/convert failure as a skip
        log.error("CAMS auto-download failed for %s: %s", out_path.name, exc)
        # Remove a partial/corrupt file so a later run can retry cleanly.
        try:
            if out_path.exists():
                out_path.unlink()
        except OSError:
            pass
        return None

    return out_path if out_path.exists() else None


def _night_dates(date_str: str) -> List[str]:
    """``['YYYY-MM-DD'(D-1), 'YYYY-MM-DD'(D)]`` so a daily file covers the whole night."""
    day = datetime.strptime(date_str[:8], "%Y%m%d").date()
    prev = day - timedelta(days=1)
    return [prev.strftime("%Y-%m-%d"), day.strftime("%Y-%m-%d")]


def _month_dates(date_str: str) -> List[str]:
    """Every ``'YYYY-MM-DD'`` in the month containing *date_str*."""
    year, month = int(date_str[:4]), int(date_str[4:6])
    start = datetime(year, month, 1).date()
    if month == 12:
        end = datetime(year + 1, 1, 1).date() - timedelta(days=1)
    else:
        end = datetime(year, month + 1, 1).date() - timedelta(days=1)
    n = (end - start).days + 1
    return [(start + timedelta(days=i)).strftime("%Y-%m-%d") for i in range(n)]
