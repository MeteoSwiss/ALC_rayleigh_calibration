# -*- coding: utf-8 -*-
"""Instrument-unit changes per stream, to mark them on the C_L time series.

A calibration constant is a property of ONE optical unit: a module swap moves it by tens of percent
(report 19: CHM15k swaps +12 % median; the Payerne CL31 optical block replacement of 2026-07-07
moved it 2.8x), so a series that hides the swap invites reading a step as drift. Three sources:

* the daily housekeeping CSV (``<key>_hk.csv``, columns ``unit_id`` / ``firmware``) written by the
  monitoring pass from the L1 global attributes -- ``optical_module_id`` (CHM15k ``TUB...``) or
  ``instrument_serial_number`` (CL61 ``X...``, Mini-MPL), and ``instrument_firmware_version``;
* the calibration NetCDFs (``optical_module_id`` per night), so CHM15k history is marked even before
  the housekeeping CSV carries the new columns;
* an operator-maintained list, ``monitoring/instrument_events.json``, for what no file reports: the
  Vaisala CL31/CL51 L1 carries the placeholder serial ``CL0``, so an optical-block replacement on
  those is invisible in the data.
"""
from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

#: Serials that identify nothing (Vaisala CL31/CL51 report "CL0"; fill values; unknown).
PLACEHOLDER_IDS = {"", "CL0", "0", "UNKNOWN", "NONE", "NAN", "-999", "-999.0", "N/A"}
EVENTS_FILE = Path(__file__).resolve().parent / "instrument_events.json"


def _clean(v) -> str:
    s = str(v).strip() if v is not None else ""
    if s.endswith(".0") and s[:-2].isdigit():        # 140016.0 read back as float
        s = s[:-2]
    return "" if s.upper() in PLACEHOLDER_IDS else s


def _norm_unit(s: str) -> str:
    """Same unit, same spelling: NetCDFs store the CHM15k module as the bare number (140016), L1 as
    TUB140016 -- compare the digits so the two sources do not report a false change."""
    s = _clean(s)
    return s[3:] if s.upper().startswith("TUB") else s


def _changes(dates, values, kind: str, what: str) -> list:
    """Run-length changes of a daily/nightly identifier, ignoring placeholders and gaps."""
    out, last, last_raw = [], None, None
    for d, v in sorted(zip(dates, values), key=lambda x: x[0]):
        raw = _clean(v)
        if not raw:
            continue
        cur = _norm_unit(raw) if kind == "unit" else raw
        if last is not None and cur != last:
            out.append({"date": pd.Timestamp(d), "kind": kind,
                        "label": f"{what} {last_raw} → {raw}"})
        last, last_raw = cur, raw
    return out


def _from_hk(fullcal_dir: Path, key: str) -> tuple[list, bool]:
    p = Path(fullcal_dir) / key / f"{key}_hk.csv"
    if not p.exists():
        return [], False
    try:
        df = pd.read_csv(p, dtype=str)
    except Exception:                                                       # noqa: BLE001
        return [], False
    if "date" not in df.columns:
        return [], False
    dates = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    out, have_unit = [], False
    if "unit_id" in df.columns and df["unit_id"].map(_clean).astype(bool).any():
        have_unit = True
        out += _changes(dates, df["unit_id"], "unit", "optical unit")
    if "firmware" in df.columns:
        out += _changes(dates, df["firmware"], "firmware", "firmware")
    return out, have_unit


def _from_netcdf(fullcal_dir: Path, key: str) -> list:
    try:
        import netCDF4
    except ImportError:                                                     # pragma: no cover
        return []
    dates, vals = [], []
    for f in sorted(Path(fullcal_dir, key).glob("**/ALC_calibration_*.nc")):
        try:
            with netCDF4.Dataset(f) as d:
                if "optical_module_id" not in d.variables:
                    continue
                t = np.asarray(d.variables["time"][:], float)
                m = np.ma.filled(d.variables["optical_module_id"][:], -999)
                dates += [datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t]
                vals += [str(int(x)) for x in np.asarray(m).ravel()]
        except Exception:                                                   # noqa: BLE001
            continue
    return _changes(dates, vals, "unit", "optical module")


def _from_events(key: str, events_file: Path = EVENTS_FILE) -> list:
    try:
        ev = json.loads(Path(events_file).read_text(encoding="utf-8"))
    except Exception:                                                       # noqa: BLE001
        return []
    out = []
    for e in ev.get(key, []):
        try:
            out.append({"date": pd.Timestamp(e["date"]), "kind": "event", "label": str(e["label"])})
        except Exception:                                                   # noqa: BLE001
            continue
    return out


def unit_changes(fullcal_dir, key: str, events_file: Path = EVENTS_FILE) -> list:
    """All known changes for one stream, oldest first: [{date, kind, label}].

    kind: 'unit' (optical module / serial changed), 'firmware', or 'event' (operator list)."""
    if not fullcal_dir:
        return _from_events(key, events_file)
    hk, have_unit = _from_hk(Path(fullcal_dir), key)
    out = list(hk)
    if not have_unit:
        out += _from_netcdf(Path(fullcal_dir), key)
    out += _from_events(key, events_file)
    return sorted(out, key=lambda c: c["date"])


def current_identity(fullcal_dir, key: str) -> dict:
    """The stream's latest {unit_id, firmware, date} from the housekeeping CSV ('' when unknown).

    A placeholder serial (Vaisala CL31/CL51 'CL0') is reported as such rather than hidden, so the
    page can say the unit is NOT identified in the data instead of showing nothing."""
    out = {"unit_id": "", "firmware": "", "date": "", "placeholder": False}
    if not fullcal_dir:
        return out
    p = Path(fullcal_dir) / key / f"{key}_hk.csv"
    if not p.exists():
        return out
    try:
        df = pd.read_csv(p, dtype=str).fillna("")
    except Exception:                                                       # noqa: BLE001
        return out
    if "date" not in df.columns:
        return out
    df = df.sort_values("date")
    if "unit_id" in df.columns:
        raw = df["unit_id"].str.strip()
        real = raw.map(_clean)
        if real.astype(bool).any():
            i = real[real.astype(bool)].index[-1]
            out.update(unit_id=real[i], date=str(df.at[i, "date"]))
        elif raw.astype(bool).any():
            out.update(unit_id=raw[raw.astype(bool)].iloc[-1], placeholder=True)
    if "firmware" in df.columns:
        fw = df["firmware"].map(_clean)
        if fw.astype(bool).any():
            out["firmware"] = fw[fw.astype(bool)].iloc[-1]
    return out


def changes_table_html(changes: list) -> str:
    """A compact, collapsible table of the instrument changes (newest first) for static pages."""
    if not changes:
        return ""
    names = {"unit": "optical unit", "firmware": "firmware", "event": "event"}
    rows = "".join(f'<tr><td>{c["date"]:%Y-%m-%d}</td><td>{names.get(c["kind"], c["kind"])}</td>'
                   f'<td>{c["label"]}</td></tr>' for c in sorted(changes, key=lambda c: c["date"],
                                                                 reverse=True))
    return (f'<details class="card tablecard"><summary>Instrument changes ({len(changes)})'
            f'</summary><table class="tbl"><thead><tr><th>date</th><th>what</th><th>change</th>'
            f'</tr></thead><tbody>{rows}</tbody></table></details>')
