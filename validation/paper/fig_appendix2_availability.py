#!/usr/bin/env python3
"""
fig_appendix2_availability.py -- Appendix 2 of the E-Profile ALC paper: data availability.

Produces, in OUT_DIR (default C:/DATA/Projects/202606_E-PROFILE_calibration/
figs_paper_validation/paper_v08/appendix2/):

  fig_A2_data_availability.png  -- per-instrument x month availability matrix (rows grouped by
                                   instrument type, then by date of entry into the network) +
                                   monthly instrument count and availability by type +
                                   per-type distribution of the per-instrument availability.
  appendix2_stats.json          -- every number quoted in the appendix text.
  appendix2_text.md             -- the appendix paragraphs and the Figure A2-1 caption, with all
                                   numbers injected from the stats (nothing hard-coded).

Why this is complementary to Figure 1 (right): Figure 1 counts the stations delivering valid data
per month and per type (network size). This figure shows, for every single instrument, the
fraction of the hours of each month with data, i.e. the continuity of each record.

Data sources
------------
* MAIN reference: D:/E-PROFILE_L2_hourly_monthly_all, the most complete local archive. One file
  per instrument and month (L2_<wigos>_<id><yyyymm>_hourly.nc) holding hourly means of the
  Level-2 product; the file only contains the hours for which the hub produced L2 data (checked
  on a random sample of 200 files: every stored hour has finite attenuated backscatter), so the
  availability of a month = number of time steps / number of hours of the month.
* SUPPLEMENT: F:/E_PROFILE_ALC (older copy of the hub archive, L1 and L2 daily files, plus
  zipped 2018 L2 files). It is used only for station-months that the main archive does not
  contain at all; for those the availability is the number of distinct hours with at least one
  profile in the daily files (one 'time' variable read per file). The zip archives are checked by
  listing their contents only.

Usage:
  python validation/paper/fig_appendix2_availability.py [--out DIR] [--rescan] [--workers 16]

The per-file reads (~30 000 hourly-monthly files + ~12 000 supplement daily files) are cached in
OUT_DIR/cache; --rescan forces a fresh read.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import date
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# %% Configuration

DEFAULT_OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/"
                   "paper_v08/appendix2")

MAIN_ROOT = "D:/E-PROFILE_L2_hourly_monthly_all"
SUPPLEMENT_ROOT = "F:/E_PROFILE_ALC"
SUPPLEMENT_SUBDIRS = ["L1_FILES", "L2_FILES", "L2_FILES_2018_TEMP"]

TYPE_ORDER = ["CHM15k", "CL31", "CL51", "CL61", "Mini-MPL", "Other"]
TYPE_COLOURS = {"CHM15k": "#d62728", "CL31": "#ff7f0e", "CL51": "#9467bd",
                "CL61": "#1f77b4", "Mini-MPL": "#2ca02c", "Other": "#7f7f7f"}

# Monthly availability classes of the matrix (fraction of the hours of the month with data).
AVAIL_BINS = [0.0, 0.5, 0.9, 1.0001]
AVAIL_LABELS = ["< 50 %", "50-90 %", "> 90 %"]
AVAIL_TINTS = [0.30, 0.62, 1.00]        # mixing weight of the type colour with white

# Identifier streams of the same type at one station sharing at most this many months are
# treated as one renamed instrument (see group_streams).
MAX_OVERLAP_MONTHS = 2

# In the matrix, types with fewer instruments than SMALL_GROUP_SIZE are drawn with rows
# SMALL_GROUP_ROW_HEIGHT times taller so that their block remains legible.
SMALL_GROUP_SIZE = 20
SMALL_GROUP_ROW_HEIGHT = 3

# Panel (c) shows a type's mean monthly availability only when at least this many instruments of
# that type are active (with fewer, single-instrument outages dominate the mean).
MIN_ACTIVE_FOR_MEAN = 5
# Paper figure = panel (a) only (operator decision 2026-10-10); --all-panels adds (b)-(d).
PANEL_A_ONLY = True

# Instrument identifier: letters (A, B, ...) or ONE digit (the DWD CHM15k use "0"). A multi-digit identifier must
# not be allowed: it would let the 12-digit timestamp of 5-minute granule files match the daily pattern.
HOURLY_RE = re.compile(r"^L2_(.+)_([A-Za-z]+|\d)(\d{6})_hourly\.nc$")
DAILY_RE = re.compile(r"(L[12])_(.+?)_([A-Za-z]+|\d)(\d{8})\.nc$")

# %% Instrument type


def normalise_type(*texts: str) -> str:
    """Map instrument_type (or, in early L1 files, 'system'/'title') onto the five classes."""
    for raw in texts:
        text = (raw or "").upper().replace(" ", "").replace("-", "").replace("_", "")
        if not text:
            continue
        if "CHM15" in text or "CHM8" in text or "NIMBUS" in text:
            return "CHM15k"
        for vaisala in ("CL31", "CL51", "CL61"):
            if vaisala in text:
                return vaisala
        if "MPL" in text:
            return "Mini-MPL"
    return "Other"


def _header_meta(ds) -> dict:
    """Instrument type, site and latitude from an open netCDF4.Dataset."""
    attrs = {a: str(getattr(ds, a)) for a in ("instrument_type", "system", "title",
                                              "site_location") if a in ds.ncattrs()}
    lat = np.nan
    if "station_latitude" in ds.variables:
        try:
            lat = float(np.ma.filled(ds.variables["station_latitude"][...], np.nan))
        except (TypeError, ValueError):
            pass
    return {"type": normalise_type(attrs.get("instrument_type", ""), attrs.get("system", ""),
                                   attrs.get("title", "")),
            "site": attrs.get("site_location", ""), "lat": lat}


# %% Main archive: hourly-monthly L2 files


def list_main_archive(root: str = MAIN_ROOT) -> list[tuple[str, str, str, str]]:
    """[(wigos, ident, yyyymm, path)] for every hourly-monthly file of the main archive."""
    found = []
    for wigos in sorted(os.listdir(root)):
        for dirpath, _dirs, files in os.walk(os.path.join(root, wigos)):
            for name in files:
                match = HOURLY_RE.match(name)
                if match and match.group(1) == wigos:
                    found.append((wigos, match.group(2), match.group(3),
                                  os.path.join(dirpath, name)))
    return found


def read_hourly_file(path: str) -> dict:
    """Hours with data, first/last time stamp and metadata of one hourly-monthly file."""
    import netCDF4

    try:
        with netCDF4.Dataset(path) as ds:
            t = np.asarray(np.ma.filled(ds.variables["time"][:], np.nan), dtype=float)
            t = t[np.isfinite(t)]
            meta = _header_meta(ds)
    except (OSError, KeyError):
        return {"hours": 0, "tmin": None, "tmax": None, "type": "Other", "site": "",
                "lat": np.nan}
    hours = np.unique(np.floor(t * 24 + 1e-6))
    meta.update({"hours": int(hours.size),
                 "tmin": float(t.min()) if t.size else None,
                 "tmax": float(t.max()) if t.size else None})
    return meta


# %% Supplement archive: daily L1/L2 files on F:


def _walk_daily(path: str) -> list[tuple[str, str, str, str, str]]:
    """[(level, wigos, ident, yyyymmdd, path)] of the daily files below one folder."""
    out = []
    for dirpath, _dirs, files in os.walk(path):
        for name in files:
            match = DAILY_RE.search(name)
            if match:
                out.append(match.groups() + (os.path.join(dirpath, name),))
    return out


def _zip_daily(path: str) -> list[tuple[str, str, str, str, str]]:
    """Same as _walk_daily for the members of one zip archive (listing only, no extraction)."""
    try:
        names = zipfile.ZipFile(path).namelist()
    except (OSError, zipfile.BadZipFile):
        return []
    return [m.groups() + ("zip:" + path,) for m in map(DAILY_RE.search, names) if m]


def list_supplement(workers: int = 16) -> list[tuple[str, str, str, str, str]]:
    """All daily files (folders and zip members) of the supplement archive."""
    entries = [os.path.join(SUPPLEMENT_ROOT, sub, name)
               for sub in SUPPLEMENT_SUBDIRS if os.path.isdir(os.path.join(SUPPLEMENT_ROOT, sub))
               for name in os.listdir(os.path.join(SUPPLEMENT_ROOT, sub))]

    def scan(entry):
        if entry.endswith(".zip"):
            return _zip_daily(entry)
        return _walk_daily(entry) if os.path.isdir(entry) else []

    with ThreadPoolExecutor(workers) as pool:  # I/O bound -> threads
        return [rec for chunk in pool.map(scan, entries) for rec in chunk]


def read_daily_hours(path: str) -> dict:
    """Distinct hours (as integer hour index since 1970) with >= 1 profile in one daily file."""
    import netCDF4

    try:
        with netCDF4.Dataset(path) as ds:
            t = np.asarray(np.ma.filled(ds.variables["time"][:], np.nan), dtype=float)
            meta = _header_meta(ds)
    except (OSError, KeyError):
        return {"hours": [], "type": "Other", "site": "", "lat": np.nan}
    t = t[np.isfinite(t)]
    meta["hours"] = np.unique(np.floor(t * 24 + 1e-6)).astype(int).tolist()
    return meta


def supplement_records(main_station_months: set, workers: int) -> tuple[list, dict]:
    """Monthly records for the station-months that only the supplement archive holds.

    Returns (records, report) where report summarises what the supplement contains and adds.
    """
    daily = list_supplement(workers)
    report = {"root": SUPPLEMENT_ROOT, "n_daily_files_listed": len(daily),
              "n_daily_files_in_zip": sum(1 for d in daily if d[4].startswith("zip:")),
              "n_station_months_listed": len({(d[1], d[3][:6]) for d in daily})}
    extra = [d for d in daily if (d[1], d[3][:6]) not in main_station_months]
    if any(d[4].startswith("zip:") for d in extra):
        # Would need extraction; does not happen with the current archives.
        print("  warning: some supplement-only days are inside zip files and are skipped")
        extra = [d for d in extra if not d[4].startswith("zip:")]
    # One file per (stream, day): prefer L1 (the instrument's own time sampling).
    by_day = {}
    for level, wigos, ident, day, path in sorted(extra):
        by_day.setdefault((wigos, ident, day), path)
    print(f"Reading {len(by_day)} supplement daily files (station-months absent from the main "
          f"archive) ...")
    keys = list(by_day)
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(read_daily_hours, [by_day[k] for k in keys], chunksize=16))

    monthly = defaultdict(lambda: {"hours": set(), "meta": None})
    for (wigos, ident, day), res in zip(keys, results):
        rec = monthly[(wigos, ident, day[:6])]
        rec["hours"].update(res["hours"])
        if rec["meta"] is None or rec["meta"]["type"] == "Other":
            rec["meta"] = res
    records = []
    for (wigos, ident, ym), rec in monthly.items():
        if not rec["hours"]:
            continue
        hours = sorted(rec["hours"])
        records.append({"wigos": wigos, "ident": ident, "ym": ym, "source": "supplement",
                        "hours": len(hours), "tmin": hours[0] / 24, "tmax": (hours[-1] + 1) / 24,
                        "type": rec["meta"]["type"], "site": rec["meta"]["site"],
                        "lat": rec["meta"]["lat"]})
    report["n_station_months_added"] = len({(r["wigos"], r["ym"]) for r in records})
    report["n_stations_with_added_months"] = len({r["wigos"] for r in records})
    added_months = sorted({r["ym"] for r in records})
    report["added_month_range"] = [added_months[0], added_months[-1]] if added_months else None
    report["added_station_months_by_year"] = dict(sorted(Counter(
        r["ym"][:4] for r in {(r["wigos"], r["ym"]): r for r in records}.values()).items()))
    return records, report


# %% Build (or load) all monthly records


def load_records(cache_dir: Path, rescan: bool, workers: int) -> tuple[list, dict]:
    """Monthly records {wigos, ident, ym, source, hours, tmin, tmax, type, site, lat}."""
    cache = cache_dir / "monthly_records.json.gz"
    if cache.is_file() and not rescan:
        with gzip.open(cache, "rt") as fh:
            blob = json.load(fh)
        return blob["records"], blob["supplement_report"]

    main_files = list_main_archive()
    print(f"Reading {len(main_files)} hourly-monthly files from {MAIN_ROOT} ...")
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(read_hourly_file, [f[3] for f in main_files], chunksize=32))
    records = [{"wigos": w, "ident": i, "ym": ym, "source": "main", **res}
               for (w, i, ym, _p), res in zip(main_files, results) if res["hours"] > 0]
    main_station_months = {(r["wigos"], r["ym"]) for r in records}

    print(f"Listing the supplement archive {SUPPLEMENT_ROOT} ...")
    extra, report = supplement_records(main_station_months, workers)
    records += extra

    cache_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(cache, "wt") as fh:
        json.dump({"records": records, "supplement_report": report}, fh)
    return records, report


# %% Instruments (streams) and identifier aliases


MIN_TYPE_MONTHS = 2   # a type reported for fewer months than this in a stream is a header glitch


def split_by_type(streams: dict) -> dict:
    """Split identifier streams whose instrument type changes into one stream per type.

    An identifier can be reused for a different instrument (e.g. Camborne 'C': a CHM15k in 2016, a
    CL61 from 2025; CL31 replaced by CL51 at a few stations). Each type reported for at least
    MIN_TYPE_MONTHS months becomes its own instrument (key ident + '/' + type for the extra types);
    months of a rarer or unknown type stay with the stream's dominant type.
    """
    out = {}
    for (wigos, ident), months in streams.items():
        counts = Counter(r["type"] for r in months.values() if r["type"] != "Other")
        kept = [t for t, n in counts.most_common() if n >= MIN_TYPE_MONTHS]
        if len(kept) <= 1:
            out[(wigos, ident)] = months
            continue
        main = kept[0]
        parts = {t: {} for t in kept}
        for ym, rec in months.items():
            parts[rec["type"] if rec["type"] in parts else main][ym] = rec
        for t, part in parts.items():
            out[(wigos, ident if t == main else f"{ident}/{t}")] = part
    return out


def group_streams(records: list) -> dict:
    """{(wigos, ident): {ym: record}}, with renamed identifiers merged into one instrument.

    Identifier letters were reassigned at some stations (e.g. all UK stations in January 2017:
    'C', 'T', 'j', ... became 'A', 'B'). Counting each letter as a separate instrument would
    inflate the instrument count and show a spurious network-wide drop. Within a station, the
    identifier streams of the same instrument type whose periods do not overlap (at most
    MAX_OVERLAP_MONTHS months in common, to allow for the transition month) are therefore
    merged, the longest stream giving its identifier to the merged row.
    """
    raw = defaultdict(dict)
    for rec in records:
        raw[(rec["wigos"], rec["ident"])][rec["ym"]] = rec
    streams = split_by_type(raw)

    def main_type(key):
        return Counter(r["type"] for r in streams[key].values()).most_common(1)[0][0]

    by_station = defaultdict(list)
    for key in streams:
        by_station[key[0]].append(key)

    merged = {}
    for keys in by_station.values():
        groups = []   # [(target_key, type, {ym: record})]
        for key in sorted(keys, key=lambda k: -len(streams[k])):
            itype = main_type(key)
            for _target, gtype, months in groups:
                overlap = len(months.keys() & streams[key].keys())
                if gtype == itype and overlap <= MAX_OVERLAP_MONTHS:
                    for ym, rec in streams[key].items():
                        if ym not in months or rec["hours"] > months[ym]["hours"]:
                            months[ym] = rec
                    break
            else:
                groups.append((key, itype, dict(streams[key])))
        for target, _gtype, months in groups:
            merged[target] = months
    return merged


def fill_unknown_types(streams: dict) -> None:
    """Months whose file carries no recognisable type inherit the nearest known type."""
    for months in streams.values():
        known = sorted((ym, r["type"]) for ym, r in months.items() if r["type"] != "Other")
        if not known:
            continue
        for ym, rec in months.items():
            if rec["type"] == "Other":
                earlier = [t for k, t in known if k <= ym]
                rec["type"] = earlier[-1] if earlier else known[0][1]


# %% Monthly availability matrix


def month_list(first: str, last: str) -> list[tuple[int, int]]:
    """All (year, month) from yyyymm `first` to `last` inclusive."""
    year, month = int(first[:4]), int(first[4:])
    months = []
    while (year, month) <= (int(last[:4]), int(last[4:])):
        months.append((year, month))
        month += 1
        if month == 13:
            year, month = year + 1, 1
    return months


def hours_in_month(year: int, month: int) -> int:
    nxt = date(year + (month == 12), month % 12 + 1, 1)
    return 24 * (nxt - date(year, month, 1)).days


def build_matrix(streams: dict, months: list) -> dict:
    """Per instrument (row) and month (column): fraction of hours with data and the type.

    Returns keys, frac (NaN = no data), type_ym, source_ym, row_type (latest type), lat, site,
    span_frac (hours with data / hours between first and last time stamp), nhours.
    """
    col = {f"{y}{m:02d}": i for i, (y, m) in enumerate(months)}
    nhours = np.array([hours_in_month(y, m) for y, m in months], dtype=float)
    keys, frac, types, sources, row_type, lats, sites, span = [], [], [], [], [], [], [], []
    for key, recs in streams.items():
        recs = {ym: r for ym, r in recs.items() if ym in col}
        if not recs:
            continue
        row = np.full(len(months), np.nan)
        type_row = np.full(len(months), "", dtype=object)
        src_row = np.full(len(months), "", dtype=object)
        for ym, rec in recs.items():
            row[col[ym]] = min(rec["hours"] / nhours[col[ym]], 1.0)
            type_row[col[ym]] = rec["type"]
            src_row[col[ym]] = rec["source"]
        # Months without data take the type of the nearest earlier month (for colouring only).
        last_type = next(t for t in type_row if t)
        for j in range(len(months)):
            if type_row[j]:
                last_type = type_row[j]
            else:
                type_row[j] = last_type
        latest = recs[max(recs)]
        t0 = min(r["tmin"] for r in recs.values())
        t1 = max(r["tmax"] for r in recs.values())
        total_hours = sum(r["hours"] for r in recs.values())
        keys.append(key)
        frac.append(row)
        types.append(type_row)
        sources.append(src_row)
        row_type.append(latest["type"])
        lats.append(latest["lat"] if latest["lat"] is not None else np.nan)
        sites.append(latest["site"])
        span.append(min(total_hours / max((t1 - t0) * 24, 1.0), 1.0))
    return {"keys": keys, "frac": np.array(frac), "type_ym": np.array(types, dtype=object),
            "source_ym": np.array(sources, dtype=object), "row_type": np.array(row_type),
            "lat": np.array(lats, dtype=float), "site": np.array(sites, dtype=object),
            "span_frac": np.array(span), "nhours": nhours}


def order_rows(mat: dict) -> np.ndarray:
    """Group rows by instrument type (TYPE_ORDER), then by first month with data, then latitude.

    Sorting by entry date inside each type turns each block into a staircase that shows when
    the instruments joined the network.
    """
    first_month = np.array([np.argmax(~np.isnan(r)) for r in mat["frac"]])
    type_rank = np.array([TYPE_ORDER.index(t) for t in mat["row_type"]])
    lat = np.nan_to_num(mat["lat"], nan=0.0)
    return np.lexsort((-lat, first_month, type_rank))


# %% Statistics quoted in the text


def monthly_by_type(mat: dict) -> dict:
    """Per type and month: number of instruments with data, and their mean hour fraction."""
    out = {}
    for itype in TYPE_ORDER:
        active = (mat["type_ym"] == itype) & ~np.isnan(mat["frac"])
        count = active.sum(axis=0)
        total = np.where(active, mat["frac"], 0.0).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean_frac = np.where(count > 0, total / np.maximum(count, 1), np.nan)
        out[itype] = {"count": count, "mean_frac": mean_frac}
    return out


def ym_label(months: list, i: int) -> str:
    return f"{months[i][0]}-{months[i][1]:02d}"


def find_network_dips(n_per_month: np.ndarray, months: list) -> list:
    """Months where the number of active instruments drops by > 10 % from the previous month."""
    return [{"month": ym_label(months, i), "from": int(n_per_month[i - 1]),
             "to": int(n_per_month[i])}
            for i in range(1, len(months)) if n_per_month[i] < 0.9 * n_per_month[i - 1]]


def country_of(site: str) -> str:
    """'PAYERNE,SWITZERLAND' -> 'Switzerland' (site_location attribute of the files)."""
    country = str(site).split(",")[-1].strip().replace("_", " ")
    if len(country) <= 3:
        return country.upper()  # UK
    return country.title() if country else "unknown"


def characterise_dip(mat: dict, months: list, i: int, recovery_months: int = 3) -> dict:
    """Describe a drop between month i-1 and month i.

    Lists the countries of the instruments lost, and whether the drop is temporary (>= 80 % of
    the lost instruments deliver again within `recovery_months`) or persistent. A temporary drop
    that hits every country at once points to an archiving/processing gap rather than to
    instrument failures; a persistent drop concentrated in one country to a national network
    leaving (or pausing) the data exchange.
    """
    lost = ~np.isnan(mat["frac"][:, i - 1]) & np.isnan(mat["frac"][:, i])
    later = (~np.isnan(mat["frac"][:, i + 1:i + 1 + recovery_months])).any(axis=1)
    back = int((lost & later).sum())
    countries = Counter(country_of(s) for s in mat["site"][lost]).most_common()
    still_missing_at_end = int((lost & np.isnan(mat["frac"][:, -1])).sum())
    # Month from which most of the lost instruments had been delivering (median entry month).
    entry = [int(np.argmax(~np.isnan(row))) for row in mat["frac"][lost]]
    entry_label = ym_label(months, int(np.median(entry))) if entry else None
    return {"countries_lost": dict(countries),
            "n_lost": int(lost.sum()),
            "n_back_within_months": back,
            "recovery_months": recovery_months,
            "temporary": bool(back >= 0.8 * lost.sum()),
            "n_lost_still_missing_last_month": still_missing_at_end,
            "lost_median_entry_month": entry_label}


def compute_stats(mat: dict, months: list, monthly: dict, supplement: dict) -> dict:
    """Collect every number used in the appendix text."""
    active = ~np.isnan(mat["frac"])
    n_per_month = active.sum(axis=0)
    labels = [ym_label(months, i) for i in range(len(months))]
    avail = mat["span_frac"]
    from_main = mat["source_ym"] == "main"

    stats = {
        "main_archive": MAIN_ROOT,
        "window": [labels[0], labels[-1]],
        "main_archive_window": [labels[int(np.argmax(from_main.any(axis=0)))],
                                labels[len(months) - 1 - int(np.argmax(
                                    from_main.any(axis=0)[::-1]))]],
        "n_months": len(months),
        "n_instruments": len(mat["keys"]),
        "n_sites_wigos": len({k[0] for k in mat["keys"]}),
        "n_instrument_months": int(active.sum()),
        "n_instrument_hours": int(np.nansum(np.nan_to_num(mat["frac"]) * mat["nhours"])),
        "share_instrument_months_from_supplement_pct": round(
            100 * float((mat["source_ym"] == "supplement").sum() / active.sum()), 2),
        "instruments_active_per_month": {
            lab: int(n_per_month[labels.index(lab)])
            for lab in ("2016-10", "2016-12", "2018-01", "2020-01", "2022-01", "2024-01",
                        "2025-11", labels[-1]) if lab in labels},
        "instruments_active_max": int(n_per_month.max()),
        "instruments_active_max_month": labels[int(n_per_month.argmax())],
        "network_median_availability_pct": round(100 * float(np.median(avail)), 1),
        "share_instruments_above_90pct": round(100 * float(np.mean(avail >= 0.9)), 1),
        "share_instruments_above_95pct": round(100 * float(np.mean(avail >= 0.95)), 1),
        "mean_monthly_hour_fraction_active_pct": round(
            100 * float(np.nanmean(mat["frac"])), 1),
        "share_active_months_above_90pct": round(
            100 * float(np.mean(mat["frac"][active] >= 0.9)), 1),
        "dips": find_network_dips(n_per_month, months),
        "supplement": supplement,
        "by_type": {},
    }
    for itype in TYPE_ORDER:
        sel = mat["row_type"] == itype
        counts = monthly[itype]["count"]
        if not sel.any():
            continue
        first = int(np.argmax(counts > 0))
        stats["by_type"][itype] = {
            "n_instruments_ever": int(sel.sum()),
            "n_active_2016-12": int(counts[labels.index("2016-12")]),
            "n_active_last_month": int(counts[-1]),
            "max_active": int(counts.max()),
            "first_month_with_data": labels[first],
            "median_availability_pct": round(100 * float(np.median(avail[sel])), 1),
            "iqr_availability_pct": [round(100 * float(q), 1)
                                     for q in np.percentile(avail[sel], [25, 75])],
        }
    # Which countries account for the larger drops (to tell outages of a whole national
    # network from scattered instrument failures).
    for dip in stats["dips"]:
        dip.update(characterise_dip(mat, months, labels.index(dip["month"])))
    return stats


# %% Figure


def availability_rgb(mat: dict, order: np.ndarray) -> np.ndarray:
    """RGB image: type colour of that month, tinted by availability class; white = no data."""
    frac = mat["frac"][order]
    types = mat["type_ym"][order]
    img = np.ones(frac.shape + (3,))
    cls = np.digitize(np.nan_to_num(frac, nan=-1.0), AVAIL_BINS) - 1   # -1 = no data
    for itype in TYPE_ORDER:
        base = np.array(to_rgb(TYPE_COLOURS[itype]))
        for k, tint in enumerate(AVAIL_TINTS):
            img[(types == itype) & (cls == k)] = 1.0 - tint * (1.0 - base)
    return img


def plot_figure(mat: dict, months: list, monthly: dict, stats: dict, out_png: Path) -> None:
    """Landscape figure: (a) matrix | (b) instrument counts, (c) monthly availability, (d) box."""
    order = order_rows(mat)
    img = availability_rgb(mat, order)
    n_rows, n_months = img.shape[:2]

    if PANEL_A_ONLY:
        fig = plt.figure(figsize=(16, 8.6), dpi=200)
        gs = GridSpec(1, 1, figure=fig, left=0.095, right=0.985, top=0.915, bottom=0.07)
        ax_m = fig.add_subplot(gs[0, 0])
    else:
        fig = plt.figure(figsize=(18, 9.2), dpi=200)
        gs = GridSpec(3, 2, figure=fig, width_ratios=[3.3, 1.0],
                      wspace=0.11, hspace=0.38, left=0.085, right=0.985, top=0.925, bottom=0.065)
        ax_m = fig.add_subplot(gs[:, 0])
        ax_n = fig.add_subplot(gs[0, 1])
        ax_f = fig.add_subplot(gs[1, 1], sharex=ax_n)
        ax_b = fig.add_subplot(gs[2, 1])

    # (a) the matrix. Types with few instruments (CL61, Mini-MPL) get taller rows so that their
    # block stays legible; the vertical axis carries no quantity, only one row per instrument.
    row_types = mat["row_type"][order]
    type_n = Counter(row_types)
    heights = np.array([SMALL_GROUP_ROW_HEIGHT if type_n[t] < SMALL_GROUP_SIZE else 1
                        for t in row_types])
    img_rows = np.repeat(img, heights, axis=0)
    n_disp = img_rows.shape[0]
    ax_m.imshow(img_rows, aspect="auto", interpolation="nearest",
                extent=(0, n_months, n_disp, 0))
    y0 = 0
    for itype in TYPE_ORDER:
        sel = row_types == itype
        n = int(sel.sum())
        if n == 0:
            continue
        h = int(heights[sel].sum())
        ax_m.axhline(y0, color="0.3", lw=0.6)
        label = f"{itype}\n(n = {n})" if h > 40 else f"{itype} (n = {n})"
        ax_m.text(-0.8, y0 + h / 2, label, ha="right", va="center",
                  fontsize=9.5, color=TYPE_COLOURS[itype], fontweight="bold")
        y0 += h
    ax_m.set_ylim(n_disp, 0)
    ax_m.set_yticks([])
    year_ticks = [i for i, (_y, m) in enumerate(months) if m == 1]
    ax_m.set_xticks(year_ticks)
    ax_m.set_xticklabels([str(months[i][0]) for i in year_ticks], fontsize=9)
    for i in year_ticks:
        ax_m.axvline(i, color="0.6", lw=0.4, ls=":", zorder=3)
    ax_m.set_xlim(0, n_months)
    ax_m.set_xlabel("Year", fontsize=10)
    ax_m.set_ylabel("Instruments (one row each)", fontsize=10, labelpad=78)
    ax_m.set_title(("" if PANEL_A_ONLY else "(a) ") + f"Fraction of hours per month with data, per instrument "
                   f"({stats['n_instruments']} instruments at {stats['n_sites_wigos']} sites, "
                   f"{stats['window'][0]} to {stats['window'][1]})", fontsize=10.5, loc="left",
                   pad=22)
    grey = np.array(to_rgb("#404040"))
    shade_handles = [Patch(facecolor=1.0 - t * (1.0 - grey), edgecolor="0.4", label=lab)
                     for t, lab in zip(AVAIL_TINTS, AVAIL_LABELS)]
    shade_handles.append(Patch(facecolor="white", edgecolor="0.4", label="no data"))
    ax_m.legend(handles=shade_handles, title=None, loc="lower right",
                bbox_to_anchor=(1.0, 1.0), fontsize=8.5, ncol=4, frameon=False,
                borderaxespad=0.1, handlelength=1.6, columnspacing=1.2)

    if not PANEL_A_ONLY:
        # (b) instruments delivering per month, per type (lines; Fig. 1 already shows the stack)
        xm = np.arange(n_months) + 0.5
        for itype in stats["by_type"]:
            ax_n.plot(xm, monthly[itype]["count"], color=TYPE_COLOURS[itype], lw=1.4, label=itype)
        ax_n.set_ylabel("Instruments", fontsize=9)
        ax_n.set_title("(b) Instruments delivering data, per month", fontsize=10, loc="left")
        ax_n.legend(fontsize=8, ncol=3, loc="upper left", frameon=False)

        # (c) mean fraction of hours with data of the active instruments of each type
        for itype in stats["by_type"]:
            series = np.where(monthly[itype]["count"] >= MIN_ACTIVE_FOR_MEAN,
                              100 * monthly[itype]["mean_frac"],
                              np.nan)  # 1-2 instruments only: too noisy to be meaningful
            ax_f.plot(xm, series, color=TYPE_COLOURS[itype], lw=1.2)
        ax_f.set_ylim(40, 101)
        ax_f.set_ylabel("Hours with data (%)", fontsize=9)
        ax_f.set_title("(c) Mean monthly availability of active instruments", fontsize=10,
                       loc="left")
        for ax in (ax_n, ax_f):
            ax.grid(alpha=0.3, lw=0.5)
            ax.set_xticks(year_ticks[::2])
            ax.set_xticklabels([str(months[i][0]) for i in year_ticks[::2]], fontsize=8)
            ax.set_xlim(0, n_months)
            ax.tick_params(labelsize=8)

        # (d) per-instrument availability over its own operating period, by type
        data, labels, colours = [], [], []
        for itype in stats["by_type"]:
            sel = mat["row_type"] == itype
            data.append(100 * mat["span_frac"][sel])
            labels.append(itype)
            colours.append(TYPE_COLOURS[itype])
        parts = ax_b.boxplot(data, vert=False, patch_artist=True, widths=0.6,
                             flierprops={"marker": ".", "markersize": 3, "alpha": 0.5},
                             medianprops={"color": "black", "lw": 1.2})
        for patch, colour in zip(parts["boxes"], colours):
            patch.set_facecolor(colour)
            patch.set_alpha(0.6)
        ax_b.set_yticks(range(1, len(labels) + 1))
        ax_b.set_yticklabels(labels, fontsize=8)
        ax_b.invert_yaxis()
        ax_b.set_xlim(0, 101)
        ax_b.set_xlabel("Hours with data between first and last record (%)", fontsize=9)
        ax_b.set_title("(d) Per-instrument availability", fontsize=10, loc="left")
        ax_b.grid(axis="x", alpha=0.3, lw=0.5)
        ax_b.tick_params(axis="x", labelsize=8)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Figure written: {out_png}")


# %% Text


def _month_name(label: str) -> str:
    """'2016-12' -> 'December 2016'."""
    names = ["January", "February", "March", "April", "May", "June", "July", "August",
             "September", "October", "November", "December"]
    return f"{names[int(label[5:]) - 1]} {label[:4]}"


def _describe_dips(stats: dict) -> str:
    """Sentences on the network-wide drops: temporary (all countries) vs persistent."""
    first = stats["window"][0]
    dips = [d for d in stats["dips"] if d["month"] != first]
    temporary = [d for d in dips if d["temporary"]]
    persistent = [d for d in dips if not d["temporary"]]
    parts = []
    if temporary:
        months = " and ".join(f"{_month_name(d['month'])} ({d['from']} to {d['to']} "
                              f"instruments)" for d in temporary)
        plural = len(temporary) > 1
        parts.append(
            f"The network-wide interruption{'s' if plural else ''} of {months} "
            f"{'affect' if plural else 'affects'} all countries at once, the instruments "
            f"delivering again in the following months; {'they are' if plural else 'it is'} "
            f"most likely {'gaps' if plural else 'a gap'} in the archive rather than "
            f"instrument outages, although this cannot be verified.")
    for d in persistent:
        country, n = next(iter(d["countries_lost"].items()))
        parts.append(
            f"The drop in {_month_name(d['month'])} ({d['from']} to {d['to']} instruments) is "
            f"due to {n} instruments in {country}, most of them delivering since "
            f"{_month_name(d['lost_median_entry_month'])}, which stopped delivering data and had "
            f"not resumed by {_month_name(stats['window'][1])}.")
    return " ".join(parts)


def write_text(stats: dict, out_md: Path) -> None:
    """Write the appendix paragraphs + caption; every number comes from `stats`."""
    bt = stats["by_type"]
    act = stats["instruments_active_per_month"]
    sup = stats["supplement"]
    first, last = stats["window"]
    main_first = stats["main_archive_window"][0]
    type_counts = ", ".join(f"{bt[t]['n_instruments_ever']} {t}" for t in bt)
    medians = ", ".join(f"{t} {bt[t]['median_availability_pct']:.0f} %" for t in bt)
    cl61 = bt["CL61"]

    para1 = (
        f"Figure A2-1 shows the availability of E-Profile ALC data from {_month_name(first)} to "
        f"{_month_name(last)}, derived from the archive of monthly files of hourly-mean Level-2 "
        f"data: an hour is available when the hub produced Level-2 data for it, i.e. when "
        f"Level-1 data were received and processed. Before {_month_name(main_first)}, and for "
        f"{sup['n_station_months_added']} station-months missing from that archive, the hourly "
        f"coverage of the Level-1 daily files is used instead. Over the period, "
        f"{stats['n_instruments']} instruments at {stats['n_sites_wigos']} sites delivered data "
        f"({type_counts}). The number of instruments delivering data grew from "
        f"{act['2016-12']} in December 2016 to {act['2020-01']} in January 2020 and peaked at "
        f"{stats['instruments_active_max']} in "
        f"{_month_name(stats['instruments_active_max_month'])} ({act[last]} in "
        f"{_month_name(last)}). The first CL61 entered the network in "
        f"{_month_name(cl61['first_month_with_data'])}; {cl61['n_active_last_month']} were "
        f"active in {_month_name(last)}."
    )
    para2 = (
        f"Once an instrument has joined the network its record is largely continuous: the "
        f"median fraction of hours with data between its first and last record is "
        f"{stats['network_median_availability_pct']:.0f} % ({medians}), and "
        f"{stats['share_active_months_above_90pct']:.0f} % of the instrument-months have data "
        f"for more than 90 % of their hours. Most gaps concern single instruments. "
        f"{_describe_dips(stats)}"
    )
    caption = (
        f"Figure A2-1. Data availability of the E-Profile ALC network, {_month_name(first)} to "
        f"{_month_name(last)}, from the archive of hourly Level-2 data (completed before "
        f"{_month_name(main_first)} with Level-1 daily files). "
        f"(a) Fraction of the hours of each month with data for each of the "
        f"{stats['n_instruments']} instruments (one row per instrument), grouped by instrument "
        f"type and sorted within each type by date of entry into the network (rows of the types "
        f"with fewer than {SMALL_GROUP_SIZE} instruments drawn {SMALL_GROUP_ROW_HEIGHT} times "
        f"taller for legibility). The colour gives "
        f"the instrument type reported in that month (CHM15k red, CL31 orange, CL51 purple, CL61 "
        f"blue, Mini-MPL green) and its shade the fraction of hours with data (< 50 %, 50-90 %, "
        f"> 90 %); white: no data. (b) Number of instruments delivering data in each month, by "
        f"type. (c) Mean monthly fraction of hours with data of the active instruments of each "
        f"type (shown when at least {MIN_ACTIVE_FOR_MEAN} are active). (d) Distribution, by type, "
        f"of the "
        f"fraction of hours with data between the first and last record of each instrument "
        f"(box: interquartile range and median; whiskers: 1.5 times the interquartile range)."
    )
    if PANEL_A_ONLY:
        caption = caption.split(" (b) Number of instruments")[0].replace(
            ". (a) Fraction of the hours", ". Fraction of the hours")
    out_md.write_text(
        "# Appendix 2: Data availability\n\n"
        f"{para1}\n\n{para2}\n\n"
        "![Figure A2-1](fig_A2_data_availability.png)\n\n"
        f"{caption}\n", encoding="utf-8")
    print(f"Text written: {out_md}")


# %% Main


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def main() -> None:
    parser = argparse.ArgumentParser(description="Appendix 2 data-availability figure")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--rescan", action="store_true", help="ignore the read cache")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--all-panels", action="store_true",
                        help="draw panels (b)-(d) next to the availability matrix (default: matrix only)")
    args = parser.parse_args()
    global PANEL_A_ONLY
    PANEL_A_ONLY = not args.all_panels

    records, supplement = load_records(args.out / "cache", args.rescan, args.workers)
    streams = group_streams(records)
    fill_unknown_types(streams)

    # Window: from the earliest record (supplement) to the last month of the MAIN archive; the
    # supplement holds a few isolated later station-months that would end the plot on an
    # almost empty month.
    first_month = min(r["ym"] for r in records)
    last_month = max(r["ym"] for r in records if r["source"] == "main")
    months = month_list(first_month, last_month)
    mat = build_matrix(streams, months)
    monthly = monthly_by_type(mat)
    stats = compute_stats(mat, months, monthly, supplement)

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "appendix2_stats.json").write_text(json.dumps(_jsonable(stats), indent=1),
                                                   encoding="utf-8")
    plot_figure(mat, months, monthly, stats, args.out / "fig_A2_data_availability.png")
    write_text(stats, args.out / "appendix2_text.md")


if __name__ == "__main__":
    main()
