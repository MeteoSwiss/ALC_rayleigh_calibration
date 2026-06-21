"""Build a single queryable SQLite index from the per-station calibration CSVs.

Two tables:
  calibrations  -- one row per (station, night): flag, lidar_constant, uncertainty, window, message
  stations      -- one row per instrument: manifest metadata + per-station aggregates

The renderer reads ONLY this index, so rebuilding it is the single 'refresh' step of the
daily dashboard. ~200 stations x ~1-3k nights is trivial for SQLite; no server needed.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from monitoring.config import DB_NAME, SUCCESS_FLAGS

# Columns written by run_all_l2monthly.py for each <key>_cl.csv.
_NUMERIC_COLS = ("flag", "lidar_constant", "uncertainty", "bottom_height", "top_height")


def _read_station_csv(path: Path, key: str) -> pd.DataFrame:
    """Load one <key>_cl.csv, coercing numerics and parsing the YYYYMMDD date."""
    df = pd.read_csv(path, dtype={"date": str})
    for col in _NUMERIC_COLS:
        # Failed nights write empty height fields and lidar_constant=-1; coerce robustly.
        df[col] = pd.to_numeric(df.get(col), errors="coerce")
    if "message" not in df.columns:
        df["message"] = ""
    df["message"] = df["message"].fillna("")
    df["key"] = key
    df["datetime"] = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    df = df.dropna(subset=["datetime"]).copy()
    df["success"] = df["flag"].isin(SUCCESS_FLAGS).astype(int)
    # Relative uncertainty (%) only meaningful on successful nights with a positive C_L.
    good = (df["success"] == 1) & (df["lidar_constant"] > 0)
    df["rel_uncertainty"] = np.where(good, 100.0 * df["uncertainty"] / df["lidar_constant"], np.nan)
    return df[["key", "date", "datetime", "flag", "success", "lidar_constant",
               "uncertainty", "rel_uncertainty", "bottom_height", "top_height", "message"]]


def _load_manifest(manifest_path: Path) -> pd.DataFrame:
    """Per-instrument metadata (lat/lon/type/...) keyed by '<wmo>_<identifier>'."""
    if not manifest_path.exists():
        return pd.DataFrame(columns=["key", "wmo", "identifier", "itype", "lat", "lon", "alt", "n_months"])
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = [
        dict(key=f"{s['wmo']}_{s['identifier']}", wmo=s.get("wmo"), identifier=s.get("identifier"),
             itype=s.get("itype"), lat=s.get("lat"), lon=s.get("lon"),
             alt=s.get("alt"), n_months=s.get("n_months"))
        for s in data
    ]
    return pd.DataFrame(rows)


def _enrich_metadata(stations: pd.DataFrame, l2_dir: Path) -> pd.DataFrame:
    """Add station name / country / institution by reading L2 NetCDF global attributes.

    E-PROFILE L2 files carry `site_location` ('NAME,COUNTRY') and `institution`. We read
    one file per instrument (header only) and cache nothing; on a missing L2 archive the
    fields stay blank and the rest of the dashboard is unaffected.
    """
    stations = stations.copy()
    stations["name"] = ""
    stations["country"] = ""
    stations["institution"] = ""
    if not Path(l2_dir).exists():
        return stations
    try:
        import netCDF4  # noqa: F401
    except ImportError:
        return stations

    for i, s in stations.iterrows():
        wmo, ident = str(s["wmo"]), str(s["identifier"])
        if not wmo or wmo == "nan":
            wmo, _, ident = str(s["key"]).rpartition("_")  # fallback: split the key
        matches = sorted(Path(l2_dir).glob(f"{wmo}/*/L2_{wmo}_{ident}*.nc"))
        if not matches:
            continue
        try:
            with netCDF4.Dataset(matches[-1]) as d:  # newest month: freshest metadata
                site = str(getattr(d, "site_location", "") or "")
                inst = str(getattr(d, "institution", "") or "")
        except OSError:
            continue
        name, _, country = site.partition(",")
        stations.at[i, "name"] = name.strip()
        stations.at[i, "country"] = country.strip()
        stations.at[i, "institution"] = inst.strip()
    return stations


def _station_aggregates(cal: pd.DataFrame) -> pd.DataFrame:
    """Per-station summary rows derived from the full calibration history."""
    records = []
    for key, g in cal.groupby("key", sort=False):
        g = g.sort_values("datetime")
        ok = g[g["success"] == 1]
        last = g.iloc[-1]
        records.append(dict(
            key=key,
            n_dates=int(len(g)),
            n_success=int(len(ok)),
            success_rate=float(100.0 * len(ok) / len(g)) if len(g) else float("nan"),
            median_cl=float(ok["lidar_constant"].median()) if len(ok) else float("nan"),
            median_rel_unc=float(ok["rel_uncertainty"].median()) if len(ok) else float("nan"),
            first_date=str(g.iloc[0]["date"]),
            last_date=str(last["date"]),
            last_flag=float(last["flag"]) if pd.notna(last["flag"]) else None,
            last_success_date=str(ok.iloc[-1]["date"]) if len(ok) else None,
            last_cl=float(ok.iloc[-1]["lidar_constant"]) if len(ok) else float("nan"),
        ))
    return pd.DataFrame(records)


def build_index(
    fullcal_dir: Path,
    manifest_path: Path,
    db_path: Path,
    limit: Optional[int] = None,
    types: Optional[Iterable[str]] = None,
    l2_dir: Optional[Path] = None,
) -> dict:
    """Scan the calibration CSVs and (re)write the SQLite index. Returns a small stats dict."""
    fullcal_dir = Path(fullcal_dir)
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = _load_manifest(Path(manifest_path))
    type_by_key = dict(zip(manifest["key"], manifest["itype"])) if len(manifest) else {}

    csvs = sorted(fullcal_dir.glob("*/*_cl.csv"))
    if not csvs:
        raise FileNotFoundError(f"No '*_cl.csv' found under {fullcal_dir}")

    want = set(types) if types else None
    selected = []
    for csv in csvs:
        key = csv.parent.name
        if want is not None and type_by_key.get(key) not in want:
            continue
        selected.append((csv, key))
    if limit:
        selected = selected[:limit]

    frames = [_read_station_csv(csv, key) for csv, key in selected]
    cal = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    agg = _station_aggregates(cal)
    # Left-join aggregates onto manifest metadata so unmatched stations still appear.
    keys = pd.DataFrame({"key": [k for _, k in selected]})
    stations = keys.merge(manifest, on="key", how="left").merge(agg, on="key", how="left")
    # itype fallback for stations missing from the manifest.
    stations["itype"] = stations["itype"].fillna("Unknown")

    # Enrich with name/country/institution from the L2 archive (best-effort).
    if l2_dir is not None:
        stations = _enrich_metadata(stations, Path(l2_dir))

    with sqlite3.connect(db_path) as con:
        cal_out = cal.drop(columns=["datetime"]) if "datetime" in cal.columns else cal
        cal_out.to_sql("calibrations", con, if_exists="replace", index=False)
        stations.to_sql("stations", con, if_exists="replace", index=False)
        con.execute("CREATE INDEX IF NOT EXISTS ix_cal_key ON calibrations(key)")
        con.execute("CREATE INDEX IF NOT EXISTS ix_cal_date ON calibrations(date)")
        con.commit()

    return dict(
        n_stations=int(len(stations)),
        n_calibrations=int(len(cal)),
        as_of=str(cal["date"].max()) if len(cal) else None,
        db_path=str(db_path),
    )
