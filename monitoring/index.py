"""Build a single queryable SQLite index from the per-station calibration CSVs.

Tables:
  calibrations -- one row per (station, method, date): flag, cal_value, uncertainty, n_profiles, window, message
  series       -- one row per (station, method): aggregates (success rate, median, last value/flag, ...)
  stations     -- one row per instrument: manifest/L2 metadata (lat/lon/type/name/country/institution)
  kalman       -- precomputed E-PROFILE Kalman best estimate per (station, method, date), if available

A "series" = one calibration method on one instrument (e.g. CL61 has a rayleigh series AND a
cloud series). The renderer reads ONLY this index, so rebuilding it is the single refresh step.

Supports both schemas:
  * unified  <key>_cal.csv : date, method, flag, cal_value, uncertainty, n_profiles, bottom_height, top_height, message
  * legacy   <key>_cl.csv  : date, flag, lidar_constant, uncertainty, bottom_height, top_height, message  (method assumed 'rayleigh')
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from monitoring.config import DB_NAME, SUCCESS_FLAGS

_CAL_COLS = ["key", "method", "date", "datetime", "flag", "success", "cal_value",
             "uncertainty", "rel_uncertainty", "n_profiles", "bottom_height", "top_height", "message"]


def _read_station_csv(path: Path, key: str) -> pd.DataFrame:
    """Load one calibration CSV (unified or legacy) into the normalized column set."""
    df = pd.read_csv(path, dtype={"date": str})
    # Legacy <key>_cl.csv: no method column, value column is 'lidar_constant'.
    if "method" not in df.columns:
        df["method"] = "rayleigh"
    if "cal_value" not in df.columns and "lidar_constant" in df.columns:
        df["cal_value"] = df["lidar_constant"]
    if "n_profiles" not in df.columns:
        df["n_profiles"] = np.nan

    for col in ("flag", "cal_value", "uncertainty", "bottom_height", "top_height", "n_profiles"):
        df[col] = pd.to_numeric(df.get(col), errors="coerce")
    df["method"] = df["method"].fillna("rayleigh").astype(str)
    df["message"] = df.get("message", "").fillna("")
    df["key"] = key
    df["datetime"] = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    df = df.dropna(subset=["datetime"]).copy()
    # A calibration counts as a success only if its flag is a success flag AND it produced a
    # positive value. This guards against rows where the flag says success but the value is
    # invalid (seen on a few CL61 cloud days: cloud_flag passes on cal_median while the derived
    # lidar_constant came out -1) -- those must not pollute the time series / median / counts.
    df["success"] = (df["flag"].isin(SUCCESS_FLAGS) & (df["cal_value"] > 0)).astype(int)
    df["rel_uncertainty"] = np.where(df["success"] == 1,
                                     100.0 * df["uncertainty"] / df["cal_value"], np.nan)
    return df[_CAL_COLS]


def _read_kalman_csv(path: Path, key: str) -> pd.DataFrame:
    """Load one <key>_kalman.csv (method, date, kalman, kalman_std) -> normalized rows."""
    if not path.exists():
        return pd.DataFrame(columns=["key", "method", "date", "kalman", "kalman_std"])
    df = pd.read_csv(path, dtype={"date": str})
    if df.empty:
        return pd.DataFrame(columns=["key", "method", "date", "kalman", "kalman_std"])
    df["key"] = key
    df["method"] = df.get("method", "rayleigh").fillna("rayleigh").astype(str)
    for col in ("kalman", "kalman_std"):
        df[col] = pd.to_numeric(df.get(col), errors="coerce")
    return df[["key", "method", "date", "kalman", "kalman_std"]]


def _scan_diagnostics(fullcal_dir: Path, keys) -> pd.DataFrame:
    """Find per-calibration diagnostic PNGs the runner emitted under <key>/plots/**/*.png.

    Filenames are '<YYYYMMDD>_<wmo>_<method>_diag_compact.png'; we key each image by
    (key, method, date) and keep its absolute source path for the renderer to copy in."""
    rows = []
    for key in keys:
        pdir = Path(fullcal_dir) / key / "plots"
        if not pdir.exists():
            continue
        for png in pdir.rglob("*_diag_compact.png"):
            name = png.name
            date = name[:8]
            if not date.isdigit():
                continue
            method = "cloud" if "_cloud_diag" in name else ("rayleigh" if "_rayleigh_diag" in name else None)
            if method is None:
                continue
            rows.append(dict(key=key, method=method, date=date, src=str(png)))
    cols = ["key", "method", "date", "src"]
    return pd.DataFrame(rows, columns=cols)


def _load_manifest(manifest_path: Path) -> pd.DataFrame:
    """Per-instrument metadata keyed by '<wmo>_<identifier>'. Handles the L2 manifest
    ({wmo, identifier, itype, ...}) and the L1 census ({wmo, ident, type, site, ...})."""
    cols = ["key", "wmo", "identifier", "itype", "lat", "lon", "alt"]
    if not manifest_path.exists():
        return pd.DataFrame(columns=cols)
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = []
    for s in data:
        ident = s.get("identifier", s.get("ident"))
        rows.append(dict(
            key=f"{s['wmo']}_{ident}", wmo=s.get("wmo"), identifier=ident,
            itype=s.get("itype", s.get("type")), lat=s.get("lat"), lon=s.get("lon"),
            alt=s.get("alt"), name=s.get("site", ""),
        ))
    return pd.DataFrame(rows)


def _enrich_metadata(stations: pd.DataFrame, l2_dir: Path) -> pd.DataFrame:
    """Add station name / country / institution from L2 NetCDF global attributes (best-effort)."""
    stations = stations.copy()
    for c in ("name", "country", "institution"):
        if c not in stations.columns:
            stations[c] = ""
        stations[c] = stations[c].fillna("")
    if not Path(l2_dir).exists():
        return stations
    try:
        import netCDF4  # noqa: F401
    except ImportError:
        return stations
    for i, s in stations.iterrows():
        if str(s.get("country") or "") and str(s.get("institution") or ""):
            continue  # already have it (e.g. from a census 'site')
        wmo, ident = str(s["wmo"]), str(s["identifier"])
        if not wmo or wmo == "nan":
            wmo, _, ident = str(s["key"]).rpartition("_")
        matches = sorted(Path(l2_dir).glob(f"{wmo}/*/L2_{wmo}_{ident}*.nc"))
        if not matches:
            continue
        try:
            with netCDF4.Dataset(matches[-1]) as d:
                site = str(getattr(d, "site_location", "") or "")
                inst = str(getattr(d, "institution", "") or "")
        except OSError:
            continue
        name, _, country = site.partition(",")
        if name.strip():
            stations.at[i, "name"] = name.strip()
        stations.at[i, "country"] = country.strip()
        stations.at[i, "institution"] = inst.strip()
    return stations


def _series_aggregates(cal: pd.DataFrame) -> pd.DataFrame:
    """Per-(station, method) summary rows derived from the calibration history in range."""
    records = []
    for (key, method), g in cal.groupby(["key", "method"], sort=False):
        g = g.sort_values("datetime")
        ok = g[g["success"] == 1]
        last = g.iloc[-1]
        # Success rate counts ATTEMPTABLE days: exclude only no-data (0) and unsuitable conditions
        # (-1 = cloudy night for Rayleigh / clear sky for cloud -- genuinely no opportunity). Every
        # other non-success, including the cloud rejections (-20..-26: dirty window, low laser energy,
        # cloud not a clean stratocumulus), counts as a FAILURE -- so the rate measures "of the days a
        # calibration could be attempted, how often it succeeded", consistent with both the headline
        # KPI (metrics.network_summary) and how Rayleigh already counts its own quality failures.
        n_suitable = int((~g["flag"].isin([0, -1])).sum())
        records.append(dict(
            key=key, method=method,
            n_dates=int(len(g)), n_success=int(len(ok)), n_suitable=n_suitable,
            success_rate=float(100.0 * len(ok) / n_suitable) if n_suitable else float("nan"),
            median_cl=float(ok["cal_value"].median()) if len(ok) else float("nan"),
            median_rel_unc=float(ok["rel_uncertainty"].median()) if len(ok) else float("nan"),
            first_date=str(g.iloc[0]["date"]), last_date=str(last["date"]),
            last_flag=float(last["flag"]) if pd.notna(last["flag"]) else None,
            last_success_date=str(ok.iloc[-1]["date"]) if len(ok) else None,
            last_cl=float(ok.iloc[-1]["cal_value"]) if len(ok) else float("nan"),
        ))
    return pd.DataFrame(records, columns=[
        "key", "method", "n_dates", "n_success", "n_suitable", "success_rate", "median_cl",
        "median_rel_unc", "first_date", "last_date", "last_flag", "last_success_date", "last_cl"])


def build_index(
    fullcal_dir: Path,
    manifest_path: Path,
    db_path: Path,
    limit: Optional[int] = None,
    types: Optional[Iterable[str]] = None,
    l2_dir: Optional[Path] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
) -> dict:
    """Scan the calibration CSVs and (re)write the SQLite index. ``start``/``end`` (YYYYMMDD)
    restrict the calibrations to a date window. Returns a small stats dict."""
    fullcal_dir = Path(fullcal_dir)
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = _load_manifest(Path(manifest_path))
    type_by_key = dict(zip(manifest["key"], manifest["itype"])) if len(manifest) else {}

    csvs = sorted(list(fullcal_dir.glob("*/*_cal.csv")) + list(fullcal_dir.glob("*/*_cl.csv")))
    if not csvs:
        raise FileNotFoundError(f"No '*_cal.csv' or '*_cl.csv' found under {fullcal_dir}")

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
    cal = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=_CAL_COLS)
    if start:
        cal = cal[cal["date"] >= start]
    if end:
        cal = cal[cal["date"] <= end]

    kal_cols = ["key", "method", "date", "kalman", "kalman_std"]
    kal_frames = [_read_kalman_csv(csv.with_name(csv.name.replace("_cal.csv", "_kalman.csv")
                                                 .replace("_cl.csv", "_kalman.csv")), key)
                  for csv, key in selected]
    kal_frames = [k for k in kal_frames if len(k)]  # drop empties (avoids concat dtype warning)
    kalman = pd.concat(kal_frames, ignore_index=True) if kal_frames else pd.DataFrame(columns=kal_cols)
    if len(kalman) and start:
        kalman = kalman[kalman["date"] >= start]
    if len(kalman) and end:
        kalman = kalman[kalman["date"] <= end]

    series = _series_aggregates(cal)
    # Stations table = per-key metadata; attach itype (also infer from series if manifest lacks it).
    keys = pd.DataFrame({"key": [k for _, k in selected]}).drop_duplicates()
    stations = keys.merge(manifest, on="key", how="left")
    stations["itype"] = stations["itype"].fillna("Unknown")
    series = series.merge(stations[["key", "itype"]], on="key", how="left")
    if l2_dir is not None:
        stations = _enrich_metadata(stations, Path(l2_dir))

    diagnostics = _scan_diagnostics(fullcal_dir, [k for _, k in selected])
    if len(diagnostics) and start:
        diagnostics = diagnostics[diagnostics["date"] >= start]
    if len(diagnostics) and end:
        diagnostics = diagnostics[diagnostics["date"] <= end]

    with sqlite3.connect(db_path) as con:
        cal.drop(columns=["datetime"]).to_sql("calibrations", con, if_exists="replace", index=False)
        series.to_sql("series", con, if_exists="replace", index=False)
        stations.to_sql("stations", con, if_exists="replace", index=False)
        kalman.to_sql("kalman", con, if_exists="replace", index=False)
        diagnostics.to_sql("diagnostics", con, if_exists="replace", index=False)
        con.execute("CREATE INDEX IF NOT EXISTS ix_cal_key ON calibrations(key, method)")
        con.execute("CREATE INDEX IF NOT EXISTS ix_kal_key ON kalman(key, method)")
        con.commit()

    return dict(
        n_stations=int(len(stations)), n_series=int(len(series)),
        n_calibrations=int(len(cal)),
        as_of=str(cal["date"].max()) if len(cal) else None,
        date_min=str(cal["date"].min()) if len(cal) else None,
        db_path=str(db_path),
    )
