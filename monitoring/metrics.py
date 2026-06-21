"""Derived metrics read back from the SQLite index: network summary + watchlist.

Kept separate from index-building so thresholds can be re-tuned and the dashboard
re-rendered without re-scanning the CSVs.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from monitoring import config


def load_frames(db_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (calibrations, stations) DataFrames with a parsed `datetime` on calibrations."""
    with sqlite3.connect(Path(db_path)) as con:
        cal = pd.read_sql("SELECT * FROM calibrations", con)
        st = pd.read_sql("SELECT * FROM stations", con)
    cal["datetime"] = pd.to_datetime(cal["date"], format="%Y%m%d", errors="coerce")
    return cal, st


def network_summary(cal: pd.DataFrame, st: pd.DataFrame) -> dict:
    """Headline KPIs for the summary page."""
    n_ok = int(cal["success"].sum())
    n_total = int(len(cal))
    as_of = str(cal["date"].max()) if n_total else None

    by_type = (
        st.groupby("itype")
        .agg(n_stations=("key", "count"),
             n_dates=("n_dates", "sum"),
             n_success=("n_success", "sum"))
        .reset_index()
    )
    by_type["success_rate"] = 100.0 * by_type["n_success"] / by_type["n_dates"].replace(0, np.nan)
    by_type = by_type.sort_values("itype")

    return dict(
        as_of=as_of,
        n_stations=int(len(st)),
        n_calibrations=n_total,
        n_success=n_ok,
        success_rate=(100.0 * n_ok / n_total) if n_total else float("nan"),
        date_min=str(cal["date"].min()) if n_total else None,
        date_max=as_of,
        by_type=by_type,
    )


def flag_distribution(cal: pd.DataFrame) -> pd.DataFrame:
    """Counts per flag value, labelled and ordered best-to-worst."""
    counts = cal["flag"].value_counts(dropna=False).reset_index()
    counts.columns = ["flag", "count"]
    counts["label"] = counts["flag"].map(config.flag_label)
    counts["color"] = counts["flag"].map(config.flag_color)
    # Order: 1, 0.5, then descending (0, -1, -2, ...) so 'good' sits on top of charts.
    counts["order"] = counts["flag"].map(lambda f: -float(f) if pd.notna(f) else 999)
    return counts.sort_values("order").drop(columns="order").reset_index(drop=True)


def _robust_sigma(x: np.ndarray) -> float:
    """MAD-based robust standard-deviation estimate (1.4826 * MAD)."""
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def watchlist(cal: pd.DataFrame, st: pd.DataFrame) -> pd.DataFrame:
    """Stations worth a look this morning: C_L drift, failure streaks, low recent success.

    Returns one row per (station, issue) with a severity-orderable `priority`.
    """
    alerts = []
    for key, g in cal.groupby("key", sort=False):
        g = g.sort_values("datetime")
        meta = st[st["key"] == key]
        itype = meta["itype"].iloc[0] if len(meta) else "Unknown"
        ok = g[g["success"] == 1]

        # 1) C_L drift: latest successful constant far from the station's robust baseline.
        if len(ok) >= 8:
            cl = ok["lidar_constant"].to_numpy(dtype=float)
            med = float(np.median(cl))
            sig = _robust_sigma(cl)
            last_cl = cl[-1]
            if np.isfinite(sig) and sig > 0:
                z = abs(last_cl - med) / sig
                if z >= config.DRIFT_SIGMA:
                    alerts.append(dict(key=key, itype=itype, issue="C_L drift",
                                       detail=f"last {last_cl:.3g} vs median {med:.3g} ({z:.1f}σ)",
                                       date=str(ok.iloc[-1]["date"]), priority=float(z)))

        # 2) Failure streak: how many most-recent nights had no successful calibration.
        streak = 0
        for s in reversed(g["success"].tolist()):
            if s == 1:
                break
            streak += 1
        if streak >= config.FAILURE_STREAK_DAYS:
            alerts.append(dict(key=key, itype=itype, issue="Failure streak",
                               detail=f"{streak} nights since last success",
                               date=str(g.iloc[-1]["date"]), priority=float(streak)))

        # 3) Low recent success over the trailing window.
        recent = g.tail(config.RECENT_WINDOW_DAYS)
        if len(recent) >= 20:
            frac = recent["success"].mean()
            if frac < config.LOW_SUCCESS_FRAC:
                alerts.append(dict(key=key, itype=itype, issue="Low recent success",
                                   detail=f"{100*frac:.0f}% over last {len(recent)} nights",
                                   date=str(g.iloc[-1]["date"]),
                                   priority=float(100 * (config.LOW_SUCCESS_FRAC - frac))))

    if not alerts:
        return pd.DataFrame(columns=["key", "itype", "issue", "detail", "date", "priority"])
    return pd.DataFrame(alerts).sort_values(["issue", "priority"], ascending=[True, False]).reset_index(drop=True)
