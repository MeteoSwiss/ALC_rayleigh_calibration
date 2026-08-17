"""Derived metrics read back from the SQLite index: network summary + watchlist.

Everything is per *series* = (station, method). Kept separate from index-building so
thresholds can be re-tuned and the dashboard re-rendered without re-scanning the CSVs.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from monitoring import config


def load_frames(db_path: Path):
    """Return (calibrations, series, stations, kalman, diagnostics) with parsed datetimes."""
    with sqlite3.connect(Path(db_path)) as con:
        cal = pd.read_sql("SELECT * FROM calibrations", con)
        series = pd.read_sql("SELECT * FROM series", con)
        st = pd.read_sql("SELECT * FROM stations", con)
        try:
            kal = pd.read_sql("SELECT * FROM kalman", con)
        except Exception:
            kal = pd.DataFrame(columns=["key", "method", "date", "kalman", "kalman_std"])
        try:
            diag = pd.read_sql("SELECT * FROM diagnostics", con)
        except Exception:
            diag = pd.DataFrame(columns=["key", "method", "date", "src"])
    cal["datetime"] = pd.to_datetime(cal["date"], format="%Y%m%d", errors="coerce")
    if len(kal):
        kal["datetime"] = pd.to_datetime(kal["date"], format="%Y%m%d", errors="coerce")
    return cal, series, st, kal, diag


def network_summary(cal: pd.DataFrame, series: pd.DataFrame, st: pd.DataFrame) -> dict:
    """Headline KPIs + per-(type, method) breakdown for the summary page."""
    n_ok = int(cal["success"].sum())
    n_total = int(len(cal))
    # Success rate = valid / ALL days: no-data (0), no-cloud/not-clear-night (-1) and the quality
    # rejections all count as non-success (true daily yield), consistent with the per-series rate.
    n_suit = int((~cal["flag"].isin([0, -1])).sum())   # attemptable days, kept for reference
    as_of = str(cal["date"].max()) if n_total else None

    by_tm = (
        series.groupby(["itype", "method"])
        .agg(n_series=("key", "count"), n_dates=("n_dates", "sum"),
             n_success=("n_success", "sum"), n_suitable=("n_suitable", "sum"))
        .reset_index()
    )
    by_tm["success_rate"] = 100.0 * by_tm["n_success"] / by_tm["n_dates"].replace(0, np.nan)
    by_tm = by_tm.sort_values(["itype", "method"])

    return dict(
        as_of=as_of,
        n_stations=int(len(st)),
        n_series=int(len(series)),
        n_calibrations=n_total,
        n_success=n_ok,
        success_rate=(100.0 * n_ok / n_total) if n_total else float("nan"),
        date_min=str(cal["date"].min()) if n_total else None,
        date_max=as_of,
        by_type_method=by_tm,
    )


def flag_distribution(cal: pd.DataFrame) -> pd.DataFrame:
    """Counts per (method, flag), labelled method-aware (so cloud 'No liquid cloud' / 'No data'
    are distinct from Rayleigh 'Not a clear night'), ordered best-to-worst."""
    counts = (cal.groupby(["method", "flag"], dropna=False).size()
              .reset_index(name="count"))
    counts["label"] = [f"{config.method_label(m)} · {config.flag_label(f, m)}"
                       for m, f in zip(counts["method"], counts["flag"])]
    counts["color"] = counts["flag"].map(config.flag_color)
    counts["order"] = counts["flag"].map(lambda f: -float(f) if pd.notna(f) else 999)
    return counts.sort_values(["order", "method"]).drop(columns="order").reset_index(drop=True)


def _robust_sigma(x: np.ndarray) -> float:
    """MAD-based robust standard-deviation estimate (1.4826 * MAD)."""
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def watchlist(cal: pd.DataFrame, st: pd.DataFrame) -> pd.DataFrame:
    """Series worth a look: value drift (both methods); failure streak / low recent success
    (Rayleigh only — for cloud a clear spell legitimately yields no calibration, so those two
    would flag every cloud series and are uninformative)."""
    itype_by_key = dict(zip(st["key"], st["itype"])) if len(st) else {}
    alerts = []
    for (key, method), g in cal.groupby(["key", "method"], sort=False):
        g = g.sort_values("datetime")
        itype = itype_by_key.get(key, "Unknown")
        ok = g[g["success"] == 1]
        mlabel = config.method_label(method)

        # 1) Value drift: latest successful value far from the series' robust baseline.
        if len(ok) >= 8:
            cl = ok["cal_value"].to_numpy(dtype=float)
            med = float(np.median(cl))
            sig = _robust_sigma(cl)
            if np.isfinite(sig) and sig > 0:
                z = abs(cl[-1] - med) / sig
                if z >= config.DRIFT_SIGMA:
                    alerts.append(dict(key=key, itype=itype, method=method, issue="Value drift",
                                       detail=f"{mlabel}: last {cl[-1]:.3g} vs median {med:.3g} ({z:.1f}σ)",
                                       date=str(ok.iloc[-1]["date"]), priority=float(z)))

        if method != "rayleigh":
            continue  # streak / low-success are not meaningful for the cloud method

        # 2) Failure streak (Rayleigh): consecutive most-recent nights with no success.
        streak = 0
        for s in reversed(g["success"].tolist()):
            if s == 1:
                break
            streak += 1
        if streak >= config.FAILURE_STREAK_DAYS:
            alerts.append(dict(key=key, itype=itype, method=method, issue="Failure streak",
                               detail=f"{mlabel}: {streak} nights since last success",
                               date=str(g.iloc[-1]["date"]), priority=float(streak)))

        # 3) Low recent success (Rayleigh) over the trailing window.
        recent = g.tail(config.RECENT_WINDOW_DAYS)
        if len(recent) >= 20:
            frac = recent["success"].mean()
            if frac < config.LOW_SUCCESS_FRAC:
                alerts.append(dict(key=key, itype=itype, method=method, issue="Low recent success",
                                   detail=f"{mlabel}: {100*frac:.0f}% over last {len(recent)} nights",
                                   date=str(g.iloc[-1]["date"]),
                                   priority=float(100 * (config.LOW_SUCCESS_FRAC - frac))))

    if not alerts:
        return pd.DataFrame(columns=["key", "itype", "method", "country", "issue", "detail", "date", "priority"])
    country_by_key = dict(zip(st["key"], st["country"])) if "country" in st.columns else {}
    df = pd.DataFrame(alerts).sort_values(["issue", "priority"], ascending=[True, False]).reset_index(drop=True)
    df["country"] = df["key"].map(lambda k: str(country_by_key.get(k, "") or ""))
    return df


#: Tile grading thresholds (absolute value of the metric, in its own unit). Display heuristics for
#: an operator's eye, NOT acceptance criteria -- the flags decide whether a night calibrates.
TILE_OK, TILE_WARN, TILE_BAD = "#1a7431", "#b7791f", "#b00020"


def _grade(v, good, warn):
    if v is None or not np.isfinite(v):
        return TILE_BAD
    a = abs(v)
    return TILE_OK if a < good else (TILE_WARN if a < warn else TILE_BAD)


def cl_headline_tiles(by_method: dict) -> list:
    """The five numbers that say whether a station's calibration is healthy right now.

    Median C_L, its p10-p90 spread, the cloud/Rayleigh ratio, the 15-day drift and how long since
    the last valid night. Colour-coded so a bad one is visible without reading the figure above it;
    the ratio in particular is the CL61 water-vapour diagnostic -- a systematic offset between two
    retrievals of the SAME constant points at the correction, not at the instrument.
    """
    ok = {m: g[g["success"] == 1].sort_values("datetime") for m, g in by_method.items()}
    ok = {m: g for m, g in ok.items() if len(g)}
    if not ok:
        return []
    ref = "cloud" if "cloud" in ok else max(ok, key=lambda m: len(ok[m]))
    g = ok[ref]
    med = float(g["cal_value"].median())
    p10, p90 = (float(g["cal_value"].quantile(q)) for q in (0.10, 0.90))
    spread = 100.0 * (p90 - p10) / med if med else None
    last = g["datetime"].max()
    win = last - pd.Timedelta(days=15)
    recent = g[g["datetime"] >= win]["cal_value"]
    prior = g[(g["datetime"] < win) & (g["datetime"] >= win - pd.Timedelta(days=45))]["cal_value"]
    drift = (100.0 * (recent.median() / prior.median() - 1.0)
             if len(recent) and len(prior) and prior.median() else None)
    ratio = None
    if "cloud" in ok and "rayleigh" in ok:
        mr = float(ok["rayleigh"]["cal_value"].median())
        ratio = float(ok["cloud"]["cal_value"].median()) / mr if mr else None
    age = int((pd.Timestamp.utcnow().tz_localize(None).normalize() - last.normalize()).days)
    return [
        {"label": f"MEDIAN {config.method_label(ref).upper()} C_L", "value": f"{med:.4g}",
         "note": f"{len(g)} calibrated nights", "color": TILE_OK},
        {"label": "SPREAD (P10-P90)", "value": "—" if spread is None else f"{spread:.0f}%",
         "note": f"{p10:.4g} … {p90:.4g}", "color": _grade(spread, 10, 20)},
        {"label": "CLOUD / RAYLEIGH RATIO", "value": "—" if ratio is None else f"{ratio:.3g}x",
         "note": ("only one method on this stream" if ratio is None else
                  ("the two retrievals agree" if abs(ratio - 1) < 0.05 else
                   "a systematic offset points at the water-vapour correction, not the instrument")),
         "color": TILE_BAD if ratio is None else _grade(ratio - 1.0, 0.05, 0.15)},
        {"label": "15-DAY DRIFT", "value": "—" if drift is None else f"{drift:+.1f}%",
         "note": "vs the preceding 45 d" if drift is not None else "not enough history",
         "color": _grade(drift, 5, 10)},
        {"label": "LAST VALID", "value": "today" if age <= 0 else f"{age} d ago",
         "note": last.strftime("%Y-%m-%d"), "color": _grade(age, 3, 8)},
    ]
