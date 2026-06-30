"""Time-period presets for the dashboard — the single source of truth.

A *period* is a labelled date window over the calibration archive. The same set drives:
  * the station-page client-side selector (``rangesync.js``, via an embedded JSON),
  * the per-period summary pages (``build_dashboard.py`` / ``render.build_site``), and
  * the per-period OmB / sensitivity images (``scripts/run_all_l1_2026.py``).

Years are *derived from the data span*, not hardcoded: every calendar year from the first
year with data to the current year (2025, 2026, then 2027, 2028 … automatically). Rolling
windows (last N days), the current calendar year and all-time are **active** (recomputed every
daily build); a past calendar year is **frozen** once it has ended *plus* a backfill guard (so
the self-healing backfill can no longer write into it) — built once, then skipped on later
builds. This keeps the daily work bounded at ~6 active windows no matter how many years pile up.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta

#: Rolling trailing windows offered alongside the years + all-time (days).
DEFAULT_ROLLING_DAYS = (30, 90, 180, 365)


@dataclass(frozen=True)
class Period:
    key: str            # id-safe slug: 'all', 'y2025', 'last90' (used in index_<key>.html / *_omb_<key>.png)
    label: str          # human label for the <select>: 'All time', '2025', 'Last 90 d'
    start: str | None   # inclusive 'YYYYMMDD', or None = open start (all-time)
    end: str | None     # inclusive 'YYYYMMDD', or None = open end
    kind: str           # 'all' | 'year' | 'rolling'
    active: bool        # True = rebuild every daily build; False = frozen (build once)

    @property
    def iso_start(self) -> str | None:
        return f"{self.start[:4]}-{self.start[4:6]}-{self.start[6:8]}" if self.start else None

    @property
    def iso_end(self) -> str | None:
        return f"{self.end[:4]}-{self.end[4:6]}-{self.end[6:8]}" if self.end else None


def _ymd(d: date) -> str:
    return d.strftime("%Y%m%d")


def _to_date(ymd: str) -> date:
    return datetime.strptime(str(ymd)[:8], "%Y%m%d").date()


def build_periods(date_min: str | None, date_max: str | None, *, today: date | None = None,
                  backfill_days: int = 5, rolling_days=DEFAULT_ROLLING_DAYS) -> list[Period]:
    """Period list for a data span ['date_min' .. 'date_max'] (YYYYMMDD strings).

    ``today`` defaults to the real current date (used only to decide which years are frozen);
    ``backfill_days`` mirrors ``ALC_BACKFILL_DAYS`` so a just-ended year stays *active* until the
    self-healing backfill window has fully passed.
    """
    if not date_min or not date_max:
        return [Period("all", "All time", None, None, "all", True)]
    today = today or date.today()
    dmin, dmax = _to_date(date_min), _to_date(date_max)
    out: list[Period] = [Period("all", "All time", None, None, "all", True)]
    for y in range(dmin.year, dmax.year + 1):
        ystart, yend = date(y, 1, 1), date(y, 12, 31)
        frozen = today > yend + timedelta(days=backfill_days)   # year ended + backfill can't touch it
        out.append(Period(f"y{y}", str(y), _ymd(ystart), _ymd(yend), "year", not frozen))
    for n in rolling_days:
        start = dmax - timedelta(days=n - 1)
        out.append(Period(f"last{n}", f"Last {n} d", _ymd(start), _ymd(dmax), "rolling", True))
    return out


def periods_to_json(periods: list[Period]) -> list[dict]:
    """JSON-friendly list for embedding in pages (consumed by ``rangesync.js``).

    ``start``/``end`` are ISO 'YYYY-MM-DD' (or null) so the browser's ``Date`` / Plotly's date
    axis parse them directly; ``ymd_start``/``ymd_end`` keep the compact form for row-date compares.
    """
    return [dict(key=p.key, label=p.label, start=p.iso_start, end=p.iso_end,
                 ymd_start=p.start, ymd_end=p.end, kind=p.kind, active=p.active)
            for p in periods]


def filter_window(df, start: str | None, end: str | None, col: str = "date"):
    """Subset a DataFrame whose ``col`` holds 'YYYYMMDD' strings to [start, end] inclusive.

    ``None`` bounds are open. Returns the (possibly unfiltered) frame; never mutates in place.
    """
    if df is None or not len(df):
        return df
    if start is not None:
        df = df[df[col] >= start]
    if end is not None:
        df = df[df[col] <= end]
    return df
