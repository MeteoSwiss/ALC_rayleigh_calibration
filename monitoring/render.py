"""Render the static dashboard: summary page (index.html) + one page per station.

Reads the SQLite index, builds Plotly figures, fills Jinja2 templates, and writes a
self-contained site (shared plotly.min.js + css under assets/). All links are relative,
so the output folder works under any URL prefix on an internal web server.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from monitoring import charts, config, metrics

_TEMPLATES = Path(__file__).parent / "templates"
_STATIC = Path(__file__).parent / "static"


def _fmt(x, spec="{:.3g}", dash="—"):
    """Format a number, returning a dash for NaN/None (keeps templates clean)."""
    try:
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return dash
        return spec.format(x)
    except (TypeError, ValueError):
        return dash


def _env() -> Environment:
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES)),
        autoescape=select_autoescape(["html"]),
    )
    env.filters["fmt"] = _fmt
    env.globals["flag_label"] = config.flag_label
    env.globals["flag_color"] = config.flag_color
    return env


def _write_assets(out_dir: Path) -> str | None:
    """Vendor the shared Plotly bundle + copy static css/js. Return the logo filename if any.

    Drop an official logo at monitoring/static/eumetnet_logo.{png,svg,jpg} to replace the
    built-in placeholder wordmark; it is copied into assets/ and used in the header.
    """
    assets = out_dir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    (assets / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")
    for name in ("style.css", "table-sort.js", "paginate.js"):
        src = _STATIC / name
        if src.exists():
            shutil.copyfile(src, assets / name)
    logo = None
    for ext in ("svg", "png", "jpg", "jpeg"):
        src = _STATIC / f"eumetnet_logo.{ext}"
        if src.exists():
            shutil.copyfile(src, assets / src.name)
            logo = src.name
            break
    return logo


def _station_table_rows(cal: pd.DataFrame, st: pd.DataFrame) -> list[dict]:
    """Build the sortable summary-table rows (one per station) with a C_L sparkline."""
    rows = []
    last_cl_by_key = {
        key: g[g["success"] == 1].sort_values("datetime")["lidar_constant"].tail(30).tolist()
        for key, g in cal.groupby("key", sort=False)
    }
    for _, s in st.iterrows():
        spark = charts.sparkline_svg(last_cl_by_key.get(s["key"], []))
        rows.append(dict(
            key=s["key"], itype=s["itype"],
            success_rate=s.get("success_rate"), n_dates=s.get("n_dates"),
            n_success=s.get("n_success"), median_cl=s.get("median_cl"),
            median_rel_unc=s.get("median_rel_unc"),
            last_date=s.get("last_date"), last_flag=s.get("last_flag"),
            spark=spark,
        ))
    return rows


def build_site(db_path: Path, out_dir: Path, limit_pages: int | None = None) -> dict:
    """Generate the full static site from the index. Returns a small stats dict."""
    out_dir = Path(out_dir)
    (out_dir / "stations").mkdir(parents=True, exist_ok=True)
    logo = _write_assets(out_dir)
    env = _env()

    cal, st = metrics.load_frames(db_path)
    summary = metrics.network_summary(cal, st)
    flags = metrics.flag_distribution(cal)
    watch = metrics.watchlist(cal, st)

    # --- Summary page -------------------------------------------------------
    summary_figs = {
        "map": charts.fig_to_div(charts.network_map(st), "fig-map"),
        "success_type": charts.fig_to_div(charts.success_by_type(summary["by_type"]), "fig-stype"),
        "flag_dist": charts.fig_to_div(charts.flag_distribution_bar(flags), "fig-flags"),
        "cl_type": charts.fig_to_div(charts.cl_by_type_box(st), "fig-cltype"),
    }
    summary_html = env.get_template("summary.html").render(
        base="", logo=logo, summary=summary,
        by_type=summary["by_type"].to_dict("records"),
        figs=summary_figs,
        watch=watch.to_dict("records"),
        rows=_station_table_rows(cal, st),
    )
    (out_dir / "index.html").write_text(summary_html, encoding="utf-8")

    # --- Per-station pages --------------------------------------------------
    keys = list(st["key"])
    if limit_pages:
        keys = keys[:limit_pages]
    station_tmpl = env.get_template("station.html")
    for key in keys:
        g = cal[cal["key"] == key].sort_values("datetime")
        meta = st[st["key"] == key].iloc[0].to_dict()
        figs = {
            "cl": charts.fig_to_div(charts.cl_timeseries(g), "fig-cl"),
            "flags": charts.fig_to_div(charts.monthly_flag_bars(g), "fig-mflags"),
            "window": charts.fig_to_div(charts.window_timeseries(g), "fig-win"),
        }
        recent = g.tail(30).iloc[::-1].to_dict("records")
        html = station_tmpl.render(base="../", logo=logo, key=key, meta=meta,
                                   figs=figs, recent=recent)
        (out_dir / "stations" / f"{key}.html").write_text(html, encoding="utf-8")

    return dict(out_dir=str(out_dir), n_pages=len(keys), as_of=summary["as_of"])
