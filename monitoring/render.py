"""Render the static dashboard: summary page (index.html) + one page per station.

Reads the SQLite index, builds Plotly figures, fills Jinja2 templates, and writes a
self-contained site (shared plotly.min.js + css under assets/). All links are relative.
Each station may carry one or two method series (Rayleigh and/or cloud).
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
    try:
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return dash
        return spec.format(x)
    except (TypeError, ValueError):
        return dash


def _env() -> Environment:
    env = Environment(loader=FileSystemLoader(str(_TEMPLATES)),
                      autoescape=select_autoescape(["html"]))
    env.filters["fmt"] = _fmt
    env.globals["flag_label"] = config.flag_label
    env.globals["flag_color"] = config.flag_color
    env.globals["method_label"] = config.method_label
    return env


def _write_assets(out_dir: Path) -> str | None:
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


def _keystats(series: pd.DataFrame, st: pd.DataFrame) -> pd.DataFrame:
    """Per-key rollup (across methods) for the network map."""
    agg = (series.groupby("key")
           .agg(n_dates=("n_dates", "sum"), n_success=("n_success", "sum"),
                methods=("method", lambda s: " + ".join(config.method_label(m) for m in sorted(set(s)))))
           .reset_index())
    agg["success_rate"] = 100.0 * agg["n_success"] / agg["n_dates"].replace(0, np.nan)
    return agg.merge(st[["key", "itype", "lat", "lon"]], on="key", how="left")


def _series_table_rows(cal: pd.DataFrame, series: pd.DataFrame) -> list[dict]:
    """One row per (station, method) with a value sparkline + method badge."""
    last_by = {(k, m): g[g["success"] == 1].sort_values("datetime")["cal_value"].tail(30).tolist()
               for (k, m), g in cal.groupby(["key", "method"], sort=False)}
    rows = []
    for _, s in series.iterrows():
        method = s["method"]
        spark = charts.sparkline_svg(last_by.get((s["key"], method), []),
                                     color=config.METHOD_COLORS.get(method, "#1f77b4"))
        rows.append(dict(
            key=s["key"], itype=s.get("itype"), method=method, method_label=config.method_label(method),
            success_rate=s.get("success_rate"), n_dates=s.get("n_dates"), n_success=s.get("n_success"),
            median_cl=s.get("median_cl"), median_rel_unc=s.get("median_rel_unc"),
            last_date=s.get("last_date"), last_flag=s.get("last_flag"), spark=spark,
        ))
    return rows


def _method_block(key, method, cal, kal, series):
    """Figures + aggregates for one method section on a station page."""
    g_m = cal[(cal["key"] == key) & (cal["method"] == method)].sort_values("datetime")
    kal_m = kal[(kal["key"] == key) & (kal["method"] == method)] if len(kal) else kal
    srow = series[(series["key"] == key) & (series["method"] == method)]
    meta = srow.iloc[0].to_dict() if len(srow) else {}
    safe = method  # 'rayleigh'/'cloud' are id-safe
    return dict(
        method=method, label=config.method_label(method), meta=meta,
        figs={
            "ts": charts.fig_to_div(charts.series_timeseries(g_m, kal_m, method), f"fig-ts-{safe}"),
            "flags": charts.fig_to_div(charts.monthly_flag_bars(g_m, method), f"fig-mf-{safe}"),
            "aux": charts.fig_to_div(charts.aux_timeseries(g_m, method), f"fig-aux-{safe}"),
        },
        recent=g_m.tail(15).iloc[::-1].to_dict("records"),
    )


def build_site(db_path: Path, out_dir: Path, limit_pages: int | None = None) -> dict:
    out_dir = Path(out_dir)
    (out_dir / "stations").mkdir(parents=True, exist_ok=True)
    logo = _write_assets(out_dir)
    env = _env()

    cal, series, st, kal = metrics.load_frames(db_path)
    summary = metrics.network_summary(cal, series, st)
    flags = metrics.flag_distribution(cal)
    watch = metrics.watchlist(cal, st)
    keystats = _keystats(series, st)

    summary_figs = {
        "map": charts.fig_to_div(charts.network_map(keystats), "fig-map"),
        "success_type": charts.fig_to_div(charts.success_by_type_method(summary["by_type_method"]), "fig-stype"),
        "flag_dist": charts.fig_to_div(charts.flag_distribution_bar(flags), "fig-flags"),
        "cl_type": charts.fig_to_div(charts.value_by_type_method_box(series), "fig-cltype"),
    }
    summary_html = env.get_template("summary.html").render(
        base="", logo=logo, summary=summary, figs=summary_figs,
        watch=watch.to_dict("records"), rows=_series_table_rows(cal, series),
    )
    (out_dir / "index.html").write_text(summary_html, encoding="utf-8")

    # --- Per-station pages (one per key; all of that key's methods) ----------
    keys = list(st["key"])
    if limit_pages:
        keys = keys[:limit_pages]
    station_tmpl = env.get_template("station.html")
    for key in keys:
        meta = st[st["key"] == key].iloc[0].to_dict()
        methods = [m for m in config.METHOD_ORDER
                   if len(cal[(cal["key"] == key) & (cal["method"] == m)])]
        blocks = [_method_block(key, m, cal, kal, series) for m in methods]
        overlay = None
        if len(methods) >= 2:
            by_method = {m: cal[(cal["key"] == key) & (cal["method"] == m)] for m in methods}
            overlay = charts.fig_to_div(charts.normalized_overlay(by_method), "fig-overlay")
        html = station_tmpl.render(base="../", logo=logo, key=key, meta=meta,
                                   blocks=blocks, overlay=overlay)
        (out_dir / "stations" / f"{key}.html").write_text(html, encoding="utf-8")

    return dict(out_dir=str(out_dir), n_pages=len(keys), n_series=int(len(series)),
                as_of=summary["as_of"])
