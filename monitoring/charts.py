"""Plotly figure builders (+ a cheap inline-SVG sparkline for tables).

Everything is per *series* = (station, method). Figures are emitted as <div> via fig_to_div();
the shared plotly.min.js is loaded once per page (see render.py), so include_plotlyjs=False.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from monitoring import config, kalman

_LAYOUT = dict(template="plotly_white", margin=dict(l=64, r=20, t=44, b=40),
               font=dict(size=12), height=340)

# Both methods report the lidar constant C_L (Wiegner); the cloud value is the O'Connor
# C_L = applied_constant / C, on the same scale as Rayleigh -- the operationally useful number.
_VALUE_NAME = {"rayleigh": "C_L", "cloud": "C_L"}


def fig_to_div(fig: go.Figure, div_id: str) -> str:
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id=div_id,
                       config={"displaylogo": False, "responsive": True})


def _value_name(method: str) -> str:
    return _VALUE_NAME.get(str(method), "value")


# --- Map helpers ------------------------------------------------------------

def _symbols_for(itypes) -> list:
    """Per-point marker symbol array keyed by instrument type. scattergeo honours a per-POINT
    marker.symbol array, so every map stays single-trace (filter.js still drives marker.size via
    customdata, and the symbols simply ride along)."""
    return [config.type_symbol(t) for t in itypes]


def _add_symbol_legend(fig: go.Figure) -> None:
    """Add a small manual symbol->instrument-type legend via dummy legend-only scattergeo traces
    (no real points). The data trace carries showlegend=False, so the legend reads purely as the
    symbol key. Placed top-left, horizontal, so it does not crowd the colorbar on the right."""
    for t in config.TYPE_ORDER:
        fig.add_trace(go.Scattergeo(
            lat=[None], lon=[None], mode="markers", name=t,
            marker=dict(size=9, symbol=config.type_symbol(t), color="#666",
                        line=dict(width=0.4, color="#555")),
            showlegend=True, hoverinfo="skip"))
    fig.update_layout(legend=dict(orientation="h", yanchor="top", y=1.0, xanchor="left", x=0.0,
                                  bgcolor="rgba(255,255,255,0.6)", font=dict(size=10),
                                  itemclick=False, itemdoubleclick=False))


# --- Summary-page figures ---------------------------------------------------

def network_map(keystats: pd.DataFrame) -> go.Figure:
    """Station map over Europe, colored by per-station success rate (max over its methods)."""
    d = keystats.dropna(subset=["lat", "lon"]).copy()
    if d.empty:
        # Empty result (e.g. a date window or type filter with no calibrations): render an
        # empty map rather than crashing Plotly's marker validator on an object-dtype Series.
        fig = go.Figure()
        fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                        showland=True, landcolor="#f5f5f5",
                        lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
        fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)},
                          title="Network — no calibrations in range")
        return fig
    n_dates = np.asarray(d["n_dates"].fillna(1), dtype=float)
    sizes = np.clip(np.sqrt(n_dates) / 2.0, 5, 18)
    country = (d["country"].fillna("") if "country" in d.columns
               else pd.Series([""] * len(d), index=d.index)).astype(str)
    name = (d["name"].fillna("") if "name" in d.columns
            else pd.Series([""] * len(d), index=d.index)).astype(str)
    # customdata = [country, type, base_size, key, name] per point: filter.js restyles
    # marker.size from [2] (0 = hidden) AND navigates to stations/<key>.html ([3]) on click.
    # A typed list (not np.column_stack, which would stringify the float size).
    customdata = [[c, t, float(sz), k, nm] for c, t, sz, k, nm in
                  zip(country.values, d["itype"].astype(str).values, sizes,
                      d["key"].astype(str).values, name.values)]
    fig = go.Figure(go.Scattergeo(
        lat=d["lat"], lon=d["lon"],
        # Hover: station name (bold) + WIGOS id (key) + details. (click opens the station page)
        text=[f"<b>{nm or k}</b><br>{k}<br>{t} · {ct} · {m}<br>{sr:.0f}% success · {nd} cal"
              for k, nm, t, ct, m, sr, nd in zip(d["key"], name, d["itype"], country, d["methods"],
                                                 d["success_rate"].fillna(0), d["n_dates"].fillna(0))],
        hoverinfo="text", customdata=customdata, showlegend=False,
        marker=dict(size=sizes, symbol=_symbols_for(d["itype"]),
                    color=np.asarray(d["success_rate"].fillna(0), dtype=float),
                    colorscale="RdYlGn", cmin=0, cmax=80,
                    colorbar=dict(title="success %"), line=dict(width=0.4, color="#555")),
    ))
    _add_symbol_legend(fig)
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)},
                      title="Network — success rate by station")
    return fig


def ratio_map(keystats: pd.DataFrame, col: str, title: str, cbar: str,
              div_id: str, cmin: float = 0.0, cmax: float = 200.0) -> go.Figure:
    """Station map colored by a percent ratio in column `col` (100 % = on target), diverging
    around 100 %. Carries the same customdata as network_map so filter.js can filter/navigate it."""
    base = keystats.dropna(subset=["lat", "lon"]).copy()
    d = base[np.isfinite(pd.to_numeric(base.get(col), errors="coerce"))].copy() if col in base else base.iloc[0:0]
    if d.empty:
        fig = go.Figure()
        fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                        showland=True, landcolor="#f5f5f5",
                        lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
        fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)},
                          title=f"{title} — no data")
        return fig
    n_dates = np.asarray(d["n_dates"].fillna(1), dtype=float)
    sizes = np.clip(np.sqrt(n_dates) / 2.0, 5, 18)
    country = (d["country"].fillna("") if "country" in d.columns
               else pd.Series([""] * len(d), index=d.index)).astype(str)
    name = (d["name"].fillna("") if "name" in d.columns
            else pd.Series([""] * len(d), index=d.index)).astype(str)
    customdata = [[c, t, float(sz), k, nm] for c, t, sz, k, nm in
                  zip(country.values, d["itype"].astype(str).values, sizes,
                      d["key"].astype(str).values, name.values)]
    vals = np.asarray(pd.to_numeric(d[col], errors="coerce"), dtype=float)
    fig = go.Figure(go.Scattergeo(
        lat=d["lat"], lon=d["lon"],
        text=[f"<b>{nm or k}</b><br>{k}<br>{t} · {ct}<br>{v:.0f}% of reference · {nd:.0f} cal"
              for k, nm, t, ct, v, nd in zip(d["key"], name, d["itype"], country, vals,
                                             d["n_dates"].fillna(0))],
        hoverinfo="text", customdata=customdata, showlegend=False,
        marker=dict(size=sizes, symbol=_symbols_for(d["itype"]),
                    color=vals, colorscale="RdBu", cmin=cmin, cmax=cmax, cmid=100.0,
                    colorbar=dict(title=cbar), line=dict(width=0.4, color="#555")),
    ))
    _add_symbol_legend(fig)
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)}, title=title)
    return fig


def _empty_map(title: str) -> go.Figure:
    """An empty Europe map with a '— no data' title (shared by the value-coloured map builders)."""
    fig = go.Figure()
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)},
                      title=f"{title} — no data")
    return fig


def omb_bias_map(keystats: pd.DataFrame) -> go.Figure:
    """Station map coloured by the mean observation-minus-background (CAMS) bias of OUR calibrated
    backscatter, ``omb_bias`` [Mm^-1 sr^-1] (= median_bias_ours * 1e6). RdBu diverging centred at 0;
    blue = we read low vs CAMS, red = high. Carries the same customdata layout + per-type symbol as
    the other maps so filter.js country/type filtering and click-to-station keep working."""
    title = "Mean OmB bias vs CAMS"
    col = "omb_bias"
    base = keystats.dropna(subset=["lat", "lon"]).copy()
    d = base[np.isfinite(pd.to_numeric(base.get(col), errors="coerce"))].copy() if col in base else base.iloc[0:0]
    if d.empty:
        return _empty_map(title)
    n_dates = np.asarray(d["n_dates"].fillna(1), dtype=float)
    sizes = np.clip(np.sqrt(n_dates) / 2.0, 5, 18)
    country = (d["country"].fillna("") if "country" in d.columns
               else pd.Series([""] * len(d), index=d.index)).astype(str)
    name = (d["name"].fillna("") if "name" in d.columns
            else pd.Series([""] * len(d), index=d.index)).astype(str)
    customdata = [[c, t, float(sz), k, nm] for c, t, sz, k, nm in
                  zip(country.values, d["itype"].astype(str).values, sizes,
                      d["key"].astype(str).values, name.values)]
    vals = np.asarray(pd.to_numeric(d[col], errors="coerce"), dtype=float)
    fig = go.Figure(go.Scattergeo(
        lat=d["lat"], lon=d["lon"],
        text=[f"<b>{nm or k}</b><br>{k}<br>{t} · {ct}<br>O-B = {v:+.2f} Mm⁻¹sr⁻¹"
              for k, nm, t, ct, v in zip(d["key"], name, d["itype"], country, vals)],
        hoverinfo="text", customdata=customdata, showlegend=False,
        marker=dict(size=sizes, symbol=_symbols_for(d["itype"]),
                    color=vals, colorscale="RdBu", cmin=-1.0, cmax=1.0, cmid=0.0,
                    colorbar=dict(title="O-B [Mm⁻¹sr⁻¹]"), line=dict(width=0.4, color="#555")),
    ))
    _add_symbol_legend(fig)
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)}, title=title)
    return fig


def icao_altitude_map(keystats: pd.DataFrame) -> go.Figure:
    """Station map coloured by the ICAO 200 ug/m3 detection altitude ``icao_alt`` [m] (= icao_alt_200):
    the highest altitude at which the instrument can still detect the ICAO low-visibility aerosol load.
    Sequential Viridis, higher = better. Stations with no detection (blank icao_alt) are omitted (they
    fall out via the finite-value filter). Same customdata + per-type symbol as the other maps."""
    title = "ICAO detection altitude"
    col = "icao_alt"
    base = keystats.dropna(subset=["lat", "lon"]).copy()
    d = base[np.isfinite(pd.to_numeric(base.get(col), errors="coerce"))].copy() if col in base else base.iloc[0:0]
    if d.empty:
        return _empty_map(title)
    n_dates = np.asarray(d["n_dates"].fillna(1), dtype=float)
    sizes = np.clip(np.sqrt(n_dates) / 2.0, 5, 18)
    country = (d["country"].fillna("") if "country" in d.columns
               else pd.Series([""] * len(d), index=d.index)).astype(str)
    name = (d["name"].fillna("") if "name" in d.columns
            else pd.Series([""] * len(d), index=d.index)).astype(str)
    customdata = [[c, t, float(sz), k, nm] for c, t, sz, k, nm in
                  zip(country.values, d["itype"].astype(str).values, sizes,
                      d["key"].astype(str).values, name.values)]
    vals = np.asarray(pd.to_numeric(d[col], errors="coerce"), dtype=float)
    fig = go.Figure(go.Scattergeo(
        lat=d["lat"], lon=d["lon"],
        text=[f"<b>{nm or k}</b><br>{k}<br>{t} · {ct}<br>ICAO 200 µg/m³ detect alt = {v:.0f} m"
              for k, nm, t, ct, v in zip(d["key"], name, d["itype"], country, vals)],
        hoverinfo="text", customdata=customdata, showlegend=False,
        marker=dict(size=sizes, symbol=_symbols_for(d["itype"]),
                    color=vals, colorscale="Viridis",
                    colorbar=dict(title=dict(text="alt [m]", side="right"), thickness=12),
                    line=dict(width=0.4, color="#555")),
    ))
    _add_symbol_legend(fig)
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)}, title=title)
    return fig


def instrument_count_over_time(activity: list) -> go.Figure:
    """Stacked monthly count of ACTIVE instruments over time, coloured by instrument type. An
    instrument (station key) is 'active' in a month if it has >=1 calibration day that month.

    ``activity`` is the list embedded for the client filter: [{key, itype, country, months:[YYYYMM,...]}].
    filter.js recomputes this same figure client-side from that list when the country/type filter
    changes (country filtering changes per-month counts, so a simple restyle won't do)."""
    from collections import defaultdict
    counts: dict = defaultdict(lambda: defaultdict(int))   # month -> itype -> distinct-key count
    months_set: set = set()
    types_present: set = set()
    for r in (activity or []):
        t = str(r.get("itype", "") or "")
        types_present.add(t)
        for m in set(r.get("months", []) or []):
            counts[m][t] += 1
            months_set.add(m)
    if not months_set:
        fig = go.Figure()
        fig.update_layout(**{**_LAYOUT, "height": 320}, title="Active instruments over time — no data")
        return fig
    months = sorted(months_set)
    x = [pd.to_datetime(m + "01", format="%Y%m%d", errors="coerce") for m in months]
    order = [t for t in config.TYPE_ORDER if t in types_present] + \
            sorted(types_present - set(config.TYPE_ORDER) - {""})
    fig = go.Figure()
    for t in order:
        y = [counts[m].get(t, 0) for m in months]
        if not any(y):
            continue
        fig.add_trace(go.Bar(x=x, y=y, name=t, marker_color=config.TYPE_COLORS.get(t, "#888"),
                             hovertemplate="%{x|%Y-%m}<br>" + t + ": %{y} instruments<extra></extra>"))
    fig.update_layout(**{**_LAYOUT, "height": 320}, barmode="stack",
                      title="Active instruments over time (by type)", yaxis_title="# instruments",
                      legend=dict(orientation="h", y=-0.2, font=dict(size=11)))
    return fig


def success_by_type_method(by_tm: pd.DataFrame) -> go.Figure:
    """Grouped bar: success rate per instrument type, one bar per method."""
    fig = go.Figure()
    for method in config.METHOD_ORDER:
        d = by_tm[by_tm["method"] == method]
        if not len(d):
            continue
        fig.add_trace(go.Bar(
            x=d["itype"], y=d["success_rate"], name=config.method_label(method),
            marker_color=config.METHOD_COLORS.get(method, "#888"),
            text=[f"{v:.0f}%" for v in d["success_rate"]], textposition="outside",
            hovertext=[f"{n} series · {s}/{t} cal" for n, s, t in
                       zip(d["n_series"], d["n_success"], d["n_dates"])],
        ))
    fig.update_layout(**{**_LAYOUT, "margin": dict(l=64, r=20, t=44, b=72)},
                      barmode="group", title="Success rate by type & method",
                      yaxis_title="success %", yaxis_range=[0, 110],
                      legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center", yanchor="top"))
    return fig


def flag_distribution_bar(flags: pd.DataFrame, method: str | None = None) -> go.Figure:
    """Horizontal bar of calibration counts per outcome. With `method`, restrict to that method and
    drop the method prefix from the labels (one chart per method)."""
    if method is not None:
        d = flags[flags["method"] == method].copy()
        title = f"Outcome distribution — {config.method_label(method)}"
    else:
        d = flags.copy()
        title = "Outcome distribution by method"
    d = d.iloc[::-1]
    ylabels = ([config.flag_label(f, method) for f in d["flag"]] if method is not None
               else list(d["label"]))
    height = max(300, 24 * len(d) + 80)
    fig = go.Figure(go.Bar(x=d["count"], y=ylabels, orientation="h",
                           marker_color=d["color"], text=d["count"], textposition="auto"))
    fig.update_layout(**{**_LAYOUT, "height": height, "margin": dict(l=210, r=20, t=44, b=40)},
                      title=title, xaxis_title="calibrations")
    return fig


def value_by_type_method_box(series: pd.DataFrame) -> go.Figure:
    """Box of per-series median lidar constant C_L, grouped by type, colored by method (log).

    C_L follows the per-instrument scale (CHM15k ~3e11, CL31/CL51 ~1e8, CL61 ~1, Mini-MPL ~5e5),
    so a log axis is essential; for CL61 the Rayleigh and cloud boxes should overlap (same C_L).
    """
    fig = go.Figure()
    for method in config.METHOD_ORDER:
        for t in config.TYPE_ORDER:
            vals = series.loc[(series["itype"] == t) & (series["method"] == method)
                              & (series["median_cl"] > 0), "median_cl"]
            if len(vals):
                fig.add_trace(go.Box(y=vals, name=t, legendgroup=method,
                                     marker_color=config.METHOD_COLORS.get(method, "#888"),
                                     boxpoints="all", jitter=0.4, pointpos=0,
                                     offsetgroup=method, showlegend=False))
    fig.update_layout(**_LAYOUT, title="Per-series median C_L by type — absolute (log)",
                      yaxis_title="C_L", yaxis_type="log", boxmode="group")
    fig.update_yaxes(exponentformat="e")
    return fig


def value_pct_theoretical_box(series: pd.DataFrame) -> go.Figure:
    """Box of per-series median C_L expressed as a PERCENT of the theoretical value, by type
    (linear). Putting every type on a common 0-?% scale makes the across-type comparison readable;
    100 % (dashed) is the nominal value."""
    fig = go.Figure()
    for method in config.METHOD_ORDER:
        for t in config.TYPE_ORDER:
            theo = config.theoretical_cl(t)
            if not theo or theo <= 0:
                continue
            vals = series.loc[(series["itype"] == t) & (series["method"] == method)
                              & (series["median_cl"] > 0), "median_cl"]
            if len(vals):
                fig.add_trace(go.Box(y=100.0 * vals / theo, name=t, legendgroup=method,
                                     marker_color=config.METHOD_COLORS.get(method, "#888"),
                                     boxpoints="all", jitter=0.4, pointpos=0,
                                     offsetgroup=method, showlegend=False))
    fig.add_hline(y=100.0, line=dict(color="#888", width=1, dash="dash"))
    fig.update_layout(**_LAYOUT, title="Per-series median C_L by type — % of theoretical (linear)",
                      yaxis_title="% of theoretical", boxmode="group")
    return fig


def cl_median_iqr_by_station(d: pd.DataFrame, itype: str) -> go.Figure:
    """Stations of ONE instrument type ranked low->high by their MEDIAN lidar constant C_L, with the
    interquartile range (Q1..Q3 of that station's successful daily C_L) drawn as an asymmetric error
    bar. Red dashed = the theoretical C_L for the type; blue dashed = the type's network median of the
    per-station medians. ``d`` has one row per station: columns key, med, q1, q3, n[, country].
    """
    d = d[d["med"] > 0].sort_values("med").reset_index(drop=True)
    if d.empty:
        fig = go.Figure()
        fig.update_layout(**_LAYOUT, title=f"{itype} — no calibrations")
        return fig
    rank = np.arange(len(d))
    color = config.TYPE_COLORS.get(itype, "#1f77b4")
    fig = go.Figure(go.Bar(
        x=rank, y=d["med"], name=itype, marker_color=color, width=1.0,
        error_y=dict(type="data", symmetric=False,
                     array=(d["q3"] - d["med"]).clip(lower=0),
                     arrayminus=(d["med"] - d["q1"]).clip(lower=0),
                     thickness=0.7, width=0, color="rgba(50,50,50,0.55)"),
        customdata=list(zip(d["key"], d["n"], d["q1"], d["q3"],
                            d["country"] if "country" in d.columns else [""] * len(d),
                            d["med"])),   # [5]=median: the client country-filter rebuilds bars from customdata
        hovertemplate=("%{customdata[0]}<br>median C_L = %{y:.3g}"
                       "<br>IQR = [%{customdata[2]:.3g}, %{customdata[3]:.3g}]  (n=%{customdata[1]})"
                       "<br><i>click to open station →</i>"
                       "<extra>" + itype + "</extra>")))
    theo = config.theoretical_cl(itype)
    if theo and theo > 0:
        fig.add_hline(y=theo, line=dict(color="#d62728", width=1.5, dash="dash"),
                      annotation_text=f"theoretical ({theo:.3g})", annotation_position="top left")
    net = float(d["med"].median())
    fig.add_hline(y=net, line=dict(color="#1f77b4", width=1.5, dash="dash"),
                  annotation_text=f"network median ({net:.3g})", annotation_position="bottom right")
    fig.update_layout(**_LAYOUT, showlegend=False,
                      title=f"{itype} — stations ranked by median C_L (n={len(d)} stations, error bar = IQR)",
                      xaxis_title="station (ranked by median C_L)", yaxis_title="C_L")
    fig.update_layout(height=400)
    fig.update_yaxes(exponentformat="e", rangemode="tozero")
    return fig


# --- Station-page figures (one set per method) ------------------------------

#: Numeric algorithm code -> the label an operator recognises (calibration/io/output.py
#: VERSION_CODES). Used to name the series from the DATA rather than from a constant.
_VERSION_NAMES = {25: "v0.25", 90: "main", 95: "improved", 100: "v1.0", 110: "v1.1", 120: "v1.2",
                  200: "v2.0", 205: "v2.0p", 220: "v2.2", 300: "earlinet", 310: "bellini",
                  1000: "O'Connor"}


def version_label(g_m: pd.DataFrame) -> str:
    """The algorithm vintage actually present in these rows, or "" when the archive does not say.

    This label used to be the literal string "v2.0", written when the pipeline ran eprof_v2 -- so a
    v2.2 archive was displayed, in the legend and in the Kalman trace, as v2.0. Deriving it from the
    data means the page cannot misreport the release; an archive predating the `version` column gets
    no version claim at all, which is honest rather than convenient.
    """
    if "version" not in g_m.columns:
        return ""
    v = pd.to_numeric(g_m["version"], errors="coerce").dropna()
    if not len(v):
        return ""
    codes = sorted(set(int(x) for x in v))
    names = [_VERSION_NAMES.get(c, f"code {c}") for c in codes]
    return " + ".join(names)


def series_timeseries(g_m: pd.DataFrame, kal_m: pd.DataFrame, method: str,
                      op_df: pd.DataFrame | None = None,
                      oldray_df: pd.DataFrame | None = None) -> go.Figure:
    """Calibration value over time for ONE method: successes + uncertainty + Kalman best estimate,
    plus (optional) the daily OPERATIONAL calibration constant from the L2 files as a black line,
    and (optional, Rayleigh only) the OLD operational Rayleigh calibration as black 'x' markers.

    The Kalman line/band is the operational E-PROFILE random-walk best estimate, preferring
    the precomputed series (kal_m) and falling back to an on-the-fly fit.
    """
    vname = _value_name(method)
    color = config.METHOD_COLORS.get(method, "#1f77b4")
    ok = g_m[g_m["success"] == 1].sort_values("datetime")
    ver = version_label(ok if len(ok) else g_m)
    fig = go.Figure()
    if op_df is not None and len(op_df):
        od = op_df.sort_values("datetime")
        fig.add_trace(go.Scatter(
            x=od["datetime"], y=od["op_coeff"], mode="lines", name="Applied in L2",
            line=dict(color="#111111", width=1.3),
            hovertemplate="%{x|%Y-%m-%d}<br>operational=%{y:.3e}<extra></extra>"))
    if oldray_df is not None and len(oldray_df):
        orr = oldray_df.sort_values("datetime")
        fig.add_trace(go.Scatter(
            x=orr["datetime"], y=orr["value"], mode="markers", name="v1.0",
            marker=dict(symbol="x", size=7, color="#000000"), visible="legendonly",
            hovertemplate="%{x|%Y-%m-%d}<br>old Rayleigh=%{y:.3e}<extra></extra>"))
    if len(ok):
        fig.add_trace(go.Scatter(
            x=ok["datetime"], y=ok["cal_value"], mode="markers",
            name=(f"{vname} · {ver}" if ver else f"{vname} (per cal)"),
            marker=dict(size=5, color=color, opacity=0.7),
            error_y=dict(type="data", array=ok["uncertainty"], visible=True,
                         thickness=0.6, width=0, color="rgba(120,120,120,0.3)"),
            hovertemplate="%{x|%Y-%m-%d}<br>" + vname + "=%{y:.3e}<extra></extra>"))
        kt, ks, kstd = _kalman_xy(ok, kal_m)
        if len(kt):
            fig.add_trace(go.Scatter(
                x=np.concatenate([kt, kt[::-1]]),
                y=np.concatenate([ks + kstd, (ks - kstd)[::-1]]),
                fill="toself", fillcolor="rgba(214,39,40,0.12)", line=dict(width=0),
                hoverinfo="skip", showlegend=False))
            fig.add_trace(go.Scatter(x=kt, y=ks, mode="lines",
                                     name=(f"Kalman estimate ({ver})" if ver
                                           else "Kalman best estimate"),
                                     line=dict(color="#d62728", width=2),
                                     hovertemplate="%{x|%Y-%m-%d}<br>Kalman=%{y:.3e}<extra></extra>"))
    fig.update_layout(**_LAYOUT, yaxis_title=vname,
                      legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center", yanchor="top"),
                      title=f"{config.method_label(method)} — {vname} over time + Kalman")
    fig.update_yaxes(exponentformat="e")
    return fig


def _kalman_xy(ok: pd.DataFrame, kal_m: pd.DataFrame):
    """Precomputed Kalman (kal_m) if present, else an on-the-fly fit of the successes."""
    if kal_m is not None and len(kal_m):
        k = kal_m.sort_values("date").copy()
        k["dt"] = pd.to_datetime(k["date"], format="%Y%m%d", errors="coerce")
        k = k.dropna(subset=["dt", "kalman"])
        if len(k):
            return (k["dt"].to_numpy(), k["kalman"].to_numpy(dtype=float),
                    k["kalman_std"].fillna(0).to_numpy(dtype=float))
    kt, ks, kstd = kalman.kalman_best_estimate(ok["datetime"].tolist(),
                                               ok["cal_value"].to_numpy())
    return kt, ks, kstd


def monthly_flag_bars(g_m: pd.DataFrame, method: str) -> go.Figure:
    """Stacked monthly outcome counts for ONE method."""
    d = g_m.copy()
    d["month"] = d["datetime"].dt.to_period("M").dt.to_timestamp()
    fig = go.Figure()
    for f in sorted(d["flag"].dropna().unique(), key=lambda f: -float(f)):
        counts = d[d["flag"] == f].groupby("month").size()
        fig.add_trace(go.Bar(x=counts.index, y=counts.values, name=config.flag_label(f, method),
                             marker_color=config.flag_color(f)))
    fig.update_layout(**{**_LAYOUT, "height": 320}, barmode="stack",
                      title=f"{config.method_label(method)} — monthly outcomes",
                      yaxis_title="calibrations", legend=dict(orientation="h", y=-0.25, font=dict(size=10)))
    return fig


def aux_timeseries(g_m: pd.DataFrame, method: str) -> go.Figure:
    """Rayleigh -> calibration window (bottom/top); cloud -> number of in-cloud profiles."""
    fig = go.Figure()
    if method == "rayleigh":
        ok = g_m[(g_m["success"] == 1) & g_m["bottom_height"].notna()].sort_values("datetime")
        if len(ok):
            # One VERTICAL SEGMENT per night, bottom to top: the window is an extent, and two
            # disconnected dot clouds ("top" green, "bottom" brown) forced the reader to pair them
            # by eye across the whole plot. None separators keep it a single cheap trace.
            xs, ys = [], []
            for t, b, tp in zip(ok["datetime"], ok["bottom_height"], ok["top_height"]):
                xs += [t, t, None]
                ys += [b, tp, None]
            col = config.METHOD_COLORS.get("rayleigh", "#1f77b4")
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name="fitted window (bottom–top)",
                                     line=dict(color=col, width=2), hoverinfo="skip"))
            mid = (ok["bottom_height"] + ok["top_height"]) / 2.0
            fig.add_trace(go.Scatter(
                x=ok["datetime"], y=mid, mode="markers", showlegend=False,
                marker=dict(size=3, color=col),
                customdata=np.stack([ok["bottom_height"], ok["top_height"]], axis=1),
                hovertemplate="%{x|%Y-%m-%d}<br>%{customdata[0]:.0f}–%{customdata[1]:.0f} m"
                              "<extra></extra>"))
        fig.update_layout(**{**_LAYOUT, "height": 300}, yaxis_title="height (m AGL)",
                          title="Calibration window", legend=dict(orientation="h", y=1.14))
    else:
        ok = g_m[g_m["success"] == 1].sort_values("datetime")
        fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["n_profiles"], mode="markers",
                                 marker=dict(size=5, color="#2ca02c"), name="in-cloud profiles"))
        fig.update_layout(**{**_LAYOUT, "height": 300}, yaxis_title="# in-cloud profiles",
                          title="Cloud profiles per calibration", showlegend=False)
    return fig


def cl_overlay(by_method: dict) -> go.Figure:
    """Both retrievals of C_L on ONE axis, with the per-night uncertainty and the median band.

    The two methods estimate the SAME Wiegner constant, so the comparison is absolute -- no
    normalisation. Bare markers made a drift or a step hard to see, so this adds three things that
    carry the reading: the per-night uncertainty as error bars, the median with a +/-10 % band (a
    stable instrument sits inside it), and a shaded 15-day window at the right so the drift quoted
    in the headline tiles can be seen rather than taken on trust.
    """
    fig = go.Figure()
    ref, all_ok = None, []
    for method, g_m in by_method.items():
        ok = g_m[g_m["success"] == 1].sort_values("datetime")
        if not len(ok):
            continue
        all_ok.append((method, ok))
        # The reference for the median band is the method with the most calibrated nights: its
        # median is the better-determined one (cloud typically runs on far more nights).
        if ref is None or len(ok) > len(ref[1]):
            ref = (method, ok)
    for method, ok in all_ok:
        col = config.METHOD_COLORS.get(method, "#888")
        unc = ok["uncertainty"] if "uncertainty" in ok.columns else None
        fig.add_trace(go.Scatter(
            x=ok["datetime"], y=ok["cal_value"], mode="markers",
            name=config.method_label(method),
            marker=dict(size=5, color=col, opacity=0.85),
            error_y=(dict(type="data", array=unc.fillna(0.0), visible=True, color=col,
                          thickness=0.9, width=0) if unc is not None else None),
            hovertemplate="%{x|%Y-%m-%d}<br>C_L=%{y:.4g}<extra></extra>"))
    if ref is not None:
        med = float(ref[1]["cal_value"].median())
        if np.isfinite(med) and med:
            fig.add_hrect(y0=med * 0.9, y1=med * 1.1, fillcolor="rgba(70,130,140,0.10)",
                          line_width=0, layer="below")
            fig.add_hline(y=med, line=dict(color="#2a6b73", width=1.2))
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                     name="±10 % of median",
                                     marker=dict(size=9, symbol="square",
                                                 color="rgba(70,130,140,0.25)")))
        last = ref[1]["datetime"].max()
        fig.add_vrect(x0=last - pd.Timedelta(days=15), x1=last,
                      fillcolor="rgba(240,173,78,0.16)", line_width=0, layer="below")
        fig.add_vline(x=last, line=dict(color="#d9534f", width=1.2))
    fig.update_layout(**_LAYOUT, title="Rayleigh vs cloud — lidar constant C_L",
                      yaxis_title="C_L", legend=dict(orientation="h", y=1.14))
    fig.update_yaxes(exponentformat="e")
    return fig


def monitoring_timeseries(hk_df: pd.DataFrame, itype: str | None = None) -> go.Figure:
    """Instrument monitoring: mean laser power/energy & window transmission (%) on the left axis and
    temperatures (degC) on the right axis, over time. Plots only the housekeeping fields the stream
    actually reports (others are blank in <key>_hk.csv and skipped); gaps are real downtime.

    The series are reindexed onto a complete uniform time grid and shipped as ``x0``/``dx`` instead of
    an explicit x array: with an explicit x, EVERY trace serialises its own full ISO-datetime string
    array, which measures 82.7 kB for five daily traces over a year (and would be 1.77 MB at hourly
    resolution). x0/dx + float32 y brings the same daily data to 19.7 kB and makes hourly resolution
    affordable (~252 kB). ``connectgaps=False`` keeps real downtime readable as a gap."""
    fig = go.Figure()
    has_temp = False
    t = pd.to_datetime(hk_df["datetime"])
    step = t.diff().dropna().min() if len(t) > 1 else pd.Timedelta(days=1)
    if pd.isna(step) or step <= pd.Timedelta(0):
        step = pd.Timedelta(days=1)
    grid = pd.date_range(t.min(), t.max(), freq=step)
    d = hk_df.set_index(t).reindex(grid)            # missing days -> NaN -> a real gap in the line
    x0 = float(grid[0].value // 10**6)              # epoch ms, as Plotly expects for a date axis
    dx = float(step.value // 10**6)
    for field, label, unit, group in config.HK_PANEL:
        if field not in d.columns:
            continue
        y = pd.to_numeric(d[field], errors="coerce")
        if not y.notna().any():
            continue
        if group == "temp":
            has_temp = True
        fig.add_trace(go.Scatter(
            x0=x0, dx=dx, y=y.to_numpy(dtype="float32"), mode="lines", name=f"{label} [{unit}]",
            line=dict(color=config.HK_COLORS.get(field), width=1.4),
            yaxis=("y2" if group == "temp" else "y"), connectgaps=False,
            hovertemplate="%{x|%Y-%m-%d}<br>" + label + "=%{y:.1f} " + unit + "<extra></extra>"))
    # No Plotly title: the station page already has an <h2> above this card.
    lay = {**_LAYOUT, "height": 300, "margin": dict(l=56, r=56, t=10, b=34),
           "legend": dict(orientation="h", y=1.16), "yaxis": dict(title="laser / window [%]"),
           "xaxis": dict(type="date")}
    if has_temp:
        lay["yaxis2"] = dict(title="temperature [degC]", overlaying="y", side="right", showgrid=False)
    fig.update_layout(**lay)
    return fig


def _discrete_colorscale(colors):
    """Piecewise-constant Plotly colorscale over z = 0..len(colors)-1 (each class gets a flat band)."""
    n = len(colors)
    out = []
    for i, c in enumerate(colors):
        out.append([i / n, c])
        out.append([(i + 1) / n, c])
    return out


def daily_availability_rows(status_df: pd.DataFrame, cal_df: pd.DataFrame | None = None,
                            methods=None) -> go.Figure | None:
    """Cloudnet-style daily strip, stacked: instrument status, mean cloud cover, and one calibration
    row per method (CL61 carries both Rayleigh and liquid-cloud, so it gets four rows).

    One cell per day over the whole record; record gaps read as 'no data'. Hover gives the decoded
    status string / the octa value / the exact calibration flag. Clicking a day drives the diagnostic
    viewer (wired in diag.js via the 'fig-avail' id) -- that contract reads only ``points[0].x``, so
    it is unaffected by the extra rows. Returns None when the stream has no decoded status history.

    Rows are separate ``go.Heatmap`` traces, NOT bars: a bar trace makes rangesync.js force
    ``yaxis.autorange`` on every period change (rangesync.js hasBars), which would scramble a fixed
    row stack, and Plotly allows only one colorscale per trace while the rows need three different
    ones. Each trace is placed with ``y0``/``dy`` on a NUMERIC axis that is then labelled by
    tickvals/ticktext: a one-row heatmap declaring a single value on a *category* axis has no
    defined cell height and Plotly draws nothing at all (verified -- the card came out blank).

    Built from ``<key>_status.csv`` (date, quality, summary, and -- once the producer ships it --
    mean_cloud_cover / cloud_cover_n / cloud_src) plus the station's ``<key>_cal.csv`` rows. The
    strip spans the full record (not the period selector) so it reads like Cloudnet's multi-year
    availability bar."""
    if status_df is None or not len(status_df):
        return None
    d = status_df.copy()
    d["dt"] = pd.to_datetime(d["date"].astype(str), format="%Y%m%d", errors="coerce")
    d = d.dropna(subset=["dt"]).sort_values("dt")
    if not len(d):
        return None

    # The x axis must span status AND calibration dates: a stream can carry calibration rows for days
    # with no decoded status (and vice versa), and clipping to the status range would silently drop
    # part of the calibration history.
    lo, hi = d["dt"].min(), d["dt"].max()
    cal_by_method = {}
    for m in (methods or []):
        if cal_df is None or not len(cal_df):
            continue
        c = cal_df[cal_df["method"] == m]
        if not len(c):
            continue
        c = c.assign(dt=pd.to_datetime(c["date"].astype(str), format="%Y%m%d", errors="coerce"))
        c = c.dropna(subset=["dt"])
        if not len(c):
            continue
        cal_by_method[m] = c
        lo, hi = min(lo, c["dt"].min()), max(hi, c["dt"].max())
    full = pd.date_range(lo, hi, freq="D")

    rows, traces = [], []

    # --- row 1: instrument status (unchanged semantics and palette) --------------------------------
    qmap = dict(zip(d["dt"], d["quality"].astype(str)))
    smap = (dict(zip(d["dt"], d["summary"].astype(str))) if "summary" in d.columns else {})
    q_order = ["pass", "warning", "error", "nodata"]
    q_idx = {q: i for i, q in enumerate(q_order)}
    quals = [qmap.get(t, "nodata") or "nodata" for t in full]
    rows.append("Instrument status")
    traces.append(go.Heatmap(
        x=list(full), y0=0, dy=1,
        z=[[q_idx.get(q, 3) for q in quals]],
        customdata=[[[config.QUALITY_LABELS.get(q, q),
                      (smap.get(t, "") if t in qmap else "No data") or ""]
                     for q, t in zip(quals, full)]],
        colorscale=_discrete_colorscale([config.QUALITY_COLORS[q] for q in q_order]),
        zmin=0, zmax=len(q_order), showscale=False,
        hovertemplate="%{x|%Y-%m-%d}<br><b>%{customdata[0]}</b>"
                      "<br>%{customdata[1]}<extra></extra>"))

    # --- row 2: mean cloud cover (octas) -----------------------------------------------------------
    # Absent until the producer writes mean_cloud_cover into <key>_status.csv; the row is simply not
    # drawn rather than drawn empty, so an old archive looks unchanged instead of looking broken.
    if "mean_cloud_cover" in d.columns:
        cc = pd.to_numeric(d["mean_cloud_cover"], errors="coerce")
        cmap = dict(zip(d["dt"], cc))
        src = (dict(zip(d["dt"], d["cloud_src"].astype(str))) if "cloud_src" in d.columns else {})
        vals = [cmap.get(t) for t in full]
        vals = [None if (v is None or not pd.notna(v)) else float(v) for v in vals]
        if any(v is not None for v in vals):
            traces.append(go.Heatmap(
                x=list(full), y0=len(rows), dy=1, z=[vals],
                customdata=[[[("cbh fraction" if src.get(t) == "cbh" else "cloud_amount")]
                             for t in full]],
                colorscale=config.CLOUD_COVER_SCALE, zmin=0, zmax=8, showscale=False,
                hovertemplate="%{x|%Y-%m-%d}<br><b>%{z:.1f} octas</b>"
                              "<br>%{customdata[0]}<extra></extra>"))
            rows.append("Mean cloud cover")

    # --- rows 3..N: one calibration outcome row per method ----------------------------------------
    cls_order = config.CAL_CLASS_ORDER
    cls_idx = {c: i for i, c in enumerate(cls_order)}
    cls_scale = _discrete_colorscale([config.CAL_CLASS_COLORS[c] for c in cls_order])
    for m, c in cal_by_method.items():
        # one row per (date, method); keep the last row if a day somehow carries duplicates
        fmap = {t: f for t, f in zip(c["dt"], c["flag"])}
        z, cd = [], []
        for t in full:
            f = fmap.get(t)
            if f is None or not pd.notna(f):
                z.append(None)
                cd.append(["—", "no calibration row"])
                continue
            k = config.cal_class(f)
            z.append(cls_idx[k])
            cd.append([config.CAL_CLASS_LABELS[k], config.flag_label(f, m)])
        traces.append(go.Heatmap(
            x=list(full), y0=len(rows), dy=1, z=[z], customdata=[cd],
            colorscale=cls_scale, zmin=0, zmax=len(cls_order), showscale=False,
            hovertemplate="%{x|%Y-%m-%d}<br><b>%{customdata[0]}</b>"
                          "<br>%{customdata[1]}<extra></extra>"))
        rows.append(f"Calibration — {config.method_label(m)}")

    fig = go.Figure(traces)
    # The HTML card already carries an <h2>; a Plotly title here would duplicate it.
    # Descending range = row 0 on top, without relying on autorange="reversed" (which the period
    # selector would fight over).
    fig.update_layout(**{**_LAYOUT, "height": 44 * len(rows) + 46,
                         "margin": dict(l=136, r=10, t=8, b=34)},
                      showlegend=False,
                      yaxis=dict(tickmode="array", tickvals=list(range(len(rows))), ticktext=rows,
                                 range=[len(rows) - 0.5, -0.5], fixedrange=True,
                                 showgrid=False, zeroline=False, ticksuffix="  "),
                      xaxis=dict(title="", type="date"))
    return fig


# --- Cheap inline sparkline -------------------------------------------------

def sparkline_svg(values, width: int = 110, height: int = 26, color: str = "#1f77b4") -> str:
    """Tiny inline-SVG line of the last values, normalized. Empty string if too few."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) < 3:
        return ""
    lo, hi = float(v.min()), float(v.max())
    span = (hi - lo) or 1.0
    n = len(v)
    pts = [f"{(i/(n-1))*(width-2)+1:.1f},{height-1-((y-lo)/span)*(height-2):.1f}"
           for i, y in enumerate(v)]
    return (f'<svg class="spark" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
            f'<polyline fill="none" stroke="{color}" stroke-width="1.2" points="{" ".join(pts)}"/></svg>')
