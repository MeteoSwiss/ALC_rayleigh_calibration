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
        hoverinfo="text", customdata=customdata,
        marker=dict(size=sizes,
                    color=np.asarray(d["success_rate"].fillna(0), dtype=float),
                    colorscale="RdYlGn", cmin=0, cmax=80,
                    colorbar=dict(title="success %"), line=dict(width=0.4, color="#555")),
    ))
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
        hoverinfo="text", customdata=customdata,
        marker=dict(size=sizes, color=vals, colorscale="RdBu", cmin=cmin, cmax=cmax, cmid=100.0,
                    colorbar=dict(title=cbar), line=dict(width=0.4, color="#555")),
    ))
    fig.update_geos(scope="europe", resolution=50, showcountries=True, countrycolor="#bbbbbb",
                    showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)}, title=title)
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
    per-station medians. ``d`` has one row per station: columns key, med, q1, q3, n.
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
        customdata=list(zip(d["key"], d["n"], d["q1"], d["q3"])),
        hovertemplate=("%{customdata[0]}<br>median C_L = %{y:.3g}"
                       "<br>IQR = [%{customdata[2]:.3g}, %{customdata[3]:.3g}]  (n=%{customdata[1]})"
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

def series_timeseries(g_m: pd.DataFrame, kal_m: pd.DataFrame, method: str,
                      op_df: pd.DataFrame | None = None) -> go.Figure:
    """Calibration value over time for ONE method: successes + uncertainty + Kalman best estimate,
    plus (optional) the daily OPERATIONAL calibration constant from the L2 files as a black line.

    The Kalman line/band is the operational E-PROFILE random-walk best estimate, preferring
    the precomputed series (kal_m) and falling back to an on-the-fly fit.
    """
    vname = _value_name(method)
    color = config.METHOD_COLORS.get(method, "#1f77b4")
    ok = g_m[g_m["success"] == 1].sort_values("datetime")
    fig = go.Figure()
    if op_df is not None and len(op_df):
        od = op_df.sort_values("datetime")
        fig.add_trace(go.Scatter(
            x=od["datetime"], y=od["op_coeff"], mode="lines", name="Operational constant (L2)",
            line=dict(color="#111111", width=1.3),
            hovertemplate="%{x|%Y-%m-%d}<br>operational=%{y:.3e}<extra></extra>"))
    if len(ok):
        fig.add_trace(go.Scatter(
            x=ok["datetime"], y=ok["cal_value"], mode="markers", name=f"{vname} (per cal)",
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
            fig.add_trace(go.Scatter(x=kt, y=ks, mode="lines", name="Kalman best estimate",
                                     line=dict(color="#d62728", width=2),
                                     hovertemplate="%{x|%Y-%m-%d}<br>Kalman=%{y:.3e}<extra></extra>"))
    fig.update_layout(**_LAYOUT, yaxis_title=vname, legend=dict(orientation="h", y=1.14),
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
            fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["top_height"], mode="markers",
                                     name="top", marker=dict(size=4, color="#2ca02c")))
            fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["bottom_height"], mode="markers",
                                     name="bottom", fill="tonexty", fillcolor="rgba(44,160,44,0.12)",
                                     marker=dict(size=4, color="#8c564b")))
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
    """Overlay the lidar constant C_L of each method (CL61) on one axis -- a direct cross-check.

    Both methods estimate the SAME C_L (Wiegner), so the Rayleigh and cloud points should
    agree; no normalization, the absolute C_L is the useful comparison.
    """
    fig = go.Figure()
    for method, g_m in by_method.items():
        ok = g_m[g_m["success"] == 1].sort_values("datetime")
        if not len(ok):
            continue
        fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["cal_value"], mode="markers",
                                 name=config.method_label(method),
                                 marker=dict(size=5, color=config.METHOD_COLORS.get(method, "#888"),
                                             opacity=0.8),
                                 hovertemplate="%{x|%Y-%m-%d}<br>C_L=%{y:.3e}<extra></extra>"))
    fig.update_layout(**_LAYOUT, title="Rayleigh vs cloud — lidar constant C_L",
                      yaxis_title="C_L", legend=dict(orientation="h", y=1.14))
    fig.update_yaxes(exponentformat="e")
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
