"""Plotly figure builders (+ a cheap inline-SVG sparkline for the station table).

Every figure is emitted as a `<div>` via fig_to_div(); the shared plotly.min.js is loaded
once per page (see render.py), so figures use include_plotlyjs=False.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from monitoring import config, kalman

# Compact, consistent styling for every figure.
_LAYOUT = dict(
    template="plotly_white",
    margin=dict(l=60, r=20, t=40, b=40),
    font=dict(size=12),
    height=340,
)


def fig_to_div(fig: go.Figure, div_id: str) -> str:
    """Render a figure to a standalone <div> that reuses the page-global Plotly."""
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id=div_id,
                       config={"displaylogo": False, "responsive": True})


# --- Summary-page figures ---------------------------------------------------

def network_map(st: pd.DataFrame) -> go.Figure:
    """Station map over Europe, colored by overall success rate, sized by activity."""
    d = st.dropna(subset=["lat", "lon"]).copy()
    fig = go.Figure(go.Scattergeo(
        lat=d["lat"], lon=d["lon"],
        text=[f"{k}<br>{t} · {sr:.0f}% success · {nd} nights"
              for k, t, sr, nd in zip(d["key"], d["itype"],
                                      d["success_rate"].fillna(0), d["n_dates"].fillna(0))],
        hoverinfo="text",
        marker=dict(
            size=np.clip(np.sqrt(d["n_dates"].fillna(1)) / 2.0, 5, 18),
            color=d["success_rate"].fillna(0), colorscale="RdYlGn", cmin=0, cmax=80,
            colorbar=dict(title="success %"), line=dict(width=0.4, color="#555"),
        ),
    ))
    fig.update_geos(scope="europe", resolution=50, showcountries=True,
                    countrycolor="#bbbbbb", showland=True, landcolor="#f5f5f5",
                    lataxis_range=config.MAP_LAT_RANGE, lonaxis_range=config.MAP_LON_RANGE)
    fig.update_layout(**{**_LAYOUT, "height": 460, "margin": dict(l=0, r=0, t=40, b=0)},
                      title="Network — success rate by station")
    return fig


def success_by_type(by_type: pd.DataFrame) -> go.Figure:
    """Bar of overall success rate per instrument type."""
    d = by_type.copy()
    fig = go.Figure(go.Bar(
        x=d["itype"], y=d["success_rate"],
        marker_color=[config.TYPE_COLORS.get(t, "#888") for t in d["itype"]],
        text=[f"{v:.0f}%" for v in d["success_rate"]], textposition="outside",
        hovertext=[f"{n} stations · {s}/{t} nights"
                   for n, s, t in zip(d["n_stations"], d["n_success"], d["n_dates"])],
    ))
    fig.update_layout(**_LAYOUT, title="Success rate by instrument type",
                      yaxis_title="success %", yaxis_range=[0, 100])
    return fig


def flag_distribution_bar(flags: pd.DataFrame) -> go.Figure:
    """Horizontal bar of night counts per flag meaning (good at top)."""
    d = flags.iloc[::-1]  # reverse so 'Successful' renders at the top
    fig = go.Figure(go.Bar(
        x=d["count"], y=d["label"], orientation="h",
        marker_color=d["color"],
        text=d["count"], textposition="auto",
    ))
    fig.update_layout(**{**_LAYOUT, "height": 380, "margin": dict(l=220, r=20, t=40, b=40)},
                      title="Outcome distribution (all nights)", xaxis_title="nights")
    return fig


def cl_by_type_box(st: pd.DataFrame) -> go.Figure:
    """Box plot of per-station median C_L grouped by instrument type (log axis)."""
    fig = go.Figure()
    for t in config.TYPE_ORDER:
        vals = st.loc[(st["itype"] == t) & (st["median_cl"] > 0), "median_cl"]
        if len(vals):
            fig.add_trace(go.Box(y=vals, name=t, marker_color=config.TYPE_COLORS.get(t, "#888"),
                                 boxpoints="all", jitter=0.4, pointpos=0))
    fig.update_layout(**_LAYOUT, title="Per-station median C_L by type",
                      yaxis_title="C_L", yaxis_type="log", showlegend=False)
    fig.update_yaxes(exponentformat="e")  # scientific notation on the log axis
    return fig


# --- Station-page figures ---------------------------------------------------

def cl_timeseries(g: pd.DataFrame) -> go.Figure:
    """Lidar constant over time: successful nights + E-PROFILE Kalman best estimate.

    The red line/band is the operational random-walk Kalman best estimate (see
    monitoring/kalman.py) -- the same routine E-PROFILE uses for the reprocessed
    best-estimate, not a plain rolling median.
    """
    ok = g[g["success"] == 1].sort_values("datetime")
    fig = go.Figure()
    if len(ok):
        fig.add_trace(go.Scatter(
            x=ok["datetime"], y=ok["lidar_constant"], mode="markers",
            name="C_L (nightly)", marker=dict(size=5, color="#1f77b4", opacity=0.7),
            error_y=dict(type="data", array=ok["uncertainty"], visible=True,
                         thickness=0.6, width=0, color="rgba(31,119,180,0.3)"),
            hovertemplate="%{x|%Y-%m-%d}<br>C_L=%{y:.3e}<extra></extra>",
        ))
        # Operational Kalman best estimate + ±1σ band.
        kt, ks, kstd = kalman.kalman_best_estimate(ok["datetime"].tolist(),
                                                   ok["lidar_constant"].to_numpy())
        if kt.size:
            fig.add_trace(go.Scatter(
                x=np.concatenate([kt, kt[::-1]]),
                y=np.concatenate([ks + kstd, (ks - kstd)[::-1]]),
                fill="toself", fillcolor="rgba(214,39,40,0.12)",
                line=dict(width=0), hoverinfo="skip", showlegend=False))
            fig.add_trace(go.Scatter(x=kt, y=ks, mode="lines", name="Kalman best estimate",
                                     line=dict(color="#d62728", width=2),
                                     hovertemplate="%{x|%Y-%m-%d}<br>Kalman=%{y:.3e}<extra></extra>"))
    fig.update_layout(**_LAYOUT, title="Lidar constant C_L over time — nightly + Kalman",
                      yaxis_title="C_L", legend=dict(orientation="h", y=1.12))
    fig.update_yaxes(exponentformat="e")  # scientific notation (1e11), not SI ("100G")
    return fig


def monthly_flag_bars(g: pd.DataFrame) -> go.Figure:
    """Stacked monthly bars of outcome counts -- shows when a station goes bad."""
    d = g.copy()
    d["month"] = d["datetime"].dt.to_period("M").dt.to_timestamp()
    fig = go.Figure()
    # One trace per flag value present, stacked, ordered good-to-bad.
    flags_present = sorted(d["flag"].dropna().unique(), key=lambda f: -float(f))
    for f in flags_present:
        sub = d[d["flag"] == f]
        counts = sub.groupby("month").size()
        fig.add_trace(go.Bar(x=counts.index, y=counts.values, name=config.flag_label(f),
                             marker_color=config.flag_color(f)))
    fig.update_layout(**{**_LAYOUT, "height": 320}, barmode="stack",
                      title="Monthly outcomes", yaxis_title="nights",
                      legend=dict(orientation="h", y=-0.25, font=dict(size=10)))
    return fig


def window_timeseries(g: pd.DataFrame) -> go.Figure:
    """Rayleigh fit window (bottom/top height) over time for successful nights."""
    ok = g[(g["success"] == 1) & g["bottom_height"].notna()].sort_values("datetime")
    fig = go.Figure()
    if len(ok):
        fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["top_height"], mode="markers",
                                 name="top", marker=dict(size=4, color="#2ca02c")))
        fig.add_trace(go.Scatter(x=ok["datetime"], y=ok["bottom_height"], mode="markers",
                                 name="bottom", fill="tonexty", fillcolor="rgba(44,160,44,0.12)",
                                 marker=dict(size=4, color="#8c564b")))
    fig.update_layout(**{**_LAYOUT, "height": 300}, title="Calibration window (m AGL)",
                      yaxis_title="height (m)", legend=dict(orientation="h", y=1.12))
    return fig


# --- Cheap inline sparkline (no Plotly, one per table row) -------------------

def sparkline_svg(values, width: int = 110, height: int = 26) -> str:
    """Tiny inline-SVG line of the last C_L values, normalized. Empty string if too few."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if len(v) < 3:
        return ""
    lo, hi = float(v.min()), float(v.max())
    span = (hi - lo) or 1.0
    n = len(v)
    pts = []
    for i, y in enumerate(v):
        px = (i / (n - 1)) * (width - 2) + 1
        py = height - 1 - ((y - lo) / span) * (height - 2)
        pts.append(f"{px:.1f},{py:.1f}")
    return (f'<svg class="spark" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">'
            f'<polyline fill="none" stroke="#1f77b4" stroke-width="1.2" '
            f'points="{" ".join(pts)}"/></svg>')
