"""Render the static dashboard: summary page (index.html) + one page per station.

Reads the SQLite index, builds Plotly figures, fills Jinja2 templates, and writes a
self-contained site (shared plotly.min.js + css under assets/). All links are relative.
Each station may carry one or two method series (Rayleigh and/or cloud).
"""
from __future__ import annotations

import json
import os
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
    env.globals["flag_anchor"] = config.flag_anchor
    env.globals["method_label"] = config.method_label
    return env


def _copy_flag_examples(flagex_dir, out_dir: Path) -> dict:
    """Copy curated per-flag example PNGs into the site. Source files are named
    '<anchor>__<caption-with-underscores>.png' (e.g. 'm1__cloud_no_liquid.png'); returns
    {flag_value: [ {rel, caption}, ... ]} for the explanation page."""
    by: dict = {}
    if not flagex_dir:
        return by
    src = Path(flagex_dir)
    if not src.exists():
        return by
    anchor_to_val = {config.flag_anchor(d["value"]): d["value"] for d in config.FLAG_DOCS}
    dst = out_dir / "flagex"
    for png in sorted(src.glob("*.png")):
        anchor = png.stem.split("__", 1)[0]
        val = anchor_to_val.get(anchor)
        if val is None:
            continue
        dst.mkdir(parents=True, exist_ok=True)
        try:
            shutil.copyfile(png, dst / png.name)
        except OSError:
            continue
        cap = png.stem.split("__", 1)[1].replace("_", " ") if "__" in png.stem else ""
        by.setdefault(val, []).append({"rel": f"flagex/{png.name}", "caption": cap})
    return by


def _write_assets(out_dir: Path) -> str | None:
    assets = out_dir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    (assets / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")
    for name in ("style.css", "table-sort.js", "paginate.js", "filter.js", "qcflag.js", "diag.js",
                 "histlink.js", "search.js"):
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
    """Per-key rollup (across methods) for the network map. Success rate excludes no-data /
    unsuitable-conditions days (uses n_suitable), matching the per-series definition."""
    agg = (series.groupby("key")
           .agg(n_dates=("n_dates", "sum"), n_success=("n_success", "sum"),
                n_suitable=("n_suitable", "sum"),
                methods=("method", lambda s: " + ".join(config.method_label(m) for m in sorted(set(s)))))
           .reset_index())
    agg["success_rate"] = 100.0 * agg["n_success"] / agg["n_suitable"].replace(0, np.nan)
    cols = [c for c in ("key", "itype", "lat", "lon", "country", "name") if c in st.columns]
    return agg.merge(st[cols], on="key", how="left")


def _load_opcoeff(csv) -> pd.DataFrame | None:
    """Load the operational calibration-constant CSV (key,date,op_coeff) produced by
    scripts/extract_l2_opcoeff.py; add a parsed datetime. Returns None if absent."""
    if not csv:
        return None
    p = Path(csv)
    if not p.exists():
        return None
    df = pd.read_csv(p, dtype={"key": str, "date": str})
    df["op_coeff"] = pd.to_numeric(df["op_coeff"], errors="coerce")
    df["datetime"] = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    return df


def _opcoeff_ratios(cal: pd.DataFrame, op_all, st: pd.DataFrame) -> pd.DataFrame:
    """Per-station median C_L expressed as a percent of (a) the theoretical value and (b) the
    operational L2 constant on the day of each calibration. Uses successful calibrations only;
    the operational ratio is a per-calibration median (robust to occasional default/garbage op)."""
    type_by_key = dict(zip(st["key"], st["itype"]))
    succ = cal[cal["success"] == 1]
    rows = []
    for key, g in succ.groupby("key"):
        cls = pd.to_numeric(g["cal_value"], errors="coerce")
        cls = cls[cls > 0]
        if not len(cls):
            continue
        median_cl = float(cls.median())
        theo = config.theoretical_cl(type_by_key.get(key))
        pct_theo = 100.0 * median_cl / theo if (theo and theo > 0) else np.nan
        pct_op = np.nan
        if op_all is not None:
            gop = g.merge(op_all[op_all["key"] == key][["date", "op_coeff"]], on="date", how="inner")
            if len(gop):
                rr = pd.to_numeric(gop["cal_value"], errors="coerce") / pd.to_numeric(gop["op_coeff"], errors="coerce")
                rr = rr[np.isfinite(rr) & (rr > 0)]
                if len(rr):
                    pct_op = 100.0 * float(rr.median())
        rows.append(dict(key=key, op_pct_theo=pct_theo, op_pct_op=pct_op, op_median_cl=median_cl))
    return pd.DataFrame(rows, columns=["key", "op_pct_theo", "op_pct_op", "op_median_cl"])


def _op_station_df(op_all, key, ref_cl):
    """Daily operational-constant series for one station's time-series overlay, with implausible
    spikes (operational defaults like 1e8 on a CL61) dropped relative to the station's C_L scale."""
    if op_all is None:
        return None
    d = op_all[op_all["key"] == key]
    if not len(d):
        return None
    op = pd.to_numeric(d["op_coeff"], errors="coerce").to_numpy()
    keep = np.isfinite(op) & (op > 0)
    ref = ref_cl if (ref_cl and np.isfinite(ref_cl) and ref_cl > 0) else (
        float(np.median(op[keep])) if keep.any() else None)
    if ref and np.isfinite(ref) and ref > 0:
        keep = keep & (op >= ref / 30.0) & (op <= ref * 30.0)
    if not keep.any():
        return None
    return pd.DataFrame({"datetime": d["datetime"].to_numpy()[keep], "op_coeff": op[keep]})


def _series_table_rows(cal: pd.DataFrame, series: pd.DataFrame, st: pd.DataFrame) -> list[dict]:
    """One row per (station, method) with a value sparkline + method badge + country (for filtering)."""
    last_by = {(k, m): g[g["success"] == 1].sort_values("datetime")["cal_value"].tail(30).tolist()
               for (k, m), g in cal.groupby(["key", "method"], sort=False)}
    country_by = dict(zip(st["key"], st["country"])) if "country" in st.columns else {}
    rows = []
    for _, s in series.iterrows():
        method = s["method"]
        spark = charts.sparkline_svg(last_by.get((s["key"], method), []),
                                     color=config.METHOD_COLORS.get(method, "#1f77b4"))
        rows.append(dict(
            key=s["key"], itype=s.get("itype"), method=method, method_label=config.method_label(method),
            country=str(country_by.get(s["key"], "") or ""),
            success_rate=s.get("success_rate"), n_dates=s.get("n_dates"), n_success=s.get("n_success"),
            median_cl=s.get("median_cl"), median_rel_unc=s.get("median_rel_unc"),
            last_date=s.get("last_date"), last_flag=s.get("last_flag"), spark=spark,
        ))
    return rows


# How per-night diagnostic PNGs are placed under the site's diag/ folder. Default 'symlink' makes the
# site reference the originals WITHOUT duplicating data -- the diagnostic set can be ~100k images /
# tens of GB and already lives next to the site on the same filesystem, so copying it is wasteful.
# 'hardlink' is similar but survives a plain rsync; 'copy' duplicates the bytes (use only when the
# site must move to a different filesystem -- e.g. rsync'd to a laptop without -L). python's
# http.server follows symlinks, so the SSH-tunnel viewer works with the default.
_DIAG_LINK_MODE = os.environ.get("ALC_DIAG_LINK", "symlink").lower()


def _materialize(src: Path, dst: Path, mode: str) -> None:
    """Make *dst* resolve to *src* as a symlink (default), hardlink, or copy. Idempotent and cheap;
    falls back to copy if linking is unsupported (cross-filesystem, or Windows without privilege)."""
    if mode == "copy":
        if dst.exists() and not dst.is_symlink() and dst.stat().st_size == src.stat().st_size:
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_symlink():
            dst.unlink()
        shutil.copyfile(src, dst)
        return
    # link modes: (re)create the link -- a metadata op, effectively free vs copying the bytes
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink() or dst.exists():
        try:
            dst.unlink()
        except OSError:
            pass
    try:
        if mode == "hardlink":
            os.link(src, dst)
        else:
            os.symlink(os.path.abspath(src), dst)
        return
    except OSError:
        shutil.copyfile(src, dst)


def _copy_diagnostics(diag: pd.DataFrame, cal: pd.DataFrame, out_dir: Path) -> dict:
    """Place per-calibration diagnostic PNGs under the site (diag/<key>/<method>_<date>.png) and
    return {(key, method): [ {date, rel, success}, ... ]} sorted by date for the station-page viewer.
    By default the PNGs are SYMLINKED, not copied, so a ~100k-image / tens-of-GB diagnostic set is
    not duplicated (set ALC_DIAG_LINK=copy for a portable site, or =hardlink). `success` (from cal)
    drives the green/grey calendar and the valid-only (left/right) navigation."""
    by: dict = {}
    if not len(diag):
        return by
    succ = set()
    if len(cal):
        s = cal[cal["success"] == 1]
        succ = set(zip(s["key"].astype(str), s["method"].astype(str), s["date"].astype(str)))
    for _, r in diag.iterrows():
        src = Path(str(r["src"]))
        if not src.exists():
            continue
        key, method, date = str(r["key"]), str(r["method"]), str(r["date"])
        fname = f"{method}_{date}.png"
        dst = out_dir / "diag" / key / fname
        try:
            _materialize(src, dst, _DIAG_LINK_MODE)
        except OSError:
            continue
        by.setdefault((key, method), []).append(
            {"date": date, "rel": f"diag/{key}/{fname}", "success": (key, method, date) in succ})
    for k in by:
        by[k].sort(key=lambda x: x["date"])
    return by


def _method_block(key, method, cal, kal, series, diags=None, op_all=None):
    """Figures + aggregates for one method section on a station page."""
    g_m = cal[(cal["key"] == key) & (cal["method"] == method)].sort_values("datetime")
    kal_m = kal[(kal["key"] == key) & (kal["method"] == method)] if len(kal) else kal
    srow = series[(series["key"] == key) & (series["method"] == method)]
    meta = srow.iloc[0].to_dict() if len(srow) else {}
    safe = method  # 'rayleigh'/'cloud' are id-safe
    ref = pd.to_numeric(g_m[g_m["success"] == 1]["cal_value"], errors="coerce")
    ref = ref[ref > 0]
    op_df = _op_station_df(op_all, key, float(ref.median()) if len(ref) else None)
    return dict(
        method=method, label=config.method_label(method), meta=meta,
        figs={
            "ts": charts.fig_to_div(charts.series_timeseries(g_m, kal_m, method, op_df), f"fig-ts-{safe}"),
            "flags": charts.fig_to_div(charts.monthly_flag_bars(g_m, method), f"fig-mf-{safe}"),
            "aux": charts.fig_to_div(charts.aux_timeseries(g_m, method), f"fig-aux-{safe}"),
        },
        recent=g_m.iloc[::-1].to_dict("records"),          # full archive, newest first (paginated)
        diag_dates=sorted({d["date"] for d in (diags or [])}),  # dates that have a diagnostic image
        diags=diags or [],
    )


def build_site(db_path: Path, out_dir: Path, limit_pages: int | None = None,
               flagex_dir=None, opcoeff_csv=None, only_keys=None) -> dict:
    out_dir = Path(out_dir)
    (out_dir / "stations").mkdir(parents=True, exist_ok=True)
    logo = _write_assets(out_dir)
    env = _env()

    cal, series, st, kal, diag = metrics.load_frames(db_path)
    summary = metrics.network_summary(cal, series, st)
    flags = metrics.flag_distribution(cal)
    watch = metrics.watchlist(cal, st)
    keystats = _keystats(series, st)
    diag_by = _copy_diagnostics(diag, cal, out_dir)

    # Operational calibration constant from the L2 files (optional): two ratio maps + per-station
    # black line on the time series.
    op_all = _load_opcoeff(opcoeff_csv)
    keystats = keystats.merge(_opcoeff_ratios(cal, op_all, st), on="key", how="left")

    summary_figs = {
        "map_theo": charts.fig_to_div(charts.ratio_map(
            keystats, "op_pct_theo", "Median C_L — % of theoretical value",
            "% of theoretical", "fig-map-theo"), "fig-map-theo"),
        "map_op": charts.fig_to_div(charts.ratio_map(
            keystats, "op_pct_op", "Median C_L — % of operational constant (L2)",
            "% of operational", "fig-map-op"), "fig-map-op"),
        "map": charts.fig_to_div(charts.network_map(keystats), "fig-map"),
        "success_type": charts.fig_to_div(charts.success_by_type_method(summary["by_type_method"]), "fig-stype"),
        "flag_dist_rayleigh": charts.fig_to_div(charts.flag_distribution_bar(flags, "rayleigh"), "fig-flags-r"),
        "flag_dist_cloud": charts.fig_to_div(charts.flag_distribution_bar(flags, "cloud"), "fig-flags-c"),
        "cl_type_abs": charts.fig_to_div(charts.value_by_type_method_box(series), "fig-cltype"),
        "cl_type_pct": charts.fig_to_div(charts.value_pct_theoretical_box(series), "fig-cltype-pct"),
    }
    # Per-station median C_L with IQR (Q1..Q3 of that station's successful daily values), one ranked
    # plot per instrument type. Pools both methods per station -- C_L is the same physical quantity.
    key_itype = dict(zip(st["key"], st["itype"]))
    okc = cal[(cal["success"] == 1) & (cal["cal_value"] > 0)].copy()
    okc["itype"] = okc["key"].map(key_itype)
    okc = okc.dropna(subset=["itype"])
    gb = okc.groupby(["itype", "key"])["cal_value"]
    sta_iqr = pd.DataFrame({"med": gb.median(), "q1": gb.quantile(0.25),
                            "q3": gb.quantile(0.75), "n": gb.size()}).reset_index()
    cl_iqr_figs = [(t, charts.fig_to_div(charts.cl_median_iqr_by_station(sta_iqr[sta_iqr["itype"] == t], t),
                                         f"fig-cliqr-{t}"))
                   for t in config.TYPE_ORDER if (sta_iqr["itype"] == t).any()]

    # Search index for the nav-bar station search (name + WIGOS id + key, all matchable).
    search_records = []
    for _, r in st.iterrows():
        k = str(r["key"])
        search_records.append({
            "key": k,
            "name": str(r.get("name", "") or ""),
            "wigos": k.rsplit("_", 1)[0] if "_" in k else k,  # drop the _A/_B/_C suffix
            "type": str(r.get("itype", "") or ""),
            "country": str(r.get("country", "") or ""),
        })
    search_json = json.dumps(search_records, ensure_ascii=False)

    countries = sorted({str(c) for c in st.get("country", pd.Series(dtype=str)).dropna()
                        if str(c).strip()})
    types = [t for t in config.TYPE_ORDER if t in set(st["itype"])] + \
            sorted(set(st["itype"]) - set(config.TYPE_ORDER) - {"Unknown"}) + \
            (["Unknown"] if "Unknown" in set(st["itype"]) else [])
    summary_html = env.get_template("summary.html").render(
        base="", logo=logo, summary=summary, figs=summary_figs, cl_iqr=cl_iqr_figs,
        watch=watch.to_dict("records"), rows=_series_table_rows(cal, series, st),
        countries=countries, types=types, search_json=search_json,
    )
    (out_dir / "index.html").write_text(summary_html, encoding="utf-8")

    # --- Flag explanation page (flags.html) ----------------------------------
    flag_examples = _copy_flag_examples(flagex_dir, out_dir)
    flags_html = env.get_template("flags.html").render(
        base="", logo=logo, flag_docs=config.FLAG_DOCS, flag_examples=flag_examples,
        search_json=search_json,
    )
    (out_dir / "flags.html").write_text(flags_html, encoding="utf-8")

    # --- Per-station pages (one per key; all of that key's methods) ----------
    keys = list(st["key"])
    if only_keys is not None:
        # incremental rebuild: re-render only the changed stations (the summary above always rebuilds);
        # unchanged station pages keep their existing HTML on disk
        only = set(only_keys)
        keys = [k for k in keys if k in only]
    if limit_pages:
        keys = keys[:limit_pages]
    station_tmpl = env.get_template("station.html")
    # Full station order (independent of only_keys/limit) so each page's up/down "previous/next
    # station" links always point at real neighbours, even on an incremental rebuild.
    all_keys = list(st["key"])
    nav_idx = {k: i for i, k in enumerate(all_keys)}
    for key in keys:
        meta = st[st["key"] == key].iloc[0].to_dict()
        methods = [m for m in config.METHOD_ORDER
                   if len(cal[(cal["key"] == key) & (cal["method"] == m)])]
        blocks = [_method_block(key, m, cal, kal, series, diag_by.get((key, m), []), op_all) for m in methods]
        overlay = None
        if len(methods) >= 2:
            by_method = {m: cal[(cal["key"] == key) & (cal["method"] == m)] for m in methods}
            overlay = charts.fig_to_div(charts.cl_overlay(by_method), "fig-overlay")
        i = nav_idx.get(key)
        prev_station = f"{all_keys[i - 1]}.html" if (i is not None and i > 0) else ""
        next_station = f"{all_keys[i + 1]}.html" if (i is not None and i < len(all_keys) - 1) else ""
        html = station_tmpl.render(base="../", logo=logo, key=key, meta=meta,
                                   blocks=blocks, overlay=overlay, search_json=search_json,
                                   prev_station=prev_station, next_station=next_station)
        (out_dir / "stations" / f"{key}.html").write_text(html, encoding="utf-8")

    return dict(out_dir=str(out_dir), n_pages=len(keys), n_series=int(len(series)),
                as_of=summary["as_of"])
