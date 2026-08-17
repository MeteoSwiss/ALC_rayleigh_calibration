"""Render the static dashboard: summary page (index.html) + one page per station.

Reads the SQLite index, builds Plotly figures, fills Jinja2 templates, and writes a
self-contained site (shared plotly.min.js + css under assets/). All links are relative.
Each station may carry one or two method series (Rayleigh and/or cloud).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from monitoring import charts, config, index, metrics, periods as periods_mod

_TEMPLATES = Path(__file__).parent / "templates"
_STATIC = Path(__file__).parent / "static"

# Project JS/CSS that get a cache-busting ?v=<hash> so browsers always pick up changes (Plotly is
# stable + large -> left unversioned). Same list is used to copy them into the site.
_VERSIONED_ASSETS = ("style.css", "table-sort.js", "paginate.js", "filter.js", "qcflag.js",
                     "diag.js", "histlink.js", "search.js", "rangesync.js",
                     "stationindex.js", "stationnav.js")


def _asset_version() -> str:
    """Short content hash of the versioned assets -> cache-busting token (?v=). Changes only when a
    JS/CSS file changes, so data-only rebuilds keep URLs cacheable but code changes are never stale."""
    h = hashlib.md5()
    for name in _VERSIONED_ASSETS:
        p = _STATIC / name
        if p.exists():
            h.update(p.read_bytes())
    return h.hexdigest()[:8]


def _write_if_changed(path: Path, text: str) -> bool:
    """Write *text* only when it differs from what is on disk. Returns True if the file was written.

    The publish step rsyncs by timestamp, so rewriting an identical file re-uploads it for nothing;
    this keeps a data file out of the daily transfer on the days it does not change."""
    try:
        if path.exists() and path.read_text(encoding="utf-8") == text:
            return False
    except OSError:
        pass
    path.write_text(text, encoding="utf-8")
    return True


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
    env.globals["asset_v"] = _asset_version()   # cache-busting token for ?v= on JS/CSS
    env.globals["has_favicon"] = (_STATIC / "favicon.png").exists()   # EUMETNET tab icon
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
        by.setdefault(val, []).append({"rel": f"{config.IMG_BASE_URL}flagex/{png.name}", "caption": cap})
    return by


def _write_assets(out_dir: Path) -> str | None:
    assets = out_dir / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    (assets / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")
    for name in _VERSIONED_ASSETS:
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
    # Dedicated square favicon (EUMETNET logo) — copied alongside so the browser tab shows it even
    # when the header logo is a wide wordmark; referenced explicitly in base.html.
    fav = _STATIC / "favicon.png"
    if fav.exists():
        shutil.copyfile(fav, assets / "favicon.png")
    return logo


def _keystats(series: pd.DataFrame, st: pd.DataFrame) -> pd.DataFrame:
    """Per-key rollup (across methods) for the network map. Success rate = valid / ALL days
    (n_dates): no-data, no-cloud/not-clear-night and rejections all count against it, matching
    the per-series definition."""
    agg = (series.groupby("key")
           .agg(n_dates=("n_dates", "sum"), n_success=("n_success", "sum"),
                n_suitable=("n_suitable", "sum"),
                methods=("method", lambda s: " + ".join(config.method_label(m) for m in sorted(set(s)))))
           .reset_index())
    agg["success_rate"] = 100.0 * agg["n_success"] / agg["n_dates"].replace(0, np.nan)
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


def _load_oldray(oldray_dir) -> dict | None:
    """Old operational Rayleigh calibrations from a tree of yearly NetCDFs
    (``ALC_calibration_<key><YYYY>.nc`` under *oldray_dir*, e.g. /scratch/mch/mhrvo/Calib_oper).
    Returns {key: DataFrame(datetime, value)} of the Rayleigh (calibration_method==0) lidar constant
    for the station-page time-series overlay, or None if the dir is absent / xarray is unavailable.
    Read with a lazy xarray import so a minimal dashboard env still builds when the overlay is off."""
    if not oldray_dir:
        return None
    d = Path(oldray_dir)
    if not d.exists():
        return None
    try:
        import xarray as xr
    except Exception:
        print(f"  oldray: xarray unavailable -> skipping old-Rayleigh overlay ({oldray_dir})", flush=True)
        return None
    import re
    by: dict = {}
    files = sorted(d.glob("**/ALC_calibration_*.nc"))
    for f in files:
        m = re.match(r"ALC_calibration_(.+)\d{4}\.nc$", f.name)
        if not m:
            continue
        key = m.group(1)
        try:
            ds = xr.open_dataset(f)
            t = ds["time"].values
            v = np.asarray(ds["lidar_constant"].values, dtype=float)
            meth = (np.asarray(ds["calibration_method"].values, dtype=float)
                    if "calibration_method" in ds.variables else np.zeros(v.shape))
            ds.close()
        except Exception:
            continue
        keep = np.isfinite(v) & (v > 0) & (meth == 0)
        if not keep.any():
            continue
        by.setdefault(key, []).append(pd.DataFrame({"datetime": pd.to_datetime(t[keep]), "value": v[keep]}))
    if not by:
        return None
    out = {k: pd.concat(v, ignore_index=True).sort_values("datetime") for k, v in by.items()}
    print(f"  oldray: old-Rayleigh overlay for {len(out)} stations (from {len(files)} files)", flush=True)
    return out


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


def _existing_diag(diag: pd.DataFrame) -> pd.DataFrame:
    """Subset of *diag* whose source PNG still exists on disk (vectorized stat). Shared by the
    materialize + index steps so the (large) 100k-row table is filtered only once.

    In IMAGES_IN_BUCKET mode the local PNGs may have been deleted (they live in the bucket and the
    site references their bucket URLs), so the existence filter is skipped and the full data-driven
    diag list is kept for the INDEX. _copy_diagnostics then stages a link only for a source that
    still exists (removing any stale one whose source was pruned), so fresh PNGs are published while
    deleted ones stay in the bucket without leaving a dangling stage link."""
    if not len(diag):
        return diag
    if config.IMAGES_IN_BUCKET:
        return diag
    src = diag["src"].astype(str).to_numpy()
    exists = np.fromiter((os.path.exists(s) for s in src), dtype=bool, count=src.shape[0])
    return diag[exists]


def _diag_index_from(d: pd.DataFrame, cal: pd.DataFrame) -> dict:
    """Build the viewer index {(key, method): [ {date, rel, success}, ... ]} from *d* (already
    filtered to existing sources). Vectorized (no per-row Python loop) -- the prior iterrows form
    cost ~23 s on the 146k-row network table; this is ~1 s."""
    by: dict = {}
    if not len(d):
        return by
    d = d[["key", "method", "date"]].astype(str).copy()
    d["rel"] = config.IMG_BASE_URL + "diag/" + d["key"] + "/" + d["method"] + "_" + d["date"] + ".png"
    if len(cal):
        s = cal[cal["success"] == 1]
        succ = set(zip(s["key"].astype(str), s["method"].astype(str), s["date"].astype(str)))
        d["success"] = [t in succ for t in zip(d["key"], d["method"], d["date"])]
    else:
        d["success"] = False
    # Classification is a product, not a pass/fail calibration -> every classified day is "valid"
    # (so the day-picker shows it available, not greyed like a rejected calibration).
    d.loc[d["method"] == "classification", "success"] = True
    d = d.sort_values("date")
    for (key, method), g in d.groupby(["key", "method"], sort=False):
        by[(key, method)] = g[["date", "rel", "success"]].to_dict("records")
    return by


def _diag_index(diag: pd.DataFrame, cal: pd.DataFrame) -> dict:
    """Read-only viewer index (no materialize), used by the parallel render workers that rely on the
    parent having already materialized the PNGs. See :func:`_diag_index_from`."""
    return _diag_index_from(_existing_diag(diag), cal)


def _copy_diagnostics(diag: pd.DataFrame, cal: pd.DataFrame, out_dir: Path) -> dict:
    """Place per-calibration diagnostic PNGs under the site (diag/<key>/<method>_<date>.png) and
    return the viewer index (see :func:`_diag_index_from`). By default the PNGs are SYMLINKED, not
    copied, so a ~100k-image / tens-of-GB diagnostic set is not duplicated (set ALC_DIAG_LINK=copy
    for a portable site, or =hardlink).

    Staging is DECOUPLED from indexing: the index keeps every diagnostic (in IMAGES_IN_BUCKET mode it
    points at the bucket URL, so history stays visible after the local PNG is pruned), but a link is
    STAGED only for a source that still exists. A source pruned after a prior upload is not re-staged,
    and any stale stage link whose source is now gone is REMOVED -- otherwise the build manufactures
    dangling symlinks that make ``aws s3 sync --follow-symlinks`` fail (exit 2). This runs over the
    whole diag table (diagnostics + classification curtains) every build, so it also clears links
    orphaned by the prune cycle."""
    d = _existing_diag(diag)
    for r in d.itertuples(index=False):
        dst = out_dir / "diag" / str(r.key) / f"{r.method}_{r.date}.png"
        if os.path.exists(str(r.src)):
            try:
                _materialize(Path(str(r.src)), dst, _DIAG_LINK_MODE)
            except OSError:
                continue
        elif dst.is_symlink():           # source pruned (image now in the bucket) -> drop the dead link
            try:
                dst.unlink()
            except OSError:
                pass
    return _diag_index_from(d, cal)


def _method_block(key, method, cal, kal, series, diags=None, op_all=None, oldray_all=None):
    """Figures + aggregates for one method section on a station page."""
    g_m = cal[(cal["key"] == key) & (cal["method"] == method)].sort_values("datetime")
    kal_m = kal[(kal["key"] == key) & (kal["method"] == method)] if len(kal) else kal
    srow = series[(series["key"] == key) & (series["method"] == method)]
    meta = srow.iloc[0].to_dict() if len(srow) else {}
    safe = method  # 'rayleigh'/'cloud' are id-safe
    ref = pd.to_numeric(g_m[g_m["success"] == 1]["cal_value"], errors="coerce")
    ref = ref[ref > 0]
    op_df = _op_station_df(op_all, key, float(ref.median()) if len(ref) else None)
    # Old operational Rayleigh overlay (black x) — Rayleigh series only.
    oldray_df = oldray_all.get(key) if (oldray_all is not None and method == "rayleigh") else None
    return dict(
        method=method, label=config.method_label(method), meta=meta,
        figs={
            "ts": charts.fig_to_div(charts.series_timeseries(g_m, kal_m, method, op_df, oldray_df), f"fig-ts-{safe}"),
            "flags": charts.fig_to_div(charts.monthly_flag_bars(g_m, method), f"fig-mf-{safe}"),
            "aux": charts.fig_to_div(charts.aux_timeseries(g_m, method), f"fig-aux-{safe}"),
        },
        recent=g_m.iloc[::-1].to_dict("records"),          # full archive, newest first (paginated)
        diag_dates=sorted({d["date"] for d in (diags or [])}),  # dates that have a diagnostic image
        diags=diags or [],
    )


def _period_row(df: pd.DataFrame, period: str):
    """Pick the summary row matching *period* from a (possibly period-keyed) OmB/sens CSV.

    Part 3 writes one row per period (a ``period`` column = 'all'/'y2025'/'last90'/…); older single-row
    CSVs have no ``period`` column → the lone row is used for every window (full-archive fallback). When
    a period has no dedicated row, fall back to the 'all' row, then to the first row."""
    if "period" not in df.columns:
        return df.iloc[0]
    pcol = df["period"].astype(str)
    for want in (str(period), "all"):
        m = df[pcol == want]
        if len(m):
            return m.iloc[0]
    return df.iloc[0]


def _ombsens_keystats(fullcal_dir, period: str = "all") -> pd.DataFrame:
    """Scan the fullcal dir for the per-stream OmB / sensitivity summaries written by the runner and
    roll them up into a per-key table for the two summary maps. Only a subset of streams have these
    files (a subset run is still populating them), so missing files are skipped. ``period`` selects the
    matching row from the period-keyed CSVs (full-archive fallback for legacy single-row files).

      <key>_omb.csv  -> omb_bias = median_bias_ours * 1e6   [Mm^-1 sr^-1]
      <key>_sens.csv -> icao_alt = icao_alt_200             [m] (blank when never detected)

    Returns columns ["key", "omb_bias", "icao_alt"] (NaN where a value is absent)."""
    rows: dict = {}
    if not fullcal_dir:
        return pd.DataFrame(columns=["key", "omb_bias", "icao_alt"])
    root = Path(fullcal_dir)
    if not root.exists():
        return pd.DataFrame(columns=["key", "omb_bias", "icao_alt"])
    for omb in root.glob("*/*_omb.csv"):
        key = omb.parent.name
        try:
            df = pd.read_csv(omb)
        except Exception:
            continue
        if not len(df):
            continue
        v = pd.to_numeric(_period_row(df, period).get("median_bias_ours"), errors="coerce")
        rows.setdefault(key, {})["omb_bias"] = (float(v) * 1e6) if np.isfinite(v) else np.nan
    for sens in root.glob("*/*_sens.csv"):
        key = sens.parent.name
        try:
            df = pd.read_csv(sens)
        except Exception:
            continue
        if not len(df):
            continue
        v = pd.to_numeric(_period_row(df, period).get("icao_alt_200"), errors="coerce")
        rows.setdefault(key, {})["icao_alt"] = float(v) if np.isfinite(v) else np.nan
    if not rows:
        return pd.DataFrame(columns=["key", "omb_bias", "icao_alt"])
    out = pd.DataFrame([{"key": k, "omb_bias": d.get("omb_bias", np.nan),
                         "icao_alt": d.get("icao_alt", np.nan)} for k, d in rows.items()])
    return out


def _stage_ombsens_pngs(fullcal_dir, key, out_dir: Path) -> dict:
    """Stage the per-station OmB / sensitivity diagnostic PNGs into the site (under ombsens/<key>/)
    the same way per-night diagnostics are staged (symlink by default, ALC_DIAG_LINK overrides), and
    return {'omb': rel, 'sens': rel} for the present files only (missing files are simply omitted)."""
    out: dict = {}
    if not fullcal_dir:
        return out
    base = Path(fullcal_dir) / key

    def _url_for(fname):
        return (f"{config.IMG_BASE_URL}ombsens/{key}/{fname}" if config.IMG_BASE_URL
                else f"../ombsens/{key}/{fname}")

    def _stage(fname):
        """Stage a local PNG if present and return its URL. In bucket mode the URL is emitted even when
        the local PNG is gone (it lives in the bucket); otherwise the file must exist on disk."""
        src = base / fname
        if src.exists():
            try:
                _materialize(src, out_dir / "ombsens" / key / fname, _DIAG_LINK_MODE)
            except OSError:
                if not config.IMAGES_IN_BUCKET:
                    return None
        elif not config.IMAGES_IN_BUCKET:
            return None
        return _url_for(fname)

    for kind in ("omb", "sens"):
        # Bucket mode gates on the persistent CSV (PNGs may have been pruned); else gate on the PNG.
        csv_path = base / f"{key}_{kind}.csv"
        if config.IMAGES_IN_BUCKET and not csv_path.exists():
            continue
        all_url = _stage(f"{key}_{kind}.png")          # all-time keeps the legacy <key>_<kind>.png name
        if all_url:
            out[kind] = all_url
        # Per-period panels: <key>_<kind>_<period>.png. In bucket mode the period list comes from the CSV
        # (PNGs may be gone); on a local site it comes from globbing the staged PNGs.
        per: dict = {}
        if config.IMAGES_IN_BUCKET and csv_path.exists():
            try:
                with open(csv_path, newline="", encoding="utf-8") as f:
                    import csv as _csv
                    for r in _csv.DictReader(f):
                        pk = (r.get("period") or "").strip()
                        if pk and pk != "all":
                            u = _stage(f"{key}_{kind}_{pk}.png")
                            if u:
                                per[pk] = u
            except OSError:
                pass
        else:
            prefix = f"{key}_{kind}_"
            for png in sorted(base.glob(f"{prefix}*.png")):
                pk = png.name[len(prefix):-4]          # strip the '<key>_<kind>_' prefix and '.png'
                if pk and pk != "all":
                    u = _stage(png.name)
                    if u:
                        per[pk] = u
        if per:
            out[f"{kind}_periods"] = per
    return out


def _load_ceda_links(path) -> dict:
    """{key: CEDA-L2 URL} from the JSON written by scripts/build_ceda_links.py; empty if absent."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as f:
            d = json.load(f)
        return {str(k): str(v) for k, v in d.items() if v}
    except (OSError, ValueError):
        return {}


def _stage_netcdfs(fullcal_dir, key, out_dir: Path) -> list:
    """Stage the per-station calibration NetCDF(s) into the site (nc/<key>/) the same way the
    diagnostics/OmB PNGs are staged (symlink by default; ALC_DIAG_LINK overrides), and return
    [{'year', 'fname', 'url'}, ...] newest year first. One NetCDF per year lives under
    <key>/<YYYY>/ALC_calibration_<WMO>_<IDENT><YYYY>.nc. Staging (not just linking a URL) means
    ops/publish.sh can sync the site's nc/ tree to the bucket."""
    out = []
    if not fullcal_dir:
        return out
    base = Path(fullcal_dir) / key
    if not base.exists():
        return out
    for nc in sorted(base.glob("*/ALC_calibration_*.nc")):
        year = nc.parent.name
        try:
            _materialize(nc, out_dir / "nc" / key / nc.name, _DIAG_LINK_MODE)
        except OSError:
            if not config.IMAGES_IN_BUCKET:
                continue
        out.append({"year": year, "fname": nc.name, "url": config.nc_url(key, nc.name)})
    out.sort(key=lambda r: r["year"], reverse=True)
    return out


def _load_status(fullcal_dir, key):
    """Per-stream decoded daily status (<key>_status.csv) for the availability bar; None if absent."""
    if not fullcal_dir:
        return None
    p = Path(fullcal_dir) / key / f"{key}_status.csv"
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p, dtype={"date": str})
    except Exception:
        return None
    if "date" not in df.columns or not len(df):
        return None
    return df


def _latest_status_by_key(fullcal_dir, keys) -> dict:
    """{key: latest daily quality class} from the last line of each ``<key>_status.csv``.

    Reads the header plus the tail of each file rather than parsing it, because this runs once for
    the whole network (436 files) purely to put a status dot on the neighbour links."""
    out = {}
    if not fullcal_dir:
        return out
    for k in keys:
        p = Path(fullcal_dir) / str(k) / f"{k}_status.csv"
        try:
            with p.open("rb") as fh:
                header = fh.readline().decode("utf-8", "replace").strip().split(",")
                if "quality" not in header:
                    continue
                idx = header.index("quality")
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(0, size - 4096))          # the tail always holds the last full row
                tail = fh.read().decode("utf-8", "replace").strip().splitlines()
        except OSError:
            continue
        for line in reversed(tail):
            parts = line.split(",")
            if len(parts) > idx and parts[0][:1].isdigit():   # skip a re-read header line
                out[str(k)] = parts[idx]
                break
    return out


def _station_index_records(st, series, fullcal_dir) -> list:
    """One compact record per station for the client-side station index.

    Feeds three things at once: the nav-bar search panel, the country/instrument filters on a
    station page, and the previous/next links (which show the neighbour's status and calibration
    constant). Emitted in the SAME order as the stations table, i.e. the order prev/next walks, so
    the client never re-sorts and the navigation order is identical to the old server-baked one.

    Keys are short because this file is fetched by every page: k=key, n=name, w=WIGOS id,
    c=country, t=instrument type, q=latest daily quality, m={method: {...}}."""
    latest_q = _latest_status_by_key(fullcal_dir, list(st["key"]))
    by_key: dict = {}
    for _, r in series.iterrows():
        by_key.setdefault(str(r["key"]), {})[str(r["method"])] = r
    out = []
    for _, r in st.iterrows():
        k = str(r["key"])
        itype = str(r.get("itype", "") or "")
        theo = config.theoretical_cl(itype)
        methods = {}
        for m, s in by_key.get(k, {}).items():
            cl = s.get("median_cl")
            cl = float(cl) if pd.notna(cl) else None
            flag = s.get("last_flag")
            methods[m] = {
                # raw C_L spans 11 orders of magnitude across types, so the percentage of the
                # type's nominal value is what makes the number readable in a one-line link
                "cl": float(f"{cl:.4g}") if cl else None,
                "pct": round(100.0 * cl / theo, 1) if (cl and theo) else None,
                "f": float(flag) if pd.notna(flag) else None,
                "d": str(s.get("last_date") or ""),
            }
        out.append({
            "k": k,
            "n": str(r.get("name", "") or ""),
            "w": k.rsplit("_", 1)[0] if "_" in k else k,
            "c": str(r.get("country", "") or ""),
            "t": itype,
            "q": latest_q.get(k, ""),
            "m": methods,
        })
    return out


def _load_hk(fullcal_dir, key):
    """Per-stream daily housekeeping (<key>_hk.csv) for the monitoring panel; None if absent/empty."""
    if not fullcal_dir:
        return None
    p = Path(fullcal_dir) / key / f"{key}_hk.csv"
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p, dtype={"date": str})
    except Exception:
        return None
    if "date" not in df.columns or not len(df):
        return None
    df["datetime"] = pd.to_datetime(df["date"], format="%Y%m%d", errors="coerce")
    hk_cols = [f for f, *_ in config.HK_PANEL if f in df.columns]
    for c in hk_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["datetime"]).sort_values("datetime")
    if not any(df[c].notna().any() for c in hk_cols):   # nothing to plot
        return None
    return df


# --- per-station page render (shared by the serial + parallel paths) ----------
def _render_one_station(key, ctx) -> str:
    """Render and write one station's HTML page from the shared context *ctx*. Independent of every
    other station (writes only stations/<key>.html and stages that key's own OmB/sens PNGs), so it is
    safe to run concurrently across stations."""
    cal, kal, series, st = ctx.cal, ctx.kal, ctx.series, ctx.st
    meta = st[st["key"] == key].iloc[0].to_dict()
    methods = [m for m in config.METHOD_ORDER
               if len(cal[(cal["key"] == key) & (cal["method"] == m)])]
    blocks = [_method_block(key, m, cal, kal, series, ctx.diag_by.get((key, m), []),
                            ctx.op_all, ctx.oldray_all) for m in methods]
    # Cloudnet target classification curtains: a per-day gallery like the calibration diagnostics, but
    # not tied to a calibration method (it has no cal rows), so it renders as its own standalone card.
    class_diags = ctx.diag_by.get((key, "classification"), [])
    overlay = None
    if len(methods) >= 2:
        by_method = {m: cal[(cal["key"] == key) & (cal["method"] == m)] for m in methods}
        overlay = charts.fig_to_div(charts.cl_overlay(by_method), "fig-overlay")
    i = ctx.nav_idx.get(key)
    prev_station = f"{ctx.all_keys[i - 1]}.html" if (i is not None and i > 0) else ""
    next_station = f"{ctx.all_keys[i + 1]}.html" if (i is not None and i < len(ctx.all_keys) - 1) else ""
    hk_df = _load_hk(ctx.fullcal_dir, key)
    monitoring = (charts.fig_to_div(charts.monitoring_timeseries(hk_df, meta.get("itype")), "fig-hk")
                  if hk_df is not None else None)
    # Cloudnet-style per-day availability/health bar + a date->status lookup for the diagnostic viewer.
    status_df = _load_status(ctx.fullcal_dir, key)
    availability, status_json = None, None
    if status_df is not None:
        avail_fig = charts.daily_availability_rows(status_df, cal[cal["key"] == key], methods)
        if avail_fig is not None:
            availability = charts.fig_to_div(avail_fig, "fig-avail")
            status_map = {str(r["date"]): {"q": str(r.get("quality", "")), "s": str(r.get("summary", ""))}
                          for _, r in status_df.iterrows()}
            status_json = json.dumps(status_map, ensure_ascii=False)
    ombsens = _stage_ombsens_pngs(ctx.fullcal_dir, key, ctx.out_dir)
    # Top-of-page links: download the calibration NetCDF(s) + the matching CEDA L2 archive page.
    nc_files = _stage_netcdfs(ctx.fullcal_dir, key, ctx.out_dir)
    ceda_url = getattr(ctx, "ceda_by", {}).get(key)
    cal_classes = [{"key": c, "label": config.CAL_CLASS_SHORT[c],
                    "title": config.CAL_CLASS_LABELS[c],
                    "color": config.CAL_CLASS_COLORS[c]} for c in config.CAL_CLASS_ORDER]
    html = ctx.tmpl.render(base="../", logo=ctx.logo, key=key, meta=meta, cal_classes=cal_classes,
                           blocks=blocks, overlay=overlay, search_json=ctx.search_json,
                           countries=getattr(ctx, "countries", []), types=getattr(ctx, "types", []),
                           prev_station=prev_station, next_station=next_station,
                           monitoring=monitoring, availability=availability, status_json=status_json,
                           ombsens=ombsens, nc_files=nc_files, ceda_url=ceda_url,
                           class_diags=class_diags,
                           periods=getattr(ctx, "periods", None),
                           periods_json=getattr(ctx, "periods_json", None))
    (ctx.out_dir / "stations" / f"{key}.html").write_text(html, encoding="utf-8")
    return key


# Per-worker context for the parallel render path. Each worker process reloads the (read-only) frames
# from the SQLite index ONCE via _render_worker_init -- cheaper and more portable (Windows spawn +
# Linux fork) than pickling the large frames to every task. The parent has already materialized the
# shared diagnostic PNGs, so workers build only the read-only diag index.
_WORKER_CTX = None


def _render_worker_init(db_path, out_dir, fullcal_dir, opcoeff_csv, oldray_dir, logo,
                        search_json, periods_json, ceda_json="", stations_v="",
                        countries_json="", types_json=""):
    global _WORKER_CTX
    if _WORKER_CTX is not None:
        return  # fork (Linux/CSCS): the worker inherited the parent's ctx -> no reload/re-index
    # spawn (Windows): rebuild the read-only context from the pickled paths
    cal, series, st, kal, diag = metrics.load_frames(Path(db_path))
    all_keys = list(st["key"])
    env = _env()
    # the parent computed this from the station-index content; a spawned worker must not emit a
    # different (or empty) cache-busting token than the serial path
    env.globals["stations_v"] = stations_v or ""
    _WORKER_CTX = SimpleNamespace(
        cal=cal, series=series, st=st, kal=kal,
        diag_by=_diag_index(diag, cal),
        op_all=_load_opcoeff(opcoeff_csv or None),
        oldray_all=_load_oldray(oldray_dir or None),
        tmpl=env.get_template("station.html"),
        out_dir=Path(out_dir), fullcal_dir=(fullcal_dir or None),
        all_keys=all_keys, nav_idx={k: i for i, k in enumerate(all_keys)},
        logo=(logo or None), search_json=search_json,
        ceda_by=(json.loads(ceda_json) if ceda_json else {}),
        countries=(json.loads(countries_json) if countries_json else []),
        types=(json.loads(types_json) if types_json else []),
        periods=(json.loads(periods_json) if periods_json else None),
        periods_json=(periods_json or None))


def _render_worker_task(key):
    return _render_one_station(key, _WORKER_CTX)


def _render_summary_page(env, logo, cal_w, st, fullcal_dir, op_all, search_json,
                         countries, types, periods, periods_json, current_key, out_path):
    """Render ONE summary page (index.html or index_<period>.html) from a date-windowed
    ``calibrations`` frame ``cal_w``. Every aggregate (KPIs, maps, success-by-type, flag distribution,
    C_L boxes, per-station IQR, watchlist, series table) is recomputed over the window. The per-station
    ``series`` aggregates are rebuilt in-memory from ``cal_w`` (not read from the prebuilt table) so the
    window is honoured. The OmB-bias / ICAO maps read the row matching ``current_key`` (period-keyed once
    the runner emits per-period summaries; full-archive fallback otherwise). Returns its summary dict."""
    series_w = index._series_aggregates(cal_w).merge(st[["key", "itype"]], on="key", how="left")
    summary = metrics.network_summary(cal_w, series_w, st)
    flags = metrics.flag_distribution(cal_w)

    # --- Reactive-filter payloads (country/type re-computed client-side; period is per-page) --------
    country_by = dict(zip(st["key"], st["country"])) if "country" in st.columns else {}
    type_by = dict(zip(st["key"], st["itype"])) if "itype" in st.columns else {}
    # (a) per-series counts -> the 4 KPIs recompute on filter change
    series_index = [
        {"key": str(r["key"]), "itype": str(r.get("itype", "") or ""),
         "country": str(country_by.get(r["key"], "") or ""), "method": str(r.get("method", "") or ""),
         "n_dates": int(r.get("n_dates") or 0), "n_success": int(r.get("n_success") or 0)}
        for _, r in series_w.iterrows()]
    series_json = json.dumps(series_index, ensure_ascii=False)
    # (b) per-instrument monthly activity -> the stacked "instruments over time" chart
    cal_m = cal_w[["key", "date"]].copy()
    cal_m["month"] = cal_m["date"].astype(str).str.slice(0, 6)
    activity = [
        {"key": str(k), "itype": str(type_by.get(k, "") or ""),
         "country": str(country_by.get(k, "") or ""),
         "months": sorted(set(g["month"].dropna().tolist()))}
        for k, g in cal_m.groupby("key", sort=False)]
    activity_json = json.dumps(activity, ensure_ascii=False)

    watch = metrics.watchlist(cal_w, st)
    keystats = _keystats(series_w, st)
    keystats = keystats.merge(_opcoeff_ratios(cal_w, op_all, st), on="key", how="left")
    keystats = keystats.merge(_ombsens_keystats(fullcal_dir, period=current_key), on="key", how="left")

    summary_figs = {
        "map_theo": charts.fig_to_div(charts.ratio_map(
            keystats, "op_pct_theo", "Median C_L — % of theoretical value",
            "% of theoretical", "fig-map-theo"), "fig-map-theo"),
        "map_op": charts.fig_to_div(charts.ratio_map(
            keystats, "op_pct_op", "Median C_L — % of operational constant (L2)",
            "% of operational", "fig-map-op"), "fig-map-op"),
        "map": charts.fig_to_div(charts.network_map(keystats), "fig-map"),
        "map_omb": charts.fig_to_div(charts.omb_bias_map(keystats), "fig-map-omb"),
        "map_icao": charts.fig_to_div(charts.icao_altitude_map(keystats), "fig-map-icao"),
        "instr": charts.fig_to_div(charts.instrument_count_over_time(activity), "fig-instr"),
        "success_type": charts.fig_to_div(charts.success_by_type_method(summary["by_type_method"]), "fig-stype"),
        "flag_dist_rayleigh": charts.fig_to_div(charts.flag_distribution_bar(flags, "rayleigh"), "fig-flags-r"),
        "flag_dist_cloud": charts.fig_to_div(charts.flag_distribution_bar(flags, "cloud"), "fig-flags-c"),
        "cl_type_abs": charts.fig_to_div(charts.value_by_type_method_box(series_w), "fig-cltype"),
        "cl_type_pct": charts.fig_to_div(charts.value_pct_theoretical_box(series_w), "fig-cltype-pct"),
    }
    # Per-station median C_L with IQR (Q1..Q3 of that station's successful daily values), one ranked
    # plot per instrument type. Pools both methods per station -- C_L is the same physical quantity.
    key_itype = dict(zip(st["key"], st["itype"]))
    okc = cal_w[(cal_w["success"] == 1) & (cal_w["cal_value"] > 0)].copy()
    okc["itype"] = okc["key"].map(key_itype)
    okc = okc.dropna(subset=["itype"])
    gb = okc.groupby(["itype", "key"])["cal_value"]
    sta_iqr = pd.DataFrame({"med": gb.median(), "q1": gb.quantile(0.25),
                            "q3": gb.quantile(0.75), "n": gb.size()}).reset_index()
    sta_iqr["country"] = sta_iqr["key"].map(
        dict(zip(st["key"], st["country"])) if "country" in st.columns else {}).fillna("")
    cl_iqr_figs = [(t, charts.fig_to_div(charts.cl_median_iqr_by_station(sta_iqr[sta_iqr["itype"] == t], t),
                                         f"fig-cliqr-{t}"))
                   for t in config.TYPE_ORDER if (sta_iqr["itype"] == t).any()]

    html = env.get_template("summary.html").render(
        base="", logo=logo, summary=summary, figs=summary_figs, cl_iqr=cl_iqr_figs,
        watch=watch.to_dict("records"), rows=_series_table_rows(cal_w, series_w, st),
        countries=countries, types=types, search_json=search_json,
        series_json=series_json, activity_json=activity_json,
        periods=periods, periods_json=periods_json, current_period=current_key)
    out_path.write_text(html, encoding="utf-8")
    return summary


def build_site(db_path: Path, out_dir: Path, limit_pages: int | None = None,
               flagex_dir=None, opcoeff_csv=None, only_keys=None, oldray_dir=None,
               fullcal_dir=None, workers: int | None = None, ceda_links=None) -> dict:
    out_dir = Path(out_dir)
    (out_dir / "stations").mkdir(parents=True, exist_ok=True)
    logo = _write_assets(out_dir)
    env = _env()

    cal, series, st, kal, diag = metrics.load_frames(db_path)
    diag_by = _copy_diagnostics(diag, cal, out_dir)

    # Operational calibration constant from the L2 files (optional): two ratio maps + per-station
    # black line on the time series. Loaded once; reused by every per-period summary page below.
    op_all = _load_opcoeff(opcoeff_csv)
    oldray_all = _load_oldray(oldray_dir)
    ceda_by = _load_ceda_links(ceda_links)   # {key: CEDA L2 URL} for the per-station link
    # NB: the summary aggregates + figures (KPIs, maps, boxes, IQR, watchlist, series table) are now
    # built per time-period in _render_summary_page(), so the page set can re-window cheaply.

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

    # Station index written ONCE to data/stations.json and fetched by every page, instead of being
    # inlined into all 436 of them. Its own content hash is the cache-busting token, so the URL only
    # changes when the data does. The inline #search-index blob stays for now as the file:// fallback
    # (a build opened by double-click cannot fetch).
    stations_json = json.dumps(_station_index_records(st, series, fullcal_dir),
                               ensure_ascii=False, separators=(",", ":"))
    (out_dir / "data").mkdir(parents=True, exist_ok=True)
    _write_if_changed(out_dir / "data" / "stations.json", stations_json)
    env.globals["stations_v"] = hashlib.md5(stations_json.encode("utf-8")).hexdigest()[:8]
    # --- Time-period set (auto-derived years; active vs frozen) ----------------
    # All-time + each calendar year (first..current) + rolling last-N-day windows. The current year,
    # the rolling windows and all-time are rebuilt every run; a past complete year is built ONCE and
    # then skipped (its data can no longer change once the backfill window has passed). See periods.py.
    as_of = str(cal["date"].max()) if len(cal) else None
    date_min = str(cal["date"].min()) if len(cal) else None
    backfill = int(os.environ.get("ALC_BACKFILL_DAYS", "5") or "5")
    period_objs = periods_mod.build_periods(date_min, as_of, backfill_days=backfill)
    period_list = periods_mod.periods_to_json(period_objs)
    periods_json = json.dumps(period_list, ensure_ascii=False)

    summary = None
    n_built = 0
    for p in period_objs:
        out_name = "index.html" if p.key == "all" else f"index_{p.key}.html"
        out_path = out_dir / out_name
        if not p.active and out_path.exists():
            continue  # frozen past year, already built -> never recomputed in real time
        cal_w = periods_mod.filter_window(cal, p.start, p.end)
        s = _render_summary_page(env, logo, cal_w, st, fullcal_dir, op_all, search_json,
                                 countries, types, period_list, periods_json, p.key, out_path)
        n_built += 1
        if p.key == "all":
            summary = s
    print(f"  summary pages: {n_built} built / {len(period_objs)} periods "
          f"({sum(1 for p in period_objs if not p.active)} frozen)", flush=True)
    if summary is None:   # 'all' is always active, but guard against an empty period set
        sa = index._series_aggregates(cal).merge(st[["key", "itype"]], on="key", how="left")
        summary = metrics.network_summary(cal, sa, st)

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
    ctx = SimpleNamespace(
        cal=cal, kal=kal, series=series, st=st, diag_by=diag_by, op_all=op_all,
        oldray_all=oldray_all, tmpl=station_tmpl, out_dir=out_dir, fullcal_dir=fullcal_dir,
        all_keys=all_keys, nav_idx=nav_idx, logo=logo, search_json=search_json,
        ceda_by=ceda_by, countries=countries, types=types,
        periods=period_list, periods_json=periods_json)

    n_workers = int(workers) if workers else 1
    if n_workers > 1 and len(keys) > 1:
        # Per-station pages are independent -> fan out across processes. Each worker reloads the
        # read-only frames once (_render_worker_init); the diagnostic/asset materialization above
        # already ran in the parent, so workers only write their own stations/<key>.html + stage
        # their own OmB/sens PNGs (distinct paths -> no races).
        print(f"  rendering {len(keys)} station pages on {n_workers} workers ...", flush=True)
        initargs = (str(db_path), str(out_dir), str(fullcal_dir) if fullcal_dir else "",
                    str(opcoeff_csv) if opcoeff_csv else "", str(oldray_dir) if oldray_dir else "",
                    logo or "", search_json, periods_json,
                    json.dumps(ceda_by, ensure_ascii=False),
                    env.globals.get("stations_v", ""),
                    json.dumps(countries, ensure_ascii=False),
                    json.dumps(types, ensure_ascii=False))
        # Expose the parent's ctx so fork()ed workers (Linux/CSCS) inherit it for free; spawn()ed
        # workers (Windows) ignore this and rebuild from initargs in the initializer.
        global _WORKER_CTX
        _WORKER_CTX = ctx
        try:
            with ProcessPoolExecutor(max_workers=n_workers, initializer=_render_worker_init,
                                     initargs=initargs) as ex:
                for _ in ex.map(_render_worker_task, keys, chunksize=4):
                    pass
        finally:
            _WORKER_CTX = None
    else:
        for key in keys:
            _render_one_station(key, ctx)

    return dict(out_dir=str(out_dir), n_pages=len(keys), n_series=int(len(series)),
                as_of=summary["as_of"])
