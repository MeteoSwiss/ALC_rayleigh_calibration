#!/usr/bin/env python3
"""Build the interactive station dashboard: one page per stream, every panel, lazy daily payloads.

This is the full version of what `mockup_station_page.py` prototyped. Differences that matter:

* **the daily panel is lazy.** Six months x two methods is ~360 payloads of ~340 KB; embedding them
  would be a 120 MB page. Instead each (day, method) is written to `data/<key>/<date>_<method>.json`
  and fetched when that day is selected, while a tiny per-day INDEX (flag, constant, has_fig) is
  embedded so the calendar, the flag chips and the day arrows are complete from the first paint.
* **station controls live in the navbar** and drive the arrows: filtering by country or instrument
  changes what up/down step through.
* **the keyboard contract is the production one** (monitoring/static/diag.js:257-271), because
  operators already have it in their fingers:
      left / right         previous / next CALIBRATED day
      Ctrl+left / right    previous / next day that has any diagnostic
      up / down            previous / next station (within the navbar filter)
      0                    open the QC flag dialog
      1 / 2 / 3            quick flag: aerosol / cloud / low signal
* **OmB, sensitivity and Cloudnet classification come from the live bucket.** They are produced by
  the operational run, not by this script, so the page references
  `<bucket>/ombsens/<key>/<key>_omb[_<period>].png` and `<bucket>/diag/<key>/classification_<d>.png`
  rather than re-deriving them. A panel whose product does not exist says so instead of vanishing.

  python scripts/build_station_dashboard.py --key 0-20000-0-06610_A --key 0-20000-0-06610_B \
      --key 0-20000-0-06610_C --db <calib_index.sqlite> --cal-dir <tree> --status-dir <tree> \
      --start 20260215 --end 20260813 --l1-root D:/E-PROFILE_L1_2026 --cams D:/CAMS_daily \
      --out dashboard_v3 --payloads
"""
from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.offline as pyo

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from monitoring import charts, config                                     # noqa: E402
from monitoring.render import _load_hk, _load_status                      # noqa: E402
from scripts.mockup_station_card import probe_cloud_cover, _month_with_data   # noqa: E402
from scripts.mockup_station_page import _cl_series, _station_index, TILE_OK   # noqa: E402
import scripts.mockup_daily_panel as PANEL                                # noqa: E402

BUCKET = "https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard"
#: Period keys the operational OmB/sensitivity renderer emits (seen on the live station pages).
PERIODS = [("all", "all"), ("last30", "last 30 d"), ("last90", "last 90 d"),
           ("last180", "last 180 d"), ("last365", "last 365 d"),
           ("y2025", "2025"), ("y2026", "2026")]
#: A CHM15k SATURATES in liquid cloud, so the cloud retrieval is not meaningful on it and is never
#: run -- not even as a control. Everything else gets both methods.
METHODS_BY_TYPE = {"CHM15k": ["rayleigh"]}


def _methods_for(itype: str) -> list:
    return METHODS_BY_TYPE.get(str(itype), ["rayleigh", "cloud"])


def _meta_from_l1(l1_root, key: str) -> dict | None:
    """Station metadata straight out of the stream's own L1 file.

    This is the AUTHORITATIVE source and the one the pipeline itself reads: site_location,
    institution, instrument_type and the station coordinates are all attributes/variables of the
    file. Sourcing them from a dashboard SQLite instead was a mistake -- that database only contains
    whatever streams the last dashboard build happened to cover, so a perfectly normal station (the
    Payerne CL31) can be absent from it while its L1 sits right there.
    """
    if not l1_root:
        return None
    wmo, _, ident = key.rpartition("_")
    files = sorted(Path(l1_root).glob(f"{wmo}/*/*/L1_{wmo}_{ident}*.nc"))
    if not files:
        return None
    try:
        import netCDF4
        with netCDF4.Dataset(files[-1]) as ds:
            site = str(getattr(ds, "site_location", "") or "")
            name, _, country = site.partition(",")
            g = lambda v: (float(np.ravel(ds.variables[v][:])[0])            # noqa: E731
                           if v in ds.variables else None)
            return {"key": key, "wmo": str(getattr(ds, "wigos_station_id", wmo) or wmo),
                    "ident": ident,
                    "itype": str(getattr(ds, "instrument_type", "") or ""),
                    "name": name.strip().title() or wmo,
                    "country": country.strip().upper(),
                    "institution": str(getattr(ds, "institution", "") or "").strip(),
                    "lat": g("station_latitude"), "lon": g("station_longitude"),
                    "alt": g("station_altitude"),
                    "constants": {}, "status": None, "status_date": None,
                    "src": f"L1 {files[-1].name}"}
    except Exception as exc:                                                 # noqa: BLE001
        print(f"  {key}: L1 metadata unreadable ({type(exc).__name__}: {exc})", flush=True)
        return None


# ============================================================================ frames from the CSVs
def _cal_frame(cal_dir: Path, key: str) -> pd.DataFrame:
    """`<key>_cal.csv` -> the frame monitoring.charts expects (it normally comes from SQLite)."""
    p = Path(cal_dir) / key / f"{key}_cal.csv"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, dtype={"date": str})
    d["datetime"] = pd.to_datetime(d["date"], format="%Y%m%d", errors="coerce")
    d["cal_value"] = pd.to_numeric(d["cal_value"], errors="coerce")
    d["flag"] = pd.to_numeric(d["flag"], errors="coerce")
    # `success` is the charts' own convention; keep it on the constant, not the flag, for the same
    # reason the panel does: flag 0.5 is a partial success that DOES yield a constant.
    d["success"] = (d["cal_value"] > 0).astype(int)
    for c in ("uncertainty", "n_profiles", "bottom_height", "top_height"):
        d[c] = pd.to_numeric(d.get(c), errors="coerce")
    d["key"] = key
    return d


def _kal_frame(cal_dir: Path, key: str) -> pd.DataFrame:
    p = Path(cal_dir) / key / f"{key}_kalman.csv"
    if not p.exists():
        return pd.DataFrame(columns=["method", "date", "kalman", "kalman_std", "datetime"])
    d = pd.read_csv(p, dtype={"date": str})
    d["datetime"] = pd.to_datetime(d["date"], format="%Y%m%d", errors="coerce")
    return d


# ================================================================================== daily payloads
def _payload_task(job: tuple) -> tuple:
    """One (day, method) in a worker process: calibrate, write the JSON, report.

    Top-level and pickle-friendly because Windows spawns rather than forks. Each job writes its OWN
    file, so the workers share nothing and need no lock; the parent only reads the directory
    afterwards. Exceptions are returned, never raised, so one bad night cannot abort the batch.
    """
    key, ds, method, ddir = job
    try:
        d = datetime.strptime(ds, "%Y%m%d")
        payload, _ = PANEL.build_payload(key, [d], [method])
        p = (payload["days"].get(ds) or {}).get(method)
        if p is None or p.get("kind") == "error":
            return ds, method, False, (p or {}).get("message", "no diagnostic")
        Path(ddir, f"{ds}_{method}.json").write_text(json.dumps(p, separators=(",", ":")),
                                                     encoding="utf-8")
        return ds, method, True, p.get("kind", "")
    except Exception as exc:                                                 # noqa: BLE001
        return ds, method, False, f"{type(exc).__name__}: {exc}"


def _write_payloads(key: str, itype: str, days: list, out: Path, force: bool,
                    workers: int) -> dict:
    """Run every (day, method) and write one JSON each. Returns the per-day index.

    The index is what the page embeds: enough to colour a calendar and step the arrows, and small
    enough that six months of it is a few kB.

    Runs in a PROCESS pool: each night is an independent read of its own L1 file, so this is
    embarrassingly parallel and there is no reason to leave 31 of 32 cores idle. Processes rather
    than threads because the work is NumPy-heavy Python, and the per-process cost of importing the
    calibration stack is paid once per worker, not once per night.
    """
    methods = _methods_for(itype)
    ddir = out / "data" / key
    ddir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    jobs = [(key, d.strftime("%Y%m%d"), m, str(ddir)) for d in days for m in methods
            if force or not (ddir / f"{d.strftime('%Y%m%d')}_{m}.json").exists()]
    print(f"  {key}: {len(days)} days x {len(methods)} method(s) -> {len(jobs)} to compute "
          f"on {workers} workers", flush=True)
    if jobs:
        done = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for ds, m, ok, note in ex.map(_payload_task, jobs, chunksize=1):
                done += 1
                if done % 25 == 0 or done == len(jobs):
                    el = time.perf_counter() - t0
                    print(f"    {done}/{len(jobs)}  {el/60:.1f} min  "
                          f"({el/max(done,1):.2f} s/night wall, {len(jobs)-done} left)", flush=True)
    # Rebuild the index from what is on disk, so a resumed run indexes earlier days too.
    index = {}
    for d in days:
        ds = d.strftime("%Y%m%d")
        for m in methods:
            f = ddir / f"{ds}_{m}.json"
            if not f.exists():
                continue
            try:
                p = json.loads(f.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            index.setdefault(ds, {})[m] = {
                "flag": p.get("flag"), "constant": p.get("constant"),
                "message": p.get("message"),
                "flag_label": p.get("flag_label"),
                "has_fig": bool((p.get("curtain") or {}).get("b64")),
            }
    return index


# ============================================================================================ page
def _fig(fig, div_id):
    return charts.fig_to_div(fig, div_id) if fig is not None else \
        '<p class="muted">not available for this stream</p>'


def _summary_rows(cal: pd.DataFrame, method: str, index: dict) -> str:
    """The calibration list. Every date is a link that loads that day's panel."""
    g = cal[cal["method"] == method].sort_values("date", ascending=False)
    out = []
    for _, r in g.iterrows():
        ds = str(r["date"])
        has = ds in index and method in index[ds]
        val = "—" if not (r["cal_value"] > 0) else f"{r['cal_value']:.4g}"
        unc = "" if not np.isfinite(r.get("uncertainty", np.nan)) else f" ± {r['uncertainty']:.3g}"
        extra = (f"{r['n_profiles']:.0f} prof" if method == "cloud" and
                 np.isfinite(r.get("n_profiles", np.nan))
                 else (f"{r['bottom_height']:.0f}–{r['top_height']:.0f} m"
                       if np.isfinite(r.get("bottom_height", np.nan)) else "—"))
        cell = (f'<a class="daylink" data-date="{ds}" data-method="{method}" href="#daily">{ds}</a>'
                if has else ds)
        out.append(f'<tr><td>{cell}</td>'
                   f'<td class="num">{r["flag"]:g}</td>'
                   f'<td>{config.flag_label(r["flag"], method)}</td>'
                   f'<td class="num">{val}{unc}</td><td class="num">{extra}</td></tr>')
    return "\n".join(out)


def build_page(key: str, rec: dict, args, index: dict, stations: list) -> str:
    itype = rec.get("itype") or ""
    methods = _methods_for(itype)
    cal = _cal_frame(args.cal_dir, key)
    kal = _kal_frame(args.cal_dir, key)
    status_df = _load_status(Path(args.status_dir), key)
    hk_df = _load_hk(Path(args.status_dir), key)

    blocks = []
    # ---- C_L over time + headline tiles
    cl_fig, tiles, cl_sub = _cl_series(key, args.cal_dir)
    tiles_html = "".join(
        f'<div class="tile"><div class="tl">{t["label"]}</div>'
        f'<div class="tv" style="color:{t["color"]}">{t["value"]}</div>'
        f'<div class="tn">{t["note"]}</div></div>' for t in tiles)
    blocks.append(f'<section class="sect" id="constant"><h2>Calibration coefficient C_L over time</h2>'
                  f'<p class="subtle">{cl_sub}</p><div class="card">{_fig(cl_fig, "fig-cl")}</div>'
                  f'<div class="tiles">{tiles_html}</div></section>')

    # ---- availability card
    avail_html = '<p class="muted">no status file for this stream</p>'
    if status_df is not None and len(cal):
        sdf = status_df
        if args.l1_root and args.l1_month:
            month = _month_with_data(Path(args.l1_root), key, args.l1_month) or args.l1_month
            cc = probe_cloud_cover(Path(args.l1_root), key, month)
            if len(cc):
                sdf = sdf.merge(cc, on="date", how="left")
        have = [m for m in config.METHOD_ORDER if len(cal[cal["method"] == m])]
        avail_html = _fig(charts.daily_availability_rows(sdf, cal, have), "fig-avail")
    legend = ('<div class="avail-legend"><b>Status</b>'
              '<span class="ak ak-pass"></span> Pass<span class="ak ak-warning"></span> Warning'
              '<span class="ak ak-error"></span> Error<span class="ak ak-nodata"></span> No data'
              '<span class="legsep"></span><b>Cloud cover</b><span class="ak ak-cc0"></span> 0'
              '<span class="ak ak-cc4"></span> 4<span class="ak ak-cc8"></span> 8 octas'
              '<span class="legsep"></span><b>Calibration</b>'
              + "".join(f'<span class="ak" style="background:{config.CAL_CLASS_COLORS[c]}" '
                        f'title="{config.CAL_CLASS_LABELS[c]}"></span> {config.CAL_CLASS_SHORT[c]}'
                        for c in config.CAL_CLASS_ORDER) + '</div>')
    blocks.append(f'<section class="sect" id="availability">'
                  f'<h2>Data availability, cloud cover &amp; calibration</h2>{legend}'
                  f'<div class="card">{avail_html}</div></section>')

    # ---- the interactive daily panel (its body comes from the panel module)
    blocks.append(f'<section class="sect" id="daily"><h2>Daily calibration</h2>'
                  f'{PANEL.PANEL_BODY}</section>')

    # ---- per-method: series, monthly outcomes, calibration window, summary table
    for m in methods:
        g = cal[cal["method"] == m]
        if not len(g):
            continue
        km = kal[kal["method"] == m] if len(kal) else pd.DataFrame()
        rows = _summary_rows(cal, m, index)
        blocks.append(
            f'<section class="sect" id="m-{m}">'
            f'<h2><span class="mtag mtag-{m}">{config.method_label(m)}</span></h2>'
            f'<div class="card">{_fig(charts.series_timeseries(g, km, m), f"fig-ts-{m}")}</div>'
            f'<div class="grid2">'
            f'<div class="card">{_fig(charts.monthly_flag_bars(g, m), f"fig-mo-{m}")}</div>'
            f'<div class="card">{_fig(charts.aux_timeseries(g, m), f"fig-aux-{m}")}</div></div>'
            f'<details class="card tablecard"><summary>All {config.method_label(m)} '
            f'calibrations ({len(g)}) — click a date to open its daily panel</summary>'
            f'<table class="caltab"><tr><th>date</th><th class="num">flag</th><th>outcome</th>'
            f'<th class="num">value</th><th class="num">window / profiles</th></tr>{rows}</table>'
            f'</details></section>')

    # ---- bucket-served products: classification, OmB, sensitivity
    blocks.append(
        '<section class="sect" id="classification"><h2>Cloudnet target classification</h2>'
        '<p class="subtle">Served from the operational bucket; follows the day selected above.</p>'
        f'<div class="card"><img id="img-class" class="bigimg" alt="classification" '
        f'data-base="{args.bucket}/diag/{key}/classification_">'
        '<p class="muted" id="class-miss" hidden>no classification image for this day</p></div>'
        '</section>')
    for kind, title, note in (
            ("omb", "Observation − Background (CAMS)",
             "calibrated backscatter against the CAMS model over the selected window"),
            ("sens", "Instrument sensitivity / detection limit",
             "ICAO detection altitude and night-time noise floor")):
        srcs = " ".join(f'data-src-{pk}="{args.bucket}/ombsens/{key}/{key}_{kind}'
                        f'{"" if pk == "all" else "_" + pk}.png"' for pk, _ in PERIODS)
        blocks.append(
            f'<section class="sect" id="{kind}"><h2>{title} <span class="muted">— {note}</span></h2>'
            f'<div class="card"><img id="img-{kind}" class="bigimg" alt="{kind}" {srcs}>'
            f'<p class="muted" id="{kind}-miss" hidden>this product has not been produced for '
            f'{key}</p></div></section>')

    # ---- instrument monitoring
    blocks.append(
        '<section class="sect" id="monitoring"><h2>Instrument monitoring '
        '<span class="muted">— daily averages: laser power, window transmission, temperatures'
        '</span></h2><div class="card">'
        + (_fig(charts.monitoring_timeseries(hk_df), "fig-hk") if hk_df is not None
           else '<p class="muted">no housekeeping file for this stream</p>')
        + '</div></section>')

    # ---- the payload the panel boots from: index for every day, full payload for none yet
    dates = sorted(index.keys())
    boot = {"station": {"key": key, "wmo": rec.get("wmo"), "ident": rec.get("ident"),
                        "type": itype, "name": rec.get("name") or ""},
            "methods": methods, "dates": dates, "days": {}, "index": index}
    meta = {"n_days": len(dates), "gz_kb": 0, "per_unit_kb": 0, "curated": False, "lazy": True}

    css = (REPO / "monitoring/static/style.css").read_text(encoding="utf-8")
    panel_js = PANEL.PANEL_JS.replace("__META__", json.dumps(meta))
    return (PAGE.replace("__STYLECSS__", css)
                .replace("__PANELCSS__", PANEL.PANEL_CSS)
                .replace("__TITLE__", f"{rec.get('name') or key} — ALC calibration")
                .replace("__HEADER__", _header(rec))
                .replace("__BLOCKS__", "\n".join(blocks))
                .replace("__STATIONS__", json.dumps(stations, separators=(",", ":")))
                .replace("__CURKEY__", json.dumps(key))
                .replace("__PERIODS__", json.dumps([{"k": k, "l": l} for k, l in PERIODS]))
                .replace("__PANELJS__", panel_js)
                .replace("__PAYLOAD__", json.dumps(boot, separators=(",", ":")))
                .replace("__PLOTLY__", pyo.get_plotlyjs()))


def _header(rec: dict) -> str:
    def f(v, spec, suffix=""):
        return "—" if v is None else format(v, spec) + suffix
    bits = [rec.get("country") or "—", rec.get("institution") or "—", rec.get("wmo", "—"),
            f"unit {rec.get('ident', '—')}",
            f"{f(rec.get('lat'), '.3f')}, {f(rec.get('lon'), '.3f')}",
            f(rec.get("alt"), ".0f", " m")]
    return (f'<h1>{rec.get("name") or rec.get("wmo")} '
            f'<span class="itype">{rec.get("itype") or ""}</span></h1>'
            f'<p class="metainline">{" · ".join(bits)}</p>')


PAGE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script>__PLOTLY__</script>
<style>
__STYLECSS__
__PANELCSS__
/* ---------- page shell + navbar-integrated station controls ---------- */
body { padding:0 0 40px; background:#f7f9fb; }
.nav { position:sticky; top:0; z-index:50; display:flex; align-items:center; gap:9px;
       flex-wrap:wrap; padding:8px 16px; background:#fff; border-bottom:1px solid #dbe3ea;
       box-shadow:0 1px 4px rgba(20,40,60,.06); }
.nav select, .nav input { font:inherit; font-size:13px; padding:4px 8px; border:1px solid #c3ceda;
       border-radius:7px; background:#fff; }
.nav .sp { flex:1; }
.nav .cnt { font-size:11.5px; color:#66707a; }
.kbd { font:11px ui-monospace,Consolas,monospace; background:#eef4fb; border:1px solid #d6e4f2;
       border-radius:5px; padding:1px 5px; color:#0b3d61; }
.wrap { padding:0 18px; }
h1 { font-size:22px; margin:14px 0 2px; }
h1 .itype { font-size:14px; font-weight:500; color:#66707a; margin-left:6px; }
.metainline { color:#66707a; font-size:13px; margin:0 0 10px; }
h2 { font-size:13px; margin:0 0 8px; text-transform:uppercase; letter-spacing:.04em; color:#66707a; }
.sect { border-top:1px solid #e3e9ef; margin-top:20px; padding-top:12px; }
.subtle { color:#66707a; font-size:12.5px; margin:0 0 8px; }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:10px;
         margin-top:11px; }
.tile { border:1px solid #dbe3ea; border-radius:9px; padding:9px 12px; background:#fff; }
.tile .tl { font:11px/1.3 ui-monospace,Consolas,monospace; color:#66707a; letter-spacing:.03em; }
.tile .tv { font-size:23px; font-weight:600; margin:3px 0 1px; font-variant-numeric:tabular-nums; }
.tile .tn { font-size:11.5px; color:#66707a; }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:12px; margin-top:12px; }
.tablecard { padding:10px 12px; margin-top:12px; }
.tablecard summary { cursor:pointer; font-size:12.5px; color:#2a5a82; }
.caltab { border-collapse:collapse; width:100%; font-size:12px; margin-top:8px; }
.caltab th { text-align:left; font-size:10.5px; text-transform:uppercase; color:#66707a;
             border-bottom:1px solid #dbe3ea; padding:3px 7px; position:sticky; top:0;
             background:#fff; }
.caltab td { padding:2px 7px; border-bottom:1px solid #f1f5f8; }
.caltab .num { text-align:right; font-variant-numeric:tabular-nums; }
.tablecard[open] { max-height:460px; overflow:auto; }
.daylink { color:#1f6feb; text-decoration:none; border-bottom:1px dotted #9dc2f0; }
.daylink:hover { background:#eef4fb; }
.bigimg { width:100%; height:auto; border-radius:7px; display:block; }
.mtag { font-size:11px; padding:1px 7px; border-radius:9px; color:#fff; }
.mtag-rayleigh { background:#1f77b4; } .mtag-cloud { background:#2ca02c; }
.qcbar { margin-top:9px; font-size:12px; color:#66707a; }
.qcchip { display:inline-block; border-radius:11px; padding:2px 9px; font-size:11.5px;
          border:1px solid #dbe3ea; background:#fff; margin-right:5px; }
@media (max-width:1100px) { .grid2 { grid-template-columns:1fr; } }
</style></head><body>

<div class="nav">
  <b style="font-size:13px">ALC calibration</b>
  <select id="f-country" title="Country filter — also limits the up/down arrows"></select>
  <select id="f-type" title="Instrument filter — also limits the up/down arrows"></select>
  <button id="st-prev" title="Previous station (Up arrow)">↑</button>
  <select id="st-sel" title="Station"></select>
  <button id="st-next" title="Next station (Down arrow)">↓</button>
  <span class="cnt" id="st-pos"></span>
  <span class="sp"></span>
  <select id="period-sel" title="Time window for the OmB and sensitivity panels"></select>
  <span class="cnt"><span class="kbd">←→</span> day
    <span class="kbd">Ctrl ←→</span> any day
    <span class="kbd">↑↓</span> station
    <span class="kbd">0</span> flag
    <span class="kbd">1 2 3</span> quick flag</span>
</div>

<div class="wrap">
__HEADER__
__BLOCKS__
</div>

<script id="stations" type="application/json">__STATIONS__</script>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>__PANELJS__</script>
<script>
// ===================================================================== station navigation (navbar)
const ST = JSON.parse(document.getElementById('stations').textContent);
const KEY = __CURKEY__;
const PERIODS = __PERIODS__;
const selC = document.getElementById('f-country'), selT = document.getElementById('f-type');
const selS = document.getElementById('st-sel'), selP = document.getElementById('period-sel');
const LS = 'alc-';                                   // sessionStorage prefix, as the live site uses

function opts(sel, vals, label) {
  sel.innerHTML = `<option value="">${label}</option>` +
    vals.map(v => `<option value="${v}">${v}</option>`).join('');
}
opts(selC, [...new Set(ST.map(s => s.country).filter(Boolean))].sort(), 'All countries');
opts(selT, [...new Set(ST.map(s => s.itype).filter(Boolean))].sort(), 'All instruments');
selP.innerHTML = PERIODS.map(p => `<option value="${p.k}">${p.l}</option>`).join('');
selC.value = sessionStorage.getItem(LS + 'country') || '';
selT.value = sessionStorage.getItem(LS + 'type') || '';
selP.value = sessionStorage.getItem(LS + 'period') || 'all';

// The arrows walk the FILTERED list. Arrows that ignored the dropdowns would jump straight out of
// the selection the operator just made, which is the whole reason the two are in one bar.
const filtered = () => ST.filter(s =>
  (!selC.value || s.country === selC.value) && (!selT.value || s.itype === selT.value));
function syncNav() {
  const l = filtered(), i = l.findIndex(s => s.key === KEY);
  selS.innerHTML = l.map(s => `<option value="${s.key}"${s.key === KEY ? ' selected' : ''}>` +
      `${s.name || s.wmo} · ${s.itype} · ${s.ident}</option>`).join('');
  if (i < 0) selS.innerHTML = `<option selected>${KEY} (outside filter)</option>` + selS.innerHTML;
  document.getElementById('st-prev').disabled = i <= 0;
  document.getElementById('st-next').disabled = i < 0 || i >= l.length - 1;
  document.getElementById('st-pos').textContent = i >= 0
    ? `${i + 1} / ${l.length}` : `outside filter (${l.length})`;
  sessionStorage.setItem(LS + 'country', selC.value);
  sessionStorage.setItem(LS + 'type', selT.value);
}
const go = k => { if (k && k !== KEY) location.href = `station_${k}.html`; };
function stepStation(step) {
  const l = filtered(), i = l.findIndex(s => s.key === KEY);
  if (i < 0) { if (l.length) go(l[0].key); return true; }
  const j = i + step;
  if (j < 0 || j >= l.length) return false;
  go(l[j].key); return true;
}
[selC, selT].forEach(s => s.addEventListener('change', syncNav));
selS.addEventListener('change', () => go(selS.value));
document.getElementById('st-prev').addEventListener('click', () => stepStation(-1));
document.getElementById('st-next').addEventListener('click', () => stepStation(1));
syncNav();

// ===================================================================== period -> OmB / sensitivity
function applyPeriod() {
  const p = selP.value || 'all';
  sessionStorage.setItem(LS + 'period', p);
  ['omb', 'sens'].forEach(kind => {
    const img = document.getElementById('img-' + kind);
    if (!img) return;
    const src = img.getAttribute('data-src-' + p) || img.getAttribute('data-src-all');
    const miss = document.getElementById(kind + '-miss');
    img.onerror = () => { img.hidden = true; if (miss) miss.hidden = false; };
    img.onload = () => { img.hidden = false; if (miss) miss.hidden = true; };
    if (src) img.src = src;
  });
}
selP.addEventListener('change', applyPeriod);
applyPeriod();

// ===================================================================== lazy daily payloads
// Only the day being looked at is fetched; the embedded index already carries every day's outcome,
// so nothing on the page has to wait for this.
const CACHE = {};
async function ensureDay(ds, m) {
  const ck = ds + '_' + m;
  if (CACHE[ck] !== undefined) return CACHE[ck];
  try {
    const r = await fetch(`data/${KEY}/${ds}_${m}.json`, { cache: 'no-cache' });
    CACHE[ck] = r.ok ? await r.json() : null;
  } catch (e) { CACHE[ck] = null; }
  if (CACHE[ck]) { D.days[ds] = D.days[ds] || {}; D.days[ds][m] = CACHE[ck]; }
  return CACHE[ck];
}
async function goDay(ds, m) {
  if (!ds) return;
  curDate = ds;
  if (m) curMethod = m;
  await Promise.all((D.methods || []).map(mm => ensureDay(ds, mm)));
  render();
  syncDayLinked();
}
// The panel calls this when it wants a day it does not have yet.
window.__onMissingDay = (ds, m) => { ensureDay(ds, m).then(p => { if (p) render(); }); };

// Cloudnet classification follows the selected day -- it is a per-day product like the panel.
function syncDayLinked() {
  const img = document.getElementById('img-class'), miss = document.getElementById('class-miss');
  if (!img) return;
  img.onerror = () => { img.hidden = true; if (miss) miss.hidden = false; };
  img.onload = () => { img.hidden = false; if (miss) miss.hidden = true; };
  img.src = img.getAttribute('data-base') + curDate + '.png';
  renderQC();
}

// ===================================================================== date links
document.addEventListener('click', e => {
  const a = e.target.closest('.daylink');
  if (!a) return;
  e.preventDefault();
  goDay(a.dataset.date, a.dataset.method);
  document.getElementById('daily').scrollIntoView({ behavior: 'smooth', block: 'start' });
});

// ===================================================================== QC flagging (0 / 1 / 2 / 3)
// Same four keys and the same three quick reasons as the operational page
// (monitoring/static/diag.js:268-271). Stored locally here; in production this posts to the QC store.
const QC_REASONS = ['aerosol contamination', 'cloud contamination', 'low signal (condensation?)'];
const qcKey = () => `alc-qc-${KEY}-${curDate}`;
function setQC(text) {
  if (text === null) localStorage.removeItem(qcKey()); else localStorage.setItem(qcKey(), text);
  renderQC();
}
function renderQC() {
  let bar = document.getElementById('qcbar');
  if (!bar) {
    const host = document.getElementById('daily');
    if (!host) return;
    bar = document.createElement('div');
    bar.id = 'qcbar'; bar.className = 'qcbar';
    host.appendChild(bar);
  }
  const v = localStorage.getItem(qcKey());
  bar.innerHTML = `<b>QC flag for ${curDate}:</b> ` + (v
    ? `<span class="qcchip" style="border-color:#f5c2c7;background:#fdecef">${v}</span>` +
      `<button id="qc-clear">clear</button>`
    : `<span class="qcchip">none</span>`) +
    ` <span class="cnt">— <span class="kbd">1</span> ${QC_REASONS[0]} · ` +
    `<span class="kbd">2</span> ${QC_REASONS[1]} · <span class="kbd">3</span> ${QC_REASONS[2]} · ` +
    `<span class="kbd">0</span> free text</span>`;
  const c = document.getElementById('qc-clear');
  if (c) c.addEventListener('click', () => setQC(null));
}

// ===================================================================== keyboard
// The production contract (monitoring/static/diag.js:257-271), preserved so the muscle memory
// carries over: left/right = calibrated days, Ctrl = every day with a diagnostic, up/down = station.
const daysWith = all => D.dates.filter(ds => (D.methods || []).some(m => {
  const s = summaryOf(ds, m);
  return s && (all ? (s.has_fig || s.constant != null) : s.constant != null);
}));
function stepDay(dir, all) {
  const list = daysWith(all);
  if (!list.length) return;
  let i = list.indexOf(curDate);
  if (i < 0) {                       // current day not in this list -> nearest in the direction
    const near = list.filter(d => dir > 0 ? d > curDate : d < curDate);
    if (!near.length) return;
    goDay(dir > 0 ? near[0] : near[near.length - 1]);
    return;
  }
  const j = i + dir;
  if (j >= 0 && j < list.length) goDay(list[j]);
}
document.addEventListener('keydown', e => {
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA')) return;
  const all = e.ctrlKey || e.metaKey;
  if (e.key === 'ArrowLeft')  { stepDay(-1, all); e.preventDefault(); }
  else if (e.key === 'ArrowRight') { stepDay(1, all); e.preventDefault(); }
  else if (e.key === 'ArrowUp')   { if (stepStation(-1)) e.preventDefault(); }
  else if (e.key === 'ArrowDown') { if (stepStation(1)) e.preventDefault(); }
  else if (e.key === '0') { const v = prompt('QC flag / comment for ' + curDate,
                              localStorage.getItem(qcKey()) || ''); if (v !== null) setQC(v || null);
                            e.preventDefault(); }
  else if (e.key === '1' || e.key === '2' || e.key === '3') {
    setQC(QC_REASONS[+e.key - 1]); e.preventDefault();
  }
});

// boot: land on the most recent day that produced a constant, else the most recent with a figure
(function () {
  const pref = daysWith(false), any = daysWith(true);
  const start = (pref.length ? pref : (any.length ? any : D.dates))[
    (pref.length ? pref : (any.length ? any : D.dates)).length - 1];
  if (start) goDay(start); else syncDayLinked();
})();
</script>
</body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", action="append", required=True)
    ap.add_argument("--db", type=Path, default=None,
                    help="optional dashboard SQLite; only used for status when L1 metadata exists")
    ap.add_argument("--cal-dir", type=Path, required=True)
    ap.add_argument("--status-dir", type=Path, default=None)
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--l1-root", default=None)
    ap.add_argument("--l1-month", default=None)
    ap.add_argument("--cams", default="")
    ap.add_argument("--bucket", default=BUCKET)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--payloads", action="store_true", help="compute the per-day JSON payloads")
    ap.add_argument("--force", action="store_true", help="recompute payloads that already exist")
    ap.add_argument("--work", type=Path, default=Path("./_dashboard_work"))
    # 32 cores here; leave a couple for the OS and for the Plotly/HTML step. Memory is
    # not the binding constraint (a night peaks well under 1 GB against 128 GB).
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    args = ap.parse_args()
    args.status_dir = args.status_dir or args.cal_dir

    if args.l1_root:
        os.environ["ALC_L1_ROOT"] = str(args.l1_root)
    if args.cams:
        os.environ["ALC_CAMS_DIR"] = args.cams
    os.environ["ALC_FULLCAL_DIR"] = str(args.work.resolve())
    os.environ["PLOTS"] = "1"
    args.work.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    d0 = datetime.strptime(args.start, "%Y%m%d")
    d1 = datetime.strptime(args.end, "%Y%m%d")
    days = [d0 + timedelta(days=i) for i in range((d1 - d0).days + 1)]

    # L1 FIRST. The dashboard SQLite only holds the streams its last build covered, so it is not a
    # station registry -- reading it as one made a normal station (the Payerne CL31) look missing.
    # The stream's own L1 file always carries its site, institution, type and coordinates.
    stations = _station_index(args.db, args.cal_dir, args.status_dir) if args.db else []
    by_key = {s["key"]: s for s in stations}
    for key in args.key:
        rec = _meta_from_l1(args.l1_root, key)
        if rec is None:
            rec = by_key.get(key)
            if rec is None:
                sys.exit(f"{key}: no L1 under {args.l1_root} and not in {args.db}")
            print(f"  {key}: no L1 found -> metadata from {args.db}", flush=True)
        else:
            db_rec = by_key.get(key)
            if db_rec:                       # keep the DB's status, prefer L1 for the identity
                rec["status"], rec["status_date"] = db_rec["status"], db_rec["status_date"]
            print(f"  {key}: {rec['itype']} '{rec['name']}' from {rec.pop('src')}", flush=True)
        cal = _cal_frame(args.cal_dir, key)
        ok = cal[cal["cal_value"] > 0] if len(cal) else cal
        for m, g in (ok.groupby("method") if len(ok) else []):
            last = g.sort_values("date").iloc[-1]
            rec["constants"][str(m)] = {"value": float(last["cal_value"]),
                                        "date": str(last["date"])}
        if key in by_key:
            stations[[s["key"] for s in stations].index(key)] = rec
        else:
            stations.append(rec)
        by_key[key] = rec
    stations.sort(key=lambda s: (s.get("country") or "", s.get("name") or "", s["key"]))

    for key in args.key:
        rec = by_key[key]
        print(f"\n=== {key} ({rec.get('itype')})", flush=True)
        idx_path = args.out / "data" / key / "_index.json"
        if args.payloads:
            index = _write_payloads(key, rec.get("itype") or "", days, args.out,
                                    args.force, args.workers)
            idx_path.parent.mkdir(parents=True, exist_ok=True)
            idx_path.write_text(json.dumps(index, separators=(",", ":")), encoding="utf-8")
        else:
            index = json.loads(idx_path.read_text(encoding="utf-8")) if idx_path.exists() else {}
            if not index:
                print("  (no payload index yet — run once with --payloads)", flush=True)
        html = build_page(key, rec, args, index, stations)
        p = args.out / f"station_{key}.html"
        p.write_text(html, encoding="utf-8")
        print(f"  -> {p}  ({p.stat().st_size / 1e6:.1f} MB, {len(index)} indexed days)", flush=True)

    (args.out / "index.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>ALC stations</title>"
        "<h1>ALC calibration — stations</h1><ul>" + "".join(
            f'<li><a href="station_{k}.html">{by_key[k].get("name") or k} · '
            f'{by_key[k].get("itype")} · unit {by_key[k].get("ident")}</a></li>'
            for k in args.key) + "</ul>", encoding="utf-8")
    print(f"\nindex: {args.out / 'index.html'}")


if __name__ == "__main__":
    main()
