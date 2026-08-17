#!/usr/bin/env python3
"""Mockup of the FULL station page: every control the operator asked for, on one page.

This is the assembly of the two half-mockups that came before it -- `mockup_station_card.py`
(availability card + instrument monitoring) and `mockup_daily_panel.py` (the unified daily
calibration panel) -- plus the station-navigation controls, so the whole page can be judged as one
thing rather than as fragments.

  compact header      NAME · COUNTRY · INSTITUTION · WMO · lat,lon · altitude, on one line
  station navigation  country + instrument dropdowns, prev/next arrows that RESPECT the filter,
                      and a station list showing each neighbour's status and calibration constant
  availability card   3 rows (instrument status / mean cloud cover / calibration), 4 for a CL61
  instrument monitor  the re-encoded x0/dx + float32 housekeeping figure
  daily calibration   the unified panel: calendar, day arrows, method switch, flags, the
                      time-height curtain and the profile selector, the diagnostics and messages

Nothing here reimplements those pieces: the card and the monitoring figure are drawn by
`monitoring.charts`, and the panel is `mockup_daily_panel`'s own CSS/body/JS fragments embedded
verbatim. What this file adds is the header, the navigation and the page shell.

  python scripts/mockup_station_page.py --key 0-20000-0-06610_C \
      --db      C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_pages/calib_index.sqlite \
      --cal-dir C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04 \
      --status-dir C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_04 \
      --dates 20260602,20260605,20260615,20260704,20260705,20260708 \
      --l1-root D:/E-PROFILE_L1_2026 --cams D:/CAMS_daily --out mock_station_page.html
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sqlite3
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import plotly.offline as pyo

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from monitoring import charts, config                                    # noqa: E402
from monitoring.render import _load_hk, _load_status                     # noqa: E402
from scripts.mockup_station_card import probe_cloud_cover, _month_with_data  # noqa: E402
import scripts.mockup_daily_panel as PANEL                               # noqa: E402


# ======================================================================= station index (for nav)
def _station_index(db: Path, cal_dir: Path, status_dir: Path) -> list[dict]:
    """Every station in the DB, with what the navigation has to show: status and constant.

    The operator's requirement was to move between stations AND see status + calibration constant
    without opening each one, so the index carries both. The constant is the LATEST per method: a
    station's headline number is what it is calibrated at now, not its mean over the archive.
    """
    with sqlite3.connect(db) as con:
        st = pd.read_sql("SELECT * FROM stations", con)
    out = []
    for _, r in st.iterrows():
        key = r["key"]
        rec = {"key": key, "wmo": r["wmo"], "ident": r["identifier"], "itype": r["itype"],
               "name": r["name"], "country": r["country"], "institution": r["institution"],
               "lat": None if pd.isna(r["lat"]) else float(r["lat"]),
               "lon": None if pd.isna(r["lon"]) else float(r["lon"]),
               "alt": None if pd.isna(r["alt"]) else float(r["alt"]),
               "constants": {}, "status": None, "status_date": None}
        cal_path = Path(cal_dir) / key / f"{key}_cal.csv"
        if cal_path.exists():
            cal = pd.read_csv(cal_path, dtype={"date": str})
            ok = cal[pd.to_numeric(cal["cal_value"], errors="coerce") > 0]
            for m, g in ok.groupby("method"):
                last = g.sort_values("date").iloc[-1]
                rec["constants"][str(m)] = {"value": float(last["cal_value"]),
                                            "date": str(last["date"])}
        sdf = _load_status(Path(status_dir), key)
        if sdf is not None and len(sdf):
            s = sdf.sort_values("date").iloc[-1]
            rec["status"] = str(s.get("quality", "") or "")
            rec["status_date"] = str(s.get("date", "") or "")
        out.append(rec)
    return out


# ================================================================================ page fragments
def _header_html(rec: dict) -> str:
    """The compact header: one title line and one meta line, nothing else.

    Replaces the old two-block header (title + a separate `.metabar` section). Every field is
    guarded -- an incomplete manifest must render as an em dash, never as the string 'None'.
    """
    def f(v, spec, suffix=""):
        return "—" if v is None else format(v, spec) + suffix
    bits = [rec.get("country") or "—", rec.get("institution") or "—", rec["wmo"],
            f"unit {rec['ident']}",
            f"{f(rec.get('lat'), '.3f')}, {f(rec.get('lon'), '.3f')}",
            f(rec.get("alt"), ".0f", " m")]
    return (f'<h1 id="st-title">{rec.get("name") or rec["wmo"]} '
            f'<span class="itype">{rec.get("itype") or ""}</span></h1>'
            f'<p class="metainline" id="st-meta">{" · ".join(bits)}</p>')


#: Tile thresholds. Each is (good_below, warn_below) on the absolute value of the metric, in the
#: metric's own unit; anything above the second bound is bad. They are display heuristics for an
#: operator's eye, NOT calibration acceptance criteria -- the flags decide that.
TILE_OK, TILE_WARN, TILE_BAD = "#1a7431", "#b7791f", "#b00020"


def _grade(v, good, warn):
    if v is None:
        return TILE_BAD
    a = abs(v)
    return TILE_OK if a < good else (TILE_WARN if a < warn else TILE_BAD)


def _cl_series(key: str, cal_dir: Path, ref_date: str | None = None):
    """C_L over time for every method, plus the five headline statistics.

    Reads the per-stream `_cal.csv` directly rather than the dashboard's SQLite, so the mockup works
    against any run tree. `charts.cl_overlay` draws the same overlay but as bare markers; this adds
    the per-night uncertainty, the median line and the +/-10 % band, which are what make a drift or
    a step visible at a glance.
    """
    import plotly.graph_objects as go
    path = Path(cal_dir) / key / f"{key}_cal.csv"
    if not path.exists():
        return None, [], ""
    cal = pd.read_csv(path, dtype={"date": str})
    cal["value"] = pd.to_numeric(cal["cal_value"], errors="coerce")
    cal["unc"] = pd.to_numeric(cal.get("uncertainty"), errors="coerce")
    cal["dt"] = pd.to_datetime(cal["date"], format="%Y%m%d", errors="coerce")
    ok = cal[(cal["value"] > 0) & cal["dt"].notna()].sort_values("dt")
    if not len(ok):
        return None, [], ""

    per = {m: g for m, g in ok.groupby("method") if len(g)}
    # "Reference" = the method the headline median describes. Prefer the cloud retrieval when it is
    # present: it runs on far more nights than Rayleigh, so its median is the better-determined one.
    ref = "cloud" if "cloud" in per else sorted(per, key=lambda m: -len(per[m]))[0]
    g = per[ref]
    med = float(g["value"].median())
    p10, p90 = (float(g["value"].quantile(q)) for q in (0.10, 0.90))
    spread = 100.0 * (p90 - p10) / med if med else None

    fig = go.Figure()
    for m, gm in per.items():
        col = config.METHOD_COLORS.get(m, "#888")
        fig.add_trace(go.Scatter(
            x=gm["dt"], y=gm["value"], mode="markers", name=config.method_label(m),
            marker=dict(size=5, color=col, opacity=0.85),
            error_y=dict(type="data", array=gm["unc"].fillna(0.0), visible=True,
                         color=col, thickness=0.9, width=0),
            hovertemplate="%{x|%Y-%m-%d}<br>C_L=%{y:.4g}<extra></extra>"))
    fig.add_hrect(y0=med * 0.9, y1=med * 1.1, fillcolor="rgba(70,130,140,0.10)",
                  line_width=0, layer="below")
    fig.add_hline(y=med, line=dict(color="#2a6b73", width=1.2))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name="±10 % of median",
                             marker=dict(size=9, color="rgba(70,130,140,0.25)", symbol="square")))

    last_dt = g["dt"].max()
    # The drift window is shaded rather than only quoted, so an operator can see WHICH nights the
    # +x % is computed from -- a drift over three nights and one over fifteen read very differently.
    win_start = last_dt - pd.Timedelta(days=15)
    fig.add_vrect(x0=win_start, x1=last_dt, fillcolor="rgba(240,173,78,0.16)", line_width=0,
                  layer="below")
    fig.add_vline(x=last_dt, line=dict(color="#d9534f", width=1.2))
    fig.update_layout(margin=dict(l=64, r=12, t=16, b=40), height=330,
                      template="plotly_white", font=dict(size=11),
                      yaxis_title="C_L", legend=dict(orientation="h", y=1.13, font=dict(size=10)))
    fig.update_yaxes(exponentformat="e")

    recent = g[g["dt"] >= win_start]["value"]
    prior = g[(g["dt"] < win_start) & (g["dt"] >= win_start - pd.Timedelta(days=45))]["value"]
    drift = (100.0 * (recent.median() / prior.median() - 1.0)
             if len(recent) and len(prior) and prior.median() else None)
    ratio = None
    if "cloud" in per and "rayleigh" in per:
        mr = float(per["rayleigh"]["value"].median())
        ratio = float(per["cloud"]["value"].median()) / mr if mr else None

    today = pd.to_datetime(ref_date, format="%Y%m%d") if ref_date else last_dt
    age = int((today - last_dt).days)
    tiles = [
        {"label": f"MEDIAN {config.method_label(ref).upper()} C_L", "value": f"{med:.4g}",
         "note": f"{len(g)} calibrated nights", "color": TILE_OK},
        {"label": "SPREAD (P10–P90)",
         "value": "—" if spread is None else f"{spread:.0f}%",
         "note": f"{p10:.4g} … {p90:.4g}", "color": _grade(spread, 10, 20)},
        {"label": "CLOUD / RAYLEIGH RATIO",
         "value": "—" if ratio is None else f"{ratio:.3g}×",
         "note": ("only one method on this stream" if ratio is None else
                  ("the two retrievals agree" if abs(ratio - 1) < 0.05
                   else "a systematic offset points at the water-vapour correction, "
                        "not at the instrument")),
         "color": TILE_BAD if ratio is None else _grade(ratio - 1.0, 0.05, 0.15)},
        {"label": "15-DAY DRIFT",
         "value": "—" if drift is None else f"{drift:+.1f}%",
         "note": "vs the preceding 45 d" if drift is not None else "not enough history",
         "color": _grade(drift, 5, 10)},
        {"label": "LAST VALID",
         "value": "today" if age <= 0 else (f"{age} d ago"),
         "note": str(g["date"].iloc[-1]), "color": _grade(age, 3, 8)},
    ]
    sub = ("Both retrievals on one axis — a systematic offset between them points at the "
           "water-vapour correction, not the instrument." if len(per) > 1
           else "One retrieval on this stream.")
    return fig, tiles, sub


def _avail_html(key: str, cal_dir: Path, status_dir: Path,
                l1_root: Path | None, l1_month: str | None) -> tuple[str, str]:
    """The availability card + instrument monitoring, drawn by the PRODUCTION chart functions."""
    status_df = _load_status(Path(status_dir), key)
    if status_df is None:
        return "", f"no {key}_status.csv under {status_dir}"
    cal_path = Path(cal_dir) / key / f"{key}_cal.csv"
    if not cal_path.exists():
        return "", f"no {cal_path.name} under {cal_dir}"
    cal = pd.read_csv(cal_path, dtype={"date": str})
    cal["key"] = key
    methods = [m for m in config.METHOD_ORDER if len(cal[cal["method"] == m])]

    note = "row 2 omitted — pass --l1-root/--l1-month to probe real octas"
    if l1_root and l1_month:
        month = _month_with_data(Path(l1_root), key, l1_month) or l1_month
        cc = probe_cloud_cover(Path(l1_root), key, month)
        if len(cc):
            status_df = status_df.merge(cc, on="date", how="left")
            note = (f"row 2 = REAL octas probed live from {len(cc)} L1 files of {month} "
                    f"(source '{cc['cloud_src'].iloc[0]}'); the rest of the record is blank because "
                    f"the operational producer does not write mean_cloud_cover yet")
        else:
            note = f"no L1 cloud_amount for {month} — row 2 omitted"

    fig = charts.daily_availability_rows(status_df, cal, methods)
    parts = ['<h2>Data availability &amp; instrument status</h2>',
             '<div class="avail-legend"><b>Status</b>'
             '<span class="ak ak-pass"></span> Pass'
             '<span class="ak ak-warning"></span> Warning'
             '<span class="ak ak-error"></span> Error'
             '<span class="ak ak-nodata"></span> No data'
             '<span class="legsep"></span><b>Cloud cover</b>'
             '<span class="ak ak-cc0"></span> 0'
             '<span class="ak ak-cc4"></span> 4'
             '<span class="ak ak-cc8"></span> 8 octas'
             '<span class="legsep"></span><b>Calibration</b>'
             + "".join(f'<span class="ak" style="background:{config.CAL_CLASS_COLORS[c]}" '
                       f'title="{config.CAL_CLASS_LABELS[c]}"></span> {config.CAL_CLASS_SHORT[c]}'
                       for c in config.CAL_CLASS_ORDER) + '</div>',
             f'<p class="muted">{note}</p>',
             f'<div class="card">{charts.fig_to_div(fig, "fig-avail")}</div>'
             if fig is not None else '<p class="muted">availability figure unavailable</p>']
    hk_df = _load_hk(Path(status_dir), key)
    if hk_df is not None:
        parts += ['<h2>Instrument monitoring</h2>',
                  f'<div class="card">'
                  f'{charts.fig_to_div(charts.monitoring_timeseries(hk_df), "fig-hk")}</div>']
    return "\n".join(parts), note


# ============================================================================================ page
PAGE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Station page — mockup</title>
<script>__PLOTLY__</script>
<style>
__STYLECSS__
__PANELCSS__
/* ---- page shell + the station-navigation controls (the part this mockup adds) ---- */
body { padding:0 0 30px; }
.wrap { padding:0 20px; }
h1 { font-size:22px; margin:14px 0 2px; }
h1 .itype { font-size:14px; font-weight:500; color:#66707a; margin-left:6px; }
.metainline { color:#66707a; font-size:13px; margin:0 0 12px; }
h2 { font-size:13px; margin:20px 0 8px; text-transform:uppercase; letter-spacing:.04em;
     color:#66707a; }
.navbar { display:flex; align-items:center; gap:10px; flex-wrap:wrap; padding:9px 11px;
          background:#fff; border:1px solid #dbe3ea; border-radius:10px; margin-bottom:14px; }
.navbar select { font:inherit; font-size:13px; padding:4px 8px; border:1px solid #c3ceda;
                 border-radius:7px; background:#fff; }
.navbar .count { font-size:12px; color:#66707a; }
.stlist { display:flex; gap:7px; flex-wrap:wrap; margin-top:9px; }
.stcard { border:1px solid #dbe3ea; border-radius:8px; padding:6px 9px; font-size:12px;
          background:#fff; cursor:pointer; min-width:150px; }
.stcard:hover { border-color:#2a5a82; background:#f4f9ff; }
.stcard.cur { border-color:#1a2530; border-width:2px; background:#eef4fb; }
.stcard .nm { font-weight:600; }
.stcard .cc { color:#66707a; }
.stcard .cl { font-variant-numeric:tabular-nums; }
.dot { display:inline-block; width:8px; height:8px; border-radius:50%; margin-right:5px; }
.sect { border-top:1px solid #e3e9ef; margin-top:22px; padding-top:4px; }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:10px;
         margin-top:11px; }
.tile { border:1px solid #dbe3ea; border-radius:9px; padding:9px 12px; background:#fff; }
.tile .tl { font:11px/1.3 ui-monospace,Consolas,monospace; color:#66707a;
            letter-spacing:.03em; }
.tile .tv { font-size:23px; font-weight:600; margin:3px 0 1px; font-variant-numeric:tabular-nums; }
.tile .tn { font-size:11.5px; color:#66707a; }
.subtle { color:#66707a; font-size:12.5px; margin:0 0 8px; }
</style></head><body>
<div class="hdr"><b>MOCKUP — full station page.</b> The availability card and the monitoring figure
are drawn by the production <code>monitoring.charts</code> functions with the real
<code>style.css</code>; the daily panel below is the same CSS/body/JS as the standalone panel
mockup, embedded rather than reimplemented. <span id="meta"></span></div>

<div class="wrap">
__HEADER__

<div class="navbar">
  <button id="st-prev">← previous station</button>
  <label class="count">Country <select id="f-country"></select></label>
  <label class="count">Instrument <select id="f-type"></select></label>
  <button id="st-next">next station →</button>
  <span class="count" id="st-pos"></span>
</div>
<div class="stlist" id="st-list"></div>

<div class="sect">
<h2>Calibration coefficient C_L over time</h2>
<p class="subtle">__CLSUB__</p>
<div class="card">__CLFIG__</div>
<div class="tiles">__CLTILES__</div>
</div>

<div class="sect">
__AVAIL__
</div>

<div class="sect">
<h2>Daily calibration</h2>
</div>
</div>
__PANELBODY__

<script id="stations" type="application/json">__STATIONS__</script>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
// ------------------------------------------------------------------ station navigation
// The arrows walk the FILTERED list, which is the whole point of pairing them with the dropdowns:
// arrows that ignored the active filter would jump out of the selection the operator just made.
const ST = JSON.parse(document.getElementById('stations').textContent);
const CUR = __CURKEY__;
const DOTC = { pass:'#2ca02c', warning:'#f0ad4e', error:'#d9534f', '':'#c3ceda' };
const selC = document.getElementById('f-country'), selT = document.getElementById('f-type');

function opts(sel, values, label) {
  sel.innerHTML = `<option value="">${label} (all)</option>` +
    values.map(v => `<option value="${v}">${v}</option>`).join('');
}
opts(selC, [...new Set(ST.map(s => s.country).filter(Boolean))].sort(), 'Country');
opts(selT, [...new Set(ST.map(s => s.itype).filter(Boolean))].sort(), 'Instrument');

const filtered = () => ST.filter(s =>
  (!selC.value || s.country === selC.value) && (!selT.value || s.itype === selT.value));

function fmtC(s) {
  const ks = Object.keys(s.constants || {});
  if (!ks.length) return '<span class="cc">no constant</span>';
  return ks.map(m => `<span class="cl">${m[0].toUpperCase()}: ` +
    `${(+s.constants[m].value).toPrecision(4)}</span>`).join(' · ');
}
function drawList() {
  const list = filtered();
  const host = document.getElementById('st-list');
  host.innerHTML = list.map(s => {
    const c = DOTC[(s.status || '').toLowerCase()] || DOTC[''];
    return `<div class="stcard${s.key === CUR ? ' cur' : ''}" data-key="${s.key}">
      <div class="nm"><span class="dot" style="background:${c}"></span>${s.name || s.wmo}</div>
      <div class="cc">${s.country || '—'} · ${s.itype || '—'} · unit ${s.ident}</div>
      <div>${fmtC(s)}</div></div>`;
  }).join('');
  const i = list.findIndex(s => s.key === CUR);
  document.getElementById('st-prev').disabled = i <= 0;
  document.getElementById('st-next').disabled = i < 0 || i >= list.length - 1;
  document.getElementById('st-pos').textContent = i >= 0
    ? `${i + 1} of ${list.length} in filter`
    // Never leave the arrows pointing nowhere without saying why: the operator filtered this
    // station out, which is a different situation from "there is nothing here".
    : `this station is outside the filter (${list.length} match${list.length === 1 ? '' : 'es'})`;
  host.querySelectorAll('.stcard').forEach(e => e.addEventListener('click', () =>
    e.dataset.key === CUR ? null
      : alert('In production this loads ' + e.dataset.key +
              '.\nThis mockup only carries the calibration payload for ' + CUR + '.')));
}
[selC, selT].forEach(s => s.addEventListener('change', drawList));
document.getElementById('st-prev').addEventListener('click', () => {
  const l = filtered(), i = l.findIndex(s => s.key === CUR);
  if (i > 0) alert('previous station: ' + l[i - 1].key);
});
document.getElementById('st-next').addEventListener('click', () => {
  const l = filtered(), i = l.findIndex(s => s.key === CUR);
  if (i >= 0 && i < l.length - 1) alert('next station: ' + l[i + 1].key);
});
drawList();
</script>
<script>__PANELJS__</script>
</body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", required=True)
    ap.add_argument("--db", type=Path, required=True, help="calib_index.sqlite (stations table)")
    ap.add_argument("--cal-dir", type=Path, required=True)
    ap.add_argument("--status-dir", type=Path, default=None)
    ap.add_argument("--dates", required=True, help="comma-separated YYYYMMDD for the daily panel")
    ap.add_argument("--methods", default="rayleigh,cloud")
    ap.add_argument("--l1-root", required=True)
    ap.add_argument("--l1-month", default=None)
    ap.add_argument("--cams", default="")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("./_station_page_work"))
    args = ap.parse_args()
    status_dir = args.status_dir or args.cal_dir

    os.environ["ALC_L1_ROOT"] = str(args.l1_root)
    os.environ["ALC_FULLCAL_DIR"] = str(args.work.resolve())
    os.environ["PLOTS"] = "1"
    if args.cams:
        os.environ["ALC_CAMS_DIR"] = args.cams
    args.work.mkdir(parents=True, exist_ok=True)

    print("station index ...", flush=True)
    index = _station_index(args.db, args.cal_dir, status_dir)
    rec = next((s for s in index if s["key"] == args.key), None)
    if rec is None:
        sys.exit(f"{args.key} is not in {args.db}")

    print("C_L series + headline tiles ...", flush=True)
    cl_fig, cl_tiles, cl_sub = _cl_series(args.key, args.cal_dir)
    cl_html = (charts.fig_to_div(cl_fig, "fig-cl") if cl_fig is not None
               else '<p class="muted">no calibrated night in this tree</p>')
    tiles_html = "".join(
        f'<div class="tile"><div class="tl">{t["label"]}</div>'
        f'<div class="tv" style="color:{t["color"]}">{t["value"]}</div>'
        f'<div class="tn">{t["note"]}</div></div>' for t in cl_tiles)

    print("availability card + monitoring ...", flush=True)
    avail, note = _avail_html(args.key, args.cal_dir, status_dir,
                              Path(args.l1_root), args.l1_month)

    print("daily calibration panel ...", flush=True)
    day_list = [datetime.strptime(s.strip(), "%Y%m%d") for s in args.dates.split(",") if s.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    payload, t0 = PANEL.build_payload(args.key, day_list, methods)
    if not payload["dates"]:
        sys.exit("no day produced a panel")
    raw = json.dumps(payload, separators=(",", ":"))
    gz_kb = len(gzip.compress(raw.encode())) / 1024
    n = len(payload["dates"])
    meta = {"n_days": n, "gz_kb": round(gz_kb, 1),
            "per_unit_kb": round(gz_kb / max(n * len(methods), 1), 1),
            "curated": True, "elapsed_s": round(time.perf_counter() - t0, 1)}

    css = (REPO / "monitoring/static/style.css").read_text(encoding="utf-8")
    # __META__ lives INSIDE the panel's JS, so it has to be substituted there before that JS is
    # spliced into the page -- doing it on the page first would run against a token that has not
    # arrived yet, and the browser would hit an undefined identifier.
    panel_js = PANEL.PANEL_JS.replace("__META__", json.dumps(meta))
    # Order matters: the panel's CSS is loaded AFTER style.css so its rules win inside the panel,
    # and the page shell is declared last so it wins over both.
    html = (PAGE
            .replace("__STYLECSS__", css)
            .replace("__PANELCSS__", PANEL.PANEL_CSS)
            .replace("__HEADER__", _header_html(rec))
            .replace("__CLSUB__", cl_sub)
            .replace("__CLFIG__", cl_html)
            .replace("__CLTILES__", tiles_html)
            .replace("__AVAIL__", avail or f'<p class="muted">{note}</p>')
            .replace("__PANELBODY__", PANEL.PANEL_BODY)
            .replace("__CURKEY__", json.dumps(args.key))
            .replace("__STATIONS__", json.dumps(index, separators=(",", ":")))
            .replace("__PANELJS__", panel_js)
            .replace("__PAYLOAD__", raw)
            .replace("__PLOTLY__", pyo.get_plotlyjs()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    print(f"\nstations in nav : {len(index)} "
          f"({len({s['country'] for s in index})} countries, "
          f"{len({s['itype'] for s in index})} instrument types)")
    print(f"panel           : {n} days x {len(methods)} methods, {gz_kb:.0f} KB gzip")
    print(f"page            : {args.out}  ({args.out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
