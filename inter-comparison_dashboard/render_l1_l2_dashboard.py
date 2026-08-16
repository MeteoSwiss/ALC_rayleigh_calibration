# -*- coding: utf-8 -*-
"""Payerne L1-vs-L2 dashboard — render stage.

Turns the JSON from build_l1_l2_dashboard.py into ONE self-contained HTML page (plotly.min.js and
the data are inlined, so the file opens anywhere with no server and no network). All the toggles are
client-side over precomputed arrays — nothing is recomputed in the browser.

Run:  python inter-comparison_dashboard/render_l1_l2_dashboard.py
Out:  inter-comparison_dashboard/index.html  (data read from C:/DATA/.../inter-comparison_dashboard)
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")                # validation.paper, monitoring
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling modules (folder name is not a valid package name)
from plotly.offline import get_plotlyjs

DATA = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard")  # inputs
HERE = Path(__file__).resolve().parent
PAGE = HERE / "index.html"   # the ONE page, kept in the project next to the code that builds it

CSS = """
:root { --bg:#f7f8fa; --card:#fff; --ink:#1c2330; --muted:#6b7280; --line:#e3e6ea;
        --brand:#0b3d61; --accent:#1f77b4; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:14px/1.45 -apple-system,"Segoe UI",Roboto,Arial,sans-serif; }
.topbar { display:flex; align-items:center; gap:12px; background:var(--brand); color:#fff;
          padding:10px 18px; position:sticky; top:0; z-index:50; }
.topbar .brand { font-weight:600; }
.topbar .spacer { flex:1; }
.topbar .asof { color:#d6e6f2; font-size:13px; }
.container { max-width:1560px; margin:0 auto; padding:18px; }
h1 { font-size:22px; margin:8px 0 4px; }
h2 { font-size:17px; margin:24px 0 10px; }
.subtitle { margin:0 0 14px; font-size:15px; color:#33405a; font-weight:500; }
.muted { color:var(--muted); font-weight:400; font-size:13px; }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
.controls { display:flex; flex-wrap:wrap; align-items:center; gap:22px; margin:14px 0 18px;
            position:sticky; top:44px; z-index:40; box-shadow:0 2px 10px rgba(0,0,0,.04); }
.ctl-group { display:flex; align-items:center; gap:10px; }
.ctl-label { font-weight:600; font-size:13px; color:#33405a; }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:8px; overflow:hidden; }
.seg button { font:inherit; font-size:13px; padding:5px 12px; background:#fff; color:var(--ink);
              border:0; border-right:1px solid var(--line); cursor:pointer; }
.seg button:last-child { border-right:0; }
.seg button.on { background:var(--accent); color:#fff; font-weight:600; }
.chk { display:inline-flex; align-items:center; gap:7px; cursor:pointer; font-size:13px; }
.chk input { width:16px; height:16px; accent-color:var(--accent); cursor:pointer; }
select { font:inherit; font-size:13px; padding:5px 8px; border:1px solid var(--line);
         border-radius:7px; background:#fff; color:var(--ink); cursor:pointer; }
select:focus { outline:none; border-color:var(--accent); }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
.panel-title { font-weight:600; font-size:14px; margin:0 0 2px; }
.panel-sub { font-size:12px; color:var(--muted); margin:0 0 6px; }
.pill { display:inline-block; font-size:11px; padding:1px 9px; border-radius:999px;
        border:1px solid var(--line); color:var(--muted); margin-left:6px; }
.pill.l1 { background:#eef4fb; border-color:#cfe0f0; color:#17537f; }
.pill.l2 { background:#f4f1ee; border-color:#e6ddd3; color:#7a5a37; }
table.stats { border-collapse:collapse; width:100%; font-size:13px; }
table.stats th, table.stats td { border-bottom:1px solid var(--line); padding:6px 10px; text-align:right; }
table.stats th:first-child, table.stats td:first-child { text-align:left; }
table.stats thead th { color:var(--muted); font-weight:600; font-size:12px; }
table.stats tbody tr:hover { background:#f6f9fc; }
.num { font-variant-numeric:tabular-nums; }
.good { color:#1a7f37; font-weight:600; }
.bad  { color:#c92a2a; font-weight:600; }
.note { font-size:12.5px; color:var(--muted); margin:8px 0 0; }
.warn { background:#fff8e6; border:1px solid #f0dca8; border-left:4px solid #e0a800;
        border-radius:8px; padding:10px 14px; margin:14px 0; font-size:13px; color:#4a3a10; }
.warn b { color:#3a2d08; }
.warn + .warn { margin-top:8px; }
.foot { color:var(--muted); font-size:12px; padding:22px 18px; text-align:center; }
.foot code { font-size:11px; }
"""

JS = r"""
const CH = D.channels, ALT = D.alt_agl, IREF = D.iref, M = D.meta, MONTHS = D.months;
const state = { cal:'v2.0', wv:true, wl:'molecular', logx:false, filter:true, r0:0, r1:MONTHS.length - 1 };
const PCFG = { displaylogo:false, responsive:true };
const BASE = { template:'plotly_white', margin:{l:64,r:16,t:34,b:46}, font:{size:12}, height:470,
               hovermode:'closest', legend:{orientation:'h', y:-0.16, x:0.5, xanchor:'center',
               yanchor:'top', font:{size:11}, traceorder:'normal'} };
const ZTOP = 4000;   /* default view; the data runs to 6 km, zoom out to see it */

const key = () => state.cal + '_wv' + (state.wv ? 1 : 0) + '_' + state.wl;
const cur = () => D.combos[key()]['r' + state.r0 + '_' + state.r1];
const rgba = (hex, a) => {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${n >> 16 & 255},${n >> 8 & 255},${n & 255},${a})`;
};
/* the noise-floor filter is applied here, not baked into the data, so it can be switched off */
const mask = (arr, keep) => state.filter ? arr.map((v, i) => keep[i] ? v : null) : arr;

/* ---------- vertical profiles: one panel per source, median + IQR band ---------- */
function profileTraces(src) {
  const P = cur()[src].profiles, out = [];
  CH.forEach((c, k) => {
    const p = P[k];
    out.push({ x:mask(p.q1, p.keep), y:ALT, mode:'lines', line:{width:0}, showlegend:false,
               hoverinfo:'skip', name:c.label + ' q1' });
    out.push({ x:mask(p.q3, p.keep), y:ALT, mode:'lines', line:{width:0}, fill:'tonextx',
               fillcolor:rgba(c.color, 0.13), showlegend:false, hoverinfo:'skip',
               name:c.label + ' IQR' });
    out.push({ x:mask(p.med, p.keep), y:ALT, mode:'lines', line:{color:c.color, width:2.2},
               name:c.label,
               hovertemplate:'%{y:.0f} m<br>' + c.label + ' = %{x:.3f}<extra></extra>' });
  });
  return out;
}

function profileLayout(title) {
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.title = { text:title, font:{size:13} };
  lay.xaxis = state.logx
    ? { title:'β<sub>att</sub> [Mm⁻¹ sr⁻¹]', type:'log', range:[-2.3, 0.3], zeroline:false }
    : { title:'β<sub>att</sub> [Mm⁻¹ sr⁻¹]', range:[0, 1.0], zeroline:false };
  lay.yaxis = { title:'Altitude a.g.l. [m]', range:[0, ZTOP] };
  lay.shapes = [zoneShape()];
  return lay;
}

/* the 500–3000 m statistics band, drawn on every altitude plot */
function zoneShape() {
  return { type:'rect', xref:'paper', yref:'y', x0:0, x1:1, y0:M.zmin, y1:M.zmax,
           fillcolor:'rgba(31,119,180,0.045)', line:{width:0}, layer:'below' };
}

/* ---------- relative-difference profiles: both sources on one axis ---------- */
function diffTraces() {
  const out = [];
  [['L2', 1.5, 'dash'], ['L1', 2.6, 'solid']].forEach(([src, w, dash]) => {
    const P = cur()[src];
    CH.forEach((c, k) => {
      if (k === IREF) return;
      out.push({ x:mask(P.diff_vs_ref[k], P.diff_keep[k]), y:ALT, mode:'lines',
                 line:{color:c.color, width:w, dash:dash},
                 name:c.label + ' — ' + src,
                 hovertemplate:'%{y:.0f} m<br>' + c.label + ' ' + src + ' %{x:+.1f}%<extra></extra>' });
    });
  });
  return out;
}

/* ---------- histograms of the per-sample relative difference over the stats band ---------- */
function histTraces() {
  const E = D.hist_edges, out = [];
  const mid = E.slice(0, -1).map((e, i) => (e + E[i + 1]) / 2);
  [['L2', 1.5, 'dash'], ['L1', 2.4, 'solid']].forEach(([src, w, dash]) => {
    const P = cur()[src];
    CH.forEach((c, k) => {
      if (k === IREF) return;
      const h = P.hist[k], tot = h.reduce((a, b) => a + b, 0);
      if (!tot) return;
      const dens = h.map(v => v / tot * 100);
      out.push({ x:mid, y:dens, mode:'lines', line:{color:c.color, width:w, dash:dash, shape:'hvh'},
                 name:c.label + ' — ' + src + ' (med ' + (P.hist_med[k] === null ? '—' :
                      (P.hist_med[k] > 0 ? '+' : '') + P.hist_med[k].toFixed(1) + '%') + ')',
                 hovertemplate:'%{x:+.1f}%<br>%{y:.2f}% of samples<extra></extra>' });
    });
  });
  return out;
}

function histLayout() {
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.title = { text:'Difference vs ' + CH[IREF].label + ' — distribution over ' +
                     M.zmin.toFixed(0) + '–' + M.zmax.toFixed(0) + ' m', font:{size:13} };
  lay.xaxis = { title:'relative difference [%]', range:[-100, 150], zeroline:true,
                zerolinewidth:1.4, zerolinecolor:'#888' };
  lay.yaxis = { title:'share of samples [%]' };
  lay.shapes = [{ type:'line', xref:'x', yref:'paper', x0:0, x1:0, y0:0, y1:1,
                  line:{color:'#888', width:1.4} }];
  return lay;
}

function diffLayout(title, xr) {
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.title = { text:title, font:{size:13} };
  lay.xaxis = { title:'relative difference [%]', range:xr, zeroline:true, zerolinewidth:1.4,
                zerolinecolor:'#888' };
  lay.yaxis = { title:'Altitude a.g.l. [m]', range:[0, ZTOP] };
  lay.shapes = [zoneShape()];
  return lay;
}

/* ---------- L1 vs L2, same instrument, same hour ---------- */
function l1l2Traces() {
  const C_ = cur();
  return CH.map((c, k) => ({ x:mask(C_.l1_vs_l2[k], C_.l1_vs_l2_keep[k]), y:ALT, mode:'lines',
    line:{color:c.color, width:2.2}, name:c.label,
    hovertemplate:'%{y:.0f} m<br>' + c.label + ': L1 %{x:+.1f}% vs L2<extra></extra>' }));
}

/* ---------- calibration coefficient time series (station-page style) ---------- */
const METHOD_COLOR = { 'Rayleigh':'#1f77b4', 'Liquid clouds':'#2ca02c' };

function calibFigure(ident, chan) {
  chan = chan || ident;
  const c = D.calib[chan];
  if (!c) return null;
  const c22 = (D.calib22 || {})[chan];
  const traces = [];
  const op = D.l2_applied[ident];
  if (op) traces.push({ x:op.date, y:op.value, mode:'lines', name:'Applied in L2',
    line:{color:'#111', width:1.3},
    hovertemplate:'%{x|%Y-%m-%d}<br>applied in L2 = %{y:.3e}<extra></extra>' });
  const k = c.kalman;
  if (k.date.length) {
    traces.push({ x:k.date.concat(k.date.slice().reverse()),
      y:k.value.map((v, i) => v + k.std[i]).concat(k.value.map((v, i) => v - k.std[i]).reverse()),
      fill:'toself', fillcolor:'rgba(214,39,40,0.12)', line:{width:0}, hoverinfo:'skip',
      showlegend:false });
    traces.push({ x:k.date, y:k.value, mode:'lines', name:'v2.0 Kalman estimate',
      line:{color:'#d62728', width:2},
      hovertemplate:'%{x|%Y-%m-%d}<br>Kalman = %{y:.3e}<extra></extra>' });
  }
  /* v2.2 alongside, only where it differs -- the cloud-calibrated channels are the SAME object
     under both variants, so drawing them twice would just be a duplicate line. */
  if (c22 && c22.kalman && c22.key === c.key && c22.created !== c.created) {
    const k2 = c22.kalman;
    if (k2.date.length) traces.push({ x:k2.date, y:k2.value, mode:'lines',
      name:'v2.2 Kalman estimate', line:{color:'#2f9e44', width:2},
      hovertemplate:'%{x|%Y-%m-%d}<br>v2.2 Kalman = %{y:.3e}<extra></extra>' });
    traces.push({ x:c22.points.map(p => p.date), y:c22.points.map(p => p.value), mode:'markers',
      name:'v2.2 — Rayleigh', marker:{ size:5, color:'#2f9e44', opacity:0.55, symbol:'triangle-up' },
      hovertemplate:'%{x|%Y-%m-%d}<br>v2.2 C_L = %{y:.3e}<extra></extra>' });
  }
  /* v2.2+dark: the measured electronic baseline subtracted inside the calibration. Same skip rule:
     cloud channels are carried over unchanged, so only genuinely different series are drawn. */
  const c22d = (D.calib22dark || {})[chan];
  if (c22d && c22d.kalman && c22d.key === c.key && c22d.created !== c.created) {
    const k3 = c22d.kalman;
    if (k3.date.length) traces.push({ x:k3.date, y:k3.value, mode:'lines',
      name:'v2.2+dark Kalman estimate', line:{color:'#7048e8', width:2},
      hovertemplate:'%{x|%Y-%m-%d}<br>v2.2+dark Kalman = %{y:.3e}<extra></extra>' });
    traces.push({ x:c22d.points.map(p => p.date), y:c22d.points.map(p => p.value), mode:'markers',
      name:'v2.2+dark — Rayleigh', marker:{ size:5, color:'#7048e8', opacity:0.5, symbol:'diamond' },
      hovertemplate:'%{x|%Y-%m-%d}<br>v2.2+dark C_L = %{y:.3e}<extra></extra>' });
  }
  const byM = {};
  c.points.forEach(p => (byM[p.method] = byM[p.method] || []).push(p));
  Object.keys(byM).sort().forEach(m => {
    const pts = byM[m];
    traces.push({ x:pts.map(p => p.date), y:pts.map(p => p.value), mode:'markers',
      name:'v2.0 — ' + m,
      marker:{ size:5.5, color:METHOD_COLOR[m] || '#1f77b4', opacity:0.8 },
      error_y:{ type:'data', array:pts.map(p => p.std || 0), visible:true, thickness:0.6, width:0,
                color:'rgba(120,120,120,0.35)' },
      hovertemplate:'%{x|%Y-%m-%d}<br>C_L = %{y:.3e}<extra>' + m + '</extra>' });
  });
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.height = 300;
  lay.title = { text:c.itype + ' (' + ident + ') — C_L over time', font:{size:13} };
  lay.yaxis = { title:'C_L', exponentformat:'e' };
  lay.xaxis = { title:'' };
  // shade the period the profile panels above are currently built from
  lay.shapes = [{ type:'rect', xref:'x', yref:'paper', x0:MONTHS[state.r0] + '-01',
                  x1:monthEnd(MONTHS[state.r1]), y0:0, y1:1,
                  fillcolor:'rgba(31,119,180,0.10)', line:{width:0}, layer:'below' }];
  return { traces, lay };
}

/* last day of a 'YYYY-MM' month, as 'YYYY-MM-DD' */
function monthEnd(ym) {
  const [y, m] = ym.split('-').map(Number);
  return new Date(Date.UTC(y, m, 0)).toISOString().slice(0, 10);
}

const MONTH_LABEL = ym => {
  const [y, m] = ym.split('-').map(Number);
  return ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][m - 1] + ' ' + y;
};

/* ---------- statistics table ---------- */
function statsTable() {
  const rows = [];
  ['L1', 'L2'].forEach(src => {
    CH.forEach((c, k) => {
      if (k === IREF) return;
      const s = cur()[src].stats[k];
      const cls = Math.abs(s.medrelbias_pct) <= 5 ? 'good' : (Math.abs(s.medrelbias_pct) > 20 ? 'bad' : '');
      rows.push(`<tr><td>${src === 'L1' ? 'L1 + v2 calibration' : 'L2 as distributed'}</td>
        <td>${c.label}</td>
        <td class="num ${cls}">${s.medrelbias_pct === null ? '—' : s.medrelbias_pct.toFixed(1) + ' %'}</td>
        <td class="num">${s.r_log === null ? '—' : s.r_log.toFixed(2)}</td>
        <td class="num">${s.relbias_pct === null ? '—' : s.relbias_pct.toFixed(1) + ' %'}</td>
        <td class="num">${s.n === null ? '—' : s.n.toLocaleString()}</td></tr>`);
    });
  });
  document.querySelector('#stats tbody').innerHTML = rows.join('');
  document.getElementById('nhours').textContent = cur().n_hours.toLocaleString();
  document.getElementById('period').textContent =
    state.r0 === state.r1 ? MONTH_LABEL(MONTHS[state.r0])
                          : MONTH_LABEL(MONTHS[state.r0]) + ' – ' + MONTH_LABEL(MONTHS[state.r1]);
}

/* ---------- draw / redraw ---------- */
function draw() {
  Plotly.react('p_l1', profileTraces('L1'), profileLayout('L1 + v2 calibration (Kalman)'), PCFG);
  Plotly.react('p_l2', profileTraces('L2'), profileLayout('L2 as distributed'), PCFG);
  Plotly.react('d_diff', diffTraces(),
               diffLayout('Difference vs ' + CH[IREF].label + ' — L1 solid, L2 dashed', [-100, 100]),
               PCFG);
  Plotly.react('p_hist', histTraces(), histLayout(), PCFG);
  const ll = diffLayout('L1 + v2 calibration relative to L2 — same instrument, same hour', [-100, 400]);
  ll.height = 500;
  Plotly.react('d_l1l2', l1l2Traces(), ll, PCFG);
  statsTable();
  CH.forEach(c => {
    const f = calibFigure(c.ident, c.chan || c.ident);
    if (f) Plotly.react('cal_' + (c.chan || c.ident), f.traces, f.lay, PCFG);
  });
}

/* ---------- controls ---------- */
function fillMonths() {
  ['r0', 'r1'].forEach(which => {
    const s = document.getElementById(which);
    s.innerHTML = MONTHS.map((m, i) =>
      `<option value="${i}">${MONTH_LABEL(m)}</option>`).join('');
    s.value = state[which];
  });
}

/* keep from <= to, snapping the other end when the user crosses over */
function onRange(which, v) {
  state[which] = v;
  if (state.r0 > state.r1) { state[which === 'r0' ? 'r1' : 'r0'] = v; }
  document.getElementById('r0').value = state.r0;
  document.getElementById('r1').value = state.r1;
  draw();
}

function bind() {
  document.getElementById('wv').addEventListener('change', e => { state.wv = e.target.checked; draw(); });
  /* Calibration variant: v2.2 changes the RAYLEIGH retrieval only, so on this page it moves the
     CHM15k and the CL61-Rayleigh channels and leaves the cloud-calibrated ones identical. */
  document.querySelectorAll('#calseg button').forEach(b => b.addEventListener('click', () => {
    state.cal = b.dataset.cal;
    document.querySelectorAll('#calseg button').forEach(x => x.classList.toggle('on', x === b));
    draw();
  }));
  document.getElementById('logx').addEventListener('change', e => { state.logx = e.target.checked; draw(); });
  document.getElementById('nfilt').addEventListener('change', e => { state.filter = e.target.checked; draw(); });
  document.getElementById('r0').addEventListener('change', e => onRange('r0', +e.target.value));
  document.getElementById('r1').addEventListener('change', e => onRange('r1', +e.target.value));
  document.querySelectorAll('#wlseg button').forEach(b => b.addEventListener('click', () => {
    state.wl = b.dataset.wl;
    document.querySelectorAll('#wlseg button').forEach(x => x.classList.toggle('on', x === b));
    draw();
  }));
}

fillMonths(); bind(); draw();
window.addEventListener('resize', () => ['p_l1','p_l2','d_diff','p_hist','d_l1l2']
  .forEach(id => Plotly.Plots.resize(document.getElementById(id))));
"""


def html(data):
    m = data["meta"]
    ch = data["channels"]
    cal_divs = "\n".join(
        f'<div class="card" style="margin-bottom:12px"><div id="cal_{c.get("chan", c["ident"])}"></div></div>'
        for c in ch)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Payerne — L1 + v2 calibration vs L2</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<div class="topbar">
  <span class="brand">E-PROFILE ALC calibration</span>
  <span class="spacer"></span>
  <span class="asof">Payerne · {m['wmo']} · {m['start']} → {m['end']}</span>
</div>
<div class="container">
  <h1>Payerne — three co-located ceilometers: L1 + v2 calibration vs the L2 product</h1>
  <p class="subtitle">Attenuated backscatter from the same hours, the same altitude grid and the
    same corrections on both sides. <span class="muted">Left = L1 <code>rcs_0</code> divided by the
    daily Kalman C<sub>L</sub> from the latest v2.0 calibration NetCDFs. Right = L2
    <code>attenuated_backscatter_0</code> exactly as distributed (CHM15k v1.0 Rayleigh, CL31 the
    1e8 default, CL61 Vaisala's internal constant). Read from the daily
    <code>E-PROFILE_L1_2026</code> / <code>E-PROFILE_L2_2026</code> archives, falling back to the
    5-minute granules on days with no concatenated file.</span></p>

  <div class="warn"><b>CL31 optical block replaced 2026-07-07 ~13:00</b> (after that day's dark
    measurement). Its per-night constant steps from ~3.1e7 (Jan–May) / 4.1e7 (Jun) to 9.1e7 (Jul–Aug).
    <b>This window spans the swap on purpose</b> — set the period to <b>Jun 2026</b> for the state
    before it and <b>Jul</b>/<b>Aug</b> for after. The smoothing is the operational dashboard's own
    filter (<code>monitoring.kalman</code>: median-normalised, 4 % day-to-day drift, 15 %/yr over
    gaps, rolling-IQR outlier rejection), which absorbs the step within days — 4.11e7 on 1 Jul,
    8.46e7 by 20 Jul, 8.00e7 on 6 Aug, identical to the station page — so the L1 panel stays valid
    either side of it. The L2 panel does not: it carries the fixed 1e8 default throughout, which
    happens to sit near the post-swap truth and far from the pre-swap one.</div>
  <div class="warn"><b>The CHM15k reference is under-constrained over this window.</b> Exactly
    <b>one</b> Rayleigh night is accepted between 8 Jul and 13 Aug (22 Jul, 8.28e11), after a
    two-month gap from 22 May (6.91e11), so the random-walk filter barely moves: it holds
    <b>6.84e11</b> — and the filter's rolling-IQR test rejects that night outright, so its series
    ends 22 May and July/August run on the clamped end value. Scale the reference to the 22 Jul night
    and CL61 moves from −14 % to about +3 %, so most of the residual sits in the reference rather
    than in CL61, whose own cloud calibration has 82 nights and tracks well. Constants come only from
    the v2.0 NetCDFs, all available years (43 CHM15k / 264 CL31 / 82 CL61 nights); the raw
    operational <code>_cal.csv</code> is deliberately not used, as its flag=0.5 nights include two
    wild CHM15k outliers (20–21 Jun, 3.4–4.1e12) that the published NetCDF drops.</div>

  <div class="card controls">
    <div class="ctl-group">
      <span class="ctl-label">Calibration</span>
      <span class="seg" id="calseg">
        <button data-cal="v2.0" class="on">v2.0 (operational)</button>
        <button data-cal="v2.2">v2.2 (noise-aware gates)</button>
        <button data-cal="v2.2dark">v2.2 + dark (electronic baseline subtracted)</button>
      </span>
    </div>
    <div class="ctl-group">
      <span class="ctl-label">Corrections</span>
      <label class="chk"><input type="checkbox" id="wv" checked> Water vapour</label>
    </div>
    <div class="ctl-group">
      <span class="ctl-label">Wavelength → {m['target']:.0f} nm</span>
      <span class="seg" id="wlseg">
        <button data-wl="none">None</button>
        <button data-wl="angstrom">Simple (Ångström α={m['alpha']:.0f})</button>
        <button data-wl="molecular" class="on">Advanced (molecular)</button>
      </span>
    </div>
    <div class="ctl-group">
      <span class="ctl-label">Period</span>
      <select id="r0"></select><span class="muted">→</span><select id="r1"></select>
    </div>
    <div class="ctl-group">
      <label class="chk"><input type="checkbox" id="nfilt" checked> Noise-floor filtering</label>
      <label class="chk"><input type="checkbox" id="logx"> log β axis</label>
    </div>
    <span class="muted">both panels identical · <b id="period">–</b> · <b id="nhours">–</b> paired hours</span>
  </div>

  <h2>Vertical profiles <span class="muted">— median (line) and inter-quartile range (band).
    Each curve stops where that instrument runs out of signal, not at a fixed height; the grid runs
    to 15 km (CHM15k 15.3, CL61 15.7, CL31 7.7 km native). Default view 0–4 km — zoom out for the
    rest.</span></h2>
  <div class="grid2">
    <div class="card"><p class="panel-title">L1 + v2 calibration <span class="pill l1">Kalman C<sub>L</sub></span></p>
      <p class="panel-sub">rcs_0 / C<sub>L</sub>(t) × 10<sup>6</sup></p><div id="p_l1"></div></div>
    <div class="card"><p class="panel-title">L2 as distributed <span class="pill l2">provider constant</span></p>
      <p class="panel-sub">attenuated_backscatter_0, no re-calibration</p><div id="p_l2"></div></div>
  </div>

  <h2>Differences vs {ch[data['iref']]['label']}
    <span class="muted">— both sources on one axis: <b>L1 solid</b>, <b>L2 dashed</b>, same colour
    per instrument. Left: the median profile of the hourly ratios. Right: the distribution every
    paired sample in the {m['zmin']:.0f}–{m['zmax']:.0f} m band falls into, whose median is the
    figure tabulated below.</span></h2>
  <div class="grid2">
    <div class="card"><div id="d_diff"></div></div>
    <div class="card"><div id="p_hist"></div></div>
  </div>

  <h2>L1 vs L2 <span class="muted">— what the re-calibration changes, per instrument</span></h2>
  <div class="card"><div id="d_l1l2"></div></div>
  <p class="note">Positive = the v2-calibrated L1 profile is higher than the L2 product at that
    altitude. CL31 is the extreme case: its L2 constant is the uncalibrated 1e8 default.</p>

  <h2>Agreement over {m['zmin']:.0f}–{m['zmax']:.0f} m a.g.l.</h2>
  <div class="card">
    <table class="stats" id="stats">
      <thead><tr><th>Source</th><th>Instrument</th><th>Median rel. bias</th><th>log r</th>
        <th>Linear rel. bias</th><th>N pairs</th></tr></thead>
      <tbody></tbody>
    </table>
    <p class="note">Median relative bias is the robust pair statistic (median of the per-sample
      ratio); log r is the Pearson correlation in log space. Both are computed against
      {ch[data['iref']]['label']} over the shaded band.</p>
  </div>

  <h2>Calibration coefficients over time</h2>
  {cal_divs}
  <p class="note">v2.0 per-night constants with their uncertainty, the operational Kalman random-walk
    best estimate (red, ±1σ band), and the constant actually applied in the L2 product (black).
    The shaded span is the period the profile panels above are built from; the v2 series itself runs
    to {max(c['points'][-1]['date'] for c in data['calib'].values())}.</p>

  <div class="foot">
    Built from <code>presentation/build_l1_l2_dashboard.py</code> +
    <code>render_l1_l2_dashboard.py</code> · WV and molecular profiles resolved <b>per day</b>
    (<code>D:/CAMS_daily</code> first, monthly 0.4° L137 archive as fallback — Jun–Aug 2026 has no
    monthly CAMS at all) · hourly medians, ≥30 min coverage, cloud/fog/quality screened ±15 min.
  </div>
</div>
<script>const D = {json.dumps(data)};</script>
<script>{JS}</script>
</body></html>
"""


def main():
    data = json.loads((DATA / "data.json").read_text(encoding="utf-8"))
    out = PAGE
    out.write_text(html(data), encoding="utf-8")
    print(f"-> {out}  ({out.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
