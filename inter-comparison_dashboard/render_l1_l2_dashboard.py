# -*- coding: utf-8 -*-
"""L1-vs-L2 inter-comparison dashboard — render stage (per SITE, see sites.py).

Turns the JSON from build_l1_l2_dashboard.py into ONE self-contained HTML page (plotly.min.js and
the data are inlined, so the file opens anywhere with no server and no network). All the toggles are
client-side over precomputed arrays — nothing is recomputed in the browser. Every L2-dependent
element (right-hand panel, L1-vs-L2 section, L2 stats rows, "Applied in L2" traces) is gated on the
payload's meta.sources, so an L1-only site renders a single-panel page.

Run:  python inter-comparison_dashboard/render_l1_l2_dashboard.py [payerne|amsterdam|lindenberg]
Out:  inter-comparison_dashboard/index_<site>.html  (payerne also writes the historic index.html;
      data read from C:/DATA/.../inter-comparison_dashboard/data_<site>.json)
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")                # validation.paper, monitoring
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling modules (folder name is not a valid package name)
from plotly.offline import get_plotlyjs

import sites

DATA = sites.DATA_ROOT       # inputs (data_<site>.json)
HERE = Path(__file__).resolve().parent

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
const SRC = M.sources, HAS_L2 = SRC.indexOf('L2') >= 0;
const VN = M.variants, VSHORT = M.variant_short || {};
const CS = D.calib_series || {};
const vname = v => VSHORT[v] || v;
const state = { cal:VN[0], wv:true, wl:'molecular', logx:false, filter:true, r0:0, r1:MONTHS.length - 1 };
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

/* ---------- relative-difference profiles: every available source on one axis ---------- */
function diffTraces() {
  const out = [];
  const specs = HAS_L2 ? [['L2', 1.5, 'dash'], ['L1', 2.6, 'solid']] : [['L1', 2.6, 'solid']];
  specs.forEach(([src, w, dash]) => {
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
  const specs = HAS_L2 ? [['L2', 1.5, 'dash'], ['L1', 2.4, 'solid']] : [['L1', 2.4, 'solid']];
  specs.forEach(([src, w, dash]) => {
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
const VAR_COLOR = { 'v2.2':'#2f9e44', 'v2.2dark':'#7048e8' };
const VAR_PALETTE = ['#2f9e44', '#7048e8', '#e8590c', '#0ca678'];
const VAR_SYMBOL = ['triangle-up', 'diamond', 'square', 'cross'];

function calibFigure(ident, chan) {
  chan = chan || ident;
  const c = (CS[VN[0]] || {})[chan];           /* the base variant's series */
  if (!c) return null;
  const traces = [];
  const op = (D.l2_applied || {})[ident];
  if (op) traces.push({ x:op.date, y:op.value, mode:'lines', name:'Applied in L2',
    line:{color:'#111', width:1.3},
    hovertemplate:'%{x|%Y-%m-%d}<br>applied in L2 = %{y:.3e}<extra></extra>' });
  const k = c.kalman;
  if (k.date.length) {
    traces.push({ x:k.date.concat(k.date.slice().reverse()),
      y:k.value.map((v, i) => v + k.std[i]).concat(k.value.map((v, i) => v - k.std[i]).reverse()),
      fill:'toself', fillcolor:'rgba(214,39,40,0.12)', line:{width:0}, hoverinfo:'skip',
      showlegend:false });
    traces.push({ x:k.date, y:k.value, mode:'lines', name:vname(VN[0]) + ' Kalman estimate',
      line:{color:'#d62728', width:2},
      hovertemplate:'%{x|%Y-%m-%d}<br>Kalman = %{y:.3e}<extra></extra>' });
  }
  /* the other variants alongside, only where they actually differ — the cloud-calibrated channels
     are CARRIED OVER identically between variants, so drawing them twice would just duplicate. */
  VN.slice(1).forEach((vn, i) => {
    const c2 = (CS[vn] || {})[chan];
    if (!c2 || !c2.kalman || c2.key !== c.key || c2.created === c.created) return;
    const col = VAR_COLOR[vn] || VAR_PALETTE[i % VAR_PALETTE.length];
    const k2 = c2.kalman;
    if (k2.date.length) traces.push({ x:k2.date, y:k2.value, mode:'lines',
      name:vname(vn) + ' Kalman estimate', line:{color:col, width:2},
      hovertemplate:'%{x|%Y-%m-%d}<br>' + vname(vn) + ' Kalman = %{y:.3e}<extra></extra>' });
    traces.push({ x:c2.points.map(p => p.date), y:c2.points.map(p => p.value), mode:'markers',
      name:vname(vn) + ' — ' + (c2.points.length ? c2.points[0].method : 'Rayleigh'),
      marker:{ size:5, color:col, opacity:0.55, symbol:VAR_SYMBOL[i % VAR_SYMBOL.length] },
      hovertemplate:'%{x|%Y-%m-%d}<br>' + vname(vn) + ' C_L = %{y:.3e}<extra></extra>' });
  });
  const byM = {};
  c.points.forEach(p => (byM[p.method] = byM[p.method] || []).push(p));
  Object.keys(byM).sort().forEach(m => {
    const pts = byM[m];
    traces.push({ x:pts.map(p => p.date), y:pts.map(p => p.value), mode:'markers',
      name:vname(VN[0]) + ' — ' + m,
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

/* ---------- time-height curtains (canonical combo, drawn ONCE — static) ---------- */
function drawPcolor() {
  const P = D.pcolor;
  if (!P) return;
  CH.forEach((c, k) => {
    const el = document.getElementById('pc_' + (c.chan || c.ident));
    if (!el || !P.ch[k]) return;
    const traces = [
      /* qf-masked-only backdrop: what the cloud/fog screening removed, in grey */
      { type:'heatmap', x:P.hours, y:P.alt, z:P.ch[k].disp, zmin:-2, zmax:1,
        colorscale:[[0, '#f1f3f5'], [1, '#868e96']], showscale:false, hoverinfo:'skip' },
      { type:'heatmap', x:P.hours, y:P.alt, z:P.ch[k].scr, zmin:-2, zmax:1,
        colorscale:'Viridis',
        colorbar:{ title:{text:'log₁₀ β<sub>att</sub>'}, thickness:14, len:0.92 },
        hovertemplate:'%{x}<br>%{y:.0f} m<br>log₁₀ β = %{z:.2f}<extra>' + c.label + '</extra>' }
    ];
    const lay = JSON.parse(JSON.stringify(BASE));
    lay.height = 320;
    lay.showlegend = false;
    lay.margin.b = 30;
    lay.title = { text:c.label + ' — attenuated backscatter (colour = clear-sky screened, ' +
                       'grey = cloud/fog-excluded hours)', font:{size:13} };
    lay.yaxis = { title:'Altitude a.g.l. [m]' };
    lay.xaxis = { title:'' };
    Plotly.react(el, traces, lay, PCFG);
  });
}

/* ---------- statistics table ---------- */
function statsTable() {
  const rows = [];
  SRC.forEach(src => {
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
  const dTitle = 'Difference vs ' + CH[IREF].label + (HAS_L2 ? ' — L1 solid, L2 dashed' : '');
  Plotly.react('d_diff', diffTraces(), diffLayout(dTitle, [-100, 100]), PCFG);
  Plotly.react('p_hist', histTraces(), histLayout(), PCFG);
  if (HAS_L2) {
    Plotly.react('p_l2', profileTraces('L2'), profileLayout('L2 as distributed'), PCFG);
    const ll = diffLayout('L1 + v2 calibration relative to L2 — same instrument, same hour', [-100, 400]);
    ll.height = 500;
    Plotly.react('d_l1l2', l1l2Traces(), ll, PCFG);
  }
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
  /* WV / wavelength controls are absent on all-1064 sites (combos collapsed to variants only) */
  const wv = document.getElementById('wv');
  if (wv) wv.addEventListener('change', e => { state.wv = e.target.checked; draw(); });
  /* Calibration variant: it changes the CONSTANT the L1 panel divides by (v2.2 moves the Rayleigh
     retrieval only, so the cloud-calibrated channels stay identical between variants). */
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

fillMonths(); bind(); draw(); drawPcolor();
window.addEventListener('resize', () => ['p_l1','p_l2','d_diff','p_hist','d_l1l2']
  .concat(CH.map(c => 'pc_' + (c.chan || c.ident)))
  .forEach(id => { const el = document.getElementById(id);
                   if (el && el.data) Plotly.Plots.resize(el); }));
"""


def html(data):
    m = data["meta"]
    ch = data["channels"]
    has_l2 = "L2" in m["sources"]
    variants = m["variants"]
    vlabels = m.get("variant_labels", {})

    warns = "\n".join(f'  <div class="warn">{w}</div>' for w in m.get("warnings", []))
    cal_buttons = "\n".join(
        '        <button data-cal="%s"%s>%s</button>'
        % (v, ' class="on"' if i == 0 else "", vlabels.get(v, v))
        for i, v in enumerate(variants))
    cal_divs = "\n".join(
        f'  <div class="card" style="margin-bottom:12px"><div id="cal_{c.get("chan", c["ident"])}"></div></div>'
        for c in ch)
    pc_divs = "\n".join(
        f'  <div class="card" style="margin-bottom:12px"><div id="pc_{c.get("chan", c["ident"])}"></div></div>'
        for c in ch)

    # correction controls: hidden entirely on all-1064 sites (both toggles are exact no-ops there)
    corr_controls = "" if m.get("collapse") else f"""
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
    </div>"""

    if has_l2:
        profile_cards = """<div class="grid2">
    <div class="card"><p class="panel-title">L1 + v2 calibration <span class="pill l1">Kalman C<sub>L</sub></span></p>
      <p class="panel-sub">rcs_0 / C<sub>L</sub>(t) × 10<sup>6</sup></p><div id="p_l1"></div></div>
    <div class="card"><p class="panel-title">L2 as distributed <span class="pill l2">provider constant</span></p>
      <p class="panel-sub">attenuated_backscatter_0, no re-calibration</p><div id="p_l2"></div></div>
  </div>"""
    else:
        profile_cards = """<div class="card"><p class="panel-title">L1 + v2 calibration <span class="pill l1">Kalman C<sub>L</sub></span></p>
      <p class="panel-sub">rcs_0 / C<sub>L</sub>(t) × 10<sup>6</sup></p><div id="p_l1"></div></div>"""

    profiles_note = m.get("profiles_note") or (
        "Each curve stops where that instrument runs out of signal, not at a fixed height; the "
        "grid runs to 15 km. Default view 0–4 km — zoom out for the rest.")
    diff_note = ("both sources on one axis: <b>L1 solid</b>, <b>L2 dashed</b>, same colour per "
                 "instrument. " if has_l2 else "")

    l1l2_section = "" if not has_l2 else """
  <h2>L1 vs L2 <span class="muted">— what the re-calibration changes, per instrument</span></h2>
  <div class="card"><div id="d_l1l2"></div></div>
  <p class="note">Positive = the v2-calibrated L1 profile is higher than the L2 product at that
    altitude. CL31 is the extreme case: its L2 constant is the uncalibrated 1e8 default.</p>
"""

    base_series = data.get("calib_series", {}).get(variants[0], {})
    last_pt = max((c["points"][-1]["date"] for c in base_series.values() if c.get("points")),
                  default=m["end"])
    applied_note = (", and the constant actually applied in the L2 product (black)"
                    if has_l2 else "")

    pcolor_note = ("" if data.get("pcolor") else
                   "<p class=\"note\">No pcolor block in this payload — rebuild the data.</p>")
    tab_title = f"{m['station']} — L1 + v2 calibration" + (" vs L2" if has_l2 else "")

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{tab_title}</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<div class="topbar">
  <span class="brand">E-PROFILE ALC calibration</span>
  <span class="spacer"></span>
  <span class="asof">{m['station']} · {m['wmo']} · {m['start']} → {m['end']}</span>
</div>
<div class="container">
  <h1>{m['title']}</h1>
  <p class="subtitle">{m['subtitle']}</p>

{warns}

  <div class="card controls">
    <div class="ctl-group">
      <span class="ctl-label">Calibration</span>
      <span class="seg" id="calseg">
{cal_buttons}
      </span>
    </div>{corr_controls}
    <div class="ctl-group">
      <span class="ctl-label">Period</span>
      <select id="r0"></select><span class="muted">→</span><select id="r1"></select>
    </div>
    <div class="ctl-group">
      <label class="chk"><input type="checkbox" id="nfilt" checked> Noise-floor filtering</label>
      <label class="chk"><input type="checkbox" id="logx"> log β axis</label>
    </div>
    <span class="muted">{'both panels identical' if has_l2 else 'strictly paired hours'} · <b id="period">–</b> · <b id="nhours">–</b> paired hours</span>
  </div>

  <h2>Vertical profiles <span class="muted">— median (line) and inter-quartile range (band).
    {profiles_note}</span></h2>
  {profile_cards}

  <h2>Differences vs {ch[data['iref']]['label']}
    <span class="muted">— {diff_note}Left: the median profile of the hourly ratios. Right: the
    distribution every paired sample in the {m['zmin']:.0f}–{m['zmax']:.0f} m band falls into,
    whose median is the figure tabulated below.</span></h2>
  <div class="grid2">
    <div class="card"><div id="d_diff"></div></div>
    <div class="card"><div id="p_hist"></div></div>
  </div>
{l1l2_section}
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
  <p class="note">{variants[0]} per-night constants with their uncertainty, the operational Kalman
    random-walk best estimate (red, ±1σ band){applied_note}.
    The shaded span is the period the profile panels above are built from; the series itself runs
    to {last_pt}.</p>

  <h2>Time–height curtains <span class="muted">— calibrated L1 attenuated backscatter,
    log₁₀ colour scale −2…1 (canonical view: {variants[0]}, WV on, molecular wavelength).
    Colour = the clear-sky screened stream the whole page is computed from; grey = hours the
    cloud/fog screening excluded. Every 2nd gate up to 8 km; static — the controls above do not
    redraw these.</span></h2>
{pc_divs}{pcolor_note}

  <div class="foot">
    Built from <code>inter-comparison_dashboard/build_l1_l2_dashboard.py</code> +
    <code>render_l1_l2_dashboard.py</code> (site table <code>sites.py</code>) · WV and molecular
    profiles resolved <b>per day</b> (<code>D:/CAMS_daily</code> first, monthly 0.4° L137 archive
    as fallback) · hourly medians, ≥30 min coverage, cloud/fog/quality screened ±15 min.
  </div>
</div>
<script>const D = {json.dumps(data)};</script>
<script>{JS}</script>
</body></html>
"""


def main(site_key=None):
    if site_key is None:
        args = [a for a in sys.argv[1:] if not a.startswith("--")]
        site_key = args[0] if args else "payerne"
    sites.get_site(site_key)                                     # validate the key early
    src = DATA / f"data_{site_key}.json"
    if not src.exists() and site_key == "payerne":
        src = DATA / "data.json"                                 # pre-generalisation payload name
    data = json.loads(src.read_text(encoding="utf-8"))
    page = html(data)
    out = HERE / f"index_{site_key}.html"
    out.write_text(page, encoding="utf-8")
    print(f"-> {out}  ({out.stat().st_size/1e6:.1f} MB)")
    if site_key == "payerne":                                    # historic name, same content
        (HERE / "index.html").write_text(page, encoding="utf-8")
        print(f"-> {HERE / 'index.html'}  (backward-compat copy)")


if __name__ == "__main__":
    main()
