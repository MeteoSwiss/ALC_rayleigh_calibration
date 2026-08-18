# -*- coding: utf-8 -*-
"""Inter-comparison dashboard v3 — render stage: ONE self-contained page for ALL sites.

Reads v3_<site>.json (build_v3.py) and writes inter-comparison_dashboard/index.html.  Everything
the controls do is arithmetic in the browser over the per-channel ingredients described in
build_v3.py — there is no precomputed combination in the payload, so a site switcher, a
per-instrument method/variant choice and a free reference channel all cost nothing.

Run:  python inter-comparison_dashboard/render_v3.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plotly.offline import get_plotlyjs

import sites
import variants_v3 as V3

DATA = sites.DATA_ROOT
HERE = Path(__file__).resolve().parent
SITE_ORDER = ["payerne", "amsterdam", "lindenberg"]

CSS = """
:root { --bg:#f7f8fa; --card:#fff; --ink:#1c2330; --muted:#6b7280; --line:#e3e6ea;
        --brand:#0b3d61; --accent:#1f77b4; --off:#f2f4f7; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:14px/1.45 -apple-system,"Segoe UI",Roboto,Arial,sans-serif; }
.topbar { display:flex; align-items:center; gap:14px; background:var(--brand); color:#fff;
          padding:8px 18px; position:sticky; top:0; z-index:60; }
.topbar .brand { font-weight:600; }
.topbar .spacer { flex:1; }
.topbar .asof { color:#d6e6f2; font-size:13px; }
.siteseg { display:inline-flex; border:1px solid rgba(255,255,255,.35); border-radius:8px;
           overflow:hidden; }
.siteseg button { font:inherit; font-size:13px; padding:5px 14px; background:transparent;
                  color:#dbeafe; border:0; border-right:1px solid rgba(255,255,255,.25);
                  cursor:pointer; }
.siteseg button:last-child { border-right:0; }
.siteseg button.on { background:#fff; color:var(--brand); font-weight:600; }
.container { max-width:1600px; margin:0 auto; padding:16px 18px 40px; }
h1 { font-size:21px; margin:6px 0 4px; }
h2 { font-size:17px; margin:26px 0 10px; }
.subtitle { margin:0 0 12px; font-size:15px; color:#33405a; font-weight:500; }
.muted { color:var(--muted); font-weight:400; font-size:13px; }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
.panel-title { font-weight:600; font-size:14px; margin:0 0 2px; }
.panel-sub { font-size:12px; color:var(--muted); margin:0 0 6px; }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
.pill { display:inline-block; font-size:11px; padding:1px 9px; border-radius:999px;
        border:1px solid var(--line); color:var(--muted); margin-left:6px; }
.pill.l1 { background:#eef4fb; border-color:#cfe0f0; color:#17537f; }
.pill.l2 { background:#f4f1ee; border-color:#e6ddd3; color:#7a5a37; }
.pill.ref { background:#fff3bf; border-color:#ffe066; color:#7a5a00; font-weight:600; }

/* ---- control block ---- */
.ctlwrap { position:sticky; top:40px; z-index:50; margin:12px 0 16px;
           box-shadow:0 2px 12px rgba(0,0,0,.06); }
/* wide tables scroll inside their own container — the page never scrolls horizontally */
.tscroll { overflow-x:auto; }
table.instr { border-collapse:collapse; width:100%; font-size:13px; }
table.instr th { color:var(--muted); font-weight:600; font-size:11.5px; text-align:left;
                 padding:2px 8px; text-transform:uppercase; letter-spacing:.03em; }
table.instr td { padding:4px 8px; border-top:1px solid var(--line); vertical-align:middle; }
table.instr tr:first-child td { border-top:0; }
.swatch { display:inline-block; width:11px; height:11px; border-radius:3px; margin-right:7px;
          vertical-align:-1px; }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:7px; overflow:hidden; }
.seg button { font:inherit; font-size:12.5px; padding:4px 11px; background:#fff; color:var(--ink);
              border:0; border-right:1px solid var(--line); cursor:pointer; white-space:nowrap; }
.seg button:last-child { border-right:0; }
.seg button.on { background:var(--accent); color:#fff; font-weight:600; }
.seg button:disabled { background:var(--off); color:#b6bcc6; cursor:not-allowed; }
select { font:inherit; font-size:12.5px; padding:4px 7px; border:1px solid var(--line);
         border-radius:6px; background:#fff; color:var(--ink); cursor:pointer; max-width:340px; }
select:focus { outline:none; border-color:var(--accent); }
select option:disabled { color:#b6bcc6; }
.nights { font-variant-numeric:tabular-nums; font-size:12.5px; color:var(--muted); }
.nights b { color:var(--ink); }
.nights.thin b { color:#c92a2a; }
.chk { display:inline-flex; align-items:center; gap:7px; cursor:pointer; font-size:13px; }
.chk input { width:15px; height:15px; accent-color:var(--accent); cursor:pointer; }
.ctl-row { display:flex; flex-wrap:wrap; align-items:center; gap:20px; padding:8px 8px 2px;
           border-top:1px solid var(--line); margin-top:6px; }
.ctl-label { font-weight:600; font-size:12.5px; color:#33405a; }
.ctl-group { display:flex; align-items:center; gap:9px; }
.ctl-group.off { opacity:.42; }
.warnbox { background:#fff8e6; border:1px solid #f0dca8; border-left:4px solid #e0a800;
        border-radius:8px; padding:9px 13px; margin:9px 0; font-size:12.8px; color:#4a3a10; }
.warnbox b { color:#3a2d08; }
.err { background:#fff5f5; border:1px solid #ffc9c9; border-left:4px solid #e03131;
       border-radius:8px; padding:9px 13px; margin:9px 0; font-size:13px; color:#7a1010; }
table.stats { border-collapse:collapse; width:100%; font-size:13px; }
table.stats th, table.stats td { border-bottom:1px solid var(--line); padding:6px 10px;
                                 text-align:right; }
table.stats th:first-child, table.stats td:first-child,
table.stats th:nth-child(2), table.stats td:nth-child(2) { text-align:left; }
table.stats thead th { color:var(--muted); font-weight:600; font-size:12px; }
table.stats tbody tr:hover { background:#f6f9fc; }
table.ladder tbody tr.sel { background:#eef4fb; }
table.ladder tbody tr.l2row { background:#faf8f5; color:#6b5d4a; }
table.ladder tbody tr.l2row td { border-bottom-width:2px; }
table.ladder tbody tr.sel td { font-weight:600; box-shadow:inset 0 0 0 9999px rgba(31,119,180,.05); }
table.ladder tbody tr.sel td:first-child { box-shadow:inset 3px 0 0 var(--accent),
                                            inset 0 0 0 9999px rgba(31,119,180,.05); }
.warnflag { color:#c92a2a; font-weight:700; margin-left:5px; cursor:help; }
.warnflag.ok { color:#1971c2; font-weight:600; }
.num { font-variant-numeric:tabular-nums; }
.good { color:#1a7f37; font-weight:600; }
.bad  { color:#c92a2a; font-weight:600; }
.note { font-size:12.5px; color:var(--muted); margin:8px 0 0; }
.foot { color:var(--muted); font-size:12px; padding:22px 18px; text-align:center; }
"""

JS = r"""
/* ===================== decoding ===================== */
function b64bytes(s) {
  const bin = atob(s), n = bin.length, out = new Uint8Array(n);
  for (let i = 0; i < n; i++) out[i] = bin.charCodeAt(i);
  return out;
}
/* log-quantised uint16 -> Float64Array, row-major (hour, gate).  Twin of build_v3.quantize. */
function decodeStream(enc) {
  const q = new Uint16Array(b64bytes(enc.b).buffer);
  const nz = enc.nz, lo = enc.lo, hi = enc.hi, out = new Float64Array(q.length);
  const ratio = new Float64Array(nz);
  for (let g = 0; g < nz; g++) ratio[g] = hi[g] / lo[g];
  for (let i = 0; i < q.length; i++) {
    const c = q[i];
    if (c === 0) { out[i] = NaN; continue; }
    const g = i % nz, sgn = c >= 32768 ? -1 : 1, k = c & 32767;
    out[i] = k === 0 ? 0 : sgn * lo[g] * Math.pow(ratio[g], (k - 1) / 32766);
  }
  return out;
}
function decodeF32(blk) { return blk ? new Float32Array(b64bytes(blk.b).buffer) : null; }

/* ===================== state ===================== */
/* boot state (operator-confirmed): WV comparison on the manufacturer spectrum, molecular
   wavelength, noise filter off, full period — per-instrument picks come from each site's
   default/default_dark in variants_v3.SITE_V3 */
const state = { site:null, wv:'ctor', wl:'molecular', filter:false, logx:false, r0:0, r1:0,
                iref:0, pick:{} };
const CACHE = {};
const PCFG = { displaylogo:false, responsive:true };
const BASE = { template:'plotly_white', margin:{l:64,r:16,t:34,b:46}, font:{size:12}, height:470,
               hovermode:'closest', legend:{orientation:'h', y:-0.16, x:0.5, xanchor:'center',
               yanchor:'top', font:{size:11}, traceorder:'normal'} };
const ZTOP = 4000;
const S = () => SITES[state.site];
const rgba = (hex, a) => { const n = parseInt(hex.slice(1), 16);
  return `rgba(${n>>16&255},${n>>8&255},${n&255},${a})`; };
const MONTH_LABEL = ym => { const [y, m] = ym.split('-').map(Number);
  return ['janv','févr','mars','avr','mai','juin','juil','août','sept','oct','nov','déc'][m-1]
         + ' ' + y; };

/* decoded arrays, per site, built once */
function siteCache() {
  const P = S();
  if (CACHE[state.site]) return CACHE[state.site];
  const c = { raw:{}, corr:{} };
  for (const k in P.streams) c.raw[k] = decodeStream(P.streams[k]);
  for (const k in P.corr) {
    const q = P.corr[k];
    c.corr[k] = { lam:q.lam, f:q.f, nz:P.z.length,
                  wv:decodeF32(q.wv), wv2:decodeF32(q.wv2),
                  bml:decodeF32(q.bml), bmt:decodeF32(q.bmt) };
  }
  CACHE[state.site] = c;
  return c;
}

/* ===================== the transform ===================== */
/* out(h,z) = molecular/Angstrom( WV( (raw - dark)/C_L(t_h) * 1e6 ) ) — the numpy twin lives in
   build_v3.transform, and check_v3.py proves the two give the v2 statistics back. */
function computeChannel(k, src) {
  const inst = S().instruments[k];
  return computeSel(inst.ident, src, state.pick[inst.ident]);
}
/* `sel` = {method, variant, dark} — explicit so the summary table can sweep the whole ladder. */
function computeSel(id, src, sel) {
  const P = S(), C = siteCache();
  const nz = P.z.length, nh = P.hours.length;
  const raw = C.raw[src + '|' + id];
  if (!raw) return null;
  const cr = C.corr[P.corr_of[src + '|' + id]];
  const rec = P.calib[id + '|' + sel.method + '|' + sel.variant];
  if (src === 'L1' && !(rec && rec.ok)) return null;
  const cl = src === 'L1' ? rec.c : null;
  /* measured hood dark b(z), in rcs_0 units, subtracted BEFORE the calibration division — the
     per-hour band values are stored per gate, so this is exact, not a band-mean approximation. */
  const dark = (src === 'L1' && sel.dark && P.dark[id]) ? P.dark[id] : null;
  const wvArr = state.wv === 'none' ? null
              : (state.wv === 'ctor' && cr.wv2 ? cr.wv2 : cr.wv);
  const needWL = cr.lam !== null && Math.abs(cr.lam - P.target) >= 1.0;
  const useMol = needWL && state.wl === 'molecular' && cr.bmt;
  const useAng = needWL && state.wl === 'angstrom';
  const out = new Float64Array(nh * nz);
  for (let h = 0; h < nh; h++) {
    const d = P.hour_day[h], o = h * nz, dof = d * nz;
    const scale = src === 'L1' ? 1e6 / cl[h] : 1.0;
    for (let g = 0; g < nz; g++) {
      let v = raw[o + g];
      if (dark) v -= dark[g];
      v *= scale;
      if (wvArr) v /= wvArr[dof + g];
      if (useMol) v = cr.bmt[dof + g] + (v - cr.bml[dof + g]) * cr.f;
      else if (useAng) v *= cr.f;
      out[o + g] = v;
    }
  }
  return out;
}

/* ===================== per-view derived quantities ===================== */
function buildView() {
  const P = S(), nz = P.z.length, nh = P.hours.length, [b0, b1] = P.band;
  const V = { nz, nh, b0, b1, ch:[], ok:[], sources:P.sources, rows:null, n_hours:0 };
  P.sources.forEach(src => { V[src] = P.instruments.map((_, k) => computeChannel(k, src)); });
  V.ok = P.instruments.map((_, k) => P.sources.every(s => V[s][k] !== null));
  /* paired hours: every ACTIVE instrument, in every source, has a finite gate in the band */
  const rows = [];
  for (let h = 0; h < nh; h++) {
    const m = P.hour_month[h];
    if (m < state.r0 || m > state.r1) continue;
    let good = true;
    for (let k = 0; k < P.instruments.length && good; k++) {
      if (!V.ok[k]) continue;
      for (const s of P.sources) {
        let any = false;
        const A = V[s][k], o = h * nz;
        for (let g = b0; g <= b1; g++) if (isFinite(A[o + g])) { any = true; break; }
        if (!any) { good = false; break; }
      }
    }
    if (good) rows.push(h);
  }
  V.rows = rows;
  V.n_hours = rows.length;
  return V;
}

function quantiles(vals) {          /* vals: sorted Float64Array slice length n */
  const n = vals.length;
  if (!n) return [NaN, NaN, NaN];
  const at = p => { const i = p * (n - 1), lo = Math.floor(i), hi = Math.ceil(i);
                    return vals[lo] + (vals[hi] - vals[lo]) * (i - lo); };
  return [at(0.5), at(0.25), at(0.75)];
}

function profileOf(V, src, k) {
  const P = S(), nz = V.nz, A = V[src][k], rows = V.rows;
  const med = new Array(nz), q1 = new Array(nz), q3 = new Array(nz), n = new Array(nz);
  const buf = new Float64Array(rows.length);
  for (let g = 0; g < nz; g++) {
    let c = 0;
    for (let i = 0; i < rows.length; i++) { const v = A[rows[i] * nz + g]; if (isFinite(v)) buf[c++] = v; }
    const sub = buf.subarray(0, c).slice().sort((a, b) => a - b);
    const [m, a, b] = quantiles(sub);
    med[g] = m; q1[g] = a; q3[g] = b; n[g] = c;
  }
  return { med, q1, q3, n, keep:keepMask(med, n, P.z, rows.length) };
}

/* noise-floor mask — port of build_l1_l2_dashboard._keep_mask (coverage, sign, sustained rise) */
function keepMask(med, n, z, nprof) {
  const nz = med.length, minN = Math.max(10, 0.20 * nprof), zmin = S().zmin;
  const good = med.map((v, i) => isFinite(v) && n[i] >= minN);
  let cut = nz;
  for (let i = 0; i < nz; i++)
    if (z[i] > zmin && (!good[i] || !(med[i] > 0))) { cut = i; break; }
  const idx = [];
  for (let i = 0; i < cut; i++) if (z[i] > zmin && good[i] && med[i] > 0) idx.push(i);
  if (idx.length) {
    /* lowest gate from which the median stays above 1.5x its running minimum ALL the way up:
       an elevated aerosol layer also lifts the median but comes back down, an r^2-amplified
       electronic offset never does. */
    const runmin = [];
    let m = Infinity;
    idx.forEach(i => { m = Math.min(m, med[i]); runmin.push(m); });
    let acc = true, first = -1;
    for (let j = idx.length - 1; j >= 0; j--) {
      acc = acc && (med[idx[j]] > 1.5 * runmin[j]);
      if (!acc) break;
      first = j;
    }
    if (first >= 0) cut = Math.min(cut, idx[first]);
  }
  return med.map((_, i) => good[i] && i < cut);
}

/* per-gate median of the per-hour relative difference (cur vs ref), in % */
function diffProfile(V, src, k, kref) {
  const nz = V.nz, A = V[src][k], B = V[src][kref], rows = V.rows;
  const out = new Array(nz).fill(null);
  const buf = new Float64Array(rows.length);
  for (let g = 0; g < nz; g++) {
    let c = 0;
    for (let i = 0; i < rows.length; i++) {
      const o = rows[i] * nz + g, a = A[o], b = B[o];
      if (isFinite(a) && isFinite(b) && b > 0) buf[c++] = (a - b) / b * 100;
    }
    if (c) { const sub = buf.subarray(0, c).slice().sort((x, y) => x - y);
             out[g] = quantiles(sub)[0]; }
  }
  return out;
}

/* the statistics table + the histogram, over the band, per sample — intercompare._stats */
const HIST_EDGES = (() => { const e = []; for (let v = -100; v <= 200.0001; v += 2.5) e.push(v);
                            return e; })();
function bandStats(V, src, k, kref) {
  return bandStatsAB(V[src][k], V[src][kref], V);
}
function bandStatsAB(A, B, V) {
  const nz = V.nz, rows = V.rows, [b0, b1] = [V.b0, V.b1];
  const a = [], b = [];
  for (let i = 0; i < rows.length; i++) {
    const o = rows[i] * nz;
    for (let g = b0; g <= b1; g++) {
      const x = A[o + g], y = B[o + g];
      if (isFinite(x) && isFinite(y)) { a.push(x); b.push(y); }
    }
  }
  const n = a.length;
  const hist = new Array(HIST_EDGES.length - 1).fill(0);
  if (n <= 2) return { n, medrel:null, mad:null, rlog:null, relbias:null, hist, histmed:null };
  let sd = 0, sb = 0;
  const rel = [];
  const la = [], lb = [];
  for (let i = 0; i < n; i++) {
    sd += a[i] - b[i]; sb += b[i];
    if (b[i] > 0) {
      const r = (a[i] - b[i]) / b[i] * 100;
      rel.push(r);
      const j = Math.floor((r + 100) / 2.5);
      if (j >= 0 && j < hist.length) hist[j]++;
    }
    if (a[i] > 0 && b[i] > 0) { la.push(Math.log10(a[i])); lb.push(Math.log10(b[i])); }
  }
  rel.sort((x, y) => x - y);
  const med = rel.length ? quantiles(rel)[0] : null;
  /* robust dispersion: 1.4826 x median absolute deviation of the same per-sample differences */
  let mad = null;
  if (rel.length) {
    const ad = rel.map(v => Math.abs(v - med)).sort((x, y) => x - y);
    mad = 1.4826 * quantiles(ad)[0];
  }
  let rlog = null;
  if (la.length > 2) {
    const ma = la.reduce((s, v) => s + v, 0) / la.length;
    const mb = lb.reduce((s, v) => s + v, 0) / lb.length;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < la.length; i++) { const dx = la[i] - ma, dy = lb[i] - mb;
                                          sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
    rlog = sxy / Math.sqrt(sxx * syy);
  }
  return { n, nrel:rel.length, medrel:med, mad, rlog,
           relbias: sb !== 0 ? 100 * (sd / n) / (sb / n) : null, hist, histmed:med };
}

/* ===================== plots ===================== */
function zoneShape() { const P = S();
  return { type:'rect', xref:'paper', yref:'y', x0:0, x1:1, y0:P.zmin, y1:P.zmax,
           fillcolor:'rgba(31,119,180,0.045)', line:{width:0}, layer:'below' }; }
const mask = (arr, keep) => state.filter ? arr.map((v, i) => keep[i] ? v : null) : arr;

function profileTraces(V, src) {
  const P = S(), out = [];
  P.instruments.forEach((c, k) => {
    if (!V.ok[k]) return;
    const p = profileOf(V, src, k);
    out.push({ x:mask(p.q1, p.keep), y:P.z, mode:'lines', line:{width:0}, showlegend:false,
               hoverinfo:'skip' });
    out.push({ x:mask(p.q3, p.keep), y:P.z, mode:'lines', line:{width:0}, fill:'tonextx',
               fillcolor:rgba(c.color, 0.13), showlegend:false, hoverinfo:'skip' });
    out.push({ x:mask(p.med, p.keep), y:P.z, mode:'lines', line:{color:c.color, width:2.2},
               name:chLabel(k),
               hovertemplate:'%{y:.0f} m<br>' + c.label + ' = %{x:.3f}<extra></extra>' });
  });
  return out;
}
function profileLayout(title) {
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.title = { text:title, font:{size:13} };
  /* Plotly v3 dropped the string shorthand for axis titles — they MUST be {text:...} objects or
     they are silently not drawn (found 2026-08-16: every axis looked titled in the source and
     none was on screen) */
  lay.xaxis = state.logx
    ? { title:{text:'β<sub>att</sub> [Mm⁻¹ sr⁻¹]'}, type:'log', range:[-2.3, 0.3],
        zeroline:false }
    : { title:{text:'β<sub>att</sub> [Mm⁻¹ sr⁻¹]'}, range:[0, 1.0], zeroline:false };
  lay.yaxis = { title:{text:'Altitude a.g.l. [m]'}, range:[0, ZTOP] };
  lay.shapes = [zoneShape()];
  return lay;
}
function diffLayout(title, xr) {
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.title = { text:title, font:{size:13} };
  lay.xaxis = { title:{text:'différence relative [%]'}, range:xr, zeroline:true,
                zerolinewidth:1.4, zerolinecolor:'#888' };
  lay.yaxis = { title:{text:'Altitude a.g.l. [m]'}, range:[0, ZTOP] };
  lay.shapes = [zoneShape()];
  return lay;
}
/* Variant naming is wavelength-aware: each type shows ITS OWN water-vapour model
   (CL31 909.7/6.0, CL51 910.0/3.4, CL61 910.74/1.0), and 1064 nm none at all. */
/* one variant id 'v2.2dark' covers two very different dark sources: the MEASURED covered-
   telescope b(z) at Payerne vs the clear-night ESTIMATED dark (demonstration-only) elsewhere —
   the label must say which one this site is showing */
const darkSuffix = v => (DARK_RUNS.indexOf(v) >= 0 && S().dark_kind === 'estimated')
  ? ' — dark ESTIMÉ ciel clair (démonstration)' : '';
const vLabel = (itype, v) =>
  ((VLBL_BY_TYPE[itype] || {})[v] || VARIANT_LABEL[v] || v) + darkSuffix(v);
const vShort = (itype, v) =>
  ((VSH_BY_TYPE[itype] || {})[v] || VSHORT[v] || v) +
  (darkSuffix(v) ? ' (dark estimé)' : '');
function chLabel(k) {
  const P = S(), inst = P.instruments[k], s = state.pick[inst.ident];
  return inst.label + ' — ' + (s.method === 'rayleigh' ? 'Ray.' : 'nuage') + ' / ' +
         vShort(inst.itype, s.variant);
}

function draw() {
  const P = S(), V = buildView();
  const kref = state.iref;
  const hasL2 = P.sources.indexOf('L2') >= 0;
  Plotly.react('p_l1', profileTraces(V, 'L1'),
               profileLayout('L1 + étalonnage v2 (Kalman)'), PCFG);
  if (hasL2) Plotly.react('p_l2', profileTraces(V, 'L2'),
                          profileLayout('L2 tel que distribué'), PCFG);

  /* differences vs the chosen reference */
  const dt = [];
  const specs = hasL2 ? [['L2', 1.5, 'dash'], ['L1', 2.6, 'solid']] : [['L1', 2.6, 'solid']];
  const keeps = {};
  specs.forEach(([src]) => { keeps[src] = P.instruments.map((_, k) =>
      V.ok[k] ? profileOf(V, src, k).keep : null); });
  specs.forEach(([src, w, dash]) => {
    P.instruments.forEach((c, k) => {
      if (k === kref || !V.ok[k] || !V.ok[kref]) return;
      const d = diffProfile(V, src, k, kref);
      const kp = keeps[src][k].map((v, i) => v && keeps[src][kref][i]);
      dt.push({ x:mask(d, kp), y:P.z, mode:'lines', line:{color:c.color, width:w, dash:dash},
                name:chLabel(k) + ' — ' + src,
                hovertemplate:'%{y:.0f} m<br>%{x:+.1f}%<extra>' + c.label + ' ' + src + '</extra>' });
    });
  });
  Plotly.react('d_diff', dt, diffLayout('Différence vs ' + chLabel(kref) +
               (hasL2 ? ' — L1 plein, L2 tireté' : ''), [-100, 100]), PCFG);

  /* histogram (the tabulated per-state agreement lives in the Tableau récapitulatif) */
  const mid = HIST_EDGES.slice(0, -1).map((e, i) => (e + HIST_EDGES[i + 1]) / 2);
  const ht = [];
  specs.forEach(([src, w, dash]) => {
    P.instruments.forEach((c, k) => {
      if (k === kref || !V.ok[k] || !V.ok[kref]) return;
      const s = bandStats(V, src, k, kref);
      /* normalise by ALL relative samples, not the in-range subset — a heavy-tailed trace
         (21% of the CL31 samples fall outside [-100,200)) must not be inflated vs the others */
      const tot = s.nrel || s.hist.reduce((a, b) => a + b, 0);
      if (tot) ht.push({ x:mid, y:s.hist.map(v => v / tot * 100), mode:'lines',
        line:{color:c.color, width:w, dash:dash, shape:'hvh'},
        name:chLabel(k) + ' — ' + src + ' (méd ' +
             (s.medrel === null ? '—' : (s.medrel > 0 ? '+' : '') + s.medrel.toFixed(1) + '%') + ')',
        hovertemplate:'%{x:+.1f}%<br>%{y:.2f}% des échantillons<extra></extra>' });
    });
  });
  const hl = JSON.parse(JSON.stringify(BASE));
  hl.title = { text:'Différence vs ' + chLabel(kref) + ' — distribution sur ' +
               P.zmin.toFixed(0) + '–' + P.zmax.toFixed(0) + ' m', font:{size:13} };
  hl.xaxis = { title:{text:'différence relative [%]'}, range:[-100, 150], zeroline:true,
               zerolinewidth:1.4, zerolinecolor:'#888' };
  hl.yaxis = { title:{text:'part des échantillons [%]'} };
  Plotly.react('p_hist', ht, hl, PCFG);

  /* L1 vs L2, same instrument */
  if (hasL2) {
    const tt = [];
    P.instruments.forEach((c, k) => {
      if (!V.ok[k]) return;
      const nz = V.nz, A = V.L1[k], B = V.L2[k], rws = V.rows;
      const out = new Array(nz).fill(null), buf = new Float64Array(rws.length);
      for (let g = 0; g < nz; g++) {
        let cc = 0;
        for (let i = 0; i < rws.length; i++) { const o = rws[i] * nz + g, x = A[o], y = B[o];
          if (isFinite(x) && isFinite(y) && y > 0) buf[cc++] = (x - y) / y * 100; }
        if (cc) { const sub = buf.subarray(0, cc).slice().sort((p, q) => p - q);
                  out[g] = quantiles(sub)[0]; }
      }
      const kp = keeps.L1[k].map((v, i) => v && keeps.L2[k][i]);
      tt.push({ x:mask(out, kp), y:P.z, mode:'lines', line:{color:c.color, width:2.2},
                name:chLabel(k),
                hovertemplate:'%{y:.0f} m<br>L1 %{x:+.1f}% vs L2<extra>' + c.label + '</extra>' });
    });
    const ll = diffLayout('L1 + étalonnage v2 relatif au L2 — même instrument, même heure',
                          [-100, 400]);
    ll.height = 500;
    Plotly.react('d_l1l2', tt, ll, PCFG);
  }

  document.getElementById('nhours').textContent = V.n_hours.toLocaleString('fr-FR');
  document.getElementById('period').textContent = state.r0 === state.r1
    ? MONTH_LABEL(P.months[state.r0])
    : MONTH_LABEL(P.months[state.r0]) + ' – ' + MONTH_LABEL(P.months[state.r1]);
  drawCalib();
  drawHopkin();
  drawLadder(V);
}

/* ---------- summary table: the WHOLE ladder at once ----------
   One row per (instrument x available variant of its selected method), so the effect of each
   hypothesis — WV spectrum, dark, pipeline version — is readable without clicking through states.
   Same per-hour arithmetic as everything else, just looped; the paired-hour set is variant
   independent (dividing by a positive constant cannot make a value non-finite), so every row is
   computed on exactly the same sample. */
function drawLadder(V) {
  const P = S(), kref = state.iref, hasL2 = P.sources.indexOf('L2') >= 0;
  const ref = V.ok[kref] ? V.L1[kref] : null;
  const rows = [];
  P.instruments.forEach((inst, k) => {
    const id = inst.ident, sel = state.pick[id];
    const list = (VARIANTS_FOR[inst.itype] || {})[sel.method] || [];
    list.forEach(vn => {
      const rec = P.calib[id + '|' + sel.method + '|' + vn];
      if (!rec || !rec.ok) return;
      /* each Rayleigh row is swept in ITS coherent dark pairing (dark on iff the variant is a
         dark-corrected run and the measured b(z) exists) — sweeping with the CURRENT selection's
         flag tabulated hybrid states ~1 pt off what clicking the variant displays (2026-08-16
         review).  Cloud rows keep the selection's profile dark: constants are dark-immune there
         and both profile states are coherent user choices. */
      const rdark = sel.method === 'rayleigh'
        ? (DARK_RUNS.indexOf(vn) >= 0 && !!P.dark[id]) : sel.dark;
      const A = computeSel(id, 'L1', { method:sel.method, variant:vn, dark:rdark });
      if (!A) return;
      let delta = '—', disp = '—';
      if (k === kref) {
        if (hasL2 && V.L2[k]) {
          const s = bandStatsAB(A, V.L2[k], V);
          delta = (s.medrel === null ? '—' : (s.medrel > 0 ? '+' : '') + s.medrel.toFixed(1) +
                   ' %<span class="muted"> vs L2</span>');
          disp = s.mad === null ? '—' : s.mad.toFixed(1) + ' %';
        } else { delta = '<span class="muted">réf.</span>'; }
      } else if (ref) {
        const s = bandStatsAB(A, ref, V);
        delta = s.medrel === null ? '—' : (s.medrel > 0 ? '+' : '') + s.medrel.toFixed(1) + ' %';
        disp = s.mad === null ? '—' : s.mad.toFixed(1) + ' %';
      }
      const hk = (P.hopkin && P.hopkin.panels[id + '|' + vn]) || null;
      const on = vn === sel.variant;
      const star = REFERENCE_VARIANTS.indexOf(vn) >= 0 ? ' ★' : '';
      rows.push(`<tr class="${on ? 'sel' : ''}">
        <td><span class="swatch" style="background:${inst.color}"></span>${inst.label}</td>
        <td>${vLabel(inst.itype, vn).replace(' ★', '')}${star}</td>
        <td class="num">${delta}</td>
        <td class="num">${disp}</td>
        <td class="num">${V.n_hours.toLocaleString('fr-FR')}</td>
        <td class="num">${rec.nights}</td>
        <td class="num">${hk ? (hk.slope > 0 ? '+' : '') + hk.slope.toFixed(1) + ' ± ' +
                              hk.se.toFixed(1) + ' %/km' : '<span class="muted">—</span>'}</td>
      </tr>`);
    });
    /* the distributed product as one more row per instrument, so the whole ladder is judged
       against what the network actually ships today (L2 vs the reference's L2, same band).
       The applied-constant provenance DIFFERS per instrument type; the tooltip appends the
       constants actually observed in the window, computed live so it cannot go stale. */
    if (hasL2 && V.ok[k] && V.L2 && V.L2[k]) {
      let delta = '<span class="muted">réf.</span>', disp = '—';
      if (k !== kref && V.ok[kref] && V.L2[kref]) {
        const s = bandStatsAB(V.L2[k], V.L2[kref], V);
        delta = s.medrel === null ? '—' : (s.medrel > 0 ? '+' : '') + s.medrel.toFixed(1) + ' %';
        disp = s.mad === null ? '—' : s.mad.toFixed(1) + ' %';
      }
      const prov = l2Provenance(inst);
      rows.push(`<tr class="l2row">
        <td><span class="swatch" style="background:${inst.color}"></span>${inst.label}</td>
        <td>L2 distribué <span class="pill l2" title="${prov.tip.replace(/"/g, '&quot;')}">${prov.pill}</span></td>
        <td class="num">${delta}</td>
        <td class="num">${disp}</td>
        <td class="num">${V.n_hours.toLocaleString('fr-FR')}</td>
        <td class="num">—</td>
        <td class="num"><span class="muted">—</span></td>
      </tr>`);
    }
  });
  document.querySelector('#ladder tbody').innerHTML = rows.join('') ||
    '<tr><td colspan="7" class="muted">aucune variante exploitable</td></tr>';
}

/* per-type provenance of the L2 product's applied constant, + the values actually observed in
   the paired window (computed from l2_applied, never hard-coded) */
function l2Provenance(inst) {
  const fmt = v => Math.abs(v) >= 1e4 ? (+v).toExponential(1).replace('e+', 'e') : (+v).toFixed(2);
  const d = (S().l2_applied || {})[inst.ident];
  let obs = '';
  if (d && d.value && d.value.length) {
    const cnt = {};
    d.value.forEach(v => { const k = fmt(v); cnt[k] = (cnt[k] || 0) + 1; });
    obs = ' Observé sur la fenêtre : ' +
          Object.keys(cnt).map(k => `${k} (${cnt[k]} j)`).join(', ') + '.';
  }
  /* Operator-stated provenance (2026-08-16, corrected same day): the distributed L2 carries the
     operational Rayleigh v1.0 calibration for the CHM15k AND — recently activated — the CL61;
     the CL31 stays on the uncalibrated default. The observed CL61 constants confirm it:
     1.00 (default) on 48 d then 2.1257 from mid-June 2026 = the v1.0 series going live.
     The appended observed values keep the tooltip honest if the files ever change. */
  const base = {
    CHM15k: { pill:'calibration Rayleigh v1.0 (opérationnelle)', short:'Rayleigh v1.0',
      tip:'Constante issue de la méthode Rayleigh opérationnelle E-PROF v1.0 (pré-v2).' },
    CL31: { pill:'non calibré (constante par défaut)', short:'constante par défaut',
      tip:'La chaîne opérationnelle n\'applique pas de calibration au CL31 — ' +
          'la constante par défaut explique l\'écart de cette ligne ; c\'est exactement le ' +
          'manque que l\'étalonnage nuage comble.' },
    CL61: { pill:'calibration Rayleigh v1.0 (opérationnelle, activée mi-2026)',
      short:'Rayleigh v1.0',
      tip:'Le CL61 est désormais calibré en opérationnel par la méthode Rayleigh v1.0 — ' +
          'la bascule est visible dans les constantes observées : 1,00 (défaut) avant, puis ' +
          'la valeur v1.0 depuis la mi-juin 2026.' },
  }[inst.itype] || { pill:'constante appliquée', short:'constante appliquée', tip:'' };
  return { pill: base.pill, short: base.short, tip: base.tip + obs };
}

/* The L2 profile card's pill: per-site summary of who ships WHICH constant (one static label
   used to claim 'constante du fournisseur' for everything — wrong since the CHM15k, then the
   CL61, went operationally Rayleigh-v1.0-calibrated). */
function fillL2Pill() {
  const el = document.getElementById('l2prov_pill');
  if (!el) return;
  const groups = {};
  S().instruments.forEach(i => {
    const p = l2Provenance(i);
    (groups[p.short] = groups[p.short] || []).push(i.itype);
  });
  el.textContent = Object.keys(groups)
    .map(sh => sh + ' : ' + [...new Set(groups[sh])].join(', ')).join(' · ');
  el.title = S().instruments.map(i => i.itype + ' — ' + l2Provenance(i).pill).join(' | ');
}

/* ---------- constants over time ---------- */
const VAR_COLOR = { 'v2.0':'#d62728', 'v2.2':'#2f9e44', 'v2.2dark':'#7048e8',
  'v2.2sansWV':'#e8590c', 'v2.2_l55s008':'#0b7285', 'v2.2_l55w10':'#ae3ec9',
  'v2.2_l55w01':'#f08c00', 'cloudWV':'#2f9e44', 'cloudNoWV':'#e8590c',
  'cloud_l55s008':'#0b7285', 'cloud_l55w10':'#ae3ec9', 'cloud_l55w01':'#f08c00',
  'cloud_cl31l910':'#5c940d', 'cloud_cl31wieg':'#a61e4d' };
/* EXCLUSIVE end of the month (first day of the next one), so the shaded stats window covers the
   last day's 24 h too — day-0 gave midnight of the last day and under-covered by one day */
function monthEnd(ym) { const [y, m] = ym.split('-').map(Number);
  return new Date(Date.UTC(y, m, 1)).toISOString().slice(0, 10); }

function drawCalib() {
  const P = S();
  P.instruments.forEach(inst => {
    const el = document.getElementById('cal_' + inst.ident);
    if (!el) return;
    const sel = state.pick[inst.ident], traces = [];
    const op = (P.l2_applied || {})[inst.ident];
    if (op) traces.push({ x:op.date, y:op.value, mode:'lines', name:'appliqué en L2',
      line:{color:'#111', width:1.3},
      hovertemplate:'%{x|%Y-%m-%d}<br>appliqué en L2 = %{y:.3e}<extra></extra>' });
    (((VARIANTS_FOR[inst.itype] || {})[sel.method]) || []).forEach(vn => {
      const rec = P.calib[inst.ident + '|' + sel.method + '|' + vn];
      if (!rec || !rec.ok) return;
      const on = vn === sel.variant, col = VAR_COLOR[vn] || '#888', nm = vShort(inst.itype, vn);
      traces.push({ x:rec.kal.d, y:rec.kal.v, mode:'lines',
        name:nm + (on ? ' — sélectionné' : ''),
        line:{color:col, width:on ? 2.6 : 1.2, dash:on ? 'solid' : 'dot'},
        opacity:on ? 1 : 0.75,
        hovertemplate:'%{x|%Y-%m-%d}<br>' + nm + ' = %{y:.3e}<extra></extra>' });
      if (on) traces.push({ x:rec.points.map(p => p.d), y:rec.points.map(p => p.v),
        mode:'markers', name:nm + ' — nuits',
        marker:{ size:5, color:col, opacity:0.7 },
        error_y:{ type:'data', array:rec.points.map(p => p.s || 0), visible:true,
                  thickness:0.6, width:0, color:'rgba(120,120,120,0.35)' },
        hovertemplate:'%{x|%Y-%m-%d}<br>C_L = %{y:.3e}<extra></extra>' });
    });
    const lay = JSON.parse(JSON.stringify(BASE));
    lay.height = 300;
    lay.title = { text:inst.itype + ' (' + inst.ident + ') — C_L, variantes « ' +
                  (sel.method === 'rayleigh' ? 'Rayleigh' : 'nuage') + ' »', font:{size:13} };
    /* the constant's unit is instrument-native (rcs_0 / beta_att; ~dimensionless for the CL61) */
    lay.yaxis = { title:{text:'C_L [unités rcs₀ de l’instrument]'}, exponentformat:'e' };
    lay.xaxis = { title:{text:'date (UTC)'} };
    lay.shapes = [{ type:'rect', xref:'x', yref:'paper', x0:P.months[state.r0] + '-01',
                    x1:monthEnd(P.months[state.r1]), y0:0, y1:1,
                    fillcolor:'rgba(31,119,180,0.10)', line:{width:0}, layer:'below' }];
    Plotly.react(el, traces, lay, PCFG);
  });
}

/* ---------- curtains (static, page default state) ---------- */
function drawPcolor() {
  const P = S(), C = P.pcolor;
  if (!C) return;
  C.ch.forEach(cc => {
    const el = document.getElementById('pc_' + cc.ident);
    const inst = P.instruments.find(i => i.ident === cc.ident);
    if (!el) return;
    const unpack = s => { const q = b64bytes(s), z = [];
      for (let y = 0; y < C.ny; y++) { const row = new Array(C.nx);
        for (let x = 0; x < C.nx; x++) { const v = q[y * C.nx + x];
          row[x] = v === 255 ? null : C.lo + v / 254 * (C.hi - C.lo); }
        z.push(row); }
      return z; };
    const traces = [
      { type:'heatmap', x:C.hours, y:C.alt, z:unpack(cc.disp), zmin:-2, zmax:1,
        colorscale:[[0, '#f1f3f5'], [1, '#868e96']], showscale:false, hoverinfo:'skip' },
      { type:'heatmap', x:C.hours, y:C.alt, z:unpack(cc.scr), zmin:-2, zmax:1,
        colorscale:'Viridis',
        colorbar:{ title:{text:'log₁₀ β<sub>att</sub> [Mm⁻¹ sr⁻¹]'}, thickness:14, len:0.92 },
        hovertemplate:'%{x}<br>%{y:.0f} m<br>log₁₀ β = %{z:.2f}<extra>' +
                      (inst ? inst.label : cc.ident) + '</extra>' }];
    const lay = JSON.parse(JSON.stringify(BASE));
    lay.height = 320; lay.showlegend = false; lay.margin.b = 30;
    lay.title = { text:(inst ? inst.label : cc.ident) +
      ' — rétrodiffusion atténuée (couleur = flux ciel clair, gris = heures écartées)',
      font:{size:13} };
    lay.yaxis = { title:{text:'Altitude a.g.l. [m]'} };
    lay.xaxis = { title:{text:'date (UTC)'} };
    lay.margin.b = 46;
    Plotly.react(el, traces, lay, PCFG);
  });
}

/* ---------- Hopkin per-profile CBH heatmap — follows the selected CLOUD variant ---------- */
function drawHopkin() {
  const P = S(), H = P.hopkin, host = document.getElementById('hkholder');
  const sec = document.getElementById('hksec');
  if (!H) { sec.style.display = 'none'; return; }
  const wanted = P.instruments.filter(i =>
    (METHOD_BY_TYPE[i.itype] || []).indexOf('cloud') >= 0);
  const cards = [];
  wanted.forEach(i => {
    const sel = state.pick[i.ident];
    const key = i.ident + '|' + (sel.method === 'cloud' ? sel.variant : 'cloudWV');
    if (H.panels[key]) cards.push([i, key, sel.method === 'cloud']);
  });
  sec.style.display = cards.length ? '' : 'none';
  if (!cards.length) return;
  /* side by side on one row (two columns as soon as there are two units, e.g. CL31 + CL61) */
  if (host.dataset.n !== String(cards.length)) {
    host.style.display = 'grid';
    host.style.gridTemplateColumns = cards.length > 1 ? '1fr 1fr' : '1fr';
    host.style.gap = '14px';
    host.innerHTML = cards.map((c, j) =>
      `<div class="card"><div id="hk_${j}"></div></div>`).join('');
  }
  host.dataset.n = String(cards.length);
  const half = cards.length > 1;
  cards.forEach(([inst, key, live], j) => {
    const p = H.panels[key], el = document.getElementById('hk_' + j);
    if (!el) return;
    const z = p.hist.map(row => row.map(v => v > 0 ? v : null));
    const xc = H.x.slice(0, -1).map((v, i) => (v + H.x[i + 1]) / 2);
    const yc = H.y.slice(0, -1).map((v, i) => (v + H.y[i + 1]) / 2);
    const traces = [
      { type:'heatmap', x:xc, y:yc, z, colorscale:'Jet', zmin:1,
        colorbar:{ title:{text:'profils'}, thickness:14, len:0.92 },
        hovertemplate:'1/C %{x:.0f} %<br>CBH %{y:.2f} km<br>%{z} profils<extra></extra>' },
      { x:p.bands.map(b => b[1]), y:p.bands.map(b => b[0]),
        error_x:{ type:'data', array:p.bands.map(b => b[2]), thickness:1, width:3, color:'#111' },
        mode:'markers', marker:{ size:6, color:'#fff', line:{color:'#111', width:1.2} },
        name:'moyenne ± écart-type par bande de 250 m',
        hovertemplate:'CBH %{y:.2f} km : %{x:.0f} %<extra></extra>' }];
    const lay = JSON.parse(JSON.stringify(BASE));
    lay.height = 430;
    lay.showlegend = false;
    /* half-width cards: the slope moves down into an annotation so the title stays readable */
    lay.title = { text:inst.label + ' — ' + vLabel(inst.itype, key.split('|')[1]) +
      (live ? '' : ' <span style="color:#888">(variante nuage non sélectionnée)</span>'),
      font:{size: half ? 12 : 13} };
    lay.annotations = [{ xref:'paper', yref:'paper', x:0.02, y:0.98, xanchor:'left',
      yanchor:'top', showarrow:false, align:'left',
      bgcolor:'rgba(255,255,255,.85)', bordercolor:'#ccc', borderwidth:1, borderpad:4,
      font:{size:11.5},
      text:'pente <b>' + (p.slope > 0 ? '+' : '') + p.slope.toFixed(1) + ' ± ' +
           p.se.toFixed(1) + ' %/km</b><br>' + p.n.toLocaleString('fr-FR') + ' profils, ' +
           p.ndays + ' jours' }];
    lay.xaxis = { title:{text:'constante par profil (1/C) / médiane de l\'unité [%]'},
                  range:[H.x[0], H.x[H.x.length - 1]] };
    lay.yaxis = { title:{text:'Hauteur de base de nuage [km]'},
                  range:[H.y[0], H.y[H.y.length - 1]] };
    lay.shapes = [{ type:'line', xref:'x', yref:'paper', x0:100, x1:100, y0:0, y1:1,
                    line:{color:'#333', width:1, dash:'dash'} }];
    Plotly.react(el, traces, lay, PCFG);
  });
}

/* ---------- PWV panel (static) — binned medians + Theil-Sen, coherent states only ----------
   The raw per-hour scatter was unreadable (678 points x 5 states on one axis); what the test
   needs is the TREND, so each coherent state is shown as its PWV-binned median (2 mm bins,
   >= 8 h per bin) with the inter-quartile band, and the slope quoted in the legend is the
   Theil-Sen estimator over the individual hours — median of pairwise slopes, insensitive to the
   heavy residual tails.  The crossed states mix two signal definitions and are not drawn. */
const PWV_BIN = 2.0, PWV_MIN_H = 8;
function theilSen(x, y) {
  const slopes = [];
  for (let i = 0; i < x.length; i++)
    for (let j = i + 1; j < x.length; j++) {
      const dx = x[j] - x[i];
      if (dx !== 0) slopes.push((y[j] - y[i]) / dx);
    }
  if (!slopes.length) return null;
  slopes.sort((a, b) => a - b);
  const s = quantiles(slopes)[0];
  const r = x.map((v, i) => y[i] - s * v).sort((a, b) => a - b);
  return { slope:s, intercept:quantiles(r)[0], n:x.length };
}
function pwvBins(x, y) {
  const bins = {};
  for (let i = 0; i < x.length; i++) {
    const b = Math.floor(x[i] / PWV_BIN);
    (bins[b] = bins[b] || []).push(y[i]);
  }
  const out = [];
  Object.keys(bins).map(Number).sort((a, b) => a - b).forEach(b => {
    const v = bins[b];
    if (v.length < PWV_MIN_H) return;
    v.sort((a, c) => a - c);
    const [m, q1, q3] = quantiles(v);
    out.push({ xc:(b + 0.5) * PWV_BIN, m, q1, q3, n:v.length });
  });
  return out;
}
function drawPwv() {
  const P = S(), W = P.pwv, el = document.getElementById('pwv_scatter');
  if (!el) return;
  if (!W) { Plotly.purge(el); return; }
  let xmin = 1e9, xmax = -1e9, ylo = 1e9, yhi = -1e9;
  const traces = [];
  W.states.filter(st => !st[2].includes('croisé')).forEach(([variant, wv, label, col, sym]) => {
    const key = variant + '|' + wv, r = W.resid[key];
    if (!r) return;
    const x = [], y = [];
    r.forEach((v, i) => { if (v !== null && W.pwv_mm[i] !== null) {
      x.push(W.pwv_mm[i]); y.push(v); } });
    if (x.length < PWV_MIN_H) return;
    const bins = pwvBins(x, y);
    if (!bins.length) return;
    const ts = theilSen(x, y);
    bins.forEach(b => { ylo = Math.min(ylo, b.q1); yhi = Math.max(yhi, b.q3); });
    xmin = Math.min(xmin, bins[0].xc - PWV_BIN); xmax = Math.max(xmax, bins[bins.length-1].xc + PWV_BIN);
    /* IQR band */
    traces.push({ x:bins.map(b => b.xc).concat(bins.map(b => b.xc).reverse()),
      y:bins.map(b => b.q3).concat(bins.map(b => b.q1).reverse()),
      fill:'toself', fillcolor:rgba(col, 0.13), line:{width:0}, hoverinfo:'skip',
      showlegend:false });
    /* binned median line */
    traces.push({ x:bins.map(b => b.xc), y:bins.map(b => b.m), mode:'lines+markers',
      line:{color:col, width:2.4}, marker:{size:7, color:col, symbol:sym},
      name:label + (ts ? ' — Theil-Sen ' + (ts.slope > 0 ? '+' : '') + ts.slope.toFixed(2) +
                         ' %/mm (n=' + ts.n + ' h)' : ''),
      text:bins.map(b => b.n),
      hovertemplate:'PWV %{x:.0f} mm<br>médiane %{y:+.1f}% (%{text} h)<extra>' + variant + '</extra>' });
    /* Theil-Sen fit line, thin */
    if (ts) traces.push({ x:[bins[0].xc - PWV_BIN/2, bins[bins.length-1].xc + PWV_BIN/2],
      y:[ts.intercept + ts.slope * (bins[0].xc - PWV_BIN/2),
         ts.intercept + ts.slope * (bins[bins.length-1].xc + PWV_BIN/2)],
      mode:'lines', line:{color:col, width:1.2, dash:'dot'}, showlegend:false,
      hoverinfo:'skip' });
  });
  const lay = JSON.parse(JSON.stringify(BASE));
  lay.height = 500;
  lay.title = { text:'Résidu CL61 (Rayleigh) vs CHM15k — médianes par classes de PWV de ' +
                      PWV_BIN.toFixed(0) + ' mm (CAMS)', font:{size:13} };
  lay.xaxis = { title:{ text:'PWV — Precipitable Water Vapour (eau précipitable intégrée) [mm]' },
                zeroline:false, range:[Math.max(0, xmin), xmax] };
  lay.yaxis = { title:{ text:'résidu relatif [%] (médiane ' + W.band[0].toFixed(0) + '–' +
                W.band[1].toFixed(0) + ' m)' }, zeroline:true, zerolinewidth:1.4,
                zerolinecolor:'#888' };
  if (yhi > ylo) lay.yaxis.range = [ylo - 4, yhi + 4];
  Plotly.react(el, traces, lay, PCFG);
}

/* ===================== control surface ===================== */
function nightsBadge(id, method, variant) {
  const rec = S().calib[id + '|' + method + '|' + variant];
  if (!rec || !rec.ok) return '<span class="nights thin"><b>—</b> nuits</span>';
  const thin = rec.nights < 20 ? ' thin' : '';
  return `<span class="nights${thin}"><b>${rec.nights}</b> nuits</span>`;
}

function buildControls() {
  const P = S();
  /* per-instrument rows */
  const rows = P.instruments.map((inst, k) => {
    const id = inst.ident, sel = state.pick[id];
    const methods = METHOD_BY_TYPE[inst.itype] || ['rayleigh'];
    const mbtn = ['rayleigh', 'cloud'].map(m => {
      const has = methods.indexOf(m) >= 0;
      const why = has ? '' : (METHOD_WHY[inst.itype + '|' + m] || 'non disponible');
      return `<button data-id="${id}" data-method="${m}" ${has ? '' : 'disabled'}
        title="${why}" class="${sel.method === m ? 'on' : ''}">${m === 'rayleigh' ? 'Rayleigh' : 'Nuages'}</button>`;
    }).join('');
    /* only entries that are physically meaningful for this instrument; a meaningful one whose run
       does not exist yet is kept but disabled, with the reason as its tooltip */
    const opts = (((VARIANTS_FOR[inst.itype] || {})[sel.method]) || []).map(v => {
      const rec = P.calib[id + '|' + sel.method + '|' + v];
      const ok = rec && rec.ok;
      const why = ok ? `${rec.nights} nuits` : (rec ? rec.why : 'variante absente du payload');
      return `<option value="${v}" ${ok ? '' : 'disabled'} ${v === sel.variant ? 'selected' : ''}
        title="${(why || '').replace(/"/g, '&quot;')}">${vLabel(inst.itype, v)}${ok ? '' : ' — indisponible'}</option>`;
    }).join('');
    const badge = REFERENCE_VARIANTS.indexOf(sel.variant) >= 0
      ? `<span class="pill ref" title="${REFERENCE_NOTE.replace(/"/g, '&quot;')}">★ hypothèse de référence</span>` : '';
    /* profile-side measured-dark checkbox: only where the covered-telescope campaign measured b(z) */
    const hasDark = !!P.dark[id];
    /* two half-states: profile dark-corrected but constants not (no twin run), and the reverse
       (a dark-run variant with the profile-side box unticked) */
    const halfProfile = sel.dark && !DARK_TWIN[sel.variant];
    const halfConst = !sel.dark && hasDark && DARK_RUNS.indexOf(sel.variant) >= 0;
    /* profile-only dark on a CLOUD variant is NOT a warning: the cloud constant is dark-immune
       (CL31 +0.08 % measured, CL61 ~0.001 %), so only the profiles need correcting */
    const warn = halfProfile
      ? (sel.method === 'cloud'
        ? `<span class="warnflag ok" title="${DARK_CLOUD_OK.replace(/"/g, '&quot;')}">ⓘ</span>`
        : `<span class="warnflag" title="${DARK_NO_TWIN_WHY.replace(/"/g, '&quot;')}">⚠</span>`)
      : (halfConst
      ? `<span class="warnflag" title="constante issue d'un run dark-corrigé mais profil non corrigé — recocher la case pour l'état cohérent">⚠</span>` : '');
    const darkCell = `<label class="chk" title="${hasDark ? 'soustrait le profil b(z) mesuré ' +
        'télescope couvert avant la division par C_L' : DARK_NO_MEAS_WHY.replace(/"/g, '&quot;')}">
      <input type="checkbox" class="dchk" data-id="${id}" ${sel.dark ? 'checked' : ''}
       ${hasDark ? '' : 'disabled'}> dark mesuré (capot)</label>${warn}`;
    return `<tr>
      <td><span class="swatch" style="background:${inst.color}"></span><b>${inst.label}</b>
          <span class="muted">${inst.itype}</span></td>
      <td><span class="seg mseg">${mbtn}</span></td>
      <td><select class="vsel" data-id="${id}">${opts}</select>${badge}</td>
      <td>${nightsBadge(id, sel.method, sel.variant)}</td>
      <td>${darkCell}</td>
      <td><label class="chk"><input type="radio" name="iref" value="${k}"
           ${state.iref === k ? 'checked' : ''}> référence</label></td></tr>`;
  }).join('');
  const anyDark = P.instruments.some(i => state.pick[i.ident].dark);
  document.getElementById('instr').innerHTML =
    `<thead><tr><th>Instrument</th><th>Méthode</th><th>Variante d'étalonnage</th>
      <th>Disponibilité</th><th>Fond électronique</th><th>Différences /</th></tr></thead>
     <tbody>${rows}</tbody>`;
  document.getElementById('darkhint').innerHTML = anyDark
    ? `<span class="muted">${DARK_HINT}</span>` : '';

  /* corrections row — inert for an all-1064 nm site */
  const has910 = P.instruments.some(i => i.itype !== 'CHM15k');
  document.getElementById('corrrow').classList.toggle('off', !has910);
  document.querySelectorAll('#corrrow button, #corrrow select').forEach(e => e.disabled = !has910);
  document.getElementById('corrwhy').textContent = has910 ? '' :
    '(tous les instruments à 1064 nm : ces corrections sont des no-ops exacts)';
  document.querySelectorAll('#wvseg button').forEach(b =>
    b.classList.toggle('on', b.dataset.wv === state.wv));
  document.querySelectorAll('#wlseg button').forEach(b =>
    b.classList.toggle('on', b.dataset.wl === state.wl));
  /* the "constructeur" WV mode only means something where a CL61 is present — and the fallback
     for the OTHER 910 nm types must be disclosed even when a CL61 sits next to them (the button
     label names the CL61 spectrum; a CL31 under it is still corrected with its own 909,7/6,0) */
  const hasCL61 = P.instruments.some(i => i.itype === 'CL61');
  const other910 = P.instruments.some(i => i.itype !== 'CHM15k' && i.itype !== 'CL61');
  const fb = 'les autres types à 910 nm (CL31/CL51) gardent leur spectre nominal propre — seul ' +
             'le spectre du CL61 change dans ce mode';
  const bctor = document.querySelector('#wvseg button[data-wv="ctor"]');
  bctor.disabled = !has910;
  bctor.title = hasCL61 ? (REFERENCE_NOTE + (other910 ? ' ' + fb.charAt(0).toUpperCase() +
                                             fb.slice(1) + '.' : ''))
    : 'aucun CL61 ici : ce mode retombe sur le spectre nominal de chaque instrument';

  /* months */
  ['r0', 'r1'].forEach(w => {
    const s = document.getElementById(w);
    s.innerHTML = P.months.map((m, i) => `<option value="${i}">${MONTH_LABEL(m)}</option>`).join('');
    s.value = state[w];
  });
  fillL2Pill();
  bindRow();
}

function bindRow() {
  document.querySelectorAll('#instr .mseg button').forEach(b => b.addEventListener('click', () => {
    if (b.disabled) return;
    const id = b.dataset.id, m = b.dataset.method, sel = state.pick[id];
    /* remember the profile-dark of the method we leave, restore it when we come back — a
       cloud->rayleigh->cloud round trip must not silently drop the cloud channel's dark */
    sel.dmem = sel.dmem || {};
    sel.dmem[sel.method] = sel.dark;
    sel.method = m;
    sel.variant = firstAvailable(id, m);
    sel.prev = null;
    if (m === 'rayleigh')  /* authoritative both ways: set for dark runs, CLEAR otherwise */
      sel.dark = DARK_RUNS.indexOf(sel.variant) >= 0 && !!S().dark[id];
    else
      sel.dark = (sel.dmem[m] !== undefined ? sel.dmem[m] : sel.dark) && !!S().dark[id];
    buildControls(); draw();
  }));
  document.querySelectorAll('#instr .vsel').forEach(s => s.addEventListener('change', e => {
    const sel = state.pick[e.target.dataset.id];
    /* defence in depth: browsers keep disabled options unselectable, but a programmatic set or
       form-state restoration could commit a not-ok variant and silently drop the instrument */
    const opt = e.target.selectedOptions && e.target.selectedOptions[0];
    if (opt && opt.disabled) { e.target.value = sel.variant; return; }
    sel.variant = e.target.value;
    sel.prev = null;
    /* Rayleigh: the dropdown is authoritative for the dark pairing — a dark-run variant ticks
       the profile-side dark (where the measured b(z) exists), anything else unticks it, so the
       constants/profile pair stays coherent without a second click.  Cloud: the constants are
       dark-immune, the profile-side dark is an independent (and recommended) choice — changing
       the variant must NOT silently untick it. */
    if (sel.method === 'rayleigh')
      sel.dark = DARK_RUNS.indexOf(sel.variant) >= 0 && !!S().dark[e.target.dataset.id];
    buildControls(); draw();
  }));
  /* Ticking the measured dark switches the CONSTANTS to their dark-corrected twin where one
     exists: a dark-corrected profile divided by a dark-free constant mixes two definitions of the
     signal and measurably worsens the comparison.  Unticking restores the previous choice. */
  document.querySelectorAll('#instr .dchk').forEach(c => c.addEventListener('change', e => {
    const id = e.target.dataset.id, sel = state.pick[id], P = S();
    sel.dark = e.target.checked;
    if (sel.dark) {
      const twin = DARK_TWIN[sel.variant];
      if (twin && (P.calib[id + '|' + sel.method + '|' + twin] || {}).ok && twin !== sel.variant) {
        sel.prev = sel.variant;
        sel.variant = twin;
      }
    } else if (sel.prev) {
      sel.variant = sel.prev;
      sel.prev = null;
    }
    buildControls(); draw();
  }));
  document.querySelectorAll('#instr input[name=iref]').forEach(r =>
    r.addEventListener('change', e => { state.iref = +e.target.value; draw(); }));
}

/* Preferred variant when the operator flips a method: the manufacturer-spectrum reference if that
   run exists for this channel, else the current operational model, else whatever is available. */
function firstAvailable(id, method) {
  const P = S(), inst = P.instruments.find(i => i.ident === id);
  const list = ((VARIANTS_FOR[inst.itype] || {})[method]) || [];
  const ok = v => (P.calib[id + '|' + method + '|' + v] || {}).ok;
  const pref = REFERENCE_VARIANTS.concat(['v2.2', 'cloudWV']).filter(v => list.indexOf(v) >= 0);
  return pref.find(ok) || list.find(ok) || list[0];
}

function setSite(key) {
  state.site = key;
  const P = S();
  state.pick = {};
  const dd = P.default_dark || {};
  P.instruments.forEach(i => {
    const d = P.default[i.ident] || ['rayleigh', 'v2.2'];
    const method = (METHOD_BY_TYPE[i.itype] || ['rayleigh']).indexOf(d[0]) >= 0 ? d[0]
                 : (METHOD_BY_TYPE[i.itype] || ['rayleigh'])[0];
    const rec = P.calib[i.ident + '|' + method + '|' + d[1]];
    const variant = (rec && rec.ok) ? d[1] : firstAvailable(i.ident, method);
    /* keep the boot dark coherent even when firstAvailable substituted the default variant:
       on a Rayleigh channel the pairing rule (dark iff dark-run) overrides default_dark */
    let dk = !!(dd[i.ident] && P.dark[i.ident]);
    if (method === 'rayleigh') dk = DARK_RUNS.indexOf(variant) >= 0 && !!P.dark[i.ident];
    state.pick[i.ident] = { method, prev:null, dark:dk, variant };
  });
  state.iref = P.iref;
  state.r0 = 0; state.r1 = P.months.length - 1;
  document.querySelectorAll('#siteseg button').forEach(b => b.classList.toggle('on', b.dataset.site === key));
  document.getElementById('h_title').innerHTML = P.title;
  document.getElementById('h_sub').innerHTML = P.subtitle;
  document.getElementById('h_asof').textContent = P.name + ' · ' + P.wmo + ' · ' + P.start + ' → ' + P.end;
  document.getElementById('warns').innerHTML =
    (P.warnings || []).map(w => `<div class="warnbox">${w}</div>`).join('');
  const hasL2 = P.sources.indexOf('L2') >= 0;
  document.getElementById('l2card').style.display = hasL2 ? '' : 'none';
  document.getElementById('l1l2sec').style.display = hasL2 ? '' : 'none';
  document.getElementById('pwvsec').style.display = P.pwv ? '' : 'none';
  document.getElementById('calholder').innerHTML = P.instruments.map(i =>
    `<div class="card" style="margin-bottom:12px"><div id="cal_${i.ident}"></div></div>`).join('');
  document.getElementById('pcholder').innerHTML = P.instruments.map(i =>
    `<div class="card" style="margin-bottom:12px"><div id="pc_${i.ident}"></div></div>`).join('');
  document.getElementById('hkholder').innerHTML = '';
  document.getElementById('hkholder').dataset.n = '';
  buildControls();
  draw();
  drawPcolor();
  drawPwv();
}

document.querySelectorAll('#siteseg button').forEach(b =>
  b.addEventListener('click', () => setSite(b.dataset.site)));
document.querySelectorAll('#wvseg button').forEach(b => b.addEventListener('click', () => {
  if (b.disabled) return; state.wv = b.dataset.wv; buildControls(); draw(); }));
document.querySelectorAll('#wlseg button').forEach(b => b.addEventListener('click', () => {
  if (b.disabled) return; state.wl = b.dataset.wl; buildControls(); draw(); }));
document.getElementById('nfilt').addEventListener('change', e => { state.filter = e.target.checked; draw(); });
document.getElementById('logx').addEventListener('change', e => { state.logx = e.target.checked; draw(); });
['r0', 'r1'].forEach(w => document.getElementById(w).addEventListener('change', e => {
  state[w] = +e.target.value;
  if (state.r0 > state.r1) state[w === 'r0' ? 'r1' : 'r0'] = +e.target.value;
  document.getElementById('r0').value = state.r0;
  document.getElementById('r1').value = state.r1;
  draw();
}));

setSite(SITE_ORDER[0]);
window.addEventListener('resize', () => document.querySelectorAll('.js-plotly-plot')
  .forEach(el => Plotly.Plots.resize(el)));
"""


def html(payloads):
    site_buttons = "\n".join(
        '      <button data-site="%s"%s>%s</button>'
        % (k, ' class="on"' if i == 0 else "", payloads[k]["name"])
        for i, k in enumerate(SITE_ORDER) if k in payloads)

    js_consts = (
        f"const SITES = {json.dumps({k: payloads[k] for k in SITE_ORDER if k in payloads})};\n"
        f"const SITE_ORDER = {json.dumps([k for k in SITE_ORDER if k in payloads])};\n"
        f"const VARIANTS = {json.dumps(V3.VARIANTS)};\n"
        f"const VARIANT_LABEL = {json.dumps(V3.VARIANT_LABEL)};\n"
        f"const VSHORT = {json.dumps(V3.VARIANT_SHORT)};\n"
        f"const VLBL_BY_TYPE = {json.dumps(V3.VARIANT_LABEL_BY_TYPE)};\n"
        f"const VSH_BY_TYPE = {json.dumps(V3.VARIANT_SHORT_BY_TYPE)};\n"
        f"const METHOD_BY_TYPE = {json.dumps(V3.METHOD_BY_TYPE)};\n"
        f"const METHOD_WHY = {json.dumps({f'{a}|{b}': w for (a, b), w in V3.METHOD_WHY.items()})};\n"
        f"const VARIANTS_FOR = {json.dumps({t: {m: V3.variants_for(t, m) for m in ms} for t, ms in V3.METHOD_BY_TYPE.items()})};\n"
        f"const DARK_TWIN = {json.dumps(V3.DARK_TWIN)};\n"
        f"const DARK_RUNS = {json.dumps(list(V3.DARK_RUNS))};\n"
        f"const DARK_NO_TWIN_WHY = {json.dumps(V3.DARK_NO_TWIN_WHY)};\n"
        f"const DARK_CLOUD_OK = {json.dumps(V3.DARK_CLOUD_OK_WHY)};\n"
        f"const DARK_NO_MEAS_WHY = {json.dumps(V3.DARK_NO_MEAS_WHY)};\n"
        f"const DARK_HINT = {json.dumps(V3.DARK_HINT)};\n"
        f"const REFERENCE_VARIANTS = {json.dumps(list(V3.REFERENCE_VARIANTS))};\n"
        f"const REFERENCE_NOTE = {json.dumps(V3.REFERENCE_NOTE)};\n")

    return f"""<!doctype html>
<html lang="fr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>E-PROFILE ALC — inter-comparaison L1/L2 (v3)</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<div class="topbar">
  <span class="brand">E-PROFILE ALC · inter-comparaison</span>
  <span class="siteseg" id="siteseg">
{site_buttons}
  </span>
  <span class="spacer"></span>
  <span class="asof" id="h_asof"></span>
</div>
<div class="container">
  <h1 id="h_title"></h1>
  <p class="subtitle" id="h_sub"></p>
  <div id="warns"></div>

  <div class="card ctlwrap">
    <div class="tscroll"><table class="instr" id="instr"></table></div>
    <div id="darkhint" style="padding:2px 8px 0"></div>
    <div class="ctl-row" id="corrrow">
      <div class="ctl-group">
        <span class="ctl-label">Correction longueur d'onde</span>
        <span class="seg" id="wlseg">
          <button data-wl="none">aucune</button>
          <button data-wl="angstrom" title="toute la rétrodiffusion multipliée par
(λ/λcible)^−α avec α=1 — moléculaire et aérosol confondus">simple (Ångström α=1)</button>
          <button data-wl="molecular" class="on" title="le moléculaire analytique (Rayleigh,
T/p CAMS) est retiré à λ source, le résidu AÉROSOL est mis à l'échelle Ångström (λ/λcible)^α
avec α=1, puis le moléculaire calculé à λ cible est ré-ajouté — moléculaire exact, α=1 ne
s'applique qu'à l'aérosol">avancée (moléculaire + aérosol α=1)</button>
        </span>
      </div>
      <div class="ctl-group">
        <span class="ctl-label">Correction vapeur d'eau (comparaison)</span>
        <span class="seg" id="wvseg">
          <button data-wv="none">aucune</button>
          <button data-wv="nom" class="on">λ910,74 / FWHM 1,0</button>
          <button data-wv="ctor">λ910,55 σ0,08 ★</button>
        </span>
        <span class="muted" id="corrwhy"></span>
      </div>
    </div>
    <div class="ctl-row">
      <div class="ctl-group">
        <span class="ctl-label">Période</span>
        <select id="r0"></select><span class="muted">→</span><select id="r1"></select>
      </div>
      <div class="ctl-group">
        <label class="chk"><input type="checkbox" id="nfilt"> filtrage du plancher de bruit</label>
        <label class="chk"><input type="checkbox" id="logx"> axe β logarithmique</label>
      </div>
      <span class="muted"><b id="period">–</b> · <b id="nhours">–</b> heures appariées</span>
    </div>
  </div>

  <h2>Tableau récapitulatif — toutes les variantes <span class="muted">— l'échelle complète des
    hypothèses d'étalonnage, d'un coup d'œil</span></h2>
  <div class="card">
    <div class="tscroll"><table class="stats ladder" id="ladder">
      <thead><tr><th>Instrument</th><th>Variante</th><th>Médiane Δ vs référence</th>
        <th>Dispersion robuste (1,4826·MAD)</th><th>Heures appariées</th>
        <th>Nuits d'étalonnage</th><th>dC/dCBH (par profil)</th></tr></thead>
      <tbody></tbody>
    </table></div>
    <p class="note">Toute l'échelle d'étalonnage à la fois : l'effet de chaque hypothèse (spectre
      de vapeur d'eau, fond électronique, version du pipeline) se lit sans passer d'un état à
      l'autre. Ligne surlignée = variante actuellement sélectionnée ; ligne « <b>L2
      distribué</b> » = le produit tel qu'expédié aujourd'hui (comparé au L2 de la référence),
      avec la provenance de sa constante par instrument sur la pastille — CHM15k et CL61 =
      calibration Rayleigh v1.0 opérationnelle (activée mi-2026 pour le CL61) ; CL31 = non
      calibré (constante par défaut, qui laisse passer l'échelle native) — et, en infobulle, les
      constantes réellement observées sur la fenêtre. Les lignes de l'instrument de
      référence donnent son écart L1–L2 lorsqu'un L2 existe, sinon « réf. ». Médiane et dispersion
      (médiane du rapport par échantillon, et 1,4826·MAD des mêmes rapports) portent sur la bande
      de statistiques et la période courantes ; la colonne dC/dCBH vient de la configuration
      par-profil correspondante (variantes nuage seulement).</p>
  </div>

  <h2>Profils verticaux <span class="muted">— médiane (trait) et intervalle interquartile
    (bande). Courbes brutes par défaut ; cocher « filtrage du plancher de bruit » pour arrêter
    chaque courbe là où l'instrument n'a plus de signal. Vue par défaut 0–4 km — dézoomer pour le
    reste.</span></h2>
  <div class="grid2">
    <div class="card"><p class="panel-title">L1 + étalonnage v2 <span class="pill l1">C<sub>L</sub>
      Kalman</span></p><p class="panel-sub">rcs_0 / C<sub>L</sub>(t) × 10<sup>6</sup></p>
      <div id="p_l1"></div></div>
    <div class="card" id="l2card"><p class="panel-title">L2 tel que distribué
      <span class="pill l2" id="l2prov_pill">constante appliquée</span></p>
      <p class="panel-sub">attenuated_backscatter_0, sans ré-étalonnage</p>
      <div id="p_l2"></div></div>
  </div>

  <h2>Différences vs la référence <span class="muted">— à gauche le profil médian des rapports
    horaires, à droite la distribution de chaque échantillon apparié dans la bande de
    statistiques.</span></h2>
  <div class="grid2">
    <div class="card"><div id="d_diff"></div></div>
    <div class="card"><div id="p_hist"></div></div>
  </div>

  <div id="l1l2sec">
    <h2>L1 vs L2 <span class="muted">— ce que le ré-étalonnage change, instrument par
      instrument</span></h2>
    <div class="card"><div id="d_l1l2"></div></div>
    <p class="note">Positif = le profil L1 ré-étalonné est plus haut que le produit L2 à cette
      altitude. Le CL31 est le cas extrême : sa constante L2 est le défaut 1e8 non étalonné.</p>
  </div>

  <div id="hksec">
    <h2>Constante par profil vs hauteur de base de nuage <span class="muted">— l'étalonnage nuage
      profil par profil, pour la variante nuage sélectionnée ci-dessus</span></h2>
    <div id="hkholder"></div>
    <p class="note">Chaque profil calibré individuellement (et non la médiane du jour) : x = sa
      constante 1/C en % de la médiane de l'unité, y = la hauteur de base du nuage utilisé. Une
      <b>pente non nulle</b> veut dire que la constante retrouvée dépend de l'épaisseur
      d'atmosphère traversée — c'est précisément ce que produit une correction de vapeur d'eau trop
      forte. La pente est un moindre-carré <b>groupé par jour</b> (les profils d'un même jour
      partagent la météo et un seul état d'étalonnage, une erreur-type par profil serait
      absurdement petite). Ces jeux viennent de
      <code>rayleigh_availability/cloud_profile_dump.py</code>, un dossier par configuration WV.</p>
  </div>

  <div id="pwvsec">
    <h2>Résidu CL61 (Rayleigh) vs PWV <span class="muted">— le test discriminant du spectre
      laser</span></h2>
    <div class="card"><div id="pwv_scatter"></div></div>
    <p class="note"><b title="Precipitable Water Vapour">PWV</b> = <i>Precipitable Water
      Vapour</i>, l'eau précipitable intégrée sur toute la colonne au-dessus de la station, en mm
      (1 mm = 1 kg m⁻²), calculée depuis les mêmes fichiers CAMS et la même fenêtre journalière que
      la correction elle-même. Chaque trait = la <b>médiane du résidu par classe de PWV de 2 mm</b>
      (classes d'au moins 8 heures ; bande = intervalle interquartile) ; la pente en légende est
      l'estimateur robuste de <b>Theil-Sen</b> sur les heures individuelles, insensible aux queues
      de distribution. Seuls les <b>trois états cohérents</b> (même hypothèse de raie des deux
      côtés) sont tracés. Une <b>pente plate</b> signale l'hypothèse de spectre correcte : si la
      correction de vapeur d'eau est du bon ordre, le résidu ne dépend plus de la quantité de
      vapeur d'eau présente. Panneau statique — il compare des états fixes et ne suit pas les
      commandes ci-dessus ; sa référence est l'état par défaut du site <b>pris dans son
      appariement cohérent</b> (à Payerne : constantes v2.2+dark ET profil corrigé du dark mesuré
      — la revue du 2026-08-16 a montré qu'une référence hybride gonflait chaque pente de
      +0,4 à +0,8 %/mm).</p>
  </div>

  <h2>Constantes d'étalonnage dans le temps</h2>
  <div id="calholder"></div>
  <p class="note">Pour chaque instrument : toutes les variantes disponibles de la méthode
    sélectionnée (la variante active en trait plein avec ses nuits, les autres en pointillé), plus
    la constante réellement appliquée dans le L2 lorsqu'elle existe. La bande grisée est la période
    dont sont tirés les panneaux ci-dessus.</p>

  <h2>Rideaux temps–altitude <span class="muted">— L1 étalonné, échelle log₁₀ −2…1, état par
    défaut de la page ; une heure sur deux et une porte sur deux jusqu'à 8 km. Statiques.</span></h2>
  <div id="pcholder"></div>

  <p class="note" style="margin-top:26px"><b>Exactitude du recalcul côté navigateur.</b> Les
    profils, les différences, l'histogramme et le tableau sont recalculés <b>heure par heure</b>
    dans la page à partir du signal brut et des facteurs de correction par (jour, porte) — aucune
    approximation « profil médian ÷ C<sub>L</sub> médian » n'est utilisée : mesurée, elle atteint
    <b>46 % sur le CL31 de Payerne</b> (dont la constante varie d'un facteur 2,1 sur la fenêtre à
    cause du remplacement du bloc optique), 34 % à Amsterdam et 12 % à Lindenberg — et elle n'est
    même pas définie quand la conversion de longueur d'onde ajoute son terme moléculaire additif,
    qui vit en unités étalonnées. Seul l'<b>ordre</b> des opérations change par rapport à la v2 : la
    v3 ré-échantillonne puis corrige, la v2 corrigeait puis ré-échantillonnait. Contrôle de non
    régression sur Payerne dans l'<b>état hérité de la v2</b> (CHM15k v2.2 / CL31 nuage λ909,7 /
    CL61 nuage λ910,74, λ avancée, WV comparaison λ910,74, dark décoché, filtrage de bruit
    désactivé — la page démarre désormais dans l'état de référence de l'opérateur, pas dans
    celui-ci ; le contrôle le reconstruit programmatiquement, <code>check_v3.py</code>) :
    <b>N = 56 952 paires identiques</b>, log r reproduit à ±0,0003 près, biais relatif médian
    <b>−11,54 % vs −11,55 %</b> (CL61 nuage), <b>+8,69 % vs +8,68 %</b> (CL61 Rayleigh v2.2) et
    <b>−0,31 % vs −0,47 %</b> (CL31, l'écart maximal : 0,16 point de pourcentage, sur l'instrument
    dont le signal est le plus proche de zéro dans la bande). Le stockage quantifié (uint16
    logarithmique, ~0,007 % de pas relatif) contribue pour moins de 0,001 point.</p>

  <div class="foot">
    Construit par <code>inter-comparison_dashboard/build_v3.py</code> +
    <code>render_v3.py</code> (matrice <code>variants_v3.py</code>, contrôle
    <code>check_v3.py</code>) · les commandes recombinent des ingrédients par canal dans le
    navigateur : rien n'est pré-calculé par combinaison · médianes horaires, ≥30 min de couverture,
    écrantage nuage/brouillard ±15 min · CAMS résolu <b>par jour</b>.
  </div>
</div>
<script>{js_consts}</script>
<script>{JS}</script>
</body></html>
"""


def main():
    payloads = {}
    for k in SITE_ORDER:
        f = DATA / f"v3_{k}.json"
        if f.exists():
            payloads[k] = json.loads(f.read_text(encoding="utf-8"))
            print(f"   {k}: {f.stat().st_size/1e6:.2f} MB")
    if not payloads:
        raise SystemExit("no v3_<site>.json — run build_v3.py first")
    out = HERE / "index.html"
    out.write_text(html(payloads), encoding="utf-8")
    print(f"-> {out}  ({out.stat().st_size/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
