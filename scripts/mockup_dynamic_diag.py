#!/usr/bin/env python3
"""Feasibility test: can the per-night diagnostic be a CLIENT-SIDE plot instead of a server PNG?

The per-night diagnostic image is the expensive product in this pipeline -- nine matplotlib axes,
five of them raster (three window heatmaps, the RCS curtain, the sensitivity grid) -- rendered once
per night per stream. Regenerating a full history of them was costed at 105-212 GB and hours of
CPU, which is why the release recompute runs with PLOTS=0.

This script answers whether the same information can be shipped as DATA and drawn in the browser:

  1. it times the calibration of one real night twice, PLOTS=0 and PLOTS=1, so the difference is
     the true marginal cost of the image (not a guess);
  2. it captures the exact arrays the plotting function was called with, and writes them as a
     compact payload (uint8 curtain, float32 profiles, small grids);
  3. it renders those arrays client-side with Plotly into a standalone mockup page, so the result
     can be compared with the PNG side by side.

Nothing is deployed: this writes a mockup and prints a comparison.

  python scripts/mockup_dynamic_diag.py --key 0-20000-0-06610_A --date 20260705 \
      --l1-root D:/E-PROFILE_L1_2026 --out C:/tmp/mock_dynamic_diag.html
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

CAPTURED: dict = {}


def _b64(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def _quant_u8(m: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """log-ish curtain -> uint8, 255 = missing. The same trick the v3 dashboard uses for its
    time-height curtains, where 1/255 of the colour range is invisible on screen."""
    q = np.full(m.shape, 255, dtype=np.uint8)
    ok = np.isfinite(m)
    if ok.any():
        v = np.clip((m[ok] - lo) / max(hi - lo, 1e-9), 0, 1)
        q[ok] = np.clip(np.round(v * 254.0), 0, 254).astype(np.uint8)
    return q


def _payload_from_kwargs(kw: dict) -> dict:
    """The captured plot arguments -> the compact payload a browser would fetch."""
    out: dict = {}
    rng = np.asarray(kw.get("range_alc", []), dtype=float)
    alt = float(kw.get("altitude", 0.0) or 0.0)
    z_km = (rng + alt) * 1e-3

    # --- 1-D profile panels: the physics an operator actually reads ---------------------------
    prof = {}
    for name in ("rcs_mean", "signal_normalized", "p_mol", "beta_att_mol"):
        v = kw.get(name)
        if v is None:
            continue
        prof[name] = _b64(np.asarray(v, dtype="float32"))
    out["z_km"] = _b64(z_km.astype("float32"))
    out["prof"] = prof
    out["n_z"] = int(z_km.size)

    # --- window grids (slope / intercept / R^2): small 2-D, ship as float32 --------------------
    grids = {}
    for name in ("slopes", "intercepts", "r_squared"):
        v = kw.get(name)
        if v is None:
            continue
        a = np.asarray(v, dtype="float32")
        grids[name] = {"b64": _b64(a), "shape": list(a.shape)}
    out["grids"] = grids
    for name in ("range_bin_m", "half_length_m"):
        v = kw.get(name)
        if v is not None:
            out[name] = _b64(np.asarray(v, dtype="float32"))

    # --- the RCS curtain: the only bulky panel -> stride + uint8 -------------------------------
    rcs = kw.get("rcs")
    if rcs is not None:
        m = np.asarray(rcs, dtype=float)
        st_t = max(1, int(np.ceil(m.shape[0] / 300)))     # 300x200 is plenty on screen; the PNG
        st_r = max(1, int(np.ceil(m.shape[1] / 200)))     # itself only draws 800x800
        sub = m[::st_t, ::st_r]
        with np.errstate(all="ignore"):
            lg = np.log10(np.where(sub > 0, sub, np.nan))
        fin = lg[np.isfinite(lg)]
        lo, hi = (float(np.percentile(fin, 5)), float(np.percentile(fin, 95))) if fin.size else (0.0, 6.0)
        out["curtain"] = {"b64": _b64(_quant_u8(lg.T, lo, hi)),
                          "shape": [int(lg.T.shape[0]), int(lg.T.shape[1])],
                          "lo": lo, "hi": hi,
                          "range_km": _b64((np.asarray(kw.get("rcs_range_alc", rng), float)[::st_r]
                                            * 1e-3).astype("float32")),
                          "n_t": int(sub.shape[0])}
    fw = (kw.get("fit_altitude_start"), kw.get("fit_altitude_end"))
    if all(v is not None for v in fw):
        out["fit_window_km"] = [float(fw[0]) * 1e-3, float(fw[1]) * 1e-3]
    out["title"] = str(kw.get("title", ""))
    return out


def _install_capture():
    """Wrap the production plotting function: time the real PNG render AND keep its arguments."""
    import calibration.plotting as P
    import calibration.rayleigh.calibration as RC

    real = P.plot_rayleigh_diagnostics_compact

    def wrapper(*a, **kw):
        t0 = time.perf_counter()
        fig = real(*a, **kw)
        dt = time.perf_counter() - t0
        CAPTURED.setdefault("renders", []).append(dt)
        CAPTURED["kwargs"] = kw
        sp = kw.get("save_path")
        if sp and Path(sp).exists():
            CAPTURED["png"] = Path(sp)
        return fig

    P.plot_rayleigh_diagnostics_compact = wrapper
    RC.plot_rayleigh_diagnostics_compact = wrapper


def _time_runner(key: str, date: str, plots: str) -> float:
    """Wall time of one night through the real runner, with images on or off."""
    env = dict(os.environ, PLOTS=plots)
    t0 = time.perf_counter()
    subprocess.run([sys.executable, str(REPO / "scripts/run_network_calibration.py"),
                    "--stream", key, "--start", date, "--end", date,
                    "--methods", "rayleigh", "--force", "--per-type", "0",
                    "--ignore-coverage", "--workers", "1"],
                   env=env, cwd=str(REPO), capture_output=True, text=True, timeout=1800)
    return time.perf_counter() - t0


HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Dynamic diagnostic — feasibility mockup</title>
<script>{plotly}</script>
<style>
 body {{ font: 14px/1.5 -apple-system,"Segoe UI",Roboto,Arial,sans-serif; margin: 22px 28px;
        color: #1a2530; }}
 .hdr {{ background:#fff3cd; border:1px solid #ffe08a; border-radius:8px; padding:10px 14px;
         margin-bottom:16px; }}
 .grid {{ display:grid; grid-template-columns:1fr 1fr 1fr; gap:14px; }}
 .card {{ border:1px solid #dbe3ea; border-radius:10px; padding:8px; }}
 table {{ border-collapse:collapse; margin:14px 0; }}
 td,th {{ border:1px solid #dbe3ea; padding:5px 10px; text-align:right; }}
 th:first-child, td:first-child {{ text-align:left; }}
 .win {{ font-weight:600; color:#b7791f; }}
</style></head><body>
<div class="hdr"><b>FEASIBILITY MOCKUP</b> — every plot below is drawn <b>in your browser</b> from a
{size_kb:.0f} KB payload ({gz_kb:.0f} KB gzipped) captured from the real calibration of
<b>{key} {date}</b>. No PNG is involved. The equivalent server-rendered image is
{png_kb:.0f} KB and took <b>{render_s:.2f} s</b> of matplotlib time to draw.</div>

<table>
<tr><th>quantity</th><th>server PNG</th><th>client-side payload</th></tr>
<tr><td>bytes on the wire</td><td>{png_kb:.0f} KB</td><td>{gz_kb:.0f} KB (gzip)</td></tr>
<tr><td>CPU to produce</td><td>{render_s:.2f} s (matplotlib)</td><td>{dump_ms:.0f} ms (array dump)</td></tr>
<tr><td>night with images vs without</td><td colspan="2">{ab}</td></tr>
</table>

<h3>RCS curtain — full night <span class="win">(gold = fit window)</span></h3>
<div class="card"><div id="curtain"></div></div>
<h3>Molecular panels</h3>
<div class="grid">
  <div class="card"><div id="p1"></div></div>
  <div class="card"><div id="p2"></div></div>
  <div class="card"><div id="p3"></div></div>
</div>
<h3>Window grids</h3>
<div class="grid">
  <div class="card"><div id="g_slopes"></div></div>
  <div class="card"><div id="g_intercepts"></div></div>
  <div class="card"><div id="g_r_squared"></div></div>
</div>

<script id="payload" type="application/json">{payload}</script>
<script>
const D = JSON.parse(document.getElementById('payload').textContent);
const dec = (b64, type) => {{
  const s = atob(b64), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  return type === 'u8' ? u : new Float32Array(u.buffer);
}};
const LAY = {{ margin:{{l:58,r:12,t:28,b:44}}, height:300, template:'plotly_white',
              font:{{size:11}} }};
const z = Array.from(dec(D.z_km, 'f32'));

function profile(div, series, title, xlab) {{
  const traces = series.filter(s => D.prof[s.k]).map(s => ({{
    x: Array.from(dec(D.prof[s.k], 'f32')), y: z, mode:'lines', name: s.n,
    line: {{ color: s.c, width: s.w || 1.1, dash: s.d }} }}));
  Plotly.newPlot(div, traces, Object.assign({{}}, LAY, {{
    title:{{text:title,font:{{size:12}}}},
    xaxis:{{title:{{text:xlab}}, type:'log', exponentformat:'e'}},
    yaxis:{{title:{{text:'Altitude a.s.l. [km]'}}}},
    showlegend:true, legend:{{orientation:'h', y:1.12, font:{{size:10}}}} }}), {{responsive:true}});
}}
profile('p1', [{{k:'rcs_mean', n:'RCS', c:'#1f77b4'}}], 'Molecular: RCS', 'RCS');
profile('p2', [{{k:'signal_normalized', n:'normalised', c:'#1f77b4'}},
               {{k:'p_mol', n:'molecular', c:'#d62728', d:'dash'}}],
        'Molecular: signal vs theory', 'signal');
profile('p3', [{{k:'beta_att_mol', n:'molecular β_att', c:'#d62728', d:'dash'}}],
        'Molecular: attenuated backscatter', 'β_att [m⁻¹ sr⁻¹]');

if (D.curtain) {{
  const c = D.curtain, q = dec(c.b64, 'u8'), [ny, nx] = c.shape;
  const rk = Array.from(dec(c.range_km, 'f32'));
  const zz = [];
  for (let j = 0; j < ny; j++) {{
    const row = new Array(nx);
    for (let i = 0; i < nx; i++) {{
      const v = q[j * nx + i];
      row[i] = v === 255 ? null : c.lo + (v / 254) * (c.hi - c.lo);
    }}
    zz.push(row);
  }}
  const shapes = D.fit_window_km ? [{{ type:'rect', xref:'paper', x0:0, x1:1,
      y0:D.fit_window_km[0], y1:D.fit_window_km[1], fillcolor:'rgba(255,193,7,0.18)',
      line:{{color:'#ffc107', width:1.3}} }}] : [];
  Plotly.newPlot('curtain', [{{ type:'heatmap', z:zz, y:rk, colorscale:'Viridis',
      zmin:c.lo, zmax:c.hi, colorbar:{{title:{{text:'log₁₀(RCS)'}}, thickness:13}},
      hovertemplate:'profile %{{x}}<br>%{{y:.2f}} km<br>log₁₀RCS %{{z:.2f}}<extra></extra>' }}],
    Object.assign({{}}, LAY, {{ height:380, shapes:shapes,
      title:{{text:'Range-corrected signal — full night', font:{{size:12}}}},
      xaxis:{{title:{{text:'profile index (strided)'}}}},
      yaxis:{{title:{{text:'Range [km]'}}}} }}), {{responsive:true}});
}}

['slopes','intercepts','r_squared'].forEach(function (g) {{
  const G = D.grids[g]; if (!G) return;
  const a = dec(G.b64, 'f32'), [n1, n2] = G.shape, zz = [];
  for (let j = 0; j < n2; j++) {{
    const row = new Array(n1);
    for (let i = 0; i < n1; i++) row[i] = a[i * n2 + j];
    zz.push(row);
  }}
  const scale = g === 'r_squared' ? 'YlGnBu' : (g === 'slopes' ? 'RdBu' : 'Magma');
  Plotly.newPlot('g_' + g, [{{ type:'heatmap', z:zz, colorscale:scale,
      colorbar:{{thickness:12}},
      hovertemplate:'%{{z:.3g}}<extra></extra>' }}],
    Object.assign({{}}, LAY, {{ title:{{text:'Window: ' + g, font:{{size:12}}}},
      xaxis:{{title:{{text:'range bin'}}}}, yaxis:{{title:{{text:'half-length bin'}}}} }}),
    {{responsive:true}});
}});
</script></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--key", required=True, help="stream key, e.g. 0-20000-0-06610_A")
    ap.add_argument("--date", required=True, help="YYYYMMDD of a night that CALIBRATES")
    ap.add_argument("--l1-root", required=True)
    ap.add_argument("--cams", default="")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("./_dyn_diag_work"))
    ap.add_argument("--skip-ab", action="store_true", help="skip the PLOTS=0/1 A/B timing")
    args = ap.parse_args()

    os.environ["ALC_L1_ROOT"] = args.l1_root
    os.environ["ALC_FULLCAL_DIR"] = str(args.work.resolve())
    if args.cams:
        os.environ["ALC_CAMS_DIR"] = args.cams
    args.work.mkdir(parents=True, exist_ok=True)

    ab = "not measured (--skip-ab)"
    if not args.skip_ab:
        print("timing the same night with images OFF then ON ...", flush=True)
        t_off = _time_runner(args.key, args.date, "0")
        t_on = _time_runner(args.key, args.date, "1")
        ab = (f"{t_off:.1f} s without images vs {t_on:.1f} s with -> "
              f"images cost {t_on - t_off:+.1f} s ({100 * (t_on - t_off) / max(t_off, 1e-9):+.0f} %)")
        print("  " + ab, flush=True)

    # in-process run with the plotting function wrapped, so the payload is EXACTLY the arrays the
    # PNG was drawn from -- not a re-derivation that could quietly differ
    os.environ["PLOTS"] = "1"
    _install_capture()
    import importlib
    from datetime import datetime
    rnc = importlib.import_module("scripts.run_network_calibration")
    census = json.loads((REPO / "validation/scope_l1_2026_census.json").read_text(encoding="utf-8"))
    rows = census if isinstance(census, list) else census.get("streams", [])
    s = next(r for r in rows if f"{r['wmo']}_{r['ident']}" == args.key)
    d = datetime.strptime(args.date, "%Y%m%d")
    print(f"calibrating {args.key} {args.date} in-process to capture the plot arrays ...", flush=True)
    rnc._do_rayleigh(s, d, d)

    if "kwargs" not in CAPTURED:
        print("\nNo diagnostic was drawn for that night -- pick a date that CALIBRATES "
              "(a rejected night draws the simpler failure plot).", file=sys.stderr)
        sys.exit(2)

    t0 = time.perf_counter()
    payload = _payload_from_kwargs(CAPTURED["kwargs"])
    raw = json.dumps(payload, separators=(",", ":"))
    dump_ms = (time.perf_counter() - t0) * 1e3
    gz_kb = len(gzip.compress(raw.encode())) / 1024
    png = CAPTURED.get("png")
    png_kb = (png.stat().st_size / 1024) if png and png.exists() else float("nan")
    render_s = max(CAPTURED.get("renders", [float("nan")]))

    import plotly.offline as pyo
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(HTML.format(plotly=pyo.get_plotlyjs(), payload=raw, key=args.key,
                                    date=args.date, size_kb=len(raw) / 1024, gz_kb=gz_kb,
                                    png_kb=png_kb, render_s=render_s, dump_ms=dump_ms, ab=ab),
                        encoding="utf-8")
    print(f"\n{'':-<64}")
    print(f"PNG           : {png_kb:8.0f} KB   rendered in {render_s:6.2f} s")
    print(f"payload       : {len(raw)/1024:8.0f} KB   built in    {dump_ms:6.0f} ms")
    print(f"payload gzip  : {gz_kb:8.0f} KB   ({png_kb / max(gz_kb, 1e-9):.1f}x smaller than the PNG)")
    print(f"speedup       : {render_s / max(dump_ms / 1e3, 1e-9):8.0f}x less CPU per night")
    print(f"mockup        : {args.out}  ({args.out.stat().st_size/1e6:.1f} MB incl. Plotly)")


if __name__ == "__main__":
    main()
