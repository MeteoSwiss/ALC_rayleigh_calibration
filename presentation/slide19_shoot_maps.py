# -*- coding: utf-8 -*-
"""Annual Meeting 2026 deck, slide 19: the three network maps of the LIVE calibration monitor
(summary page): success rate, mean OmB bias vs CAMS, ICAO detection altitude.

Each map is the dashboard's own Plotly figure (monitoring/charts.py: network_map, omb_bias_map,
icao_altitude_map; divs fig-map, fig-map-omb, fig-map-icao). After it renders, only the LAYOUT is
changed (size, fonts, title, horizontal colorbar, symbol legend hidden -- the slide carries one
shared legend) via Plotly.relayout/restyle; the data are untouched. The three shots are stitched
side by side into one PNG, and the per-type statistics of the plotted values are dumped to JSON.

Run from the repo root:  python presentation/slide19_shoot_maps.py
Output: presentation/figs_annual2026/fig_maps_slide.png (+ slide19_maps_stats.json)
"""
import json
from pathlib import Path

import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright

BASE = "https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host"
OUTDIR = Path(__file__).resolve().parent / "figs_annual2026"
OUT = OUTDIR / "fig_maps_slide.png"
STATS = OUTDIR / "slide19_maps_stats.json"
W, H = 576, 610          # CSS px of one panel = its size on the 1920x1080 slide
MAPS = [("fig-map", "Calibration success rate (%)", 1.0),
        ("fig-map-omb", "Mean O−B vs CAMS (Mm⁻¹ sr⁻¹)", 1.0),
        ("fig-map-icao", "ICAO detection altitude (km)", 0.001)]

RELAYOUT = """
([id, title, W, H, scale]) => {
  // Re-plot the dashboard's own data (decoded arrays from _fullData) in a fresh, fixed-size div:
  // same points, colours, colour scale and map extent; only the layout is slide-sized.
  const gd = document.getElementById(id);
  if (!gd || !window.Plotly || !gd._fullData) return {err: 'no-plot ' + id};
  const fd = gd._fullData[0], fl = gd._fullLayout, g = fl.geo;
  const A = (v) => (v && v.length !== undefined && typeof v !== 'string') ? Array.from(v) : v;
  const col = A(fd.marker.color).map(v => (v === null || v === undefined) ? null : v * scale);
  const out = {id: id, n: (A(fd.lat) || []).length, color: A(fd.marker.color), symbol: A(fd.marker.symbol),
               cd: (A(fd.customdata) || []).map(r => [r[0], r[1], r[3], r[4]])};
  const div = document.createElement('div');
  div.id = 'shot-' + id;
  div.style.cssText = `position:fixed;left:0;top:0;width:${W}px;height:${H}px;z-index:99999;background:#fff`;
  document.body.appendChild(div);
  const mk = fd.marker;
  const trace = {type: 'scattergeo', mode: 'markers', lat: A(fd.lat), lon: A(fd.lon), hoverinfo: 'skip',
    showlegend: false,
    marker: {size: A(mk.size), symbol: A(mk.symbol), color: col, colorscale: mk.colorscale,
             cmin: mk.cmin * scale, cmax: mk.cmax * scale, opacity: mk.opacity,
             line: {width: 0.5, color: '#555'}, showscale: true,
             colorbar: {orientation: 'h', x: 0.5, xanchor: 'center', y: -0.02, yanchor: 'top',
                        len: 0.9, thickness: 16, tickfont: {size: 18}, outlinewidth: 0}}};
  const layout = {width: W, height: H, autosize: false, paper_bgcolor: '#fff',
    font: {family: 'Arial, sans-serif', size: 18, color: '#1C1F23'},
    title: {text: title, x: 0.02, xanchor: 'left', y: 0.985, yanchor: 'top', font: {size: 23}},
    margin: {t: 46, l: 4, r: 4, b: 70},
    geo: {scope: g.scope, resolution: 50, showcountries: true, countrycolor: '#bbbbbb',
          showland: true, landcolor: '#f5f5f5', showframe: false, bgcolor: '#fff',
          showcoastlines: true, coastlinecolor: '#bbbbbb',
          projection: {type: g.projection.type},
          lataxis: {range: A(g.lataxis.range)}, lonaxis: {range: A(g.lonaxis.range)}}};
  return Plotly.newPlot(div, [trace], layout, {staticPlot: true, displayModeBar: false})
    .then(() => out);
}
"""

shots, stats = [], {}
with sync_playwright() as p:
    br = p.chromium.launch()
    pg = br.new_page(viewport={"width": 1900, "height": 1400}, device_scale_factor=2)
    try:
        pg.goto(BASE + "/", wait_until="networkidle", timeout=90000)
    except Exception as e:
        print("goto-warn:", e)
    pg.wait_for_timeout(8000)
    print("page text:", pg.evaluate("() => document.body.innerText.slice(0, 900)").replace("\n", " / "))
    for div, title, scale in MAPS:
        r = pg.evaluate(RELAYOUT, [div, title, W, H, scale])
        if "err" in r:
            print(r["err"]); continue
        pg.wait_for_timeout(2500)
        f = OUTDIR / f"_map_{div}.png"
        pg.locator("#shot-" + div).screenshot(path=str(f))
        shots.append(f)
        types = np.array([c[1] for c in r["cd"]])
        vals = np.array([np.nan if v is None else float(v) for v in r["color"]])
        st = {"n": int(r["n"])}
        for t in ["all"] + sorted(set(types)):
            m = np.isfinite(vals) & ((types == t) if t != "all" else True)
            if m.sum():
                v = vals[m]
                st[t] = {"n": int(m.sum()), "median": float(np.median(v)), "p25": float(np.percentile(v, 25)),
                         "p75": float(np.percentile(v, 75)), "mean": float(np.mean(v))}
        # per-country medians of the success map, for the conclusion
        if div == "fig-map":
            ctry = np.array([c[0] for c in r["cd"]])
            st["by_country_type"] = {}
            for c in sorted(set(ctry)):
                for t in sorted(set(types)):
                    m = (ctry == c) & (types == t) & np.isfinite(vals)
                    if m.sum() >= 3:
                        st["by_country_type"][f"{c}|{t}"] = [int(m.sum()), float(np.median(vals[m]))]
        stats[div] = st
    br.close()

ims = [Image.open(f).convert("RGBA") for f in shots]
h = max(i.height for i in ims)
canvas = Image.new("RGBA", (sum(i.width for i in ims), h), (255, 255, 255, 255))
x = 0
for i in ims:
    canvas.paste(i, (x, 0)); x += i.width
canvas.convert("RGB").save(OUT)
for f in shots:
    f.unlink()
STATS.write_text(json.dumps(stats, indent=1))
print("saved", OUT, canvas.size)
print(json.dumps({k: {t: v for t, v in s.items() if t != "by_country_type"} for k, s in stats.items()}, indent=1))
print("SHOOT_DONE")
