# -*- coding: utf-8 -*-
"""Slide 10b of the Annual Meeting 2026 deck: the paper Fig. 3 Rayleigh night (Hohenpeissenberg
CHM15k, 12-13 Sep 2026) re-shot from the LIVE dashboard daily panel with slide-sized fonts.

The panel is the dashboard's own Plotly figure (monitoring/panel.py); after it renders, only its
layout is changed (font sizes, height, legend position) via Plotly.relayout -- the data are untouched.

Run from the repo root:  python presentation/slide10_shoot_rayleigh_example.py
Output: presentation/figs_annual2026/fig_rayleigh_example_slide.png
"""
from pathlib import Path
from playwright.sync_api import sync_playwright

BASE = "https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host"
KEY, DAY = "0-20000-0-10962_0", "20260913"   # panel date = morning of the night (12-13 Sep)
OUT = Path(__file__).resolve().parent / "figs_annual2026" / "fig_rayleigh_example_slide.png"

RELAYOUT = """
() => {
  const gd = document.getElementById('panel');
  if (!gd || !window.Plotly || !gd.data) return 'no-plot';
  const names = gd.data.map(t => t.name || '').join(' | ');
  return Plotly.relayout(gd, {
    'font.size': 20, 'showlegend': false,
    'height': 410, 'margin.t': 14, 'margin.l': 80, 'margin.b': 92, 'margin.r': 16,
    'xaxis2.dtick': 1, 'xaxis2.exponentformat': 'none',
    'xaxis.title.font.size': 21, 'yaxis.title.font.size': 21, 'xaxis2.title.font.size': 21,
    'yaxis.title.text': 'Range (km)',
  }).then(() => {
    const cb = gd.data.findIndex(t => t.type === 'heatmap');
    if (cb >= 0) return Plotly.restyle(gd, {'colorbar.thickness': 18, 'colorbar.tickfont.size': 18,
                                            'colorbar.title.text': ''}, [cb]).then(() => names);
    return names;
  });
}
"""

with sync_playwright() as p:
    br = p.chromium.launch()
    pg = br.new_page(viewport={"width": 2050, "height": 1200}, device_scale_factor=2)
    url = f"{BASE}/stations/{KEY}.html#d={DAY}"
    try:
        pg.goto(url, wait_until="networkidle", timeout=60000)
    except Exception as e:
        print("goto-warn:", e)
    pg.wait_for_timeout(6000)
    # make sure the Rayleigh method is shown
    btn = pg.locator("#mswitch button", has_text="Rayleigh")
    if btn.count():
        btn.first.click()
        pg.wait_for_timeout(3000)
    print("date label:", pg.locator("#datelbl").inner_text(timeout=5000))
    for sel in ("#dp-flag", "#msg", "#why", "#curtitle"):
        try:
            print(sel, "::", pg.locator(sel).first.inner_text(timeout=3000).replace("\n", " / ")[:1500])
        except Exception as e:
            print(sel, ":: n/a", e)
    print("traces:", pg.evaluate(RELAYOUT))
    pg.wait_for_timeout(2500)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pg.locator("#panel").screenshot(path=str(OUT))
    print("saved", OUT)
    br.close()
print("SHOOT_DONE")
