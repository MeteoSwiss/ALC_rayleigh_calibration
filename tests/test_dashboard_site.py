"""Regression tests for the station dashboard — one test per way it has actually broken.

Two tiers:

  UNIT tier (always runs)      -- the panel fragments, chart labels and pipeline invariants,
                                  straight from the source tree.
  SITE tier (needs a build)    -- structural checks of a BUILT site. Point ALC_SITE_DIR at the
                                  output of scripts/build_dashboard.py; skipped when unset, so CI
                                  without a build still runs the unit tier.

        ALC_SITE_DIR=/path/to/dash_prod python -m pytest tests/test_dashboard_site.py -q

  BROWSER tier (not automated here, run before a release) -- serve the site, open each station
  page in a real browser and check: zero console errors; the daily panel draws (>= 1 trace);
  the calendar shows ONE month and its arrows step months; #period-sel changes the x-range of
  fig-ts-*; arrows/0/1/2/3 act exactly once per keypress; a .diaglink click moves BOTH the image
  viewer and the panel; a file:// open shows the explicit fetch warning, not a blank panel.

Every test below is annotated with the regression that motivated it. If one of these fails,
something that already broke an operator once is broken again.
"""
from __future__ import annotations

import json
import os
import re
import sys
from base64 import b64decode
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

SITE = Path(os.environ["ALC_SITE_DIR"]) if os.environ.get("ALC_SITE_DIR") else None
needs_site = pytest.mark.skipif(SITE is None, reason="ALC_SITE_DIR not set (no built site)")


def _panel():
    import monitoring.panel as M
    return M


def _css_classes(css: str) -> set:
    out = set()
    for rule in re.split(r"}", css):
        sel = rule.split("{")[0]
        out |= set(re.findall(r"\.([a-zA-Z][\w-]*)", sel))
    return out


# =====================================================================================
# UNIT TIER
# =====================================================================================

def test_panel_css_is_fully_scoped():
    """Regression lineage: .msg collided with style.css (the 360px verdict box); then the audit
    found bare h2/:root/body/* rules restyling every host heading and hijacking the host CSS
    tokens -- and the old class-intersection test was DEFEATED by an allowlist and by matching
    classes only. The invariant now: every selector in the embedded sheet is scoped under
    .dp-wrap, and every custom property is namespaced --dp-*."""
    import re
    css = _panel().PANEL_CSS
    i = 0
    while True:
        j = css.find("{", i)
        if j < 0:
            break
        sel = css[i:j]
        sel = sel.split("}")[-1]                 # drop a preceding @media close
        sel = sel[sel.rfind("*/") + 2 if "*/" in sel else 0:].strip()
        if sel.startswith("@media"):
            i = j + 1
            continue
        if sel:
            for one in sel.split(","):
                assert one.strip().startswith(".dp-wrap"), f"unscoped selector: {one.strip()[:60]}"
        i = css.find("}", j) + 1
    assert ":root" not in css, "the embedded sheet must not touch :root"
    assert not re.search(r"--(line|ink|dim|ok|bad)\s*:", css), "un-namespaced custom property"
    assert "body" in _panel().STANDALONE_CSS, "standalone chrome must live in STANDALONE_CSS"


def test_panel_keyboard_stands_down_when_diag_owns_it():
    """Regression: PANEL_JS and diag.js both bound document keydown, so every arrow press
    navigated twice and the image viewer drifted apart from the panel."""
    js = _panel().PANEL_JS
    assert "window.__panelEmbedded" in js.split("keydown", 1)[1][:900], \
        "PANEL_JS keydown must stand down via the explicit embed handshake"
    tpl = (REPO / "monitoring/templates/station.html").read_text(encoding="utf-8")
    assert "window.__panelEmbedded = true" in tpl, "host must set the handshake before PANEL_JS"


def test_no_stale_wording_in_panel_fragments():
    """Regression: PNG-era jargon ('no figure drawn') and superseded labels survived in served
    pages long after being 'fixed', because the fix landed in one generation of the code and the
    page was built from another."""
    m = _panel()
    blob = m.PANEL_JS + m.PANEL_CSS + m.PANEL_BODY
    for phrase in ("no figure drawn", "before anything could be plotted",
                   "screened / not used", "where B accumulates",
                   "open it to see why", "no usable measurement", "unflagged profiles",
                   "Signal / range", "Calibrated signal"):
        assert phrase not in blob, f"stale wording resurfaced: {phrase!r}"


def test_series_labels_are_derived_not_asserted():
    """Regression: charts.series_timeseries hard-coded the string 'v2.0', so a v2.2 archive was
    presented as v2.0 — a provenance error. Labels must come from the data's version column."""
    src = (REPO / "monitoring/charts.py").read_text(encoding="utf-8")
    assert 'name=("v2.0"' not in src and '"v2.0 Kalman' not in src
    assert "def version_label" in src
    runner = (REPO / "scripts/run_network_calibration.py").read_text(encoding="utf-8")
    assert '"version", "message"]' in runner, "CSV_FIELDS must carry the algorithm version"


def test_clear_day_is_not_blamed_on_peak_shape():
    """Regression: a day with zero clouds was reported as flag -22 'peak not sharp' because the
    peak-shape filters run on every profile and their counter dominated. If the instrument saw no
    cloud base all day, the day is -1."""
    runner = (REPO / "scripts/run_network_calibration.py").read_text(encoding="utf-8")
    assert "if flag <= -20:" in runner and "flag = -1.0" in runner, \
        "the no-CBH override on cloud rejection attribution is gone"


def test_dailypanel_js_is_a_staged_asset():
    """Regression: dailypanel.js was written but never added to _VERSIONED_ASSETS, so it was never
    copied into the site and the panel silently did not load."""
    from monitoring import render
    assert "dailypanel.js" in render._VERSIONED_ASSETS
    assert (REPO / "monitoring/static/dailypanel.js").exists()


def test_panel_fragments_have_no_unfilled_tokens_after_render():
    """Regression: __META__ was substituted on the page after PANEL_JS had already been spliced in,
    so the browser hit a literal undefined __META__ identifier."""
    m = _panel()
    assert m.PANEL_JS.count("__META__") == 1, "render code substitutes exactly one __META__ token"
    filled = m.PANEL_JS.replace("__META__", "{}")
    assert "__" + "META" + "__" not in filled
    assert m.PANEL_JS.count("{") == m.PANEL_JS.count("}"), "unbalanced braces in PANEL_JS"


def test_payload_tool_offers_data_only_mode():
    """Regression: the payload generator wrote its own index.html/station_*.html into the
    production site directory, clobbering the summary page with a stub."""
    src = (REPO / "scripts/build_station_dashboard.py").read_text(encoding="utf-8")
    assert "--no-pages" in src and "args.no_pages" in src
    assert "index.pop(d.strftime" in src,         "_write_payloads must MERGE into the existing index, not rebuild from this run's window"


def test_flag_labels_agree_between_pipeline_and_dashboard():
    """The dashboard must not invent its own meaning for a pipeline flag."""
    from calibration.flags import FLAG_MEANINGS
    from monitoring import config as mcfg
    for f, label in FLAG_MEANINGS.items():
        got = mcfg.flag_label(f, "rayleigh")
        assert isinstance(got, str) and got, f"dashboard has no label for pipeline flag {f}"


# =====================================================================================
# SITE TIER — structural checks of a built site
# =====================================================================================

def _station_pages():
    return sorted((SITE / "stations").glob("*.html"))


@needs_site
def test_site_has_production_summary_not_a_stub():
    """Regression: a stub '<h2>Stations</h2>' list overwrote the production summary page."""
    idx = (SITE / "index.html").read_text(encoding="utf-8", errors="replace")
    assert "Calibration monitor" in idx, "index.html is not the production summary page"
    stray = list(SITE.glob("station_*.html"))
    assert not stray, f"stray generator pages in the site root: {[p.name for p in stray]}"


@needs_site
def test_payload_index_matches_files_on_disk():
    """Regression: an 11-day test run overwrote a 180-day _index.json, shrinking the calendar to
    one month and marking every other day 'no data'. Index and files must agree BOTH ways."""
    data = SITE / "data"
    if not data.exists():
        pytest.skip("site has no daily payloads")
    for idx_file in data.glob("*/_index.json"):
        key = idx_file.parent.name
        index = json.loads(idx_file.read_text(encoding="utf-8"))
        files = {f.stem for f in idx_file.parent.glob("2*.json")}
        indexed = {f"{ds}_{m}" for ds, by in index.items() for m in by}
        missing = sorted(indexed - files)[:5]
        orphans = sorted(files - indexed)[:5]
        assert not missing, f"{key}: indexed but no file: {missing}"
        assert not orphans, f"{key}: on disk but not indexed (stale index?): {orphans}"
        # "a daily plot in all conditions": every payload parses, none is an error, and an entry
        # without a figure is only legitimate when the night truly has nothing to draw (kind none).
        for ds, by in index.items():
            for meth, summ in by.items():
                payload = json.loads((idx_file.parent / f"{ds}_{meth}.json")
                                     .read_text(encoding="utf-8"))
                assert payload.get("kind") != "error", f"{key} {ds} {meth}: error payload"
                if not summ.get("has_fig"):
                    assert payload.get("kind") == "none",                         f"{key} {ds} {meth}: no figure but kind={payload.get('kind')}"


@needs_site
def test_embedded_payload_agrees_with_index_file():
    """The page embeds a boot payload built FROM _index.json; if they diverge the calendar and the
    arrows disagree with what fetches can actually return."""
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        m = re.search(r'<script id="payload" type="application/json">(.*?)</script>', html, re.S)
        if not m:
            continue                       # station without payloads: no panel, by design
        boot = json.loads(m.group(1))
        idx_file = SITE / "data" / boot["station"]["key"] / "_index.json"
        assert idx_file.exists(), f"{page.name}: panel embedded but no _index.json"
        index = json.loads(idx_file.read_text(encoding="utf-8"))
        assert set(boot["index"]) == set(index), \
            f"{page.name}: embedded index ({len(boot['index'])} days) != disk ({len(index)} days)"


@needs_site
def test_sampled_payloads_decode():
    """Regression class: quantised curtains whose b64 does not match the declared shape draw
    garbage or nothing. Sample a few per station rather than trusting all of them."""
    for idx_file in (SITE / "data").glob("*/_index.json"):
        index = json.loads(idx_file.read_text(encoding="utf-8"))
        picked = 0
        for ds in sorted(index, reverse=True):
            for meth, s in index[ds].items():
                if not s.get("has_fig") or picked >= 3:
                    continue
                f = idx_file.parent / f"{ds}_{meth}.json"
                p = json.loads(f.read_text(encoding="utf-8"))
                cur = p.get("curtain") or {}
                ny, nx = cur["shape"]
                assert len(b64decode(cur["b64"])) == ny * nx, f"{f.name}: b64 != shape"
                assert cur["lo"] < cur["hi"], f"{f.name}: degenerate colour limits"
                picked += 1
        assert picked, f"{idx_file.parent.name}: no decodable payload found"


@needs_site
def test_every_referenced_asset_exists():
    """Regression: dailypanel.js was referenced by the template but absent from assets/."""
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        for ref in re.findall(r'(?:src|href)="\.\./(assets/[^"?]+)', html):
            assert (SITE / ref).exists(), f"{page.name} references missing {ref}"


@needs_site
def test_no_template_or_wording_leftovers_in_built_pages():
    """Regression: served pages kept superseded wording because the build that carried the fix was
    never run against them. The forbidden list IS the changelog of past wording bugs."""
    forbidden = ("no figure drawn", "before anything could be plotted", "MOCKUP —",
                 "open it to see why", "no usable measurement", "unflagged profiles",
                 "no PNG involved", "Signal / range",
                 "{%", "{{ ")
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        for phrase in forbidden:
            assert phrase not in html, f"{page.name} still contains {phrase!r}"


@needs_site
def test_calendar_is_single_month():
    """Regression: six stacked month grids made the rail taller than the plots; the calendar must
    ship the one-month version with month navigation."""
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        if 'id="payload"' in html:
            assert "calnav" in html, f"{page.name}: month-navigation calendar missing"


def test_cloudless_profiles_do_not_vote_on_the_rejection_reason():
    """Regression: a cloudless day read 'flag -22 peak not sharp' / '-24 aerosol below cloud'
    because the funnel runs on every profile's PSEUDO-peak (strongest aerosol/noise gate) and those
    rejections out-voted the real ones. Cloudless profiles are normal sky: they tally under
    no_cloud_rejected and never decide the flag; only cloudy profiles' reasons do."""
    from types import SimpleNamespace
    import numpy as np
    from calibration.cloud._filters import apply_cloud_filters
    from calibration.flags import dominant_cloud_reject_flag

    rng = np.arange(0.0, 3100.0, 100.0)              # 31 gates, 100 m
    n = 6
    beta = np.full((rng.size, n), 1.0)               # flat "clear sky with aerosol" profiles
    cbh = np.full(n, np.nan)
    # profile 5: a genuine cloud at 1.5 km with heavy aerosol below -> must fail the RATIO filter
    beta[:, 5] = 0.001
    beta[1:11, 5] = 5.0                              # aerosol layer below cloud
    beta[15, 5] = 100.0                              # sharp cloud peak
    cbh[5] = 800.0
    data = SimpleNamespace(range=rng, cbh=cbh)
    cfg = SimpleNamespace(cal_minheight=100.0, cal_maxheight=2400.0, attenuation_factor=2.0,
                          ratio_filter=0.1, cbh_minheight=100.0, cbh_maxheight=2400.0)
    _, stats = apply_cloud_filters(np.ones(n), beta, data, cfg)
    assert stats["no_cloud_rejected"] == 5, stats
    assert stats["ratio_rejected"] == 1, stats
    assert stats["above_rejected"] == stats["below_rejected"] == stats["cbh_rejected"] == 0, stats

    flag, reason, _ = dominant_cloud_reject_flag(None, stats, None)
    assert flag == -24.0 and reason == "ratio_rejected", (flag, reason)
    # ... and with ONLY cloudless rejections, the day is simply "no liquid cloud"
    flag2, reason2, _ = dominant_cloud_reject_flag(None, {"no_cloud_rejected": 5}, None)
    assert flag2 == -1.0 and reason2 == "no liquid cloud", (flag2, reason2)


@needs_site
def test_availability_card_present_with_data():
    """Regression: the Payerne CL31 page skipped the 'Data availability, cloud cover &
    calibration' card entirely because its status.csv was missing — the section must always
    render, and on this site every stream must carry real heatmap rows."""
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        assert "Data availability, cloud cover" in html, f"{page.name}: availability card missing"
        # A stream whose status has not been recorded shows the explicit placeholder -- the card
        # may never silently vanish, and a page with data must carry real heatmap rows.
        if 'id="fig-avail"' in html:
            assert '"heatmap"' in html, f"{page.name}: availability card has no heatmap data"
            assert 'id="status-index"' in html, f"{page.name}: status index missing"
        else:
            assert "No instrument status has been recorded" in html,                 f"{page.name}: neither availability data nor the placeholder"


@needs_site
def test_bottom_products_present_and_ordered():
    """The three operational products must close every page, in the order OmB -> Cloudnet
    classification -> sensitivity, after the method blocks — present even when a product is
    absent for the stream (a section that vanishes reads as a broken page)."""
    for page in _station_pages():
        html = page.read_text(encoding="utf-8", errors="replace")
        i_omb = html.find('id="sec-omb"')
        i_cls = html.find('id="sec-classification"')
        i_sen = html.find('id="sec-sens"')
        assert -1 not in (i_omb, i_cls, i_sen), f"{page.name}: a product section is missing"
        assert i_omb < i_cls < i_sen, f"{page.name}: product sections out of order"
        last_method = html.rfind('class="methodblock"')
        assert i_omb > last_method, f"{page.name}: products are not at the bottom"
        for sec, nxt in (("sec-omb", i_cls), ("sec-classification", i_sen), ("sec-sens", None)):
            start = html.find(f'id="{sec}"')
            chunk = html[start: nxt if nxt else len(html)]
            assert ("data-src-all" in chunk or "diag-data" in chunk
                    or 'class="muted"' in chunk), f"{page.name}: {sec} has neither data nor notice"


def test_profile_views_are_raw_s_rcs_beta_in_that_order():
    """The selector reads Raw S -> Raw RCS -> Att. backscatter (the order of derivation); the
    default is beta when the night has it, Raw RCS otherwise. 'Signal / range^2' named the
    operation, not the quantity -- dividing RCS by r^2 UNDOES the range correction."""
    js = _panel().PANEL_JS
    assert "const VIEW_ORDER = ['raws', 'rcs', 'beta'];" in js
    assert "raws:'Raw S'" in js and "beta:'Att. backscatter'" in js
    assert "have.includes('beta') ? 'beta' : (have.includes('rcs')" in js


def test_ratio_tile_absent_on_single_method_streams():
    """A red em-dash 'CLOUD / RAYLEIGH RATIO' tile on a one-method stream reads as a fault; the
    tile must simply not exist there."""
    import pandas as pd
    from monitoring import metrics
    g = pd.DataFrame({"success": [1, 1, 1], "cal_value": [1.0, 1.1, 0.9],
                      "datetime": pd.to_datetime(["2026-08-01", "2026-08-02", "2026-08-03"])})
    tiles = metrics.cl_headline_tiles({"cloud": g})
    assert all("RATIO" not in t["label"] for t in tiles), tiles


def test_ctrl_arrows_jump_to_valid_calibrations():
    """Operator request: Ctrl+Left/Right must move to the previous/next night that produced a
    constant. Registered in the capture phase so it wins over diag.js's step-every-day binding."""
    js = (REPO / "monitoring/static/dailypanel.js").read_text(encoding="utf-8")
    assert "e.ctrlKey || e.metaKey" in js and "constant !== null" in js
    assert js.rstrip().count("}, true);") >= 1, "Ctrl handler must be a capture-phase listener"


def test_calendar_colour_classes_are_all_defined():
    """Regression: 'the reject days are gone' — merging the two rejected classes removed
    CAL_COL.none while a lookup still referenced it, so rejected days resolved undefined and
    rendered unstyled grey. Every CAL_COL.<x> reference must be a defined key, and hasFig must
    understand index summaries (has_fig), not only full payloads."""
    import re
    js = _panel().PANEL_JS
    block = js[js.index("const CAL_COL"): js.index("};", js.index("const CAL_COL"))]
    defined = set(re.findall(r"^\s*(\w+):\s*\{", block, re.M))
    used = set(re.findall(r"CAL_COL\.(\w+)", js))
    assert used <= defined, f"dangling CAL_COL references: {sorted(used - defined)}"
    assert "p.has_fig" in js, "hasFig must read the index summary field"


def test_calendar_rail_beside_the_card():
    """Operator request: the calendar stands OUTSIDE the daily card, collapsible, with the date and
    day arrows on the card's title row."""
    m = _panel()
    assert 'id="cal"' in m.PANEL_RAIL and 'id="cal"' not in m.PANEL_BODY
    assert "dp-title" in m.PANEL_BODY.split("</h2>")[0], "title must open the toolbar row"
    assert "cal-toggle" in m.PANEL_BODY and "dp-wrap" in m.PANEL_CSS
    assert "right:100%" in m.PANEL_CSS, "rail must float in the page margin, not steal card width"
    assert 'id="cst"' not in m.PANEL_BODY, "C_L belongs to the verdict line, not the title row"
    tpl = (REPO / "monitoring/templates/station.html").read_text(encoding="utf-8")
    assert "dp-rail" in tpl and "daily_panel.rail" in tpl
    assert "no PNG involved" not in tpl


def test_calendar_month_navigation_does_not_snap_back():
    """Regression: the month arrows flipped the calendar and it bounced straight back -- markCal
    (called at the end of every buildCal) snapped calMonth to the selected day's month
    unconditionally, so browsing any other month was impossible. The follow-the-day behaviour must
    trigger only when the DAY changed. Day-click wiring is asserted alongside; the click-through
    itself is the browser tier's job."""
    js = _panel().PANEL_JS
    assert "lastMarkedDate" in js, "markCal needs the day-changed guard"
    guard = js.index("curDate !== lastMarkedDate")
    snap = js.index("calMonth = curDate.slice(0, 6)")
    assert guard < snap, "the month snap-back must sit INSIDE the day-changed guard"
    assert "#cal-prev" in js and "#cal-next" in js, "month arrows missing"
    assert "cell.addEventListener('click'" in js, "day cells must be clickable"


def test_rangesync_never_touches_the_panel():
    """Audit-confirmed: #period-sel relayouted the panel's one-night figures with station-history
    date ranges, corrupting their axes. rangesync must exclude figures inside #daily."""
    js = (REPO / "monitoring/static/rangesync.js").read_text(encoding="utf-8")
    assert "#daily" in js and "filter" in js, "rangesync must skip the daily panel's figures"


def test_panel_module_is_production_code():
    """Audit-confirmed: production render.py imported panel fragments from scripts/mockup_*. The
    panel lives in monitoring/panel.py; the mockup path is a thin shim."""
    r = (REPO / "monitoring/render.py").read_text(encoding="utf-8")
    assert "from monitoring import panel as PANEL" in r
    shim = (REPO / "scripts/mockup_daily_panel.py").read_text(encoding="utf-8")
    assert "from monitoring.panel import" in shim and len(shim) < 1000
