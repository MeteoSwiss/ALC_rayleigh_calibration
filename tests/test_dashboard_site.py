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
    import scripts.mockup_daily_panel as M
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

def test_panel_css_does_not_collide_with_production_css():
    """Regression: the panel reused ``.msg``, which style.css caps at 360 px — the night's
    verdict box rendered as a narrow strip. Any selector defined by BOTH sheets is a landmine;
    the two below are shared ON PURPOSE (compatible box styling). Adding a name to this list
    requires checking both definitions side by side, not just appending it."""
    allowed_shared = {"card", "empty"}
    style = (REPO / "monitoring/static/style.css").read_text(encoding="utf-8")
    inter = _css_classes(_panel().PANEL_CSS) & _css_classes(style)
    assert inter <= allowed_shared, f"new CSS collision(s) with style.css: {sorted(inter - allowed_shared)}"


def test_panel_keyboard_stands_down_when_diag_owns_it():
    """Regression: PANEL_JS and diag.js both bound document keydown, so every arrow press
    navigated twice and the image viewer drifted apart from the panel."""
    js = _panel().PANEL_JS
    assert ".diag-data" in js.split("keydown", 1)[1][:600], \
        "PANEL_JS keydown must stand down when the production diag viewer is on the page"


def test_no_stale_wording_in_panel_fragments():
    """Regression: PNG-era jargon ('no figure drawn') and superseded labels survived in served
    pages long after being 'fixed', because the fix landed in one generation of the code and the
    page was built from another."""
    m = _panel()
    blob = m.PANEL_JS + m.PANEL_CSS + m.PANEL_BODY
    for phrase in ("no figure drawn", "before anything could be plotted",
                   "screened / not used", "where B accumulates"):
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
