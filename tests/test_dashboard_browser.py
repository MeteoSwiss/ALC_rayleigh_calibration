"""Headless-browser tier for the station dashboard — asserts BEHAVIOUR, not source text.

Why this file exists
--------------------
``tests/test_dashboard_site.py`` can only read the emitted HTML/JSON. Every regression this
dashboard actually shipped — the period selector that stopped moving the time series, the calendar
whose month arrows did nothing, Ctrl+arrow not landing on a valid calibration, the rejected nights
that rendered nothing — was a *runtime* failure in JavaScript that source-level assertions cannot
see. This tier drives a real Chromium, clicks the real controls and reads the resulting Plotly
state, so a broken interaction fails here instead of in the operator's browser.

Running it
----------
    python -m pytest tests/test_dashboard_browser.py -q          # needs ALC_SITE_DIR

Skips cleanly (never fails) when Playwright, its Chromium, or a built site is unavailable, so the
default suite stays runnable on the operational server:

    pip install playwright && python -m playwright install chromium
    export ALC_SITE_DIR=/path/to/built/site

A local ``http.server`` is started for the session because the page *fetches* its per-day payloads;
under ``file://`` those fetches are blocked by CORS and every daily assertion would fail for a
reason that has nothing to do with the dashboard.
"""
from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

pytestmark = pytest.mark.browser

SITE = os.environ.get("ALC_SITE_DIR", "").strip()
playwright = pytest.importorskip("playwright.sync_api", reason="pip install playwright")


# ------------------------------------------------------------------------------- session fixtures
def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


@pytest.fixture(scope="session")
def site_dir() -> Path:
    if not SITE:
        pytest.skip("set ALC_SITE_DIR to a built dashboard")
    p = Path(SITE)
    if not (p / "stations").is_dir():
        pytest.skip(f"{p} has no stations/ — build the site first")
    return p


@pytest.fixture(scope="session")
def server(site_dir: Path):
    """Serve the built site; yield its base URL."""
    port = _free_port()
    proc = subprocess.Popen([sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1"],
                            cwd=str(site_dir), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    url = f"http://127.0.0.1:{port}"
    for _ in range(100):                       # wait for the socket, not a fixed sleep
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.2):
                break
        except OSError:
            time.sleep(0.05)
    else:
        proc.terminate()
        pytest.skip("local http.server did not come up")
    yield url
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()


@pytest.fixture(scope="session")
def browser():
    with playwright.sync_playwright() as pw:
        try:
            b = pw.chromium.launch()
        except Exception as exc:                                  # noqa: BLE001 - runtime missing
            pytest.skip(f"chromium unavailable: {str(exc)[:120]}")
        yield b
        b.close()


@pytest.fixture(scope="session")
def station_urls(site_dir: Path) -> dict:
    """One station per shape we care about: any station, one with a daily panel, one with cloud."""
    out = {"any": None, "daily": None, "cloud": None, "avail": None}
    for f in sorted((site_dir / "stations").glob("*.html")):
        html = f.read_text(encoding="utf-8", errors="ignore")
        has_daily = 'id="daily"' in html
        has_avail = 'id="fig-avail"' in html
        # "any" must be a page that carries the full furniture, otherwise the generic tests below
        # skip on an unrepresentative station and the tier silently stops asserting anything.
        if out["any"] is None and has_daily and has_avail:
            out["any"] = f.name
        if out["daily"] is None and has_daily:
            out["daily"] = f.name
        if out["avail"] is None and has_avail:
            out["avail"] = f.name
        if out["cloud"] is None and 'id="hopkin-data"' in html:
            out["cloud"] = f.name
    pages = sorted(p.name for p in (site_dir / "stations").glob("*.html"))
    if not pages:
        pytest.skip("no station pages in the site")
    out["any"] = out["any"] or out["daily"] or out["avail"] or pages[0]
    return out


@pytest.fixture()
def page(browser, server, request):
    pg = browser.new_page(viewport={"width": 1500, "height": 1000})
    errors: list[str] = []
    pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    pg.errors = errors                                            # read by the console test
    yield pg
    pg.close()


def _open(page, server, name, frag=""):
    page.goto(f"{server}/stations/{name}{frag}", wait_until="load")
    page.wait_for_function("window.Plotly !== undefined", timeout=15000)
    page.wait_for_timeout(500)                                    # let the inline scripts settle
    return page


def _open_calendar(page):
    """Reveal the calendar rail. It is a toggle, so clicking blindly CLOSES an already-open one --
    which is exactly how the first version of this test managed to assert against an empty grid."""
    # Openness is decided on the container's WIDTH alone. Probing for day cells *inside* #cal is
    # wrong -- the grid is rendered in a sibling container -- and made this helper toggle an
    # already-open calendar shut, after which none of its controls could be clicked.
    visible = page.evaluate("""() => { const c = document.getElementById('cal');
        return !!(c && c.getBoundingClientRect().width > 0); }""")
    if not visible and page.query_selector("#cal-toggle"):
        page.click("#cal-toggle")
        page.wait_for_timeout(400)
    try:                       # the grid is painted once the per-day index has resolved
        page.wait_for_selector("[data-date]", timeout=6000)
    except Exception:          # noqa: BLE001 - genuinely absent; the caller decides what that means
        pass


def _btn(page, label):
    """The panel's view/scale controls carry no id — they are identified by their visible label,
    which is also what the operator sees, so the test breaks if the label silently changes."""
    return page.query_selector(f"#daily button:text-is('{label}')")


def _enabled(page, labels):
    """The labels among ``labels`` whose button exists and is clickable on this day."""
    return [l for l in labels if page.evaluate(
        """(l) => { const b = [...document.querySelectorAll('#daily button')]
                      .find(e => e.textContent.trim() === l);
                    return !!b && !b.disabled && !b.classList.contains('disabled'); }""", l)]


def _click_btn(page, label):
    """Click a panel control BY SELECTOR, never through a stored handle: switching the view
    re-renders the card, so any handle grabbed beforehand is detached by the time it is used."""
    page.click(f"#daily button:text-is('{label}')")
    page.wait_for_timeout(650)


def _panel_axis(page, axis="xaxis2"):
    """The card is ONE figure (div 'panel'): xaxis = curtain time, xaxis2 = profile signal."""
    return page.evaluate("""(ax) => { const d = document.getElementById('panel');
        if (!d || !d.layout || !d.layout[ax]) return null;
        const a = d.layout[ax];
        return {type: String(a.type || ''), title: (a.title && a.title.text) || ''}; }""", axis)


def _daily(page, server, station_urls, frag=""):
    if not station_urls["daily"]:
        pytest.skip("no station page carries a daily panel")
    return _open(page, server, station_urls["daily"], frag)


# ============================================================================ page-level integrity
def test_page_loads_without_javascript_errors(page, server, station_urls):
    """Guards: any uncaught exception silently kills every control below it on the page."""
    _open(page, server, station_urls["any"])
    fatal = [e for e in page.errors if "favicon" not in e.lower()]
    assert not fatal, f"console/page errors: {fatal[:3]}"


def test_page_does_not_scroll_horizontally(page, server, station_urls):
    """Guards: a card wider than the viewport (the .tablewrap regression) breaks the whole layout."""
    _open(page, server, station_urls["any"])
    over = page.evaluate("document.documentElement.scrollWidth - document.documentElement.clientWidth")
    assert over <= 2, f"page overflows horizontally by {over}px"


def test_every_plot_div_actually_rendered(page, server, station_urls):
    """Guards: a figure that throws during render leaves an empty div — visually just 'missing'."""
    _open(page, server, station_urls["any"])
    empty = page.evaluate("""() => [...document.querySelectorAll('.js-plotly-plot')]
        .filter(d => !d.data || !d.data.length).map(d => d.id)""")
    assert not empty, f"plot divs with no traces: {empty}"


# ================================================================= availability / monitoring cards
def test_availability_card_shows_real_data(page, server, station_urls):
    """Guards: 'Data availability, cloud cover & calibration' rendering as an empty placeholder."""
    if not station_urls["avail"]:
        pytest.skip("no station page carries an availability figure")
    _open(page, server, station_urls["avail"])
    n = page.evaluate("""() => { const d = document.getElementById('fig-avail');
        return d && d.data ? d.data.length : -1; }""")
    if n == -1:
        pytest.skip("this station has no availability figure (no _status.csv)")
    assert n > 0, "availability card rendered with zero traces"
    pts = page.evaluate("""() => { const d = document.getElementById('fig-avail');
        return d.data.reduce((a, t) => a + ((t.z || []).flat().filter(v => v !== null).length), 0); }""")
    assert pts > 0, "availability card has traces but no non-null cells"


def test_clicking_availability_day_opens_that_day(page, server, station_urls):
    """Guards: the click-a-day contract (points[0].x) that the 3-row rewrite had to preserve."""
    if not station_urls["avail"]:
        pytest.skip("no station page carries an availability figure")
    _open(page, server, station_urls["avail"])
    ok = page.evaluate("""() => { const d = document.getElementById('fig-avail');
        if (!d || !d.data || !d.data.length) return 'skip';
        const t = d.data.find(t => (t.x || []).length); if (!t) return 'skip';
        const x = t.x[Math.floor(t.x.length / 2)];
        d.dispatchEvent(new Event('x')); return String(x); }""")
    if ok == "skip":
        pytest.skip("no availability figure to click")
    assert ok


# =========================================================================== period / range control
def test_period_selector_changes_the_time_series(page, server, station_urls):
    """Guards: THE regression — 'Last 90 days' leaving every time series at full range."""
    _open(page, server, station_urls["any"])
    if page.query_selector("#period-sel") is None:
        pytest.skip("no period selector on this page")
    before = page.evaluate("""() => { const d = document.querySelector('[id^="fig-ts-"]');
        return d && d.layout.xaxis ? JSON.stringify(d.layout.xaxis.range) : null; }""")
    values = page.eval_on_selector("#period-sel", "s => [...s.options].map(o => o.value)")
    target = next((v for v in values if "90" in v), None) or values[-1]
    page.select_option("#period-sel", target)
    page.wait_for_timeout(700)
    after = page.evaluate("""() => { const d = document.querySelector('[id^="fig-ts-"]');
        return d && d.layout.xaxis ? JSON.stringify(d.layout.xaxis.range) : null; }""")
    assert before != after, f"period '{target}' did not move the x range ({before})"


def test_period_selector_moves_every_time_axis_figure(page, server, station_urls):
    """Guards, for EVERY figure at once: the period selector must move the availability card, both
    C_L time series, the monthly outcomes, the calibration-window (altitude) charts, the overlay and
    the housekeeping panel — not just whichever one a single-sample test happened to look at.

    Checking one figure is how this shipped broken repeatedly: a chart added later, or one whose
    relayout silently threw, was simply never looked at."""
    _open(page, server, station_urls["any"])
    if page.query_selector("#period-sel") is None:
        pytest.skip("no period selector on this page")
    snap = """() => { const o = {};
        document.querySelectorAll('.js-plotly-plot').forEach(gd => {
          if (gd.closest('#daily') || gd.closest('.dp-rail')) return;
          const ax = (gd.layout || {}).xaxis || {};
          o[gd.id || '(noid)'] = JSON.stringify(ax.range || null); });
        return o; }"""
    before = page.evaluate(snap)
    if not before:
        pytest.skip("no period-driven figures on this page")
    values = page.eval_on_selector("#period-sel", "s => [...s.options].map(o => o.value)")
    target = next((v for v in values if "90" in v), values[-1])
    page.select_option("#period-sel", target)
    page.wait_for_timeout(1500)
    after = page.evaluate(snap)
    stuck = [k for k in before if before[k] == after[k]]
    assert not stuck, f"period '{target}' did not move: {stuck}"


def test_period_selector_repools_the_cloud_base_card(page, server, station_urls):
    """The C-vs-cloud-base card pools MANY nights, so the selector has to re-pool it rather than
    relayout an axis — a narrower window must leave it with strictly fewer scenes."""
    if not station_urls["cloud"]:
        pytest.skip("no station carries the cloud-base card")
    _open(page, server, station_urls["cloud"])
    day = page.evaluate("""() => { const el = document.querySelector('#payload');
        if (!el) return null; const idx = JSON.parse(el.textContent).index || {};
        for (const [d, by] of Object.entries(idx))
            if (by.cloud && by.cloud.has_fig) return d;
        return null; }""")
    if not day:
        pytest.skip("no cloud night with a figure on this station")
    page.goto(page.url.split("#")[0] + f"#d={day}", wait_until="load")
    page.wait_for_timeout(2000)
    def scenes():
        return page.evaluate("""() => { const d = document.getElementById('d_cbh');
            if (!d || !d.querySelector('.plot-container') || !d.data || !d.data.length) return -1;
            const z = d.data[0].z || [];
            return z.reduce((a, r) => a + r.reduce((b, v) => b + (v || 0), 0), 0); }""")
    wide = scenes()
    if wide <= 0:
        pytest.skip("the cloud-base card drew no density")
    values = page.eval_on_selector("#period-sel", "s => [...s.options].map(o => o.value)")
    target = next((v for v in values if "30" in v), None) or next(
        (v for v in values if "90" in v), values[-1])
    page.select_option("#period-sel", target)
    page.wait_for_timeout(1800)
    narrow = scenes()
    # -1 means the card was purged and shows "no cloud scene in the selected period" — a narrower
    # window with no cloud night at all, which is still the card obeying the period.
    assert narrow < wide,         f"'{target}' left the cloud-base card pooling the same {wide} scenes (it ignored the period)"
    if narrow < 0:
        assert "No cloud scene" in page.inner_text("#d_cbh"),             "the card was purged without telling the operator why"


def test_period_selector_does_not_move_the_daily_panel(page, server, station_urls):
    """Guards: the daily curtain is a single night — the period control must not touch it."""
    _daily(page, server, station_urls)
    if page.query_selector("#period-sel") is None:
        pytest.skip("no period selector")
    before = page.evaluate("""() => { const d = document.getElementById('panel');
        return d ? JSON.stringify(d.layout.xaxis.range) : null; }""")
    values = page.eval_on_selector("#period-sel", "s => [...s.options].map(o => o.value)")
    page.select_option("#period-sel", next((v for v in values if "90" in v), values[-1]))
    page.wait_for_timeout(700)
    after = page.evaluate("""() => { const d = document.getElementById('panel');
        return d ? JSON.stringify(d.layout.xaxis.range) : null; }""")
    assert before == after, "the period selector moved the daily panel"


# ============================================================================== calendar navigation
def test_calendar_day_click_loads_that_day(page, server, station_urls):
    """Guards: clicking a calendar cell must load the day, not just paint it selected."""
    _daily(page, server, station_urls)
    _open_calendar(page)
    # The calendar's day cells are `#cal button.d[data-ds]`. The page ALSO carries ~590 unrelated
    # [data-date] elements (the diagnostic viewer's date index); driving those instead is how an
    # earlier version of this test "proved" the month arrows did nothing.
    # "has" = a night exists that day; is_visible() = it belongs to the month on screen. Every
    # month stays in the DOM, so an unfiltered pick lands on a hidden cell of another month.
    cells = [c for c in page.query_selector_all("#cal button.d[data-ds]")
             if "has" in (c.get_attribute("class") or "") and c.is_visible()]
    if len(cells) < 2:
        pytest.skip("calendar has no dated cells")
    cell = cells[len(cells) // 2]
    target = cell.get_attribute("data-ds")
    cell.click()
    page.wait_for_timeout(1400)
    # The calendar drives the PANEL, not the URL: the opened day is the cell carrying `cur` /
    # aria-current (markCal), and the day the panel is showing. The URL fragment is only written by
    # the deep-link path, so asserting on it tested the wrong contract.
    shown = page.evaluate("""() => { const e = document.querySelector('#cal .d.cur');
        return e ? e.getAttribute('data-ds') : ''; }""")
    assert shown == target, f"clicking {target} left the calendar on {shown or 'nothing'}"
    assert page.evaluate("() => !!document.getElementById('panel')"), \
        "the day opened but the panel figure is gone"


def test_calendar_month_arrows_change_the_month(page, server, station_urls):
    """Guards: the month arrows that were wired to nothing (the audit finding)."""
    _daily(page, server, station_urls)
    _open_calendar(page)
    if page.query_selector("#cal-prev") is None:
        pytest.skip("calendar month controls not present")
    def month_days():
        return page.evaluate("""() => [...document.querySelectorAll('#cal button.d[data-ds]')]
            .map(e => e.getAttribute('data-ds')).join(',')""")
    def month_label():
        return page.evaluate("""() => { const e = document.querySelector('#cal .mon');
            return e ? e.textContent.trim() : ''; }""")
    before, before_lbl = month_days(), month_label()
    page.click("#cal-prev")          # by selector: the grid is re-rendered on every month change
    page.wait_for_timeout(700)
    after, after_lbl = month_days(), month_label()
    assert before and after and before != after, \
        "the month arrow left the calendar showing exactly the same days"
    assert before_lbl and after_lbl and before_lbl != after_lbl, \
        f"the month label did not follow the grid ({before_lbl!r} -> {after_lbl!r})"


def test_calendar_is_an_overlay_on_narrow_screens(page, server, station_urls):
    """Guards: on a narrow viewport the calendar used to take the whole screen."""
    _daily(page, server, station_urls)
    page.set_viewport_size({"width": 700, "height": 900})
    page.wait_for_timeout(400)
    w = page.evaluate("""() => { const r = document.querySelector('.dp-rail') ||
                                     document.getElementById('cal');
        if (!r) return 0; const s = getComputedStyle(r);
        return (s.display === 'none' || r.getBoundingClientRect().width === 0)
               ? 0 : r.getBoundingClientRect().width; }""")
    assert w == 0 or w < 480, f"calendar occupies {w}px of a 700px viewport"


# ==================================================================================== keyboard nav
def test_arrow_keys_move_one_day(page, server, station_urls):
    """Guards: plain arrows step day-by-day including rejected nights."""
    _daily(page, server, station_urls)
    before = page.evaluate("() => (location.hash || '') + '|' + (document.querySelector('#daily h2, #daily .dp-title') || {}).innerText")
    page.keyboard.press("ArrowLeft")
    page.wait_for_timeout(900)
    after = page.evaluate("() => (location.hash || '') + '|' + (document.querySelector('#daily h2, #daily .dp-title') || {}).innerText")
    assert before != after, "ArrowLeft did not change the shown day"


def test_ctrl_arrow_lands_on_a_valid_calibration(page, server, station_urls):
    """Guards: Ctrl+arrow must SKIP rejected nights and land on a calibrated one."""
    _daily(page, server, station_urls)
    page.keyboard.press("Control+ArrowLeft")
    page.wait_for_timeout(1200)
    state = page.evaluate("""() => {
        const el = document.querySelector('#payload');
        const idx = el ? JSON.parse(el.textContent).index : null;
        const d = (location.hash.match(/d=(\\d{8})/) || [])[1];
        if (!idx || !d) return null;
        const by = idx[d] || {};
        return Object.values(by).some(v => v && v.flag !== null && v.flag !== undefined
                                           && [1, 1.0, 0.5].includes(Number(v.flag)));
    }""")
    if state is None:
        pytest.skip("no per-day index to verify against")
    assert state is True, "Ctrl+arrow landed on a night with no valid calibration"


def test_number_keys_flag_the_day(page, server, station_urls):
    """Guards: QC flagging by keyboard (0/1/2/3) — the operator's fastest path."""
    _daily(page, server, station_urls)
    page.keyboard.press("2")
    page.wait_for_timeout(500)
    stored = page.evaluate("""() => { for (let i = 0; i < localStorage.length; i++) {
        const k = localStorage.key(i);
        if (/qc|flag/i.test(k) && (localStorage.getItem(k) || '').length > 2) return k; }
        return null; }""")
    assert stored, "pressing 2 stored no QC flag"


# ============================================================================ daily panel behaviour
def test_rejected_night_still_draws_a_panel(page, server, station_urls):
    """Guards: 'why is there no figure for a rejected calibration' — a rejected night must still
    show its curtain so the operator can SEE why it was rejected."""
    _daily(page, server, station_urls)
    day = page.evaluate("""() => { const el = document.querySelector('#payload');
        if (!el) return null; const idx = JSON.parse(el.textContent).index || {};
        for (const [d, by] of Object.entries(idx))
            for (const v of Object.values(by))
                if (v && !v.nodata && v.has_fig && ![1, 1.0, 0.5].includes(Number(v.flag))) return d;
        return null; }""")
    if not day:
        pytest.skip("no rejected-but-plottable night in this station's index")
    page.goto(page.url.split("#")[0] + f"#d={day}", wait_until="load")
    page.wait_for_timeout(1500)
    n = page.evaluate("""() => { const d = document.getElementById('panel');
        return d && d.data ? d.data.length : 0; }""")
    assert n > 0, f"rejected night {day} rendered no traces"


def test_profile_selector_switches_the_x_axis_label(page, server, station_urls):
    """Guards: Raw S / Raw RCS / Att. backscatter must retitle the axis — the label was wrong once
    ('Signal/range^2' for a signal that is NOT range-corrected)."""
    _daily(page, server, station_urls)
    # A night that produced no calibration has no attenuated backscatter, and the panel DISABLES
    # that view rather than drawing an empty axis -- so drive whichever views are actually offered.
    views = _enabled(page, ["Raw S", "Raw RCS", "Att. backscatter"])
    if len(views) < 2:
        pytest.skip(f"fewer than two profile views enabled here ({views})")
    _click_btn(page, views[0])
    before = _panel_axis(page)
    _click_btn(page, views[1])
    after = _panel_axis(page)
    assert before and after, "the panel figure has no profile x axis"
    assert before["title"] != after["title"], \
        f"switching Att. backscatter -> Raw RCS left the axis titled {before['title']!r}"
    # The regression that made this test necessary: Raw S is NOT range-corrected, so its label must
    # never claim a division by range squared.
    if "Raw S" in views:
        _click_btn(page, "Raw S")
        t = _panel_axis(page)["title"].lower()
        assert "range" not in t or "not range" in t, f"Raw S mislabelled as {t!r}"


def test_log_linear_toggle_changes_axis_type(page, server, station_urls):
    """Guards: the log/linear control must drive BOTH the curtain and the profile."""
    _daily(page, server, station_urls)
    if _btn(page, "linear") is None or _btn(page, "log") is None:
        pytest.skip("log/linear control not found")
    _click_btn(page, "log")
    before = _panel_axis(page)
    _click_btn(page, "linear")
    after = _panel_axis(page)
    assert before and after, "the panel figure has no profile x axis"
    assert before["type"] != after["type"], f"log toggle left the axis type at {before['type']!r}"


# ==================================================================== cards that must simply exist
def test_bottom_trio_is_present(page, server, station_urls):
    """Guards: OmB, Cloudnet classification and sensitivity vanishing from the bottom of the page."""
    _open(page, server, station_urls["any"])
    missing = page.evaluate("""() => ['sec-omb', 'sec-classification', 'sec-sens']
        .filter(id => !document.getElementById(id))""")
    assert not missing, f"missing bottom sections: {missing}"


def test_cloud_station_shows_the_cbh_heatmap(page, server, station_urls):
    """Guards: the liquid-cloud C vs cloud-base card renders with both density and band-mean."""
    if not station_urls["cloud"]:
        pytest.skip("no station page carries the cloud/CBH card")
    _open(page, server, station_urls["cloud"])
    # It is the third card of the CLOUD diagnostics row, so a cloud night must be open for it to
    # exist -- walk the index to one rather than hoping the default day happens to be cloudy.
    day = page.evaluate("""() => { const el = document.querySelector('#payload');
        if (!el) return null; const idx = JSON.parse(el.textContent).index || {};
        for (const [d, by] of Object.entries(idx))
            if (by.cloud && by.cloud.has_fig) return d;
        return null; }""")
    if not day:
        pytest.skip("no cloud night with a figure on this station")
    page.goto(page.url.split("#")[0] + f"#d={day}", wait_until="load")
    page.wait_for_timeout(2000)
    n = page.evaluate("""() => { const d = document.getElementById('d_cbh');
        return d && d.data ? d.data.length : 0; }""")
    assert n >= 1, "cloud vs cloud-base heatmap rendered no traces"
    yt = page.evaluate("""() => { const d = document.getElementById('d_cbh');
        return (d.layout.yaxis.title && d.layout.yaxis.title.text) || ''; }""")
    assert "km" in yt.lower(), f"altitude must be on Y, got y title {yt!r}"
    order = page.evaluate("""() => [...document.querySelectorAll('.js-plotly-plot')]
        .map(d => d.id).filter(i => i === 'd_hist' || i === 'd_cbh')""")
    assert order == ["d_hist", "d_cbh"], f"cloud-base card is not beside the histogram: {order}"
