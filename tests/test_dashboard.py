"""Regression tests for the monitoring dashboard.

Builds a tiny dashboard from synthetic fixtures (a couple of stations with cal /
omb / sens CSVs + dummy figures + an operational-constant CSV) and asserts the
key functionalities keep working across versions:

* the SQLite index builds from the per-stream CSVs + manifest,
* the full site builds (index page + one page per station),
* every summary map is present and the value maps carry data,
* the new OmB-bias and ICAO-altitude maps render with per-instrument-type symbols,
* the "% of operational constant" map is populated when an opcoeff CSV is given,
* per-station pages embed the OmB + sensitivity figures and the diagnostic calendar,
* the sensitivity period aggregation uses the MEDIAN, and the OmB figure labels the
  calibration "v2" (not "ours"), and the ICAO colorbar is the short label.

The tests avoid the heavy calibration/CAMS pipeline: all inputs are tiny and
synthetic, so the suite runs in a second or two.
"""
from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from monitoring import charts, config, index, render


# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------
def _png_bytes() -> bytes:
    """A minimal valid 1x1 PNG (so _stage_ombsens_pngs has a real file to copy)."""
    def chunk(typ, data):
        c = typ + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
    sig = b"\x89PNG\r\n\x1a\n"
    ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    idat = zlib.compress(b"\x00\xff\xff\xff")
    return sig + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")


_STATIONS = [
    dict(wmo="0-20000-0-00001", ident="A", type="CHM15k", site="ALPHA", lat=47.5, lon=8.5, alt=500.0),
    dict(wmo="0-20000-0-00002", ident="C", type="CL61", site="BETA", lat=46.8, lon=6.9, alt=490.0),
    dict(wmo="0-20000-0-00003", ident="B", type="CL31", site="GAMMA", lat=52.2, lon=13.4, alt=80.0),
]
_DATES = ["20260110", "20260115", "20260220", "20260315", "20260410", "20260512"]


@pytest.fixture
def dash(tmp_path):
    """Build a tiny dashboard; return paths + the loaded frames for assertions."""
    fullcal = tmp_path / "fullcal"
    for s in _STATIONS:
        key = f"{s['wmo']}_{s['ident']}"
        d = fullcal / key
        d.mkdir(parents=True)
        method = "rayleigh" if s["type"] in ("CHM15k", "CL61") else "cloud"
        theo = config.theoretical_cl(s["type"]) or 1.0
        cal_rows = ["date,method,flag,cal_value,uncertainty,n_profiles,bottom_height,top_height,message"]
        for i, dt in enumerate(_DATES):
            cl = theo * (0.9 + 0.02 * i)            # near the theoretical value
            cal_rows.append(f"{dt},{method},1,{cl:.6e},{cl*0.05:.3e},60,2000,6000,OK")
        (d / f"{key}_cal.csv").write_text("\n".join(cal_rows) + "\n")
        # per-night diagnostic PNGs so the station-page calendar (b.diags) renders
        pdir = d / "plots" / s["wmo"] / "2026"
        pdir.mkdir(parents=True, exist_ok=True)
        for dt in _DATES:
            (pdir / f"{dt}_{s['wmo']}_{method}_diag_compact.png").write_bytes(_png_bytes())
        # OmB + sensitivity summaries (one row each) + dummy figures
        (d / f"{key}_omb.csv").write_text(
            "date_start,date_end,wavelength,median_bias_ours,median_bias_op,"
            "median_bias_ours_wv,rms_ours,n_obs\n"
            f"20260101,20260531,910.5,-1.1e-07,-2.6e-07,4.0e-08,1.7e-06,2000\n")
        (d / f"{key}_sens.csv").write_text(
            "date_start,date_end,wavelength,icao_alt_200,icao_alt_2000,icao_alt_4000,"
            "sigma_night_3000,n_days_night,n_days_day\n"
            f"20250101,20261231,910.5,8500,12000,15000,5.5e-02,300,300\n")
        (d / f"{key}_omb.png").write_bytes(_png_bytes())
        (d / f"{key}_sens.png").write_bytes(_png_bytes())
        # decoded per-day status (for the Cloudnet availability bar) — mix of quality classes
        quals = ["pass", "warning", "error", "pass", "pass", "warning"]
        st_rows = ["date,quality,n_profiles,coverage_pct,gap_hours,n_alarm_h,n_warn_h,flags_json,summary"]
        for dt, q in zip(_DATES, quals):
            summ = {"pass": "No warning or error recorded", "warning": "'Window contamination (W)' 6 h",
                    "error": "'Transmitter failure (A)' 3 h"}[q]
            na, nw = (3, 0) if q == "error" else ((0, 6) if q == "warning" else (0, 0))
            st_rows.append(f'{dt},{q},400,100.0,0,{na},{nw},"{{}}","{summ}"')
        (d / f"{key}_status.csv").write_text("\n".join(st_rows) + "\n")
        # a per-year calibration NetCDF (for the download link + nc/ staging); content irrelevant here
        ncd = d / "2026"
        ncd.mkdir(parents=True, exist_ok=True)
        (ncd / f"ALC_calibration_{s['wmo']}_{s['ident']}2026.nc").write_bytes(b"CDF\x01dummy")

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(_STATIONS))

    # CEDA L2 links map {key: url} for the per-station "CEDA L2 data" link
    ceda = tmp_path / "ceda_links.json"
    ceda.write_text(json.dumps({f"{s['wmo']}_{s['ident']}":
                                f"https://data.ceda.ac.uk/badc/eprofile/data/daily_files/x/y/z_{s['ident']}"
                                for s in _STATIONS}))

    # operational-constant CSV (key,date,op_coeff) for the % of operational map
    op_rows = ["key,date,op_coeff"]
    for s in _STATIONS:
        key = f"{s['wmo']}_{s['ident']}"
        theo = config.theoretical_cl(s["type"]) or 1.0
        for dt in _DATES:
            op_rows.append(f"{key},{dt},{theo*0.8:.6e}")   # op constant 80% of theo
    opcoeff = tmp_path / "opcoeff.csv"
    opcoeff.write_text("\n".join(op_rows) + "\n")

    db = tmp_path / "index.sqlite"
    stats = index.build_index(fullcal, manifest, db)
    out = tmp_path / "site"
    render.build_site(db, out, opcoeff_csv=opcoeff, fullcal_dir=fullcal, ceda_links=ceda)
    return dict(out=out, db=db, fullcal=fullcal, opcoeff=opcoeff, ceda=ceda, stats=stats)


# ----------------------------------------------------------------------------
# End-to-end build
# ----------------------------------------------------------------------------
def test_index_built(dash):
    assert dash["stats"]["n_stations"] == len(_STATIONS)
    assert dash["db"].exists()


def test_site_pages_exist(dash):
    out = dash["out"]
    assert (out / "index.html").is_file()
    for s in _STATIONS:
        assert (out / "stations" / f"{s['wmo']}_{s['ident']}.html").is_file()


def test_summary_has_all_maps(dash):
    html = (dash["out"] / "index.html").read_text(encoding="utf-8")
    for div in ("fig-map", "fig-map-theo", "fig-map-op", "fig-map-omb", "fig-map-icao"):
        assert div in html, f"missing map {div} on summary page"


def test_op_map_populated(dash):
    """% of operational constant map must carry data when an opcoeff CSV is given
    (regression for the 'op map does not work' bug = missing opcoeff)."""
    html = (dash["out"] / "index.html").read_text(encoding="utf-8")
    # the op map's "no data" fallback title must NOT appear for fig-map-op
    assert "operational constant" in html
    # at least one finite op ratio reached keystats
    cal, series, st, kal, diag = _frames(dash)
    ratios = render._opcoeff_ratios(cal, render._load_opcoeff(dash["opcoeff"]), st)
    assert ratios["op_pct_op"].notna().any(), "op_pct_op all NaN -> op map empty"


def test_station_page_has_availability_bar_and_status(dash):
    """Cloudnet-style daily availability bar + per-day status field + status index are on the page."""
    page = (dash["out"] / "stations" / f"{_STATIONS[0]['wmo']}_{_STATIONS[0]['ident']}.html").read_text(encoding="utf-8")
    assert "fig-avail" in page                 # the availability bar figure
    assert "availcard" in page                 # its card (class="card availcard")
    assert 'id="status-index"' in page         # per-day status lookup for the viewer
    assert "diag-status" in page               # the reporting-status field under the image
    assert "'Transmitter failure (A)' 3 h" in page or "Window contamination" in page  # decoded summary


def test_station_page_has_download_and_ceda_links(dash):
    """Top-of-page links: a calibration-NetCDF download + the CEDA L2 archive link; nc/ is staged."""
    key = f"{_STATIONS[0]['wmo']}_{_STATIONS[0]['ident']}"
    page = (dash["out"] / "stations" / f"{key}.html").read_text(encoding="utf-8")
    assert "Calibration NetCDF" in page
    assert "data.ceda.ac.uk" in page and "CEDA L2 data" in page
    assert (dash["out"] / "nc" / key / f"ALC_calibration_{_STATIONS[0]['wmo']}_{_STATIONS[0]['ident']}2026.nc").exists()


def test_summary_reactive_kpis_and_instrument_chart(dash):
    """Summary embeds the per-series + activity JSON and the stacked instrument chart, and the KPI
    divs carry ids so filter.js can update them on country/type change."""
    html = (dash["out"] / "index.html").read_text(encoding="utf-8")
    assert "fig-instr" in html                                  # instruments-over-time chart
    assert 'id="series-index"' in html and 'id="instr-activity"' in html
    # "method series" was dropped from the KPI row (operator request 2026-08-17)
    assert 'id="kpi-series"' not in html
    for kid in ("kpi-instruments", "kpi-cals", "kpi-success"):
        assert f'id="{kid}"' in html, f"missing KPI id {kid}"
    # filter.js recomputes those ids
    js = (Path(charts.__file__).resolve().parent / "static" / "filter.js").read_text(encoding="utf-8")
    assert "recomputeKpis" in js and "recomputeInstrumentChart" in js


def test_favicon_present(dash):
    """The EUMETNET favicon is staged and referenced in the page head."""
    assert (dash["out"] / "assets" / "favicon.png").exists()
    page = (dash["out"] / "index.html").read_text(encoding="utf-8")
    assert 'rel="icon"' in page and "favicon.png" in page


def _status_frame():
    return pd.DataFrame({"date": ["20260101", "20260103", "20260110"],  # a gap between 03 and 10
                         "quality": ["pass", "error", "warning"],
                         "summary": ["No warning or error recorded",
                                     "'Transmitter failure (A)' 2 h",
                                     "'Window contamination (W)' 5 h"]})


def test_availability_rows_status_only():
    """With no calibration rows the card is the status row alone, spanning the whole record."""
    fig = charts.daily_availability_rows(_status_frame())
    assert fig is not None and len(fig.data) == 1
    assert len(fig.data[0].x) == 10                    # 20260101..20260110 inclusive
    assert list(fig.layout.yaxis.ticktext) == ["Instrument status"]
    assert fig.layout.title.text is None               # the HTML card carries the title
    q = list(fig.data[0].z[0])
    assert q[0] == 0 and q[2] == 2 and q[3] == 3       # pass, error, then gap days -> nodata


def test_availability_rows_adds_cloud_and_calibration_rows():
    """A CL61-like stream gets four rows: status, cloud cover, and one row per method."""
    status = _status_frame()
    status["mean_cloud_cover"] = [0.0, 4.0, 8.0]
    status["cloud_src"] = ["cloud_amount"] * 3
    cal = pd.DataFrame({"date": ["20260101", "20260103", "20260110"] * 2,
                        "method": ["rayleigh"] * 3 + ["cloud"] * 3,
                        "flag": [1.0, -1.0, -2.0, 1.0, -24.0, -21.0]})
    fig = charts.daily_availability_rows(status, cal, ["rayleigh", "cloud"])
    assert list(fig.layout.yaxis.ticktext) == [
        "Instrument status", "Mean cloud cover", "Calibration — Rayleigh",
        "Calibration — Liquid-cloud"]
    assert len(fig.data) == 4
    # rows are stacked top-down: y0 increases, and the axis range is descending
    assert [t.y0 for t in fig.data] == [0, 1, 2, 3]
    assert fig.layout.yaxis.range[0] > fig.layout.yaxis.range[1]
    order = config.CAL_CLASS_ORDER
    ray = list(fig.data[2].z[0])
    assert ray[0] == order.index("ok") and ray[2] == order.index("noscene")
    cld = list(fig.data[3].z[0])
    assert cld[9] == order.index("instrument")         # flag -21 = laser energy too low


def test_station_index_payload(dash):
    """data/stations.json carries one record per station, in the navigation order, with the
    status and the constant the neighbour links display."""
    p = dash["out"] / "data" / "stations.json"
    assert p.exists()
    recs = json.loads(p.read_text(encoding="utf-8"))
    assert [r["k"] for r in recs] == [f"{s['wmo']}_{s['ident']}" for s in _STATIONS]
    a = recs[0]
    assert a["t"] == "CHM15k" and a["n"] == "ALPHA" and a["w"] == "0-20000-0-00001"
    assert a["q"] == "warning"                       # last row of the synthetic status CSV
    m = a["m"]["rayleigh"]
    assert m["d"] == _DATES[-1] and m["f"] == 1.0
    # the constant is also expressed as a percent of the type's nominal value, which is what the
    # links show (raw C_L spans 11 orders of magnitude across types)
    assert 90 <= m["pct"] <= 110
    # the CL31 station calibrates by cloud, so its record carries that method instead
    assert "cloud" in recs[2]["m"]


def test_station_pages_fetch_the_index_and_defer_navigation(dash):
    """Station pages ship the filters + an EMPTY prev/next (stationnav.js fills them in), and the
    index URL carries a cache-busting token."""
    page = (dash["out"] / "stations" / f"{_STATIONS[1]['wmo']}_{_STATIONS[1]['ident']}.html").read_text(
        encoding="utf-8")
    assert 'id="f-country"' in page and 'id="f-type"' in page
    assert 'data-key="0-20000-0-00002_C"' in page
    assert 'data-prev=""' in page and 'data-next=""' in page      # filled client-side
    assert "data/stations.json?v=" in page
    for asset in ("stationindex.js", "stationnav.js"):
        assert asset in page and (dash["out"] / "assets" / asset).exists()


def test_station_index_url_token_tracks_content(dash, tmp_path):
    """The ?v= token is the index's own content hash, so an unchanged index keeps its cached URL."""
    page = (dash["out"] / "index.html").read_text(encoding="utf-8")
    tok = page.split("data/stations.json?v=")[1].split('"')[0]
    import hashlib
    body = (dash["out"] / "data" / "stations.json").read_text(encoding="utf-8")
    assert tok == hashlib.md5(body.encode("utf-8")).hexdigest()[:8]


def test_cal_class_groups_flags():
    """Every documented flag maps to a class, and an unknown flag never reads as a success."""
    assert config.cal_class(1) == "ok" and config.cal_class(0.5) == "ok"
    assert config.cal_class(-1) == "noscene"
    assert config.cal_class(-24) == "atmos" and config.cal_class(-2) == "atmos"
    assert config.cal_class(-20) == "instrument"
    assert config.cal_class(-3) == "retrieval"
    assert config.cal_class(0) == "nodata"
    assert config.cal_class(-77) == "retrieval"        # unknown -> never "ok"
    assert config.cal_class(None) == "nodata"


def test_monitoring_timeseries_uses_compact_x_encoding():
    """The HK figure ships x0/dx (not a per-trace datetime array) and carries no Plotly title."""
    hk = pd.DataFrame({"datetime": pd.to_datetime(["20260101", "20260102", "20260104"],
                                                  format="%Y%m%d"),
                       "laser": [100.0, 99.0, 98.0], "window": [80.0, 81.0, 82.0]})
    fig = charts.monitoring_timeseries(hk)
    assert fig.data and fig.data[0].x is None and fig.data[0].x0 is not None
    assert fig.data[0].dx == 86400000.0                # one day in ms
    assert len(fig.data[0].y) == 4                     # reindexed onto the full grid (a gap on 03)
    assert fig.layout.title.text is None


def test_instrument_count_chart_builder():
    """instrument_count_over_time stacks per-type monthly instrument counts."""
    activity = [
        {"key": "a", "itype": "CHM15k", "country": "CH", "months": ["202601", "202602"]},
        {"key": "b", "itype": "CHM15k", "country": "DE", "months": ["202601"]},
        {"key": "c", "itype": "CL61", "country": "CH", "months": ["202602"]},
    ]
    fig = charts.instrument_count_over_time(activity)
    names = [t.name for t in fig.data]
    assert "CHM15k" in names and "CL61" in names
    chm = next(t for t in fig.data if t.name == "CHM15k")
    assert list(chm.y) == [2, 1]   # Jan: a+b, Feb: a


def test_station_page_has_ombsens_and_calendar(dash):
    page = (dash["out"] / "stations" / f"{_STATIONS[0]['wmo']}_{_STATIONS[0]['ident']}.html").read_text(encoding="utf-8")
    assert "Observation" in page and "Background" in page          # OmB section
    assert "sensitivity" in page.lower()                            # sensitivity section
    assert 'class="diag-cal"' in page                               # calendar mount point
    assert "diag.js" in page                                        # calendar script loaded


# ----------------------------------------------------------------------------
# Chart builders (unit)
# ----------------------------------------------------------------------------
def _keystats():
    return pd.DataFrame([
        dict(key="0-20000-0-00001_A", lat=47.5, lon=8.5, itype="CHM15k", country="CH",
             name="ALPHA", n_dates=100, success_rate=70.0, methods="Rayleigh",
             op_pct_theo=95.0, op_pct_op=120.0, omb_bias=-0.11, icao_alt=8500.0),
        dict(key="0-20000-0-00002_C", lat=46.8, lon=6.9, itype="CL61", country="CH",
             name="BETA", n_dates=80, success_rate=60.0, methods="Rayleigh + Liquid-cloud",
             op_pct_theo=90.0, op_pct_op=110.0, omb_bias=0.05, icao_alt=15000.0),
        dict(key="0-20000-0-00003_B", lat=52.2, lon=13.4, itype="CL31", country="DE",
             name="GAMMA", n_dates=120, success_rate=50.0, methods="Liquid-cloud",
             op_pct_theo=100.0, op_pct_op=100.0, omb_bias=-0.30, icao_alt=5500.0),
    ])


@pytest.mark.parametrize("fn,kwargs", [
    (charts.network_map, {}),
    (charts.omb_bias_map, {}),
    (charts.icao_altitude_map, {}),
])
def test_maps_render_with_data(fn, kwargs):
    fig = fn(_keystats(), **kwargs)
    assert fig.data, f"{fn.__name__} produced no traces"
    geo = [t for t in fig.data if t.type == "scattergeo"]
    assert geo and any(len(t.lat) for t in geo), f"{fn.__name__} has no map points"


def test_ratio_op_map_has_points():
    fig = charts.ratio_map(_keystats(), "op_pct_op", "Median C_L — % of operational constant (L2)",
                           "% of operational", "fig-map-op")
    pts = [t for t in fig.data if t.type == "scattergeo" and len(t.lat)]
    assert pts, "op ratio map has no points despite finite op_pct_op"


def test_maps_use_per_type_symbols():
    """Every map distinguishes instrument type by marker symbol (per-point symbol
    array or per-type legend traces)."""
    fig = charts.icao_altitude_map(_keystats())
    syms = set()
    for t in fig.data:
        sym = getattr(t.marker, "symbol", None)
        if sym is None:
            continue
        if isinstance(sym, (list, tuple, np.ndarray)):
            syms.update(map(str, sym))
        else:
            syms.add(str(sym))
    assert len(syms) >= 2, f"expected multiple instrument symbols, got {syms}"


def test_icao_colorbar_label_short():
    """ICAO colorbar must be the short label (regression: 'remove ICAO from legend')."""
    fig = charts.icao_altitude_map(_keystats())
    titles = []
    for t in fig.data:
        cb = getattr(t.marker, "colorbar", None)
        ttl = getattr(getattr(cb, "title", None), "text", None) if cb is not None else None
        if ttl:
            titles.append(ttl)
    assert titles, "ICAO map has no colorbar title"
    assert all("ICAO" not in ttl for ttl in titles), f"colorbar still says ICAO: {titles}"


# ----------------------------------------------------------------------------
# Sensitivity aggregation = MEDIAN; OmB figure says "v2"; calendar JS present
# ----------------------------------------------------------------------------
def test_combine_sens_uses_median():
    from calibration.sensitivity.network import SensResult, combine_sens_results
    z = np.array([1000.0, 3000.0])
    # one clean month + one bad (high-noise) month: median << mean
    clean = np.array([[1e-3, 1e-3, 1e-3, 1e-3], [1e-2, 1e-2, 1e-2, 1e-2]])
    bad = np.array([[1.0], [1.0]])
    r1 = SensResult(wavelength=910.0, z_ctr=z, dates=np.array(["20250101", "20250102",
                    "20250103", "20250104"], dtype="datetime64[D]"),
                    bmin_night=clean, bmin_day=clean)
    r2 = SensResult(wavelength=910.0, z_ctr=z, dates=np.array(["20250201"], dtype="datetime64[D]"),
                    bmin_night=bad, bmin_day=bad)
    res = combine_sens_results([r1, r2])
    # median of [1e-3 x4, 1.0] = 1e-3 (mean would be ~0.2): mass stays tiny -> detect high
    assert res.sigma_probe[1000] < 1e-2, "aggregation is not the median (bad day dominates)"


def test_omb_figure_labels_v2():
    """OmB plot labels the v2 calibration 'v2', not 'ours' (regression)."""
    src = Path(charts.__file__).resolve().parent.parent / "calibration" / "omb" / "figures.py"
    text = src.read_text(encoding="utf-8")
    assert "obs (v2)" in text and "v2 - CAMS" in text
    assert "obs (ours)" not in text and "ours - CAMS" not in text


def test_calendar_is_three_month_with_dropdown():
    """diag.js renders a 3-month window with arrows + a month dropdown (regression)."""
    js = (Path(charts.__file__).resolve().parent / "static" / "diag.js").read_text(encoding="utf-8")
    assert "CAL_WIN" in js and "cal-select" in js and "cal-prev" in js and "cal-next" in js


def test_calendar_window_is_vertical():
    """The 3 station-page calendars stack VERTICALLY (regression)."""
    import re
    css = (Path(charts.__file__).resolve().parent / "static" / "style.css").read_text(encoding="utf-8")
    m = re.search(r"\.cal-window\s*\{[^}]*\}", css)
    assert m and "column" in m.group(0), f"cal-window not vertical: {m.group(0) if m else 'absent'}"


def test_success_rate_counts_all_nonsuccess():
    """Success rate = valid / ALL days: no-data, no-liquid-cloud and rejection flags all count
    against it (regression for the cloud rate reading ~100%). Identical formula for both methods."""
    from monitoring import index as idx
    # 10 cloud days for one station: 3 valid, 5 no-liquid-cloud(-1), 1 no-data(0), 1 rejection(-20)
    flags = [1, 1, 1, -1, -1, -1, -1, -1, 0, -20]
    cal = pd.DataFrame({
        "key": "K", "method": "cloud",
        "date": [f"2026011{i}" for i in range(10)],
        "flag": flags,
        "cal_value": [1.0 if f == 1 else -1 for f in flags],
    })
    cal["datetime"] = pd.to_datetime(cal["date"], format="%Y%m%d")
    cal["success"] = ((cal["flag"] == 1) & (cal["cal_value"] > 0)).astype(int)
    cal["rel_uncertainty"] = 0.0
    agg = idx._series_aggregates(cal)
    sr = float(agg.iloc[0]["success_rate"])
    assert abs(sr - 30.0) < 0.01, f"expected 30% (3/10), got {sr}"   # NOT 3/(3+1)=75% (old) or ~100%


def test_rayleigh_timeseries_legend_below_and_renamed():
    """Rayleigh C_L time series: legend below the plot + v1.0/v2.0 line names (regression)."""
    dt = pd.to_datetime(["2026-01-01", "2026-01-05", "2026-01-09"])
    g = pd.DataFrame({"datetime": dt, "success": [1, 1, 1],
                      "cal_value": [3e11, 3.1e11, 2.9e11], "uncertainty": [1e10] * 3})
    op = pd.DataFrame({"datetime": dt, "op_coeff": [3e11] * 3})
    ora = pd.DataFrame({"datetime": dt, "value": [2.8e11] * 3})
    fig = charts.series_timeseries(g, pd.DataFrame(), "rayleigh", op_df=op, oldray_df=ora)
    assert fig.layout.legend.y is not None and fig.layout.legend.y < 0, "legend not moved below plot"
    names = [t.name for t in fig.data if t.name]
    assert {"Applied in L2", "v1.0", "v2.0"} <= set(names), names
    # cloud must NOT be mislabeled with the rayleigh version tags
    namesc = [t.name for t in charts.series_timeseries(g, pd.DataFrame(), "cloud", op_df=op).data if t.name]
    assert "v2.0" not in namesc and "v1.0" not in namesc, namesc


# ----------------------------------------------------------------------------
# Image base URL -> absolute image references (EWC object-storage hosting)
# ----------------------------------------------------------------------------
def test_img_base_url_makes_image_refs_absolute(dash, tmp_path, monkeypatch):
    """With ALC_IMG_BASE_URL set, the per-night diagnostic, OmB/sensitivity and flag-example image
    references become ABSOLUTE URLs under that base (so the bulky images can be hosted in an S3 bucket
    while the HTML is served from a small web server); without it they stay site-relative. Regression
    for the EWC port."""
    base = "https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/"
    monkeypatch.setattr(config, "IMG_BASE_URL", base)
    flagex = tmp_path / "flagex"
    flagex.mkdir()
    (flagex / "1__successful.png").write_bytes(_png_bytes())
    out2 = tmp_path / "site_abs"
    render.build_site(dash["db"], out2, opcoeff_csv=dash["opcoeff"], fullcal_dir=dash["fullcal"],
                      flagex_dir=flagex)
    key = f"{_STATIONS[0]['wmo']}_{_STATIONS[0]['ident']}"
    page = (out2 / "stations" / f"{key}.html").read_text(encoding="utf-8")
    assert f"{base}diag/{key}/" in page, "diagnostic image refs are not absolute under the base URL"
    assert f"{base}ombsens/{key}/" in page, "OmB/sensitivity image refs are not absolute under the base URL"
    flags_html = (out2 / "flags.html").read_text(encoding="utf-8")
    assert f"{base}flagex/1__successful.png" in flags_html, "flag-example refs are not absolute"
    # back-compat: the default (no base) build keeps relative paths
    rel_page = (dash["out"] / "stations" / f"{key}.html").read_text(encoding="utf-8")
    assert base not in rel_page, "relative build unexpectedly carries an absolute base URL"
    assert f"../ombsens/{key}/" in rel_page, "relative build lost its relative ombsens path"


def test_diag_js_handles_absolute_urls():
    """diag.js must use an absolute it.rel verbatim (only prefixing '../' for relative paths), so the
    diagnostic viewer can load images straight from the bucket (regression for the EWC port)."""
    js = (Path(charts.__file__).resolve().parent / "static" / "diag.js").read_text(encoding="utf-8")
    assert "^https?" in js, "diag.js does not pass through absolute image URLs"


# ----------------------------------------------------------------------------
def _frames(dash):
    from monitoring import metrics
    return metrics.load_frames(dash["db"])
