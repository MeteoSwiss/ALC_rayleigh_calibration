// Nav-bar filters: filter every filterable table (series + watchlist) AND the network map by
// country / instrument type, in sync. Filtered-out rows get data-fhidden="1"; paginated tables
// are re-paged over the matching rows (paginate.js). The map's scattergeo trace carries per-point
// customdata = [country, type, base_size, key, name]; we restyle marker.size (0 = hidden) and
// navigate to the station page on click. No dependencies beyond the page-global Plotly.
(function () {
  function pad(n) { return (n < 10 ? "0" : "") + n; }

  // All scattergeo maps share the customdata layout [country, type, size, key, name] on data[0].
  // (Maps also carry data[1..5] = dummy symbol-legend traces; we only ever touch trace [0].)
  var MAP_IDS = ["fig-map", "fig-map-theo", "fig-map-op", "fig-map-omb", "fig-map-icao"];

  // Instrument palette (mirror of monitoring/config.TYPE_ORDER / TYPE_COLORS) for the client-side
  // rebuild of the "instruments over time" stacked chart on filter change.
  var TYPE_ORDER = ["CHM15k", "CL31", "CL51", "CL61", "Mini-MPL"];
  var TYPE_COLORS = { "CHM15k": "#1f77b4", "CL31": "#ff7f0e", "CL51": "#2ca02c",
                      "CL61": "#d62728", "Mini-MPL": "#9467bd" };

  function readJson(id) {
    var el = document.getElementById(id);
    if (!el) return null;
    try { return JSON.parse(el.textContent); } catch (e) { return null; }
  }
  var seriesIndex = readJson("series-index");   // [{key,itype,country,method,n_dates,n_success}]
  var activity = readJson("instr-activity");     // [{key,itype,country,months:[YYYYMM,...]}]

  // Recompute the 4 headline KPIs from the filtered per-series rows (period is already baked per page).
  function recomputeKpis(c, t) {
    if (!seriesIndex) return;
    var nSeries = 0, cals = 0, succ = 0, keys = {};
    seriesIndex.forEach(function (r) {
      if ((c && r.country !== c) || (t && r.itype !== t)) return;
      nSeries++; cals += (r.n_dates || 0); succ += (r.n_success || 0); keys[r.key] = 1;
    });
    function set(id, v) { var e = document.getElementById(id); if (e) e.textContent = v; }
    set("kpi-instruments", Object.keys(keys).length);
    set("kpi-series", nSeries);
    set("kpi-cals", cals);
    set("kpi-success", (cals ? (100 * succ / cals).toFixed(1) : "0.0") + "%");
  }

  // Rebuild the stacked "active instruments over time" chart for the current filter. Country filtering
  // changes per-month counts (not a simple restyle), so recompute from the embedded activity list.
  function recomputeInstrumentChart(c, t) {
    var gd = document.getElementById("fig-instr");
    if (!gd || !window.Plotly || !activity) return;
    var monthsSet = {}, typesPresent = {};
    activity.forEach(function (r) {
      typesPresent[r.itype] = 1;
      (r.months || []).forEach(function (m) { monthsSet[m] = 1; });
    });
    var months = Object.keys(monthsSet).sort();
    if (!months.length) return;
    var mIdx = {}; months.forEach(function (m, i) { mIdx[m] = i; });
    var x = months.map(function (m) { return m.slice(0, 4) + "-" + m.slice(4, 6) + "-01"; });
    var order = TYPE_ORDER.filter(function (tp) { return typesPresent[tp]; })
      .concat(Object.keys(typesPresent).filter(function (tp) { return tp && TYPE_ORDER.indexOf(tp) < 0; }).sort());
    if (t) order = order.filter(function (tp) { return tp === t; });
    var counts = {};
    order.forEach(function (tp) { counts[tp] = months.map(function () { return 0; }); });
    activity.forEach(function (r) {
      if ((c && r.country !== c) || (t && r.itype !== t) || !counts[r.itype]) return;
      var seen = {};
      (r.months || []).forEach(function (m) {
        if (m in mIdx && !seen[m]) { seen[m] = 1; counts[r.itype][mIdx[m]]++; }
      });
    });
    var data = order.filter(function (tp) { return counts[tp].some(function (v) { return v > 0; }); })
      .map(function (tp) { return { type: "bar", x: x, y: counts[tp], name: tp,
                                    marker: { color: TYPE_COLORS[tp] || "#888" } }; });
    window.Plotly.react(gd, data, gd.layout || {});
  }

  // Click a map marker -> open that station's page (retry until Plotly has initialised).
  function wireMapClick(id, tries) {
    var gd = document.getElementById(id);
    if (!gd) return;
    if (typeof gd.on !== "function") {
      if ((tries || 0) < 40) setTimeout(function () { wireMapClick(id, (tries || 0) + 1); }, 150);
      return;
    }
    gd.on("plotly_click", function (data) {
      if (!data || !data.points || !data.points.length) return;
      var cd = data.points[0].customdata;  // [country, type, size, key, name]
      if (cd && cd[3]) window.location.href = "stations/" + cd[3] + ".html";
    });
  }
  MAP_IDS.forEach(function (id) { wireMapClick(id, 0); });

  var fc = document.getElementById("f-country");
  var ft = document.getElementById("f-type");
  if (!fc && !ft) return;
  // Station pages carry the same two selects, but there they scope the previous/next navigation
  // (stationnav.js) and there is nothing filterable on the page. Without this guard every change
  // would run five no-op passes over tables/maps/histograms that do not exist here.
  if (!document.querySelector("table.filterable")) return;
  var tables = Array.prototype.slice.call(document.querySelectorAll("table.filterable"));
  var count = document.getElementById("filter-count");
  // Per-instrument-type ranked-histogram sections ("Median C_L per station, by type"): when a type
  // is selected, show only its histogram (and the section header); otherwise show all of them.
  var cliqrSections = Array.prototype.slice.call(document.querySelectorAll("section.cliqr"));
  var cliqrHeader = document.querySelector("h2.sec-h");

  function val(sel) { return sel ? sel.value : ""; }

  function filterTables(c, t) {
    var seriesShown = 0;
    tables.forEach(function (table) {
      var tbody = table.tBodies[0];
      if (!tbody) return;
      var shown = 0;
      Array.prototype.forEach.call(tbody.rows, function (r) {
        var ok = (!c || r.getAttribute("data-country") === c) &&
                 (!t || r.getAttribute("data-type") === t);
        if (ok) { r.removeAttribute("data-fhidden"); shown++; }
        else { r.setAttribute("data-fhidden", "1"); r.style.display = "none"; }
      });
      if (typeof table._repaginate === "function") table._repaginate();
      else Array.prototype.forEach.call(tbody.rows, function (r) {
        r.style.display = r.getAttribute("data-fhidden") === "1" ? "none" : "";
      });
      if (table.id === "stations") seriesShown = shown;
    });
    if (count) count.textContent = seriesShown + " series match";
  }

  function filterMap(c, t) {
    if (!window.Plotly) return;
    MAP_IDS.forEach(function (id) {
      var gd = document.getElementById(id);
      if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return;
      var sizes = gd.data[0].customdata.map(function (p) {
        var ok = (!c || p[0] === c) && (!t || p[1] === t);
        return ok ? p[2] : 0;  // 0 hides the marker
      });
      window.Plotly.restyle(gd, { "marker.size": [sizes] }, [0]);
    });
  }

  // Filter one ranked-histogram figure's BARS to a country (re-ranking the survivors 0..M-1) and
  // update its station-count title; returns the number of bars shown. The full bar arrays are cached
  // on first call so every filter is taken from the complete set. The theoretical / network-median
  // reference lines are layout shapes and intentionally stay at their full-network values.
  function filterBars(gd, c) {
    if (!gd.data || !gd.data[0] || !gd.data[0].customdata) return 1;   // not ready -> don't hide it
    if (!gd._fullbars) {
      // Read EVERYTHING from customdata = [key, n, q1, q3, country, median]. Plotly stores y and
      // error_y.array as base64-encoded typed arrays (Array.from gives []), so reconstruct the bar
      // height from the median (cd[5]) and the IQR error bars from q1/q3 (cd[2]/cd[3]) -- otherwise
      // the bars get undefined heights and vanish even though the count/axis update.
      gd._fullbars = { cd: Array.from(gd.data[0].customdata || []), itype: gd.data[0].name || "" };
    }
    var f = gd._fullbars, idx = [];
    for (var i = 0; i < f.cd.length; i++) {
      if (!c || (f.cd[i] && f.cd[i][4] === c)) idx.push(i);   // [4] = country
    }
    window.Plotly.restyle(gd, {
      x: [idx.map(function (_, j) { return j; })],
      y: [idx.map(function (i) { return f.cd[i][5]; })],                                  // median
      "error_y.array": [idx.map(function (i) { return Math.max(0, f.cd[i][3] - f.cd[i][5]); })],     // q3 - median
      "error_y.arrayminus": [idx.map(function (i) { return Math.max(0, f.cd[i][5] - f.cd[i][2]); })], // median - q1
      customdata: [idx.map(function (i) { return f.cd[i]; })]
    }, [0]);
    window.Plotly.relayout(gd, { "title.text": f.itype +
      " — stations ranked by median C_L (n=" + idx.length + " stations, error bar = IQR)" });
    return idx.length;
  }

  function filterHistograms(c, t) {
    var anyShown = 0;
    cliqrSections.forEach(function (s) {
      var gd = s.querySelector(".plotly-graph-div");
      var n = (gd && window.Plotly) ? filterBars(gd, c) : 1;   // subset bars by country (all sections)
      var show = (!t || s.getAttribute("data-type") === t) && n > 0;   // hide empty / off-type ones
      s.style.display = show ? "" : "none";
      if (show) { anyShown++; try { window.Plotly.Plots.resize(gd); } catch (e) {} }
    });
    if (cliqrHeader) cliqrHeader.style.display = anyShown ? "" : "none";
  }

  function apply() {
    var c = val(fc), t = val(ft);
    // Same sessionStorage key stationnav.js uses, so one filter follows the operator in BOTH
    // directions between the summary and the station pages instead of resetting at each hop.
    try { sessionStorage.setItem("alc-filter", JSON.stringify({ c: c, t: t })); } catch (e) {}
    filterTables(c, t);
    filterMap(c, t);
    filterHistograms(c, t);
    recomputeKpis(c, t);
    recomputeInstrumentChart(c, t);
  }
  if (fc) fc.addEventListener("change", apply);
  if (ft) ft.addEventListener("change", apply);
  // Restore a filter carried over from a station page. Deferred to the next tick so the Plotly
  // figure divs the map/histogram passes need are already in the DOM.
  setTimeout(function () {
    var f;
    try { f = JSON.parse(sessionStorage.getItem("alc-filter") || "{}"); } catch (e) { f = {}; }
    if (!f || (!f.c && !f.t)) return;
    if (fc && f.c) fc.value = f.c;
    if (ft && f.t) ft.value = f.t;
    if (fc && fc.selectedIndex < 0) fc.value = "";     // a value this build no longer offers
    if (ft && ft.selectedIndex < 0) ft.value = "";
    if (val(fc) || val(ft)) apply();
  }, 0);
})();
