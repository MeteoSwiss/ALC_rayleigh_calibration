// Nav-bar filters: filter every filterable table (series + watchlist) AND the network map by
// country / instrument type, in sync. Filtered-out rows get data-fhidden="1"; paginated tables
// are re-paged over the matching rows (paginate.js). The map's scattergeo trace carries per-point
// customdata = [country, type, base_size, key, name]; we restyle marker.size (0 = hidden) and
// navigate to the station page on click. No dependencies beyond the page-global Plotly.
(function () {
  function pad(n) { return (n < 10 ? "0" : "") + n; }

  // Click a map marker -> open that station's page (retry until Plotly has initialised).
  function wireMapClick(tries) {
    var gd = document.getElementById("fig-map");
    if (!gd) return;
    if (typeof gd.on !== "function") {
      if ((tries || 0) < 40) setTimeout(function () { wireMapClick((tries || 0) + 1); }, 150);
      return;
    }
    gd.on("plotly_click", function (data) {
      if (!data || !data.points || !data.points.length) return;
      var cd = data.points[0].customdata;  // [country, type, size, key, name]
      if (cd && cd[3]) window.location.href = "stations/" + cd[3] + ".html";
    });
  }
  wireMapClick(0);

  var fc = document.getElementById("f-country");
  var ft = document.getElementById("f-type");
  if (!fc && !ft) return;
  var tables = Array.prototype.slice.call(document.querySelectorAll("table.filterable"));
  var count = document.getElementById("filter-count");

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
    var gd = document.getElementById("fig-map");
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata || !window.Plotly) return;
    var cd = gd.data[0].customdata;
    var sizes = cd.map(function (p) {
      var ok = (!c || p[0] === c) && (!t || p[1] === t);
      return ok ? p[2] : 0;  // 0 hides the marker
    });
    window.Plotly.restyle(gd, { "marker.size": [sizes] }, [0]);
  }

  function apply() {
    var c = val(fc), t = val(ft);
    filterTables(c, t);
    filterMap(c, t);
  }
  if (fc) fc.addEventListener("change", apply);
  if (ft) ft.addEventListener("change", apply);
})();
