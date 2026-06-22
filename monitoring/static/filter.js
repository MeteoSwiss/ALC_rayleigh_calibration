// Nav-bar filters: filter the series table AND the network map by country / instrument type.
// The map's scattergeo trace carries per-point customdata = [country, type, base_size]; we
// restyle marker.size (0 = hidden) so the map stays in sync with the table. No dependencies
// beyond the page-global Plotly.
(function () {
  // Click a map marker -> open that station's page (customdata[3] = key). Retry until the
  // Plotly map has initialised (gd.on appears only after Plotly.newPlot, which may run after
  // this script). Wired regardless of whether the filter controls are present.
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
  var rows = Array.prototype.slice.call(document.querySelectorAll("#stations tbody tr"));
  var count = document.getElementById("filter-count");

  function val(sel) { return sel ? sel.value : ""; }

  function filterTable(c, t) {
    var shown = 0;
    rows.forEach(function (r) {
      var ok = (!c || r.getAttribute("data-country") === c) &&
               (!t || r.getAttribute("data-type") === t);
      r.style.display = ok ? "" : "none";
      if (ok) shown++;
    });
    if (count) count.textContent = shown + " of " + rows.length + " series";
  }

  function filterMap(c, t) {
    var gd = document.getElementById("fig-map");
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata || !window.Plotly) return;
    var cd = gd.data[0].customdata;            // [country, type, baseSize] per point
    var sizes = cd.map(function (p) {
      var ok = (!c || p[0] === c) && (!t || p[1] === t);
      return ok ? p[2] : 0;                    // 0 hides the marker
    });
    window.Plotly.restyle(gd, { "marker.size": [sizes] }, [0]);
  }

  function apply() {
    var c = val(fc), t = val(ft);
    filterTable(c, t);
    filterMap(c, t);
  }
  if (fc) fc.addEventListener("change", apply);
  if (ft) ft.addEventListener("change", apply);
})();
