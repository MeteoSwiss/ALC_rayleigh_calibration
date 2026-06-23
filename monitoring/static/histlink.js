// On the summary page, clicking a bar in the per-type "stations ranked by median C_L" figures opens
// that station's page. Each bar carries the station key as customdata[0]. Plotly attaches its event
// API asynchronously, so we retry until gd.on is available (same approach as diag.js::wireTimeSeries).
(function () {
  function wire(gd, tries) {
    if (typeof gd.on !== "function") {
      if ((tries || 0) < 40) { setTimeout(function () { wire(gd, (tries || 0) + 1); }, 150); }
      return;
    }
    gd.style.cursor = "pointer";
    gd.on("plotly_click", function (data) {
      if (!data || !data.points || !data.points.length) { return; }
      var cd = data.points[0].customdata;
      var key = Array.isArray(cd) ? cd[0] : cd;
      if (key) { window.location.href = "stations/" + key + ".html"; }
    });
  }
  Array.prototype.forEach.call(document.querySelectorAll('[id^="fig-cliqr-"]'), function (gd) { wire(gd, 0); });
})();
