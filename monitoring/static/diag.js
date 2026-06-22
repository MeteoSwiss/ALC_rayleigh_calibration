// Per-station diagnostic viewer: shows the latest successful calibration image, a month
// calendar (non-calibration days greyed), prev/next + arrow-key navigation, and syncs with a
// click on the Plotly time series. One viewer per method section (.diag). No dependencies
// beyond the page-global Plotly. Station pages live at stations/<key>.html, so image paths
// (stored site-root-relative as "diag/<key>/<file>") are prefixed with "../".
(function () {
  var active = null;  // viewer that arrow keys control (last interacted with)

  function pad(n) { return (n < 10 ? "0" : "") + n; }

  function buildCalendar(calEl, dates, onPick) {
    var avail = {};
    dates.forEach(function (d) { avail[d] = 1; });
    var first = dates[0], last = dates[dates.length - 1];
    var y = +first.slice(0, 4), m = +first.slice(4, 6) - 1;
    var y1 = +last.slice(0, 4), m1 = +last.slice(4, 6) - 1;
    var WD = ["M", "T", "W", "T", "F", "S", "S"];
    var html = "";
    while (y < y1 || (y === y1 && m <= m1)) {
      html += '<div class="cal-month"><div class="cal-title">' + y + "-" + pad(m + 1) + "</div><div class=\"cal-grid\">";
      WD.forEach(function (w) { html += '<div class="cal-wd">' + w + "</div>"; });
      var offset = (new Date(y, m, 1).getDay() + 6) % 7;   // Monday-first
      for (var i = 0; i < offset; i++) html += '<div class="cal-day empty"></div>';
      var nd = new Date(y, m + 1, 0).getDate();
      for (var dd = 1; dd <= nd; dd++) {
        var ds = "" + y + pad(m + 1) + pad(dd);
        html += '<div class="cal-day ' + (avail[ds] ? "avail" : "none") + '" data-date="' + ds + '">' + dd + "</div>";
      }
      html += "</div></div>";
      m++; if (m > 11) { m = 0; y++; }
    }
    calEl.innerHTML = html;
    Array.prototype.forEach.call(calEl.querySelectorAll(".cal-day.avail"), function (c) {
      c.addEventListener("click", function () { onPick(c.getAttribute("data-date")); });
    });
  }

  function wireTimeSeries(tsId, viewer, tries) {
    var gd = document.getElementById(tsId);
    if (!gd) return;
    if (typeof gd.on !== "function") {
      if ((tries || 0) < 40) setTimeout(function () { wireTimeSeries(tsId, viewer, (tries || 0) + 1); }, 150);
      return;
    }
    gd.on("plotly_click", function (data) {
      if (!data || !data.points || !data.points.length) return;
      var x = data.points[0].x, ds;
      if (typeof x === "string" && /^\d{4}-\d{2}-\d{2}/.test(x)) ds = x.slice(0, 10).replace(/-/g, "");
      else { var d = new Date(x); ds = "" + d.getUTCFullYear() + pad(d.getUTCMonth() + 1) + pad(d.getUTCDate()); }
      viewer.jumpNearest(ds);
      active = viewer;
    });
  }

  function buildViewer(section) {
    var dataEl = section.querySelector(".diag-data");
    var items;
    try { items = JSON.parse(dataEl.textContent); } catch (e) { return; }
    if (!items || !items.length) return;
    items.sort(function (a, b) { return a.date < b.date ? -1 : 1; });
    var dates = items.map(function (it) { return it.date; });
    var pos = {}; dates.forEach(function (d, i) { pos[d] = i; });
    var img = section.querySelector(".diag-img");
    var link = section.querySelector(".diag-imglink");
    var label = section.querySelector(".diag-date");
    var calEl = section.querySelector(".diag-cal");
    var idx = items.length - 1;  // default: latest successful

    function show(i) {
      if (i < 0 || i >= items.length) return;
      idx = i;
      var it = items[idx], src = "../" + it.rel;
      img.src = src; link.href = src;
      label.textContent = it.date.slice(0, 4) + "-" + it.date.slice(4, 6) + "-" + it.date.slice(6, 8) +
        "  (" + (idx + 1) + " of " + items.length + ")";
      Array.prototype.forEach.call(calEl.querySelectorAll(".cal-day.sel"), function (c) { c.classList.remove("sel"); });
      var cell = calEl.querySelector('.cal-day[data-date="' + it.date + '"]');
      if (cell) cell.classList.add("sel");
      active = viewer;
    }
    var viewer = {
      prev: function () { show(idx - 1); },
      next: function () { show(idx + 1); },
      jump: function (d) { if (d in pos) show(pos[d]); },
      jumpNearest: function (d) {
        var best = null, bd = Infinity;
        dates.forEach(function (x) { var diff = Math.abs(+x - +d); if (diff < bd) { bd = diff; best = x; } });
        if (best) show(pos[best]);
      },
    };

    buildCalendar(calEl, dates, function (d) { viewer.jump(d); active = viewer; });
    section.querySelector(".diag-prev").addEventListener("click", function () { viewer.prev(); });
    section.querySelector(".diag-next").addEventListener("click", function () { viewer.next(); });
    section.addEventListener("mouseenter", function () { active = viewer; });
    section.addEventListener("focus", function () { active = viewer; });
    wireTimeSeries(section.getAttribute("data-ts"), viewer);
    show(idx);
  }

  document.addEventListener("keydown", function (e) {
    if (!active || e.target.tagName === "SELECT" || e.target.tagName === "INPUT") return;
    if (e.key === "ArrowLeft") { active.prev(); e.preventDefault(); }
    else if (e.key === "ArrowRight") { active.next(); e.preventDefault(); }
  });

  Array.prototype.forEach.call(document.querySelectorAll(".diag"), buildViewer);
})();
