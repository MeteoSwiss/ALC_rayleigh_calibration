// Per-station diagnostic viewer: shows the diagnostic image for a chosen day, a month calendar
// (good calibrations green, rejected-but-imaged days grey, days with no image blank), and two
// navigation axes:
//   * Left / Right  (and the prev/next buttons): step between VALID (successful) calibrations.
//   * Up / Down     (and the up/down buttons):   step through ALL imaged days (success + rejected).
// Also syncs with a click on the Plotly time series and on a date in the "all calibrations" table.
// One viewer per method section (.diag). Station pages live at stations/<key>.html, so image paths
// (stored site-root-relative as "diag/<key>/<file>") are prefixed with "../".
(function () {
  var active = null;  // viewer that the arrow keys control (last interacted with)

  function pad(n) { return (n < 10 ? "0" : "") + n; }

  function buildCalendar(calEl, items, onPick) {
    var info = {};
    items.forEach(function (it) { info[it.date] = !!it.success; });
    var first = items[0].date, last = items[items.length - 1].date;
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
        var cls = (ds in info) ? ("avail " + (info[ds] ? "good" : "rejected")) : "none";
        html += '<div class="cal-day ' + cls + '" data-date="' + ds + '">' + dd + "</div>";
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
    var validIdx = [];
    items.forEach(function (it, i) { if (it.success) validIdx.push(i); });
    var img = section.querySelector(".diag-img");
    var link = section.querySelector(".diag-imglink");
    var label = section.querySelector(".diag-date");
    var calEl = section.querySelector(".diag-cal");
    var idx = validIdx.length ? validIdx[validIdx.length - 1] : items.length - 1;  // default: latest valid

    function show(i) {
      if (i < 0 || i >= items.length) return;
      idx = i;
      var it = items[idx], src = "../" + it.rel;
      img.src = src; link.href = src;
      label.textContent = it.date.slice(0, 4) + "-" + it.date.slice(4, 6) + "-" + it.date.slice(6, 8) +
        (it.success ? "  ✓ valid" : "  ✗ rejected") + "  (day " + (idx + 1) + " of " + items.length + ")";
      Array.prototype.forEach.call(calEl.querySelectorAll(".cal-day.sel"), function (c) { c.classList.remove("sel"); });
      var cell = calEl.querySelector('.cal-day[data-date="' + it.date + '"]');
      if (cell) cell.classList.add("sel");
      active = viewer;
    }
    var viewer = {
      prevValid: function () { var c = null; validIdx.forEach(function (p) { if (p < idx) c = p; }); if (c !== null) show(c); },
      nextValid: function () { for (var k = 0; k < validIdx.length; k++) { if (validIdx[k] > idx) { show(validIdx[k]); return; } } },
      prevAny: function () { show(idx - 1); },
      nextAny: function () { show(idx + 1); },
      jump: function (d) { if (d in pos) show(pos[d]); },
      jumpNearest: function (d) {
        var best = null, bd = Infinity;
        dates.forEach(function (x) { var diff = Math.abs(+x - +d); if (diff < bd) { bd = diff; best = x; } });
        if (best) show(pos[best]);
      },
    };

    buildCalendar(calEl, items, function (d) { viewer.jump(d); active = viewer; });
    section.querySelector(".diag-prev").addEventListener("click", function () { viewer.prevValid(); });
    section.querySelector(".diag-next").addEventListener("click", function () { viewer.nextValid(); });
    var up = section.querySelector(".diag-up"), dn = section.querySelector(".diag-down");
    if (up) up.addEventListener("click", function () { viewer.prevAny(); });
    if (dn) dn.addEventListener("click", function () { viewer.nextAny(); });
    section.addEventListener("mouseenter", function () { active = viewer; });
    section.addEventListener("focus", function () { active = viewer; });
    wireTimeSeries(section.getAttribute("data-ts"), viewer);

    // Clicking a date in this method's "all calibrations" table loads that diagnostic + scrolls up.
    var block = section.closest(".methodblock");
    if (block) {
      Array.prototype.forEach.call(block.querySelectorAll(".diaglink[data-date]"), function (a) {
        a.addEventListener("click", function (e) {
          e.preventDefault();
          viewer.jump(a.getAttribute("data-date"));
          active = viewer;
          section.scrollIntoView({ behavior: "smooth", block: "start" });
        });
      });
    }
    show(idx);
  }

  document.addEventListener("keydown", function (e) {
    if (!active || e.target.tagName === "SELECT" || e.target.tagName === "INPUT") return;
    if (e.key === "ArrowLeft") { active.prevValid(); e.preventDefault(); }
    else if (e.key === "ArrowRight") { active.nextValid(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { active.prevAny(); e.preventDefault(); }
    else if (e.key === "ArrowDown") { active.nextAny(); e.preventDefault(); }
  });

  Array.prototype.forEach.call(document.querySelectorAll(".diag"), buildViewer);
})();
