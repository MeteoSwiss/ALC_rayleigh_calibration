// Per-station diagnostic viewer: shows the diagnostic image for a chosen day, a month calendar
// (good calibrations green, rejected-but-imaged days grey, days with no image blank). When a viewer
// is hovered/focused, three keyboard axes (mirrored by on-card buttons / header links) navigate:
//   * Left / Right       (prev/next-cal buttons): step between VALID (successful) calibrations.
//   * Ctrl+Left / Right  (the day buttons):       step through ALL imaged days (success + rejected).
//   * Up / Down          (header links):          go to the previous / next STATION page.
// Also syncs with a click on the Plotly time series and on a date in the "all calibrations" table.
// One viewer per method section (.diag). Station pages live at stations/<key>.html, so image paths
// (stored site-root-relative as "diag/<key>/<file>") are prefixed with "../". Prev/next station URLs
// come from #station-nav (data-prev / data-next).
(function () {
  var active = null;  // viewer that the arrow keys control (last interacted with)

  function pad(n) { return (n < 10 ? "0" : "") + n; }

  // Calendar showing a WINDOW of 3 months at a time, with prev/next-month arrows and a
  // month dropdown. Returns a controller: .ensureVisible(date) scrolls the window to a
  // date's month and highlights it (called by the viewer's show()).
  var CAL_WIN = 3;
  function buildCalendar(calEl, items, onPick) {
    var info = {};
    items.forEach(function (it) { info[it.date] = !!it.success; });
    var WD = ["M", "T", "W", "T", "F", "S", "S"];
    // full ordered list of YYYYMM from the first to the last imaged day
    var months = [];
    var y = +items[0].date.slice(0, 4), m = +items[0].date.slice(4, 6) - 1;
    var y1 = +items[items.length - 1].date.slice(0, 4), m1 = +items[items.length - 1].date.slice(4, 6) - 1;
    while (y < y1 || (y === y1 && m <= m1)) { months.push("" + y + pad(m + 1)); m++; if (m > 11) { m = 0; y++; } }

    calEl.innerHTML =
      '<div class="cal-nav">' +
      '<button type="button" class="cal-prev" title="Earlier months">‹</button>' +
      '<select class="cal-select" title="Jump to month"></select>' +
      '<button type="button" class="cal-next" title="Later months">›</button>' +
      '</div><div class="cal-window"></div>';
    var selEl = calEl.querySelector(".cal-select");
    var winEl = calEl.querySelector(".cal-window");
    var prevB = calEl.querySelector(".cal-prev");
    var nextB = calEl.querySelector(".cal-next");
    months.forEach(function (ym, i) {
      var o = document.createElement("option");
      o.value = i; o.textContent = ym.slice(0, 4) + "-" + ym.slice(4, 6);
      selEl.appendChild(o);
    });
    var start = Math.max(0, months.length - CAL_WIN);   // default: the latest 3 months
    var ctl = { selDate: null };

    function monthHtml(ym) {
      var yy = +ym.slice(0, 4), mm = +ym.slice(4, 6) - 1;
      var h = '<div class="cal-month"><div class="cal-title">' + ym.slice(0, 4) + "-" + ym.slice(4, 6) +
        '</div><div class="cal-grid">';
      WD.forEach(function (w) { h += '<div class="cal-wd">' + w + "</div>"; });
      var offset = (new Date(yy, mm, 1).getDay() + 6) % 7;   // Monday-first
      for (var i = 0; i < offset; i++) h += '<div class="cal-day empty"></div>';
      var nd = new Date(yy, mm + 1, 0).getDate();
      for (var dd = 1; dd <= nd; dd++) {
        var ds = "" + yy + pad(mm + 1) + pad(dd);
        var cls = (ds in info) ? ("avail " + (info[ds] ? "good" : "rejected")) : "none";
        h += '<div class="cal-day ' + cls + '" data-date="' + ds + '">' + dd + "</div>";
      }
      return h + "</div></div>";
    }
    function render() {
      var maxStart = Math.max(0, months.length - CAL_WIN);
      start = Math.min(Math.max(0, start), maxStart);
      winEl.innerHTML = months.slice(start, start + CAL_WIN).map(monthHtml).join("");
      selEl.value = start;
      prevB.disabled = (start <= 0);
      nextB.disabled = (start >= maxStart);
      Array.prototype.forEach.call(winEl.querySelectorAll(".cal-day.avail"), function (c) {
        c.addEventListener("click", function () { onPick(c.getAttribute("data-date")); });
      });
      if (ctl.selDate) {
        var cell = winEl.querySelector('.cal-day[data-date="' + ctl.selDate + '"]');
        if (cell) cell.classList.add("sel");
      }
    }
    prevB.addEventListener("click", function () { start -= 1; render(); });
    nextB.addEventListener("click", function () { start += 1; render(); });
    selEl.addEventListener("change", function () { start = +selEl.value; render(); });

    ctl.ensureVisible = function (date) {
      ctl.selDate = date;
      var mi = months.indexOf(date.slice(0, 6));
      if (mi >= 0 && (mi < start || mi >= start + CAL_WIN)) {
        start = mi - (CAL_WIN - 1);   // show the selected month as the most recent of the 3
      }
      render();
    };
    render();
    return ctl;
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
    var cal = null;   // calendar controller (assigned by buildCalendar below)
    var flagBtn = section.querySelector(".diag-flag");
    var idx = validIdx.length ? validIdx[validIdx.length - 1] : items.length - 1;  // default: latest valid

    function curRec() {
      var k = section.getAttribute("data-key") || "";
      var li = k.lastIndexOf("_");
      return {
        key: k,
        wmo: li >= 0 ? k.slice(0, li) : k,
        identifier: li >= 0 ? k.slice(li + 1) : "",
        itype: section.getAttribute("data-itype") || "",
        method: section.getAttribute("data-method") || "",
        date: items[idx] ? items[idx].date : ""
      };
    }
    function syncFlag() {
      if (flagBtn && window.QCFlags) { flagBtn.classList.toggle("flagged", window.QCFlags.has(curRec())); }
    }

    function show(i) {
      if (i < 0 || i >= items.length) return;
      idx = i;
      var it = items[idx], src = /^https?:\/\//.test(it.rel) ? it.rel : "../" + it.rel;
      img.src = src; link.href = src;
      label.textContent = it.date.slice(0, 4) + "-" + it.date.slice(4, 6) + "-" + it.date.slice(6, 8) +
        (it.success ? "  ✓ valid" : "  ✗ rejected") + "  (day " + (idx + 1) + " of " + items.length + ")";
      if (cal) cal.ensureVisible(it.date);   // scroll the 3-month window to this date + highlight
      syncFlag();
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
      openFlag: function () { if (window.QCFlags) window.QCFlags.openDialog(curRec(), syncFlag); },
      quickFlag: function (i) {
        if (window.QCFlags && window.QCFlags.PRESETS && window.QCFlags.PRESETS[i] != null) {
          window.QCFlags.set(curRec(), window.QCFlags.PRESETS[i]); syncFlag();
        }
      },
    };

    cal = buildCalendar(calEl, items, function (d) { viewer.jump(d); active = viewer; });
    section.querySelector(".diag-prev").addEventListener("click", function () { viewer.prevValid(); });
    section.querySelector(".diag-next").addEventListener("click", function () { viewer.nextValid(); });
    var up = section.querySelector(".diag-up"), dn = section.querySelector(".diag-down");
    if (up) up.addEventListener("click", function () { viewer.prevAny(); });
    if (dn) dn.addEventListener("click", function () { viewer.nextAny(); });
    if (flagBtn) flagBtn.addEventListener("click", function () {
      if (window.QCFlags) { window.QCFlags.openDialog(curRec(), syncFlag); }
    });
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

  function navStation(dir) {
    var nav = document.getElementById("station-nav");
    if (!nav) return false;
    var url = dir < 0 ? nav.getAttribute("data-prev") : nav.getAttribute("data-next");
    if (url) { window.location.href = url; return true; }
    return false;  // already at the first / last station
  }

  document.addEventListener("keydown", function (e) {
    if (window.QCFlags && window.QCFlags.isDialogOpen && window.QCFlags.isDialogOpen()) return;  // dialog owns the keyboard
    if (!active || e.target.tagName === "SELECT" || e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA") return;
    var allDays = e.ctrlKey || e.metaKey;   // Ctrl (Cmd on macOS) -> step through ALL imaged days
    if (e.key === "ArrowLeft") { (allDays ? active.prevAny() : active.prevValid()); e.preventDefault(); }
    else if (e.key === "ArrowRight") { (allDays ? active.nextAny() : active.nextValid()); e.preventDefault(); }
    else if (e.key === "ArrowUp") { if (navStation(-1)) e.preventDefault(); }
    else if (e.key === "ArrowDown") { if (navStation(1)) e.preventDefault(); }
    else if (e.key === "0") { active.openFlag(); e.preventDefault(); }            // open flag + comment dialog
    else if (e.key === "1") { active.quickFlag(0); e.preventDefault(); }          // aerosol contamination
    else if (e.key === "2") { active.quickFlag(1); e.preventDefault(); }          // cloud contamination
    else if (e.key === "3") { active.quickFlag(2); e.preventDefault(); }          // low signal (condensation?)
  });

  Array.prototype.forEach.call(document.querySelectorAll(".diag"), buildViewer);
})();
