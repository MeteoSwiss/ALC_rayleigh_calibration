// Period selector — one header control drives a whole page's time window client-side.
//
// Station pages embed every figure's full history + the full calibration table + the diagnostic
// calendar, so picking a window is pure relayout/hide: no rebuild, no server round-trip.
//   #period-sel  (station pages) -> Plotly.relayout every figure to the window (+ rescale Y to the
//                                   visible points for line/marker plots), hide out-of-window table
//                                   rows, swap the OmB/sensitivity images, and tell the calendar.
//   #period-nav  (summary page)  -> navigate to the chosen prebuilt per-period summary page.
//
// The period list is embedded as JSON in #period-index (produced by monitoring/periods.py).
(function () {
  var idxEl = document.getElementById("period-index");
  var periods = [];
  if (idxEl) { try { periods = JSON.parse(idxEl.textContent); } catch (e) { periods = []; } }
  var byKey = {};
  periods.forEach(function (p) { byKey[p.key] = p; });

  // ---- summary page: the control just navigates between prebuilt pages ----
  var nav = document.getElementById("period-nav");
  if (nav) nav.addEventListener("change", function () { if (nav.value) location.href = nav.value; });

  // ---- station page: client-side relayout ----
  var sel = document.getElementById("period-sel");
  if (!sel) return;

  function plots() {
    return Array.prototype.slice.call(document.querySelectorAll(".js-plotly-plot"));
  }

  function hasBars(gd) {
    return (gd.data || []).some(function (tr) { return tr.type === "bar"; });
  }

  function visibleYRange(gd, x0, x1) {
    // min/max over y of points whose x falls in [x0,x1], grouped by the trace's y-axis ('y'/'y2').
    var lo = {}, hi = {};
    (gd.data || []).forEach(function (tr) {
      if (tr.visible === "legendonly" || tr.visible === false) return;
      var xs = tr.x, ys = tr.y;
      if (!xs || !ys) return;
      var ax = tr.yaxis || "y";
      for (var i = 0; i < xs.length; i++) {
        var xv = +new Date(xs[i]);
        if (isNaN(xv) || xv < x0 || xv > x1) continue;
        var yv = +ys[i];
        if (ys[i] == null || isNaN(yv) || !isFinite(yv)) continue;
        if (lo[ax] === undefined || yv < lo[ax]) lo[ax] = yv;
        if (hi[ax] === undefined || yv > hi[ax]) hi[ax] = yv;
      }
    });
    return { lo: lo, hi: hi };
  }

  function applyToPlots(p) {
    var allTime = p.kind === "all" || (!p.start && !p.end);
    var lo = p.start ? p.start + " 00:00:00" : undefined;
    var hi = p.end ? p.end + " 23:59:59" : undefined;
    plots().forEach(function (gd) {
      try {
        var rl = {};
        if (allTime) {
          rl["xaxis.autorange"] = true; rl["yaxis.autorange"] = true; rl["yaxis2.autorange"] = true;
          Plotly.relayout(gd, rl); return;
        }
        rl["xaxis.range"] = [lo, hi];
        rl["xaxis.autorange"] = false;
        // Rescale Y to the visible points for line/marker plots; bar charts keep their from-zero
        // autorange so stacked monthly bars are never clipped.
        if (hasBars(gd)) {
          rl["yaxis.autorange"] = true; rl["yaxis2.autorange"] = true;
        } else {
          var yr = visibleYRange(gd, +new Date(lo), +new Date(hi));
          ["y", "y2"].forEach(function (ax) {
            var k = ax === "y" ? "yaxis" : "yaxis2";
            if (yr.lo[ax] !== undefined && yr.hi[ax] > yr.lo[ax]) {
              var pad = (yr.hi[ax] - yr.lo[ax]) * 0.08 || Math.abs(yr.hi[ax]) * 0.08 || 1;
              rl[k + ".range"] = [yr.lo[ax] - pad, yr.hi[ax] + pad];
              rl[k + ".autorange"] = false;
            } else {
              rl[k + ".autorange"] = true;   // nothing visible on this axis -> leave it
            }
          });
        }
        Plotly.relayout(gd, rl);
      } catch (e) { /* leave this figure as-is */ }
    });
  }

  function applyToTables(p) {
    document.querySelectorAll("table.stationtable").forEach(function (table) {
      var tb = table.tBodies[0]; if (!tb) return;
      var changed = false;
      Array.prototype.forEach.call(tb.rows, function (r) {
        var d = r.getAttribute("data-date"); if (!d) return;
        var out = (p.ymd_start && d < p.ymd_start) || (p.ymd_end && d > p.ymd_end);
        var now = out ? "1" : "0";
        if (r.getAttribute("data-phidden") !== now) { r.setAttribute("data-phidden", now); changed = true; }
      });
      if (changed && typeof table._repaginate === "function") table._repaginate();
    });
  }

  function applyToImages(p) {
    // <img data-period-img> carries per-period URLs as data-src-<periodkey>; swap to the chosen one
    // (falls back to the all-time image when a given window has no dedicated render).
    document.querySelectorAll("img[data-period-img]").forEach(function (img) {
      var url = img.getAttribute("data-src-" + p.key) || img.getAttribute("data-src-all");
      if (!url) return;
      img.src = url;
      if (img.parentElement && img.parentElement.tagName === "A") img.parentElement.href = url;
    });
  }

  function applyToCalendar(p) {
    // diag.js may expose window.__diagSetWindow(ymd_start, ymd_end) to constrain its calendar.
    if (typeof window.__diagSetWindow === "function") {
      try { window.__diagSetWindow(p.ymd_start || null, p.ymd_end || null); } catch (e) {}
    }
  }

  function apply(key) {
    var p = byKey[key]; if (!p) return;
    try { sessionStorage.setItem("alc-period", key); } catch (e) {}
    applyToPlots(p); applyToTables(p); applyToImages(p); applyToCalendar(p);
  }

  sel.addEventListener("change", function () { apply(sel.value); });

  // Restore a previously chosen window when paging station -> station.
  var saved = null; try { saved = sessionStorage.getItem("alc-period"); } catch (e) {}
  if (saved && byKey[saved] && saved !== sel.value) { sel.value = saved; apply(saved); }
})();
