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
    return Array.prototype.slice.call(Array.prototype.filter.call(document.querySelectorAll(".js-plotly-plot"), function (gd) { return !gd.closest("#daily") && !gd.closest(".dp-rail"); }));
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
      // Compact encoding: a uniform-grid trace ships x0/dx instead of an explicit x array (see
      // charts.monitoring_timeseries). Without this branch the whole trace is skipped and Y stops
      // rescaling on period change -- silently, which is why it must land with that encoding.
      if (!xs && ys && tr.x0 !== undefined && tr.dx) {
        xs = new Array(ys.length);
        for (var j = 0; j < ys.length; j++) xs[j] = +tr.x0 + j * (+tr.dx);
      }
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

  // Console tracing. Every control change prints one grouped line under the [ALC] prefix so a
  // report of "it did not move" can be checked in the operator's own browser: which figures were
  // seen, which were moved, which refused and why. Silent per-figure try/catch is what made the
  // earlier failures invisible.
  var LOG = "[ALC]";

  function applyToPlots(p) {
    var allTime = p.kind === "all" || (!p.start && !p.end);
    var lo = p.start ? p.start + " 00:00:00" : undefined;
    var hi = p.end ? p.end + " 23:59:59" : undefined;
    var moved = [], failed = [];
    plots().forEach(function (gd) {
      var id = gd.id || "(noid)";
      var was = gd.layout && gd.layout.xaxis ? JSON.stringify(gd.layout.xaxis.range) : null;
      try {
        // THE X WINDOW GOES FIRST, ON ITS OWN. Plotly.relayout is atomic: one bad key throws and
        // the whole call is discarded. Bundling the Y rescale with the window meant a figure whose
        // Y could not be touched (the availability heatmap's reversed category axis, and the C_L
        // series, both raise "_inputDomain") silently kept its full time range -- the control said
        // "Last 180 d" while the plot still showed everything. The window is what the operator
        // asked for, so it is applied alone and can no longer be cancelled by the cosmetic part.
        if (allTime) {
          Plotly.relayout(gd, {"xaxis.autorange": true});
        } else {
          Plotly.relayout(gd, {"xaxis.range": [lo, hi], "xaxis.autorange": false});
        }
      } catch (e) {
        failed.push(id + " (x: " + (e && e.message ? e.message : e) + ")");
        return;
      }
      try {
        // Y rescale: best effort, never fatal. Category axes (the availability rows) must keep
        // their own ordering, so they are left alone entirely.
        var yCat = gd.layout && gd.layout.yaxis && gd.layout.yaxis.type === "category";
        var hasY2b = !!(gd.layout && gd.layout.yaxis2);
        var ry = {};
        if (!yCat) {
          if (allTime || hasBars(gd)) {
            ry["yaxis.autorange"] = true;
            if (hasY2b) ry["yaxis2.autorange"] = true;
          } else {
            var yr = visibleYRange(gd, +new Date(lo), +new Date(hi));
            ["y", "y2"].forEach(function (ax) {
              var k = ax === "y" ? "yaxis" : "yaxis2";
              if (k === "yaxis2" && !hasY2b) return;
              if (yr.lo[ax] !== undefined && yr.hi[ax] > yr.lo[ax]) {
                var pad = (yr.hi[ax] - yr.lo[ax]) * 0.08 || Math.abs(yr.hi[ax]) * 0.08 || 1;
                ry[k + ".range"] = [yr.lo[ax] - pad, yr.hi[ax] + pad];
                ry[k + ".autorange"] = false;
              } else {
                ry[k + ".autorange"] = true;
              }
            });
          }
        }
        if (Object.keys(ry).length) Plotly.relayout(gd, ry);
      } catch (e2) { /* the window is already applied; the Y scale is cosmetic */ }
      var now = gd.layout && gd.layout.xaxis ? JSON.stringify(gd.layout.xaxis.range) : null;
      (was === now && !allTime ? failed : moved).push(id);
    });
    console.info(LOG, "period ->", p.key,
                 allTime ? "(all time)" : (p.ymd_start || "?") + ".." + (p.ymd_end || "?"),
                 "| moved:", moved.length ? moved.join(", ") : "none",
                 "| unchanged:", failed.length ? failed.join(", ") : "none");
  }

  function applyToTables(p) {
    var nTables = 0, nHidden = 0;
    document.querySelectorAll("table.stationtable").forEach(function (table) {
      nTables++;
      var tb = table.tBodies[0]; if (!tb) return;
      var changed = false;
      Array.prototype.forEach.call(tb.rows, function (r) {
        var d = r.getAttribute("data-date"); if (!d) return;
        var out = (p.ymd_start && d < p.ymd_start) || (p.ymd_end && d > p.ymd_end);
        if (out) nHidden++;
        var now = out ? "1" : "0";
        if (r.getAttribute("data-phidden") !== now) { r.setAttribute("data-phidden", now); changed = true; }
      });
      if (changed && typeof table._repaginate === "function") table._repaginate();
    });
    console.info(LOG, "tables:", nTables, "| rows hidden by the window:", nHidden);
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

  function applyToHopkin(p) {
    // The pooled C-vs-cloud-base card is not a relayout: it re-pools the nights inside the window,
    // so it is handed the window explicitly (the panel exposes the hook).
    if (typeof window.__hopkinSetWindow === "function") {
      try {
        window.__hopkinSetWindow(p.ymd_start || null, p.ymd_end || null);
        console.info(LOG, "cloud-base card re-pooled for", p.key);
      } catch (e) { console.warn(LOG, "cloud-base card failed:", e); }
    } else {
      console.info(LOG, "cloud-base card: not on this page");
    }
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
    applyToHopkin(p);
  }

  sel.addEventListener("change", function () {
    console.info(LOG, "selector changed by hand ->", sel.value);
    apply(sel.value);
  });

  // Restore a previously chosen window when paging station -> station.
  var saved = null; try { saved = sessionStorage.getItem("alc-period"); } catch (e) {}
  if (saved && byKey[saved] && saved !== sel.value) {
    sel.value = saved;
    // The figures are instantiated by inline scripts further down the page, so applying the
    // remembered window immediately can run against zero (or half) of them -- and since each
    // figure is relayouted inside its own try/catch, the miss is silent and the control ends up
    // reading a window the plots never took. Wait for a stable figure count first.
    var tries = 0, lastN = -1;
    (function settle() {
      var n = document.querySelectorAll(".js-plotly-plot").length;
      if (n > 0 && n === lastN) {
        console.info(LOG, "restoring remembered period", saved, "over", n, "figures");
        apply(saved);
        return;
      }
      lastN = n;
      if (++tries < 40) setTimeout(settle, 150);
      else { console.warn(LOG, "figures never settled; applying", saved, "anyway"); apply(saved); }
    })();
  }
  console.info(LOG, "rangesync ready ·", periods.length, "periods · current", sel.value);
})();
