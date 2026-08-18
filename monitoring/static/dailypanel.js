/* Lazy loader + host-page glue for the interactive daily calibration panel.
 *
 * The panel itself (PANEL_JS, emitted inline by the station template) knows how to DRAW a night.
 * This file is what connects it to the rest of the station page:
 *
 *   - it fetches one day's payload on demand, because six months x two methods is far too much to
 *     embed; the embedded index already carries every day's outcome, so the calendar, the flag
 *     chips and the day arrows are complete before any fetch happens;
 *   - it FOLLOWS the existing per-day diagnostic viewer instead of competing with it. diag.js owns
 *     the keyboard (left/right days, Ctrl for all days, up/down stations, 0/1/2/3 flags) and the
 *     .diaglink dates in the calibration tables; adding a second keydown handler here would give
 *     every arrow press two different meanings. Instead we watch the viewer's own date label and
 *     move the panel to match, so one set of controls drives both.
 */
(function () {
  "use strict";
  if (typeof D === "undefined" || typeof render !== "function") return;   // panel not on this page

  var KEY = (D.station && D.station.key) || "";
  var CACHE = {};
  var pending = 0;

  function url(ds, m) { return "../data/" + KEY + "/" + ds + "_" + m + ".json"; }
  // The public site serves payloads from the image bucket (publish.sh data leg); the base is the
  // same one every diag <img> already uses, recovered from any of them.
  function bucketUrl(ds, m) {
    // Explicit base first (station.html bakes it in whenever the site is built in bucket mode).
    if (window.__payloadBase) return window.__payloadBase + KEY + "/" + ds + "_" + m + ".json";
    // Fallback for older builds: recover the base from any bucket-hosted image on the page.
    var img = document.querySelector('img[src*="/diag/"], img[data-src-all*="/ombsens/"]');
    var src = img ? (img.getAttribute("src") || img.getAttribute("data-src-all") || "") : "";
    var i = src.indexOf("/diag/") >= 0 ? src.indexOf("/diag/") : src.indexOf("/ombsens/");
    if (i < 0) return null;
    return src.slice(0, i) + "/data/" + KEY + "/" + ds + "_" + m + ".json";
  }

  // A page opened by double-click runs on file://, where fetch() of a sibling JSON is blocked as a
  // cross-origin request. The panel then draws nothing at all, which looks exactly like a broken
  // build. Say so instead of failing silently -- this is the most likely way anyone first opens it.
  function warn(why) {
    // A sibling element, NOT #panel: drawPanel rewrites #panel on every render, which both erased
    // this warning an instant after it appeared and (via the old `warned` latch) prevented it from
    // ever coming back.
    var host = document.getElementById("panel");
    if (!host) return;
    var box = document.getElementById("panel-warn");
    if (!box) {
      box = document.createElement("div");
      box.id = "panel-warn";
      host.parentNode.insertBefore(box, host);
    }
    host = box;
    var fileUrl = location.protocol === "file:";
    host.innerHTML =
      '<div style="border:1px solid #f5c2c7;background:#fdecef;border-left:5px solid #b00020;' +
      'border-radius:8px;padding:12px 14px;margin:8px 0">' +
      '<b>The daily panel could not load its data.</b><div style="margin-top:5px;font-size:13px">' +
      (fileUrl
        ? "This page is open from the filesystem (<code>file://</code>), and browsers block a page " +
          "from reading sibling files that way. Serve the folder over HTTP and open it from there — " +
          'e.g. <code>python -m http.server 8000</code> in the dashboard directory, then ' +
          "<code>http://localhost:8000/stations/</code>. Every other panel on this page works " +
          "offline; only this one needs the per-day files."
        : "Could not fetch <code>" + why + "</code>. The per-day payloads may not have been " +
          "generated for this station yet.") +
      "</div></div>";
  }

  function ensureDay(ds, m) {
    var ck = ds + "_" + m;
    if (CACHE[ck] !== undefined) return Promise.resolve(CACHE[ck]);
    pending++;
    return fetch(url(ds, m), { cache: "no-cache" })
      .then(function (r) { return r.ok ? r.json() : Promise.reject(new Error(r.status)); })
      .catch(function () {
        var alt = bucketUrl(ds, m);
        if (!alt) { warn(url(ds, m)); return null; }
        return fetch(alt, { cache: "no-cache" })
          .then(function (r) { return r.ok ? r.json() : Promise.reject(new Error(r.status)); })
          .catch(function () { warn(alt); return null; });
      })
      .then(function (p) {
        CACHE[ck] = p;
        if (p) { D.days[ds] = D.days[ds] || {}; D.days[ds][m] = p; }
        pending--;
        return p;
      });
  }

  function goDay(ds) {
    if (!ds || !/^\d{8}$/.test(ds)) return;
    // Linkable nights: the selected day rides the URL fragment, so back/forward and copied links
    // land on the same night. replaceState -- stepping through days must not spam history.
    try { history.replaceState(null, "", "#d=" + ds); } catch (e) { /* file:// */ }
    if (D.index && !D.index[ds]) {
      var w = document.getElementById("panel-warn");
      if (w) w.textContent = "";
      return;                                  // outside the panel window: nothing to fetch
    }
    curDate = ds;
    Promise.all((D.methods || []).map(function (m) { return ensureDay(ds, m); }))
      .then(function () { if (pending === 0) render(); });
  }

  // The panel asks for a day it does not hold yet (calendar click, its own arrows).
  window.__onMissingDay = function (ds, m) {
    ensureDay(ds, m).then(function (p) { if (p && pending === 0) render(); });
  };

  // ---- follow the production diagnostic viewer -------------------------------------------------
  // diag.js writes the day it is showing into .diag-date. Watching that element means every route
  // into it -- the arrow keys, the calendar, a .diaglink in a table -- moves the panel too, with no
  // second implementation of any of them.
  function watchViewer() {
    var labels = document.querySelectorAll(".diag-date");
    if (!labels.length) return false;
    var obs = new MutationObserver(function () {
      for (var i = 0; i < labels.length; i++) {
        var m = (labels[i].textContent || "").match(/(\d{8})|(\d{4})-(\d{2})-(\d{2})/);
        if (!m) continue;
        var ds = m[1] || (m[2] + m[3] + m[4]);
        if (ds !== curDate) { goDay(ds); return; }
      }
    });
    for (var i = 0; i < labels.length; i++) {
      obs.observe(labels[i], { childList: true, characterData: true, subtree: true });
    }
    return true;
  }

  // A .diaglink click is handled by diag.js for the image viewer; mirror it for the panel so the
  // two stay on the same night even before the observer sees the label change.
  document.addEventListener("click", function (e) {
    var a = e.target.closest && e.target.closest(".diaglink");
    if (a && a.dataset && a.dataset.date) goDay(a.dataset.date);
  });

  function onReady(fn) {
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", fn);
    else fn();
  }

  function qcRec() {
    var st = D.station || {};
    var wmo = st.wmo || (st.key || "").split("_")[0];
    return { key: st.key, wmo: wmo, identifier: st.ident || (st.key || "").split("_").pop(),
             itype: st.type || "", method: (typeof curMethod !== "undefined" && curMethod) || "",
             date: curDate };
  }

  onReady(function () {
    // QC flagging goes through the ONE store (qcflag.js). The visible button always works; the
    // 0/1/2/3 keys are bound here only when no calibration diag viewer owns them (diag.js binds
    // the same keys to the same store when its sections exist -- two writers would double-flag).
    var flagBtn = document.getElementById("dp-flag");
    if (flagBtn && window.QCFlags) {
      flagBtn.hidden = false;
      flagBtn.addEventListener("click", function () { window.QCFlags.openDialog(qcRec()); });
    }
    // The include sits above the method blocks, so at parse time the diag sections do not exist
    // yet -- this is why the panel<->viewer sync never attached before.
    watchViewer();

  // ---- hourly housekeeping on zoom --------------------------------------------------------------
  // The daily chart stays the default; when the operator zooms fig-hk below ~45 days the hourly
  // sidecar is fetched ONCE and the visible traces swap to hourly arrays, restoring the stored
  // daily arrays on zoom-out. Nothing loads for operators who never zoom.
  (function hourlyHk() {
    var gd = document.getElementById("fig-hk");
    if (!gd || !gd.on) return;
    var daily = null, hourly = null, mode = "daily", loading = false;
    function spanDays(e) {
      var r0 = e["xaxis.range[0]"], r1 = e["xaxis.range[1]"];
      if (r0 === undefined || r1 === undefined) {
        var rr = gd.layout && gd.layout.xaxis && gd.layout.xaxis.range;
        if (!rr) return Infinity;
        r0 = rr[0]; r1 = rr[1];
      }
      return (new Date(r1) - new Date(r0)) / 86400000;
    }
    function snapshotDaily() {
      if (daily) return;
      daily = (gd.data || []).map(function (tr) {
        return { x: tr.x, x0: tr.x0, dx: tr.dx, y: tr.y };
      });
    }
    function apply(which) {
      if (which === mode) return;
      if (which === "hourly" && hourly) {
        snapshotDaily();
        (gd.data || []).forEach(function (tr, i) {
          var ys = hourly.series[tr.name];
          if (!ys) return;
          Plotly.restyle(gd, { x: [hourly.t], x0: [null], dx: [null], y: [ys] }, [i]);
        });
        mode = "hourly";
      } else if (which === "daily" && daily) {
        (gd.data || []).forEach(function (tr, i) {
          var d0 = daily[i];
          if (!d0) return;
          Plotly.restyle(gd, { x: [d0.x || null], x0: [d0.x0 || null], dx: [d0.dx || null],
                               y: [d0.y] }, [i]);
        });
        mode = "daily";
      }
    }
    gd.on("plotly_relayout", function (e) {
      var span = spanDays(e);
      if (span < 45) {
        if (hourly) { apply("hourly"); return; }
        if (loading) return;
        loading = true;
        fetch("../data/" + KEY + "/hk_hourly.json", { cache: "no-cache" })
          .then(function (r) { return r.ok ? r.json() : null; })
          .catch(function () { return null; })
          .then(function (j) {
            loading = false;
            if (j && j.t && j.t.length) { hourly = j; apply("hourly"); }
          });
      } else {
        apply("daily");
      }
    });
  })();
    // The availability card promises "click a day to load its diagnostic". With a diag viewer the
    // click reaches the panel through the viewer's date label; without one, nothing happened --
    // so the panel takes the click itself.
    var avail = document.getElementById("fig-avail");
    if (avail && avail.on
        && !document.querySelector('section.diag[data-method="rayleigh"], '
                                   + 'section.diag[data-method="cloud"]')) {
      avail.on("plotly_click", function (ev) {
        var pt = ev && ev.points && ev.points[0];
        if (!pt || pt.x === undefined) return;
        var dt = new Date(pt.x);
        if (isNaN(+dt)) return;
        var ds = dt.toISOString().slice(0, 10).replace(/-/g, "");
        goDay(ds);
      });
    }

    // Clicking a point on any dated station figure opens that night below. The availability card
    // always advertised this; the C_L time series, the calibration-window chart and the Rayleigh
    // overlay carry one marker PER NIGHT, so a click on them means exactly one day and it is the
    // most direct way to inspect an outlier the operator has just spotted.
    function dayFromPoint(pt) {
      if (!pt || pt.x === undefined || pt.x === null) return null;
      var dt = new Date(pt.x);
      if (isNaN(+dt)) return null;
      return dt.toISOString().slice(0, 10).replace(/-/g, "");
    }
    function wireDayClick(gd, tries) {
      if (!gd) return;
      if (typeof gd.on !== "function") {         // Plotly may not have taken the div over yet
        if ((tries || 0) < 40) setTimeout(function () { wireDayClick(gd, (tries || 0) + 1); }, 150);
        return;
      }
      if (gd.__dayClickWired) return;
      gd.__dayClickWired = true;
      gd.on("plotly_click", function (ev) {
        var ds = dayFromPoint(ev && ev.points && ev.points[0]);
        if (!ds) return;
        console.info("[ALC] time-series click ->", ds, "(" + (gd.id || "?") + ")");
        goDay(ds);
        // Keep the per-day diagnostic images in step with the panel, exactly as a click on the
        // availability card does.
        if (typeof window.__diagJump === "function") {
          try { window.__diagJump(ds); } catch (e) { /* viewers are optional */ }
        }
      });
    }
    document.querySelectorAll('[id^="fig-ts-"], [id^="fig-aux-"], #fig-overlay')
      .forEach(function (gd) { wireDayClick(gd, 0); });

    // Dates in the "all calibrations" tables open that night in the panel (the per-calibration
    // diagnostic-PNG viewer they used to load was retired in favour of the panel).
    document.addEventListener("click", function (e) {
      var a = e.target && e.target.closest ? e.target.closest("a.paneldate[data-date]") : null;
      if (!a) return;
      e.preventDefault();
      var ds = a.getAttribute("data-date");
      goDay(ds);
      if (typeof window.__diagJump === "function") {
        try { window.__diagJump(ds); } catch (err) { /* classification card is optional */ }
      }
      var host = document.querySelector(".dp-wrap");
      if (host) host.scrollIntoView({ behavior: "smooth", block: "start" });
    });

    // Plain arrows: diag.js owns them when a CALIBRATION viewer exists; otherwise they would do
    // nothing at all, so the panel steps its own indexed days.
    if (!document.querySelector('section.diag[data-method="rayleigh"], '
                                + 'section.diag[data-method="cloud"]')) {
      document.addEventListener("keydown", function (e) {
        if (e.ctrlKey || e.metaKey) return;
        var t = e.target;
        if (t && (t.tagName === "INPUT" || t.tagName === "SELECT" || t.tagName === "TEXTAREA")) return;
        if (window.QCFlags && window.QCFlags.isDialogOpen && window.QCFlags.isDialogOpen()) return;
        if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
          var dir = e.key === "ArrowRight" ? 1 : -1;
          var i = (D.dates || []).indexOf(curDate);
          var j2 = i + dir;
          if (i < 0 || j2 < 0 || j2 >= D.dates.length) return;
          e.preventDefault();
          goDay(D.dates[j2]);
          return;
        }
        if (!window.QCFlags) return;
        if (e.key === "0") { window.QCFlags.openDialog(qcRec()); e.preventDefault(); }
        else if (e.key === "1" || e.key === "2" || e.key === "3") {
          window.QCFlags.set(qcRec(), window.QCFlags.PRESETS[+e.key - 1]);
          e.preventDefault();
        }
      });
    }
  });

  // Ctrl+Left / Ctrl+Right jump to the previous / next night that PRODUCED A CONSTANT, straight
  // from the embedded index. Capture phase, so it wins over diag.js's own Ctrl binding (step every
  // imaged day): one key must not mean two things depending on handler order.
  document.addEventListener("keydown", function (e) {
    if (!(e.ctrlKey || e.metaKey)) return;
    if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
    var t = e.target;
    if (t && (t.tagName === "INPUT" || t.tagName === "SELECT" || t.tagName === "TEXTAREA")) return;
    var dir = e.key === "ArrowRight" ? 1 : -1;
    var valid = (D.dates || []).filter(function (ds) {
      return (D.methods || []).some(function (m) {
        var s2 = (D.index[ds] || {})[m];
        return s2 && s2.constant !== null && s2.constant !== undefined;
      });
    });
    if (!valid.length) return;
    var i = valid.indexOf(curDate), target;
    if (i < 0) {
      var side = valid.filter(function (d2) { return dir > 0 ? d2 > curDate : d2 < curDate; });
      if (!side.length) return;
      target = dir > 0 ? side[0] : side[side.length - 1];
    } else {
      var k = i + dir;
      if (k < 0 || k >= valid.length) return;
      target = valid[k];
    }
    e.preventDefault();
    e.stopPropagation();
    goDay(target);
  }, true);

  // ---- boot: the most recent night that actually produced a constant ---------------------------
  var withCal = (D.dates || []).filter(function (ds) {
    return (D.methods || []).some(function (m) {
      var s = (D.index[ds] || {})[m];
      return s && s.constant !== null && s.constant !== undefined;
    });
  });
  var start = (withCal.length ? withCal : (D.dates || []))[
    (withCal.length ? withCal : (D.dates || [])).length - 1];
  var frag = (location.hash.match(/[#&]d=(\d{8})/) || [])[1];
  if (frag && D.index && D.index[frag]) start = frag;
  if (start) onReady(function () { goDay(start); });
})();
