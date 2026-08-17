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
      .catch(function () { warn(url(ds, m)); return null; })
      .then(function (p) {
        CACHE[ck] = p;
        if (p) { D.days[ds] = D.days[ds] || {}; D.days[ds][m] = p; }
        pending--;
        return p;
      });
  }

  function goDay(ds) {
    if (!ds || !/^\d{8}$/.test(ds)) return;
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
  if (start) onReady(function () { goDay(start); });
})();
