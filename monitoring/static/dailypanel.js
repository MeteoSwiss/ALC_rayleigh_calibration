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

  function ensureDay(ds, m) {
    var ck = ds + "_" + m;
    if (CACHE[ck] !== undefined) return Promise.resolve(CACHE[ck]);
    pending++;
    return fetch(url(ds, m), { cache: "no-cache" })
      .then(function (r) { return r.ok ? r.json() : null; })
      .catch(function () { return null; })
      .then(function (p) {
        CACHE[ck] = p;
        if (p) { D.days[ds] = D.days[ds] || {}; D.days[ds][m] = p; }
        pending--;
        return p;
      });
  }

  function goDay(ds) {
    if (!ds || !/^\d{8}$/.test(ds)) return;
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

  watchViewer();

  // ---- boot: the most recent night that actually produced a constant ---------------------------
  var withCal = (D.dates || []).filter(function (ds) {
    return (D.methods || []).some(function (m) {
      var s = (D.index[ds] || {})[m];
      return s && s.constant !== null && s.constant !== undefined;
    });
  });
  var start = (withCal.length ? withCal : (D.dates || []))[
    (withCal.length ? withCal : (D.dates || [])).length - 1];
  if (start) goDay(start);
})();
