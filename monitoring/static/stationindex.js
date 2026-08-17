// Shared station index: one fetch of data/stations.json per browser session, reused by the nav-bar
// search (search.js) and by the filtered previous/next navigation (stationnav.js).
//
// Why a fetch and not an inline blob: the index is the same ~100 kB for every page, so inlining it
// duplicates it into all 436 station pages. Fetched once, it is browser-cached across the whole
// site; its URL carries a content hash so a changed index busts that cache and an unchanged one
// does not.
//
// file:// fallback: an operator opening a local build by double-click cannot fetch (the browser
// blocks it for file:// origins), so we fall back to the inline #search-index blob, which carries
// the same records minus the status/constant fields. Callers must therefore treat q/m as optional.
(function () {
  var URL_EL = document.querySelector(".navsearch");
  var url = URL_EL && URL_EL.getAttribute("data-index");
  var promise = null;

  function fromInline() {
    var el = document.getElementById("search-index");
    if (!el) return [];
    try {
      // the inline blob uses long field names; normalise it to the compact schema
      return JSON.parse(el.textContent).map(function (r) {
        return { k: r.key, n: r.name, w: r.wigos, c: r.country, t: r.type, q: "", m: {} };
      });
    } catch (e) { return []; }
  }

  function load() {
    if (promise) return promise;
    if (!url || !window.fetch || location.protocol === "file:") {
      promise = Promise.resolve(fromInline());
      return promise;
    }
    promise = fetch(url)
      .then(function (r) { if (!r.ok) throw new Error("HTTP " + r.status); return r.json(); })
      .catch(function () { return fromInline(); });
    return promise;
  }

  // Percent of the instrument type's nominal constant is what makes C_L readable at a glance --
  // the raw value spans 11 orders of magnitude between a CL61 (~1) and a CHM15k (~3e11).
  function clLabel(rec, method) {
    var m = rec && rec.m && rec.m[method];
    if (!m || m.pct == null) return "";
    return Math.round(m.pct) + "% nom.";
  }

  function bestMethod(rec) {
    if (!rec || !rec.m) return null;
    if (rec.m.rayleigh) return "rayleigh";
    var ks = Object.keys(rec.m);
    return ks.length ? ks[0] : null;
  }

  window.ALCStations = { load: load, clLabel: clLabel, bestMethod: bestMethod };
})();
