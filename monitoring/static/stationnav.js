// Station-page navigation: the country/instrument filters in the navbar scope the previous/next
// station links, and those links show the neighbour's status and calibration constant.
//
// The links used to be baked into each page server-side. They are now computed here from
// data/stations.json, so a filter can change what "next station" means without rebuilding 436
// pages. The mechanism is deliberately minimal: this file WRITES data-prev/data-next back onto the
// existing #station-nav div, so diag.js::navStation (which reads exactly those two attributes) and
// the Up/Down keyboard binding need no change at all.
//
// The filter is kept in sessionStorage under "alc-filter", mirroring how rangesync.js persists
// "alc-period", so it survives the page load that prev/next triggers.
(function () {
  var nav = document.getElementById("station-nav");
  if (!nav || !window.ALCStations) return;
  var here = nav.getAttribute("data-key");
  var base = (document.querySelector(".navsearch") || {}).getAttribute
    ? document.querySelector(".navsearch").getAttribute("data-base") || ""
    : "";
  var fc = document.getElementById("f-country");
  var ft = document.getElementById("f-type");
  var prevA = document.getElementById("nav-prev");
  var nextA = document.getElementById("nav-next");
  var hint = document.getElementById("nav-hint");
  var records = [];

  function saved() {
    try { return JSON.parse(sessionStorage.getItem("alc-filter") || "{}"); } catch (e) { return {}; }
  }
  function save(f) {
    try { sessionStorage.setItem("alc-filter", JSON.stringify(f)); } catch (e) {}
  }

  function label(rec, arrow) {
    if (!rec) return "";
    var dot = rec.q ? '<span class="ndot ndot-' + rec.q + '"></span>' : "";
    var cl = window.ALCStations.clLabel(rec, window.ALCStations.bestMethod(rec));
    return arrow + " " + dot + '<span class="nname">' + (rec.n || rec.k) + "</span>" +
      (cl ? '<span class="ncl">' + cl + "</span>" : "");
  }

  function apply() {
    var c = fc ? fc.value : "";
    var t = ft ? ft.value : "";
    save({ c: c, t: t });
    var subset = records.filter(function (r) {
      return (!c || r.c === c) && (!t || r.t === t);
    });
    var i = subset.findIndex(function (r) { return r.k === here; });
    var scoped = true;
    if (i < 0) {
      // The station being viewed is outside its own filter (e.g. the operator filtered to Germany
      // while looking at a Swiss station). Dead arrows would be worse than useless, so fall back to
      // the unfiltered neighbours and say so.
      subset = records;
      i = subset.findIndex(function (r) { return r.k === here; });
      scoped = false;
    }
    var prev = i > 0 ? subset[i - 1] : null;
    var next = i >= 0 && i < subset.length - 1 ? subset[i + 1] : null;

    nav.setAttribute("data-prev", prev ? base + "stations/" + prev.k + ".html" : "");
    nav.setAttribute("data-next", next ? base + "stations/" + next.k + ".html" : "");
    [[prevA, prev, "&uarr;"], [nextA, next, "&darr;"]].forEach(function (t2) {
      var a = t2[0], rec = t2[1];
      if (!a) return;
      if (!rec) { a.hidden = true; return; }
      a.hidden = false;
      a.href = base + "stations/" + rec.k + ".html";
      a.innerHTML = label(rec, t2[2]);
      a.title = rec.k + (rec.c ? " · " + rec.c : "") + (rec.t ? " · " + rec.t : "");
    });
    if (hint) {
      if (!scoped && (c || t)) {
        hint.hidden = false;
        hint.textContent = "this station is outside the filter — arrows walk the full network";
      } else if (c || t) {
        hint.hidden = false;
        hint.textContent = subset.length + " stations in filter";
      } else {
        hint.hidden = true;
      }
    }
  }

  window.ALCStations.load().then(function (recs) {
    records = recs || [];
    var f = saved();
    if (fc && f.c) fc.value = f.c;
    if (ft && f.t) ft.value = f.t;
    // a stored value that no longer exists in the option list would silently filter everything out
    if (fc && fc.selectedIndex < 0) fc.value = "";
    if (ft && ft.selectedIndex < 0) ft.value = "";
    if (fc) fc.addEventListener("change", apply);
    if (ft) ft.addEventListener("change", apply);
    apply();
  });
})();
