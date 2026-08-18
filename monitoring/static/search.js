// Nav-bar station search. Matches a free-text query against each station's name, WIGOS id and
// key (substring, case-insensitive), so "Payerne", "Pay", "0-20000-0-06610" and "06610" all find
// the 0-20000-0-06610_A/_B/_C pages. Arrow keys move the selection, Enter opens it, click opens it.
// Station index is embedded as JSON (#search-index); navigation is relative to data-base.
(function () {
  var box = document.getElementById("station-search");
  var panel = document.getElementById("search-results");
  var wrap = document.querySelector(".navsearch");
  if (!box || !panel || !wrap || !window.ALCStations) return;
  var base = wrap.getAttribute("data-base") || "";
  var records = [];
  window.ALCStations.load().then(function (recs) {
    records = (recs || []).map(function (r) {
      r._s = (r.n + " " + r.w + " " + r.k + " " + r.c + " " + r.t).toLowerCase();
      return r;
    });
    // the operator may have typed OR just clicked into the (still empty) box before the fetch
    // landed; in the empty case update() renders the neighbourhood, which needs these records
    if (box.value.trim() || document.activeElement === box) update();
  });
  var matches = [], sel = -1;

  function go(rec) { if (rec) window.location.href = base + "stations/" + rec.k + ".html"; }

  function rank(r, q) {  // prefix matches first, then name-contains, then anything
    var n = (r.n || "").toLowerCase(), w = (r.w || "").toLowerCase();
    if (n.indexOf(q) === 0) return 0;
    if (w.indexOf(q) === 0) return 1;
    if (n.indexOf(q) >= 0) return 2;
    return 3;
  }

  function render() {
    if (!matches.length) { panel.hidden = true; panel.innerHTML = ""; return; }
    panel.innerHTML = matches.map(function (r, i) {
      // status dot + constant as a percent of the type's nominal value, so the panel doubles as a
      // "browse the network" view instead of being search-only
      var cl = window.ALCStations.clLabel(r, window.ALCStations.bestMethod(r));
      return '<div class="sr' + (i === sel ? " sel" : "") + '" data-i="' + i + '">' +
        window.ALCStations.statusDot(r) + '<span class="sr-name">' + (r.n || r.w) + '</span>' +
        '<span class="sr-key">' + r.k + '</span>' +
        '<span class="sr-meta">' + [r.t, r.c].filter(Boolean).join(" · ") + '</span>' +
        (cl ? '<span class="sr-cl">' + cl + '</span>' : "") + '</div>';
    }).join("");
    panel.hidden = false;
    Array.prototype.forEach.call(panel.querySelectorAll(".sr"), function (el) {
      el.addEventListener("mousedown", function (e) { e.preventDefault(); go(matches[+el.getAttribute("data-i")]); });
    });
  }

  function update() {
    var q = box.value.trim().toLowerCase();
    sel = -1;
    if (!q) {
      // Empty box on focus = "browse the neighbours": show the stations around this one, honouring
      // the navbar filter if there is one, so browsing and searching are the same panel.
      matches = neighbourhood();
      render();
      return;
    }
    matches = records.filter(function (r) { return r._s.indexOf(q) >= 0; });
    matches.sort(function (a, b) {
      var ra = rank(a, q), rb = rank(b, q);
      if (ra !== rb) return ra - rb;
      return (a.n || a.k).localeCompare(b.n || b.k);
    });
    matches = matches.slice(0, 15);
    render();
  }

  function neighbourhood() {
    var navEl = document.getElementById("station-nav");
    var here = navEl && navEl.getAttribute("data-key");
    var fc = document.getElementById("f-country"), ft = document.getElementById("f-type");
    var c = fc ? fc.value : "", t = ft ? ft.value : "";
    var subset = records.filter(function (r) { return (!c || r.c === c) && (!t || r.t === t); });
    if (!subset.length) return [];
    var i = here ? subset.findIndex(function (r) { return r.k === here; }) : -1;
    if (i < 0) return subset.slice(0, 15);
    return subset.slice(Math.max(0, i - 7), Math.max(0, i - 7) + 15);
  }

  box.addEventListener("input", update);
  box.addEventListener("focus", update);
  box.addEventListener("keydown", function (e) {
    if (e.key === "ArrowDown") { sel = Math.min(sel + 1, matches.length - 1); render(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); render(); e.preventDefault(); }
    // Enter opens an explicit selection, or the best hit of a TYPED query. It must not open
    // anything from the empty-box browse list — before the browse list existed, Enter on an empty
    // box did nothing, and silently jumping to a station would be a nasty surprise.
    else if (e.key === "Enter") {
      if (sel >= 0) go(matches[sel]);
      else if (box.value.trim()) go(matches[0]);
      e.preventDefault();
    }
    else if (e.key === "Escape") { matches = []; render(); box.blur(); }
  });
  document.addEventListener("click", function (e) { if (!wrap.contains(e.target)) { matches = []; render(); } });
})();
