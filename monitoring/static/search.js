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
  var matches = [], sel = -1, hereKey = null;   // hereKey: current station, highlighted in browse mode

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
      return '<div class="sr' + (i === sel ? " sel" : "") + (r.k === hereKey ? " here" : "") +
        '" data-i="' + i + '">' +
        window.ALCStations.statusDot(r) + '<span class="sr-name">' + (r.n || r.w) + '</span>' +
        '<span class="sr-key">' + r.k + '</span>' +
        '<span class="sr-meta">' + [r.t, r.c].filter(Boolean).join(" · ") + '</span>' +
        (cl ? '<span class="sr-cl">' + cl + '</span>' : "") + '</div>';
    }).join("");
    panel.hidden = false;
    Array.prototype.forEach.call(panel.querySelectorAll(".sr"), function (el) {
      el.addEventListener("mousedown", function (e) { e.preventDefault(); go(matches[+el.getAttribute("data-i")]); });
    });
    // Keep the interesting row in view of the (wheel-scrollable) panel. Manual scrollTop, not
    // scrollIntoView: the latter also scrolls every ancestor, i.e. the page itself.
    var tgt = panel.querySelector(".sr.sel");
    if (tgt) {   // arrow keys: minimal scroll to keep the selection visible
      if (tgt.offsetTop < panel.scrollTop) panel.scrollTop = tgt.offsetTop;
      else if (tgt.offsetTop + tgt.offsetHeight > panel.scrollTop + panel.clientHeight)
        panel.scrollTop = tgt.offsetTop + tgt.offsetHeight - panel.clientHeight;
    } else if ((tgt = panel.querySelector(".sr.here"))) {   // browse mode: centre this station
      panel.scrollTop = Math.max(0, tgt.offsetTop - (panel.clientHeight - tgt.offsetHeight) / 2);
    } else {
      panel.scrollTop = 0;
    }
  }

  function update() {
    var q = box.value.trim().toLowerCase();
    sel = -1;
    if (!q) {
      // Empty box on focus = "browse the network": the WHOLE scoped station list, wheel-scrollable,
      // centred on (and highlighting) the station being viewed. The old 15-row neighbourhood window
      // made the scrollbar pointless -- there was nothing beyond it to scroll to.
      matches = browseList();
      render();
      return;
    }
    hereKey = null;
    matches = records.filter(function (r) { return r._s.indexOf(q) >= 0; });
    matches.sort(function (a, b) {
      var ra = rank(a, q), rb = rank(b, q);
      if (ra !== rb) return ra - rb;
      return (a.n || a.k).localeCompare(b.n || b.k);
    });
    render();
  }

  function browseList() {
    var navEl = document.getElementById("station-nav");
    var here = navEl && navEl.getAttribute("data-key");
    var fc = document.getElementById("f-country"), ft = document.getElementById("f-type");
    var c = fc ? fc.value : "", t = ft ? ft.value : "";
    var subset = records.filter(function (r) { return (!c || r.c === c) && (!t || r.t === t); });
    hereKey = (here && subset.some(function (r) { return r.k === here; })) ? here : null;
    return subset;
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
