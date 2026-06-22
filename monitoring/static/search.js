// Nav-bar station search. Matches a free-text query against each station's name, WIGOS id and
// key (substring, case-insensitive), so "Payerne", "Pay", "0-20000-0-06610" and "06610" all find
// the 0-20000-0-06610_A/_B/_C pages. Arrow keys move the selection, Enter opens it, click opens it.
// Station index is embedded as JSON (#search-index); navigation is relative to data-base.
(function () {
  var box = document.getElementById("station-search");
  var panel = document.getElementById("search-results");
  var idxEl = document.getElementById("search-index");
  var wrap = document.querySelector(".navsearch");
  if (!box || !panel || !idxEl || !wrap) return;
  var base = wrap.getAttribute("data-base") || "";
  var records;
  try { records = JSON.parse(idxEl.textContent); } catch (e) { return; }
  records.forEach(function (r) {
    r._s = (r.name + " " + r.wigos + " " + r.key + " " + r.country + " " + r.type).toLowerCase();
  });
  var matches = [], sel = -1;

  function go(rec) { if (rec) window.location.href = base + "stations/" + rec.key + ".html"; }

  function rank(r, q) {  // prefix matches first, then name-contains, then anything
    var n = r.name.toLowerCase(), w = r.wigos.toLowerCase();
    if (n.indexOf(q) === 0) return 0;
    if (w.indexOf(q) === 0) return 1;
    if (n.indexOf(q) >= 0) return 2;
    return 3;
  }

  function render() {
    if (!matches.length) { panel.hidden = true; panel.innerHTML = ""; return; }
    panel.innerHTML = matches.map(function (r, i) {
      return '<div class="sr' + (i === sel ? " sel" : "") + '" data-i="' + i + '">' +
        '<span class="sr-name">' + (r.name || r.wigos) + '</span>' +
        '<span class="sr-key">' + r.key + '</span>' +
        '<span class="sr-meta">' + [r.type, r.country].filter(Boolean).join(" · ") + '</span></div>';
    }).join("");
    panel.hidden = false;
    Array.prototype.forEach.call(panel.querySelectorAll(".sr"), function (el) {
      el.addEventListener("mousedown", function (e) { e.preventDefault(); go(matches[+el.getAttribute("data-i")]); });
    });
  }

  function update() {
    var q = box.value.trim().toLowerCase();
    sel = -1;
    if (!q) { matches = []; render(); return; }
    matches = records.filter(function (r) { return r._s.indexOf(q) >= 0; });
    matches.sort(function (a, b) {
      var ra = rank(a, q), rb = rank(b, q);
      if (ra !== rb) return ra - rb;
      return (a.name || a.key).localeCompare(b.name || b.key);
    });
    matches = matches.slice(0, 15);
    render();
  }

  box.addEventListener("input", update);
  box.addEventListener("focus", update);
  box.addEventListener("keydown", function (e) {
    if (e.key === "ArrowDown") { sel = Math.min(sel + 1, matches.length - 1); render(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); render(); e.preventDefault(); }
    else if (e.key === "Enter") { go(matches[sel >= 0 ? sel : 0]); e.preventDefault(); }
    else if (e.key === "Escape") { matches = []; render(); box.blur(); }
  });
  document.addEventListener("click", function (e) { if (!wrap.contains(e.target)) { matches = []; render(); } });
})();
