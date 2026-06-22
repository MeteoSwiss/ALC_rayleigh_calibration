// Client-side pagination for long tables. A table with class "paginated" and
// data-page-size="N" shows N rows at a time with Prev/Next controls.
// Filter-aware: rows marked data-fhidden="1" (by filter.js) are excluded from the page set,
// and sorting (table-sort.js) re-reads DOM order. Each table gets a table._repaginate() hook.
(function () {
  function paginate(table) {
    var size = parseInt(table.getAttribute("data-page-size") || "25", 10);
    var tbody = table.tBodies[0];
    if (!tbody) return;
    var cur = 0;

    var nav = document.createElement("div");
    nav.className = "pager";
    var prev = document.createElement("button"); prev.type = "button"; prev.textContent = "‹ Prev";
    var next = document.createElement("button"); next.type = "button"; next.textContent = "Next ›";
    var label = document.createElement("span"); label.className = "pager-label";
    nav.appendChild(prev); nav.appendChild(label); nav.appendChild(next);
    table.parentNode.insertBefore(nav, table.nextSibling);

    function active() {  // rows in current DOM order (respects sort), minus filter-hidden ones
      return Array.prototype.filter.call(tbody.rows, function (r) {
        return r.getAttribute("data-fhidden") !== "1";
      });
    }
    function render() {
      var rows = active();
      var pages = Math.max(1, Math.ceil(rows.length / size));
      if (cur >= pages) cur = pages - 1;
      if (cur < 0) cur = 0;
      Array.prototype.forEach.call(tbody.rows, function (r) {
        if (r.getAttribute("data-fhidden") !== "1") r.style.display = "none";
      });
      var start = cur * size;
      rows.slice(start, start + size).forEach(function (r) { r.style.display = ""; });
      label.textContent = "Page " + (cur + 1) + " of " + pages + " · " + rows.length +
        (rows.length === 1 ? " row" : " rows");
      prev.disabled = cur === 0;
      next.disabled = cur >= pages - 1;
      nav.style.display = rows.length > size ? "" : "none";
    }
    prev.addEventListener("click", function () { if (cur > 0) { cur--; render(); } });
    next.addEventListener("click", function () { cur++; render(); });
    table._repaginate = function () { cur = 0; render(); };  // reset to page 1 + redraw
    table._render = render;
    render();
  }

  document.querySelectorAll("table.paginated").forEach(paginate);
})();
