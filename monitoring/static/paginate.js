// Client-side pagination for long tables. A table with class "paginated" and
// data-page-size="N" shows N rows at a time with Prev/Next controls. No dependencies.
(function () {
  function paginate(table) {
    var size = parseInt(table.getAttribute("data-page-size") || "25", 10);
    var tbody = table.tBodies[0];
    if (!tbody) return;
    var rows = Array.prototype.slice.call(tbody.rows);
    if (rows.length <= size) return;                 // nothing to paginate
    var pages = Math.ceil(rows.length / size);
    var cur = 0;

    var nav = document.createElement("div");
    nav.className = "pager";
    var prev = document.createElement("button");
    prev.type = "button"; prev.textContent = "‹ Prev";
    var next = document.createElement("button");
    next.type = "button"; next.textContent = "Next ›";
    var label = document.createElement("span");
    label.className = "pager-label";
    nav.appendChild(prev); nav.appendChild(label); nav.appendChild(next);
    table.parentNode.insertBefore(nav, table.nextSibling);

    function render() {
      var start = cur * size, end = start + size;
      rows.forEach(function (r, i) { r.style.display = (i >= start && i < end) ? "" : "none"; });
      label.textContent = "Page " + (cur + 1) + " of " + pages + " · " + rows.length + " total";
      prev.disabled = cur === 0;
      next.disabled = cur === pages - 1;
    }
    prev.addEventListener("click", function () { if (cur > 0) { cur--; render(); } });
    next.addEventListener("click", function () { if (cur < pages - 1) { cur++; render(); } });
    render();
  }

  document.querySelectorAll("table.paginated").forEach(paginate);
})();
