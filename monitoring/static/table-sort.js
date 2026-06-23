// Minimal click-to-sort for the station table. No dependencies.
// A header cell with data-type="num" sorts on each row's matching cell, preferring a
// `data-val` attribute (set server-side for correct numeric ordering of formatted text).
(function () {
  function cellValue(row, idx, numeric) {
    var td = row.children[idx];
    if (!td) return numeric ? -Infinity : "";
    var raw = td.getAttribute("data-val");
    if (raw === null) raw = td.textContent.trim();
    return numeric ? (parseFloat(raw) || 0) : raw.toLowerCase();
  }

  function makeSortable(table) {
    var headers = table.tHead ? table.tHead.rows[0].cells : [];
    Array.prototype.forEach.call(headers, function (th, idx) {
      var type = th.getAttribute("data-type");
      if (!type) return;                       // columns without data-type are not sortable
      th.addEventListener("click", function () {
        var numeric = type === "num";
        var asc = !th.classList.contains("sort-asc");
        Array.prototype.forEach.call(headers, function (h) {
          h.classList.remove("sort-asc", "sort-desc");
        });
        th.classList.add(asc ? "sort-asc" : "sort-desc");
        var tbody = table.tBodies[0];
        var rows = Array.prototype.slice.call(tbody.rows);
        rows.sort(function (a, b) {
          var va = cellValue(a, idx, numeric), vb = cellValue(b, idx, numeric);
          if (va < vb) return asc ? -1 : 1;
          if (va > vb) return asc ? 1 : -1;
          return 0;
        });
        rows.forEach(function (r) { tbody.appendChild(r); });
        if (typeof table._repaginate === "function") table._repaginate();  // re-page in new order
      });
    });
  }

  document.querySelectorAll("table.sortable").forEach(makeSortable);
})();
