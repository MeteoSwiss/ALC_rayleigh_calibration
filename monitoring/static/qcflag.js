// Client-side QC flagging for the static dashboard (no server). Clicking the red "flag" button in a
// diagnostic viewer records {key, wmo, identifier, itype, method, date} for the currently shown
// calibration into the browser's localStorage. A floating widget shows the running count and exports
// the whole list as CSV, to later refine the QC on those dates. diag.js calls window.QCFlags.toggle().
(function () {
  var KEY = "alc_qc_flags";

  function load() {
    try { return JSON.parse(localStorage.getItem(KEY)) || []; } catch (e) { return []; }
  }
  function store(a) {
    try { localStorage.setItem(KEY, JSON.stringify(a)); } catch (e) { /* quota / private mode */ }
    render();
  }
  function rid(r) { return r.key + "|" + r.method + "|" + r.date; }   // identity = station + method + day

  var API = {
    has: function (r) { var k = rid(r); return load().some(function (x) { return rid(x) === k; }); },
    toggle: function (r) {
      var a = load(), k = rid(r), i = -1;
      a.forEach(function (x, j) { if (rid(x) === k) i = j; });
      if (i >= 0) { a.splice(i, 1); } else { a.push(r); }
      store(a);
      return i < 0;   // true => now flagged
    },
    list: load,
    clear: function () { if (confirm("Clear all " + load().length + " QC flags?")) store([]); },
    exportCsv: function () {
      var a = load();
      if (!a.length) { return; }
      var cols = ["key", "wmo", "identifier", "itype", "method", "date"];
      var esc = function (v) {
        v = (v == null ? "" : "" + v);
        return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
      };
      var lines = [cols.join(",")].concat(a.map(function (r) {
        return cols.map(function (c) { return esc(r[c]); }).join(",");
      }));
      var blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
      var url = URL.createObjectURL(blob), link = document.createElement("a");
      link.href = url; link.download = "qc_flags.csv";
      document.body.appendChild(link); link.click(); document.body.removeChild(link);
      URL.revokeObjectURL(url);
    }
  };

  function render() {
    var n = load().length, w = document.getElementById("qcflag-widget");
    if (!w) { w = document.createElement("div"); w.id = "qcflag-widget"; document.body.appendChild(w); }
    if (!n) { w.style.display = "none"; return; }
    w.style.display = "flex";
    w.innerHTML =
      '<span class="qcw-n" title="dates flagged for QC review (saved in this browser)">⚑ ' + n + ' flagged</span>' +
      '<button type="button" class="qcw-exp">export CSV</button>' +
      '<button type="button" class="qcw-clr" title="clear all flags">clear</button>';
    w.querySelector(".qcw-exp").addEventListener("click", API.exportCsv);
    w.querySelector(".qcw-clr").addEventListener("click", API.clear);
  }

  window.QCFlags = API;
  if (document.readyState !== "loading") { render(); }
  else { document.addEventListener("DOMContentLoaded", render); }
})();
