// Client-side QC flagging for the static dashboard (no server). Clicking the red "flag" button in a
// diagnostic viewer records {key, wmo, identifier, itype, method, date} for the currently shown
// calibration. The list is kept in localStorage (primary) AND mirrored to a site-wide cookie, so it
// survives navigation between pages even in browsers/contexts where localStorage is not shared (e.g.
// some file:// setups). A floating widget shows the running count and exports the whole list as CSV.
// diag.js calls window.QCFlags.toggle().
(function () {
  var KEY = "alc_qc_flags";

  function readCookie() {
    var m = document.cookie.match(/(?:^|;\s*)alc_qc_flags=([^;]*)/);
    if (!m) { return null; }
    try { return JSON.parse(decodeURIComponent(m[1])); } catch (e) { return null; }
  }
  function writeCookie(a) {
    try {
      var v = encodeURIComponent(JSON.stringify(a));
      // cookies are capped near 4 KB; only mirror when it fits, else drop it (localStorage stays
      // authoritative and the widget still exports everything).
      if (a.length && v.length < 3800) {
        document.cookie = KEY + "=" + v + "; path=/; max-age=31536000; SameSite=Lax";
      } else {
        document.cookie = KEY + "=; path=/; max-age=0; SameSite=Lax";
      }
    } catch (e) { /* cookies disabled */ }
  }
  function load() {
    var ls = null;
    try { ls = JSON.parse(localStorage.getItem(KEY)); } catch (e) { ls = null; }
    if (ls && ls.length) { return ls; }
    // localStorage empty/unavailable on this page -> restore from the cookie mirror.
    var ck = readCookie();
    if (ck && ck.length) {
      try { localStorage.setItem(KEY, JSON.stringify(ck)); } catch (e) { /* ignore */ }
      return ck;
    }
    return ls || [];
  }
  function store(a) {
    try { localStorage.setItem(KEY, JSON.stringify(a)); } catch (e) { /* quota / private mode */ }
    writeCookie(a);
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
