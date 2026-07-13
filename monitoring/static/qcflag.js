// Client-side QC flagging for the static dashboard (no server). The red "flag" button in a diagnostic
// viewer opens a comment dialog (or diag.js quick-flags via shortcuts) and records
// {key, wmo, identifier, itype, method, date, comment} for the currently shown calibration. The list
// is kept in localStorage (primary) AND mirrored to a site-wide cookie, so it survives navigation
// between pages even where localStorage is not shared (e.g. some file:// setups). A floating widget
// shows the running count and exports the whole list (with comments) as CSV. diag.js calls
// window.QCFlags.openDialog() / .set(); the three preset comments live in window.QCFlags.PRESETS.
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
  var FIELDS = ["key", "wmo", "identifier", "itype", "method", "date"];

  // Preset QC comments offered in the dialog and bound to the 1/2/3 quick-flag shortcuts (diag.js).
  var PRESETS = ["Aerosol contamination", "Cloud contamination", "Low signal (condensation?)"];

  function _attr(s) { return ("" + s).replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;"); }

  function closeDialog() {
    var ov = document.getElementById("qcm-overlay");
    if (ov && ov.parentNode) { ov.parentNode.removeChild(ov); }
    document.removeEventListener("keydown", _dialogKeys, true);
  }
  function _dialogKeys(e) {
    if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); closeDialog(); }
  }
  // Modal: pick a preset (one click = flag with that comment), type a custom note, remove, or cancel.
  function openDialog(rec, onChange) {
    closeDialog();   // never stack dialogs
    var existing = API.get(rec);
    var ds = rec.date ? (rec.date.slice(0, 4) + "-" + rec.date.slice(4, 6) + "-" + rec.date.slice(6, 8)) : "";
    var presetBtns = PRESETS.map(function (p, i) {
      return '<button type="button" class="qcm-preset" data-c="' + _attr(p) + '">' +
             '<b>' + (i + 1) + '</b> · ' + _attr(p) + '</button>';
    }).join("");
    var ov = document.createElement("div");
    ov.id = "qcm-overlay"; ov.className = "qcm-overlay";
    ov.innerHTML =
      '<div class="qcm-dialog" role="dialog" aria-modal="true" aria-label="Flag calibration">' +
        '<div class="qcm-head">⚑ Flag calibration for QC review</div>' +
        '<div class="qcm-sub">' + _attr(rec.key) + ' · ' + _attr(rec.method) + ' · ' + ds + '</div>' +
        '<div class="qcm-presets">' + presetBtns + '</div>' +
        '<label class="qcm-label" for="qcm-text">Or a custom note</label>' +
        '<textarea id="qcm-text" class="qcm-text" rows="3" placeholder="Type a comment…"></textarea>' +
        '<div class="qcm-actions">' +
          (existing ? '<button type="button" class="qcm-remove">Remove flag</button>' : '<span></span>') +
          '<span class="qcm-right">' +
            '<button type="button" class="qcm-cancel">Cancel</button>' +
            '<button type="button" class="qcm-save">Save flag</button>' +
          '</span>' +
        '</div>' +
      '</div>';
    document.body.appendChild(ov);
    var ta = ov.querySelector(".qcm-text");
    if (existing && existing.comment) { ta.value = existing.comment; }
    function done() { closeDialog(); if (onChange) { onChange(); } }
    Array.prototype.forEach.call(ov.querySelectorAll(".qcm-preset"), function (b) {
      b.addEventListener("click", function () { API.set(rec, b.getAttribute("data-c")); done(); });
    });
    ov.querySelector(".qcm-save").addEventListener("click", function () { API.set(rec, ta.value.trim()); done(); });
    ov.querySelector(".qcm-cancel").addEventListener("click", closeDialog);
    var rm = ov.querySelector(".qcm-remove");
    if (rm) { rm.addEventListener("click", function () { API.remove(rec); done(); }); }
    ov.addEventListener("mousedown", function (e) { if (e.target === ov) { closeDialog(); } });  // click backdrop
    document.addEventListener("keydown", _dialogKeys, true);
    ta.focus();
  }

  var API = {
    has: function (r) { var k = rid(r); return load().some(function (x) { return rid(x) === k; }); },
    get: function (r) { var k = rid(r); var m = load().filter(function (x) { return rid(x) === k; }); return m.length ? m[0] : null; },
    set: function (r, comment) {
      var a = load(), k = rid(r), found = false;
      a.forEach(function (x) { if (rid(x) === k) { x.comment = comment || ""; found = true; } });
      if (!found) { var rec = {}; FIELDS.forEach(function (c) { rec[c] = r[c]; }); rec.comment = comment || ""; a.push(rec); }
      store(a);
    },
    remove: function (r) { var k = rid(r); store(load().filter(function (x) { return rid(x) !== k; })); },
    toggle: function (r) {   // kept for completeness; the UI now flags through the comment dialog
      if (API.has(r)) { API.remove(r); return false; }
      API.set(r, ""); return true;
    },
    list: load,
    PRESETS: PRESETS,
    openDialog: openDialog,
    isDialogOpen: function () { return !!document.getElementById("qcm-overlay"); },
    clear: function () { if (confirm("Clear all " + load().length + " QC flags?")) store([]); },
    exportCsv: function () {
      var a = load();
      if (!a.length) { return; }
      var cols = ["key", "wmo", "identifier", "itype", "method", "date", "comment"];
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
