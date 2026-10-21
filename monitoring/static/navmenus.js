// Station-topbar panels (design 2B): the ⇩ Data and ⚙ Filters dropdowns, the active-filter count
// badge, and the period pills.
//
// The pills are a PROXY for the hidden #period-sel select — rangesync.js owns that select (window
// relayout of every figure + the "alc-period" sessionStorage), so a pill click just writes the
// select and fires its change event; nothing about the period machinery is duplicated here. The
// country/type scope selects likewise stay owned by stationnav.js; this file only counts them for
// the badge. No-op on pages without these elements (summary, flags).
(function () {
  "use strict";

  // ---- dropdown open/close (one open at a time; outside click / Escape closes) ----
  var drops = [];
  function bindDrop(btnId, panelId) {
    var btn = document.getElementById(btnId);
    var panel = document.getElementById(panelId);
    if (!btn || !panel) return;
    drops.push({ btn: btn, panel: panel });
    btn.addEventListener("click", function (e) {
      e.stopPropagation();
      var open = panel.hidden;
      closeAll();
      if (open) { panel.hidden = false; btn.setAttribute("aria-expanded", "true"); }
    });
    panel.addEventListener("click", function (e) { e.stopPropagation(); });
  }
  function closeAll() {
    drops.forEach(function (d) {
      d.panel.hidden = true;
      d.btn.setAttribute("aria-expanded", "false");
    });
  }
  document.addEventListener("click", closeAll);
  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closeAll(); });
  bindDrop("nv-data-btn", "nv-data-panel");
  bindDrop("nv-filters-btn", "nv-filters-panel");

  // ---- period pills -> hidden #period-sel ----
  var sel = document.getElementById("period-sel");
  var pills = document.getElementById("period-pills");
  function markPill(key) {
    if (!pills) return;
    Array.prototype.forEach.call(pills.querySelectorAll(".nvpill"), function (b) {
      b.classList.toggle("on", b.getAttribute("data-period") === key);
    });
  }
  if (sel && pills) {
    pills.addEventListener("click", function (e) {
      var b = e.target.closest ? e.target.closest("[data-period]") : null;
      if (!b) return;
      sel.value = b.getAttribute("data-period");
      sel.dispatchEvent(new Event("change"));   // rangesync.js applies + persists the window
      markPill(sel.value);
      updateCount();
    });
    // Restore the remembered window's pill. Read sessionStorage directly rather than the select:
    // rangesync.js loads AFTER this file, so at this point the select still says its baked default.
    var saved = null;
    try { saved = sessionStorage.getItem("alc-period"); } catch (e) {}
    if (saved && pills.querySelector('[data-period="' + saved + '"]')) markPill(saved);
    else markPill(sel.value);
  }

  // ---- active-filter count badge on the ⚙ button ----
  var badge = document.getElementById("nv-filter-count");
  var fc = document.getElementById("f-country");
  var ft = document.getElementById("f-type");
  function periodKey() {
    var saved = null;
    try { saved = sessionStorage.getItem("alc-period"); } catch (e) {}
    return saved || (sel ? sel.value : "all");
  }
  function updateCount() {
    if (!badge) return;
    var n = 0;
    if (fc && fc.value) n++;
    if (ft && ft.value) n++;
    if (sel && periodKey() !== "all") n++;
    badge.hidden = n === 0;
    badge.textContent = String(n);
  }
  if (fc) fc.addEventListener("change", updateCount);
  if (ft) ft.addEventListener("change", updateCount);
  if (sel) sel.addEventListener("change", updateCount);
  // stationnav.js restores the stored scope into the selects synchronously before this runs
  // (script order in base.html), so counting on load sees the restored values.
  updateCount();
})();
