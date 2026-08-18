/* Headline tiles follow the C_L overlay's legend.
 *
 * The tiles default to the POOLED statistics (both retrievals estimate the same constant); when
 * the operator isolates one method by clicking the other out of the overlay's legend, the tiles
 * re-render from that method's own set. Server renders the combined set into the DOM, so the page
 * is complete without JS; this script only swaps between the precomputed sets.
 */
(function () {
  "use strict";
  var dataEl = document.getElementById("cl-tiles-data");
  var host = document.getElementById("cl-tiles");
  if (!dataEl || !host) return;
  var S;
  try { S = JSON.parse(dataEl.textContent); } catch (e) { return; }

  function esc(x) { return String(x).replace(/&/g, "&amp;").replace(/</g, "&lt;"); }
  function paint(key) {
    var tiles = S[key] || S.combined;
    if (!tiles) return;
    host.innerHTML = tiles.map(function (t) {
      return '<div class="tile"><div class="tl">' + esc(t.label) + '</div>' +
             '<div class="tv" style="color:' + esc(t.color) + '">' + esc(t.value) + '</div>' +
             '<div class="tn">' + esc(t.note) + '</div></div>';
    }).join("");
  }

  var gd = document.getElementById("fig-overlay");
  if (!gd || !gd.on) return;                 // one-method stream: combined already rendered
  gd.on("plotly_restyle", function () {
    var vis = [];
    (gd.data || []).forEach(function (tr) {
      if (tr.visible === "legendonly") return;
      if (tr.name === "Rayleigh") vis.push("rayleigh");
      if (tr.name === "Liquid-cloud") vis.push("cloud");
    });
    paint(vis.length === 1 ? vis[0] : "combined");
  });
})();
