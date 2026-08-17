# 12 — Dashboard professionalization: audit findings and plan

**Date** 2026-08-17 · **Scope** the station dashboard (production `monitoring/` renderer + the
interactive daily-calibration panel embedded into it) · **Method** six-lens multi-agent audit
(CSS collisions, JS controls, data/labels, duplication, UX, pipeline), 87 raw findings deduplicated,
the top findings adversarially verified (8 CONFIRMED, 2 REFUTED), plus a completeness critic.
Verification transcripts: session workflow `wf_ec93a10b-05f`.

The context that made the audit necessary: the panel was developed as a standalone mockup
(`scripts/mockup_daily_panel.py`) and then embedded into the production page. Most findings are
consequences of exactly that — two generations of code sharing one page and one output directory.
An operator-visible regression treadmill was the symptom; `tests/test_dashboard_site.py`
(25 tests, each annotated with the regression it guards) is the countermeasure already in place.

## Already fixed during the review (for the record)

`.msg` 360 px CSS collision → `.dp-verdict`; clear-day cloud rejections misattributed to
−22/−23/−24 (cloudless profiles no longer vote; `no_cloud_rejected`); hard-coded "v2.0" label →
`version` column end-to-end; "Signal / range²" → "Signal (not range-corrected)" with the view
set renamed/ordered Raw S → Raw RCS → Att. backscatter; rejected days rendered grey
(`CAL_COL.none` dangling reference + `hasFig` blind to index summaries); month navigation
snap-back; payload worker race (per-PID work dirs) and corrupt-payload self-healing; calendar
rail in the page margin with narrow-screen overlay; one-line header; availability card always
present (real octas incl. `_B`); OmB → Classification → Sensitivity always at the bottom with
bucket fallback; vertical-segment calibration-window chart; Ctrl+←/→ = valid calibrations.

## P0 — the embedding is subtly broken (fix before anything else)

| finding | file | why it matters |
|---|---|---|
| `watchViewer()` runs before the diag sections exist and is never retried | `dailypanel.js:105` | the panel↔image-viewer sync NEVER attaches — the script is included mid-page, above the method blocks |
| PANEL_JS still ships its own arrow keydown into the production page; the `.diag-data` stand-down guard is defeated by the classification section (which carries `.diag-data` with no calibration viewer) | `mockup_daily_panel.py` | arrows dead on some stations, double-stepping on others — consolidate ALL embedded keyboard into `dailypanel.js`, keep PANEL_JS keys for the standalone mockup only (explicit handshake flag, not DOM sniffing) |
| `#period-sel` relayouts the panel's Plotly figures with station-history date ranges | `rangesync.js:27` | corrupts the one-night curtain axes — skip figures inside `#daily` |
| keydown guards miss SELECT/TEXTAREA | `mockup_daily_panel.py` | arrows inside the period dropdown step the panel |
| fetch-failure warning overwritten by the next render, then latched off | `dailypanel.js:29` | render the warning outside `#panel` |
| `_write_payloads` rebuilds `_index.json` from the RUN's date range only | `build_station_dashboard.py:184` | an incremental daily run would wipe the panel's history — merge into the existing index (the 11-day-calendar incident was this bug wearing another hat) |
| unscoped `PANEL_CSS` (`h2`, `body`, `:root`, `*`, `button`, `.card`, `.empty`) restyles the host page | `mockup_daily_panel.py:812-880` | every heading on panel-bearing pages is 13 px uppercase grey; `:root` hijacks `--line/--ink`; `.empty` will balloon the diag calendar's 17 px cells the day both render together — scope everything under `#daily`, tokens → `--dp-*`, split standalone-only chrome; **extend the collision test to element selectors and custom properties (it was defeated by an allowlist and by matching classes only)** |

## P1 — single source of truth

* **Promote the panel out of `scripts/mockup_*`**: `monitoring/panel.py` (fragments + payload
  builder), `scripts/` keeps thin CLIs. Production code importing from a mockup module is the
  root enabler of the two-generations problem.
* **Reduce `build_station_dashboard.py` to a payloads-only tool** (its `--no-pages` mode becomes
  the only mode): its page generator, hand-copied topbar (dead search box), `stations.json`
  writer (incompatible schema), prompt-based QC store and hard-coded PERIODS list are all
  duplicate implementations the audit flagged; `monitoring/render.py` is canonical.
* **One palette, one wording**: three legends disagree on the colours of "calibrated"/"rejected"
  (panel calendar vs availability card vs flag dots) — drive all from `monitoring/config.py`;
  panel flag labels through `config.flag_label(flag, method)`; one C_L precision (4 sig figs)
  via a shared formatter.
* **Scientific label bug (verified)**: the calibration-window chart plots `bottom/top_height`
  (ASL) under an axis that says **m AGL** (`charts.py:519`) — relabel or subtract the station
  altitude so it matches the panel's AGL axis. *(Y axis stays altitude/range, per house rule.)*
* **Tiles**: `LAST VALID` ages against wall-clock `utcnow` — every station of an archived build
  reads red; age against the archive's max date passed in as `as_of`. Single implementation in
  `metrics.py` (the `_cl_series` fork in the old generator dies with P1).
* **Flag documentation**: −10 and −11 absent from `FLAG_DOCS` (dead anchors from station tables)
  and −11 from `FLAG_COLORS`; `flags.html` claims −1/0 are excluded from success denominators
  while every computed rate includes them — fix the text to the implemented definition.

## P2 — operational integration (blocking for going live)

* **`ops_daily.py` never regenerates payloads**: the panel freezes at the last manual run. Add a
  post-calibration step: compute D-1's `<ds>_<method>.json` for panel-enabled stations and
  merge-reindex. Costs ~0.2 s/station-night measured.
* **`publish.sh`**: `data/<key>/` payloads (tens of MB per stream) currently ride the ssh HTML
  leg and `--delete` semantics need auditing — publish payloads via the S3 bucket leg (fetch from
  `ALC_IMG_BASE_URL/data/...`), or generate them server-side only.
* **`--changed-only` is blind to payload updates** and `.last_build` is stamped at build END
  (missing CSVs written during the render) — mark keys changed when `_index.json` is newer, stamp
  the marker with the pre-index time.
* `_write_assets` rewrites every asset (incl. 4.7 MB plotly) with fresh mtimes → route through
  `_write_if_changed`. Daily rebuild silently drops the curated flag-example images (`--flagex`
  never passed by ops) → `ALC_FLAGEX_DIR` in `config.sh`.

## P3 — UX polish and accessibility

Dead affordance texts ("click a day…" where nothing is wired) → wire `fig-avail` clicks and table
dates to `goDay()`; URL-hash for the selected night (linkable, back/forward works); day cells as
`<button>` with `aria-current` (keyboard/screen-reader access); wide tables in `overflow-x:auto`
containers; rail offset under the sticky topbar; explicit empty-state for `_B`-style streams'
monitoring section; **the promised hourly instrument monitoring (`<key>_hk_hourly.csv`, separate
file, `_hour_of` binning, fetch-on-zoom) — confirmed absent, still owed**; browser-tier tests
automated (playwright: console-clean, one-keypress-one-action, period selector, `file://` warning).

## Refuted by verification (do not "fix")

* "Two parallel QC-flag stores in the embedded panel" — the embedded panel has no flag UI;
  `qcflag.js` is the single store. (The duplicate lived in the old generator's page, which P1 deletes.)
* One duplicate-implementation claim counted the same code path twice.

## Sequencing

P0 in one sitting (all are small, all are tested afterwards by extending the existing suite);
P1 as one refactor commit (move + delete, no behaviour change, suite green before/after);
P2 alongside the v2.2 cutover work — the ops-payload step and publish routing must exist before
the panel goes live network-wide; P3 continuously. Every fix lands with its regression test —
that is the actual professionalization.
