# Plan — Cloudnet target classification for the ALC network (updated)

Run `ceiloclass` (CloudnetPy-style aerosol / cloud / ice / drizzle classification) for every
E-PROFILE stream daily, driven from **CAMS** temperature, riding on the **read-once** shared load,
with the curtains shown on the dashboard.

**Status: the hard prerequisite (read-once) is done and exceeded. What's left is wiring the
already-validated classifier into the operational pipeline.**

---

## Done (the foundation)

- **Read-once shared load** — `calibration/io/instrument_day.py::load_instrument_day` reads L1 once
  + the CAMS closest cell once, builds the coarse **working grid** (30 s / 10 m) and the WV
  transmission, and every calibration pass (Rayleigh / cloud / OmB / sensitivity) now consumes it.
  → classification's inputs (β on the working grid, CAMS temperature) are **already in memory**, so
  a classification consumer is near-free. This was the ~1–1.5 week bulk of the original estimate.
- **Resolution normalization** — the working grid already applies the ≥10 m / ≥30 s coarsening the
  classifier wants (only CL61 4.8 m range actually coarsens).
- **The port + adapter exist and are validated** — `Python/ceiloclass/payerne_eprofile_test/`:
  `eprofile_l1.py` (E-PROFILE L1 → `ceilopyter`-style β via rcs_0 × per-instrument factor),
  `cams_model.py::cams_to_model` (CAMS `Model` from the L137 temperature via the repo's
  `water_vapor.cams_temperature_pressure_profile` — **still works after the WV unification**),
  run on 5 Payerne days for CL31 / CHM15k / CL61 with cross-instrument agreement.
- **Q5 value confirmed** — on CL61 Payerne 2026-03-06 the classifier flags the persistent
  depol-ice layer (2.5–3.7 km, 0.2 depol) that slips the Rayleigh scattering-ratio gate.

## Left (the integration)

| # | Task | Detail | Effort |
|---|------|--------|--------|
| 1 | ✅ **Dependencies** | `ceiloclass` + `ceilopyter` (+ atmoslib, cloudnet-api-client) installed into `.venv`; recorded as the `classify` optional extra in `pyproject.toml`. Install: `pip install "ceiloclass @ git+https://github.com/actris-cloudnet/ceiloclass.git"`. **Still needed:** the same install in the ops `$ALC_VENV`. | done |
| 2 | ✅ **Adapter in the repo** | `calibration/classify/` — `eprofile_l1.read_eprofile_l1` (file) + `ceilo_from_shared(idd)` (read-once working grid, incl. CL61 depol) + `cams_model.cams_to_model` (CAMS→`Model` via the unified WV reader). The loader now carries `depol`. Validated: shared path == file path to ~0.02%/class and both detect the Q5 depol-ice contamination (`tests/test_classify.py`). | done (commit 01b6a0c) |
| 3 | ✅ **`_do_classification(shared)` consumer** | `run_network_calibration.py --classify`: per stream-day `ceilo_from_shared(idd.slice_to_date(d))` → `cams_to_model(...)` → `classify(..., use_wet_bulb=False)` → writes a classification NetCDF (+ curtain PNG when PLOTS=1) under `<key>/classification/<wmo>/<year>/`. Flag-gated; composes with `--no-cal`/`--omb`/`--sens`; optional-dep + CAMS-less days degrade gracefully. Validated on the real runner (Payerne CL61 2026-03-06, ICE ~2.1%). | done (commit 71a6c4a) |
| 4 | ✅ **Storage / publish** | The NetCDF + PNG land under `<key>/classification/`; `--classify` skips days whose NetCDF exists (resume, `--force` overrides); `ops/publish.sh` prunes the curtain PNGs to the bucket (they stage into `diag/`, the NetCDF stays as the marker). | done (180df44, 6244aee) |
| 5 | ✅ **Dashboard** | Per-day "Cloudnet target classification" curtain card on the station page, reusing the diagnostic day-viewer (calendar + prev/next-day). Verified in-browser (card renders, curtain loads, no JS errors). | done (commit 6244aee) |
| 6 | ✅ **Q5 — Rayleigh screening** | `calibrate_rayleigh(contam_profile=…)` post-fit gate → new flag **-11** when >30% of the selected molecular window is classified cloud/ice over the night (catches the CL61 depol-ice the scattering-ratio gate misses); `_do_rayleigh` reads the d-1+d classification when `--classify`. No-op without the classification product. Validated (fires over the window, not outside/none). | done (commit 3eedf45) |
| 7 | ✅ **Ops enablement (code)** | `ops_daily.py` appends `--classify` when `ALC_CLASSIFY=1` (`ops/config.sh`, default 0); the runner degrades gracefully if the deps/CAMS are absent. Smoke-tested end-to-end (Payerne CL61 2026-03-16 → valid classification NetCDF + curtain PNG under `<key>/classification/`). | done |
| — | **Deploy (server, manual)** | On `zueub434`: `source ops/config.sh; source "$ALC_VENV/bin/activate"; pip install "ceiloclass @ git+https://github.com/actris-cloudnet/ceiloclass.git"` (needs `https_proxy`), then set `ALC_CLASSIFY=1` in `ops/config.sh`. The 15:00 cron then classifies daily. | ~0.5 d |

**All code is done.** Only the server-side dependency install + flipping `ALC_CLASSIFY=1` remain (manual, on the ops VM).

### CAMS temperature — time-varying (fixed 2026-07-11)

The classifier's temperature (melting layer / ice-vs-liquid) now uses **CAMS at its native 3-hourly
resolution, interpolated to each observation profile** (`cams_model.cams_to_model`), instead of a
day-mean broadcast to every step. On a frontal day (Payerne 2026-01-27) the real 0 °C level moves
~1.5 km through the day; a constant line mislabels the evening precip (snow read as rain). Separately,
`water_vapor._cams_levels` now averages CAMS on **geometric height** rather than model-level index
(the latter is only exact at constant surface pressure); impact on the Rayleigh constant is negligible
(≤0.25 % WV, 0.002 % molecular density) but it is the correct aggregation. Cloud WV and OmB were
already per-time.

## Notes / gotchas

- **910 nm nights without usable CAMS** — classification needs temperature only for the melting
  layer; a WV-free night can still classify (unlike the mandatory WV correction for the 910 nm cloud
  calibration). Degrade gracefully rather than skip.
- **Compute** — classification is ~12–18 s/instrument-day, the new bottleneck (calibration is 1–4 s).
  For the full network that's the dominant cost; keep it flag-gated and parallel per stream.
- **Near-surface CAMS temperature** — the coarse-orography cold bias (see the WV-resolution study)
  is irrelevant to the Rayleigh fit (3.5–5 km) but matters for boundary-layer classification; the ISA
  lapse fill below the lowest level is already handled in `cams_to_model`.
- The classifier only reads the shared β / temperature view — it **cannot shift any calibration
  number** (except deliberately, via task 6).

## Recommendation

The shared loader makes task 3 cheap, so the original "standalone Phase A batch first" is no longer
needed — go straight to the runner consumer. Order: **1 → 2** (unblocks everything) → **3 → 4**
(classification products daily) → **5** (dashboard) → **6** (Q5 screening) as a separate,
output-preserving change.
