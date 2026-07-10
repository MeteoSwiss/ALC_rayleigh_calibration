# Phase B — read-once / share-many (detailed design)

*Load each instrument-day's L1 (once) and CAMS (once, closest cell), coarsen to a shared
working grid, compute the WV transmission once, and fan everything out to classification +
Rayleigh + cloud + housekeeping + OmB + sensitivity. Companion to
`ceiloclass_integration_plan.md`.*

---

## 0. The two efficiency checks (answered)

| Question | Verdict | Evidence |
|:--|:--|:--|
| **Compute WV two-way transmission once?** | **Yes** — done 3× today, not cached | `rayleigh/calibration.py:531`, `cloud/calibration.py:877`(→1807), `omb/omb.py:321` all call it independently. Pure fn of (CAMS WV profile, λ, range grid, CAMS times) → compute **once per 910 nm day** on the shared grid, share. CHM15k (1064 nm) skips WV. |
| **Lazy-load L1?** | **No — not the lever.** Read-once is. | L1 opens with `netCDF4.Dataset` (`l1_window.py:33`, `cloud:396`); every consumer needs the **full** `rcs_0`. Nothing to skip lazily. Keep netCDF4 (fast, no xarray/pandas). |
| **Lazy-load CAMS?** | **Yes — big lever** | 0.4° monthly file is **8.5 GB**, chunked. Closest-cell lazy read = **2.7 s** vs full-domain **41 s** (**15×**). Read the one cell once, share. |

---

## 1. Architecture: one read, in-memory views

```
load_instrument_day(stream, dates)
├─ L1  : netCDF4, read the day/night span ONCE → native arrays
├─ CAMS: xarray .sel(nearest).values ONCE → single-cell columns (t,q,z,lnsp,aer532/1064)
└─ build InstrumentDayData:
     • native rcs/cbh/HK/range/time            (for sensitivity & OmB noise fidelity)
     • working-grid view = coarsen(native)      (≥10 m / ≥15 s; Rayleigh, cloud, classification, HK)
     • CamsProfile: heights (a/b+lnsp), T, P, n_wv, aerosol β   (closest cell)
     • wv_t2(z, cams_time): two-way WV transmission, computed ONCE on the working grid (910 nm)
```

Consumers **never re-open a file** — they select the view they need. The expensive things
(L1 full read, CAMS 8.5 GB open, WV HITRAN convolution) happen exactly once per instrument-day.

---

## 2. Resolution normalization — dynamic, simple, shared

**Config-agnostic** (resolution varies by station and over time, so never hardcode per type):
derive the factor from *each file's own* grid.

```python
def coarsen_factor(native, target, tol=0.9):
    # 1 unless the file is meaningfully finer than the target; then nearest whole factor.
    return 1 if native >= tol * target else max(1, round(target / native))

def block_average(a, ft, fr):                 # simple mean of consecutive bins, no interpolation
    nt, nr = a.shape
    a = a[:(nt // ft) * ft, :(nr // fr) * fr]
    return ma.mean(a.reshape(nt // ft, ft, nr // fr, fr), axis=(1, 3))
```

Only the **config target** (default 10 m / 15 s) is a knob; the factor is per-file. Measured:

| Type | native | factor | working grid |
|:--|:--|:--|:--|
| CHM15k | 15 m / 15 s | r÷1 t÷1 | unchanged |
| CL31 / CL51 | 10 m / ~30 s | ÷1 ÷1 | unchanged |
| **CL61** | **4.8 m** / 30 s | **r÷2** t÷1 | 9.6 m (**50 % fewer cells**) |

Range averaging preserves the layers (9.6 m still resolves thin liquid); time stays ≥15 s.
Coordinates (`range`, `time`) are block-averaged too; per-profile fields (`cbh`, HK) are
time-averaged only (median of finite for CBH — it's a height, not a signal).

**Two guards:**
- **No double-averaging.** Cloud currently hardcodes `average_time_s=30, average_range_m=10`
  (`_do_cloud`, ~line 307). Replace with the shared working grid; cloud may *further* time-average
  to 30 s **on top** of the working grid (cheap, on the already-small array) if it still wants 30 s.
- **Sensitivity noise fidelity.** Block-averaging N gates drops per-gate noise by √N, so feeding
  the *coarsened* stack to the **sensitivity** step would report a better detection limit than the
  instrument achieves natively. Sensitivity (and OmB's SNR screen) therefore read the **native
  view** — which is free, it's already in memory from the single read. Everything else uses the
  working grid. This honours "read once, one config" without misstating the native noise.

---

## 3. CAMS — closest cell, lazy, once; WV transmission once; correct ground

- **One lazy cell read:** `ds = xr.open_dataset(f); sub = ds.sel(latitude=lat, longitude=lon,
  method="nearest"); {v: sub[v].values ...}` for `t,q,z,lnsp,aerbackscatgnd532,aerbackscatgnd1064`.
  2.7 s vs 41 s. Use `ALC_CAMS_DIR` (0.4° `CAMS_Monthly_04`), 1° fallback.
- **Heights:** `z`/`lnsp` are *surface* fields → per-level pressure from L137 a/b + `exp(lnsp)`,
  hydrostatic height (reuse `water_vapor._cams_levels`).
- **Ground extrapolation (temperature):** below the lowest CAMS level use the **standard ISA lapse
  (6.5 K/km) anchored at the lowest level**, *not* a two-point slope (the ~20 m-spaced, often
  surface-inverted lowest levels give a wild slope → −110 °C at the surface). With the 0.4° cell the
  extrapolation span is only ~300 m (station 490 m, cell surface ~789 m) and T then matches
  Cloudnet ≤ 1 °C from 1 km up. `n_wv` keeps the existing constant-fill below surface.
- **WV transmission once (910 nm):** compute `two_way_wv_transmission` on the **working grid** at the
  CAMS time steps a single time; store in `InstrumentDayData`. Rayleigh, cloud and OmB each
  time-interpolate it to their grid. This only works because they now share one range grid.
- **Cross-stream (future):** ~400 European streams share the *same* box file. A pre-pass that
  extracts all station columns in one open (2.7 s × 400 → one pass) is a further win beyond Phase B's
  per-stream read-once — note it, don't block on it.

---

## 4. Concrete refactor steps (sequenced, output-preserving, flag-gated)

All behind `ALC_READONCE=1`; each step keeps the old path until its diff check passes.

| # | Change | Touches |
|:--|:--|:--|
| **B1** | New `InstrumentDayData` + `load_instrument_day()` — L1 once (netCDF4), CAMS once (lazy cell). No behaviour change yet. | new `calibration/io/instrument_day.py` |
| **B2** | Single CAMS cell-reader; make WV transmission a **compute-once** shared field (pass into Rayleigh/cloud/OmB instead of each recomputing). | `water_vapor.py`, `cloud/calibration.py:877`, `omb.py:321`, `rayleigh/calibration.py:531` |
| **B3** | Coarsening → **working grid** producer; wire Rayleigh + cloud + classification to it; drop cloud's internal 30 s/10 m double-average. | `io/instrument_day.py`, `_do_cloud`, `_do_rayleigh` |
| **B4** | Make Rayleigh window search **metre-based** (`increment_bins=8` is a gate count → convert to metres) so coarsening can't change search granularity. | `rayleigh/rayleigh_fit.py:104` |
| **B5** | Thread the shared object into `_do_rayleigh/_do_cloud/_do_monitoring/_do_omb/_do_sens`; each stops re-opening. Sensitivity/OmB take the **native view**. | `scripts/run_network_calibration.py` |
| **B6** | Add `_do_classification(shared)` consumer (β = rcs×factor on the working grid + CAMS T). | `scripts/run_network_calibration.py`, `calibration/classify/…` |

Because all steps already run in one subprocess per stream (`_process_stream`), this is threading
one object through — **not** a process re-architecture.

---

## 5. Validation

- **Output-preserving diffs.** Old vs new on a sample of streams/days; the calibration/OmB/sensitivity
  CSVs must match. The *only* intended numeric change is CL61 range 4.8 → 9.6 m for Rayleigh/cloud;
  quantify ΔC_L on ~5 CL61 nights (expect ≪ 1 %). CHM15k/CL31/CL51 have factor 1 → bit-identical.
  Sensitivity/OmB use the native view → unchanged.
- **Perf.** Per-stream wall-clock before/after. Expect L1 3–5×→1×, CAMS 2–3×→1× (and 41 s→2.7 s if any
  path was reading full fields), WV 3×→1×.

## 6. Effort & risks

- **~1–1.5 weeks.** Risks & mitigations: shared-load regression → flag + diff harness; the
  `increment_bins` metre conversion (small, isolated); sensitivity/OmB must keep the native view
  (explicit); memory holds native + working view (CL61 ~0.5 GB/worker — fine at 6 workers).
- **Order:** B1→B2 give the read-once + WV-once wins with the least surface area; B3–B4 add the shared
  grid; B5–B6 complete the fan-out and add classification. Ship B1–B2 first if a quick I/O win is wanted.
