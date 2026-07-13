# 11 — Efficient batch reprocessing: chunk-outer read-once design

**Status:** M1 + M2 + M4 implemented on branch `batch-readonce` (2026-07-13); M3 deferred ·
**Scope:** `scripts/run_network_calibration.py` batch path ·
**Goal:** reprocess ~425 streams × 1.5 years reading the L1 archive ~once, on any of the three
hosts (home Windows PC, zueub434, CSCS balfrin), without giving up the measured advantage of
**multi-day batched reads** (a batched read of several days is much faster than day-by-day).

## 1. Problem — measured baseline

The 2026-07-13 home reprocess (24 workers, L1 on a SATA QLC disk) stalled: disk pinned at 100 %
serving ~55 MB/s of random reads, CPU ~7 %, and after 1.5 h only **1 of 48 started streams had a
`_cal.csv`** (47 had classification output only). The identical code ran perfectly as an isolated
CSCS task (Payerne A, June, all products: 5 min 51 s, correct June 20–21 rejection). The failure is
not the science code — it is **I/O amplification** in the batch driver.

Measured anchors:

| Anchor | Value |
|---|---|
| Home, 1 stream, 16 d, Rayleigh+plots | 20 s (~1.25 s/day incl. reads) |
| Home, 1 stream, 545 d, Rayleigh-only | 684 s |
| CSCS postproc, 1 stream, 30 d, ALL products | 5 min 51 s (~12 s/day, 4 EPYC cores) |
| CSCS 2026-06 array, per-stream mean | cal 2013 s, sens 1630 s (max 133 / 182 min) |
| Home D: (Samsung 870 QVO SATA) | ~55 MB/s effective under 24-worker random load; ~530 MB/s sequential cap |
| L1 archive 2025-01 → 2026-07, ~425 streams | ≈ 1.2–1.5 TB (sample-based estimate) |

## 2. Root cause — the batch driver multiplies reads ×10

Two multipliers stack in `_process_stream`:

1. **Night pairing ×2.** `build_file_paths` (data_loader.py) returns *previous + current day*
   for every night, so within one pass each daily L1 file is read twice (as day D of night D and
   as day D−1 of night D+1). The per-night reader's cache is 3 days — sized for the daily cron,
   useless across a 545-day sweep.
2. **Pass-outer ×5.** The batch path used to run each product over the whole window sequentially —
   classify → Rayleigh → cloud → OmB → sens — so by the time pass N+1 revisited day 1, the cache
   was evicted ~544 days earlier. Every pass re-read the whole period.

Net: **~10× the archive bytes per stream**, plus a monitoring pass re-opening every file for the
1-D housekeeping vars. On Lustre (CSCS) this is absorbed; on a QLC SATA disk (home) or a slow NAS
(zueub434 backfills) it is the bottleneck — ~13 TB of reads against a ~1.3 TB archive.

## 3. Design — stream × chunk, read once, calibrate late (implemented)

The parallelism unit stays the **stream** (Windows pool worker ≡ zueub434 `--workers N` ≡ CSCS
array task). Inside a stream the loop is inverted **per chunk** — not per day, so the measured
multi-day batched-read advantage is preserved:

```
reader = _make_chunk_reader(s, start, end)             # ALC_CHUNK_DAYS per batched read (+1 head day)
for (cs, ce) in _chunks(start, end, CHUNK_DAYS):       # sweep 1: calibration trio, chunk-outer
    classify(cs..ce, reader)                           #   curtain .nc per day (contam for Rayleigh)
    rayleigh(cs..ce, reader)                           #   nights = slice_night(D) views of the chunk
    cloud(cs..ce, reader)                              #   days   = slice_to_date(D) views (memoised)
kalman <- cal rows                                     # needs the WHOLE window -> after sweep 1
_do_omb_sens(start, end, kalman, fresh reader)         # sweep 2: OmB + sens INTERLEAVED, one day
                                                       #   sweep sharing each day's coarse dict
```

Mechanics (all in `calibration/io/instrument_day.py` + `scripts/run_network_calibration.py`):

- `load_instrument_day` already accepts a **file list** → the chunk read is the existing batched
  reader over `[chunk_start − 1 … chunk_end]` (head day = the first night's D−1; ~3 % overlap).
- **`slice_night(D)`** (new) serves a night as a row-subset **copy** of the chunk — exactly what
  `load_instrument_day([D−1, D])` returns today (`working` aliased to native, CAMS/WV dropped), so
  consumers cannot tell the difference and cannot corrupt the chunk. The **day-slice memo is
  shared** with the chunk, so each day is re-coarsened once per chunk (better than today's
  once per pass).
- **Parity fallbacks:** a night with no file for D−1 nor D returns `None` (as today); a failed
  chunk concatenation (e.g. mid-chunk range-grid change) degrades that chunk to per-night reads.
- **Classification order:** curtains are written chunk-by-chunk *before* that chunk's Rayleigh, so
  the flag −11 screen finds D−1's curtain on disk (previous chunk or this one) — the read-back is
  a few KB/day and stays as-is.
- **The Kalman dependency:** OmB and sensitivity consume the full-window Kalman C_L, which only
  exists after all nights are calibrated (the filter's measurement-noise/outlier statistics are
  full-series, so it cannot be finalised chunk-incrementally). They therefore run as **one second
  interleaved sweep** (`_do_omb_sens`): both products share each day's
  `slice_to_date(D).to_omb_dict()` and the same chunk reader — one batched read per chunk, half
  the sweeps of the old back-to-back pair, per-day cache updates and aggregation verbatim.

`ALC_CHUNK_DAYS` (default **7**) sets the batch size: a 1-day window (the daily cron) is a single
1-day chunk = exactly today's `[D−1, D]` pair read once — **operations are unaffected by
construction**. Reprocessing sets 31 (14 for CL61 if memory-tight: CL61 native ≈ 6 GB/31 d,
CHM15k ≈ 0.7 GB).

## 4. Milestones

| # | What | State |
|---|---|---|
| M0 | `validation/_run_with_iocount.py` — wall time + OS read/write bytes around a runner call | done |
| M1 | Chunk loop in `_process_stream` + `slice_night()` + `_make_chunk_reader` (parity fallbacks) | done |
| M2 | `_do_omb_sens`: OmB + sens as ONE interleaved chunked sweep (shared day dict) | done |
| M3 | C-agnostic OmB/sens caches (apply C_L at aggregation) → drops sweep 2 entirely (≈1.03×); schema bump ⇒ one-off network cache rebuild | **deferred** — do not land the same week the reprocessed caches ship |
| M4 | `tests/test_chunk_readonce.py` (slice_night/boundary parity — 3/3 pass) + `validation/diff_fullcal.py` golden diff on 06610 A/B/C × 2 months | done — see §7 |

## 5. Per-machine behaviour

| Host | Setup | Expected |
|---|---|---|
| **Home PC** (i9-14900KF, 128 GB; L1 on QLC SATA; CAMS moved to NVMe) | 10–16 workers, `ALC_CHUNK_DAYS=31` | reads ~sequential, ~2 sweeps ≈ 2.6 TB ⇒ ~2–3 h disk, run CPU-bound **~10–18 h overnight** (vs infeasible today) |
| **zueub434** (slow NAS, no sudo) | cron unchanged (1-day window ≡ today); backfills `ALC_CHUNK_DAYS=31` | reprocess/backfill ≈ **×4–5 faster** |
| **CSCS balfrin** (postproc array) | sbatch + `export ALC_CHUNK_DAYS=31` | **−30–40 % core-hours** (decode-once); wall 5–8 h → ~3–5 h; big CL61 tail streams shrink most |

## 6. Estimated gains (per stream, S = its share of the archive)

| Stage | L1 bytes read | vs today |
|---|---|---|
| Today (pass-outer) | ≈ 10 S (+ HK opens) | — |
| M1 (chunked calibration trio; OmB+sens as before) | ≈ 1 S + 4 S = 5 S | ×2 |
| **M1+M2 (this branch: + interleaved OmB/sens sweep)** | ≈ **2.06 S** | **×5** |
| M3 (C-agnostic caches, deferred) | ≈ 1.03 S | ×10 |

Compute drops too (netCDF decode/coarsen that happened 5× happens ~2×; day re-coarsening now once
per chunk instead of once per pass).

## 7. Golden validation (M4)

Setup: Payerne 0-20000-0-06610 **A** (CHM15k, Rayleigh+classify, June aerosol-rejection
regression), **B** (CL31, cloud), **C** (CL61, both + depol classify, stream starts mid-window) ·
window **2026-05-01 → 2026-06-30** · `--methods rayleigh,cloud --sens --omb --classify`,
`PLOTS=0` · baseline = `development` @ fc2b512 (pass-outer) vs this branch with
`ALC_CHUNK_DAYS=31`. Compared: `_cal/_kalman/_hk/_status/_sens/_omb` CSVs (float tol 1e-9),
yearly `ALC_calibration_*.nc` arrays, classification NetCDF sets + sampled
`target_classification` equality (`validation/diff_fullcal.py`).

Results (2026-07-13, home i9-14900KF, single worker, warm page cache):

| Stream | Baseline (pass-outer) | Chunked (31 d) | Wall |
|---|---|---|---|
| A · CHM15k (ray+classify) | 118.6 s | 83.6 s | −30 % |
| B · CL31 (cloud) | 117.6 s | 121.6 s | ≈ |
| C · CL61 (ray+cloud+depol classify) | 324.2 s | 216.7 s | −33 % |
| **total** | **560.4 s** | **421.9 s** | **−25 %** |

Correctness: **golden diff 3/3 IDENTICAL.** Bit-identical across A/B/C: `_cal.csv`, `_kalman.csv`,
`_hk.csv`, `_status.csv`, the 51 classification NetCDFs (sampled `target_classification` equal),
and the yearly `ALC_calibration_*.nc` (identical on time / lidar_constant / calibration_method once
compared as a (method, time)-sorted set — chunked writes rows in append-order interleaved by chunk,
the daily cron already appends per-day, so row order was never a semantic invariant; the golden diff
sorts before comparing). **OmB** compared identical where produced (B, `_omb.csv` 645 B).
`test_chunk_readonce.py` 3/3 pass.

Coverage honesty — **`_sens.csv` was NOT compared** in this window:
- A/C got **0 successful rayleigh nights** in this aerosol-heavy May–June window → empty rayleigh-Kalman
  → no rayleigh-method sens/omb by design (A/C use the rayleigh method).
- B (CL31, cloud method) has a 40-row Kalman and *should* produce sens, but its `_sens_cache.npz`
  write hit the Windows `os.replace` flake repeatedly **on the baseline** (main has no retry), so
  there was no baseline sens to diff against.
Sens is therefore covered here by (i) **OmB parity** — sens rides the identical `_do_omb_sens`
day-sweep + incremental-cache pattern as OmB, which *is* identical; (ii) the read-path unit tests.
A headline-sens golden diff needs a clear-night window (≥5 rayleigh successes) and the retry fix on
both checkouts — worth a follow-up run before the branch merges.

Caveat on the numbers: this ran with the files warm in the Windows page cache (the baseline run
touched them minutes earlier) and the M0 `read_bytes` probe reported ~0 GB on Windows — so the
**I/O reduction is NOT what this benchmark measured**; even the −25 % wall is mostly the fewer
decode/coarsen passes, and the real ×5 read win only appears on a **cold / disk-bound** run (the
actual reprocess). The value proven here is *correctness parity* + a compute win; the I/O win is
by construction (§6) and will be re-measured on the first cold CSCS/home batch.

A Windows-only robustness fix rode along: `_savez_atomic` (incremental.py) now retries the cache
`os.replace` (WinError 5/32 when AV/indexer scans the fresh `.tmp.npz`) — without it a network
reprocess drops the occasional station's OmB/sensitivity product (hit twice in this very run).

## 8. Risks / notes

- **Chunk-boundary night**: +1 head day per chunk (~3 % extra reads) — verified by
  `test_chunk_reader_serves_nights_and_boundaries`.
- **Missing days**: per-day file checks unchanged; a night with no files stays `None`.
- **Mid-chunk grid change**: chunk concatenation failure degrades that chunk to per-night reads
  (today's exact behaviour), never fails the stream.
- **Chunk metadata**: `rcs_units`/`l0_wavelength` are read from the chunk's first file (was: each
  night's D−1 file). Only matters across a mid-chunk instrument swap — negligible and self-heals
  at the next chunk.
- **Daily cron**: a 1-day window is one 1-day chunk — same two files, read once, as today.
