# v2.2 network release → interactive dashboard → EWC publish

Plan only, written 2026-08-17. Nothing below has been executed. It covers the three things asked
for: recompute every site on balfrin, rebuild the dashboard **with the interactive daily panel**,
and publish to EWC — plus the one decision that has to be taken before any of it starts.

---

## 0. Where we already are

| piece | state |
|---|---|
| v2.2 calibration, 434 streams, 2025-01-01 → 2026-08-13 | **done** on balfrin, 3 h05, `/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel` |
| OmB + sensitivity add-ons (`--no-cal --sens --omb`) | **running**, job 5123654, 354/866 CSVs, 0 failures, ~3 h left |
| acceptance gates G1–G5 | script written; G2/G4 passed on the earlier run, G1 re-run needed, G3 blocked on the add-ons |
| interactive daily panel + tests | **done and committed** (`3de57ca` + follow-up), 60 tests green |
| per-day payloads for the network | **not started** — this is the new, and only expensive, step |

The calibration itself is therefore *already recomputed*. What is missing for an interactive
dashboard is the per-day payload store, and that is where the whole plan hinges.

---

## Scope and resolution: full period, full resolution

**Settled: the full time series, 2025-01-01 → 2026-08-13, at the resolution the panel already
produces.** No averaging, no shortened window.

The sizing question that looked alarming in isolation answers itself once the payload is compared
to the thing it REPLACES rather than to zero. Measured on 2 874 real diagnostic PNGs from a
dashboard build, and on the 25-station payload build:

| store, for the same 234 794 (day, method) | size |
|---|---|
| diagnostic PNGs — today's product | median **1 190 KB** each → **266 GB** |
| interactive payloads, full resolution | median **357 KB** each → **80 GB** |
| same, served gzipped | → **53 GB** |

The interactive panel is **3.3× smaller than the PNG store at full resolution, 5× gzipped**, and it
is interactive. Reducing the curtain grid would have been optimising the wrong side of a 3× win.

For the record, since the analysis is done: base64 size is exactly `ceil(nt·nz/3)·4`, the curtain is
720 × ~332 (the `_curtain` target is 800 × 800, so `st_t = st_r = 1` for a normal night), and ~46 KB
of each payload is diagnostics (28 KB) + profile (16.5 KB) rather than curtain. Those are the knobs
**if** a future constraint ever demands them. Nothing here needs them.

### What follows from keeping full resolution

1. **Serve gzipped** (`Content-Encoding: gzip` on upload, or `gzip_static` on nginx). It is free,
   costs one flag, and takes 80 GB to 53 GB.
2. **Retire the diagnostic PNGs for every day that has a payload.** This is the point: the panel
   renders the same night from data, so keeping both means paying 266 GB to duplicate what the
   payload already shows. Net storage change for the release is then **−213 GB**, not +80 GB.
   `ops/publish.sh` already prunes local PNGs after upload; the change is to stop *producing and
   uploading* them for covered days, and to drop the `section.diag` viewer from pages whose panel
   covers the same period (the panel already owns the keyboard when a viewer is absent).
3. **Transfer becomes the long pole, not storage.** 235 000 objects / 80 GB across
   balfrin → zueub434 → EWC. Plan it as hours, run it per-station with `--partial --inplace`, and
   keep it inside the cron's `flock` (§6).

**One thing to confirm, because it is not mine to assume:** the EWC bucket has to hold ~53 GB and
235 000 objects. That is well below the PNG store it replaces, so the quota is very likely already
adequate — but it is worth one `aws s3 ls --summarize` on the existing bucket before the fleet run,
since the two stores coexist until the PNG retirement completes.

---

## OPEN DEFECT found during execution: the sensitivity product is not being produced

Measured on the release tree while the add-ons job ran:

| product | streams |
|---|---|
| `_omb_cache.npz` | **429** / 434 |
| `_omb.csv` | 387 |
| `_sens_cache.npz` | **15** |
| `_sens.csv` | 15 |

The live bucket carries `_sens.png` for **431** stations, so this is a REGRESSION of this run, not a
pre-existing gap. What has been ruled out, with evidence:

* **not a crash** — zero `failed:` lines in the job log;
* **not the Kalman gate** (`if not kmap: return None`) — the failing streams have 556 and 463
  rayleigh Kalman rows, *more* than a stream that succeeds (41);
* **not the regression guard** — it fires only when a non-empty `_sens.csv` exists without its
  cache, and only 15 `_sens.csv` exist at all;
* **not the kernel** — driven directly on a real day for a failing stream,
  `sensitivity_over_period` returns a `SensResult` with `dates=1`;
* **not OmB consuming the shared per-day read first** — all 15 streams with sensitivity also have
  OmB, so the two are not mutually exclusive.

What is left, and what the evidence points at: the streams that succeed are the SHORT ones
(0-20000-0-00202_A has 52 days of L1), while 500-plus-day streams produce nothing — so the suspect
is the day loop in `_do_sens` over a long window (it accumulates one `SensResult` per day in
`parts` before writing, unlike OmB which updates its cache per day).

**Impact and decision.** It blocks nothing in this release: the payloads, the constants, the Kalman,
OmB, classification and every time series are complete. It costs one of the three cards at the
bottom of the station page. The release therefore proceeds, and:

* **OmB images ARE re-uploaded** (complete and v2.2);
* **sensitivity images are NOT** — the bucket keeps the operational v2.0 ones rather than a
  15-station v2.2 patchwork;
* the fix is a targeted `--sens`-only re-run once the day loop is understood, which touches no
  calibration value.

---

## 1. Prerequisites (must be true before step 2 starts)

1. **Add-ons job finished** — `ALC_V22_ADDONS_DONE rc=0`, `count(_omb_cache.npz) == count(*_omb.csv)`.
   The `.npz` caches are not optional: the daily runner's regression guard compares cache to CSV
   and *skips* rather than rebuilds, so an archive transferred without them freezes OmB and
   sensitivity permanently.
2. **Acceptance gates G1–G5 re-run** on the finished tree, in particular
   * G1 no mixed archive (every rayleigh row version 220, no duplicate (key, method, date));
   * G3 completeness — `_cal/_kalman/_hk/_status/_sens/_omb` + annual NetCDF + `classification/`
     + **both** `.npz` for every census stream.
3. **`classification/` copied** from the old tree — it is method-independent, and regenerating it
   is not part of this release. Forgetting the copy removes the gallery from the whole network.
4. **Code on balfrin matches the commit being released** — `git -C ~/alc_v22_code fetch && reset
   --hard <sha>`; the payload generator and the panel must be the versions the tests ran against.

## 2. Payload generation on balfrin (new)

New `ops/cscs/alc_v22_payloads.sbatch`, modelled on the two existing ones:

* `--partition=postproc` (**CPU only** — the partition assertion at the top of the existing scripts
  is copied verbatim), 1 node, `--cpus-per-task=256`, `--mem=0`, `--chdir=/scratch/mch/mhrvo`,
  no `--account`.
* Same environment as the release job: `ALC_MOLECULAR_METHOD=eprof_v2.2`,
  `ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'`, `ALC_FULLCAL_DIR=.../E_PROFILE_calout_v22_rel`,
  `STREAM_TIMEOUT=28800`. **These are not optional under a payload run**: the payload generator
  *re-runs the retrieval in-process* to capture the curtain and the diagnostics, so a missing
  switch produces v2.0 panels sitting under v2.2 numbers — a mixed-vintage page.
* `scripts/build_station_dashboard.py --key ... --start 20250101 --end <D_END>
  --cal-dir $ALC_FULLCAL_DIR --l1-root ... --out /scratch/mch/mhrvo/dash_payloads
  --payloads --no-pages --workers 200`
* Stream list = **union** of the live census and the directories present in the tree, same rule as
  the release run, so a station commissioned since the census snapshot is not dropped.
* Per-PID work dirs are already implemented; concurrent writers to a shared per-ident NetCDF
  corrupt it, which is why this must not be run any other way.

**Before the full run**: one 3-stream × 7-day canary to confirm wall-time per night and mean
payload size, then extrapolate. Payload generation re-runs the retrieval, so expect the same order
as the release recompute itself: **≈ 3–4 h** for 590 days × 434 streams on one 256-core node.

No curtain change is needed — full resolution is the decision (see above). Scratch must have room
for ~80 GB before the run starts.

## 3. Build the site

On balfrin (or on zueub434 after transfer — it needs only the CSV tree + the payload dir):

```
python scripts/build_dashboard.py --fullcal <v22 tree> --out <staging site> --workers 32
```

Full build, **no `--changed-only`**, into a directory that is not the live site.

Then the gates that exist for this:

* `ALC_SITE_DIR=<staging> python -m pytest tests/test_dashboard_site.py -q` — 40 structural checks;
* the browser tier — `pytest tests/test_dashboard_browser.py` — needs Playwright + Chromium, which
  is unlikely to be installable on the server. **Run it on the workstation against a synced
  sample** (10–20 stations is enough; the tier picks representative pages by feature). It is the
  only thing that proves the interactive panel actually works, so it must run somewhere.
* page count ≥ 434 and zero `REGRESSION-GUARD` lines in the build log.

## 4. Transfer

balfrin → zueub434, over the existing route:

* the calibration tree **including both `.npz` per stream** (see §1.1);
* the payload store `data/<key>/*.json[.gz]` — ~235 000 objects, ~80 GB. rsync one directory per station
  rather than one flat sweep, and `--partial --inplace`; a stalled transfer then resumes per
  station instead of restarting.

Nothing is deleted at this stage. The live site stays exactly as it is.

## 5. Publish to EWC

`ops/publish.sh` already splits the site: images → S3 bucket `eprofile-alc-dashboard`,
HTML → web VM `hem@136.156.139.31:/var/www/alc`. The payloads follow the **images** path, not the
HTML path — `dailypanel.js` already falls back to the bucket when a payload is not same-origin, so
this needs no client change.

Additions to `publish.sh`:

1. `aws s3 sync <site>/data s3://eprofile-alc-dashboard/data --size-only` with
   `--content-encoding gzip` if the gzip decision is taken. Parallelism matters at 235 000 objects;
   set `max_concurrent_requests` in the AWS config rather than looping.
2. A **prune** rule mirroring the existing diag-PNG prune: payloads older than the window are
   removed from the bucket after a successful sync, if a horizon is adopted (§0 decision 2);
   otherwise the store grows by ~40 MB/day, ~15 GB/year.
3. HTML rsync unchanged — and remember it intermittently returns **rc=2** while the HTML still
   lands; re-run `bash ops/publish.sh` rather than debugging it.

Verify live: open three stations of different type (CHM15k / CL31 / CL61) and check the console
shows `[ALC] rangesync ready`, the panel draws, and the period selector reports
`unchanged: none`. That console tracing was added for exactly this.

## 6. Protecting it from the daily cron

This is the part that interacts with tomorrow morning's server change.

* The 15:00 cron runs `ops/run_daily.sh` → a 434-stream calibration + `--changed-only` build +
  publish. If it fires mid-transfer, two 433-stream runs collide.
  **Take the same `flock` the cron uses for the whole transfer + build window**, or disable the
  cron entry for the duration. Do not rely on timing.
* The cutover itself is **one edit of `ops/config.sh`** (method + spectrum + `ALC_FULLCAL_DIR`),
  never a partial one: with the method flipped but the tree unchanged, the cron merges v2.2 rows
  into the v2.0 archive and `_preserve_existing_rows` keeps both per (method, window) —
  irreversibly.
* The daily flow must gain a **payload step** (D-1 only, ~456 payloads) and the matching prune,
  otherwise the interactive panel freezes on the release date while every other chart advances —
  the most confusing possible failure mode.
* **Back up the live census first.** `ALC_CENSUS` points *inside* the checkout and is rewritten
  daily; a `git checkout` on the deploy silently reverts every station commissioned since.

## 7. Rollback

Restore `ops/config.sh.pre_v22_*`, full rebuild, publish. Never a partial rebuild with
`rsync --delete`. The payload store can stay in the bucket — an old HTML simply does not reference
it, and the prune will retire it.

---

## Order of work, and what needs you

1. *(automatic)* add-ons job finishes — ~3 h.
2. Gates G1–G5 on the finished tree — ~30 min.
3. Payload canary (3 streams × 7 days), then the fleet run on balfrin — **~3–4 h**, full period,
   full resolution.
4. Staging build + site tier on the server, browser tier on the workstation — ~30 min.
5. Transfer + publish, under the cron's flock — **the long pole**: 235 000 objects / 80 GB.
6. Retire the diagnostic PNGs for covered days (the −213 GB), once the panel is verified live.

Steps 3–6 can all happen before tomorrow morning; only step 7 touches the operational server, which
is where your change comes in.
