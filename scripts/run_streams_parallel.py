# -*- coding: utf-8 -*-
"""Run a FEW streams over a LONG window on many cores, by chunking the date range.

Why this exists
---------------
``run_network_calibration.py`` parallelises over STREAMS: one single-threaded process each. That is
the right unit for a network run (hundreds of streams), but it caps a small-site diagnostic run at
one core per stream -- 7 streams on a 32-core box is a 22 % ceiling, and the real figure is lower
still because ~40 % of the wall time is I/O (L1 read, PNG write) rather than CPU.

This launcher splits each stream's window into date CHUNKS and runs (stream x chunk) as independent
processes, each with its OWN output directory -- separate directories are mandatory, not tidiness:
the per-stream CSV and the yearly NetCDF are written without locking, so two processes sharing an
output tree corrupt each other (observed: successful calibrations silently lost as EXC OSError).
The per-chunk outputs are then merged back into one directory per stream.

  7 streams x 4 chunks = 28 processes -> ~3.5x faster than 7 processes on this machine.

What the merge produces (identical in content to a sequential run):
  <key>_cal.csv      all rows, deduped on (method, date), sorted
  <key>_kalman.csv   RECOMPUTED from the merged rows with the runner's own _kalman_rows()
  <key>_hk.csv / _sens.csv / _status.csv   concatenated, deduped on the first column
  plots/**           hardlinked into one tree (no copy, no extra disk)

NOT merged: the yearly E-PROFILE NetCDF. Each chunk writes its own under parts_nc/<chunk>/ because
merging the exchange format needs the writer's own record semantics. A chunked run is therefore for
DIAGNOSTICS and dashboards (which read the CSVs); use the sequential runner when the NetCDF is the
product.

Example
-------
  python scripts/run_streams_parallel.py \\
      --streams 0-20000-0-06610_A,0-20000-0-06610_C,0-20000-0-10393_0 \\
      --start 20250101 --end 20260814 --chunks 4 --methods rayleigh,cloud --plots --sens \\
      --out C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v20
"""
from __future__ import annotations
import argparse
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def _chunks(start: str, end: str, n: int):
    """Split [start, end] into n contiguous YYYYMMDD sub-windows (last one absorbs the remainder)."""
    d0 = datetime.strptime(start, "%Y%m%d")
    d1 = datetime.strptime(end, "%Y%m%d")
    total = (d1 - d0).days + 1
    n = max(1, min(int(n), total))
    step = total // n
    out = []
    for k in range(n):
        a = d0 + timedelta(days=k * step)
        b = d1 if k == n - 1 else d0 + timedelta(days=(k + 1) * step - 1)
        out.append((a.strftime("%Y%m%d"), b.strftime("%Y%m%d")))
    return out


def _run_one(job, args, parts_dir):
    key, (c0, c1), idx = job
    outdir = parts_dir / f"{key}__{idx}"
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["ALC_FULLCAL_DIR"] = str(outdir)
    env["PLOTS"] = "1" if args.plots else "0"
    if args.cams:
        env["ALC_CAMS_DIR"] = args.cams
    cmd = [sys.executable, str(REPO / "scripts" / "run_network_calibration.py"),
           "--stream", key, "--start", c0, "--end", c1, "--methods", args.methods, "--force"]
    if args.sens:
        cmd.append("--sens")
    if args.omb:
        cmd.append("--omb")
    log = parts_dir / f"{key}__{idx}.log"
    with open(log, "w", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT)
    return key, idx, c0, c1, rc


def _concat_csv(paths, out_path, dedupe_on):
    """Concatenate CSVs sharing a header; keep the FIRST row per dedupe key; sort by that key."""
    rows, header = {}, None
    for p in paths:
        if not p.exists():
            continue
        with open(p, newline="", encoding="utf-8") as fh:
            rd = csv.DictReader(fh)
            if rd.fieldnames:
                header = header or rd.fieldnames
                for r in rd:
                    rows.setdefault(tuple(r.get(k, "") for k in dedupe_on), r)
    if not header or not rows:
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])
    return len(rows)


def _link_tree(src: Path, dst: Path):
    """Hardlink every file of src into dst (same volume); fall back to copy across volumes."""
    n = 0
    for f in src.rglob("*"):
        if not f.is_file():
            continue
        t = dst / f.relative_to(src)
        t.parent.mkdir(parents=True, exist_ok=True)
        if t.exists():
            continue
        try:
            os.link(f, t)
        except OSError:
            shutil.copy2(f, t)
        n += 1
    return n


def merge(key, parts_dir, out_dir, n_chunks):
    """Merge every chunk of one stream into <out_dir>/<key>/, as a sequential run would have left it."""
    from run_network_calibration import _kalman_rows  # the runner's own Kalman, unchanged

    dirs = [parts_dir / f"{key}__{i}" / key for i in range(n_chunks)]
    sdir = out_dir / key
    sdir.mkdir(parents=True, exist_ok=True)

    n_cal = _concat_csv([d / f"{key}_cal.csv" for d in dirs], sdir / f"{key}_cal.csv",
                        dedupe_on=("method", "date"))
    for name, keys in ((f"{key}_hk.csv", ("date",)), (f"{key}_sens.csv", ("date",)),
                       (f"{key}_status.csv", ("date",)), (f"{key}_omb.csv", ("date",))):
        _concat_csv([d / name for d in dirs], sdir / name, dedupe_on=keys)

    # Kalman must be recomputed on the MERGED rows: each chunk only saw its own slice of history.
    with open(sdir / f"{key}_cal.csv", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    krows = _kalman_rows(rows)
    with open(sdir / f"{key}_kalman.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["method", "date", "kalman", "kalman_std"])
        w.writeheader()
        w.writerows(krows)

    n_png = sum(_link_tree(d / "plots", sdir / "plots") for d in dirs if (d / "plots").is_dir())
    for i, d in enumerate(dirs):                      # keep the per-chunk NetCDFs, unmerged
        for y in d.glob("[0-9][0-9][0-9][0-9]"):
            if y.is_dir():
                _link_tree(y, sdir / "parts_nc" / f"chunk{i}" / y.name)
    return n_cal, len(krows), n_png


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--streams", required=True, help="comma list of '<wmo>_<ident>' keys")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--chunks", type=int, default=4, help="date chunks per stream (default 4)")
    ap.add_argument("--methods", default="rayleigh,cloud")
    ap.add_argument("--out", required=True, help="merged output dir (one subdir per stream)")
    ap.add_argument("--parts", default="", help="dir for the per-chunk outputs (default <out>/_parts)")
    ap.add_argument("--workers", type=int, default=0, help="parallel processes (default: cores-4)")
    ap.add_argument("--cams", default=os.environ.get("ALC_CAMS_DIR", ""), help="ALC_CAMS_DIR override")
    ap.add_argument("--plots", action="store_true", help="per-night diagnostic PNGs (PLOTS=1)")
    ap.add_argument("--sens", action="store_true")
    ap.add_argument("--omb", action="store_true")
    ap.add_argument("--merge-only", action="store_true", help="skip compute, just merge existing parts")
    args = ap.parse_args()

    keys = [s.strip() for s in args.streams.split(",") if s.strip()]
    out_dir = Path(args.out)
    parts_dir = Path(args.parts) if args.parts else out_dir / "_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    windows = _chunks(args.start, args.end, args.chunks)
    workers = args.workers or max(1, (os.cpu_count() or 8) - 4)

    if not args.merge_only:
        jobs = [(k, w, i) for k in keys for i, w in enumerate(windows)]
        print(f"{len(keys)} streams x {len(windows)} chunks = {len(jobs)} processes "
              f"over {workers} workers | {args.start}..{args.end} | plots={int(args.plots)}",
              flush=True)
        done = 0
        with ThreadPoolExecutor(max_workers=workers) as ex:      # threads: each job is a subprocess
            futs = [ex.submit(_run_one, j, args, parts_dir) for j in jobs]
            for fut in as_completed(futs):
                key, idx, c0, c1, rc = fut.result()
                done += 1
                print(f"  [{done}/{len(jobs)}] {key} chunk{idx} {c0}..{c1} rc={rc}", flush=True)

    print("merging ...", flush=True)
    for k in keys:
        n_cal, n_kal, n_png = merge(k, parts_dir, out_dir, len(windows))
        print(f"  {k}: {n_cal} cal rows, {n_kal} kalman rows, {n_png} files linked", flush=True)
    print(f"-> {out_dir}\nPARALLEL_RUN_DONE", flush=True)


if __name__ == "__main__":
    main()
