"""Extract the OPERATIONAL calibration constant from E-PROFILE L2 daily files.

Each L2 file (e.g. D:/E-PROFILE_L2_2026/<wigos>/2026/<MM>/L2_<key><YYYYMMDD>.nc) carries
`calibration_constant_0` (units m^3*sr*counts/s, dim time) -- the constant actually applied
operationally to produce the L2 attenuated backscatter. It is constant within a day, so we store
one value per (key, date). This is the same physical quantity as our calibrated C_L, so the ratio
C_L / operational is meaningful.

Writes a CSV (key,date,op_coeff) that the dashboard build joins in. Resumable: existing (key,date)
rows are kept and skipped. Restricted to the stations in our census by default.

Run:  python scripts/extract_l2_opcoeff.py <l2dir> <out_csv> [--start YYYYMMDD --end YYYYMMDD] [--all-keys]
"""
import argparse
import csv
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import netCDF4  # noqa: E402
import run_all_l1_2026 as R  # noqa: E402

_FNAME = re.compile(r"L2_(.+?)(\d{8})\.nc$", re.IGNORECASE)


def _read_one(fp):
    """Return (key, date, op_coeff) for one L2 file, or None."""
    m = _FNAME.search(os.path.basename(fp))
    if not m:
        return None
    key, date = m.group(1).rstrip("_"), m.group(2)
    try:
        ds = netCDF4.Dataset(fp)
        try:
            v = np.asarray(ds.variables["calibration_constant_0"][:], dtype=float)
        finally:
            ds.close()
    except Exception:  # noqa: BLE001
        return None
    fin = v[np.isfinite(v)]
    if not fin.size:
        return None
    return key, date, float(np.median(fin))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("l2dir", type=Path)
    ap.add_argument("out_csv", type=Path)
    ap.add_argument("--start", default="20260301")
    ap.add_argument("--end", default="20260531")
    ap.add_argument("--all-keys", action="store_true",
                    help="extract every station found, not just the census set")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    census = json.loads(R.CENSUS.read_text(encoding="utf-8"))
    wigos_wanted = {R._key(s).rsplit("_", 1)[0] for s in census}  # dir names are the bare WIGOS

    # gather candidate files within the date window, restricted to wanted stations
    files = []
    for wig_dir in sorted(args.l2dir.iterdir()):
        if not wig_dir.is_dir():
            continue
        if not args.all_keys and wig_dir.name not in wigos_wanted:
            continue
        for fp in wig_dir.glob("20*/*/L2_*.nc"):  # any year subdir (2025, 2026, ...)
            m = _FNAME.search(fp.name)
            if m and args.start <= m.group(2) <= args.end:
                files.append(str(fp))
    print(f"{len(files)} L2 files to read (window {args.start}..{args.end})", flush=True)

    # resume: keep existing rows
    have = set()
    rows = []
    if args.out_csv.exists():
        with args.out_csv.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                have.add((r["key"], r["date"]))
                rows.append((r["key"], r["date"], r["op_coeff"]))
    todo = [fp for fp in files if (lambda m: (m.group(1).rstrip("_"), m.group(2)) not in have)
            (_FNAME.search(os.path.basename(fp)))]
    print(f"{len(have)} cached, {len(todo)} to read", flush=True)

    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(_read_one, fp) for fp in todo]):
            r = fut.result()
            done += 1
            if r:
                rows.append((r[0], r[1], f"{r[2]:.8e}"))
            if done % 2000 == 0:
                print(f"  {done}/{len(todo)}", flush=True)

    rows.sort(key=lambda x: (x[0], x[1]))
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "date", "op_coeff"])
        w.writerows(rows)
    n_keys = len({r[0] for r in rows})
    print(f"wrote {len(rows)} rows ({n_keys} stations) -> {args.out_csv}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
