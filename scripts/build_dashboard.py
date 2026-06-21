#!/usr/bin/env python3
"""Build the calibration monitoring dashboard (static site) from the full-network outputs.

  1. Index the per-station <key>_cl.csv files into a single SQLite file.
  2. Render a summary page (index.html) + one page per instrument.

Run after run_all_l2monthly.py (or on a schedule). The output is a folder of static
files; copy it to an internal web dir to share with the team.

Examples
--------
  # Full network, default paths:
  python scripts/build_dashboard.py

  # Quick look: index everything but render only the first 5 station pages:
  python scripts/build_dashboard.py --limit-pages 5 --open

  # A subset of instrument types, custom output folder:
  python scripts/build_dashboard.py --types CHM15k,CL61 --out D:/tmp/caldash
"""
from __future__ import annotations

import argparse
import sys
import time
import webbrowser
from pathlib import Path

# Allow running as a plain script (python scripts/build_dashboard.py) by putting the
# repo root on the path so `import monitoring` resolves.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from monitoring import config, index, render  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fullcal", type=Path, default=config.DEFAULT_FULLCAL_DIR,
                    help="folder containing <key>/<key>_cl.csv (default: %(default)s)")
    ap.add_argument("--manifest", type=Path, default=config.DEFAULT_MANIFEST,
                    help="stations manifest JSON (lat/lon/type) (default: %(default)s)")
    ap.add_argument("--l2dir", type=Path, default=config.DEFAULT_L2_DIR,
                    help="L2 archive for station name/country/institution (default: %(default)s)")
    ap.add_argument("--out", type=Path, default=config.DEFAULT_OUT_DIR,
                    help="output site folder (default: %(default)s)")
    ap.add_argument("--types", default=None,
                    help="comma-separated instrument types to include (default: all)")
    ap.add_argument("--limit", type=int, default=None,
                    help="index only the first N stations (debug)")
    ap.add_argument("--limit-pages", type=int, default=None,
                    help="render only the first N station pages (summary still covers all indexed)")
    ap.add_argument("--open", action="store_true", help="open the result in a browser when done")
    args = ap.parse_args()

    types = [t.strip() for t in args.types.split(",")] if args.types else None
    db_path = args.out / config.DB_NAME

    t0 = time.perf_counter()
    print(f"Indexing {args.fullcal} ...", flush=True)
    stats = index.build_index(args.fullcal, args.manifest, db_path, limit=args.limit,
                              types=types, l2_dir=args.l2dir)
    print(f"  {stats['n_stations']} stations, {stats['n_calibrations']} nights, "
          f"as of {stats['as_of']}  ({time.perf_counter() - t0:.1f}s)", flush=True)

    print("Rendering site ...", flush=True)
    site = render.build_site(db_path, args.out, limit_pages=args.limit_pages)
    print(f"  {site['n_pages']} station pages -> {site['out_dir']}  "
          f"({time.perf_counter() - t0:.1f}s total)", flush=True)

    index_html = args.out / "index.html"
    print(f"\nDone. Open: {index_html}", flush=True)
    if args.open:
        webbrowser.open(index_html.resolve().as_uri())


if __name__ == "__main__":
    main()
