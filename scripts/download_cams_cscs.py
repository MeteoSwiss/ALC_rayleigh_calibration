#!/usr/bin/env python3
"""Bulk-download CAMS_Beta monthly files (0.4 deg, aerosol + T/RH) for OmB + WV.

Re-downloads the CAMS archive for a date range at the **native 0.4 deg** ADS grid
(``REGRID_TO_1DEG = False``) with the lean OmB variable set (aerosol backscatter
532/1064 nm + t/q/z/lnsp). The finer grid markedly improves the water-vapour
correction and the OmB near-surface comparison, especially at elevated/Alpine and
coastal sites where the 1 deg cell orography is unrepresentative.

Domain: the Europe+Arctic box ``AREA`` in download_cams_beta (80 N .. 27 N,
-30 E .. 45 E) covers 421/427 E-PROFILE census stations. The 6 far-flung
affiliates (Canada, Bonaire, New Zealand) are outside it and need a separate
regional download (or use sensitivity-only, which needs no CAMS aerosol).
Override the domain with the ALC_CAMS_AREA env var ("N,W,S,E").

Resumable: a month whose CAMS_Beta_YYYYMM.nc already exists AND carries aerosol
backscatter is skipped, so the job can be re-run after an interruption.

WHERE TO RUN (CSCS): the ADS download needs outbound HTTPS, which CSCS *compute*
nodes usually lack. Run on a login / data-mover node inside tmux/screen (the job
is network-bound, not compute-bound), exporting the proxy if required, e.g.:

    export ALC_CAMS_DIR=/capstor/scratch/<user>/CAMS
    export ADS_API_KEY=<your-ADS-key>          # or have ~/.cdsapirc
    export HTTPS_PROXY=http://proxy.cscs.ch:8080   # if CSCS requires a proxy
    conda activate alc
    python scripts/download_cams_cscs.py --start 202501 --end 202612

Storage: each 0.4 deg monthly file is ~5-8x the 1 deg one (~several GB); 24 months
is ~100-200 GB. The ADS is queued/rate-limited, so the full run can take days of
wall-clock; the per-month resume makes that safe.
"""
from __future__ import annotations

import argparse
import calendar
import os
import sys
import traceback
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from calibration.io import download_cams_beta as D  # noqa: E402
from calibration.io.cams import _has_backscatter  # noqa: E402


def month_list(start: str, end: str) -> list[str]:
    """Inclusive list of 'YYYYMM' from start to end."""
    y, m = int(start[:4]), int(start[4:6])
    ey, em = int(end[:4]), int(end[4:6])
    out = []
    while (y, m) <= (ey, em):
        out.append(f"{y:04d}{m:02d}")
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return out


def month_dates(ym: str, last_allowed: date) -> list[str]:
    """'YYYY-MM-DD' for every day of the month, capped at last_allowed.

    CAMS forecasts publish ~next-day, so for the current/most-recent month we cap
    at yesterday; requesting days the ADS does not have yet would error.
    """
    y, m = int(ym[:4]), int(ym[4:6])
    n = calendar.monthrange(y, m)[1]
    out = []
    for d in range(1, n + 1):
        day = date(y, m, d)
        if day <= last_allowed:
            out.append(day.strftime("%Y-%m-%d"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="202501", help="first month YYYYMM")
    ap.add_argument("--end", default="202612", help="last month YYYYMM")
    ap.add_argument("--out", default=os.environ.get("ALC_CAMS_DIR", "."),
                    help="output folder (default $ALC_CAMS_DIR or .)")
    ap.add_argument("--full", action="store_true",
                    help="download the full variable set (aerosol 355/532/1064 + "
                         "extinction + t/q/z/lnsp) instead of the lean OmB set")
    ap.add_argument("--lean", action="store_true",
                    help="download ONLY the calibration fields (t/q/z/lnsp, no aerosol) -> "
                         "enables the WV correction + molecular profile but NOT OmB. Use for "
                         "out-of-domain stations where the aerosol archive is not wanted.")
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the monthly file already has backscatter")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.lean:
        variables = D.CALIBRATION_VARIABLES   # t/q/z/lnsp only: calibration, no OmB
    elif args.full:
        variables = D.VARIABLES
    else:
        variables = D.OMB_VARIABLES
    yesterday = datetime.now(timezone.utc).date() - timedelta(days=1)   # UTC day convention

    print(f"[cams] domain AREA (N,W,S,E) = {D.AREA}  grid = native 0.4 deg "
          f"(REGRID_TO_1DEG={D.REGRID_TO_1DEG})")
    print(f"[cams] {len(variables)} variables -> {out_dir}")
    months = month_list(args.start, args.end)
    n_done = n_skip = n_fail = 0
    for ym in months:
        out_path = out_dir / f"CAMS_Beta_{ym}.nc"
        # A lean file has no backscatter by design, so "done" = exists; otherwise require backscatter.
        already_done = out_path.exists() and (args.lean or _has_backscatter(out_path))
        if already_done and not args.force:
            print(f"[skip] {out_path.name} (exists{'' if args.lean else ', has backscatter'})")
            n_skip += 1
            continue
        dates = month_dates(ym, yesterday)
        if not dates:
            print(f"[skip] {ym}: no days on/before {yesterday} yet")
            continue
        print(f"[get ] {out_path.name}  ({dates[0]} .. {dates[-1]}, {len(dates)} days)")
        try:
            D.download_to_netcdf(dates, out_path, variables=variables)
            ok = _has_backscatter(out_path)
            print(f"[ok  ] {out_path.name}  backscatter={'yes' if ok else 'NO'}")
            n_done += 1
        except Exception as exc:  # noqa: BLE001 — keep going, the run is resumable
            print(f"[FAIL] {ym}: {exc}")
            traceback.print_exc()
            n_fail += 1
    print(f"[cams] done: {n_done} downloaded, {n_skip} skipped, {n_fail} failed "
          f"of {len(months)} months")


if __name__ == "__main__":
    main()
