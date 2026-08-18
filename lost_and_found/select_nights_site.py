#!/usr/bin/env python3
"""
Select up to 50 random nights successful in BOTH configs, for a given site.

Reads C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/<site>/<site>_cl_{withfix,nofix}.csv, takes the
intersection of nights with flag in {1, 0.5}, and writes a reproducible random sample
(seed 42) to <site>/selected_nights.txt.

Usage:  python select_nights_site.py --site granada
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh")
SUCCESS_FLAGS = {1.0, 0.5}
N_SELECT = 50
SEED = 42


def _successful(csv_path: Path) -> set[str]:
    ok = set()
    if not csv_path.exists():
        return ok
    for row in csv.DictReader(open(csv_path, newline="", encoding="utf-8")):
        try:
            if float(row["flag"]) in SUCCESS_FLAGS:
                ok.add(row["date"])
        except (ValueError, KeyError):
            pass
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True)
    args = ap.parse_args()
    site = args.site

    wf = _successful(BASE / site / f"{site}_cl_withfix.csv")
    nf = _successful(BASE / site / f"{site}_cl_nofix.csv")
    common = sorted(wf & nf)
    print(f"[{site}] successful withfix={len(wf)} nofix={len(nf)} both={len(common)}")

    n = min(N_SELECT, len(common))
    rng = np.random.default_rng(SEED)
    chosen = sorted(rng.choice(common, size=n, replace=False).tolist()) if common else []

    out = BASE / site / "selected_nights.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(chosen) + ("\n" if chosen else ""))
    print(f"[{site}] selected {len(chosen)} nights (seed={SEED}) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
