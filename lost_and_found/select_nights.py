#!/usr/bin/env python3
"""
Select 50 random nights that calibrate successfully in BOTH configs.

Reads the two Phase-A CSVs (withfix / nofix), takes the intersection of nights with
flag in {1, 0.5} (successful / partially-clear), and writes a reproducible random
selection of up to 50 dates (YYYYMMDD, one per line) to selected_nights.txt.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh")
SUCCESS_FLAGS = {1.0, 0.5}
N_SELECT = 50
SEED = 42


def _successful_dates(csv_path: Path) -> set[str]:
    ok = set()
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                flag = float(row["flag"])
            except (ValueError, KeyError):
                continue
            if flag in SUCCESS_FLAGS:
                ok.add(row["date"])
    return ok


def main() -> int:
    wf = _successful_dates(BASE / "payerne_cl_withfix.csv")
    nf = _successful_dates(BASE / "payerne_cl_nofix.csv")
    common = sorted(wf & nf)

    print(f"successful withfix : {len(wf)}")
    print(f"successful nofix   : {len(nf)}")
    print(f"successful in both : {len(common)}")

    n = min(N_SELECT, len(common))
    rng = np.random.default_rng(SEED)
    chosen = sorted(rng.choice(common, size=n, replace=False).tolist()) if common else []

    out = BASE / "selected_nights.txt"
    out.write_text("\n".join(chosen) + ("\n" if chosen else ""))
    print(f"selected {len(chosen)} nights (seed={SEED}) -> {out}")
    if chosen:
        print("  " + ", ".join(chosen))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
