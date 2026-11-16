# -*- coding: utf-8 -*-
"""q16g (task 1c, history): does the L1 conversion (raw2l1, all versions used 2016-2024) change the CHM15k
signal? Payerne raw Lufft daily files (A:/CHM15k/PAY) vs the Swiss E-PROFILE L1 archive
(D:/E-PROFILE_L1_Switzerland_2017_2024), the 2nd of every month 2016-2024: ratio rcs_0 / beta_raw on matched
profiles, raw2l1 version and firmware. Output: q16g_raw_vs_l1_swiss.csv
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q16g_l1_provenance as prov  # noqa: E402

SWISS = Path("D:/E-PROFILE_L1_Switzerland_2017_2024")


def swiss_file(wmo, ident, d):
    return SWISS / wmo / f"{d.year}" / f"{d.month:02d}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"


def one(day):
    prov.l1_file = swiss_file          # same comparison, other L1 root
    return prov.raw_vs_l1(day)


def main():
    days = [date(y, m, 2) for y in range(2016, 2025) for m in range(1, 13)]
    with ProcessPoolExecutor(4) as ex:
        rows = [r for r in ex.map(one, days) if r]
    prov.write_csv(rows, prov.Q16 / "q16g_raw_vs_l1_swiss.csv")
    print(f"{len(rows)} days")


if __name__ == "__main__":
    main()
