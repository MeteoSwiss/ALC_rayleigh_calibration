"""method_test_flag_counts.py - Section 6.1 numbers not in method_comparison_summary.json.

  * Mini-MPL: "on 68-91 % of the days with data the signal does not drop by the required factor of 20 within
    300 m above the cloud peak" = cloud rows flagged 'peak not sharp above' / cloud rows (days with data).
    2026-10-09: 07014_A 90/113, 07110_A 104/147, 07617_A 102/149, 07774_A 125/138.
  * Lanzhot CL61 (0-203-10-LNG_A): cloud rows rejected with 'laser energy too low' (121 of 133).
Input: the forced-method run DATA/method_test_2026 (validation/paper/run_method_test.ps1).
Run:  python -m validation.paper.method_test_flag_counts
History: inline grep/awk counts of 2026-10-09, rewritten in Python.
"""
import csv
from pathlib import Path

RUN = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/method_test_2026")
COUNTS = {
    "peak not sharp above": ("0-20000-0-07014_A", "0-20000-0-07110_A", "0-20000-0-07617_A", "0-20000-0-07774_A"),
    "laser energy too low": ("0-203-10-LNG_A",),
}


def main():
    for pattern, keys in COUNTS.items():
        for key in keys:
            with open(RUN / key / f"{key}_cal.csv", newline="", encoding="utf-8", errors="replace") as fh:
                rows = [r for r in csv.reader(fh) if len(r) > 1 and r[1] == "cloud"]
            n = sum(1 for r in rows if pattern.lower() in ",".join(r).lower())
            print(f"{key}: '{pattern}' {n}/{len(rows)} cloud rows ({100 * n / max(len(rows), 1):.1f} %)")


if __name__ == "__main__":
    main()
