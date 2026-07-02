"""dashboard_to_calib.py — feed the operational DASHBOARD calibration (fullcal_l1_2026) into the
paper-validation report.

The dashboard (scripts/run_all_l1_2026.py) calibrates every stream from the native L1 archive over
2025-2026 and writes, per stream <WMO>_<ident>:
    <key>_cal.csv     date, method, flag, cal_value, uncertainty, ...   (raw per-night/day)
    <key>_kalman.csv  method, date, kalman, kalman_std                  (E-PROFILE best estimate)
In the dashboard, BOTH methods store the absolute Wiegner lidar constant C_L (cloud C_L =
applied_default / O'Connor_C, with applied_default = INSTRUMENT_CAL_DEFAULT).

The paper report (validation/paper/intercompare.py) reads, per benchmark channel
<WMO>_<ident>_<calib>, a CSV `time, C_daily, C_daily_std, C_kalman, C_kalman_std` whose value is the
**absolute Wiegner lidar constant C_L for BOTH methods** (single physical constant everywhere —
displayed as C_L in every figure). The consumers derive their multipliers from it:
  - Rayleigh : corr = calibration_constant_0 / C_L
  - Cloud    : corr = INSTRUMENT_CAL_DEFAULT / C_L   (the O'Connor multiplier, computed at use time)

We write the SAME dashboard-derived series to `<key>.csv`, `<key>_L1.csv` and `<key>_L2.csv` for every
benchmark + EARLINET channel, so whichever suffix a reader picks (intercompare calibLevel, the
calibration-time-series figure, earlinet.load_calib_series) it gets the dashboard numbers.

Usage:  python -m validation.paper.dashboard_to_calib
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper.calib_benchmark import BENCHMARK, key_of

FC = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/fullcal_l1_2026")
CALIB = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python/calib")
CALIB.mkdir(parents=True, exist_ok=True)

_SUCCESS = ("1", "1.0", "0.5")


def _to_report_value(itype, calib, cl):
    """Dashboard value -> report value: the lidar constant C_L, unchanged, for BOTH methods."""
    if cl in (None, 0) or cl != cl:
        return None
    return cl


def _read_daily(stream_key, method):
    """Raw per-night/day successful cal_value (+uncertainty) from <key>_cal.csv, date->(val, unc)."""
    f = FC / stream_key / f"{stream_key}_cal.csv"
    out = {}
    if not f.is_file():
        return out
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != method or str(r.get("flag")) not in _SUCCESS:
            continue
        try:
            v = float(r["cal_value"]); u = float(r.get("uncertainty") or 0.0)
        except (TypeError, ValueError):
            continue
        if v > 0 and v == v:
            out[r["date"]] = (v, u)
    return out


def _read_kalman(stream_key, method):
    """E-PROFILE Kalman best estimate from <key>_kalman.csv, date->(kalman, kalman_std)."""
    f = FC / stream_key / f"{stream_key}_kalman.csv"
    out = {}
    if not f.is_file():
        return out
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r.get("method") != method:
            continue
        try:
            v = float(r["kalman"]); s = float(r.get("kalman_std") or 0.0)
        except (TypeError, ValueError):
            continue
        if v > 0 and v == v:
            out[r["date"]] = (v, s)
    return out


def _fmt_date(yyyymmdd):
    s = str(yyyymmdd)
    return f"{s[:4]}-{s[4:6]}-{s[6:8]}" if len(s) >= 8 else s


def write_channel(ch):
    itype, calib = ch["itype"], ch["calib"]
    stream_key = f"{ch['wmo']}_{ch['ident']}"
    key = key_of(ch)
    daily = _read_daily(stream_key, calib)
    kal = _read_kalman(stream_key, calib)
    if not kal and not daily:
        print(f"  {key:42s} -> no dashboard {calib} series (skipped)")
        return False
    dates = sorted(set(daily) | set(kal))
    rows = []
    for d in dates:
        cd = cds = ck = cks = ""
        if d in daily:
            v, u = daily[d]
            rv = _to_report_value(itype, calib, v)
            if rv is not None:
                cd = f"{rv:.6e}"
                # preserve relative uncertainty through the C_L<->C reciprocal
                cds = f"{abs(rv) * (u / v):.6e}" if v else ""
        if d in kal:
            v, s = kal[d]
            rv = _to_report_value(itype, calib, v)
            if rv is not None:
                ck = f"{rv:.6e}"
                cks = f"{abs(rv) * (s / v):.6e}" if v else ""
        rows.append((_fmt_date(d), cd, cds, ck, cks))
    body = ["time,C_daily,C_daily_std,C_kalman,C_kalman_std"]
    body += [",".join(r) for r in rows]
    txt = "\n".join(body) + "\n"
    for suffix in ("", "_L1", "_L2"):
        (CALIB / f"{key}{suffix}.csv").write_text(txt, encoding="utf-8")
    nok_d = sum(1 for r in rows if r[1])
    nok_k = sum(1 for r in rows if r[3])
    print(f"  {key:42s} -> {len(rows):4d} rows ({nok_d} daily, {nok_k} kalman)  [{itype} {calib}]")
    return True


def main():
    seen = set()
    n = 0
    for st in BENCHMARK.values():
        for ch in st["channels"]:
            k = key_of(ch)
            if k in seen:
                continue
            seen.add(k)
            if write_channel(ch):
                n += 1
    print(f"DASHBOARD_TO_CALIB_DONE — wrote {n} channels into {CALIB}")


if __name__ == "__main__":
    main()
