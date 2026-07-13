"""Golden diff of two fullcal output dirs (chunk read-once validation, M4).

Compares, per stream key: the per-stream CSVs (_cal, _kalman, _hk, _status, _sens, _omb) as
text first, then field-by-field with a float tolerance; the yearly ALC_calibration NetCDFs
(time / lidar_constant / calibration_method); and the classification NetCDFs (count + exact
target_classification equality on up to 3 sampled days).

Run:  python validation/diff_fullcal.py DIR_A DIR_B KEY [KEY ...]
Exit 0 when everything matches (within tolerance), 1 otherwise.
"""
import csv
import sys
from pathlib import Path

import numpy as np

RTOL = 1e-9          # float fields must agree to this relative tolerance
CSVS = ("_cal.csv", "_kalman.csv", "_hk.csv", "_status.csv", "_sens.csv", "_omb.csv")


def _read_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def _cmp_cell(a, b):
    if a == b:
        return True
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return False
    if np.isnan(fa) and np.isnan(fb):
        return True
    return np.isclose(fa, fb, rtol=RTOL, atol=0.0)


def _cmp_csv(pa, pb, label, report):
    ea, eb = pa.exists(), pb.exists()
    if not ea and not eb:
        return True
    if ea != eb:
        report.append(f"  {label}: present only in {'A' if ea else 'B'}")
        return False
    if pa.read_bytes() == pb.read_bytes():
        return True
    ra, rb = _read_rows(pa), _read_rows(pb)
    if len(ra) != len(rb):
        report.append(f"  {label}: row count {len(ra)} vs {len(rb)}")
        return False
    for i, (rowa, rowb) in enumerate(zip(ra, rb)):
        if len(rowa) != len(rowb) or not all(_cmp_cell(x, y) for x, y in zip(rowa, rowb)):
            report.append(f"  {label}: first differing row {i}:\n    A: {rowa}\n    B: {rowb}")
            return False
    return True   # differed as text but every field within tolerance


def _cmp_yearly_nc(da, db, key, report):
    import netCDF4
    ok = True
    for nca in sorted(da.glob("*/ALC_calibration_*.nc")):
        ncb = db / nca.relative_to(da)
        if not ncb.exists():
            report.append(f"  {nca.name}: missing in B")
            ok = False
            continue
        with netCDF4.Dataset(nca) as a, netCDF4.Dataset(ncb) as b:
            for var in ("time", "lidar_constant", "calibration_method"):
                if var not in a.variables or var not in b.variables:
                    continue
                va = np.ma.filled(a.variables[var][:].astype("f8"), np.nan)
                vb = np.ma.filled(b.variables[var][:].astype("f8"), np.nan)
                if va.shape != vb.shape or not np.allclose(va, vb, rtol=RTOL, equal_nan=True):
                    report.append(f"  {nca.name}:{var}: arrays differ "
                                  f"(shapes {va.shape} vs {vb.shape})")
                    ok = False
    return ok


def _cmp_classification(da, db, report):
    import netCDF4
    fa = sorted((da / "classification").rglob("*_classification.nc")) if (da / "classification").is_dir() else []
    fb = sorted((db / "classification").rglob("*_classification.nc")) if (db / "classification").is_dir() else []
    na = {p.name for p in fa}
    nb = {p.name for p in fb}
    if na != nb:
        report.append(f"  classification: file sets differ (A-only {sorted(na - nb)[:3]}, "
                      f"B-only {sorted(nb - na)[:3]})")
        return False
    ok = True
    for pa in fa[:: max(1, len(fa) // 3)][:3]:      # sample up to 3 days spread over the window
        pb = next(p for p in fb if p.name == pa.name)
        with netCDF4.Dataset(pa) as a, netCDF4.Dataset(pb) as b:
            ta = np.asarray(a.variables["target_classification"][:])
            tb = np.asarray(b.variables["target_classification"][:])
            if ta.shape != tb.shape or not np.array_equal(ta, tb):
                report.append(f"  classification {pa.name}: target_classification differs")
                ok = False
    return ok


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    dir_a, dir_b, keys = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    failures = 0
    for key in keys:
        da, db = dir_a / key, dir_b / key
        report = []
        ok = True
        for suffix in CSVS:
            ok &= _cmp_csv(da / f"{key}{suffix}", db / f"{key}{suffix}", f"{key}{suffix}", report)
        ok &= _cmp_yearly_nc(da, db, key, report)
        ok &= _cmp_classification(da, db, report)
        print(f"{key}: {'IDENTICAL (within tol)' if ok else 'DIFFERS'}")
        for line in report:
            print(line)
        failures += 0 if ok else 1
    print(f"\n{len(keys) - failures}/{len(keys)} streams identical")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
