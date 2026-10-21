# -*- coding: utf-8 -*-
"""q16g (task 1b/1c): provenance of the CHM15k L1 rcs_0 over 2025-2026.

Two questions about the processing UPSTREAM of the calibration code:

1. Did anything in the L1 conversion change over time (raw2l1 / hermes version, firmware, overlap
   flag, rcs_0 units or long_name, optical module) and does it coincide with the network-wide
   same-month jumps 2026 vs 2025 found in q16b, or with the season?
   -> scan the global and rcs_0 attributes of one L1 file per stream every ~7 days.

2. Does raw2l1 change the signal? raw2l1 (reader lufft_chm15k_nimbus.calc_pr2) writes
   rcs_0 = beta_raw for firmware >= 0.7 (pass-through; beta_att/c_cal for netcdf_mode 1).
   Test it on data: Payerne has the raw Lufft daily files (A:/CHM15k/PAY, same module TUB140016)
   for 2025-2026 -> ratio rcs_0(L1) / beta_raw(raw) on matched profiles, 1 day every ~5 days.

Outputs (q16 folder): q16g_l1_attrs.csv (one row per sampled file), q16g_l1_attr_changes.csv
(per stream: every attribute change with its date), q16g_raw_vs_l1.csv.
"""
from __future__ import annotations

import csv
import glob
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import date, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import L1_ROOT, chm15k_streams, stream_key  # noqa: E402

Q16 = Path("C:/DATA/Projects/202610_CHM15k_stability/q16")
N_WORKERS = 4
GLOBAL_ATTRS = ("history", "hermes_history", "instrument_firmware_version", "instrument_serial_number",
                "optical_module_id", "overlap_function", "overlap_is_corrected", "title", "institution",
                "source", "Conventions", "instrument_type", "comment", "references")
RCS_ATTRS = ("units", "long_name", "overlap", "detection_mode")
SAMPLE_STEP_D = 7


def l1_file(wmo, ident, d):
    return L1_ROOT / wmo / f"{d.year}" / f"{d.month:02d}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"


def scan_stream(args):
    """Attributes of one L1 file every SAMPLE_STEP_D days (first existing file within the step)."""
    import netCDF4
    wmo, ident = args
    rows = []
    d = date(2025, 1, 1)
    end = date(2026, 9, 30)
    while d <= end:
        f = None
        for k in range(SAMPLE_STEP_D):
            cand = l1_file(wmo, ident, d + timedelta(days=k))
            if cand.exists():
                f = cand
                break
        if f is not None:
            try:
                with netCDF4.Dataset(f) as ds:
                    r = dict(stream=stream_key(wmo, ident), file_date=f.stem[-8:])
                    for a in GLOBAL_ATTRS:
                        r[a] = str(ds.getncattr(a)) if a in ds.ncattrs() else ""
                    # raw2l1 writes its version after the date in history: keep the version only
                    r["raw2l1_version"] = r["history"].split("raw2l1")[-1].strip() if "raw2l1" in r["history"] else ""
                    v = ds["rcs_0"]
                    for a in RCS_ATTRS:
                        r[f"rcs_{a}"] = str(v.getncattr(a)) if a in v.ncattrs() else ""
                    for name in ("l0_wavelength", "range_resol", "time_resol", "tilt_angle"):
                        if name in ds.variables:
                            x = np.ma.filled(np.asarray(ds[name][:], dtype="f8").ravel(), np.nan)
                            r[name] = float(np.nanmedian(x)) if x.size else np.nan
                        else:
                            r[name] = np.nan
                    r["n_range"] = int(ds.dimensions["range"].size)
                    rng = np.asarray(ds["range"][:2], dtype="f8")
                    r["range0"] = float(rng[0])
            except Exception as e:  # noqa: BLE001
                r = dict(stream=stream_key(wmo, ident), file_date=f.stem[-8:], history=f"ERROR {e}")
            rows.append(r)
        d += timedelta(days=SAMPLE_STEP_D)
    return rows


def attr_changes(rows):
    """Per stream: dates at which any tracked attribute changes value."""
    keys = ("raw2l1_version", "hermes_history", "instrument_firmware_version", "optical_module_id",
            "overlap_function", "overlap_is_corrected", "rcs_units", "rcs_long_name", "l0_wavelength",
            "range_resol", "n_range", "range0", "title", "instrument_serial_number")
    out = []
    by = {}
    for r in rows:
        by.setdefault(r["stream"], []).append(r)
    for st, rr in by.items():
        rr = sorted(rr, key=lambda x: x["file_date"])
        for k in keys:
            prev = None
            for r in rr:
                v = r.get(k, "")
                v = f"{v:.4g}" if isinstance(v, float) else str(v)
                if prev is not None and v != prev:
                    out.append(dict(stream=st, attribute=k, date=r["file_date"], old=prev, new=v))
                prev = v
    return out


# ------------------------------------------------------------------------------------------------
# raw Lufft vs L1 at Payerne
# ------------------------------------------------------------------------------------------------
def raw_vs_l1(day):
    import netCDF4
    raw = sorted(glob.glob(f"A:/CHM15k/PAY/{day:%Y}/{day:%m}/{day:%Y%m%d}_pay_CHM*.nc"))
    l1 = l1_file("0-20000-0-06610", "A", day)
    if not raw or not l1.exists():
        return None
    try:
        with netCDF4.Dataset(raw[0]) as d:
            tr = np.asarray(d["time"][:], dtype="f8") / 86400.0 - (date(1970, 1, 1) - date(1904, 1, 1)).days
            br = np.ma.filled(np.asarray(d["beta_raw"][:], dtype="f8"), np.nan)
            rr = np.asarray(d["range"][:], dtype="f8")
            sw = str(d.getncattr("software_version"))
            serlom = str(d.getncattr("serlom"))
            scaling = float(d["scaling"][:])
            mode = int(d.getncattr("netcdf_mode")) if "netcdf_mode" in d.ncattrs() else -1
            ovl = str(d.getncattr("overlap_file")) if "overlap_file" in d.ncattrs() else ""
        with netCDF4.Dataset(l1) as d:
            tl = np.asarray(d["time"][:], dtype="f8")
            rl = np.ma.filled(np.asarray(d["rcs_0"][:], dtype="f8"), np.nan)
            rgl = np.asarray(d["range"][:], dtype="f8")
            hist = str(d.getncattr("history"))
    except Exception as e:  # noqa: BLE001
        return dict(day=f"{day:%Y%m%d}", error=str(e))
    # match profiles within 1 s
    il = np.searchsorted(tl, tr)
    il = np.clip(il, 1, len(tl) - 1)
    il = np.where(np.abs(tl[il - 1] - tr) < np.abs(tl[il] - tr), il - 1, il)
    ok = np.abs(tl[il] - tr) < 1.0 / 86400.0
    if ok.sum() < 10 or rr.size != rgl.size or np.max(np.abs(rr - rgl)) > 0.1:
        return dict(day=f"{day:%Y%m%d}", n_match=int(ok.sum()), error="no match / range mismatch")
    a = rl[il[ok]]
    b = br[ok]
    out = dict(day=f"{day:%Y%m%d}", n_match=int(ok.sum()), software=sw, serlom=serlom, scaling=scaling,
               netcdf_mode=mode, overlap_file=ovl, l1_history=hist)
    for tag, (z0, z1) in {"0.2-1.5km": (200, 1500), "1.5-8km": (1500, 8000), "8-15km": (8000, 15000)}.items():
        m = (rr >= z0) & (rr < z1)
        aa, bb = a[:, m], b[:, m]
        good = np.isfinite(aa) & np.isfinite(bb) & (np.abs(bb) > 0)
        ratio = aa[good] / bb[good]
        out[f"ratio_med_{tag}"] = float(np.median(ratio)) if ratio.size else np.nan
        out[f"ratio_p01_{tag}"] = float(np.percentile(ratio, 1)) if ratio.size else np.nan
        out[f"ratio_p99_{tag}"] = float(np.percentile(ratio, 99)) if ratio.size else np.nan
        out[f"maxabsdiff_{tag}"] = float(np.max(np.abs(aa[good] - bb[good]))) if ratio.size else np.nan
        out[f"sum_ratio_{tag}"] = float(np.sum(aa[good]) / np.sum(bb[good])) if ratio.size else np.nan
    return out


def write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    streams = [(w, i) for w, i, *_ in chm15k_streams()]
    print(f"scanning {len(streams)} CHM15k streams", flush=True)
    rows = []
    with ProcessPoolExecutor(N_WORKERS) as ex:
        for k, rr in enumerate(ex.map(scan_stream, streams)):
            rows.extend(rr)
            if k % 20 == 0:
                print(f"  {k}/{len(streams)} streams, {len(rows)} files", flush=True)
    write_csv(rows, Q16 / "q16g_l1_attrs.csv")
    write_csv(attr_changes(rows), Q16 / "q16g_l1_attr_changes.csv")
    days = [date(2025, 1, 2) + timedelta(days=k) for k in range(0, 640, 5)]
    with ProcessPoolExecutor(N_WORKERS) as ex:
        rv = [r for r in ex.map(raw_vs_l1, days) if r]
    write_csv(rv, Q16 / "q16g_raw_vs_l1.csv")
    print(f"raw vs L1: {len(rv)} days", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
