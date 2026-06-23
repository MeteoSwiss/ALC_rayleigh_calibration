#!/usr/bin/env python3
"""Unified L1 calibration runner — Rayleigh AND liquid-cloud, one instrument stream at a time.

Terminology
-----------
- **stream** (the unit this runner loops over) = one physical INSTRUMENT, identified by its
  (WMO, identifier) pair; a single WMO site can host several (A/B/C...). It is independent of
  technique: a CL61 stream runs BOTH calibrations, a CL31 stream runs cloud only. [FR: « flux
  instrument » — un instrument unique (WMO, identifiant).]
- **series** (how the dashboard groups results) = one calibration METHOD on one instrument =
  (stream, method). A CL61 stream -> 2 series (Rayleigh + cloud); a CL31 stream -> 1 series.
  [FR: « série (instrument x méthode) » — une calibration d'un instrument par une technique.]

For the E-PROFILE L1 2026 archive (D:/E-PROFILE_L1_2026), this runs:
  * Rayleigh calibration  for CL61 / CHM15k / Mini-MPL   (per night, daily L1 files)
  * Liquid-cloud (O'Connor/Hopkin) calibration for CL31 / CL51 / CL61   (per day)
CL61 gets BOTH, which lets the dashboard cross-check the two methods.

Both methods write the SAME outputs per instrument <WMO>_<ident>:
  <key>/<key>_cal.csv      homogenized rows: date, method, flag, cal_value, uncertainty,
                           n_profiles, bottom_height, top_height, message
  <key>/<key>_kalman.csv   E-PROFILE Kalman best estimate per method (date, kalman, std)
  <key>/2026/ALC_calibration_<key>2026.nc   standard NetCDF; calibration_method tags each
                           row 0=Rayleigh / 1=Liquid_water_clouds (Rayleigh rows are written
                           by calibrate_rayleigh itself, cloud rows are written here).

Flags are the homogenized cloud/Rayleigh taxonomy (calibration.flags). The cloud headline
value is the O'Connor coefficient C = cal_median (~1 when well-calibrated): L1 carries no
applied calibration constant, so the absolute Wiegner C_L is undefined and C is the
monitorable quantity (its drift = calibration drift). The WV correction is applied by the
cloud core, gated by wavelength (CL31/CL51/CL61 only).

Date semantics: for a calendar date D, the Rayleigh row is the *preceding night* (D-1 -> D,
night-filtered; its NetCDF central time is floor(max kept-profile time), usually D) while the
cloud row is the *full daytime of D*. Both are labelled date=D and share calendar day D, so the
dashboard/NetCDF treat them as the same day -- a method cross-check on a fast-changing day
compares slightly different windows.

Usage:
  python scripts/run_all_l1_2026.py [--start 20260301] [--end 20260531]
        [--types CL31,CL51,CL61,CHM15k,Mini-MPL] [--per-type 3] [--limit N] [--workers 6]

--per-type N keeps the N best-covered streams per instrument type (a balanced test subset);
--per-type 0 runs every stream of the selected types. Resumable: a stream whose _cal.csv
exists is skipped.
"""
from __future__ import annotations

# Single-threaded BLAS BEFORE numpy (main + every spawned worker re-imports this on Windows).
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import json
import logging
import math
import subprocess
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent

# Each stream runs in its OWN subprocess with a hard timeout, so one hang-prone cloud day
# (a known issue -- see validation/_cloud_probe.py) cannot freeze the whole run. A normal
# stream finishes in well under a minute; this only fires on a true hang. Override via env.
STREAM_TIMEOUT = int(os.environ.get("STREAM_TIMEOUT", "900"))  # seconds
sys.path.insert(0, str(REPO))

from calibration import (  # noqa: E402
    calibrate_rayleigh, CalibrationOptions, CalibrationResult,
    InstrumentInfo, DataLevel,
)
from calibration.config import InstrumentType  # noqa: E402
from calibration.cloud import CloudCalConfig  # noqa: E402
from calibration.cloud.calibration import (  # noqa: E402
    read_ceilometer_data, liquid_cloud_calibration_from_data, set_defaults)
from calibration.flags import cloud_flag, flag_label, dominant_cloud_reject_flag  # noqa: E402
from calibration.io.output import write_calibration_result, strip_calibration_method  # noqa: E402
from calibration.plotting import plot_cloud_diagnostics_compact  # noqa: E402
from monitoring.kalman import kalman_best_estimate  # noqa: E402  (self-contained leaf)

# PLOTS=1 emits a diagnostic PNG per SUCCESSFUL calibration (Rayleigh via plot_main, cloud via
# plot_cloud_diagnostics_compact). Env-controlled so it propagates to every per-stream subprocess.
PLOT_ENABLED = os.environ.get("PLOTS", "0") == "1"

# --- Paths / configuration --------------------------------------------------
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/fullcal_l1_2026")
CAMS = Path("D:/CAMS")
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"   # bundled 910 nm WV LUT
EPOCH = datetime(1970, 1, 1)

ITYPE = {
    "CL31": InstrumentType.CL31, "CL51": InstrumentType.CL51, "CL61": InstrumentType.CL61,
    "CHM15k": InstrumentType.CHM15k, "Mini-MPL": InstrumentType.MINI_MPL,
}
RAYLEIGH_TYPES = {"CL61", "CHM15k", "Mini-MPL"}
CLOUD_TYPES = {"CL31", "CL51", "CL61"}

CSV_FIELDS = ["date", "method", "flag", "cal_value", "uncertainty",
              "n_profiles", "bottom_height", "top_height", "message"]
_SUCCESS = (1, 1.0, 0.5)
_HK_NAN = {k: float("nan") for k in ("laser_life_time", "status_detector", "status_laser",
                                     "temperature_optical_module", "window_transmission",
                                     "optical_module_id")}


# --- Helpers ----------------------------------------------------------------
def _key(s):
    return f"{s['wmo']}_{s['ident']}"


def _info(s):
    return InstrumentInfo(
        site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
        instrument_type=ITYPE[s["type"]], latitude=s["lat"], longitude=s["lon"],
        altitude=s.get("alt", 0.0),
    )


def _l1_file(wmo, ident, d):
    ds = d.strftime("%Y%m%d")
    # The L1 archive is laid out as <root>/<wmo>/<YYYY>/<MM>/ -- use the date's own year (the tree
    # under E-PROFILE_L1_2026 holds both 2025 and 2026), not a hard-coded year.
    return L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"


def _days(start, end):
    d = start
    while d <= end:
        yield d
        d += timedelta(days=1)


def _epoch_days(d):
    return float((datetime(d.year, d.month, d.day) - EPOCH).days)


def _is_success(flag):
    try:
        return float(flag) in _SUCCESS
    except (TypeError, ValueError):
        return False


def _write_csv_atomic(path, fieldnames, rows):
    """Write a CSV via temp-file + os.replace so a killed run never leaves a partial file
    under the final name (which existence-only resume would wrongly treat as 'done')."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


# --- Rayleigh (per night) ---------------------------------------------------
def _do_rayleigh(s, start, end):
    info = _info(s)
    key = _key(s)
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = OUT / key          # calibrate_rayleigh writes its NetCDF here (method 0)
    o.cams_folder = CAMS
    o.plot_all = False
    o.plot_main = PLOT_ENABLED           # PLOTS=1 -> Rayleigh diagnostic PNG per success
    rows = []
    for d in _days(start, end):
        if not _l1_file(s["wmo"], s["ident"], d).exists():
            continue
        ds = d.strftime("%Y%m%d")
        try:
            r = calibrate_rayleigh(ds, info, o)
            rows.append(dict(date=ds, method="rayleigh", flag=r.flag, cal_value=r.lidar_constant,
                             uncertainty=r.uncertainty, n_profiles="",
                             bottom_height=r.calibration_bottom_height,
                             top_height=r.calibration_top_height, message=r.message))
        except Exception as exc:  # noqa: BLE001 - one bad night must not kill the stream
            rows.append(dict(date=ds, method="rayleigh", flag=-99, cal_value=-1, uncertainty=0,
                             n_profiles="", bottom_height=None, top_height=None,
                             message=f"{type(exc).__name__}: {exc}"))
    return rows


# --- Cloud (per day) --------------------------------------------------------
def _do_cloud(s, start, end):
    info = _info(s)
    key = _key(s)
    # Drop any cloud (method=1) rows already in the NetCDF for THIS date window, so a re-run leaves no
    # stale row on days that no longer calibrate (the writer only overwrites on success). Scoped to the
    # window so other periods in the same yearly file are untouched; Rayleigh (method=0) is untouched.
    try:
        strip_calibration_method(OUT / key, info, method=1,
                                 t_start=_epoch_days(start), t_end=_epoch_days(end) + 1.0)
    except Exception:  # noqa: BLE001 - a cleanup failure must not abort the calibration
        pass
    rows = []
    for d in _days(start, end):
        fp = _l1_file(s["wmo"], s["ident"], d)
        if not fp.exists():
            continue
        ds = d.strftime("%Y%m%d")
        try:
            cfg = set_defaults(CloudCalConfig(
                nc_file=str(fp), instrument=s["type"], apply_wv_correction=True,
                apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
                cams_folder=str(CAMS), abs_cs_lookup_table=str(WV_LUT),
                station_latitude=s["lat"], station_longitude=s["lon"],
                average_time_s=30.0, average_range_m=10.0,   # finer cadence -> more valid cloud cals
            ))
            # Read FIRST so we can tell NO DATA (file present but no usable signal -> flag 0)
            # apart from NO CLOUD (data fine, but clear sky / no liquid cloud -> flag -1).
            data, status = read_ceilometer_data(cfg.nc_file, cfg)
            beta = getattr(data, "beta", None) if data is not None else None
            if status != 0 or beta is None or not np.any(np.isfinite(np.asarray(beta, dtype=float))):
                rows.append(dict(date=ds, method="cloud", flag=0, cal_value=-1, uncertainty=0,
                                 n_profiles=0, bottom_height=None, top_height=None,
                                 message="No data (no usable signal)"))
                continue
            res = liquid_cloud_calibration_from_data(data, cfg)
            n = int(getattr(res, "n_profiles", 0))
            coef = float(res.cal_median)
            std = float(res.cal_std)
            flag = cloud_flag(n, coef, std)  # n==0 with data present -> -1 (no liquid cloud)
            if flag == -1:
                # No usable profile: report WHY -- the dominant filter rejection (window / energy /
                # above / below / ratio / cbh / consistency), or genuine -1 if nothing was rejected.
                flag, _reason, _rej = dominant_cloud_reject_flag(
                    getattr(res, "filter_stats", None), getattr(res, "cloud_stats", None),
                    getattr(res, "consistency_stats", None))
            # Headline value is the lidar constant C_L = applied_constant / C (Wiegner) -- the
            # operationally useful quantity, on the SAME scale as Rayleigh -- NOT the O'Connor
            # coefficient C. C_L's relative uncertainty equals the coefficient's (C_L = const / C).
            cl = float(res.lidar_constant)
            ok = _is_success(flag) and math.isfinite(cl) and cl > 0
            cl_unc = (cl * std / coef) if (ok and math.isfinite(coef) and coef != 0) else 0.0
            msg = (f"OK ({n} profiles)" if ok else flag_label(flag))
            rows.append(dict(date=ds, method="cloud", flag=flag,
                             cal_value=(cl if ok else -1),
                             uncertainty=(cl_unc if ok else 0), n_profiles=n,
                             bottom_height=None, top_height=None, message=msg))
            if ok:
                # Same standard NetCDF as Rayleigh, tagged method=1 (cloud, per-day window).
                rr = CalibrationResult(lidar_constant=cl, flag=flag, uncertainty=cl_unc,
                                       calibration_bottom_height=None, calibration_top_height=None,
                                       message="liquid-cloud (O'Connor) calibration")
                write_calibration_result(
                    output_dir=OUT / key, info=info, result=rr,
                    date_epoch=_epoch_days(d) + 0.5, time_start=_epoch_days(d),
                    time_end=_epoch_days(d) + 1.0,
                    wavelength_nm=info.instrument_type.wavelength_nm,
                    housekeeping=_HK_NAN, method=1,
                )
            # Diagnostic image for successes AND informative rejections (a cloud was present but a
            # filter rejected it: flags -20..-26). Genuine clear sky (-1) / no data (0) get none --
            # there is nothing to diagnose, and over a multi-year run that would be a huge image flood.
            if PLOT_ENABLED and (ok or flag <= -20):
                pdir = OUT / key / "plots" / info.wmo_id / ds[:4]
                pdir.mkdir(parents=True, exist_ok=True)
                try:
                    plot_cloud_diagnostics_compact(
                        data, res,
                        title=f"{s.get('site', key)} ({info.wmo_id}) — {ds} — Cloud diagnostics "
                              f"[{flag_label(flag) if not ok else 'OK'}]",
                        save_path=pdir / f"{ds}_{info.wmo_id}_cloud_diag_compact.png")
                except Exception:  # noqa: BLE001 - a plot failure must not lose the calibration
                    pass
        except Exception as exc:  # noqa: BLE001
            msg = str(exc)
            low = msg.lower()
            # Only a genuinely MISSING CAMS file is -4 "missing aux data"; other WV failures
            # (e.g. a missing LUT, or a real bug in the WV chain) must NOT hide behind -4.
            flag = -4 if ("no cams file" in low or "cams_beta" in low) else -99
            rows.append(dict(date=ds, method="cloud", flag=flag, cal_value=-1, uncertainty=0,
                             n_profiles=0, bottom_height=None, top_height=None,
                             message=f"{type(exc).__name__}: {msg[:140]}"))
    return rows


# --- Kalman best estimate (per method) --------------------------------------
def _kalman_rows(rows):
    out = []
    for method in ("rayleigh", "cloud"):
        pts = []
        for r in rows:
            if r["method"] != method or not _is_success(r["flag"]):
                continue
            try:
                v = float(r["cal_value"])
            except (TypeError, ValueError):
                continue
            if math.isfinite(v) and v > 0:
                pts.append((datetime.strptime(r["date"], "%Y%m%d"), v))
        if len(pts) < 5:
            continue
        pts.sort()
        grid, state, std = kalman_best_estimate([p[0] for p in pts], [p[1] for p in pts])
        for g, sv, sd in zip(grid, state, std):
            out.append(dict(method=method,
                            date=np.datetime_as_string(g, unit="D").replace("-", ""),
                            kalman=f"{sv:.6e}", kalman_std=f"{sd:.6e}"))
    return out


def _preserve_other_methods(csv_path, methods):
    """Return rows from an existing per-stream CSV whose method was NOT recomputed this run, so a
    partial run (e.g. --methods cloud) keeps the other method's rows instead of dropping them."""
    if not csv_path.exists():
        return []
    keep = []
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("method") not in methods:
                    keep.append({k: r.get(k, "") for k in CSV_FIELDS})
    except (OSError, csv.Error):
        return []
    return keep


# --- One instrument stream --------------------------------------------------
def _process_stream(payload):
    s, start, end, methods = payload
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    key = _key(s)
    sdir = OUT / key
    sdir.mkdir(parents=True, exist_ok=True)

    rows = []
    if "rayleigh" in methods and s["type"] in RAYLEIGH_TYPES:
        rows += _do_rayleigh(s, start, end)
    if "cloud" in methods and s["type"] in CLOUD_TYPES:
        rows += _do_cloud(s, start, end)
    # preserve rows for any method we did NOT recompute (keeps Rayleigh during a cloud-only rerun)
    rows += _preserve_other_methods(sdir / f"{key}_cal.csv", methods)
    rows.sort(key=lambda r: (r["method"], r["date"]))

    _write_csv_atomic(sdir / f"{key}_cal.csv", CSV_FIELDS, rows)
    _write_csv_atomic(sdir / f"{key}_kalman.csv",
                      ["method", "date", "kalman", "kalman_std"], _kalman_rows(rows))

    n_ok = sum(1 for r in rows if _is_success(r["flag"]))
    return key, s["type"], len(rows), n_ok


# --- Stream selection -------------------------------------------------------
def select(census, types, per_type, limit, start, end):
    s0, s1 = start.strftime("%Y%m%d"), end.strftime("%Y%m%d")
    streams = [s for s in census
               if s.get("type") in types
               and s.get("first", "") <= s1 and s.get("last", "") >= s0]
    if per_type:
        sub = []
        for t in types:
            same = sorted([s for s in streams if s["type"] == t], key=lambda s: -s.get("n_days", 0))
            sub += same[:per_type]
        streams = sub
    streams.sort(key=lambda s: _key(s))
    if limit:
        streams = streams[:limit]
    return streams


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", default="20260301")
    ap.add_argument("--end", default="20260531")
    ap.add_argument("--types", default="CL31,CL51,CL61,CHM15k,Mini-MPL")
    ap.add_argument("--per-type", type=int, default=3,
                    help="keep the N best-covered streams per type (0 = all)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--stream", default=None,
                    help="internal: process ONE stream by '<wmo>_<ident>' key and exit")
    ap.add_argument("--methods", default="rayleigh,cloud",
                    help="comma list of methods to (re)compute: rayleigh, cloud. Rows for a method NOT "
                         "listed are preserved from the existing per-stream CSV, so e.g. "
                         "'--methods cloud --force' updates only the cloud rows and keeps Rayleigh.")
    ap.add_argument("--force", action="store_true",
                    help="re-process streams even if a *_cal.csv already exists (needed for in-place "
                         "reruns such as a cloud-only update that keeps the previous Rayleigh output)")
    args = ap.parse_args()

    start = datetime.strptime(args.start, "%Y%m%d")
    end = datetime.strptime(args.end, "%Y%m%d")
    methods = {m.strip().lower() for m in args.methods.split(",") if m.strip()}
    if not methods or not methods <= {"rayleigh", "cloud"}:
        ap.error("--methods must be a comma list drawn from: rayleigh, cloud")
    census = json.loads(CENSUS.read_text(encoding="utf-8"))

    # --- single-stream worker mode: spawned as an isolated, killable subprocess ----------
    if args.stream:
        warnings.filterwarnings("ignore")
        logging.getLogger().setLevel(logging.CRITICAL)
        s = next((x for x in census if _key(x) == args.stream), None)
        if s is None:
            print(f"stream {args.stream} not in census", flush=True)
            return
        _process_stream((s, start, end, methods))
        return

    types = [t.strip() for t in args.types.split(",")]
    streams = select(census, types, args.per_type, args.limit, start, end)
    # drop streams that have no work for the requested methods (e.g. CHM15k/Mini-MPL under --methods cloud)
    relevant = (RAYLEIGH_TYPES if "rayleigh" in methods else set()) | (CLOUD_TYPES if "cloud" in methods else set())
    streams = [s for s in streams if s["type"] in relevant]
    OUT.mkdir(parents=True, exist_ok=True)
    done = {p.parent.name for p in OUT.glob("*/*_cal.csv")}
    todo = streams if args.force else [s for s in streams if _key(s) not in done]
    print(f"window {args.start}..{args.end} | methods={','.join(sorted(methods))} | "
          f"{len(streams)} streams selected; {len(done)} with output; {len(todo)} to do"
          f"{' (forced)' if args.force else ''}; {args.workers} workers "
          f"(per-stream timeout {STREAM_TIMEOUT}s)", flush=True)

    def _run_stream(s):
        """Run ONE stream in a separate process; kill it if it hangs past the timeout."""
        key = _key(s)
        cmd = [sys.executable, str(Path(__file__).resolve()), "--stream", key,
               "--start", args.start, "--end", args.end, "--methods", args.methods]
        try:
            subprocess.run(cmd, timeout=STREAM_TIMEOUT, check=False,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            return key, "TIMEOUT (killed)"
        return key, ("ok" if (OUT / key / f"{key}_cal.csv").exists() else "no-output")

    # The thread pool only SUPERVISES the per-stream subprocesses (the heavy numpy/netCDF
    # work runs in the children). A hung child is killed at the timeout and the run goes on
    # -- unlike a ProcessPoolExecutor worker, which hangs the whole pool indefinitely.
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_run_stream, s): s for s in todo}
        for k, fut in enumerate(as_completed(futs), 1):
            s = futs[fut]
            try:
                key, status = fut.result()
                print(f"[{k}/{len(todo)}] {key} ({s['type']}): {status}", flush=True)
            except Exception as exc:  # noqa: BLE001
                print(f"[{k}/{len(todo)}] FAILED {_key(s)}: {type(exc).__name__}: {exc}", flush=True)

    Path(OUT / "ALL_DONE.flag").write_text("done")
    print("L1_2026_DONE", flush=True)


if __name__ == "__main__":
    main()
