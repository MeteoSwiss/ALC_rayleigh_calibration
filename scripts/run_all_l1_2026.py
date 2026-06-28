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
import time
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
# Paths are env-overridable (set them in ops/config.sh on the server); the defaults are the local
# Windows dev locations, so nothing changes unless the ALC_* vars are exported.
L1_ROOT = Path(os.environ.get("ALC_L1_ROOT", "D:/E-PROFILE_L1_2026"))
CENSUS = Path(os.environ.get("ALC_CENSUS", str(REPO / "validation" / "scope_l1_2026_census.json")))
OUT = Path(os.environ.get("ALC_FULLCAL_DIR", "C:/DATA/Projects/202606_E-PROFILE_calibration/fullcal_l1_2026"))
CAMS = Path(os.environ.get("ALC_CAMS_DIR", "D:/CAMS"))
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"   # bundled 910 nm WV LUT
EPOCH = datetime(1970, 1, 1)

ITYPE = {
    "CL31": InstrumentType.CL31, "CL51": InstrumentType.CL51, "CL61": InstrumentType.CL61,
    "CHM15k": InstrumentType.CHM15k, "Mini-MPL": InstrumentType.MINI_MPL,
}
RAYLEIGH_TYPES = {"CL61", "CHM15k", "Mini-MPL"}
CLOUD_TYPES = {"CL31", "CL51", "CL61"}

# Operational L2 archive (for the OmB operational-constant overlay) and the
# authoritative Kalman method per instrument type for the OmB/sensitivity C_L
# (Rayleigh where the molecular return is usable; cloud for the weak 910 nm
# CL31/CL51 that have no Rayleigh calibration).
L2_ROOT = Path(os.environ.get("ALC_L2_DIR", "D:/E-PROFILE_L2_2026"))
SENS_OMB_METHOD = {"CHM15k": "rayleigh", "CL61": "rayleigh", "Mini-MPL": "rayleigh",
                   "CL31": "cloud", "CL51": "cloud"}
OMB_FIELDS = ["date_start", "date_end", "wavelength", "median_bias_ours",
              "median_bias_op", "median_bias_ours_wv", "rms_ours", "n_obs"]
SENS_FIELDS = ["date_start", "date_end", "wavelength", "icao_alt_200", "icao_alt_2000",
               "icao_alt_4000", "sigma_night_3000", "n_days_night", "n_days_day"]

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


def _fmt_hms(seconds):
    """Format a duration in seconds as H:MM:SS for progress / ETA lines."""
    seconds = int(max(0.0, seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


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
            # Station outside the CAMS domain (nearest grid point too far) -> -10, distinct from a
            # genuinely MISSING CAMS file (-4). Other WV failures (missing LUT, real bug) -> -99.
            if "cams data too far" in low or "outside cams domain" in low:
                flag = -10
            elif "no cams file" in low or "cams_beta" in low:
                flag = -4
            else:
                flag = -99
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


def _preserve_existing_rows(csv_path, methods, start, end):
    """Rows from an existing per-stream CSV to KEEP unchanged: those of a method we did NOT recompute,
    OR of a recomputed method but with a date OUTSIDE the processed [start, end] window. This lets a
    daily / partial run replace only the dates it actually processed and ACCUMULATE history, instead
    of overwriting the whole file with just the processed window (which would wipe prior days)."""
    if not csv_path.exists():
        return []
    s, e = start.strftime("%Y%m%d"), end.strftime("%Y%m%d")   # YYYYMMDD sorts chronologically
    keep = []
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                method, date = r.get("method"), str(r.get("date", ""))
                if (method not in methods) or not (s <= date <= e):
                    keep.append({k: r.get(k, "") for k in CSV_FIELDS})
    except (OSError, csv.Error):
        return []
    return keep


# --- Instrument monitoring (housekeeping) -----------------------------------
# Daily means of the laser / optics / temperature essentials for the station-page monitoring panel.
# Canonical field -> candidate L1 variable names (first present wins); manufacturers differ -- Lufft
# CHM15k uses status_laser / temperature_optical_module / temperature_detector, Vaisala CL31/51/61
# use laser_energy / temperature_laser. Temperatures are stored in degC (the L1 files store K).
HK_FIELDS = ["date", "laser", "window", "temp_optics", "temp_internal", "temp_detector"]
_HK_SOURCES = {
    "laser":         ("status_laser", "laser_energy"),                    # laser power / pulse energy (%)
    "window":        ("window_transmission",),                           # window transmission (%)
    "temp_optics":   ("temperature_optical_module", "temperature_laser"),
    "temp_internal": ("temp_int",),
    "temp_detector": ("temperature_detector",),
}
_HK_TEMP = {"temp_optics", "temp_internal", "temp_detector"}             # K -> degC


def _do_hk(s, start, end):
    """Per-day daily-mean housekeeping (laser/optics/temperature) for one stream. Cheap: opens each
    L1 file but reads ONLY the small 1-D HK variables (never the rcs matrix). Runs for every
    instrument type and every day with data, independent of calibration -> the monitoring panel."""
    import netCDF4  # lazy: only this leaf needs it
    rows = []
    for d in _days(start, end):
        fp = _l1_file(s["wmo"], s["ident"], d)
        if not fp.exists():
            continue
        row = {"date": d.strftime("%Y%m%d")}
        try:
            nc = netCDF4.Dataset(str(fp))
            try:
                for field, cands in _HK_SOURCES.items():
                    val = float("nan")
                    for nm in cands:
                        if nm in nc.variables:
                            a = np.asarray(nc.variables[nm][:], dtype=float).ravel()
                            a = np.where(np.isfinite(a) & (a > -990.0), a, np.nan)
                            if np.isfinite(a).any():
                                m = float(np.nanmean(a))
                                val = (m - 273.15) if field in _HK_TEMP else m
                            break
                    row[field] = "" if val != val else f"{val:.3f}"   # val!=val -> NaN
            finally:
                nc.close()
        except Exception:  # noqa: BLE001 - one unreadable file must not kill the stream
            continue
        if any(row[f] for f in HK_FIELDS if f != "date"):
            rows.append(row)
    return rows


def _preserve_existing_hk(csv_path, start, end):
    """Existing _hk.csv rows OUTSIDE the processed window, so partial/daily runs accumulate history."""
    if not csv_path.exists():
        return []
    s, e = start.strftime("%Y%m%d"), end.strftime("%Y%m%d")
    keep = []
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if not (s <= str(r.get("date", "")) <= e):
                    keep.append({k: r.get(k, "") for k in HK_FIELDS})
    except (OSError, csv.Error):
        return []
    return keep


# --- OmB + sensitivity (consume the Kalman C_L) -----------------------------
def _read_kalman_csv(path):
    """Existing <key>_kalman.csv rows (method, date, kalman) — for --no-cal reuse."""
    if not path.exists():
        return []
    try:
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except (OSError, csv.Error):
        return []


def _kalman_map(kalman_rows, method):
    """date(YYYYMMDD) -> Kalman best-estimate C_L for the given method."""
    out = {}
    for r in kalman_rows:
        if r.get("method") == method:
            try:
                out[str(r["date"])] = float(r["kalman"])
            except (TypeError, ValueError, KeyError):
                pass
    return out


def _op_map(s, start, end):
    """date(YYYYMMDD) -> operational L2 calibration_constant_0 (median over the day)."""
    import netCDF4
    key = _key(s)
    out = {}
    for d in _days(start, end):
        ds = d.strftime("%Y%m%d")
        # L2 archives differ in depth: 3-level (<wmo>/<YYYY>/<MM>/) locally, but balfrin
        # nests 2-level (<wmo>/<YYYY>/). Try both before giving up.
        cands = [L2_ROOT / s["wmo"] / ds[:4] / ds[4:6] / f"L2_{key}{ds}.nc",
                 L2_ROOT / s["wmo"] / ds[:4] / f"L2_{key}{ds}.nc"]
        fp = next((c for c in cands if c.exists()), None)
        if fp is None:
            continue
        try:
            nc = netCDF4.Dataset(str(fp))
            try:
                v = np.asarray(nc.variables["calibration_constant_0"][:], dtype=float)
            finally:
                nc.close()
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                out[ds] = float(np.median(v))
        except Exception:  # noqa: BLE001
            pass
    return out


def _const_per_profile(time, cmap, fallback):
    """Per-profile calibration constant from a date->constant map (Kalman has all
    days, so gaps are rare; fill with the series median, else the default)."""
    med = float(np.median(list(cmap.values()))) if cmap else fallback
    dates = np.char.replace(
        np.datetime_as_string(time.astype("datetime64[D]")).astype("U10"), "-", "")
    return np.array([cmap.get(d, med) for d in dates], dtype="float64")


def _load_l1_window(s, start, end):
    from calibration.io.l1_window import load_l1_window
    paths = [str(_l1_file(s["wmo"], s["ident"], d)) for d in _days(start, end)
             if _l1_file(s["wmo"], s["ident"], d).exists()]
    return load_l1_window(paths) if paths else None


def _cache_coverage_regression(key, sdir, cache_name, output_name):
    """Regression guard for the historic per-station sens/omb caches.

    The dashboard panels read the aggregated ``<key>_sens.csv`` / ``<key>_omb.csv``,
    which are produced by aggregating the incremental ``_sens_cache.npz`` /
    ``_omb_cache.npz`` over the full 2025-2026 window. Those caches are built offline
    (chunked monthly jobs). A daily run that processes a SINGLE day would, for a station
    whose cache is missing, create a brand-new 1-day cache and then aggregate+overwrite
    the rich historic output with a 1-day value -- silently regressing the panel.

    Return True (=> caller must SKIP, leaving the output untouched and logging a WARNING)
    only when BOTH hold:
      * the product cache file does NOT exist (so any update would start a fresh cache), AND
      * a NON-EMPTY historic output already exists for this station+product.
    A genuinely new station (neither cache nor output) returns False so it proceeds
    normally and starts accumulating its cache. Caches that exist already get appended to
    (de-dup by date) as before, so they also return False. No calibration math is touched.
    """
    cache_path = sdir / cache_name
    if cache_path.exists():
        return False  # cache present -> normal incremental append, no regression risk
    out_path = sdir / output_name
    try:
        # "non-empty historic output" = at least one data row beyond the CSV header
        has_history = out_path.exists() and sum(1 for _ in open(out_path, encoding="utf-8")) > 1
    except OSError:
        has_history = False
    if has_history:
        # Visible on the worker's stdout (logging is forced to CRITICAL in _process_stream).
        print(f"REGRESSION-GUARD: {key}: missing {cache_name} but {output_name} has history "
              f"-> SKIP (not overwriting historic output with a partial-window cache)",
              flush=True)
        return True
    return False  # brand-new station (no cache, no output) -> proceed normally


def _do_omb(s, start, end, kalman_rows):
    """Observation-minus-Background vs CAMS (operational + our-calibrated). Writes
    <key>_omb.png and a one-row <key>_omb.csv. Returns the summary row or None."""
    from calibration.cloud.calibration import INSTRUMENT_CAL_DEFAULT
    from calibration.io.cams import find_cams_file, _has_backscatter
    from calibration.omb.omb import compute_omb
    from calibration.omb.figures import plot_omb_station
    key, itype = _key(s), s["type"]
    # Regression guard: never let a partial-window daily run replace a non-empty historic
    # OmB output for a station whose cache hasn't been built yet (see helper docstring).
    if _cache_coverage_regression(key, OUT / key, "_omb_cache.npz", f"{key}_omb.csv"):
        return None
    method = SENS_OMB_METHOD.get(itype, "rayleigh")
    kmap = _kalman_map(kalman_rows, method)
    if not kmap:
        return None  # no calibrated C_L -> skip (don't fabricate 'ours' from the default constant)
    cams = find_cams_file(CAMS, end.strftime("%Y%m%d"))
    if cams is None or not _has_backscatter(cams):
        return None  # needs the 0.4 deg CAMS-with-backscatter download first
    data = _load_l1_window(s, start, end)
    if data is None:
        return None
    default = INSTRUMENT_CAL_DEFAULT.get(itype, 1.0)
    c_ours = _const_per_profile(data["time"], kmap, default)
    c_op = _const_per_profile(data["time"], _op_map(s, start, end), default)
    rcs = data["rcs"].astype("float64")
    res = compute_omb(
        time=data["time"], range_agl=data["range"],
        beta_sources={"op": rcs / c_op[:, None], "ours": rcs / c_ours[:, None]},
        station_lat=data["lat"], station_lon=data["lon"], station_alt=data["alt"],
        wavelength=data["wl"], cams_file=str(cams), instrument=itype,
        cloud_base_height=data["cbh"], abs_cs_lookup_table=str(WV_LUT),
    )
    sdir = OUT / key
    from calibration.incremental import omb_cache_update, omb_cache_aggregate
    omb_cache_update(sdir, res)
    res = omb_cache_aggregate(sdir)
    plot_omb_station(res, itype, sdir / f"{key}_omb.png",
                     title=f"{s.get('site', key)} ({s['wmo']}) — OmB "
                           f"{start:%Y-%m-%d}..{end:%Y-%m-%d}")
    sc = res.scalar
    wv = sc.get("ours_wv", {}).get("median_bias", float("nan"))
    row = dict(date_start=start.strftime("%Y%m%d"), date_end=end.strftime("%Y%m%d"),
               wavelength=f"{data['wl']:.1f}",
               median_bias_ours=f"{sc['ours']['median_bias']:.6e}",
               median_bias_op=f"{sc['op']['median_bias']:.6e}",
               median_bias_ours_wv=f"{wv:.6e}",
               rms_ours=f"{sc['ours']['rms']:.6e}", n_obs=sc["ours"]["n_obs"])
    _write_csv_atomic(sdir / f"{key}_omb.csv", OMB_FIELDS, [row])
    return row


def _do_sens(s, start, end, kalman_rows):
    """Per-day noise -> detection thresholds over the window. Writes <key>_sens.png
    and a one-row <key>_sens.csv (headline ICAO detection altitude)."""
    from calibration.cloud.calibration import INSTRUMENT_CAL_DEFAULT
    from calibration.sensitivity.network import (
        sensitivity_over_period, combine_sens_results, plot_sensitivity_station)
    key, itype = _key(s), s["type"]
    # Regression guard: never let a partial-window daily run replace a non-empty historic
    # sensitivity output for a station whose cache hasn't been built yet (see helper docstring).
    if _cache_coverage_regression(key, OUT / key, "_sens_cache.npz", f"{key}_sens.csv"):
        return None
    method = SENS_OMB_METHOD.get(itype, "rayleigh")
    kmap = _kalman_map(kalman_rows, method)
    if not kmap:
        return None  # no calibrated C_L -> skip (sensitivity scale depends on it)
    default = INSTRUMENT_CAL_DEFAULT.get(itype, 1.0)
    # Process MONTH BY MONTH and stitch: loading a whole multi-month window at once OOMs
    # (CL61 17 months > 64 GB even in float32, mostly the concatenate peak). Each month
    # (~1 GB) is loaded, reduced to its daily beta_min columns, then freed.
    parts = []
    m0 = datetime(start.year, start.month, 1)
    while m0 <= end:
        m_end = (datetime(m0.year + (m0.month == 12), (m0.month % 12) + 1, 1)
                 - timedelta(days=1))
        ws, we = max(m0, start), min(m_end, end)
        data = _load_l1_window(s, ws, we)
        if data is not None:
            c = _const_per_profile(data["time"], kmap, default).astype("float32")
            beta = (data["rcs"] / c[:, None]) * np.float32(1e6)  # Mm^-1 sr^-1
            parts.append(sensitivity_over_period(
                time=data["time"], beta=beta, range_agl=data["range"], cbh=data["cbh"],
                lat=data["lat"], lon=data["lon"], wavelength=data["wl"]))
            del data, beta
        m0 = m_end + timedelta(days=1)
    from calibration.incremental import sens_cache_update, sens_cache_aggregate
    sdir = OUT / key
    for _p in parts:
        if _p is not None and getattr(_p, "dates", None) is not None and _p.dates.size:
            sens_cache_update(sdir, _p)
    res = sens_cache_aggregate(sdir)
    if res is None:
        return None
    plot_sensitivity_station(res, itype, sdir / f"{key}_sens.png",
                             title=f"{s.get('site', key)} ({s['wmo']}) — Sensitivity "
                                   f"{start:%Y-%m-%d}..{end:%Y-%m-%d}")

    def _fa(v):
        return "" if (v is None or v != v) else f"{v:.0f}"

    a = res.icao_alt
    row = dict(date_start=start.strftime("%Y%m%d"), date_end=end.strftime("%Y%m%d"),
               wavelength=f"{res.wavelength:.1f}",
               icao_alt_200=_fa(a.get(200.0)), icao_alt_2000=_fa(a.get(2000.0)),
               icao_alt_4000=_fa(a.get(4000.0)),
               sigma_night_3000=f"{res.sigma_probe.get(3000, float('nan')):.6e}",
               n_days_night=res.n_days_night, n_days_day=res.n_days_day)
    _write_csv_atomic(sdir / f"{key}_sens.csv", SENS_FIELDS, [row])
    return row


# --- One instrument stream --------------------------------------------------
def _process_stream(payload):
    s, start, end, methods, hk_only, do_sens, do_omb, no_cal = payload
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.CRITICAL)
    key = _key(s)
    sdir = OUT / key
    sdir.mkdir(parents=True, exist_ok=True)

    def _write_hk():
        hk = _do_hk(s, start, end)
        if not hk:
            return
        hk += _preserve_existing_hk(sdir / f"{key}_hk.csv", start, end)
        hk.sort(key=lambda r: r["date"])
        _write_csv_atomic(sdir / f"{key}_hk.csv", HK_FIELDS, hk)

    # Backfill / refresh the monitoring CSV only (skip everything else) -- cheap.
    if hk_only:
        _write_hk()
        return key, s["type"], 0, 0

    rows = []
    n_ok = 0
    # --no-cal reuses the EXISTING calibration + Kalman (for a sens/omb-only pass);
    # otherwise (re)compute the per-night calibration, Kalman best estimate and HK.
    if not no_cal:
        if "rayleigh" in methods and s["type"] in RAYLEIGH_TYPES:
            rows += _do_rayleigh(s, start, end)
        if "cloud" in methods and s["type"] in CLOUD_TYPES:
            rows += _do_cloud(s, start, end)
        if rows:
            # keep prior rows outside the processed window (and other methods) so daily/partial
            # runs accumulate history instead of overwriting with just the processed dates
            rows += _preserve_existing_rows(sdir / f"{key}_cal.csv", methods, start, end)
            rows.sort(key=lambda r: (r["method"], r["date"]))
            _write_csv_atomic(sdir / f"{key}_cal.csv", CSV_FIELDS, rows)
            _write_csv_atomic(sdir / f"{key}_kalman.csv",
                              ["method", "date", "kalman", "kalman_std"], _kalman_rows(rows))
            _write_hk()   # monitoring panel: daily HK means for every processed day
            n_ok = sum(1 for r in rows if _is_success(r["flag"]))

    # OmB / sensitivity add-ons consume the Kalman C_L: the fresh one if we just
    # calibrated, else the existing <key>_kalman.csv (the --no-cal reuse path).
    if do_sens or do_omb:
        kalman_rows = _kalman_rows(rows) if rows else _read_kalman_csv(sdir / f"{key}_kalman.csv")
        try:
            if do_omb:
                _do_omb(s, start, end, kalman_rows)
            if do_sens:
                _do_sens(s, start, end, kalman_rows)
        except Exception as exc:  # noqa: BLE001 - an add-on failure must not lose the calibration
            print(f"{key}: sens/omb failed: {type(exc).__name__}: {exc}", flush=True)

    return key, s["type"], len(rows), n_ok


# --- Stream selection -------------------------------------------------------
def select(census, types, per_type, limit, start, end, ignore_coverage=False):
    s0, s1 = start.strftime("%Y%m%d"), end.strftime("%Y%m%d")
    # The census first/last coverage is a static snapshot; for a daily run (date past the snapshot) it
    # would wrongly drop streams. ignore_coverage selects every census stream of the type and lets the
    # per-stream L1-file check decide whether there is data for the day.
    streams = [s for s in census
               if s.get("type") in types
               and (ignore_coverage or (s.get("first", "") <= s1 and s.get("last", "") >= s0))]
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
    ap.add_argument("--ignore-coverage", action="store_true",
                    help="select every census stream of the type, ignoring its first/last coverage "
                         "snapshot (use for daily runs whose date is past the snapshot)")
    ap.add_argument("--hk-only", action="store_true",
                    help="only (re)extract the per-day housekeeping <key>_hk.csv (laser/optics/temperature "
                         "daily means) for the dashboard monitoring panel; skip calibration. Cheap backfill.")
    ap.add_argument("--omb", action="store_true",
                    help="also produce Observation-minus-Background vs CAMS (<key>_omb.png/.csv): "
                         "operational + our-calibrated backscatter. Needs a CAMS_Beta file WITH aerosol "
                         "backscatter for the window's month.")
    ap.add_argument("--sens", action="store_true",
                    help="also produce the instrument sensitivity / detection-threshold product "
                         "(<key>_sens.png/.csv): per-day noise -> ICAO detection altitude.")
    ap.add_argument("--no-cal", action="store_true",
                    help="skip (re)calibration and REUSE the existing <key>_kalman.csv for --omb/--sens "
                         "(the decoupled add-on pass / operational D-1 step). Requires --omb and/or --sens.")
    args = ap.parse_args()
    if args.no_cal and not (args.omb or args.sens):
        ap.error("--no-cal requires --omb and/or --sens (nothing to do otherwise)")

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
        _process_stream((s, start, end, methods, args.hk_only,
                         args.sens, args.omb, args.no_cal))
        return

    types = [t.strip() for t in args.types.split(",")]
    streams = select(census, types, args.per_type, args.limit, start, end,
                     ignore_coverage=args.ignore_coverage)
    # drop streams that have no work for the requested methods (e.g. CHM15k/Mini-MPL under --methods cloud);
    # --hk-only extracts housekeeping for EVERY instrument type, so keep all selected streams.
    if args.hk_only:
        relevant = set(ITYPE)
    elif args.no_cal:
        # sens/omb-only: every calibratable stream (it reuses its existing Kalman)
        relevant = RAYLEIGH_TYPES | CLOUD_TYPES
    else:
        relevant = (RAYLEIGH_TYPES if "rayleigh" in methods else set()) | (CLOUD_TYPES if "cloud" in methods else set())
    streams = [s for s in streams if s["type"] in relevant]
    OUT.mkdir(parents=True, exist_ok=True)
    # Per-run output marker drives resume (skip streams that already produced it). The marker
    # is the run's PRODUCT regardless of --no-cal, so adding --omb/--sens to a normal run does
    # not get skipped on a pre-existing _cal.csv. sens before omb: sensitivity has no CAMS
    # dependency, so its CSV is the reliable "this stream ran" marker when both are requested.
    out_marker = ("_hk.csv" if args.hk_only else
                  "_sens.csv" if args.sens else
                  "_omb.csv" if args.omb else
                  "_cal.csv")
    done = {p.parent.name for p in OUT.glob(f"*/*{out_marker}")}
    todo = streams if (args.force or args.hk_only or args.no_cal) else [s for s in streams if _key(s) not in done]
    print(f"window {args.start}..{args.end} | methods={','.join(sorted(methods))} | "
          f"{len(streams)} streams selected; {len(done)} with output; {len(todo)} to do"
          f"{' (forced)' if args.force else ''}; {args.workers} workers "
          f"(per-stream timeout {STREAM_TIMEOUT}s)", flush=True)

    def _run_stream(s):
        """Run ONE stream in a separate process; kill it if it hangs past the timeout."""
        key = _key(s)
        cmd = [sys.executable, str(Path(__file__).resolve()), "--stream", key,
               "--start", args.start, "--end", args.end, "--methods", args.methods]
        if args.hk_only:
            cmd.append("--hk-only")
        if args.omb:
            cmd.append("--omb")
        if args.sens:
            cmd.append("--sens")
        if args.no_cal:
            cmd.append("--no-cal")
        try:
            subprocess.run(cmd, timeout=STREAM_TIMEOUT, check=False,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            return key, "TIMEOUT (killed)"
        return key, ("ok" if (OUT / key / f"{key}{out_marker}").exists() else "no-output")

    # The thread pool only SUPERVISES the per-stream subprocesses (the heavy numpy/netCDF
    # work runs in the children). A hung child is killed at the timeout and the run goes on
    # -- unlike a ProcessPoolExecutor worker, which hangs the whole pool indefinitely.
    n_todo = len(todo)
    t0 = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_run_stream, s): s for s in todo}
        for k, fut in enumerate(as_completed(futs), 1):
            s = futs[fut]
            try:
                key, status = fut.result()
                msg = f"{key} ({s['type']}): {status}"
            except Exception as exc:  # noqa: BLE001
                msg = f"FAILED {_key(s)}: {type(exc).__name__}: {exc}"
            # ETA: linear extrapolation from the average per-stream time so far (refines as it runs)
            elapsed = time.monotonic() - t0
            remaining = (n_todo - k) * elapsed / k if k else 0.0
            eta_clock = (datetime.now() + timedelta(seconds=remaining)).strftime("%a %H:%M")
            print(f"[{k}/{n_todo}] {msg} | elapsed {_fmt_hms(elapsed)} | "
                  f"ETA {_fmt_hms(remaining)} (~{eta_clock})", flush=True)

    Path(OUT / "ALL_DONE.flag").write_text("done")
    print(f"L1_2026_DONE — {n_todo} streams in {_fmt_hms(time.monotonic() - t0)}", flush=True)


if __name__ == "__main__":
    main()
