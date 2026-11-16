# -*- coding: utf-8 -*-
"""q16a - recalibration experiments with the real Rayleigh code (report 19 follow-up).

Questions
  (1) Does the annual cycle of the CHM15k Rayleigh constant C_L depend on the altitude of the
      calibration window?
  (2) Does remaining sunlight (twilight) contaminate the night calibration?

Method
  A harness around calibration.rayleigh.calibrate_rayleigh (the operational function, unmodified):
  every night is recalibrated from the native L1 (D:/E-PROFILE_L1_2026) with the operational options
  (options.json + molecular_method eprof_v2.2, which is what produced the published constants:
  calibration_version 220). The package is NEVER edited; three module-level names of
  calibration.rayleigh.calibration are monkeypatched inside this process only:
    * write_calibration_result -> no-op (nothing is written to any operational folder);
    * find_optimal_molecular_window -> for a FORCED window run, an OLS fit of the time-cleaned
      signal (same eprof_v2.2 contaminated-cell flagging) on p_mol over a fixed band [z1, z2] m AGL;
      the rest of the pipeline (Klett LR 52 sr below the window, the 125-member LR x shift x time
      ensemble, the median C_L, every QC gate) is the operational one;
    * filter_time_range -> for a SOLAR-DEPRESSION run, the operational SZA night selection followed
      by an extra SZA band / dusk-dawn restriction (same SZA formula as the operational code);
    * average_ceilometer_data -> pass-through that records the times of the cloud-screened profiles
      actually used (for mean SZA and mean T_int of each subset).

Runs per night
  op                operational (SZA > 100, free window)                      -> validation + reference
  f2_3 .. f6_8      forced windows 2-3, 3-4, 4-5, 5-6, 6-8 km AGL, free intercept -> question (1)
  z2_3 .. z6_8      the same windows, intercept forced to 0                     -> question (1)
  thr108            sza_night_threshold 108 instead of 100                      -> question (2)
  tw                SZA in [100, 106)  (sun 10-16 deg below the horizon)         -> question (2)
  dusk / dawn       the evening / morning halves of tw                          -> T_int confound check
  deep110           SZA > 110                                                    -> question (2)
  nt103 / core103   SZA in [100, 103) / SZA > 103 (exists also near the solstice at 52-55 N)
  Sub-night SZA runs use min_time_range = 1 h (it only matters on partly cloudy nights; 3 h would
  reject every twilight subset by construction).

Sub-commands
  validate  : reproduce the published nightly constants on 3 units x ~30 nights
  run       : the experiment on ~20 units, all nights with a published constant 2025-01..2026-07
  analyse   : tables + figures + summary inputs

Outputs only under C:/DATA/Projects/202610_CHM15k_stability/q16/q16a_*.
"""
from __future__ import annotations

import argparse
import copy
import csv
import dataclasses
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import OUT, L1_ROOT, chm15k_streams, load_cache, load_constants  # noqa: E402

Q16 = OUT / "q16"
RUN_DIR = Q16 / "q16a_runs"
EPOCH = datetime(1970, 1, 1)
CENSUS = REPO / "validation/scope_l1_2026_census.json"

# Forced calibration windows, m AGL (range; the CHM15k points vertically).
_B = {"2_3": (2000.0, 3000.0), "3_4": (3000.0, 4000.0), "4_5": (4000.0, 5000.0),
      "5_6": (5000.0, 6000.0), "6_8": (6000.0, 8000.0)}
# f*: free-intercept OLS (operational fit form, operational rel-error gate);
# z*: intercept forced to 0 (least squares through the origin, as E-PROF v0.25/MATLAB) - the
#     free intercept makes the slope of a fixed 1-km band very noisy at low SNR.
BANDS = {**{f"f{k}": (v, False) for k, v in _B.items()}, **{f"z{k}": (v, True) for k, v in _B.items()}}
# Solar-depression subsets: (sza_low_exclusive, sza_high_exclusive_or_None, part, min_time_range_h)
SZA_RUNS = {"thr108": (108.0, None, None, None),
            "tw": (100.0, 106.0, None, 1.0),
            "dusk": (100.0, 106.0, "dusk", 1.0),
            "dawn": (100.0, 106.0, "dawn", 1.0),
            "deep110": (110.0, None, None, 1.0),
            "nt103": (100.0, 103.0, None, 1.0),
            "core103": (103.0, None, None, 1.0)}
RUN_ORDER = ["op"] + list(BANDS) + list(SZA_RUNS)
SUCCESS = (1.0, 0.5)

# Units: spread over countries/latitudes, Schiphol A-D + Budapest A/B + Magurele A/B pairs,
# Lindenberg, southern (Granada, Evora, Potenza) and northern (Schleswig, Oslo) sites. Payerne A
# (abnormal receiver) is excluded.
UNITS = [
    ("0-20000-0-06240", "A"), ("0-20000-0-06240", "B"), ("0-20000-0-06240", "C"),
    ("0-20000-0-06240", "D"), ("0-20000-0-10393", "0"), ("0-20008-0-UGR", "A"),
    ("0-620-3704-5480", "B"), ("0-20008-0-POT", "A"), ("0-100-20000-0000001", "A"),
    ("0-20008-0-INO", "A"), ("0-20008-0-INO", "B"), ("0-380-5-1", "0"),
    ("0-20000-0-12843", "A"), ("0-20000-0-12843", "B"), ("0-20000-0-10962", "0"),
    ("0-250-1001-07151", "B"), ("0-20000-0-06736", "A"), ("0-20000-0-03590", "A"),
    ("0-20000-0-10035", "0"), ("0-20000-0-01492", "A"),
]
VALIDATION_UNITS = [("0-20000-0-10393", "0"), ("0-20000-0-06240", "A"), ("0-20008-0-UGR", "A")]

FIELDS = ["wmo", "ident", "date", "run", "flag", "c", "u", "zb_agl", "zt_agl", "n_native",
          "t_mean", "t_first", "t_last", "sza_mean", "sza_min", "sza_max", "tint_mean",
          "bg_mean", "snr_2_3", "snr_3_4", "snr_4_5", "snr_5_6", "snr_6_8", "snr_8_10",
          "slope", "r2", "rel_err", "cpu_s", "message"]

# --------------------------------------------------------------------------------------------
# Monkeypatching of calibration.rayleigh.calibration (this process only)
# --------------------------------------------------------------------------------------------
STATE = {"band": None, "sza": None, "capture": {}}
_PATCHED = {"done": False}


def _sza_and_solar_hours(times_days, lat, lon, date_str):
    """Solar zenith angle exactly as calibration.io.data_loader.filter_time_range computes it."""
    t = np.asarray(times_days, float)
    utc_h = (t - np.floor(t)) * 24.0
    solar_h = utc_h + lon / 15.0
    solar_h = np.mod(solar_h, 24.0)
    current = datetime.strptime(date_str, "%Y%m%d").date()
    doy = current.timetuple().tm_yday
    decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
    la = np.radians(lat)
    ha = np.radians((solar_h - 12.0) * 15.0)
    cos_sza = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)
    return np.degrees(np.arccos(np.clip(cos_sza, -1.0, 1.0))), solar_h


def _subset(data, mask):
    """Index every time-dimension field of a CeilometerData by a boolean mask."""
    n = len(data.time)
    kw = {}
    for f in dataclasses.fields(data):
        v = getattr(data, f.name)
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n and f.name not in (
                "range_alc", "altitude_grid"):
            kw[f.name] = v[mask]
        elif isinstance(v, list) and len(v) == n:
            kw[f.name] = [x for x, k in zip(v, mask) if k]
    return dataclasses.replace(data, **kw)


def _install_patches():
    if _PATCHED["done"]:
        return
    import calibration.rayleigh.calibration as RC
    from calibration.rayleigh import molecular_methods as MM
    from calibration.rayleigh.rayleigh_fit import RayleighFitResult
    from scipy.stats import linregress

    orig_find = RC.find_optimal_molecular_window
    orig_filter = RC.filter_time_range
    orig_avg = RC.average_ceilometer_data

    def no_write(*a, **k):
        return None

    def find(**kw):
        res = _find(**kw)
        STATE["capture"]["fit"] = (float(res.slope), float(res.r_squared), float(res.relative_error))
        return res

    def _find(**kw):
        band = STATE["band"]
        if band is None:
            return orig_find(**kw)
        (z1, z2), zero_icpt = band
        signal = np.asarray(kw["signal"], float)
        p_mol = np.asarray(kw["p_mol"], float)
        rng = np.asarray(kw["range_alc"], float)
        stack = kw.get("signal_stack")
        prm = MM.DEFAULT_PARAMS["eprof_v2.2"]
        clean = signal
        if stack is not None and np.ndim(stack) == 2 and np.shape(stack)[0] >= 5:
            # the same time-resolved contaminated-cell flagging as the operational eprof_v2.2
            flag = MM.flag_contaminated_cells(stack, p_mol, rng, nmad=prm["flag_nmad"],
                                              min_excess=prm["flag_min_excess"])
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                clean = np.nanmean(np.where(flag, np.nan, np.asarray(stack, float)), axis=0)
        sel = (rng >= z1) & (rng <= z2) & np.isfinite(clean) & np.isfinite(p_mol)
        nan_fit = RayleighFitResult(slope=np.nan, intercept=np.nan, r_squared=np.nan,
                                    std_error=np.nan, p_value=np.nan, center_range_m=np.nan,
                                    half_length_m=np.nan, range_start_m=np.nan,
                                    range_end_m=np.nan, altitude_start=0, altitude_end=0,
                                    relative_error=np.inf, search_diagnostics=None)
        if sel.sum() < 5:
            return nan_fit
        x, y = p_mol[sel], clean[sel]
        if zero_icpt:
            a = float(np.sum(x * y) / np.sum(x * x))
            b, p = 0.0, np.nan
            res = y - a * x
            ss_tot = float(np.sum((y - y.mean()) ** 2))
            r = np.sqrt(max(0.0, 1.0 - float(np.sum(res ** 2)) / ss_tot)) if ss_tot > 0 else 0.0
            se = float(np.sqrt(np.sum(res ** 2) / max(1, x.size - 1) / np.sum(x * x)))
        else:
            a, b, r, p, se = linregress(x, y)
        ratio = y / x
        med = np.median(ratio[np.isfinite(ratio)])
        rel = abs((a - med) / med * 100.0) if med != 0 else np.inf
        return RayleighFitResult(slope=a, intercept=b, r_squared=r ** 2, std_error=se, p_value=p,
                                 center_range_m=0.5 * (z1 + z2), half_length_m=0.5 * (z2 - z1),
                                 range_start_m=z1, range_end_m=z2, altitude_start=0,
                                 altitude_end=0, relative_error=rel, search_diagnostics=None)

    def filt(data, date_str, options):
        out = orig_filter(data, date_str, options)
        spec = STATE["sza"]
        if spec is None or len(out.time) == 0:
            return out
        lo, hi, part = spec
        sza, sh = _sza_and_solar_hours(out.time, out.latitude, out.longitude, date_str)
        keep = sza > lo
        if hi is not None:
            keep &= sza < hi
        if part == "dusk":
            keep &= sh >= 12.0
        elif part == "dawn":
            keep &= sh < 12.0
        return _subset(out, keep)

    def avg(data, *a, **k):
        STATE["capture"]["times"] = np.asarray(data.time, float).copy()
        STATE["capture"]["lat"] = float(data.latitude)
        STATE["capture"]["lon"] = float(data.longitude)
        return orig_avg(data, *a, **k)

    RC.write_calibration_result = no_write
    RC.find_optimal_molecular_window = find
    RC.filter_time_range = filt
    RC.average_ceilometer_data = avg
    _PATCHED["done"] = True


# --------------------------------------------------------------------------------------------
# Per-unit set-up
# --------------------------------------------------------------------------------------------
def _census_entry(wmo, ident):
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["wmo"] == wmo and s["ident"] == ident:
            return s
    raise KeyError(f"{wmo} {ident} not in census")


def _options():
    from calibration import CalibrationOptions, DataLevel
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT
    o.data_level = DataLevel.L1
    o.folder_output = Q16 / "q16a_nowrite"     # never written (write is patched out)
    o.plot_all = False
    o.plot_main = False
    o.molecular_method = "eprof_v2.2"          # calibration_version 220 of the published files
    return o


def _info(s):
    from calibration import InstrumentInfo
    from calibration.config import InstrumentType
    return InstrumentInfo(site_name=s.get("site", s["wmo"]), wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=InstrumentType.CHM15k, latitude=s["lat"],
                          longitude=s["lon"], altitude=s.get("alt", 0.0))


def night_dates_from_constants(con, lon):
    """Night key D (YYYYMMDD) of each published constant: the solar date of start_time + 12 h."""
    t0 = np.asarray(con["t0"], float)
    d = np.floor(t0 + lon / 360.0 + 0.5)
    return [(EPOCH + timedelta(days=float(x))).strftime("%Y%m%d") for x in d]


class HK:
    """10-min housekeeping cache lookup (temp_int, bckgrd_rcs_0) for a list of profile times."""

    def __init__(self, wmo, ident):
        c = load_cache(wmo, ident)
        self.ok = c is not None
        if not self.ok:
            return
        self.day = np.asarray(c["day"], float)
        names = [str(x) for x in c["hk_vars"]]
        self.tint = c["hk"][:, :, names.index("temp_int")].astype(float)
        self.bg = c["hk"][:, :, names.index("bckgrd_rcs_0")].astype(float)
        self.idx = {int(d): i for i, d in enumerate(self.day)}

    def mean(self, times):
        if not self.ok or times is None or len(times) == 0:
            return np.nan, np.nan
        t = np.asarray(times, float)
        di = np.floor(t).astype(int)
        bi = np.clip(((t - di) * 144).astype(int), 0, 143)
        rows = np.array([self.idx.get(int(d), -1) for d in di])
        ok = rows >= 0
        if not ok.any():
            return np.nan, np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return (float(np.nanmean(self.tint[rows[ok], bi[ok]])),
                    float(np.nanmean(self.bg[rows[ok], bi[ok]])))


def _snr_bands(fi):
    """Median per-bin SNR (night-mean signal / its photon noise) in fixed bands, op run only."""
    out = {}
    try:
        sig = np.asarray(fi["signal"], float)
        sg = np.asarray(fi["sigma_signal"], float)
        rng = np.asarray(fi["range_alc"], float)
    except (KeyError, TypeError):
        return out
    for name, (a, b) in {"snr_2_3": (2000, 3000), "snr_3_4": (3000, 4000), "snr_4_5": (4000, 5000),
                         "snr_5_6": (5000, 6000), "snr_6_8": (6000, 8000),
                         "snr_8_10": (8000, 10000)}.items():
        m = (rng >= a) & (rng < b) & np.isfinite(sig) & np.isfinite(sg) & (sg > 0)
        out[name] = float(np.median(sig[m] / sg[m])) if m.sum() >= 3 else np.nan
    return out


def calibrate_night(s, date_str, runs, hk, native=None):
    """Run the requested variants on one night; return a list of row dicts."""
    _install_patches()
    from calibration import calibrate_rayleigh
    from calibration.io.data_loader import build_file_paths, load_data
    info = _info(s)
    o_base = _options()
    rows = []
    t_load = time.process_time()
    if native is None:
        files = [f for f in build_file_paths(date_str, info, o_base) if f.exists()]
        if not files:
            return rows
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            native = load_data(files, info.instrument_type, o_base.data_level)
    if native is None:
        return rows
    load_cpu = time.process_time() - t_load
    op_flag = None
    for run in runs:
        if run != "op" and op_flag is not None and op_flag in (0, -1):
            break                                  # no data / not a clear night: nothing to vary
        o = copy.copy(o_base)
        STATE["band"] = BANDS.get(run)
        STATE["sza"] = None
        if run in SZA_RUNS:
            lo, hi, part, mtr = SZA_RUNS[run]
            o.sza_night_threshold = lo
            STATE["sza"] = (lo, hi, part)
            if mtr is not None:
                o.min_time_range = mtr
        STATE["capture"] = {}
        fi = {}
        t0 = time.process_time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = calibrate_rayleigh(date_str, info, o, preloaded_data=native,
                                       fit_inputs_out=fi)
            flag, c, u = float(r.flag), float(r.lidar_constant), float(r.uncertainty)
            zb = (r.calibration_bottom_height - info.altitude
                  if r.calibration_bottom_height is not None else np.nan)
            zt = (r.calibration_top_height - info.altitude
                  if r.calibration_top_height is not None else np.nan)
            msg = r.message
        except Exception as exc:  # noqa: BLE001 - one bad variant must not lose the night
            flag, c, u, zb, zt = -99.0, np.nan, np.nan, np.nan, np.nan
            msg = f"{type(exc).__name__}: {exc}"
        cpu = time.process_time() - t0 + (load_cpu if run == "op" else 0.0)
        if run == "op":
            op_flag = flag
        times = STATE["capture"].get("times")
        row = dict(wmo=s["wmo"], ident=s["ident"], date=date_str, run=run, flag=flag,
                   c=c if flag in SUCCESS else np.nan, u=u if flag in SUCCESS else np.nan,
                   zb_agl=zb, zt_agl=zt, cpu_s=round(cpu, 3), message=str(msg)[:120])
        if times is not None and len(times):
            sza, _ = _sza_and_solar_hours(times, info.latitude, info.longitude, date_str)
            row.update(n_native=len(times), t_mean=float(np.mean(times)), t_first=float(times.min()),
                       t_last=float(times.max()), sza_mean=float(np.mean(sza)),
                       sza_min=float(sza.min()), sza_max=float(sza.max()))
            tm, bg = hk.mean(times)
            row.update(tint_mean=tm, bg_mean=bg)
        if "fit" in STATE["capture"]:
            row.update(dict(zip(("slope", "r2", "rel_err"), STATE["capture"]["fit"])))
        if run == "op" and fi:
            row.update(_snr_bands(fi))
        rows.append(row)
    STATE["band"] = None
    STATE["sza"] = None
    return rows


def _write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    os.replace(tmp, path)


def work_chunk(wmo, ident, dates, runs, out_path):
    """Worker: one unit, a list of nights -> one CSV (resume = skip existing CSV)."""
    if out_path.exists():
        return str(out_path), 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    s = _census_entry(wmo, ident)
    hk = HK(wmo, ident)
    rows = []
    t0 = time.time()
    for d in dates:
        try:
            rows.extend(calibrate_night(s, d, runs, hk))
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(wmo=wmo, ident=ident, date=d, run="op", flag=-98,
                             message=f"{type(exc).__name__}: {exc}"[:120]))
    _write_rows(out_path, rows)
    return str(out_path), len(rows), time.time() - t0


def _pool_run(tasks, workers):
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(work_chunk, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(tasks)}] {Path(p).name}: {n} rows in {dt:.0f} s "
                  f"(elapsed {time.time() - t0:.0f} s)", flush=True)


# --------------------------------------------------------------------------------------------
# validate / run
# --------------------------------------------------------------------------------------------
def cmd_validate(args):
    tasks = []
    for wmo, ident in VALIDATION_UNITS:
        s = _census_entry(wmo, ident)
        con = load_constants(wmo, ident)
        dates = night_dates_from_constants(con, s["lon"])
        dates = sorted(set(dates))
        # 30 nights spread over the whole record (all seasons) + 10 consecutive calendar nights
        # (published success AND failure) to check the success/failure agreement too.
        pick = [dates[int(i)] for i in np.linspace(0, len(dates) - 1, 30)]
        start = datetime.strptime("20250901", "%Y%m%d")
        cal = [(start + timedelta(days=k)).strftime("%Y%m%d") for k in range(10)]
        for k, chunk in enumerate(np.array_split(sorted(set(pick) | set(cal)), 4)):
            tasks.append((wmo, ident, list(chunk), ["op"],
                          RUN_DIR.parent / "q16a_validate" / f"{wmo}_{ident}_{k}.csv"))
    _pool_run(tasks, args.workers)


def cmd_run(args):
    runs = RUN_ORDER if not args.runs else args.runs.split(",")
    tasks = []
    units = UNITS if not args.units else [tuple(u.split(":")) for u in args.units.split(",")]
    for wmo, ident in units:
        s = _census_entry(wmo, ident)
        con = load_constants(wmo, ident)
        if con is None:
            continue
        # ALL calendar nights 2025-01-01..2026-07-31 with an L1 file for day D (the variants are
        # skipped automatically when the operational run finds no data / no clear night).
        d0 = datetime(2025, 1, 1)
        dates = []
        for k in range((datetime(2026, 7, 31) - d0).days + 1):
            d = d0 + timedelta(days=k)
            f = L1_ROOT / wmo / f"{d:%Y}" / f"{d:%m}" / f"L1_{wmo}_{ident}{d:%Y%m%d}.nc"
            if f.exists():
                dates.append(f"{d:%Y%m%d}")
        if args.stride > 1:
            dates = dates[::args.stride]
        # one task per unit-month (load balance + resume granularity)
        months = sorted({d[:6] for d in dates})
        for m in months:
            ds = [d for d in dates if d[:6] == m]
            tasks.append((wmo, ident, ds, runs, RUN_DIR / f"{wmo}_{ident}" / f"{m}.csv"))
    if args.limit:
        tasks = tasks[:args.limit]
    print(f"{len(tasks)} tasks, {sum(len(t[2]) for t in tasks)} unit-nights", flush=True)
    _pool_run(tasks, args.workers)


# --------------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------------
W = 2 * np.pi / 365.25
ZRUNS = ["z2_3", "z3_4", "z4_5", "z5_6", "z6_8"]
FRUNS = ["f2_3", "f3_4", "f4_5", "f5_6", "f6_8"]
ZCENTRE = {"2_3": 2.5, "3_4": 3.5, "4_5": 4.5, "5_6": 5.5, "6_8": 7.0}
COMMON = ["op", "z3_4", "z4_5", "z5_6", "z6_8"]   # same-night set (2-3 km is aerosol-laden)


def _date_days(dates):
    return np.array([(datetime.strptime(str(d), "%Y%m%d") - EPOCH).days for d in dates], float)


def load_runs():
    import pandas as pd
    frames = []
    for f in sorted(RUN_DIR.glob("*/*.csv")):
        try:
            frames.append(pd.read_csv(f, dtype={"date": str, "ident": str}))
        except Exception:  # noqa: BLE001 - empty chunk
            continue
    df = pd.concat(frames, ignore_index=True)
    df["key"] = df["wmo"] + "_" + df["ident"]
    df["ok"] = df["flag"].isin(SUCCESS) & np.isfinite(df["c"]) & (df["c"] > 0)
    df["lnc"] = np.where(df["ok"], np.log(df["c"].where(df["c"] > 0)), np.nan)
    df["tday"] = _date_days(df["date"])
    return df


def main_module_nights(wmo, ident, dates):
    """Boolean per night: day D belongs to the unit's main optical module (from the 10-min cache)."""
    c = load_cache(wmo, ident)
    if c is None:
        return np.ones(len(dates), bool), ""
    mods = {int(d): str(m) for d, m in zip(c["day"], c["module"])}
    nm = [mods.get(int(x), "") for x in _date_days(dates)]
    vals, cnt = np.unique([m for m in nm if m], return_counts=True)
    if not len(vals):
        return np.ones(len(dates), bool), ""
    main = vals[np.argmax(cnt)]
    return np.array([(m == main) or (m == "") for m in nm]), main


def kalman_annual(tday, y):
    """Paper Kalman prediction model (q13 method): (amp %, peak DOY, scatter %) of C."""
    from q13_kalman_network import annual
    t = [EPOCH + timedelta(days=float(x)) for x in tday]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return annual(t, np.asarray(y, float))
        except Exception:  # noqa: BLE001 - curve_fit failure
            return np.nan, np.nan, np.nan


def harmonic(tday, y, peak_ref=None, clip=True):
    """OLS annual harmonic of y (ln units). Returns amp [%], peak DOY, in-phase amp vs peak_ref [%],
    residual sd [%], n. In-phase = component along cos(w (doy - peak_ref))."""
    from q1_temperature import robust_mask
    t = np.asarray(tday, float)
    y = np.asarray(y, float)
    m = robust_mask(y) if clip else np.isfinite(y)
    if m.sum() < 30:
        return dict(amp=np.nan, peak=np.nan, inph=np.nan, sd=np.nan, n=int(m.sum()))
    doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t[m]], float)
    X = np.column_stack([np.ones(m.sum()), np.cos(W * doy), np.sin(W * doy)])
    coef, *_ = np.linalg.lstsq(X, y[m], rcond=None)
    a, b = coef[1], coef[2]
    amp = np.hypot(a, b)
    peak = (np.arctan2(b, a) / W) % 365.25
    inph = np.nan
    if peak_ref is not None and np.isfinite(peak_ref):
        inph = amp * np.cos(W * (peak - peak_ref))
    sd = np.std(y[m] - X @ coef)
    return dict(amp=100 * amp, peak=peak, inph=100 * inph, sd=100 * sd, n=int(m.sum()))


def n2n_scatter(tday, y):
    """Night-to-night scatter [%]: robust sd of ln C around its 31-day running median."""
    from q1_temperature import running_median, robust_mask
    t = np.asarray(tday, float)
    y = np.asarray(y, float)
    o = np.argsort(t)
    t, y = t[o], y[o]
    m = robust_mask(y)
    rm = running_median(t[m], y[m], 31.0)
    an = y[m] - rm
    an = an[np.isfinite(an)]
    if an.size < 20:
        return np.nan
    return 100 * 1.4826 * np.median(np.abs(an - np.median(an)))


def _bootstrap_median(x, nboot=2000, seed=1):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    b = np.median(rng.choice(x, size=(nboot, x.size), replace=True), axis=1)
    return float(np.median(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def _q(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return "–"
    return f"{np.median(x):+.1f} [{np.percentile(x, 25):+.1f}, {np.percentile(x, 75):+.1f}]"


def _q_abs(x, fmt="{:.1f}"):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return "–"
    return (fmt.format(np.median(x)) + " [" + fmt.format(np.percentile(x, 25)) + ", "
            + fmt.format(np.percentile(x, 75)) + "]")


def bg_vs_sza(keys, census):
    """Night background (bckgrd_rcs_0, 10-min cache) vs solar zenith angle, normalised per night by
    the median of the blocks with SZA >= 110 of the same night. Returns a DataFrame of SZA-bin
    medians per unit (direct measurement of solar leakage into the background)."""
    import pandas as pd
    edges = np.arange(90.0, 126.0, 2.0)
    rows = []
    for key in keys:
        wmo, ident = key.rsplit("_", 1)
        s = census[(wmo, ident)]
        c = load_cache(wmo, ident)
        if c is None:
            continue
        names = [str(x) for x in c["hk_vars"]]
        bg = c["hk"][:, :, names.index("bckgrd_rcs_0")].astype(float).ravel()
        t = (np.asarray(c["day"], float)[:, None] + (np.arange(144)[None, :] + 0.5) / 144.0).ravel()
        lon, lat = s["lon"], s["lat"]
        # solar date of the night (noon-to-noon in solar time) and SZA, same formula as the code
        tsol = t + lon / 360.0
        night = np.floor(tsol + 0.5)
        sh = np.mod((tsol - np.floor(tsol)) * 24.0, 24.0)
        doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in night])
        decl = np.radians(-23.44) * np.cos(2.0 * np.pi * (doy + 10) / 365.25)
        la = np.radians(lat)
        cz = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(np.radians((sh - 12) * 15))
        sza = np.degrees(np.arccos(np.clip(cz, -1, 1)))
        ok = np.isfinite(bg) & (bg > 0)
        df = pd.DataFrame(dict(night=night[ok], sza=sza[ok], bg=bg[ok]))
        ref = df[df.sza >= 110].groupby("night")["bg"].median()
        df = df[df.night.isin(ref.index)].copy()
        df["r"] = df.bg.values / ref.loc[df.night].values
        df["bin"] = np.digitize(df.sza, edges)
        for b, g in df.groupby("bin"):
            if 1 <= b < len(edges) and len(g) >= 20:
                rows.append(dict(key=key, lat=lat, sza_lo=edges[b - 1], sza_hi=edges[b],
                                 ratio_median=float(np.median(g.r)), n=len(g)))
    return pd.DataFrame(rows)


def analyse():
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    warnings.filterwarnings("ignore")
    df = load_runs()
    census = {(s["wmo"], s["ident"]): s for s in json.load(open(CENSUS, encoding="utf-8"))}
    kin = {}
    for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8")):
        try:
            kin[r["stream"]] = float(r["b_T_ctrl"])
        except ValueError:
            pass
    out_lines = []

    # ---------------- 0. validation against the published constants ----------------
    val_rows = []
    keep_main = []
    for (wmo, ident), g in df[df.run == "op"].groupby(["wmo", "ident"]):
        s = census[(wmo, ident)]
        con = load_constants(wmo, ident)
        pub = dict(zip(night_dates_from_constants(con, s["lon"]), con["c"]))
        p = np.array([pub.get(d, np.nan) for d in g.date])
        both = g.ok.values & np.isfinite(p)
        r = g.c.values[both] / p[both]
        val_rows.append(dict(
            key=f"{wmo}_{ident}", site=s["site"], lat=s["lat"], n_nights=len(g),
            n_pub=int(np.isfinite(p).sum()), n_ours=int(g.ok.sum()), n_both=int(both.sum()),
            pub_only=int((np.isfinite(p) & ~g.ok.values).sum()),
            ours_only=int((~np.isfinite(p) & g.ok.values).sum()),
            ratio_median=float(np.median(r)) if r.size else np.nan,
            frac_identical=float(np.mean(np.abs(r - 1) < 1e-9)) if r.size else np.nan,
            ratio_p05=float(np.percentile(r, 5)) if r.size else np.nan,
            ratio_p95=float(np.percentile(r, 95)) if r.size else np.nan,
            max_abs_dev_pct=float(100 * np.max(np.abs(r - 1))) if r.size else np.nan))
    val = pd.DataFrame(val_rows)
    val.to_csv(Q16 / "q16a_validation.csv", index=False)

    # main-module night mask
    mm = {}
    for (wmo, ident), g in df.groupby(["wmo", "ident"]):
        dates = sorted(g.date.unique())
        m, main = main_module_nights(wmo, ident, dates)
        mm[f"{wmo}_{ident}"] = (dict(zip(dates, m)), main)
    df["main_mod"] = [mm[k][0].get(d, True) for k, d in zip(df.key, df.date)]
    d = df[df.main_mod].copy()

    # wide table: one row per unit-night, ln C per run
    wide = d.pivot_table(index=["key", "date"], columns="run", values="lnc", aggfunc="first")
    wide = wide.reset_index()
    wide["tday"] = _date_days(wide["date"])
    cpu = df.groupby("run")["cpu_s"].sum()

    # ---------------- 1. window test ----------------
    win_rows, ratio_rows = [], []
    snr_rows = []
    for key, g in wide.groupby("key"):
        wmo, ident = key.rsplit("_", 1)
        s = census[(wmo, ident)]
        opok = np.isfinite(g["op"].values)
        common = np.all(np.isfinite(g[COMMON].values), axis=1)
        h_op = harmonic(g.tday[opok], g["op"][opok])
        for run in ["op"] + ZRUNS + FRUNS:
            if run not in g:
                continue
            y = g[run].values
            for setname, msk in (("own", opok & np.isfinite(y)), ("common", common & np.isfinite(y))):
                if msk.sum() < 40:
                    continue
                ka, kp, ks = kalman_annual(g.tday.values[msk], np.exp(y[msk]))
                h = harmonic(g.tday.values[msk], y[msk], h_op["peak"])
                win_rows.append(dict(key=key, site=s["site"], lat=s["lat"], run=run, set=setname,
                                     n=int(msk.sum()), frac_of_op=msk.sum() / max(1, opok.sum()),
                                     amp_kalman=ka, peak_kalman=kp, scat_kalman=ks,
                                     amp_harm=h["amp"], peak_harm=h["peak"], inphase_harm=h["inph"],
                                     n2n=n2n_scatter(g.tday.values[msk], y[msk])))
        # window ratios on the common nights
        for a, b in (("z6_8", "z4_5"), ("z5_6", "z3_4"), ("z6_8", "z3_4"), ("z4_5", "z3_4"),
                     ("z3_4", "op"), ("z4_5", "op"), ("z5_6", "op"), ("z6_8", "op"),
                     ("z2_3", "op")):
            msk = common & np.isfinite(g[a].values) & np.isfinite(g[b].values)
            if msk.sum() < 40:
                continue
            lr = g[a].values[msk] - g[b].values[msk]
            h = harmonic(g.tday.values[msk], lr, h_op["peak"])
            ratio_rows.append(dict(key=key, site=s["site"], lat=s["lat"], pair=f"{a}/{b}",
                                   n=int(msk.sum()), median_pct=100 * np.median(lr),
                                   amp=h["amp"], peak=h["peak"], inphase=h["inph"],
                                   se_amp=h["sd"] * np.sqrt(2.0 / max(h["n"], 1)),
                                   n2n=n2n_scatter(g.tday.values[msk], lr)))
    win = pd.DataFrame(win_rows)
    rat = pd.DataFrame(ratio_rows)
    win.to_csv(Q16 / "q16a_windows_units.csv", index=False)
    rat.to_csv(Q16 / "q16a_window_ratios.csv", index=False)

    # SNR / availability per band from the op rows (per-bin SNR of the night-mean, 30 m bins)
    op = d[(d.run == "op") & d.ok].copy()
    op["month"] = op.date.str[4:6].astype(int)
    for band in ("snr_2_3", "snr_3_4", "snr_4_5", "snr_5_6", "snr_6_8", "snr_8_10"):
        for season, months in (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)), ("JJA", (6, 7, 8)),
                               ("SON", (9, 10, 11)), ("all", tuple(range(1, 13)))):
            x = op[op.month.isin(months)].groupby("key")[band].median()
            snr_rows.append(dict(band=band, season=season, median_unit_snr=float(np.nanmedian(x)),
                                 p25=float(np.nanpercentile(x, 25)), p75=float(np.nanpercentile(x, 75))))
    snr = pd.DataFrame(snr_rows)
    snr.to_csv(Q16 / "q16a_snr.csv", index=False)
    # success of the forced windows on op-success nights
    succ = []
    for run in ZRUNS + FRUNS:
        x = d[(d.run == run)]
        ops = d[(d.run == "op") & d.ok][["key", "date"]]
        mrg = ops.merge(x[["key", "date", "ok"]], on=["key", "date"], how="left")
        succ.append(dict(run=run, frac_success_on_op_nights=float(mrg.ok.fillna(False).mean())))
    succ = pd.DataFrame(succ)
    succ.to_csv(Q16 / "q16a_window_success.csv", index=False)

    # monthly network medians of ln(C_run / C_op) on the common nights (heatmap)
    wc = wide[np.all(np.isfinite(wide[COMMON].values), axis=1)].copy()
    wc["month"] = wc.date.str[4:6].astype(int)
    heat = np.full((len(ZRUNS), 12), np.nan)
    for i, run in enumerate(ZRUNS):
        x = wc[np.isfinite(wc[run])].copy()
        x["lr"] = x[run] - x["op"]
        # per unit monthly median, then network median (each unit weighted once)
        um = x.groupby(["key", "month"])["lr"].median().reset_index()
        for mth, gg in um.groupby("month"):
            if len(gg) >= 5:
                heat[i, mth - 1] = 100 * np.median(gg.lr)
    # monthly network medians of ln C per run, normalised by the unit's annual mean (cycle shape)
    cyc = np.full((len(COMMON) + 1, 12), np.nan)
    for i, run in enumerate(["op"] + ZRUNS[1:] + ["z2_3"]):
        x = wc[np.isfinite(wc[run])].copy()
        x["y"] = x[run] - x.groupby("key")[run].transform("median")
        um = x.groupby(["key", "month"])["y"].median().reset_index()
        for mth, gg in um.groupby("month"):
            if len(gg) >= 5:
                cyc[i, mth - 1] = 100 * np.median(gg.y)

    # ---------------- 2. solar-depression test ----------------
    sza_pairs = [("tw", "deep110"), ("nt103", "core103"), ("dusk", "deep110"), ("dawn", "deep110"),
                 ("dusk", "dawn"), ("thr108", "op"), ("deep110", "op"), ("tw", "op")]
    side = d.pivot_table(index=["key", "date"], columns="run",
                         values=["tint_mean", "bg_mean", "sza_mean", "zb_agl", "zt_agl", "n_native"],
                         aggfunc="first")
    tw_rows = []
    for a, b in sza_pairs:
        if a not in wide or b not in wide:
            continue
        msk = np.isfinite(wide[a].values) & np.isfinite(wide[b].values) & np.isfinite(wide["op"].values)
        sub = wide[msk]
        for _, r in sub.iterrows():
            ix = (r.key, r.date)
            sd = side.loc[ix]
            wmo, ident = r.key.rsplit("_", 1)
            s = census[(wmo, ident)]
            tw_rows.append(dict(pair=f"{a}/{b}", key=r.key, site=s["site"], lat=s["lat"],
                                date=r.date, tday=r.tday, month=int(r.date[4:6]),
                                lnr=r[a] - r[b], dT=sd[("tint_mean", a)] - sd[("tint_mean", b)],
                                bg_ratio=sd[("bg_mean", a)] / sd[("bg_mean", b)],
                                sza_a=sd[("sza_mean", a)], sza_b=sd[("sza_mean", b)],
                                dzb=sd[("zb_agl", a)] - sd[("zb_agl", b)],
                                n_a=sd[("n_native", a)], n_b=sd[("n_native", b)],
                                kappa_unit=kin.get(r.key, np.nan)))
    tw = pd.DataFrame(tw_rows)
    # Symmetric estimator: mean of ln(dusk/deep) and ln(dawn/deep) on nights where both halves
    # exist. Any drift that is linear in time through the night (T_int cooling, elapsed time,
    # aerosol evolution) cancels to first order; a sunlight excess (same sign at dusk and dawn)
    # does not.
    du = tw[tw.pair == "dusk/deep110"].set_index(["key", "date"])
    da = tw[tw.pair == "dawn/deep110"].set_index(["key", "date"])
    both = du.index.intersection(da.index)
    if len(both):
        sym = du.loc[both].copy()
        sym["lnr"] = 0.5 * (du.loc[both, "lnr"] + da.loc[both, "lnr"])
        sym["dT"] = 0.5 * (du.loc[both, "dT"] + da.loc[both, "dT"])
        sym["bg_ratio"] = 0.5 * (du.loc[both, "bg_ratio"] + da.loc[both, "bg_ratio"])
        sym["sza_a"] = 0.5 * (du.loc[both, "sza_a"] + da.loc[both, "sza_a"])
        sym["dzb"] = 0.5 * (du.loc[both, "dzb"] + da.loc[both, "dzb"])
        sym["pair"] = "sym(dusk,dawn)/deep110"
        tw = pd.concat([tw, sym.reset_index()], ignore_index=True)
    tw.to_csv(Q16 / "q16a_twilight_nights.csv", index=False)
    from q1_temperature import robust_mask
    tw_sum = []
    for pair, g in tw.groupby("pair"):
        y = 100 * g.lnr.values
        m = robust_mask(y) & np.isfinite(g.dT.values)
        X = np.column_stack([np.ones(m.sum()), g.dT.values[m]])
        coef, *_ = np.linalg.lstsq(X, y[m], rcond=None)
        res = y[m] - X @ coef
        cov = np.linalg.inv(X.T @ X) * (res @ res) / max(1, m.sum() - 2)
        ycorr = y - g.kappa_unit.fillna(1.5).values * g.dT.values      # unit intra-night kappa
        med, lo, hi = _bootstrap_median(y)
        cmed, clo, chi = _bootstrap_median(ycorr)
        # per-unit medians first, then network median (each unit counted once)
        per_unit = pd.Series(ycorr).groupby(g.key.values).median()
        tw_sum.append(dict(pair=pair, n_nights=len(g), n_units=g.key.nunique(),
                           median_pct=med, ci_lo=lo, ci_hi=hi,
                           median_Tcorr_pct=cmed, ci_lo_Tcorr=clo, ci_hi_Tcorr=chi,
                           unit_median_Tcorr=float(per_unit.median()),
                           intercept_dT0_pct=coef[0], se_intercept=float(np.sqrt(cov[0, 0])),
                           kappa_fit=coef[1], se_kappa=float(np.sqrt(cov[1, 1])),
                           dT_median=float(np.nanmedian(g.dT)), bg_ratio_median=float(np.nanmedian(g.bg_ratio)),
                           dzb_median=float(np.nanmedian(g.dzb)),
                           sza_a=float(np.nanmedian(g.sza_a)), sza_b=float(np.nanmedian(g.sza_b)),
                           robust_sd_pct=float(1.4826 * np.median(np.abs(y - np.median(y))))))
    tws = pd.DataFrame(tw_sum)
    tws.to_csv(Q16 / "q16a_twilight_summary.csv", index=False)
    # seasonal / latitude split of the T-corrected twilight ratios
    tw["ycorr"] = 100 * tw.lnr - tw.kappa_unit.fillna(1.5) * tw.dT
    tw["season"] = tw.month.map({12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
                                 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"})
    tw["latband"] = pd.cut(tw.lat, [0, 45, 50, 53, 90], labels=["<45", "45-50", "50-53", ">53"])
    split = (tw.groupby(["pair", "season"])["ycorr"].agg(["median", "count"]).reset_index())
    split2 = (tw.groupby(["pair", "latband"], observed=True)["ycorr"].agg(["median", "count"]).reset_index())
    split.to_csv(Q16 / "q16a_twilight_season.csv", index=False)
    split2.to_csv(Q16 / "q16a_twilight_latitude.csv", index=False)

    # annual cycle with sza_night_threshold 108 vs 100 (same nights and own nights)
    thr_rows = []
    for key, g in wide.groupby("key"):
        wmo, ident = key.rsplit("_", 1)
        s = census[(wmo, ident)]
        opok = np.isfinite(g["op"].values)
        t8 = np.isfinite(g["thr108"].values) if "thr108" in g else np.zeros(len(g), bool)
        both = opok & t8
        row = dict(key=key, site=s["site"], lat=s["lat"], n_op=int(opok.sum()), n_108=int(t8.sum()),
                   n_both=int(both.sum()))
        for name, run, msk in (("op_all", "op", opok), ("op_same", "op", both), ("thr108", "thr108", both)):
            if msk.sum() >= 40:
                ka, kp, ks = kalman_annual(g.tday.values[msk], np.exp(g[run].values[msk]))
                h = harmonic(g.tday.values[msk], g[run].values[msk])
                row.update({f"{name}_amp_kalman": ka, f"{name}_peak_kalman": kp,
                            f"{name}_amp_harm": h["amp"], f"{name}_peak_harm": h["peak"],
                            f"{name}_n2n": n2n_scatter(g.tday.values[msk], g[run].values[msk])})
        thr_rows.append(row)
    thr = pd.DataFrame(thr_rows)
    # annual harmonic of the SZA-subset ratios per unit (in phase with the unit's operational C)
    sr_rows = []
    peak_op = {k: harmonic(g.tday[np.isfinite(g["op"])], g["op"][np.isfinite(g["op"])])["peak"]
               for k, g in wide.groupby("key")}
    for (pair, key), g in tw.groupby(["pair", "key"]):
        y = g.lnr.values
        yc = g.ycorr.values / 100.0
        h = harmonic(g.tday.values, y, peak_op.get(key))
        hc = harmonic(g.tday.values, yc, peak_op.get(key))
        sr_rows.append(dict(pair=pair, key=key, n=h["n"], inphase=h["inph"], inphase_Tcorr=hc["inph"],
                            se_amp=h["sd"] * np.sqrt(2.0 / max(h["n"], 1))))
    srat = pd.DataFrame(sr_rows)
    srat.to_csv(Q16 / "q16a_twilight_ratio_harmonics.csv", index=False)
    thr.to_csv(Q16 / "q16a_thr108_units.csv", index=False)

    bgs = bg_vs_sza(sorted(df.key.unique()), census)
    bgs.to_csv(Q16 / "q16a_background_vs_sza.csv", index=False)
    bgnet = bgs.groupby("sza_lo")["ratio_median"].agg(
        median="median", p25=lambda v: np.percentile(v, 25), p75=lambda v: np.percentile(v, 75),
        n_units="count")

    # ---------------- figures ----------------
    fig, ax = plt.subplots(1, 4, figsize=(18, 5.5), sharey=True)
    ycen = {r: ZCENTRE[r[1:]] for r in ZRUNS}
    for j, (col, lab) in enumerate((("amp_kalman", "annual amplitude of C_L [%]"),
                                    ("peak_kalman", "peak DOY"),
                                    ("n2n", "night-to-night scatter [%]"))):
        a = ax[j]
        for run in ZRUNS:
            x = win[(win.run == run) & (win.set == "common")][col].values
            yy = np.full(x.size, ycen[run]) + np.random.default_rng(0).uniform(-0.12, 0.12, x.size)
            a.plot(x, yy, ".", color="0.6", ms=5)
            if x.size:
                a.errorbar(np.nanmedian(x), ycen[run],
                           xerr=[[np.nanmedian(x) - np.nanpercentile(x, 25)],
                                 [np.nanpercentile(x, 75) - np.nanmedian(x)]],
                           fmt="o", color="k", ms=8, capsize=3)
        x = win[(win.run == "op") & (win.set == "common")][col].values
        zop = (d[(d.run == "op") & d.ok].zb_agl.median() + d[(d.run == "op") & d.ok].zt_agl.median()) / 2e3
        if x.size:
            a.errorbar(np.nanmedian(x), zop, xerr=[[np.nanmedian(x) - np.nanpercentile(x, 25)],
                                                   [np.nanpercentile(x, 75) - np.nanmedian(x)]],
                       fmt="s", color="tab:red", ms=9, capsize=3, label="operational free window")
        a.set_xlabel(lab)
        a.grid(alpha=0.3)
    ax[0].set_ylabel("calibration window centre [km AGL]")
    ax[0].legend(loc="upper right", fontsize=8)
    a = ax[3]
    sx = snr[snr.season.isin(["DJF", "JJA", "all"])]
    for season, col in (("DJF", "tab:blue"), ("JJA", "tab:orange"), ("all", "k")):
        g = sx[sx.season == season]
        zz = [ZCENTRE[b[4:]] if b[4:] in ZCENTRE else 9.0 for b in g.band]
        a.plot(g.median_unit_snr, zz, "o-", color=col, label=season)
    a.axvline(1, color="0.5", ls="--")
    a.set_xscale("log")
    a.set_xlabel("per-bin SNR of the night mean (30 m bins)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    fig.suptitle("q16a - forced calibration windows (intercept 0), same nights, main optical module; "
                 f"{win[(win.run == 'op') & (win.set == 'common')].key.nunique()} CHM15k units")
    fig.tight_layout()
    fig.savefig(Q16 / "q16a_windows.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[1.2, 1]))
    a = ax[0]
    im = a.pcolormesh(np.arange(0.5, 13), [2, 3, 4, 5, 6, 8], heat, cmap="RdBu_r", vmin=-15, vmax=15,
                      shading="flat")
    plt.colorbar(im, ax=a, label="network median ln(C_window / C_op) [%]")
    a.set_xlabel("month")
    a.set_ylabel("forced window [km AGL]")
    a.set_title("window / operational, same nights")
    a = ax[1]
    labels = ["op"] + ZRUNS[1:] + ["z2_3"]
    for i, run in enumerate(labels):
        a.plot(np.arange(1, 13), cyc[i], "o-", label=run,
               color=("tab:red" if run == "op" else plt.cm.viridis(i / len(labels))))
    a.set_xlabel("month")
    a.set_ylabel("network median ln C - unit median [%]")
    a.set_title("seasonal shape of C_L per window (same nights)")
    a.grid(alpha=0.3)
    a.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(Q16 / "q16a_window_months.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 4, figsize=(22, 5.2))
    g = tw[tw.pair == "tw/deep110"]
    a = ax[0]
    a.plot(g.dT, 100 * g.lnr, ".", color="0.5", ms=4, alpha=0.6)
    rr = tws.set_index("pair")
    if "tw/deep110" in rr.index:
        xx = np.linspace(np.nanpercentile(g.dT, 1), np.nanpercentile(g.dT, 99), 10)
        a.plot(xx, rr.loc["tw/deep110", "intercept_dT0_pct"] + rr.loc["tw/deep110", "kappa_fit"] * xx,
               "r-", label=(f"fit: {rr.loc['tw/deep110', 'intercept_dT0_pct']:+.1f} % at dT=0, "
                            f"{rr.loc['tw/deep110', 'kappa_fit']:+.2f} %/K"))
    a.set_ylim(-60, 60)
    a.set_xlabel("T_int(twilight) - T_int(deep night) [K]")
    a.set_ylabel("ln(C_twilight / C_deep) [%]")
    a.set_title("SZA 100-106 vs SZA > 110, per night")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a = ax[1]
    for pair, col in (("tw/deep110", "k"), ("nt103/core103", "tab:blue"), ("dusk/deep110", "tab:orange"),
                      ("dawn/deep110", "tab:green"), ("sym(dusk,dawn)/deep110", "tab:red"),
                      ("thr108/op", "tab:purple")):
        gg = tw[tw.pair == pair]
        mo = gg.groupby("month")["ycorr"].median()
        cn = gg.groupby("month")["ycorr"].count()
        mo = mo[cn >= 15]
        a.plot(mo.index, mo.values, "o-", color=col, label=pair)
    a.axhline(0, color="0.5")
    a.set_xlabel("month")
    a.set_ylabel("median ln ratio, T_int-corrected [%]")
    a.set_title("twilight / deep-night ratio by month (>=15 nights)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a = ax[2]
    for pair, col in (("tw/deep110", "k"), ("nt103/core103", "tab:blue")):
        gg = tw[tw.pair == pair]
        um = gg.groupby("key").agg(lat=("lat", "first"), y=("ycorr", "median"), n=("ycorr", "count"))
        um = um[um.n >= 10]
        a.plot(um.lat, um.y, "o", color=col, label=pair)
    a.axhline(0, color="0.5")
    a.set_xlabel("latitude [deg N]")
    a.set_ylabel("unit median ln ratio, T_int-corrected [%]")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a = ax[3]
    for key, gg in bgs.groupby("key"):
        a.plot(gg.sza_lo + 1, gg.ratio_median, "-", color="0.75", lw=0.8)
    a.plot(bgnet.index + 1, bgnet["median"], "ko-", label="network median")
    a.axvspan(100, 106, color="tab:orange", alpha=0.15, label="twilight subset 100-106")
    a.axvline(100, color="tab:red", ls="--", label="operational threshold 100")
    a.set_ylim(0.995, 1.03)   # SZA 90-92 (x1.21) is off-scale on purpose: the twilight band is what matters
    a.set_xlabel("solar zenith angle [deg]")
    a.set_ylabel("background / deep-night background (SZA>=110)")
    a.set_title("measured solar leakage into bckgrd_rcs_0")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(Q16 / "q16a_twilight.png", dpi=130)
    plt.close(fig)

    # ---------------- console digest (used to write the summary) ----------------
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 40)
    print("VALIDATION\n", val.to_string())
    print("\nCPU seconds per run (all unit-nights):\n", cpu.to_string())
    print("\nWINDOWS (common nights) network median [IQR]:")
    for run in ["op"] + ZRUNS + FRUNS:
        for setname in ("common", "own"):
            x = win[(win.run == run) & (win.set == setname)]
            if len(x):
                print(f"{run:6s} {setname:6s} units {x.key.nunique():2d} n/unit {x.n.median():5.0f} "
                      f"amp_k {_q_abs(x.amp_kalman)} peak_k {_q_abs(x.peak_kalman, '{:.0f}')} "
                      f"scat_k {_q_abs(x.scat_kalman)} amp_h {_q_abs(x.amp_harm)} inph {_q(x.inphase_harm)} "
                      f"n2n {_q_abs(x.n2n)} frac {x.frac_of_op.median():.2f}")
    print("\nRATIOS (common nights):")
    for pair, x in rat.groupby("pair"):
        bm, blo, bhi = _bootstrap_median(x.inphase.values)
        print(f"{pair:10s} units {len(x):2d} median {_q(x.median_pct)} in-phase amp {_q(x.inphase)} "
              f"(95% CI of network median {blo:+.1f}..{bhi:+.1f}) amp {_q_abs(x.amp)} "
              f"noise-se {_q_abs(x.se_amp)} peak {_q_abs(x.peak, '{:.0f}')} n2n {_q_abs(x.n2n)}")
    print("\nSUCCESS\n", succ.to_string())
    print("\nSNR\n", snr.to_string())
    print("\nHEAT (month cols) window rows z2_3..z6_8\n", np.round(heat, 1))
    print("\nCYCLE rows", labels, "\n", np.round(cyc, 1))
    print("\nBACKGROUND vs SZA (network)\n", bgnet.to_string())
    print("\nTWILIGHT\n", tws.to_string())
    print("\nSZA-RATIO HARMONICS (per unit, in phase with C_op):")
    for pair, x in srat.groupby("pair"):
        x = x[x.n >= 30]
        if len(x) < 3:
            continue
        bm, blo, bhi = _bootstrap_median(x.inphase.values)
        cm, clo, chi = _bootstrap_median(x.inphase_Tcorr.values)
        print(f"{pair:24s} units {len(x):2d} in-phase {bm:+.2f} ({blo:+.2f}..{bhi:+.2f}) "
              f"T-corr {cm:+.2f} ({clo:+.2f}..{chi:+.2f}) noise-se {np.median(x.se_amp):.2f}")
    print("\nSEASON\n", split.to_string())
    print("\nLATITUDE\n", split2.to_string())
    print("\nTHR108\n", thr.to_string())
    for col in ("op_all_amp_kalman", "op_same_amp_kalman", "thr108_amp_kalman", "op_same_peak_kalman",
                "thr108_peak_kalman", "op_same_n2n", "thr108_n2n"):
        if col in thr:
            print(col, _q_abs(thr[col]))
    return df


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--workers", type=int, default=12)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=14)
    r.add_argument("--runs", default="")
    r.add_argument("--units", default="")
    r.add_argument("--stride", type=int, default=1)
    r.add_argument("--limit", type=int, default=0)
    t = sub.add_parser("time")
    t.add_argument("--unit", default="0-20000-0-10393:0")
    t.add_argument("--date", default="20250705")
    a = sub.add_parser("analyse")
    args = ap.parse_args()
    if args.cmd == "validate":
        cmd_validate(args)
    elif args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "time":
        wmo, ident = args.unit.split(":")
        s = _census_entry(wmo, ident)
        hk = HK(wmo, ident)
        t0 = time.time()
        rows = calibrate_night(s, args.date, RUN_ORDER, hk)
        for rw in rows:
            print({k: rw.get(k) for k in ("run", "flag", "c", "zb_agl", "zt_agl", "n_native",
                                          "sza_mean", "tint_mean", "snr_6_8", "cpu_s", "message")})
        print(f"wall {time.time() - t0:.1f} s")
    elif args.cmd == "analyse":
        analyse()


if __name__ == "__main__":
    main()
