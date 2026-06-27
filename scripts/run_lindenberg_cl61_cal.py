"""
run_lindenberg_cl61_cal.py
==========================

Rayleigh (optimal molecular method) AND liquid-cloud calibration of the Vaisala
CL61 at Lindenberg, over the period 2024-05-01 .. 2026-05-31 (inclusive).

Both calibrations are run twice — *with* and *without* the 910 nm water-vapor
absorption correction — so the WV impact can be quantified directly.

The per-day calibration coefficients are then turned into a smooth daily best
estimate with the operational **E-PROFILE Kalman filter** (random-walk state +
seasonal predict model) taken from

    C:\\Users\\hervo\\OneDrive\\Documents\\Python\\improve_alc_calib
        src/improve_alc_calib/cal_best_estimate.py
            kalman_predict(), kalman_update(), constant()

— exactly the same routine used in the E-PROFILE reprocessing.

Raw data
--------
Daily Cloudnet raw files (~454 MB each):
    A:\\CL61_Cloudnet\\Lindenberg\\YYYYMMDD.nc
read through the package RAW reader (``DataLevel.RAW`` with
``folder_root = A:\\CL61_Cloudnet`` and ``wmo_id = "Lindenberg"``; the reader uses the
previous + current day to build each night).

Outputs (under ``OUT``)
-----------------------
- ``rayleigh_lindenberg_cl61.csv``  : per-night Rayleigh lidar_constant (wv on/off)
- ``cloud_lindenberg_cl61.csv``     : per-day cloud calibration coefficient (wv on/off)
- ``*_kalman.csv``                  : daily Kalman best estimate for each series
- ``lindenberg_cl61_rayleigh_diag.png``
- ``lindenberg_cl61_cloud_diag.png``
- ``lindenberg_cl61_wv_impact.png``

Run
---
    .venv\\Scripts\\python.exe run_lindenberg_cl61_cal.py

Re-running skips days already present in the CSVs unless FORCE_RECAL=1 is set in
the environment (the heavy per-day NetCDF reads are checkpointed to the CSVs).
"""
from __future__ import annotations

import os
import sys

# Keep BLAS single-threaded: the work is parallel across days, not within numpy.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import csv
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta, date
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from calibration import (
    calibrate_rayleigh,
    CalibrationOptions,
    InstrumentInfo,
    DataLevel,
)
from calibration.config import InstrumentType

# --- E-PROFILE operational Kalman filter -----------------------------------
# Imported from the separate improve_alc_calib project (NOT vendored: we use the
# exact operational routine so the best estimate matches E-PROFILE reprocessing).
_KALMAN_SRC = Path(r"C:/Users/hervo/OneDrive/Documents/Python/improve_alc_calib/src/improve_alc_calib")
if not (_KALMAN_SRC / "cal_best_estimate.py").is_file():
    raise FileNotFoundError(
        f"E-PROFILE Kalman source not found at {_KALMAN_SRC}. "
        "Update _KALMAN_SRC to point at improve_alc_calib/src/improve_alc_calib."
    )
sys.path.insert(0, str(_KALMAN_SRC))
from cal_best_estimate import kalman_predict, kalman_update, constant  # noqa: E402


# ===========================================================================
#  Configuration
# ===========================================================================
DATA_ROOT = Path("A:/CL61_Cloudnet")          # folder_root; site sub-folder = wmo_id
WMO = "Lindenberg"                              # sub-folder name under DATA_ROOT
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/lindenberg_cl61_cal")
OUT.mkdir(parents=True, exist_ok=True)

DATE_START = date(2024, 5, 1)
DATE_END = date(2026, 6, 15)                     # inclusive — full available archive

# Lindenberg (Meteorological Observatory) coordinates.
SITE = dict(label="Lindenberg_CL61", lat=52.21, lon=14.12, alt=123.0)

# WV correction needs CAMS monthly files + the HITRAN cross-section LUT.
CAMS_FOLDER = "D:/CAMS/"
ABS_CS_LUT = "C:/Users/hervo/OneDrive/Documents/MATLAB/MDA/monitoring_alc_monthly/abs_cross_647_full_levels_1000.nc"

FORCE = 0

# Parallelism: 30 CPUs available + fast disk. Each worker reads ~454 MB/day, so we
# leave a couple of cores free for I/O and the OS. Override with N_WORKERS env var.
N_WORKERS = int(os.environ.get("N_WORKERS", "28"))

# Kalman variance parameters.
# Both Rayleigh and Cloud are expressed in the Wiegner (2014) C_L convention
# (see below), so both are O(1–3) dimensionless — same process-noise scale.
#   var_eps_const   : long-term drift variance (random walk per day)
#   var_eps_temp    : seasonal amplitude / total days
KALMAN_PARAMS = {
    "cloud":   dict(var_eps_const=0.04 ** 2, var_eps_temp=(0.15 / 365.0) ** 2),
    "rayleigh": dict(var_eps_const=0.04 ** 2, var_eps_temp=(0.15 / 365.0) ** 2),
}


# ===========================================================================
#  Helpers
# ===========================================================================
def daterange(d0: date, d1: date):
    """Yield every calendar date from d0 to d1 inclusive."""
    d = d0
    while d <= d1:
        yield d
        d += timedelta(days=1)


def available_days() -> set[str]:
    """YYYYMMDD stems for which a non-empty daily raw file exists."""
    root = DATA_ROOT / WMO
    days = set()
    for p in root.glob("????????.nc"):
        try:
            if p.stat().st_size > 0:
                days.add(p.stem)
        except OSError:
            pass
    return days


def base_rayleigh_options(apply_wv: bool) -> CalibrationOptions:
    """Rayleigh options: RAW CL61, *optimal* molecular method, WV on/off."""
    o = CalibrationOptions()
    o.folder_root = DATA_ROOT
    o.data_level = DataLevel.RAW
    o.folder_output = OUT
    o.molecular_method = "eprof_v2"
    o.average_time_s = 60.0
    o.average_range_m = 30.0
    o.molecular_source = "standard"     # US Std 1976 (matches the MATLAB reference)
    o.use_std_atm = True
    o.use_sza_night = True
    o.apply_wv_correction = apply_wv
    o.cams_folder = Path(CAMS_FOLDER)
    o.abs_cs_lookup_table = Path(ABS_CS_LUT)
    o.plot_main = True
    return o


def make_instrument() -> InstrumentInfo:
    return InstrumentInfo(
        site_name=SITE["label"],
        wmo_id=WMO,
        identifier="",
        instrument_type=InstrumentType.CL61,
        latitude=SITE["lat"],
        longitude=SITE["lon"],
        altitude=SITE["alt"],
    )


# ---------------------------------------------------------------------------
#  CSV checkpoint I/O
# ---------------------------------------------------------------------------
def load_csv(path: Path) -> dict[str, dict]:
    """Read a checkpoint CSV keyed by the YYYYMMDD 'date' column."""
    rows: dict[str, dict] = {}
    if path.is_file():
        with path.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                rows[row["date"]] = row
    return rows


def write_csv(path: Path, header: list[str], rows: dict[str, dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


# ===========================================================================
#  Rayleigh calibration (optimal method, WV on + off)
# ===========================================================================
def _rayleigh_one(ds: str) -> dict:
    """Worker: calibrate one night (WV on + off). Runs in a separate process."""
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    info = make_instrument()
    rec = {"date": ds}
    for tag, apply_wv in (("wv", True), ("nowv", False)):
        opt = base_rayleigh_options(apply_wv=apply_wv)
        try:
            res = calibrate_rayleigh(ds, info, opt)
            rec[f"lc_{tag}"] = f"{res.lidar_constant:.6e}"
            rec[f"flag_{tag}"] = str(res.flag)
            rec[f"unc_{tag}"] = f"{res.uncertainty:.6e}"
        except Exception as exc:  # noqa: BLE001 - log + continue on bad nights
            rec[f"lc_{tag}"] = "nan"
            rec[f"flag_{tag}"] = "0"
            rec[f"unc_{tag}"] = "nan"
            logging.debug("Rayleigh %s (%s) failed: %s", ds, tag, exc)
    return rec


def run_rayleigh(days: set[str]) -> dict[str, dict]:
    """Per-night Rayleigh lidar_constant with the optimal method, WV on and off.

    The per-night calibrations are independent, so they run across a process pool
    (``N_WORKERS``). Results are checkpointed to CSV as workers complete.
    """
    csv_path = OUT / "rayleigh_lindenberg_cl61.csv"
    header = ["date", "lc_wv", "flag_wv", "unc_wv", "lc_nowv", "flag_nowv", "unc_nowv"]
    rows = {} if FORCE else load_csv(csv_path)

    nights = sorted(
        d.strftime("%Y%m%d")
        for d in daterange(DATE_START, DATE_END)
        if d.strftime("%Y%m%d") in days
        and (d - timedelta(days=1)).strftime("%Y%m%d") in days  # need prev day too
    )

    todo = [n for n in nights if n not in rows]
    print(f"[Rayleigh] {len(nights)} candidate nights, {len(todo)} to compute "
          f"({len(nights) - len(todo)} cached); {N_WORKERS} workers.")

    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(_rayleigh_one, ds): ds for ds in todo}
        for fut in as_completed(futures):
            ds = futures[fut]
            try:
                rows[ds] = fut.result()
            except Exception as exc:  # noqa: BLE001
                logging.debug("Rayleigh worker %s crashed: %s", ds, exc)
                continue
            done += 1
            if done % 20 == 0 or done == len(todo):
                print(f"  [Rayleigh] {done}/{len(todo)} done")
                write_csv(csv_path, header, rows)  # periodic checkpoint

    write_csv(csv_path, header, rows)
    return rows


# ===========================================================================
#  Cloud calibration (liquid-water clouds, WV on + off)
# ===========================================================================
def _load_cloud_csv_migrated(path: Path) -> dict[str, dict]:
    """Load the cloud CSV, converting old C-convention rows to C_L = 1/C on the fly.

    Old format used columns  c_wv / c_nowv  (C = beta_true / beta_measured, O'Connor
    multiplier).  New format uses  cl_wv / cl_nowv  (C_L = 1/C, Wiegner 2014 lidar
    constant).  Rows that already carry cl_* are passed through unchanged.
    """
    rows = load_csv(path)
    if not rows:
        return rows
    first = next(iter(rows.values()))
    if "cl_wv" in first:
        return rows   # already in new format
    # Migrate: convert C -> C_L = 1/C for every row
    migrated: dict[str, dict] = {}
    for ds, row in rows.items():
        new_row: dict[str, str] = {"date": ds}
        for tag in ("wv", "nowv"):
            c_str = row.get(f"c_{tag}", "nan")
            n_str = row.get(f"n_{tag}", "0")
            std_str = row.get(f"std_{tag}", "nan")
            try:
                c = float(c_str)
                n = float(n_str)
                std = float(std_str)
                if c > 0 and n > 0:
                    cl = 1.0 / c
                    cl_std = std / (c ** 2)
                    new_row[f"cl_{tag}"] = f"{cl:.6f}"
                    new_row[f"n_{tag}"] = str(int(n))
                    new_row[f"std_{tag}"] = f"{cl_std:.6f}"
                else:
                    new_row[f"cl_{tag}"] = "nan"
                    new_row[f"n_{tag}"] = "0"
                    new_row[f"std_{tag}"] = "nan"
            except (ValueError, ZeroDivisionError):
                new_row[f"cl_{tag}"] = "nan"
                new_row[f"n_{tag}"] = "0"
                new_row[f"std_{tag}"] = "nan"
        migrated[ds] = new_row
    print(f"  [Cloud] migrated {len(migrated)} rows from C -> C_L = 1/C convention")
    return migrated


def _cloud_one(ds: str) -> dict:
    """Worker: cloud-calibrate one day (WV on + off). Runs in a separate process."""
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    from calibration.cloud import (
        liquid_cloud_calibration,
        CloudCalConfig,
    )

    daily_nc = DATA_ROOT / WMO / f"{ds}.nc"
    rec = {"date": ds}
    for tag, wv in (("wv", True), ("nowv", False)):
        try:
            cfg = CloudCalConfig(
                nc_file=str(daily_nc),
                instrument="CL61",
                apply_wv_correction=wv,
                cams_folder=CAMS_FOLDER,
                abs_cs_lookup_table=ABS_CS_LUT,
                station_latitude=SITE["lat"],
                station_longitude=SITE["lon"],
                aerosol_lidar_ratio=50.0,
            )
            res = liquid_cloud_calibration(cfg)
            ok = res.n_profiles > 0 and np.isfinite(res.lidar_constant) and res.lidar_constant > 0
            # Wiegner & Geiss (2012) convention: report the absolute lidar constant
            # C_L = calibration_constant_0 / C, where C is the O'Connor cloud multiplier
            # (C = beta_true / beta_file). This puts cloud and Rayleigh C_L on the same axis.
            cl = res.lidar_constant if ok else float("nan")
            rec[f"cl_{tag}"] = f"{cl:.6f}" if ok else "nan"
            rec[f"n_{tag}"] = str(int(res.n_profiles))
            # std propagated through C_L = calConst / C:  std(C_L) = C_L * std(C) / median(C)
            cl_std = res.lidar_constant * res.cal_std / res.cal_median if ok else float("nan")
            rec[f"std_{tag}"] = f"{cl_std:.6f}" if ok else "nan"
        except Exception as exc:  # noqa: BLE001
            rec[f"cl_{tag}"] = "nan"
            rec[f"n_{tag}"] = "0"
            rec[f"std_{tag}"] = "nan"
            logging.debug("Cloud %s (%s) failed: %s", ds, tag, exc)
    return rec


def run_cloud(days: set[str]) -> dict[str, dict]:
    """Per-day liquid-cloud calibration coefficient, WV on and off.

    Independent per day, so calibrations run across a process pool (``N_WORKERS``),
    checkpointing to CSV as workers complete.
    """
    csv_path = OUT / "cloud_lindenberg_cl61.csv"
    header = ["date", "cl_wv", "n_wv", "std_wv", "cl_nowv", "n_nowv", "std_nowv"]
    rows = {} if FORCE else _load_cloud_csv_migrated(csv_path)

    all_days = sorted(
        d.strftime("%Y%m%d")
        for d in daterange(DATE_START, DATE_END)
        if d.strftime("%Y%m%d") in days
    )
    todo = [d for d in all_days if d not in rows]
    print(f"[Cloud] {len(all_days)} candidate days, {len(todo)} to compute "
          f"({len(all_days) - len(todo)} cached); {N_WORKERS} workers.")

    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(_cloud_one, ds): ds for ds in todo}
        for fut in as_completed(futures):
            ds = futures[fut]
            try:
                rows[ds] = fut.result()
            except Exception as exc:  # noqa: BLE001
                logging.debug("Cloud worker %s crashed: %s", ds, exc)
                continue
            done += 1
            if done % 20 == 0 or done == len(todo):
                print(f"  [Cloud] {done}/{len(todo)} done")
                write_csv(csv_path, header, rows)

    write_csv(csv_path, header, rows)
    return rows


# ===========================================================================
#  E-PROFILE Kalman best estimate
# ===========================================================================
def kalman_best_estimate(times, values, kind: str):
    """Run the E-PROFILE Kalman filter on an irregular (time, value) series.

    Mirrors run_kalman_from_matlab.py: daily-median aggregation, rolling-IQR
    outlier rejection, measurement-noise variance from rolling-mean residuals,
    predict-only on gap days. ``kind`` selects the process-noise scale
    ('cloud' or 'rayleigh').

    Returns (daily_grid_dates, kalman_state, kalman_std) on a contiguous daily grid.
    """
    params = KALMAN_PARAMS[kind]
    times = np.asarray(times)
    values = np.asarray(values, dtype=float)

    good = np.isfinite(values)
    times, values = times[good], values[good]
    if values.size < 5:
        return np.array([]), np.array([]), np.array([])

    order = np.argsort(times)
    times, values = times[order], values[order]

    # --- daily-median aggregation -----------------------------------------
    day_keys = np.array([t.date() for t in times])
    uniq_days = sorted(set(day_keys))
    daily_t = np.array([datetime(d.year, d.month, d.day) for d in uniq_days])
    daily_v = np.array([np.median(values[day_keys == d]) for d in uniq_days])

    n_days = daily_v.size
    if n_days < 5:
        return np.array([]), np.array([]), np.array([])

    # --- rolling-IQR outlier flag (window in days) ------------------------
    win = 30 if n_days > 100 else 10
    half = win // 2
    is_outlier = np.zeros(n_days, dtype=bool)
    for k in range(n_days):
        lo, hi = max(0, k - half), min(n_days, k + half + 1)
        w = daily_v[lo:hi]
        med = np.median(w)
        q25, q75 = np.percentile(w, [25, 75])
        iqr = q75 - q25
        if daily_v[k] < med - 1.5 * iqr or daily_v[k] > med + 1.5 * iqr:
            is_outlier[k] = True

    clean_t = daily_t[~is_outlier]
    clean_v = daily_v[~is_outlier]
    if clean_v.size < 5:
        clean_t, clean_v = daily_t, daily_v

    # --- predict model + measurement-noise variance ----------------------
    clean_t_list = [t.to_pydatetime() if hasattr(t, "to_pydatetime")
                    else datetime(t.year, t.month, t.day) for t in clean_t]
    predict_func, _, _ = constant(clean_t_list, clean_v)

    # rolling mean (10-day) residuals -> measurement variance
    roll = np.copy(clean_v)
    for k in range(clean_v.size):
        lo, hi = max(0, k - 5), min(clean_v.size, k + 6)
        roll[k] = np.mean(clean_v[lo:hi])
    var_meas = float(np.mean((clean_v - roll) ** 2))
    if not np.isfinite(var_meas) or var_meas <= 0:
        var_meas = float(np.var(clean_v)) or 1.0

    # --- contiguous daily grid -------------------------------------------
    grid = [clean_t_list[0] + timedelta(days=k)
            for k in range((clean_t_list[-1] - clean_t_list[0]).days + 1)]
    obs = {t.date(): v for t, v in zip(clean_t_list, clean_v)}

    # --- Kalman loop ------------------------------------------------------
    t_prev = clean_t_list[0] - timedelta(days=10)
    x_est = float(roll[0])
    var_est = float(var_meas)

    state, variance = [], []
    for t_cur in grid:
        y = obs.get(t_cur.date(), np.nan)
        if np.isfinite(y):
            x_a, var_a = kalman_predict(t_cur, t_prev, x_est, var_est,
                                        predict_func, **params)
            x_est, var_est = kalman_update(float(y), x_a, var_meas, var_a)
            t_prev = t_cur
            if not np.isfinite(x_est):
                x_est, var_est = x_a, var_a
        else:
            x_est, var_est = kalman_predict(t_cur, t_prev, x_est, var_est,
                                            predict_func, **params)
        state.append(x_est)
        variance.append(var_est)

    return (np.array(grid, dtype="datetime64[ns]"),
            np.array(state),
            np.sqrt(np.array(variance)))


# ---------------------------------------------------------------------------
#  Parsing helpers for the diagnostics
# ---------------------------------------------------------------------------
def _to_dt(ymd: str) -> datetime:
    return datetime.strptime(ymd, "%Y%m%d")


def _f(row: dict, key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return float("nan")


def extract_rayleigh(rows: dict[str, dict]):
    """-> dict tag -> (times, values) keeping only flag==1 / 0.5 nights."""
    out = {}
    for tag in ("wv", "nowv"):
        t, v = [], []
        for ds in sorted(rows):
            flag = _f(rows[ds], f"flag_{tag}")
            lc = _f(rows[ds], f"lc_{tag}")
            if flag in (1.0, 0.5) and np.isfinite(lc) and lc > 0:
                t.append(_to_dt(ds))
                v.append(lc)
        out[tag] = (np.array(t), np.array(v))
    return out


def extract_cloud(rows: dict[str, dict]):
    """Extract cloud C_L values (Wiegner 2014 convention: C_L = 1/C)."""
    out = {}
    for tag in ("wv", "nowv"):
        t, v = [], []
        for ds in sorted(rows):
            cl = _f(rows[ds], f"cl_{tag}")
            n = _f(rows[ds], f"n_{tag}")
            if np.isfinite(cl) and cl > 0 and n > 0:
                t.append(_to_dt(ds))
                v.append(cl)
        out[tag] = (np.array(t), np.array(v))
    return out


# ===========================================================================
#  Diagnostic plots
# ===========================================================================
_COL = {"wv": "#1f77b4", "nowv": "#d62728"}
_LBL = {"wv": "WV correction ON", "nowv": "WV correction OFF"}


def plot_series(series: dict, kind: str, ylabel: str, title: str, outfile: Path):
    """Scatter of daily values + Kalman best estimate, WV on vs off."""
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    kalman_curves = {}
    for tag in ("wv", "nowv"):
        t, v = series[tag]
        if t.size == 0:
            continue
        ax.plot(t, v, ".", color=_COL[tag], ms=4, alpha=0.35,
                label=f"{_LBL[tag]} (daily, N={v.size})")
        gt, gs, gstd = kalman_best_estimate(t, v, kind)
        if gt.size:
            kalman_curves[tag] = (gt, gs, gstd)
            ax.fill_between(gt, gs - gstd, gs + gstd, color=_COL[tag], alpha=0.15)
            ax.plot(gt, gs, "-", color=_COL[tag], lw=2.2,
                    label=f"{_LBL[tag]} Kalman (mean={np.nanmean(gs):.3g})")

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="best")

    # Bottom panel: relative WV impact on the Kalman estimate (on shared dates).
    if "wv" in kalman_curves and "nowv" in kalman_curves:
        gt_w, gs_w, _ = kalman_curves["wv"]
        gt_n, gs_n, _ = kalman_curves["nowv"]
        # align to common dates
        common = np.intersect1d(gt_w, gt_n)
        if common.size:
            sw = gs_w[np.isin(gt_w, common)]
            sn = gs_n[np.isin(gt_n, common)]
            rel = 100.0 * (sw - sn) / sn
            axr.plot(common, rel, "-", color="#444444", lw=1.5)
            axr.axhline(0, color="k", lw=0.8, ls=":")
            axr.set_ylabel("WV impact\n(ON-OFF) [%]")
            axr.text(0.01, 0.85, f"mean {np.nanmean(rel):+.2f}%",
                     transform=axr.transAxes, fontsize=8, va="top")
    axr.grid(True, alpha=0.3)
    axr.set_xlabel("Date")
    axr.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  saved {outfile}")
    return kalman_curves


def plot_wv_impact(ray_k: dict, cloud_k: dict, outfile: Path):
    """Summary: Rayleigh vs cloud, WV-on Kalman normalised to their own mean."""
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    panels = [
        ("Rayleigh — $C_L$", ray_k, axes[0]),
        ("Liquid-cloud — $C_L = 1/C$", cloud_k, axes[1]),
    ]
    for name, curves, ax in panels:
        for tag in ("wv", "nowv"):
            if tag in curves:
                gt, gs, _ = curves[tag]
                m = np.nanmean(gs)
                ax.plot(gt, gs / m, "-", color=_COL[tag], lw=2, label=_LBL[tag])
        ax.set_title(f"{name} — Kalman estimate (normalised)")
        ax.set_ylabel("$C_L$ / period mean")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(30)
            lbl.set_ha("right")

    fig.suptitle("Lindenberg CL61 — WV-correction impact on $C_L$ (Wiegner 2014 convention)")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  saved {outfile}")


def _plot_rayleigh_vs_cloud(ray_k: dict, cloud_k: dict, outfile: Path):
    """Overlay Rayleigh C_L and Cloud C_L (both WV-on) on the same axis.

    Both series should agree within ~10–20 % if the calibration is consistent
    (Rayleigh uses molecular reference at ~5–8 km; cloud uses liquid-water clouds
    at ~0.5–2.4 km).  The ratio C_L(Rayleigh)/C_L(Cloud) = C_Rayleigh × C_cloud
    should be ~1.
    """
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    colors = {"rayleigh": "#1f77b4", "cloud": "#2ca02c"}
    labels = {
        "rayleigh": r"Rayleigh $C_L$ (WV on)",
        "cloud":    r"Cloud $C_L = 1/C$ (WV on)",
    }
    kalmans = {}
    for key, curves in (("rayleigh", ray_k), ("cloud", cloud_k)):
        if "wv" not in curves:
            continue
        gt, gs, gstd = curves["wv"]
        if gt.size == 0:
            continue
        kalmans[key] = (gt, gs, gstd)
        ax.fill_between(gt, gs - gstd, gs + gstd, color=colors[key], alpha=0.12)
        ax.plot(gt, gs, "-", color=colors[key], lw=2.2, label=labels[key])

    ax.set_ylabel(r"$C_L$ (Wiegner 2014)  [—]")
    ax.set_title(
        r"Lindenberg CL61 — Rayleigh vs Cloud $C_L$ (WV on); both in Wiegner (2014) convention"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    if "rayleigh" in kalmans and "cloud" in kalmans:
        gt_r, gs_r, _ = kalmans["rayleigh"]
        gt_c, gs_c, _ = kalmans["cloud"]
        common = np.intersect1d(gt_r, gt_c)
        if common.size:
            sr = gs_r[np.isin(gt_r, common)]
            sc = gs_c[np.isin(gt_c, common)]
            ratio = sr / sc
            axr.plot(common, ratio, "-", color="#444444", lw=1.5)
            axr.axhline(1.0, color="k", lw=0.8, ls=":")
            axr.set_ylabel(r"$C_L^\mathrm{Ray} / C_L^\mathrm{Cloud}$")
            axr.text(0.01, 0.85, f"mean ratio = {np.nanmean(ratio):.3f}",
                     transform=axr.transAxes, fontsize=8, va="top")
    axr.grid(True, alpha=0.3)
    axr.set_xlabel("Date")
    axr.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  saved {outfile}")


def save_kalman_csv(curves: dict, outfile: Path):
    """Write the daily Kalman best estimate (both WV variants) to CSV."""
    dates = sorted({d for tag in curves for d in curves[tag][0]})
    if not dates:
        return
    idx = {tag: {curves[tag][0][i]: (curves[tag][1][i], curves[tag][2][i])
                 for i in range(curves[tag][0].size)} for tag in curves}
    with outfile.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["date", "kalman_wv", "kalman_std_wv", "kalman_nowv", "kalman_std_nowv"])
        for d in dates:
            kw, sw = idx.get("wv", {}).get(d, (np.nan, np.nan))
            kn, sn = idx.get("nowv", {}).get(d, (np.nan, np.nan))
            w.writerow([np.datetime_as_string(d, unit="D"),
                        f"{kw:.6e}", f"{sw:.6e}", f"{kn:.6e}", f"{sn:.6e}"])
    print(f"  saved {outfile}")


# ===========================================================================
#  Main
# ===========================================================================
def main():
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)

    print(f"Lindenberg CL61 calibration  {DATE_START} .. {DATE_END}")
    print(f"  data : {DATA_ROOT / WMO}")
    print(f"  out  : {OUT.resolve()}")

    days = available_days()
    in_range = {d.strftime("%Y%m%d") for d in daterange(DATE_START, DATE_END)} & days
    print(f"  {len(in_range)} daily raw files present in the requested range.")
    if not in_range:
        print("No data found in range — nothing to do.")
        return

    # --- Rayleigh (optimal) -----------------------------------------------
    ray_rows = run_rayleigh(days)
    ray_series = extract_rayleigh(ray_rows)

    # --- Cloud ------------------------------------------------------------
    cloud_rows = run_cloud(days)
    cloud_series = extract_cloud(cloud_rows)

    # --- Diagnostics ------------------------------------------------------
    print("Building diagnostic plots...")
    ray_k = plot_series(
        ray_series, "rayleigh",
        ylabel=r"Rayleigh lidar constant  $C_L$  (Wiegner 2014)  [—]",
        title="Lindenberg CL61 — Rayleigh calibration (improved method) + E-PROFILE Kalman",
        outfile=OUT / "lindenberg_cl61_rayleigh_diag.png",
    )
    cloud_k = plot_series(
        cloud_series, "cloud",
        ylabel=r"Cloud lidar constant  $C_L = 1/C$  (Wiegner 2014)  [—]",
        title="Lindenberg CL61 — Liquid-cloud calibration ($C_L = 1/C$) + E-PROFILE Kalman",
        outfile=OUT / "lindenberg_cl61_cloud_diag.png",
    )
    plot_wv_impact(ray_k, cloud_k, OUT / "lindenberg_cl61_wv_impact.png")
    _plot_rayleigh_vs_cloud(ray_k, cloud_k, OUT / "lindenberg_cl61_ray_vs_cloud.png")

    save_kalman_csv(ray_k, OUT / "rayleigh_lindenberg_cl61_kalman.csv")
    save_kalman_csv(cloud_k, OUT / "cloud_lindenberg_cl61_kalman.csv")

    print("DONE.")


if __name__ == "__main__":
    main()
