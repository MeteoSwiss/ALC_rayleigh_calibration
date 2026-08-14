# -*- coding: utf-8 -*-
"""v2.0 calibration coefficients for the Payerne L1-vs-L2 dashboard.

Input is ONLY the v2.0 calibration NetCDFs (ALC_calibration_<wigos>_<ident><YYYY>.nc, the E-PROFILE
exchange format written by this pipeline). Nothing else -- in particular NOT the raw operational
`<key>_cal.csv`, which additionally carries flag=0.5 (degraded) nights; two of those are wild at
Payerne CHM15k (2026-06-20 3.38e12, 2026-06-21 4.13e12, ~6x the median) and the published NetCDF
drops them (45 CSV successes -> 43 NetCDF nights).

The NetCDFs are per YEAR and one year is not enough history: the 2026 file alone gives CHM15k just
13 Rayleigh nights. Every available year is read and concatenated (2025+2026 -> 43 CHM15k nights,
264 CL31, 75 CL61 -- matching the per-night traces on the live dashboard station pages).

Smoothing uses the SAME filter and parameters as the operational dashboard,
monitoring.kalman.kalman_best_estimate: median-normalised so one set of RELATIVE noise parameters
fits any instrument magnitude, 4 % day-to-day random-walk drift, 15 %/yr seasonal accumulation over
gaps, daily-median aggregation, rolling-IQR outlier rejection, measurement noise from rolling-mean
residuals, predict-only on gap days.

That choice matters at Payerne: the CL31 optical block was replaced on 2026-07-07 ~13:00 and its
constant steps ~2.8x (4.1e7 -> 9.1e7). With 4 % daily drift the operational filter absorbs the step
within days (4.11e7 on 1 Jul -> 8.46e7 by 20 Jul, exactly as the dashboard shows), so no manual
break-point is needed. The alternative bridge (calib_benchmark.kalman -> run_kalman_from_matlab.py)
estimates its measurement noise from residuals over the WHOLE series, so a step inflates that
variance and the filter crawls -- it reached only 3.6e7 by August, ~2.4x too small.

Writes <OUT>/<wmo>_<ident>_<calib>_L1.csv (time, C_daily, C_daily_std, C_kalman, C_kalman_std) so
intercompare.load_calib_series picks them up unchanged, plus a points.json for the charts.
"""
from __future__ import annotations
import csv
import json
import sys
import warnings
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")                # validation.paper, monitoring
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling modules (folder name is not a valid package name)
import numpy as np
from netCDF4 import Dataset

from monitoring.kalman import kalman_best_estimate

WMO = "0-20000-0-06610"
NC_DIR = Path(r"C:/Users/hervo/Downloads")
YEARS = (2025, 2026)
OUT = Path(r"C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/calib")
# eprof_v2.2 (the noise-aware molecular-window gates) for comparison. It changes the RAYLEIGH
# retrieval ONLY, so at Payerne it moves the CHM15k (A) and leaves the cloud-calibrated CL31 (B)
# and CL61 (C) bit-identical -- those two are copied across so both variants are complete series.
OUT_V22 = OUT.parent / "calib_v22"
RAW = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/calib_raw")
# Per-channel v2.2 source. The CHM15k comes from the availability corpus; the CL61's RAYLEIGH
# series had to be produced separately (run_payerne_v22.py) because the CL61 was only ever a
# cloud-calibrated reference in that study -- and it is a Rayleigh retrieval, so v2.2 moves it too.
V22_SRC = {
    "A":  Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/candidates"
               "/cand_N2.5_PAYERNE_CHM15k_A.json"),
    "Cr": RAW / "payerne_C_v2.2.json",
}

# ident -> (instrument type, the calib tag the channel key uses)
# The CL61 carries BOTH calibrations -- Rayleigh on clear nights and liquid-cloud on cloudy ones --
# and they are separate retrievals of the same constant, so they get separate series rather than
# being merged. (Merging was harmless for a single reference line but hides the disagreement between
# the two methods, and only the Rayleigh one responds to the v2.2 gate change.)
INSTR = {
    "A":  dict(itype="CHM15k", calib="rayleigh", methods=(0,)),
    "B":  dict(itype="CL31",   calib="cloud",    methods=(1,)),
    "C":  dict(itype="CL61",   calib="cloud",    methods=(1,)),
    "Cr": dict(itype="CL61",   calib="rayleigh", methods=(0,), ident="C"),
}
METHOD_NAME = {0: "Rayleigh", 1: "Liquid clouds", 2: "Ground-based lidar", 3: "Satellite lidar"}


def read_calib_nc(ident, year):
    """Per-night v2.0 constants from one yearly file: (dates, C, Cstd, method, history, path)."""
    f = NC_DIR / f"ALC_calibration_{WMO}_{ident}{year}.nc"
    with Dataset(f) as nc:
        t = np.ma.filled(nc.variables["time"][:].astype("f8"), np.nan)
        lc = np.ma.filled(nc.variables["lidar_constant"][:].astype("f8"), np.nan)
        u = np.ma.filled(nc.variables["lidar_constant_uncertainty"][:].astype("f8"), np.nan)
        m = np.ma.filled(nc.variables["calibration_method"][:], -1).astype(int).ravel()
        created = getattr(nc, "history", "")
    ok = np.isfinite(t) & np.isfinite(lc) & (lc > 0)
    dates = [datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t[ok]]
    order = np.argsort([d.toordinal() for d in dates])
    return ([dates[i] for i in order], lc[ok][order], u[ok][order], m[ok][order], created, f)


def all_nights(ident):
    """Every available year concatenated. A night carrying BOTH methods (CL61 does) keeps both --
    they are true lidar constants on the same scale, and the filter aggregates per day anyway."""
    dates, C, U, M, created, src = [], [], [], [], "", None
    per_year = []
    for year in YEARS:
        if not (NC_DIR / f"ALC_calibration_{WMO}_{ident}{year}.nc").exists():
            continue
        d, c, u, m, created, src = read_calib_nc(ident, year)
        per_year.append(f"{year}:{len(d)}")
        dates += d
        C.append(c); U.append(u); M.append(m)
    if not dates:
        raise FileNotFoundError(f"no ALC_calibration_{WMO}_{ident}<year>.nc in {NC_DIR}")
    print(f"     nights per year: {', '.join(per_year)}  -> {len(dates)} total")
    return dates, np.concatenate(C), np.concatenate(U), np.concatenate(M), created, src


def build(chan, outdir=OUT):
    """Smooth one channel's v2.0 constants with the operational filter; write the CSV.

    ``chan`` keys INSTR, which may map several channels onto one physical unit (the CL61's two
    calibration methods); ``methods`` selects which calibration_method codes belong to it.
    """
    spec = INSTR[chan]
    ident = spec.get("ident", chan)
    dates, C, Cstd, method, created, src = all_nights(ident)
    keep = np.isin(method, np.asarray(spec["methods"]))
    if keep.sum() < 5:
        print(f"  {chan} ({spec['itype']}, {spec['calib']}): only {int(keep.sum())} nights "
              f"of method {spec['methods']} -- skipped")
        return None
    dates = [d for d, k in zip(dates, keep) if k]
    C, Cstd, method = C[keep], Cstd[keep], method[keep]
    key = f"{WMO}_{ident}_{spec['calib']}"
    outdir.mkdir(parents=True, exist_ok=True)

    kt, ks, kstd = kalman_best_estimate(dates, C)
    if not len(kt):
        print(f"  {chan} ({spec['itype']}): too few nights ({len(C)}) for the operational filter")
        return None

    # per-night values on the daily grid, for the C_daily columns / the chart
    daily = defaultdict(list)
    for d, c in zip(dates, C):
        daily[d.date()].append(float(c))
    rows = []
    for t, v, s in zip(kt, ks, kstd):
        day = t.astype("datetime64[D]").astype(object)
        obs = daily.get(day, [])
        rows.append([str(day), (np.median(obs) if obs else ""),
                     (np.std(obs) if len(obs) > 1 else (0.0 if obs else "")), v, s])
    with open(outdir / f"{key}_L1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
        w.writerows(rows)
    print(f"  {chan} ({spec['itype']}, {spec['calib']}): {len(C)} nights {dates[0]:%Y-%m-%d}..{dates[-1]:%Y-%m-%d}  "
          f"median C_L={np.median(C):.4g}  -> Kalman {len(kt)} days "
          f"{str(kt[0])[:10]}..{str(kt[-1])[:10]}, {len(np.unique(np.round(ks, 6)))} distinct, "
          f"last {ks[-1]:.4g}")

    return dict(
        ident=ident, itype=spec["itype"], calib=spec["calib"], key=key,
        source=src.name, created=created,
        points=[dict(date=d.strftime("%Y-%m-%d"), value=float(c),
                     std=(float(s) if np.isfinite(s) else None),
                     method=METHOD_NAME.get(int(mm), str(mm)))
                for d, c, s, mm in zip(dates, C, Cstd, method)],
        kalman=dict(date=[str(t)[:10] for t in kt],
                    value=[float(v) for v in ks],
                    std=[float(s) if np.isfinite(s) else 0.0 for s in kstd]),
    )


def v22_nights(chan="A"):
    """Per-night eprof_v2.2 Rayleigh constants for the Payerne CHM15k, from the availability study.

    Record format is [flag, C_L, uncertainty, bottom, top, message]; only nights the QC accepted
    (flag 1 or 0.5) carry a constant. Returns (dates, C, uncertainty) or None if the study output
    is not present -- the dashboard then simply offers the v2.0 variant alone.
    """
    src = V22_SRC.get(chan)
    if src is None or not src.exists():
        return None
    rec = json.loads(src.read_text())
    d, c, u = [], [], []
    for day in sorted(rec):
        v = rec[day]
        if v[0] in (1.0, 0.5) and v[1] and v[1] > 0:
            d.append(datetime.strptime(day, "%Y%m%d"))
            c.append(float(v[1]))
            u.append(float(v[2]) if v[2] else np.nan)
    return (d, np.array(c), np.array(u)) if len(d) >= 5 else None


def build_v22_channel(chan, outdir=OUT_V22):
    """One channel's v2.2 series from its own per-night constants."""
    nights = v22_nights(chan)
    if nights is None:
        print(f"  {chan}: no eprof_v2.2 output ({V22_SRC.get(chan)}) -> skipped")
        return None
    spec = INSTR[chan]
    ident = spec.get("ident", chan)
    outdir.mkdir(parents=True, exist_ok=True)
    dates, C, Cstd = nights
    key = f"{WMO}_{ident}_{spec['calib']}"
    kt, ks, kstd = kalman_best_estimate(dates, C)
    daily = defaultdict(list)
    for d, c in zip(dates, C):
        daily[d.date()].append(float(c))
    rows = []
    for t, v, s in zip(kt, ks, kstd):
        day = t.astype("datetime64[D]").astype(object)
        obs = daily.get(day, [])
        rows.append([str(day), (np.median(obs) if obs else ""),
                     (np.std(obs) if len(obs) > 1 else (0.0 if obs else "")), v, s])
    with open(outdir / f"{key}_L1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time", "C_daily", "C_daily_std", "C_kalman", "C_kalman_std"])
        w.writerows(rows)
    print(f"  {chan} ({spec['itype']}, {spec['calib']}) v2.2: {len(C)} nights "
          f"{dates[0]:%Y-%m-%d}..{dates[-1]:%Y-%m-%d}  median C_L={np.median(C):.4g}  "
          f"-> Kalman {len(kt)} days, last {ks[-1]:.4g}")
    return dict(
        ident=ident, itype=spec["itype"], calib=spec["calib"], key=key,
        source=V22_SRC[chan].name, created="eprof_v2.2 / N2.5",
        points=[dict(date=d.strftime("%Y-%m-%d"), value=float(c),
                     std=(float(s) if np.isfinite(s) else None), method="Rayleigh")
                for d, c, s in zip(dates, C, Cstd)],
        kalman=dict(date=[str(t)[:10] for t in kt], value=[float(v) for v in ks],
                    std=[float(s) if np.isfinite(s) else 0.0 for s in kstd]),
    )


def main():
    print("== Smoothing the v2.0 calibration NetCDFs (operational dashboard filter) ==")
    out = {}
    for chan in ("A", "B", "C", "Cr"):
        r = build(chan)
        if r is not None:
            out[chan] = r
    (OUT / "points.json").write_text(json.dumps(out), encoding="utf-8")
    print(f"-> {OUT}")

    # v2.2 changes the RAYLEIGH retrieval only, so the two Rayleigh channels (CHM15k A, CL61 Cr)
    # are rebuilt from their own v2.2 runs and the cloud-calibrated ones are carried over unchanged.
    print("== eprof_v2.2 variant (Rayleigh only: CHM15k A and CL61 Rayleigh move) ==")
    import shutil
    out22, any_v22 = {}, False
    for chan in ("A", "B", "C", "Cr"):
        if INSTR[chan]["calib"] == "rayleigh":
            r = build_v22_channel(chan)
            if r is not None:
                out22[chan] = r
                any_v22 = True
                continue
        r0 = out.get(chan)
        if r0 is not None:                       # cloud channels: identical under both variants
            out22[chan] = r0
            src = OUT / f"{r0['key']}_L1.csv"
            if src.exists():
                OUT_V22.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, OUT_V22 / src.name)
    if any_v22:
        (OUT_V22 / "points.json").write_text(json.dumps(out22), encoding="utf-8")
        print(f"-> {OUT_V22}")
    return {"v2.0": out, "v2.2": out22 if any_v22 else {}}


if __name__ == "__main__":
    main()
