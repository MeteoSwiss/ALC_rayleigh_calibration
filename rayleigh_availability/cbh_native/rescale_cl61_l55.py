"""Re-express the CL61 cloud constants of calout_v22_04 at the correct WV line (910.55 nm / 0.188 nm).

Context: doc/reports/wv_cl61_91074_audit.md. The network run calout_v22_04 (2025-01-01 .. 2026-08-13)
calibrated the CL61 with the retired WV model 910.74 nm / FWHM 1.0 nm. The CBH-slope studies
(phase4 section 3, cbh_slope_heterogeneity, altitude_independence_audit, cbh_hopkin_study) are built
on the per-day scenes extracted from that run (cbh_native/cbh_scenes.csv).

Method (the same first-order replay as rayleigh_availability/wv_lambda_sweep.py): in the cloud
calibration the WV correction is applied before the integral (beta /= T2_wv), and over the ~300 m
cloud return T2_wv is smooth, so C_L carries a factor 1/T2_wv(CBH). Each CL61 stream-day is
therefore re-expressed with
    ln C_L(910.55/0.188) = ln C_L(910.74/1.0) - [ln T2_l55(CBH) - ln T2_l74(CBH)]
using the pipeline's own two_way_wv_transmission, the bundled LUT and CAMS 0.4 deg (monthly, then
daily) at the scene's median CBH. Scene acceptance is assumed unchanged. The approximation is checked
against the full l55 reruns where they exist (diag_v22_l55s008: Payerne C, Lindenberg C).

Run:    python rayleigh_availability/cbh_native/rescale_cl61_l55.py [--workers 16]
Output: C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_native_l55/cbh_scenes.csv (CL61 rows
        rescaled, other types copied unchanged, extra columns lnT2_l74 / lnT2_l55 / cal_value_l74),
        + validation printout. The original cbh_scenes.csv is not modified.
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from calibration.io.cams import find_cams_file                                  # noqa: E402
from calibration.water_vapor_correction.water_vapor import (                    # noqa: E402
    cams_water_vapor_profile, two_way_wv_transmission)

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
OUT_DIR = DATA / "cbh_native_l55"
CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
SPEC_L74 = (910.74, 1.0)    # retired model used by calout_v22_04
SPEC_L55 = (910.55, 0.188)  # correct CL61 line


def scene_lnT2(job):
    """(ln T2_l74, ln T2_l55) at the scene CBH, or None without CAMS."""
    date, lat, lon, alt, cbh = job
    cams = find_cams_file(CAMS, date)
    if cams is None:
        return None
    t0 = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}T00:00:00")
    prof = cams_water_vapor_profile(cams, lat, lon, t0, t0 + np.timedelta64(24, "h"))
    if prof is None:
        return None
    h_wv, n_wv = prof
    grid = alt + np.arange(0.0, cbh + 400.0, 30.0)
    i = int(np.argmin(np.abs((grid - alt) - cbh)))
    out = []
    for lam0, fwhm in (SPEC_L74, SPEC_L55):
        v = float(two_way_wv_transmission(grid, alt, h_wv, n_wv, LUT, lam0, fwhm)[i])
        out.append(np.log(v) if np.isfinite(v) and v > 0 else np.nan)
    return out


def validate(df):
    """Compare the rescaled C_L with the full l55 reruns (same days)."""
    import pandas as pd
    for key in ("0-20000-0-06610_C", "0-20000-0-10393_C"):
        f = DATA / "diag_v22_l55s008" / key / f"{key}_cal.csv"
        if not f.exists():
            print(f"  [validate] {key}: no full l55 rerun")
            continue
        full = pd.read_csv(f, dtype={"date": str})
        full = full[(full.method == "cloud") & full.flag.isin([1, 0.5])][["date", "cal_value"]]
        sub = df[(df.key == key)][["date", "cal_value", "cal_value_l74"]].copy()
        sub["date"] = sub["date"].astype(str)
        j = sub.merge(full, on="date", suffixes=("_resc", "_full"))
        if len(j):
            r_resc = j.cal_value_resc / j.cal_value_full
            r_old = j.cal_value_l74 / j.cal_value_full
            print(f"  [validate] {key}: {len(j)} common days; rescaled/full-l55 median "
                  f"{r_resc.median():.4f} (IQR {r_resc.quantile(.25):.4f}-{r_resc.quantile(.75):.4f}); "
                  f"retired(910.74)/full-l55 median {r_old.median():.4f}")


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    out_csv = OUT_DIR / "cbh_scenes.csv"
    if out_csv.exists():
        print(f"{out_csv} exists -> validation only")
        validate(pd.read_csv(out_csv, dtype={"date": str}, low_memory=False))
        return
    src = REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv"
    df = pd.read_csv(src, dtype={"date": str}, low_memory=False)
    census = {f"{s['wmo']}_{s['ident']}": s for s in
              json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))}
    df["cal_value_l74"] = pd.to_numeric(df["cal_value"], errors="coerce")
    df["lnT2_l74"] = np.nan
    df["lnT2_l55"] = np.nan
    cbh = pd.to_numeric(df["cbh_med"], errors="coerce")
    idx = df.index[(df["type"] == "CL61") & cbh.notna() & (df["cal_value_l74"] > 0)]
    jobs, keep = [], []
    for i in idx:
        c = census.get(df.at[i, "key"])
        if c is None:
            continue
        jobs.append((df.at[i, "date"], float(c["lat"]), float(c["lon"]), float(c["alt"]),
                     float(cbh[i])))
        keep.append(i)
    print(f"{len(jobs)} CL61 stream-days to re-express", flush=True)
    with ProcessPoolExecutor(max_workers=min(args.workers, 16)) as ex:
        res = list(ex.map(scene_lnT2, jobs, chunksize=8))
    n_ok = 0
    for i, r in zip(keep, res):
        if r is None or not np.all(np.isfinite(r)):
            continue
        df.at[i, "lnT2_l74"], df.at[i, "lnT2_l55"] = r
        df.at[i, "cal_value"] = df.at[i, "cal_value_l74"] * np.exp(-(r[1] - r[0]))
        n_ok += 1
    # CL61 days without a usable CAMS profile cannot be re-expressed: drop them (NaN constant)
    bad = df.index.isin(idx) & df["lnT2_l55"].isna()
    df.loc[bad, "cal_value"] = np.nan
    print(f"{n_ok} re-expressed, {int(bad.sum())} CL61 days dropped (no CAMS)", flush=True)
    m = df["lnT2_l55"].notna()
    ratio = np.exp(-(df.loc[m, "lnT2_l55"] - df.loc[m, "lnT2_l74"]))
    print(f"C_L(910.55)/C_L(910.74): median {ratio.median():.4f}, "
          f"P05-P95 {ratio.quantile(.05):.4f}-{ratio.quantile(.95):.4f}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"-> {out_csv}")
    validate(df)


if __name__ == "__main__":
    main()
