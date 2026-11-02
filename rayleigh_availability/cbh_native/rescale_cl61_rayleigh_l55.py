"""Re-express the CL61 RAYLEIGH constants of calout_v22_04 at the correct WV line (910.55 / 0.188).

Companion of rescale_cl61_l55.py (cloud). In the Rayleigh fit the range-corrected signal is divided by
T2_wv(z) before the lidar-constant fit, so C_L carries 1/<T2_wv> over the fit window
[bottom_height, top_height] (stored per night). Each successful CL61 night is re-expressed with
    C_L(910.55/0.188) = C_L(910.74/1.0) * <T2_l74>_window / <T2_l55>_window
(same CAMS 0.4 deg / LUT / transmission function as the pipeline; window selection assumed unchanged),
then the per-stream Rayleigh/cloud ratio (phase4 section 1) is recomputed with the cloud constants of
cbh_native_l55/cbh_scenes.csv. Checked against the full l55 reruns diag_v22_l55s008.

Run:    python rayleigh_availability/cbh_native/rescale_cl61_rayleigh_l55.py
Output: C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_native_l55/cl61_rayleigh_l55.csv + printout.
"""
from __future__ import annotations

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
RUN = DATA / "calout_v22_04"
OUT = DATA / "cbh_native_l55"
CAMS = "D:/CAMS_Monthly_04;D:/CAMS_daily"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"


def window_t2(job):
    date, lat, lon, alt, zb, zt = job
    cams = find_cams_file(CAMS, date)
    if cams is None:
        return None
    t0 = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}T00:00:00")
    prof = cams_water_vapor_profile(cams, lat, lon, t0, t0 + np.timedelta64(24, "h"))
    if prof is None:
        return None
    rng = np.arange(0.0, zt + 60.0, 30.0)
    m = (rng >= zb) & (rng <= zt)
    out = []
    for lam0, fwhm in ((910.74, 1.0), (910.55, 0.188)):
        t2 = two_way_wv_transmission(alt + rng, alt, prof[0], prof[1], LUT, lam0, fwhm)
        out.append(float(np.nanmean(t2[m])))
    return out


def main():
    import pandas as pd
    census = {f"{s['wmo']}_{s['ident']}": s for s in
              json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))}
    cl61 = [k for k, s in census.items() if s.get("type") == "CL61" and (RUN / k / f"{k}_cal.csv").exists()]
    rows, jobs = [], []
    for k in cl61:
        d = pd.read_csv(RUN / k / f"{k}_cal.csv", dtype={"date": str})
        d = d[(d.method == "rayleigh") & d.flag.isin([1, 0.5])]
        c = census[k]
        for _, r in d.iterrows():
            if not (np.isfinite(r.bottom_height) and np.isfinite(r.top_height)):
                continue
            rows.append(dict(key=k, date=r.date, cal_value_l74=r.cal_value))
            jobs.append((r.date, float(c["lat"]), float(c["lon"]), float(c["alt"]),
                         float(r.bottom_height), float(r.top_height)))
    print(f"{len(jobs)} CL61 Rayleigh nights in {len(cl61)} streams", flush=True)
    with ProcessPoolExecutor(max_workers=16) as ex:
        res = list(ex.map(window_t2, jobs, chunksize=8))
    for row, r in zip(rows, res):
        row["t2_l74"], row["t2_l55"] = (r if r else (np.nan, np.nan))
        row["cal_value"] = row["cal_value_l74"] * row["t2_l74"] / row["t2_l55"]
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "cl61_rayleigh_l55.csv", index=False)
    rr = (df.cal_value / df.cal_value_l74).dropna()
    print(f"Rayleigh C_L(910.55)/C_L(910.74): median {rr.median():.4f} (P05-P95 {rr.quantile(.05):.4f}-{rr.quantile(.95):.4f})")
    # validation against full reruns
    for key in ("0-20000-0-06610_C", "0-20000-0-10393_C"):
        f = DATA / "diag_v22_l55s008" / key / f"{key}_cal.csv"
        full = pd.read_csv(f, dtype={"date": str})
        full = full[(full.method == "rayleigh") & full.flag.isin([1, 0.5])][["date", "cal_value"]]
        j = df[df.key == key].merge(full, on="date", suffixes=("_resc", "_full"))
        if len(j):
            print(f"  [validate] {key}: {len(j)} nights; rescaled/full {np.median(j.cal_value_resc / j.cal_value_full):.4f}; "
                  f"retired/full {np.median(j.cal_value_l74 / j.cal_value_full):.4f}")
    # per-stream Rayleigh/cloud ratio, old vs new (medians over each stream's record)
    sc = pd.read_csv(OUT / "cbh_scenes.csv", dtype={"date": str}, low_memory=False)
    sc = sc[sc.type == "CL61"]
    print("\nstream                 ray/cloud (910.74)  ray/cloud (910.55)  n_ray n_cloud")
    rat_old, rat_new = [], []
    for k in sorted(set(df.key)):
        a = df[df.key == k]
        b = sc[sc.key == k]
        if len(a) < 3 or len(b) < 3:
            continue
        ro = a.cal_value_l74.median() / b.cal_value_l74.median()
        rn = a.cal_value.median() / b.cal_value.median()
        rat_old.append(ro)
        rat_new.append(rn)
        print(f"{k:22s} {ro:18.3f} {rn:19.3f} {len(a):6d} {len(b):6d}")
    print(f"median over streams: {np.median(rat_old):.3f} (910.74) -> {np.median(rat_new):.3f} (910.55); "
          f"IQR new {np.percentile(rat_new, 25):.3f}-{np.percentile(rat_new, 75):.3f}")


if __name__ == "__main__":
    main()
