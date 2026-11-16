# -*- coding: utf-8 -*-
"""ref23 step 2: per-night tables from the harness pickles (no L1 is read here).

  rays.pkl      one row per (site, ident, date, variant v22|v23): flag, constant, tier, window and
                every v2.3 budget component (1 sigma, relative) -- the totals are recomputed from
                them in the analysis -- with the REVIEWED budget (runs_rayleigh_reviewed/:
                window-mean noise, vertical cap 15 %, molecular climatology ASL)
  rays_prereview.pkl  the same with the pre-review budget (the first pass, runs/)
  clouds.pkl    one row per (site, ident, UTC day): the O'Connor cloud constant, its reported
                uncertainty (C_L x std/median of the per-profile coefficients), n, the robust
                spread and the number of independent cloud scenes of the day
  cloud_profiles.pkl  per-profile coefficients / times / cloud bases of every successful day
  published.pkl the published per-night constants (dashboard NetCDFs fetched 2026-10-09):
                Rayleigh (v2.2, version 220, 2 sigma ensemble) and cloud (version 1000, 1 sigma)
  ingredients.pkl   the per-night inter-comparison ingredients (hourly rcs, admission masks,
                CAMS corrections), keyed (site, date)

Usage: python validation/v23_reference/ref23_tables.py
"""
from __future__ import annotations

import pickle
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402

TAB = RC.OUT / "tables"
RAY_FIELDS = ["flag", "c", "c_all", "tier", "zb", "zt", "u_ens", "u_total", "sig_ens", "sig_vert",
              "sig_mol", "sig_aer", "sig_noise", "g", "g_sigma", "resid", "chi2", "mol_bias",
              "n_sub", "z_lo", "z_hi", "exceeds", "budget_model", "n_prof", "t_first", "t_last",
              "t_mean", "msg", "window_msg", "noise_estimator", "mol_href", "vertical_cap",
              "sig_vert_uncapped"]
CLOUD_FIELDS = ["flag", "c", "u", "n", "coef", "std", "cc_applied", "wl", "fwhm", "msg"]


def load_records(root=None):
    recs = []
    for f in sorted((root or RC.RUN_DIR).glob("*/*.pkl")):
        with open(f, "rb") as fh:
            recs.extend(pickle.load(fh))
    return recs


def n_scenes(t_days, gap_s=600.0):
    """Number of separate cloud scenes: runs of accepted profiles separated by > gap_s."""
    t = np.sort(np.asarray(t_days, float))
    if t.size == 0:
        return 0
    return int(1 + np.sum(np.diff(t) * 86400.0 > gap_s))


def build(recs):
    ray_rows, cloud_rows, cprof, ing = [], [], {}, {}
    fatal = []
    for r in recs:
        if "fatal" in r:
            fatal.append((r["site"], r["date"], r["fatal"]))
            continue
        site, ds = r["site"], r["date"]
        for ident, vv in r.get("ray", {}).items():
            for v, row in vv.items():
                d = dict(site=site, ident=ident, date=ds, variant=v)
                d.update({k: row.get(k, np.nan) for k in RAY_FIELDS})
                d["sub_z"] = row.get("sub_z")
                d["sub_c"] = row.get("sub_c")
                d["sub_s"] = row.get("sub_s")
                ray_rows.append(d)
        for ident, row in r.get("cloud", {}).items():
            d = dict(site=site, ident=ident, date=ds)
            d.update({k: row.get(k, np.nan) for k in CLOUD_FIELDS})
            vc = row.get("vc")
            if vc is not None and np.size(vc):
                vc = np.asarray(vc, float)
                lv = np.log(vc[vc > 0])
                d.update(p_med=float(np.median(vc)),
                         p_mad_rel=float(1.4826 * np.median(np.abs(lv - np.median(lv)))),
                         p_iqr_rel=float(np.subtract(*np.percentile(lv, [75, 25]))),
                         n_scenes=n_scenes(row.get("vt", [])),
                         cbh_med=float(np.nanmedian(row.get("vcbh", [np.nan]))),
                         t_med=float(np.median(row.get("vt", [np.nan]))))
                cprof[(site, ident, ds)] = dict(vc=row.get("vc"), vt=row.get("vt"),
                                                vcbh=row.get("vcbh"))
            cloud_rows.append(d)
        if r.get("ing") is not None:
            ing[(site, ds)] = r["ing"]
    rays = pd.DataFrame(ray_rows)
    clouds = pd.DataFrame(cloud_rows)
    return rays, clouds, cprof, ing, fatal


def published():
    """Published per-night constants of every study stream (both methods)."""
    import netCDF4
    rows = []
    cen = RC.census()
    for sk, s in RC.SITES.items():
        for ident in s["streams"]:
            lon = cen[(s["wmo"], ident)]["lon"]
            for y in (2025, 2026):
                f = RC.PUBLISHED_NC / f"ALC_calibration_{s['wmo']}_{ident}{y}.nc"
                if not f.exists():
                    continue
                with netCDF4.Dataset(f) as ds:
                    m = np.ma.filled(ds["calibration_method"][:], -1).astype(int)
                    t = np.ma.filled(ds["time"][:], np.nan)
                    t0 = np.ma.filled(ds["start_time"][:], np.nan)
                    c = np.ma.filled(ds["lidar_constant"][:].astype("f8"), np.nan)
                    u = np.ma.filled(ds["lidar_constant_uncertainty"][:].astype("f8"), np.nan)
                    zb = np.ma.filled(ds["calibration_bottom_height"][:].astype("f8"), np.nan)
                    zt = np.ma.filled(ds["calibration_top_height"][:].astype("f8"), np.nan)
                    ver = (np.ma.filled(ds["calibration_version"][:], -1).astype(int)
                           if "calibration_version" in ds.variables else np.full(m.size, -1))
                for k in range(m.size):
                    if m[k] == 0:
                        # night key: solar date of start_time + 12 h (q16a convention)
                        dd = np.floor(t0[k] + lon / 360.0 + 0.5)
                        method = "ray"
                    else:
                        dd = np.floor(t[k])            # cloud: UTC day, time = day + 0.5
                        method = "cloud"
                    ds_ = (RC.EPOCH + timedelta(days=float(dd))).strftime("%Y%m%d")
                    rows.append(dict(site=sk, ident=ident, method=method, date=ds_, c=c[k], u=u[k],
                                     zb=zb[k], zt=zt[k], version=ver[k]))
    pub = pd.DataFrame(rows)
    # a night key can repeat when a file was appended twice: keep the last row
    pub = pub.drop_duplicates(["site", "ident", "method", "date"], keep="last")
    return pub


def main():
    warnings.filterwarnings("ignore")
    TAB.mkdir(parents=True, exist_ok=True)
    recs = load_records()
    rays_pre, clouds, cprof, ing, fatal = build(recs)
    # the reviewed budget (Rayleigh-only re-run); the cloud / ingredient part is budget-free
    recs_rev = load_records(RC.RUN_DIR_REVIEWED)
    rays, _c, _p, _i, fatal_rev = build(recs_rev)
    fatal += fatal_rev
    pub = published()
    rays.to_pickle(TAB / "rays.pkl")
    rays_pre.to_pickle(TAB / "rays_prereview.pkl")
    clouds.to_pickle(TAB / "clouds.pkl")
    pub.to_pickle(TAB / "published.pkl")
    with open(TAB / "cloud_profiles.pkl", "wb") as f:
        pickle.dump(cprof, f, protocol=pickle.HIGHEST_PROTOCOL)
    with open(TAB / "ingredients.pkl", "wb") as f:
        pickle.dump(ing, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"records {len(recs)} (+{len(recs_rev)} reviewed), fatal {len(fatal)}, ray rows "
          f"{len(rays)} reviewed / {len(rays_pre)} pre-review, cloud rows {len(clouds)}, "
          f"published rows {len(pub)}, ingredient nights {len(ing)}")
    for f in fatal[:10]:
        print("FATAL", f)


if __name__ == "__main__":
    main()
