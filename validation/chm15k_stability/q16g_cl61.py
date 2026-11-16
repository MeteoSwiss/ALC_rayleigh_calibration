# -*- coding: utf-8 -*-
"""q16g (task 3): the CL61 negative control, redone with the most reliable CL61 constants available locally.

1. Vintage check. The CL61 constants used by q15f come from common.NC_DIR (paper_v08 calib_dashboard/nc).
   Are they the OPERATIONAL 910.55 nm / sigma 0.08 ("l55s008") constants or the retired 910.74 / FWHM 1.0?
   Compare night-for-night with the l55s008 reruns (diag_v22_l55s008: Payerne C, Lindenberg C) and with the
   pre-switch local copy (C:/Users/hervo/Downloads, Payerne C).
2. Per site (Payerne, Lindenberg, Camborne, Aosta): annual amplitude / peak DOY of
   CL61 Rayleigh (method 0), CL61 cloud (method 1; liquid-cloud, no molecular reference at all),
   CHM15k raw and kappa_in-normalised on the SAME period, and the CHM15k/CL61 ratio on common nights
   (Rayleigh/Rayleigh: same night; CHM15k/cloud: CL61 cloud within +-1 day). Harmonic + trend when the span
   exceeds 400 d; monthly-median harmonics as a robustness check.
3. Clear-sky molecular-normalised signal X(7-8 km) of both instruments on the SAME CHM15k calibration windows
   (10-min cache: 4 CL61 cached) -> amplitude / phase of X_CHM, X_CL61 and of their ratio (an aerosol-free,
   calibration-code-free test).

Outputs: q16g_cl61_vintage.csv, q16g_cl61_sites.csv, q16g_cl61.png
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from common import NC_DIR, load_cache, load_constants  # noqa: E402
from q16b_insolation import Q16, W, load_kappas  # noqa: E402
from validation.paper.kalman_improve_alc import flag_outliers  # noqa: E402

EPOCH = datetime(1970, 1, 1)
SITES = {"Payerne": ("0-20000-0-06610", "A", "C", 491.0),
         "Lindenberg": ("0-20000-0-10393", "0", "C", 98.0),
         "Camborne": ("0-20000-0-03808", "A", "C", 88.0),
         "Aosta": ("0-380-5-1", "0", "B", 545.0)}
L55 = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_l55s008")
OLD = Path("C:/Users/hervo/Downloads")


def read_nc(path, method):
    out = {k: [] for k in ("t", "t0", "t1", "c", "zb", "zt", "ver")}
    for f in sorted(Path(path).parent.glob(Path(path).name)):
        with netCDF4.Dataset(f) as ds:
            m = np.ma.filled(ds["calibration_method"][:], -1) == method
            out["t"].append(np.ma.filled(ds["time"][:], np.nan)[m])
            out["t0"].append(np.ma.filled(ds["start_time"][:], np.nan)[m])
            out["t1"].append(np.ma.filled(ds["end_time"][:], np.nan)[m])
            out["c"].append(np.ma.filled(ds["lidar_constant"][:].astype("f8"), np.nan)[m])
            out["zb"].append(np.ma.filled(ds["calibration_bottom_height"][:].astype("f8"), np.nan)[m])
            out["zt"].append(np.ma.filled(ds["calibration_top_height"][:].astype("f8"), np.nan)[m])
            v = (np.ma.filled(ds["calibration_version"][:], -1) if "calibration_version" in ds.variables
                 else np.full(m.size, -1))
            out["ver"].append(np.asarray(v)[m])
    if not out["t"]:
        return None
    o = {k: np.concatenate(v) for k, v in out.items()}
    s = np.argsort(o["t"])
    o = {k: v[s] for k, v in o.items()}
    ok = np.isfinite(o["c"]) & (o["c"] > 0)
    return {k: v[ok] for k, v in o.items()}


def read_l55_csv(key, method):
    f = L55 / key / f"{key}_cal.csv"
    if not f.exists():
        return None
    t, c = [], []
    for r in csv.DictReader(open(f, encoding="utf-8")):
        if r["method"] != method or r["cal_value"] in ("", "nan"):
            continue
        if r.get("flag", "1") not in ("1", "1.0"):
            continue
        t.append((datetime.strptime(r["date"], "%Y%m%d") - EPOCH).days)
        c.append(float(r["cal_value"]))
    return (np.array(t, float), np.array(c, float)) if t else None


def match_days(ta, tb, tol=0.75):
    """Index pairs of nearest neighbours within tol days."""
    if len(ta) == 0 or len(tb) == 0:
        return np.array([], int), np.array([], int)
    j = np.clip(np.searchsorted(tb, ta), 1, len(tb) - 1)
    j = np.where(np.abs(tb[j - 1] - ta) < np.abs(tb[j] - ta), j - 1, j)
    ok = np.abs(tb[j] - ta) <= tol
    return np.where(ok)[0], j[ok]


def vintage():
    rows = []
    for site, (wmo, ich, icl, _) in SITES.items():
        key = f"{wmo}_{icl}"
        for meth, mname in ((0, "rayleigh"), (1, "cloud")):
            cur = read_nc(NC_DIR / f"ALC_calibration_{key}20*.nc", meth)
            if cur is None:
                continue
            row = dict(site=site, stream=key, method=mname, n=len(cur["t"]),
                       versions=";".join(str(v) for v in np.unique(cur["ver"])))
            l55 = read_l55_csv(key, mname)
            if l55 is not None:
                ia, ib = match_days(np.floor(cur["t"]), l55[0], 0.1)
                r = cur["c"][ia] / l55[1][ib]
                row.update(n_vs_l55=len(r), med_ratio_vs_l55=float(np.median(r)) if r.size else np.nan,
                           frac_within2pct_l55=float(np.mean(np.abs(r - 1) < 0.02)) if r.size else np.nan)
            old = read_nc(OLD / f"ALC_calibration_{key}20*.nc", meth)
            if old is not None:
                ia, ib = match_days(cur["t"], old["t"], 0.1)
                r = cur["c"][ia] / old["c"][ib]
                row.update(n_vs_old=len(r), med_ratio_vs_old=float(np.median(r)) if r.size else np.nan,
                           frac_within2pct_old=float(np.mean(np.abs(r - 1) < 0.02)) if r.size else np.nan)
            rows.append(row)
            print("vintage", row, flush=True)
    write_csv(rows, Q16 / "q16g_cl61_vintage.csv")


def doy_phase(t):
    return W * np.mod(t, 365.25)


def harm(t, y, trend=None, min_n=25):
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if len(y) < min_n or t.max() - t.min() < 200:
        return np.nan, np.nan, len(y)
    if trend is None:
        trend = t.max() - t.min() > 400
    cols = [np.ones_like(t), np.cos(doy_phase(t)), np.sin(doy_phase(t))]
    if trend:
        cols.append((t - t.mean()) / 365.25)
    c, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)
    p = c[1] + 1j * c[2]
    return 100 * abs(p), float(np.mod(np.angle(p) / W, 365.25)), len(y)


def monthly_harm(t, y):
    mon = np.array([(EPOCH + timedelta(days=float(x))).strftime("%Y%m") for x in t])
    tm, ym = [], []
    for mm in np.unique(mon):
        s = (mon == mm) & np.isfinite(y)
        if s.sum() >= 3:
            tm.append(np.median(t[s]))
            ym.append(np.median(y[s]))
    tm, ym = np.array(tm), np.array(ym)
    if len(tm) < 8:
        return np.nan, np.nan, len(tm)
    a, p, _ = harm(tm, ym, min_n=8)
    return a, p, len(tm)


def screened(t, c):
    lc = np.log(c)
    bad = flag_outliers(lc)
    return t[~bad], lc[~bad]


def x78(cache, t0, t1, alt, wavelength_scale=1.0):
    """Clear-sky rcs / molecular attenuated backscatter (US-1976 shape) in 7-8 km over [t0, t1] windows."""
    from q15e_noise_profile import molecular_att
    tilt = float(cache.get("tilt", 0.0))
    tilt = 0.0 if not np.isfinite(tilt) else tilt
    edges = cache["range_edges"]
    rmid = 0.5 * (edges[:-1] + edges[1:])
    fine = np.arange(5.0, 8000.0, 10.0)
    bm = molecular_att(alt, fine * np.cos(np.deg2rad(tilt)))
    mol = np.array([bm[(fine >= edges[k]) & (fine < edges[k + 1])].mean() for k in range(len(rmid))])
    band = (rmid >= 7000) & (rmid < 8000)
    tb = (cache["day"][:, None] + (np.arange(144)[None, :] + 0.5) / 144).ravel()
    o = np.argsort(tb)
    rcs = cache["rcs"].reshape(-1, len(rmid))[o]
    clr = (cache["n_clear"].ravel()[o] >= 0.9 * np.maximum(cache["n_all"].ravel()[o], 1)) & (cache["n_clear"].ravel()[o] > 0)
    tb = tb[o]
    out = np.full(len(t0), np.nan)
    i0 = np.searchsorted(tb, t0)
    i1 = np.searchsorted(tb, t1, side="right")
    for k in range(len(t0)):
        seg = rcs[i0[k]:i1[k]][clr[i0[k]:i1[k]]]
        if seg.shape[0] < 6:
            continue
        v = np.nanmean(np.nanmean(seg, axis=0)[band] / mol[band])
        out[k] = np.log(v) if v > 0 else np.nan
    return out


def sites():
    k_in, _, _ = load_kappas()
    rows, plot = [], {}
    for site, (wmo, ich, icl, alt) in SITES.items():
        ch = load_constants(wmo, ich)
        mod = max(set(ch["module"]), key=list(ch["module"]).count)
        mc = (ch["module"] == mod) & (ch["c"] > 0)
        ch = {k: v[mc] for k, v in ch.items()}
        good = ~flag_outliers(np.log(ch["c"]))
        ch = {k: v[good] for k, v in ch.items()}
        lch = np.log(ch["c"])
        kin = k_in.get(f"{wmo}_{ich}", np.nan)
        cache = load_cache(wmo, ich)
        from q1_temperature import night_hk
        hk = night_hk(cache, ch["t0"], ch["t1"])
        T = hk[:, list(cache["hk_vars"]).index("temp_int")]
        lchk = lch - (kin / 100.0) * (T - np.nanmedian(T)) if np.isfinite(kin) else np.full(len(lch), np.nan)
        clr = read_nc(NC_DIR / f"ALC_calibration_{wmo}_{icl}20*.nc", 0)
        cld = read_nc(NC_DIR / f"ALC_calibration_{wmo}_{icl}20*.nc", 1)
        row = dict(site=site, chm=f"{wmo}_{ich}", cl61=f"{wmo}_{icl}", module=mod, kappa_in=kin)
        series = {}
        for tag, d in (("cl61_ray", clr), ("cl61_cloud", cld)):
            if d is None or len(d["t"]) < 25:
                continue
            t, y = screened(0.5 * (d["t0"] + d["t1"]) if tag == "cl61_ray" else d["t"], d["c"])
            series[tag] = (t, y)
        tmin = min(v[0].min() for v in series.values()) if series else ch["t"].min()
        tmax = max(v[0].max() for v in series.values()) if series else ch["t"].max()
        tc = 0.5 * (ch["t0"] + ch["t1"])
        inper = (tc >= tmin - 15) & (tc <= tmax + 15)
        series["chm_raw"] = (tc[inper], lch[inper])
        series["chm_kin"] = (tc[inper], lchk[inper])
        series["chm_raw_all"] = (tc, lch)
        for tag, (t, y) in series.items():
            a, p, n = harm(t, y)
            am, pm, nm = monthly_harm(t, y)
            row.update({f"amp_{tag}": a, f"peak_{tag}": p, f"n_{tag}": n, f"span_{tag}": float(t.max() - t.min()),
                        f"amp_mon_{tag}": am, f"peak_mon_{tag}": pm, f"nmon_{tag}": nm})
        # ratios on common nights
        for cl in ("cl61_ray", "cl61_cloud"):
            if cl not in series:
                continue
            tl, yl = series[cl]
            for chtag in ("chm_raw", "chm_kin"):
                t, y = series[chtag]
                ia, ib = match_days(t, tl, 3.0 if cl == "cl61_cloud" else 0.5)
                if len(ia) < 25:
                    continue
                rr = y[ia] - yl[ib]
                a, p, n = harm(t[ia], rr)
                am, pm, nm = monthly_harm(t[ia], rr)
                row.update({f"amp_{chtag}_over_{cl}": a, f"peak_{chtag}_over_{cl}": p, f"n_{chtag}_over_{cl}": n,
                            f"amp_mon_{chtag}_over_{cl}": am, f"peak_mon_{chtag}_over_{cl}": pm})
                if chtag == "chm_raw":
                    plot.setdefault(site, {})[f"ratio_{cl}"] = (t[ia], rr)
        # clear-sky X(7-8 km) of both instruments on the CHM15k windows
        cl_cache = load_cache(wmo, icl)
        if cl_cache is not None:
            xa = x78(cache, ch["t0"], ch["t1"], alt)
            xb = x78(cl_cache, ch["t0"], ch["t1"], alt)
            ok = np.isfinite(xa) & np.isfinite(xb)
            for tag, y in (("x78_chm", xa[ok]), ("x78_cl61", xb[ok]), ("x78_ratio", xa[ok] - xb[ok])):
                y = y.copy()
                t = tc[ok]
                bad = flag_outliers(y)
                a, p, n = harm(t[~bad], y[~bad])
                am, pm, nm = monthly_harm(t[~bad], y[~bad])
                row.update({f"amp_{tag}": a, f"peak_{tag}": p, f"n_{tag}": n, f"amp_mon_{tag}": am, f"peak_mon_{tag}": pm})
                plot.setdefault(site, {})[tag] = (t[~bad], y[~bad])
        for tag in ("cl61_ray", "cl61_cloud", "chm_raw", "chm_kin"):
            if tag in series:
                plot.setdefault(site, {})[tag] = series[tag]
        rows.append(row)
        print({k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    write_csv(rows, Q16 / "q16g_cl61_sites.csv")
    figure(plot)


def figure(plot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(20, 8), sharex=True)
    styles = {"chm_raw": ("CHM15k raw ln C", "0.4"), "cl61_ray": ("CL61 Rayleigh", "C0"),
              "cl61_cloud": ("CL61 cloud", "C1"), "x78_chm": ("CHM15k X(7-8 km)", "0.4"),
              "x78_cl61": ("CL61 X(7-8 km)", "C0"), "x78_ratio": ("CHM15k/CL61 X(7-8 km)", "C3"),
              "ratio_cl61_ray": ("CHM15k/CL61 Rayleigh", "C3"), "ratio_cl61_cloud": ("CHM15k/CL61 cloud", "C4")}
    for j, (site, d) in enumerate(plot.items()):
        for i, keys in enumerate((("chm_raw", "cl61_ray", "cl61_cloud", "ratio_cl61_cloud"),
                                  ("x78_chm", "x78_cl61", "x78_ratio"))):
            ax = axs[i, j]
            for k in keys:
                if k not in d:
                    continue
                t, y = d[k]
                lab, c = styles[k]
                dt = [EPOCH + timedelta(days=float(x)) for x in t]
                y = 100 * (y - np.nanmedian(y))
                ax.plot(dt, y, ".", ms=2, color=c, alpha=0.35)
                # 30-day running median
                o = np.argsort(t)
                tt, yy = t[o], y[o]
                rm = np.array([np.nanmedian(yy[(tt >= x - 15) & (tt <= x + 15)]) for x in tt])
                ax.plot([EPOCH + timedelta(days=float(x)) for x in tt], rm, "-", color=c, lw=1.6, label=lab)
            ax.axhline(0, color="k", lw=0.5)
            ax.set_ylim(-60, 60)
            ax.set_title(f"{site}: " + ("constants" if i == 0 else "clear-sky X(7-8 km), CHM15k windows"), fontsize=9)
            if j == 0:
                ax.set_ylabel("ln ratio to own median [%]")
            ax.legend(fontsize=7, loc="lower left")
    fig.suptitle("q16g - CL61 negative control (operational 910.55 nm constants; 30-day running medians)", fontsize=10)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(Q16 / "q16g_cl61.png", dpi=120)


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


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    vintage()
    sites()
