# -*- coding: utf-8 -*-
"""Where does the cloud calibration land for alpha = 10, 20, 30 /km (and a_G 3.3 / 3.7 um)?

Exact per-profile rescaling on the stored in-cloud windows of the 43 study units (22 CL31, 10 CL51, 11 CL61,
1.0 M accepted profiles, `incloud/<key>.npz`): for the same profile, C_L,v = C_L,oper x I_v / I_oper with
I_x = sum over the window gates of (beta_before / T2_wv) x eta_x(range) (cbh_slope_windowed.py; matches exact
re-runs to a median 0.998). "oper" = the 2026-07 tables (PVC a_G 5.5 um, alpha 10 /km, `_ETA_PVC55_*`).

Outputs (C:/DATA/Projects/202611_eta_study/alpha_ladder/):
  level_by_cbh.csv      median C_L,v / C_L,oper per type in CBH bins around 0.5, 1 and 2 km
  cbh_slope.csv         per-type CBH slope of ln C_L (unit x month fixed effects, 95 % CI clustered by unit-day)
  cl61_ratios.csv       cloud / Rayleigh ratio of medians per CL61 unit and variant (Payerne: hood dark v3 and v1)
  daily_lnC.csv.gz      daily medians per unit and variant
  summary.txt

Run:  python validation/eta_study/alpha_ladder.py [--workers 8]
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

STUDY = Path("C:/DATA/Projects/202611_eta_study")
DIAG = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
SRC = STUDY / "incloud"
OUT = STUDY / "alpha_ladder"
VARIANTS = {"oper": None}
for ag in (3.3, 3.7):
    for al in (10.0, 20.0, 30.0):
        VARIANTS[f"ag{ag:.1f}_a{al:.0f}"] = (ag, al)
RAY = {
    "0-20000-0-06610_C": [("Rayleigh+dark v3", DIAG / "diag_v22_dark_v3"), ("Rayleigh+dark v1", DIAG / "diag_v22_l55s008dark"),
                          ("Rayleigh no dark", DIAG / "diag_v22_l55s008")],
    "0-20000-0-10393_C": [("Rayleigh no dark", DIAG / "diag_v22_l55s008")],
}
RAY_OTHERS = STUDY / "step4" / "rayleigh_v22_l55_cl61"
CBH_BINS = {"0.5 km": (0.25, 0.75), "1 km": (0.75, 1.25), "2 km": (1.75, 2.25)}


def tables(itype):
    import calibration.cloud._filters as CF
    from validation.multiple_scattering_eta import CBH_GRID_KM, eta_table
    oper = getattr(CF, f"_ETA_PVC55_{itype}", None)
    if oper is None:                      # code without the 2026-11 tables: _ETA_<type> IS the 2026-07 table
        oper = getattr(CF, f"_ETA_{itype}")
    n = oper.shape[0]
    out = {"oper": oper}
    for name, spec in VARIANTS.items():
        if spec is None:
            continue
        ag, al = spec
        out[name] = np.column_stack([CBH_GRID_KM[:13], eta_table(itype, ag, alpha=al * 1e-3)[:13]])
    return out


def unit_stats(job):
    key, workers_unused = job
    from calibration.cloud._filters import _interp1_linear_extrap
    z = np.load(SRC / f"{key}.npz", allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    itype = meta["itype"]
    tab = tables(itype)
    B = z["beta_before"].astype("float64")
    T2 = z["t2_wv"].astype("float64")
    rr = z["range_rel_m"].astype("float64")
    r_abs = (z["peak_range_m"].astype("float64")[:, None] + rr[None, :]) / 1000.0
    w = np.where(np.isfinite(B) & np.isfinite(T2) & (T2 > 0), B / np.where(T2 > 0, T2, 1.0), 0.0)
    I = {v: (w * _interp1_linear_extrap(t[:, 0], t[:, 1], r_abs)).sum(1) for v, t in tab.items()}
    c = z["c_profile"].astype("float64")
    cbh = z["cbh_m"].astype("float64") / 1000.0
    ok = np.isfinite(c) & (c > 0) & np.isfinite(cbh) & (I["oper"] > 0)
    for v in I:
        ok &= np.isfinite(I[v]) & (I[v] > 0)
    date = z["date"].astype(int)[ok]
    month = date // 100
    x = cbh[ok]
    lnC = {v: -np.log(c[ok]) + np.log(I[v][ok] / I["oper"][ok]) for v in I}
    ratio = {v: I[v][ok] / I["oper"][ok] for v in I}
    # per-profile ratio by CBH bin (level)
    lev = []
    for bname, (lo, hi) in CBH_BINS.items():
        m = (x >= lo) & (x < hi)
        for v in I:
            if m.sum():
                lev.append(dict(key=key, itype=itype, bin=bname, variant=v, n=int(m.sum()),
                                med_ratio=float(np.median(ratio[v][m]))))
    daily, stats = [], []
    days = np.unique(date)
    for d in days:
        m = date == d
        daily.append(dict(key=key, itype=itype, date=int(d), n=int(m.sum()), cbh=float(np.median(x[m])),
                          **{f"lnC_{v}": float(np.median(lnC[v][m])) for v in I}))
    xt = x.copy()
    yt = {v: lnC[v].copy() for v in I}
    for g in np.unique(month):
        m = month == g
        xt[m] -= x[m].mean()
        for v in I:
            yt[v][m] -= lnC[v][m].mean()
    for d in days:
        m = date == d
        stats.append(dict(key=key, itype=itype, date=int(d), n=int(m.sum()), sxx=float((xt[m] ** 2).sum()),
                          **{f"sxy_{v}": float((xt[m] * yt[v][m]).sum()) for v in I}))
    return daily, stats, lev


def fe_slope(s, v):
    sxx = s.sxx.sum()
    b = s[f"sxy_{v}"].sum() / sxx
    g = (s[f"sxy_{v}"] - b * s.sxx).to_numpy()
    G = len(g)
    se = np.sqrt(G / (G - 1) * (g ** 2).sum()) / sxx
    return 100 * b, 100 * (b - 1.96 * se), 100 * (b + 1.96 * se)


def rayleigh_series(path):
    import pandas as pd
    f = Path(path)
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d = d[(d.method == "rayleigh") & (d.flag > 0) & (d.cal_value > 0)]
    return d[["date", "cal_value"]].astype({"date": int}) if len(d) else None


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    keys = sorted(p.stem for p in SRC.glob("*.npz"))
    D, S, L = [], [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for d, s, lev in ex.map(unit_stats, [(k, a.workers) for k in keys], chunksize=1):
            D += d
            S += s
            L += lev
    D, S, L = pd.DataFrame(D), pd.DataFrame(S), pd.DataFrame(L)
    D.to_csv(OUT / "daily_lnC.csv.gz", index=False, float_format="%.6g")
    lines = [f"{len(keys)} units, {int(D.n.sum())} profiles, {len(D)} unit-days; variants {list(VARIANTS)}"]

    # level per type and CBH bin: median over units of the unit median ratio (weights units equally)
    lev_rows = []
    lines.append("\n== C_L change against the 2026-07 tables (PVC 5.5 um, alpha 10), median over units of the per-unit "
                 "median per-profile ratio, % ==")
    hdr = f"{'variant':14s}" + "".join(f"{it + ' ' + b:>14s}" for it in ("CL31", "CL51", "CL61") for b in CBH_BINS)
    lines.append(hdr)
    for v in VARIANTS:
        if v == "oper":
            continue
        row = f"{v:14s}"
        for it in ("CL31", "CL51", "CL61"):
            for b in CBH_BINS:
                s = L[(L.itype == it) & (L.bin == b) & (L.variant == v)]
                val = 100 * (s.med_ratio.median() - 1) if len(s) else np.nan
                lev_rows.append(dict(itype=it, bin=b, variant=v, n_units=len(s), change_pct=val))
                row += f"{val:>13.1f} "
        lines.append(row)
    pd.DataFrame(lev_rows).to_csv(OUT / "level_by_cbh.csv", index=False, float_format="%.3f")

    # CBH slopes
    sl_rows = []
    lines.append("\n== CBH slope of ln C_L, %/km, unit x month fixed effects [95 % CI clustered by unit-day] ==")
    lines.append(f"{'variant':14s}" + "".join(f"{it:>26s}" for it in ("CL31", "CL51", "CL61")))
    for v in VARIANTS:
        row = f"{v:14s}"
        for it in ("CL31", "CL51", "CL61"):
            s = S[S.itype == it]
            b, lo, hi = fe_slope(s, v)
            sl_rows.append(dict(itype=it, variant=v, n_units=s.key.nunique(), slope=b, lo=lo, hi=hi))
            row += f"{b:>+8.2f} [{lo:+.2f},{hi:+.2f}] "
        lines.append(row)
    pd.DataFrame(sl_rows).to_csv(OUT / "cbh_slope.csv", index=False, float_format="%.3f")

    # CL61 cloud / Rayleigh ratios
    rat_rows = []
    lines.append("\n== CL61: cloud C_L / Rayleigh C_L, ratio of medians over the common period ==")
    for key in sorted(D[D.itype == "CL61"].key.unique()):
        dk = D[D.key == key]
        refs = RAY.get(key, [("Rayleigh no dark", RAY_OTHERS)])
        for rname, folder in refs:
            ray = rayleigh_series(Path(folder) / key / f"{key}_cal.csv")
            if ray is None:
                continue
            t0 = max(dk.date.min(), ray.date.min())
            t1 = min(dk.date.max(), ray.date.max())
            c2 = dk[(dk.date >= t0) & (dk.date <= t1)]
            r2 = ray[(ray.date >= t0) & (ray.date <= t1)]
            if len(c2) < 5 or len(r2) < 5:
                continue
            rmed = float(np.median(r2.cal_value))
            row = dict(key=key, reference=rname, period=f"{t0}-{t1}", n_cloud=len(c2), n_ray=len(r2))
            for v in VARIANTS:
                row[v] = float(np.exp(np.median(c2[f"lnC_{v}"])) / rmed)
            rat_rows.append(row)
    R = pd.DataFrame(rat_rows)
    R.to_csv(OUT / "cl61_ratios.csv", index=False, float_format="%.4f")
    show = R.copy()
    show["unit"] = show.key.str.replace("0-20000-0-", "").str.replace("0-20008-0-", "").str.replace("0-203-10-", "").str.replace("0-380-5-", "").str.replace("0-756-4-", "")
    cols = ["unit", "reference", "n_cloud", "n_ray"] + list(VARIANTS)
    lines.append(show[cols].round(3).to_string(index=False))
    oth = R[(R.reference == "Rayleigh no dark") & (~R.key.isin(["0-20000-0-06610_C"]))]
    lines.append("median over the dark-free CL61 units (" + str(len(oth)) + "): " +
                 "  ".join(f"{v} {oth[v].median():.3f}" for v in VARIANTS))
    lines.append("median |ln ratio| over them: " + "  ".join(f"{v} {np.median(np.abs(np.log(oth[v]))):.3f}" for v in VARIANTS))
    txt = "\n".join(lines)
    print(txt)
    (OUT / "summary.txt").write_text(txt, encoding="utf-8")


if __name__ == "__main__":
    main()
