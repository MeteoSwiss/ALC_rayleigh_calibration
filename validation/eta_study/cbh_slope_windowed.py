"""CBH flatness of the per-profile cloud C_L for eta variants, network-wide, from the stored
in-cloud windows (no re-calibration).

eta is a per-range multiplier, so for the SAME accepted profile the lidar constant of a variant v is
  C_L,v = C_L,oper * I_v / I_a,   I_x = sum over the window gates of (beta_before / T2_wv) * eta_x(range)
(the window peak-60..peak+300 m holds nearly all of the integrated cloud signal; the rest -- below-cloud
aerosol, <= 5 % by the ratio filter -- carries eta ratios within a few % of the window's, so the error is
<~0.2 %; checked against the exact CL61 re-runs of step 4). ln C_L,oper = ln(cc_applied) - ln(c_profile)
and the unit constant drops out of every slope below.

Variants per type: a = operational table (_filters._ETA_<type>); b = legacy Hopkin ladder (CL31: own;
CL51 and CL61: the CL51 ladder, as before 2026-07); c = PVC (type optics, a_G 3.3 um, alpha = the type's
network alpha_k, --alpha-k CL31=..,CL51=..,CL61=..).

Statistics (slope of ln C_L against CBH, in %/km):
  (a) daily: Theil-Sen on day medians, per-unit normalised (ln C_L - unit median, CBH - unit median);
  (b) per profile, unit fixed effects (within-unit demeaned OLS), CI clustered by unit-day;
  (c) per profile, unit x month fixed effects (removes season), CI clustered by unit-day.
Per-unit sufficient statistics are written by `--stats` (small CSVs, runnable on CSCS next to the data)
and pooled by `--pool` (also for subsets: per type, with / without Payerne CL31).

Run:  python validation/eta_study/cbh_slope_windowed.py --stats --alpha-k CL31=27,CL51=30,CL61=31 [--src DIR] [--out DIR]
      python validation/eta_study/cbh_slope_windowed.py --pool --stats-dirs DIR1;DIR2 [--validate]
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

A_G = 3.3
VARIANTS = ("a", "b", "c")
PAYERNE_CL31 = "0-20000-0-06610_B"


def _interp_extrap(x, y, xi):
    """MATLAB interp1 'linear','extrap' -- the operational eta interpolation."""
    from calibration.cloud._filters import _interp1_linear_extrap
    return _interp1_linear_extrap(x, y, xi)


def tables(itype, alpha_k):
    import calibration.cloud._filters as CF
    from validation.multiple_scattering_eta import CBH_GRID_KM, eta_table
    # "operational" of this study = the 2026-07 tables; since 2026-11 _ETA_<type> holds the v2 tables
    oper = {"CL31": getattr(CF, "_ETA_PVC55_CL31", CF._ETA_CL31), "CL51": getattr(CF, "_ETA_PVC55_CL51", CF._ETA_CL51),
            "CL61": getattr(CF, "_ETA_PVC55_CL61", CF._ETA_CL61)}[itype]
    legacy = CF._ETA_LEGACY_CL31 if itype == "CL31" else CF._ETA_LEGACY_CL51
    n = oper.shape[0]
    pvc = np.column_stack([CBH_GRID_KM[:n], eta_table(itype, A_G, alpha=alpha_k * 1e-3)[:n]])
    return {"a": oper, "b": legacy, "c": pvc}


def unit_stats(job):
    """One unit -> (daily medians rows, sufficient-statistics rows)."""
    key, src, alpha_k = job
    z = np.load(Path(src) / f"{key}.npz", allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    itype = meta["itype"]
    tab = tables(itype, alpha_k[itype])
    B = z["beta_before"].astype("float64")
    T2 = z["t2_wv"].astype("float64")
    rr = z["range_rel_m"].astype("float64")
    r_abs = (z["peak_range_m"].astype("float64")[:, None] + rr[None, :]) / 1000.0     # km
    w = np.where(np.isfinite(B) & np.isfinite(T2) & (T2 > 0), B / np.where(T2 > 0, T2, 1.0), 0.0)
    I = {v: (w * _interp_extrap(t[:, 0], t[:, 1], r_abs)).sum(1) for v, t in tab.items()}
    c = z["c_profile"].astype("float64")
    cbh = z["cbh_m"].astype("float64") / 1000.0
    ok = np.isfinite(c) & (c > 0) & np.isfinite(cbh) & (I["a"] > 0)
    for v in VARIANTS:
        ok &= np.isfinite(I[v]) & (I[v] > 0)
    date = z["date"].astype(int)[ok]
    month = date // 100
    x = cbh[ok]
    y = {v: -np.log(c[ok]) + np.log(I[v][ok] / I["a"][ok]) for v in VARIANTS}
    daily, stats = [], []
    days = np.unique(date)
    for d in days:
        m = date == d
        daily.append(dict(key=key, itype=itype, date=int(d), n=int(m.sum()), cbh=float(np.median(x[m])),
                          **{f"lnC_{v}": float(np.median(y[v][m])) for v in VARIANTS}))
    for scheme, grp in (("unit", np.zeros_like(date)), ("unit_month", month)):
        xt = x.copy()
        yt = {v: y[v].copy() for v in VARIANTS}
        for g in np.unique(grp):
            m = grp == g
            xt[m] -= x[m].mean()
            for v in VARIANTS:
                yt[v][m] -= y[v][m].mean()
        for d in days:
            m = date == d
            stats.append(dict(key=key, itype=itype, date=int(d), scheme=scheme, n=int(m.sum()),
                              sxx=float((xt[m] ** 2).sum()),
                              **{f"sxy_{v}": float((xt[m] * yt[v][m]).sum()) for v in VARIANTS}))
    return daily, stats


def run_stats(a):
    import pandas as pd
    alpha_k = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in a.alpha_k.split(",")}
    src = Path(a.src)
    keys = sorted(p.stem for p in src.glob("*.npz"))
    if a.keys:
        want = {l.strip() for l in open(a.keys, encoding="utf-8") if l.strip()}
        keys = [k for k in keys if k in want]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    D, S = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for d, s in ex.map(unit_stats, [(k, str(src), alpha_k) for k in keys], chunksize=1):
            D += d
            S += s
    pd.DataFrame(D).to_csv(out / "daily_lnC.csv.gz", index=False, float_format="%.6g")
    pd.DataFrame(S).to_csv(out / "suffstats.csv.gz", index=False, float_format="%.8g")
    print(f"{len(keys)} units -> {out}", flush=True)


def fe_slope(s, v):
    """Pooled within-group slope and day-clustered 95 % CI from per-day sufficient statistics."""
    sxx = s.sxx.sum()
    b = s[f"sxy_{v}"].sum() / sxx
    g = (s[f"sxy_{v}"] - b * s.sxx).to_numpy()            # per-cluster score
    G = len(g)
    se = np.sqrt(G / (G - 1) * (g ** 2).sum()) / sxx
    return 100 * b, 100 * (b - 1.96 * se), 100 * (b + 1.96 * se)


def run_pool(a):
    import pandas as pd
    from scipy.stats import theilslopes
    dirs = [Path(d) for d in a.stats_dirs.split(";")]
    D = pd.concat([pd.read_csv(d / "daily_lnC.csv.gz") for d in dirs], ignore_index=True)
    S = pd.concat([pd.read_csv(d / "suffstats.csv.gz") for d in dirs], ignore_index=True)
    D = D.drop_duplicates(["key", "date"])
    S = S.drop_duplicates(["key", "date", "scheme"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    groups = [("CL31 all", D.itype == "CL31"), ("CL31 excl. Payerne CL31", (D.itype == "CL31") & (D.key != PAYERNE_CL31)),
              ("CL51 all", D.itype == "CL51"), ("CL61 all", D.itype == "CL61"),
              ("Payerne CL31 only", D.key == PAYERNE_CL31)]
    rows = []
    for gname, mask in groups:
        keys = set(D.loc[mask, "key"])
        d = D[D.key.isin(keys)].copy()
        if d.empty:
            continue
        for v in VARIANTS:
            yv = d[f"lnC_{v}"] - d.groupby("key")[f"lnC_{v}"].transform("median")
            xv = d.cbh - d.groupby("key").cbh.transform("median")
            ts, _, lo, hi = theilslopes(yv, xv)
            r = dict(group=gname, variant=v, n_units=len(keys), n_days=len(d),
                     daily_theilsen=100 * ts, daily_lo=100 * lo, daily_hi=100 * hi)
            for scheme in ("unit", "unit_month"):
                s = S[S.key.isin(keys) & (S.scheme == scheme)]
                b, lo2, hi2 = fe_slope(s, v)
                r[f"{scheme}_FE"] = b
                r[f"{scheme}_lo"] = lo2
                r[f"{scheme}_hi"] = hi2
            r["n_profiles"] = int(S[S.key.isin(keys) & (S.scheme == "unit")].n.sum())
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(out / "cbh_slope_network.csv", index=False, float_format="%.3f")
    with open(out / "cbh_slope_network.txt", "w", encoding="utf-8") as f:
        f.write(t.round(2).to_string(index=False) + "\n")
    print(t.round(2).to_string(index=False), flush=True)
    if a.validate:
        validate(D, out)


def validate(D, out):
    """Windowed variant ratios vs the exact step-4 re-runs (CL61 daily C_L, variants c/a and b/a)."""
    import pandas as pd
    ex = pd.read_csv(Path("C:/DATA/Projects/202611_eta_study/step4/daily_cloud_variants.csv"))
    ex = ex.pivot_table(index=["key", "date"], columns="variant", values="C_L")
    d = D[D.itype == "CL61"].set_index(["key", "date"])
    j = d.join(ex, how="inner")
    lines = []
    for v in ("b", "c"):
        win = np.exp(j[f"lnC_{v}"] - j["lnC_a"])
        exact = j[v] / j["a"]
        r = (win / exact).dropna()
        lines.append(f"windowed vs exact daily C_L ratio {v}/a: n={len(r)} days, median {r.median():.4f}, "
                     f"P5 {r.quantile(0.05):.4f}, P95 {r.quantile(0.95):.4f}")
    print("\n".join(lines), flush=True)
    (out / "validation_windowed_vs_exact.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", action="store_true")
    ap.add_argument("--pool", action="store_true")
    ap.add_argument("--alpha-k", default="CL31=27,CL51=30,CL61=31")
    ap.add_argument("--src", default=os.environ.get("ETA_OUT", "C:/DATA/Projects/202611_eta_study/incloud"))
    ap.add_argument("--keys", help="optional key list (one per line) to restrict --stats")
    ap.add_argument("--out", default="C:/DATA/Projects/202611_eta_study/slope")
    ap.add_argument("--stats-dirs", default="")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()
    if a.stats:
        run_stats(a)
    if a.pool:
        run_pool(a)


if __name__ == "__main__":
    main()
