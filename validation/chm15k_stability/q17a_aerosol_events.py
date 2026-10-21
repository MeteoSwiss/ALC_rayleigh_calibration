# -*- coding: utf-8 -*-
"""q17a - Aerosol event studies and geography: is the residual CHM15k annual cycle free-tropospheric aerosol?

Hypothesis under test: the solstice-phased ~+-9 % (after kappa; ~+-17 % raw) annual cycle of the CHM15k
Rayleigh constant C_L is aerosol contaminating the 1064 nm 'molecular' reference (C_op = C_true * R_window).
Predictions (operator): (1) lower C after strong rain (washout) and recovery over days, (2) higher C during /
after dust and smoke advection, (3) geographic structure following an aerosol climatology, (4) an amplitude
decreasing with height towards an aerosol-free value.

Inputs (read-only):
  q15/q15f_nightly.csv            per-night lnC, T_int, window, AERONET aod/ang (136 units, main module)
  q15/q15e_series.npz             per-night ln X (clear-sky molecular-normalised signal) in 12 range bands
  q16/q16b_units.csv              kappa_in per unit (intra-night), lat/lon
  q16/q16a_runs/<unit>/*.csv      forced fixed-window recalibrations (20 units)
  q17/q17a_era5_tp_daily.csv      ERA5 daily precipitation (q17a_precip_era5.py)
  q17/q17a_cams_points.npz        CAMS AOD by species + aerosol extinction 1064 nm on 9 levels (q17a_cams_ads.py)
Outputs: q17/q17a_*.csv, q17a_*.png

Quantities per night (all on ln scale, x100 = %):
  lnC      raw operational constant
  res      lnC - kappa_in * (T_int - T_int_31d): kappa applied to synoptic / sub-monthly T deviations only
  resf     lnC - kappa_in * T_int (full kappa, for the residual-amplitude geography)
  anomalies = value - centred 31-day running median (per unit)
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from q1_temperature import running_median  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path("C:/DATA/Projects/202610_CHM15k_stability")
Q15, Q16, Q17 = ROOT / "q15", ROOT / "q16", ROOT / "q17"
REPO = Path(__file__).resolve().parents[2]
CENSUS = REPO / "validation/scope_l1_2026_census.json"
PAYERNE = "0-20000-0-06610_A"
W = 2 * np.pi / 365.25
PEAK_REF = 176.0                       # network C_L peak DOY (in-phase projection)
LR = 50.0                              # sr, lidar ratio used to turn CAMS extinction into backscatter
BMOL0 = 0.089                          # Mm-1 sr-1, molecular backscatter at 1064 nm, sea level
HMOL = 8.0                             # km, molecular scale height
BANDS = [(0.2, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0),
         (3.0, 4.0), (4.0, 5.0), (5.0, 6.0), (6.0, 7.0), (7.0, 8.0)]
ML_ALT_KM = {89: 7.77, 92: 6.89, 95: 6.00, 98: 5.13, 101: 4.30, 104: 3.54, 107: 2.87, 110: 2.30, 113: 1.81}
RNG = np.random.default_rng(17)
NBOOT = 1000


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def anomaly(t, x, window=31.0, clip=0.6):
    """31-day running-median anomaly; |anomaly| > clip (ln units, ~5 sigma) set to NaN (gross outliers)."""
    a = x - running_median(t, x, window)
    if clip is not None:
        a = np.where(np.abs(a) > clip, np.nan, a)
    return a


def harmonic(t, y):
    """OLS annual harmonic; returns (amp [same units as y], peak DOY, in-phase amp on PEAK_REF, n)."""
    ok = np.isfinite(t) & np.isfinite(y)
    if ok.sum() < 40 or np.ptp(t[ok]) < 300:
        return np.nan, np.nan, np.nan, int(ok.sum())
    X = np.column_stack([np.ones(ok.sum()), np.cos(W * t[ok]), np.sin(W * t[ok])])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    amp = np.hypot(c[1], c[2])
    ph = np.arctan2(c[2], c[1])                       # cos(W t - ph)
    peak = (ph / W) % 365.25
    # t is days since 1970-01-01, DOY 1 of 1970 is t=0 -> peak DOY ~ peak+1 (drift < 1 d over the period)
    peak_doy = (peak + 1.0) % 365.25
    inph = amp * np.cos(W * (peak_doy - PEAK_REF))
    return amp, peak_doy, inph, int(ok.sum())


def doy(t):
    return (pd.to_datetime(np.asarray(t), unit="D").dayofyear.values).astype(float)


def bmol(z_km):
    return BMOL0 * np.exp(-np.asarray(z_km) / HMOL)


def cluster_mean(df, ycol, gcol, groups, ucol="stream", nboot=NBOOT):
    """Pooled mean of ycol per group with a unit-cluster bootstrap. Returns DataFrame(group, n, n_units, mean, lo, hi, se)."""
    d = df[[ucol, gcol, ycol]].dropna()
    units = d[ucol].unique()
    ui = {u: i for i, u in enumerate(units)}
    S = np.zeros((len(units), len(groups)))
    N = np.zeros((len(units), len(groups)))
    gi = {g: j for j, g in enumerate(groups)}
    for (u, g), s in d.groupby([ucol, gcol])[ycol]:
        if g in gi:
            S[ui[u], gi[g]] += s.sum()
            N[ui[u], gi[g]] += len(s)
    mean = S.sum(0) / np.maximum(N.sum(0), 1)
    w = RNG.multinomial(len(units), np.ones(len(units)) / len(units), size=nboot)
    bs = (w @ S) / np.maximum(w @ N, 1)
    rows = []
    for j, g in enumerate(groups):
        rows.append(dict(group=g, n=int(N[:, j].sum()), n_units=int((N[:, j] > 0).sum()), mean=mean[j],
                         lo=np.percentile(bs[:, j], 2.5), hi=np.percentile(bs[:, j], 97.5), se=bs[:, j].std()))
    return pd.DataFrame(rows), bs, mean


def cluster_diff(df, ycol, mask_a, mask_b, ucol="stream", nboot=NBOOT):
    """mean(y | a) - mean(y | b), pooled, unit-cluster bootstrap; also the within-unit (paired) version."""
    d = df.assign(_g=np.where(mask_a, "a", np.where(mask_b, "b", "x")))
    tab, bs, mean = cluster_mean(d, ycol, "_g", ["a", "b"], ucol, nboot)
    diff = mean[0] - mean[1]
    bd = bs[:, 0] - bs[:, 1]
    if tab.n[0] == 0 or tab.n[1] == 0:
        return dict(diff=np.nan, lo=np.nan, hi=np.nan, se=np.nan, n_a=int(tab.n[0]), n_b=int(tab.n[1]),
                    units_a=int(tab.n_units[0]), units_b=int(tab.n_units[1]), z=np.nan)
    return dict(diff=diff, lo=np.percentile(bd, 2.5), hi=np.percentile(bd, 97.5), se=bd.std(),
                n_a=int(tab.n[0]), n_b=int(tab.n[1]), units_a=int(tab.n_units[0]), units_b=int(tab.n_units[1]),
                z=diff / bd.std() if bd.std() > 0 else np.nan)


# ---------------------------------------------------------------------------------------------
# data assembly
# ---------------------------------------------------------------------------------------------
def load_nights():
    # full-precision per-night table (the CSV copy is rounded to 6 digits, too coarse to match q15e by lnC)
    frames = []
    for meta, rows in np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"]:
        f = pd.DataFrame({k: np.asarray(v) for k, v in rows.items()})
        f.insert(0, "stream", meta["stream"])
        f.insert(1, "module", str(meta["module"]))
        frames.append(f)
    n = pd.concat(frames, ignore_index=True)
    u = pd.read_csv(Q16 / "q16b_units.csv")
    cen = {f"{s['wmo']}_{s['ident']}": float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    kap = dict(zip(u.stream, u.kappa_in))
    n = n.sort_values(["stream", "t"]).reset_index(drop=True)
    n["alt_km"] = n.stream.map(cen) / 1000.0
    n["lat"] = n.stream.map(dict(zip(u.stream, u.lat)))
    n["lon"] = n.stream.map(dict(zip(u.stream, u.lon)))
    n["site"] = n.stream.map(dict(zip(u.stream, u.site)))
    n["dref"] = np.floor(n.t - 0.5).astype(int)          # the evening day of the night (UTC)
    ser = np.load(Q15 / "q15e_series.npz")
    skeys = set(k.split("__")[0] for k in ser.files)
    out = []
    for s, g in n.groupby("stream", sort=False):
        g = g.copy()
        t = g.t.values
        k = kap.get(s, np.nan)
        k = 1.5 if not np.isfinite(k) else k
        g["kappa"] = k
        T31 = running_median(t, g.T_int.values, 31.0)
        g["res"] = g.lnC.values - k / 100.0 * (g.T_int.values - T31)
        g["resf"] = g.lnC.values - k / 100.0 * g.T_int.values
        g["aC"] = anomaly(t, g.lnC.values)
        g["aR"] = anomaly(t, g.res.values)
        g["aT"] = g.T_int.values - T31
        with np.errstate(divide="ignore", invalid="ignore"):
            g["lnaod"] = np.log(g.aod.values)
        g["a_lnaod"] = anomaly(t, g.lnaod.values, clip=None) if np.isfinite(g.lnaod).sum() > 20 else np.nan
        g["a_ang"] = anomaly(t, g.ang.values, clip=None) if np.isfinite(g.ang).sum() > 20 else np.nan
        g["a_aer"] = anomaly(t, g.aer.values, clip=None)
        lx = np.full((len(g), len(BANDS)), np.nan)
        if s in skeys:
            idx = {round(float(v), 9): i for i, v in enumerate(ser[f"{s}__lnc"])}
            j = np.array([idx.get(round(float(v), 9), -1) for v in g.lnC.values])
            lx = np.where(j[:, None] >= 0, ser[f"{s}__lnx"][np.clip(j, 0, None)], np.nan)
        for jb, (a, b) in enumerate(BANDS):
            g[f"x{a:.1f}_{b:.1f}"] = lx[:, jb]
            g[f"ax{a:.1f}_{b:.1f}"] = anomaly(t, lx[:, jb], clip=0.6 if a >= 4.0 else 1.5) if np.isfinite(lx[:, jb]).sum() > 20 else np.nan
        out.append(g)
    return pd.concat(out, ignore_index=True)


def add_precip(n):
    f = Q17 / "q17a_era5_tp_daily.csv"
    pts = pd.read_csv(Q17 / "q17a_era5_points.csv")
    tp = pd.read_csv(f)
    tp["day"] = (pd.to_datetime(tp.date) - pd.Timestamp("1970-01-01")).dt.days
    gp = dict(zip(pts.stream, zip(pts.glat.round(2), pts.glon.round(2))))
    series = {}
    for (la, lo), g in tp.groupby([tp.glat.round(2), tp.glon.round(2)]):
        s = pd.Series(g.tp_mm.values, index=g.day.values).sort_index()
        c = pd.Series(g.cp_mm.values, index=g.day.values).sort_index()
        series[(la, lo)] = (s, c)
    cols = {k: np.full(len(n), np.nan) for k in ("P0", "P_prev3", "L5", "L10", "L20", "dry", "cpfrac_last10")}
    for s, ix in n.groupby("stream").groups.items():
        key = gp.get(s)
        if key not in series:
            continue
        P, C = series[key]
        d0, d1 = P.index.min(), P.index.max()
        full = P.reindex(np.arange(d0, d1 + 1))
        fullc = C.reindex(np.arange(d0, d1 + 1))
        arr, arrc = full.values, fullc.values
        for i in ix:
            d = int(n.at[i, "dref"])
            if d > d1 or d - 20 < d0:
                continue
            k = d - d0
            win = arr[k - 20:k + 1][::-1]               # lag 0..20
            cols["P0"][i] = win[0]
            cols["P_prev3"][i] = win[:3].sum()
            for thr, name in ((5, "L5"), (10, "L10"), (20, "L20")):
                hit = np.where(win > thr)[0]
                cols[name][i] = hit[0] if hit.size else 21
            wet = np.where(win >= 1.0)[0]
            cols["dry"][i] = wet[0] if wet.size else 21
            h10 = np.where(win > 10)[0]
            if h10.size:
                kk = h10[0]
                cols["cpfrac_last10"][i] = arrc[k - kk] / max(arr[k - kk], 1e-6)
    for k, v in cols.items():
        n[k] = v
    return n


def add_cams(n):
    f = Q17 / "q17a_cams_points.npz"
    if not f.exists():
        print("no CAMS extract yet")
        return n
    z = np.load(f, allow_pickle=True)
    keys = list(z["keys"])
    slt = np.round(z["sl_t"]).astype(int) if "sl_t" in z.files else None
    names = list(z["sl_names"]) if "sl_names" in z.files else []
    mlt = np.round(z["ml_t"]).astype(int) if "ml_t" in z.files else None
    out = {k: np.full(len(n), np.nan) for k in ("c_dust", "c_om", "c_bc", "c_su", "c_ss", "c_tot", "c_tot1240",
                                                  "c_pred_win", "c_pred_ft", "c_ext_ft")}
    vmap = {"duaod550": "c_dust", "omaod550": "c_om", "bcaod550": "c_bc", "suaod550": "c_su", "ssaod550": "c_ss",
            "aod550": "c_tot", "aod1240": "c_tot1240"}
    if mlt is not None:
        lev = np.asarray(z["ml_levels"]).astype(int) if "ml_levels" in z.files else np.array(sorted(ML_ALT_KM))
        zlev = np.array([ML_ALT_KM[int(l)] for l in lev])
        mld = z["ml_data"][..., 0]                       # time, station, level
        if mld.ndim == 3 and mld.shape[-1] != len(lev):
            mld = np.moveaxis(mld, 1, -1)
    for s, ix in n.groupby("stream").groups.items():
        if s not in keys:
            continue
        j = keys.index(s)
        tnext = n.loc[ix, "dref"].values + 1           # 00 UTC of the night
        if slt is not None:
            pos = {t: i for i, t in enumerate(slt)}
            ii = np.array([pos.get(t, -1) for t in tnext])
            for vi, nm in enumerate(names):
                col = vmap.get(nm)
                if col:
                    v = z["sl_data"][np.clip(ii, 0, None), j, vi]
                    out[col][ix] = np.where(ii >= 0, v, np.nan)
        if mlt is not None:
            pos = {t: i for i, t in enumerate(mlt)}
            ii = np.array([pos.get(t, -1) for t in tnext])
            ext = np.where(ii[:, None] >= 0, mld[np.clip(ii, 0, None), j, :], np.nan) * 1e6   # Mm-1
            rm1 = ext / LR / bmol(zlev)[None, :]                                              # predicted R-1
            zb = n.loc[ix, "zb_km"].values + n.loc[ix, "alt_km"].values
            zt = n.loc[ix, "zt_km"].values + n.loc[ix, "alt_km"].values
            inwin = (zlev[None, :] >= zb[:, None] - 0.3) & (zlev[None, :] <= zt[:, None] + 0.3)
            pw = np.where(inwin, rm1, np.nan)
            with np.errstate(invalid="ignore"):
                out["c_pred_win"][ix] = np.nanmean(pw, axis=1)
                ft = (zlev >= 3.4) & (zlev <= 7.0)
                out["c_pred_ft"][ix] = np.nanmean(rm1[:, ft], axis=1)
                out["c_ext_ft"][ix] = np.nanmean(ext[:, ft], axis=1)
    for k, v in out.items():
        n[k] = v
    for s, g in n.groupby("stream"):
        t = g.t.values
        for k in ("c_dust", "c_om", "c_tot", "c_pred_win", "c_pred_ft"):
            if np.isfinite(g[k]).sum() > 20:
                n.loc[g.index, "a_" + k] = anomaly(t, g[k].values, clip=None)
    return n


# ---------------------------------------------------------------------------------------------
# 1. washout
# ---------------------------------------------------------------------------------------------
YVARS = [("aC", "ln C raw"), ("aR", "ln C, synoptic-kappa residual"), ("ax2.0_2.5", "ln X 2-2.5 km"),
         ("ax2.5_3.0", "ln X 2.5-3 km"), ("ax4.0_5.0", "ln X 4-5 km"), ("ax7.0_8.0", "ln X 7-8 km"),
         ("a_lnaod", "ln AOD1064 (AERONET)"), ("a_aer", "lidar aerosol contrast"), ("a_c_tot", "CAMS AOD550"),
         ("a_c_pred_win", "CAMS-predicted R-1 in window")]
LAGBINS = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 5), (6, 8), (9, 12), (13, 21)]


def lag_label(lo, hi):
    return f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 21 else f">={lo}")


def washout(n):
    rows, curves = [], {}
    d = n[n.stream != PAYERNE].copy()
    for thr in ("L10", "L20", "L5"):
        lab = np.full(len(d), None, dtype=object)
        for lo, hi in LAGBINS:
            m = (d[thr] >= lo) & (d[thr] <= hi)
            lab[m.values] = lag_label(lo, hi)
        d["_lag"] = lab
        groups = [lag_label(lo, hi) for lo, hi in LAGBINS]
        for y, name in YVARS:
            if y not in d or not np.isfinite(d[y]).any():
                continue
            tab, _, _ = cluster_mean(d.assign(_y=100 * d[y]), "_y", "_lag", groups)
            tab["var"], tab["thr"] = y, thr
            curves[(thr, y)] = tab
            rows.append(tab)
    curves_df = pd.concat(rows, ignore_index=True)
    # contrasts: first clear night(s) after heavy rain vs long dry spell
    con = []
    for thr, mlab in (("L10", ">10 mm"), ("L20", ">20 mm")):
        a = (d[thr] <= 1).values
        b = ((d.dry >= 5) & (d.L10 >= 10)).values
        for y, name in YVARS:
            if y not in d or not np.isfinite(d[y]).any():
                continue
            r = cluster_diff(d.assign(_y=100 * d[y]), "_y", a, b)
            con.append(dict(test=f"after {mlab} (lag<=1) - dry spell (>=5 dry d, no >10 mm in 10 d)", var=y,
                            name=name, **r))
        # large-scale (frontal) vs convective washout
        for kind, mk in (("frontal (cp/tp<0.5)", d.cpfrac_last10 < 0.5), ("convective (cp/tp>=0.5)", d.cpfrac_last10 >= 0.5)):
            if thr != "L10":
                continue
            a2 = ((d.L10 <= 1) & mk).values
            for y in ("aC", "aR", "ax2.5_3.0", "ax7.0_8.0", "a_lnaod"):
                if y in d and np.isfinite(d[y]).any():
                    r = cluster_diff(d.assign(_y=100 * d[y]), "_y", a2, b)
                    con.append(dict(test=f"after >10 mm {kind} (lag<=1) - dry spell", var=y,
                                    name=dict(YVARS).get(y, y), **r))
    # dose-response: anomaly vs dry-spell length and vs log precipitation of the last 3 days (pooled, unit fixed effects)
    for y in ("aC", "aR", "ax2.5_3.0", "ax7.0_8.0", "a_lnaod", "a_c_tot"):
        if y not in d:
            continue
        sub = d[["stream", y, "P_prev3", "dry"]].dropna()
        if len(sub) < 200:
            continue
        x = np.log1p(sub.P_prev3.values)
        yy = 100 * sub[y].values
        xs = x - sub.groupby("stream").P_prev3.transform(lambda v: np.log1p(v).mean()).values
        ys = yy - sub.groupby("stream")[y].transform("mean").values * 100
        b = (xs @ ys) / (xs @ xs)
        units = sub.stream.unique()
        bs = []
        gidx = sub.groupby("stream").indices
        for _ in range(500):
            pick = RNG.choice(units, len(units))
            ii = np.concatenate([gidx[u] for u in pick])
            bs.append((xs[ii] @ ys[ii]) / (xs[ii] @ xs[ii]))
        con.append(dict(test="slope per unit ln(1+P_3d [mm])", var=y, name=dict(YVARS).get(y, y), diff=b,
                        lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), se=np.std(bs), n_a=len(sub),
                        z=b / np.std(bs)))
    con = pd.DataFrame(con)
    curves_df.to_csv(Q17 / "q17a_washout_curves.csv", index=False)
    con.to_csv(Q17 / "q17a_washout_contrasts.csv", index=False)
    return curves_df, con


# ---------------------------------------------------------------------------------------------
# 2. advection (dust / smoke)
# ---------------------------------------------------------------------------------------------
def unit_pct(d, col, q):
    return d.groupby("stream")[col].transform(lambda v: np.nanpercentile(v, q) if np.isfinite(v).sum() > 10 else np.nan)


def advection(n):
    d = n[n.stream != PAYERNE].copy()
    res = []
    # AERONET classes (column AOD at 1064, Angstrom 440-870)
    if np.isfinite(d.aod).sum() > 200:
        p85, p50 = unit_pct(d, "aod", 85), unit_pct(d, "aod", 50)
        clean = (d.aod < p50).values
        cls = {"AERONET dust (AOD>P85, AE<0.7)": ((d.aod > p85) & (d.ang < 0.7)).values,
               "AERONET fine/smoke (AOD>P85, AE>1.4)": ((d.aod > p85) & (d.ang > 1.4)).values,
               "AERONET any high (AOD>P85)": (d.aod > p85).values}
        for k, m in cls.items():
            for y, name in YVARS:
                if y in d and np.isfinite(d[y]).any():
                    r = cluster_diff(d.assign(_y=100 * d[y]), "_y", m, clean)
                    res.append(dict(source="AERONET", event=k, var=y, name=name, **r))
    if "c_dust" in d and np.isfinite(d.c_dust).sum() > 200:
        p50t = unit_pct(d, "c_tot", 50)
        clean = (d.c_tot < p50t).values
        cls = {"CAMS dust (dust AOD550>P90 and >0.05)": ((d.c_dust > unit_pct(d, "c_dust", 90)) & (d.c_dust > 0.05)).values,
               "CAMS smoke (OM AOD550>P90 and >0.05)": ((d.c_om > unit_pct(d, "c_om", 90)) & (d.c_om > 0.05)).values,
               "CAMS FT 3.4-7 km extinction > P90": (d.c_ext_ft > unit_pct(d, "c_ext_ft", 90)).values,
               "CAMS predicted window R-1 > P90": (d.c_pred_win > unit_pct(d, "c_pred_win", 90)).values}
        for k, m in cls.items():
            for y, name in YVARS:
                if y in d and np.isfinite(d[y]).any():
                    r = cluster_diff(d.assign(_y=100 * d[y]), "_y", m, clean)
                    res.append(dict(source="CAMS", event=k, var=y, name=name, **r))
        # continuous transfer: ln C anomaly on CAMS predicted R-1 anomaly (expected slope 1 if CAMS is right)
        for y in ("aC", "aR", "ax2.5_3.0", "ax4.0_5.0", "ax7.0_8.0"):
            for xcol in ("a_c_pred_win", "a_c_pred_ft"):
                if xcol not in d or y not in d:
                    continue
                sub = d[["stream", y, xcol]].dropna()
                if len(sub) < 200:
                    continue
                x, yy = sub[xcol].values, 100 * sub[y].values
                b = (x @ yy) / (x @ x) / 100.0
                gidx = sub.groupby("stream").indices
                units = list(gidx)
                bs = []
                for _ in range(500):
                    ii = np.concatenate([gidx[u] for u in RNG.choice(units, len(units))])
                    bs.append((x[ii] @ yy[ii]) / (x[ii] @ x[ii]) / 100.0)
                res.append(dict(source="CAMS", event=f"transfer d lnC / d(R-1) [{xcol}] (1 = aerosol fully transferred)",
                                var=y, name=dict(YVARS).get(y, y), diff=b, lo=np.percentile(bs, 2.5),
                                hi=np.percentile(bs, 97.5), se=np.std(bs), n_a=len(sub), z=b / np.std(bs),
                                r=np.corrcoef(x, yy)[0, 1]))
    res = pd.DataFrame(res)
    res.to_csv(Q17 / "q17a_advection.csv", index=False)
    return res


def network_daily(n):
    """Network-coherent daily anomalies vs CAMS tracers; known events."""
    d = n[n.stream != PAYERNE]
    cols = [c for c in ("aC", "aR", "ax7.0_8.0", "ax2.5_3.0", "c_dust", "c_om", "c_tot", "c_pred_ft", "c_pred_win",
                        "a_c_pred_win", "a_lnaod", "P_prev3") if c in d]
    g = d.groupby("dref")
    daily = g[cols].median()
    daily["n_units"] = g.size()
    daily = daily[daily.n_units >= 15]
    daily.index.name = "day"
    daily["date"] = pd.to_datetime(daily.index, unit="D").strftime("%Y-%m-%d")
    daily.to_csv(Q17 / "q17a_network_daily.csv")
    out = {}
    for x in ("c_dust", "c_om", "c_tot", "c_pred_win", "a_c_pred_win", "P_prev3"):
        if x in daily and np.isfinite(daily[x]).sum() > 50:
            ok = np.isfinite(daily[x]) & np.isfinite(daily.aC)
            from scipy.stats import spearmanr
            rho, p = spearmanr(daily[x][ok], daily.aC[ok])
            ok2 = np.isfinite(daily[x]) & np.isfinite(daily["ax7.0_8.0"])
            rho2, p2 = spearmanr(daily[x][ok2], daily["ax7.0_8.0"][ok2]) if ok2.sum() > 20 else (np.nan, np.nan)
            out[x] = dict(rho_aC=rho, p_aC=p, rho_x78=rho2, p_x78=p2, n_days=int(ok.sum()))
    # known events (literature / CAMS news): Canadian smoke over NW/central Europe 2025-06-01..06-10,
    # first plume S/SE Europe 2025-05-18..05-22
    ev = {"smoke 2025-06-01..06-10 (NW/central Europe)": ("2025-06-01", "2025-06-10"),
          "smoke 2025-05-18..05-22 (Mediterranean)": ("2025-05-18", "2025-05-22")}
    evrows = []
    for k, (a, b) in ev.items():
        ta = (pd.Timestamp(a) - pd.Timestamp("1970-01-01")).days
        tb = (pd.Timestamp(b) - pd.Timestamp("1970-01-01")).days
        m = (d.dref >= ta) & (d.dref <= tb)
        ref = (d.dref >= ta - 30) & (d.dref <= tb + 30) & ~m
        for y in ("aC", "aR", "ax2.5_3.0", "ax7.0_8.0", "a_lnaod", "a_c_om", "a_c_pred_win"):
            if y in d and np.isfinite(d[y]).any():
                r = cluster_diff(d.assign(_y=100 * d[y]), "_y", m.values, ref.values)
                evrows.append(dict(event=k, var=y, **r))
    evdf = pd.DataFrame(evrows)
    evdf.to_csv(Q17 / "q17a_known_events.csv", index=False)
    return daily, out, evdf


# ---------------------------------------------------------------------------------------------
# 3. geography
# ---------------------------------------------------------------------------------------------
def haversine(lat1, lon1, lat2, lon2):
    p = np.pi / 180
    a = np.sin((lat2 - lat1) * p / 2) ** 2 + np.cos(lat1 * p) * np.cos(lat2 * p) * np.sin((lon2 - lon1) * p / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(a))


def geography(n):
    from scipy.stats import spearmanr
    rows = []
    for s, g in n.groupby("stream"):
        t = g.t.values
        r = dict(stream=s, site=g.site.iloc[0], lat=g.lat.iloc[0], lon=g.lon.iloc[0], alt_km=g.alt_km.iloc[0],
                 n=len(g), kappa=g.kappa.iloc[0])
        for y, lab in (("lnC", "C"), ("res", "res"), ("resf", "resf"), ("x7.0_8.0", "x78"), ("x2.5_3.0", "x253"),
                       ("x4.0_5.0", "x45")):
            amp, pk, ip, nn = harmonic(t, 100 * g[y].values)
            r[f"amp_{lab}"], r[f"peak_{lab}"], r[f"ip_{lab}"] = amp, pk, ip
        dy = doy(t)
        jja, djf = (dy >= 152) & (dy < 244), (dy >= 335) | (dy < 60)
        for y, lab in (("lnC", "C"), ("resf", "resf")):
            r[f"sw_{lab}"] = 100 * (np.nanmedian(g[y][jja]) - np.nanmedian(g[y][djf])) if jja.sum() > 5 and djf.sum() > 5 else np.nan
        # aerosol predictors
        for y, lab in (("c_pred_win", "pred_win"), ("c_pred_ft", "pred_ft"), ("aod", "aod")):
            if y in g and np.isfinite(g[y]).sum() > 40:
                amp, pk, ip, nn = harmonic(t, 100 * g[y].values if y != "aod" else g[y].values)
                r[f"amp_{lab}"], r[f"peak_{lab}"], r[f"ip_{lab}"] = amp, pk, ip
                r[f"mean_{lab}"] = np.nanmean(g[y]) * (100 if y != "aod" else 1)
                r[f"sw_{lab}"] = (np.nanmedian(g[y][jja]) - np.nanmedian(g[y][djf])) * (100 if y != "aod" else 1)
        if "c_tot" in g and np.isfinite(g.c_tot).sum() > 40:
            r["ip_ln_ctot"] = harmonic(t, 100 * np.log(g.c_tot.clip(lower=1e-4).values))[2]
            r["ip_c_dust"] = harmonic(t, g.c_dust.values)[2]
            r["ip_c_om"] = harmonic(t, g.c_om.values)[2]
        for y in ("c_dust", "c_om", "c_tot"):
            if y in g:
                r[f"mean_{y}"] = np.nanmean(g[y])
                r[f"jja_{y}"] = np.nanmean(g[y][jja]) if jja.sum() else np.nan
        r["dist_sahara_km"] = haversine(r["lat"], r["lon"], 30.0, 5.0)
        rows.append(r)
    u = pd.DataFrame(rows)
    u.to_csv(Q17 / "q17a_units.csv", index=False)
    v = u[(u.stream != PAYERNE) & (u.lon > -20)]
    preds = ["lat", "lon", "alt_km", "dist_sahara_km", "mean_c_dust", "jja_c_om", "mean_c_tot", "ip_pred_win",
             "sw_pred_win", "ip_pred_ft", "mean_pred_win", "ip_aod", "sw_aod", "mean_aod", "ip_ln_ctot", "ip_c_dust",
             "ip_c_om"]
    targets = ["ip_C", "ip_resf", "ip_res", "ip_x78", "sw_C", "sw_resf", "amp_C", "amp_resf", "peak_C"]
    cor = []
    for tg in targets:
        for p in preds:
            if p not in v or tg not in v:
                continue
            ok = np.isfinite(v[p]) & np.isfinite(v[tg])
            if ok.sum() < 8:
                continue
            rho, pv = spearmanr(v[p][ok], v[tg][ok])
            X = np.column_stack([np.ones(ok.sum()), v[p][ok]])
            c, *_ = np.linalg.lstsq(X, v[tg][ok], rcond=None)
            cor.append(dict(target=tg, predictor=p, n=int(ok.sum()), rho=rho, p=pv, slope=c[1]))
    cor = pd.DataFrame(cor)
    cor.to_csv(Q17 / "q17a_geography_correlations.csv", index=False)
    return u, cor


# ---------------------------------------------------------------------------------------------
# 4. height dependence
# ---------------------------------------------------------------------------------------------
def height(n, u):
    d = u[u.stream != PAYERNE]
    rows = []
    for a, b in BANDS:
        lab = f"x{a:.1f}_{b:.1f}"
        amps, ips = [], []
        for s, g in n[n.stream != PAYERNE].groupby("stream"):
            amp, pk, ip, nn = harmonic(g.t.values, 100 * g[lab].values)
            if np.isfinite(amp):
                amps.append(amp)
                ips.append(ip)
        rows.append(dict(source="q15e clear-sky ln X band (all units)", z_lo=a, z_hi=b, z_mid=(a + b) / 2,
                         n_units=len(ips), amp_med=np.median(amps), ip_med=np.median(ips),
                         ip_q25=np.percentile(ips, 25), ip_q75=np.percentile(ips, 75)))
    # CAMS-predicted R-1 annual in-phase amplitude per level (station nearest point), at the units' altitude
    z = np.load(Q17 / "q17a_cams_points.npz", allow_pickle=True) if (Q17 / "q17a_cams_points.npz").exists() else None
    if z is not None and "ml_t" in z.files:
        lev = np.asarray(z["ml_levels"]).astype(int) if "ml_levels" in z.files else np.array(sorted(ML_ALT_KM))
        mld = z["ml_data"][..., 0]
        tt = z["ml_t"]
        keys = list(z["keys"])
        alts = dict(zip(n.stream, n.alt_km))
        for jl, L in enumerate(lev):
            zl = ML_ALT_KM[int(L)]
            ips, means, sws = [], [], []
            for s in d.stream:
                if s not in keys:
                    continue
                e = mld[:, keys.index(s), jl] * 1e6 / LR / bmol(zl) * 100      # % (R-1)
                if not np.isfinite(e).any() or zl < alts.get(s, 0) + 0.3:
                    continue
                amp, pk, ip, nn = harmonic(tt, e)
                ips.append(ip)
                means.append(np.nanmean(e))
                dy = doy(tt)
                sws.append(np.nanmean(e[(dy >= 152) & (dy < 244)]) - np.nanmean(e[(dy >= 335) | (dy < 60)]))
            if ips:
                rows.append(dict(source="CAMS predicted R-1 (ext1064/50 sr / beta_mol), all days 00 UTC", z_lo=zl, z_hi=zl,
                                 z_mid=zl, n_units=len(ips), amp_med=np.median(means), ip_med=np.median(ips),
                                 ip_q25=np.percentile(ips, 25), ip_q75=np.percentile(ips, 75), sw_med=np.median(sws)))
    # q16a fixed windows (20 units)
    runs = load_q16a()
    if runs is not None:
        for run, g0 in runs.groupby("run"):
            ips = []
            for s, g in g0.groupby("stream"):
                amp, pk, ip, nn = harmonic(g.t.values, 100 * g.lnc.values)
                if np.isfinite(ip):
                    ips.append(ip)
            zz = {"op": (np.nan, np.nan), "z2_3": (2, 3), "z3_4": (3, 4), "z4_5": (4, 5), "z5_6": (5, 6), "z6_8": (6, 8)}[run]
            rows.append(dict(source=f"q16a forced window {run} (20 units)", z_lo=zz[0], z_hi=zz[1],
                             z_mid=np.nanmean(zz) if np.isfinite(zz[0]) else np.nan, n_units=len(ips),
                             ip_med=np.median(ips) if ips else np.nan,
                             ip_q25=np.percentile(ips, 25) if ips else np.nan,
                             ip_q75=np.percentile(ips, 75) if ips else np.nan))
        # same 20 units: clear-sky ln X bands (to compare X with the forced-window C like for like)
        su = set(runs.stream)
        for a, b in BANDS[6:]:
            lab = f"x{a:.1f}_{b:.1f}"
            ips = [harmonic(g.t.values, 100 * g[lab].values)[2] for s, g in n[n.stream.isin(su)].groupby("stream")]
            ips = [x for x in ips if np.isfinite(x)]
            rows.append(dict(source="q15e ln X band, the 20 q16a units", z_lo=a, z_hi=b, z_mid=(a + b) / 2,
                             n_units=len(ips), ip_med=np.median(ips), ip_q25=np.percentile(ips, 25),
                             ip_q75=np.percentile(ips, 75)))
    h = pd.DataFrame(rows)
    h.to_csv(Q17 / "q17a_height_amplitude.csv", index=False)
    return h, runs


def load_q16a():
    root = Q16 / "q16a_runs"
    frames = []
    u = pd.read_csv(Q16 / "q16b_units.csv")
    kap = dict(zip(u.stream, u.kappa_in))
    for dd in sorted(root.iterdir()):
        if not dd.is_dir():
            continue
        for f in sorted(dd.glob("*.csv")):
            x = pd.read_csv(f, usecols=["wmo", "ident", "date", "run", "flag", "c", "t_mean", "tint_mean", "zb_agl", "zt_agl"])
            frames.append(x[x.run.isin(["op", "z2_3", "z3_4", "z4_5", "z5_6", "z6_8"]) & (x.flag > 0) & (x.c > 0)])
    if not frames:
        return None
    r = pd.concat(frames, ignore_index=True)
    r["stream"] = r.wmo + "_" + r.ident.astype(str)
    r["t"] = r.t_mean
    r["lnc"] = np.log(r.c)
    r["dref"] = (pd.to_datetime(r.date.astype(str), format="%Y%m%d") - pd.Timestamp("1970-01-01")).dt.days
    out = []
    for (s, run), g in r.groupby(["stream", "run"]):
        g = g.sort_values("t").copy()
        # the Kalman global IQR screen per series (as operations) - simple 4-MAD clip here
        med = np.median(g.lnc)
        mad = 1.4826 * np.median(np.abs(g.lnc - med))
        g = g[np.abs(g.lnc - med) < 4 * mad]
        if len(g) < 30:
            continue
        k = kap.get(s, 1.5)
        k = 1.5 if not np.isfinite(k) else k
        T31 = running_median(g.t.values, g.tint_mean.values, 31.0)
        g["aC"] = anomaly(g.t.values, g.lnc.values)
        g["aR"] = anomaly(g.t.values, g.lnc.values - k / 100 * (g.tint_mean.values - T31))
        out.append(g)
    return pd.concat(out, ignore_index=True)


def events_by_height(n, runs):
    """Event anomalies (washout / dust / smoke / CAMS FT) on every height: q15e ln X bands and q16a windows."""
    d = n[n.stream != PAYERNE].copy()
    tests = {}
    tests["washout: lag<=1 after >10 mm - dry spell"] = ((d.L10 <= 1).values, ((d.dry >= 5) & (d.L10 >= 10)).values)
    if "c_dust" in d and np.isfinite(d.c_dust).any():
        cl = (d.c_tot < unit_pct(d, "c_tot", 50)).values
        tests["CAMS dust > P90 - clean"] = (((d.c_dust > unit_pct(d, "c_dust", 90)) & (d.c_dust > 0.05)).values, cl)
        tests["CAMS smoke (OM) > P90 - clean"] = (((d.c_om > unit_pct(d, "c_om", 90)) & (d.c_om > 0.05)).values, cl)
        tests["CAMS FT ext > P90 - clean"] = ((d.c_ext_ft > unit_pct(d, "c_ext_ft", 90)).values, cl)
    if np.isfinite(d.aod).sum() > 200:
        tests["AERONET AOD > P85 - < P50"] = ((d.aod > unit_pct(d, "aod", 85)).values, (d.aod < unit_pct(d, "aod", 50)).values)
    rows = []
    for k, (a, b) in tests.items():
        for (lo, hi) in BANDS:
            y = f"ax{lo:.1f}_{hi:.1f}"
            r = cluster_diff(d.assign(_y=100 * d[y]), "_y", a, b, nboot=400)
            rows.append(dict(test=k, level=f"ln X {lo:.1f}-{hi:.1f} km", z_mid=(lo + hi) / 2, **r))
        for y in ("aC", "aR"):
            r = cluster_diff(d.assign(_y=100 * d[y]), "_y", a, b, nboot=400)
            rows.append(dict(test=k, level=f"ln C ({y})", z_mid=np.nan, **r))
    # q16a forced windows: join event flags by (stream, dref)
    if runs is not None:
        flags = d[["stream", "dref", "L10", "dry", "aod"] + [c for c in ("c_dust", "c_om", "c_tot", "c_ext_ft") if c in d]]
        rr = runs.merge(flags, on=["stream", "dref"], how="left")
        for run, g in rr.groupby("run"):
            tw = {"washout: lag<=1 after >10 mm - dry spell": ((g.L10 <= 1).values, ((g.dry >= 5) & (g.L10 >= 10)).values)}
            if "c_dust" in g and np.isfinite(g.c_dust).any():
                cl = (g.c_tot < unit_pct(g, "c_tot", 50)).values
                tw["CAMS dust > P90 - clean"] = (((g.c_dust > unit_pct(g, "c_dust", 90)) & (g.c_dust > 0.05)).values, cl)
                tw["CAMS smoke (OM) > P90 - clean"] = (((g.c_om > unit_pct(g, "c_om", 90)) & (g.c_om > 0.05)).values, cl)
                tw["CAMS FT ext > P90 - clean"] = ((g.c_ext_ft > unit_pct(g, "c_ext_ft", 90)).values, cl)
            for k, (a, b) in tw.items():
                r = cluster_diff(g.assign(_y=100 * g.aC), "_y", a, b, nboot=400)
                zz = {"op": np.nan, "z2_3": 2.5, "z3_4": 3.5, "z4_5": 4.5, "z5_6": 5.5, "z6_8": 7.0}[run]
                rows.append(dict(test=k, level=f"q16a C {run}", z_mid=zz, **r))
    e = pd.DataFrame(rows)
    e.to_csv(Q17 / "q17a_events_by_height.csv", index=False)
    return e


# ---------------------------------------------------------------------------------------------
# 5. attribution: event transfer coefficient x seasonal aerosol proxy -> aerosol share of the cycle
# ---------------------------------------------------------------------------------------------
def joint_diff(d, a, b, cols, nboot=NBOOT):
    """Event-minus-reference differences of several columns on a COMMON sample, same bootstrap draws."""
    m = np.all([np.isfinite(d[c].values) for c in cols], axis=0) & (a | b)
    dd = d[m]
    ga = a[m]
    units = dd.stream.unique()
    ui = dd.stream.map({u: i for i, u in enumerate(units)}).values
    S = np.zeros((len(units), len(cols), 2))
    N = np.zeros((len(units), 2))
    for j, c in enumerate(cols):
        v = 100 * dd[c].values
        np.add.at(S[:, j, 0], ui[ga], v[ga])
        np.add.at(S[:, j, 1], ui[~ga], v[~ga])
    np.add.at(N[:, 0], ui[ga], 1)
    np.add.at(N[:, 1], ui[~ga], 1)
    w = RNG.multinomial(len(units), np.ones(len(units)) / len(units), size=nboot)
    full = S.sum(0)[:, 0] / max(N[:, 0].sum(), 1) - S.sum(0)[:, 1] / max(N[:, 1].sum(), 1)
    bs = np.einsum("bu,ujg->bjg", w, S) / np.maximum(np.einsum("bu,ug->bg", w, N)[:, None, :], 1)
    bsd = bs[:, :, 0] - bs[:, :, 1]
    return dict(zip(cols, full)), {c: bsd[:, j] for j, c in enumerate(cols)}, int(N[:, 0].sum()), int(N[:, 1].sum())


def attribution(n, u):
    d = n[n.stream != PAYERNE].copy()
    d["P253"] = d["ax2.5_3.0"] - d["ax7.0_8.0"]
    d["P45"] = d["ax4.0_5.0"] - d["ax7.0_8.0"]
    ev = {"washout (lag<=1 after >10 mm) vs dry spell": ((d.L10 <= 1).values, ((d.dry >= 5) & (d.L10 >= 10)).values),
          "aerosol build-up: 9+ d after >10 mm vs 1-3 d after": (((d.L10 >= 9) & (d.L10 <= 21)).values,
                                                                 ((d.L10 >= 1) & (d.L10 <= 3)).values)}
    if np.isfinite(d.aod).sum() > 200:
        p85, p50 = unit_pct(d, "aod", 85), unit_pct(d, "aod", 50)
        ev["AERONET AOD>P85 vs <P50"] = ((d.aod > p85).values, (d.aod < p50).values)
        ev["AERONET dust (AE<0.7) vs <P50"] = (((d.aod > p85) & (d.ang < 0.7)).values, (d.aod < p50).values)
        ev["AERONET fine/smoke (AE>1.4) vs <P50"] = (((d.aod > p85) & (d.ang > 1.4)).values, (d.aod < p50).values)
    if "c_dust" in d and np.isfinite(d.c_dust).sum() > 200:
        cl = (d.c_tot < unit_pct(d, "c_tot", 50)).values
        ev["CAMS dust>P90 vs clean"] = (((d.c_dust > unit_pct(d, "c_dust", 90)) & (d.c_dust > 0.05)).values, cl)
        ev["CAMS smoke (OM)>P90 vs clean"] = (((d.c_om > unit_pct(d, "c_om", 90)) & (d.c_om > 0.05)).values, cl)
        if "c_ext_ft" in d and np.isfinite(d.c_ext_ft).sum() > 200:
            ev["CAMS FT ext (3.4-7 km)>P90 vs clean"] = ((d.c_ext_ft > unit_pct(d, "c_ext_ft", 90)).values, cl)
    # seasonal in-phase amplitudes of the lidar aerosol proxies and of C, per unit (raw series)
    seas = {}
    for s, g in d.groupby("stream"):
        t = g.t.values
        r = {}
        for lab, y in (("P253", g["x2.5_3.0"] - g["x7.0_8.0"]), ("P45", g["x4.0_5.0"] - g["x7.0_8.0"]),
                       ("res", g.res), ("resf", g.resf), ("lnC", g.lnC), ("x78", g["x7.0_8.0"])):
            r[lab] = harmonic(t, 100 * y.values)[2]
        seas[s] = r
    S = pd.DataFrame(seas).T
    med = S.median()
    rows = []
    for k, (a, b) in ev.items():
        cols = ["aR", "aC", "P253", "P45", "ax7.0_8.0"]
        full, bs, na, nb = joint_diff(d, a, b, cols)
        for prox in ("P253", "P45"):
            kk = full["aR"] / full[prox] if full[prox] != 0 else np.nan
            kb = bs["aR"] / bs[prox]
            pred = kk * med[prox]
            pb = kb * med[prox]
            rows.append(dict(event=k, proxy=prox, n_a=na, n_b=nb, d_aR=full["aR"],
                             d_aR_lo=np.percentile(bs["aR"], 2.5), d_aR_hi=np.percentile(bs["aR"], 97.5),
                             d_aC=full["aC"], d_proxy=full[prox], d_proxy_lo=np.percentile(bs[prox], 2.5),
                             d_proxy_hi=np.percentile(bs[prox], 97.5), d_x78=full["ax7.0_8.0"],
                             k=kk, k_lo=np.nanpercentile(kb, 2.5), k_hi=np.nanpercentile(kb, 97.5),
                             seasonal_ip_proxy=med[prox], pred_ip_C_aer=pred,
                             pred_lo=np.nanpercentile(pb, 2.5), pred_hi=np.nanpercentile(pb, 97.5),
                             obs_ip_res=med["res"], obs_ip_resf=med["resf"], obs_ip_lnC=med["lnC"],
                             frac_of_res=pred / med["res"], frac_of_resf=pred / med["resf"]))
    # AERONET continuous: synoptic slope of aR on ln AOD anomaly (unit FE) x seasonal ip of ln AOD
    sub = d[["stream", "aR", "a_lnaod", "t", "lnaod"]].dropna()
    if len(sub) > 300:
        x, y = sub.a_lnaod.values, 100 * sub.aR.values
        gidx = sub.groupby("stream").indices
        kk = (x @ y) / (x @ x)
        kb = []
        for _ in range(NBOOT):
            ii = np.concatenate([gidx[uu] for uu in RNG.choice(list(gidx), len(gidx))])
            kb.append((x[ii] @ y[ii]) / (x[ii] @ x[ii]))
        ips = []
        for s, g in d.groupby("stream"):
            if np.isfinite(g.lnaod).sum() > 60:
                ips.append(harmonic(g.t.values, g.lnaod.values)[2])
        ipa = np.nanmedian(ips)
        rows.append(dict(event="AERONET synoptic regression aR ~ ln AOD anomaly (24 units)", proxy="ln AOD1064",
                         n_a=len(sub), k=kk, k_lo=np.percentile(kb, 2.5), k_hi=np.percentile(kb, 97.5),
                         seasonal_ip_proxy=100 * ipa, pred_ip_C_aer=kk * ipa,
                         pred_lo=np.percentile(kb, 2.5) * ipa, pred_hi=np.percentile(kb, 97.5) * ipa,
                         obs_ip_res=med["res"], obs_ip_resf=med["resf"], obs_ip_lnC=med["lnC"],
                         frac_of_res=kk * ipa / med["res"], frac_of_resf=kk * ipa / med["resf"]))
    for lab, preds, logx in (("CAMS ln total AOD550", ["c_tot"], True),
                             ("CAMS species AOD550 (dust, OM, BC, sulphate, sea salt)", ["c_dust", "c_om", "c_bc", "c_su", "c_ss"], False),
                             ("CAMS predicted window R-1 (ext1064/50 sr)", ["c_pred_win"], False),
                             ("CAMS predicted FT 3.4-7 km R-1", ["c_pred_ft"], False)):
        r = regression_attribution(d, preds, logx, lab, med)
        if r:
            rows.append(r)
    # clean-half test: the annual cycle fitted only on the cleaner half of the nights (proxy below its own
    # 61-day running median) vs the dirtier half. An aerosol-driven cycle must shrink on the clean half.
    for lab, proxy in (("CAMS total AOD550", "c_tot"), ("CAMS dust+OM AOD550", "c_dustom"),
                       ("lidar ln X(2.5-3) - ln X(7-8)", "P253raw"), ("AERONET AOD1064", "aod")):
        dd = d.copy()
        dd["c_dustom"] = dd.get("c_dust", np.nan) + dd.get("c_om", np.nan)
        dd["P253raw"] = dd["x2.5_3.0"] - dd["x7.0_8.0"]
        ipc, ipd, ipa = [], [], []
        for s, g in dd.groupby("stream"):
            if np.isfinite(g[proxy]).sum() < 80:
                continue
            rm = running_median(g.t.values, g[proxy].values, 61.0)
            cl = (g[proxy].values < rm)
            dr = (g[proxy].values >= rm)
            a1 = harmonic(g.t.values[cl], 100 * g.res.values[cl])[2]
            a2 = harmonic(g.t.values[dr], 100 * g.res.values[dr])[2]
            a3 = harmonic(g.t.values[np.isfinite(g[proxy].values)], 100 * g.res.values[np.isfinite(g[proxy].values)])[2]
            if np.isfinite(a1) and np.isfinite(a2):
                ipc.append(a1)
                ipd.append(a2)
                ipa.append(a3)
        if ipc:
            dlt = np.array(ipd) - np.array(ipc)
            bsm = [np.median(RNG.choice(dlt, len(dlt))) for _ in range(NBOOT)]
            rows.append(dict(event=f"clean-half test: annual cycle on nights with {lab} below vs above its 61-d median",
                             proxy=proxy, n_a=len(ipc), obs_ip_res=np.median(ipa), pred_ip_C_aer=np.median(dlt),
                             pred_lo=np.percentile(bsm, 2.5), pred_hi=np.percentile(bsm, 97.5),
                             k=f"clean {np.median(ipc):.2f} / dirty {np.median(ipd):.2f}"))
    A = pd.DataFrame(rows)
    A.to_csv(Q17 / "q17a_attribution.csv", index=False)
    S.to_csv(Q17 / "q17a_seasonal_proxies_units.csv")
    return A, S


def regression_attribution(d, preds, logx, lab, med):
    """Synoptic (31-day anomaly) multiple regression of aR on CAMS predictors with unit fixed effects, then the
    predicted in-phase annual C amplitude = sum_i k_i x median-over-units in-phase amplitude of predictor i."""
    if not all(p in d and np.isfinite(d[p]).sum() > 500 for p in preds):
        return None
    X, Xr = [], []
    for p in preds:
        raw = np.log(d[p].clip(lower=1e-4)) if logx else d[p]
        an = np.full(len(d), np.nan)
        for s, ix in d.groupby("stream").indices.items():
            g = d.iloc[ix]
            an[ix] = anomaly(g.t.values, raw.values[ix], clip=None)
        X.append(an)
        Xr.append(raw.values)
    X = np.column_stack(X)
    Xr = np.column_stack(Xr)
    y = 100 * d.aR.values
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    st = d.stream.values[ok]
    Xo, yo = X[ok], y[ok]
    k, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    gidx = pd.Series(np.arange(ok.sum())).groupby(st).indices
    kb = []
    for _ in range(500):
        ii = np.concatenate([gidx[u] for u in RNG.choice(list(gidx), len(gidx))])
        kb.append(np.linalg.lstsq(Xo[ii], yo[ii], rcond=None)[0])
    kb = np.array(kb)
    # per-unit seasonal in-phase amplitude of each predictor (raw, on the calibration nights)
    ips = []
    for s, ix in d.groupby("stream").indices.items():
        tt = d.t.values[ix]
        ips.append([harmonic(tt, Xr[ix, j])[2] for j in range(len(preds))])
    ip_med = np.nanmedian(np.array(ips), axis=0)
    pred = float(k @ ip_med)
    pb = kb @ ip_med
    r2 = 1 - np.sum((yo - Xo @ k) ** 2) / np.sum((yo - yo.mean()) ** 2)
    return dict(event=f"synoptic regression aR ~ {lab}", proxy="+".join(preds), n_a=int(ok.sum()),
                k=";".join(f"{v:.3g}" for v in k), k_lo=";".join(f"{v:.3g}" for v in np.percentile(kb, 2.5, axis=0)),
                k_hi=";".join(f"{v:.3g}" for v in np.percentile(kb, 97.5, axis=0)),
                seasonal_ip_proxy=";".join(f"{v:.4g}" for v in ip_med), pred_ip_C_aer=pred,
                pred_lo=np.percentile(pb, 2.5), pred_hi=np.percentile(pb, 97.5), r2=r2,
                obs_ip_res=med["res"], obs_ip_resf=med["resf"], obs_ip_lnC=med["lnC"],
                frac_of_res=pred / med["res"], frac_of_resf=pred / med["resf"])


def interannual(n):
    """Same calendar month 2026 - 2025 per unit (identical astronomy): does the C change follow aerosol?"""
    d = n[n.stream != PAYERNE].copy()
    dt = pd.to_datetime(d.dref, unit="D")
    d["year"], d["month"] = dt.dt.year, dt.dt.month
    d["ln_ctot"] = np.log(d.c_tot) if "c_tot" in d else np.nan
    d["ln_cdust"] = np.log(d.c_dust) if "c_dust" in d else np.nan
    cols = [c for c in ("res", "lnC", "x7.0_8.0", "x2.5_3.0", "lnaod", "ln_ctot", "ln_cdust", "c_pred_win", "P_prev3")
            if c in d]
    g = d.groupby(["stream", "month", "year"])[cols].agg(["mean", "count"])
    rows = []
    for (s, m), gg in g.groupby(level=[0, 1]):
        if (s, m, 2025) not in gg.index or (s, m, 2026) not in gg.index:
            continue
        r = dict(stream=s, month=m)
        for c in cols:
            a, b = gg.loc[(s, m, 2026), (c, "mean")], gg.loc[(s, m, 2025), (c, "mean")]
            na, nb = gg.loc[(s, m, 2026), (c, "count")], gg.loc[(s, m, 2025), (c, "count")]
            r[f"d_{c}"] = (a - b) if (na >= 3 and nb >= 3) else np.nan
        rows.append(r)
    D = pd.DataFrame(rows)
    D["d_P253"] = D["d_x2.5_3.0"] - D["d_x7.0_8.0"]          # aerosol at 2.5-3 km relative to the 7-8 km reference
    D.to_csv(Q17 / "q17a_interannual_units.csv", index=False)
    mon = D.groupby("month").median(numeric_only=True)
    mon["n_units"] = D.groupby("month").size()
    mon.to_csv(Q17 / "q17a_interannual_months.csv")
    out = []
    for x in ("d_lnaod", "d_ln_ctot", "d_ln_cdust", "d_c_pred_win", "d_P253", "d_x7.0_8.0", "d_P_prev3"):
        if x not in D:
            continue
        sub = D[["stream", "month", "d_res", x]].dropna()
        if len(sub) < 20:
            continue
        xx = sub[x].values - sub.groupby("month")[x].transform("mean").values      # within-month (unit contrast)
        yy = 100 * (sub.d_res.values - sub.groupby("month").d_res.transform("mean").values)
        b = (xx @ yy) / (xx @ xx)
        # month-level (network) association
        mm = sub.groupby("month")[["d_res", x]].median()
        rho = np.corrcoef(mm[x], mm.d_res)[0, 1] if len(mm) > 3 else np.nan
        bs = []
        gidx = sub.groupby("stream").indices
        for _ in range(500):
            ii = np.concatenate([gidx[u] for u in RNG.choice(list(gidx), len(gidx))])
            bs.append((xx[ii] @ yy[ii]) / (xx[ii] @ xx[ii]))
        out.append(dict(predictor=x, n=len(sub), slope_within_month=b, lo=np.percentile(bs, 2.5),
                        hi=np.percentile(bs, 97.5), r_month_medians=rho, n_months=len(mm)))
    O = pd.DataFrame(out)
    O.to_csv(Q17 / "q17a_interannual_tests.csv", index=False)
    return mon, O


def profile_fit(h):
    """Fit A(z) = A_inf + A1 exp(-(z - 2.75)/H) to the in-phase amplitude of ln X for z >= 2.5 km."""
    from scipy.optimize import curve_fit
    a = h[(h.source == "q15e clear-sky ln X band (all units)") & (h.z_lo >= 2.5)]
    z, y = a.z_mid.values, a.ip_med.values
    f = lambda zz, ainf, a1, H: ainf + a1 * np.exp(-(zz - 2.75) / H)  # noqa: E731
    p, cov = curve_fit(f, z, y, p0=(12, 25, 1.5), maxfev=20000)
    return dict(A_inf=p[0], A1=p[1], H_km=p[2], se=np.sqrt(np.diag(cov)).tolist(),
                aerosol_at_7p5=p[1] * np.exp(-(7.5 - 2.75) / p[2]))


# ---------------------------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------------------------
def fig_washout(curves, con):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
    groups = [lag_label(lo, hi) for lo, hi in LAGBINS]
    xx = np.arange(len(groups))
    for ax, thr in zip(axs[:2], ("L10", "L20")):
        for y, col in (("aC", "k"), ("aR", "tab:blue"), ("ax2.5_3.0", "tab:orange"), ("ax7.0_8.0", "tab:green"),
                       ("a_aer", "tab:red"), ("a_c_tot", "tab:purple")):
            c = curves[(curves.thr == thr) & (curves["var"] == y)]
            if c.empty:
                continue
            c = c.set_index("group").reindex(groups)
            ax.errorbar(xx + 0.05 * list(dict(aC=0, aR=1, **{"ax2.5_3.0": 2, "ax7.0_8.0": 3, "a_aer": 4, "a_c_tot": 5}).keys()).index(y),
                        c["mean"], yerr=[c["mean"] - c.lo, c.hi - c["mean"]], color=col, marker="o", capsize=2,
                        label=dict(YVARS)[y])
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xticks(xx)
        ax.set_xticklabels(groups)
        ax.set_xlabel(f"days since last day with ERA5 precipitation > {thr[1:]} mm")
        ax.set_ylabel("31-day anomaly [%] (pooled, 95 % unit-bootstrap CI)")
        ax.set_title(f"Superposed epoch after > {thr[1:]} mm/day")
        ax.set_ylim(-25, 15)
    axs[0].legend(fontsize=8, loc="lower right")
    c = con[con.test.str.startswith("after >10 mm (lag")]
    ax = axs[2]
    yy = np.arange(len(c))
    ax.errorbar(c["diff"], yy, xerr=[c["diff"] - c.lo, c.hi - c["diff"]], fmt="o", color="k", capsize=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(c.name, fontsize=8)
    ax.axvline(0, color="grey")
    ax.set_xlabel("first night(s) after > 10 mm minus dry spell [%]")
    ax.set_title("Washout contrast (positive controls vs C)")
    fig.tight_layout()
    fig.savefig(Q17 / "q17a_washout.png", dpi=130)
    plt.close(fig)


def fig_advection(adv, daily, A=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(18, 9))
    ax = fig.add_axes([0.05, 0.56, 0.9, 0.38])
    dd = daily.copy()
    t = pd.to_datetime(dd.index, unit="D")
    ax.plot(t, 100 * dd.aC.rolling(5, center=True, min_periods=3).mean(), "k", label="network median ln C anomaly [%] (5-d mean)")
    if "ax7.0_8.0" in dd:
        ax.plot(t, 100 * dd["ax7.0_8.0"].rolling(5, center=True, min_periods=3).mean(), color="tab:green", lw=1,
                label="network median ln X 7-8 km anomaly [%]")
    ax.set_ylabel("anomaly [%]")
    ax.set_ylim(-20, 20)
    ax2 = ax.twinx()
    for c, col, lab in (("c_dust", "tab:brown", "CAMS dust AOD550"), ("c_om", "tab:red", "CAMS OM AOD550"),
                        ("c_pred_win", "tab:purple", "CAMS predicted R-1 window /10")):
        if c in dd:
            ax2.plot(t, dd[c] / (10 if c == "c_pred_win" else 1), color=col, lw=0.9, alpha=0.8, label=lab)
    ax2.set_ylabel("network median CAMS (calibration nights)")
    for a, b in (("2025-06-01", "2025-06-10"), ("2025-05-18", "2025-05-22")):
        ax.axvspan(pd.Timestamp(a), pd.Timestamp(b), color="orange", alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)
    ax2.legend(loc="upper right", fontsize=8)
    ax.set_title("Network-coherent daily anomalies vs CAMS aerosol tracers (orange = Canadian smoke over Europe, CAMS)")
    ax3 = fig.add_axes([0.21, 0.07, 0.30, 0.38])
    ev = adv[(adv["var"].isin(["aC", "aR", "ax2.5_3.0", "ax7.0_8.0", "a_lnaod"])) & ~adv.event.str.startswith("transfer")]
    ev = ev[np.isfinite(ev["diff"])]
    evs = list(dict.fromkeys(ev.event))
    vs = ["aC", "aR", "ax2.5_3.0", "ax7.0_8.0", "a_lnaod"]
    cols = ["k", "tab:blue", "tab:orange", "tab:green", "tab:red"]
    for iv, (v, col) in enumerate(zip(vs, cols)):
        e = ev[ev["var"] == v].set_index("event").reindex(evs)
        yy = np.arange(len(evs)) + 0.15 * (iv - 2)
        ax3.errorbar(e["diff"], yy, xerr=[e["diff"] - e.lo, e.hi - e["diff"]], fmt="o", color=col, capsize=2,
                     label=dict(YVARS)[v], ms=4)
    ax3.set_yticks(np.arange(len(evs)))
    ax3.set_yticklabels(evs, fontsize=7)
    ax3.axvline(0, color="grey")
    ax3.set_xlabel("event minus clean nights, 31-day anomaly [%] (95 % CI)")
    ax3.legend(fontsize=7, loc="lower right")
    ax3.set_xlim(-15, 130)
    ax3.set_title("Advection events: C rises by +4...+6 %, X(7-8 km) barely moves")
    ax4 = fig.add_axes([0.70, 0.07, 0.27, 0.38])
    if A is not None and not A.empty:
        a = A[~A.event.str.startswith("clean") & ~((A.proxy == "P45") & (A.event.str.contains("smoke (OM)", regex=False)))]
        a = a[np.isfinite(a.pred_ip_C_aer) & (a.pred_ip_C_aer.abs() < 50)]
        lab = [f"{e[:48]} [{p}]" for e, p in zip(a.event, a.proxy)]
        yy = np.arange(len(a))
        ax4.errorbar(a.pred_ip_C_aer, yy, xerr=[a.pred_ip_C_aer - a.pred_lo, a.pred_hi - a.pred_ip_C_aer], fmt="o",
                     color="k", capsize=3)
        ax4.set_yticks(yy)
        ax4.set_yticklabels(lab, fontsize=6)
        med = A.obs_ip_res.dropna().iloc[0]
        medf = A.obs_ip_resf.dropna().iloc[0]
        ax4.axvline(med, color="tab:blue", ls="--", label=f"observed in-phase amplitude, synoptic-kappa residual ({med:.1f} %)")
        ax4.axvline(medf, color="tab:purple", ls=":", label=f"observed, full-kappa residual ({medf:.1f} %)")
        ax4.axvline(0, color="grey")
        ax4.set_xlim(-5, 25)
        ax4.legend(fontsize=7, loc="lower right")
        ax4.set_xlabel("aerosol-predicted in-phase annual amplitude of ln C [%] (95 % CI)")
        ax4.set_title("Attribution: event sensitivity x seasonal aerosol cycle")
    fig.savefig(Q17 / "q17a_advection.png", dpi=130)
    plt.close(fig)


def fig_geography(u, cor):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    v = u[(u.stream != PAYERNE) & (u.lon > -20)].copy()
    v["pred_cams"] = 0.0534 * v.get("ip_ln_ctot", np.nan)
    cor = pd.concat([cor, cor[cor.predictor == "ip_ln_ctot"].assign(predictor="pred_cams")], ignore_index=True)
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    sc = axs[0].scatter(v.lon, v.lat, c=v.ip_resf, cmap="viridis", s=35, vmin=0, vmax=25, edgecolor="k", lw=0.3)
    plt.colorbar(sc, ax=axs[0], label="in-phase annual amplitude, full-kappa residual [%]")
    axs[0].set_xlabel("longitude")
    axs[0].set_ylabel("latitude")
    axs[0].set_title("Residual amplitude per unit")
    for ax, p, lab in ((axs[1], "lon", "longitude [deg E]"), (axs[2], "dist_sahara_km", "distance to 30N 5E [km]"),
                       (axs[3], "pred_cams", "aerosol-predicted in-phase amplitude of ln C [%] (synoptic k = 5.3 %/ln AOD x CAMS ln AOD cycle)")):
        if p not in v:
            continue
        ax.scatter(v[p], v.ip_resf, s=20, label="residual (full kappa)", color="tab:blue")
        ax.scatter(v[p], v.ip_C, s=12, label="raw ln C", color="k", marker="x")
        for tg in ("ip_resf", "ip_C"):
            c = cor[(cor.target == tg) & (cor.predictor == p)]
            if not c.empty:
                ax.text(0.02, 0.95 - 0.07 * (tg == "ip_C"), f"{tg}: rho={c.rho.iloc[0]:+.2f} (p={c.p.iloc[0]:.2f}, n={c.n.iloc[0]})",
                        transform=ax.transAxes, fontsize=8)
        ax.set_xlabel(lab)
        ax.set_ylabel("in-phase annual amplitude [%]")
        if p == "pred_cams":
            ax.plot([-2, 30], [-2, 30], "r--", lw=0.8, label="1:1 (aerosol explains all)")
            ax.set_xlim(-2, 8)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(Q17 / "q17a_geography.png", dpi=130)
    plt.close(fig)


def fig_height(h, e):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    ax = axs[0]
    a = h[h.source == "q15e clear-sky ln X band (all units)"]
    a = a[a.z_lo >= 1.0]
    try:
        pf = profile_fit(h)
        zz = np.linspace(2.5, 8.0, 50)
        ax.plot(pf["A_inf"] + pf["A1"] * np.exp(-(zz - 2.75) / pf["H_km"]), zz, "k:", lw=1.2,
                label=f"fit A_inf + A1 exp(-z/H): A_inf={pf['A_inf']:.1f} %, H={pf['H_km']:.2f} km")
        ax.axvline(pf["A_inf"], color="k", ls="--", lw=0.6)
    except Exception:  # noqa: BLE001
        pass
    s20 = h[h.source == "q15e ln X band, the 20 q16a units"]
    if not s20.empty:
        ax.plot(s20.ip_med, s20.z_mid, "o--", color="tab:blue", mfc="none", label="clear-sky ln X, same 20 units")
    ax.errorbar(a.ip_med, a.z_mid, xerr=[a.ip_med - a.ip_q25, a.ip_q75 - a.ip_med], fmt="o-", color="k", capsize=3,
                label="observed: clear-sky ln X per band (all units, median, IQR)")
    q = h[h.source.str.startswith("q16a") & np.isfinite(h.z_mid)]
    ax.errorbar(q.ip_med, q.z_mid, xerr=[q.ip_med - q.ip_q25, q.ip_q75 - q.ip_med], fmt="s-", color="tab:blue",
                capsize=3, label="observed: forced-window C (q16a, 20 units)")
    c = h[h.source.str.startswith("CAMS")]
    if not c.empty:
        ax.errorbar(c.ip_med, c.z_mid, xerr=[c.ip_med - c.ip_q25, c.ip_q75 - c.ip_med], fmt="^-", color="tab:red",
                    capsize=3, label="CAMS-predicted aerosol R-1 at that height (in-phase)")
    ax.axvline(0, color="grey")
    ax.set_xlabel("in-phase annual amplitude [%] (projected on DOY 176)")
    ax.set_ylabel("height [km] (bands AGL; CAMS levels ASL)")
    ax.set_title("Annual amplitude vs height: observed vs aerosol-predicted")
    ax.legend(fontsize=8)
    ax = axs[1]
    cols = {"washout: lag<=1 after >10 mm - dry spell": "tab:blue", "CAMS dust > P90 - clean": "tab:brown",
            "CAMS smoke (OM) > P90 - clean": "tab:red", "CAMS FT ext > P90 - clean": "tab:purple",
            "AERONET AOD > P85 - < P50": "tab:orange"}
    for k, col in cols.items():
        s = e[(e.test == k) & e.level.str.startswith("ln X")]
        s = s[s.z_mid >= 1.0]
        if s.empty:
            continue
        ax.errorbar(s["diff"], s.z_mid, xerr=[s["diff"] - s.lo, s.hi - s["diff"]], fmt="o-", color=col, capsize=2, label=k)
        cc = e[(e.test == k) & (e.level == "ln C (aC)")]
        if not cc.empty:
            ax.errorbar(cc["diff"], [8.6], xerr=[cc["diff"] - cc.lo, cc.hi - cc["diff"]], fmt="*", color=col, ms=12, capsize=2)
    ax.axvline(0, color="grey")
    ax.set_xlabel("event minus reference, 31-day anomaly of ln X [%] (stars at 8.6 km = ln C)")
    ax.set_ylabel("height AGL [km]")
    ax.set_title("Event anomalies by height (clear-sky molecular-normalised signal)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(Q17 / "q17a_height.png", dpi=130)
    plt.close(fig)


def main():
    n = load_nights()
    n = add_precip(n)
    n = add_cams(n)
    n.to_csv(Q17 / "q17a_nights.csv", index=False)
    print("nights", len(n), "units", n.stream.nunique(), "with precip", np.isfinite(n.L10).sum(),
          "with CAMS", np.isfinite(n.get("c_tot", pd.Series(np.nan))).sum(), flush=True)
    curves, con = washout(n)
    print(con.to_string(max_colwidth=60), flush=True)
    adv = advection(n)
    print(adv.to_string(max_colwidth=70), flush=True)
    daily, netcor, evdf = network_daily(n)
    print(json.dumps(netcor, indent=1, default=float))
    print(evdf.to_string(), flush=True)
    u, cor = geography(n)
    print(cor.sort_values("p").head(40).to_string(), flush=True)
    h, runs = height(n, u)
    print(h.to_string(), flush=True)
    e = events_by_height(n, runs)
    print(e.to_string(max_colwidth=50), flush=True)
    pf = profile_fit(h)
    print("profile fit", pf, flush=True)
    json.dump(pf, open(Q17 / "q17a_profile_fit.json", "w"), indent=1, default=float)
    A, S = attribution(n, u)
    print(A.round(3).to_string(), flush=True)
    print("seasonal proxies (median over units):", S.median().round(2).to_dict(), flush=True)
    mon, O = interannual(n)
    print((mon * 100).round(1).to_string(), flush=True)
    print(O.round(3).to_string(), flush=True)
    fig_washout(curves, con)
    fig_advection(adv, daily, A)
    fig_geography(u, cor)
    fig_height(h, e)


if __name__ == "__main__":
    sys.exit(main())
