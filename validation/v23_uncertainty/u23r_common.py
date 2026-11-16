# -*- coding: utf-8 -*-
"""Shared loaders and statistics of the independent review of the v2.3 uncertainty budget (u23r_*).

Truths (as in doc/v2.3/uncertainty.md section 3.2, with the leakage controls the review adds):

* time series: anomaly of ln C_L(night) - kappa_u (T_int - its 31-day median) from the running median
  of the unit's PUBLISHED series. Options: exclusion half-width around the night (the u23 truth
  excludes only the night itself), causal (past-only) window, kappa on/off, and the standard error
  of the running median itself (the truth is not error-free: 1.2533 sigma / sqrt(n));
* same-night twins: ln(C_A / C_B) anomaly from the running median of the published pair ratio.

The recalibrated constant of the night is used (not the published one), so the same truth serves
eprof_v2.2 and eprof_v2.3 runs.
"""
from __future__ import annotations

import glob
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))
sys.path.insert(0, str(HERE))

import u23_analysis as A  # noqa: E402

OUT = Path("C:/DATA/Projects/202610_v23/uncertainty_review")
U23 = Path("C:/DATA/Projects/202610_v23/uncertainty")
EPOCH = datetime(1970, 1, 1)
PAYERNE = "0-20000-0-06610_A"
SCHLEIZ = "0-20000-0-10564_0"
MESSINA = "0-20000-0-00203_A"
TWINS = A.TWINS
# Fresh co-located CHM15k / CL61 pair (not in the u23 design set)
FRESH_COLOC = [("0-20000-0-03808_A", "0-20000-0-03808_C")]


def pull_stats(z):
    return A.pull_stats(z)


def rstd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size >= 5 else np.nan


# --------------------------------------------------------------------------------------------
# Loading the review runs
# --------------------------------------------------------------------------------------------
def _parse_list(s):
    if not isinstance(s, str) or not s:
        return np.array([])
    return np.array([float(x) for x in s.split(",")])


def load_runs(set_name, method):
    """All rows of one review run (one row per unit-night)."""
    import pandas as pd
    fs = sorted(glob.glob(str(OUT / "runs" / set_name / method / "*" / "*.csv")))
    frames = []
    for f in fs:
        try:
            frames.append(pd.read_csv(f, dtype={"date": str, "ident": str, "sub": str,
                                                "sub_odd": str, "sub_even": str, "sub_h1": str,
                                                "sub_h2": str, "ens_shift": str, "ens_lr": str,
                                                "ens_time": str, "sel_msg": str}))
        except Exception:  # noqa: BLE001 - empty block
            continue
    df = pd.concat(frames, ignore_index=True)
    df["key"] = df["wmo"] + "_" + df["ident"]
    df = df.drop_duplicates(subset=["key", "date"]).reset_index(drop=True)
    parsed = [A.parse_sub(s) for s in df["sub"]]
    df["_z"] = [p[0] for p in parsed]
    df["_r"] = [p[1] for p in parsed]
    df["_s"] = [p[2] for p in parsed]
    for c in ("sub_odd", "sub_even", "sub_h1", "sub_h2", "ens_shift", "ens_lr", "ens_time"):
        if c in df:
            df["_" + c] = [_parse_list(s) for s in df[c]]
    df["ok"] = df["flag"].isin([1.0, 0.5]) & np.isfinite(df["c"]) & (df["c"] > 0)
    df["tday"] = [(datetime.strptime(d, "%Y%m%d") - EPOCH).days for d in df["date"]]
    # shipped-default budget, 1 sigma relative (the run stores the absolute 2-sigma total)
    df["sig23"] = df["def_u"] / (2.0 * df["c_all"])
    df["zc_agl"] = 0.5 * (df["zb_agl"] + df["zt_agl"])
    df["month"] = [int(d[4:6]) for d in df["date"]]
    return df


# --------------------------------------------------------------------------------------------
# Truth (i): time series
# --------------------------------------------------------------------------------------------
def running_median_ex(t, x, t_eval, window=31.0, exclude_days=0.5, causal=False, min_n=5):
    """Running median of x(t) at t_eval excluding |t - t_eval| <= exclude_days.

    Returns (median, standard error of that median, n used). The standard error is the large-n
    value 1.2533 x robust sigma / sqrt(n): the truth carries it into every anomaly."""
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    te = np.atleast_1d(np.asarray(t_eval, float))
    med = np.full(te.size, np.nan)
    se = np.full(te.size, np.nan)
    nn = np.zeros(te.size, int)
    fin = np.isfinite(x)
    for i, tt in enumerate(te):
        if causal:
            m = fin & (t >= tt - window) & (t < tt - exclude_days)
        else:
            m = fin & (np.abs(t - tt) <= window / 2.0) & (np.abs(t - tt) > exclude_days)
        n = int(m.sum())
        nn[i] = n
        if n >= min_n:
            v = x[m]
            med[i] = np.median(v)
            se[i] = 1.2533 * 1.4826 * np.median(np.abs(v - med[i])) / np.sqrt(n)
    return med, se, nn


class Truth:
    """Published series of every unit (lazy), kappa table, census."""

    def __init__(self):
        self.cen = A.census()
        self.ktab = A.kappa_table()
        self.pubs = {}

    def pub(self, key):
        if key not in self.pubs:
            wmo, ident = key.rsplit("_", 1)
            s = self.cen[(wmo, ident)]
            self.pubs[key] = A.Published(wmo, ident, s["lon"])
        return self.pubs[key]

    def kappa(self, key, typ):
        return A.unit_kappa(key, typ, self.ktab)


def add_truth(df, T: Truth, exclude_days=0.5, window=31.0, causal=False, use_kappa=True,
              col="anom", own_series=False):
    """Time-series anomaly of every successful night of ``df`` (recalibrated constant).

    own_series=True builds the running median from the recalibrated (sampled) series itself
    instead of the published one (no published series: the Mini-MPLs)."""
    df[col] = np.nan
    df[col + "_se"] = np.nan
    df[col + "_n"] = 0
    df[col + "_inpub"] = False
    for key, g in df[df.ok].groupby("key"):
        typ = g["type"].iloc[0]
        P = T.pub(key)
        kap = T.kappa(key, typ) if use_kappa else 0.0
        if not P.ok or own_series:
            if len(g) < 10:
                continue
            t = g.tday.values.astype(float)
            lnc = np.log(g.c.values)
            if kap and "tint_mean" in g and np.any(np.isfinite(g.tint_mean.values)):
                tm, _, _ = running_median_ex(t, g.tint_mean.values, t, window=window, exclude_days=-1)
                dT = np.where(np.isfinite(g.tint_mean.values) & np.isfinite(tm),
                              g.tint_mean.values - tm, 0.0)
                lnc = lnc - kap / 100.0 * dT
            base, se, n = running_median_ex(t, lnc, t, window, exclude_days, causal)
            df.loc[g.index, col] = lnc - base
            df.loc[g.index, col + "_se"] = se
            df.loc[g.index, col + "_n"] = n
            df.loc[g.index, col + "_inpub"] = True          # own-series truth: usable
            continue
        lncc = P.corrected(kap)
        dT_pub = np.zeros(P.t.size)
        if kap and np.any(np.isfinite(P.T)):
            trm = A.running_median(P.t, P.T)
            dT_pub = np.where(np.isfinite(P.T) & np.isfinite(trm), P.T - trm, 0.0)
        t_eval = g.tday.values.astype(float)
        base, se, n = running_median_ex(P.t, lncc, t_eval, window, exclude_days, causal)
        corr = np.zeros(len(g))
        inpub = np.zeros(len(g), bool)
        for k, d in enumerate(g.date.values):
            j = P.by_date.get(d)
            if j is not None:
                corr[k] = kap / 100.0 * dT_pub[j]
                inpub[k] = True
        df.loc[g.index, col] = np.log(g.c.values) - corr - base
        df.loc[g.index, col + "_se"] = se
        df.loc[g.index, col + "_n"] = n
        df.loc[g.index, col + "_inpub"] = inpub
    return df


# --------------------------------------------------------------------------------------------
# Truth (ii): same-night twins
# --------------------------------------------------------------------------------------------
def twin_frame(df, T: Truth, exclude_days=0.5, window=31.0, use_kappa=True, pairs=None):
    """One row per twin pair-night where both units calibrated: anomaly of ln(C_A/C_B) (recalibrated
    constants, kappa-corrected) from the running median of the PUBLISHED pair ratio."""
    import pandas as pd
    rows = []
    idx = df[df.ok].set_index(["key", "date"])
    if pairs is None:
        pairs = []
        for wmo, idents in TWINS.items():
            for i, a in enumerate(idents):
                for b in idents[i + 1:]:
                    pairs.append((f"{wmo}_{a}", f"{wmo}_{b}"))
    for ka, kb in pairs:
        Pa, Pb = T.pub(ka), T.pub(kb)
        if not (Pa.ok and Pb.ok):
            continue
        ta = df.loc[df.key == ka, "type"]
        tb = df.loc[df.key == kb, "type"]
        if ta.empty or tb.empty:
            continue
        kapa = T.kappa(ka, ta.iloc[0]) if use_kappa else 0.0
        kapb = T.kappa(kb, tb.iloc[0]) if use_kappa else 0.0
        la, lb = Pa.corrected(kapa), Pb.corrected(kapb)
        common = sorted(set(Pa.dates) & set(Pb.dates))
        if len(common) < 20:
            continue
        ia = np.array([Pa.by_date[d] for d in common])
        ib = np.array([Pb.by_date[d] for d in common])
        tc = Pa.t[ia]
        r_pub = la[ia] - lb[ib]
        corr_a = la - Pa.lnc                      # -kappa dT per published night
        corr_b = lb - Pb.lnc
        da = set(df.date[df.ok & (df.key == ka)])
        db = set(df.date[df.ok & (df.key == kb)])
        for d in sorted(da & db):
            ra, rb = idx.loc[(ka, d)], idx.loc[(kb, d)]
            tt = float(ra["tday"])
            base, se, n = running_median_ex(tc, r_pub, [tt], window, exclude_days)
            ja, jb = Pa.by_date.get(d), Pb.by_date.get(d)
            ca = corr_a[ja] if ja is not None else 0.0
            cb = corr_b[jb] if jb is not None else 0.0
            r = (np.log(ra["c"]) + ca) - (np.log(rb["c"]) + cb)
            rows.append(dict(site=ka.rsplit("_", 1)[0], ka=ka, kb=kb, date=d, tday=tt,
                             anom=r - base[0], anom_se=se[0], anom_n=n[0]))
    return pd.DataFrame(rows)


def twin_sigma(tw, df, sig_col, drop_cols=()):
    """sqrt(sigma_A^2 + sigma_B^2) of every twin pair-night, sigma from ``sig_col`` minus the
    columns ``drop_cols`` (the terms common to co-located units) removed in quadrature."""
    s = df.set_index(["key", "date"])
    v = s[sig_col] ** 2
    for c in drop_cols:
        v = v - s[c].fillna(0.0) ** 2
    v = np.sqrt(np.clip(v, 0.0, None))
    sa = np.array([v.get((a, d), np.nan) for a, d in zip(tw.ka, tw.date)])
    sb = np.array([v.get((b, d), np.nan) for b, d in zip(tw.kb, tw.date)])
    return np.hypot(sa, sb)


# --------------------------------------------------------------------------------------------
# Budget variants from the stored components
# --------------------------------------------------------------------------------------------
def sigma_vert_variant(df, D_km=1.0, est="abs", resid=True, g="def_g", gs="def_gs",
                       res="def_res", nsub="def_nsub"):
    """sigma_vert as the module computes it, for another D / gradient estimator (shipped range)."""
    gg = df[g].values.astype(float)
    gsv = df[gs].values.astype(float)
    rr = df[res].values.astype(float)
    ns = df[nsub].values.astype(float)
    if est == "abs":
        a = np.abs(gg)
    elif est == "debias":
        a = np.sqrt(np.clip(gg ** 2 - np.nan_to_num(gsv) ** 2, 0.0, None))
    else:
        a = np.hypot(gg, np.nan_to_num(gsv))
    out = a * D_km
    if resid:
        out = np.where((ns >= 3) & np.isfinite(rr), np.hypot(out, np.nan_to_num(rr)), out)
    fb = ~np.isfinite(gg)                     # module fallback: g ~ 2 sigma_ens / 0.3 per km
    out = np.where(fb, 2.0 * df["sig_ens"].values / 0.3 * D_km, out)
    return out


def total(df, vert, noise_col="def_sig_noise", mol_col="def_sig_mol", aer_col="def_sig_aer",
          extra=None):
    parts = (df["sig_ens"].values ** 2 + np.nan_to_num(vert) ** 2
             + np.nan_to_num(df[noise_col].values) ** 2 + np.nan_to_num(df[mol_col].values) ** 2
             + np.nan_to_num(df[aer_col].values) ** 2)
    if extra is not None:
        parts = parts + np.nan_to_num(extra) ** 2
    return np.sqrt(parts)
