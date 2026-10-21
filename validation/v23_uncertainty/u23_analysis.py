# -*- coding: utf-8 -*-
"""Offline analysis of the v2.3 uncertainty validation runs (see u23_validation.py).

Truth references (doc/v2.3/uncertainty.md):
  (i)  time series: anomaly of ln C_L from the 31-day running median of the PUBLISHED series of the
       same unit (all nights; the harness reproduces them bit-identically), after removing
       kappa x (T_int - its 31-day running median), i.e. the sub-daily instrument-temperature term;
  (ii) same-night twins: ln(C_A / C_B) anomaly from its own 31-day running median; the atmosphere,
       the molecular model and the common-mode gain cancel, so sigma_pair^2 = sigma_A^2 + sigma_B^2
       without the (common) molecular term.
A realistic uncertainty gives pulls with a robust std of ~1 and ~5 % of the nights outside 2 sigma.
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
STAB = REPO / "validation" / "chm15k_stability"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAB))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from u23_validation import (OUT, RUN_DIR, TWINS, UNITS, census, calibrate_night,  # noqa: E402
                            _options)

EPOCH = datetime(1970, 1, 1)
FIG = OUT
PAYERNE = ("0-20000-0-06610", "A")
SPECIAL = {"0-20000-0-06610_A": "Payerne", "0-20000-0-00203_A": "Messina",
           "0-20000-0-10564_0": "Schleiz"}
W_SUB = 500.0                 # recorded sub-window width (m)
KAPPA_FALLBACK = 1.5          # %/K, network intra-night median (report 19)

# --------------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------------


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


def parse_sub(s):
    if not isinstance(s, str) or "|" not in s:
        return np.array([]), np.array([]), np.array([])
    a, b, c = s.split("|")
    return (np.array([float(x) for x in a.split(",")]), np.array([float(x) for x in b.split(",")]),
            np.array([float(x) for x in c.split(",")]))


def load_runs():
    import pandas as pd
    frames = []
    for f in sorted(RUN_DIR.glob("*/*.csv")):
        try:
            frames.append(pd.read_csv(f, dtype={"date": str, "ident": str, "sub": str}))
        except Exception:  # noqa: BLE001 - empty block
            continue
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates(subset=["wmo", "ident", "date"]).reset_index(drop=True)
    parsed = [parse_sub(s) for s in df["sub"]]
    df["_z"] = [p[0] for p in parsed]
    df["_r"] = [p[1] for p in parsed]
    df["_s"] = [p[2] for p in parsed]
    df["key"] = df["wmo"] + "_" + df["ident"]
    df["ok"] = df["flag"].isin([1.0, 0.5]) & np.isfinite(df["c"]) & (df["c"] > 0)
    # v2.2 success also requires the old 'u > C' gate (not applied in the v2.3 recording)
    df["ok22"] = df["ok"] & ~(df["u_ens"] > df["c"])
    df["tday"] = [(datetime.strptime(d, "%Y%m%d") - EPOCH).days for d in df["date"]]
    return df.drop_duplicates(subset=["key", "date"]).reset_index(drop=True)


def kappa_table():
    out = {}
    p = Path("C:/DATA/Projects/202610_CHM15k_stability/q2a_intranight.csv")
    if p.exists():
        for r in csv.DictReader(open(p, encoding="utf-8")):
            k, se = _f(r.get("b_T_ctrl")), _f(r.get("se_T_ctrl"))
            if np.isfinite(k) and np.isfinite(se) and se < 0.6 and -1.0 < k < 5.0:
                out[r["stream"]] = k
    return out


def unit_kappa(key, typ, ktab):
    if key == "0-20000-0-06610_A":
        return 0.0            # report 19: Payerne's night kappa is unreliable, never apply it
    if typ != "CHM15k":
        return 0.0            # kappa is a CHM15k property (CL61 / Mini-MPL: no correction)
    return ktab.get(key, KAPPA_FALLBACK)


def running_median(t, x, window=31.0, exclude_self=False, t_eval=None):
    """Centred running median of x(t) evaluated at t_eval (default t); >= 5 points required."""
    t0 = np.asarray(t, float)
    x0 = np.asarray(x, float)
    o = np.argsort(t0, kind="stable")
    t, x = t0[o], x0[o]
    te = t0 if t_eval is None else np.asarray(t_eval, float)   # evaluated in the CALLER's order
    lo = np.searchsorted(t, te - window / 2.0)
    hi = np.searchsorted(t, te + window / 2.0, side="right")
    out = np.full(te.size, np.nan)
    for i in range(te.size):
        seg_t = t[lo[i]:hi[i]]
        seg = x[lo[i]:hi[i]]
        m = np.isfinite(seg)
        if exclude_self:
            m &= np.abs(seg_t - te[i]) > 0.5
        if m.sum() >= 5:
            out[i] = np.median(seg[m])
    return out


class Published:
    """Published nightly constants of one unit (NC_DIR), keyed by night date D, with T_int."""

    def __init__(self, wmo, ident, lon):
        from common import load_constants
        import q16a_recalibration as Q
        con = load_constants(wmo, ident)
        self.ok = con is not None and len(con["c"]) >= 20
        if not self.ok:
            return
        m = np.isfinite(con["c"]) & (con["c"] > 0)
        con = {k: v[m] for k, v in con.items()}
        mods = [str(x) for x in con["module"]]
        valid = [x for x in mods if x not in ("", "-999", "--", "nan")]
        self.main = max(set(valid), key=valid.count) if valid else ""
        keep = np.array([(x == self.main) or (self.main == "") or (x in ("", "-999", "--", "nan"))
                         for x in mods])
        self.dates = np.array(Q.night_dates_from_constants(con, lon))[keep]
        self.t = np.array([(datetime.strptime(d, "%Y%m%d") - EPOCH).days for d in self.dates], float)
        self.lnc = np.log(con["c"][keep])
        self.c = con["c"][keep]
        self.u = con["u"][keep]
        self.t0, self.t1 = con["t0"][keep], con["t1"][keep]
        hk = Q.HK(wmo, ident)
        self.T = np.full(self.t.size, np.nan)
        if hk.ok:
            for i, (a, b) in enumerate(zip(self.t0, self.t1)):
                if np.isfinite(a) and np.isfinite(b) and b > a:
                    self.T[i] = hk.mean(np.linspace(a, b, 40))[0]
        self.by_date = {d: i for i, d in enumerate(self.dates)}

    def corrected(self, kappa):
        """ln C - kappa (T - T_31d) (kappa in %/K), on every published night."""
        if kappa == 0 or not np.any(np.isfinite(self.T)):
            return self.lnc.copy()
        trm = running_median(self.t, self.T)
        dT = np.where(np.isfinite(self.T) & np.isfinite(trm), self.T - trm, 0.0)
        return self.lnc - kappa / 100.0 * dT


# --------------------------------------------------------------------------------------------
# Offline vertical term (mirrors calibration.rayleigh.uncertainty)
# --------------------------------------------------------------------------------------------
BASE_CFG = dict(below_m=1000.0, above_m=1000.0, min_alt_m=1500.0, snr_min=3.0, D_km=1.0,
                est="abs", resid=True, floor=0.01)


def select_subwindows(z, r, s, zb, zt, cfg):
    """Emulate the module's sub-window range on the recorded sub-windows."""
    if z.size == 0:
        return z, r, s
    lo = min(max(zb - cfg["below_m"], cfg["min_alt_m"]), zb)
    keep = (z >= lo - 1.0)                                  # centre at/above the range bottom
    above = (z - 0.5 * W_SUB) >= zt - 1e-6
    top = max(zt + cfg["above_m"], zt + W_SUB)              # module: sub-window END <= top
    keep &= (z + 0.5 * W_SUB) <= top + 1e-6
    # SNR stop above the window: equivalent per-gate mean SNR of a 500 m / 30 m-gate sub-window
    if np.isfinite(cfg["snr_min"]) and cfg["snr_min"] > 0:
        smax = 1.0 / (cfg["snr_min"] * np.sqrt(W_SUB / 30.0))
        bad = above & ~(np.isfinite(s) & (s <= smax))
        if np.any(bad):
            first_bad = z[bad].min()
            keep &= ~(above & (z >= first_bad))
    return z[keep], r[keep], s[keep]


def line_fit(z, r, s, zb, zt, cfg):
    """(g, g_sigma, resid, n) of the module's line fit on the emulated sub-window range."""
    from calibration.rayleigh.uncertainty import fit_log_gradient
    zz, rr, ss = select_subwindows(z, r, s, zb, zt, cfg)
    fit = fit_log_gradient(zz, rr, ss, 0.5 * (zb + zt), cfg["floor"])
    g, gs, res = fit["g"], fit["g_sigma"], fit["resid_rel"]
    if not np.isfinite(g) and rr.size == 2:
        g = (np.log(rr[1]) - np.log(rr[0])) / ((zz[1] - zz[0]) / 1000.0)
    return g, gs, res, rr.size


def vertical_from_fit(g, gs, res, D_km, est="abs", resid=True):
    """Vertical term (1 sigma, relative) from fitted gradient arrays (vectorised)."""
    g = np.asarray(g, float)
    gs = np.asarray(gs, float)
    if est == "abs":
        gg = np.abs(g)
    elif est == "debias":
        gg = np.sqrt(np.clip(g * g - np.nan_to_num(gs) ** 2, 0.0, None))
    else:                                                    # "quad"
        gg = np.hypot(g, np.nan_to_num(gs))
    sig = gg * D_km
    if resid:
        sig = np.hypot(sig, np.nan_to_num(np.asarray(res, float)))
    return sig


def add_fit(df, cfg, tag=""):
    """Line-fit columns g/gs/res/nsub for every successful night with an offline range choice."""
    out = np.full((len(df), 4), np.nan)
    for i, (z, r, s, zb, zt, ok) in enumerate(zip(df["_z"], df["_r"], df["_s"], df["zb_agl"],
                                                  df["zt_agl"], df["ok"])):
        if ok:
            out[i] = line_fit(z, r, s, float(zb), float(zt), cfg)
    for j, name in enumerate(("g", "gs", "res", "nsub")):
        df[name + tag] = out[:, j]
    return df


def add_vertical(df, cfg, col="sig_vert"):
    """Recompute sigma_vert for every successful night with an offline configuration."""
    df = add_fit(df, cfg)
    df[col] = vertical_from_fit(df["g"], df["gs"], df["res"], cfg["D_km"], cfg["est"], cfg["resid"])
    return df


def total_sigma(df, vert_col="sig_vert", with_mol=True):
    parts = [df["sig_ens"] ** 2, df[vert_col].fillna(0.0) ** 2, df["sig_noise"].fillna(0.0) ** 2]
    if with_mol:
        parts.append(df["sig_mol"].fillna(0.0) ** 2)
    return np.sqrt(sum(parts))


# --------------------------------------------------------------------------------------------
# Truth (i): anomalies from the published running median
# --------------------------------------------------------------------------------------------
def anomalies(df, cen):
    ktab = kappa_table()
    df = df.copy()
    df["anom"] = np.nan
    df["kappa"] = np.nan
    df["in_pub"] = False
    pubs = {}
    for key, g in df.groupby("key"):
        wmo, ident = key.rsplit("_", 1)
        s = cen[(wmo, ident)]
        P = Published(wmo, ident, s["lon"])
        pubs[key] = P
        if not P.ok:
            continue
        kap = unit_kappa(key, s["type"], ktab)
        lncc = P.corrected(kap)
        trm = running_median(P.t, P.T) if np.any(np.isfinite(P.T)) else np.full(P.t.size, np.nan)
        for idx, row in g.iterrows():
            if not row["ok"]:
                continue
            j = P.by_date.get(row["date"])
            if j is None:
                continue
            df.at[idx, "in_pub"] = True
            base = running_median(P.t, lncc, exclude_self=True, t_eval=[P.t[j]])[0]
            df.at[idx, "anom"] = lncc[j] - base
            df.at[idx, "kappa"] = kap
            df.at[idx, "pub_ratio"] = row["c"] / P.c[j]
            df.at[idx, "T_dev"] = (P.T[j] - trm[j]) if np.isfinite(P.T[j]) else np.nan
    return df, pubs


def pull_stats(z):
    z = np.asarray(z, float)
    z = z[np.isfinite(z)]
    if z.size < 5:
        return dict(n=int(z.size), rstd=np.nan, out2=np.nan, std=np.nan)
    return dict(n=int(z.size), rstd=float(1.4826 * np.median(np.abs(z - np.median(z)))),
                out2=float(np.mean(np.abs(z) > 2.0)), std=float(np.std(z)))


# --------------------------------------------------------------------------------------------
# Truth (ii): same-night twins
# --------------------------------------------------------------------------------------------
def twin_anomalies(df, pubs, cen):
    """Per twin pair and sampled night: ln(C_A/C_B) anomaly from the published pair's running
    median (kappa-corrected units). Returns DataFrame(site, pair, ka, kb, date, anom)."""
    import pandas as pd
    ktab = kappa_table()
    rows = []
    okset = set(zip(df.key[df.ok], df.date[df.ok]))
    for wmo, idents in TWINS.items():
        for i, a in enumerate(idents):
            for b in idents[i + 1:]:
                ka, kb = f"{wmo}_{a}", f"{wmo}_{b}"
                Pa, Pb = pubs.get(ka), pubs.get(kb)
                if Pa is None or Pb is None or not (Pa.ok and Pb.ok):
                    continue
                common = sorted(set(Pa.dates) & set(Pb.dates))
                if len(common) < 20:
                    continue
                la = Pa.corrected(unit_kappa(ka, cen[(wmo, a)]["type"], ktab))
                lb = Pb.corrected(unit_kappa(kb, cen[(wmo, b)]["type"], ktab))
                ia = np.array([Pa.by_date[d] for d in common])
                ib = np.array([Pb.by_date[d] for d in common])
                t = Pa.t[ia]
                r = la[ia] - lb[ib]
                for k, d in enumerate(common):
                    if (ka, d) in okset and (kb, d) in okset:
                        base = running_median(t, r, exclude_self=True, t_eval=[t[k]])[0]
                        rows.append(dict(site=wmo, pair=f"{a}/{b}", ka=ka, kb=kb, date=d,
                                         anom=r[k] - base))
    return pd.DataFrame(rows)


def twin_pull_stats(tw, sig_by_key_date):
    """Pull stats of the twin anomalies for a {(key, date): sigma} map (sigma without the
    common molecular term)."""
    sa = np.array([sig_by_key_date.get((a, d), np.nan) for a, d in zip(tw.ka, tw.date)])
    sb = np.array([sig_by_key_date.get((b, d), np.nan) for b, d in zip(tw.kb, tw.date)])
    z = tw.anom.values / np.hypot(sa, sb)
    return pull_stats(z), z


# --------------------------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------------------------
def _prepare():
    import pandas as pd
    cen = census()
    df = load_runs()
    df, pubs = anomalies(df, cen)
    df["site"] = [cen[(w, i)].get("site", w) for w, i in zip(df.wmo, df.ident)]
    df["lat"] = [cen[(w, i)]["lat"] for w, i in zip(df.wmo, df.ident)]
    df["sig22"] = df["u_ens"] / (2.0 * df["c"])
    return df, pubs, cen, pd


def _sweep_range(args):
    """Worker: one sub-window range choice -> per-night line-fit arrays."""
    below, above, snr, recs = args
    cfg = dict(BASE_CFG, below_m=below, above_m=above, snr_min=snr)
    out = np.full((len(recs), 4), np.nan)
    for i, (z, r, s, zb, zt) in enumerate(recs):
        out[i] = line_fit(z, r, s, zb, zt, cfg)
    return below, above, snr, out


def sweep(df, tw, workers=16):
    """Grid over the vertical-term choices: pooled CHM15k pulls (Payerne excluded) and twins."""
    import pandas as pd
    from concurrent.futures import ProcessPoolExecutor
    base = df[df.ok & (df.type == "CHM15k") & (df.key != "0-20000-0-06610_A")].copy()
    recs = [(z, r, s, float(zb), float(zt)) for z, r, s, zb, zt in
            zip(base["_z"], base["_r"], base["_s"], base["zb_agl"], base["zt_agl"])]
    ranges = [(b, a, s, recs) for b in (0.0, 500.0, 1000.0, 1500.0)
              for a in (0.0, 500.0, 1000.0, 1e9) for s in (3.0, 1.0)]
    fits = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for below, above, snr, arr in ex.map(_sweep_range, ranges):
            fits[(below, above, snr)] = arr
    has_anom = (base.in_pub & np.isfinite(base.anom)).values
    keys = list(zip(base.key, base.date))
    rows = []
    for (below, above, snr), arr in fits.items():
        g, gs, res = arr[:, 0], arr[:, 1], arr[:, 2]
        for est in ("abs", "debias", "quad"):
            for resid in (False, True):
                for D in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
                    sv = vertical_from_fit(g, gs, res, D, est, resid)
                    sv = np.where(np.isfinite(sv), sv, 0.0)
                    nomol = (base.sig_ens.values ** 2 + sv ** 2
                             + np.nan_to_num(base.sig_noise.values) ** 2)
                    sig = np.sqrt(nomol + np.nan_to_num(base.sig_mol.values) ** 2)
                    z = base.anom.values[has_anom] / sig[has_anom]
                    st = pull_stats(z)
                    unit_r = pd.Series(z).groupby(base.key.values[has_anom]).apply(
                        lambda v: pull_stats(v)["rstd"])
                    tst, _ = twin_pull_stats(tw, dict(zip(keys, np.sqrt(nomol))))
                    rows.append(dict(below=below, above=above, snr=snr, est=est, resid=resid,
                                     D=D, n=st["n"], rstd=st["rstd"], out2=st["out2"],
                                     unit_rstd_med=float(unit_r.median()),
                                     unit_rstd_p25=float(unit_r.quantile(0.25)),
                                     unit_rstd_p75=float(unit_r.quantile(0.75)),
                                     twin_n=tst["n"], twin_rstd=tst["rstd"], twin_out2=tst["out2"],
                                     sig_med=float(np.median(sig[has_anom]))))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "u23_sweep.csv", index=False)
    return out, fits, base


def cmd_sweep(args):
    import pickle
    df, pubs, cen, pd = _prepare()
    df.drop(columns=["sub", "_z", "_r", "_s"]).to_csv(OUT / "u23_nights.csv", index=False)
    ok = df[df.ok]
    print(f"{len(df)} unit-nights, {int(df.ok.sum())} successes, {df.key.nunique()} units")
    r = ok.pub_ratio.dropna()
    print(f"reproduction vs published: n={r.size}, median ratio {np.median(r):.6f}, "
          f"bit-identical {np.mean(np.abs(r - 1) < 1e-9) * 100:.1f} %, "
          f"within 1e-6 {np.mean(np.abs(r - 1) < 1e-6) * 100:.1f} %")
    tw = twin_anomalies(df, pubs, cen)
    tw.to_csv(OUT / "u23_twin_anomalies.csv", index=False)
    print(f"twin pair-nights: {len(tw)}")
    sw, fits, base = sweep(df, tw)
    with open(OUT / "u23_prepared.pkl", "wb") as fh:
        pickle.dump(dict(df=df, tw=tw), fh)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    for D in (0.75, 1.0, 1.5):
        print(f"\n--- D = {D} km, sorted by |rstd - 1| ---")
        x = sw[sw.D == D].copy()
        x["dev"] = (x.rstd - 1).abs() + (x.twin_rstd - 1).abs()
        print(x.sort_values("dev").head(25).to_string())


# --------------------------------------------------------------------------------------------
# Payerne 2026-09-07 (window 6.06-6.54 km ASL; published C_L 1.02e12 vs unit median ~5.3e11)
# --------------------------------------------------------------------------------------------
PAY_NIGHTS = ("20260907", "20260906", "20260912", "20260923")   # the case + neighbours


def _classify_payerne(days):
    """Cloudnet target classification of Payerne A (the runner's own _do_classification), so the
    night is calibrated WITH the classification screen, as the published series was."""
    import os
    os.environ["ALC_FULLCAL_DIR"] = str(OUT / "payerne_cls")
    os.environ["ALC_CAMS_DIR"] = "D:/CAMS_daily"
    os.environ["ALC_L1_ROOT"] = "D:/E-PROFILE_L1_2026"
    sys.path.insert(0, str(REPO / "scripts"))
    import run_network_calibration as RNC
    s = census()[PAYERNE]
    for d in days:
        dd = datetime.strptime(d, "%Y%m%d")
        RNC._do_classification(s, dd - timedelta(days=1), dd)
    return OUT / "payerne_cls"


def cmd_payerne(args):
    import logging
    import pickle
    from calibration import calibrate_rayleigh
    from calibration.io.classification import (classification_files_for_night,
                                               read_classification_curtain)
    from calibration.io.data_loader import build_file_paths, load_data
    import u23_validation as V
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    V._install()
    root = _classify_payerne(PAY_NIGHTS)
    s = census()[PAYERNE]
    info = V._info(s)
    key = f"{s['wmo']}_{s['ident']}"
    out = {}
    for d in PAY_NIGHTS:
        files = [f for f in build_file_paths(d, info, _options()) if f.exists()]
        native = load_data(files, info.instrument_type, _options().data_level)
        cls = read_classification_curtain(classification_files_for_night(root, key, s["wmo"], d))
        res = {}
        for label, model, params in (("ensemble", "ensemble", {}),
                                     ("record", "v2.3", None),
                                     ("default", "v2.3", {})):
            o = _options(params, model)
            for screen in (True, False):
                fi = {}
                r = calibrate_rayleigh(d, info, o, preloaded_data=native, fit_inputs_out=fi,
                                       classification=(cls if screen else None))
                res[(label, screen)] = dict(flag=r.flag, c=r.lidar_constant, u=r.uncertainty,
                                            zb=r.calibration_bottom_height,
                                            zt=r.calibration_top_height, msg=r.message,
                                            budget=fi.get("uncertainty_budget"))
                b = fi.get("uncertainty_budget")
                txt = ""
                if b is not None:
                    txt = (f" | v2.3 2s {100 * b.rel_total:.0f}% (ens {200 * b.sigma_ensemble_rel:.0f}, "
                           f"vert {200 * b.sigma_vertical_rel:.0f}, mol {200 * b.sigma_molecular_rel:.1f}, "
                           f"noise {200 * b.sigma_noise_rel:.1f}) g {100 * b.gradient_per_km:+.0f}%/km "
                           f"n {b.n_subwindows} {b.z_lo_m:.0f}-{b.z_hi_m:.0f} m")
                print(f"{d} {label:8s} screen={screen!s:5s} flag {r.flag:5} C {r.lidar_constant:.3e} "
                      f"u {r.uncertainty:.3e}{txt}", flush=True)
        out[d] = res
    with open(OUT / "u23_payerne.pkl", "wb") as fh:
        pickle.dump(out, fh)


# --------------------------------------------------------------------------------------------
# verify: the shipped defaults, rerun in the code, against the offline emulation
# --------------------------------------------------------------------------------------------
VERIFY_UNITS = [("0-20000-0-06240", "B"), ("0-20000-0-10564", "0"), ("0-20000-0-06610", "A"),
                ("0-20008-0-UGR", "A"), ("0-20000-0-07014", "A")]


def cmd_verify(args):
    """Rerun 2 blocks of 5 units with the SHIPPED defaults (code path) and compare with the
    offline recomputation from the recorded sub-windows."""
    import u23_validation as V
    import pandas as pd
    from calibration.rayleigh.uncertainty import DEFAULT_PARAMS
    tasks = []
    for wmo, ident in VERIFY_UNITS:
        files = sorted((RUN_DIR / f"{wmo}_{ident}").glob("*.csv"))
        for f in files[len(files) // 3: len(files) // 3 + 2]:
            dates = pd.read_csv(f, dtype={"date": str})["date"].tolist()
            tasks.append((wmo, ident, dates, OUT / "verify" / f"{wmo}_{ident}" / f.name,
                          dict(DEFAULT_PARAMS, max_rel_uncertainty=1e9)))
    with __import__("concurrent.futures").futures.ProcessPoolExecutor(max_workers=8) as ex:
        list(ex.map(V.work_chunk, *zip(*tasks)))
    code = pd.concat([pd.read_csv(f, dtype={"date": str, "ident": str})
                      for f in (OUT / "verify").glob("*/*.csv")])
    code["key"] = code.wmo + "_" + code.ident
    rec = load_runs()
    cfg = dict(BASE_CFG, below_m=DEFAULT_PARAMS["extend_below_m"],
               above_m=DEFAULT_PARAMS.get("extend_above_m", 1e9),
               min_alt_m=DEFAULT_PARAMS["min_altitude_m"], snr_min=DEFAULT_PARAMS["snr_min"],
               D_km=DEFAULT_PARAMS["window_displacement_m"] / 1000.0,
               resid=DEFAULT_PARAMS["include_residual"])
    rec = add_vertical(rec[rec.ok].copy(), cfg, "sig_vert_off")
    m = code.merge(rec[["key", "date", "sig_vert_off", "c"]], on=["key", "date"],
                   suffixes=("", "_rec"))
    m = m[np.isfinite(m.c) & np.isfinite(m.sig_vert_off)]
    d = m.sig_vert_rec - m.sig_vert_off
    print(f"verify: {len(m)} nights; same C {np.mean(m.c == m.c_rec) * 100:.1f} %; "
          f"code vs offline sigma_vert: median |diff| {np.median(np.abs(d)) * 100:.3f} %-pt, "
          f"p95 {np.percentile(np.abs(d), 95) * 100:.3f} %-pt, "
          f"corr {np.corrcoef(m.sig_vert_rec, m.sig_vert_off)[0, 1]:.4f}")
    m.to_csv(OUT / "u23_verify.csv", index=False)


# --------------------------------------------------------------------------------------------
# cost: runtime of the budget with the shipped defaults
# --------------------------------------------------------------------------------------------
def cmd_cost(args):
    import time as _t
    import logging
    from calibration import calibrate_rayleigh
    from calibration.io.data_loader import build_file_paths, load_data
    import u23_validation as V
    import pandas as pd
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    V._install()
    rec = load_runs()
    rec = rec[rec.ok & (rec.type == "CHM15k")]
    pick = rec.sample(n=min(args.n, len(rec)), random_state=1)
    cen = census()
    rows = []
    for _, r in pick.iterrows():
        s = cen[(r.wmo, r.ident)]
        info = V._info(s)
        o_e = _options({}, "ensemble")
        o_v = _options({}, "v2.3")
        files = [f for f in build_file_paths(r.date, info, o_e) if f.exists()]
        t0 = _t.perf_counter()
        native = load_data(files, info.instrument_type, o_e.data_level)
        t_read = _t.perf_counter() - t0
        t = {}
        for lab, o in (("ensemble", o_e), ("v2.3", o_v), ("ensemble2", o_e), ("v2.32", o_v)):
            t0 = _t.perf_counter()
            calibrate_rayleigh(r.date, info, o, preloaded_data=native)
            t[lab] = _t.perf_counter() - t0
        rows.append(dict(key=r.key, date=r.date, read=t_read,
                         ens=min(t["ensemble"], t["ensemble2"]), v23=min(t["v2.3"], t["v2.32"])))
    c = pd.DataFrame(rows)
    c.to_csv(OUT / "u23_cost.csv", index=False)
    print(f"{len(c)} CHM15k nights: L1 read median {c.read.median():.2f} s; calibration "
          f"ensemble {c.ens.median():.3f} s, v2.3 {c.v23.median():.3f} s "
          f"(+{(c.v23 - c.ens).median():.3f} s, +{100 * ((c.v23 - c.ens) / c.ens).median():.0f} % "
          f"of the calibration, +{100 * ((c.v23 - c.ens) / (c.ens + c.read)).median():.0f} % incl. "
          f"the read)")


# --------------------------------------------------------------------------------------------
# report: final configuration -> tables + figures
# --------------------------------------------------------------------------------------------
def shipped_cfg():
    from calibration.rayleigh.uncertainty import DEFAULT_PARAMS as P
    return dict(BASE_CFG, below_m=P["extend_below_m"], above_m=P.get("extend_above_m", 1e9),
                min_alt_m=P["min_altitude_m"], snr_min=P["snr_min"],
                D_km=P["window_displacement_m"] / 1000.0, resid=P["include_residual"],
                floor=P["noise_floor_rel"])


def _load_prepared():
    import pickle
    with open(OUT / "u23_prepared.pkl", "rb") as fh:
        d = pickle.load(fh)
    return d["df"], d["tw"]


def final_frame(cfg=None):
    """Every night with the v2.2 and v2.3 sigmas (relative, 1 sigma) and its pulls."""
    df, tw = _load_prepared()
    cfg = cfg or shipped_cfg()
    df = add_vertical(df, cfg)
    df.loc[~df.ok, "sig_vert"] = np.nan
    df["sig23"] = total_sigma(df)
    df["sig23_nomol"] = total_sigma(df, with_mol=False)
    df["sig22"] = df.sig_ens
    df["z22"] = df.anom / df.sig22
    df["z23"] = df.anom / df.sig23
    df["u23_rel"] = 2 * df.sig23
    return df, tw, cfg


def unit_table(df):
    import pandas as pd
    rows = []
    for key, g in df[df.ok & np.isfinite(df.anom)].groupby("key"):
        s22, s23 = pull_stats(g.z22), pull_stats(g.z23)
        rows.append(dict(key=key, site=g.site.iloc[0], type=g.type.iloc[0], n=len(g),
                         obs_sigma=100 * pull_stats(g.anom)["rstd"],
                         rep22=100 * g.sig22.median(), rep23=100 * g.sig23.median(),
                         vert=100 * g.sig_vert.median(), noise=100 * g.sig_noise.median(),
                         mol=100 * g.sig_mol.median(), ens=100 * g.sig_ens.median(),
                         g_med=100 * g.g.median(), g_abs=100 * g.g.abs().median(),
                         pull22=s22["rstd"], pull23=s23["rstd"], out22=s22["out2"],
                         out23=s23["out2"]))
    return pd.DataFrame(rows).sort_values("obs_sigma")


def cmd_report(args):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    df, tw, cfg = final_frame()
    print("configuration:", cfg)
    chm = df.ok & np.isfinite(df.anom) & (df.type == "CHM15k")
    pay = df.key == "0-20000-0-06610_A"
    lines = []

    def report(label, m):
        a, b = pull_stats(df.z22[m]), pull_stats(df.z23[m])
        txt = (f"{label:34s} n {a['n']:5d} | v2.2 pull rstd {a['rstd']:.2f} out2 {100 * a['out2']:.0f}% "
               f"| v2.3 pull rstd {b['rstd']:.2f} out2 {100 * b['out2']:.1f}% | reported sigma "
               f"{100 * df.sig22[m].median():.1f} -> {100 * df.sig23[m].median():.1f} %, observed "
               f"{100 * pull_stats(df.anom[m])['rstd']:.1f} %")
        print(txt)
        lines.append(txt)

    report("CHM15k (Payerne A excluded)", chm & ~pay)
    report("CHM15k Payerne A", chm & pay)
    for typ in ("CL61", "Mini-MPL"):
        report(typ, df.ok & np.isfinite(df.anom) & (df.type == typ))
    report("all units", df.ok & np.isfinite(df.anom))
    # twins (no molecular term: it cancels in the ratio)
    sig_nomol = dict(zip(zip(df.key, df.date), df.sig23_nomol))
    sig22 = dict(zip(zip(df.key, df.date), df.sig22))
    t22, z22t = twin_pull_stats(tw, sig22)
    t23, z23t = twin_pull_stats(tw, sig_nomol)
    txt = (f"twins ({len(tw)} pair-nights, {tw.pair.nunique()} pair labels, 3 sites): v2.2 pull rstd "
           f"{t22['rstd']:.2f} out2 {100 * t22['out2']:.0f}% | v2.3 rstd {t23['rstd']:.2f} out2 "
           f"{100 * t23['out2']:.1f}% | observed pair sigma {100 * pull_stats(tw.anom)['rstd']:.1f} %")
    print(txt)
    lines.append(txt)
    for site, g in tw.groupby("site"):
        _, za = twin_pull_stats(g, sig22)
        _, zb = twin_pull_stats(g, sig_nomol)
        a, b = pull_stats(za), pull_stats(zb)
        txt = (f"   {site}: n {a['n']} v2.2 {a['rstd']:.2f}/{100 * a['out2']:.0f}% -> v2.3 "
               f"{b['rstd']:.2f}/{100 * b['out2']:.1f}%; observed {100 * pull_stats(g.anom)['rstd']:.1f} %")
        print(txt)
        lines.append(txt)
    # how much of the vertical term is common to the twins (atmosphere) and must cancel in the ratio
    idx = df.set_index(["key", "date"])
    ga = idx.reindex(list(zip(tw.ka, tw.date)))
    gb = idx.reindex(list(zip(tw.kb, tw.date)))
    okg = np.isfinite(ga.g.values) & np.isfinite(gb.g.values)
    zca = 0.5 * (ga.zb_agl.values + ga.zt_agl.values)
    zcb = 0.5 * (gb.zb_agl.values + gb.zt_agl.values)
    okz = np.isfinite(zca) & np.isfinite(zcb)
    txt = (f"twin same-night correlation: gradient g r = "
           f"{np.corrcoef(ga.g.values[okg], gb.g.values[okg])[0, 1]:+.2f} (n {okg.sum()}), window "
           f"centre r = {np.corrcoef(zca[okz], zcb[okz])[0, 1]:+.2f}; median |g| "
           f"{100 * np.median(np.abs(ga.g.values[okg])):.1f} %/km, median |g_A - g_B| "
           f"{100 * np.median(np.abs(ga.g.values[okg] - gb.g.values[okg])):.1f} %/km")
    print(txt)
    lines.append(txt)
    ut = unit_table(df)
    ut.to_csv(OUT / "u23_units.csv", index=False)
    rho22 = spearmanr(ut.rep22, ut.obs_sigma)[0]
    rho23 = spearmanr(ut.rep23, ut.obs_sigma)[0]
    txt = f"Spearman(reported, observed) across {len(ut)} units: v2.2 {rho22:+.2f}, v2.3 {rho23:+.2f}"
    print(txt)
    lines.append(txt)
    ut["rank23"] = ut.rep23.rank(ascending=False).astype(int)
    ut["rank_obs"] = ut.obs_sigma.rank(ascending=False).astype(int)
    pd.set_option("display.width", 250)
    print(ut.sort_values("rep23", ascending=False).round(2).to_string(index=False))
    (OUT / "u23_report.txt").write_text("\n".join(lines), encoding="utf-8")
    df.drop(columns=["sub", "_z", "_r", "_s"]).to_csv(OUT / "u23_final_nights.csv", index=False)
    _figures(df, tw, ut, sig22, sig_nomol, plt)


def _figures(df, tw, ut, sig22, sig_nomol, plt):
    chm = df.ok & np.isfinite(df.anom) & (df.type == "CHM15k") & (df.key != "0-20000-0-06610_A")
    # 1 -- pulls ------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    bins = np.linspace(-8, 8, 81)
    xx = np.linspace(-8, 8, 400)
    gauss = np.exp(-0.5 * xx ** 2) / np.sqrt(2 * np.pi)
    for a, (lab, z22, z23) in zip(ax[:2], (
            ("CHM15k, 31-day running median (Payerne excluded)", df.z22[chm], df.z23[chm]),
            ("same-night twins (Schiphol, Budapest, Magurele)",
             twin_pull_stats(tw, sig22)[1], twin_pull_stats(tw, sig_nomol)[1]))):
        s22, s23 = pull_stats(z22), pull_stats(z23)
        a.hist(np.clip(z22, -8, 8), bins=bins, density=True, histtype="step", color="0.45", lw=1.6,
               label=f"v2.2 ensemble: robust sd {s22['rstd']:.2f}, {100 * s22['out2']:.0f} % > 2 sigma")
        a.hist(np.clip(z23, -8, 8), bins=bins, density=True, histtype="stepfilled", alpha=0.45,
               color="tab:blue",
               label=f"v2.3 budget: robust sd {s23['rstd']:.2f}, {100 * s23['out2']:.1f} % > 2 sigma")
        a.plot(xx, gauss, "k--", lw=1, label="N(0, 1)")
        a.set_xlabel("pull = anomaly / reported sigma")
        a.set_ylabel("density")
        a.set_title(lab, fontsize=10)
        a.legend(fontsize=8, loc="upper left")
        a.grid(alpha=0.3)
    a = ax[2]
    for typ, mk in (("CHM15k", "o"), ("CL61", "s"), ("Mini-MPL", "^")):
        u = ut[ut.type == typ]
        a.plot(u.rep22, u.obs_sigma, mk, mfc="none", color="0.5", ms=6)
        a.plot(u.rep23, u.obs_sigma, mk, color="tab:blue", ms=6, label=f"{typ} v2.3")
        for _, r in u.iterrows():
            a.plot([r.rep22, r.rep23], [r.obs_sigma, r.obs_sigma], "-", color="0.8", lw=0.8, zorder=0)
            if r.key in SPECIAL:
                a.annotate(SPECIAL[r.key], (r.rep23, r.obs_sigma), xytext=(4, 3),
                           textcoords="offset points", fontsize=8)
    lim = [0.5, 60]
    a.plot(lim, lim, "k--", lw=1, label="1:1")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel("reported sigma, unit median [%] (grey: v2.2, blue: v2.3)")
    a.set_ylabel("observed night-to-night sigma [%]")
    a.set_title("per unit: reported vs observed", fontsize=10)
    a.legend(fontsize=8, loc="upper left")
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "u23_pulls.png", dpi=130)
    plt.close(fig)

    # 2 -- component breakdown per unit ---------------------------------------------------------
    u = ut.sort_values("rep23")
    fig, ax = plt.subplots(figsize=(17, 6.5))
    x = np.arange(len(u))
    bottom = np.zeros(len(u))
    for col, c, lab in (("ens", "0.6", "ensemble (v2.2)"), ("vert", "tab:red", "vertical"),
                        ("mol", "tab:green", "molecular"), ("noise", "tab:purple", "photon noise")):
        v = u[col].values ** 2
        ax.bar(x, v, bottom=bottom, color=c, label=lab, width=0.8)
        bottom += v
    ax.plot(x, u.obs_sigma.values ** 2, "k_", ms=14, mew=2, label="observed variance")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.site[:14]} {r.key.rsplit('_', 1)[1]} ({r.type[:4]})"
                        for _, r in u.iterrows()], rotation=70, ha="right", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylabel("variance, unit median [%^2] (stacked)")
    ax.set_title("v2.3 budget per unit (median of the nightly components) and the observed "
                 "night-to-night variance")
    ax.legend(fontsize=8, ncol=5)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "u23_components.png", dpi=130)
    plt.close(fig)

    # 3 -- QC: anomaly vs relative uncertainty --------------------------------------------------
    m = df.ok & np.isfinite(df.anom)
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.2))
    a = ax[0]
    a.plot(100 * df.u23_rel[m], 100 * np.abs(df.anom[m]), ".", ms=2, alpha=0.3, color="tab:blue")
    a.plot([0, 300], [0, 300], "k--", lw=1, label="|anomaly| = 2 sigma (reported)")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(2, 300)
    a.set_ylim(0.01, 300)
    a.set_xlabel("v2.3 reported u / C_L (2 sigma) [%]")
    a.set_ylabel("|anomaly| from the 31-day running median [%]")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    a = ax[1]
    thr = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0])
    frac_rej, out_kept, out_rej = [], [], []
    for t in thr:
        rej = df.u23_rel[m] > t
        frac_rej.append(100 * rej.mean())
        big = np.abs(df.anom[m]) > np.log(1.25)
        out_kept.append(100 * big[~rej].mean())
        out_rej.append(100 * big[rej].mean() if rej.any() else np.nan)
    a.plot(100 * thr, frac_rej, "o-", label="nights rejected")
    a.plot(100 * thr, out_kept, "s-", label="kept nights off by > 25 %")
    a.plot(100 * thr, out_rej, "^-", label="rejected nights off by > 25 %")
    a.set_xscale("log")
    a.set_xlabel("QC limit on u / C_L (2 sigma) [%]")
    a.set_ylabel("[%]")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "u23_qc.png", dpi=130)
    plt.close(fig)
    qc = pd_qc(df[m], thr)
    qc.to_csv(OUT / "u23_qc.csv", index=False)
    print(qc.round(3).to_string(index=False))


def pd_qc(d, thr):
    import pandas as pd
    rows = []
    for t in thr:
        rej = d.u23_rel > t
        big = np.abs(d.anom) > np.log(1.25)
        kept = d[~rej]
        rows.append(dict(limit=t, frac_rejected=rej.mean(), n_rejected=int(rej.sum()),
                         outlier_rate_kept=big[~rej].mean(),
                         outlier_rate_rejected=big[rej].mean() if rej.any() else np.nan,
                         obs_sigma_kept=pull_stats(kept.anom)["rstd"],
                         pull23_kept=pull_stats(kept.z23)["rstd"],
                         frac_payerne_rejected=float(rej[d.key == "0-20000-0-06610_A"].mean())))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------
# Window displacement statistics (justification of D)
# --------------------------------------------------------------------------------------------
def cmd_windows(args):
    """How far does the molecular window move between nights? Per unit, from the PUBLISHED
    series (all nights; window centre = middle of calibration_bottom/top_height)."""
    import pandas as pd
    from common import load_constants
    cen = census()
    rows = []
    for wmo, ident in UNITS:
        s = cen[(wmo, ident)]
        con = load_constants(wmo, ident)
        if con is None:
            continue
        ok = np.isfinite(con["zb"]) & np.isfinite(con["zt"]) & (con["c"] > 0)
        t = con["t"][ok]
        zc = (0.5 * (con["zb"] + con["zt"])[ok] - s.get("alt", 0.0)) / 1000.0
        if zc.size < 30:
            continue
        o = np.argsort(t)
        t, zc = t[o], zc[o]
        rm = running_median(t, zc)
        doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t])
        summer = (doy >= 152) & (doy <= 243)
        winter = (doy <= 59) | (doy >= 335)
        rows.append(dict(
            key=f"{wmo}_{ident}", site=s["site"], type=s["type"], n=int(zc.size),
            zc_median=float(np.median(zc)),
            zc_sd=float(np.std(zc)),
            zc_rsd=float(1.4826 * np.median(np.abs(zc - np.median(zc)))),
            within_month_sd=float(np.nanstd(zc - rm)),
            consecutive_abs=float(np.median(np.abs(np.diff(zc)))),
            summer_minus_winter=(float(np.median(zc[summer]) - np.median(zc[winter]))
                                 if summer.sum() > 5 and winter.sum() > 5 else np.nan),
            p10_p90=float(np.percentile(zc, 90) - np.percentile(zc, 10))))
    w = pd.DataFrame(rows)
    w.to_csv(OUT / "u23_window_displacement.csv", index=False)
    pd.set_option("display.width", 220)
    print(w.round(2).to_string(index=False))
    for col in ("zc_sd", "zc_rsd", "within_month_sd", "consecutive_abs", "summer_minus_winter",
                "p10_p90"):
        print(f"{col:22s} network median {w[col].median():.2f} km [IQR {w[col].quantile(0.25):.2f}"
              f"-{w[col].quantile(0.75):.2f}]")


# --------------------------------------------------------------------------------------------
# Kalman: weighting the best estimate by the uncertainty
# --------------------------------------------------------------------------------------------
def _kalman_series(times, values, unc):
    from monitoring.kalman import kalman_best_estimate
    g, s, sd = kalman_best_estimate(times, values, uncertainties=unc)
    days = np.array([(np.datetime64(x, "D") - np.datetime64("1970-01-01", "D")).astype(int)
                     for x in g])
    return days, np.asarray(s, float)


def _innovations(days, state, t_obs, lnc):
    """One-step-ahead innovations: ln C(t) - ln state(t - 1 day) (forward filter: the state of the
    previous day has not seen night t)."""
    idx = {int(d): i for i, d in enumerate(days)}
    e = np.full(len(t_obs), np.nan)
    for k, (t, y) in enumerate(zip(t_obs, lnc)):
        i = idx.get(int(t) - 1)
        if i is not None and i > 10 and np.isfinite(state[i]) and state[i] > 0:
            e[k] = y - np.log(state[i])
    return e


def cmd_kalman(args):
    import pandas as pd
    df, tw, cfg = final_frame()
    ok = df[df.ok].copy()
    rows = []
    series = {}
    for key, g in ok.groupby("key"):
        g = g.sort_values("tday")
        if len(g) < 30:
            continue
        times = [EPOCH + timedelta(days=int(x)) for x in g.tday]
        vals = g.c.values
        out = {}
        for lab, unc in (("flat", None), ("w22", g.u_ens.values),
                         ("w23", 2 * g.sig23.values * vals)):
            days, st = _kalman_series(times, vals, unc)
            out[lab] = (days, st)
            e = _innovations(days, st, g.tday.values, np.log(vals))
            good = (g.sig23.values <= np.nanmedian(g.sig23.values))   # the more certain half
            rows.append(dict(key=key, site=g.site.iloc[0], type=g.type.iloc[0], model=lab,
                             n=len(g), innov_all=pull_stats(e)["rstd"],
                             innov_good=pull_stats(e[good])["rstd"],
                             innov_rms_good=float(np.sqrt(np.nanmean(e[good] ** 2))),
                             jitter=float(np.nanstd(np.diff(np.log(st))))))
        series[key] = out
    k = pd.DataFrame(rows)
    k.to_csv(OUT / "u23_kalman_units.csv", index=False)
    piv = k.pivot_table(index="key", columns="model", values=["innov_good", "innov_all", "jitter"])
    for met in ("innov_good", "innov_all", "jitter"):
        d23 = 100 * (piv[(met, "w23")] / piv[(met, "flat")] - 1)
        d22 = 100 * (piv[(met, "w22")] / piv[(met, "flat")] - 1)
        print(f"{met:11s}: weighted(v2.3) vs flat median {d23.median():+.1f} % "
              f"[IQR {d23.quantile(0.25):+.1f}..{d23.quantile(0.75):+.1f}], better on "
              f"{100 * np.mean(d23 < 0):.0f} % of units; weighted(v2.2) {d22.median():+.1f} %")
    pay = piv.loc["0-20000-0-06610_A"] if "0-20000-0-06610_A" in piv.index else None
    if pay is not None:
        print("Payerne A:", pay.round(4).to_dict())
    # twins: stability of the Kalman ratio (deviation from its own 31-day running median)
    trows = []
    for wmo, idents in TWINS.items():
        for i, a in enumerate(idents):
            for b in idents[i + 1:]:
                ka, kb = f"{wmo}_{a}", f"{wmo}_{b}"
                if ka not in series or kb not in series:
                    continue
                for lab in ("flat", "w22", "w23"):
                    da, sa = series[ka][lab]
                    db, sb = series[kb][lab]
                    common, ia, ib = np.intersect1d(da, db, return_indices=True)
                    if common.size < 60:
                        continue
                    r = np.log(sa[ia] / sb[ib])
                    dev = r - running_median(common.astype(float), r)
                    trows.append(dict(pair=f"{ka}/{kb}", model=lab,
                                      ratio_dev_rstd=pull_stats(dev)["rstd"],
                                      ratio_dev_p95=float(np.nanpercentile(np.abs(dev), 95))))
    t = pd.DataFrame(trows)
    t.to_csv(OUT / "u23_kalman_twins.csv", index=False)
    tp = t.pivot_table(index="pair", columns="model", values=["ratio_dev_rstd", "ratio_dev_p95"])
    print((100 * tp).round(2).to_string())
    for met in ("ratio_dev_rstd", "ratio_dev_p95"):
        d = 100 * (tp[(met, "w23")] / tp[(met, "flat")] - 1)
        print(f"twins {met}: v2.3-weighted vs flat median {d.median():+.1f} %, better on "
              f"{100 * np.mean(d < 0):.0f} % of pairs")
    _kalman_figure(series, ok)


def _kalman_figure(series, ok):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [k for k in ("0-20000-0-06610_A", "0-20000-0-10564_0", "0-20000-0-06240_B")
            if k in series]
    fig, ax = plt.subplots(len(keys), 1, figsize=(17, 3.6 * len(keys)), squeeze=False)
    for a, key in zip(ax[:, 0], keys):
        g = ok[ok.key == key].sort_values("tday")
        t = np.array([np.datetime64(EPOCH + timedelta(days=int(x))) for x in g.tday])
        scale = np.median(g.c)
        a.errorbar(t, g.c / scale, yerr=g.sig23 * g.c / scale, fmt="o", ms=3, color="0.55",
                   ecolor="0.8", lw=0.8, label="nightly C_L +- 1 sigma (v2.3)")
        for lab, col, txt in (("flat", "k", "Kalman, equal weights (operational)"),
                              ("w23", "tab:blue", "Kalman, weighted by the v2.3 uncertainty")):
            d, s = series[key][lab]
            a.plot(d.astype("datetime64[D]"), s / scale, "-", color=col, lw=1.6, label=txt)
        a.set_ylim(0.3, 2.0)
        a.set_ylabel("C_L / unit median")
        a.set_title(f"{g.site.iloc[0]} ({key})", fontsize=10)
        a.grid(alpha=0.3)
        a.legend(fontsize=8, ncol=3, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "u23_kalman.png", dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------------------------
# Example sub-window profiles (altitude on the Y axis)
# --------------------------------------------------------------------------------------------
def cmd_profiles(args):
    import pickle
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df, tw, cfg = final_frame()
    ok = df[df.ok & np.isfinite(df.sig23)]
    panels = []
    # (1) a clean night of a twin unit: small gradient, many sub-windows
    cand = ok[(ok.key == "0-20000-0-06240_B") & (ok.nsub >= 6)]
    if len(cand):
        panels.append(("Schiphol B, a typical night", cand.iloc[(cand.g.abs()).argsort().iloc[len(cand) // 2]]))
    # (2) a typical Payerne A night (median |g|)
    cand = ok[(ok.key == "0-20000-0-06610_A") & (ok.nsub >= 4)]
    if len(cand):
        panels.append(("Payerne A, a typical night",
                       cand.iloc[(cand.g.abs() - cand.g.abs().median()).abs().argsort().iloc[0]]))
    # (3) Schleiz, its largest-gradient night
    cand = ok[(ok.key == "0-20000-0-10564_0") & (ok.nsub >= 4)]
    if len(cand):
        panels.append(("Schleiz, steepest night", cand.loc[cand.g.abs().idxmax()]))
    pay = None
    pf = OUT / "u23_payerne.pkl"
    if pf.exists():
        with open(pf, "rb") as fh:
            pay = pickle.load(fh).get("20260907")
    fig, ax = plt.subplots(1, len(panels) + (1 if pay else 0), figsize=(18, 5.6), sharex=False)
    ax = np.atleast_1d(ax)
    for a, (title, r) in zip(ax, panels):
        z, rr, s = r["_z"], r["_r"], r["_s"]
        zz, r2, s2 = select_subwindows(z, rr, s, r.zb_agl, r.zt_agl, cfg)
        a.errorbar(rr, z / 1000, xerr=np.nan_to_num(s), fmt="o", color="0.75", ms=4,
                   label="all recorded sub-windows")
        a.errorbar(r2, zz / 1000, xerr=np.nan_to_num(s2), fmt="o", color="tab:blue", ms=5,
                   label="used for the vertical term")
        zc = 0.5 * (r.zb_agl + r.zt_agl)
        zline = np.linspace(zz.min(), zz.max(), 20) if zz.size else np.array([zc])
        if np.isfinite(r.g):
            a.plot(np.exp(r.g * (zline - zc) / 1000.0) * np.exp(np.mean(np.log(r2)) - r.g *
                                                                 (np.mean(zz) - zc) / 1000.0),
                   zline / 1000, "-", color="tab:red", label=f"fit: {100 * r.g:+.1f} %/km")
        a.axhspan(r.zb_agl / 1000, r.zt_agl / 1000, color="tab:green", alpha=0.15,
                  label="operational window")
        a.axvline(1.0, color="k", lw=0.8)
        a.set_xlabel("C_L(sub-window) / C_L(operational)")
        a.set_ylabel("range [km]")
        a.set_title(f"{title}\n{r.date}: v2.2 2s {200 * r.sig22:.0f} % -> v2.3 2s {200 * r.sig23:.0f} %",
                    fontsize=9)
        a.legend(fontsize=7, loc="lower right")
        a.grid(alpha=0.3)
    if pay:
        a = ax[-1]
        rec = pay.get(("record", True))
        dfl = pay.get(("default", True))
        b = rec["budget"] if rec else None
        if b is not None and b.sub_cl.size:
            c0 = b.lidar_constant
            a.errorbar(b.sub_cl / c0, b.sub_centre_m / 1000, xerr=b.sub_sigma_rel, fmt="o",
                       color="0.75", ms=4, label="all recorded sub-windows")
        bd = dfl["budget"] if dfl else None
        if bd is not None and bd.sub_cl.size:
            c0 = bd.lidar_constant
            a.errorbar(bd.sub_cl / c0, bd.sub_centre_m / 1000, xerr=bd.sub_sigma_rel, fmt="o",
                       color="tab:blue", ms=5, label="used (shipped defaults)")
            zc = 0.5 * (dfl["zb"] + dfl["zt"]) - 490.0
            zline = np.linspace(bd.sub_centre_m.min(), bd.sub_centre_m.max(), 20)
            a.plot(np.exp(bd.gradient_per_km * (zline - zc) / 1000.0), zline / 1000, "-",
                   color="tab:red", label=f"fit: {100 * bd.gradient_per_km:+.0f} %/km")
            a.axhspan((dfl["zb"] - 490) / 1000, (dfl["zt"] - 490) / 1000, color="tab:green",
                      alpha=0.15, label="operational window")
            a.set_title(f"Payerne A, 2026-09-07 (classification screen on)\nv2.2 2s "
                        f"{100 * bd.u_ensemble / c0:.0f} % -> v2.3 2s {100 * bd.rel_total:.0f} %",
                        fontsize=9)
        a.axvline(1.0, color="k", lw=0.8)
        a.set_xlabel("C_L(sub-window) / C_L(operational)")
        a.set_ylabel("range [km]")
        a.legend(fontsize=7, loc="lower right")
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "u23_profiles.png", dpi=130)
    plt.close(fig)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    a = ap.parse_args()
    globals()[f"cmd_{a.cmd}"](a)


if __name__ == "__main__":
    main()
