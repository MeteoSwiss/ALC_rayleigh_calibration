# -*- coding: utf-8 -*-
"""Review attack (claim 2): do "the two truths bracket the answer", and is the "thirds" argument sound?

1. Variance shares from the PUBLISHED series (kappa-corrected 31-day anomalies), with bootstrap CIs:
   twin pairs (shared = atmosphere + CHM15k gain) and co-located CHM15k / CL61 pairs (shared =
   atmosphere only, scaled by the CL61's relative sensitivity w in [0.63, 1]: an aerosol error scales
   as (910/1064)^3, a molecular-model error does not). The fresh Camborne pair is out of sample.
2. The pulls a perfectly calibrated budget would give under those shares (time series and twins).
3. A sharper twin test that needs no variance decomposition: a ratio-specific budget built only from
   the terms that do NOT cancel in ln(C_A/C_B) -- the photon noise and the time-subset spread of
   both units, the DIFFERENTIAL gradient (g_A - g_B) x D, the common gradient times the difference of
   the two windows ((g_A + g_B)/2 x (z_A - z_B)), and optionally the residual scatter. If the
   budget's unit-specific part is right this pull is ~1.

Usage: python u23r_bracket.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
import u23_analysis as A  # noqa: E402
from u23r_truth import prepare  # noqa: E402
from u23r_components import rsd5  # noqa: E402

warnings.filterwarnings("ignore")

TWIN_PAIRS = [(f"{w}_{a}", f"{w}_{b}") for w, ids in C.TWINS.items()
              for i, a in enumerate(ids) for b in ids[i + 1:]]
CL61_PAIRS = [("0-20000-0-10393_0", "0-20000-0-10393_C"), ("0-380-5-1_0", "0-380-5-1_B"),
              ("0-20000-0-06610_A", "0-20000-0-06610_C"), ("0-20000-0-03808_A", "0-20000-0-03808_C")]


def pub_anoms(T, key):
    P = T.pub(key)
    if not P.ok:
        return {}
    typ = T.cen[tuple(key.rsplit("_", 1))]["type"]
    lncc = P.corrected(T.kappa(key, typ))
    base = A.running_median(P.t, lncc, exclude_self=True)
    return dict(zip(P.dates, lncc - base))


def shares(xa, xb, n_boot=1000, seed=0):
    """Robust variance of a, b and their covariance (4-sigma clipped), with bootstrap CIs."""
    def stats(x, y):
        def clip(v):
            s = C.rstd(v)
            return np.abs(v - np.median(v)) <= 4 * s
        m = clip(x) & clip(y)
        x, y = x[m], y[m]
        va, vb = np.var(x), np.var(y)
        cov = np.cov(x, y)[0, 1]
        return va, vb, cov, cov / np.sqrt(va * vb)
    est = stats(xa, xb)
    rs = np.random.default_rng(seed)
    bs = np.array([stats(xa[ii], xb[ii]) for ii in (rs.integers(0, xa.size, xa.size)
                                                     for _ in range(n_boot))])
    return est, bs


def main():
    import pandas as pd
    T = C.Truth()
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    # ---- 1. variance shares ------------------------------------------------------------------
    rows = []
    cache = {}
    for kind, pairs in (("twin", TWIN_PAIRS), ("CHM15k/CL61", CL61_PAIRS)):
        for a, b in pairs:
            for k in (a, b):
                if k not in cache:
                    cache[k] = pub_anoms(T, k)
            common = sorted(set(cache[a]) & set(cache[b]))
            xa = np.array([cache[a][d] for d in common])
            xb = np.array([cache[b][d] for d in common])
            ok = np.isfinite(xa) & np.isfinite(xb)
            xa, xb = xa[ok], xb[ok]
            if xa.size < 20:
                continue
            (va, vb, cov, r), bs = shares(xa, xb)
            rows.append(dict(kind=kind, a=a, b=b, n=xa.size, sigma_a=100 * np.sqrt(va),
                             sigma_b=100 * np.sqrt(vb), r=r, r_lo=np.percentile(bs[:, 3], 2.5),
                             r_hi=np.percentile(bs[:, 3], 97.5), share_a=cov / va,
                             share_a_lo=np.percentile(bs[:, 2] / bs[:, 0], 2.5),
                             share_a_hi=np.percentile(bs[:, 2] / bs[:, 0], 97.5),
                             fresh=(a == "0-20000-0-03808_A")))
    sh = pd.DataFrame(rows)
    sh.to_csv(C.OUT / "u23r_bracket_shares.csv", index=False)
    pd.set_option("display.width", 250)
    say(sh.round(3).to_string(index=False))

    # ---- 2. pulls a perfect budget would give -------------------------------------------------
    st = sh[sh.kind == "twin"].share_a.values
    sc = sh[(sh.kind == "CHM15k/CL61") & (sh.a != "0-20000-0-06610_A")].share_a.values
    X2 = 0.0135 ** 2 + 0.010 ** 2 + 0.02 ** 2            # molecular + aerosol terms (cancel in twins)
    V = 0.110 ** 2                                         # CHM15k observed night-to-night variance
    grid = []
    for s_t in (np.min(st), np.median(st), np.max(st)):
        for s_c in (np.min(sc), np.median(sc), np.max(sc)):
            for w in (0.63, 1.0):
                v_atm = min(s_c / w, s_t)                  # cannot exceed what twins share
                g = s_t - v_atm
                v_unit = 1 - s_t
                e2 = v_atm + v_unit
                ts = np.sqrt(1.0 / e2)
                tw = np.sqrt(v_unit / max(e2 - X2 / V, 1e-6))
                grid.append(dict(share_twin=s_t, share_cl61=s_c, w=w, atm=v_atm, gain=g,
                                 unit=v_unit, ts_pull=ts, twin_pull=tw))
    gr = pd.DataFrame(grid)
    gr.to_csv(C.OUT / "u23r_bracket_expected.csv", index=False)
    say(gr.round(3).to_string(index=False))
    say(f"perfect-budget expectation over the share CIs: time-series pull {gr.ts_pull.min():.2f}-"
        f"{gr.ts_pull.max():.2f}, twin pull {gr.twin_pull.min():.2f}-{gr.twin_pull.max():.2f}")

    # ---- 3. ratio-specific twin budget ---------------------------------------------------------
    d = prepare("orig", "eprof_v2.2")
    df, tw = d["df"], d["tw"]
    s = df.set_index(["key", "date"])
    s_time = {k: rsd5(v) for k, v in zip(zip(df.key, df.date), df._ens_time)}
    out = []
    for D in (0.75, 1.0, 1.25):
        for with_res in (True, False):
            sig = []
            for r in tw.itertuples():
                A_, B_ = s.loc[(r.ka, r.date)], s.loc[(r.kb, r.date)]
                v = (A_.def_sig_noise ** 2 + B_.def_sig_noise ** 2
                     + np.nan_to_num(s_time.get((r.ka, r.date), np.nan)) ** 2
                     + np.nan_to_num(s_time.get((r.kb, r.date), np.nan)) ** 2)
                ga, gb = A_.def_g, B_.def_g
                if np.isfinite(ga) and np.isfinite(gb):
                    v += (D * (ga - gb)) ** 2
                    v += (0.5 * (ga + gb) * (A_.zc_agl - B_.zc_agl) / 1000.0) ** 2
                if with_res:
                    v += np.nan_to_num(A_.def_res) ** 2 + np.nan_to_num(B_.def_res) ** 2
                sig.append(np.sqrt(v))
            sig = np.array(sig)
            z = tw.anom.values / sig
            ps = C.pull_stats(z)
            zc = tw.anom.values / np.sqrt(sig ** 2 + tw.anom_se.values ** 2)
            out.append(dict(D_km=D, residual=with_res, n=ps["n"], twin_pull=ps["rstd"],
                            out2=ps["out2"], pull_truth_err_removed=C.pull_stats(zc)["rstd"],
                            med_sig_ratio=100 * np.nanmedian(sig)))
    rt = pd.DataFrame(out)
    rt.to_csv(C.OUT / "u23r_bracket_ratio_budget.csv", index=False)
    say(rt.round(3).to_string(index=False))
    s23 = C.twin_sigma(tw, df, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
    say(f"shipped twin denominator: median {100 * np.nanmedian(s23):.2f} %, pull "
        f"{C.pull_stats(tw.anom / s23)['rstd']:.3f}; observed twin anomaly robust sigma "
        f"{100 * C.rstd(tw.anom):.2f} %, median truth error {100 * np.nanmedian(tw.anom_se):.2f} %")
    (C.OUT / "u23r_bracket.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
