# -*- coding: utf-8 -*-
"""Review attack: WHY is the budget not calibrated as a function of the reported sigma?

(1) Twin truth in bins of the reported pair sigma (real gain changes common to co-located units cancel
    there): if the quiet end is fine on the twins, the time-series floor is a real gain change.
(2) Time-series truth in bins of |g| (the vertical term): does the observed scatter grow like |g| x D?
(3) Is the vertical error random (|g| x D) or set by the night's own window displacement? The anomaly
    is regressed on g x dz, dz = the night's window centre minus the 31-day running median of the
    unit's PUBLISHED window centres (the windows the truth's median is built from).

Usage: python u23r_reliability2.py [orig|fresh] [eprof_v2.2|eprof_v2.3]
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

warnings.filterwarnings("ignore")


def window_dz(d, T):
    """Night window centre minus the 31-day running median of the published window centres (km)."""
    out = np.full(len(d), np.nan)
    for key, g in d.groupby("key"):
        P = T.pub(key)
        if not P.ok:
            continue
        wmo, ident = key.rsplit("_", 1)
        from common import load_constants
        con = load_constants(wmo, ident)
        alt = T.cen[(wmo, ident)].get("alt", 0.0)
        zc = 0.5 * (con["zb"] + con["zt"]) - alt
        dates = A.Published.__dict__  # noqa: F841 - (keeps the import local)
        import q16a_recalibration as Q
        nd = np.array(Q.night_dates_from_constants(con, T.cen[(wmo, ident)]["lon"]))
        t = np.array([(np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:]}") - np.datetime64("1970-01-01"))
                      .astype(int) for x in nd], float)
        ok = np.isfinite(zc)
        med, _, _ = C.running_median_ex(t[ok], zc[ok], g.tday.values.astype(float), exclude_days=0.5)
        out[d.index.get_indexer(g.index)] = (g.zc_agl.values - med) / 1000.0
    return out


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    p = prepare(set_name, method)
    df = p["df"]
    T = C.Truth()
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    # (1) twins by pair sigma
    tw = p.get("tw")
    if tw is not None and len(tw):
        s23 = C.twin_sigma(tw, df, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
        t = tw.assign(s=s23)
        t = t[np.isfinite(t.s)]
        t["q"] = pd.qcut(t.s, 5, labels=False)
        say("twins, quintiles of the reported pair sigma (v2.3, molecular+aerosol removed):")
        for q, g in t.groupby("q"):
            z = g.anom / g.s
            say(f"   q{q}: n {len(g):3d} reported pair sigma {100 * g.s.median():5.1f} % observed "
                f"{100 * C.rstd(g.anom):5.1f} % (truth error {100 * g.anom_se.median():4.1f} %) pull "
                f"{C.pull_stats(z)['rstd']:.2f} outside 2 sigma {100 * C.pull_stats(z)['out2']:4.1f} %")

    d = df[df.ok & np.isfinite(df.anom) & np.isfinite(df.sig23) & (df.type == "CHM15k")
           & (df.key != C.PAYERNE)].copy()
    # (2) by |g|
    d["gq"] = pd.qcut(np.abs(d.def_g), 8, labels=False)
    say("time series, octiles of |g| (CHM15k, Payerne excl.):")
    for q, g in d.groupby("gq"):
        say(f"   q{q}: n {len(g):3d} median |g| {100 * np.abs(g.def_g).median():5.1f} %/km vertical "
            f"{100 * g.def_sig_vert.median():5.1f} % total {100 * g.sig23.median():5.1f} % observed "
            f"{100 * C.rstd(g.anom):5.1f} % pull {C.pull_stats(g.anom / g.sig23)['rstd']:.2f}")
    # (3) displacement
    d["dz"] = window_dz(d, T)
    ok = np.isfinite(d.dz) & np.isfinite(d.def_g)
    x = (d.def_g * d.dz)[ok].values
    y = d.anom[ok].values
    clip = np.abs(y - np.median(y)) < 5 * C.rstd(y)
    slope = np.polyfit(x[clip], y[clip], 1)[0]
    from scipy.stats import spearmanr
    say(f"anomaly vs g x dz (n {int(ok.sum())}): OLS slope {slope:+.2f} (1 = the window displacement "
        f"explains the error fully), Spearman {spearmanr(x, y)[0]:+.3f}; corr(anomaly, dz) "
        f"{spearmanr(d.dz[ok], y)[0]:+.3f}; median |dz| {np.median(np.abs(d.dz[ok])):.2f} km")
    for lab, m in (("|g| x |dz| small", np.abs(x) < np.percentile(np.abs(x), 50)),
                   ("|g| x |dz| large", np.abs(x) >= np.percentile(np.abs(x), 80))):
        dd = d[ok][m]
        say(f"   {lab}: n {len(dd)} observed {100 * C.rstd(dd.anom):.1f} % reported {100 * dd.sig23.median():.1f} %"
            f" pull {C.pull_stats(dd.anom / dd.sig23)['rstd']:.2f}")
    d[["key", "date", "anom", "sig23", "def_g", "dz", "def_sig_vert"]].to_csv(
        C.OUT / f"u23r_reliability2_{set_name}_{method}.csv", index=False)
    (C.OUT / f"u23r_reliability2_{set_name}_{method}.txt").write_text("\n".join(lines),
                                                                         encoding="utf-8")


if __name__ == "__main__":
    main()
