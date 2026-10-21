# -*- coding: utf-8 -*-
"""ref23: what is the uncertainty the O'Connor cloud calibration reports, and is it meaningful?

The cloud runner (scripts/run_network_calibration._do_cloud) publishes for a successful day
    lidar_constant_uncertainty = C_L x std(C_i) / median(C_i)
where C_i are the per-profile O'Connor coefficients of the day's accepted liquid-cloud profiles
(calibration.cloud.calibration: cal_std = np.std(valid_C, ddof=1)).  So it is the 1-sigma SPREAD of
single-profile values -- not the uncertainty of the day's median (which is what is used), and,
unlike the Rayleigh product, a 1-sigma (not 2-sigma) number.  It contains no systematic term
(multiple-scattering eta, S = 18.8 sr, water vapour, aerosol transmission).

The test here: the day-to-day reproducibility of the cloud constant.  For a channel the structure
function D(L) = robust variance of ln C(d + L) - ln C(d) over pairs of successful days L days apart
is 2 sigma_w^2 + drift(L); sigma_w = sqrt(D(1)/2) is the actual random error of one day value (an
upper bound: it also contains real day-to-day instrument change).  It is compared with
  * sigma_rep  = std / median                         (the published number)
  * sigma_sem  = 1.2533 x robust spread / sqrt(n_scenes)   (standard error of the median, counting
                independent cloud scenes -- consecutive accepted profiles are not independent)

Writes cloud_uncertainty.csv, cloud_structure.csv, cloud_vs_published.csv and the figure
ref23_cloud_uncertainty.png; returns the per-channel model used by the pair analysis.
"""
from __future__ import annotations

import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402
import ref23_pairs as P  # noqa: E402

CHANNELS = [("payerne", "C"), ("lindenberg", "C"), ("aosta", "B"), ("payerne", "B")]
LABEL = {("payerne", "C"): "Payerne CL61", ("lindenberg", "C"): "Lindenberg CL61",
         ("aosta", "B"): "Aosta CL61", ("payerne", "B"): "Payerne CL31"}
MAX_LAG = 21


def _robust_var(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    return float((1.4826 * np.median(np.abs(x - np.median(x)))) ** 2)


def structure_function(cd, max_lag=MAX_LAG):
    """Robust D(L) of ln C between successful days L apart (same era only)."""
    dn = np.array([(datetime.strptime(d, "%Y%m%d") - RC.EPOCH).days for d in cd.index])
    y = cd.lnc.values
    era = cd.era.values
    out = []
    for L in range(1, max_lag + 1):
        i, j = np.nonzero((dn[None, :] - dn[:, None] == L) & (era[None, :] == era[:, None]))
        dif = y[j] - y[i]
        out.append(dict(lag=L, n=int(dif.size), D=_robust_var(dif)))
    return pd.DataFrame(out)


FIT_LAGS = 7          # the bridging range: D(L) is modelled over L = 1..7 days


def channel_model(cd, sf):
    """Noise + drift model of one channel's cloud constant: D(L) = 2 sigma_w^2 + q L, fitted by
    weighted least squares (weights = number of day pairs) over L = 1..FIT_LAGS, with q >= 0.
    sigma_w = the random error of one day value; q L = the variance of the real change of the
    constant over L days (what a bridged value misses).  drift(L) = sqrt(q L)."""
    g = sf[(sf.lag <= FIT_LAGS) & np.isfinite(sf.D) & (sf.n >= 5)]
    if len(g) < 3:
        d1 = float(np.nanmin(sf.D))
        return float(np.sqrt(d1 / 2.0)), 0.0
    L, D, w = g.lag.values.astype(float), g.D.values, g.n.values.astype(float)
    X = np.column_stack([np.full(L.size, 2.0), L])
    W = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * W[:, None], D * W, rcond=None)
    s2, q = float(coef[0]), float(coef[1])
    if q < 0:                                  # no measurable drift: constant + noise
        q = 0.0
        s2 = float(np.sum(w * D) / np.sum(w) / 2.0)
    return float(np.sqrt(max(s2, 0.0))), q


def run(make_figure=True):
    warnings.filterwarnings("ignore")
    rays, clouds, pub, _ing = P.load_tables()
    rows, sfs, models = [], [], {}
    cmp_rows = []
    for site, ident in CHANNELS:
        cd = P.cloud_days(clouds, site, ident)
        if len(cd) < 10:
            continue
        sf = structure_function(cd)
        sf["site"], sf["ident"] = site, ident
        sfs.append(sf)
        sig_w, q = channel_model(cd, sf)

        def drift(L, q=q):
            return float(np.sqrt(q * max(float(L), 0.0)))
        # does the reported spread predict the day-to-day difference? consecutive days only
        dn = np.array([(datetime.strptime(d, "%Y%m%d") - RC.EPOCH).days for d in cd.index])
        nxt = np.nonzero(np.diff(dn) == 1)[0]
        dd = np.abs(np.diff(cd.lnc.values)[nxt])
        pr = np.hypot(cd.sig_rep.values[nxt], cd.sig_rep.values[nxt + 1])
        ps = np.hypot(cd.sig_sem.values[nxt], cd.sig_sem.values[nxt + 1])
        from scipy.stats import spearmanr
        rs = spearmanr(dd, pr)[0] if dd.size > 10 else np.nan
        rows.append(dict(
            site=site, ident=ident, label=LABEL[(site, ident)], days=len(cd),
            n_prof_med=float(cd.n.median()), n_scenes_med=float(cd.n_scenes.median()),
            sig_rep_med=100 * float(cd.sig_rep.median()),
            sig_sem_med=100 * float(cd.sig_sem.median()),
            spread_robust_med=100 * float(cd.p_mad_rel.median()),
            sig_w_fit=100 * sig_w,
            sig_w_lag1=(100 * float(np.sqrt(sf.loc[sf.lag == 1, "D"].iloc[0] / 2.0))
                        if (sf.lag == 1).any() else np.nan),
            n_pairs_lag1=int(sf.loc[sf.lag == 1, "n"].iloc[0]) if (sf.lag == 1).any() else 0,
            q_per_day=q, drift_3d=100 * drift(3), drift_7d=100 * drift(7),
            pull_rep_consec=float(1.4826 * np.median(dd / pr)) if dd.size else np.nan,
            pull_sem_consec=float(1.4826 * np.median(dd / ps)) if dd.size else np.nan,
            spearman_rep_vs_dd=float(rs),
        ))
        models[(site, ident)] = dict(sig_w=sig_w, q=q, drift=drift)
        # my cloud rerun against the published cloud constants (calibration_version 1000)
        pb = pub[(pub.site == site) & (pub.ident == ident) & (pub.method == "cloud")]
        pb = pb.set_index("date")
        both = cd.index.intersection(pb.index)
        if len(both):
            r = np.exp(cd.loc[both, "lnc"].values) / pb.loc[both, "c"].values
            cmp_rows.append(dict(site=site, ident=ident, n_rerun=len(cd), n_pub=len(pb),
                                 n_both=len(both), ratio_med=float(np.median(r)),
                                 within_1pct=100 * float(np.mean(np.abs(r - 1) < 0.01)),
                                 within_1e6=100 * float(np.mean(np.abs(r - 1) < 1e-6)),
                                 u_pub_vs_rerun=float(np.median(
                                     (pb.loc[both, "u"].values / pb.loc[both, "c"].values)
                                     / cd.loc[both, "sig_rep"].values))))
    tab = pd.DataFrame(rows)
    sfa = pd.concat(sfs, ignore_index=True) if sfs else pd.DataFrame()
    cmp = pd.DataFrame(cmp_rows)
    RC.OUT.mkdir(parents=True, exist_ok=True)
    tab.to_csv(RC.OUT / "cloud_uncertainty.csv", index=False)
    sfa.to_csv(RC.OUT / "cloud_structure.csv", index=False)
    cmp.to_csv(RC.OUT / "cloud_vs_published.csv", index=False)
    if make_figure and len(tab):
        figure(clouds, tab, sfa, models)
    return tab, sfa, cmp, models


def figure(clouds, tab, sfa, models=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {("payerne", "C"): "#2ca02c", ("lindenberg", "C"): "#1f77b4",
            ("aosta", "B"): "#9467bd", ("payerne", "B"): "#ff7f0e"}
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.4))
    a = ax[0]
    for (site, ident), g in sfa.groupby(["site", "ident"]):
        a.plot(g.lag, 100 * np.sqrt(g.D / 2.0), "o", ms=3, color=cols[(site, ident)],
               label=LABEL[(site, ident)])
        if models and (site, ident) in models:
            m = models[(site, ident)]
            LL = np.arange(1, MAX_LAG + 1)
            a.plot(LL, 100 * np.sqrt((2 * m["sig_w"] ** 2 + m["q"] * LL) / 2.0), "-",
                   color=cols[(site, ident)], lw=1.2)
        t = tab[(tab.site == site) & (tab.ident == ident)].iloc[0]
        a.axhline(t.sig_rep_med, color=cols[(site, ident)], ls="--", lw=1)
        a.axhline(t.sig_sem_med, color=cols[(site, ident)], ls=":", lw=1)
    a.axvline(FIT_LAGS, color="0.5", lw=0.8, ls="-.")
    a.set_xlabel("lag between cloud-calibration days [days]")
    a.set_ylabel("sqrt(D(L) / 2) of ln C  [%]")
    a.set_title("day-to-day reproducibility (dots; line: noise + drift fit, L <= 7 d)\n"
                "vs reported spread (dashed) and standard error of the median (dotted)",
                fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
    a = ax[1]
    for (site, ident) in cols:
        cd = P.cloud_days(clouds, site, ident)
        if len(cd) < 10:
            continue
        a.hist(100 * cd.sig_rep, bins=np.arange(0, 30.5, 1), histtype="step", lw=1.5,
               color=cols[(site, ident)], label=LABEL[(site, ident)], density=True)
    a.set_xlabel("reported 1-sigma (std / median of the per-profile coefficients) [%]")
    a.set_ylabel("density")
    a.set_title("published cloud uncertainty, per day", fontsize=10)
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
    a = ax[2]
    for (site, ident) in cols:
        cd = P.cloud_days(clouds, site, ident)
        if len(cd) < 10:
            continue
        a.plot(cd.n_scenes, 100 * cd.sig_rep, ".", ms=3, alpha=0.5, color=cols[(site, ident)],
               label=LABEL[(site, ident)])
    a.set_xscale("log")
    a.set_xlabel("independent cloud scenes in the day")
    a.set_ylabel("reported 1-sigma [%]")
    a.set_title("reported spread vs number of scenes", fontsize=10)
    a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(RC.FIG_DIR / "ref23_cloud_uncertainty.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    t, s, c, m = run()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(t.to_string())
    print(c.to_string())
