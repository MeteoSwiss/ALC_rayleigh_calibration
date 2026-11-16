# -*- coding: utf-8 -*-
"""q19 - per-unit characteristic gradient of R(z) from the unit's own clean tier-1 PROFILES.

The availability study (doc/reports/rayleigh_availability.md 4.1) showed that a unit's C_L-vs-height
gradient cannot be estimated from v2's accepted constants: their windows span too little height. Here the
gradient is taken from the night-mean PROFILE instead, which spans 2.5-8 km on every night regardless of
where the window sits: d ln R / dz over [max(window bottom, 2.5 km), min(z_hi, 8 km)] (WLS, measured noise)
on winter (Nov-Mar), fully clear tier-1 nights whose window passes the vertical aerosol test. Compared
with the q16a fixed-window recalibrations (6-8 km vs 3-4 km windows, all nights).
Writes q19_unit_gradient.csv and q19_unit_gradient.png.
"""
from __future__ import annotations

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q = L.Q
warnings.filterwarnings("ignore")


def band_gradient(P_arr, z0, z1):
    z, sig, p, sg, _ = P_arr
    m = (z >= z0) & (z <= z1) & np.isfinite(sig) & np.isfinite(p) & (p > 0) & np.isfinite(sg) & (sg > 0)
    if m.sum() < 30:
        return np.nan, np.nan
    zc = 0.5 * (z0 + z1) / 1e3
    g, s, chi, _, _ = Q._grad(sig[m], p[m], z[m] / 1e3, 1.0 / sg[m] ** 2, zc)
    return g, s


def main():
    W = pickle.load(open(Q.Q19 / "q19_blocks_windows.pkl", "rb"))
    Wd = {(w["key"], w["date"]): w for w in W if w["variant"] == "op" and w["which"] == "sel"}
    G = Q.load_gate_variants()[0][1]["gate"]
    df = L.load(("CHM15k",))
    op = df[(df.variant == "op") & (df.tier == 1) & df.ok & (df.flag == 1.0)
            & df.month.isin([11, 12, 1, 2, 3])]
    rows = []
    for r in op.itertuples():
        b = Wd.get((r.key, r.date))
        if b is None or L.vc_verdict(b, k=G["k"], eps_tol=G["eps_tol"], tol_lo=np.inf) != "accept":
            continue
        P = L.profile_of(r.key, r.date)
        if P is None:
            continue
        z0 = max(r.sel_zs, 2500.0)
        z1 = min(r.sel_z_hi, 8000.0)
        if z1 - z0 < 2500:
            continue
        g, s = band_gradient(P, z0, z1)
        rows.append(dict(key=r.key, date=r.date, g=g, sg=s, z0=z0, z1=z1))
    N = pd.DataFrame(rows)
    U = N.groupby("key").agg(g_med=("g", "median"), g_q25=("g", lambda x: np.percentile(x, 25)),
                            g_q75=("g", lambda x: np.percentile(x, 75)), n=("g", "size"),
                            sg_med=("sg", "median")).reset_index()
    U["g_se"] = 1.2533 * N.groupby("key").g.std().values / np.sqrt(U.n.values)
    q16 = pd.read_csv(Q.Q19.parent / "q16" / "q16a_window_ratios.csv")
    pv = q16.pivot_table(index="key", columns="pair", values="median_pct")
    U["g_q16a"] = U.key.map(pv["z6_8/z3_4"] / 3.5)
    U.to_csv(Q.Q19 / "q19_unit_gradient.csv", index=False)
    print(U.round(2).to_string())
    m = np.isfinite(U.g_q16a)
    print("r(profile, q16a fixed windows) =", np.corrcoef(U.g_med[m], U.g_q16a[m])[0, 1].round(2),
          "n =", int(m.sum()))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw=dict(width_ratios=[1, 1.6]))
    a = ax[0]
    a.errorbar(U.g_q16a[m], U.g_med[m], yerr=U.g_se[m], fmt="o", color="k", ms=5, capsize=2)
    for _, r in U[m].iterrows():
        a.annotate(r.key.replace("0-20000-0-", "").replace("0-20008-0-", ""), (r.g_q16a, r.g_med),
                   fontsize=6, xytext=(3, 2), textcoords="offset points")
    lim = [-6, 14]
    a.plot(lim, lim, color="0.6", lw=0.8)
    a.set_xlabel("q16a: ln(C 6-8 km / C 3-4 km) / 3.5 km  [%/km], all nights")
    a.set_ylabel("q19: median d ln R/dz of clean winter tier-1 profiles [%/km]")
    a.set_title("unit characteristic gradient: profile vs fixed-window recalibration")
    a.grid(alpha=0.3)
    a = ax[1]
    o = U.sort_values("g_med")
    a.errorbar(np.arange(len(o)), o.g_med, yerr=[o.g_med - o.g_q25, o.g_q75 - o.g_med], fmt="o", color="k",
               capsize=2, label="median [IQR] of nights")
    a.axhline(0, color="0.5")
    a.set_xticks(np.arange(len(o)))
    a.set_xticklabels([k.replace("0-20000-0-", "").replace("0-20008-0-", "").replace("0-724-60000-0000000", "") for k in o.key], rotation=70,
                      fontsize=7)
    a.set_ylabel("d ln R/dz, 2.5-8 km [%/km]")
    a.set_title("per-unit profile gradient (clean winter tier-1 nights)")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q.Q19 / "q19_unit_gradient.png", dpi=120)


if __name__ == "__main__":
    main()
