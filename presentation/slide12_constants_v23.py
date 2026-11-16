"""Slide 12 of the Annual Meeting 2026 deck, v2.3 edition ("Per-instrument constants matter"): paper Figs. 4
and 7 of manuscript v1.0 as ONE two-row figure with slide-sized fonts, the layout of presentation/
slide12_constants.py. Data: the calibration v2.3 release archive (validation/paper/fig_v10_chm15k_network.collect:
median and interquartile range of the successful nightly Rayleigh constants of every CHM15k, flags 1 / 0.5 /
0.25) and the Kleine Scheidegg nightly constants smoothed with the improve_alc Kalman filter weighted by the v2.3
uncertainty (the rule of validation/paper/fig_v10_kalman_kleine_scheidegg.py).

Run from the repo root:  python presentation/slide12_constants_v23.py
Output: presentation/figs_annual2026/fig_constants_v23_slide.png (+ STATS line)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation.paper import fig_v10_chm15k_network as NET  # noqa: E402
from validation.paper.build_calib_v23 import nightly  # noqa: E402
from validation.paper.fig_v10_kalman_kleine_scheidegg import CAL as KS_CAL  # noqa: E402
from validation.paper.kalman_improve_alc import fitfunc_sin, kalman_improve_alc  # noqa: E402

OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_constants_v23_slide.png"
LABEL = {"0-20000-0-06736_A": "Kleine Scheidegg", "0-20000-0-06610_A": "Payerne",
         "0-20000-0-10962_0": "Hohenpeissenberg"}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 17, "axes.labelsize": 17, "xtick.labelsize": 15.5,
                         "ytick.labelsize": 15.5, "legend.fontsize": 15})
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(17.28, 6.3), gridspec_kw=dict(hspace=0.36))
    fig.subplots_adjust(left=0.06, right=0.995, top=0.985, bottom=0.085)

    # ---- (a) sorted network constants, v2.3 release
    df = NET.collect()
    x = np.arange(len(df))
    a1.bar(x, df.med, width=0.85, color="#c0392b")
    a1.errorbar(x, df.med, yerr=[df.med - df.q25, df.q75 - df.med], fmt="none", ecolor="0.3",
                elinewidth=0.7, capsize=0)
    net_med = float(np.median(df.med))
    dflt = NET.DEFAULT_CL / NET.UNIT
    a1.axhline(dflt, color="k", ls="--", lw=1.8, label="default constant 3 × 10$^{11}$")
    a1.set_ylim(0, 10)
    a1.set_xlim(-1, len(df))
    for i, r in df.iterrows():
        if r.key in LABEL:
            y = min(r.q75, 9) + 0.3
            # a unit near the top of the axis gets its label to the left, inside the axes
            if y + 2.4 > 9.6:
                a1.annotate(LABEL[r.key], (i, min(r.med, 9.5)), xytext=(i - 14, 8.6), ha="right", va="center",
                            fontsize=15, arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
            else:
                a1.annotate(LABEL[r.key], (i, y), xytext=(i, y + 2.4), ha="center", fontsize=15,
                            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
    a1.set_ylabel("C$_L$ (10$^{11}$)")
    a1.set_xlabel(f"{len(df)} CHM15k, sorted by their median Rayleigh constant (v2.3)", labelpad=2)
    a1.set_xticks([])
    a1.grid(axis="y", alpha=0.3)
    a1.legend(loc="upper left", frameon=False)
    n_off = int((df.q75 > 10).sum())

    # ---- (b) Kalman best estimate, Kleine Scheidegg, weighted by the v2.3 uncertainty
    d, v, u = nightly(KS_CAL, "rayleigh")
    t = np.array(d)
    r = kalman_improve_alc(list(d), v, "CHM15k", uncertainties=u)
    od = set(r["outlier_dates"])
    out = np.array([pd.Timestamp(z).date() in od for z in t])
    g = pd.to_datetime(r["grid"])
    st, sd = r["state"], r["std"]
    a2.fill_between(g, (st - sd) / 1e11, (st + sd) / 1e11, color="tab:green", alpha=0.25, lw=0)
    a2.plot(g, st / 1e11, color="tab:green", lw=2.2, label="Kalman best estimate ±1σ")
    if r["model"] == "annual cycle":
        doy = np.array([z.timetuple().tm_yday for z in g])
        a2.plot(g, fitfunc_sin(doy, *r["params"]) / 1e11, ":", color="k", lw=1.8, label="annual-cycle model")
        amp, _ph, c0 = r["params"]
    else:
        amp, c0 = np.nan, r["params"][0]
    a2.plot(t[~out], v[~out] / 1e11, "o", color="tab:red", ms=4.5, label=f"nightly C$_L$ ({int((~out).sum())})")
    a2.plot(t[out], v[out] / 1e11, "x", color="0.4", ms=8, mew=1.6, label=f"outlier ({int(out.sum())})")
    a2.set_ylim(0.2, 3.5)
    a2.set_xlim(pd.Timestamp("2024-12-15"), pd.Timestamp("2026-11-15"))
    a2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11]))
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
    a2.set_ylabel("C$_L$ (10$^{11}$)")
    a2.set_yticks([1, 2, 3])
    a2.grid(alpha=0.3)
    a2.legend(loc="lower center", ncol=4, frameon=True, framealpha=0.95, borderpad=0.3, columnspacing=1.6)
    a2.text(0.005, 0.95, "Kleine Scheidegg CHM15k, v2.3", transform=a2.transAxes, va="top", fontsize=16,
            fontweight="bold")
    fig.align_ylabels([a1, a2])
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    rel_u = u / v / 2.0
    stats = dict(n_units=len(df), net_median=round(net_med, 2), n_clipped=n_off,
                 off20=round(float((abs(df.med / dflt - 1) > 0.2).mean()), 3),
                 range=[round(float(df.med.min()), 2), round(float(df.med.max()), 2)],
                 n_lowsun_units=int((df.n_lowsun > 0).sum()),
                 ks_nights=int(v.size), ks_outliers=int(out.sum()), ks_model=r["model"],
                 ks_p2p_pct=round(float(2 * amp / c0 * 100), 1) if amp == amp else None,
                 ks_sigma_m_pct=round(float(np.sqrt(r["var_meas"]) / np.median(v) * 100), 1),
                 ks_median_1sigma_pct=round(float(100 * np.median(rel_u)), 1),
                 ks_last=str(pd.Timestamp(t.max()).date()))
    print("wrote", OUT_PNG)
    print("STATS", json.dumps(stats))


if __name__ == "__main__":
    main()
