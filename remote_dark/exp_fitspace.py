# -*- coding: utf-8 -*-
"""Experiment (operator's proposal): is the Rayleigh amplitude more STABLE when fitted in
S-view (signal, log scale) than in RCS view?

Three estimators of the same nightly constant, same nights, same window:

  E1  WLS in RCS view (the current one):     C1 = sum(w S M) / sum(w M^2),  w = 1/sigma^2
  E2  plain fit in S-view (P = rcs/z^2):     C2 = sum(P Mp) / sum(Mp^2)    (no sigma weights --
      the view IS the weighting: 1/z^4 relative to RCS)
  E3  median log-ratio in S-view:            C3 = exp( median_z log(P/Mp) ) over gates P,Mp > 0
      -- "fit the molecular as an offset in log-S", the estimator the operator's eye performs
  E4  median ratio, linear:                  C4 = median_z (S/M)  with NO positivity requirement
      -- same robust statistic as E3 but negative-noise gates are kept, so it cannot carry the
      log-estimator's positivity selection bias at low SNR (the control that separates
      "robust median helps" from "dropping negative gates fakes stability")

The metric is the one that matters operationally: the night-to-night robust CV of the constant
(MAD/median over clear nights), per window altitude. Real T2_aerosol variability is common to all
estimators, so DIFFERENCES between them isolate their noise/misfit handling. The companion bias
table (median C_k/C_WLS per window) is the honesty gauge: stability bought by bias is worthless.
"""
from __future__ import annotations

import numpy as np

from remote_dark.common import V4_DIR, ensure_out

WIN_M = 1000.0
CENTRES = np.arange(2500.0, 9501.0, 500.0)

STREAMS = [("0-20000-0-06610_A_20250402_20260813_220.npz", "Payerne A"),
           ("0-20000-0-06610_C_20250101_20260813_220.npz", "Payerne C"),
           ("0-20000-0-06240_A_20250101_20260813_220.npz", "Schiphol A"),
           ("0-20000-0-06240_B_20250101_20260813_220.npz", "Schiphol B"),
           ("0-20000-0-06240_C_20250101_20260813_220.npz", "Schiphol C"),
           ("0-20000-0-06240_D_20250101_20260813_220.npz", "Schiphol D")]


def estimators(S, sig, M, rng, zc):
    m = (rng >= zc - WIN_M / 2) & (rng <= zc + WIN_M / 2)
    if m.sum() < 10:
        return np.full((S.shape[0], 3), np.nan)
    z2 = rng[m] ** 2
    Sw, Mw, gw = S[:, m], M[:, m], 1.0 / np.maximum(sig[:, m], 1e-300) ** 2
    ok = np.isfinite(Sw) & np.isfinite(Mw)
    gw = np.where(ok, gw, 0.0)
    S0, M0 = np.nan_to_num(Sw), np.nan_to_num(Mw)

    c1 = np.nansum(gw * S0 * M0, axis=1) / np.maximum(np.nansum(gw * M0 * M0, axis=1), 1e-300)

    P, Mp = S0 / z2[None, :], M0 / z2[None, :]
    okp = ok
    c2 = (np.nansum(np.where(okp, P * Mp, 0.0), axis=1)
          / np.maximum(np.nansum(np.where(okp, Mp * Mp, 0.0), axis=1), 1e-300))

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(ok & (Sw > 0) & (Mw > 0), Sw / Mw, np.nan)
        c3 = np.exp(np.nanmedian(np.log(ratio), axis=1))
        # E4: same median, linear space, negatives KEPT -> no positivity selection possible
        ratio_lin = np.where(ok & (Mw > 0), Sw / Mw, np.nan)
        c4 = np.nanmedian(ratio_lin, axis=1)
    return np.column_stack([c1, c2, c3, c4])


def robust_cv(x):
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < 20:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med)) / med


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG = ensure_out("figs")

    fig, axes = plt.subplots(2, len(STREAMS), figsize=(15.5, 8.6), sharey=True)
    print(f"{'stream':12s} {'window':>7s}  {'CV WLS-RCS':>10s} {'CV S-view':>10s} "
          f"{'CV S-log-med':>12s} {'CV S-lin-med':>12s}")
    for icol, (name, lab) in enumerate(STREAMS):
        ax, axb = axes[0, icol], axes[1, icol]
        z = np.load(V4_DIR / "cache" / name, allow_pickle=True)
        S, sig, M, rng = z["S"], z["sigma"], z["M"], z["rng"]
        cvs = np.full((CENTRES.size, 4), np.nan)
        bias = np.full((CENTRES.size, 4), np.nan)
        for j, zc in enumerate(CENTRES):
            C = estimators(S, sig, M, rng, zc)
            cvs[j] = [robust_cv(C[:, k]) for k in range(4)]
            ok = np.isfinite(C[:, 0]) & (C[:, 0] > 0)
            if ok.sum() >= 20:
                with np.errstate(invalid="ignore"):
                    bias[j] = [np.nanmedian(C[ok, k] / C[ok, 0]) for k in range(4)]
        styles = (("WLS in RCS (current)", "k", "-"),
                  ("plain fit in S", "steelblue", "--"),
                  ("median log-ratio in S", "crimson", "-"),
                  ("median ratio, linear", "seagreen", "-"))
        for k, (lbl, colr, ls) in enumerate(styles):
            ax.plot(cvs[:, k] * 100, CENTRES / 1e3, ls, color=colr, lw=1.6, label=lbl)
            if k >= 2:
                axb.plot((bias[:, k] - 1) * 100, CENTRES / 1e3, ls, color=colr, lw=1.6,
                         label=lbl)
        axb.axvline(0, color="k", lw=0.6)
        ax.set_xlim(0, 80)
        axb.set_xlim(-30, 210)
        ax.set_title(f"{lab}\n{S.shape[0]} nights", fontsize=9)
        ax.set_xlabel("night-to-night CV of Ĉ (%)")
        axb.set_xlabel("median bias vs WLS (%)")
        ax.grid(alpha=0.25)
        axb.grid(alpha=0.25)
        if icol == 0:
            ax.set_ylabel("window-centre altitude (km)")
            axb.set_ylabel("window-centre altitude (km)")
            ax.legend(fontsize=7, loc="upper right")
            axb.legend(fontsize=7, loc="upper right")
        jref = int(np.argmin(np.abs(CENTRES - 4500)))
        print(f"{lab:12s} {'4.5 km':>7s}  {cvs[jref,0]*100:10.1f} {cvs[jref,1]*100:10.1f} "
              f"{cvs[jref,2]*100:12.1f} {cvs[jref,3]*100:12.1f}")
    fig.suptitle("Fit-space experiment — same nights, same windows: which estimator of the "
                 "Rayleigh constant is most stable night-to-night?\n"
                 "Top: stability (lower = better).  Bottom: the honesty gauge — median bias vs "
                 "WLS (stability bought by bias is worthless).", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig11_fitspace.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
