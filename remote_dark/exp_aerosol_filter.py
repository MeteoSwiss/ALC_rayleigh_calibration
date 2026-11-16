# -*- coding: utf-8 -*-
"""Fig 14 — operator's question: the aerosol concentration CHANGES over time while the
electronic dark does not. Can that difference filter the aerosol out of the sky evidence?

Setup (Payerne A = CHM15k, C = CL61; hood truth available for both). Design lesson re-learned
while writing this: the per-gate intercept is identifiable ONLY because the atmospheric part of
the signal VARIES night to night (transmission +-40 %) while the dark does not — so the
regressor must be each night's PREDICTED atmospheric signal x_n(z) = C_n·M_n(z) (C_n = that
night's constant from the stable 4.5-6.5 km window), never a night-normalised signal (dividing
by C_n destroys the leverage and makes the dark a varying quantity; first attempt returned
theta = -100). Per gate, four intercepts are compared:

  d_all    S = c·(C_n M) + d            over ALL nights     (the v4-evidence design)
  d_clean  same, cleanest third of nights by CAMS aerosol load
  d_dirty  same, dirtiest third
  d_zero   S = c·(C_n M) + e·(C_n Aer) + d — each night's CAMS aerosol profile is its own
           regressor, so d is the ZERO-AEROSOL extrapolation (the operator's idea in
           regression form: what varies with the aerosol is aerosol, what stays is dark+rest)

If the intercept were pure dark, clean and dirty nights would give the SAME d. The spread
clean->dirty measures how much aerosol sits in the evidence; d_zero shows how much of it the
CAMS-shaped time variability can remove; the distance d_zero -> hood is the part that CANNOT
be filtered this way (mean load never sampled at zero + shapes CAMS does not model).

Judged over 2-8 km in rcs view: theta = <d,hood>/<hood,hood> (amplitude of hood contained in
d; target 1) and RMS(d - hood).
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import DATA, V4_DIR, PAYERNE, ensure_out

FIG = ensure_out("figs")
M1_DIR = DATA / "remote_dark" / "m1"
CWIN = (4500.0, 6500.0)            #: per-night constant window (the fig10 stable band)
ABAND = (2000.0, 6000.0)           #: aerosol-load index band
JBAND = (2000.0, 8000.0)           #: judgement band (rcs view)

UNITS = [("A", "CHM15k", "0-20000-0-06610_A_20250402_20260813_220.npz"),
         ("C", "CL61", "0-20000-0-06610_C_20250101_20260813_220.npz")]


def nightly_c(S, sig, M, rng):
    m = (rng >= CWIN[0]) & (rng <= CWIN[1])
    w = 1.0 / np.maximum(sig[:, m], 1e-300) ** 2
    ok = np.isfinite(S[:, m]) & np.isfinite(M[:, m])
    w = np.where(ok, w, 0.0)
    S0, M0 = np.nan_to_num(S[:, m]), np.nan_to_num(M[:, m])
    return np.nansum(w * S0 * M0, axis=1) / np.maximum(np.nansum(w * M0 * M0, axis=1), 1e-300)


def gate_fit(Sp, M, Aer=None):
    """Per-gate least squares over nights: Sp = c*M (+ e*Aer) + d. Returns d(z)."""
    n = np.isfinite(Sp) & np.isfinite(M) & (np.isfinite(Aer) if Aer is not None else True)
    cnt = n.sum(axis=0).astype(float)
    X = [np.where(n, M, 0.0), np.ones_like(M) * n]
    if Aer is not None:
        X.insert(1, np.where(n, Aer, 0.0))
    Y = np.where(n, Sp, 0.0)
    d = np.full(M.shape[1], np.nan)
    for j in range(M.shape[1]):
        if cnt[j] < 8:
            continue
        A = np.column_stack([x[:, j] for x in X])
        try:
            beta = np.linalg.lstsq(A, Y[:, j], rcond=None)[0]
            d[j] = beta[-1]
        except np.linalg.LinAlgError:
            pass
    return d


def smooth(x, k=7):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def judge(d_rcs, hood_b, rng):
    m = (rng >= JBAND[0]) & (rng <= JBAND[1]) & np.isfinite(d_rcs) & np.isfinite(hood_b)
    th = float(np.nansum(d_rcs[m] * hood_b[m]) / np.nansum(hood_b[m] ** 2))
    rms = float(np.sqrt(np.nanmean((d_rcs[m] - hood_b[m]) ** 2)))
    return th, rms


def main():
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.2))
    print(f"{'unit':8s} {'variant':>16s} {'theta(vs hood)':>14s} {'RMS vs hood':>12s}")
    for ax, (ident, itype, cache_name) in zip(axes, UNITS):
        z = np.load(V4_DIR / "cache" / cache_name, allow_pickle=True)
        S, sig, M, Aer, rng = z["S"], z["sigma"], z["M"], z["Aer"], z["rng"]
        m1 = np.load(M1_DIR / f"{PAYERNE['wmo']}_{ident}_m1.npz", allow_pickle=True)
        hood_b = np.asarray(m1["hood_b"], float)

        C_n = nightly_c(S, sig, M, rng)
        okn = np.isfinite(C_n) & (C_n > 0)
        Y = S[okn]
        Xm = C_n[okn, None] * M[okn]        # each night's PREDICTED molecular signal
        Xa = C_n[okn, None] * Aer[okn]      # ... and its predicted aerosol signal

        # per-night aerosol-load index from the CAMS profile itself (relative to molecular)
        ab = (rng >= ABAND[0]) & (rng <= ABAND[1])
        with np.errstate(invalid="ignore", divide="ignore"):
            a_idx = np.nanmean(Aer[okn][:, ab] / M[okn][:, ab], axis=1)
        t1, t2 = np.nanpercentile(a_idx, [33.3, 66.7])

        variants = [
            ("all nights", gate_fit(Y, Xm), "steelblue", "-"),
            (f"cleanest 1/3 (n={int((a_idx <= t1).sum())})",
             gate_fit(Y[a_idx <= t1], Xm[a_idx <= t1]), "seagreen", "-"),
            (f"dirtiest 1/3 (n={int((a_idx >= t2).sum())})",
             gate_fit(Y[a_idx >= t2], Xm[a_idx >= t2]), "darkorange", "-"),
            ("zero-aerosol extrapolation\n(CAMS profile as regressor)",
             gate_fit(Y, Xm, Xa), "crimson", "-"),
        ]

        z2 = np.where(rng > 0, rng, np.nan) ** 2
        mm = (rng >= 1000) & (rng <= 8000)
        ax.plot(smooth(hood_b / z2)[mm], rng[mm] / 1e3, "k-", lw=2.2, label="hood truth")
        for lab, d_rcs, colr, ls in variants:
            th, rms = judge(d_rcs, hood_b, rng)
            ax.plot(smooth(d_rcs / z2)[mm], rng[mm] / 1e3, ls, color=colr, lw=1.5,
                    label=f"{lab}  [θ={th:+.1f}]")
            print(f"{itype:8s} {lab.splitlines()[0]:>16.16s} {th:14.2f} {rms:12.3e}")
        ax.axvline(0, color="k", lw=0.5)
        ref = smooth(hood_b / z2)[mm]
        span = np.nanmax(np.abs(ref))
        ax.set_xlim(-6 * span, 6 * span)
        ax.set_title(f"{itype} (Payerne {ident}) — {int(okn.sum())} nights\n"
                     "does the aerosol's time variability filter it out of the evidence?",
                     fontsize=10)
        ax.set_xlabel("intercept d, P-view (rcs/z²)")
        ax.set_ylabel("altitude (km)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Aerosol filtering by time variability — per-gate intercepts from night "
                 "subsets and from zero-aerosol extrapolation, vs the hood truth "
                 "(θ = amplitude of hood contained in the curve; target 1)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig14_aerosol_filter.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
