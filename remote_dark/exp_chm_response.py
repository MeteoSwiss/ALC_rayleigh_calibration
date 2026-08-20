# -*- coding: utf-8 -*-
"""Fig 19 — CHM15k above-cloud response with cloud brightness Q as a CONTINUOUS input
(operator: "use the different cloud intensities to better characterize the detector physics").

Four questions, one per panel, on the 36k opaque-cloud night profiles:

  P1  how does the above-cloud profile deform as Q grows?      (median profile per Q sextile)
  P2  is the response LINEAR in Q at fixed height?             (y vs Q at three z' slices,
      hood value plotted at Q -> 0 as the quiescent anchor)
  P3  WHERE does the response live?                            (per-gate kernel k(z') from a
      line through the six sextile medians - robust to tails)
  P4  the rescue test: does the per-gate extrapolation Q -> 0  (intercept d0(z') vs hood)
      recover the quiescent dark the terciles could not?

Interpretation guide: if P2 is linear over the observed Q range, the Q->0 intercept is
trustworthy and the CHM15k natural hood is RESCUED by regression; if it curves (saturation),
the intercept inherits the curvature and P4's disagreement with the hood measures that bias.
The kernel P3 doubles as the prediction for the HOOD's own internal-reflection response
(a "cloud" at 0 m): fig20 compares it against the hood session's drift shape.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import DATA, ensure_out
from remote_dark.exp_cloud_impact import collect, ZREL
from remote_dark import hood

FIG = ensure_out("figs")
CACHE = DATA / "remote_dark" / "chm_cloud_stack.npz"
NBIN = 6
SLICES = (1500.0, 2500.0, 4000.0)
JBAND = (500.0, 3000.0)


def load_stack():
    if CACHE.exists():
        z = np.load(CACHE)
        return z["Y"], z["Q"], z["CB"], int(z["n_days"])
    Y, Q, CB, n_days = collect("A")
    np.savez_compressed(CACHE, Y=Y, Q=Q, CB=CB, n_days=n_days)
    return Y, Q, CB, n_days


def smooth(x, k=5):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    Y, Q, CB, n_days = load_stack()
    edges = np.nanpercentile(Q, np.linspace(0, 100, NBIN + 1))
    qmid, prof = [], []
    for i in range(NBIN):
        m = (Q >= edges[i]) & (Q <= edges[i + 1])
        qmid.append(float(np.nanmedian(Q[m])))
        prof.append(np.nanmedian(Y[m], axis=0))
    qmid = np.array(qmid)
    prof = np.vstack(prof)                       # (NBIN, ZREL)

    trng, b_rcs, _sem, _bp = hood.truth("A")
    med_cb = float(np.nanmedian(CB))
    zt = np.where(trng > 0, trng, np.nan)
    hood_rel = np.interp(med_cb + ZREL, trng, b_rcs / zt ** 2, left=np.nan, right=np.nan)

    # per-gate line through the six sextile medians: kernel k and Q->0 intercept d0
    k = np.full(ZREL.size, np.nan)
    d0 = np.full(ZREL.size, np.nan)
    for j in range(ZREL.size):
        y = prof[:, j]
        ok = np.isfinite(y)
        if ok.sum() >= 4:
            p = np.polyfit(qmid[ok], y[ok], 1)
            k[j], d0[j] = p[0], p[1]

    fig, (a1, a2, a3, a4) = plt.subplots(1, 4, figsize=(16.5, 6.2))
    sc = 1e3
    cmap = plt.cm.viridis(np.linspace(0.1, 0.95, NBIN))

    for i in range(NBIN):
        a1.plot(smooth(prof[i]) * sc, ZREL / 1e3, "-", color=cmap[i], lw=1.4,
                label=f"Q sextile {i + 1}")
    a1.plot(smooth(hood_rel) * sc, ZREL / 1e3, "k--", lw=1.8, label="hood dark")
    a1.axvline(0, color="k", lw=0.5)
    vm = ZREL >= 1000
    span = float(np.nanpercentile(np.abs(prof[:, vm]) * sc, 97))
    a1.set_xlim(-1.4 * span, 0.6 * span)
    a1.set_ylabel("height above cloud base (km)")
    a1.set_xlabel("offset, P-view (×1e-3)")
    a1.set_title("P1 — profile per Q sextile\n(the bowl deepens with charge)", fontsize=10)
    a1.legend(fontsize=7, loc="lower left")
    a1.grid(alpha=0.25)

    for zc, colr in zip(SLICES, ("steelblue", "crimson", "darkorange")):
        j = int(np.argmin(np.abs(ZREL - zc)))
        ys = np.array([np.nanmedian(smooth(prof[i])[max(0, j - 2):j + 3])
                       for i in range(NBIN)])
        a2.plot(qmid, ys * sc, "o-", color=colr, lw=1.5, ms=5,
                label=f"z' = {zc / 1e3:.1f} km")
        jh = smooth(hood_rel)[max(0, j - 2):j + 3]
        a2.plot(0.0, float(np.nanmedian(jh)) * sc, "s", color=colr, ms=8, mfc="none",
                mew=1.8)
    a2.axhline(0, color="k", lw=0.5)
    a2.set_xlabel("cloud-return charge Q (recorded)")
    a2.set_ylabel("offset, P-view (×1e-3)")
    a2.set_title("P2 — response vs Q at fixed height\n(open squares: hood value at Q=0)",
                 fontsize=10)
    a2.legend(fontsize=8)
    a2.grid(alpha=0.25)

    a3.plot(smooth(k * np.nanmedian(np.diff(qmid)) * NBIN) * sc, ZREL / 1e3, "-",
            color="purple", lw=1.8)
    a3.axvline(0, color="k", lw=0.5)
    a3.set_xlabel("kernel k(z') × Q-range, P-view (×1e-3)")
    a3.set_title("P3 — response kernel: where the\ncloud's fingerprint lives", fontsize=10)
    a3.grid(alpha=0.25)

    a4.plot(smooth(d0) * sc, ZREL / 1e3, "-", color="crimson", lw=1.9,
            label="Q→0 extrapolation d0(z')")
    a4.plot(smooth(hood_rel) * sc, ZREL / 1e3, "k--", lw=1.8, label="hood dark")
    a4.axvline(0, color="k", lw=0.5)
    jm = (ZREL >= JBAND[0]) & (ZREL <= JBAND[1]) & np.isfinite(d0) & np.isfinite(hood_rel)
    th = float(np.nansum(d0[jm] * hood_rel[jm]) / np.nansum(hood_rel[jm] ** 2))
    rms = float(np.sqrt(np.nanmean((d0[jm] - hood_rel[jm]) ** 2))
                / np.sqrt(np.nanmean(hood_rel[jm] ** 2)))
    span4 = float(np.nanpercentile(np.abs(np.concatenate(
        [d0[vm], hood_rel[vm]])) * sc, 97))
    a4.set_xlim(-3 * span4, 3 * span4)
    a4.set_xlabel("offset, P-view (×1e-3)")
    a4.set_title(f"P4 — the rescue test\nθ vs hood = {th:+.2f}, RMS/|hood| = {rms:.2f}",
                 fontsize=10)
    a4.legend(fontsize=8)
    a4.grid(alpha=0.25)

    # linearity metric: curvature of P2 at the middle slice (quad term vs linear term)
    j = int(np.argmin(np.abs(ZREL - SLICES[1])))
    ys = prof[:, j]
    c2 = np.polyfit(qmid, ys, 2)
    curv = float(abs(c2[0]) * (qmid[-1] - qmid[0]) / max(abs(c2[1]), 1e-300))
    print(f"{Y.shape[0]} profiles, {n_days} days; Q sextile medians "
          + " ".join(f"{q:.2e}" for q in qmid))
    print(f"linearity at 2.5 km: |quad/lin| over the Q range = {curv:.2f} "
          f"(0 = perfectly linear)")
    print(f"rescue test: theta = {th:+.2f}, RMS/|hood| = {rms:.2f}")
    fig.suptitle("CHM15k above-cloud response with cloud brightness as a continuous input — "
                 "profile deformation, linearity, kernel, and the Q→0 rescue test",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    f = FIG / "fig19_chm_response.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
