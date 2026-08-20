# -*- coding: utf-8 -*-
"""Fig 8 (v2: successive-difference noise) — the simple question, asked simply: is the measured clear-night noise just
molecular shot noise + a background floor, or is there instrument noise on top?

One panel per Schiphol unit. Everything in P-view (raw signal, rcs/z²) where the physics is
plain: shot noise follows sqrt(signal), the background/detector floor is FLAT, and any excess —
especially the near-range rise where electronic noise is amplified by 1/overlap — is instrument.

    measured  NOISE(z) = median over clear nights of the per-profile noise  sigma_n·sqrt(n_prof)/z²
    model     sqrt( α·P_signal(z) + β )   α fitted at 4–8 km (shot regime), β from the top gates

Altitude on Y, noise on X (log). The gap between black and red IS the instrument noise.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import V4_DIR, ensure_out

WMO = "0-20000-0-06240"
IDENTS = "ABCD"
FIG = ensure_out("figs")


def load(ident):
    z = np.load(V4_DIR / "cache" / f"{WMO}_{ident}_20250101_20260813_220.npz", allow_pickle=True)
    return {k: z[k] for k in z.files}


def direct_noise(ident, dates, n_sample=25):
    """Per-profile noise from RAW L1 via successive differences — immune to real atmospheric
    variability, which the cache's spread-around-the-median is not (operator's catch: at
    500-1000 m the cache metric runs ~5x the true fast noise because the nocturnal boundary
    layer MOVES; successive differencing kills everything slow)."""
    from datetime import datetime
    from remote_dark.common import read_day, solar_elevation_deg
    per = []
    rng = None
    step = max(1, len(dates) // n_sample)
    for ds in list(dates)[::step]:
        d = read_day(WMO, ident, datetime.strptime(str(int(ds)), "%Y%m%d"))
        if d is None:
            continue
        el = solar_elevation_deg(52.317, 4.804, d["times"])
        m = el < -6.0
        if m.sum() < 60:
            continue
        rng = d["rng"]
        dif = np.diff(d["rcs"][m], axis=0)
        per.append(1.4826 * np.nanmedian(np.abs(dif - np.nanmedian(dif, axis=0)), axis=0)
                   / np.sqrt(2))
    return rng, (np.nanmedian(np.vstack(per), axis=0) if per else None), len(per)


def main():
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.2), sharey=True)
    for ax, ident in zip(axes, IDENTS):
        c = load(ident)
        rng0, nd, n_used = direct_noise(ident, c["dates"])
        rng = c["rng"] if rng0 is None else rng0
        z2 = np.where(rng > 0, rng, np.nan) ** 2
        noise = nd / z2
        psig = np.nanmedian(c["S"], axis=0) / z2

        # the two-parameter physical model: alpha from the shot regime, beta from the top gates
        shot_band = (rng >= 4000) & (rng <= 8000) & np.isfinite(noise) & (psig > 0)
        top = rng >= np.nanpercentile(rng, 92)
        beta = float(np.nanmedian(noise[top] ** 2))
        alpha = float(np.nanmedian((noise[shot_band] ** 2 - beta) / psig[shot_band]))
        alpha = max(alpha, 0.0)
        model = np.sqrt(np.maximum(alpha * psig + beta, 0.0))

        m = (rng >= 100) & (rng <= 15000) & np.isfinite(noise)
        ax.plot(noise[m], rng[m] / 1e3, "k-", lw=1.6,
                label=f"measured fast noise ({n_used} nights, succ. diff.)")
        ax.plot(model[m], rng[m] / 1e3, "r--", lw=1.6,
                label="model: signal shot + background floor")
        ax.plot(np.sqrt(np.maximum(alpha * psig, 0))[m], rng[m] / 1e3, "-",
                color="steelblue", lw=1.0, alpha=0.8, label="signal shot alone")
        ax.axvline(np.sqrt(beta), color="0.55", lw=1.0, ls=":",
                   label="background/detector floor")
        ax.set_xscale("log")
        ax.set_xlabel("noise per profile, P-view (rcs/z²)")
        if ident == "A":
            ax.set_ylabel("altitude (km)")
        # the number the panel exists to show: the near-range excess factor
        near = (rng >= 100) & (rng <= 500)
        exc = float(np.nanmedian(noise[near] / model[near]))
        ax.set_title(f"unit {ident}\nnear-range excess ×{exc:.1f} (100–500 m)", fontsize=10)
        ax.set_ylim(0, 15)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("Schiphol quadruple — is the clear-night noise explained by signal shot noise + "
                 "a flat background floor?  The gap is instrument noise.", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig8_noise_budget.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
