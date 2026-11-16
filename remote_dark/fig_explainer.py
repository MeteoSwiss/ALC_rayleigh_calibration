# -*- coding: utf-8 -*-
"""Fig 13 — the campaign explained on one page, per instrument type, in ABSOLUTE units.

Two different "electronic" quantities, one row each, judged against the Payerne hood truth:

  row 1  the DARK  — the fixed structured offset the electronics ADD to every profile
         (a bias; it shifts the Rayleigh fit). Hard to retrieve from sky.
  row 2  the NOISE — the random fluctuation of each profile around its mean (a scatter;
         it blurs but does not shift). Easy to retrieve from sky, and validated here.

Row 1 shows, per type: the hood truth (black), the assumption-free sky evidence (blue,
the v4 per-gate regression intercept), and the M1 retrieval (red) with its honest label —
for the CL61 the red uses the HOOD'S OWN SHAPE and only the amplitude comes from the sky.

Row 2 is a genuine like-for-like closure: the same successive-difference noise estimator on
(a) the covered-telescope hood frames and (b) ordinary clear nights. Wherever the sky signal
is negligible the two must coincide if the sky retrieval is right — no shape is imposed.
CL31 era note: the sky nights predate the 2026-07-07 optic swap, so its hood-noise reference
uses the PRE-swap sessions (the pooled hood truth of row 1 is post-swap by construction).
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import DATA, V4_DIR, PAYERNE, ensure_out, read_day, solar_elevation_deg
from remote_dark import hood

FIG = ensure_out("figs")
M1_DIR = DATA / "remote_dark" / "m1"

UNITS = [
    ("B", "CL31", "0-20000-0-06610_B_borrowA_20250101_20260707_220.npz",
     "fit output NOT trustable: injection\nself-test gain −0.03 (blind at 910 nm)",
     "pre_swap", (300.0, 7700.0)),
    ("A", "CHM15k", "0-20000-0-06610_A_20250402_20260813_220.npz",
     "family fit: near shape OK but\namplitude ×3.9 too big (aerosol leaks in)",
     "single", (1000.0, 8000.0)),
    ("C", "CL61", "0-20000-0-06610_C_20250101_20260813_220.npz",
     "red = HOOD's shape, only the\namplitude is from sky: ×1.13 (target 1)",
     "single", (1000.0, 8000.0)),
]
N_SKY_NIGHTS = 14           #: nights sampled for the sky-noise estimate (median across them)


def smooth(x, k=7):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def sky_noise(ident, dates):
    """Per-profile fast noise from raw L1 nights (successive differences, night mask)."""
    per, rng = [], None
    step = max(1, len(dates) // N_SKY_NIGHTS)
    for ds in list(dates)[::step]:
        d = read_day(PAYERNE["wmo"], ident, datetime.strptime(str(int(ds)), "%Y%m%d"))
        if d is None:
            continue
        el = solar_elevation_deg(PAYERNE["lat"], PAYERNE["lon"], d["times"])
        m = el < -6.0
        if m.sum() < 60:
            continue
        rng = d["rng"]
        dif = np.diff(d["rcs"][m], axis=0)
        per.append(1.4826 * np.nanmedian(np.abs(dif - np.nanmedian(dif, axis=0)), axis=0)
                   / np.sqrt(2))
    return rng, (np.nanmedian(np.vstack(per), axis=0) if per else None), len(per)


def hood_noise(ident, era):
    """The SAME estimator on the covered-telescope frames — the truth for row 2."""
    per, rng = [], None
    for s in hood.session_frames(ident, min_hours=2.0):
        if era != "single" and s["era"] != era:
            continue
        dif = np.diff(s["dark"], axis=0)
        rng = s["rng"]
        per.append(1.4826 * np.nanmedian(np.abs(dif - np.nanmedian(dif, axis=0)), axis=0)
                   / np.sqrt(2))
    return rng, (np.nanmedian(np.vstack(per), axis=0) if per else None), len(per)


def main():
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.2))
    for col, (ident, itype, cache_name, verdict, era, band) in enumerate(UNITS):
        m1 = np.load(M1_DIR / f"{PAYERNE['wmo']}_{ident}_m1.npz", allow_pickle=True)
        rng = m1["rng"]
        z2 = np.where(rng > 0, rng, np.nan) ** 2

        # ---- row 1: the DARK, absolute P-view ------------------------------------------------
        ax = axes[0, col]
        mm = (rng >= band[0]) & (rng <= band[1])
        hood_p = smooth(m1["hood_b"] / z2)
        hat_p = smooth(m1["b_hat"] / z2)
        free_p = smooth(m1["b_free"] / z2)
        ax.plot(hood_p[mm], rng[mm] / 1e3, "k-", lw=2.0,
                label="hood truth (covered telescope)")
        ax.plot(free_p[mm], rng[mm] / 1e3, "-", color="steelblue", lw=1.0,
                alpha=0.9, label="sky, no assumptions (v4 evidence)")
        ax.plot(hat_p[mm], rng[mm] / 1e3, "r-", lw=1.8,
                label="sky, M1 retrieval")
        ax.axvline(0, color="k", lw=0.5)
        # x-limits from truth + retrieval in the display band; the near-range free-evidence
        # tail (overlap territory) is allowed to clip
        ref = np.concatenate([hood_p[mm], hat_p[mm]])
        lo, hi = np.nanpercentile(ref, [0.5, 99.5])
        pad = 0.35 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_title(f"{itype} (Payerne {ident})\n{verdict}", fontsize=10)
        ax.set_xlabel("dark offset, P-view (rcs/z²)")
        if col == 0:
            ax.set_ylabel("THE DARK (fixed bias)\naltitude (km)")
            ax.legend(fontsize=8, loc="best")
        ax.set_ylim(0, 8)
        ax.grid(alpha=0.25)

        # ---- row 2: the NOISE, same estimator sky vs hood ------------------------------------
        ax = axes[1, col]
        cache = np.load(V4_DIR / "cache" / cache_name, allow_pickle=True)
        r_s, n_sky, k_sky = sky_noise(ident, cache["dates"])
        r_h, n_hood, k_hood = hood_noise(ident, era)
        for r, n, lab, col_, lw in ((r_h, n_hood, f"under the hood ({k_hood} sessions)",
                                     "k", 2.0),
                                    (r_s, n_sky, f"ordinary clear nights ({k_sky})",
                                     "crimson", 1.6)):
            if n is None:
                continue
            z2n = np.where(r > 0, r, np.nan) ** 2
            m = (r >= 100) & (r <= 15000)
            ax.plot((n / z2n)[m], r[m] / 1e3, "-", color=col_, lw=lw, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel("noise per profile, P-view (rcs/z²)")
        if col == 0:
            ax.set_ylabel("THE NOISE (random scatter)\naltitude (km)")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(0, 15)
        ax.grid(alpha=0.25, which="both")
        # closure number where the sky signal is negligible
        if n_sky is not None and n_hood is not None and r_s.size == r_h.size:
            hi = r_s >= 0.75 * r_s.max()
            ratio = float(np.nanmedian(n_sky[hi] / n_hood[hi]))
            ax.set_title(f"sky/hood noise ratio (top gates): ×{ratio:.2f}", fontsize=10)

    fig.suptitle("The two 'electronic' quantities, per instrument type — the DARK (a fixed "
                 "bias, row 1: hard from sky) and the NOISE (a random scatter, row 2: retrieved "
                 "from sky and validated against the hood)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    f = FIG / "fig13_explainer.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
