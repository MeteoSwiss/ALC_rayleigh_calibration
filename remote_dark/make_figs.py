# -*- coding: utf-8 -*-
"""Key figures of the M1+M3 session. Landscape, range on the Y axis everywhere."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark import hood, models
from remote_dark.common import OUT_DIR, ensure_out

FIG = ensure_out("figs")


def _norm(x, m):
    s = float(np.nanmax(np.abs(x[m]))) if np.isfinite(x[m]).any() else 1.0
    return x / (s if s > 0 else 1.0)


# ---------------------------------------------------------------- fig 1: families vs hood truth --
def fig_families():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    for ax, ident, itype in zip(axes, "BAC", ("CL31", "CHM15k", "CL61")):
        rng, b, sem, braw = hood.truth(ident)
        fam = models.fit_family(itype, rng, braw, sem)
        lo, hi = models.FIT_BAND[itype]
        m = (rng >= lo) & (rng <= hi) & (rng > 0)
        # P-view everywhere: it is where the electronics live and where each structure is visible.
        # x-limits from the STRUCTURE band, not from the near-range spike that dwarfs everything.
        x_t = braw / np.where(rng > 0, rng, np.nan) ** 2
        x_f = np.where(m, fam.p_view(rng), np.nan)
        if itype == "CL31":
            band, ylim = (rng >= 150) & (rng <= 3000), (0, 3.0)
            score = f"R2(P, 60-1500 m) = {fam.meta['r2_p_near']:.3f}"
        elif itype == "CHM15k":
            band, ylim = (rng >= 300) & (rng <= 12000), (0, 12.0)
            score = f"R2(b,struct) = {fam.meta['r2_b_struct']:.3f} - ray-band {fam.meta['ray_pct']:.0f}%"
        else:
            band, ylim = (rng >= 10) & (rng <= 1500), (0, 1.5)
            score = "template (exact by construction)"
        sc = float(np.nanmax(np.abs(x_f[band & m]))) if np.isfinite(x_f[band & m]).any() else 1.0
        ax.plot(x_t / sc, rng / 1e3, color="0.55", lw=1.0, label="hood truth (raw)")
        ax.plot(x_f / sc, rng / 1e3, color="crimson", lw=1.8, label=f"family ({fam.kind})")
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(*ylim)
        ax.set_xlabel("P = b/z$^2$, normalised on the structure band")
        ax.set_title(f"{itype}  (Payerne {ident})\n{score}", fontsize=10)
        ax.set_ylabel("range (km)")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.25)
    fig.suptitle("M1 step 1 - parametric families fitted to the Payerne hood truth "
                 "(the ceiling of what the sky fit can recover)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    f = FIG / "fig1_families_vs_hood.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    return f


# ------------------------------------------------------------------- fig 2: M1 sky vs hood shape --
def fig_m1():
    runs = [("0-20000-0-06610_A_m1.npz", "CHM15k (Payerne A)"),
            ("0-20000-0-06610_C_m1.npz", "CL61 (Payerne C)")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax, (name, lab) in zip(axes, runs):
        z = np.load(OUT_DIR / "m1" / name, allow_pickle=True)
        rng, bh, bf, hb = z["rng"], z["b_hat"], z["b_free"], z["hood_b"]
        band = (rng >= 2000) & (rng <= 8000)
        m = band & np.isfinite(bh) & np.isfinite(hb)
        cc = float(np.corrcoef(bh[m], hb[m])[0, 1])
        th = float(z["theta_corr"]) if "theta_corr" in z.files else float("nan")
        ax.plot(_norm(hb, m)[band], rng[band] / 1e3, color="0.2", lw=1.8, label="hood truth")
        ax.plot(_norm(bh, m)[band], rng[band] / 1e3, color="crimson", lw=1.6,
                label="M1 anchored family (sky)")
        mf = band & np.isfinite(bf)
        ax.plot(_norm(bf, mf)[band], rng[band] / 1e3, color="steelblue", lw=0.9, alpha=0.7,
                label="free b (v4 evidence)")
        ax.set_title(f"{lab} — {int(z['n_nights'])} nights\n"
                     f"shape corr vs hood (2–8 km): {cc:+.2f}   "
                     f"[amplitude: open units bug]", fontsize=10)
        ax.set_xlabel("dark, normalised on 2–8 km")
        ax.set_ylabel("range (km)")
        ax.set_ylim(2, 8)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("M1 — sky-retrieved dark vs hood truth, Rayleigh band (shape comparison)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    f = FIG / "fig2_m1_sky_vs_hood.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    return f


# -------------------------------------------------------- fig 3: M3 twilight — the honest result --
def fig_twilight():
    z = np.load(OUT_DIR / "twilight" / "0-20000-0-06610_B.npz", allow_pickle=True)
    rng = z["rng"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax = axes[0]
    m = (rng >= 60) & (rng <= 2000)
    ax.plot(z["sky_D_dawn"][m], rng[m] / 1e3, lw=1.2, label=f"sky dawn (n={int(z['sky_n_events'])//2})")
    ax.plot(z["sky_D_dusk"][m], rng[m] / 1e3, lw=1.2, label="sky dusk")
    ax.plot(z["sky_D_med"][m], rng[m] / 1e3, "k-", lw=1.8, label="sky median")
    ax.set_title("sky twilight D_bg(z) — dawn/dusk consistency", fontsize=10)
    ax.set_xlabel("D_bg (rcs_0 per B-unit)")
    ax.set_ylabel("range (km)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax = axes[1]
    for k in z.files:
        if k.startswith("hood_") and k.endswith("_D_med"):
            lab = k.replace("hood_", "").replace("_D_med", "")
            ax.plot(z[k][m], rng[m] / 1e3, lw=1.4, label=f"hood {lab}")
    ax.plot(z["sky_D_med"][m], rng[m] / 1e3, "k-", lw=1.8, label="sky median")
    ax.set_title("sky vs hood truth — corr ≈ −0.1..−0.2: NOT validated\n"
                 "(the negative result: terminator is BL-locked; fast-B truth weak under hood)",
                 fontsize=10)
    ax.set_xlabel("D_bg (rcs_0 per B-unit)")
    ax.set_ylabel("range (km)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.suptitle("M3 — twilight background sweep, CL31 (state: refuted as designed, "
                 "joint B+T redesign next)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig3_twilight_state.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    return f


if __name__ == "__main__":
    for fn in (fig_families, fig_m1, fig_twilight):
        print(fn())
