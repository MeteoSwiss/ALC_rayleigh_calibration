# -*- coding: utf-8 -*-
"""Fig 9 — preview: does the D optical-module swap (TUB150037 -> TUB160055, 2026-07-11) show in
the noise? Local data ends 2026-07-12, so the post-swap side is 1-2 nights — a preview with its
sample size printed on it, not a verdict. The full test runs once the balfrin extension lands."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import V4_DIR, ensure_out

FIG = ensure_out("figs")


def main():
    z = np.load(V4_DIR / "cache" / "0-20000-0-06240_D_20250101_20260813_220.npz",
                allow_pickle=True)
    rng = z["rng"]
    z2 = np.where(rng > 0, rng, np.nan) ** 2
    d = z["dates"]
    n_prof = np.maximum(z["n_prof"].astype(float), 1.0)
    noise = z["sigma"] * np.sqrt(n_prof)[:, None] / z2[None, :]
    pre = d < 20260711
    post = d >= 20260711

    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    m = (rng >= 100) & (rng <= 15000)
    med_pre = np.nanmedian(noise[pre], axis=0)
    lo_pre = np.nanpercentile(noise[pre], 25, axis=0)
    hi_pre = np.nanpercentile(noise[pre], 75, axis=0)
    ax.fill_betweenx(rng[m] / 1e3, lo_pre[m], hi_pre[m], color="0.8", alpha=0.6,
                     label=f"TUB150037 p25–p75 ({int(pre.sum())} nights)")
    ax.plot(med_pre[m], rng[m] / 1e3, "k-", lw=1.6,
            label="TUB150037 median (pre-swap)")
    for k in np.where(post)[0]:
        ax.plot(noise[k][m], rng[m] / 1e3, "-", color="crimson", lw=1.1, alpha=0.85,
                label=f"TUB160055 night {int(d[k])}")
    ax.set_xscale("log")
    ax.set_xlabel("noise per profile, P-view (rcs/z²)")
    ax.set_ylabel("altitude (km)")
    ax.set_ylim(0, 15)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, loc="upper right")
    top = rng >= np.nanpercentile(rng, 92)
    f_pre = float(np.nanmedian(med_pre[top]))
    f_post = float(np.nanmedian(np.nanmedian(noise[post], axis=0)[top]))
    ax.set_title("Schiphol D — optical module swap 2026-07-11, noise view\n"
                 f"floor: {f_pre:.2e} (old) vs {f_post:.2e} (new, n={int(post.sum())} nights) "
                 f"= ×{f_post / f_pre:.2f}   [PREVIEW — full test awaits the balfrin extension]",
                 fontsize=10)
    fig.tight_layout()
    f = FIG / "fig9_d_swap_preview.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f, f"floor ratio x{f_post / f_pre:.2f}")


if __name__ == "__main__":
    main()
