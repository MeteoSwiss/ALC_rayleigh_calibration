# -*- coding: utf-8 -*-
"""Fig 6 — the M3 joint decomposition: D_slow and D_fast per instrument-era, from the hood.

Landscape, range on Y. P-view (coupling/z²) so the near-range structure is visible; the CL61
panel keeps its own scale (its couplings are ~12 orders below the CL31's — that IS the result)."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import OUT_DIR, ensure_out

FIG = ensure_out("figs")


def main():
    z = np.load(OUT_DIR / "m3" / "hood_joint.npz", allow_pickle=True)
    keys = sorted({k.rsplit("_", 1)[0] for k in z.files if k.endswith("_rng")})
    fig, axes = plt.subplots(1, len(keys), figsize=(4.2 * len(keys), 5.0))
    for ax, key in zip(np.atleast_1d(axes), keys):
        rng = z[f"{key}_rng"]
        with np.errstate(invalid="ignore", divide="ignore"):
            ds = z[f"{key}_D_slow"] / rng ** 2
            df = z[f"{key}_D_fast"] / rng ** 2
        m = (rng >= 60) & (rng <= 3000)
        sc = float(np.nanmedian(np.abs(ds[m & (rng <= 500)]))) * 6 + 1e-300
        ax.plot(ds[m] / sc, rng[m] / 1e3, color="crimson", lw=1.4, label="D_slow (per °C)")
        ax.plot(df[m] / sc, rng[m] / 1e3, color="steelblue", lw=1.1, label="D_fast (per B)")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(0, 3)
        near = m & (rng <= 500)
        ax.set_title(f"{key}\n|D_slow| 60–500 m = "
                     f"{np.nanmedian(np.abs(z[f'{key}_D_slow'][near])):.3g} rcs/°C", fontsize=9)
        ax.set_xlabel("coupling / z² (norm.)")
        ax.set_ylabel("range (km)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("M3 — the dark's modulation couplings under the hood: slow (thermal) and fast "
                 "(background), per instrument-era", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig6_m3_joint_couplings.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
