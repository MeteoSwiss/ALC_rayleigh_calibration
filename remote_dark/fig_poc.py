# -*- coding: utf-8 -*-
"""Fig 7 — the Schiphol-quadruple POC. Landscape; range on Y for the profile panels."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import OUT_DIR, ensure_out

FIG = ensure_out("figs")


def main():
    z = np.load(OUT_DIR / "m2" / "poc_amsterdam.npz", allow_pickle=True)
    rng = z["rng"]
    node = z["node"]
    idents = [str(x) for x in z["idents"]]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))

    # (a) node dark-difference profiles
    ax = axes[0]
    m = (rng >= 1000) & (rng <= 10000)
    for k, i in enumerate(idents):
        ax.plot(node[k][m] * 1e9, rng[m] / 1e3, lw=1.2, label=f"unit {i}")
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("dark offset vs quad mean (×10⁻⁹, normalised units)")
    ax.set_ylabel("range (km)")
    ax.set_title(f"node dark differences — {int(z['n_shared'])} shared clear nights", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    # (b) the multiplicative channel vs prior knowledge
    ax = axes[1]
    mm = (rng >= 100) & (rng <= 3000)
    for pair, ls in (("AB", "-"), ("AC", "--"), ("AD", ":")):
        ax.plot(z[f"rho_{pair}"][mm], rng[mm] / 1e3, ls, lw=1.3, label=f"ρ {pair[0]}/{pair[1]}")
    ax.axvline(1.0, color="k", lw=0.6)
    ax.axvspan(0.75, 0.85, ymin=0.5 / 3, ymax=1.0 / 3, color="orange", alpha=0.18)
    ax.text(0.76, 1.05, "known B overlap deficit\n(~0.8 at 0.5–1 km)", fontsize=7, color="peru")
    ax.set_xlim(0.6, 1.4)
    ax.set_xlabel("ρ(z) — response / overlap ratio")
    ax.set_ylabel("range (km)")
    ax.set_title("multiplicative channel: A/B reproduces the\nknown Schiphol-B overlap deficit "
                 "(blind)", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")

    # (c) delta_s summary (leak-controlled)
    ax = axes[2]
    s = z["s_node"]
    ypos = np.arange(len(idents))
    ax.barh(ypos, s, height=0.5, color=["#888" if abs(v) < 0.45 else "#c0392b" for v in s])
    ax.set_yticks(ypos, [f"unit {i}" for i in idents])
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Δs vs quad mean (Payerne-dark units, two-basis projection)")
    ax.set_title("per-unit dark amplitude offsets\n(common-mode absolute NOT included — "
                 "Payerne anchors that)", fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    fig.suptitle("M2 POC — four co-located CHM15k at Schiphol: the atmosphere common-modes away, "
                 "instrument differences remain", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig7_poc_amsterdam.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
