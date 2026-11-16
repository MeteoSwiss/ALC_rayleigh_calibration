# -*- coding: utf-8 -*-
"""Jenoptik slide: a 1 nm passband drifting at the JENOFILT spec (0.02 nm/K) against a fixed laser line.

Left: Gaussian passband (FWHM 1 nm, the Lufft datasheet value) at night and after a +6 K rise (summer
noon median of temp_int), shifted by 0.02 nm/K; the laser line sits at the fitted detuning
Delta0/W = 0.19 (q11 median). Right: the filter drift the fit requires for W = 1 nm (alpha/W IQR
0.015-0.024 /K over the units) against the Jenoptik JENOFILT datasheet coefficient.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_illustrations import C, FIG, _math_rc, new_figure, save  # noqa: E402


def jenoptik_filter():
    import matplotlib.pyplot as plt
    _math_rc()
    fig0, _ = new_figure("figure")
    plt.close(fig0)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.8, 3.3), gridspec_kw={"width_ratios": [1.15, 1]})
    W, alpha, d0, dT = 1.0, 0.020, 0.19, 6.0          # nm, nm/K, fitted detuning (W units), summer noon rise
    x = np.linspace(-1.3, 1.3, 400)

    def trans(c):
        return np.exp(-4 * np.log(2) * ((x - c) / W) ** 2)
    t0 = np.exp(-4 * np.log(2) * d0 ** 2)
    t1 = np.exp(-4 * np.log(2) * (d0 - alpha * dT / W) ** 2)
    ax.plot(x, trans(0), color=C["blue_dark"], lw=2.2, label="filter at night")
    ax.plot(x, trans(alpha * dT), color=C["brown"], lw=2.2, ls="--", label=f"filter +{dT:.0f} K (noon)")
    ax.axvline(d0, color=C["red"], lw=2.5)
    ax.text(d0 + 0.05, 0.25, "laser line\n(fixed)", color=C["red"], fontsize=10)
    ax.annotate("", xy=(d0, t1), xytext=(d0, t0), arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
    ax.text(0.42, 1.06, f"+{100 * (t1 / t0 - 1):.0f} % at the laser line", fontsize=10.5, fontweight="bold")
    ax.annotate("", xy=(alpha * dT, 1.04), xytext=(0, 1.04),
                arrowprops=dict(arrowstyle="->", color=C["brown"], lw=1.5))
    ax.text(-0.04, 1.06, f"{1000 * alpha * dT:.0f} pm", fontsize=9.5, color=C["brown"], ha="right")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("wavelength − filter centre at night [nm]", fontsize=10)
    ax.set_ylabel("filter transmission", fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(-0.01, 0.62))
    ax.tick_params(labelsize=9)
    # right: drift required by the fit vs the Jenoptik specification
    ax2.axvspan(15, 24, color=C["red"], alpha=0.18, lw=0)
    ax2.axvline(20, color=C["red"], lw=2)
    ax2.text(14.3, 1.55, "required by the fit\n(W = 1 nm, IQR of units)", color=C["red"], fontsize=9.5,
             va="center", ha="right")
    ax2.plot([20], [0.6], "o", color=C["brown"], ms=11, zorder=3)
    ax2.text(25, 0.6, "Jenoptik\nJENOFILT spec\n0.02 nm/K", color=C["brown"], fontsize=10, ha="left", va="center")
    ax2.set_xlim(0, 32)
    ax2.set_ylim(-0.6, 2.1)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.set_xlabel("filter drift dλ/dT [pm/K]", fontsize=10)
    top = ax2.secondary_xaxis("top", functions=(lambda a: 0.1054 * a, lambda k: k / 0.1054))
    top.set_xlabel("predicted κ for a 1 nm filter [%/K]", fontsize=9.5)
    top.tick_params(labelsize=9)
    ax2.tick_params(labelsize=9)
    fig.tight_layout(w_pad=2.5)
    save(fig, FIG / "fig_jenoptik_filter.png")


if __name__ == "__main__":
    jenoptik_filter()
