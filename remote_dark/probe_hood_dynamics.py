# -*- coding: utf-8 -*-
"""What actually modulates the dark? Empirical probe of the two long CL31 hood sessions.

Decides the M3 design before any regression is written:
  * does the solar background reach the detector under the hood (is bckgrd swept, or flat)?
  * how strongly does the dark respond to internal temperature / background, and WHERE in range?
Prints numbers; writes one landscape figure per session.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark import hood
from remote_dark.common import ensure_out


def probe(ident: str = "B"):
    out = ensure_out("probe")
    for s in hood.session_frames(ident, min_hours=10):
        t = s["times"]
        hrs = np.array([(x - t[0]).total_seconds() / 3600.0 for x in t])
        dark, rng = s["dark"], s["rng"]
        bg, ti = s["hk"]["bckgrd"], s["hk"]["t_int"]
        lab = f"{ident} {s['t0']:%Y%m%d} ({s['era']})"
        print(f"== {lab}: {dark.shape[0]} profiles over {hrs[-1]:.1f} h")
        print(f"   bckgrd: min {np.nanmin(bg):.4g}  max {np.nanmax(bg):.4g}  "
              f"(x{np.nanmax(bg)/max(np.nanmin(bg),1e-30):.1f})")
        print(f"   t_int : min {np.nanmin(ti):.1f}  max {np.nanmax(ti):.1f} degC")

        # per-gate correlation of the dark with each modulator (5-min averages kill shot noise)
        n_avg = max(1, int(round(300.0 / max(np.median(np.diff(hrs)) * 3600.0, 1e-6))))
        nb = dark.shape[0] // n_avg
        dsub = dark[:nb * n_avg].reshape(nb, n_avg, -1).mean(axis=1)
        bsub = bg[:nb * n_avg].reshape(nb, n_avg).mean(axis=1)
        tsub = ti[:nb * n_avg].reshape(nb, n_avg).mean(axis=1)
        hsub = hrs[:nb * n_avg].reshape(nb, n_avg).mean(axis=1)

        def corr_profile(x):
            xc = x - np.nanmean(x)
            d = dsub - np.nanmean(dsub, axis=0)
            sx = np.nanstd(x)
            sd = np.nanstd(dsub, axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                return np.nanmean(d * xc[:, None], axis=0) / (sx * sd)
        c_bg, c_ti = corr_profile(bsub), corr_profile(tsub)
        for name, c in (("bckgrd", c_bg), ("t_int", c_ti)):
            band = lambda lo, hi: float(np.nanmedian(np.abs(c[(rng >= lo) & (rng <= hi)])))
            print(f"   |corr(dark, {name:6s})| median: 60-500m {band(60,500):.2f}   "
                  f"0.5-2km {band(500,2000):.2f}   2-7km {band(2000,7000):.2f}")

        # figure: modulators + the dark's drift at three heights (landscape)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        ax = axes[0]
        ax.plot(hsub, bsub, "k-", lw=1)
        ax.set_xlabel("hours"); ax.set_ylabel("bckgrd_rcs_0"); ax.set_title(f"{lab} — background")
        ax2 = ax.twinx(); ax2.plot(hsub, tsub, "r-", lw=1); ax2.set_ylabel("t_int (degC)", color="r")
        ax = axes[1]
        for lo, hi in ((60, 300), (500, 1000), (2000, 4000)):
            m = (rng >= lo) & (rng <= hi)
            ax.plot(hsub, np.nanmean(dsub[:, m] / rng[m] ** 2, axis=1) * 1e9,
                    lw=1, label=f"{lo}-{hi} m")
        ax.set_xlabel("hours"); ax.set_ylabel("P = dark/z$^2$ (x1e-9)"); ax.legend(fontsize=8)
        ax.set_title("dark drift by band (P-view)")
        ax = axes[2]
        ax.plot(c_bg, rng / 1e3, lw=1, label="corr vs bckgrd")
        ax.plot(c_ti, rng / 1e3, lw=1, label="corr vs t_int")
        ax.axvline(0, color="k", lw=0.5); ax.set_ylim(0, 7.7)
        ax.set_xlabel("correlation"); ax.set_ylabel("range (km)"); ax.legend(fontsize=8)
        ax.set_title("per-gate correlation")
        fig.tight_layout()
        f = out / f"hood_dynamics_{ident}_{s['t0']:%Y%m%d}.png"
        fig.savefig(f, dpi=130); plt.close(fig)
        print(f"   fig -> {f}")


if __name__ == "__main__":
    probe("B")
