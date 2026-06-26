"""Per-station OmB figure (2x3, faithful to E_PROFILE_ALC_Monthly_OB.m).

Top: observation pcolor (our C_L) | median profiles (op / ours / ours-WV / CAMS).
Bottom: CAMS aerosol-backscatter pcolor | bias (O - B) profiles.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .omb import OmBResult  # noqa: E402

_COL = {"ours": (0.85, 0.10, 0.10), "op": (0.25, 0.25, 0.25),
        "ours_wv": (0.10, 0.30, 0.85), "cams": (0.0, 0.0, 0.0)}


def plot_omb_station(res: OmBResult, instrument: str, save_path, *, title: str = "",
                     unit: float = 1e6, range_top: float = 15000.0) -> str:
    """Render the per-station OmB figure to *save_path*. ``unit`` converts the
    internal m^-1 sr^-1 to the display unit (1e6 -> Mm^-1 sr^-1)."""
    have_wv = "ours_wv" in res.obs_interp
    z = res.z_cams
    t = res.time_cams
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3)

    # (0,0:2) observation pcolor (our-calibrated), log10
    ax = fig.add_subplot(gs[0, 0:2])
    obs = res.obs_mean["ours"].T * unit
    pcm = ax.pcolormesh(t, res.range_mean, np.log10(np.clip(obs, 1e-2, None)),
                        cmap="jet", vmin=-2, vmax=2, shading="auto")
    ax.set_ylim(0, range_top); ax.set_ylabel("Range AGL [m]")
    ax.set_title(f"{instrument} observation (v2 C_L), {res.wavelength:.0f} nm")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    fig.colorbar(pcm, ax=ax, label=r"log$_{10}\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]")

    # (0,2) median profiles
    ax = fig.add_subplot(gs[0, 2])
    ax.fill_betweenx(z, res.prof["ours"]["obs_p25"] * unit,
                     res.prof["ours"]["obs_p75"] * unit, color=_COL["ours"], alpha=0.15)
    ax.plot(res.prof["ours"]["obs_med"] * unit, z, "-", color=_COL["ours"], label="obs (v2)")
    ax.plot(res.prof["op"]["obs_med"] * unit, z, "-", color=_COL["op"], label="obs (op)")
    if have_wv:
        ax.plot(res.prof["ours_wv"]["obs_med"] * unit, z, "--", color=_COL["ours_wv"],
                label="obs (v2, WV)")
    ax.plot(res.cams_med * unit, z, "-", color=_COL["cams"], lw=2, label="CAMS")
    ax.set_ylim(0, range_top); ax.set_xlim(left=0)
    ax.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]"); ax.set_ylabel("Altitude ASL [m]")
    ax.set_title("Median profiles"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (1,0:2) CAMS pcolor
    ax = fig.add_subplot(gs[1, 0:2])
    pcm = ax.pcolormesh(t, z, np.log10(np.clip(res.cams_beta * unit, 1e-2, None)),
                        cmap="jet", vmin=-2, vmax=2, shading="auto")
    ax.set_ylim(0, range_top); ax.set_ylabel("Altitude ASL [m]"); ax.set_xlabel("Day of month")
    ax.set_title(f"CAMS aerosol backscatter at {res.wavelength:.0f} nm")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    fig.colorbar(pcm, ax=ax, label=r"log$_{10}\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]")

    # (1,2) bias profiles
    ax = fig.add_subplot(gs[1, 2])
    for src, lab in (("op", "op - CAMS"), ("ours", "v2 - CAMS"),
                     ("ours_wv", "v2 WV - CAMS")):
        if src in res.prof:
            ax.plot(res.prof[src]["bias_med"] * unit, z, color=_COL[src], label=lab,
                    ls="--" if src == "ours_wv" else "-")
    ax.fill_betweenx(z, res.prof["ours"]["bias_p25"] * unit,
                     res.prof["ours"]["bias_p75"] * unit, color=_COL["ours"], alpha=0.12)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_ylim(0, range_top); ax.set_xlabel(r"O - B [Mm$^{-1}$sr$^{-1}$]")
    ax.set_ylabel("Altitude ASL [m]"); ax.set_title("Bias profile"); ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(title or f"OmB - {instrument}", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150, facecolor="w")
    plt.close(fig)
    return str(save_path)
