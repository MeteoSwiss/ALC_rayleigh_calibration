# -*- coding: utf-8 -*-
"""Fig 16 — operator's idea: use completely obscuring low clouds at night as a NATURAL
TERMINATION HOOD. Above an opaque cloud the two-way transmission is ~zero, so whatever the
instrument records up there is its own additive background — measured per unit, per night,
network-wide, with no site visit and NO atmospheric model (the atmosphere is switched OFF, not
modelled — this walks around the aerosol identifiability wall instead of fighting it).

Known contaminants the hood does not have: multiple-scattering afterglow just above cloud top,
and signal-induced afterpulse / baseline recovery from the strong cloud return (dead time in
the cloud itself does not matter — we only use gates above). Both live just above the cloud,
so each profile is masked below its own CBH + BUFFER and the buffer sensitivity is the test.

Payerne validation: stack winter overcast nights for A (CHM15k), B (CL31), C (CL61) and compare
the cloud-derived offset with the hood truth. Selection per profile: night (sun < -6 deg),
first cloud base 100-1500 m, no second layer below 6 km (broken-sky guard).
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import PAYERNE, l1_file, solar_elevation_deg, ensure_out
from remote_dark import hood

FIG = ensure_out("figs")
D0, D1 = datetime(2025, 10, 1), datetime(2026, 4, 30)
BUFFER_M = 1200.0
CBH_MAX = 1500.0
MIN_PROF_DAY = 100
MAX_DAYS = 30
ERA_B = "pre_swap"          #: CL31 winter 2025/26 = pre-swap optic block


def day_offsets(wmo, ident, day):
    """One day's opaque-overcast night profiles, each masked below its own CBH+buffer.
    Returns (rng, per-day median offset profile, n_profiles) or None."""
    f = l1_file(wmo, ident, day)
    if not f.exists():
        return None
    try:
        with netCDF4.Dataset(f) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
            cbh = np.ma.filled(ds.variables["cloud_base_height"][:].astype("f8"), np.nan)
            times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
            el = solar_elevation_deg(PAYERNE["lat"], PAYERNE["lon"], times)
            sel = ((el < -6.0) & (cbh[:, 0] > 100.0) & (cbh[:, 0] < CBH_MAX)
                   & ~((cbh[:, 1] > 0) & (cbh[:, 1] < 6000.0)))
            if sel.sum() < MIN_PROF_DAY:
                return None
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][sel, :].astype("f8"), np.nan)
    except Exception:                                                       # noqa: BLE001
        return None
    rcs = np.where(rng[None, :] > 0, rcs, np.nan)
    lim = cbh[sel, 0][:, None] + BUFFER_M
    rcs = np.where(rng[None, :] >= lim, rcs, np.nan)
    return rng, np.nanmedian(rcs, axis=0), int(sel.sum())


def collect(ident):
    days, rng = [], None
    day = D0
    while day <= D1 and len(days) < MAX_DAYS:
        r = day_offsets(PAYERNE["wmo"], ident, day)
        day += timedelta(days=1)
        if r is None:
            continue
        rng = r[0]
        days.append((r[1], r[2]))
    if not days:
        return None, None, 0, 0
    prof = np.vstack([d[0] for d in days])
    return rng, prof, len(days), sum(d[1] for d in days)


def smooth(x, k=7):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.0))
    for ax, ident, itype in zip(axes, "ABC", ("CHM15k", "CL31", "CL61")):
        rng, prof, n_days, n_prof = collect(ident)
        trng, b_rcs, _sem, _bp = hood.truth(ident)
        if ident == "B":                      # winter nights are the PRE-swap optic block
            sess = [s for s in hood.session_frames("B", min_hours=2.0) if s["era"] == ERA_B]
            if sess:
                b_rcs = np.nanmedian(np.vstack(
                    [np.nanmedian(s["dark"], axis=0) for s in sess]), axis=0)
                trng = sess[0]["rng"]
        bh = np.interp(rng, trng, b_rcs, left=np.nan, right=np.nan) if rng is not None else None

        if rng is None:
            ax.text(0.5, 0.5, "no qualifying overcast nights", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(f"{itype} (Payerne {ident})", fontsize=10)
            continue
        z2 = np.where(rng > 0, rng, np.nan) ** 2
        med = np.nanmedian(prof, axis=0)
        lo, hi = np.nanpercentile(prof, [25, 75], axis=0)
        top = 7.7 if ident == "B" else 12.0
        mm = (rng >= 2000) & (rng <= top * 1e3)
        ax.fill_betweenx(rng[mm] / 1e3, smooth(lo / z2)[mm], smooth(hi / z2)[mm],
                         color="mistyrose", label="p25-p75 (overcast days)")
        ax.plot(smooth(med / z2)[mm], rng[mm] / 1e3, "-", color="crimson", lw=1.8,
                label=f"above opaque cloud ({n_days} d, {n_prof} prof)")
        ax.plot(smooth(bh / z2)[mm], rng[mm] / 1e3, "k-", lw=2.0, label="hood truth")
        ax.axvline(0, color="k", lw=0.5)
        # judge: amplitude of hood contained in the cloud-derived curve + RMS, 2.5-7 km
        jm = (rng >= 2500) & (rng <= 7000) & np.isfinite(med) & np.isfinite(bh)
        th = float(np.nansum(med[jm] * bh[jm]) / np.nansum(bh[jm] ** 2))
        rms = float(np.sqrt(np.nanmean((med[jm] - bh[jm]) ** 2)))
        rms_h = float(np.sqrt(np.nanmean(bh[jm] ** 2)))
        ax.set_title(f"{itype} (Payerne {ident})\n"
                     f"theta vs hood (2.5-7 km): {th:+.2f}   RMS/|hood|: {rms / rms_h:.2f}",
                     fontsize=10)
        ax.set_xlabel("offset, P-view (rcs/z²)")
        if ident == "A":
            ax.set_ylabel("altitude (km)")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.25)
        ref = smooth(bh / z2)[mm]
        span = np.nanmax(np.abs(ref[np.isfinite(ref)])) if np.isfinite(ref).any() else 1.0
        ax.set_xlim(-5 * span, 5 * span)
        print(f"{itype}: {n_days} days, {n_prof} profiles, theta={th:+.2f}, "
              f"RMS/|hood|={rms / rms_h:.2f}")
    fig.suptitle("Opaque low cloud as a natural termination hood (operator's idea) — signal "
                 "above fully-attenuating night cloud vs the real hood, Payerne winter "
                 f"2025/26 (each profile masked below its own CBH + {BUFFER_M:.0f} m)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig16_cloud_dark.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
