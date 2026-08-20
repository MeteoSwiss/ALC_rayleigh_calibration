# -*- coding: utf-8 -*-
"""Fig 20 — what drives the SLOW INCREASE of the CHM15k signal during a hood (dark) session?
(operator: "is it the recovery of the internal reflection, or other physical phenomena?")

The hood sessions carry time-resolved housekeeping alongside the dark frames (internal
temperature, laser power, background), so the driver is identifiable instead of guessable:

  - thermal:      the drift tracks t_int (M3 already measured D_slow ~ 21 rcs/degC at 60-500 m)
  - laser power:  the drift tracks the laser monitor
  - settling:     the drift is a relaxation from the pre-covering sky exposure (decays with
                  hours-since-start regardless of housekeeping)
  - internal reflection: the hood/window reflection is an "opaque cloud at 0 m" - if its
                  response drives the drift, the drift's ALTITUDE SHAPE should look like the
                  fig19 cloud-response kernel k(z') (taken at z' = z), and its amplitude should
                  follow the laser (the reflection's charge is proportional to emitted power).

Panels: the two longest A sessions as time series (5-min bins, z-scored bands + drivers), and
the altitude shape of the drift (per-gate regression on t_int) against the two candidate
shapes (hood dark itself vs the fig19 kernel).
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import DATA, ensure_out
from remote_dark import hood

FIG = ensure_out("figs")
BANDS = ((60.0, 500.0), (1000.0, 3000.0), (5000.0, 10000.0))
BIN_MIN = 5.0


def zscore(x):
    s = np.nanstd(x)
    return (x - np.nanmean(x)) / (s if s > 0 else 1.0)


def binned(hours, x, width_min=BIN_MIN):
    edges = np.arange(0, np.nanmax(hours) + 1e-9, width_min / 60.0)
    idx = np.digitize(hours, edges) - 1
    out = np.full(edges.size, np.nan)
    for i in range(edges.size):
        m = idx == i
        if m.sum() >= 3:
            out[i] = np.nanmedian(x[m])
    return edges + width_min / 120.0, out


def smooth(x, k=7):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sessions = sorted(hood.session_frames("A", min_hours=3.0),
                      key=lambda s: -(s["t1"] - s["t0"]).total_seconds())[:2]
    fig, axes = plt.subplots(1, len(sessions) + 1, figsize=(15.5, 6.2))

    beta_T = None
    rng_ref = None
    for ax, s in zip(axes, sessions):
        rng = s["rng"]
        z2 = np.where(rng > 0, rng, np.nan) ** 2
        P = s["dark"] / z2[None, :]
        hours = np.array([(t - s["times"][0]).total_seconds() / 3600.0
                          for t in s["times"]])
        drivers = {"t_int": s["hk"]["t_int"], "laser": s["hk"]["laser"]}
        corrs = {}
        for (b0, b1), colr in zip(BANDS, ("crimson", "darkorange", "steelblue")):
            m = (rng >= b0) & (rng <= b1)
            band = np.nanmedian(P[:, m], axis=1)
            hb, xb = binned(hours, band)
            lab = f"{b0 / 1e3:.2g}-{b1 / 1e3:.2g} km"
            ax.plot(hb, zscore(xb), "-", color=colr, lw=1.5, label=lab)
            for dn, dv in drivers.items():
                _, db = binned(hours, dv)
                ok = np.isfinite(xb) & np.isfinite(db)
                if ok.sum() > 10:
                    corrs[(lab, dn)] = float(np.corrcoef(xb[ok], db[ok])[0, 1])
            # absolute drift over the session, in % of the hood-dark scale in that band
            first = np.nanmedian(xb[:max(3, int(60 / BIN_MIN))])
            last = np.nanmedian(xb[-max(3, int(60 / BIN_MIN)):])
            print(f"  {s['t0']:%Y-%m-%d} band {lab}: drift {last - first:+.2e} P-view "
                  f"({100 * (last - first) / max(abs(first), 1e-300):+.0f} % of its level)")
        for dn, ls in (("t_int", "--"), ("laser", ":")):
            hb, db = binned(hours, drivers[dn])
            ax.plot(hb, zscore(db), ls, color="k", lw=1.3, label=dn)
        dur = (s["t1"] - s["t0"]).total_seconds() / 3600.0
        cstr = "  ".join(f"r({k[0]},{k[1]})={v:+.2f}" for k, v in corrs.items()
                         if k[1] == "t_int")
        ax.set_title(f"session {s['t0']:%Y-%m-%d} ({dur:.0f} h)\n{cstr}", fontsize=9)
        ax.set_xlabel("hours under the hood")
        ax.set_ylabel("z-score")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, loc="best")
        print(f"session {s['t0']:%Y-%m-%d}: " +
              "  ".join(f"r({k[0]}|{k[1]})={v:+.2f}" for k, v in sorted(corrs.items())))

        if beta_T is None:
            # altitude shape of the drift: per-gate slope of P on t_int (5-min bins)
            hb, tb = binned(hours, drivers["t_int"])
            Pb = np.full((hb.size, rng.size), np.nan)
            idx = np.digitize(hours, np.arange(0, np.nanmax(hours) + 1e-9, BIN_MIN / 60.0)) - 1
            for i in range(hb.size):
                m = idx == i
                if m.sum() >= 3:
                    Pb[i] = np.nanmedian(P[m], axis=0)
            ok = np.isfinite(tb)
            tc = tb[ok] - np.nanmean(tb[ok])
            Yc = Pb[ok] - np.nanmean(Pb[ok], axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                beta_T = np.nansum(tc[:, None] * Yc, axis=0) / np.nansum(tc ** 2)
            rng_ref = rng

    # ---- shape panel: what does the drift look like in altitude? -----------------------------
    ax = axes[-1]
    z2 = np.where(rng_ref > 0, rng_ref, np.nan) ** 2
    mm = (rng_ref >= 100) & (rng_ref <= 10000)

    def norm(x):
        s = np.nanmax(np.abs(smooth(x)[mm]))
        return smooth(x) / (s if s > 0 else 1.0)

    trng, b_rcs, _sem, _bp = hood.truth("A")
    hood_p = np.interp(rng_ref, trng, b_rcs / np.where(trng > 0, trng, np.nan) ** 2)
    kz = None
    try:
        from remote_dark.exp_chm_response import CACHE
        from remote_dark.exp_cloud_impact import ZREL
        if CACHE.exists():
            zc = np.load(CACHE)
            Y, Q = zc["Y"], zc["Q"]
            t1, t2 = np.nanpercentile(Q, [33.3, 66.7])
            kern = np.nanmedian(Y[Q >= t2], axis=0) - np.nanmedian(Y[Q <= t1], axis=0)
            kz = np.interp(rng_ref, ZREL, kern, left=np.nan, right=np.nan)
    except Exception:                                                       # noqa: BLE001
        pass
    ax.plot(norm(beta_T)[mm], rng_ref[mm] / 1e3, "-", color="seagreen", lw=2.0,
            label="drift shape β_T(z) (dP/dT_int)")
    ax.plot(norm(hood_p)[mm], rng_ref[mm] / 1e3, "k--", lw=1.6, label="hood dark shape")
    if kz is not None:
        ax.plot(norm(kz)[mm], rng_ref[mm] / 1e3, "-", color="purple", lw=1.4, alpha=0.8,
                label="cloud-response kernel (fig19, z'=z)")
        for other, lab in ((hood_p, "hood"), (kz, "kernel")):
            o = np.isfinite(beta_T) & np.isfinite(other) & mm
            r = float(np.corrcoef(beta_T[o], other[o])[0, 1])
            print(f"shape corr(beta_T, {lab}) over 0.1-10 km: {r:+.2f}")
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("shape (normalised)")
    ax.set_ylabel("altitude (km)")
    ax.set_title("the drift's altitude shape vs the\ntwo candidate origins", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best")

    fig.suptitle("CHM15k hood-session slow drift — driver identification (time series vs "
                 "housekeeping) and altitude-shape arbitration (thermal dark vs internal-"
                 "reflection response)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig20_hood_drift.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
