# -*- coding: utf-8 -*-
"""Fig 17 — CL61 undershoot deconvolution: turn opaque-cloud scenes into quiescent-dark
measurements by removing the signal-induced part of the AC-coupling response.

Physics. The CL61 receiver is AC-coupled. For a single-pole high-pass with time constant
tau_u, output y and input x obey  y = x - w  with the coupling state  dw/dt = (x - w)/tau_u
= y/tau_u  —  so the state is EXACTLY the plain integral of the RECORDED signal:

    w(z) = (1/L_u) * integral_0^z y(z') dz'          L_u = c*tau_u/2 = 4570 m  (hood fit)

After the huge cloud return, w is charged high and relaxes: that is the deepened undershoot
fig16 measured (x2.1 the hood's quiescent state). The inversion is per profile and closed-form:
x_hat = y + w. The positive fast lobe (hood fit: L_p = 711 m) is modelled alongside as an
exponentially-forgetting low-pass of y, coefficient fitted.

Estimators compared against the hood truth (all in the 2.5-7 km judgement band):
  raw       median above-cloud signal              (fig16's curve, theta +2.1)
  deconv    median of  y + w                       (slow pole inverted, coefficient forced 1)
  regress   per-gate fit over profiles  rcs = b(z) + a1*U1 + a2*U2,  U1 = w*z^2 (slow state),
            U2 = v*z^2 (fast state) — cloud-brightness variation between profiles is the
            leverage; b(z) is the quiescent dark with BOTH signal-induced terms projected out,
            and a1 ~ 1 is the physics check (the state model predicts unit coefficient).

Selection: night, CBH1 500-1500 m (>=500 m keeps the charge integral out of the incomplete-
overlap zone), no second layer below 6 km; regression targets only at z >= CBH + 800 m.
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
L_U = 4570.0                #: slow undershoot decay length (hood physical-model fit, 30.5 us)
L_P = 711.0                 #: fast lobe decay length (hood fit, 4.7 us)
CBH_MIN, CBH_MAX = 500.0, 1500.0
BUFFER_M = 800.0
MIN_PROF_DAY = 100
MAX_DAYS = 30
JBAND = (2500.0, 7000.0)


def collect_profiles():
    """All qualifying overcast night profiles: full rcs (for the state integrals) + CBH."""
    rcs_all, cbh_all, rng = [], [], None
    day, n_days = D0, 0
    while day <= D1 and n_days < MAX_DAYS:
        f = l1_file(PAYERNE["wmo"], "C", day)
        day += timedelta(days=1)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
                cbh = np.ma.filled(ds.variables["cloud_base_height"][:].astype("f8"), np.nan)
                times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
                el = solar_elevation_deg(PAYERNE["lat"], PAYERNE["lon"], times)
                sel = ((el < -6.0) & (cbh[:, 0] >= CBH_MIN) & (cbh[:, 0] < CBH_MAX)
                       & ~((cbh[:, 1] > 0) & (cbh[:, 1] < 6000.0)))
                if sel.sum() < MIN_PROF_DAY:
                    continue
                rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
                rcs_all.append(np.ma.filled(
                    ds.variables["rcs_0"][sel, :].astype("f4"), np.nan))
                cbh_all.append(cbh[sel, 0])
                n_days += 1
        except Exception:                                                   # noqa: BLE001
            continue
    return rng, np.vstack(rcs_all), np.concatenate(cbh_all), n_days


def states(rcs, rng):
    """Per-profile AC-coupling states from the recorded signal, vectorised over profiles.

    slow: w(z) = (1/L_u) * cumulative integral of y   (exact single-pole state)
    fast: v(z) = (1/L_p) * integral of y * exp(-(z-z')/L_p)   (forgetting integral)
    """
    z = np.where(rng > 0, rng, np.nan)
    y = np.nan_to_num(rcs / z[None, :] ** 2)
    dz = np.median(np.diff(rng))
    w = np.cumsum(y, axis=1) * dz / L_U
    k = np.exp(-dz / L_P)
    v = np.empty_like(y)
    acc = np.zeros(y.shape[0])
    for j in range(y.shape[1]):
        acc = acc * k + y[:, j] * dz / L_P
        v[:, j] = acc
    return w, v


def smooth(x, k=9):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def judge(est, bh, rng):
    m = ((rng >= JBAND[0]) & (rng <= JBAND[1]) & np.isfinite(est) & np.isfinite(bh))
    th = float(np.nansum(est[m] * bh[m]) / np.nansum(bh[m] ** 2))
    rms = float(np.sqrt(np.nanmean((est[m] - bh[m]) ** 2)))
    rms_h = float(np.sqrt(np.nanmean(bh[m] ** 2)))
    return th, rms / rms_h


def main():
    rng, rcs, cbh, n_days = collect_profiles()
    print(f"{rcs.shape[0]} profiles / {n_days} days")
    z2 = np.where(rng > 0, rng, np.nan) ** 2
    w, v = states(rcs, rng)
    above = rng[None, :] >= (cbh[:, None] + BUFFER_M)
    R = np.where(above, rcs, np.nan)

    # --- estimator 1: raw above-cloud median (fig16) -----------------------------------------
    raw = np.nanmedian(R, axis=0)
    # --- estimator 2: closed-form slow-pole inversion, coefficient forced to 1 ---------------
    deconv = np.nanmedian(np.where(above, rcs + w * z2[None, :], np.nan), axis=0)
    # --- estimator 3: per-gate regression, both states, cloud-brightness leverage ------------
    U1 = np.where(above, w * z2[None, :], np.nan)
    U2 = np.where(above, v * z2[None, :], np.nan)
    ok = np.isfinite(R) & np.isfinite(U1) & np.isfinite(U2)
    b_fit = np.full(rng.size, np.nan)
    a1_fit = np.full(rng.size, np.nan)
    for j in range(rng.size):
        m = ok[:, j]
        if m.sum() < 200:
            continue
        A = np.column_stack([U1[m, j], U2[m, j], np.ones(int(m.sum()))])
        try:
            beta = np.linalg.lstsq(A, R[m, j], rcond=None)[0]
            a1_fit[j], b_fit[j] = beta[0], beta[2]
        except np.linalg.LinAlgError:
            pass

    trng, b_rcs, _sem, _bp = hood.truth("C")
    bh = np.interp(rng, trng, b_rcs, left=np.nan, right=np.nan)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.2),
                                  gridspec_kw={"width_ratios": [1.6, 1.0]})
    mm = (rng >= 2000) & (rng <= 12000)
    curves = [(raw, "raw above-cloud (fig16)", "0.55", 1.4),
              (deconv, "deconvolved:  y + w  (slow pole, coeff = 1)", "darkorange", 1.7),
              (b_fit, "regression intercept b(z)  (both states)", "crimson", 1.9)]
    print(f"{'estimator':>42s} {'theta':>7s} {'RMS/|hood|':>11s}")
    for est, lab, colr, lw in curves:
        th, rr = judge(est, bh, rng)
        ax.plot(smooth(est / z2)[mm], rng[mm] / 1e3, "-", color=colr, lw=lw,
                label=f"{lab}  [θ={th:+.2f}, RMS {rr:.2f}]")
        print(f"{lab:>42.42s} {th:7.2f} {rr:11.2f}")
    ax.plot(smooth(bh / z2)[mm], rng[mm] / 1e3, "k-", lw=2.2, label="hood truth")
    ax.axvline(0, color="k", lw=0.5)
    ref = smooth(bh / z2)[mm]
    span = np.nanmax(np.abs(ref[np.isfinite(ref)]))
    ax.set_xlim(-4 * span, 4 * span)
    ax.set_xlabel("offset, P-view (rcs/z²)")
    ax.set_ylabel("altitude (km)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(f"CL61 (Payerne C) — {rcs.shape[0]} profiles, {n_days} overcast days\n"
                 "removing the signal-induced undershoot recovers the quiescent dark",
                 fontsize=10)

    mm2 = (rng >= 2500) & (rng <= 10000)
    ax2.plot(a1_fit[mm2], rng[mm2] / 1e3, "-", color="seagreen", lw=1.2)
    ax2.axvline(1.0, color="k", lw=0.8, ls="--", label="state model predicts a1 = 1")
    med_a1 = float(np.nanmedian(a1_fit[mm2]))
    ax2.axvline(med_a1, color="seagreen", lw=0.8, ls=":",
                label=f"median a1 = {med_a1:+.2f}")
    ax2.set_xlim(-1, 3)
    ax2.set_xlabel("fitted slow-state coefficient a1(z)")
    ax2.grid(alpha=0.25)
    ax2.legend(fontsize=8, loc="upper right")
    ax2.set_title("physics check: the fitted coefficient of the\n"
                  "exact AC-coupling state (target 1)", fontsize=10)

    fig.suptitle("CL61 undershoot deconvolution — the cloud return charges the AC coupling; "
                 "its state w(z) = ∫y/L_u is computed per profile from the data itself and "
                 "removed", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig17_cl61_deconv.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
