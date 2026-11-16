# -*- coding: utf-8 -*-
"""Fig 18 — measure the cloud's impact on the rest of the profile EMPIRICALLY, and use it as
detector physics (operator's questions, both at once).

Frame: height ABOVE the cloud base (z' = z - CBH), because any signal-induced response is
locked to the cloud, not to absolute altitude. Profiles are split by the cloud-return charge
Q = integral of the recorded signal up to CBH+400 m (brightness terciles). Then:

    bright tercile - dim tercile  =  the pure signal-induced pulse-response tail
                                     (quiescent dark AND internal-pulse response cancel)

That difference is detector physics measured from the sky: its decay length (CL61) or ripple
wavelength (CL31) can be compared with the hood-fitted circuit constants (L_u = 4570 m;
fast ring Lambda = 1053 m), and its magnitude quantifies how much the cloud contaminates the
above-cloud dark estimate of fig16 — including the CL31 mismatch below ~4 km absolute.

Lesson from fig17 carried in: the IDEAL single-pole state (w = int y / L_u, coefficient 1 in
recorded units) over-predicts the coupling by ~1e4 and its fitted coefficient is 0.00 — the
recorded, calibrated signal does not charge the coupling at face value (firmware baseline
handling / electrical dynamic-range compression in between). Hence this empirical route.
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
CBH_MIN, CBH_MAX = 500.0, 1500.0
MIN_PROF_DAY = 100
MAX_DAYS = 30
ZREL = np.arange(0.0, 5501.0, 25.0)          #: height above cloud base (m)


def collect(ident):
    """Per-profile y-view profiles on the cloud-relative grid + charge Q per profile."""
    ys, qs, cbhs = [], [], []
    day, n_days = D0, 0
    while day <= D1 and n_days < MAX_DAYS:
        f = l1_file(PAYERNE["wmo"], ident, day)
        day += timedelta(days=1)
        if not f.exists():
            continue
        try:
            with netCDF4.Dataset(f) as ds:
                t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
                cbh = np.ma.filled(ds.variables["cloud_base_height"][:].astype("f8"), np.nan)
                times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x))
                                  for x in t])
                el = solar_elevation_deg(PAYERNE["lat"], PAYERNE["lon"], times)
                sel = ((el < -6.0) & (cbh[:, 0] >= CBH_MIN) & (cbh[:, 0] < CBH_MAX)
                       & ~((cbh[:, 1] > 0) & (cbh[:, 1] < 6000.0)))
                if sel.sum() < MIN_PROF_DAY:
                    continue
                rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
                rcs = np.ma.filled(ds.variables["rcs_0"][sel, :].astype("f8"), np.nan)
        except Exception:                                                   # noqa: BLE001
            continue
        z = np.where(rng > 0, rng, np.nan)
        y = rcs / z[None, :] ** 2
        cb = cbh[sel, 0]
        dz = float(np.median(np.diff(rng)))
        for i in range(y.shape[0]):
            # charge = the CLOUD return only — integrating from 0 would fold the CL31's own
            # near-range dark lobe (negative!) into Q and corrupt the brightness split
            qm = (rng >= cb[i] - 100.0) & (rng <= cb[i] + 400.0) & np.isfinite(y[i])
            q = float(np.nansum(y[i][qm]) * dz)
            yi = np.interp(cb[i] + ZREL, rng, np.nan_to_num(y[i]),
                           left=np.nan, right=np.nan)
            ys.append(yi.astype("f4"))
            qs.append(q)
            cbhs.append(cb[i])
        n_days += 1
    return np.vstack(ys), np.array(qs), np.array(cbhs), n_days


def smooth(x, k=5):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    from scipy.optimize import curve_fit
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.4), sharey=True)
    scales = {"A": (1e3, "(×1e-3)"), "B": (1e5, "(×1e-5)"), "C": (1e15, "(×1e-15)")}
    for ax, ident, itype in zip(axes, "ABC", ("CHM15k", "CL31", "CL61")):
        Y, Q, CB, n_days = collect(ident)
        t1, t2 = np.nanpercentile(Q, [33.3, 66.7])
        dim = np.nanmedian(Y[Q <= t1], axis=0)
        bright = np.nanmedian(Y[Q >= t2], axis=0)
        diff = bright - dim
        lever = float(np.nanmedian(Q[Q >= t2]) / np.nanmedian(Q[Q <= t1]))

        # hood reference shifted into the relative frame by the median CBH
        trng, b_rcs, _sem, b_p = hood.truth(ident)
        if ident == "B":
            sess = [s for s in hood.session_frames("B", min_hours=2.0)
                    if s["era"] == "pre_swap"]
            if sess:
                b_rcs = np.nanmedian(np.vstack(
                    [np.nanmedian(s["dark"], axis=0) for s in sess]), axis=0)
                trng = sess[0]["rng"]
        med_cb = float(np.nanmedian(CB))
        zt = np.where(trng > 0, trng, np.nan)
        hood_rel = np.interp(med_cb + ZREL, trng, b_rcs / zt ** 2,
                             left=np.nan, right=np.nan)

        sc0 = scales[ident][0]
        ax.plot(smooth(dim) * sc0, ZREL / 1e3, "-", color="steelblue", lw=1.4,
                label=f"dim tercile (n={int((Q <= t1).sum())})")
        ax.plot(smooth(bright) * sc0, ZREL / 1e3, "-", color="darkorange", lw=1.4,
                label=f"bright tercile (n={int((Q >= t2).sum())}, Q ×{lever:.1f})")
        ax.plot(smooth(diff) * sc0, ZREL / 1e3, "-", color="crimson", lw=2.0,
                label="bright − dim = pulse-response tail")
        ax.plot(smooth(hood_rel) * sc0, ZREL / 1e3, "k--", lw=1.6,
                label="hood dark (shifted by median CBH)")
        ax.axvline(0, color="k", lw=0.5)

        # quantify: impact size vs the hood-dark scale in the same window
        jm = (ZREL >= 500) & (ZREL <= 3000)
        imp = float(np.sqrt(np.nanmean(diff[jm] ** 2)))
        hs = float(np.sqrt(np.nanmean(hood_rel[jm] ** 2)))
        # detector constant from the sky: decay length (C) / ripple wavelength (B)
        note = ""
        d_s = smooth(diff)
        try:
            if ident in "AC":
                sgn = -1.0 if np.nanmedian(d_s[jm]) < 0 else 1.0
                m = jm & (sgn * d_s > 0)
                p = np.polyfit(ZREL[m], np.log(sgn * d_s[m]), 1)
                ref_note = " (hood slow pole: 4570 m)" if ident == "C" else                     " (CHM15k: no single hood pole)"
                note = f"decay L = {-1 / p[0]:.0f} m{ref_note}"
            else:
                fm = (ZREL >= 200) & (ZREL <= 3000)

                def ring(z, A, L, lam, ph):
                    return A * np.exp(-z / L) * np.cos(2 * np.pi * z / lam + ph)
                p, _ = curve_fit(ring, ZREL[fm], d_s[fm],
                                 p0=[np.nanmax(np.abs(d_s[fm])), 1500.0, 5080.0, 0.0],
                                 maxfev=20000)
                note = f"ripple Λ = {abs(p[2]):.0f} m (hood slow mode: 5080 m)"
        except Exception:                                                   # noqa: BLE001
            note = "constant fit failed"
        ax.set_title(f"{itype} (Payerne {ident}) — {Y.shape[0]} profiles, {n_days} days\n"
                     f"impact RMS(0.5-3 km) = {imp / hs:.2f} × hood   |   {note}",
                     fontsize=9)
        lin = float(np.sqrt(np.nanmean(diff[jm] ** 2))
                    / max(np.sqrt(np.nanmean((dim[jm] - hood_rel[jm]) ** 2)), 1e-300))
        ax.set_xlabel("offset, P-view " + scales[ident][1])
        # x-limits from the ABOVE-cloud structure, robust (p95, z' >= 1 km) so neither the
        # cloud peak nor the CHM15k first-km saturation spike sets the scale — they clip,
        # and the off-scale factor is annotated instead
        vm = (ZREL >= 1000)
        ref = np.concatenate([smooth(diff)[vm], smooth(hood_rel)[vm]]) * sc0
        span = float(np.nanpercentile(np.abs(ref[np.isfinite(ref)]), 95))
        ax.set_xlim(-3.0 * span, 3.0 * span)
        pk = float(np.nanmax(np.abs(smooth(diff)[(ZREL >= 100) & (ZREL < 1000)])) * sc0)
        if pk > 3.0 * span:
            ax.text(0.03, 0.03, f"0–1 km spike off-scale (×{pk / span:.0f})",
                    transform=ax.transAxes, fontsize=8, color="crimson")
        if ident == "A":
            ax.set_ylabel("height above cloud base (km)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="best")
        print(f"{itype}: {Y.shape[0]} prof / {n_days} d, Q lever ×{lever:.1f}, "
              f"impact/hood = {imp / hs:.2f}, charge-linearity = {lin:.2f}, {note}")
    fig.suptitle("The cloud's fingerprint on the profile above it — brightness terciles in the "
                 "cloud-relative frame: (bright − dim) isolates the signal-induced detector "
                 "response, and its shape measures the circuit constants from the sky",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig18_cloud_impact.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
