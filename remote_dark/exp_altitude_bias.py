# -*- coding: utf-8 -*-
"""Fig 15 — operator's proposal: fit the Rayleigh constant at SEVERAL altitudes and use the
altitude dependence to estimate (and remove) the calibration bias, instead of trying to
retrieve the dark profile itself.

The reasoning being tested: every additive contaminant (electronic dark AND residual aerosol)
decays away with altitude, while the true constant does not. So the fitted constant vs window
altitude, C-hat(zc), should tilt at low windows and flatten toward the unbiased value aloft —
and the high-altitude asymptote is a bias-corrected constant that needs NO dark knowledge, NO
hood, NO aerosol model.

The referee is the Payerne hood: compute the same window scan on the dark-SUBTRACTED signal
(S - b_hood). If the raw scan and the dark-corrected scan CONVERGE at high windows, the
altitude route removes the dark automatically — the operator's idea works for the dark part.
Whatever common level both converge to below zero is the residual-aerosol share of the
reference window, which the altitude route also sheds but a hood never could.

Payerne A (CHM15k, the big-dark case) and Payerne C (CL61, small-dark control). All curves are
shown relative to C_darkcorr = the dark-corrected constant at the operational-like reference
window (4.5-6.5 km), so "0" means "what a hood correction would give".
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import V4_DIR, ensure_out
from remote_dark import hood

FIG = ensure_out("figs")
WIN_M = 1000.0
CENTRES = np.arange(2500.0, 9501.0, 250.0)
REF_BAND = (4500.0, 6500.0)
ASYM_BAND = (8000.0, 9500.0)

UNITS = [("A", "CHM15k", "0-20000-0-06610_A_20250402_20260813_220.npz"),
         ("C", "CL61", "0-20000-0-06610_C_20250101_20260813_220.npz")]


def window_fit(S, sig, M, rng, zc):
    m = (rng >= zc - WIN_M / 2) & (rng <= zc + WIN_M / 2)
    if m.sum() < 10:
        return np.full(S.shape[0], np.nan)
    w = 1.0 / np.maximum(sig[:, m], 1e-300) ** 2
    ok = np.isfinite(S[:, m]) & np.isfinite(M[:, m])
    w = np.where(ok, w, 0.0)
    S0, M0 = np.nan_to_num(S[:, m]), np.nan_to_num(M[:, m])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(w * S0 * M0, axis=1) / np.nansum(w * M0 * M0, axis=1)


def scan(S, sig, M, rng):
    """median over nights of C-hat(zc), plus the median ref-window constant."""
    C = np.column_stack([window_fit(S, sig, M, rng, zc) for zc in CENTRES])
    refm = (CENTRES >= REF_BAND[0]) & (CENTRES <= REF_BAND[1])
    c_ref_n = np.nanmedian(C[:, refm], axis=1)
    return np.nanmedian(C, axis=0), np.nanpercentile(C, [25, 75], axis=0), \
        float(np.nanmedian(c_ref_n))


def main():
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4), sharey=True)
    for ax, (ident, itype, cache_name) in zip(axes, UNITS):
        z = np.load(V4_DIR / "cache" / cache_name, allow_pickle=True)
        S, sig, M, rng = z["S"], z["sigma"], z["M"], z["rng"]
        trng, b_rcs, _sem, _bp = hood.truth(ident)
        b = np.interp(rng, trng, b_rcs, left=np.nan, right=np.nan)
        okb = np.isfinite(b)
        Sc = S.copy()
        Sc[:, okb] = Sc[:, okb] - b[None, okb]
        Sc[:, ~okb] = np.nan

        med_raw, (lo, hi), cref_raw = scan(S, sig, M, rng)
        med_cor, _, cref_cor = scan(Sc, sig, M, rng)

        # everything relative to the dark-corrected reference-window constant
        x_raw = (med_raw / cref_cor - 1.0) * 100.0
        x_cor = (med_cor / cref_cor - 1.0) * 100.0
        am = (CENTRES >= ASYM_BAND[0]) & (CENTRES <= ASYM_BAND[1])
        c_asym = float(np.nanmedian(med_raw[am]))
        bias_ref = (cref_raw / cref_cor - 1.0) * 100.0
        bias_asym = (c_asym / cref_cor - 1.0) * 100.0

        ax.fill_betweenx(CENTRES / 1e3, (lo / cref_cor - 1) * 100, (hi / cref_cor - 1) * 100,
                         color="0.88", label="p25-p75 nights (raw)")
        ax.plot(x_raw, CENTRES / 1e3, "k-", lw=1.8, label="raw scan  Ĉ(z_win)")
        ax.plot(x_cor, CENTRES / 1e3, "-", color="seagreen", lw=1.8,
                label="scan after hood-dark subtraction")
        ax.axvline(0, color="seagreen", lw=0.8, ls=":")
        ax.axvline(bias_ref, color="k", lw=0.8, ls=":")
        ax.axhspan(REF_BAND[0] / 1e3, REF_BAND[1] / 1e3, color="steelblue", alpha=0.10,
                   label="reference window 4.5-6.5 km")
        ax.axhspan(ASYM_BAND[0] / 1e3, ASYM_BAND[1] / 1e3, color="crimson", alpha=0.10,
                   label="asymptote band 8-9.5 km")
        ax.set_title(f"{itype} (Payerne {ident}) — {S.shape[0]} nights\n"
                     f"fixed-window bias (no dark corr): {bias_ref:+.1f} %   →   "
                     f"altitude-asymptote residual: {bias_asym:+.1f} %", fontsize=10)
        ax.set_xlabel("Ĉ / C_dark-corrected(4.5-6.5 km) − 1   (%)")
        ax.grid(alpha=0.25)
        if ident == "A":
            ax.set_ylabel("window-centre altitude (km)")
            ax.legend(fontsize=8, loc="upper left")
        ax.set_xlim(-30, 30)
        print(f"{itype}: C_ref(raw) vs dark-corrected: {bias_ref:+.1f} %   "
              f"asymptote 8-9.5 km vs dark-corrected: {bias_asym:+.1f} %   "
              f"(raw-corrected convergence aloft: "
              f"{np.nanmedian(np.abs(x_raw[am] - x_cor[am])):.1f} % apart)")
    fig.suptitle("Multi-altitude Rayleigh fit as a bias estimator (operator's idea) — does the "
                 "high-window asymptote reach the dark-corrected constant without knowing the "
                 "dark?  Black→green gap = the dark's share; both curves' common tilt = aerosol.",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    f = FIG / "fig15_altitude_bias.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
