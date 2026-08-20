# -*- coding: utf-8 -*-
"""M1-bis — the window scan (operator's idea, 2026-08-20).

Slide the Rayleigh fit window in altitude and watch the fitted constant move. Within ONE night
the true amplitude (C·T²) is common to every window, so the ratio

    r_n(zc) = Ĉ_n(window at zc) / Ĉ_n(reference window)

is SELF-REFERENCED: the night's unknown amplitude cancels, and what remains is the additive
baseline's window-integrated bias,

    Ĉ = A·(1 + <M·b>_w / (A·<M²>_w))     =>     r(zc) maps b's altitude structure,

measured per night with no cross-night amplitude assumption at all — a different projection of
the identifiability problem than the v4 per-gate line, and one the operator can read directly:
"which fit altitudes are biased by how much".

Validation: at Payerne the hood b predicts r(zc) exactly (same formula with b_hood inserted);
measured-vs-predicted over the scan is the test. At Schiphol the four co-located units give the
per-unit contrast on one sky.
"""
from __future__ import annotations

import numpy as np

from remote_dark.common import V4_DIR, ensure_out
from remote_dark import hood

WIN_M = 1000.0                     #: window length (the operational Rayleigh window scale)
CENTRES = np.arange(2500.0, 9501.0, 250.0)
REF_BAND = (4000.0, 6000.0)        #: reference = the mid scan, where the hood dark is moderate


def scan_stream(cache_file, ref_band=REF_BAND):
    z = np.load(cache_file, allow_pickle=True)
    S, sig, M, rng = z["S"], z["sigma"], z["M"], z["rng"]
    n_nights = S.shape[0]
    A_hat = np.full((n_nights, CENTRES.size), np.nan)
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    for j, zc in enumerate(CENTRES):
        m = (rng >= zc - WIN_M / 2) & (rng <= zc + WIN_M / 2)
        if m.sum() < 10:
            continue
        num = np.nansum(np.where(m, w_all * S * M, 0.0), axis=1)
        den = np.nansum(np.where(m, w_all * M * M, 0.0), axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            A_hat[:, j] = num / den
    refm = (CENTRES >= ref_band[0]) & (CENTRES <= ref_band[1])
    ref = np.nanmedian(A_hat[:, refm], axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = A_hat / ref[:, None]
    r[(ref <= 0) | ~np.isfinite(ref)] = np.nan
    return {"centres": CENTRES, "r_med": np.nanmedian(r, axis=0),
            "r_lo": np.nanpercentile(r, 25, axis=0), "r_hi": np.nanpercentile(r, 75, axis=0),
            "n_nights": n_nights, "rng": rng, "M": M, "sig": sig, "ref": ref}


def predict_from_hood(scanres, ident):
    """The same ratio with the hood b inserted: what the scan SHOULD show if the hood dark is
    the whole story. Per night (its own M and amplitude), then median — so the prediction
    carries the real nights' weighting, not an idealised atmosphere."""
    trng, t_brcs, _sem, _tp = hood.truth(ident)
    rng = scanres["rng"]
    b = np.interp(rng, trng, t_brcs, left=np.nan, right=np.nan)
    M, sig, ref = scanres["M"], scanres["sig"], scanres["ref"]
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    pred = np.full((M.shape[0], CENTRES.size), np.nan)
    for j, zc in enumerate(CENTRES):
        m = (rng >= zc - WIN_M / 2) & (rng <= zc + WIN_M / 2) & np.isfinite(b)
        if m.sum() < 10:
            continue
        num = np.nansum(np.where(m, w_all * M * b[None, :], 0.0), axis=1)
        den = np.nansum(np.where(m, w_all * M * M, 0.0), axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            pred[:, j] = 1.0 + (num / den) / ref
    refm = (CENTRES >= REF_BAND[0]) & (CENTRES <= REF_BAND[1])
    pr = pred / np.nanmedian(pred[:, refm], axis=1)[:, None]
    return np.nanmedian(pr, axis=0)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG = ensure_out("figs")

    # --- Payerne validation panels + Schiphol contrast panels, one simple row ------------------
    jobs = [("0-20000-0-06610_A_20250402_20260813_220.npz", "Payerne A (CHM15k)", "A"),
            ("0-20000-0-06610_C_20250101_20260813_220.npz", "Payerne C (CL61)", "C"),
            ("0-20000-0-06240_A_20250101_20260813_220.npz", "Schiphol A", None),
            ("0-20000-0-06240_B_20250101_20260813_220.npz", "Schiphol B", None),
            ("0-20000-0-06240_C_20250101_20260813_220.npz", "Schiphol C", None),
            ("0-20000-0-06240_D_20250101_20260813_220.npz", "Schiphol D", None)]
    fig, axes = plt.subplots(1, len(jobs), figsize=(15.5, 5.0), sharey=True)
    for ax, (name, lab, hood_ident) in zip(axes, jobs):
        res = scan_stream(V4_DIR / "cache" / name)
        x = (res["r_med"] - 1.0) * 100.0
        ax.fill_betweenx(CENTRES / 1e3, (res["r_lo"] - 1) * 100, (res["r_hi"] - 1) * 100,
                         color="0.85", label="p25–p75 (nights)")
        ax.plot(x, CENTRES / 1e3, "k-", lw=1.7, label="measured")
        if hood_ident:
            pr = predict_from_hood(res, hood_ident)
            ax.plot((pr - 1.0) * 100.0, CENTRES / 1e3, "r--", lw=1.6, label="hood-predicted")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlim(-25, 25)
        ax.set_title(f"{lab}\n{res['n_nights']} nights", fontsize=9)
        ax.set_xlabel("Ĉ(window) / Ĉ(4–6 km) − 1  (%)")
        ax.grid(alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel("window-centre altitude (km)")
            ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Window scan — slide the Rayleigh fit window and watch the constant move: "
                 "the shape IS the additive baseline's bias, self-referenced per night",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    f = FIG / "fig10_window_scan.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
