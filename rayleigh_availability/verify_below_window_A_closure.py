# -*- coding: utf-8 -*-
"""ADVERSARIAL CHECK A -- closure of the shipped chain + an INDEPENDENT forward/inverse model.

Two things are tested here, both prerequisites for believing anything downstream:

 1. CLOSURE. Run the repo's own ``run_closure_tests`` (pure molecular, noiseless) and confirm the
    shipped retrieval returns C_true. The claim under review reported "CLOSURE: not reported"; the
    stored JSON of the analyst's run DOES contain it, so this re-run checks the number.

 2. INDEPENDENCE. The whole scan under review rests on ONE quantity: the error the retrieval makes
    on the two-way transmission of an aerosol layer that sits BELOW the molecular fit window, when
    the assumed lidar ratio (52 sr) differs from the true one. I re-derive that quantity with my
    OWN forward model and my OWN Fernald/Klett backward inversion, sharing nothing with the repo
    except the molecular profile, and compare.

FORWARD EQUATION (mine, written out; SI throughout)
    beta_tot(z) = beta_mol(z) + beta_aer(z)                [m^-1 sr^-1]
    ext_tot(z)  = S_mol*beta_mol(z) + S_true*beta_aer(z)   [m^-1],  S_mol = 8*pi/3 sr
    tau(z)      = INT_0^z ext_tot dz'                      [-]
    rcs(z)      = C_true * beta_tot(z) * exp(-2*tau(z))    [instrument units]

INVERSE (mine): Fernald backward inversion from a reference gate in the window where the aerosol
is assumed zero, with the ASSUMED ratio S_a = 52 sr, then
    C_L(z) = rcs(z) / beta_tot_ret(z) * exp(+2*tau_ret(z))
and C_L is the median over the window -- the same estimator as
``rayleigh_fit.calculate_lidar_constant``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

from forward_model import (  # noqa: E402
    PAYERNE_CHM15K, default_grid, default_options, forward_signal, make_atmosphere,
    retrieve, run_closure_tests,
)
from calibration.rayleigh.atmosphere import MOLECULAR_LIDAR_RATIO  # noqa: E402

OUT = Path(__file__).parent / "verif_A.json"


# ------------------------------------------------------------------------------------------------
# my own inversion -- shares nothing with calibration/rayleigh
# ------------------------------------------------------------------------------------------------
def my_fernald_cl(z, rcs, beta_mol, s_assumed, i_ref, i_lo, i_hi):
    """Independent Fernald(1984) backward inversion + lidar constant over [i_lo, i_hi].

    Backward (downward) recursion from ``i_ref`` where beta_aer is assumed 0, in the standard
    two-component form.  Working variable X(z) = rcs(z) (already range-corrected).

    Discrete backward step (Fernald 1984 eq. 4, trapezoid form), from level k+1 down to k::

        A       = (S_a - S_m) * (beta_mol[k] + beta_mol[k+1]) * dz
        beta[k] = X[k]*exp(A) / ( X[k+1]/beta[k+1] + S_a*(X[k]*exp(A) + X[k+1])*dz )

    with beta = beta_tot.  Then ext_aer = S_a*(beta_tot - beta_mol), and the constant is
    median over the window of rcs/beta_tot*exp(2*tau_ret), tau_ret integrated from the ground
    with the trapezoid rule (leading rectangle included, matching the forward convention).
    """
    n = z.size
    dz = float(z[1] - z[0])
    s_m = float(MOLECULAR_LIDAR_RATIO)
    beta_tot = np.full(n, np.nan)
    beta_tot[i_ref] = beta_mol[i_ref]                      # reference: pure molecular
    for k in range(i_ref - 1, -1, -1):
        A = (s_assumed - s_m) * (beta_mol[k] + beta_mol[k + 1]) * dz
        num = rcs[k] * np.exp(A)
        den = rcs[k + 1] / beta_tot[k + 1] + s_assumed * (rcs[k] * np.exp(A) + rcs[k + 1]) * dz
        beta_tot[k] = num / den
    # above the reference we simply keep it molecular (production does the same: the Klett is
    # run over [i_start, reference] and the profile above is treated as molecular)
    beta_tot[i_ref:] = beta_mol[i_ref:]

    ext_aer = np.maximum(beta_tot - beta_mol, 0.0) * s_assumed
    ext_tot = beta_mol * s_m + ext_aer
    tau = np.empty(n)
    tau[0] = ext_tot[0] * z[0]
    tau[1:] = tau[0] + np.cumsum(0.5 * np.diff(z) * (ext_tot[:-1] + ext_tot[1:]))
    cl = rcs / beta_tot * np.exp(2.0 * tau)
    return float(np.median(cl[i_lo:i_hi + 1])), beta_tot, ext_aer, tau


def my_forward(z, C_true, beta_mol, beta_aer, s_true):
    s_m = float(MOLECULAR_LIDAR_RATIO)
    ext_tot = beta_mol * s_m + beta_aer * s_true
    tau = np.empty(z.size)
    tau[0] = ext_tot[0] * z[0]
    tau[1:] = tau[0] + np.cumsum(0.5 * np.diff(z) * (ext_tot[:-1] + ext_tot[1:]))
    return C_true * (beta_mol + beta_aer) * np.exp(-2.0 * tau)


def main():
    res = {}

    # ---- 1. closure of the SHIPPED chain -------------------------------------------------------
    clos = run_closure_tests(fig_path=None, n_seeds=50, verbose=False)
    res["closure_shipped"] = dict(
        max_abs_dev_pct=clos["clean_max_abs_dev_pct"],
        slope_pct_per_km=clos["clean_slope_pct_per_km"],
        linearity_max_abs_err=clos["linearity_max_abs_err"],
        dev_by_centre_pct=dict(zip([str(int(c)) for c in clos["centres_m"]],
                                   clos["clean_dev_pct"])),
    )
    print(f"[A1] SHIPPED closure: max|dev| = {clos['clean_max_abs_dev_pct']:.5f} % , "
          f"slope = {clos['clean_slope_pct_per_km']:+.6f} %/km , "
          f"linearity err = {clos['linearity_max_abs_err']:.1e}")

    # ---- 2. closure of MY chain (must also return C_true on pure molecular) --------------------
    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    beta_mol = atm["beta_mol"]
    zeros = np.zeros_like(z)

    half = 490.0
    centres = np.arange(2500.0, 6501.0, 500.0)
    my_clean = []
    for c in centres:
        i_lo = int(np.argmin(np.abs(z - (c - half))))
        i_hi = int(np.argmin(np.abs(z - (c + half))))
        i_ref = i_hi
        rcs = my_forward(z, C_true, beta_mol, zeros, 52.0)
        cl, *_ = my_fernald_cl(z, rcs, beta_mol, 52.0, i_ref, i_lo, i_hi)
        my_clean.append((cl / C_true - 1.0) * 100.0)
    res["closure_mine"] = dict(centres_m=centres.tolist(), dev_pct=my_clean,
                               max_abs_dev_pct=float(np.max(np.abs(my_clean))))
    print(f"[A2] MY closure:      max|dev| = {np.max(np.abs(my_clean)):.5f} %")

    # ---- 3. head-to-head on the below-window layer ---------------------------------------------
    # Scene under review: aerosol confined to 100-1500 m AGL, window centred at 3.0 and 6.0 km.
    from forward_scan_transmission import scene  # noqa: E402
    options = default_options()
    rows = []
    for aod in (0.05, 0.20, 0.35, 0.50):
        for s_true in (20.0, 30.0, 52.0, 80.0, 90.0):
            ba, ea = scene(z, "below", aod=aod, lidar_ratio_true=s_true, top_m=1500.0)
            aod_check = float(np.trapezoid(ea, z))
            # --- shipped
            rcs_ship = forward_signal(z, C_true, beta_mol, ba, ea)
            r3 = retrieve(z, rcs_ship, atm, options=options, window=(3000. - half, 3000. + half))
            r6 = retrieve(z, rcs_ship, atm, options=options, window=(6000. - half, 6000. + half))
            ship3 = (r3["C_L"] / C_true - 1.0) * 100.0 if r3["ok"] else np.nan
            ship6 = (r6["C_L"] / C_true - 1.0) * 100.0 if r6["ok"] else np.nan
            # --- mine (independent forward AND inverse)
            rcs_mine = my_forward(z, C_true, beta_mol, ba, s_true)
            fwd_rel = float(np.max(np.abs(rcs_mine / rcs_ship - 1.0)))
            mine = {}
            for c in (3000.0, 6000.0):
                i_lo = int(np.argmin(np.abs(z - (c - half))))
                i_hi = int(np.argmin(np.abs(z - (c + half))))
                cl, *_ = my_fernald_cl(z, rcs_mine, beta_mol, 52.0, i_hi, i_lo, i_hi)
                mine[c] = (cl / C_true - 1.0) * 100.0
            # --- pure-analytic expectation: window above the layer, so beta_tot is molecular
            #     there and the ONLY error is the optical depth of the layer:
            #     C_L/C_true = exp(2*(tau_ret - tau_true)) with tau_ret ~ AOD*S_a/S_true
            analytic = (np.exp(2.0 * aod * (52.0 / s_true - 1.0)) - 1.0) * 100.0
            rows.append(dict(aod=aod, aod_check=aod_check, s_true=s_true,
                             shipped_3km_pct=ship3, shipped_6km_pct=ship6,
                             mine_3km_pct=mine[3000.0], mine_6km_pct=mine[6000.0],
                             fwd_max_rel_diff=fwd_rel,
                             analytic_linear_pct=analytic))
            print(f"  AOD={aod:.2f} S_true={s_true:4.0f} sr | shipped 3km {ship3:+8.3f} % "
                  f"6km {ship6:+8.3f} % | mine 3km {mine[3000.]:+8.3f} % 6km {mine[6000.]:+8.3f} % "
                  f"| linear-Klett analytic {analytic:+8.2f} % | fwd diff {fwd_rel:.2e}")
    res["head_to_head_below"] = rows

    d = np.array([r["mine_3km_pct"] - r["shipped_3km_pct"] for r in rows])
    print(f"[A3] independent-vs-shipped max |difference| at 3 km = {np.nanmax(np.abs(d)):.3f} "
          f"percentage points (n={len(rows)})")
    res["indep_vs_shipped_max_abs_pp"] = float(np.nanmax(np.abs(d)))

    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
