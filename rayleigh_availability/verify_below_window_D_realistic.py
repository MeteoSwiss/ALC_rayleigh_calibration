# -*- coding: utf-8 -*-
"""ADVERSARIAL CHECK D -- the required AOD in a scene that also matches the OBSERVED R(z).

The claim isolates one term: a layer entirely BELOW the 2-6 km fit window, with a lidar-ratio
error. But the very nights it is trying to explain are not like that. The measured observable on
the v2.2-recovered Payerne nights is

    d ln(signal / p_mol) / dz = -19.5 %/km  between 2 and 6 km        (kept nights: -3.7 %/km)

which says there IS aerosol inside the window. So the honest question is not "how much
below-window AOD alone reproduces -22.5 %" but "how much below-window AOD reproduces -22.5 % in a
profile that ALSO reproduces -19.5 %/km".

Scene: exponential haze beta_aer(z) = beta_s * exp(-z/H), S = 52 sr (no ratio error of its own, so
it contributes no term-T of its own), with beta_s tuned so the simulated R(z) slope over 2-6 km
matches the target; PLUS the claim's 100-1500 m AGL layer of optical depth AOD and ratio S_true.
Bias is read at the window centre a recovered night actually uses (4500 m AGL) and, for
comparison, at 3000 m AGL, half-length 490 m in both cases.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

REPO = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

from forward_model import (  # noqa: E402
    PAYERNE_CHM15K, default_grid, default_options, exponential_aerosol, forward_signal,
    make_atmosphere, retrieve,
)
from forward_scan_transmission import scene, signal_ratio_slope  # noqa: E402

OUT = Path(__file__).parent / "verif_D.json"
HALF = 490.0
R_TARGET_RECOVERED = -19.5   # %/km, measured, 2-6 km
R_TARGET_KEPT = -3.7         # %/km, measured, 2-6 km


def build(z, atm, C_true, beta_s, H, aod_below, s_true):
    ba_h, ea_h = exponential_aerosol(z, beta_s, H, 52.0)
    ba_b, eb_b = scene(z, "below", aod=aod_below, lidar_ratio_true=s_true, top_m=1500.0)
    ba, ea = ba_h + ba_b, ea_h + eb_b
    return forward_signal(z, C_true, atm["beta_mol"], ba, ea), ba, ea


def main():
    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    options = default_options()
    out = {}

    # ---- 1. tune the in-window haze to the measured R(z) slope ---------------------------------
    H = 3000.0
    def r_slope(beta_s, aod_below=0.0, s_true=52.0):
        rcs, _, _ = build(z, atm, C_true, beta_s, H, aod_below, s_true)
        return signal_ratio_slope(z, rcs, atm)

    beta_s_rec = brentq(lambda b: r_slope(b) - R_TARGET_RECOVERED, 1e-12, 5e-6, xtol=1e-14)
    beta_s_kept = brentq(lambda b: r_slope(b) - R_TARGET_KEPT, 1e-14, 5e-6, xtol=1e-16)
    j45 = int(np.argmin(np.abs(z - 4500.0)))
    sr45 = 1.0 + exponential_aerosol(z, beta_s_rec, H, 52.0)[0][j45] / atm["beta_mol"][j45]
    print(f"D1  in-window haze, H = {H:.0f} m, S = 52 sr")
    print(f"    recovered-night target R' = {R_TARGET_RECOVERED:+.1f} %/km -> "
          f"beta_aer(0) = {beta_s_rec:.3e} m^-1 sr^-1, scattering ratio at 4.5 km = {sr45:.2f}, "
          f"column AOD_1064 of the haze = {np.trapezoid(exponential_aerosol(z, beta_s_rec, H, 52.0)[1], z):.4f}")
    print(f"    kept-night     target R' = {R_TARGET_KEPT:+.1f} %/km -> "
          f"beta_aer(0) = {beta_s_kept:.3e} m^-1 sr^-1")
    out["haze"] = dict(scale_height_m=H, beta_s_recovered=beta_s_rec, beta_s_kept=beta_s_kept,
                       sr_4500m_recovered=float(sr45),
                       haze_aod_recovered=float(np.trapezoid(
                           exponential_aerosol(z, beta_s_rec, H, 52.0)[1], z)))

    # ---- 2. bias vs below-window AOD, with and without that haze -------------------------------
    def bias(beta_s, aod_below, s_true, centre):
        rcs, _, _ = build(z, atm, C_true, beta_s, H, aod_below, s_true)
        r = retrieve(z, rcs, atm, options=options, window=(centre - HALF, centre + HALF))
        return ((r["C_L"] / C_true - 1.0) * 100.0) if r["ok"] else np.nan

    print("\nD2  C_L bias [%] vs below-window AOD_1064, S_true = 90 sr")
    print(f"{'AOD':>6} | {'no haze @3km':>13} {'no haze @4.5km':>15} | "
          f"{'+haze @3km':>11} {'+haze @4.5km':>13}")
    rows = []
    for aod in (0.0, 0.05, 0.10, 0.20, 0.38, 0.60, 0.90, 1.20):
        a3 = bias(0.0, aod, 90.0, 3000.0)
        a45 = bias(0.0, aod, 90.0, 4500.0)
        h3 = bias(beta_s_rec, aod, 90.0, 3000.0)
        h45 = bias(beta_s_rec, aod, 90.0, 4500.0)
        rows.append(dict(aod=aod, nohaze_3km=a3, nohaze_45km=a45, haze_3km=h3, haze_45km=h45))
        print(f"{aod:6.2f} | {a3:13.2f} {a45:15.2f} | {h3:11.2f} {h45:13.2f}")
    out["bias_vs_aod_S90"] = rows

    # ---- 3. what the observation actually asks for --------------------------------------------
    # The published -22.5 % is (median C on RECOVERED nights, fitted high) minus (median C on
    # KEPT nights, fitted low). So the simulated equivalent is
    #     bias(recovered scene, window 4500 m)  -  bias(kept scene, window 3000 m).
    kept_ref = bias(beta_s_kept, 0.0, 52.0, 3000.0)
    print(f"\nD3  reference: kept-night scene at 3000 m -> {kept_ref:+.2f} %")
    print(f"{'AOD':>6} {'recovered @4.5km':>17} {'recovered - kept':>18}")
    solve = []
    for aod in (0.0, 0.10, 0.20, 0.38, 0.60, 0.90, 1.20, 1.60):
        rec = bias(beta_s_rec, aod, 90.0, 4500.0)
        solve.append(dict(aod=aod, recovered_45km=rec, difference=rec - kept_ref))
        print(f"{aod:6.2f} {rec:17.2f} {rec - kept_ref:18.2f}")
    out["kept_reference_pct"] = kept_ref
    out["recovered_minus_kept"] = solve

    d = np.array([r["difference"] for r in solve])
    a = np.array([r["aod"] for r in solve])
    if np.nanmin(d) <= -22.5:
        k = int(np.argmax(d <= -22.5))
        req = a[k] if k == 0 else float(a[k - 1] + (-22.5 - d[k - 1]) * (a[k] - a[k - 1])
                                        / (d[k] - d[k - 1]))
        print(f"\n    -> required below-window AOD_1064 in the REALISTIC scene: {req:.3f}")
        out["required_aod_realistic"] = req
    else:
        print(f"\n    -> -22.5 % is NOT reachable even at AOD_1064 = {a[-1]:.2f} "
              f"(best {np.nanmin(d):+.2f} %)")
        out["required_aod_realistic"] = None
        out["best_reachable_pct"] = float(np.nanmin(d))

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
