# -*- coding: utf-8 -*-
"""ADVERSARIAL CHECK C -- refine the required AOD, test the interaction, test the sign story.

Three attacks on the claim "AOD_1064 = 0.38 at S_true = 90 sr is needed for -22.5 %":

 C1. REFINED SOLVE. The claim interpolated a coarse (AOD, S) grid. Re-solve on a dense grid, for
     the target of EVERY site with a published recovered-minus-kept offset, and report the
     (AOD_1064, S_true) each site would need. Sites with a POSITIVE offset need S_true < 52 sr;
     that requirement is part of the verdict and has to be quantified, not asserted.

 C2. THE ZERO-MISMATCH FLOOR. At S_true = 52 sr the mechanism should be exactly zero. The shipped
     chain returns -0.65 % at AOD 0.35 and -1.31 % at AOD 0.50 while an independent Fernald
     returns -0.008 %. That residual is an ESTIMATOR artefact of the shipped Klett, not the
     mechanism; it is quantified here because it sets a floor the claim never mentions.

 C3. INTERACTION. The claim was measured with the below-window layer ALONE. The real recovered
     nights also carry (i) aerosol INSIDE the 2-6 km window -- that is what the measured
     signal/p_mol slope of -19.5 %/km says -- and (ii) an additive residual b (the free-intercept
     test measured about -6.5 % at Payerne). If the below-window term is not additive-in-log with
     those, the claim's amplitude does not transfer to a real night.

Everything is evaluated with the SHIPPED retrieval on forced windows; bias is always quoted as
C_L/C_true - 1 in %, at named window centres (m AGL) with the production half-length 490 m.
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
    PAYERNE_CHM15K, default_grid, default_options, forward_signal, make_atmosphere, retrieve,
)
from forward_scan_transmission import scene  # noqa: E402

OUT = Path(__file__).parent / "verif_C.json"
HALF = 490.0
# Published recovered-minus-kept medians (task brief).
TARGETS = {"Payerne CHM15k": -22.5, "Payerne CL61": -17.5, "Lindenberg CHM15k": +23.2,
           "Aosta CHM15k": +31.7, "Aosta CL61": +32.1, "SIRTA CHM15k": +7.4}


def bias_below(z, atm, C_true, options, aod, s_true, centre=3000.0,
               extra_haze=None, additive=0.0):
    """C_L/C_true - 1 [%] for a 100-1500 m AGL layer of optical depth ``aod`` and ratio ``s_true``.

    ``extra_haze`` = (beta, ext) of a second, in-window aerosol field added on top;
    ``additive`` = constant added to rcs (the free-intercept residual), in rcs units.
    """
    ba, ea = scene(z, "below", aod=aod, lidar_ratio_true=s_true, top_m=1500.0)
    if extra_haze is not None:
        ba = ba + extra_haze[0]
        ea = ea + extra_haze[1]
    rcs = forward_signal(z, C_true, atm["beta_mol"], ba, ea, additive_residual=additive)
    r = retrieve(z, rcs, atm, options=options, window=(centre - HALF, centre + HALF))
    return ((r["C_L"] / C_true - 1.0) * 100.0) if r["ok"] else np.nan


def main():
    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    options = default_options()
    out = {}

    # ---- C1: dense (AOD, S) map, then solve each site's target --------------------------------
    aods = [0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20, 0.28, 0.35, 0.45, 0.60, 0.80]
    ss = [15.0, 20.0, 25.0, 30.0, 40.0, 52.0, 65.0, 80.0, 90.0, 110.0]
    grid = {}
    print("C1  bias [%] at window centre 3000 m AGL, layer 100-1500 m AGL")
    print("       AOD |" + "".join(f"{s:8.0f}" for s in ss))
    for aod in aods:
        row = [bias_below(z, atm, C_true, options, aod, s) for s in ss]
        grid[f"{aod:.2f}"] = row
        print(f"    {aod:6.2f} |" + "".join(f"{v:8.2f}" for v in row))
    out["grid"] = dict(aods=aods, lidar_ratios_sr=ss, bias_pct_at_3km=grid)

    def required_aod(target_pct, s_true):
        """Smallest AOD whose bias reaches ``target_pct`` at ``s_true`` (linear interp), or None."""
        col = np.array([grid[f'{a:.2f}'][ss.index(s_true)] for a in aods])
        a = np.array(aods)
        if target_pct < 0:
            ok = col <= target_pct
        else:
            ok = col >= target_pct
        if not ok.any():
            return None
        k = int(np.argmax(ok))
        if k == 0:
            return float(a[0])
        x0, x1, y0, y1 = a[k - 1], a[k], col[k - 1], col[k]
        return float(x0 + (target_pct - y0) * (x1 - x0) / (y1 - y0))

    req = {}
    print("\nC1  AOD_1064 required to reach each published offset")
    for site, tgt in TARGETS.items():
        cand = {}
        for s in ss:
            if s == 52.0:
                continue
            if (tgt < 0) == (s > 52.0):
                r = required_aod(tgt, s)
                if r is not None:
                    cand[s] = r
        req[site] = dict(target_pct=tgt, required_aod_by_S=cand,
                         min_required_aod=(min(cand.values()) if cand else None),
                         S_at_min=(min(cand, key=cand.get) if cand else None))
        if cand:
            s_best = min(cand, key=cand.get)
            print(f"    {site:20s} target {tgt:+6.1f} %  -> cheapest is S_true = {s_best:5.0f} sr "
                  f"with AOD_1064 = {cand[s_best]:.3f}   "
                  f"(needs S_true {'>' if tgt < 0 else '<'} 52 sr)")
        else:
            print(f"    {site:20s} target {tgt:+6.1f} %  -> UNREACHABLE on this grid")
    out["required"] = req

    # ---- C2: the zero-mismatch floor of the shipped chain --------------------------------------
    floor = {f"{a:.2f}": bias_below(z, atm, C_true, options, a, 52.0) for a in aods}
    print("\nC2  shipped-chain bias with NO lidar-ratio error (S_true = S_assumed = 52 sr):")
    print("    " + "  ".join(f"AOD {k}: {v:+.3f} %" for k, v in list(floor.items())[::3]))
    out["zero_mismatch_floor_pct"] = floor

    # ---- C3: interaction with in-window haze and with the additive residual --------------------
    from forward_model import exponential_aerosol
    # In-window haze tuned to the measured signal/p_mol slope regime; S fixed at 52 so the haze
    # contributes NO lidar-ratio error of its own and any change is pure interaction.
    haze = exponential_aerosol(z, 3.0e-8, 3000.0, 52.0)
    b_amp = PAYERNE_CHM15K["additive_residual_signal"] * (3000.0 ** 2)   # rcs units at 3 km scale
    inter = []
    print("\nC3  interaction (window centre 3000 m AGL)")
    print(f"{'AOD':>6} {'S':>5} {'alone':>9} {'+haze':>9} {'+b':>9} {'+both':>9} "
          f"{'sum-of-parts':>13} {'non-add':>8}")
    for aod in (0.05, 0.20, 0.38):
        for s_true in (30.0, 90.0):
            a0 = bias_below(z, atm, C_true, options, aod, s_true)
            ah = bias_below(z, atm, C_true, options, aod, s_true, extra_haze=haze)
            ab = bias_below(z, atm, C_true, options, aod, s_true, additive=b_amp)
            abo = bias_below(z, atm, C_true, options, aod, s_true, extra_haze=haze,
                             additive=b_amp)
            h_only = bias_below(z, atm, C_true, options, 1e-9, 52.0, extra_haze=haze)
            b_only = bias_below(z, atm, C_true, options, 1e-9, 52.0, additive=b_amp)
            # log-additive prediction (biases are multiplicative transmission factors)
            pred = (np.exp(np.log1p(a0 / 100) + np.log1p(h_only / 100)
                           + np.log1p(b_only / 100)) - 1.0) * 100.0
            inter.append(dict(aod=aod, s_true=s_true, alone=a0, with_haze=ah, with_b=ab,
                              with_both=abo, haze_only=h_only, b_only=b_only,
                              log_additive_prediction=pred, non_additivity_pp=abo - pred))
            print(f"{aod:6.2f} {s_true:5.0f} {a0:9.2f} {ah:9.2f} {ab:9.2f} {abo:9.2f} "
                  f"{pred:13.2f} {abo - pred:8.2f}")
    out["interaction"] = inter
    out["interaction_note"] = dict(haze_only_pct=inter[0]["haze_only"],
                                   b_only_pct=inter[0]["b_only"],
                                   b_amplitude_rcs=b_amp)

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
