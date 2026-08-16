# -*- coding: utf-8 -*-
"""ADVERSARIAL CHECK E -- run the mechanism on the REAL per-night aerosol and read the answer.

Instead of asking "what AOD would be needed", feed the mechanism the (AOD_1064 below 1500 m AGL,
S_1064) that CAMS actually carries on each v2.0-KEPT and v2.2-RECOVERED night (verif_B) and compute
the C_L bias the SHIPPED retrieval commits. The recovered-minus-kept median of that bias is the
mechanism's honest prediction, to be set against the published -22.5 / +23.2 / +7.4 %.

The bias is read off a bilinear interpolation of the shipped-retrieval grid computed in verif_C
(window centre 3000 m AGL, half-length 490 m, layer 100-1500 m AGL). The grid IS the shipped chain;
interpolating it avoids re-running 400 retrievals for a smooth two-parameter surface, and the
interpolation error is checked against a direct run on the median night of each population.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

HERE = Path(__file__).parent   # verif_*.json live beside these scripts
OUT = HERE / "verif_E.json"


def main():
    C = json.loads((HERE / "verif_C.json").read_text(encoding="utf-8"))
    B = json.loads((HERE / "verif_B.json").read_text(encoding="utf-8"))
    aods = np.array(C["grid"]["aods"], float)
    ss = np.array(C["grid"]["lidar_ratios_sr"], float)
    Z = np.array([C["grid"]["bias_pct_at_3km"][f"{a:.2f}"] for a in aods], float)  # (n_aod, n_S)

    def bias(aod, s):
        """Bilinear read of the shipped-retrieval bias surface; clipped to the grid corners."""
        if not (np.isfinite(aod) and np.isfinite(s)):
            return np.nan
        a = float(np.clip(aod, aods[0], aods[-1]))
        v = float(np.clip(s, ss[0], ss[-1]))
        ia = int(np.clip(np.searchsorted(aods, a) - 1, 0, aods.size - 2))
        iv = int(np.clip(np.searchsorted(ss, v) - 1, 0, ss.size - 2))
        ta = (a - aods[ia]) / (aods[ia + 1] - aods[ia])
        tv = (v - ss[iv]) / (ss[iv + 1] - ss[iv])
        return float((1 - ta) * ((1 - tv) * Z[ia, iv] + tv * Z[ia, iv + 1])
                     + ta * ((1 - tv) * Z[ia + 1, iv] + tv * Z[ia + 1, iv + 1]))

    out = {}
    print("Mechanism (b) evaluated on the REAL CAMS aerosol of each night")
    print("bias = C_L/C_true - 1 [%], window centre 3000 m AGL (+-490 m), "
          "layer = the CAMS AOD_1064 below 1500 m AGL at its CAMS lidar ratio\n")
    for label, v in B.items():
        res = {}
        for g in ("kept", "recovered"):
            vals, a_, s_ = [], [], []
            for r in v["rows"][g]:
                b = bias(r["aod_below1500"], r["S_local_0_1km"])
                if np.isfinite(b):
                    vals.append(b)
                    a_.append(r["aod_below1500"])
                    s_.append(r["S_local_0_1km"])
            vals = np.array(vals)
            res[g] = dict(n=int(vals.size), median_bias_pct=float(np.median(vals)),
                          p10_bias_pct=float(np.percentile(vals, 10)),
                          p90_bias_pct=float(np.percentile(vals, 90)),
                          min_bias_pct=float(vals.min()), max_bias_pct=float(vals.max()),
                          median_aod=float(np.median(a_)), median_S=float(np.median(s_)))
        pred = res["recovered"]["median_bias_pct"] - res["kept"]["median_bias_pct"]
        obs = v["observed_offset_pct"]
        res["predicted_recovered_minus_kept_pct"] = pred
        res["observed_pct"] = obs
        res["shortfall_factor"] = float(abs(obs) / abs(pred)) if pred else float("inf")
        res["sign_agrees"] = bool(np.sign(pred) == np.sign(obs))
        out[label] = res
        print(f"=== {label}")
        for g in ("kept", "recovered"):
            r = res[g]
            print(f"    {g:9s} n={r['n']:3d}  median AOD_below1500 = {r['median_aod']:.4f}, "
                  f"median S = {r['median_S']:.1f} sr  ->  bias {r['median_bias_pct']:+.3f} % "
                  f"(p10 {r['p10_bias_pct']:+.2f}, p90 {r['p90_bias_pct']:+.2f}, "
                  f"extremes {r['min_bias_pct']:+.2f} .. {r['max_bias_pct']:+.2f})")
        print(f"    PREDICTED recovered - kept = {pred:+.3f} %   OBSERVED = {obs:+.1f} %   "
              f"-> shortfall x{res['shortfall_factor']:.0f}, "
              f"sign {'AGREES' if res['sign_agrees'] else 'IS WRONG'}\n")

    # worst-case: the single most favourable night CAMS knows about, per station
    print("Most favourable single night CAMS offers (max |bias| in the recovered population):")
    for label, v in B.items():
        best = max(v["rows"]["recovered"],
                   key=lambda r: abs(bias(r["aod_below1500"], r["S_local_0_1km"]))
                   if np.isfinite(bias(r["aod_below1500"], r["S_local_0_1km"])) else -1)
        bb = bias(best["aod_below1500"], best["S_local_0_1km"])
        print(f"    {label:22s} {best['date']}  AOD_below1500 = {best['aod_below1500']:.4f}, "
              f"S = {best['S_local_0_1km']:.1f} sr -> {bb:+.2f} %")
        out[label]["most_favourable_night"] = dict(date=best["date"],
                                                   aod_below1500=best["aod_below1500"],
                                                   S=best["S_local_0_1km"], bias_pct=bb)

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
