# -*- coding: utf-8 -*-
"""Downward extension of the sub-window range limited by an altitude floor.

Probing 1 km under the window measures the boundary-layer aerosol for LOW windows (the shipped
default therefore starts at the window bottom), but for a HIGH window the layer underneath is free
troposphere, where a strong gradient is a real warning (Payerne 2026-09-07). Variant: extend 1 km
below the window, but not below a floor z_floor (sub-windows centred under it are dropped).
Evaluated on the recorded sub-windows: time-series and twin pulls, unit ranking, Payerne 2026-09-07.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from u23_analysis import (OUT, BASE_CFG, final_frame, pull_stats, twin_pull_stats,  # noqa: E402
                          unit_table, select_subwindows, vertical_from_fit)


def main():
    from scipy.stats import spearmanr
    from calibration.rayleigh.uncertainty import fit_log_gradient
    pay = pickle.load(open(OUT / "u23_payerne.pkl", "rb"))
    rows = []
    for floor in (None, 2500.0, 3000.0, 3500.0, 4000.0):
        if floor is None:
            cfg = dict(BASE_CFG, below_m=0.0, above_m=1000.0, snr_min=3.0, est="abs", resid=True,
                       D_km=1.0)
        else:
            cfg = dict(BASE_CFG, below_m=1000.0, min_alt_m=floor, above_m=1000.0, snr_min=3.0,
                       est="abs", resid=True, D_km=1.0)
        df, tw, _ = final_frame(cfg)
        df["z22"] = df.anom / df.sig22
        m = df.ok & np.isfinite(df.anom) & (df.type == "CHM15k") & (df.key != "0-20000-0-06610_A")
        s = pull_stats(df.z23[m])
        t, _ = twin_pull_stats(tw, dict(zip(zip(df.key, df.date), df.sig23_nomol)))
        ut = unit_table(df)
        uc = ut[ut.type == "CHM15k"]
        rho = spearmanr(uc.rep23, uc.obs_sigma)[0]
        rank = uc.rep23.rank(ascending=False)
        ranks = {k: int(rank[uc.key == k].iloc[0]) for k in
                 ("0-20000-0-00203_A", "0-20000-0-06610_A", "0-20000-0-10564_0")}
        # Payerne 2026-09-07 from its recorded (wide) sub-windows
        rec = pay["20260907"][("record", True)]
        b = rec["budget"]
        zb = rec["zb"] - 490.0
        zt = rec["zt"] - 490.0
        z, r, sg = b.sub_centre_m, b.sub_cl / b.lidar_constant, b.sub_sigma_rel
        zz, rr, ss = select_subwindows(z, r, sg, zb, zt, cfg)
        f = fit_log_gradient(zz, rr, ss, 0.5 * (zb + zt), cfg["floor"])
        sv = vertical_from_fit([f["g"]], [f["g_sigma"]], [f["resid_rel"]], 1.0, "abs", True)[0]
        tot = 2 * np.sqrt(b.sigma_ensemble_rel ** 2 + sv ** 2 + b.sigma_molecular_rel ** 2
                          + 0.02 ** 2 + b.sigma_noise_rel ** 2)
        rows.append(dict(floor_m=floor if floor else "window bottom", ts_rstd=s["rstd"],
                         ts_out2=s["out2"], twin_rstd=t["rstd"], twin_out2=t["out2"],
                         spearman_chm=rho, sig_med=float(np.median(df.sig23[m])),
                         rank_messina=ranks["0-20000-0-00203_A"],
                         rank_payerne=ranks["0-20000-0-06610_A"],
                         rank_schleiz=ranks["0-20000-0-10564_0"],
                         payerne_0907_n=len(zz), payerne_0907_g=f["g"],
                         payerne_0907_u2s=tot))
    import pandas as pd
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "u23_floor_test.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
