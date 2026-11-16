"""Step 4 complement: CBH flatness of the daily CL61 cloud C_L for each eta variant.

The eta variants differ less in level than in their CBH dependence (PVC a_G 3.3 / alpha 31 /km falls
from 0.88 at 0.5 km to 0.73 at 2 km; the operational 5.5 um / 10 /km table from 0.92 to 0.81). A
correct eta makes the calibration independent of the cloud height, so the slope of ln C_L against
the day's median CBH is a level-free test (no Rayleigh reference, no dark question).

Daily median CBH = median over the accepted profiles of that unit-day (incloud/<key>.npz, the
operational variant's scenes). Slope by Theil-Sen (robust), in %/km, with its 95 % interval;
also per-unit and pooled (unit-normalised: ln C_L minus the unit median).

Run:  python validation/eta_study/step4_cbh_slope.py
Output: C:/DATA/Projects/202611_eta_study/step4/cbh_slope_by_variant.csv (+ printout)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

BASE = Path("C:/DATA/Projects/202611_eta_study")
VARIANTS = ["a", "b", "c", "d", "e", "f_lo", "f_hi", "f_ag25", "f_ag50"]


def main():
    import pandas as pd
    from scipy.stats import theilslopes
    df = pd.read_csv(BASE / "step4" / "daily_cloud_variants.csv")
    cbh = []
    for k in df.key.unique():
        z = np.load(BASE / "incloud" / f"{k}.npz", allow_pickle=True)
        d = pd.DataFrame(dict(date=z["date"], cbh=z["cbh_m"].astype(float)))
        g = d.groupby("date").cbh.median().reset_index()
        g["key"] = k
        cbh.append(g)
    cb = pd.concat(cbh)
    df = df.merge(cb, on=["key", "date"], how="inner")
    df = df[np.isfinite(df.C_L) & (df.C_L > 0) & np.isfinite(df.cbh)]
    df["lnC"] = np.log(df.C_L)
    df["lnC_n"] = df.lnC - df.groupby(["key", "variant"]).lnC.transform("median")
    rows = []
    for var in VARIANTS:
        for scope, g in [("pooled 11 CL61", df[df.variant == var])] + \
                        [(k, df[(df.variant == var) & (df.key == k)]) for k in
                         ("0-20000-0-06610_C", "0-20000-0-10393_C")]:
            if len(g) < 15:
                continue
            y = g.lnC_n if scope.startswith("pooled") else g.lnC
            s, i, lo, hi = theilslopes(y, g.cbh / 1000.0)
            rows.append(dict(variant=var, scope=scope, n_days=len(g), slope_pct_per_km=100 * s,
                             ci_lo=100 * lo, ci_hi=100 * hi, cbh_med_km=float(g.cbh.median() / 1000)))
    t = pd.DataFrame(rows)
    t.to_csv(BASE / "step4" / "cbh_slope_by_variant.csv", index=False, float_format="%.3f")
    print(t.pivot_table(index="variant", columns="scope", values="slope_pct_per_km").reindex(VARIANTS).round(2).to_string())
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
