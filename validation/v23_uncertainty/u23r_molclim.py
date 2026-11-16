# -*- coding: utf-8 -*-
"""Review (code): reproduce the molecular-model climatology of calibration/rayleigh/uncertainty.py.

The coefficients _MOL_COEF were fitted on the q15d per-night CAMS-vs-US-1976 ln-biases
(C:/DATA/Projects/202610_CHM15k_stability/q15/q15d_nights.csv, windows in m ASL) by a script that is
not in the repository. This refits the published model with the window centre taken
  (a) above sea level and (b) above ground level,
to find which convention the shipped coefficients correspond to -- the module evaluates it with the
window centre in km AGL. Also reports the residual scatter, the realised RMS, and the bias the
latitude clip (35-62 N) implies for stations outside the fit range is NOT testable here (no nights).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
OUT = Path("C:/DATA/Projects/202610_v23/uncertainty_review")
Q15D = Path("C:/DATA/Projects/202610_CHM15k_stability/q15/q15d_nights.csv")
SHIPPED = (-0.421, 1.516, 0.799, 0.080, -0.017, -0.493, -0.210, 0.008)


def design(lat, doy, zc, h):
    w = 2 * np.pi / 365.25
    c, s = np.cos(w * doy), np.sin(w * doy)
    return np.column_stack([np.ones_like(lat), c, s, lat - 47, (lat - 47) * c, zc - 4.5,
                            (zc - 4.5) * c, h])


def main():
    import pandas as pd
    cen = {f"{s['wmo']}_{s['ident']}": s for s in
           json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))}
    d = pd.read_csv(Q15D)
    d = d[np.isfinite(d.b_mol) & np.isfinite(d.zb) & np.isfinite(d.zt)]
    d = d[d.stream.isin(cen)]
    lat = np.array([cen[s]["lat"] for s in d.stream])
    h = np.array([cen[s].get("alt", 0.0) for s in d.stream]) / 1000.0
    doy = np.array([(datetime(1970, 1, 1) + timedelta(days=float(t))).timetuple().tm_yday
                    for t in d.t])
    y = 100.0 * d.b_mol.values
    zc_asl = 0.5 * (d.zb.values + d.zt.values) / 1000.0
    zc_agl = zc_asl - h
    lines = [f"{len(d)} nights, {d.stream.nunique()} units, latitude {lat.min():.1f}-{lat.max():.1f} N, "
             f"realised RMS of b {np.sqrt(np.mean(y ** 2)):.2f} %"]
    for lab, zc in (("window centre ASL", zc_asl), ("window centre AGL", zc_agl)):
        X = design(lat, doy, zc, h)
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ coef
        pred_ship = X @ np.array(SHIPPED)
        lines.append(f"{lab}: coef " + ", ".join(f"{c:+.3f}" for c in coef)
                     + f" | residual sd {np.std(res):.3f} % | shipped coefficients applied with this "
                     f"convention: residual sd {np.std(y - pred_ship):.3f} %, mean {np.mean(y - pred_ship):+.3f} %")
    # what the convention mismatch does at high stations (shipped coefficients)
    for hs in (0.5, 1.0, 1.6, 2.1):
        dz = -SHIPPED[5] * hs
        lines.append(f"station at {hs:.1f} km: window centre AGL instead of ASL shifts b by "
                     f"{dz:+.2f} % (+/- {abs(SHIPPED[6]) * hs:.2f} % seasonal)")
    print("\n".join(lines))
    (OUT / "u23r_molclim.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
