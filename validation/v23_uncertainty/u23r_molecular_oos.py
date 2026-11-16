# -*- coding: utf-8 -*-
"""Review (claim: molecular term): the US-1976 climatology OUT OF SAMPLE.

The climatology was fitted on 127 CHM15k units at 37-55 N (q15d). Here the realised per-night ln-bias
b (CAMS T/p vs US-1976, on the night's own published window, q15d machinery) is computed for CHM15k
units NOT in that fit -- the fresh review units and every census CHM15k north of the fit range
(55-77 N, where the module clips the latitude at 62 N) -- and compared with the climatology evaluated
(a) as the module does (window centre AGL) and (b) with the fit's own convention (ASL).
Stations outside the European CAMS box (Bonaire, Ontario) are skipped.

Output: u23r_molecular_oos.csv / .txt
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))
sys.path.insert(0, str(HERE))

warnings.filterwarnings("ignore")
OUT = Path("C:/DATA/Projects/202610_v23/uncertainty_review")


def main():
    import pandas as pd
    import q15d_atmosphere as Q
    from common import load_constants
    from calibration.rayleigh import uncertainty as U
    from u23r_harness import FRESH_CHM15K
    import u23_validation as V
    Q.STD = Q.StdAtm()
    cen = {f"{s['wmo']}_{s['ident']}": s for s in
           json.load(open(REPO / "validation" / "scope_l1_2026_census.json", encoding="utf-8"))}
    fitted = set(pd.read_csv(Q.Q15 / "q15d_nights.csv", usecols=["stream"]).stream.unique())
    design = {f"{w}_{i}" for w, i in V.UNITS}
    keys = [f"{w}_{i}" for w, i in FRESH_CHM15K]
    keys += [k for k, s in cen.items() if s["type"] == "CHM15k" and s["lat"] > 55.0 and k not in keys]
    cams = Q.load_cams()
    kidx = {k: i for i, k in enumerate(cams["keys"])}
    rows = []
    for key in keys:
        s = cen[key]
        if not (27.0 <= s["lat"] <= 80.0 and -30.0 <= s["lon"] <= 45.0):
            continue                                   # outside the European CAMS box
        ip = kidx.get(key)
        con = load_constants(s["wmo"], s["ident"])
        if ip is None or con is None:
            continue
        alt = float(s.get("alt", 0.0))
        for t, t0, t1, zb, zt in zip(con["t"], con["t0"], con["t1"], con["zb"], con["zt"]):
            if not (np.isfinite(zb) and np.isfinite(zt) and np.isfinite(t0) and np.isfinite(t1)):
                continue
            prof = Q.cams_night(cams, ip, t0, t1)
            if prof is None:
                continue
            r = Q.night_cams_terms(prof, alt, zb, zt)
            if r is None:
                continue
            doy = (Q.EPOCH + Q.timedelta(days=float(t))).timetuple().tm_yday
            zc_asl = 0.5 * (zb + zt) / 1000.0
            rows.append(dict(key=key, site=s["site"], lat=s["lat"], alt=alt, t=t, doy=doy,
                             zc_asl=zc_asl, b=100 * r["b_mol"],
                             clim_agl=100 * U.molecular_bias_climatology(s["lat"], doy, zc_asl - alt / 1000.0, alt / 1000.0),
                             clim_asl=100 * U.molecular_bias_climatology(s["lat"], doy, zc_asl, alt / 1000.0),
                             in_fit=key in fitted, in_design=key in design))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "u23r_molecular_oos.csv", index=False)
    lines = []
    for lab, m in (("fitted units (in sample)", d.in_fit),
                   ("not in the fit, <= 55 N", ~d.in_fit & (d.lat <= 55)),
                   ("not in the fit, 55-62 N", ~d.in_fit & (d.lat > 55) & (d.lat <= 62)),
                   ("not in the fit, > 62 N (clipped)", ~d.in_fit & (d.lat > 62)),
                   ("stations > 1 km ASL", d.alt > 1000)):
        g = d[m]
        if len(g) < 10:
            continue
        for conv in ("clim_agl", "clim_asl"):
            e = g.b - g[conv]
            # the budget's sigma_mol = sqrt(b_clim^2 + 1.35^2): does it cover the realised |b|?
            sig = np.sqrt(g[conv] ** 2 + 1.35 ** 2)
            cover = np.mean(np.abs(g.b) <= 2 * sig)
            lines.append(f"{lab:34s} {conv}: n {len(g):5d} units {g.key.nunique():3d} | realised RMS b "
                         f"{np.sqrt(np.mean(g.b ** 2)):.2f} % | error of the climatology: mean {e.mean():+.2f} %, "
                         f"sd {e.std():.2f} % | |b| <= 2 sigma_mol on {100 * cover:.1f} % of the nights; "
                         f"pull of b: {np.std(g.b / sig):.2f}")
    per = d[~d.in_fit].groupby(["key", "site", "lat", "alt"]).apply(
        lambda g: pd.Series(dict(n=len(g), rms_b=np.sqrt(np.mean(g.b ** 2)),
                                 err_agl=np.mean(g.b - g.clim_agl), err_asl=np.mean(g.b - g.clim_asl),
                                 sd_resid_asl=np.std(g.b - g.clim_asl)))).reset_index()
    lines.append(per.round(2).sort_values("lat").to_string(index=False))
    print("\n".join(lines))
    (OUT / "u23r_molecular_oos.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
