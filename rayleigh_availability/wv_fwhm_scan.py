# -*- coding: utf-8 -*-
"""Scan FWHM a lambda0 = 910.74 nm (la valeur MESUREE au spectrometre, appendice A1) :
l'absorption WV effective du CL61 s'effondre-t-elle pour une emission sub-nm ?

Le spectrometre borne le FWHM vrai a < 1.5 nm sans le resoudre ; Vaisala revendique une
"very narrow bandwidth" concue pour un creux d'absorption. On calcule, sur les scenes nuage
CL61 reelles, la dependance CBH du terme WV (la "pente injectee") pour FWHM 0.05..1.5 nm :
si elle tombe vers ~8 % de la valeur a 1.0 nm, le modele spectral etroit EXPLIQUE le deficit
observe ; sinon, le creux n'existe pas a la resolution de la LUT (8.3 pm) et l'explication
est ailleurs (raies sub-LUT, spectre local different).

Run : python rayleigh_availability/wv_fwhm_scan.py [--workers 10]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration.io.cams import find_cams_file                                  # noqa: E402
from calibration.water_vapor_correction.water_vapor import (                    # noqa: E402
    cams_water_vapor_profile, two_way_wv_transmission)

FIG = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"
CAMS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
LAM0 = 910.74                                   # nm, MESURE (appendice A1, +-0.10)
FWHMS = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 1.5)


def scene_t2(job):
    date, lat, lon, alt, cbh = job
    cams = find_cams_file(CAMS, date)
    if cams is None:
        return None
    t0 = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}T00:00:00")
    prof = cams_water_vapor_profile(cams, lat, lon, t0, t0 + np.timedelta64(24, "h"))
    if prof is None:
        return None
    h_wv, n_wv = prof
    grid = alt + np.arange(0.0, cbh + 400.0, 30.0)
    i = int(np.argmin(np.abs((grid - alt) - cbh)))
    out = []
    for fw in FWHMS:
        t2 = two_way_wv_transmission(grid, alt, h_wv, n_wv, LUT, LAM0, fw)
        v = float(t2[i])
        out.append(np.log(v) if np.isfinite(v) and v > 0 else np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()

    import pandas as pd
    scenes = pd.read_csv(REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv")
    scenes = scenes[(scenes["type"] == "CL61") & np.isfinite(scenes["cbh_med"])
                    & (scenes["cal_value"] > 0)].reset_index(drop=True)
    census = {f"{s['wmo']}_{s['ident']}": s for s in
              json.load(open(REPO / "validation" / "scope_l1_2026_census.json",
                             encoding="utf-8"))}
    jobs, keep = [], []
    for i, r in scenes.iterrows():
        c = census.get(r["key"])
        if c:
            jobs.append((str(int(r["date"])), float(c["lat"]), float(c["lon"]),
                         float(c["alt"]), float(r["cbh_med"])))
            keep.append(i)
    scenes = scenes.loc[keep].reset_index(drop=True)
    print(f"{len(jobs)} scenes CL61, lambda0 = {LAM0} nm (mesure), FWHM = {FWHMS}")

    res = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, out in enumerate(ex.map(scene_t2, jobs, chunksize=16)):
            res[i] = out
    lnT = np.full((len(jobs), len(FWHMS)), np.nan)
    for i, r in enumerate(res):
        if r is not None:
            lnT[i] = r

    cbh = scenes["cbh_med"].to_numpy(float) / 1000.0
    unit = scenes["key"].to_numpy()
    print("\nFWHM [nm] | pente CBH du terme WV [%/km] | <T2> a la CBH mediane")
    slopes = []
    for j, fw in enumerate(FWHMS):
        xd, yd = [], []
        for u in np.unique(unit):
            m = (unit == u) & np.isfinite(lnT[:, j])
            if m.sum() >= 8:
                xd.append(cbh[m] - cbh[m].mean())
                yd.append(lnT[m][:, j] - lnT[m][:, j].mean())
        x = np.concatenate(xd)
        y = np.concatenate(yd)
        s = -100.0 * np.sum(x * y) / np.sum(x ** 2)     # terme dans ln C = -ln T2
        slopes.append(s)
        print(f"  {fw:5.2f}   |        {s:+7.2f}          |   {np.exp(np.nanmedian(lnT[:, j])):.3f}")
    print("\nRappel : le brut CL61 ne montre que ~8 % du terme a FWHM 1.0 "
          f"({0.08 * slopes[FWHMS.index(1.0)]:+.1f} %/km equivalent).")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    ax.plot(FWHMS, slopes, "o-", color="#d62728", lw=1.8)
    ax.axhline(0.08 * slopes[FWHMS.index(1.0)], color="#444", ls="--", lw=1.2,
               label="niveau vu par le signal brut CL61 (~8 % du terme a 1.0 nm)")
    ax.axvspan(1.5, max(FWHMS), color="#999", alpha=0.15)
    ax.set_xlabel("FWHM suppose de l'emission CL61 [nm]  (mesure : < 1,5 nm, non resolu)")
    ax.set_ylabel("pente CBH injectee par le terme WV [%/km]")
    ax.set_title(f"Scan FWHM a lambda0 = {LAM0} nm (mesure au Qmini) — "
                 "un sub-nm explique-t-il le deficit d'absorption du CL61 ?")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    p = FIG / "wv_fwhm_scan.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
