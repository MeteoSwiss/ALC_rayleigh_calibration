# -*- coding: utf-8 -*-
"""Test lambda0 +-0.3 nm et SANS correction WV : la pente dC/dCBH des CL61 bouge-t-elle ?

PRINCIPE (semi-analytique, fidele au pipeline). Dans la calibration nuage, la correction WV est
appliquee AVANT l'integrale : beta /= T2_wv(z). Sur l'epaisseur du retour nuage (~300 m) T2 est
lisse, donc l'effet sur l'integrale est le facteur scalaire T2_wv(CBH) : ln C_L contient
-ln T2_wv(CBH; lambda0, FWHM). On peut donc rejouer chaque scene EXISTANTE du run reseau avec un
lambda0 decale (ou sans correction) en recalculant seulement T2_wv(CBH) — avec les MEMES
fonctions (two_way_wv_transmission), la MEME LUT et les MEMES CAMS 0.4 deg que le pipeline —
sans relancer une seule calibration :

    ln C(variante) = ln C(nominal) - [ln T2_var(CBH) - ln T2_nom(CBH)]
    ln C(sans WV)  = ln C(nominal) + ln T2_nom(CBH)

Pente testee : dC/dCBH en pool par type avec demoyennage par unite (la meme statistique que le
rapport d'heterogeneite). Controles : CL31 (FWHM 6 nm) et CL51 (3.4 nm), predits insensibles au
decalage de lambda0. Approximation assumee : T2 evalue a la CBH (pas re-integre sur le retour),
acceptance des scenes inchangee — exacte au premier ordre pour une pente.

Run : python rayleigh_availability/wv_lambda_sweep.py [--nmax-ctl 2000] [--workers 8]
Sorties : <DATA>/rayleigh_availability/wv_sweep_scenes.csv + table imprimee
          + doc/reports/figs_cbh_heterogeneity/wv_lambda_sweep.png
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
    cams_water_vapor_profile, laser_spectrum_for, two_way_wv_transmission)

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
FIG = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"
CAMS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
OFFSETS = (-0.30, -0.15, 0.0, +0.15, +0.30)          # nm autour du lambda0 nominal du type


def scene_t2(job):
    """ln T2_wv(CBH) pour chaque decalage de lambda0, pour UNE scene. None si CAMS absent."""
    date, lat, lon, alt, cbh, itype = job
    cams = find_cams_file(CAMS, date)
    if cams is None:
        return None
    t0 = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}T00:00:00")
    t1 = t0 + np.timedelta64(24, "h")
    prof = cams_water_vapor_profile(cams, lat, lon, t0, t1)
    if prof is None:
        return None
    h_wv, n_wv = prof
    lam0, fwhm = laser_spectrum_for(itype, 910.0)
    grid = alt + np.arange(0.0, cbh + 400.0, 30.0)   # ASL, du sol a la CBH
    out = []
    for off in OFFSETS:
        t2 = two_way_wv_transmission(grid, alt, h_wv, n_wv, LUT, lam0 + off, fwhm)
        i = int(np.argmin(np.abs((grid - alt) - cbh)))
        v = float(t2[i])
        out.append(np.log(v) if np.isfinite(v) and v > 0 else np.nan)
    return out


def fe_slope(lnc, cbh_km, unit):
    """Pente poolée avec demoyennage par unite (>= 8 scenes), en %/km, + SE OLS."""
    lnc = np.asarray(lnc, float)
    cbh_km = np.asarray(cbh_km, float)
    unit = np.asarray(unit)
    xd, yd = [], []
    for u in np.unique(unit):
        m = (unit == u) & np.isfinite(lnc) & np.isfinite(cbh_km)
        if m.sum() >= 8:
            xd.append(cbh_km[m] - cbh_km[m].mean())
            yd.append(lnc[m] - lnc[m].mean())
    if not xd:
        return np.nan, np.nan, 0
    x = np.concatenate(xd)
    y = np.concatenate(yd)
    vx = np.sum(x ** 2)
    if vx <= 0:
        return np.nan, np.nan, x.size
    s = np.sum(x * y) / vx
    r = y - s * x
    se = np.sqrt(np.sum(r ** 2) / max(x.size - 2, 1) / vx)
    return 100.0 * s, 100.0 * se, x.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax-ctl", type=int, default=2000,
                    help="scenes max par type de CONTROLE (CL31/CL51) ; CL61 = toutes")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    import pandas as pd
    scenes = pd.read_csv(REPO / "rayleigh_availability" / "cbh_native" / "cbh_scenes.csv")
    census = {f"{s['wmo']}_{s['ident']}": s for s in
              json.load(open(REPO / "validation" / "scope_l1_2026_census.json",
                             encoding="utf-8"))}
    scenes = scenes[np.isfinite(scenes["cbh_med"]) & (scenes["cal_value"] > 0)]

    # Controles CL31/CL51 : echantillonner des UNITES ENTIERES (le demoyennage par unite exige
    # >= 8 scenes par unite ; un tirage aleatoire de scenes les eparpille et tue le fit).
    rng = np.random.default_rng(0)
    parts = [scenes[scenes["type"] == "CL61"]]
    for t in ("CL51", "CL31"):
        st = scenes[scenes["type"] == t]
        counts = st.groupby("key").size()
        units = counts[counts >= 50].index.to_numpy()
        n_units = max(1, args.nmax_ctl // 130)
        if len(units) > n_units:
            units = rng.choice(units, n_units, replace=False)
        parts.append(st[st["key"].isin(units)])
    scenes = pd.concat(parts, ignore_index=True)

    jobs, keep = [], []
    for i, r in scenes.iterrows():
        c = census.get(r["key"])
        if c is None:
            continue
        jobs.append((str(int(r["date"])), float(c["lat"]), float(c["lon"]),
                     float(c["alt"]), float(r["cbh_med"]), r["type"]))
        keep.append(i)
    scenes = scenes.loc[keep].reset_index(drop=True)
    print(f"{len(jobs)} scenes a rejouer "
          f"({(scenes['type'] == 'CL61').sum()} CL61, le reste en controle)")

    res = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, out in enumerate(ex.map(scene_t2, jobs, chunksize=16)):
            res[i] = out
    ok = np.array([r is not None for r in res])
    print(f"{ok.sum()} scenes avec CAMS/T2 valides")
    lnT = np.full((len(jobs), len(OFFSETS)), np.nan)
    for i, r in enumerate(res):
        if r is not None:
            lnT[i] = r
    scenes = scenes.assign(**{f"lnT2_{k}": lnT[:, j] for j, k in
                              enumerate(["m30", "m15", "nom", "p15", "p30"])})
    scenes.to_csv(DATA / "wv_sweep_scenes.csv", index=False)

    inom = OFFSETS.index(0.0)
    print("\n=== pente dC/dCBH poolée (demoyennage par unite) [%/km] ===")
    header = (f"{'type':6s} {'nominal':>12s} " +
              " ".join(f"{f'l0{o:+.2f}':>12s}" for o in OFFSETS if o != 0.0) +
              f" {'SANS WV':>12s}")
    print(header)
    table = {}
    for t in ("CL61", "CL51", "CL31"):
        st = scenes[scenes["type"] == t]
        lnc = np.log(st["cal_value"].to_numpy(float))
        cbh = st["cbh_med"].to_numpy(float) / 1000.0
        unit = st["key"].to_numpy()
        lnt_nom = st["lnT2_nom"].to_numpy(float)
        row = {}
        s0, e0, n0 = fe_slope(lnc, cbh, unit)
        row["nominal"] = (s0, e0, n0)
        cells = [f"{s0:+7.2f}±{e0:4.2f}"]
        for j, off in enumerate(OFFSETS):
            if off == 0.0:
                continue
            lnc_v = lnc - (lnT[scenes["type"] == t][:, j] - lnt_nom)
            s, e, n = fe_slope(lnc_v, cbh, unit)
            row[f"l0{off:+.2f}"] = (s, e, n)
            cells.append(f"{s:+7.2f}±{e:4.2f}")
        lnc_no = lnc + lnt_nom
        s, e, n = fe_slope(lnc_no, cbh, unit)
        row["noWV"] = (s, e, n)
        cells.append(f"{s:+7.2f}±{e:4.2f}")
        table[t] = row
        print(f"{t:6s} " + " ".join(f"{c:>12s}" for c in cells) + f"   (n={n0})")

    # figure paysage : pente vs decalage lambda0 par type + point SANS WV
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12.5, 6.0))
    cols = {"CL61": "#d62728", "CL51": "#7d3c98", "CL31": "#1f77b4"}
    for t, col in cols.items():
        xs, ys, es = [], [], []
        for off in OFFSETS:
            k = "nominal" if off == 0.0 else f"l0{off:+.2f}"
            s, e, _ = table[t][k]
            xs.append(off)
            ys.append(s)
            es.append(e)
        ax.errorbar(xs, ys, yerr=es, fmt="o-", color=col, capsize=3, lw=1.6,
                    label=f"{t} (λ0 nominal {laser_spectrum_for(t, 910.0)[0]:.2f} nm, "
                          f"FWHM {laser_spectrum_for(t, 910.0)[1]:.1f})")
        s, e, _ = table[t]["noWV"]
        ax.errorbar([0.42], [s], yerr=[e], fmt="s", color=col, ms=9, capsize=3)
    ax.axvline(0.38, color="#666", lw=0.8, ls=":")
    ax.text(0.42, ax.get_ylim()[1], "SANS\ncorrection WV", fontsize=9, ha="center", va="top")
    ax.axhline(0, color="#444", ls=":", lw=1)
    ax.set_xlabel("décalage de λ0 par rapport au nominal [nm]")
    ax.set_ylabel("pente dC/dCBH poolée [%/km]")
    ax.set_title("Balayage λ0 et suppression de la correction WV — effet sur la pente dC/dCBH\n"
                 "(scènes réelles rejouées semi-analytiquement, T2 via la LUT et CAMS 0,4° du pipeline)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / "wv_lambda_sweep.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
