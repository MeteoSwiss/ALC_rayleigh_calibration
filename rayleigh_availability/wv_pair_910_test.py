# -*- coding: utf-8 -*-
"""TEST DECISIF v2 : CL61 vs CL31 a Payerne — DEUX instruments a 910 nm, meme aerosol.

POURQUOI CETTE VERSION. Le test CL61-vs-CHM15k (wv_shape_test.py) est ruine par un confondant :
910 vs 1064 nm ne voient pas le meme aerosol, et la croissance hygroscopique de l'aerosol est
correlee au PWV -> Payerne et Lindenberg donnent des reponses incompatibles (+344 % vs +20 %).

Ici on compare deux instruments a la MEME longueur d'onde nominale (910 nm), co-localises a
Payerne : le CL61 (C, raie etroite ~910,6-910,7) et le CL31 (B, multimode large 5-7 nm centre
909,7). Ils voient le meme aerosol, la meme croissance hygroscopique, le meme moleculaire : tout
cela s'annule dans le rapport. Il ne reste que la difference d'absorption VAPEUR D'EAU due a
leurs spectres d'emission differents.

    ln[S_CL61(z)/S_CL31(z)] = const(z, statique : overlap, etalonnage) + [ln T2_C(z) - ln T2_B(z)]

En notant f_X la FRACTION de l'absorption modelisee que l'instrument X subit reellement :
    d/dz de la partie WV = f_C * dlnT2_C/dz - f_B * dlnT2_B/dz
et comme tau_C ~ tau_B (spectres differents mais absorption modelisee comparable), on mesure
essentiellement (f_C - f_B) en unites de la pente modelisee du CL61 :

    R = [d(pente observee)/dPWV] / [d(pente modelisee CL61)/dPWV]  ~  f_C - f_B

  * analyse NUAGE vraie (f_C = 0,08 ; f_B = 0,80)  ->  R ~ -0,7  (le CL31 est plus attenue)
  * les deux instruments subissent la meme chose   ->  R ~  0
Le signe est discriminant et l'amplitude aussi ; aucun etalonnage n'intervient.

Run : python rayleigh_availability/wv_pair_910_test.py
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

from calibration.io.cams import find_cams_file                              # noqa: E402
from calibration.water_vapor_correction.water_vapor import (                # noqa: E402
    cams_water_vapor_profile, laser_spectrum_for, two_way_wv_transmission)
from dark_from_clearsky import process_night                                # noqa: E402

CACHE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/"
             "dark_clearsky/cache")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
FIG = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
CAMS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
L1 = "D:/E-PROFILE_L1_2026"

WMO, LAT, LON, ALT = "0-20000-0-06610", 46.8137, 6.9425, 491.0
GRID = np.arange(300.0, 3600.0, 30.0)
BANDS = {"0,5-3 km": (500.0, 3000.0), "1-2,5 km": (1000.0, 2500.0)}


def load_cache(ident):
    hits = sorted(CACHE.glob(f"{WMO}_{ident}_*.npz"))
    if not hits:
        return None
    z = np.load(hits[-1])
    return dict(rng=z["rng"].astype(float), S=z["S"].astype(float),
                dates=[str(int(d)) for d in z["dates"]])


def slope(x_km, y, m):
    if m.sum() < 8:
        return np.nan
    xx = x_km[m] - x_km[m].mean()
    yy = y[m] - y[m].mean()
    v = np.sum(xx * xx)
    return float(np.sum(xx * yy) / v) if v > 0 else np.nan


def main():
    c = load_cache("C")
    if c is None:
        print("cache CL61 manquant")
        return
    print(f"CL61 : {len(c['dates'])} nuits claires en cache")

    # --- CL31 (B) : pas de nuits claires dans le run reseau -> on extrait sur les MEMES nuits
    b_cache = OUT / "dark_clearsky" / "cache" / f"{WMO}_B_pairtest.npz"
    if b_cache.exists():
        z = np.load(b_cache)
        b = dict(rng=z["rng"].astype(float), S=z["S"].astype(float),
                 dates=[str(int(d)) for d in z["dates"]])
        print(f"CL31 : {len(b['dates'])} nuits depuis le cache dedie")
    else:
        jobs = [dict(date=d, wmo=WMO, ident="B", l1_root=L1, cams_folder=CAMS)
                for d in c["dates"]]
        got = []
        with ProcessPoolExecutor(max_workers=10) as ex:
            for r in ex.map(process_night, jobs, chunksize=4):
                if r is not None:
                    got.append(r)
        if len(got) < 20:
            print(f"CL31 : seulement {len(got)} nuits exploitables, abandon")
            return
        rng0 = got[0]["rng"]
        got = [g for g in got if g["rng"].size == rng0.size]
        b = dict(rng=rng0, S=np.vstack([g["S"] for g in got]),
                 dates=[g["date"] for g in got])
        np.savez_compressed(b_cache, rng=b["rng"], S=b["S"],
                            dates=np.array([int(d) for d in b["dates"]]))
        print(f"CL31 : {len(b['dates'])} nuits extraites -> {b_cache.name}")

    common = sorted(set(b["dates"]) & set(c["dates"]))
    print(f"nuits communes CL31/CL61 : {len(common)}")
    ib = {d: i for i, d in enumerate(b["dates"])}
    ic = {d: i for i, d in enumerate(c["dates"])}

    lam_c, fw_c = laser_spectrum_for("CL61", 910.0)
    lam_b, fw_b = laser_spectrum_for("CL31", 910.0)
    print(f"spectres modelises : CL61 ({lam_c}, {fw_c})  CL31 ({lam_b}, {fw_b})")

    rows = []
    for d in common:
        sb = np.interp(GRID, b["rng"], b["S"][ib[d]], left=np.nan, right=np.nan)
        sc = np.interp(GRID, c["rng"], c["S"][ic[d]], left=np.nan, right=np.nan)
        ok = np.isfinite(sb) & np.isfinite(sc) & (sb > 0) & (sc > 0)
        if ok.sum() < 40:
            continue
        cams = find_cams_file(CAMS, d)
        if cams is None:
            continue
        t0 = np.datetime64(f"{d[:4]}-{d[4:6]}-{d[6:8]}T00:00:00")
        prof = cams_water_vapor_profile(cams, LAT, LON, t0, t0 + np.timedelta64(24, "h"))
        if prof is None:
            continue
        h_wv, n_wv = prof
        m = h_wv >= ALT
        pwv_mm = float(np.trapezoid(n_wv[m], h_wv[m]) * 18.015 / 6.022e23 * 1e-3)
        grid_asl = ALT + GRID
        t2c = np.log(np.maximum(two_way_wv_transmission(
            grid_asl, ALT, h_wv, n_wv, LUT, lam_c, fw_c), 1e-300))
        t2b = np.log(np.maximum(two_way_wv_transmission(
            grid_asl, ALT, h_wv, n_wv, LUT, lam_b, fw_b), 1e-300))
        lr = np.where(ok, np.log(sc / sb), np.nan)
        row = dict(date=d, pwv=pwv_mm)
        for name, (z0, z1) in BANDS.items():
            mm = ok & (GRID >= z0) & (GRID <= z1) & np.isfinite(lr)
            row[f"obs_{name}"] = slope(GRID / 1000.0, lr, mm)
            row[f"mdif_{name}"] = slope(GRID / 1000.0, t2c - t2b, mm)   # difference modelisee
            row[f"mc_{name}"] = slope(GRID / 1000.0, t2c, mm)           # echelle : le CL61 seul
        rows.append(row)

    print(f"{len(rows)} nuits exploitables")
    if len(rows) < 15:
        return
    pwv = np.array([r["pwv"] for r in rows])
    print(f"PWV {pwv.min():.1f}-{pwv.max():.1f} mm")

    def reg(y):
        good = np.isfinite(y)
        x = pwv[good] - pwv[good].mean()
        yy = y[good] - y[good].mean()
        v = np.sum(x * x)
        a = np.sum(x * yy) / v
        r = yy - a * x
        se = np.sqrt(np.sum(r * r) / max(good.sum() - 2, 1) / v)
        return a, se, int(good.sum())

    print("\n=== d(pente du log-rapport)/dPWV  [1/km par mm] ===")
    for name in BANDS:
        obs = np.array([r[f"obs_{name}"] for r in rows])
        mdif = np.array([r[f"mdif_{name}"] for r in rows])
        mc = np.array([r[f"mc_{name}"] for r in rows])
        a_o, se_o, n = reg(obs)
        a_d, _, _ = reg(mdif)
        a_c, _, _ = reg(mc)
        R = a_o / a_c if abs(a_c) > 1e-12 else np.nan
        seR = se_o / abs(a_c) if abs(a_c) > 1e-12 else np.nan
        print(f"  {name:9s} observe {a_o:+.3e} ± {se_o:.1e} (n={n}) | "
              f"modele CL61-CL31 {a_d:+.3e} | echelle CL61 {a_c:+.3e}")
        print(f"            -> R = (f_CL61 - f_CL31) = {R:+.2f} ± {seR:.2f}"
              f"   [nuage predit -0,7 ; instruments identiques 0]")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12.5, 6.0))
    name = "0,5-3 km"
    obs = np.array([r[f"obs_{name}"] for r in rows])
    mc = np.array([r[f"mc_{name}"] for r in rows])
    a_o, se_o, n = reg(obs)
    a_c, _, _ = reg(mc)
    ax.plot(pwv, obs, "o", ms=6, color="#d62728", label="observe : pente de ln(CL61/CL31)")
    xx = np.linspace(pwv.min(), pwv.max(), 10)
    ax.plot(xx, np.median(obs) + a_o * (xx - np.median(pwv)), "-", color="#d62728", lw=2,
            label=f"ajustement observe ({a_o:+.2e} /km/mm)")
    ax.plot(xx, np.median(obs) - 0.72 * a_c * (xx - np.median(pwv)), "--", color="#1f77b4", lw=2,
            label="attendu si l'analyse nuage est vraie (f_C−f_B = −0,72)")
    ax.plot(xx, np.median(obs) + 0 * xx, ":", color="#2ca02c", lw=2,
            label="attendu si les deux subissent la meme absorption")
    ax.set_xlabel("PWV de la nuit [mm]")
    ax.set_ylabel(f"pente de ln(S$_{{CL61}}$/S$_{{CL31}}$) sur {name} [1/km]")
    ax.set_title("Test decisif : CL61 vs CL31 a Payerne — memes 910 nm, meme aerosol\n"
                 f"R = f_CL61 − f_CL31 = {a_o/a_c:+.2f} ± {se_o/abs(a_c):.2f}  (n={n})")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    p = FIG / "wv_pair_910_test.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
