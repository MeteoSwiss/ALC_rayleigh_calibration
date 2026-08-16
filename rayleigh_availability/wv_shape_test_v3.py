# -*- coding: utf-8 -*-
"""TEST v3 : CL61 vs CHM15k, avec l'aerosol EN COVARIABLE (le confondant, enfin controle).

HISTORIQUE DES DEUX ECHECS (a ne pas repeter).
  v1 (wv_shape_test.py) : pente de ln(S_910/S_1064) regressee sur le seul PWV -> Payerne +344 %,
      Lindenberg +20 %. Incompatibles : la croissance hygroscopique de l'aerosol est correlee au
      PWV et le rapport 910/1064 y est tres sensible.
  v2 (wv_pair_910_test.py) : CL61 vs CL31 (meme longueur d'onde, l'aerosol s'annule) -> R = -7,9,
      hors de toute plage physique (|R| <= ~1,5). Le CL31 de Payerne a son propre dark (-14 a
      -60 % en proche portee) et un faible SNR : sa forme depend du NIVEAU de signal, donc de la
      saison, donc du PWV.

CE QUE FAIT v3. Regression MULTIPLE de la pente observee sur DEUX predicteurs modelises, tous
deux calcules nuit par nuit a partir des memes CAMS que le pipeline :
    pente_obs(nuit) = a * pente_WV(nuit) + b * pente_AER(nuit) + c
  * pente_WV  : pente en portee de ln T2_wv (modele complet de l'instrument teste) ;
  * pente_AER : pente en portee de ln(beta_aer_910 / beta_aer_1064) (les deux colonnes CAMS
                deja stockees dans les caches) -> le confondant spectral aerosol, explicite.
Le coefficient **a** est alors la FRACTION de l'absorption WV modelisee que l'instrument subit
reellement, nette de l'aerosol : a ~ 1 -> il voit tout (la correction est justifiee) ;
a ~ 0 -> il est aveugle. Les biais STATIQUES (overlap, etalonnage, dark constant) sont dans c.

Run : python rayleigh_availability/wv_shape_test_v3.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration.io.cams import find_cams_file                              # noqa: E402
from calibration.water_vapor_correction.water_vapor import (                # noqa: E402
    cams_water_vapor_profile, laser_spectrum_for, two_way_wv_transmission)

CACHE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/"
             "dark_clearsky/cache")
FIG = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"
LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"
CAMS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
DARK_NPZ = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/"
                "dark_profiles_payerne.npz")

SITES = {
    "Payerne": dict(wmo="0-20000-0-06610", ref="A", tst="C", lat=46.8137, lon=6.9425,
                    alt=491.0, tst_type="CL61", dark=True),
    "Lindenberg": dict(wmo="0-20000-0-10393", ref="0", tst="C", lat=52.21, lon=14.12,
                       alt=123.0, tst_type="CL61", dark=False),
}
GRID = np.arange(400.0, 3600.0, 30.0)
BAND = (500.0, 3000.0)


def load_cache(wmo, ident):
    hits = sorted(CACHE.glob(f"{wmo}_{ident}_*.npz"))
    if not hits:
        return None
    z = np.load(hits[-1])
    return dict(rng=z["rng"].astype(float), S=z["S"].astype(float), Aer=z["Aer"].astype(float),
                dates=[str(int(d)) for d in z["dates"]])


def dark_profile(ident, rng):
    """Dark MESURE sous capot (Payerne uniquement), interpole sur la grille, en rcs_0."""
    if not DARK_NPZ.exists():
        return np.zeros_like(rng)
    z = np.load(DARK_NPZ)
    key = f"{ident}_b_rcs"
    if key not in z:
        return np.zeros_like(rng)
    return np.interp(rng, z[f"{ident}_range"], np.nan_to_num(z[key]), left=0.0, right=0.0)


def slope(y, m):
    if m.sum() < 8:
        return np.nan
    x = GRID[m] / 1000.0
    xx = x - x.mean()
    yy = y[m] - np.nanmean(y[m])
    v = np.sum(xx * xx)
    return float(np.sum(xx * yy) / v) if v > 0 else np.nan


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.0))
    for si, (site, sp) in enumerate(SITES.items()):
        ref, tst = load_cache(sp["wmo"], sp["ref"]), load_cache(sp["wmo"], sp["tst"])
        if ref is None or tst is None:
            print(f"[{site}] cache manquant")
            continue
        common = sorted(set(ref["dates"]) & set(tst["dates"]))
        ir = {d: i for i, d in enumerate(ref["dates"])}
        it = {d: i for i, d in enumerate(tst["dates"])}
        # dark mesure (Payerne) : on le retire des DEUX profils avant tout
        d_ref = dark_profile(sp["ref"], ref["rng"]) if sp["dark"] else 0.0
        d_tst = dark_profile(sp["tst"], tst["rng"]) if sp["dark"] else 0.0
        lam0, fwhm = laser_spectrum_for(sp["tst_type"], 910.0)

        rows = []
        for d in common:
            sr = np.interp(GRID, ref["rng"], ref["S"][ir[d]] - d_ref, left=np.nan, right=np.nan)
            st = np.interp(GRID, tst["rng"], tst["S"][it[d]] - d_tst, left=np.nan, right=np.nan)
            ar = np.interp(GRID, ref["rng"], ref["Aer"][ir[d]], left=np.nan, right=np.nan)
            at = np.interp(GRID, tst["rng"], tst["Aer"][it[d]], left=np.nan, right=np.nan)
            ok = np.isfinite(sr) & np.isfinite(st) & (sr > 0) & (st > 0)
            m = ok & (GRID >= BAND[0]) & (GRID <= BAND[1])
            if m.sum() < 30:
                continue
            cams = find_cams_file(CAMS, d)
            if cams is None:
                continue
            t0 = np.datetime64(f"{d[:4]}-{d[4:6]}-{d[6:8]}T00:00:00")
            prof = cams_water_vapor_profile(cams, sp["lat"], sp["lon"], t0,
                                            t0 + np.timedelta64(24, "h"))
            if prof is None:
                continue
            h_wv, n_wv = prof
            t2 = two_way_wv_transmission(sp["alt"] + GRID, sp["alt"], h_wv, n_wv, LUT, lam0, fwhm)
            s_wv = slope(np.log(np.maximum(t2, 1e-300)), m)
            # confondant aerosol : rapport spectral CAMS 910/1064, meme grille
            with np.errstate(divide="ignore", invalid="ignore"):
                la = np.log(np.where((at > 0) & (ar > 0), at / ar, np.nan))
            s_aer = slope(la, m & np.isfinite(la))
            s_obs = slope(np.log(st / sr), m)
            hm = h_wv >= sp["alt"]
            pwv = float(np.trapezoid(n_wv[hm], h_wv[hm]) * 18.015 / 6.022e23 * 1e-3)
            if np.isfinite(s_obs) and np.isfinite(s_wv):
                rows.append((s_obs, s_wv, s_aer if np.isfinite(s_aer) else 0.0, pwv))

        if len(rows) < 20:
            print(f"[{site}] {len(rows)} nuits — insuffisant")
            continue
        y = np.array([r[0] for r in rows])
        xw = np.array([r[1] for r in rows])
        xa = np.array([r[2] for r in rows])
        pwv = np.array([r[3] for r in rows])
        print(f"\n[{site}] {len(rows)} nuits, PWV {pwv.min():.1f}-{pwv.max():.1f} mm")

        def fit(X, names):
            A = np.column_stack(X + [np.ones(len(y))])
            beta, *_ = np.linalg.lstsq(A, y, rcond=None)
            r = y - A @ beta
            s2 = float(r @ r) / max(len(y) - A.shape[1], 1)
            cov = s2 * np.linalg.pinv(A.T @ A)
            for k, nm in enumerate(names):
                print(f"   {nm:24s} = {beta[k]:+.3f} ± {np.sqrt(cov[k, k]):.3f}")
            return beta, np.sqrt(np.diag(cov))

        print("  -- sans controle aerosol (= le test v1, pour memoire) --")
        fit([xw], ["fraction WV vue (a)"])
        print("  -- AVEC controle aerosol --")
        beta, se = fit([xw, xa], ["fraction WV vue (a)", "coef aerosol (b)"])

        ax = axes[si]
        # residus partiels : y nettoye de l'aerosol, contre le predicteur WV
        y_adj = y - beta[1] * xa
        ax.plot(xw, y_adj, "o", ms=5, color="#d62728", label="observe (net d'aerosol)")
        xx = np.linspace(np.nanmin(xw), np.nanmax(xw), 10)
        ax.plot(xx, beta[0] * xx + beta[2], "-", color="#d62728", lw=2,
                label=f"ajustement : a = {beta[0]:+.2f} ± {se[0]:.2f}")
        ax.plot(xx, 1.0 * xx + beta[2], "--", color="#1f77b4", lw=2,
                label="a = 1 (l'instrument voit toute la WV)")
        ax.plot(xx, 0.0 * xx + beta[2], ":", color="#2ca02c", lw=2,
                label="a = 0 (aveugle)")
        ax.set_xlabel("pente modelisee de ln T2$_{WV}$ [1/km]")
        ax.set_ylabel("pente observee de ln(S$_{910}$/S$_{1064}$) [1/km]")
        ax.set_title(f"{site} — fraction WV vue = {beta[0]:+.2f} ± {se[0]:.2f}  (n={len(y)})")
        ax.legend(fontsize=8.5)
        ax.grid(alpha=0.3)

    fig.suptitle("Test v3 : fraction d'absorption WV reellement subie par le CL61, "
                 "aerosol en covariable", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = FIG / "wv_shape_test_v3.png"
    fig.savefig(p, dpi=140)
    print(f"\nfigure -> {p}")


if __name__ == "__main__":
    main()
