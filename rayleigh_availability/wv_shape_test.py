# -*- coding: utf-8 -*-
"""TEST DECISIF, SANS AUCUNE CONSTANTE DE CALIBRATION : le CL61 subit-il l'absorption WV ?

POURQUOI CE TEST. Deux analyses se contredisent :
  * NUAGE  : la pente dC/dCBH s'annule quand on retire la correction WV -> le brut CL61 ne
             porterait que ~8 % de l'absorption modelisee ;
  * RAYLEIGH (dashboard) : retirer la correction des DEUX etages degrade l'accord avec le
             CHM15k (-1,9 % -> +19 %) -> la correction ferait du vrai travail.
Les deux passent par une CONSTANTE de calibration, donc par la fenetre de fit, le dark, le
niveau... Ce test-ci n'en utilise aucune : il ne regarde que la FORME du rapport de profils
bruts entre le CL61 (910 nm) et le CHM15k (1064 nm) co-localises, et sa dependance a l'humidite.

PRINCIPE. rcs_0 ~ beta_att x C_L, donc
    ln[S_910(z)/S_1064(z)] = ln(beta_910/beta_1064) + ln T2_mol_ratio + ln T2_aer_ratio
                             + ln T2_wv(z) + const
La constante C_L disparait dans la derivee en z. Seul le terme WV depend de l'humidite du jour :
    d/dz ln T2_wv = -2 sigma n_wv(z)
On regresse donc la PENTE du log-rapport (0,5-3 km) contre le PWV de la nuit, et on compare le
coefficient obtenu a celui que PREDIT notre modele WV (meme regression sur les T2 modelises).
    ratio = coef_observe / coef_modelise = fraction de l'absorption modelisee reellement subie.
  ~1 -> l'instrument voit toute la WV (la correction est justifiee, l'analyse nuage a un defaut)
  ~0 -> l'instrument est aveugle (mitigation par conception confirmee independamment)

Confondants assumes et traites : l'overlap et les differences d'etalonnage sont STATIQUES -> ils
entrent dans l'ordonnee a l'origine de la regression, pas dans la pente vs PWV. Reste la
croissance hygroscopique de l'aerosol (correlee au PWV) : on la borne en refaisant le test sur
une bande haute (1,5-3,5 km) ou l'aerosol de couche limite pese beaucoup moins.

Entrees : les caches de nuits claires de dark_from_clearsky.py (S = profil median de nuit en
rcs_0), deja calcules pour Payerne A/C et Lindenberg 0/C.

Run : python rayleigh_availability/wv_shape_test.py
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

SITES = {
    "Payerne": dict(wmo="0-20000-0-06610", ref="A", tst="C", lat=46.8137, lon=6.9425, alt=491.0,
                    tst_type="CL61"),
    "Lindenberg": dict(wmo="0-20000-0-10393", ref="0", tst="C", lat=52.21, lon=14.12, alt=123.0,
                       tst_type="CL61"),
}
BANDS = {"0,5-3 km": (500.0, 3000.0), "1,5-3,5 km": (1500.0, 3500.0)}
GRID = np.arange(300.0, 3600.0, 30.0)


def load_cache(wmo, ident):
    hits = sorted(CACHE.glob(f"{wmo}_{ident}_*.npz"))
    if not hits:
        return None
    z = np.load(hits[-1])
    return dict(rng=z["rng"].astype(float), S=z["S"].astype(float),
                dates=[str(int(d)) for d in z["dates"]])


def pwv_and_t2(date, lat, lon, alt, itype):
    """(PWV [mm], pente modelisee de ln T2_wv par bande [1/km]) pour la nuit."""
    cams = find_cams_file(CAMS, date)
    if cams is None:
        return None
    t0 = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}T00:00:00")
    prof = cams_water_vapor_profile(cams, lat, lon, t0, t0 + np.timedelta64(24, "h"))
    if prof is None:
        return None
    h_wv, n_wv = prof
    m = h_wv >= alt
    if m.sum() < 5:
        return None
    # PWV = colonne de vapeur au-dessus de la station, en mm d'eau liquide
    pwv = float(np.trapezoid(n_wv[m], h_wv[m]) * 18.015 / 6.022e23 * 1000.0 / 1000.0)
    lam0, fwhm = laser_spectrum_for(itype, 910.0)
    grid_asl = alt + GRID
    t2 = two_way_wv_transmission(grid_asl, alt, h_wv, n_wv, LUT, lam0, fwhm)
    return pwv, np.log(np.maximum(t2, 1e-300))


def slope(x_km, y, m):
    """pente OLS de y vs x sur le masque m, en unites de y par km."""
    if m.sum() < 8:
        return np.nan
    xx = x_km[m] - x_km[m].mean()
    yy = y[m] - y[m].mean()
    v = np.sum(xx * xx)
    return float(np.sum(xx * yy) / v) if v > 0 else np.nan


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.0))
    results = {}
    for si, (site, spec) in enumerate(SITES.items()):
        ref = load_cache(spec["wmo"], spec["ref"])
        tst = load_cache(spec["wmo"], spec["tst"])
        if ref is None or tst is None:
            print(f"[{site}] cache manquant (ref={ref is not None}, tst={tst is not None})")
            continue
        common = sorted(set(ref["dates"]) & set(tst["dates"]))
        print(f"[{site}] {len(common)} nuits claires communes")
        if len(common) < 15:
            continue
        ir = {d: i for i, d in enumerate(ref["dates"])}
        it = {d: i for i, d in enumerate(tst["dates"])}

        rows = []
        for d in common:
            sr = np.interp(GRID, ref["rng"], ref["S"][ir[d]], left=np.nan, right=np.nan)
            st = np.interp(GRID, tst["rng"], tst["S"][it[d]], left=np.nan, right=np.nan)
            ok = np.isfinite(sr) & np.isfinite(st) & (sr > 0) & (st > 0)
            if ok.sum() < 40:
                continue
            lr = np.where(ok, np.log(np.abs(st / sr)), np.nan)
            got = pwv_and_t2(d, spec["lat"], spec["lon"], spec["alt"], spec["tst_type"])
            if got is None:
                continue
            pwv, ln_t2 = got
            row = dict(date=d, pwv=pwv)
            for name, (z0, z1) in BANDS.items():
                m = ok & (GRID >= z0) & (GRID <= z1) & np.isfinite(lr)
                row[f"obs_{name}"] = slope(GRID / 1000.0, lr, m)
                row[f"mod_{name}"] = slope(GRID / 1000.0, ln_t2, m)
            rows.append(row)

        if len(rows) < 15:
            print(f"[{site}] trop peu de nuits exploitables ({len(rows)})")
            continue
        pwv = np.array([r["pwv"] for r in rows])
        print(f"  PWV {pwv.min():.1f}-{pwv.max():.1f} mm sur {len(rows)} nuits")
        res_site = {}
        for name in BANDS:
            obs = np.array([r[f"obs_{name}"] for r in rows])
            mod = np.array([r[f"mod_{name}"] for r in rows])
            good = np.isfinite(obs) & np.isfinite(mod)
            # regression robuste (Theil-Sen) de la pente observee et modelisee contre le PWV
            def ts(y):
                x = pwv[good]
                yy = y[good]
                sl = []
                for i in range(len(x)):
                    dx = x[i + 1:] - x[i]
                    dy = yy[i + 1:] - yy[i]
                    k = np.abs(dx) > 1e-6
                    sl.extend((dy[k] / dx[k]).tolist())
                return float(np.median(sl)) if sl else np.nan
            a_obs, a_mod = ts(obs), ts(mod)
            frac = a_obs / a_mod if np.isfinite(a_mod) and abs(a_mod) > 1e-9 else np.nan
            # incertitude : bootstrap des nuits
            rs = np.random.default_rng(0)
            bs = []
            idx = np.where(good)[0]
            for _ in range(300):
                p = rs.choice(idx, idx.size, replace=True)
                x, yo, ym = pwv[p], obs[p], mod[p]
                xc = x - x.mean()
                vo = np.sum(xc * (yo - yo.mean())) / np.sum(xc * xc)
                vm = np.sum(xc * (ym - ym.mean())) / np.sum(xc * xc)
                if abs(vm) > 1e-9:
                    bs.append(vo / vm)
            se = float(np.std(bs)) if len(bs) > 50 else np.nan
            res_site[name] = (a_obs, a_mod, frac, se, int(good.sum()))
            print(f"  bande {name:10s} : d(pente)/dPWV observe {a_obs:+.5f} /km/mm | "
                  f"modele {a_mod:+.5f} | FRACTION VUE = {100*frac:+.0f} % +- {100*se:.0f} "
                  f"(n={good.sum()})")
        results[site] = res_site

        ax = axes[si]
        name = "0,5-3 km"
        obs = np.array([r[f"obs_{name}"] for r in rows])
        mod = np.array([r[f"mod_{name}"] for r in rows])
        ax.plot(pwv, obs, "o", ms=5, color="#d62728", label="observe (rapport de profils bruts)")
        ax.plot(pwv, mod - np.nanmedian(mod) + np.nanmedian(obs), "s", ms=4, color="#1f77b4",
                alpha=0.7, label="modele WV (recale en niveau)")
        a_obs, a_mod, frac, se, n = results[site][name]
        xx = np.linspace(pwv.min(), pwv.max(), 10)
        ax.plot(xx, np.nanmedian(obs) + a_obs * (xx - np.median(pwv)), "-", color="#d62728", lw=2)
        ax.plot(xx, np.nanmedian(obs) + a_mod * (xx - np.median(pwv)), "--", color="#1f77b4", lw=2)
        ax.set_xlabel("PWV de la nuit [mm]")
        ax.set_ylabel(f"pente de ln(S$_{{910}}$/S$_{{1064}}$) sur {name} [1/km]")
        ax.set_title(f"{site} — fraction de l'absorption WV\nreellement subie : "
                     f"{100*frac:+.0f} % ± {100*se:.0f} (n={n})")
        ax.legend(fontsize=8.5)
        ax.grid(alpha=0.3)

    fig.suptitle("Test SANS constante de calibration : le CL61 subit-il l'absorption vapeur d'eau ?",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = FIG / "wv_shape_test.png"
    fig.savefig(p, dpi=140)
    print(f"\nfigure -> {p}")


if __name__ == "__main__":
    main()
