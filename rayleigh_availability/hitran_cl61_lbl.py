# -*- coding: utf-8 -*-
"""Simulation HITRAN line-by-line de l'absorption WV vue par le CL61.

POURQUOI. Notre LUT operationnelle (`abs_cross_wv_910nm.nc`) est echantillonnee a 8,3 pm en
MOYENNE PAR BIN, alors que les raies elargies par pression font ~17 pm FWHM au sol : elle ne
peut pas resoudre les fenetres inter-raies etroites. Or Vaisala annonce que l'absorption WV du
CL61 est « mitigated by selecting a different wavelength and a very narrow bandwidth » (User
Guide M212475EN-E, relaye par le refere RC1 du preprint Le & O'Connor), et nos scenes reelles
montrent que son signal brut ne porte que ~8 % de l'absorption que notre modele lui applique.
Ce script tranche avec un calcul line-by-line (HITRAN via hapi, profils de Voigt) a resolution
fine, pour trois familles de longueurs d'onde :

  * 910,55 nm  : valeur NOMINALE Vaisala (spec constructeur)
  * 910,74 nm  : valeur MESUREE au spectrometre Qmini a Payerne (+-0,10 nm, annexe A1)
  * derive thermique : lambda0 +- (0,05..0,30 nm/K x Delta T) autour des deux precedentes

et pour des largeurs d'emission de 0,01 a 1,5 nm (la mesure borne le FWHM vrai a < 1,5 nm,
non resolu par le spectrometre).

SORTIE. Pour chaque (lambda0, FWHM) : la section efficace effective, la transmission
deux-voies a la CBH mediane des scenes, et surtout la PENTE INJECTEE dC/dCBH (%/km) -- la
grandeur directement comparable au +9,4 %/km observe et au niveau ~8 % qu'un signal brut
« mitige par conception » exhiberait.

Run :  python rayleigh_availability/hitran_cl61_lbl.py [--dnu 0.0005] [--workers 8]
       (premiere execution : telechargement des raies HITRAN, ~qq Mo, mis en cache)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
CACHE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/hitran_cache")
FIG = REPO / "doc" / "reports" / "figs_cbh_heterogeneity"

# Fenetre spectrale : large pour couvrir toute derive thermique plausible
LAM_MIN, LAM_MAX = 906.0, 914.0          # nm
NU_MIN, NU_MAX = 1e7 / LAM_MAX, 1e7 / LAM_MIN     # cm-1  (~10941 .. 11038)

# Longueurs d'onde testees
LAM_NOMINAL = 910.55        # spec Vaisala
LAM_MESUREE = 910.74        # Qmini, annexe A1 (+-0,10)
DRIFTS_NM = (-0.60, -0.30, -0.15, 0.0, +0.15, +0.30, +0.60)
FWHMS_NM = (0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 1.50)

# Atmosphere de reference (Payerne, ete moyen) : niveaux du sol a 3 km, la zone qui compte
# pour l'absorption sous-nuage (CBH des scenes 0,5-2,4 km).
Z_KM = np.array([0.49, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0])
T_K = 288.0 - 6.5 * (Z_KM - 0.49)                       # gradient standard
P_ATM = np.exp(-(Z_KM - 0.49) / 8.4) * 0.955            # atm (Payerne ~955 hPa au sol)
# vapeur d'eau : echelle de hauteur 2 km, 8 g/kg au sol -> densite numerique [molec/cm3]
RH_SCALE_KM = 2.0
N_WV_SURF = 4.0e17                                       # molec/cm3 (~12 g/m3)
N_WV = N_WV_SURF * np.exp(-(Z_KM - 0.49) / RH_SCALE_KM)


def fetch_lines(dnu_cm: float):
    """Sections efficaces HITRAN par niveau, sur une grille fine. Cache disque."""
    import hapi
    CACHE.mkdir(parents=True, exist_ok=True)
    hapi.db_begin(str(CACHE))
    tbl = "H2O_910"
    if tbl not in hapi.getTableList():
        print(f"telechargement HITRAN H2O {NU_MIN:.0f}-{NU_MAX:.0f} cm-1 ...")
        hapi.fetch_by_ids(tbl, [1, 2, 3, 4, 5, 6, 7], NU_MIN - 5, NU_MAX + 5)
    n_lines = len(hapi.getColumn(tbl, "nu"))
    print(f"raies H2O disponibles : {n_lines}")

    sigmas = []
    for k, (p, t) in enumerate(zip(P_ATM, T_K)):
        nu, coef = hapi.absorptionCoefficient_Voigt(
            SourceTables=tbl, HITRAN_units=True,            # cm2/molecule
            Environment={"p": float(p), "T": float(t)},
            WavenumberRange=[NU_MIN, NU_MAX], WavenumberStep=dnu_cm,
            Diluent={"air": 1.0},
        )
        sigmas.append(np.asarray(coef, float))
        print(f"  niveau {k+1}/{len(P_ATM)}  z={Z_KM[k]:.2f} km  p={p:.3f} atm  T={t:.1f} K  "
              f"sigma max={np.max(coef):.2e} cm2")
    nu = np.asarray(nu, float)
    lam = 1e7 / nu                                          # nm, decroissant
    order = np.argsort(lam)
    return lam[order], np.vstack(sigmas)[:, order]


def sigma_eff(lam, sig_lvl, lam0, fwhm):
    """Section efficace effective par niveau pour une raie gaussienne (lam0, fwhm)."""
    s = fwhm / 2.3548
    g = np.exp(-0.5 * ((lam - lam0) / s) ** 2)
    g /= g.sum()
    return sig_lvl @ g


def cbh_slope(lam, sig_lvl):
    """(T2 a 1,25 km, pente dC/dCBH injectee sur 0,5-2,4 km) pour une famille (lam0, fwhm)."""
    def t2_of(sig_eff):
        # profondeur optique une-voie cumulee du sol a z (integration trapeze en cm)
        dz_cm = np.diff(Z_KM) * 1e5
        tau = np.concatenate([[0.0], np.cumsum(
            0.5 * (sig_eff[1:] * N_WV[1:] + sig_eff[:-1] * N_WV[:-1]) * dz_cm)])
        return np.exp(-2.0 * tau)
    return t2_of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dnu", type=float, default=0.0005,
                    help="pas spectral en cm-1 (0.0005 ~ 0.04 pm a 910 nm)")
    args = ap.parse_args()

    lam, sig = fetch_lines(args.dnu)
    dlam_pm = 1000 * np.median(np.diff(lam))
    print(f"grille : {lam.size} points, pas median {dlam_pm:.3f} pm "
          f"({lam.min():.2f}-{lam.max():.2f} nm)")

    t2_of = cbh_slope(lam, sig)
    z_fit = (Z_KM >= 0.5) & (Z_KM <= 2.4)

    def metrics(lam0, fwhm):
        se = sigma_eff(lam, sig, lam0, fwhm)
        t2 = t2_of(se)
        ln_t2 = np.log(np.maximum(t2, 1e-300))
        # la correction ajoute -ln T2 a ln C : pente = -d(ln T2)/dz
        A = np.vstack([Z_KM[z_fit], np.ones(z_fit.sum())]).T
        slope = -100.0 * np.linalg.lstsq(A, ln_t2[z_fit], rcond=None)[0][0]
        i125 = int(np.argmin(np.abs(Z_KM - 1.25)))
        return se[0], t2[i125], slope

    # reference : la meme metrique avec NOTRE LUT (8,3 pm) pour le meme couple
    from calibration.water_vapor_correction.water_vapor import load_abs_cross_section
    wl_lut, h_lut, cs_lut = load_abs_cross_section(
        REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc")
    lut_lvl = np.vstack([cs_lut[:, int(np.argmin(np.abs(h_lut - z * 1000)))] for z in Z_KM])

    def metrics_lut(lam0, fwhm):
        s = fwhm / 2.3548
        g = np.exp(-0.5 * ((wl_lut - lam0) / s) ** 2)
        g /= g.sum()
        se = lut_lvl @ g
        t2 = t2_of(se)
        ln_t2 = np.log(np.maximum(t2, 1e-300))
        A = np.vstack([Z_KM[z_fit], np.ones(z_fit.sum())]).T
        return -100.0 * np.linalg.lstsq(A, ln_t2[z_fit], rcond=None)[0][0]

    out = {"grid_pm": float(dlam_pm), "cases": []}
    print("\n=== pente dC/dCBH injectee [%/km] : HITRAN line-by-line vs notre LUT (8,3 pm) ===")
    print(f"{'lambda0':>9s} {'FWHM':>6s} | {'sigma_eff sol':>14s} {'T2(1.25km)':>11s} "
          f"{'pente LBL':>10s} {'pente LUT':>10s} {'LBL/LUT':>8s}")
    for base, tag in ((LAM_NOMINAL, "nominal"), (LAM_MESUREE, "mesure")):
        for d in DRIFTS_NM:
            lam0 = base + d
            for fw in FWHMS_NM:
                se0, t2, sl = metrics(lam0, fw)
                sl_lut = metrics_lut(lam0, fw)
                ratio = sl / sl_lut if abs(sl_lut) > 1e-9 else np.nan
                out["cases"].append(dict(base=tag, lam0=lam0, drift=d, fwhm=fw,
                                         sigma_surface=float(se0), t2_125=float(t2),
                                         slope_lbl=float(sl), slope_lut=float(sl_lut),
                                         ratio=float(ratio)))
                if d in (0.0, -0.30, +0.30) and fw in (0.02, 0.10, 0.50, 1.00):
                    print(f"{lam0:9.3f} {fw:6.2f} | {se0:14.3e} {t2:11.4f} "
                          f"{sl:10.2f} {sl_lut:10.2f} {ratio:8.2f}")

    (CACHE / "hitran_cl61_results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"\n-> {CACHE / 'hitran_cl61_results.json'}")

    # ---------------- figures (paysage) ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.0))
    m = (lam >= 909.8) & (lam <= 911.5)
    ax = axes[0]
    ax.semilogy(lam[m], sig[0][m], lw=0.5, color="#1f77b4",
                label=f"HITRAN line-by-line, sol ({dlam_pm:.2f} pm)")
    ml = (wl_lut >= 909.8) & (wl_lut <= 911.5)
    ax.semilogy(wl_lut[ml], lut_lvl[0][ml], lw=1.4, color="#e8871a",
                label="notre LUT (moyenne par bin de 8,3 pm)")
    ax.axvline(LAM_NOMINAL, color="#2ca02c", lw=1.6, ls="--", label=f"nominal {LAM_NOMINAL}")
    ax.axvline(LAM_MESUREE, color="#d62728", lw=1.6, label=f"mesure {LAM_MESUREE}")
    ax.set_xlabel("longueur d'onde [nm]")
    ax.set_ylabel(r"section efficace H$_2$O [cm$^2$/molecule]")
    ax.set_title("Le spectre vrai vs celui que voit notre LUT")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    for fw, c in zip((0.02, 0.10, 0.50, 1.00), ("#d62728", "#e8871a", "#2ca02c", "#1f77b4")):
        lams = np.linspace(910.0, 911.3, 261)
        sl = [metrics(l0, fw)[2] for l0 in lams]
        ax.plot(lams, sl, lw=1.5, color=c, label=f"FWHM {fw:.2f} nm")
    ax.axvline(LAM_NOMINAL, color="#2ca02c", lw=1.4, ls="--")
    ax.axvline(LAM_MESUREE, color="#d62728", lw=1.4)
    ax.axhline(0.8, color="#444", ls=":", lw=1.4,
               label="niveau vu par le signal brut CL61 (~8 %)")
    ax.set_xlabel(r"$\lambda_0$ [nm]")
    ax.set_ylabel("pente dC/dCBH injectee [%/km]")
    ax.set_title("Balayage fin de $\\lambda_0$ (line-by-line) :\n"
                 "existe-t-il une fenetre qui explique le deficit du CL61 ?")
    ax.legend(fontsize=8.5)
    ax.grid(alpha=0.3)
    fig.suptitle("CL61 : absorption vapeur d'eau line-by-line (HITRAN) — "
                 "nominal 910,55 vs mesure 910,74 nm", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / "hitran_cl61_lbl.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
