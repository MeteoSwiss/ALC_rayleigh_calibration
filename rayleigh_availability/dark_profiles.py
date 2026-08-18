# -*- coding: utf-8 -*-
"""Profils de bruit d'obscurite (dark) mesures a Payerne, par instrument.

POURQUOI. La campagne dark de Payerne (mai-juillet 2026) couvre le telescope de chaque instrument
et mesure donc directement la LIGNE DE BASE electronique b(z) : ce que le detecteur lit en l'absence
totale de signal atmospherique. Cette ligne de base n'est pas nulle et, surtout, elle n'est pas
plate en portee -- or la restitution suppose implicitement qu'elle l'est (le signal est utilise tel
quel, `signal = rcs_mean / range^2`, sans soustraction d'un profil de fond).

Consequences mesurees par l'audit (doc/reports/dark_aeronet_sonde_audit.md) :
  * CL61   : b(z) CROIT en portee -> fabrique une pente. La soustraire retire +8,89 +- 0,57 %/km
             du gradient intra-nuit (15 nuits claires, t = 15,5), soit 66 % du gradient.
  * CHM15k : b(z) suit la decroissance moleculaire (forme en "cuvette") -> ne fabrique quasi aucune
             pente, mais un piedestal de -2,557e-4 +- 4,73e-5 counts/s qui biaise la constante de
             -17,2 % a 4,5 km. C'est une erreur d'ECHELLE, pas de forme.

PIEGE A NE PAS REPRODUIRE. `options.subtract_background` soustrait l'ORDONNEE A L'ORIGINE AJUSTEE de
la droite signal-vs-p_mol, pas une mesure. Le modele direct a montre qu'avec b_vrai = 0 et une simple
brume, cet intercept ajuste vaut deja -2,2e-4 : le soustraire revient a soustraire de l'atmosphere,
et cela AGGRAVE le gradient sur les vraies nuits. Ce module fournit au contraire un profil MESURE,
sous capot, ou aucune atmosphere ne peut entrer.

UNITES. b(z) est exprime dans les unites de `signal` du pipeline, c'est-a-dire rcs_0 / z^2, avec le
rcs_0 tel qu'il est stocke dans les L1 E-PROFILE (chaque type d'instrument a sa propre echelle : le
facteur vers beta_att est CHM15k ~3e-12, CL31 1e-8, CL61 1,0). On ne convertit donc PAS : le profil
est utilisable directement par soustraction sur `signal`.

CRENEAUX. Les fenetres du logbook sont ROGNEES de MARGIN_S de chaque cote pour ecarter la pose et la
depose du capot (le script MATLAB de reference fait de meme : 09:35-14:50 la ou le logbook dit
09:24-14:58 le 12 mai).

Run :  python rayleigh_availability/dark_profiles.py [--plot]
Sortie : <DATA>/dark_profiles_payerne.npz  (une entree par instrument : range, b, sem, n)
         + doc/reports/figs_altitude_audit/dark_profiles_measured.png avec --plot
"""
from __future__ import annotations
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
L1 = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610")
FIG = REPO / "doc" / "reports" / "figs_altitude_audit"
OUT = DATA / "dark_profiles_payerne.npz"

# Marge rognee de chaque cote d'un creneau (pose/depose du capot, stabilisation).
MARGIN_S = 360.0

# Lissage en PORTEE du profil poole. Indispensable : mesure sur mesure, le bruit porte-a-porte vaut
# 1,28x l'amplitude du profil CHM15k (0,57x pour le CL61) et la SEM par porte est du meme ordre que
# |b| -- soustraire le profil brut revient a injecter du bruit, ce qui fait echouer l'ajustement
# moleculaire. La ligne de base electronique est physiquement LISSE en portee (reponse du detecteur,
# queue d'afterpulse) : le script MATLAB de reference lisse de meme (medfilt2 [9 5] puis
# smoothdata "rlowess" 50). Fenetres exprimees en metres pour rester independantes de la grille.
MEDIAN_M = 150.0     # filtre median : tue les portes aberrantes isolees
MEAN_M = 450.0       # moyenne glissante : ne garde que la forme lentement variable

# Logbook officiel de la campagne (UTC). (debut, fin) par instrument et par date.
# ident -> liste de (debut, fin) en datetime naive UTC.
WINDOWS = {
    "A": [                                                     # CHM15k, carton ~5 mm
        ("2026-05-12 09:24", "2026-05-12 14:53"),
        ("2026-05-26 11:42", "2026-05-27 13:21"),              # cycle 24 h (dome de chaleur)
        ("2026-06-09 09:16", "2026-06-09 11:55"),
        ("2026-06-23 12:32", "2026-06-23 15:09"),
        ("2026-07-07 08:47", "2026-07-07 11:32"),
        ("2026-07-21 07:58", "2026-07-21 12:10"),
    ],
    "B": [                                                     # CL31, capot constructeur
        ("2026-05-12 09:33", "2026-05-12 14:55"),
        ("2026-05-26 11:35", "2026-05-27 13:16"),
        ("2026-06-09 09:16", "2026-06-09 11:57"),
        ("2026-06-23 10:12", "2026-06-23 12:13"),
        ("2026-07-07 08:46", "2026-07-07 11:19"),              # bloc optique remplace APRES, ~13:00
        ("2026-07-09 09:45", "2026-07-10 06:45"),
        ("2026-07-21 08:06", "2026-07-21 12:18"),
    ],
    "C": [                                                     # CL61, capot constructeur
        ("2026-05-12 09:17", "2026-05-12 14:58"),
        ("2026-05-26 11:37", "2026-05-27 13:26"),
        ("2026-06-09 09:10", "2026-06-09 11:51"),
        ("2026-06-23 10:00", "2026-06-23 12:25"),
        ("2026-07-07 11:35", "2026-07-07 14:04"),
        ("2026-07-21 07:54", "2026-07-21 12:06"),
    ],
}
# Le bloc optique du CL31 a ete remplace le 07/07/2026 ~13:00 : le gain change (x2,255 mesure sur
# C_daily). Les creneaux d'avant et d'apres ne decrivent donc PAS le meme instrument.
CL31_SWAP = datetime(2026, 7, 7, 13, 0)


def _parse(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M")


def _read_day(ident, day):
    """(times, range, signal=rcs_0/z^2) d'un fichier L1, ou None."""
    f = L1 / f"{day:%Y}" / f"{day:%m}" / f"L1_0-20000-0-06610_{ident}{day:%Y%m%d}.nc"
    if not f.exists():
        return None
    try:
        with netCDF4.Dataset(f) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)   # jours depuis 1970
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f8"), np.nan)
    except Exception as exc:
        print(f"    lecture impossible {f.name}: {exc}")
        return None
    if rcs.ndim != 2 or rcs.shape[1] != rng.size:
        return None
    times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
    # On travaille en unites rcs_0 (counts range-corriges), PAS en rcs_0/z^2 : c'est en counts que
    # la ligne de base electronique est une fonction LISSE de la portee. Divisee par z^2 elle
    # explose quand z -> 0 (le CL61 a une premiere porte a 0 m), le lissage en portee traine alors
    # ces valeurs enormes vers le haut, et comme le Klett integre l'epaisseur optique DEPUIS LE SOL
    # la corruption en proche portee contamine toutes les altitudes.
    rcs = np.where(rng[None, :] > 0, rcs, np.nan)
    return times, rng, rcs


def collect(ident, verbose=True):
    """Profils moyens par creneau + profil poole, pour un instrument."""
    per_window, rng_ref = [], None
    for w0, w1 in WINDOWS[ident]:
        t0, t1 = _parse(w0) + timedelta(seconds=MARGIN_S), _parse(w1) - timedelta(seconds=MARGIN_S)
        chunks = []
        day = t0.date()
        while day <= t1.date():                       # un creneau peut franchir minuit (26-27 mai)
            got = _read_day(ident, datetime.combine(day, datetime.min.time()))
            if got is not None:
                times, rng, sig = got
                if rng_ref is None:
                    rng_ref = rng
                if rng.size == rng_ref.size:
                    m = (times >= t0) & (times <= t1)
                    if m.any():
                        chunks.append(sig[m])
            day += timedelta(days=1)
        if not chunks:
            if verbose:
                print(f"    {w0[:16]} : aucun profil")
            continue
        block = np.vstack(chunks)
        # Mediane en temps : robuste aux rares pics (le MATLAB de reference filtre de meme).
        prof = np.nanmedian(block, axis=0)
        n = int(np.isfinite(block).sum(axis=0).max())
        per_window.append(dict(start=w0, n=n, profile=prof,
                               swap_before=(_parse(w0) < CL31_SWAP)))
        if verbose:
            band = (rng_ref >= 2500) & (rng_ref <= 6500)
            print(f"    {w0[:16]} : n={n:5d} profils, <b> 2,5-6,5 km = {np.nanmean(prof[band]):+.3e}")
    return rng_ref, per_window


def _smooth_range(rng, b):
    """Filtre median puis moyenne glissante, en fenetres exprimees en metres (voir MEDIAN_M/MEAN_M)."""
    dz = float(np.median(np.diff(rng))) if rng.size > 1 else 1.0
    if not np.isfinite(dz) or dz <= 0:
        return b
    k_med = max(3, int(round(MEDIAN_M / dz)) | 1)      # impair
    k_mean = max(3, int(round(MEAN_M / dz)) | 1)
    x = np.asarray(b, float).copy()
    good = np.isfinite(x)
    if good.sum() < k_mean:
        return x
    x[~good] = np.interp(rng[~good], rng[good], x[good])
    half = k_med // 2
    pad = np.pad(x, half, mode="edge")
    med = np.array([np.median(pad[i:i + k_med]) for i in range(x.size)])
    half = k_mean // 2
    pad = np.pad(med, half, mode="edge")
    ker = np.ones(k_mean) / k_mean
    out = np.convolve(pad, ker, mode="valid")
    out[~good] = np.nan                                 # ne pas inventer la ou rien n'a ete mesure
    return out


def pool(ident, rng, per_window):
    """Profil poole + erreur-type inter-creneaux (la dispersion honnete, cf. audit)."""
    if ident == "B":
        # CL31 : ne pooler QUE l'ere en cours (apres le remplacement du bloc optique), sinon on
        # moyenne deux instruments de gain different (facteur ~2,255 mesure).
        sel = [w for w in per_window if not w["swap_before"]]
        if len(sel) < 2:
            sel = [w for w in per_window if w["swap_before"]]
            print("    (CL31 : trop peu de creneaux post-swap, pool sur l'ere PRE-swap)")
    else:
        sel = per_window
    if not sel:
        return None
    P = np.vstack([_smooth_range(rng, w["profile"]) for w in sel])   # lisser AVANT de pooler
    b = np.nanmedian(P, axis=0)
    sem = np.nanstd(P, axis=0) / np.sqrt(max(P.shape[0], 1))
    return dict(range=rng, b=b, sem=sem, n_windows=P.shape[0],
                n_profiles=int(sum(w["n"] for w in sel)))


def main():
    out = {}
    for ident, name in (("A", "CHM15k"), ("B", "CL31"), ("C", "CL61")):
        print(f"\n== {name} (ident {ident})")
        rng, per_window = collect(ident)
        if rng is None or not per_window:
            print("   aucune donnee")
            continue
        pooled = pool(ident, rng, per_window)
        if pooled is None:
            continue
        band = (rng >= 2500) & (rng <= 6500)
        b = pooled["b"]
        # Pente de b(z) sur la bande de fit : c'est elle qui fabrique (ou non) un dC_L/dz.
        ok = band & np.isfinite(b)
        slope = np.polyfit(rng[ok] / 1000.0, b[ok], 1)[0] if ok.sum() > 10 else np.nan
        print(f"   poole sur {pooled['n_windows']} creneaux, {pooled['n_profiles']} profils")
        print(f"   <b> 2,5-6,5 km = {np.nanmean(b[band]):+.4e}   pente = {slope:+.3e} /km")
        for z in (1000, 2000, 3000, 4500, 6000, 8000):
            i = int(np.argmin(np.abs(rng - z)))
            print(f"     b({rng[i]:6.0f} m) = {b[i]:+.4e}")
        out[f"{ident}_range"] = rng
        out[f"{ident}_b_rcs"] = b                       # unites rcs_0 : ce que le pipeline soustrait
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"{ident}_b"] = np.where(rng > 0, b / rng ** 2, np.nan)   # unites `signal`, diagnostic
        out[f"{ident}_sem"] = pooled["sem"]
        out[f"{ident}_nwin"] = pooled["n_windows"]
        out[f"{ident}_nprof"] = pooled["n_profiles"]
    if out:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        np.savez(OUT, **out)
        print(f"\n-> {OUT}")
    if "--plot" in sys.argv and out:
        _figure(out)


def _figure(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = {"A": ("CHM15k (1064 nm)", "#2ca02c"), "B": ("CL31 (910 nm)", "#1f77b4"),
             "C": ("CL61 (910,55 nm)", "#d62728")}
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2))
    for ident, (lab, col) in names.items():
        if f"{ident}_b" not in out:
            continue
        rng, b, sem = out[f"{ident}_range"], out[f"{ident}_b"], out[f"{ident}_sem"]
        m = (rng > 100) & (rng <= 12000) & np.isfinite(b)
        # 1) profil brut, chaque instrument dans ses propres unites -> normalise pour la forme
        sc = np.nanmax(np.abs(b[m])) or 1.0
        axes[0].plot(b[m] / sc, rng[m] / 1000.0, "-", color=col, lw=1.5, label=lab)
        axes[0].fill_betweenx(rng[m] / 1000.0, (b[m] - sem[m]) / sc, (b[m] + sem[m]) / sc,
                              color=col, alpha=0.18, lw=0)
        # 2) le meme, en % du signal moleculaire typique a cette altitude
        axes[1].plot(b[m] / sc, rng[m] / 1000.0, "-", color=col, lw=1.5)
    axes[0].axvline(0, color="#444", ls=":", lw=1)
    axes[0].set_xlabel("b(z) normalise par son maximum")
    axes[0].set_ylabel("portee [km]")
    axes[0].set_title("Ligne de base mesuree sous capot")
    axes[0].legend(fontsize=8.5)
    axes[0].grid(alpha=0.3)
    axes[1].axvline(0, color="#444", ls=":", lw=1)
    axes[1].set_ylim(0, 3)
    axes[1].set_xlabel("b(z) normalise par son maximum")
    axes[1].set_ylabel("portee [km]")
    axes[1].set_title("Zoom courte portee (0-3 km)")
    axes[1].grid(alpha=0.3)
    # 3) pente locale de b(z) : c'est elle qui fabrique un gradient de C_L
    for ident, (lab, col) in names.items():
        if f"{ident}_b" not in out:
            continue
        rng, b = out[f"{ident}_range"], out[f"{ident}_b"]
        m = (rng > 500) & (rng <= 10000) & np.isfinite(b)
        z, bb = rng[m] / 1000.0, b[m]
        sc = np.nanmax(np.abs(bb)) or 1.0
        w = 40
        d = np.full(z.size, np.nan)
        for k in range(w, z.size - w):
            d[k] = np.polyfit(z[k - w:k + w], bb[k - w:k + w], 1)[0] / sc
        axes[2].plot(d, z, "-", color=col, lw=1.5, label=lab)
    axes[2].axvline(0, color="#444", ls=":", lw=1)
    axes[2].set_xlabel("db/dz local (normalise) [1/km]")
    axes[2].set_ylabel("portee [km]")
    axes[2].set_title("Pente de la ligne de base\n(>0 = fabrique un gradient de $C_L$)")
    axes[2].legend(fontsize=8.5)
    axes[2].grid(alpha=0.3)
    fig.suptitle("Profils de bruit d'obscurite mesures — campagne Payerne mai–juillet 2026",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / "dark_profiles_measured.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
