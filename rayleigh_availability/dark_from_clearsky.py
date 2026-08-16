# -*- coding: utf-8 -*-
"""Estimation de la ligne de base electronique b(z) depuis les nuits claires (sans capot).

POURQUOI. La campagne capot de Payerne (dark_profiles.py) mesure b(z) directement, mais seul
Payerne a un capot. L'extraction nuit-claire historique (rapport 08 par 3.8) ne recupere que la
composante PERIODIQUE (ripple du numeriseur, gate-fold) : elle est structurellement aveugle aux
composantes LISSES -- le piedestal du CHM15k (-17 % du moleculaire a 4,5 km) et la ligne de base
croissante du CL61 -- precisement celles qui biaisent le fit Rayleigh. Ce module tente de les
estimer quand meme, en exploitant la seule information qui distingue le moleculaire de
l'electronique : LE MOLECULAIRE VARIE D'UNE NUIT A L'AUTRE d'une facon CONNUE (CAMS T/p + vapeur
d'eau), l'electronique non.

MODELE. Par nuit claire n et porte z, en unites rcs_0 (counts range-corriges des L1) :

    S_n(z) = A_n * m_n(z) + b(z) + d_n * z^2 + bruit

  * m_n(z)   : backscatter moleculaire attenue calcule de CAMS pour la nuit n (Bucholtz + T2_mol,
               x T2_wv pour les instruments a 910 nm). Forme quasi exponentielle, CONNUE.
  * A_n      : amplitude libre par nuit = C * T2_aerosol (transmission aerosol + derive du laser).
  * b(z)     : ligne de base electronique STATIQUE -- la cible. Libre par porte, lissee a la fin
               avec les memes fenetres que la mesure capot (median 150 m + moyenne 450 m).
  * d_n*z^2  : residu de la soustraction de fond du firmware. Un fond residuel est plat en counts
               BRUTS, donc croit en z^2 dans rcs_0 ; il VARIE par nuit (fond de ciel, lune,
               temperature). Sa part STATIQUE est reversee dans b (contrainte mediane(d_n)=0),
               comme dans la mesure capot qui la contient aussi.

Les trois regresseurs ont des formes tres differentes (exponentielle decroissante, libre lisse,
parabole croissante) : le fit alterne robuste (ALS) les separe si A_n varie assez d'une nuit a
l'autre. La degenerescence residuelle (une composante de b parallele a la forme moleculaire
moyenne) est QUANTIFIEE par le self-test d'injection (--selftest), pas supposee nulle.

PORTEE HONNETE. Le fit travaille au-dessus de Z0 = 2 km : en dessous, aerosol de couche limite et
overlap sont indissociables d'un b lisse sans capot. Le npz produit vaut donc pour la zone des
fenetres Rayleigh (2-8 km) ; sous 2 km il vaut NaN -> le pipeline (_dark_on_grid) y applique 0.
La mesure capot reste la seule verite en proche portee (rapport 08 par 3.9).

VALIDATION. A Payerne, b_est est confronte au b mesure sous capot (dark_profiles_payerne.npz)
pour les trois instruments -- c'est le critere de robustesse demande avant tout deploiement
reseau. Diagnostics par flux : reproductibilite split-half (nuits paires/impaires), SEM bootstrap,
et biais Rayleigh predit (b integre sur les fenetres reellement ajustees, en % du moleculaire).

Run (un flux) :
    python rayleigh_availability/dark_from_clearsky.py --wmo 0-20000-0-06610 --ident A \
        --start 20250601 [--selftest] [--plot]
Sortie : <DATA>/dark_clearsky/<wmo>_<ident>.npz  (+ <...>_diag.json, + figure avec --plot)
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration.io.cams import find_cams_file                              # noqa: E402
from calibration.omb.cams_aerosol import cams_aerosol_backscatter           # noqa: E402
from calibration.rayleigh.atmosphere import calculate_molecular_properties  # noqa: E402
from calibration.water_vapor_correction.water_vapor import (                # noqa: E402
    cams_temperature_pressure_profile,
    cams_water_vapor_profile,
    in_water_vapor_band,
    laser_spectrum_for,
    two_way_wv_transmission,
)

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
OUT_DIR = DATA / "dark_clearsky"
DEFAULT_L1 = Path("D:/E-PROFILE_L1_2026")
DEFAULT_CSV = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04")
DEFAULT_CAMS = "A:/CAMS_Monthly_04;D:/CAMS_daily"
WV_LUT = REPO / "calibration" / "data" / "abs_cross_wv_910nm.nc"

# Nuits retenues : tout flag ou l'ecran nuit-claire est passe (le fit a ete TENTE). -1 (pas une
# nuit claire), 0 (pas de profils), -4 (pas de CAMS), -10/-11 sont exclus.
CLEAR_FLAGS = {1.0, 0.5, -2.0, -3.0, -9.0}

# Fenetre du fit en portee. Sous Z0 l'aerosol de couche limite et l'overlap dominent ; au-dessus
# de Z1 beaucoup de L1 tronquent ou saturent en bruit.
Z0_M = 2000.0
Z1_M = 15000.0

# Soleil sous -8 deg = nuit exploitable (meme esprit que l'ecran clear-night du pipeline) ;
# segment du MATIN du fichier L1 du jour d (avant 14 h UTC) pour rester dans la nuit que la ligne
# CSV du jour d qualifie de claire.
SUN_ELEV_MAX = -8.0
MORNING_CUT_H = 14.0
MIN_PROFILES = 100

# Lissage en portee du b final : memes fenetres que la mesure capot (dark_profiles.py), pour que
# les deux profils soient comparables terme a terme.
MEDIAN_M = 150.0
MEAN_M = 450.0

# Lissage de la reponse multiplicative g(z) : au-dessus de l'overlap complet, la PHYSIQUE impose
# une reponse quasi constante -- toute structure plus fine que ~2,5 km dans g est de l'atmosphere
# (desaccord aerosol CAMS/reel) qu'on refuse de laisser fuir dans le canal multiplicatif.
G_SMOOTH_M = 2500.0

# Ecran physique par nuit : A_n = C * T2_aerosol, donc A_n / C_ref doit rester dans l'enveloppe
# [T2_aerosol plausible en nuit claire] x [derive saisonniere de C]. A 1064 nm T2_aer ~ 0.6-1.0,
# et C varie de +-30 % dans l'annee (cycle mesure a Payerne) ; C_ref est une MEDIANE, donc les
# pics d'hiver legitimes montent a ~1.3-1.4. Une nuit hors [0.35, 1.45] viole le modele (nuage
# residuel, aerosol aliase dans la pente) et est ecartee -- au premier reglage [0.25, 1.15]
# l'ecran coupait les pics d'hiver et tuait le levier (CV(A) 0.31, lambda 0.26).
T2_MIN, T2_MAX = 0.35, 1.45

# Bande de la projection theta : sous 3 km l'aerosol de couche limite contamine l'ordonnee,
# au-dessus de ~9,5 km les erreurs de d_n (x z^2) dominent.
THETA_Z = (3000.0, 9500.0)

N_ALS_ITER = 12
ALS_TOL = 1e-4


# ---------------------------------------------------------------------------
# Petites briques
# ---------------------------------------------------------------------------

def solar_elevation_deg(lat: float, lon: float, times: np.ndarray) -> np.ndarray:
    """Elevation solaire (deg) pour des datetimes UTC naifs -- formule NOAA approchee (~0,3 deg),
    largement suffisante pour un seuil de nuit a -8 deg."""
    j2000 = datetime(2000, 1, 1, 12)
    ndays = np.array([(t - j2000).total_seconds() / 86400.0 for t in times])
    g = np.radians((357.528 + 0.9856003 * ndays) % 360.0)
    lam = np.radians((280.460 + 0.9856474 * ndays) % 360.0
                     + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    eps = np.radians(23.439 - 4.0e-7 * ndays)
    delta = np.arcsin(np.sin(eps) * np.sin(lam))
    alpha = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    ut_h = np.array([t.hour + t.minute / 60.0 + t.second / 3600.0 for t in times])
    gmst_h = (6.697375 + 0.0657098242 * ndays + 1.0027379 * ut_h) % 24.0
    ha = np.radians(((gmst_h + lon / 15.0) * 15.0) % 360.0) - alpha
    sin_el = (np.sin(np.radians(lat)) * np.sin(delta)
              + np.cos(np.radians(lat)) * np.cos(delta) * np.cos(ha))
    return np.degrees(np.arcsin(np.clip(sin_el, -1, 1)))


def smooth_range(rng: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Median 150 m puis moyenne glissante 450 m -- copie de dark_profiles._smooth_range pour que
    l'estime et la mesure capot subissent exactement le meme traitement."""
    dz = float(np.median(np.diff(rng))) if rng.size > 1 else 1.0
    if not np.isfinite(dz) or dz <= 0:
        return b
    k_med = max(3, int(round(MEDIAN_M / dz)) | 1)
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
    out = np.convolve(pad, np.ones(k_mean) / k_mean, mode="valid")
    out[~good] = np.nan
    return out


def clear_nights_from_csv(csv_root: Path, wmo: str, ident: str,
                          start: str, end: str) -> list[str]:
    """Dates (YYYYMMDD) des nuits que le run reseau a qualifiees de claires (fit tente)."""
    key = f"{wmo}_{ident}"
    csv = csv_root / key / f"{key}_cal.csv"
    if not csv.exists():
        raise FileNotFoundError(f"CSV reseau introuvable : {csv}")
    dates = []
    with open(csv, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip().split(",")
        i_date, i_meth, i_flag = (header.index(c) for c in ("date", "method", "flag"))
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) <= max(i_date, i_meth, i_flag):
                continue
            if parts[i_meth] != "rayleigh":
                continue
            try:
                flag = float(parts[i_flag])
            except ValueError:
                continue
            d = parts[i_date]
            if flag in CLEAR_FLAGS and start <= d <= end:
                dates.append(d)
    return sorted(set(dates))


def fitted_windows_from_csv(csv_root: Path, wmo: str, ident: str) -> list[tuple[str, float, float]]:
    """(date, bottom_m, top_m) des nuits Rayleigh REUSSIES -- pour le biais predit."""
    key = f"{wmo}_{ident}"
    csv = csv_root / key / f"{key}_cal.csv"
    rows = []
    with open(csv, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip().split(",")
        idx = {c: header.index(c) for c in ("date", "method", "flag", "bottom_height", "top_height")}
        for line in fh:
            p = line.rstrip("\n").split(",")
            if len(p) <= idx["top_height"] or p[idx["method"]] != "rayleigh":
                continue
            try:
                if float(p[idx["flag"]]) not in (1.0, 0.5):
                    continue
                rows.append((p[idx["date"]], float(p[idx["bottom_height"]]), float(p[idx["top_height"]])))
            except ValueError:
                continue
    return rows


# ---------------------------------------------------------------------------
# Une nuit -> (S, sigma, m)   [execute dans les workers]
# ---------------------------------------------------------------------------

def process_night(job: dict):
    """Lit le L1 du jour d, garde le segment de nuit du matin, mediane par porte, et construit le
    predicteur moleculaire CAMS de la nuit. Retourne None si la nuit est inutilisable."""
    import netCDF4  # import dans le worker (spawn Windows)

    date, wmo, ident = job["date"], job["wmo"], job["ident"]
    l1 = (Path(job["l1_root"]) / wmo / date[:4] / date[4:6]
          / f"L1_{wmo}_{ident}{date}.nc")
    if not l1.exists():
        return None
    try:
        with netCDF4.Dataset(l1) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][:].astype("f8"), np.nan)
            lat = float(np.ma.filled(ds.variables["station_latitude"][:], np.nan).flat[0])
            lon = float(np.ma.filled(ds.variables["station_longitude"][:], np.nan).flat[0])
            alt = float(np.ma.filled(ds.variables["station_altitude"][:], np.nan).flat[0])
            try:
                wl_nm = float(np.ma.filled(ds.variables["l0_wavelength"][:], np.nan).flat[0])
            except Exception:  # noqa: BLE001
                wl_nm = np.nan
            itype = str(getattr(ds, "instrument_type", "")).strip()
    except Exception:  # noqa: BLE001 -- fichier corrompu -> nuit sautee
        return None
    if rcs.ndim != 2 or rcs.shape[1] != rng.size or not np.isfinite(lat + lon + alt):
        return None
    if not np.isfinite(wl_nm) or wl_nm <= 0:
        wl_nm = {"CHM15k": 1064.0, "CL61": 910.55, "CL31": 910.0, "CL51": 910.0}.get(itype, np.nan)
        if not np.isfinite(wl_nm):
            return None

    times = np.array([datetime(1970, 1, 1) + timedelta(days=float(x)) for x in t])
    hours = np.array([tt.hour + tt.minute / 60.0 for tt in times])
    elev = solar_elevation_deg(lat, lon, times)
    keep = (elev < SUN_ELEV_MAX) & (hours < MORNING_CUT_H)
    if keep.sum() < MIN_PROFILES:
        return None
    block = rcs[keep]
    signal = np.nanmedian(block, axis=0)
    n_eff = np.maximum(np.isfinite(block).sum(axis=0), 1)
    sigma = 1.4826 * np.nanmedian(np.abs(block - signal[None, :]), axis=0) / np.sqrt(n_eff)
    t0 = np.datetime64(times[keep][0].strftime("%Y-%m-%dT%H:%M:%S"))
    t1 = np.datetime64(times[keep][-1].strftime("%Y-%m-%dT%H:%M:%S"))

    cams = find_cams_file(job["cams_folder"], date)
    if cams is None:
        return None
    prof = cams_temperature_pressure_profile(cams, lat, lon, t0, t1)
    if prof is None:
        return None
    h_c, t_c, p_c = prof
    alt_grid = alt + rng
    temp = np.interp(alt_grid, h_c, t_c, left=np.nan, right=np.nan)
    pres = np.interp(alt_grid, h_c, p_c, left=np.nan, right=np.nan)
    ok = np.where(np.isfinite(temp))[0]
    if ok.size == 0:
        return None
    temp[: ok[0]] = temp[ok[0]]
    pres[: ok[0]] = pres[ok[0]]

    mol = calculate_molecular_properties(temp, pres, rng, wl_nm * 1e-9)
    m = mol.beta_att_mol
    if in_water_vapor_band(wl_nm):
        lam0, fwhm = laser_spectrum_for(itype, wl_nm)
        wv = cams_water_vapor_profile(cams, lat, lon, t0, t1)
        if wv is None:
            return None
        t2 = two_way_wv_transmission(alt_grid, alt, wv[0], wv[1], WV_LUT, lam0, fwhm)
        m = m * t2

    # Regresseur aerosol CAMS (prevision assimilee AOD) : absorbe l'aerosol persistant de la
    # troposphere libre qui sinon se projette dans b via la direction quasi degeneree ~m(z).
    # Le CAMS beta_aer est "attenuated from ground" -- la meme forme que ce que voit le lidar.
    a = np.zeros(rng.size)
    try:
        t_num, z_mod, b_aer = cams_aerosol_backscatter(str(cams), lat, lon, wl_nm)
        # time_num est un datenum MATLAB ; fenetre de nuit elargie aux pas 3 h adjacents
        dn0 = (t0 - np.datetime64("1970-01-01")) / np.timedelta64(1, "D") + 719529.0 - 0.09
        dn1 = (t1 - np.datetime64("1970-01-01")) / np.timedelta64(1, "D") + 719529.0 + 0.09
        sel = (t_num >= dn0) & (t_num <= dn1)
        if sel.any():
            zcol = np.nanmean(z_mod[:, sel], axis=1)
            acol = np.nanmean(b_aer[:, sel], axis=1)
            order = np.argsort(zcol)
            good = np.isfinite(zcol[order]) & np.isfinite(acol[order])
            if good.sum() > 5:
                a = np.interp(alt_grid, zcol[order][good], acol[order][good],
                              left=np.nan, right=0.0)
                first = np.where(np.isfinite(a))[0]
                if first.size:
                    a[: first[0]] = a[first[0]]
                a = np.nan_to_num(a, nan=0.0)
    except Exception:  # noqa: BLE001 -- pas d'aerosol CAMS -> regresseur nul (fit a 2 termes)
        a = np.zeros(rng.size)

    return dict(date=date, rng=rng, S=signal, sigma=sigma, m=np.asarray(m, float),
                a=np.asarray(a, float),
                n_prof=int(keep.sum()), lat=lat, lon=lon, alt=alt, wl_nm=wl_nm, itype=itype)


# ---------------------------------------------------------------------------
# Le fit alterne robuste
# ---------------------------------------------------------------------------

def _wls_night(y, regs, w, min_gates=30):
    """WLS robuste (2 passes, ecretage 4xMAD) de y sur les colonnes de regs (liste de 1D).
    Retourne les coefficients ou None."""
    good0 = np.isfinite(y) & np.isfinite(w)
    for r in regs:
        good0 &= np.isfinite(r)
    if good0.sum() < min_gates:
        return None
    gm = good0.copy()
    coef = None
    for _pass in range(2):
        X = np.column_stack([r[gm] for r in regs])
        sw = np.sqrt(w[gm])
        try:
            coef, *_ = np.linalg.lstsq(X * sw[:, None], y[gm] * sw, rcond=None)
        except np.linalg.LinAlgError:
            return None
        resid = y - sum(c * r for c, r in zip(coef, regs))
        mad = 1.4826 * np.nanmedian(np.abs(resid[gm]))
        if not np.isfinite(mad) or mad <= 0:
            break
        gm = good0 & (np.abs(resid) < 4 * mad)
        if gm.sum() < min_gates:
            break
    return coef


def als_fit(S: np.ndarray, sig: np.ndarray, M: np.ndarray, Aer: np.ndarray,
            rng: np.ndarray, fit_mask: np.ndarray, n_iter: int = N_ALS_ITER):
    """S, sig, M, Aer : (n_nuits, n_portes). Retourne A, B, d (n) et b (n_portes).

    Alternance : (1) par nuit, WLS robuste sur S_n - b avec regresseurs [m_n, a_n, z^2] --
    a_n est l'aerosol CAMS de la nuit, contraint B_n >= 0 (un chargement negatif signifierait
    utiliser la forme aerosol pour ajuster l'electronique : refit sans aerosol) ;
    (2) b(z) = mediane_nuits du residu ; la part statique du terme z^2 est reversee dans b
    (la mesure capot la contient aussi), la part statique de l'aerosol reste ATMOSPHERE
    (elle est portee par B_n * a_n, jamais par b).
    """
    z2 = (rng / 1000.0) ** 2                       # km^2 : conditionne mieux le systeme
    n_nights, n_gates = S.shape
    b = np.zeros(n_gates)
    A = np.full(n_nights, np.nan)
    B = np.zeros(n_nights)
    d = np.zeros(n_nights)
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    scale = np.nanmedian(np.abs(S[:, fit_mask])) or 1.0

    for it in range(n_iter):
        for n in range(n_nights):
            y0 = np.where(fit_mask, S[n] - b, np.nan)
            has_aer = np.nanmax(np.abs(Aer[n])) > 0
            coef = _wls_night(y0, [M[n], Aer[n], z2], w_all[n]) if has_aer else None
            if coef is not None and coef[1] < 0:
                coef = None                        # chargement aerosol negatif -> sans aerosol
            if coef is None:
                c2 = _wls_night(y0, [M[n], z2], w_all[n])
                if c2 is None:
                    A[n] = np.nan
                    continue
                A[n], d[n] = c2
                B[n] = 0.0
            else:
                A[n], B[n], d[n] = coef
        ok_n = np.isfinite(A) & (A > 0)
        if ok_n.sum() < 5:
            return A, B, d, np.full(n_gates, np.nan), ok_n
        # la part statique du terme z^2 appartient a b (la mesure capot la contient aussi)
        d_med = np.median(d[ok_n])
        d[ok_n] -= d_med
        R = (S[ok_n] - A[ok_n, None] * M[ok_n] - B[ok_n, None] * Aer[ok_n]
             - d[ok_n, None] * z2[None, :])
        b_new = np.where(fit_mask, np.nanmedian(R, axis=0) + d_med * z2, np.nan)
        delta = np.nanmax(np.abs(np.nan_to_num(b_new) - np.nan_to_num(b))) / scale
        b = b_new
        if it >= 2 and delta < ALS_TOL:
            break
    return A, B, d, b, ok_n


def read_cref(csv_root: Path, wmo: str, ident: str) -> float:
    """C de reference du flux (mediane des nuits Rayleigh reussies) pour l'ecran physique
    T2_aer = A_n / C_ref. NaN si aucun succes (l'ecran est alors desactive)."""
    key = f"{wmo}_{ident}"
    csv = csv_root / key / f"{key}_cal.csv"
    vals = []
    try:
        with open(csv, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().strip().split(",")
            idx = {c: header.index(c) for c in ("method", "flag", "cal_value")}
            for line in fh:
                p = line.rstrip("\n").split(",")
                if len(p) <= idx["cal_value"] or p[idx["method"]] != "rayleigh":
                    continue
                try:
                    if float(p[idx["flag"]]) in (1.0, 0.5):
                        vals.append(float(p[idx["cal_value"]]))
                except ValueError:
                    continue
    except OSError:
        return np.nan
    return float(np.median(vals)) if len(vals) >= 5 else np.nan


def load_template(itype: str, rng: np.ndarray, fit_mask: np.ndarray) -> np.ndarray | None:
    """Forme de dark par TYPE d'instrument, tiree de la campagne capot de Payerne
    (A=CHM15k, B=CL31, C=CL61). La physique de la ligne de base est materielle et commune au
    type (rapport 08 par 3.3) : on estime alors UNE amplitude theta par unite au lieu d'un profil
    libre -- c'est la parade standard a l'identifiabilite partielle. theta = 1 signifie "meme
    dark que l'unite Payerne du meme type" (les unites rcs_0 d'un type sont homogenes reseau).

    La forme est ORTHOGONALISEE contre {z^2, 1} sur la bande de fit : la part z^2-statique du
    dark est portee par d_med (bien identifiee nuit par nuit), pas par theta -- sinon les deux
    se disputent la meme direction et l'ALS oscille."""
    ident = {"CHM15k": "A", "CL31": "B", "CL61": "C"}.get(itype)
    npz = DATA / "dark_profiles_payerne.npz"
    if ident is None or not npz.exists():
        return None
    z = np.load(npz)
    if f"{ident}_b_rcs" not in z:
        return None
    f_full = np.interp(rng, z[f"{ident}_range"], np.nan_to_num(z[f"{ident}_b_rcs"]),
                       left=0.0, right=0.0)
    band = fit_mask & np.isfinite(f_full)
    if band.sum() < 30:
        return None
    z2 = (rng / 1000.0) ** 2
    # seul z^2 est retire (il est porte par d_n dans le fit) ; la composante CONSTANTE du dark
    # est identifiable (rien dans le modele par nuit n'est plat en rcs_0) et reste dans f.
    # f_full (la forme capot complete) sert a RECONSTRUIRE le dark soustrayable theta * f_full.
    coef = float(np.sum(z2[band] * f_full[band]) / np.sum(z2[band] ** 2))
    f = f_full - coef * z2
    f[~fit_mask] = 0.0
    return f, f_full


def template_fit(S, sig, M, Aer, rng, fit_mask, f, n_iter: int = N_ALS_ITER):
    """Variante rigide de l'ALS : b(z) est remplace par theta * f(z), theta scalaire commun a
    toutes les nuits. Retourne A, B, d (par nuit), theta, sem_theta (jackknife nuits), d_med.

    theta est mis a jour par WLS global robuste des residus sur f ; la degenerescence residuelle
    (contenu de f colineaire a la forme moleculaire moyenne) se mesure au selftest : l'injection
    capot avec theta_vrai = 1 donne le facteur de recuperation lambda, qui CORRIGE theta."""
    z2 = (rng / 1000.0) ** 2
    n_nights = S.shape[0]
    A = np.full(n_nights, np.nan)
    B = np.zeros(n_nights)
    d = np.zeros(n_nights)
    theta = 0.0
    d_med = 0.0
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    u = np.zeros(n_nights)
    v = np.zeros(n_nights)

    for it in range(n_iter):
        bcur = theta * f + d_med * z2
        for n in range(n_nights):
            y0 = np.where(fit_mask, S[n] - bcur, np.nan)
            has_aer = np.nanmax(np.abs(Aer[n])) > 0
            coef = _wls_night(y0, [M[n], Aer[n], z2], w_all[n]) if has_aer else None
            if coef is not None and coef[1] < 0:
                coef = None
            if coef is None:
                c2 = _wls_night(y0, [M[n], z2], w_all[n])
                if c2 is None:
                    A[n] = np.nan
                    continue
                A[n], dn = c2
                B[n] = 0.0
            else:
                A[n], B[n], dn = coef
            d[n] = dn + d_med                     # d porte le TOTAL z^2 de la nuit
        ok_n = np.isfinite(A) & (A > 0)
        if ok_n.sum() < 5:
            return A, B, d, np.nan, np.nan, np.nan, ok_n
        d_med = float(np.median(d[ok_n]))
        # theta : WLS global robuste des residus (hors z^2 statique) sur f
        u[:] = 0.0
        v[:] = 0.0
        for n in np.where(ok_n)[0]:
            r = S[n] - A[n] * M[n] - B[n] * Aer[n] - d[n] * z2
            good = fit_mask & np.isfinite(r) & np.isfinite(w_all[n]) & (np.abs(f) > 0)
            if good.sum() < 30:
                continue
            resid_f = r - theta * f
            mad = 1.4826 * np.nanmedian(np.abs(resid_f[good]))
            if np.isfinite(mad) and mad > 0:
                good &= np.abs(resid_f) < 4 * mad
            u[n] = np.sum(w_all[n][good] * r[good] * f[good])
            v[n] = np.sum(w_all[n][good] * f[good] ** 2)
        if v.sum() <= 0:
            return A, B, d, np.nan, np.nan, d_med, ok_n
        theta_new = float(u.sum() / v.sum())
        conv = abs(theta_new - theta) < 1e-3 * (abs(theta_new) + 1e-9)
        theta = theta_new
        if it >= 2 and conv:
            break
    # SE jackknife sur les nuits (retire une nuit du ratio somme(u)/somme(v))
    okI = np.where(ok_n)[0]
    jk = np.array([(u.sum() - u[n]) / (v.sum() - v[n])
                   for n in okI if (v.sum() - v[n]) > 0])
    sem_theta = float(np.sqrt((jk.size - 1) * np.var(jk))) if jk.size > 3 else np.nan
    return A, B, d, theta, sem_theta, d_med, ok_n


def pergate_estimate(S, sig, M, Aer, rng, fit_mask, n_outer: int = 4, cref: float = np.nan):
    """L'estimateur DEFINITIF (v4) : pente ET ordonnee libres PAR PORTE.

        S_n(z) = A_n * m_n(z) * g(z) + B_n * a_n(z) + d_n * z^2 + b(z) + bruit

    g(z) statique capte toute erreur multiplicative de portee (residu d'overlap, normalisation
    proche portee) qui, dans le modele a pente commune (als_fit), se deversait dans b via le gros
    signal moleculaire de 2-4 km. b(z) est l'ordonnee a l'origine de la droite S vs A*m PAR PORTE
    a travers les nuits : avec CV(A) ~ 0.9 la droite est bien conditionnee et b est identifiable
    MEME dans la direction moleculaire (la degenerescence du modele a pente commune disparait --
    physiquement : les nuits a faible transmission revelent le plancher additif).

    Alternance : (a) par porte, droite ponderee robuste (S - Ba - dz^2) vs (A*m) -> g(z), b(z) ;
    (b) par nuit, WLS robuste [m*g, a, z^2] sur (S - b) -> A_n, B_n, d_n ;
    (c) normalisation <g> = 1 sur 2.5-6 km (l'echelle est portee par A_n).
    Retourne A, B, d, g, b_raw, ok_n.
    """
    z2 = (rng / 1000.0) ** 2
    n_nights, n_gates = S.shape
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    w_all[~np.isfinite(w_all)] = 0.0

    # init : fit a pente commune (il fournit des A_n de depart raisonnables)
    A, B, d, _b0, ok_n = als_fit(S, sig, M, Aer, rng, fit_mask, n_iter=4)
    if ok_n.sum() < 8:
        return A, B, d, np.ones(n_gates), np.full(n_gates, np.nan), ok_n
    g = np.ones(n_gates)
    b = np.zeros(n_gates)
    norm_band = fit_mask & (rng >= 2500) & (rng <= 6000)

    for _outer in range(n_outer):
        okI = np.where(ok_n)[0]
        # --- (a) droites par porte, vectorisees sur les portes ; 2 passes d'ecretage
        X = A[okI, None] * M[okI]                       # (n_ok, g)
        Y = S[okI] - B[okI, None] * Aer[okI] - d[okI, None] * z2[None, :]
        W = w_all[okI].copy()
        W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
        Xs = np.nan_to_num(X)
        Ys = np.nan_to_num(Y)
        for _pass in range(3):
            sw = W.sum(axis=0)
            sw[sw <= 0] = np.nan
            xbar = (W * Xs).sum(axis=0) / sw
            ybar = (W * Ys).sum(axis=0) / sw
            covxy = (W * (Xs - xbar) * (Ys - ybar)).sum(axis=0) / sw
            varx = (W * (Xs - xbar) ** 2).sum(axis=0) / sw
            with np.errstate(divide="ignore", invalid="ignore"):
                slope = covxy / varx
            slope = np.clip(slope, 0.5, 2.0)            # garde-fou portes sans levier
            # portes sans variation de x (tres haut, m minuscule) : droite degeneree -> g = 1
            # (on attribue le moleculaire predit tel quel, le reste va dans l'ordonnee)
            no_lev = ~np.isfinite(varx) | (varx < 1e-6 * np.maximum(xbar, 1e-300) ** 2)
            slope[no_lev] = 1.0
            # au-dessus de l'overlap la reponse est physiquement quasi constante : toute
            # structure < G_SMOOTH_M dans la pente est de l'atmosphere -> lissage fort
            dz = float(np.median(np.diff(rng))) if rng.size > 1 else 1.0
            if np.isfinite(dz) and dz > 0:
                k = max(3, int(round(G_SMOOTH_M / dz)) | 1)
                half = k // 2
                pad = np.pad(slope, half, mode="edge")
                slope = np.convolve(pad, np.ones(k) / k, mode="valid")
            inter = ybar - slope * xbar
            resid = Ys - slope[None, :] * Xs - inter[None, :]
            mad = 1.4826 * np.nanmedian(np.where(W > 0, np.abs(resid), np.nan), axis=0)
            bad = np.abs(resid) > 4 * np.where(np.isfinite(mad) & (mad > 0), mad, np.inf)[None, :]
            W = np.where(bad, 0.0, w_all[okI])
            W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
        g = np.where(fit_mask, slope, 1.0)
        b = np.where(fit_mask, inter, np.nan)
        # --- (c) normalisation d'echelle
        nb = norm_band & np.isfinite(g)
        c = float(np.nanmedian(g[nb])) if nb.sum() > 10 else 1.0
        if np.isfinite(c) and c > 0:
            g = np.where(fit_mask, g / c, 1.0)
            A[ok_n] *= c
        # --- (b) re-fit par nuit avec la reponse g
        Mg = M * g[None, :]
        for n in range(n_nights):
            y0 = np.where(fit_mask, S[n] - b, np.nan)
            has_aer = np.nanmax(np.abs(Aer[n])) > 0
            coef = _wls_night(y0, [Mg[n], Aer[n], z2], w_all[n]) if has_aer else None
            if coef is not None and coef[1] < 0:
                coef = None
            if coef is None:
                c2 = _wls_night(y0, [Mg[n], z2], w_all[n])
                if c2 is None:
                    A[n] = np.nan
                    continue
                A[n], d[n] = c2
                B[n] = 0.0
            else:
                A[n], B[n], d[n] = coef
        ok_n = np.isfinite(A) & (A > 0)
        # PAS d'ecran d'amplitude sur A_n : teste sous deux formes (borne au C_ref du run, puis
        # ecremage auto-reference [0.4, 2.2] x mediane), l'ecran DEGRADE le detecteur dans les
        # deux cas -- les nuits a fort A ancrent les droites par porte, et les retirer amplifie
        # l'attenuation erreurs-dans-les-variables (a Payerne A : lambda 0.74 -> 0.24, faux
        # positif theta0 0.00 -> -0.33, significativite 6 sigma -> 2 sigma). L'ecretage robuste
        # 4xMAD par porte suffit a contenir les nuits aberrantes.
        if ok_n.sum() < 8:
            break
    return A, B, d, g, b, ok_n


def decompose_b(rng, b, m_mean, a_mean, fit_mask):
    """Scinde b en sa part COLINEAIRE au sous-espace {m_moyen, a_moyen} (aliasee avec les
    amplitudes moleculaire/aerosol moyennes -> a ne pas croire) et sa part ORTHOGONALE (la
    signature identifiable ; le z^2 statique et la constante SONT identifiables : rien dans le
    modele par nuit ne les imite). Projection LS sur la bande de fit."""
    band = fit_mask & np.isfinite(b)
    cols = [m_mean, a_mean]
    G = np.column_stack([c[band] for c in cols])
    keep = np.all(np.isfinite(G), axis=1)
    coll = np.full(rng.size, np.nan)
    orth = np.full(rng.size, np.nan)
    if keep.sum() > 20:
        coef, *_ = np.linalg.lstsq(G[keep], b[band][keep], rcond=None)
        proj = sum(c * co for c, co in zip(cols, coef))
        coll[band] = proj[band]
        orth[band] = b[band] - proj[band]
    return coll, orth


def theta_projection(rng, b, f, fit_mask):
    """Projection de b sur la forme capot f (bande 2-9 km) : le resume scalaire 'theta' --
    theta = 1 signifie 'meme dark que l'unite Payerne du meme type'."""
    if f is None or b is None:
        return np.nan
    band = (fit_mask & (rng >= THETA_Z[0]) & (rng <= THETA_Z[1])
            & np.isfinite(b) & np.isfinite(f) & (np.abs(f) > 0))
    if band.sum() < 30 or np.sum(f[band] ** 2) <= 0:
        return np.nan
    return float(np.sum(b[band] * f[band]) / np.sum(f[band] ** 2))


def run_estimator(nights: list[dict], rng: np.ndarray, cref: float = np.nan):
    """Empile les nuits, applique l'estimateur par porte (v4), et calcule les diagnostics."""
    S = np.vstack([n["S"] for n in nights])
    sig = np.vstack([n["sigma"] for n in nights])
    M = np.vstack([n["m"] for n in nights])
    Aer = np.vstack([n["a"] for n in nights])
    fit_mask = (rng >= Z0_M) & (rng <= Z1_M)

    A, B, d, g, b_raw, ok_n = pergate_estimate(S, sig, M, Aer, rng, fit_mask, cref=cref)
    b = smooth_range(rng, b_raw)
    okI = np.where(ok_n)[0]

    # Decomposition (desormais un simple diagnostic : la v4 identifie aussi la direction ~m)
    m_mean = np.nanmean(M[okI], axis=0) if okI.size else np.full(rng.size, np.nan)
    a_mean = np.nanmean(Aer[okI], axis=0) if okI.size else np.zeros(rng.size)
    b_coll, b_orth = decompose_b(rng, b, m_mean, a_mean, fit_mask)

    # SEM bootstrap des ordonnees par porte (conditionnel aux A_n, B_n, d_n, g finals)
    z2 = (rng / 1000.0) ** 2
    X = A[okI, None] * M[okI] * g[None, :]
    Y = S[okI] - B[okI, None] * Aer[okI] - d[okI, None] * z2[None, :]
    W = 1.0 / np.maximum(sig[okI], 1e-300) ** 2
    W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
    Xs, Ys = np.nan_to_num(X), np.nan_to_num(Y)
    rs = np.random.default_rng(0)
    boots = np.full((200, rng.size), np.nan)
    for k in range(200):
        pick = rs.choice(okI.size, okI.size, replace=True)
        w, x, y = W[pick], Xs[pick], Ys[pick]
        sw = w.sum(axis=0)
        sw[sw <= 0] = np.nan
        xbar = (w * x).sum(axis=0) / sw
        ybar = (w * y).sum(axis=0) / sw
        covxy = (w * (x - xbar) * (y - ybar)).sum(axis=0) / sw
        varx = (w * (x - xbar) ** 2).sum(axis=0) / sw
        with np.errstate(divide="ignore", invalid="ignore"):
            sl = np.clip(covxy / varx, 0.2, 3.0)
        boots[k] = ybar - sl * xbar
    sem = np.where(fit_mask, np.nanstd(boots, axis=0), np.nan)

    # Split-half : re-fit v4 complet sur nuits paires / impaires (reproductibilite, rapport 08)
    halves = []
    for par in (0, 1):
        sel = okI[par::2]
        if sel.size >= 8:
            _, _, _, _, bh, _ = pergate_estimate(S[sel], sig[sel], M[sel], Aer[sel],
                                                 rng, fit_mask, n_outer=3, cref=cref)
            halves.append(smooth_range(rng, bh))
    split_r = split_amp = np.nan
    if len(halves) == 2:
        band = fit_mask & (rng <= 9000) & np.isfinite(halves[0]) & np.isfinite(halves[1])
        if band.sum() > 30:
            h0, h1 = halves[0][band], halves[1][band]
            if np.nanstd(h0) > 0 and np.nanstd(h1) > 0:
                split_r = float(np.corrcoef(h0, h1)[0, 1])
                split_amp = float(np.nanstd(h1) / np.nanstd(h0))

    # Resume scalaire : projection sur la forme capot du type (+ SE via les memes bootstraps).
    # La reconstruction soustrayable utilise la forme capot COMPLETE : theta * f_full.
    itype = nights[0]["itype"]
    tpl = load_template(itype, rng, fit_mask)
    f_tpl, f_full = tpl if tpl is not None else (None, None)
    theta = theta_projection(rng, b, f_tpl, fit_mask)
    th_boot = [theta_projection(rng, boots[k], f_tpl, fit_mask) for k in range(boots.shape[0])]
    th_boot = [t for t in th_boot if np.isfinite(t)]
    sem_theta = float(np.std(th_boot)) if len(th_boot) > 20 else np.nan
    b_tpl = theta * f_full if (f_full is not None and np.isfinite(theta)) \
        else np.full(rng.size, np.nan)

    return dict(A=A, B=B, d=d, g=g, b=b, b_raw=b_raw, b_coll=b_coll, b_orth=b_orth,
                sem=sem, ok_n=ok_n, m_mean=m_mean, a_mean=a_mean, cref=cref,
                theta=theta, sem_theta=sem_theta, f_tpl=f_tpl, f_full=f_full, b_tpl=b_tpl,
                split_r=split_r, split_amp=split_amp, fit_mask=fit_mask)


def predicted_bias_pct(rng, b, M, A, ok_n, windows, dates):
    """Biais Rayleigh predit : b integre sur les fenetres reellement ajustees, en % du signal
    moleculaire de la nuit. Negatif = la constante actuelle est biaisee BASSE (la correction
    l'augmentera), comme le -17 % / +24 % de Payerne A."""
    date_to_i = {d: i for i, d in enumerate(dates)}
    vals = []
    for date, bot, top in windows:
        i = date_to_i.get(date)
        if i is None or not ok_n[i]:
            continue
        w = (rng >= bot) & (rng <= top) & np.isfinite(b)
        if w.sum() < 5:
            continue
        mol = A[i] * M[i][w]
        if np.nansum(mol) <= 0:
            continue
        vals.append(100.0 * np.nansum(b[w]) / np.nansum(mol))
    if not vals:
        return np.nan, 0
    return float(np.median(vals)), len(vals)


# ---------------------------------------------------------------------------
# Self-test d'injection : la reponse chiffree a "est-ce identifiable ?"
# ---------------------------------------------------------------------------

def selftest(nights, rng, est, dark_npz: Path | None, ident: str):
    """Reconstruit des nuits synthetiques S' = A_n m_n + B_n a_n + b_vrai + d_n z^2 + bruit avec
    les A_n, B_n, d_n et sigma REELS du flux, pour trois b_vrai : (0) zero -- plancher de faux
    positif ; (1) la mesure capot Payerne (si dispo pour cet ident) ; (2) le pire cas degenere,
    un profil proportionnel a la forme moleculaire moyenne. La recuperation est rapportee POUR LE
    TOTAL et POUR CHAQUE COMPOSANTE (orthogonale = identifiable, colineaire = degeneree)."""
    S = np.vstack([n["S"] for n in nights])
    sig = np.vstack([n["sigma"] for n in nights])
    M = np.vstack([n["m"] for n in nights])
    Aer = np.vstack([n["a"] for n in nights])
    A, B, d, g = est["A"], est["B"], est["d"], est["g"]
    ok_n = est["ok_n"]
    fit_mask = est["fit_mask"]
    okI = np.where(ok_n)[0]
    z2 = (rng / 1000.0) ** 2
    rs = np.random.default_rng(1)

    m_mean, a_mean = est["m_mean"], est["a_mean"]
    cases = {"zero": np.zeros(rng.size)}
    # Injection "capot" = la forme capot du TYPE (theta_vrai = 1) pour TOUS les flux -- c'est
    # elle qui calibre la reponse affine (theta0, pente) du detecteur sur chaque flux ; a
    # Payerne elle coincide avec le capot de l'unite elle-meme.
    if est.get("f_full") is not None:
        cases["capot"] = np.nan_to_num(est["f_full"])
    elif dark_npz is not None and dark_npz.exists():
        z = np.load(dark_npz)
        if f"{ident}_b_rcs" in z:
            cases["capot"] = np.interp(rng, z[f"{ident}_range"],
                                       np.nan_to_num(z[f"{ident}_b_rcs"]), left=0, right=0)
    amp = np.nanmedian(np.abs(S[okI][:, fit_mask])) * 0.15
    cases["molecul"] = amp * m_mean / (np.nanmax(m_mean[fit_mask]) or 1.0)

    f_tpl = est.get("f_tpl")
    scale = float(np.nanmedian(np.abs(S[okI][:, fit_mask]))) or 1.0

    out = {}
    for name, b_true in cases.items():
        Ssyn = (A[okI, None] * M[okI] * g[None, :] + B[okI, None] * Aer[okI]
                + b_true[None, :] + d[okI, None] * z2[None, :]
                + rs.standard_normal((okI.size, rng.size)) * sig[okI])
        _, _, _, _, b_hat, _ = pergate_estimate(Ssyn, sig[okI], M[okI], Aer[okI],
                                                rng, fit_mask, n_outer=3,
                                                cref=est.get("cref", np.nan))
        b_hat_s = smooth_range(rng, b_hat)
        hc, ho = decompose_b(rng, b_hat_s, m_mean, a_mean, fit_mask)
        tc, to = decompose_b(rng, np.where(fit_mask, b_true, np.nan), m_mean, a_mean, fit_mask)
        band = fit_mask & (rng <= 9000) & np.isfinite(b_hat_s)

        def _rec(bh, bt):
            err = bh[band] - bt[band]
            rms_err = float(np.sqrt(np.nanmean(err ** 2)))
            denom = float(np.sqrt(np.nanmean(bt[band] ** 2)))
            # un denominateur < 1e-4 du signal ne definit pas une "recuperation"
            rec = float(1.0 - rms_err / denom) if denom > 1e-4 * scale else np.nan
            return dict(rms_err=rms_err, rms_true=denom, recovery=rec)

        theta_hat = theta_projection(rng, b_hat_s, f_tpl, fit_mask)

        out[name] = dict(
            total=_rec(b_hat_s, np.where(np.isfinite(b_hat_s), b_true, np.nan)),
            orth=_rec(ho, to),
            coll=_rec(hc, tc),
            theta_hat=float(theta_hat),
            corr=float(np.corrcoef(b_hat_s[band], b_true[band])[0, 1])
            if np.nanstd(b_true[band]) > 0 and np.nanstd(b_hat_s[band]) > 0 else np.nan,
        )
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _night_cache_path(wmo, ident, args) -> Path:
    return OUT_DIR / "cache" / f"{wmo}_{ident}_{args.start}_{args.end}_{args.max_nights}.npz"


def run_stream(wmo: str, ident: str, args) -> dict | None:
    dates = clear_nights_from_csv(Path(args.csv_root), wmo, ident, args.start, args.end)
    if len(dates) > args.max_nights:
        idx = np.linspace(0, len(dates) - 1, args.max_nights).astype(int)
        dates = [dates[i] for i in idx]
    print(f"[{wmo}_{ident}] {len(dates)} nuits claires candidates ({args.start}..{args.end})")
    if len(dates) < 10:
        print("  -> trop peu de nuits, abandon")
        return None

    # Cache des extractions par nuit : les lectures L1+CAMS dominent le cout ; en cache, une
    # variante d'estimateur se re-evalue en secondes.
    cache = _night_cache_path(wmo, ident, args)
    nights = None
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        rng = z["rng"]
        nights = [dict(date=str(int(dd)), rng=rng, S=z["S"][i], sigma=z["sigma"][i],
                       m=z["M"][i], a=z["Aer"][i], n_prof=int(z["n_prof"][i]),
                       itype=str(z["itype"]), wl_nm=float(z["wl_nm"]))
                  for i, dd in enumerate(z["dates"])]
        print(f"  {len(nights)} nuits depuis le cache ({cache.name})")
    if nights is None:
        jobs = [dict(date=d, wmo=wmo, ident=ident, l1_root=args.l1_root,
                     cams_folder=args.cams) for d in dates]
        nights = []
        if args.workers > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                for res in ex.map(process_night, jobs, chunksize=4):
                    if res is not None:
                        nights.append(res)
        else:
            for j in jobs:
                res = process_night(j)
                if res is not None:
                    nights.append(res)
        if len(nights) < 10:
            print(f"  -> {len(nights)} nuits exploitables, abandon")
            return None
        rng = nights[0]["rng"]
        nights = [n for n in nights if n["rng"].size == rng.size]
        print(f"  {len(nights)} nuits exploitables")
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            cache, rng=rng,
            S=np.vstack([n["S"] for n in nights]),
            sigma=np.vstack([n["sigma"] for n in nights]),
            M=np.vstack([n["m"] for n in nights]),
            Aer=np.vstack([n["a"] for n in nights]),
            dates=np.array([int(n["date"]) for n in nights]),
            n_prof=np.array([n["n_prof"] for n in nights]),
            itype=np.array(nights[0]["itype"]), wl_nm=np.array(nights[0]["wl_nm"]))
    if len(nights) < 10:
        print(f"  -> {len(nights)} nuits exploitables, abandon")
        return None

    if getattr(args, "no_aerosol", False):
        # test de sensibilite : sans le regresseur aerosol CAMS (B_n force a 0 partout)
        for n in nights:
            n["a"] = np.zeros_like(n["a"])
    cref = read_cref(Path(args.csv_root), wmo, ident)
    est = run_estimator(nights, rng, cref=cref)
    n_used = int(est["ok_n"].sum())
    A_ok = est["A"][est["ok_n"]]
    cv_A = float(np.nanstd(A_ok) / np.nanmean(A_ok)) if n_used else np.nan

    windows = fitted_windows_from_csv(Path(args.csv_root), wmo, ident)
    M = np.vstack([n["m"] for n in nights])
    bias_pct, n_bias = predicted_bias_pct(
        rng, est["b"], M, est["A"], est["ok_n"], windows, [n["date"] for n in nights])
    bias_orth, _ = predicted_bias_pct(
        rng, est["b_orth"], M, est["A"], est["ok_n"], windows, [n["date"] for n in nights])
    bias_tpl, _ = predicted_bias_pct(
        rng, est["b_tpl"], M, est["A"], est["ok_n"], windows, [n["date"] for n in nights])

    frac_aer = float(np.mean(est["B"][est["ok_n"]] > 0)) if n_used else np.nan
    band = est["fit_mask"] & (rng <= 9000)
    rms_orth = float(np.sqrt(np.nanmean(est["b_orth"][band] ** 2)))
    rms_coll = float(np.sqrt(np.nanmean(est["b_coll"][band] ** 2)))

    print(f"  ALS: {n_used} nuits retenues, CV(A) = {cv_A:.3f} (le levier de separation), "
          f"aerosol CAMS charge sur {100 * frac_aer:.0f} % des nuits")
    print(f"  split-half r = {est['split_r']:.3f}, ratio amplitude = {est['split_amp']:.3f}")
    print(f"  b: RMS orthogonal (identifiable) = {rms_orth:.3e}, "
          f"RMS colineaire (degenere) = {rms_coll:.3e}")
    print(f"  theta (forme capot {est.get('f_tpl') is not None and 'du type' or 'ABSENTE'}) = "
          f"{est['theta']:+.3f} +- {est['sem_theta']:.3f}")
    print(f"  biais Rayleigh predit = {bias_pct:+.2f} % (profil libre, contamine) / "
          f"{bias_orth:+.2f} % (orthogonal) / {bias_tpl:+.2f} % (template theta) "
          f"sur {n_bias} fenetres")

    st = None
    theta_corr = np.nan
    if args.selftest:
        dark_npz = DATA / "dark_profiles_payerne.npz" if wmo == "0-20000-0-06610" else None
        st = selftest(nights, rng, est, dark_npz, ident)
        for name, r in st.items():
            print(f"  selftest[{name}]: recovery total = {r['total']['recovery']:.2f} "
                  f"orth = {r['orth']['recovery']:.2f} coll = {r['coll']['recovery']:.2f} "
                  f"theta_hat = {r['theta_hat']:+.3f} corr = {r['corr']:.2f}")
        # Calibration AFFINE du detecteur par le selftest : l'injection zero mesure le biais
        # erreurs-dans-les-variables theta0 (A_n est estime du meme signal que y), l'injection
        # capot mesure la pente de reponse. theta_corr est l'estime debiaisee en unites
        # "dark Payerne du type".
        th0 = st.get("zero", {}).get("theta_hat", np.nan)
        th1 = st.get("capot", {}).get("theta_hat", np.nan)
        if np.isfinite(th0) and np.isfinite(th1) and abs(th1 - th0) > 0.05 \
                and np.isfinite(est["theta"]):
            theta_corr = float((est["theta"] - th0) / (th1 - th0))
            print(f"  calibration affine du selftest (theta0={th0:+.3f}, pente="
                  f"{th1 - th0:.3f}) : theta_corr = {theta_corr:+.3f}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sfx = "_noaer" if getattr(args, "no_aerosol", False) else ""
    out_npz = OUT_DIR / f"{wmo}_{ident}{sfx}.npz"
    np.savez(out_npz, **{
        f"{ident}_range": rng,
        f"{ident}_b_rcs": est["b"],          # directement consommable par ALC_DARK_PROFILE
        f"{ident}_b_raw": est["b_raw"],
        f"{ident}_b_coll": est["b_coll"],
        f"{ident}_b_orth": est["b_orth"],
        f"{ident}_b_tpl": est["b_tpl"],
        f"{ident}_sem": est["sem"],
        f"{ident}_g": est["g"],
        "A": est["A"], "B": est["B"], "d": est["d"],
        "theta": est["theta"], "sem_theta": est["sem_theta"],
        "dates": np.array([int(n["date"]) for n in nights]),
        "ok_n": est["ok_n"],
    })
    diag = dict(wmo=wmo, ident=ident, itype=nights[0]["itype"], wl_nm=nights[0]["wl_nm"],
                n_candidates=len(dates), n_usable=len(nights), n_used=n_used,
                cv_A=cv_A, frac_aer=frac_aer, split_r=est["split_r"], split_amp=est["split_amp"],
                rms_orth=rms_orth, rms_coll=rms_coll,
                theta=float(est["theta"]), sem_theta=float(est["sem_theta"]),
                theta_corr=theta_corr,
                bias_pct=bias_pct, bias_orth=bias_orth, bias_tpl=bias_tpl, n_bias=n_bias,
                start=args.start, end=args.end, selftest=st)
    with open(OUT_DIR / f"{wmo}_{ident}{sfx}_diag.json", "w", encoding="utf-8") as fh:
        json.dump(diag, fh, indent=1, ensure_ascii=False, default=float)
    print(f"  -> {out_npz}")

    if args.plot:
        _figure(wmo, ident, rng, est, nights)
    return diag


def _figure(wmo, ident, rng, est, nights):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.4))
    km = rng / 1000.0
    m = np.isfinite(est["b"])
    axes[0].plot(est["b_raw"][m], km[m], color="#bbb", lw=0.6, label="b brut")
    axes[0].plot(est["b"][m], km[m], color="#d62728", lw=1.8, label="b lisse")
    axes[0].plot(est["b_orth"][m], km[m], color="#1f77b4", lw=1.4, ls="--",
                 label="part orthogonale (identifiable)")
    if np.any(np.isfinite(est["b_tpl"])):
        mt = np.isfinite(est["b_tpl"])
        axes[0].plot(est["b_tpl"][mt], km[mt], color="#e8871a", lw=1.6,
                     label=f"template $\\theta$={est['theta']:+.2f}")
    axes[0].fill_betweenx(km[m], (est["b"] - est["sem"])[m], (est["b"] + est["sem"])[m],
                          color="#d62728", alpha=0.15, lw=0)
    dark_npz = DATA / "dark_profiles_payerne.npz"
    if wmo == "0-20000-0-06610" and dark_npz.exists():
        z = np.load(dark_npz)
        if f"{ident}_b_rcs" in z:
            axes[0].plot(z[f"{ident}_b_rcs"], z[f"{ident}_range"] / 1000.0,
                         color="#2ca02c", lw=1.8, label="mesure capot")
    axes[0].axvline(0, color="#444", ls=":", lw=1)
    axes[0].set_ylim(0, 15)
    axes[0].set_xlabel("b(z) [unites rcs_0]")
    axes[0].set_ylabel("portee [km]")
    axes[0].set_title("Ligne de base estimee (nuits claires)")
    axes[0].legend(fontsize=8.5)
    axes[0].grid(alpha=0.3)

    okI = np.where(est["ok_n"])[0]
    dts = [datetime.strptime(nights[i]["date"], "%Y%m%d") for i in okI]
    axes[1].plot(dts, est["A"][okI], ".", ms=4, color="#1f77b4")
    axes[1].set_title("Amplitude moleculaire A_n (levier du fit)")
    axes[1].set_ylabel("A_n [rcs_0 / beta_att]")
    axes[1].grid(alpha=0.3)
    for lab in axes[1].get_xticklabels():
        lab.set_rotation(30)

    mg = np.isfinite(est["g"]) & (rng >= Z0_M) & (rng <= Z1_M)
    axes[2].plot(est["g"][mg], km[mg], "-", color="#7048e8", lw=1.3)
    axes[2].axvline(1, color="#444", ls=":", lw=1)
    axes[2].set_xlim(0.5, 1.5)
    axes[2].set_ylim(0, 15)
    axes[2].set_xlabel("g(z) (reponse multiplicative statique)")
    axes[2].set_ylabel("portee [km]")
    axes[2].set_title("Erreur multiplicative capturee\n(residu overlap / normalisation)")
    axes[2].grid(alpha=0.3)

    fig.suptitle(f"Dark estime des nuits claires — {wmo} {ident} "
                 f"(split-half r={est['split_r']:.2f})", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = OUT_DIR / f"{wmo}_{ident}_dark_est.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    print(f"  figure -> {p}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wmo")
    ap.add_argument("--ident")
    ap.add_argument("--list", dest="stream_list",
                    help="fichier texte: une ligne par flux 'wmo ident [start]' (# = commentaire)")
    ap.add_argument("--start", default="20250101")
    ap.add_argument("--end", default="20260813")
    ap.add_argument("--l1-root", dest="l1_root", default=str(DEFAULT_L1))
    ap.add_argument("--csv-root", dest="csv_root", default=str(DEFAULT_CSV))
    ap.add_argument("--cams", default=DEFAULT_CAMS)
    ap.add_argument("--max-nights", dest="max_nights", type=int, default=220)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--no-aerosol", dest="no_aerosol", action="store_true",
                    help="test de sensibilite : sans regresseur aerosol CAMS (sorties _noaer)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    if args.stream_list:
        default_start = args.start
        for line in Path(args.stream_list).read_text(encoding="utf-8").splitlines():
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            args.start = parts[2] if len(parts) > 2 else default_start
            try:
                run_stream(parts[0], parts[1], args)
            except Exception as exc:  # noqa: BLE001 -- un flux casse n'arrete pas le scan
                print(f"[{parts[0]}_{parts[1]}] ECHEC : {exc}")
    else:
        if not (args.wmo and args.ident):
            ap.error("--wmo/--ident ou --list requis")
        run_stream(args.wmo, args.ident, args)


if __name__ == "__main__":
    main()
