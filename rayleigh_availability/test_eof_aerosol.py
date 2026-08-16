# -*- coding: utf-8 -*-
"""Test : regresseur aerosol EMPIRIQUE par modes temporels (EOF) au lieu de CAMS.

Idee operateur : aerosols, bruit et moleculaire n'ont pas le meme temps de vie. La version
regressive de cette idee : les EOF de la matrice des residus inter-nuits (centree par porte,
donc aveugle a tout ce qui est STATIQUE, dark compris) capturent les formes verticales des
composantes qui FLUCTUENT — l'aerosol au premier chef. On reinjecte ces formes comme
regresseurs a chargement libre par nuit, a la place du profil CAMS dont la localisation
verticale est discutable. Benchmark : Payerne A (cache), verite capot theta = 1.

Run : python rayleigh_availability/test_eof_aerosol.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[1]))

from dark_from_clearsky import (  # noqa: E402
    _wls_night, load_template, smooth_range, theta_projection,
    Z0_M, Z1_M, G_SMOOTH_M, DATA,
)

CACHE = DATA / "dark_clearsky" / "cache" / "0-20000-0-06610_A_20250601_20260813_220.npz"
N_MODES = 2


def pergate_fit_eof(S, sig, M, V, rng, fit_mask, n_outer=4):
    """Variante v4 avec K regresseurs STATIQUES en forme (V: K x n_gates), chargements libres
    par nuit (signe libre : le signe d'une EOF est arbitraire)."""
    z2 = (rng / 1000.0) ** 2
    n_nights, n_gates = S.shape
    K = V.shape[0]
    w_all = 1.0 / np.maximum(sig, 1e-300) ** 2
    w_all[~np.isfinite(w_all)] = 0.0
    A = np.full(n_nights, np.nan)
    L = np.zeros((n_nights, K))
    d = np.zeros(n_nights)
    g = np.ones(n_gates)
    b = np.zeros(n_gates)
    norm_band = fit_mask & (rng >= 2500) & (rng <= 6000)
    dz = float(np.median(np.diff(rng)))
    k_g = max(3, int(round(G_SMOOTH_M / dz)) | 1)

    for _outer in range(n_outer):
        Mg = M * g[None, :]
        for n in range(n_nights):
            y0 = np.where(fit_mask, S[n] - b, np.nan)
            coef = _wls_night(y0, [Mg[n]] + [V[k] for k in range(K)] + [z2], w_all[n])
            if coef is None:
                A[n] = np.nan
                continue
            A[n] = coef[0]
            L[n] = coef[1:1 + K]
            d[n] = coef[-1]
        ok_n = np.isfinite(A) & (A > 0)
        if ok_n.sum() < 8:
            return A, L, d, g, np.full(n_gates, np.nan), ok_n
        okI = np.where(ok_n)[0]
        X = A[okI, None] * M[okI] * g[None, :]
        Y = (S[okI] - np.einsum("nk,kz->nz", L[okI], V) - d[okI, None] * z2[None, :])
        W = w_all[okI].copy()
        W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
        Xs, Ys = np.nan_to_num(X), np.nan_to_num(Y)
        for _pass in range(3):
            sw = W.sum(axis=0)
            sw[sw <= 0] = np.nan
            xbar = (W * Xs).sum(axis=0) / sw
            ybar = (W * Ys).sum(axis=0) / sw
            covxy = (W * (Xs - xbar) * (Ys - ybar)).sum(axis=0) / sw
            varx = (W * (Xs - xbar) ** 2).sum(axis=0) / sw
            with np.errstate(divide="ignore", invalid="ignore"):
                slope = np.clip(covxy / varx, 0.5, 2.0)
            no_lev = ~np.isfinite(varx) | (varx < 1e-6 * np.maximum(xbar, 1e-300) ** 2)
            slope[no_lev] = 1.0
            half = k_g // 2
            pad = np.pad(slope, half, mode="edge")
            slope = np.convolve(pad, np.ones(k_g) / k_g, mode="valid")
            inter = ybar - slope * xbar
            resid = Ys - slope[None, :] * Xs - inter[None, :]
            mad = 1.4826 * np.nanmedian(np.where(W > 0, np.abs(resid), np.nan), axis=0)
            bad = np.abs(resid) > 4 * np.where(np.isfinite(mad) & (mad > 0), mad, np.inf)[None, :]
            W = np.where(bad, 0.0, w_all[okI])
            W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
        gnew = np.where(fit_mask, slope / (np.nanmedian(slope[norm_band]) or 1.0), 1.0)
        A[ok_n] *= float(np.nanmedian(slope[norm_band]) or 1.0)
        g = gnew
        b = np.where(fit_mask, inter, np.nan)
    return A, L, d, g, b, ok_n


def eof_shapes(S, sig, M, rng, fit_mask, n_modes=N_MODES):
    """Fit sans aerosol -> EOF des residus centres par porte (aveugles au statique)."""
    zero = np.zeros_like(S)
    A, L, d, g, b, ok_n = pergate_fit_eof(S, sig, M, np.zeros((1, rng.size)), rng, fit_mask)
    okI = np.where(ok_n)[0]
    z2 = (rng / 1000.0) ** 2
    R = (S[okI] - A[okI, None] * M[okI] * g[None, :] - d[okI, None] * z2[None, :]
         - np.nan_to_num(b)[None, :])
    band = fit_mask & (rng <= 12000)
    Rb = np.nan_to_num(R[:, band])
    Rb = Rb - Rb.mean(axis=0, keepdims=True)      # centrage par porte : tue le statique
    # ponderation par la dispersion de porte pour eviter que le champ lointain bruyant domine
    scale = np.median(np.abs(Rb), axis=0) + 1e-300
    U, s, Vt = np.linalg.svd(Rb / scale[None, :], full_matrices=False)
    V = np.zeros((n_modes, rng.size))
    for k in range(n_modes):
        V[k, band] = Vt[k] * scale               # re-dimensionne la forme
        V[k] = smooth_range(rng, np.where(band, V[k], np.nan))
        V[k] = np.nan_to_num(V[k])
    ev = (s ** 2 / np.sum(s ** 2))[:n_modes]
    return V, ev, (A, d, g, b, ok_n)


def run_case(S, sig, M, V, rng, fit_mask, f_tpl, label):
    A, L, d, g, b, ok_n = pergate_fit_eof(S, sig, M, V, rng, fit_mask)
    b_s = smooth_range(rng, b)
    th = theta_projection(rng, b_s, f_tpl, fit_mask)
    print(f"  [{label}] nuits={int(ok_n.sum())}  theta = {th:+.3f}")
    return th, (A, L, d, g, ok_n)


def main():
    z = np.load(CACHE)
    rng = z["rng"]
    S, sig, M = z["S"], z["sigma"], z["M"]
    fit_mask = (rng >= Z0_M) & (rng <= Z1_M)
    tpl = load_template("CHM15k", rng, fit_mask)
    f_tpl, f_full = tpl

    print("== extraction des modes EOF temporels ==")
    V, ev, _ = eof_shapes(S, sig, M, rng, fit_mask)
    print(f"  variance expliquee par les {N_MODES} modes : {100 * ev.sum():.1f} % "
          f"({', '.join(f'{100 * e:.1f}%' for e in ev)})")

    print("== fit avec regresseurs EOF (donnees reelles) ==")
    th_real, (A, L, d, g, ok_n) = run_case(S, sig, M, V, rng, fit_mask, f_tpl, "reel")

    # calibration affine par injection, avec le MEME pipeline EOF
    okI = np.where(ok_n)[0]
    z2 = (rng / 1000.0) ** 2
    rs = np.random.default_rng(1)
    ths = {}
    for name, b_true in (("zero", np.zeros(rng.size)),
                         ("capot", np.nan_to_num(f_full))):
        Ssyn = (A[okI, None] * M[okI] * g[None, :]
                + np.einsum("nk,kz->nz", L[okI], V)
                + d[okI, None] * z2[None, :] + b_true[None, :]
                + rs.standard_normal((okI.size, rng.size)) * sig[okI])
        # re-extraction EOF sur le synthetique (pipeline complet, honnete)
        Vs, _, _ = eof_shapes(Ssyn, sig[okI], M[okI], rng, fit_mask)
        th, _ = run_case(Ssyn, sig[okI], M[okI], Vs, rng, fit_mask, f_tpl, f"inj:{name}")
        ths[name] = th
    lam = ths["capot"] - ths["zero"]
    print(f"== calibration affine : theta0 = {ths['zero']:+.3f}, pente = {lam:.3f} ==")
    if abs(lam) > 0.05:
        print(f"  theta_corr (verite = 1) : {(th_real - ths['zero']) / lam:+.3f}")
    print(f"  Rappel CAMS : theta = +1.05 +- 0.39 (corr +1.49) ; sans aerosol : -0.12")


if __name__ == "__main__":
    main()
