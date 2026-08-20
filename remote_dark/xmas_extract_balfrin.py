# -*- coding: utf-8 -*-
"""Runs ON balfrin (login node, light) — extract the 2025-12-25/26 clean-airmass residuals.

Per stream whose P0 cache certifies those nights clear: the night's Rayleigh constant C_n
(stable 4.5-6.5 km window), and the RELATIVE residual per night

    delta_n(z) = (S_n - C_n * M_n) / (C_n * M_n)

i.e. the additive part expressed as a fraction of the molecular signal — dimensionless,
directly calibration-relevant, comparable across units regardless of firmware scale. The
CAMS-predicted aerosol fraction Aer/M is stored alongside (what the model itself expects of
the airmass). Output: one small npz to pull back.
"""
from __future__ import annotations

import numpy as np
from pathlib import Path

OUT = Path("/scratch/mch/mhrvo/remote_dark_net")
DATES = (20251225, 20251226)
CWIN = (4500.0, 6500.0)


def nightly_c(S, sig, M, rng):
    m = (rng >= CWIN[0]) & (rng <= CWIN[1])
    w = 1.0 / np.maximum(sig[:, m], 1e-300) ** 2
    ok = np.isfinite(S[:, m]) & np.isfinite(M[:, m])
    w = np.where(ok, w, 0.0)
    S0, M0 = np.nan_to_num(S[:, m]), np.nan_to_num(M[:, m])
    with np.errstate(all="ignore"):
        return (np.nansum(w * S0 * M0, axis=1)
                / np.maximum(np.nansum(w * M0 * M0, axis=1), 1e-300))


res = {}
n_streams = 0
for f in sorted(OUT.glob("0-*.npz")):
    try:
        z = np.load(f, allow_pickle=True)
        dates = z["dates"]
    except Exception:
        continue
    sel = np.isin(dates, DATES)
    if not sel.any():
        continue
    S, sig, M, Aer, rng = z["S"][sel], z["sigma"][sel], z["M"][sel], z["Aer"][sel], z["rng"]
    c = nightly_c(S, sig, M, rng)
    ok = np.isfinite(c) & (c > 0)
    if not ok.any():
        continue
    with np.errstate(all="ignore"):
        delta = (S[ok] - c[ok, None] * M[ok]) / (c[ok, None] * M[ok])
        aer_frac = np.nanmedian(np.where(M[ok] > 0, Aer[ok] / M[ok], np.nan), axis=0)
    key = f.stem
    res[f"rng_{key}"] = rng.astype("f4")
    res[f"delta_{key}"] = delta.astype("f4")
    res[f"aer_{key}"] = aer_frac.astype("f4")
    res[f"dates_{key}"] = dates[sel][ok]
    res[f"c_{key}"] = c[ok]
    n_streams += 1

np.savez_compressed(OUT / "xmas_residuals.npz", **res)
print(f"XMAS_DONE {n_streams} streams")
