# -*- coding: utf-8 -*-
"""M1 — physics-anchored parametric retrieval of the dark, from the v4 night caches.

What changes relative to v4. The v4 estimator fits a FREE per-gate baseline b(z) and is therefore
blind, by proven principle, to any component collinear with the mean molecular shape (worst case:
the CHM15k cuvette). M1 replaces the free b with the type's parametric FAMILY (models.py — a
handful of parameters whose shapes come from the electronics), plus one absolute anchor:

  * family        b(z) = fam(z; θ), the flat/z² parts excluded (they ride the per-night d_n·z²
                  regressor, exactly as in v4);
  * anchor        the nightly amplitudes A_n are no longer fully free: their median is pulled
                  toward C_ref·T̄² with a soft prior (T̄² = 0.93 ± 0.07, the climatological clear-
                  night aerosol two-way transmission). C_ref comes from the v2.2 archive; for
                  CL31/CL51/CL61 it is CLOUD-anchored and dark-immune, which is what makes the
                  prior legitimate. It converts the fatal "uniform A_n shift" freedom into a
                  bounded ±7 % leak — quantified per stream by the same injection self-test.

Alternation (reusing the v4 building blocks):
  (a) v4 pergate_estimate provides the initial state (A, B, d, g, ok_n) and the free b̂ for
      comparison;
  (b) family step: fit θ to the per-gate weighted intercepts of (S − B·a − d·z²) vs A·m·g,
      i.e. to the same free-b profile v4 would extract, but projected on the family;
  (c) night step: v4-style WLS per night on (S − fam), with the A-prior added as one synthetic
      observation per night (weight = prior);
  (d) iterate (b)-(c); the anchor is what keeps the family's molecular-collinear component from
      sliding.

Outputs one npz per stream under remote_dark/m1/: the family parameters, the retrieved b(z) in
rcs view, the free-b for reference, and the hood comparison when the stream is a Payerne one.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from rayleigh_availability import dark_from_clearsky as v4
from remote_dark import hood, models
from remote_dark.common import V4_DIR, ensure_out

#: The clear-night two-way aerosol transmission prior: median and softness of the anchor.
T2_PRIOR, T2_SD = 0.93, 0.07


def _load_cache(cache_file: Path):
    z = np.load(cache_file, allow_pickle=True)
    return {k: z[k] for k in z.files}


def _family_step(itype, rng, fit_mask, S, sig, M, Aer, A, B, d, g, ok_n, fam0):
    """Fit the family parameters to the per-gate intercept profile (the 'free b' of this state).

    The intercept and its SEM per gate come from the same weighted per-gate line v4 uses; the
    family is then fitted to THAT profile with SEM weights — so the family sees exactly the
    evidence a free b would see, but can only bend along its own shapes.
    """
    z2 = (rng / 1000.0) ** 2
    okI = np.where(ok_n)[0]
    X = A[okI, None] * M[okI] * g[None, :]
    Y = S[okI] - B[okI, None] * Aer[okI] - d[okI, None] * z2[None, :]
    W = 1.0 / np.maximum(sig[okI], 1e-300) ** 2
    W[~(np.isfinite(X) & np.isfinite(Y))] = 0.0
    Xs, Ys = np.nan_to_num(X), np.nan_to_num(Y)
    sw = W.sum(axis=0)
    sw[sw <= 0] = np.nan
    xbar = (W * Xs).sum(axis=0) / sw
    ybar = (W * Ys).sum(axis=0) / sw
    covxy = (W * (Xs - xbar) * (Ys - ybar)).sum(axis=0) / sw
    varx = (W * (Xs - xbar) ** 2).sum(axis=0) / sw
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = np.clip(covxy / varx, 0.5, 2.0)
    b_free = ybar - slope * xbar                       # per-gate intercept = the free-b evidence
    sem = np.sqrt(1.0 / np.maximum(sw, 1e-300))

    m = fit_mask & np.isfinite(b_free) & np.isfinite(sem) & (rng > 0)
    zf, yf, wf = rng[m], b_free[m], 1.0 / np.maximum(sem[m], 1e-300)
    wf /= np.nanmedian(wf)

    if fam0.kind == "damped2":
        # structured part only: c frozen at 0 (flat rides d_n); start from the hood-prior params
        p0 = fam0.params.copy()
        p0[0] = 0.0

        def resid(p):
            pp = p.copy()
            pp[0] = 0.0
            return (models.damped_sinusoids(zf, pp, 2) * zf ** 0 - yf / zf ** 2) * wf * zf ** 2
        r = least_squares(resid, p0, max_nfev=3000)
        pars = r.x
        pars[0] = 0.0
        fam = models.Family(fam0.itype, "damped2", pars, dict(fam0.meta))
    elif fam0.kind == "negexp":
        p0 = fam0.params.copy()
        p0[2] = 0.0

        def resid(p):
            return (models.negexp(zf, [p[0], p[1], 0.0]) - yf / zf ** 2) * wf * zf ** 2
        r = least_squares(resid, p0[:2].tolist() + [0.0],
                          bounds=([0, 100, -1e-30], [np.inf, 5e4, 1e-30]), max_nfev=2000)
        fam = models.Family(fam0.itype, "negexp", [r.x[0], r.x[1], 0.0], dict(fam0.meta))
    else:                                             # template: one amplitude, linear WLS
        t = fam0.p_view(zf) / max(float(fam0.params[0]), 1e-30) * zf ** 2   # template in rcs view
        num = float(np.nansum(wf * t * yf))
        den = float(np.nansum(wf * t * t)) + 1e-30
        fam = models.Family(fam0.itype, "template", np.array([num / den]), dict(fam0.meta))
    return fam, b_free, sem


def _night_step(rng, fit_mask, S, sig, M, Aer, g, b_fam, ok_n, a_prior):
    """v4-style per-night WLS on (S − b_fam), with the amplitude anchor as a synthetic point."""
    z2 = (rng / 1000.0) ** 2
    n = S.shape[0]
    A = np.full(n, np.nan)
    B = np.zeros(n)
    d = np.zeros(n)
    mu, sd = a_prior
    for i in range(n):
        if not ok_n[i]:
            continue
        y = S[i] - b_fam
        regs = [M[i] * g, Aer[i], z2]
        w = 1.0 / np.maximum(sig[i], 1e-300) ** 2
        m = fit_mask & np.isfinite(y) & np.all([np.isfinite(r) for r in regs], axis=0)
        if m.sum() < 30:
            ok_n[i] = False
            continue
        X = np.column_stack([r[m] for r in regs])
        # the anchor: one synthetic observation A_i ~ N(mu, sd) — a row [1, 0, 0] with weight
        # 1/sd² pulling the amplitude; it is the whole difference between "free" and "anchored"
        Xa = np.vstack([X * np.sqrt(w[m])[:, None],
                        np.array([[1.0, 0.0, 0.0]]) * (1.0 / sd)])
        ya = np.concatenate([y[m] * np.sqrt(w[m]), [mu / sd]])
        try:
            beta, *_ = np.linalg.lstsq(Xa, ya, rcond=None)
        except np.linalg.LinAlgError:
            ok_n[i] = False
            continue
        A[i], B[i], d[i] = beta
        B[i] = max(B[i], 0.0)
    return A, B, d, ok_n


def run_stream(cache_file: Path, ident_hint: str | None = None, n_outer: int = 3):
    c = _load_cache(cache_file)
    S, sig, M, Aer = c["S"], c["sigma"], c["M"], c["Aer"]
    rng = c["rng"]
    itype = str(c["itype"])
    dates = c["dates"]
    print(f"== {cache_file.name}: {S.shape[0]} nights, {itype}")

    fit_mask = np.isfinite(rng) & (rng >= 2000.0)
    if itype in ("CL31", "CL51"):
        fit_mask = np.isfinite(rng) & (rng >= 300.0)   # the CL31 dark lives near range

    # hood-prior family for the type (fit on the matching Payerne instrument)
    pay_ident = {"CHM15k": "A", "CL31": "B", "CL61": "C"}.get(itype)
    trng, tb, tsem, tbraw = hood.truth(pay_ident)
    fam0 = models.fit_family(itype, trng, tbraw, tsem)

    # C_ref from the v2.2 archive -> the anchor
    wmo, ident = cache_file.name.split("_")[0], cache_file.name.split("_")[1]
    cref = v4.read_cref(v4.DEFAULT_CSV, wmo, ident)
    if not np.isfinite(cref):
        print("   no C_ref -> anchor disabled (fit stays v4-free, family-only)")
        a_prior = (np.nan, np.inf)
    else:
        a_prior = (T2_PRIOR * cref, T2_SD * cref)

    # v4 initial state
    A, B, d, b0, ok_n = v4.als_fit(S, sig, M, Aer, rng, fit_mask, n_iter=4)
    g = np.ones(rng.size)
    fam = fam0
    b_free = sem = None
    for _ in range(n_outer):
        fam, b_free, sem = _family_step(itype, rng, fit_mask, S, sig, M, Aer,
                                        A, B, d, g, ok_n, fam)
        if np.isfinite(a_prior[0]):
            A, B, d, ok_n = _night_step(rng, fit_mask, S, sig, M, Aer, g,
                                        fam.rcs_view(rng), ok_n, a_prior)
    b_hat = fam.rcs_view(rng)

    out = {"rng": rng, "b_hat": b_hat, "b_free": b_free, "sem": sem,
           "params": fam.params, "kind": np.array(fam.kind), "itype": np.array(itype),
           "cref": np.array(cref), "n_nights": np.array(int(ok_n.sum()))}

    # hood comparison for the Payerne streams — the verdict
    if wmo == "0-20000-0-06610":
        hb = np.interp(rng, trng, tb)
        m = fit_mask & np.isfinite(b_hat) & np.isfinite(hb)
        band = (rng >= 2000) & (rng <= 8000) & m
        nearb = (rng >= 300) & (rng <= 1500) & m
        for name, mm in (("ray 2-8km", band), ("near 0.3-1.5km", nearb)):
            if mm.sum() < 10:
                continue
            cc = float(np.corrcoef(b_hat[mm], hb[mm])[0, 1])
            gain = float(np.nansum(b_hat[mm] * hb[mm]) / max(np.nansum(hb[mm] ** 2), 1e-30))
            print(f"   vs hood [{name}]: corr {cc:+.2f}  amplitude {gain:+.2f} (1 = perfect)")
        out["hood_b"] = hb
    d_out = ensure_out("m1")
    np.savez(d_out / f"{wmo}_{ident}_m1.npz", **out)
    print(f"   -> {d_out / f'{wmo}_{ident}_m1.npz'}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", help="explicit cache npz; default: every Payerne cache present")
    a = ap.parse_args()
    files = ([Path(a.cache)] if a.cache else
             sorted((V4_DIR / "cache").glob("0-20000-0-06610_*_2*.npz")))
    for f in files:
        try:
            run_stream(f)
        except Exception as exc:                                            # noqa: BLE001
            print(f"   FAILED {f.name}: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
