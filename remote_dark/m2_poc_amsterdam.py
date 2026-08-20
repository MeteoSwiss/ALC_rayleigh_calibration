# -*- coding: utf-8 -*-
"""M2 proof-of-concept — the Schiphol quadruple: four CHM15k, one sky.

Why this cluster is the cleanest possible test of M2's core claim. The fleet-scale design must
separate the shared dark D(z) from the shared contamination K(z) (aerosol/misfit projecting onto
molecular). Here that whole problem VANISHES: units A–D at 0-20000-0-06240 sit at the same
coordinates, so on a shared night unit v IS the atmospheric reference for unit u — no CAMS, no
molecular model, no aerosol regressor. Per gate z, across shared clear nights n:

    S_u,n(z) = rho_uv(z) * S_v,n(z) + db_uv(z) + noise

  * rho_uv(z)  multiplicative: response/overlap/calibration ratio — must reproduce the KNOWN
               Schiphol-B static overlap deficit (~+20 % signal at 500–1000 m, back to 1 by
               ~1.6 km) as a free cross-validation;
  * db_uv(z)   additive: the DARK DIFFERENCE b_u − rho·b_v — the M2 quantity, contamination-free.

What co-location cannot give, stated up front: the common-mode absolute (all four share the
invisible component along the common sky) — that stays with the Payerne anchor. The POC tests
everything else M2 relies on: that unit dark DIFFERENCES are measurable from sky data, that they
live in the type family (gexp projection), and that the pairwise machinery is self-consistent
(triangle closure over the 6 pairs, 4 nodes, 1 gauge freedom).

Caveat carried, not hidden: plain WLS of S_u on S_v has errors-in-variables in the regressor;
with slope ~1 and night-to-night CV ~0.9 the attenuation is second-order, and the closure test
would expose it if it mattered.
"""
from __future__ import annotations

import numpy as np

from remote_dark.common import V4_DIR, ensure_out
from remote_dark import hood, models
from remote_dark.m1_fit import _nightly_cref

WMO = "0-20000-0-06240"
IDENTS = "ABCD"
CACHE = V4_DIR / "cache"


def _load(ident: str):
    z = np.load(CACHE / f"{WMO}_{ident}_20250101_20260813_220.npz", allow_pickle=True)
    return {k: z[k] for k in z.files}


def _pair_fit(Su, sig_u, Sv, sig_v, rng):
    """Per-gate robust line S_u = rho*S_v + db across shared nights (2 IRLS passes, 4-MAD)."""
    W = 1.0 / np.maximum(sig_u ** 2 + sig_v ** 2, 1e-300)      # both sides noisy
    W[~(np.isfinite(Su) & np.isfinite(Sv))] = 0.0
    X = np.nan_to_num(Sv)
    Y = np.nan_to_num(Su)
    rho = np.full(rng.size, np.nan)
    db = np.full(rng.size, np.nan)
    for _ in range(2):
        sw = W.sum(axis=0)
        sw[sw <= 0] = np.nan
        xb = (W * X).sum(axis=0) / sw
        yb = (W * Y).sum(axis=0) / sw
        cxy = (W * (X - xb) * (Y - yb)).sum(axis=0) / sw
        vx = (W * (X - xb) ** 2).sum(axis=0) / sw
        with np.errstate(invalid="ignore", divide="ignore"):
            rho = np.clip(cxy / vx, 0.3, 3.0)
        db = yb - rho * xb
        resid = Y - rho[None, :] * X - db[None, :]
        mad = 1.4826 * np.nanmedian(np.where(W > 0, np.abs(resid), np.nan), axis=0) + 1e-300
        W = np.where(np.abs(resid) < 4.0 * mad, W, 0.0)
    sem = np.sqrt(1.0 / np.maximum(np.where(np.isfinite(sw), sw, 0.0), 1e-300))
    return rho, db, sem


def main():
    data = {i: _load(i) for i in IDENTS}
    rng = data["A"]["rng"]

    # Per-night CALIBRATION normalisation: without it the multiplicative channel mixes overlap
    # with the units' calibration scales (measured: C runs ~20 % hot aloft), additive closure is
    # broken by rho != 1, and db across units is in incomparable counts. S/C_n puts every unit in
    # attenuated-backscatter-like units; per-NIGHT so drifts and undocumented swaps self-correct.
    for i in IDENTS:
        cmap = _nightly_cref(WMO, i, "CHM15k")
        ds8 = [str(int(d)) for d in data[i]["dates"]]
        C_n = np.array([cmap.get(d, np.nan) for d in ds8])
        med = np.nanmedian(C_n)
        C_n = np.where(np.isfinite(C_n), C_n, med)
        data[i]["S"] = data[i]["S"] / C_n[:, None]
        data[i]["sigma"] = data[i]["sigma"] / C_n[:, None]
        print(f"  {i}: era-median C = {med:.4g}  ({np.isfinite(C_n).sum()} nightly constants)")
    # shared nights: the sky is one, but each unit has its own good-night subset
    sets = {i: set(int(d) for d in data[i]["dates"]) for i in IDENTS}
    shared = sorted(set.intersection(*sets.values()))
    print(f"shared clear nights across the four units: {len(shared)}")
    idx = {i: {int(d): k for k, d in enumerate(data[i]["dates"])} for i in IDENTS}

    def nights_of(i):
        rows = [idx[i][d] for d in shared]
        return data[i]["S"][rows], data[i]["sigma"][rows]

    S = {i: nights_of(i)[0] for i in IDENTS}
    G = {i: nights_of(i)[1] for i in IDENTS}

    pairs = [(u, v) for ki, u in enumerate(IDENTS) for v in IDENTS[ki + 1:]]
    fits = {}
    for u, v in pairs:
        rho, db, sem = _pair_fit(S[u], G[u], S[v], G[v], rng)
        fits[u + v] = {"rho": rho, "db": db, "sem": sem}

    # ---- triangle closure: db_uv + db_vw - db_uw ~ 0 (rho ~ 1) --------------------------------
    band = (rng >= 1000) & (rng <= 10000)
    tri = []
    for a, b, c in (("A", "B", "C"), ("A", "B", "D"), ("A", "C", "D"), ("B", "C", "D")):
        d1 = fits[a + b]["db"]
        d2 = fits[b + c]["db"]
        d3 = fits[a + c]["db"]
        tri.append(np.nanmedian(np.abs((d1 + d2 - d3)[band])))
    typ = np.nanmedian([np.nanmedian(np.abs(f["db"][band])) for f in fits.values()])
    print(f"closure |db_uv+db_vw-db_uw| median (1-10 km): {np.nanmedian(tri):.3g} rcs  "
          f"vs typical |db| {typ:.3g}  -> ratio {np.nanmedian(tri)/max(typ,1e-30):.2f}")

    # ---- node solution: 6 pairwise db -> per-unit b offsets (gauge: mean = 0) -----------------
    P = np.zeros((len(pairs), 4))
    for r, (u, v) in enumerate(pairs):
        P[r, IDENTS.index(u)] = 1.0
        P[r, IDENTS.index(v)] = -1.0
    Pg = np.vstack([P, np.ones((1, 4)) * 0.5])                  # gauge row: sum b = 0
    node = np.full((4, rng.size), np.nan)
    for z in range(rng.size):
        y = np.array([fits[u + v]["db"][z] for u, v in pairs] + [0.0])
        if not np.all(np.isfinite(y[:-1])):
            continue
        sol, *_ = np.linalg.lstsq(Pg, y, rcond=None)
        node[:, z] = sol
    resid = np.array([fits[u + v]["db"] - (node[IDENTS.index(u)] - node[IDENTS.index(v)])
                      for u, v in pairs])
    print(f"node-solution residual median (1-10 km): "
          f"{np.nanmedian(np.abs(resid[:, band])):.3g} rcs")

    # ---- family projection: do the node differences live in the CHM15k dark family? ------------
    trng, t_brcs, tsem, _tp = hood.truth("A")
    fam0 = models.fit_family("CHM15k", trng, t_brcs, tsem)
    # the template in the SAME normalised units: the Payerne-A dark divided by Payerne-A's C.
    # delta_s is then "dark, in units of the Payerne dark" — comparable across any unit.
    cmap_pay = _nightly_cref("0-20000-0-06610", "A", "CHM15k")
    c_pay = float(np.nanmedian(list(cmap_pay.values()))) if cmap_pay else float("nan")
    D = fam0.rcs_view(rng) / c_pay
    print(f"\nPayerne A era-median C = {c_pay:.4g} (template scale)")
    m = band & np.isfinite(D)

    # TWO-basis projection, the same structure the fleet fit will use. Even with the atmosphere
    # common-moded away, the pairwise machinery has one leakage channel of its own: the noise of
    # rho-hat(z) times the mean signal -- smooth, molecular-shaped structure entering db. So each
    # node profile is projected on [D, L] with L = the quad-mean signal shape; delta_s is the
    # coefficient on D with the leakage direction absorbed by L, never silently mixed in.
    L = np.nanmean(np.stack([np.nanmean(S[i], axis=0) for i in IDENTS]), axis=0)
    L = L / max(float(np.nanmax(np.abs(L[m]))), 1e-30)
    print("\nper-unit dark offset vs the quad mean "
          "(two-basis projection [hood family, mean-signal leakage], 1-10 km):")
    s_node = {}
    for k, i in enumerate(IDENTS):
        bi = node[k]
        mm = m & np.isfinite(bi) & np.isfinite(L)
        Xb = np.column_stack([D[mm], L[mm]])
        try:
            beta, *_ = np.linalg.lstsq(Xb, bi[mm], rcond=None)
        except np.linalg.LinAlgError:
            continue
        pred = Xb @ beta
        r2 = 1.0 - float(np.nansum((bi[mm] - pred) ** 2)
                         / max(np.nansum((bi[mm] - np.nanmean(bi[mm])) ** 2), 1e-30))
        r2_d = 1.0 - float(np.nansum((bi[mm] - beta[0] * D[mm]) ** 2)
                           / max(np.nansum((bi[mm] - np.nanmean(bi[mm])) ** 2), 1e-30))
        s_node[i] = float(beta[0])
        print(f"  {i}: delta_s = {beta[0]:+.3f}  (leak coef {beta[1]:+.3g})   "
              f"R2 two-basis = {r2:+.2f}  (family alone {r2_d:+.2f})")

    # ---- the known cross-check: Schiphol B's static overlap deficit in rho --------------------
    print("\nrho(z) medians (multiplicative channel — the overlap/response ratio):")
    for u, v in (("A", "B"), ("A", "C"), ("A", "D")):
        r = fits[u + v]["rho"]
        b1 = float(np.nanmedian(r[(rng >= 500) & (rng <= 1000)]))
        b2 = float(np.nanmedian(r[(rng >= 1600) & (rng <= 3000)]))
        print(f"  {u}/{v}: 500-1000 m = {b1:.3f}   1.6-3 km = {b2:.3f}")
    print("  (prior knowledge: B's applied overlap ~25 % low at 500-1000 m -> S_B high there, "
          "back to ~1 by 1.6 km; A/B should sit ~0.8 in that band and ~1 above)")

    out = ensure_out("m2")
    np.savez(out / "poc_amsterdam.npz",
             rng=rng, node=node, idents=np.array(list(IDENTS)),
             **{f"db_{u}{v}": fits[u + v]["db"] for u, v in pairs},
             **{f"rho_{u}{v}": fits[u + v]["rho"] for u, v in pairs},
             n_shared=np.array(len(shared)),
             s_node=np.array([s_node[i] for i in IDENTS]))
    print(f"\n-> {out / 'poc_amsterdam.npz'}")


if __name__ == "__main__":
    main()
