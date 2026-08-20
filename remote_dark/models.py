# -*- coding: utf-8 -*-
"""Parametric dark families per instrument type, and their fit to the Payerne hood truth.

Why parametric. The v4 campaign proved that a FREE per-gate baseline is degenerate with the mean
molecular shape (the identifiability wall). The families here replace ~2000 free gates with a
handful of parameters whose SHAPES come from the electronics:

  CL31 / CL51  P(z) = c + Σ_k A_k · exp(−z/τ_k) · sin(2π z/λ_k + φ_k)    (K under-damped resonances
               + a range-flat offset; the hood physical-model work found K=2: AC-coupling ring +
               transmitter ripple, and the network-offset work showed the flat part is real)
  CHM15k       P(z) = −A · exp(−z/τ) + c                                 (firmware background
               over-subtraction; the "cuvette" in rcs view is this times z²)
  CL61         b(z) = s · template(z)                                    (the hood shape itself with
               a free amplitude — the v4 θ-projection formalised; the CL61 dark is small and
               type-common, Looschelders 2025)

All families are defined in P-view (P = b/z², raw-signal units) where the electronics are smooth,
and evaluated in rcs view (b = P·z²) where the fits and the truth live. Fitting the families to
the hood profiles does two jobs at once: (1) it measures how well each family can represent the
truth at all (the ceiling of what M1 can ever achieve), and (2) it produces the per-type parameter
PRIORS the anchored sky fit starts from.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from remote_dark.common import ensure_out
from remote_dark import hood

#: Fit band per type (m). Below the floor the hood itself is transient-contaminated; above the
#: ceiling the truth is pure noise. CL31 range ends at 7.7 km; CL61 structure lives <1.5 km.
FIT_BAND = {"CL31": (60.0, 7000.0), "CL51": (60.0, 7000.0),
            "CHM15k": (200.0, 12000.0), "CL61": (10.0, 3000.0)}


# ------------------------------------------------------------------------------------ families ---
def damped_sinusoids(z, p, K):
    """c + Σ_k A_k exp(−z/τ_k) sin(2π z/λ_k + φ_k) — p = [c] + [A, τ, λ, φ] × K, P-view.

    The flat c matters twice over: in b view it is the dominant far-range shape (c·z²), and in the
    SKY fit it is exactly the shape the per-night d_n·z² regressor absorbs — so the sky fit uses
    the family WITHOUT c and the truth comparison separates the structured part accordingly.
    """
    out = np.full_like(z, p[0], dtype=float)
    for k in range(K):
        A, tau, lam, phi = p[1 + 4 * k:1 + 4 * k + 4]
        out += A * np.exp(-z / tau) * np.sin(2 * np.pi * z / lam + phi)
    return out


def negexp(z, p):
    """−A exp(−z/τ) + c — the CHM15k over-subtraction family, P-view."""
    A, tau, c = p
    return -A * np.exp(-z / tau) + c


class Family:
    """One fitted family: evaluate in P or rcs view, serialise, and report its truth ceiling."""

    def __init__(self, itype: str, kind: str, params: np.ndarray, meta: dict):
        self.itype, self.kind, self.params, self.meta = itype, kind, np.asarray(params), meta

    def p_view(self, z):
        if self.kind == "damped2":
            return damped_sinusoids(z, self.params, 2)
        if self.kind == "negexp":
            return negexp(z, self.params)
        if self.kind == "template":
            t = np.interp(z, self.meta["tpl_rng"], self.meta["tpl_p"], left=0.0, right=0.0)
            return self.params[0] * t
        raise ValueError(self.kind)

    def rcs_view(self, z):
        return self.p_view(z) * z ** 2

    def to_dict(self):
        return {"itype": self.itype, "kind": self.kind, "params": self.params.tolist(),
                "r2_b": self.meta.get("r2_b"), "r2_p": self.meta.get("r2_p"),
                "r2_p_near": self.meta.get("r2_p_near"),
                "r2_b_struct": self.meta.get("r2_b_struct"),
                "ray_pct": self.meta.get("ray_pct"), "band_m": self.meta.get("band_m")}


# ------------------------------------------------------------------------------- hood fitting ---
def _r2(y, yhat, w):
    m = np.isfinite(y) & np.isfinite(yhat) & np.isfinite(w)
    if m.sum() < 10:
        return float("nan")
    resid = (y[m] - yhat[m]) * w[m]
    base = (y[m] - np.average(y[m], weights=w[m])) * w[m]
    return 1.0 - float(np.sum(resid ** 2) / np.sum(base ** 2))


def fit_family(itype: str, rng: np.ndarray, b_raw: np.ndarray, sem: np.ndarray) -> Family:
    """Fit the type's family to a hood profile (rcs view in, P-view fit, both R² out).

    The fit weights are 1/sem·z² — i.e. uniform in P-view where the electronics live; without the
    z² the far gates (huge in rcs view) would drown the near-range structure that IS the dark.
    """
    lo, hi = FIT_BAND.get(itype, (60.0, 8000.0))
    m = (rng >= lo) & (rng <= hi) & np.isfinite(b_raw) & (rng > 0)
    z, y = rng[m], b_raw[m] / rng[m] ** 2                       # -> P-view
    w = np.ones_like(z)
    s = np.where(np.isfinite(sem[m]) & (sem[m] > 0), sem[m] / rng[m] ** 2, np.nan)
    if np.isfinite(s).sum() > 10:
        w = 1.0 / np.where(np.isfinite(s), s, np.nanmedian(s))
        w /= np.nanmedian(w)

    if itype in ("CL31", "CL51"):
        # Multi-start on the two range periods: the loss is very multi-modal in λ. Periods from
        # the hood physics: the ring ~1 km-scale, the transmitter ripple slower; scan generously.
        best = None
        amp0 = np.nanstd(y)
        c0 = float(np.nanmedian(y[z > 0.6 * z.max()]))
        for lam1 in (300.0, 600.0, 1000.0, 1500.0):
            for lam2 in (2000.0, 3500.0, 5000.0):
                p0 = [c0, amp0, 800.0, lam1, 0.0, amp0 / 3, 3000.0, lam2, 0.0]
                try:
                    r = least_squares(
                        lambda p: (damped_sinusoids(z, p, 2) - y) * w, p0,
                        bounds=([-np.inf, 0, 50, 100, -np.pi, 0, 50, 500, -np.pi],
                                [np.inf, np.inf, 2e4, 4000, np.pi, np.inf, 5e4, 2e4, np.pi]),
                        max_nfev=4000)
                except Exception:                                           # noqa: BLE001
                    continue
                if best is None or r.cost < best.cost:
                    best = r
        fam = Family(itype, "damped2", best.x, {})
    elif itype == "CHM15k":
        p0 = [max(np.nanmax(-y), 1e-12), 2000.0, float(np.nanmedian(y[-50:]))]
        r = least_squares(lambda p: (negexp(z, p) - y) * w, p0,
                          bounds=([0, 100, -np.inf], [np.inf, 5e4, np.inf]), max_nfev=2000)
        fam = Family(itype, "negexp", r.x, {})
    else:                                                        # CL61 and anything template-like
        tpl = y / max(float(np.nanmax(np.abs(y))), 1e-30)
        fam = Family(itype, "template", np.array([float(np.nanmax(np.abs(y)))]),
                     {"tpl_rng": z.copy(), "tpl_p": tpl.copy()})

    fam.meta["band_m"] = [lo, hi]
    fam.meta["r2_p"] = _r2(y, fam.p_view(z), w)
    bm = b_raw[m]
    fam.meta["r2_b"] = _r2(bm, fam.rcs_view(z), np.ones_like(bm))

    # The metrics that decide anything, separated by what consumes them:
    #  * r2_p_near   — P-view, 60–1500 m: the band where the CL31 correction acts on raw signal;
    #  * r2_b_struct — b view with the flat c·z² removed from BOTH sides: the structured part,
    #                  which is what the sky fit must recover (its z² part rides d_n);
    #  * ray_pct     — the family-vs-truth disagreement integrated over 2–8 km in b view, as % of
    #                  the truth integral: what a Rayleigh window would feel if the family were
    #                  subtracted instead of the hood profile.
    flat = (fam.params[0] if fam.kind == "damped2" else
            fam.params[2] if fam.kind == "negexp" else 0.0)
    y_st, f_st = y - flat, fam.p_view(z) - flat
    near = z <= 1500.0
    fam.meta["r2_p_near"] = _r2(y[near], fam.p_view(z[near]), w[near])
    fam.meta["r2_b_struct"] = _r2(y_st * z ** 2, f_st * z ** 2, np.ones_like(z))
    ray = (z >= 2000.0) & (z <= 8000.0)
    ti = float(np.nansum(np.abs(y_st[ray] * z[ray] ** 2)))
    fam.meta["ray_pct"] = (100.0 * float(np.nansum(np.abs((y_st - f_st)[ray] * z[ray] ** 2)))
                           / ti) if ti > 0 else float("nan")
    return fam


def fit_payerne(save: bool = True) -> dict:
    """Fit all three Payerne instruments; save params + a JSON summary. Returns {ident: Family}."""
    from remote_dark.common import PAYERNE
    out = {}
    summary = {}
    for ident, itype in PAYERNE["types"].items():
        rng, b, sem, b_raw = hood.truth(ident)
        fam = fit_family(itype, rng, b_raw, sem)
        out[ident] = fam
        summary[ident] = fam.to_dict()
        print(f"  {ident} ({itype:6s}) kind={fam.kind:9s} "
              f"R2(P)={fam.meta['r2_p']:.3f}  R2(P,near)={fam.meta['r2_p_near']:.3f}  "
              f"R2(b,struct)={fam.meta['r2_b_struct']:.3f}  ray-band misfit="
              f"{fam.meta['ray_pct']:.0f}%")
    if save:
        d = ensure_out("models")
        (d / "payerne_families.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
        np.savez(d / "payerne_families.npz", **{
            f"{ident}_params": fam.params for ident, fam in out.items()
        } | {f"{ident}_kind": np.array(fam.kind) for ident, fam in out.items()})
    return out


if __name__ == "__main__":
    fit_payerne()
