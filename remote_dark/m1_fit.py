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


def _nightly_cref(wmo: str, ident: str, itype: str) -> dict:
    """date(YYYYMMDD) -> raw successful-night constant, from the v2.2 archive (per-night, never
    the Kalman: it smooths across instrument swaps)."""
    import csv as _csv
    key = f"{wmo}_{ident}"
    f = v4.DEFAULT_CSV / key / f"{key}_cal.csv"
    out = {}
    if not f.exists():
        return out
    method = "cloud" if itype in ("CL31", "CL51", "CL61") else "rayleigh"
    try:
        with f.open(newline="", encoding="utf-8") as fh:
            for r in _csv.DictReader(fh):
                if r.get("method") != method:
                    continue
                try:
                    if float(r.get("flag", -9)) in (1.0, 0.5) and float(r["cal_value"]) > 0:
                        out[str(r.get("date", ""))[:8]] = float(r["cal_value"])
                except (TypeError, ValueError):
                    pass
    except OSError:
        return {}
    return out


def _era_cref(wmo: str, ident: str, itype: str, dates) -> tuple[float, float]:
    """C_ref for the anchor, scoped to THE CACHE'S OWN DATE SPAN — never the whole archive.

    Session-2 lesson: a whole-archive median mixes instrument eras (the Payerne CHM15k unit swap
    changes C by ~2x), and an anchor built on the wrong era drags every amplitude toward the wrong
    unit — the family then absorbs the missing molecular signal and the fit degrades with each
    iteration.

    Method choice and honesty: CL31/CL51/CL61 use the CLOUD Kalman — dark-immune, which is what
    makes the anchor legitimate. The CHM15k has no cloud calibration (it saturates in liquid
    cloud); its only constant is the Rayleigh one, which is itself dark-AFFECTED — weakly circular.
    It is still used, with the prior width DOUBLED, and the circularity is stated rather than
    hidden: a few-percent bias on C moves the anchor by a few percent, far less than the 2x an
    era-mixed median caused.

    Returns (cref, sd_scale) — sd_scale multiplies T2_SD.
    """
    import csv as _csv
    key = f"{wmo}_{ident}"
    # RAW per-night constants, not the Kalman: the Kalman smooths ACROSS instrument swaps (the
    # documented CL31 failure mode — a whole-year filter is 2.5x wrong at a gain break), so even an
    # era-scoped Kalman median is contaminated near a swap. Raw successful-night medians are not.
    f = v4.DEFAULT_CSV / key / f"{key}_cal.csv"
    if not f.exists():
        return float("nan"), 1.0
    method = "cloud" if itype in ("CL31", "CL51", "CL61") else "rayleigh"
    ds = sorted(str(d) for d in np.asarray(dates).ravel())
    lo, hi = ds[0][:8], ds[-1][:8]
    vals = []
    try:
        with f.open(newline="", encoding="utf-8") as fh:
            for r in _csv.DictReader(fh):
                if r.get("method") != method:
                    continue
                d8 = str(r.get("date", ""))[:8]
                if not (lo <= d8 <= hi):
                    continue
                try:
                    if float(r.get("flag", -9)) in (1.0, 0.5):
                        v = float(r["cal_value"])
                        if v > 0:
                            vals.append(v)
                except (TypeError, ValueError):
                    pass
    except OSError:
        return float("nan"), 1.0
    if len(vals) < 5:
        return float("nan"), 1.0
    return float(np.nanmedian(vals)), (2.0 if method == "rayleigh" else 1.0)


def _load_cache(cache_file: Path):
    z = np.load(cache_file, allow_pickle=True)
    return {k: z[k] for k in z.files}


def _family_step(itype, rng, fit_mask, S, sig, M, Aer, A, B, d, g, ok_n, fam0,
                 fam_prior=None):
    """Fit the family parameters to the per-gate intercept profile (the 'free b' of this state).

    The intercept and its SEM per gate come from the same weighted per-gate line v4 uses; the
    family is then fitted to THAT profile with SEM weights — so the family sees exactly the
    evidence a free b would see, but can only bend along its own shapes.
    """
    z2 = (rng / 1000.0) ** 2
    okI = np.where(ok_n)[0]
    X = A[okI, None] * M[okI] * g[None, :]
    # v4's median(d)=0 convention, reinstated: the STATIC part of the firmware z^2 term belongs to
    # the baseline evidence, only the per-night VARIATION rides d_n. Without this the night step
    # absorbs any z^2-collinear family component (the CL61 template died of exactly that).
    d_ctr = d - np.nanmedian(d[okI])
    Y = S[okI] - B[okI, None] * Aer[okI] - d_ctr[okI, None] * z2[None, :]
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
        # Structured part only: c frozen at 0 (flat rides d_n). The SHAPE parameters (decay
        # lengths, range periods) are bounded to the hood prior: that is the whole point of having
        # hood physics -- the sky data adjusts amplitudes and phases, it does not reinvent the
        # electronics. Unbounded, the resonances also over-ring above 1.2 km (fig 1).
        pr = (fam_prior or fam0).params.copy()          # bounds NEVER walk with the iterate
        p0 = fam0.params.copy()
        p0[0] = 0.0
        lob = np.full(p0.size, -np.inf)
        hib = np.full(p0.size, np.inf)
        lob[0], hib[0] = -1e-30, 1e-30                       # c == 0
        for k in range(2):
            iA, itau, ilam, iphi = 1 + 4 * k, 2 + 4 * k, 3 + 4 * k, 4 + 4 * k
            lob[iA], hib[iA] = 0.0, np.inf                   # amplitude free (>= 0)
            lob[itau], hib[itau] = 0.6 * pr[itau], 1.6 * pr[itau]
            lob[ilam], hib[ilam] = 0.85 * pr[ilam], 1.18 * pr[ilam]
            lob[iphi], hib[iphi] = pr[iphi] - 0.7, pr[iphi] + 0.7
        def resid(p):
            return (models.damped_sinusoids(zf, p, 2) - yf / zf ** 2) * wf * zf ** 2
        r = least_squares(resid, np.clip(p0, lob + 1e-12, hib - 1e-12),
                          bounds=(lob, hib), max_nfev=3000)
        fam = models.Family(fam0.itype, "damped2", r.x, dict(fam0.meta))
    elif fam0.kind == "gexp":
        # Shape (k, tau) bounded to the FIXED hood prior; amplitude free; c = 0 (its static part
        # rides the evidence now that median(d)=0 is reinstated). Unbounded tau turned the family
        # into a pure z^2 shape and exploded the amplitude (session 1); walking bounds collapsed
        # tau to the floor every iteration (session-2 diagnostics).
        pr = (fam_prior or fam0).params
        k0, tau0 = float(pr[1]), float(pr[2])
        p0 = fam0.params.copy()
        def resid(p):
            return (models.gexp(zf, [p[0], p[1], p[2], 0.0]) - yf / zf ** 2) * wf * zf ** 2
        r = least_squares(resid, [max(float(p0[0]), 1e-12),
                                  float(np.clip(float(p0[1]), 0.7 * k0, 1.4 * k0)),
                                  float(np.clip(float(p0[2]), 0.7 * tau0, 1.4 * tau0))],
                          bounds=([0, 0.7 * k0, 0.7 * tau0], [np.inf, 1.4 * k0, 1.4 * tau0]),
                          max_nfev=3000)
        fam = models.Family(fam0.itype, "gexp", [r.x[0], r.x[1], r.x[2], 0.0], dict(fam0.meta))
    else:                                             # template: one amplitude, linear WLS
        t = fam0.p_view(zf) / max(float(fam0.params[0]), 1e-30) * zf ** 2   # template in rcs view
        num = float(np.nansum(wf * t * yf))
        den = float(np.nansum(wf * t * t)) + 1e-30
        fam = models.Family(fam0.itype, "template", np.array([num / den]), dict(fam0.meta))
    return fam, b_free, sem


def _night_step(rng, fit_mask, S, sig, M, Aer, g, b_fam, ok_n, mu_arr, sd_rel):
    """v4-style per-night WLS on (S − b_fam), with a PER-NIGHT amplitude prior (NaN = no prior)."""
    z2 = (rng / 1000.0) ** 2
    n = S.shape[0]
    A = np.full(n, np.nan)
    B = np.zeros(n)
    d = np.zeros(n)
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
        # the anchor: one synthetic observation A_i ~ N(mu_i, sd_rel·mu_i) — nights without an
        # archive constant get no prior row and stay free
        mu = mu_arr[i]
        if np.isfinite(mu):
            sd = sd_rel * mu
            Xa = np.vstack([X * np.sqrt(w[m])[:, None],
                            np.array([[1.0, 0.0, 0.0]]) * (1.0 / sd)])
            ya = np.concatenate([y[m] * np.sqrt(w[m]), [mu / sd]])
        else:
            Xa = X * np.sqrt(w[m])[:, None]
            ya = y[m] * np.sqrt(w[m])
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
    elif itype == "CHM15k":
        fit_mask = np.isfinite(rng) & (rng >= 1000.0)  # below 1 km the hood shows overlap, not dark

    # hood-prior family for the type (fit on the matching Payerne instrument)
    pay_ident = {"CHM15k": "A", "CL31": "B", "CL61": "C"}.get(itype)
    trng, t_brcs, tsem, _tbp = hood.truth(pay_ident)
    fam0 = models.fit_family(itype, trng, t_brcs, tsem)

    # per-night archive constants -> RELATIVE anchor with a one-shot units factor k_hat
    wmo, ident = cache_file.name.split("_")[0], cache_file.name.split("_")[1]
    cmap = _nightly_cref(wmo, ident, itype)
    cref, sd_scale = _era_cref(wmo, ident, itype, dates)   # kept for the diagnostic print only

    hb_inj = np.interp(rng, trng, t_brcs)              # the injection shape = the hood truth
    hb_inj = np.where(np.isfinite(hb_inj), hb_inj, 0.0)

    def _fit_once(S_in):
        """The whole estimator on one S matrix -> (fam, b_free, sem, theta). Used twice: on the
        real nights, and on the same nights + 1x the hood dark (the self-test)."""
        A_, B_, d_, _b0, okn_ = v4.als_fit(S_in, sig, M, Aer, rng, fit_mask, n_iter=4)
        g_ = np.ones(rng.size)
        ds8_ = np.array([str(x)[:8] for x in np.asarray(dates).ravel()])
        C_n_ = np.array([cmap.get(x, np.nan) for x in ds8_])
        with np.errstate(invalid="ignore", divide="ignore"):
            rat_ = A_ / C_n_
        kh_ = float(np.nanmedian(rat_[okn_])) if np.isfinite(rat_[okn_]).any() else np.nan
        mu_ = kh_ * C_n_ if np.isfinite(kh_) else np.full(C_n_.size, np.nan)
        fam_ = fam0
        bf_ = sm_ = None
        n_eff_ = 1 if fam0.kind in ("template", "gexp") else n_outer
        for _ in range(n_eff_):
            fam_, bf_, sm_ = _family_step(itype, rng, fit_mask, S_in, sig, M, Aer,
                                          A_, B_, d_, g_, okn_, fam_, fam_prior=fam0)
            if fam0.kind not in ("template", "gexp") and np.isfinite(kh_):
                A_, B_, d_, okn_ = _night_step(rng, fit_mask, S_in, sig, M, Aer, g_,
                                               fam_.rcs_view(rng), okn_, mu_, sd_rel)
        bh_ = fam_.rcs_view(rng)
        tb_ = (rng >= 3000) & (rng <= 9500) & fit_mask & np.isfinite(bh_) & np.isfinite(hb_inj)
        th_ = (float(np.nansum(bh_[tb_] * hb_inj[tb_])
                     / max(np.nansum(hb_inj[tb_] ** 2), 1e-30)) if tb_.sum() > 10 else np.nan)
        return fam_, bf_, sm_, A_, okn_, th_

    # v4 initial state
    A, B, d, b0, ok_n = v4.als_fit(S, sig, M, Aer, rng, fit_mask, n_iter=4)
    g = np.ones(rng.size)
    fam = fam0
    b_free = sem = None
    ds8 = np.array([str(d)[:8] for d in np.asarray(dates).ravel()])
    C_n = np.array([cmap.get(d, np.nan) for d in ds8])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratios = A / C_n
    k_hat = float(np.nanmedian(ratios[ok_n])) if np.isfinite(ratios[ok_n]).any() else np.nan
    mu_arr = k_hat * C_n if np.isfinite(k_hat) else np.full(C_n.size, np.nan)
    n_pri = int(np.isfinite(mu_arr).sum())
    sd_rel = T2_SD * (sd_scale if np.isfinite(cref) else 1.0)
    print(f"   era cref={cref:.4g}  units factor k_hat={k_hat:.3f}  "
          f"nights with prior: {n_pri}/{S.shape[0]}  "
          f"A_init median={np.nanmedian(A[ok_n]):.4g}")
    # The template family does NOT alternate. Its far-range shape is nearly z^2/molecular-
    # collinear (the CL61 dark is a small quasi-flat P), so every anchored night step drains a bit
    # of it into A_n and d_n — session-2 traces: amplitude 2.8e-9 -> 3.6e-10 -> -3.5e-11. The v4
    # configuration (free amplitudes, single per-gate pass) is exactly what validated theta=0.96,
    # so the template takes its ONE family step on the v4-initialised state and stops.
    # Session-2 conclusion, generalised from the CL61: the anchored alternation DRIFTS for every
    # family whose far-range shape brushes the molecular/z^2 directions -- each night step
    # re-partitions the collinear component and the family walks (CHM15k: A -35 %, corr 0.46->0.24
    # over two extra iterations, with per-night priors in place). The single family step on the
    # v4-initialised state is the stable configuration; the alternation stays available behind
    # --alternate for the injection study.
    n_eff = 1 if fam0.kind in ("template", "gexp") else n_outer
    for it in range(n_eff):
        fam, b_free, sem = _family_step(itype, rng, fit_mask, S, sig, M, Aer,
                                        A, B, d, g, ok_n, fam, fam_prior=fam0)
        if fam0.kind != "template" and np.isfinite(k_hat):
            A, B, d, ok_n = _night_step(rng, fit_mask, S, sig, M, Aer, g,
                                        fam.rcs_view(rng), ok_n, mu_arr, sd_rel)
        bmag = float(np.nanmedian(np.abs(fam.rcs_view(rng)[fit_mask])))
        print(f"   iter {it}: A median={np.nanmedian(A[ok_n]):.4g}  "
              f"|b_fam| median={bmag:.4g}  params={np.array2string(np.asarray(fam.params)[:4], precision=3)}")
    b_hat = fam.rcs_view(rng)

    # --- the self-test: same estimator, same nights, + 1x the hood dark injected ---------------
    fam_i, _bfi, _smi, _Ai, _oki, th_inj = _fit_once(S + hb_inj[None, :])
    _f0, _bf0, _sm0, _A0, _ok0, th_raw = _fit_once(S)
    gain = th_inj - th_raw                    # what a TRUE unit dark adds to the estimate
    th_corr = th_raw / gain if np.isfinite(gain) and abs(gain) > 0.05 else np.nan
    print(f"   selftest: theta_raw={th_raw:+.2f}  injection gain={gain:+.2f}  "
          f"theta_CORR={th_corr:+.2f}  (target 1)")

    out = {"rng": rng, "b_hat": b_hat, "b_free": b_free, "sem": sem,
           "theta_raw": np.array(th_raw), "inj_gain": np.array(gain),
           "theta_corr": np.array(th_corr),
           "params": fam.params, "kind": np.array(fam.kind), "itype": np.array(itype),
           "cref": np.array(cref), "n_nights": np.array(int(ok_n.sum()))}

    # hood comparison for the Payerne streams — the verdict, IN RCS VIEW on b_rcs (the
    # subtractable product; the P-view array cost a z² mix-up once)
    if wmo == "0-20000-0-06610":
        hb = np.interp(rng, trng, t_brcs)
        m = fit_mask & np.isfinite(b_hat) & np.isfinite(hb)
        band = (rng >= 2000) & (rng <= 8000) & m
        nearb = (rng >= 300) & (rng <= 1500) & m
        for name, mm in (("ray 2-8km", band), ("near 0.3-1.5km", nearb)):
            if mm.sum() < 10:
                continue
            cc = float(np.corrcoef(b_hat[mm], hb[mm])[0, 1])
            gain = float(np.nansum(b_hat[mm] * hb[mm]) / max(np.nansum(hb[mm] ** 2), 1e-30))
            print(f"   vs hood [{name}]: corr {cc:+.2f}  amplitude {gain:+.2f} (1 = perfect)")
        # the v4-equivalent theta: project the retrieved profile on the hood shape over the v4
        # band (3-9.5 km). theta = 1 is a perfect amplitude on the hood template.
        tband = (rng >= 3000) & (rng <= 9500) & m
        if tband.sum() > 10:
            th = float(np.nansum(b_hat[tband] * hb[tband])
                       / max(np.nansum(hb[tband] ** 2), 1e-30))
            print(f"   theta (3-9.5 km, v4 convention): {th:+.2f}  (target 1)")
            out["theta"] = np.array(th)
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
