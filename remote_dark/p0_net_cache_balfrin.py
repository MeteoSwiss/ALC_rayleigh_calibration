# -*- coding: utf-8 -*-
"""Runs ON balfrin — P0 of the CHM15k dark plan (report 18): build v4 night caches for every
network CHM15k over 2025-01-01..2026-08-13 and compute the per-unit dark-impact BOUND.

Per stream:
  1. clear nights from the v2.2 release CSVs, cache S/sigma/M/Aer per night (restartable:
     existing npz is reused).
  2. the evidence intercept d(z): per-gate OLS of the raw night signal S_n(z) against that
     night's PREDICTED atmospheric signal x_n(z) = C_n * M_n(z)  (C_n = the night's constant
     in the stable 4.5-6.5 km window). The transmission variation between nights is the
     leverage; night-normalising instead destroys identifiability (report 16 lesson).
  3. the conservative dark-impact bound: refit the nightly constants with d(z) subtracted;
     bound = |median C' / median C - 1|. d(z) contains dark + mean aerosol residue (~x4 the
     dark at Payerne), so the bound OVERSTATES, never hides, a large dark — that is what a
     triage bound must do. Calibration against Payerne A's known +24 % happens locally.

Outputs on /scratch/mch/mhrvo/remote_dark_net/: per-stream cache npz, p0_bounds.csv,
p0_evidence.npz (per-stream d(z) + rng, for local shape analysis).
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "alc_v22_code"))

from rayleigh_availability import dark_from_clearsky as v4          # noqa: E402

L1 = "/scratch/mch/mhrvo/E-PROFILE_L1"
CAMS = "/scratch/mch/mhrvo/CAMS_Monthly_04;/scratch/mch/mhrvo/CAMS"
CSV = Path("/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel")
OUT = Path("/scratch/mch/mhrvo/remote_dark_net")
ERA = ("20250101", "20260813")
CWIN = (4500.0, 6500.0)
WORKERS = 32
MIN_NIGHTS = 20


def cache_stream(wmo, ident):
    f = OUT / f"{wmo}_{ident}.npz"
    if f.exists():
        return np.load(f, allow_pickle=True)
    dates = v4.clear_nights_from_csv(CSV, wmo, ident, *ERA)
    if len(dates) < MIN_NIGHTS:
        return None
    jobs = [dict(date=d, wmo=wmo, ident=ident, l1_root=L1, cams_folder=CAMS) for d in dates]
    nights = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for res in ex.map(v4.process_night, jobs, chunksize=1):
            if res is not None:
                nights.append(res)
    if len(nights) < MIN_NIGHTS:
        return None
    rng = nights[0]["rng"]
    keep = [n for n in nights if n["rng"].size == rng.size]
    np.savez(f,
             S=np.vstack([n["S"] for n in keep]),
             sigma=np.vstack([n["sigma"] for n in keep]),
             M=np.vstack([n["m"] for n in keep]),
             Aer=np.vstack([n["a"] for n in keep]),
             dates=np.array([int(n["date"]) for n in keep]),
             n_prof=np.array([n["n_prof"] for n in keep]),
             rng=rng, itype=np.array("CHM15k"))
    return np.load(f, allow_pickle=True)


def nightly_c(S, sig, M, rng):
    m = (rng >= CWIN[0]) & (rng <= CWIN[1])
    w = 1.0 / np.maximum(sig[:, m], 1e-300) ** 2
    ok = np.isfinite(S[:, m]) & np.isfinite(M[:, m])
    w = np.where(ok, w, 0.0)
    S0, M0 = np.nan_to_num(S[:, m]), np.nan_to_num(M[:, m])
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.nansum(w * S0 * M0, axis=1)
                / np.maximum(np.nansum(w * M0 * M0, axis=1), 1e-300))


def evidence_and_bound(z):
    S, sig, M, rng = z["S"], z["sigma"], z["M"], z["rng"]
    c = nightly_c(S, sig, M, rng)
    ok = np.isfinite(c) & (c > 0)
    if ok.sum() < MIN_NIGHTS:
        return None
    S, sig, M = S[ok], sig[ok], M[ok]
    c = c[ok]
    x = c[:, None] * M
    n = np.isfinite(S) & np.isfinite(x)
    cnt = n.sum(axis=0).astype(float)
    Sm, Xm = np.where(n, S, 0.0), np.where(n, x, 0.0)
    sx, ss = Xm.sum(axis=0), Sm.sum(axis=0)
    sxx, sxs = (Xm * Xm).sum(axis=0), (Xm * Sm).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        den = cnt * sxx - sx * sx
        d = (sxx * ss - sx * sxs) / den
    d[cnt < 8] = np.nan
    d0 = np.nan_to_num(d)
    c2 = nightly_c(S - d0[None, :], sig, M, rng)
    ok2 = np.isfinite(c2) & (c2 > 0)
    # the reference is the raw constant of the SAME nights
    bound = (float(np.median(c2[ok2]) / np.median(c[ok2]) - 1.0)
             if ok2.sum() >= MIN_NIGHTS else np.nan)
    return dict(d=d, rng=rng, n_nights=int(ok.sum()), c_med=float(np.median(c)),
                bound=bound)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    streams = json.load(open(Path.home() / "chm_streams.json"))
    rows = []
    ev = {}
    for i, s in enumerate(streams):
        wmo, ident = s["wmo"], s["ident"]
        key = f"{wmo}_{ident}"
        try:
            z = cache_stream(wmo, ident)
        except Exception as e:                                              # noqa: BLE001
            print(f"[{i + 1}/{len(streams)}] {key}: CACHE ERROR {e}", flush=True)
            continue
        if z is None:
            print(f"[{i + 1}/{len(streams)}] {key}: insufficient nights", flush=True)
            rows.append(f"{wmo},{ident},{s['site']},0,,")
            continue
        r = evidence_and_bound(z)
        if r is None:
            print(f"[{i + 1}/{len(streams)}] {key}: bound not computable", flush=True)
            rows.append(f"{wmo},{ident},{s['site']},{z['S'].shape[0]},,")
            continue
        ev[f"d_{key}"] = r["d"]
        ev[f"rng_{key}"] = r["rng"]
        rows.append(f"{wmo},{ident},{s['site']},{r['n_nights']},{r['c_med']:.6e},"
                    f"{100 * r['bound']:+.2f}")
        print(f"[{i + 1}/{len(streams)}] {key}: {r['n_nights']} nights, "
              f"bound {100 * r['bound']:+.1f} %", flush=True)
    with open(OUT / "p0_bounds.csv", "w") as fh:
        fh.write("wmo,ident,site,n_nights,c_med,bound_pct\n")
        fh.write("\n".join(rows) + "\n")
    np.savez_compressed(OUT / "p0_evidence.npz", **ev)
    print("P0_DONE", flush=True)


if __name__ == "__main__":
    main()
