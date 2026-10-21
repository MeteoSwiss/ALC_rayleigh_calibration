# -*- coding: utf-8 -*-
"""Top-down decomposition of the CHM15k calibration-constant variability (report 19, q15 agent f).

The other q15 agents build the budget bottom-up, one physical factor at a time. This script tests the
budget top-down, empirically:

1. Common-mode vs instrument split with co-located units. Atmosphere and calibration method (aerosol in
   the window, molecular model, transmission) are common to co-located units and cancel in the
   same-night ratio C_A/C_B; instrumental factors do not. Annual harmonic of ln C per unit and of
   ln(C_A/C_B) on common nights, before and after the kappa normalisation; correlation of the synoptic
   anomalies between co-located units (= common-mode share of the night-to-night variance).
   Same with CHM15k vs CL61 (different wavelength: weaker test).
2. Network-wide multivariate regression per unit era (main optical module), ln C nightly against
   housekeeping, calibration-window and aerosol covariates, in two regimes:
   (a) synoptic anomalies (31-day running median removed): group dominance analysis (Shapley R^2);
   (b) long term: (b1) annual harmonic of each covariate x synoptic coefficient (transfer model),
       projected on the C_L phasor; (b2) pooled monthly means with unit fixed effects, dominance analysis
       with the calendar harmonic as a competing group.
3. Budget of the annual amplitude and of the night-to-night scatter.
4. Is the residual annual cycle after kappa normalisation structured by geography/climate (atmosphere)
   or by module generation (instrument)?

Run:  .venv/Scripts/python.exe validation/chm15k_stability/q15f_decomposition.py [--rebuild]
Outputs: C:/DATA/Projects/202610_CHM15k_stability/q15/q15f_*
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
import warnings
from datetime import datetime, timedelta
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
from common import (BLOCKS_PER_DAY, CENSUS, OUT, PAIR_SITES, chm15k_streams, load_cache,  # noqa: E402
                    load_constants, stream_key)
from q1_temperature import ANOM_WINDOW, night_aerosol_contrast, running_median  # noqa: E402
from validation.paper.kalman_improve_alc import flag_outliers  # noqa: E402

warnings.filterwarnings("ignore")
Q15 = OUT / "q15"
EPOCH = datetime(1970, 1, 1)
W = 2 * np.pi / 365.25
PEAK_REF = 176.0                       # network peak DOY of C_L (brief)
E_REF = np.exp(1j * W * PEAK_REF)      # unit phasor of the reference phase
PAYERNE = "0-20000-0-06610_A"
MIN_NIGHTS = 100
AERONET_DIR = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-"
                   "calibration/e6120125-e8c9-4e0b-8352-26777f3067de/scratchpad/aeronet_raw")

# Covariates per night and the factor groups they belong to.
GROUPS = {
    "instrument temperature": ["T_int", "T_lom", "T_det"],
    "receiver (p_calc, dark, noise)": ["ln_pcal", "ln_bg", "ln_sd"],
    "laser (pulse count)": ["ln_pulses"],
    "window transmission": ["wt"],
    "calibration window / sampling": ["zb_km", "zt_km", "dur_h", "ln_nclear"],
    "aerosol in/below window": ["aer"],
}
COVS = [v for g in GROUPS.values() for v in g]
HK_MAP = {"T_int": "temp_int", "T_lom": "temperature_optical_module", "T_det": "temperature_detector",
          "pcal": "calibration_pulse", "wt": "window_transmission", "bg": "bckgrd_rcs_0", "sd": "stddev",
          "pulses": "laser_pulses", "life": "laser_life_time", "T_ext": "temp_ext"}


# ---------------------------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------------------------
def night_key(t0):
    """Integer night index: the night starting on calendar day D (noon D to noon D+1) -> D."""
    return np.floor(np.asarray(t0) - 0.5).astype(int)


def harmonic(t, y):
    """Fit y = m + a cos(w doy) + b sin(w doy). Returns (phasor a+ib, mean, residuals)."""
    ok = np.isfinite(t) & np.isfinite(y)
    d = np.mod(t[ok], 365.25)
    X = np.column_stack([np.ones(ok.sum()), np.cos(W * d), np.sin(W * d)])
    c, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    res = np.full(len(y), np.nan)
    d_all = np.mod(t, 365.25)
    res[ok] = y[ok] - X @ c
    fit = c[0] + c[1] * np.cos(W * d_all) + c[2] * np.sin(W * d_all)
    return c[1] + 1j * c[2], c[0], res, fit


def peak_doy(p):
    return float(np.mod(np.angle(p) / W, 365.25))


def inphase(p):
    """Component of a phasor in phase with the network C_L peak (DOY 176)."""
    return float(np.real(p * np.conj(E_REF)))


def r2_ols(X, y):
    if X.shape[1] == 0:
        return 0.0
    A = np.column_stack([np.ones(len(y)), X])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ c
    return 1.0 - res @ res / np.sum((y - y.mean()) ** 2)


def dominance(blocks, y):
    """General dominance (Shapley value of R^2) of each block of predictors. Sum = full R^2."""
    G = len(blocks)
    r2 = {}
    for mask in range(2 ** G):
        cols = [blocks[g] for g in range(G) if mask >> g & 1]
        r2[mask] = r2_ols(np.column_stack(cols) if cols else np.empty((len(y), 0)), y)
    dom = np.zeros(G)
    for g in range(G):
        by_size = []
        for k in range(G):
            inc = [r2[m | (1 << g)] - r2[m] for m in range(2 ** G)
                   if not (m >> g & 1) and bin(m).count("1") == k]
            by_size.append(np.mean(inc))
        dom[g] = np.mean(by_size)
    full = r2[2 ** G - 1]
    unique = np.array([full - r2[(2 ** G - 1) & ~(1 << g)] for g in range(G)])
    alone = np.array([r2[1 << g] for g in range(G)])
    return dom, full, unique, alone


def anomaly(t, x):
    return x - running_median(t, x, ANOM_WINDOW)


# ---------------------------------------------------------------------------------------------
# Stage 1: nightly table (one row per calibration night of each unit's main optical module)
# ---------------------------------------------------------------------------------------------
def window_stats(cache, t0, t1):
    """Median of every HK variable and the sum of clear / all profiles over each night's window."""
    tb = (cache["day"][:, None] + (np.arange(BLOCKS_PER_DAY)[None, :] + 0.5) / BLOCKS_PER_DAY).ravel()
    nv = cache["hk"].shape[-1]
    hk = cache["hk"].reshape(-1, nv)
    nc = cache["n_clear"].ravel().astype(float)
    na = cache["n_all"].ravel().astype(float)
    order = np.argsort(tb)
    tb, hk, nc, na = tb[order], hk[order], nc[order], na[order]
    med = np.full((len(t0), nv), np.nan)
    sc = np.full(len(t0), np.nan)
    sa = np.full(len(t0), np.nan)
    i0 = np.searchsorted(tb, t0)
    i1 = np.searchsorted(tb, t1, side="right")
    for k in range(len(t0)):
        if i1[k] > i0[k]:
            med[k] = np.nanmedian(hk[i0[k]:i1[k]], axis=0)
            sc[k] = np.nansum(nc[i0[k]:i1[k]])
            sa[k] = np.nansum(na[i0[k]:i1[k]])
    return med, sc, sa


def build_unit(args):
    wmo, ident, site, lat, lon, alt, aer_daily = args
    key = stream_key(wmo, ident)
    cons, cache = load_constants(wmo, ident), load_cache(wmo, ident)
    if cons is None or cache is None:
        return None
    mod = max(set(cons["module"]), key=list(cons["module"]).count)
    m = (cons["module"] == mod) & (cons["c"] > 0) & np.isfinite(cons["c"])
    if m.sum() < MIN_NIGHTS:
        return None
    cons = {k: v[m] for k, v in cons.items()}
    out_flag = flag_outliers(cons["c"])           # the Kalman's global IQR screen (per module era)
    cons = {k: v[~out_flag] for k, v in cons.items()}
    med, sc, sa = window_stats(cache, cons["t0"], cons["t1"])
    names = list(cache["hk_vars"])

    def hk(v):
        return med[:, names.index(HK_MAP[v])] if HK_MAP[v] in names else np.full(len(cons["t"]), np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        rows = dict(
            t=cons["t0"] * 0.5 + cons["t1"] * 0.5, night=night_key(cons["t0"]), lnC=np.log(cons["c"]),
            rel_u=cons["u"] / cons["c"],
            T_int=hk("T_int"), T_lom=hk("T_lom"), T_det=hk("T_det"), T_ext=hk("T_ext"),
            ln_pcal=np.log(hk("pcal")), wt=hk("wt"), ln_bg=np.log(hk("bg")), ln_sd=np.log(hk("sd")),
            ln_pulses=np.log(hk("pulses")), life_kh=hk("life") / 1000.0,
            zb_km=(cons["zb"] - alt) / 1000.0, zt_km=(cons["zt"] - alt) / 1000.0,
            dur_h=(cons["t1"] - cons["t0"]) * 24.0, ln_nclear=np.log(np.where(sc > 0, sc, np.nan)),
            clear_frac=sc / np.where(sa > 0, sa, np.nan),
            aer=night_aerosol_contrast(cache, cons, alt))
    aod = np.full(len(cons["t"]), np.nan)
    ang = np.full(len(cons["t"]), np.nan)
    if aer_daily:
        for k, (t0, t1) in enumerate(zip(cons["t0"], cons["t1"])):
            vals = [aer_daily[x] for x in {(EPOCH + timedelta(days=float(t0))).date(),
                                           (EPOCH + timedelta(days=float(t1))).date()} if x in aer_daily]
            if vals:
                aod[k] = np.mean([v[0] for v in vals])
                ang[k] = np.mean([v[2] for v in vals])
    rows["aod"] = aod
    rows["ang"] = ang
    order = np.argsort(rows["t"])
    rows = {k: v[order] for k, v in rows.items()}
    meta = dict(stream=key, site=site, module=mod, lat=lat, lon=lon, alt=alt)
    return meta, rows


def build_table(rebuild=False):
    f = Q15 / "q15f_nightly.npz"
    if f.exists() and not rebuild:
        z = np.load(f, allow_pickle=True)
        return list(z["units"])
    alts = {(s["wmo"], s["ident"]): float(s["alt"]) for s in json.load(open(CENSUS, encoding="utf-8"))}
    from q5_aeronet import read_aeronet_daily
    match = {k: ae for k, _, ae, _ in json.load(open(AERONET_DIR.parent / "aeronet_match.json", encoding="utf-8"))}
    jobs = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = stream_key(wmo, ident)
        daily = None
        if key in match and (AERONET_DIR / f"AOD15_daily_{match[key]}.txt").exists():
            daily = read_aeronet_daily(AERONET_DIR / f"AOD15_daily_{match[key]}.txt")
        jobs.append((wmo, ident, site, lat, lon, alts.get((wmo, ident), 0.0), daily))
    with Pool(4) as pool:
        res = [r for r in pool.map(build_unit, jobs) if r is not None]
    np.savez_compressed(f, units=np.array(res, dtype=object))
    # flat CSV of every night (for other agents)
    with open(Q15 / "q15f_nightly.csv", "w", newline="", encoding="utf-8") as fh:
        cols = list(res[0][1].keys())
        w = csv.writer(fh)
        w.writerow(["stream", "module"] + cols)
        for meta, rows in res:
            for i in range(len(rows["t"])):
                w.writerow([meta["stream"], meta["module"]] + [f"{rows[c][i]:.6g}" for c in cols])
    print(f"nightly table: {len(res)} units -> {f}")
    return res


# ---------------------------------------------------------------------------------------------
# kappa per unit
# ---------------------------------------------------------------------------------------------
def load_kappas():
    k_in = {r["stream"]: float(r["b_T_ctrl"]) for r in csv.DictReader(open(OUT / "q2a_intranight.csv", encoding="utf-8"))}
    k_m2 = {}
    for r in csv.DictReader(open(OUT / "q1_slopes.csv", encoding="utf-8")):
        if r["arm"] == "A" and r["tvar"] == "temp_int" and r["b2"] not in ("", "nan"):
            k_m2[(r["stream"], r["era"])] = float(r["b2"])
    return k_in, k_m2


def kappa_norm(rows, kappa):
    """ln C normalised to the unit's median T_int with kappa [%/K]."""
    T = rows["T_int"]
    return rows["lnC"] - kappa / 100.0 * (T - np.nanmedian(T))


# ---------------------------------------------------------------------------------------------
# Task 1: co-located pairs (same wavelength) and CHM15k vs CL61
# ---------------------------------------------------------------------------------------------
def pair_regress(ra, rb, ia, ib):
    """Synoptic anomalies on common nights: common part (A+B)/2 on the mean covariate anomalies and
    differential part A-B on the covariate differences. Adjusted R^2 and group dominance."""
    t = 0.5 * (ra["t"][ia] + rb["t"][ib])
    o = np.argsort(t)
    t = t[o]
    la = anomaly(t, ra["lnC"][ia][o])
    lb = anomaly(t, rb["lnC"][ib][o])
    out = {}
    for lab, y, comb in (("common", 0.5 * (la + lb), lambda a, b: 0.5 * (a + b)), ("diff", la - lb, lambda a, b: a - b)):
        cols, names = [], []
        for v in COVS:
            xa = anomaly(t, np.where(np.isfinite(ra[v][ia][o]), ra[v][ia][o], np.nan))
            xb = anomaly(t, np.where(np.isfinite(rb[v][ib][o]), rb[v][ib][o], np.nan))
            x = comb(xa, xb)
            if np.mean(np.isfinite(x)) < 0.7:
                continue
            cols.append(x)
            names.append(v)
        X = np.column_stack(cols)
        ok = np.isfinite(y) & np.all(np.isfinite(X), 1)
        X, yy = X[ok], y[ok]
        keep = X.std(0) > 1e-9
        X = X[:, keep]
        names = [v for v, k in zip(names, keep) if k]
        Z = (X - X.mean(0)) / X.std(0)
        groups = [(g, [names.index(v) for v in vs if v in names]) for g, vs in GROUPS.items()]
        groups = [(g, ix) for g, ix in groups if ix]
        dom, full, _, _ = dominance([Z[:, ix] for _, ix in groups], yy)
        n, k = Z.shape
        adj = 1 - (1 - full) * (n - 1) / max(n - k - 1, 1)
        out[f"{lab}_sd"] = 100 * np.std(yy)
        out[f"{lab}_r2adj"] = adj
        for (g, _), d in zip(groups, dom):
            out[f"{lab}_dom|{g}"] = d * max(adj, 0) / full if full > 0 else 0.0
    return out


def task1_pairs(units, k_in, k_m2):
    by = {m["stream"]: (m, r) for m, r in units}
    unit_rows, pair_rows, site_rows = [], [], []
    for wmo, idents in PAIR_SITES.items():
        members = [f"{wmo}_{i}" for i in idents if f"{wmo}_{i}" in by]
        phasors = {}
        for s in members:
            m, r = by[s]
            ok = np.isfinite(r["T_int"])
            t = r["t"][ok]
            series = {"op": r["lnC"][ok]}
            km = k_m2.get((s, m["module"]), np.nan)
            if np.isfinite(km):
                series["k_m2"] = kappa_norm({k: v[ok] for k, v in r.items()}, km)
            if s in k_in:
                series["k_in"] = kappa_norm({k: v[ok] for k, v in r.items()}, k_in[s])
            pT, _, _, _ = harmonic(t, r["T_int"][ok])
            phasors[s] = {}
            for lab, y in series.items():
                p, _, res, _ = harmonic(t, y)
                phasors[s][lab] = p
                unit_rows.append(dict(site_wmo=wmo, stream=s, module=m["module"], variant=lab, n=int(ok.sum()),
                                      amp=100 * abs(p), peak=peak_doy(p), amp_inphase=100 * inphase(p),
                                      scat=100 * np.nanstd(res), T_amp=abs(pT), T_peak=peak_doy(pT),
                                      kappa_m2=km, kappa_in=k_in.get(s, np.nan)))
        # pairs on common nights
        for a, b in itertools.combinations(members, 2):
            (ma, ra), (mb, rb) = by[a], by[b]
            common, ia, ib = np.intersect1d(ra["night"], rb["night"], return_indices=True)
            if len(common) < 60:
                continue
            t = 0.5 * (ra["t"][ia] + rb["t"][ib])
            for lab in ("op", "k_m2", "k_in"):
                if lab == "op":
                    ya, yb = ra["lnC"][ia], rb["lnC"][ib]
                else:
                    ka = k_m2.get((a, ma["module"]), np.nan) if lab == "k_m2" else k_in.get(a, np.nan)
                    kb = k_m2.get((b, mb["module"]), np.nan) if lab == "k_m2" else k_in.get(b, np.nan)
                    if not (np.isfinite(ka) and np.isfinite(kb)):
                        continue
                    ya, yb = kappa_norm(ra, ka)[ia], kappa_norm(rb, kb)[ib]
                ok = np.isfinite(ya) & np.isfinite(yb)
                tt, ya, yb = t[ok], ya[ok], yb[ok]
                pa, _, resa, _ = harmonic(tt, ya)
                pb, _, resb, _ = harmonic(tt, yb)
                pr, _, resr, _ = harmonic(tt, ya - yb)
                la, lb = anomaly(tt, ya), anomaly(tt, yb)
                g = np.isfinite(la) & np.isfinite(lb)
                cov_res = np.cov(resa, resb)
                cov_an = np.cov(la[g], lb[g])
                pair_rows.append(dict(
                    pair=f"{a} / {b.rsplit('_', 1)[1]}", variant=lab, n_common=int(ok.sum()),
                    amp_A=100 * abs(pa), amp_B=100 * abs(pb), peak_A=peak_doy(pa), peak_B=peak_doy(pb),
                    amp_ratio=100 * abs(pr), peak_ratio=peak_doy(pr), inphase_ratio=100 * inphase(pr),
                    # share of the annual amplitude that is common (shared phasor = (pa+pb)/2 projected)
                    common_amp=100 * abs(0.5 * (pa + pb)),
                    frac_amp_differential=abs(pr) / (0.5 * (abs(pa) + abs(pb))),
                    scat_A=100 * np.std(resa), scat_B=100 * np.std(resb), scat_ratio=100 * np.std(resr),
                    r_resid=float(np.corrcoef(resa, resb)[0, 1]),
                    common_var_share_resid=float(cov_res[0, 1] / np.sqrt(cov_res[0, 0] * cov_res[1, 1])),
                    sd_anom_A=100 * np.sqrt(cov_an[0, 0]), sd_anom_B=100 * np.sqrt(cov_an[1, 1]),
                    sd_anom_ratio=100 * np.nanstd(la[g] - lb[g]),
                    r_anom=float(cov_an[0, 1] / np.sqrt(cov_an[0, 0] * cov_an[1, 1])),
                    rel_u_A=100 * np.nanmedian(ra["rel_u"][ia][ok]), rel_u_B=100 * np.nanmedian(rb["rel_u"][ib][ok]),
                    var_common=1e4 * cov_an[0, 1], var_A=1e4 * cov_an[0, 0], var_B=1e4 * cov_an[1, 1]))
                if lab == "op":
                    pair_rows[-1].update(pair_regress(ra, rb, ia, ib))
        # site-level common phasor (mean over members) vs unit-specific spread
        for lab in ("op", "k_m2", "k_in"):
            ps = np.array([phasors[s][lab] for s in members if lab in phasors[s]])
            if len(ps) < 2:
                continue
            mean = ps.mean()
            spread = np.sqrt(np.mean(np.abs(ps - mean) ** 2))
            site_rows.append(dict(site_wmo=wmo, variant=lab, n_units=len(ps), common_amp=100 * abs(mean),
                                  common_peak=peak_doy(mean), common_inphase=100 * inphase(mean),
                                  specific_rms=100 * spread, mean_unit_amp=100 * np.mean(np.abs(ps))))
    return unit_rows, pair_rows, site_rows


def task1_cl61(by=None, k_m2=None):
    """CHM15k vs CL61 constants on common nights (different wavelength: weak common-mode test)."""
    sites = {"Payerne": ("0-20000-0-06610", "A", "C"), "Lindenberg": ("0-20000-0-10393", "0", "C"),
             "Camborne": ("0-20000-0-03808", "A", "C"), "Aosta": ("0-380-5-1", "0", "B")}
    rows = []
    for site, (wmo, ich, icl) in sites.items():
        ch, cl = load_constants(wmo, ich), load_constants(wmo, icl)
        if ch is None or cl is None:
            rows.append(dict(site=site, n_common=0, note="no constants"))
            continue
        mod = max(set(ch["module"]), key=list(ch["module"]).count)
        mc = (ch["module"] == mod) & (ch["c"] > 0)
        ch = {k: v[mc] for k, v in ch.items()}
        cl = {k: v[cl["c"] > 0] for k, v in cl.items()}
        ch = {k: v[~flag_outliers(ch["c"])] for k, v in ch.items()}
        cl = {k: v[~flag_outliers(cl["c"])] for k, v in cl.items()}
        common, ia, ib = np.intersect1d(night_key(ch["t0"]), night_key(cl["t0"]), return_indices=True)
        if len(common) < 40:
            rows.append(dict(site=site, n_common=len(common), note="too few common nights"))
            continue
        t = 0.5 * (ch["t0"][ia] + ch["t1"][ia])
        ya, yb = np.log(ch["c"][ia]), np.log(cl["c"][ib])
        o = np.argsort(t)
        t, ya, yb = t[o], ya[o], yb[o]
        pa, _, resa, _ = harmonic(t, ya)
        pb, _, resb, _ = harmonic(t, yb)
        pr, _, resr, _ = harmonic(t, ya - yb)
        la, lb = anomaly(t, ya), anomaly(t, yb)
        g = np.isfinite(la) & np.isfinite(lb)
        span = t.max() - t.min()
        # CHM15k kappa-normalised (multi-day kappa) on the same nights, from the nightly table
        extra = {}
        key = f"{wmo}_{ich}"
        if by is not None and key in by:
            meta, r = by[key]
            km = k_m2.get((key, meta["module"]), np.nan)
            nk = night_key(cl["t0"][ib][o])
            common2, j1, j2 = np.intersect1d(r["night"], nk, return_indices=True)
            if np.isfinite(km) and len(common2) >= 40:
                yk = kappa_norm(r, km)[j1]
                ok2 = np.isfinite(yk)
                t2 = r["t"][j1][ok2]
                pk, _, _, _ = harmonic(t2, yk[ok2])
                prk, _, _, _ = harmonic(t2, yk[ok2] - yb[j2][ok2])
                extra = dict(kappa_m2=km, amp_CHM_k=100 * abs(pk), peak_CHM_k=peak_doy(pk),
                             amp_ratio_k=100 * abs(prk), peak_ratio_k=peak_doy(prk))
        rows.append(dict(site=site, n_common=len(common), span_d=span, amp_CHM=100 * abs(pa), peak_CHM=peak_doy(pa),
                         amp_CL61=100 * abs(pb), peak_CL61=peak_doy(pb), amp_ratio=100 * abs(pr),
                         peak_ratio=peak_doy(pr), r_resid=float(np.corrcoef(resa, resb)[0, 1]),
                         r_anom=float(np.corrcoef(la[g], lb[g])[0, 1]) if g.sum() > 20 else np.nan,
                         sd_anom_CHM=100 * np.std(la[g]), sd_anom_CL61=100 * np.std(lb[g]),
                         note="span < 330 d: harmonic poorly constrained" if span < 330 else "", **extra))
    return rows


# ---------------------------------------------------------------------------------------------
# Task 2/3: network regression, regime (a) synoptic and (b1) transfer to the annual harmonic
# ---------------------------------------------------------------------------------------------
def prepare_design(rows, covs, extra=()):
    """Anomaly design matrix (31-day running median removed) for the requested covariates."""
    t = rows["t"]
    y = anomaly(t, rows["lnC"])
    use, X = [], []
    for v in list(covs) + list(extra):
        x = np.where(np.isfinite(rows[v]), rows[v], np.nan)
        if np.mean(np.isfinite(x)) < 0.7:
            continue
        xa = anomaly(t, x)
        if np.nanstd(xa) < 1e-9:
            continue
        use.append(v)
        X.append(xa)
    if not use:
        return None
    X = np.column_stack(X)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[ok], y[ok]
    keep = X.std(0) > 1e-9 if len(y) else np.zeros(X.shape[1], bool)
    use = [v for v, k in zip(use, keep) if k]
    return use, X[:, keep], y, ok


def unit_regression(meta, rows, kappa_m2):
    t = rows["t"]
    lnC = rows["lnC"]
    d = prepare_design(rows, COVS)
    if d is None or d[1].shape[0] < 60:
        return None
    use, X, y, ok = d
    n, p = X.shape
    mu, sd = X.mean(0), X.std(0)
    Z = (X - mu) / sd
    groups = [(g, [use.index(v) for v in vs if v in use]) for g, vs in GROUPS.items()]
    groups = [(g, ix) for g, ix in groups if ix]
    dom, full, uniq, alone = dominance([Z[:, ix] for _, ix in groups], y)
    adj = 1 - (1 - full) * (n - 1) / max(n - p - 1, 1)
    scale = max(adj, 0) / full if full > 0 else 0.0
    A = np.column_stack([np.ones(n), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    beta = dict(zip(use, coef[1:]))                       # per original unit (anomaly regression)
    beta_std = dict(zip(use, coef[1:] * sd / np.std(y)))  # standardised coefficients

    # Variance budget of the residual around the annual harmonic
    p_c, m_c, res_h, fit_h = harmonic(t, lnC)
    slow = running_median(t, lnC, ANOM_WINDOW) - fit_h
    an = anomaly(t, lnC)
    out = dict(stream=meta["stream"], site=meta["site"], module=meta["module"], lat=meta["lat"], lon=meta["lon"],
               alt=meta["alt"], n=int(len(lnC)), n_reg=n, n_cov=p,
               amp=100 * abs(p_c), peak=peak_doy(p_c), amp_inphase=100 * inphase(p_c),
               scat_tot=100 * np.nanstd(res_h), scat_slow=100 * np.nanstd(slow), scat_syn=100 * np.nanstd(an),
               rel_u=100 * np.nanmedian(rows["rel_u"]), r2_syn=full, r2adj_syn=adj,
               kappa_syn=100 * beta.get("T_int", np.nan), kappa_m2=kappa_m2)
    for v in COVS:
        out[f"bstd|{v}"] = beta_std.get(v, np.nan)
        out[f"b|{v}"] = beta.get(v, np.nan)
    for (g, _), dv, uv, av in zip(groups, dom, uniq, alone):
        out[f"dom|{g}"] = dv * scale
        out[f"uniq|{g}"] = uv
        out[f"alone|{g}"] = av
    # (b1) transfer: annual harmonic of each covariate x synoptic coefficient
    pred_total = 0j
    for g, vs in GROUPS.items():
        pg = 0j
        for v in vs:
            if v in beta and np.mean(np.isfinite(rows[v])) > 0.7:
                pv, _, _, _ = harmonic(t, rows[v])
                pg += beta[v] * pv
                out[f"trv|{v}"] = 100 * float(np.real(beta[v] * pv * np.conj(p_c)) / abs(p_c))
                out[f"harm_inphase|{v}"] = inphase(pv)
        out[f"tr_amp|{g}"] = 100 * abs(pg)
        out[f"tr_proj|{g}"] = 100 * float(np.real(pg * np.conj(p_c)) / abs(p_c))   # in phase with this unit's C
        out[f"tr_peak|{g}"] = peak_doy(pg) if abs(pg) > 0 else np.nan
        pred_total += pg
    out["tr_amp_total"] = 100 * abs(pred_total)
    out["tr_proj_total"] = 100 * float(np.real(pred_total * np.conj(p_c)) / abs(p_c))
    out["amp_resid_tr"] = 100 * abs(p_c - pred_total)
    out["peak_resid_tr"] = peak_doy(p_c - pred_total)
    out["inphase_resid_tr"] = 100 * inphase(p_c - pred_total)
    # kappa normalisation (multi-day kappa) -> residual cycle (task 4)
    T = rows["T_int"]
    if np.isfinite(kappa_m2):
        okT = np.isfinite(T)
        pk, _, resk, _ = harmonic(t[okT], kappa_norm(rows, kappa_m2)[okT])
        out.update(amp_k=100 * abs(pk), peak_k=peak_doy(pk), inphase_k=100 * inphase(pk), scat_k=100 * np.nanstd(resk))
    # trend (ageing) with the harmonic, and the change it makes to the amplitude
    dd = np.mod(t, 365.25)
    Xt = np.column_stack([np.ones(len(t)), np.cos(W * dd), np.sin(W * dd), (t - t.mean()) / 365.25])
    ct, *_ = np.linalg.lstsq(Xt, lnC, rcond=None)
    out.update(trend_pct_yr=100 * ct[3], amp_with_trend=100 * abs(ct[1] + 1j * ct[2]), span_d=float(t.max() - t.min()))
    # climate descriptors for task 4: annual amplitudes of night T_ext and T_int, mean T_ext
    for v in ("T_ext", "T_int"):
        okv = np.isfinite(rows[v])
        if okv.sum() > 60:
            pv, mv, _, _ = harmonic(t[okv], rows[v][okv])
            out[f"{v}_amp"] = abs(pv)
            out[f"{v}_mean"] = mv
            out[f"{v}_peak"] = peak_doy(pv)
    out["zb_med"] = float(np.nanmedian(rows["zb_km"]))
    out["aer_amp_inphase"] = inphase(harmonic(t, rows["aer"])[0]) if np.mean(np.isfinite(rows["aer"])) > 0.7 else np.nan
    # AERONET extension: does AOD add explained synoptic variance?
    if np.mean(np.isfinite(rows["aod"])) > 0.3:
        sub = {k: v[np.isfinite(rows["aod"])] for k, v in rows.items()}
        d2 = prepare_design(sub, COVS, extra=("aod", "ang"))
        if d2 is not None and d2[1].shape[0] >= 40 and "aod" in d2[0]:
            use2, X2, y2, _ = d2
            Z2 = (X2 - X2.mean(0)) / X2.std(0)
            gr2 = [(g, [use2.index(v) for v in vs if v in use2]) for g, vs in GROUPS.items()]
            gr2 = [(g, ix) for g, ix in gr2 if ix] + [("AERONET AOD", [use2.index(v) for v in ("aod", "ang") if v in use2])]
            dom2, full2, uniq2, _ = dominance([Z2[:, ix] for _, ix in gr2], y2)
            n2, p2 = X2.shape
            adj2 = 1 - (1 - full2) * (n2 - 1) / max(n2 - p2 - 1, 1)
            out.update(n_aod=n2, r2adj_aod=adj2, dom_aod=dom2[-1] * max(adj2, 0) / full2, uniq_aod=uniq2[-1],
                       dom_aer_with_aod=dict(zip([g for g, _ in gr2], dom2))["aerosol in/below window"] * max(adj2, 0) / full2
                       if "aerosol in/below window" in [g for g, _ in gr2] else np.nan)
            # transfer of the AOD seasonal cycle with its synoptic coefficient
            A2 = np.column_stack([np.ones(n2), X2])
            c2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
            b2 = dict(zip(use2, c2[1:]))
            pa, _, _, _ = harmonic(sub["t"], sub["aod"])
            out["tr_proj_aod"] = 100 * float(np.real(b2["aod"] * pa * np.conj(p_c)) / abs(p_c))
            out["aod_amp"] = abs(pa)
            out["aod_peak"] = peak_doy(pa)
            out["dlnC_per_0p1aod"] = 100 * b2["aod"] * 0.1
    return out


# ---------------------------------------------------------------------------------------------
# (b2) pooled monthly means with unit fixed effects
# ---------------------------------------------------------------------------------------------
def monthly_pooled(units, k_m2, exclude):
    Y, Yk, Xs, CAL, ids = [], [], [], [], []
    for meta, r in units:
        if meta["stream"] in exclude:
            continue
        dts = [EPOCH + timedelta(days=float(x)) for x in r["t"]]
        mon = np.array([d.year * 12 + d.month - 1 for d in dts])
        km = k_m2.get((meta["stream"], meta["module"]), np.nan)
        lk = kappa_norm(r, km) if np.isfinite(km) else np.full(len(mon), np.nan)
        rows_y, rows_yk, rows_x, rows_cal = [], [], [], []
        for mm in np.unique(mon):
            s = mon == mm
            if s.sum() < 4:
                continue
            rows_y.append(np.nanmean(r["lnC"][s]))
            rows_yk.append(np.nanmean(lk[s]))
            rows_x.append([np.nanmean(r[v][s]) for v in COVS])
            dmid = 2 * np.pi * ((mm % 12) + 0.5) / 12
            rows_cal.append([np.cos(dmid), np.sin(dmid)])
        if len(rows_y) < 8:
            continue
        y = np.array(rows_y)
        yk = np.array(rows_yk)
        x = np.array(rows_x)
        Y.append(y - y.mean())
        Yk.append(yk - np.nanmean(yk))
        Xs.append(x - np.nanmean(x, axis=0))
        CAL.append(np.array(rows_cal))
        ids += [meta["stream"]] * len(y)
    Y, Yk, X, CAL = np.concatenate(Y), np.concatenate(Yk), np.concatenate(Xs), np.concatenate(CAL)
    out = {}
    for lab, yy in (("op", Y), ("k_m2", Yk)):
        ok = np.isfinite(yy) & np.all(np.isfinite(X), axis=1)
        Z = X[ok] / X[ok].std(0)
        blocks, names = [], []
        for g, vs in GROUPS.items():
            ix = [COVS.index(v) for v in vs if np.std(Z[:, COVS.index(v)]) > 0]
            blocks.append(Z[:, ix])
            names.append(g)
        blocks.append(CAL[ok])
        names.append("calendar (season, unexplained)")
        dom, full, uniq, alone = dominance(blocks, yy[ok])
        out[lab] = dict(n=int(ok.sum()), sd=100 * np.std(yy[ok]), full=full,
                        groups={n_: (d, u, a) for n_, d, u, a in zip(names, dom, uniq, alone)},
                        r2_cal_alone=alone[-1])
    return out


# ---------------------------------------------------------------------------------------------
# Task 4: is the residual cycle geography- or generation-structured?
# ---------------------------------------------------------------------------------------------
def spearman(a, b):
    from scipy.stats import spearmanr
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    r, p = spearmanr(a[ok], b[ok])
    return float(r), float(p), int(ok.sum())


def task4(reg_rows):
    gen = {}
    for r in csv.DictReader(open(OUT / "q14_generation.csv", encoding="utf-8")):
        gen[r["stream"]] = (float(r["year"]) if r["year"] not in ("", "nan") else np.nan)
    R = [r for r in reg_rows if r["stream"] != PAYERNE]
    year = np.array([gen.get(r["stream"], np.nan) for r in R])
    cols = {"latitude": "lat", "longitude": "lon", "altitude": "alt", "continentality (night T_ext amp)": "T_ext_amp",
            "mean night T_ext": "T_ext_mean", "T_int annual amp": "T_int_amp", "window bottom (median, km AGL)": "zb_med",
            "kappa (multi-day)": "kappa_m2"}
    targets = {"amp_op (in phase)": "amp_inphase", "residual after kappa (in phase)": "inphase_k",
               "residual after all covariates (in phase)": "inphase_resid_tr"}
    out = []
    for tl, tk in targets.items():
        y = np.array([r.get(tk, np.nan) for r in R], float)
        for cl, ck in cols.items():
            x = np.array([r.get(ck, np.nan) for r in R], float)
            rho, p, n = spearman(x, y)
            out.append(dict(target=tl, predictor=cl, rho=rho, p=p, n=n))
        rho, p, n = spearman(year, y)
        out.append(dict(target=tl, predictor="module production year", rho=rho, p=p, n=n))
    # multiple regression of the kappa-residual in-phase amplitude: geography/climate block vs generation block
    y = np.array([r.get("inphase_k", np.nan) for r in R], float)
    geo = np.column_stack([np.array([r.get(k, np.nan) for r in R], float) for k in ("lat", "alt", "T_ext_amp", "T_ext_mean")])
    genb = np.column_stack([year, (year >= 2016).astype(float), (year >= 2020).astype(float)])
    ok = np.isfinite(y) & np.all(np.isfinite(geo), 1) & np.all(np.isfinite(genb), 1)
    Zg = (geo[ok] - geo[ok].mean(0)) / geo[ok].std(0)
    Zn = (genb[ok] - genb[ok].mean(0)) / np.where(genb[ok].std(0) > 0, genb[ok].std(0), 1)
    dom, full, uniq, alone = dominance([Zg, Zn], y[ok])
    from scipy.stats import f as fdist
    n_ = int(ok.sum())
    F = (full / 7) / ((1 - full) / (n_ - 8))
    blocks = dict(p_full=float(fdist.sf(F, 7, n_ - 8)), n=int(ok.sum()), r2_full=full, dom_geo=dom[0], dom_gen=dom[1], uniq_geo=uniq[0], uniq_gen=uniq[1],
                  alone_geo=alone[0], alone_gen=alone[1])
    return out, blocks, year


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------
def figures(reg_rows, pair_rows, site_rows, unit_rows, year, monthly):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = [r for r in reg_rows if r["stream"] != PAYERNE]
    gnames = list(GROUPS.keys())

    # Figure 1: budgets
    fig, axs = plt.subplots(1, 3, figsize=(21, 6.5))
    ax = axs[0]
    data = [[r.get(f"dom|{g}", 0) * 100 for r in R] for g in gnames]
    data.append([(1 - r["r2adj_syn"]) * 100 for r in R])
    ax.boxplot(data, vert=False, showfliers=False)
    ax.set_yticks(range(1, len(gnames) + 2))
    ax.set_yticklabels(gnames + ["unexplained"], fontsize=8)
    ax.set_xlabel("share of synoptic variance of ln C_L [%] (Shapley R², adjusted)")
    ax.set_title(f"(a) night-to-night (31-day anomalies), {len(R)} units", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[1]
    data = [[r.get(f"tr_proj|{g}", np.nan) for r in R] for g in gnames]
    data.append([r["amp"] for r in R])
    data.append([r["amp_resid_tr"] for r in R])
    ax.boxplot([np.array(d)[np.isfinite(d)] for d in data], vert=False, showfliers=False)
    ax.set_yticks(range(1, len(data) + 1))
    ax.set_yticklabels(gnames + ["observed |C_L| amplitude", "residual amplitude"], fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("annual amplitude of ln C_L, in phase with C_L [%]")
    ax.set_title("(b1) annual cycle: covariate harmonic × synoptic sensitivity", fontsize=10)
    ax.grid(alpha=0.3)
    ax = axs[2]
    labs = list(monthly["op"]["groups"].keys())
    x = np.arange(len(labs))
    for k, (lab, col) in enumerate((("op", "#1f77b4"), ("k_m2", "#ff7f0e"))):
        dv = [monthly[lab]["groups"][g][0] * 100 for g in labs]
        uv = [monthly[lab]["groups"][g][1] * 100 for g in labs]
        ax.barh(x + 0.38 * k - 0.19, dv, height=0.36, color=col,
                label=f"{'operational' if lab == 'op' else 'kappa-normalised'} (R² {monthly[lab]['full'] * 100:.0f} %)")
        ax.plot(uv, x + 0.38 * k - 0.19, "k|", ms=10)
    ax.set_yticks(x)
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xlabel("share of monthly-mean variance [%] (bar: Shapley; tick: unique)")
    ax.set_title("(b2) monthly means, unit fixed effects, pooled", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15f_budget.png", dpi=130)
    plt.close(fig)

    # Figure 2: pairs
    fig, axs = plt.subplots(1, 3, figsize=(21, 6.5))
    ax = axs[0]
    for v, mk in (("op", "o"), ("k_m2", "s")):
        rr = [r for r in unit_rows if r["variant"] == v]
        amp = np.array([r["amp"] for r in rr])
        pk = np.array([r["peak"] for r in rr])
        ax.scatter(pk, amp, marker=mk, s=50, label="operational" if v == "op" else "kappa (multi-day) normalised")
        for r, a_, p_ in zip(rr, amp, pk):
            ax.annotate(r["stream"].split("-")[-1], (p_, a_), fontsize=7)
    ax.set_xlabel("peak day of year")
    ax.set_ylabel("annual amplitude of ln C_L [%]")
    ax.set_title("(a) co-located units: annual harmonic", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axs[1]
    pr = [r for r in pair_rows if r["variant"] == "op"]
    lab = [r["pair"].replace("0-20000-0-", "").replace("0-20008-0-", "") for r in pr]
    yy = np.arange(len(pr))
    ax.barh(yy - 0.2, [r["common_amp"] for r in pr], height=0.4, label="common (mean phasor) amplitude")
    ax.barh(yy + 0.2, [r["amp_ratio"] for r in pr], height=0.4, label="ratio C_A/C_B amplitude (instrument-specific)")
    prk = {r["pair"]: r for r in pair_rows if r["variant"] == "k_m2"}
    ax.plot([prk[r["pair"]]["amp_ratio"] if r["pair"] in prk else np.nan for r in pr], yy + 0.2, "kD", ms=5,
            label="ratio amplitude after kappa")
    ax.plot([prk[r["pair"]]["common_amp"] if r["pair"] in prk else np.nan for r in pr], yy - 0.2, "k|", ms=10,
            label="common amplitude after kappa")
    ax.set_yticks(yy)
    ax.set_yticklabels(lab, fontsize=8)
    ax.set_xlabel("annual amplitude [%]")
    ax.set_title("(b) common-mode vs differential annual amplitude", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax = axs[2]
    ax.barh(yy - 0.2, [r["r_anom"] * 100 for r in pr], height=0.4, label="operational")
    ax.barh(yy + 0.2, [prk[r["pair"]]["r_anom"] * 100 if r["pair"] in prk else np.nan for r in pr], height=0.4,
            label="kappa-normalised")
    ax.set_yticks(yy)
    ax.set_yticklabels(lab, fontsize=8)
    ax.set_xlabel("correlation of synoptic anomalies A vs B [%] = common share of variance")
    ax.set_title("(c) night-to-night: common-mode share", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15f_pairs.png", dpi=130)
    plt.close(fig)

    # Figure 3: task 4
    fig, axs = plt.subplots(1, 4, figsize=(24, 6))
    ink = np.array([r.get("inphase_k", np.nan) for r in R])
    lat = np.array([r["lat"] for r in R])
    lon = np.array([r["lon"] for r in R])
    cont = np.array([r.get("T_ext_amp", np.nan) for r in R])
    ax = axs[0]
    sc = ax.scatter(lon, lat, c=ink, cmap="RdBu_r", vmin=-20, vmax=20, s=40, edgecolor="k", lw=0.3)
    plt.colorbar(sc, ax=ax, label="residual in-phase amplitude after kappa [%]")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_xlim(-12, 32)
    ax.set_ylim(35, 71)
    ax.set_title("(a) residual annual cycle after kappa, map", fontsize=10)
    for a, xv, xl in ((axs[1], lat, "latitude [°N]"), (axs[2], cont, "night T_ext annual amplitude [K] (continentality)"),
                      (axs[3], year, "module production year")):
        a.scatter(xv, ink, s=25)
        rho, p, n = spearman(xv, ink)
        a.set_xlabel(xl)
        a.set_ylabel("residual in-phase amplitude after kappa [%]")
        a.set_title(f"Spearman ρ = {rho:+.2f} (p {p:.1g}, n {n})", fontsize=10)
        a.axhline(0, color="k", lw=0.8)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(Q15 / "q15f_residual_structure.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
def write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def med_iqr(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return "n/a"
    return f"{np.median(x):+.1f} [{np.percentile(x, 25):+.1f}, {np.percentile(x, 75):+.1f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    Q15.mkdir(parents=True, exist_ok=True)
    units = build_table(a.rebuild)
    k_in, k_m2 = load_kappas()

    # Task 1
    unit_rows, pair_rows, site_rows = task1_pairs(units, k_in, k_m2)
    write_csv(unit_rows, Q15 / "q15f_pairs_units.csv")
    write_csv(pair_rows, Q15 / "q15f_pairs.csv")
    write_csv(site_rows, Q15 / "q15f_pairs_sites.csv")
    cl61_rows = task1_cl61({m["stream"]: (m, r) for m, r in units}, k_m2)
    write_csv(cl61_rows, Q15 / "q15f_cl61.csv")
    print("\n=== Task 1: co-located pairs ===")
    for r in pair_rows:
        print(f"  {r['pair']:28s} {r['variant']:5s} n {r['n_common']:3d}  amp A {r['amp_A']:5.1f} B {r['amp_B']:5.1f}  "
              f"common {r['common_amp']:5.1f}  ratio {r['amp_ratio']:5.1f} (peak {r['peak_ratio']:5.0f})  "
              f"scat A {r['scat_A']:4.1f} B {r['scat_B']:4.1f} ratio {r['scat_ratio']:4.1f}  r_res {r['r_resid']:+.2f}  "
              f"r_anom {r['r_anom']:+.2f}  u {r['rel_u_A']:.1f}/{r['rel_u_B']:.1f}")
    for r in pair_rows:
        if "common_sd" in r:
            print(f"  {r['pair']:28s} synoptic var A {r['var_A']:.0f} B {r['var_B']:.0f} common(cov) {r['var_common']:.0f} %^2 | "
                  f"common sd {r['common_sd']:.1f} R2adj {r['common_r2adj']:.2f} | diff sd {r['diff_sd']:.1f} R2adj {r['diff_r2adj']:.2f} | "
                  + " ".join(f"{g[:8]} c{100 * r.get('common_dom|' + g, 0):.0f}/d{100 * r.get('diff_dom|' + g, 0):.0f}" for g in GROUPS))
    for r in site_rows:
        print(f"  site {r['site_wmo']:18s} {r['variant']:5s} units {r['n_units']} common amp {r['common_amp']:5.1f} "
              f"(peak {r['common_peak']:4.0f})  unit-specific rms {r['specific_rms']:5.1f}  mean |amp| {r['mean_unit_amp']:5.1f}")
    for r in cl61_rows:
        print("  CL61", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})

    # Task 2/3
    reg_rows = []
    for meta, rows in units:
        km = k_m2.get((meta["stream"], meta["module"]), np.nan)
        r = unit_regression(meta, rows, km)
        if r is not None:
            reg_rows.append(r)
    write_csv(reg_rows, Q15 / "q15f_regression_units.csv")
    R = [r for r in reg_rows if r["stream"] != PAYERNE]
    print(f"\n=== Task 2/3: {len(R)} units (Payerne excluded) ===")
    print(f"  ln-harmonic amplitude {med_iqr([r['amp'] for r in R])} %, peak {np.median([r['peak'] for r in R]):.0f}; "
          f"scatter around harmonic {med_iqr([r['scat_tot'] for r in R])} %")
    print(f"  slow (31-d median - harmonic) {med_iqr([r['scat_slow'] for r in R])}  synoptic {med_iqr([r['scat_syn'] for r in R])}"
          f"  retrieval u/c {med_iqr([r['rel_u'] for r in R])}")
    print(f"  synoptic R2 adj {med_iqr([100 * r['r2adj_syn'] for r in R])} %   kappa_syn {med_iqr([r['kappa_syn'] for r in R])} vs "
          f"kappa_m2 {med_iqr([r['kappa_m2'] for r in R])}")
    for g in GROUPS:
        print(f"   {g:34s} dom {med_iqr([100 * r.get(f'dom|{g}', 0) for r in R])} %  unique "
              f"{med_iqr([100 * r.get(f'uniq|{g}', 0) for r in R])}  | annual transfer proj "
              f"{med_iqr([r.get(f'tr_proj|{g}', np.nan) for r in R])} %  amp {med_iqr([r.get(f'tr_amp|{g}', np.nan) for r in R])}"
              f"  peak {np.nanmedian([r.get(f'tr_peak|{g}', np.nan) for r in R]):.0f}")
    print(f"   transfer total proj {med_iqr([r['tr_proj_total'] for r in R])}  residual amp {med_iqr([r['amp_resid_tr'] for r in R])}"
          f" (peak {np.nanmedian([r['peak_resid_tr'] for r in R]):.0f})  after kappa {med_iqr([r.get('amp_k', np.nan) for r in R])}")
    print(f"   trend {med_iqr([r['trend_pct_yr'] for r in R])} %/yr; amp with trend {med_iqr([r['amp_with_trend'] for r in R])}")
    print("   per covariate: median standardised synoptic coefficient | annual transfer projection [%] | covariate in-phase amp")
    for v in COVS:
        print(f"     {v:10s} {med_iqr([r.get(f'bstd|{v}', np.nan) for r in R])} | {med_iqr([r.get(f'trv|{v}', np.nan) for r in R])}"
              f" | {med_iqr([r.get(f'harm_inphase|{v}', np.nan) for r in R])}")
    ip = np.array([r["amp_inphase"] for r in R])
    ipk = np.array([r.get("inphase_k", np.nan) for r in R])
    print(f"   in-phase amp op {med_iqr(ip)} -> after kappa {med_iqr(ipk)}; removed {med_iqr(ip - ipk)}; "
          f"peak after kappa {np.nanmedian([r.get('peak_k', np.nan) for r in R]):.0f}; units with residual in-phase > 0: "
          f"{np.mean(ipk[np.isfinite(ipk)] > 0) * 100:.0f} %")
    A = [r for r in R if "dom_aod" in r]
    print(f"   AERONET subset {len(A)} units: dom AOD {med_iqr([100 * r['dom_aod'] for r in A])} %  unique "
          f"{med_iqr([100 * r['uniq_aod'] for r in A])}  aer-contrast dom with AOD {med_iqr([100 * r['dom_aer_with_aod'] for r in A])}"
          f"  annual proj AOD {med_iqr([r['tr_proj_aod'] for r in A])} %  dlnC/0.1AOD {med_iqr([r['dlnC_per_0p1aod'] for r in A])}"
          f"  AOD peak {np.nanmedian([r['aod_peak'] for r in A]):.0f}")
    pay = [r for r in reg_rows if r["stream"] == PAYERNE]
    if pay:
        p = pay[0]
        print(f"  Payerne A: amp {p['amp']:.1f} peak {p['peak']:.0f} scat {p['scat_tot']:.1f} R2adj {p['r2adj_syn']:.2f}")

    monthly = monthly_pooled(units, k_m2, exclude={PAYERNE})
    mrows = []
    for lab, d in monthly.items():
        print(f"\n  monthly pooled [{lab}] n {d['n']} sd {d['sd']:.1f} % R2 full {d['full']:.2f}  calendar alone {d['r2_cal_alone']:.2f}")
        for g, (dv, uv, av) in d["groups"].items():
            print(f"     {g:34s} dom {100 * dv:5.1f}  unique {100 * uv:5.1f}  alone {100 * av:5.1f}")
            mrows.append(dict(variant=lab, group=g, n=d["n"], sd_pct=d["sd"], r2_full=d["full"], dominance=dv,
                              unique=uv, alone=av))
    write_csv(mrows, Q15 / "q15f_monthly_pooled.csv")

    # Task 4
    t4, blocks, year = task4(reg_rows)
    write_csv(t4, Q15 / "q15f_residual_structure.csv")
    print("\n=== Task 4 ===")
    for r in t4:
        print(f"  {r['target']:42s} vs {r['predictor']:36s} rho {r['rho']:+.2f} p {r['p']:.2g} n {r['n']}")
    print("  blocks:", {k: round(v, 3) if isinstance(v, float) else v for k, v in blocks.items()})
    json.dump(blocks, open(Q15 / "q15f_residual_blocks.json", "w"), indent=1)

    figures(reg_rows, pair_rows, site_rows, unit_rows, year, monthly)
    print(f"-> figures in {Q15}")


if __name__ == "__main__":
    main()
