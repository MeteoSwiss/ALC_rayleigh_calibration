# -*- coding: utf-8 -*-
"""q20 stage 2 - apply every quantified correction to the nightly CHM15k Rayleigh constant, one at a time and
cumulatively, and measure what each removes (report 19, q20).

Input:  q20_nights.csv (q20_build.py), q2a_intranight.csv (kappa), the q19 record run (tier per night).
Sample: the published nights 2025-01-01 .. 2026-07-31 that have housekeeping (needed for steps 1, 3, 5).
Metrics per series version (Kalman global-IQR screen re-applied to every version, as the operational
best-estimate filter would):
  * annual harmonic of ln C (OLS on day of year): amplitude [%], peak DOY, in-phase component along the
    operational phasor of the same unit [%];
  * night-to-night scatter: 1.4826 MAD of ln C about its 31-day running median [%];
  * scatter about the annual model (RMS) [%];
  * paper Kalman (kalman_improve_alc): annual-model amplitude and peak, and the daily best estimate.

Run:     .venv/Scripts/python.exe validation/chm15k_stability/q20_corrections.py
Outputs: q20_series.csv (every version per night), q20_unit_steps.csv (metrics per unit and version),
         q20_terms.csv (amplitude / scatter content of each correction term), q20_network.csv (medians),
         q20_other_test.csv (step-8 regression), q20_twins.csv (co-located ratios), q20_kalman.npz
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q20_lib import (AER_AERONET_LN, AER_CAMS_LN, AER_SPECIES, ALL_UNITS, K_PCAL, K_WIN, PAYERNE,  # noqa: E402
                     PCAL_MIN_TREND, Q20, STEPS, TIER2_BIAS, TWIN_PAIRS, UNITS, W, centred, doy, harmonic,
                     kalman_metrics, n2n_scatter, running_median, screen, series_metrics)
from common import OUT  # noqa: E402

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(2020)
NBOOT = 1000
Q19_RECORD = OUT / "q19" / "record"
OTHER_VARS = {"zb_km": "window bottom [km AGL]", "zt_km": "window top [km AGL]", "dur_h": "window duration [h]",
              "ln_nclear": "ln clear-profile count"}


# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------
def load_tiers(keys):
    """Tier of every successful 'op' night of the q19 record run: {(key, date): (tier, c)}."""
    out = {}
    for key in keys:
        for f in sorted((Q19_RECORD / key).glob("*.csv")):
            try:
                r = pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str})
            except Exception:  # noqa: BLE001 - a file being written by the record run
                continue
            r = r[(r.variant == "op") & r.flag.isin([1.0, 0.5])]
            for dt, tier, c in zip(r.date, r.tier, r.c):
                out[(key, dt)] = (tier, c)
    return out


def load_kappa():
    k = pd.read_csv(OUT / "q2a_intranight.csv").set_index("stream")
    return {s: (float(k.loc[s, "b_T_ctrl"]), float(k.loc[s, "se_T_ctrl"])) for s in k.index}


def prepare(d, tiers, kap):
    d = d.copy()
    d["tmid"] = 0.5 * (d.t0 + d.t1)
    d["lnc"] = np.log(d.c)
    d["zb_km"] = (d.zb_asl - d.alt) / 1000.0
    d["zt_km"] = (d.zt_asl - d.alt) / 1000.0
    d["dur_h"] = (d.t1 - d.t0) * 24.0
    d["ln_nclear"] = np.log(d.n_clear.where(d.n_clear > 0))
    tr = [tiers.get((s, dt), (np.nan, np.nan)) for s, dt in zip(d.stream, d.date)]
    d["tier"] = [x[0] for x in tr]
    d["c_q19"] = [x[1] for x in tr]
    # trust the label only where the record run reproduces the published constant exactly
    exact = np.abs(d.c_q19 / d.c - 1) < 1e-6
    d.loc[~exact, "tier"] = np.nan
    d["kappa"] = [kap.get(s, (np.nan, np.nan))[0] for s in d.stream]
    d["kappa_se"] = [kap.get(s, (np.nan, np.nan))[1] for s in d.stream]
    return d


# ---------------------------------------------------------------------------------------------
# Correction terms (ln units, centred on the unit's median; ln C_corrected = ln C - term)
# ---------------------------------------------------------------------------------------------
def term_kappa(g, full=False):
    T = g.T_int.values
    k = g.kappa.values[0] / 100.0
    if full:
        return k * centred(T)
    order = np.argsort(g.tmid.values)
    t_s, T_s = g.tmid.values[order], T[order]
    rm = np.empty(len(T))
    rm[order] = running_median(t_s, T_s, 31.0)
    dev = T - rm
    return np.where(np.isfinite(dev), k * dev, 0.0)


def term_aerosol(g, how="species"):
    if how == "species":
        out = np.zeros(len(g))
        for col, k in AER_SPECIES.items():
            out += k / 100.0 * centred(g[col].values)
        return out
    if how == "cams_ln":
        return AER_CAMS_LN / 100.0 * centred(np.log(g.aod550.clip(lower=1e-4).values))
    if how == "aeronet_ln":
        x = np.log(g.aod1064_aeronet.clip(lower=1e-4).values)
        return AER_AERONET_LN / 100.0 * centred(x) if np.isfinite(x).sum() >= 30 else np.zeros(len(g))
    raise ValueError(how)


def term_window(g):
    return K_WIN * centred(np.log(g.wt.where(g.wt > 0).values))


def term_module(g, lnc):
    """Era offsets (optical_module_id) fitted jointly with the annual harmonic on the IQR-screened series;
    the main era (most nights) is the reference. Returns (term, {module: offset %})."""
    mods = g.module.values
    if len(set(mods)) == 1:
        return np.zeros(len(g)), {}
    main = max(set(mods), key=list(mods).count)
    others = [m for m in sorted(set(mods)) if m != main]
    keep = screen(np.exp(lnc))
    d = doy(g.tmid.values)
    X = np.column_stack([np.ones(len(g)), np.cos(W * d), np.sin(W * d)] + [(mods == m).astype(float) for m in others])
    c, *_ = np.linalg.lstsq(X[keep], lnc[keep], rcond=None)
    off = {m: float(c[3 + i]) for i, m in enumerate(others)}
    term = np.array([off.get(m, 0.0) for m in mods])
    return term, {m: 100 * (np.exp(v) - 1) for m, v in off.items()}      # multiplicative step vs the main era


def term_pcal(g):
    """Imported p_calc ageing: per module era, linear trend of ln p_calc (with the annual harmonic) and,
    where |trend| >= PCAL_MIN_TREND %/yr, term = -K_PCAL x trend part (C moves opposite to p_calc)."""
    term = np.zeros(len(g))
    info = {}
    lnp = np.log(g.pcal.where(g.pcal > 0).values)
    for m in sorted(set(g.module.values)):
        e = (g.module.values == m) & np.isfinite(lnp)
        if e.sum() < 40:
            continue
        t = g.tmid.values[e]
        span = t.max() - t.min()
        if span < 180:
            continue
        tt = (t - t.mean()) / 365.25
        if span >= 300:
            d = doy(t)
            X = np.column_stack([np.ones(e.sum()), np.cos(W * d), np.sin(W * d), tt])
        else:
            X = np.column_stack([np.ones(e.sum()), tt])
        # robust: two passes, dropping 4-MAD outliers of ln p_calc
        c, *_ = np.linalg.lstsq(X, lnp[e], rcond=None)
        r = lnp[e] - X @ c
        good = np.abs(r - np.median(r)) <= 4 * 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-12
        c, *_ = np.linalg.lstsq(X[good], lnp[e][good], rcond=None)
        trend = 100 * c[-1]
        info[m] = trend
        if abs(trend) >= PCAL_MIN_TREND:
            term[e] = -K_PCAL * c[-1] * tt
    return centred_nonzero(term), info


def centred_nonzero(x):
    return x - np.median(x) if np.any(x != 0) else x


def term_molecular(g):
    return centred(g.b_mol.values)


# ---------------------------------------------------------------------------------------------
# Step 8: pooled synoptic regression of the step-7 series on window geometry / sampling
# ---------------------------------------------------------------------------------------------
def anomaly_by_unit(df, col):
    out = np.full(len(df), np.nan)
    for s, ix in df.groupby("stream").indices.items():
        g = df.iloc[ix]
        order = np.argsort(g.tmid.values)
        x = g[col].values.astype(float)
        rm = np.empty(len(ix))
        rm[order] = running_median(g.tmid.values[order], x[order], 31.0)
        out[ix] = x - rm
    return out


def other_test(series):
    """Coefficients (% per unit) of the step-7 ln C anomalies on the anomalies of OTHER_VARS, unit-cluster
    bootstrap. Returns a DataFrame and the dict of kept coefficients (95 % CI excluding 0 and |z| >= 3)."""
    df = series[series.stream.isin([k for k, *_ in UNITS]) & series.kept_tier2].copy()
    df["y"] = anomaly_by_unit(df, "lnc_tier2")
    cols = list(OTHER_VARS)
    for c in cols:
        df["a_" + c] = anomaly_by_unit(df, c)
    ok = np.isfinite(df.y) & np.all(np.isfinite(df[["a_" + c for c in cols]]), axis=1)
    # 4-MAD clip of the response, as an outlier guard
    yy = df.y[ok]
    ok &= np.abs(df.y - yy.median()) <= 4 * 1.4826 * np.median(np.abs(yy - yy.median()))
    X = df.loc[ok, ["a_" + c for c in cols]].values
    y = 100 * df.loc[ok, "y"].values
    k, *_ = np.linalg.lstsq(X, y, rcond=None)
    gidx = pd.Series(np.arange(ok.sum())).groupby(df.stream.values[ok]).indices
    units = list(gidx)
    kb = []
    for _ in range(NBOOT):
        ii = np.concatenate([gidx[u] for u in RNG.choice(units, len(units))])
        kb.append(np.linalg.lstsq(X[ii], y[ii], rcond=None)[0])
    kb = np.array(kb)
    se = kb.std(axis=0)
    lo, hi = np.percentile(kb, 2.5, axis=0), np.percentile(kb, 97.5, axis=0)
    rows, kept = [], {}
    sd_x = X.std(axis=0)
    for i, c in enumerate(cols):
        z = k[i] / se[i] if se[i] > 0 else np.nan
        keep = (lo[i] > 0 or hi[i] < 0) and abs(z) >= 3
        rows.append(dict(variable=c, meaning=OTHER_VARS[c], coef_pct_per_unit=k[i], se=se[i], ci_lo=lo[i], ci_hi=hi[i],
                         z=z, effect_per_sd_pct=k[i] * sd_x[i], kept=keep, n=int(ok.sum()), units=len(units)))
        if keep:
            kept[c] = k[i]
    return pd.DataFrame(rows), kept


def term_other(g, kept):
    out = np.zeros(len(g))
    for c, k in kept.items():
        out += k / 100.0 * centred(g[c].values)
    return out


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------
def build_series(d, kept_other=None):
    """Every version of ln C per night. Columns lnc_<step> (cumulative), lnc_only_<step> (one at a time),
    term_<step> (the correction term), plus sensitivity variants."""
    parts = []
    module_info, pcal_info = {}, {}
    for key, *_ in ALL_UNITS:
        g = d[d.stream == key].sort_values("tmid").copy()
        g = g[np.isfinite(g.T_int)]                 # common sample: nights with housekeeping
        lnc = g.lnc.values
        terms = {"kappa": term_kappa(g), "aerosol": term_aerosol(g, "species"), "window": term_window(g)}
        terms["module"], module_info[key] = term_module(g, lnc)
        terms["pcal"], pcal_info[key] = term_pcal(g)
        terms["molecular"] = term_molecular(g)
        terms["other"] = term_other(g, kept_other) if kept_other else np.zeros(len(g))
        tier2 = (g.tier.values == 2)
        g["is_tier2"] = tier2
        g["tier_known"] = np.isfinite(g.tier.values)
        cur = lnc.copy()
        g["lnc_raw"] = lnc
        for st in STEPS[1:]:
            if st == "tier2":
                g["lnc_tier2"] = cur.copy()          # values unchanged; the nights are dropped below
                g["term_tier2"] = np.where(tier2, np.log(1 + TIER2_BIAS / 100.0), 0.0)
                g["lnc_only_tier2"] = lnc.copy()
                continue
            cur = cur - terms[st]
            g[f"lnc_{st}"] = cur.copy()
            g[f"term_{st}"] = terms[st]
            g[f"lnc_only_{st}"] = lnc - terms[st]
        g["kept_tier2"] = ~tier2
        # sensitivity variants
        g["term_kappa_full"] = term_kappa(g, full=True)
        g["lnc_only_kappa_full"] = lnc - g.term_kappa_full.values
        g["term_aer_cams_ln"] = term_aerosol(g, "cams_ln") if "aod550" in g else 0.0
        g["lnc_only_aer_cams_ln"] = lnc - g.term_aer_cams_ln.values
        g["term_aer_aeronet"] = term_aerosol(g, "aeronet_ln")
        g["lnc_only_aer_aeronet"] = lnc - g.term_aer_aeronet.values
        # final with the full kappa instead of the synoptic one (for the twin test / sensitivity)
        g["term_kappa_seasonal"] = g.term_kappa_full.values - terms["kappa"]
        g["lnc_final_kfull"] = g["lnc_other"].values - g.term_kappa_seasonal.values
        # final with tier-2 nights scaled by -9.4 % instead of dropped
        g["lnc_final_t2scale"] = g["lnc_other"].values - g.term_tier2.values
        parts.append(g)
    return pd.concat(parts, ignore_index=True), module_info, pcal_info


VERSIONS = (["raw"] + [f"cum_{s}" for s in STEPS[1:]] + [f"only_{s}" for s in STEPS[1:]]
            + ["only_kappa_full", "only_aer_cams_ln", "only_aer_aeronet", "final_kfull", "final_ktwin",
               "final_t2scale"])


def version_values(g, v):
    """(times, ln C) of one series version of a unit (tier-2 nights dropped from step 7 on)."""
    t = g.tmid.values
    drop = g.is_tier2.values
    if v == "raw":
        return t, g.lnc_raw.values
    if v.startswith("cum_"):
        st = v[4:]
        y = g[f"lnc_{st}"].values
        if STEPS.index(st) >= STEPS.index("tier2"):
            return t[~drop], y[~drop]
        return t, y
    if v == "only_tier2":
        return t[~drop], g.lnc_raw.values[~drop]
    if v.startswith("only_"):
        return t, g[f"lnc_{v}"].values
    if v == "final_kfull":
        return t[~drop], g.lnc_final_kfull.values[~drop]
    if v == "final_ktwin":
        return t[~drop], g.lnc_final_ktwin.values[~drop]
    if v == "final_t2scale":
        return t, g.lnc_final_t2scale.values
    raise ValueError(v)


def unit_metrics(series):
    rows, terms, kal = [], [], {}
    for key, label, country, setting in ALL_UNITS:
        g = series[series.stream == key]
        if len(g) < 30:
            continue
        p_raw = series_metrics(*version_values(g, "raw"))["phasor"]
        for v in VERSIONS:
            t, y = version_values(g, v)
            m = series_metrics(t, y, pref=p_raw)
            rows.append(dict(stream=key, label=label, version=v, n=m["n"], amp=m["amp"], peak=m["peak"],
                             inphase=m["inphase"], n2n=m["n2n"], scat_model=m["scat_model"]))
        # content of each correction term on the raw (screened) sample
        keep = screen(g.c.values)
        for st in STEPS[1:]:
            x = g[f"term_{st}"].values
            p, *_ = harmonic(g.tmid.values[keep], x[keep])
            terms.append(dict(stream=key, label=label, term=st, amp=100 * abs(p) if np.isfinite(p) else np.nan,
                              peak=np.mod(np.angle(p) / W, 365.25) if np.isfinite(p) and abs(p) > 0 else np.nan,
                              inphase=100 * np.real(p * np.conj(p_raw)) / abs(p_raw) if np.isfinite(p) else np.nan,
                              sd=100 * np.std(x[keep]), n2n=n2n_scatter(g.tmid.values[keep], x[keep]) if np.std(x) > 0 else 0.0,
                              n_nonzero=int(np.sum(x != 0)), n_tier2=int(g.is_tier2.sum()) if st == "tier2" else np.nan))
        for v in ("raw", "cum_other", "final_kfull", "final_ktwin"):
            t, y = version_values(g, v)
            k = kalman_metrics(t, np.exp(y))
            if k is not None:
                kal[f"{key}__{v}"] = k
                rows.append(dict(stream=key, label=label, version=f"kalman_{v}", n=len(y), amp=k["amp"], peak=k["peak"],
                                 inphase=np.nan, n2n=np.nan, scat_model=100 * np.sqrt(k["var_meas"]) / np.median(np.exp(y))))
    return pd.DataFrame(rows), pd.DataFrame(terms), kal


def twins(series):
    rows = []
    for a, b in TWIN_PAIRS:
        ga = series[series.stream == a].drop_duplicates("night").set_index("night")
        gb = series[series.stream == b].drop_duplicates("night").set_index("night")
        common = ga.index.intersection(gb.index)
        for v, col in (("raw", "lnc_raw"), ("cum_kappa", "lnc_kappa"), ("final", "lnc_other"),
                       ("final_ktwin", "lnc_final_ktwin"), ("final_kfull", "lnc_final_kfull")):
            ok = common
            if v.startswith("final"):
                ok = [n for n in common if not ga.loc[n, "is_tier2"] and not gb.loc[n, "is_tier2"]]
            r = ga.loc[ok, col].values - gb.loc[ok, col].values
            t = ga.loc[ok, "tmid"].values
            keep = screen(np.exp(r))
            p, *_ = harmonic(t[keep], r[keep])
            rows.append(dict(pair=f"{a[-1]}/{b[-1]}" if a[:-2] == b[:-2] else f"{a}/{b}", site=a.rsplit("_", 1)[0],
                             version=v, n=int(keep.sum()), amp=100 * abs(p), peak=np.mod(np.angle(p) / W, 365.25),
                             n2n=n2n_scatter(t[keep], r[keep]), mean=100 * np.median(r[keep])))
    return pd.DataFrame(rows)


def twoway_demean(df, cols, n_iter=30):
    """Remove night and unit fixed effects (alternating projections on an unbalanced panel)."""
    out = df[cols].astype(float).copy()
    for _ in range(n_iter):
        out = out - out.groupby(df.night.values).transform("mean")
        out = out - out.groupby(df.stream.values).transform("mean")
    return out


def kappa_fraction(series):
    """Fraction f of the SEASONAL kappa x T_int term expressed between co-located twins.

    Panel regression on the twin sites: ln C_u,n (all corrections, synoptic kappa) = g_n + d_u + f s_u,n + e,
    with s_u = kappa_u x (31-day running median of T_int - unit median) (= term_kappa_seasonal), night fixed
    effects g_n (the atmosphere and any gain common to the co-located units cancel) and unit fixed effects
    d_u. f = 0: kappa acts on synoptic deviations only; f = 1: kappa acts on the full T_int. CI: bootstrap
    over calendar months (blocks)."""
    sites = {"Schiphol A-D": [k for k, *_ in UNITS if k.startswith("0-20000-0-06240")],
             "Budapest A/B": [k for k, *_ in UNITS if k.startswith("0-20000-0-12843")]}
    sites["Schiphol + Budapest"] = sites["Schiphol A-D"] + sites["Budapest A/B"]
    rows = []
    for name, units in sites.items():
        df = series[series.stream.isin(units) & series.kept_tier2][
            ["stream", "night", "tmid", "lnc_other", "term_kappa_seasonal"]].copy()
        if name == "Schiphol + Budapest":
            df["night"] = df["night"].astype(str) + np.where(df.stream.str.startswith("0-20000-0-06240"), "S", "B")
        df = df[df.groupby("night").stream.transform("count") >= 2]
        dm = twoway_demean(df, ["lnc_other", "term_kappa_seasonal"])
        y, x = dm.lnc_other.values, dm.term_kappa_seasonal.values
        med = np.median(y)
        mad = 1.4826 * np.median(np.abs(y - med))
        ok = np.abs(y - med) <= 4 * mad
        f = float(x[ok] @ y[ok] / (x[ok] @ x[ok]))
        months = pd.to_datetime(df.tmid.values, unit="D").strftime("%Y-%m").values
        gidx = pd.Series(np.arange(len(df))).groupby(months).indices
        keys = list(gidx)
        fb = []
        for _ in range(NBOOT):
            ii = np.concatenate([gidx[m] for m in RNG.choice(keys, len(keys))])
            ii = ii[ok[ii]]
            if (x[ii] @ x[ii]) > 0:
                fb.append(x[ii] @ y[ii] / (x[ii] @ x[ii]))
        rows.append(dict(sites=name, units=len(units), nights=int(df.night.nunique()), obs=int(ok.sum()), f=f,
                         ci_lo=float(np.percentile(fb, 2.5)), ci_hi=float(np.percentile(fb, 97.5)),
                         sd_seasonal_term=100 * float(np.std(x[ok]))))
    return pd.DataFrame(rows)


def between_unit_fraction(series):
    """Cross-unit test: annual phasor of the fully corrected series (synoptic kappa) of each unit, P_u, against
    the annual phasor of its seasonal kappa term, S_u (kappa_u x seasonal T_int), over the 20 units:
        P_u = P_G + f S_u + e_u          (P_G complex = gain common to all units, f real)
    P_G is the part of the cycle that does not scale with the unit's kappa x T swing: its modulus is the
    common CHM15k gain amplitude, its argument the peak day. Unit bootstrap (95 % CI)."""
    P, S, keys = [], [], []
    for key, *_ in UNITS:
        g = series[(series.stream == key) & series.kept_tier2]
        keep = screen(np.exp(g.lnc_other.values))
        p, *_ = harmonic(g.tmid.values[keep], g.lnc_other.values[keep])
        s_, *_ = harmonic(g.tmid.values[keep], g.term_kappa_seasonal.values[keep])
        if np.isfinite(p) and np.isfinite(s_):
            P.append(p)
            S.append(s_)
            keys.append(key)
    P, S = np.array(P), np.array(S)

    def fit(ii):
        y = np.r_[P[ii].real, P[ii].imag]
        n = len(ii)
        X = np.zeros((2 * n, 3))
        X[:n, 0] = 1
        X[n:, 1] = 1
        X[:n, 2] = S[ii].real
        X[n:, 2] = S[ii].imag
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        return c[0] + 1j * c[1], c[2]

    idx = np.arange(len(P))
    pg, f = fit(idx)
    boot = [fit(RNG.choice(idx, len(idx))) for _ in range(NBOOT)]
    fb = np.array([b[1] for b in boot])
    ab = np.array([abs(b[0]) for b in boot])
    pk = np.array([np.mod(np.angle(b[0]) / W, 365.25) for b in boot])
    # fixed f = 0 and f = 1 for reference (mean phasor of P - f S)
    rows = [dict(model="free f", f=f, f_lo=np.percentile(fb, 2.5), f_hi=np.percentile(fb, 97.5),
                 common_amp=100 * abs(pg), common_amp_lo=100 * np.percentile(ab, 2.5),
                 common_amp_hi=100 * np.percentile(ab, 97.5), common_peak=np.mod(np.angle(pg) / W, 365.25),
                 common_peak_iqr=np.percentile(pk, 75) - np.percentile(pk, 25), units=len(P),
                 rms_misfit=100 * np.sqrt(np.mean(np.abs(P - pg - f * S) ** 2)))]
    for f0 in (0.0, 1.0):
        pg0 = np.mean(P - f0 * S)
        rows.append(dict(model=f"f fixed {f0:.0f}", f=f0, common_amp=100 * abs(pg0),
                         common_peak=np.mod(np.angle(pg0) / W, 365.25), units=len(P),
                         rms_misfit=100 * np.sqrt(np.mean(np.abs(P - pg0 - f0 * S) ** 2))))
    per_unit = pd.DataFrame(dict(stream=keys, P_amp=100 * np.abs(P), P_peak=np.mod(np.angle(P) / W, 365.25),
                                 S_amp=100 * np.abs(S), S_peak=np.mod(np.angle(S) / W, 365.25)))
    return pd.DataFrame(rows), per_unit


def mirror_test(series):
    """Solstice-mirror test as in q16b: mean of spring (DOY 60-141) minus mean of autumn (DOY 203-284), per unit,
    on the 2025 nights only (one complete year; avoids mixing two springs with one autumn). A solstice-phased
    driver gives 0; a driver following T_int (peak ~DOY 206) gives spring < autumn."""
    rows = []
    y2025 = (series.tmid >= T_START_2025) & (series.tmid < T_START_2025 + 365)
    d = doy(series.tmid.values)
    spring = (d >= 60) & (d <= 141)
    autumn = (d >= 203) & (d <= 284)
    for key, label, *_ in ALL_UNITS:
        g = (series.stream == key).values & y2025.values
        out = dict(stream=key, label=label)
        for name, col, drop in (("lnC_raw", "lnc_raw", False), ("final_kappa_synoptic", "lnc_other", True),
                                ("final_kappa_full", "lnc_final_kfull", True), ("T_int", "T_int", False)):
            m = g & (~series.is_tier2.values if drop else True)
            v = series[col].values
            sp, au = v[m & spring], v[m & autumn]
            scale = 1.0 if name == "T_int" else 100.0
            out[name] = scale * (np.nanmean(sp) - np.nanmean(au)) if (np.isfinite(sp).sum() >= 8 and
                                                                      np.isfinite(au).sum() >= 8) else np.nan
            out[f"n_{name}"] = f"{np.isfinite(sp).sum()}/{np.isfinite(au).sum()}"
        rows.append(out)
    return pd.DataFrame(rows)


T_START_2025 = float((pd.Timestamp("2025-01-01") - pd.Timestamp("1970-01-01")).days)


def tier2_offset(g, days=45.0, min_n=3):
    """Median offset [%] of the unit's tier-2 nights from its tier-1 nights within +-days (q16e definition)."""
    t1 = g[(g.tier == 1)]
    off = []
    for t, y in zip(g.tmid[g.is_tier2], g.lnc_raw[g.is_tier2]):
        m = np.abs(t1.tmid.values - t) <= days
        if m.sum() >= min_n:
            off.append(y - np.median(t1.lnc_raw.values[m]))
    return (100 * (np.exp(np.median(off)) - 1), len(off)) if off else (np.nan, 0)


def budget_table(steps, terms):
    """Network budget per correction step (20 units): content of the term, and amplitude / phase / scatter
    after the cumulative step, with per-unit changes (median [IQR]) cumulative and one-at-a-time."""
    keys = [k for k, *_ in UNITS]
    st = steps[steps.stream.isin(keys)]
    te = terms[terms.stream.isin(keys)]
    P = st.pivot(index="stream", columns="version", values="amp")
    N = st.pivot(index="stream", columns="version", values="n2n")
    K = st.pivot(index="stream", columns="version", values="peak")

    def q3(x):
        x = np.asarray(x, float)
        return np.nanmedian(x), np.nanpercentile(x, 25), np.nanpercentile(x, 75)

    rows = []
    prev = "raw"
    for stp in STEPS:
        v = "raw" if stp == "raw" else f"cum_{stp}"
        r = dict(step=stp, n_units_affected=int((te[te.term == stp].n_nonzero > 0).sum()) if stp != "raw" else 20)
        r["amp_med"], r["amp_p25"], r["amp_p75"] = q3(P[v])
        r["peak_circ"] = circ_mean(K[v].values)
        r["n2n_med"], r["n2n_p25"], r["n2n_p75"] = q3(N[v])
        if stp != "raw":
            tt = te[te.term == stp]
            r["term_inphase_med"], r["term_inphase_p25"], r["term_inphase_p75"] = q3(tt.inphase)
            r["term_n2n_med"] = np.nanmedian(tt.n2n)
            r["d_amp_cum_med"], r["d_amp_cum_p25"], r["d_amp_cum_p75"] = q3(P[v] - P[prev])
            r["d_n2n_cum_med"], r["d_n2n_cum_p25"], r["d_n2n_cum_p75"] = q3(N[v] - N[prev])
            r["d_amp_only_med"], _, _ = q3(P[f"only_{stp}"] - P["raw"])
            r["d_n2n_only_med"], _, _ = q3(N[f"only_{stp}"] - N["raw"])
        rows.append(r)
        prev = v
    for v in ("final_ktwin", "final_kfull", "final_t2scale", "only_kappa_full", "only_aer_cams_ln", "only_aer_aeronet"):
        r = dict(step=v)
        r["amp_med"], r["amp_p25"], r["amp_p75"] = q3(P[v])
        r["peak_circ"] = circ_mean(K[v].values)
        r["n2n_med"], r["n2n_p25"], r["n2n_p75"] = q3(N[v])
        r["d_amp_only_med"], _, _ = q3(P[v] - P["raw"])
        r["d_n2n_only_med"], _, _ = q3(N[v] - N["raw"])
        rows.append(r)
    return pd.DataFrame(rows)


def network(steps):
    sub = steps[steps.stream.isin([k for k, *_ in UNITS])]
    rows = []
    for v in VERSIONS + ["kalman_raw", "kalman_cum_other", "kalman_final_kfull", "kalman_final_ktwin"]:
        x = sub[sub.version == v]
        if x.empty:
            continue
        rows.append(dict(version=v, units=int(np.isfinite(x.amp).sum()), amp_med=np.nanmedian(x.amp),
                         amp_p25=np.nanpercentile(x.amp, 25), amp_p75=np.nanpercentile(x.amp, 75),
                         inphase_med=np.nanmedian(x.inphase), peak_circmed=circ_mean(x.peak.values),
                         n2n_med=np.nanmedian(x.n2n), n2n_p25=np.nanpercentile(x.n2n, 25) if x.n2n.notna().any() else np.nan,
                         n2n_p75=np.nanpercentile(x.n2n, 75) if x.n2n.notna().any() else np.nan,
                         scat_model_med=np.nanmedian(x.scat_model)))
    return pd.DataFrame(rows)


def circ_mean(peaks):
    p = np.asarray(peaks, float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan
    z = np.mean(np.exp(1j * W * p))
    return float(np.mod(np.angle(z) / W, 365.25))


def main():
    d = pd.read_csv(Q20 / "q20_nights.csv", dtype={"date": str, "module": str})
    keys = [k for k, *_ in ALL_UNITS]
    tiers = load_tiers(keys)
    d = prepare(d, tiers, load_kappa())
    # tier coverage and reproduction check
    ok = np.isfinite(d.c_q19)
    ratio = (d.c_q19 / d.c)[ok]
    print(f"tier labels: {ok.sum()}/{len(d)} nights matched; |C_q19/C_pub - 1| < 1e-6 on "
          f"{(np.abs(ratio - 1) < 1e-6).sum()} of them; tier-2 share {np.mean(d.tier[ok] == 2) * 100:.1f} %")
    # pass 1 without step 8, regression test, pass 2 with the kept terms
    series, _, _ = build_series(d)
    other, kept = other_test(series)
    other.to_csv(Q20 / "q20_other_test.csv", index=False)
    print(other.to_string(index=False))
    print("step 8 kept:", kept or "none")
    series, module_info, pcal_info = build_series(d, kept)
    kf = kappa_fraction(series)
    kf.to_csv(Q20 / "q20_kappa_fraction.csv", index=False)
    print(kf.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    bu, bu_units = between_unit_fraction(series)
    bu.to_csv(Q20 / "q20_between_units.csv", index=False)
    bu_units.to_csv(Q20 / "q20_between_units_phasors.csv", index=False)
    print(bu.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    f_twin = float(kf.loc[kf.sites == "Schiphol + Budapest", "f"].iloc[0])
    series["f_twin"] = f_twin
    series["lnc_final_ktwin"] = series.lnc_other - f_twin * series.term_kappa_seasonal
    series.to_csv(Q20 / "q20_series.csv", index=False)
    steps, terms, kal = unit_metrics(series)
    steps.to_csv(Q20 / "q20_unit_steps.csv", index=False)
    terms.to_csv(Q20 / "q20_terms.csv", index=False)
    np.savez_compressed(Q20 / "q20_kalman.npz", **{f"{k}__{f}": v[f] for k, v in kal.items()
                                                   for f in ("grid", "state", "std")})
    net = network(steps)
    net.to_csv(Q20 / "q20_network.csv", index=False)
    bt = budget_table(steps, terms)
    bt.to_csv(Q20 / "q20_budget.csv", index=False)
    print(bt.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    tw = twins(series)
    tw.to_csv(Q20 / "q20_twins.csv", index=False)
    mt = mirror_test(series)
    mt.to_csv(Q20 / "q20_mirror.csv", index=False)
    sub = mt[mt.stream.isin([k for k, *_ in UNITS])]
    for c in ("lnC_raw", "final_kappa_synoptic", "final_kappa_full", "T_int"):
        x = sub[c].dropna().values
        bs = [np.median(RNG.choice(x, len(x))) for _ in range(NBOOT)]
        print(f"mirror spring-autumn 2025 {c:22s}: median {np.median(x):+.2f} (95 % CI {np.percentile(bs, 2.5):+.2f} .. "
              f"{np.percentile(bs, 97.5):+.2f}), n = {len(x)}")
    # per-unit info table
    info = []
    for key, label, country, setting in ALL_UNITS:
        g = series[series.stream == key]
        if g.empty:
            continue
        mods_in_order = list(dict.fromkeys(g.sort_values("tmid").module))
        t2_off, t2_n = tier2_offset(g)
        info.append(dict(stream=key, label=label, country=country, setting=setting, lat=g.lat.iloc[0], lon=g.lon.iloc[0],
                         alt=g.alt.iloc[0], modules="/".join(mods_in_order),
                         module_years="/".join(str(2000 + int(m[:2])) for m in mods_in_order),
                         rel_u_med=100 * float(np.median(g.u / g.c)), n_nights=len(g),
                         n_tier2=int(g.is_tier2.sum()), tier_known=int(g.tier_known.sum()),
                         tier2_offset_pct=t2_off, tier2_offset_n=t2_n,
                         n_aeronet=int(np.isfinite(g.aod1064_aeronet).sum()),
                         kappa=g.kappa.iloc[0], kappa_se=g.kappa_se.iloc[0],
                         module_offsets="; ".join(f"{m}: {v:+.1f} %" for m, v in module_info[key].items()),
                         pcal_trend="; ".join(f"{m}: {v:+.1f} %/yr" for m, v in pcal_info[key].items()),
                         aeronet=g.aeronet_site.iloc[0] if isinstance(g.aeronet_site.iloc[0], str) else ""))
    pd.DataFrame(info).to_csv(Q20 / "q20_units.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(info)[["label", "modules", "n_nights", "n_tier2", "tier_known", "kappa", "module_offsets",
                              "pcal_trend"]].to_string(index=False))
    print(net.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(tw.to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()
