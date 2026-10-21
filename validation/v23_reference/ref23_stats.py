# -*- coding: utf-8 -*-
"""ref23 step 4: statistics of ln rho against the predicted uncertainty.

Per site x pair x calibration version (v2.2 published, v2.2 constant + v2.3 budget, v2.3):
  n             nights with a value
  bias          median ln rho [%], 95 % CI from a weekly block bootstrap (nights of one calendar
                week resampled together: consecutive nights share their aerosol and their gain state)
  scatter       robust std of ln rho = 1.4826 x MAD [%]
  sig_pred      median predicted 1 sigma of ln rho [%]
  pull_u        ln rho / sig_pred (tests bias + scatter): robust std about ZERO = 1.4826 x median|pull|,
                and the fraction |pull| < 2 (Gaussian 95.4 %)
  pull_c        (ln rho - median) / sig_pred (tests the scatter alone): same two numbers
  spearman      rank correlation of |ln rho - median| with sig_pred (does the uncertainty predict
                the disagreement?), with its p-value
Predicted sigma (1 sigma, relative):  sqrt(sig_X^2 + sig_R^2 + sig_conv^2) where
  * Rayleigh-vs-Rayleigh pairs drop the molecular-model term of both units (common: US-1976 vs the
    real T/p in similar windows; doc/v2.3/uncertainty.md drops it for twins) -- "with mol" is the
    sensitivity;
  * Rayleigh-vs-cloud pairs keep it (the cloud calibration has no molecular model);
  * sig_conv = the lambda-conversion term of the night, |ln rho(alpha=1.5) - ln rho(alpha=0.5)| / 2
    (Angstrom exponent alpha = 1 +- 0.5), only for 910 nm vs 1064 nm pairs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

VERSIONS = {            # version -> (ln rho column state, sigma state of the Rayleigh channels)
    "v22pub": ("v22pub", "v22pub"),
    "v22": ("v22", "v22"),
    "v22b": ("v22", "v22b"),
    "v23": ("v23", "v23"),             # eprof_v2.3 + the REVIEWED v2.3 budget (main)
    "v23pre": ("v23pre", "v23pre"),    # eprof_v2.3 + the pre-review budget (secondary)
    "cloud": ("cloud", "cloud"),       # cloud / cloud pairs (no Rayleigh member)
}
VERSION_LABEL = {"v22pub": "v2.2 published", "v22": "v2.2 rerun",
                 "v22b": "v2.2 const. + reviewed v2.3 budget", "v23": "v2.3 (reviewed budget)",
                 "v23pre": "v2.3 (pre-review budget)", "cloud": "cloud / cloud"}


def lnr_col(df, version, band="", alpha="", w="", ref=""):
    st = VERSIONS[version][0]
    return f"lnr_{st}{band}{alpha}{w}{ref}"


def _conv_cols(df, st_l, band, w):
    """The alpha = 1.5 / 0.5 columns of the same observable (falls back to the band's default
    5-min observable when the variant has no alpha sweep, e.g. the 60-min window)."""
    hi, lo = f"lnr_{st_l}{band}{w}_a1.5", f"lnr_{st_l}{band}{w}_a0.5"
    if hi not in df:
        hi, lo = f"lnr_{st_l}{band}_a1.5", f"lnr_{st_l}{band}_a0.5"
    nan = pd.Series(np.nan, index=df.index)
    return df.get(hi, nan).values, df.get(lo, nan).values


def sigma_pred(df, version, with_mol=False, with_conv=True, band="", cloud_model=None, ref="",
               w=""):
    """Predicted 1-sigma of ln rho per row for one version (see the module docstring).
    cloud_model: None = the primary cloud-reference sigma (cr_sig), else cr_sig_<cloud_model>;
    ref "_cl" = the local (+-7 d median) cloud reference and its sigma (cr_sig_cl)."""
    kind = df["kind"].values
    st_l, st_s = VERSIONS[version]
    n = len(df)
    sx = np.full(n, np.nan)
    sr = np.full(n, np.nan)
    rr = np.isin(kind, ("RR_twin", "RR_mixed"))
    col = (lambda pre, s_: f"{pre}_{s_}")
    sxa = df.get(col("sx", st_s), pd.Series(np.nan, index=df.index)).values
    sxn = df.get(col("sxnm", st_s), pd.Series(np.nan, index=df.index)).values
    sra = df.get(col("sr", st_s), pd.Series(np.nan, index=df.index)).values
    srn = df.get(col("srnm", st_s), pd.Series(np.nan, index=df.index)).values
    sfx = "" if cloud_model is None else f"_{cloud_model}"
    if ref == "_cl":
        sfx = "_cl"
    crs = df.get(f"cr_sig{sfx}", pd.Series(np.nan, index=df.index)).values
    cxs = df.get(f"cx_sig{sfx}", pd.Series(np.nan, index=df.index)).values
    # Rayleigh X
    xr = np.isin(kind, ("RR_twin", "RR_mixed", "RC", "self"))
    sx[xr] = np.where(rr[xr] & (not with_mol), sxn[xr], sxa[xr])
    sx[kind == "CC"] = cxs[kind == "CC"]
    sr[rr] = srn[rr] if not with_mol else sra[rr]
    cref = np.isin(kind, ("RC", "self", "CC"))
    sr[cref] = crs[cref]
    tot = sx ** 2 + sr ** 2
    if with_conv:
        sc = sigma_conv(df, version, band, w)
        conv = np.isin(kind, ("RC", "RR_mixed"))
        sc = np.where(conv & np.isfinite(sc), sc, 0.0)
        tot = tot + sc ** 2
    return np.sqrt(tot)


def sigma_conv(df, version, band="", w=""):
    """Per-night lambda-conversion term: half the change of ln rho for alpha 0.5 -> 1.5."""
    st_l = VERSIONS[version][0]
    a_hi, a_lo = _conv_cols(df, st_l, band, w)
    return np.abs(a_hi - a_lo) / 2.0


def block_bootstrap_median(dates, x, nboot=4000, seed=7):
    """95 % CI of the median with calendar weeks as resampling blocks."""
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    x = x[ok]
    if x.size < 5:
        return np.nan, np.nan
    wk = pd.to_datetime(pd.Series(np.asarray(dates)[ok])).dt.to_period("W").astype(str).values
    groups = [x[wk == w] for w in np.unique(wk)]
    rng = np.random.default_rng(seed)
    meds = np.empty(nboot)
    for b in range(nboot):
        pick = rng.integers(0, len(groups), len(groups))
        meds[b] = np.median(np.concatenate([groups[i] for i in pick]))
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def robust_std(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return np.nan
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def robust_std0(x):
    """Robust std about ZERO (for the uncentred pull: bias and scatter together)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return np.nan
    return float(1.4826 * np.median(np.abs(x)))


def summarize(df, version, band="", w="", with_mol=False, with_conv=True, min_n=8,
              cloud_model=None, ref="", center_by=None):
    """One dict of statistics for the rows of df (one site x pair, or a pooled group).
    center_by: a column of df whose groups are centred separately (pooled groups: each pair on its
    own median) for the centred pull, the observed scatter and the Spearman test."""
    col = lnr_col(df, version, band, "", w, ref)
    if col not in df:
        return None
    y = df[col].values.astype(float)
    sp = sigma_pred(df, version, with_mol=with_mol, with_conv=with_conv, band=band,
                    cloud_model=cloud_model, ref=ref, w=w)
    ok = np.isfinite(y) & np.isfinite(sp) & (sp > 0)
    n = int(ok.sum())
    if n < min_n:
        return dict(n=n)
    has_conv = bool(np.isin(df["kind"].values, ("RC", "RR_mixed")).any())
    scv = sigma_conv(df, version, band, w)[ok]
    y, sp = y[ok], sp[ok]
    dates = df["date"].values[ok]
    med = float(np.median(y))
    lo, hi = block_bootstrap_median(dates, y)
    if center_by is not None:
        grp = df[center_by].values[ok]
        cen = pd.Series(y).groupby(grp).transform("median").values
    else:
        cen = np.full(y.size, med)
    pu = y / sp
    pc = (y - cen) / sp
    rho_s, p_s = spearmanr(np.abs(y - cen), sp)
    return dict(
        n=n, bias=100 * med, ci_lo=100 * lo, ci_hi=100 * hi,
        bias_sig=bool(np.isfinite(lo) and (lo > 0 or hi < 0)),
        scatter=100 * robust_std(y - cen), rms=100 * float(np.sqrt(np.mean((y - cen) ** 2))),
        sig_pred=100 * float(np.median(sp)),
        pull_u=robust_std0(pu), in2_u=100 * float(np.mean(np.abs(pu) < 2)),
        pull_c=robust_std0(pc), in2_c=100 * float(np.mean(np.abs(pc) < 2)),
        spearman=float(rho_s), spearman_p=float(p_s),
        sig_conv=(100 * float(np.nanmedian(scv)) if has_conv and np.any(np.isfinite(scv))
                  else 0.0),
    )


PULL_LO, PULL_HI = 0.75, 1.33     # centred-pull robust std accepted as "scatter predicted"
COVER_OK = 90.0                   # % of nights within 2 sigma (uncentred) accepted as "covered"


def verdict(s):
    """Short verdict from the statistics of one pair (rules stated in the report):
    consistent        no significant bias, >= 90 % of the nights within 2 sigma (uncentred) and a
                      centred pull robust std within [0.75, 1.33];
    significant bias  the 95 % CI of the median ln rho excludes 0;
    u under-/over-stated  centred pull robust std > 1.33 / < 0.75 (the night-to-night scatter is
                      larger / smaller than the predicted sigma)."""
    if s is None or s.get("n", 0) < 8:
        return "too few nights"
    sc_ok = PULL_LO <= s["pull_c"] <= PULL_HI
    if not s["bias_sig"] and s["in2_u"] >= COVER_OK and sc_ok:
        return "consistent"
    parts = []
    if s["bias_sig"]:
        cov = "within" if s["in2_u"] >= COVER_OK else "beyond"
        parts.append(f"significant bias ({s['bias']:+.0f} %, {cov} the stated u)")
    if s["pull_c"] > PULL_HI:
        parts.append("u under-stated")
    elif s["pull_c"] < PULL_LO:
        parts.append("u over-stated")
    elif not s["bias_sig"]:
        parts.append(f"no significant bias, scatter as predicted, {s['in2_u']:.0f} % within 2 sigma")
    else:
        parts.append("scatter as predicted")
    return "; ".join(parts)
