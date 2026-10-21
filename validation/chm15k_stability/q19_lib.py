# -*- coding: utf-8 -*-
"""q19 - shared analysis helpers (offsets vs tier-1 neighbours, offline gate evaluation, re-selection)."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
import q19_gradient_gate as Q  # noqa: E402

PAYERNE_KEY = "_".join(Q.PAYERNE)
CELL_COLS = list(Q.CELL_COLS)


def load(types=("CHM15k", "CL61", "CL31")):
    """Record run, one row per (unit, night, variant); adds t, ok, lnc, month."""
    df = Q.load_record()
    df = df[df.type.isin(types)].copy()
    df["t"] = pd.to_datetime(df.date)
    df["month"] = df.t.dt.month
    df["ok"] = df.flag.isin(Q.SUCCESS) & (df.c > 0)
    df["lnc"] = np.where(df.ok, np.log(df.c.where(df.c > 0)), np.nan)
    return df


def neighbour_offset(t_eval, lnc_eval, t_ref, lnc_ref, days=45, min_n=3):
    """ln C of each evaluated night minus the median ln C of the reference nights within +-days."""
    t_ref = np.asarray(t_ref, "datetime64[D]").astype(float)
    t_eval = np.asarray(t_eval, "datetime64[D]").astype(float)
    out = np.full(len(t_eval), np.nan)
    nref = np.zeros(len(t_eval), int)
    for k, (t, y) in enumerate(zip(t_eval, lnc_eval)):
        m = np.abs(t_ref - t) <= days
        nref[k] = int(m.sum())
        if m.sum() >= min_n and np.isfinite(y):
            out[k] = y - np.median(lnc_ref[m])
    return out, nref


def tier1_offsets(df_unit_variant, t_eval, lnc_eval, days=45):
    """Offset of evaluated nights vs the unit's TIER-1 nights (op variant) within +-days."""
    r = df_unit_variant[(df_unit_variant.tier == 1) & df_unit_variant.ok]
    return neighbour_offset(t_eval, lnc_eval, r.t.values, r.lnc.values, days)


# ------------------------------------------------------------------------------------------------
# gate rules evaluated OFFLINE on stored statistics (dict-like rows with the STAT_KEYS columns)
# ------------------------------------------------------------------------------------------------
def rule_verdict(row, rule, g_ref=0.0):
    """Vectorisable verdict of one rule on a pandas row/dict of statistics (prefix already removed).

    rule = dict(stat='up'|'ext', g_tol, k, s_max, min_above_m). Returns one of
    'accept', 'reject_aerosol', 'reject_rising', 'reject_unverifiable'."""
    st = rule.get("stat", "up")
    g = row.get(f"g_{st}", np.nan)
    sg = row.get(f"sg_{st}", np.nan)
    above = row.get("z_hi", np.nan) - row.get("ze", np.nan)
    if not (np.isfinite(g) and np.isfinite(sg)):
        return "reject_unverifiable"
    if rule.get("min_above_m", 0) > 0 and not (np.isfinite(above) and above >= rule["min_above_m"]):
        return "reject_unverifiable"
    if sg > rule.get("s_max", np.inf):
        return "reject_unverifiable"
    d = g - g_ref
    lim = rule["g_tol"] + rule["k"] * sg
    if d < -lim:
        return "reject_aerosol"
    if d > lim:
        return "reject_rising" if rule.get("reject_rising", True) else "accept_rising"
    return "accept"


def accepted(v):
    return np.char.startswith(np.asarray(v, str), "accept")


def verdicts(frame, prefix, rule, g_ref=None):
    """Apply rule_verdict row by row on columns '<prefix><stat>'. g_ref: scalar or per-row array."""
    cols = [c for c in frame.columns if c.startswith(prefix)]
    sub = frame[cols].rename(columns=lambda c: c[len(prefix):])
    gr = np.zeros(len(frame)) if g_ref is None else np.broadcast_to(np.asarray(g_ref, float), (len(frame),))
    return np.array([rule_verdict(r, rule, float(x)) for r, x in zip(sub.to_dict("records"), gr)])


def cells_of(key, date, variant="op"):
    """Stored relaxed (tier-2) cells of one night as a DataFrame (CELL_COLS), or None."""
    wmo, ident = key.rsplit("_", 1)
    f = Q.REC_DIR / key / f"{date[:6]}.npz"
    if not f.exists():
        return None
    with np.load(f) as z:
        k = f"c_{variant}_{date}"
        if k not in z.files:
            return None
        a = z[k]
    return pd.DataFrame(a.astype(float), columns=CELL_COLS)


def profile_of(key, date):
    """Stored profile (z, signal, p_mol, sigma_emp, sigma_op) of one night, or None."""
    f = Q.REC_DIR / key / f"{date[:6]}.npz"
    if not f.exists():
        return None
    with np.load(f) as z:
        k = f"p_{date}"
        if k not in z.files:
            return None
        return z[k]


def reselect(cells, rule, g_ref=0.0):
    """Offline gated re-selection: highest-quality relaxed cell that passes the rule.

    Returns (index into cells or -1, verdict of the operational first choice, n_tested)."""
    if cells is None or len(cells) == 0:
        return -1, "none", 0
    order = np.argsort(-np.nan_to_num(cells["quality"].values, nan=-np.inf), kind="stable")
    first = None
    for n, i in enumerate(order, 1):
        v = rule_verdict(cells.iloc[i].to_dict(), rule, g_ref)
        if first is None:
            first = v
        if v.startswith("accept"):
            return int(i), first, n
    return -1, first, len(order)


def robust_sd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def qfmt(x, pct=True, digits=1):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not x.size:
        return "-"
    f = 100.0 if pct else 1.0
    return (f"{f * np.median(x):+.{digits}f} [{f * np.percentile(x, 25):+.{digits}f}, "
            f"{f * np.percentile(x, 75):+.{digits}f}]")


def load_unit(key):
    """Record rows of one unit (both variants)."""
    fr = []
    for f in sorted((Q.REC_DIR / key).glob("*.csv")):
        try:
            fr.append(pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str}))
        except Exception:  # noqa: BLE001
            continue
    df = pd.concat(fr, ignore_index=True)
    df["key"] = df["wmo"] + "_" + df["ident"]
    return df


def block_verdict(b, k=2.0, eps_tol=0.10, s_max=np.inf, include_below=False, rising="blocks",
                  min_blocks=1):
    """Cleanest-block verdict from q19_offline.block_arrays output.

    aerosol : some block (above the window; also below it if include_below) is darker than the window
              by more than eps_tol + k sigma (lower confidence bound of the excess > eps_tol);
    rising  : 'blocks' -> some block ABOVE is brighter than the window by more than eps_tol + k sigma;
              'none'   -> never rejected for rising R.
    unverifiable : no block above the window, or (s_max finite) no block >= 1 km above the window top
              with sigma(excess) <= s_max.
    Blocks below the window never count as 'rising': a boundary layer there is normal.
    """
    if b is None:
        return "reject_unverifiable"
    up = b["is_up"]
    if up.sum() < min_blocks:
        return "reject_unverifiable"
    if np.isfinite(s_max):
        far = up & (b["zc"] >= b["ze"] + 1000.0)
        if not np.any(b["sig"][far] <= s_max):
            return "reject_unverifiable"
    use = up | (~up if include_below else np.zeros_like(up))
    exc, sig = b["exc"], b["sig"]
    if np.max(exc[use] - k * sig[use]) > eps_tol:
        return "reject_aerosol"
    if rising == "blocks" and np.min(exc[up] + k * sig[up]) < -eps_tol:
        return "reject_rising"
    return "accept"


def local_ref(t_ref, y_ref, t_eval, days=45, min_n=6, min_side=2):
    """Expected tier-1 level at each evaluation date: OLS line through the reference nights within +-days
    (one 2.5-MAD clipping pass), evaluated at the date. A plain +-45 d median is biased around the
    solstice-peaked annual cycle when the evaluated nights cluster in June and the reference nights in
    May/July; the local line removes that first-order bias. Falls back to NaN when the reference nights
    are too few or one-sided. Returns (ref, n_used)."""
    tr = np.asarray(t_ref, "datetime64[D]").astype(float)
    yr = np.asarray(y_ref, float)
    te = np.asarray(t_eval, "datetime64[D]").astype(float)
    out = np.full(te.size, np.nan)
    nn = np.zeros(te.size, int)
    for k, t in enumerate(te):
        m = (np.abs(tr - t) <= days) & np.isfinite(yr)
        if m.sum() < min_n or (tr[m] < t).sum() < min_side or (tr[m] > t).sum() < min_side:
            continue
        x, y = tr[m] - t, yr[m]
        for _ in range(2):
            A = np.column_stack([np.ones(x.size), x])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            r = y - A @ c
            s = 1.4826 * np.median(np.abs(r - np.median(r)))
            keep = np.abs(r) <= max(2.5 * s, 1e-6)
            if keep.all() or keep.sum() < min_n:
                break
            x, y = x[keep], y[keep]
        out[k] = c[0]
        nn[k] = x.size
    return out, nn


def unit_offsets(g_unit, t_eval, y_eval, ref_mask=None, days=45):
    """Offsets of y_eval (ln C) vs the unit's tier-1 nights: (plain median, local-line reference)."""
    r = g_unit[(g_unit.tier == 1) & g_unit.ok]
    if ref_mask is not None:
        r = r[ref_mask.reindex(r.index).fillna(False).values]
    med, nmed = neighbour_offset(t_eval, y_eval, r.t.values, r.lnc.values, days)
    ref, nref = local_ref(r.t.values, r.lnc.values, t_eval, days)
    return med, y_eval - ref, nref


def twin_offsets(series, tier1_ok, eval_rows, days=45, min_n=3):
    """Same-night twin test.

    series   : dict key -> DataFrame(date, lnc) of the VERIFIED reference nights of each unit (tier 1)
    tier1_ok : same structure, the tier-1 nights used for the pair's reference ratio
    eval_rows: DataFrame(key, date, lnc_eval, label) of the nights to evaluate (e.g. tier-2 nights)
    For an evaluated night of unit X and every twin Y with a tier-1 constant the SAME night:
      offset = [ln C_X(night) - ln C_Y(night)] - median over +-days of [ln C_X - ln C_Y] on nights where
      both are tier 1. The annual gain cycle common to co-located CHM15k cancels in the ratio.
    Returns a DataFrame (key, twin, date, label, offset, n_ref)."""
    out = []
    for r in eval_rows.itertuples():
        wmo, ident = r.key.rsplit("_", 1)
        if wmo not in Q.TWINS:
            continue
        for tid in Q.TWINS[wmo]:
            if tid == ident:
                continue
            ky = f"{wmo}_{tid}"
            if ky not in tier1_ok or r.key not in tier1_ok:
                continue
            ty = tier1_ok[ky]
            same = ty[ty.date == r.date]
            if not len(same):
                continue
            tx = tier1_ok[r.key]
            both = tx.merge(ty, on="date", suffixes=("_x", "_y"))
            tb = pd.to_datetime(both.date)
            m = np.abs((tb - pd.Timestamp(r.date)).dt.days.values) <= days
            if m.sum() < min_n:
                continue
            ref = np.median((both.lnc_x - both.lnc_y).values[m])
            out.append(dict(key=r.key, twin=ky, date=r.date, label=r.label,
                            offset=(r.lnc_eval - same.lnc.values[0]) - ref, n_ref=int(m.sum())))
    return pd.DataFrame(out)


def series_stats(t, lnc):
    """Annual-harmonic amplitude [%], night-to-night robust scatter [%] (q16a definitions), n."""
    sys.path.insert(0, str(HERE))
    from q16a_recalibration import harmonic, n2n_scatter
    t = pd.to_datetime(pd.Series(t))
    tday = (t - pd.Timestamp("1970-01-01")).dt.days.values.astype(float)
    y = np.asarray(lnc, float)
    m = np.isfinite(y)
    if m.sum() < 40:
        return dict(amp=np.nan, peak=np.nan, n2n=np.nan, n=int(m.sum()))
    h = harmonic(tday[m], y[m])
    return dict(amp=h["amp"], peak=h["peak"], n2n=n2n_scatter(tday[m], y[m]), n=int(m.sum()))


def vc_verdict(b, k=2.0, eps_tol=0.10, tol_lo=0.10, s_max=np.inf):
    """FINAL vertical-consistency rule (cleanest block + consistency with the 1 km below the window).

    b = q19_offline.block_arrays(...): window excess R_win/R_block - 1 and its noise sigma for every
    750-m block above the window (to the SNR limit, <= 9 km) and below it (to max(1.5 km, zb - 1 km)).
      'unverifiable' : no block above the window (or, s_max finite, none >= 1 km above measured to s_max)
      'aerosol'      : some block (above, or within 1 km below) DARKER than the window by > eps_tol + k sigma
      'tail'         : some block within 1 km below BRIGHTER than the window by > tol_lo + k sigma
                       (the window sits in an aerosol tail still decaying into it)
      'accept'       : otherwise. Blocks ABOVE that are brighter (cirrus, elevated layers) never reject.
    """
    if b is None:
        return "unverifiable"
    up = b["is_up"]
    if up.sum() < 1:
        return "unverifiable"
    if np.isfinite(s_max):
        far = up & (b["zc"] >= b["ze"] + 1000.0)
        if not np.any(b["sig"][far] <= s_max):
            return "unverifiable"
    exc, sig = b["exc"], b["sig"]
    if np.max(exc - k * sig) > eps_tol:
        return "aerosol"
    lo = ~up
    if lo.any() and np.isfinite(tol_lo) and np.min(exc[lo] + k * sig[lo]) < -tol_lo:
        return "tail"
    return "accept"


def vc_verdict2(b, k=2.0, eps_tol=0.10, tol_lo=0.10):
    """Revised rule (q19 final): as vc_verdict, but a window with NO measurable block above it (it reaches
    the SNR limit) is not rejected for that alone -- its flatness inside the window is already tested by
    the tier-2 chi2red gate once the noise is right -- and is judged on the blocks within 1 km below it.
    'unverifiable' only when there is no block at all (window from the 1.5 km floor to the SNR limit)."""
    if b is None or b["exc"].size == 0:
        return "unverifiable"
    exc, sig, up = b["exc"], b["sig"], b["is_up"]
    if np.max(exc - k * sig) > eps_tol:
        return "aerosol"
    lo = ~up
    if lo.any() and np.isfinite(tol_lo) and np.min(exc[lo] + k * sig[lo]) < -tol_lo:
        return "tail"
    return "accept"
