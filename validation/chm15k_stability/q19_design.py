# -*- coding: utf-8 -*-
"""q19 design step: evaluate candidate gate rules OFFLINE on the record run (no recalibration).

Inputs: the record run (q19_gradient_gate.py record) and the offline block statistics
(q19_offline.py -> q19_blocks_windows.pkl, q19_blocks_cells.pkl).

For every CHM15k unit (Payerne A apart) and both noise bases ("op" = operational v2.2 with its
over-estimated sigma feeding the existing tier-2 gates; "fix" = the same with the measured noise):
  * tier-2 nights and their offset vs the unit's tier-1 nights within +-45 d, both as the plain median
    (q16e definition) and against a local line through the tier-1 nights (removes the bias of the
    solstice-peaked annual cycle on nights clustered in June);
  * the verdict of each candidate rule on the operational first-choice window and after the offline
    gated re-selection over the stored relaxed cells (C of a re-selected window estimated as
    C x ratio_med(new) / ratio_med(old); exact constants come from the `gate` run);
  * the false-rejection rate on CLEAN tier-1 nights (Nov-Mar, fully clear, window bottom <= 3 km,
    constant within +-10 % of the local level): the window the noise-aware tier WOULD choose there.
Writes q19_design_rules.csv (one row per rule x base) and q19_design_tier2_<base>.csv.
"""
from __future__ import annotations

import itertools
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q19 = L.Q.Q19
warnings.filterwarnings("ignore")


def clean_mask(op):
    """Clean tier-1 nights (same definition as the injection set; no gate statistic involved)."""
    out = pd.Series(False, index=op.index)
    for key, g in op.groupby("key"):
        ok = (g.ok & (g.tier == 1)).values
        t = g.t.values.astype("datetime64[D]").astype(float)
        lc = g.lnc.values
        loc = np.full(len(g), np.nan)
        for k in range(len(g)):
            m = ok & (np.abs(t - t[k]) <= 30)
            if m.sum() >= 3:
                loc[k] = np.median(lc[m])
        m = (ok & g.month.isin([11, 12, 1, 2, 3]).values & (g.flag.values == 1.0)
             & (g.zb.values <= 3000) & (np.abs(lc - loc) <= 0.10))
        out.loc[g.index[m]] = True
    return out


def tier2_frame(df, base):
    """Tier-2 nights of one noise base with offsets vs the unit's tier-1 nights (op rows = tier 1 is
    identical in both bases)."""
    rows = []
    for key, g in df.groupby("key"):
        gb = g[g.variant == base]
        gop = g[g.variant == "op"]
        t2 = gb[(gb.tier == 2) & gb.ok]
        if not len(t2):
            continue
        med, loc, nref = L.unit_offsets(gop, t2.t.values, t2.lnc.values)
        rows.append(t2.assign(off_med=med, off_loc=loc, nref=nref))
    return pd.concat(rows) if rows else pd.DataFrame()


# ------------------------------------------------------------------------------------------------
# rule families
# ------------------------------------------------------------------------------------------------
def grad_rules():
    out = []
    for stat, tol, k, s_max in itertools.product(("up", "ext"), (3.0, 5.0, 8.0), (2.0, 3.0),
                                                 (np.inf, 4.0, 2.5)):
        out.append(dict(family="grad", stat=stat, g_tol=tol, k=k, s_max=s_max, min_above_m=0.0,
                        reject_rising=True,
                        name=f"grad_{stat}_tol{tol:g}_k{k:g}_s{s_max:g}"))
    return out


def block_rules():
    out = []
    for k, tol, below, s_max, rising in itertools.product((2.0, 3.0), (0.05, 0.10, 0.15), (False, True),
                                                          (np.inf, 0.15, 0.10, 0.07), ("none", "blocks")):
        out.append(dict(family="block", k=k, eps_tol=tol, include_below=below, s_max=s_max, rising=rising,
                        name=f"blk_k{k:g}_tol{tol:g}_{'below' if below else 'above'}_s{s_max:g}_{rising}"))
    return out


def verdict_of(rule, b=None, stats=None):
    if rule["family"] == "block":
        return L.block_verdict(b, k=rule["k"], eps_tol=rule["eps_tol"], s_max=rule["s_max"],
                               include_below=rule["include_below"], rising=rule["rising"])
    return L.rule_verdict(stats, rule)


def summarise(off, acc):
    fin = np.isfinite(off)
    a = fin & acc
    r = fin & ~acc
    return dict(n=int(fin.sum()), acc=float(a.sum() / max(fin.sum(), 1)),
                off_acc=float(np.median(off[a])) if a.any() else np.nan,
                mean_acc=float(np.mean(off[a])) if a.any() else np.nan,
                sd_acc=L.robust_sd(off[a]),
                p90_acc=float(np.percentile(np.abs(off[a]), 90)) if a.any() else np.nan,
                off_rej=float(np.median(off[r])) if r.any() else np.nan)


def main():
    df = L.load(("CHM15k",))
    df = df[df.key != L.PAYERNE_KEY]
    op = df[df.variant == "op"]
    clean = clean_mask(op)
    print(f"units {op.key.nunique()}, nights {len(op)}, tier-1 {int(((op.tier == 1) & op.ok).sum())}, "
          f"clean tier-1 {int(clean.sum())}", flush=True)
    W = pickle.load(open(Q19 / "q19_blocks_windows.pkl", "rb"))
    C = pickle.load(open(Q19 / "q19_blocks_cells.pkl", "rb"))
    Wd = {(w["key"], w["date"], w["variant"], w["which"]): w for w in W}
    T2 = {b: tier2_frame(df, b) for b in ("op", "fix")}
    for b, t2 in T2.items():
        t2.to_csv(Q19 / f"q19_design_tier2_{b}.csv", index=False)
        for col in ("off_med", "off_loc"):
            o = t2[col].values
            print(f"tier-2 {b}: n {len(t2)}, {col} n {np.isfinite(o).sum()}: median {100 * np.nanmedian(o):+.1f} %, "
                  f"mean {100 * np.nanmean(o):+.1f} %, IQR {100 * np.nanpercentile(o, 25):+.1f}..{100 * np.nanpercentile(o, 75):+.1f}")
    ck = set(zip(op[clean].key, op[clean].date))
    res = []
    for rule in grad_rules() + block_rules():
        for base in ("op", "fix"):
            out = dict(rule=rule["name"], family=rule["family"], base=base,
                       **{k: v for k, v in rule.items() if k not in ("name", "family")})
            # false rejection on clean tier-1 nights: window the noise-aware tier would choose
            fb = df[(df.variant == base)]
            fb = fb[[(k, d) in ck for k, d in zip(fb.key, fb.date)]]
            v = []
            for r in fb.itertuples():
                if rule["family"] == "block":
                    v.append(verdict_of(rule, b=Wd.get((r.key, r.date, base, "t2"))))
                else:
                    st = {c[3:]: getattr(r, c) for c in fb.columns if c.startswith("t2_")}
                    v.append(verdict_of(rule, stats=st))
            v = np.array(v)
            fin = np.array([isinstance(x, str) for x in v])
            out["frr_clean"] = float(np.mean(~L.accepted(v)))
            out["frr_clean_aer"] = float(np.mean(v == "reject_aerosol"))
            out["frr_clean_unv"] = float(np.mean(v == "reject_unverifiable"))
            t2 = T2[base]
            vf, kept, newoff_loc, newoff_med = [], [], [], []
            for r in t2.itertuples():
                if rule["family"] == "block":
                    vv = verdict_of(rule, b=Wd.get((r.key, r.date, base, "sel")))
                else:
                    st = {c[4:]: getattr(r, c) for c in t2.columns if c.startswith("sel_")}
                    vv = verdict_of(rule, stats=st)
                vf.append(vv)
                # offline re-selection over the relaxed cells
                cc = C.get((r.key, r.date, base))
                ok, dl = False, np.nan
                if cc is not None:
                    cells, blist = cc
                    order = np.argsort(-np.nan_to_num(cells["quality"].values, nan=-np.inf), kind="stable")
                    for i in order:
                        if rule["family"] == "block":
                            vv2 = verdict_of(rule, b=blist[i])
                        else:
                            vv2 = verdict_of(rule, stats=cells.iloc[i].to_dict())
                        if vv2.startswith("accept"):
                            rm_new, rm_old = cells["ratio_med"].values[i], r.sel_ratio_med
                            dl = np.log(rm_new / rm_old) if (rm_new > 0 and rm_old > 0) else 0.0
                            ok = True
                            break
                kept.append(ok)
                newoff_loc.append(r.off_loc + dl if ok else np.nan)
                newoff_med.append(r.off_med + dl if ok else np.nan)
            vf = np.array(vf)
            acc = L.accepted(vf)
            for col, o in (("med", t2.off_med.values), ("loc", t2.off_loc.values)):
                s = summarise(o, acc)
                out.update({f"first_{col}_{k}": v for k, v in s.items()})
            for cat in ("reject_aerosol", "reject_rising", "reject_unverifiable"):
                m = (vf == cat) & np.isfinite(t2.off_loc.values)
                out[f"n_{cat}"] = int(m.sum())
                out[f"off_loc_{cat}"] = float(np.median(t2.off_loc.values[m])) if m.any() else np.nan
            kept = np.array(kept)
            for col, o in (("med", np.array(newoff_med)), ("loc", np.array(newoff_loc))):
                s = summarise(o, kept)
                out.update({f"resel_{col}_{k}": v for k, v in s.items()})
            out["resel_kept"] = float(np.mean(kept))
            res.append(out)
    R = pd.DataFrame(res)
    R.to_csv(Q19 / "q19_design_rules.csv", index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_rows", 400)
    cols = ["rule", "base", "frr_clean", "frr_clean_unv", "first_loc_acc", "first_loc_off_acc",
            "first_loc_off_rej", "first_loc_sd_acc", "resel_kept", "resel_loc_off_acc", "resel_loc_mean_acc",
            "resel_loc_sd_acc", "resel_loc_p90_acc", "resel_med_off_acc"]
    print(R[cols].round(3).to_string())


if __name__ == "__main__":
    main()
