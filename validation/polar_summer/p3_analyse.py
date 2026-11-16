# -*- coding: utf-8 -*-
"""p3 - polar summer: analysis of the p2 recalibration runs (validation, availability, bias versus
solar elevation, solar-background residual, continuity with the dark-night series, QC threshold).

Reads C:/DATA/Projects/202610_v23/polar/runs/*/*.csv (p2_recalibrate.py) and writes tables and
figures to C:/DATA/Projects/202610_v23/polar/ (p3_*.csv, p3_*.png, p3_digest.txt).

Conventions
  * ln ratios are reported in % (100 * ln), medians with bootstrap 95 % CI.
  * Temperature: the CHM15k constant follows the internal temperature (report 19): same-night ratios
    are corrected with each unit's intra-night kappa (q2a_intranight.csv, b_T_ctrl, %/K; fallback
    1.5 %/K); the CL61 is not corrected (no kappa established).
  * 'dark' = SZA > 100 (the operational night); 'low-sun' = the v2.3 fallback window.
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import load_constants, chm15k_streams  # noqa: E402

OUT = Path("C:/DATA/Projects/202610_v23/polar")
RUN_DIR = OUT / "runs"
CENSUS = REPO / "validation/scope_l1_2026_census.json"
KAPPA_CSV = Path("C:/DATA/Projects/202610_CHM15k_stability/q2a_intranight.csv")
EPOCH = datetime(1970, 1, 1)
OK_FLAGS = (1.0, 0.5)
LOW_SUN = 0.25
BANDS = ["b80_85", "b85_90", "b90_94", "b94_98", "b98_100"]
BAND_MID = {"b80_85": 82.5, "b85_90": 87.5, "b90_94": 92.0, "b94_98": 96.0, "b98_100": 99.0}
THR_RUNS = ["t98", "t96", "t94", "t92", "t90"]
XCOLS = ["x_2_3", "x_3_4", "x_4_5", "x_5_6", "x_6_8", "x_8_10", "x_10_12", "x_12_14"]
XMID = [2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.0]
POLAR_KEYS = ["0-20000-0-01062_A", "0-20000-0-01028_A", "0-20000-0-01001_A", "0-20000-0-01010_A",
              "0-208-0-06011_A", "0-20000-0-01311_A", "0-20000-0-03005_A", "0-20000-0-01492_A",
              "0-20000-0-03018_A", "0-20000-0-03023_A", "0-20008-0-BIR_A"]

warnings.filterwarnings("ignore")
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def boot_median(x, nboot=2000, seed=1):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return (float(np.median(x)) if x.size else np.nan), np.nan, np.nan
    rng = np.random.default_rng(seed)
    b = np.median(rng.choice(x, size=(nboot, x.size), replace=True), axis=1)
    return float(np.median(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def robust_sd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


FINAL_DIR = OUT / "runs_final"
LOWSUN_RUN = "final"          # the recommended configuration (CalibrationOptions defaults, option on)


def load_runs():
    frames = []
    for f in sorted(RUN_DIR.glob("*/*.csv")):
        try:
            frames.append(pd.read_csv(f, dtype={"date": str, "ident": str, "module": str}))
        except Exception:  # noqa: BLE001 - empty chunk
            continue
    df = pd.concat(frames, ignore_index=True)
    # 'final' = the recommended defaults (6 h centred on solar midnight, SZA >= 78, window >=
    # min_time_range). It selects exactly the profiles of the evaluated 'pol_h6' variant, except on
    # the Hopen / Bjornoya solstice nights where the SZA floor shortens the window: pol_h6 rejected
    # those (window < 6 h) and they were re-run with the final code (p2_recalibrate.py final).
    fin = df[df.run == "pol_h6"].copy()
    fin["run"] = "final"
    re = [pd.read_csv(f, dtype={"date": str, "ident": str, "module": str})
          for f in sorted(FINAL_DIR.glob("*.csv"))] if FINAL_DIR.exists() else []
    if re:
        re = pd.concat(re, ignore_index=True)
        k_re = set(zip(re.wmo, re.ident, re.date))
        fin = fin[[(w, i, d) not in k_re for w, i, d in zip(fin.wmo, fin.ident, fin.date)]]
        fin = pd.concat([fin, re], ignore_index=True)
    df = pd.concat([df, fin], ignore_index=True)
    df["key"] = df["wmo"] + "_" + df["ident"]
    df["ok"] = df["flag"].isin(OK_FLAGS + (LOW_SUN,)) & np.isfinite(df["c"]) & (df["c"] > 0)
    df["lnc"] = np.where(df["ok"], np.log(df["c"].where(df["c"] > 0)), np.nan)
    df["tday"] = [(datetime.strptime(d, "%Y%m%d") - EPOCH).days for d in df["date"]]
    df["doy"] = [datetime.strptime(d, "%Y%m%d").timetuple().tm_yday for d in df["date"]]
    df["polar_unit"] = df["key"].isin(POLAR_KEYS)
    return df


def census():
    return {f"{s['wmo']}_{s['ident']}": s for s in json.load(open(CENSUS, encoding="utf-8"))}


def kappas():
    k = {}
    if KAPPA_CSV.exists():
        for r in csv.DictReader(open(KAPPA_CSV, encoding="utf-8")):
            try:
                k[r["stream"]] = float(r["b_T_ctrl"])
            except ValueError:
                pass
    return k


def night_dates_from_constants(con, lon):
    """Night key D of each published constant: the solar date of start_time + 12 h (q16a)."""
    t0 = np.asarray(con["t0"], float)
    d = np.floor(t0 + lon / 360.0 + 0.5)
    return [(EPOCH + timedelta(days=float(x))).strftime("%Y%m%d") for x in d]


# --------------------------------------------------------------------------------------------
# A. validation
# --------------------------------------------------------------------------------------------
def validation(df, cen):
    rows = []
    op = df[df.run == "op"]
    for key, g in op.groupby("key"):
        s = cen[key]
        con = load_constants(s["wmo"], s["ident"])
        if con is None:
            continue
        pub = dict(zip(night_dates_from_constants(con, s["lon"]), con["c"]))
        p = np.array([pub.get(d, np.nan) for d in g.date])
        both = g.ok.values & np.isfinite(p)
        r = g.c.values[both] / p[both]
        rows.append(dict(key=key, site=s["site"], type=s["type"], lat=s["lat"], n_nights=len(g),
                         n_pub=int(np.isfinite(p).sum()), n_ours=int(g.ok.sum()), n_both=int(both.sum()),
                         frac_identical=float(np.mean(np.abs(r - 1) < 1e-9)) if r.size else np.nan,
                         ratio_median=float(np.median(r)) if r.size else np.nan,
                         max_abs_dev_pct=float(100 * np.max(np.abs(r - 1))) if r.size else np.nan))
    val = pd.DataFrame(rows)
    val.to_csv(OUT / "p3_validation.csv", index=False)
    say("\n== A. validation: op vs published v2.2 constants ==")
    say(val.to_string())
    say(f"all units: {int(val.n_both.sum())} common nights, "
        f"{100 * np.nansum(val.frac_identical * val.n_both) / max(1, val.n_both.sum()):.1f} % bit-identical")
    # option ON with a real night must equal the option OFF (real-data bit identity)
    w = df[df.run.isin(["op", "polar"])].pivot_table(index=["key", "date"], columns="run",
                                                      values=["flag", "c", "dark_h"], aggfunc="first")
    dark = w[("dark_h", "op")] > 0
    same_flag = w[("flag", "op")][dark] == w[("flag", "polar")][dark]
    c_op, c_po = w[("c", "op")][dark], w[("c", "polar")][dark]
    same_c = (c_op == c_po) | (c_op.isna() & c_po.isna())
    say(f"option ON vs OFF on {int(dark.sum())} nights with darkness: flags identical "
        f"{int(same_flag.sum())}/{int(dark.sum())}, constants identical {int(same_c.sum())}/{int(dark.sum())}")
    return val


# --------------------------------------------------------------------------------------------
# B. availability
# --------------------------------------------------------------------------------------------
def availability(df, cen):
    say("\n== B. availability on nights without darkness (no profile with SZA > 100) ==")
    nn = df[(df.dark_h <= 0) & df.polar_unit]
    rows = []
    for key, g in nn.groupby("key"):
        s = cen[key]
        nights = g.date.unique()
        row = dict(key=key, site=s["site"], type=s["type"], lat=s["lat"], n_no_night=len(nights))
        for run in ["op"] + THR_RUNS + ["polar", "pol_h2", "pol_h6", "final"]:
            x = g[g.run == run]
            row[run] = int(x.ok.sum())
        # nights where ANY variant could calibrate = an upper bound set by the clouds
        row["any"] = int(g[g.ok].date.nunique())
        rows.append(row)
    av = pd.DataFrame(rows).sort_values("lat", ascending=False)
    tot = av.drop(columns=["key", "site", "type", "lat"]).sum()
    av.to_csv(OUT / "p3_availability.csv", index=False)
    say(av.to_string())
    say("total:", tot.to_dict())
    # whole-period gain: successful nights per unit with and without the fallback
    g2 = []
    for key, g in df[df.polar_unit].groupby("key"):
        op = g[(g.run == "op") & g.ok]
        # with the option on: the operational result where a night exists, the low-sun one elsewhere
        po = pd.concat([op[op.dark_h > 0], g[(g.run == "final") & g.ok & (g.dark_h <= 0)]])
        row = dict(key=key, site=cen[key]["site"], nights=g.date.nunique(), op=len(op),
                   polar_total=len(po), gain=len(po) - len(op),
                   gain_pct=100.0 * (len(po) - len(op)) / max(1, len(op)))
        # longest stretch without any calibration that overlaps the summer (Apr-Sep), per year
        for name, src in (("gap_op", op), ("gap_polar", po)):
            gaps = []
            for yr in (2025, 2026):
                t = np.sort(src.tday.values)
                lo = (datetime(yr, 3, 1) - EPOCH).days
                hi = (datetime(yr, 10, 31) - EPOCH).days
                last = (datetime(2026, 8, 31) - EPOCH).days
                have = np.unique(g.tday.values)
                inyr = have[(have >= lo) & (have <= hi)]
                if inyr.size < 100:                                # no summer record that year
                    gaps.append(np.nan)
                    continue
                t = np.r_[lo, t[(t >= lo) & (t <= hi)], min(hi, last, inyr.max())]
                gaps.append(int(np.max(np.diff(t))) if t.size > 1 else np.nan)
            row[name + "_2025"], row[name + "_2026"] = gaps
        g2.append(row)
    gain = pd.DataFrame(g2)
    gain.to_csv(OUT / "p3_gain.csv", index=False)
    say(gain.to_string())
    return av, gain


# --------------------------------------------------------------------------------------------
# C. bias versus solar elevation (same night, SZA bands vs SZA > 100)
# --------------------------------------------------------------------------------------------
def band_pairs(df, kap):
    typ_of = df.groupby("key")["type"].first().to_dict()
    side = df[df.run.isin(BANDS + ["b100"])]
    wide = side.pivot_table(index=["key", "date"], columns="run",
                            values=["lnc", "tint_mean", "bg_mean", "resid_pct", "resid_q30",
                                    "resid_tmed", "sza_mean", "n_used", "u", "c", "err_pct"]
                            + XCOLS, aggfunc="first")
    rows = []
    for (key, date), r in wide.iterrows():
        ref = r[("lnc", "b100")] if ("lnc", "b100") in r.index else np.nan
        if not np.isfinite(ref):
            continue
        typ = typ_of.get(key, "CHM15k")
        k = kap.get(key, 1.5) if typ == "CHM15k" else 0.0
        for b in BANDS:
            if ("lnc", b) not in r.index or not np.isfinite(r[("lnc", b)]):
                continue
            dT = r[("tint_mean", b)] - r[("tint_mean", "b100")]
            lnr = 100.0 * (r[("lnc", b)] - ref)
            row = dict(key=key, date=date, band=b, sza=r[("sza_mean", b)], typ=typ,
                       lnr=lnr, dT=dT, kappa=k,
                       lnr_T=lnr - k * dT if np.isfinite(dT) else np.nan,
                       bg_ratio=r[("bg_mean", b)] / r[("bg_mean", "b100")],
                       d_resid=r[("resid_pct", b)] - r[("resid_pct", "b100")],
                       d_resid_q30=r[("resid_q30", b)] - r[("resid_q30", "b100")],
                       d_resid_tmed=r[("resid_tmed", b)] - r[("resid_tmed", "b100")],
                       rel_u=r[("u", b)] / r[("c", b)], rel_u_ref=r[("u", "b100")] / r[("c", "b100")],
                       err=r[("err_pct", b)], err_ref=r[("err_pct", "b100")])
            for xc in XCOLS:
                row[f"d{xc}"] = r[(xc, b)] - r[(xc, "b100")]
            rows.append(row)
    return pd.DataFrame(rows)


def bias_vs_sza(df, kap, cen):
    say("\n== C. same-night bias versus SZA band (relative to SZA > 100), T_int-corrected ==")
    bp = band_pairs(df, kap)
    bp["polar_unit"] = bp.key.isin(POLAR_KEYS)
    bp["lat"] = [cen[k]["lat"] for k in bp.key]
    bp.to_csv(OUT / "p3_band_pairs.csv", index=False)
    rows = []
    for (grp, typ), g0 in bp.groupby([np.where(bp.polar_unit, "polar", "control"), "typ"]):
        for b in BANDS:
            g = g0[g0.band == b]
            if len(g) < 5:
                continue
            m, lo, hi = boot_median(g.lnr_T)
            mr, _, _ = boot_median(g.lnr)
            per_unit = g.groupby("key").lnr_T.median()
            rows.append(dict(group=grp, typ=typ, band=b, n_nights=len(g), n_units=g.key.nunique(),
                             median_lnr=mr, median_lnr_T=m, ci_lo=lo, ci_hi=hi,
                             unit_median_T=float(per_unit.median()) if len(per_unit) else np.nan,
                             sd_lnr_T=robust_sd(g.lnr_T), dT=float(np.nanmedian(g.dT)),
                             bg_ratio=float(np.nanmedian(g.bg_ratio)),
                             d_resid=float(np.nanmedian(g.d_resid)),
                             d_resid_q30=float(np.nanmedian(g.d_resid_q30)),
                             d_resid_tmed=float(np.nanmedian(g.d_resid_tmed)),
                             rel_u=float(np.nanmedian(g.rel_u)), rel_u_ref=float(np.nanmedian(g.rel_u_ref)),
                             err=float(np.nanmedian(g.err)), err_ref=float(np.nanmedian(g.err_ref))))
    bs = pd.DataFrame(rows)
    bs.to_csv(OUT / "p3_bias_sza.csv", index=False)
    say(bs.to_string())
    # Does the far-range residual predict the constant's change? (same night, subset - dark)
    from scipy.stats import spearmanr
    for est in ("d_resid", "d_resid_tmed"):
        g = bp[np.isfinite(bp[est]) & np.isfinite(bp.lnr_T)]
        if len(g) < 20:
            continue
        rho, p = spearmanr(g[est], g.lnr_T)
        q = pd.qcut(g[est], 5, duplicates="drop")
        tab = g.groupby(q, observed=True).lnr_T.agg(["median", "count"])
        say(f"{est} vs ln(C_subset/C_dark) T-corr: Spearman {rho:+.2f} (p={p:.2g}, n={len(g)})")
        say(tab.to_string())
    return bp, bs


def threshold_bias(df, kap):
    """Same-night C(threshold t) / C(op) on nights with a short dark period (both succeed)."""
    say("\n== C2. lowered threshold vs operational, same nights (dark < 6 h), T-corrected ==")
    w = df[df.run.isin(["op"] + THR_RUNS + ["pol_force"])].pivot_table(
        index=["key", "date"], columns="run", values=["lnc", "tint_mean", "dark_h"], aggfunc="first")
    rows = []
    for run in THR_RUNS + ["pol_force"]:
        if ("lnc", run) not in w.columns:
            continue
        m = np.isfinite(w[("lnc", run)]) & np.isfinite(w[("lnc", "op")])
        g = w[m]
        k = np.array([kap.get(key, 1.5) for key, _ in g.index])
        lnr = 100 * (g[("lnc", run)] - g[("lnc", "op")]).values
        dT = (g[("tint_mean", run)] - g[("tint_mean", "op")]).values
        lnr_T = lnr - k * np.nan_to_num(dT)
        med, lo, hi = boot_median(lnr_T)
        rows.append(dict(run=run, n=int(m.sum()), median_lnr=float(np.median(lnr)) if lnr.size else np.nan,
                         median_lnr_T=med, ci_lo=lo, ci_hi=hi, sd=robust_sd(lnr_T),
                         dark_h_median=float(np.nanmedian(g[("dark_h", "op")]))))
    tb = pd.DataFrame(rows)
    tb.to_csv(OUT / "p3_threshold_bias.csv", index=False)
    say(tb.to_string())
    return tb


def window_variants(df, kap):
    """Low-sun window length: 2 h (pol_h2) / 4 h or 3 deg (polar) / 6 h (final = recommended),
    on the same nights without darkness, relative to the recommended one."""
    say("\n== C4. low-sun window variants on nights without darkness (same night, vs final 6 h) ==")
    ref = LOWSUN_RUN
    nn = df[(df.dark_h <= 0) & df.run.isin(["polar", "pol_h2", ref])]
    w = nn.pivot_table(index=["key", "date"], columns="run",
                       values=["lnc", "tint_mean", "u", "c", "sza_min", "bg_mean"], aggfunc="first")
    rows = []
    for run in ("pol_h2", "polar", ref):
        okr = np.isfinite(w[("lnc", run)])
        row = dict(run=run, n_success=int(okr.sum()),
                   relu_median=float(np.nanmedian((w[("u", run)] / w[("c", run)])[okr])),
                   sza_min_median=float(np.nanmedian(w[("sza_min", run)][okr])))
        if run != ref:
            m = okr & np.isfinite(w[("lnc", ref)])
            k = np.array([kap.get(key, 1.5) for key, _ in w[m].index])
            lnr = 100 * (w[("lnc", run)][m] - w[("lnc", ref)][m]).values
            dT = (w[("tint_mean", run)][m] - w[("tint_mean", ref)][m]).values
            med, lo, hi = boot_median(lnr - k * np.nan_to_num(dT))
            row.update(n_both=int(m.sum()), median_vs_final_T=med, ci_lo=lo, ci_hi=hi,
                       rsd_vs_final=robust_sd(lnr),
                       bg_ratio_vs_final=float(np.nanmedian((w[("bg_mean", run)] / w[("bg_mean", ref)])[m])))
        rows.append(row)
    wv = pd.DataFrame(rows)
    wv.to_csv(OUT / "p3_window_variants.csv", index=False)
    say(wv.to_string())
    return wv


def trigger_study(df, kap):
    """Nights with a SHORT dark period (0 < dark < 3 h): operational (dark profiles only) versus the
    forced low-sun window (pol_force). Residuals are taken against a +-15-day running median of the
    unit's ordinary dark-night constants (dark >= 3 h), so both selections are judged against the
    same reference."""
    say("\n== C3. trigger: short-dark nights, operational vs forced low-sun window ==")
    rows = []
    out = []
    for key, g in df.groupby("key"):
        ref = g[(g.run == "op") & g.ok & (g.dark_h >= 3.0)][["tday", "lnc"]].sort_values("tday")
        sh = g[(g.dark_h > 0) & (g.dark_h < 3.0)]
        if sh.empty:
            continue
        for run in ("op", "pol_force"):
            x = sh[sh.run == run]
            for _, r in x.iterrows():
                m = np.abs(ref.tday.values - r.tday) <= 15
                base = np.median(ref.lnc.values[m]) if m.sum() >= 3 else np.nan
                out.append(dict(key=key, date=r.date, run=run, dark_h=r.dark_h, ok=bool(r.ok),
                                resid=100 * (r.lnc - base) if r.ok else np.nan,
                                relu=r.u / r.c if r.ok else np.nan))
    tr = pd.DataFrame(out)
    if tr.empty:
        say("no short-dark nights")
        return tr
    tr["dbin"] = pd.cut(tr.dark_h, [0, 0.5, 1.0, 2.0, 3.0], right=False)
    for (dbin, run), g in tr.groupby(["dbin", "run"], observed=True):
        rows.append(dict(dark_h=str(dbin), run=run, n_nights=len(g), success=float(g.ok.mean()),
                         n_ok=int(g.ok.sum()), resid_median=float(np.nanmedian(g.resid)),
                         resid_rsd=robust_sd(g.resid), relu_median=float(np.nanmedian(g.relu))))
    ts = pd.DataFrame(rows)
    ts.to_csv(OUT / "p3_trigger.csv", index=False)
    say(ts.to_string())
    return ts


# --------------------------------------------------------------------------------------------
# D. continuity with the dark-night series (spring / autumn) and the network seasonal course
# --------------------------------------------------------------------------------------------
def network_seasonal_expectation(pairs_of_windows):
    """Network median of ln(C in DOY window S / C in DOY window R) over mid-latitude CHM15k units
    (published v2.2 constants, lat < 55 N, main optical module), per (S, R) window pair."""
    out = {}
    units = [(w, i, la) for (w, i, site, la, lo) in chm15k_streams() if la < 55.0]
    cache = {}
    for wmo, ident, la in units:
        con = load_constants(wmo, ident)
        if con is None or len(con["c"]) < 100:
            continue
        mod = max(set(con["module"]), key=list(con["module"]).count)
        m = (con["module"] == mod) & (con["c"] > 0)
        t = con["t"][m]
        doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in t])
        cache[(wmo, ident)] = (doy, np.log(con["c"][m]))
    for name, (s_win, r_win) in pairs_of_windows.items():
        vals = []
        for (doy, lnc) in cache.values():
            ms = _in_window(doy, s_win)
            mr = _in_window(doy, r_win)
            if ms.sum() >= 5 and mr.sum() >= 10:
                vals.append(100 * (np.median(lnc[ms]) - np.median(lnc[mr])))
        out[name] = (float(np.median(vals)) if vals else np.nan,
                     float(np.percentile(vals, 25)) if vals else np.nan,
                     float(np.percentile(vals, 75)) if vals else np.nan, len(vals))
    return out


def _in_window(doy, win):
    """win = (lo, hi) DOY window, or a list of such windows (union)."""
    doy = np.asarray(doy)
    if isinstance(win, list):
        return np.any([_in_window(doy, w) for w in win], axis=0)
    lo, hi = win
    return (doy >= lo) & (doy <= hi)


def continuity(df, cen, kap):
    """Per polar unit: low-sun constants (no-night season) vs the dark-night series around it."""
    from monitoring.kalman import kalman_best_estimate
    say("\n== D. continuity: low-sun constants vs the dark-night series (same optical module) ==")
    rows = []
    windows = {}
    series = {}
    for key in POLAR_KEYS:
        g = df[df.key == key]
        if g.empty:
            continue
        op = g[(g.run == "op") & g.ok].sort_values("tday")
        po = g[(g.run == LOWSUN_RUN) & g.ok & (g.flag == LOW_SUN)].sort_values("tday")
        series[key] = (op, po, np.sort(g[g.dark_h <= 0].tday.unique()))
        if len(po) < 3 or len(op) < 10:
            rows.append(dict(key=key, site=cen[key]["site"], n_op=len(op), n_low_sun=len(po)))
            continue
        # main optical module of the low-sun nights, compare within it
        mod = po.module.mode().iloc[0] if po.module.notna().any() else None
        if mod is not None and op.module.notna().any():
            op_m = op[op.module == mod]
        else:
            op_m = op
        nn_doys = g[g.dark_h <= 0].doy
        d_first, d_last = int(nn_doys.min()), int(nn_doys.max())
        spring = op_m[(op_m.doy < d_first) & (op_m.doy >= d_first - 45)]
        autumn = op_m[(op_m.doy > d_last) & (op_m.doy <= d_last + 45)]
        shoulder = pd.concat([spring, autumn])
        lvl = np.median(shoulder.lnc) if len(shoulder) >= 5 else np.nan
        windows[key] = ((d_first, d_last), [(d_first - 45, d_first - 1), (d_last + 1, d_last + 45)])
        # Kalman (dashboard filter) on the dark-night series, evaluated on the low-sun dates
        kal = np.nan
        try:
            grid, state, std = kalman_best_estimate(
                [EPOCH + timedelta(days=float(x)) for x in op_m.tday], op_m.c.values)
            if len(grid):
                gd = (grid.astype("datetime64[D]").astype(int)).astype(float)
                ks = np.interp(po.tday.values, gd, np.log(state))
                kal = float(np.median(100 * (po.lnc.values - ks)))
        except Exception:  # noqa: BLE001
            pass
        k = kap.get(key, 1.5) if cen[key]["type"] == "CHM15k" else 0.0
        rows.append(dict(
            key=key, site=cen[key]["site"], type=cen[key]["type"], lat=cen[key]["lat"], module=mod,
            n_op=len(op_m), n_low_sun=len(po), n_shoulder=len(shoulder),
            low_sun_vs_shoulder=100 * (np.median(po.lnc) - lvl) if np.isfinite(lvl) else np.nan,
            low_sun_vs_kalman=kal,
            n2n_low_sun=n2n(po.tday.values, po.lnc.values),
            n2n_op=n2n(op_m.tday.values, op_m.lnc.values),
            relu_low_sun=float(np.median(po.u / po.c)), relu_op=float(np.median(op_m.u / op_m.c)),
            tint_low_sun=float(np.nanmedian(po.tint_mean)), tint_shoulder=float(np.nanmedian(shoulder.tint_mean))
            if len(shoulder) else np.nan, kappa=k,
            d_first=d_first, d_last=d_last))
    cont = pd.DataFrame(rows)
    # network expectation of the same DOY contrast (seasonal CHM15k gain cycle)
    exp = network_seasonal_expectation(windows) if windows else {}
    cont["network_expected"] = [exp.get(k, (np.nan,))[0] for k in cont.key]
    cont["network_expected_iqr"] = [f"{exp[k][1]:+.1f}..{exp[k][2]:+.1f} ({exp[k][3]} units)"
                                    if k in exp else "" for k in cont.key]
    cont["excess_vs_network"] = cont.low_sun_vs_shoulder - cont.network_expected
    cont.to_csv(OUT / "p3_continuity.csv", index=False)
    say(cont.to_string())
    return cont, series


W_YEAR = 2 * np.pi / 365.25


def network_harmonic():
    """Median annual harmonic of ln C over mid-latitude CHM15k units (published v2.2 constants,
    lat < 55 N, main optical module): s(doy) = b cos(w doy) + c sin(w doy), in ln units. This is the
    instrumental CHM15k gain cycle of report 19 (peak near the solstice), used to remove the
    seasonal course when low-sun constants are compared with dark nights a few weeks away."""
    bs, cs = [], []
    for wmo, ident, site, la, lo in chm15k_streams():
        if la >= 55.0:
            continue
        con = load_constants(wmo, ident)
        if con is None or len(con["c"]) < 100:
            continue
        mod = max(set(con["module"]), key=list(con["module"]).count)
        m = (con["module"] == mod) & (con["c"] > 0)
        if m.sum() < 100:
            continue
        y = np.log(con["c"][m])
        doy = np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in con["t"][m]])
        med = np.median(y)
        keep = np.abs(y - med) < 4 * 1.4826 * np.median(np.abs(y - med))
        X = np.column_stack([np.ones(keep.sum()), np.cos(W_YEAR * doy[keep]), np.sin(W_YEAR * doy[keep])])
        coef, *_ = np.linalg.lstsq(X, y[keep], rcond=None)
        bs.append(coef[1])
        cs.append(coef[2])
    b, c = float(np.median(bs)), float(np.median(cs))
    amp = np.hypot(b, c)
    peak = (np.arctan2(c, b) / W_YEAR) % 365.25
    say(f"network CHM15k annual harmonic ({len(bs)} units, lat < 55): amplitude {100 * amp:.1f} %, "
        f"peak DOY {peak:.0f}")
    return lambda doy: b * np.cos(W_YEAR * np.asarray(doy, float)) + c * np.sin(W_YEAR * np.asarray(doy, float))


def local_continuity(df, cen, season, kap, window_days=45.0):
    """Each low-sun constant against the dark-night constants of the same unit and optical module
    within +-window_days (the edges of the no-night season), after removing the network seasonal
    course of C (CHM15k gain cycle, report 19). The same statistic for every dark night against its
    dark neighbours (leave-one-out) is the control: its median must be ~0 and its spread is the
    ordinary night-to-night scatter. lnr_T additionally removes kappa * (T_int - T_int of the
    neighbours): an upper bound of the temperature part, since the network course already contains
    the average seasonal temperature effect."""
    say(f"\n== D2. local continuity: constant vs dark nights within +-{window_days:.0f} d, "
        f"network seasonal cycle removed ==")
    rows = []
    for key in POLAR_KEYS:
        g = df[df.key == key]
        op = g[(g.run == "op") & g.ok]
        po = g[(g.run == LOWSUN_RUN) & g.ok & (g.flag == LOW_SUN)]
        if op.empty:
            continue
        k = kap.get(key, 1.5) if cen[key]["type"] == "CHM15k" else 0.0
        op_y = op.lnc.values - season(op.doy.values)
        for kind, src in (("low-sun", po), ("dark", op)):
            for idx, r in src.iterrows():
                m = (np.abs(op.tday.values - r.tday) <= window_days) & (op.module.values == r.module)
                if kind == "dark":
                    m &= op.index.values != idx
                if m.sum() < 2:
                    continue
                y = r.lnc - season(r.doy)
                lnr = 100 * (y - np.median(op_y[m]))
                dT = r.tint_mean - np.nanmedian(op.tint_mean.values[m])
                rows.append(dict(key=key, site=cen[key]["site"], kind=kind, date=r.date,
                                 n_ref=int(m.sum()), lnr=lnr,
                                 lnr_T=lnr - k * dT if np.isfinite(dT) else np.nan,
                                 dt_days=float(np.median(np.abs(op.tday.values[m] - r.tday))),
                                 relu=r.u / r.c, tint=r.tint_mean, dT=dT))
    lc = pd.DataFrame(rows)
    lc.to_csv(OUT / "p3_local_continuity.csv", index=False)
    out = []
    for (key, kind), g in lc.groupby(["key", "kind"]):
        m, lo, hi = boot_median(g.lnr)
        out.append(dict(key=key, site=cen[key]["site"], kind=kind, n=len(g), median=m, ci_lo=lo,
                        ci_hi=hi, median_T=float(np.nanmedian(g.lnr_T)), rsd=robust_sd(g.lnr),
                        relu=float(np.median(g.relu)), dT=float(np.nanmedian(g.dT))))
    for kind, g in lc.groupby("kind"):
        m, lo, hi = boot_median(g.lnr)
        mt, lot, hit = boot_median(g.lnr_T)
        out.append(dict(key="ALL", site="all polar units", kind=kind, n=len(g), median=m,
                        ci_lo=lo, ci_hi=hi, median_T=mt, ci_lo_T=lot, ci_hi_T=hit,
                        rsd=robust_sd(g.lnr), relu=float(np.median(g.relu)),
                        dT=float(np.nanmedian(g.dT))))
    ls = pd.DataFrame(out)
    ls.to_csv(OUT / "p3_local_continuity_summary.csv", index=False)
    say(ls.to_string())
    return lc, ls


def n2n(tday, lnc, window=31.0):
    """Night-to-night scatter [%]: robust sd of ln C around a running median (+-15 d)."""
    t = np.asarray(tday, float)
    y = np.asarray(lnc, float)
    if t.size < 8:
        return np.nan
    o = np.argsort(t)
    t, y = t[o], y[o]
    an = []
    for i in range(t.size):
        m = np.abs(t - t[i]) <= window / 2
        if m.sum() >= 4:
            an.append(y[i] - np.median(y[m]))
    return 100 * robust_sd(an) if len(an) >= 6 else np.nan


# --------------------------------------------------------------------------------------------
# E. residual QC threshold
# --------------------------------------------------------------------------------------------
def residual_qc(df):
    say("\n== E. far-range residual (implied window bias, %) on dark nights vs low-sun nights ==")
    dark = df[(df.run == "op") & df.ok]
    low = df[(df.run == LOWSUN_RUN) & df.ok & (df.flag == LOW_SUN)]
    rows = []
    for est in ("resid_pct", "resid_q30", "resid_tmed"):
        for name, g in (("dark (op)", dark), ("low-sun (polar)", low)):
            x = g[est].values
            x = x[np.isfinite(x)]
            if not x.size:
                continue
            rows.append(dict(estimator=est, sample=name, n=x.size, median=np.median(x),
                             rsd=robust_sd(x), p01=np.percentile(x, 1), p99=np.percentile(x, 99),
                             p995=np.percentile(np.abs(x), 99.5), frac_abs_gt10=np.mean(np.abs(x) > 10),
                             frac_abs_gt15=np.mean(np.abs(x) > 15), frac_abs_gt20=np.mean(np.abs(x) > 20)))
    qc = pd.DataFrame(rows)
    qc.to_csv(OUT / "p3_residual_qc.csv", index=False)
    say(qc.to_string())
    # rejection rate versus threshold (|implied bias| > thr), dark vs low-sun, per estimator
    rr = []
    for est in ("resid_pct", "resid_tmed"):
        for thr in (5, 10, 15, 20, 30, 50, 100):
            row = dict(estimator=est, threshold=thr)
            for name, g in (("dark", dark), ("low_sun", low)):
                x = g[est].values
                x = x[np.isfinite(x)]
                row[f"reject_{name}"] = float(np.mean(np.abs(x) > thr)) if x.size else np.nan
            rr.append(row)
    rr = pd.DataFrame(rr)
    rr.to_csv(OUT / "p3_residual_thresholds.csv", index=False)
    say(rr.to_string())
    # does the residual predict the constant? low-sun constants vs the local dark level
    return qc


# --------------------------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------------------------
def _site(s):
    return s.title().replace("_", " ")


def fig_availability(av):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(16, 5.6))
    av = av.sort_values("lat", ascending=False)
    x = np.arange(len(av))
    runs = [("t98", "threshold SZA > 98"), ("t96", "SZA > 96"), ("t94", "SZA > 94"),
            ("t92", "SZA > 92"), ("t90", "SZA > 90"), ("polar", "low-sun window 4 h"),
            (LOWSUN_RUN, "low-sun window 6 h (v2.3 default)")]
    cols = ["#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#9ecae1", "#2171b5"]
    w = 0.11
    for i, ((run, lab), c) in enumerate(zip(runs, cols)):
        ax.bar(x + (i - 3) * w, av[run], w, color=c, label=lab, edgecolor="k", lw=0.3)
    ax.plot(x, av["any"], "k_", ms=28, mew=2, label="nights calibratable by any variant (cloud limit)")
    for xi, (_, r) in zip(x, av.iterrows()):
        ax.text(xi, max(r["any"], r[LOWSUN_RUN]) + 1.5, f"{int(r.n_no_night)} nights\nwithout darkness",
                ha="center", fontsize=7, color="0.3")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_site(s)}\n{la:.1f} N {t}" for s, la, t in zip(av.site, av.lat, av.type)],
                       fontsize=8)
    ax.set_ylabel("calibrated nights (no profile with SZA > 100)")
    ax.set_title("Polar summer, 2025-01..2026-08: nights calibrated when there is no night "
                 "(operational: 0 everywhere)")
    ax.legend(fontsize=8, ncol=4, loc="upper right")
    ax.grid(alpha=0.3, axis="y")
    ax.set_ylim(0, max(av["any"].max(), av[LOWSUN_RUN].max()) * 1.35 + 5)
    fig.tight_layout()
    fig.savefig(OUT / "p3_availability.png", dpi=130)
    plt.close(fig)


def fig_bias(bp, bs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(19, 6), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    a = ax[0]
    groups = [("control", "CHM15k", "tab:blue", "CHM15k, controls 37-55 N"),
              ("polar", "CHM15k", "tab:red", "CHM15k, polar stations"),
              ("control", "CL61", "tab:green", "CL61, controls"),
              ("polar", "CL61", "tab:olive", "CL61, Birkenes")]
    for k, (grp, typ, col, lab) in enumerate(groups):
        g = bs[(bs.group == grp) & (bs.typ == typ)]
        if g.empty:
            continue
        xm = np.array([BAND_MID[b] for b in g.band]) + (k - 1.5) * 0.35
        a.errorbar(xm, g.median_lnr_T, yerr=[g.median_lnr_T - g.ci_lo, g.ci_hi - g.median_lnr_T],
                   fmt="o-", color=col, capsize=3, label=f"{lab} (n={int(g.n_nights.sum())})")
        sub = bp[(np.where(bp.polar_unit, "polar", "control") == grp) & (bp.typ == typ)]
        xs = np.array([BAND_MID[b] for b in sub.band]) + (k - 1.5) * 0.35
        a.plot(xs + np.random.default_rng(k).uniform(-0.12, 0.12, len(xs)), sub.lnr_T, ".",
               color=col, alpha=0.12, ms=3)
    a.axhline(0, color="0.4", lw=1)
    a.axvline(100, color="tab:red", lw=1)
    a.set_ylim(-30, 30)
    a.set_xlim(79, 101)
    a.set_xlabel("solar zenith angle of the subset [deg]")
    a.set_ylabel("ln(C_subset / C_dark) [%], T_int-corrected (kappa)")
    a.set_title("Same night, SZA subset vs SZA > 100 (median, 95 % CI)")
    a.legend(fontsize=8, loc="lower right")
    a.grid(alpha=0.3)

    a = ax[1]
    for k, (grp, typ, col, lab) in enumerate(groups[:2]):
        g = bs[(bs.group == grp) & (bs.typ == typ)]
        if g.empty:
            continue
        xm = np.array([BAND_MID[b] for b in g.band])
        a.plot(xm, g.bg_ratio, "o-", color=col, label=f"background ratio, {lab}")
    a.set_xlabel("solar zenith angle of the subset [deg]")
    a.set_ylabel("bckgrd_rcs_0 (subset) / bckgrd_rcs_0 (SZA > 100), same night")
    a.set_ylim(0.9, 3.0)
    a.axhline(1, color="0.4", lw=1)
    a2 = a.twinx()
    for k, (grp, typ, col, lab) in enumerate(groups[:2]):
        g = bs[(bs.group == grp) & (bs.typ == typ)]
        if g.empty:
            continue
        xm = np.array([BAND_MID[b] for b in g.band])
        a2.plot(xm, g.d_resid_tmed, "s--", color=col, alpha=0.7,
                label=f"far-range residual (window bias), {lab}")
    a2.set_ylabel("far-range residual, subset - dark [% at the window]")
    a.set_title("Solar background and far-range residual (medians)")
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper right")
    a.grid(alpha=0.3, which="both")
    a.set_xlim(79, 101)

    a = ax[2]
    cmap = plt.cm.plasma
    chm = bp[bp.typ == "CHM15k"]
    for i, b in enumerate(BANDS):
        g = chm[chm.band == b]
        if len(g) < 10:
            continue
        prof = [float(np.nanmedian(g[f"d{xc}"])) for xc in XCOLS]
        a.plot(100 * np.array(prof), XMID, "o-", color=cmap(i / len(BANDS)),
               label=f"SZA {b[1:].replace('_', '-')} (n={len(g)})")
    a.axvline(0, color="0.4", lw=1)
    a.axhspan(2.5, 5.5, color="0.9", zorder=0)
    a.text(0.98, 4.0, "typical\ncalibration\nwindow", transform=a.get_yaxis_transform(), ha="right",
           va="center", fontsize=7, color="0.4")
    a.set_xlabel("X(subset) - X(SZA > 100) [%-pt], X = signal / (slope p_mol)")
    a.set_ylabel("range [km]")
    a.set_title("Excess of the molecular-normalised profile (CHM15k)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "p3_bias_sza.png", dpi=130)
    plt.close(fig)


def fig_continuity(series, ls, cen):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from monitoring.kalman import kalman_best_estimate
    show = [k for k in ("0-20000-0-01062_A", "0-20000-0-01010_A", "0-20000-0-01311_A",
                        "0-20000-0-01492_A") if k in series]
    fig = plt.figure(figsize=(19, 8.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.8])
    for i, key in enumerate(show):
        a = fig.add_subplot(gs[i // 2, i % 2])
        op, po, nn = series[key]
        # normalise each point by the median dark-night constant of ITS optical module (Oslo swapped
        # modules in March 2026, a x3 level step)
        mref = op.groupby("module").c.median().to_dict()
        op = op.assign(c=op.c / op.module.map(mref))
        po = po.assign(c=po.c / po.module.map(mref).fillna(np.median(list(mref.values()) or [1.0])))
        ref = 1.0
        if len(nn):                              # shade the nights without darkness (SZA <= 100)
            brk = np.where(np.diff(nn) > 1)[0]
            starts = np.r_[nn[0], nn[brk + 1]]
            ends = np.r_[nn[brk], nn[-1]]
            for s0, s1 in zip(starts, ends):
                a.axvspan(EPOCH + timedelta(days=float(s0)), EPOCH + timedelta(days=float(s1) + 1),
                          color="gold", alpha=0.15, lw=0)
        a.plot([EPOCH + timedelta(days=float(t)) for t in op.tday], op.c / ref, "o", ms=4,
               color="k", label=f"dark night (SZA > 100), n={len(op)}")
        a.plot([EPOCH + timedelta(days=float(t)) for t in po.tday], po.c / ref, "o", ms=5,
               color="tab:orange", mec="k", mew=0.4, label=f"low-sun fallback, n={len(po)}")
        try:
            grid, state, std = kalman_best_estimate(
                [EPOCH + timedelta(days=float(x)) for x in op.tday], op.c.values)
            if len(grid):
                a.plot(grid.astype("datetime64[s]").astype(datetime), state / ref, "-",
                       color="0.5", lw=1, label="Kalman (dashboard filter) of the dark nights")
        except Exception:  # noqa: BLE001
            pass
        import matplotlib.dates as mdates
        a.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 4, 7, 10)))
        a.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
        a.set_xlim(datetime(2025, 1, 1), datetime(2026, 9, 1))
        a.set_ylim(0.4, 2.5)
        a.set_ylabel("C_L / median of the module's dark nights")
        if len(nn):
            a.text(0.99, 0.02, "shaded: nights without darkness (SZA <= 100)", transform=a.transAxes,
                   ha="right", fontsize=7, color="0.4")
        a.set_title(f"{_site(cen[key]['site'])} {cen[key]['lat']:.1f} N ({cen[key]['type']})")
        a.legend(fontsize=7.5, loc="upper left")
        a.grid(alpha=0.3)
    a = fig.add_subplot(gs[:, 2])
    low = ls[(ls.kind == "low-sun")].copy()
    dark = ls[(ls.kind == "dark")].set_index("key")
    low["lat"] = [cen[k]["lat"] if k in cen else 90.0 for k in low.key]
    low = low.sort_values("lat")
    y = np.arange(len(low))
    for yi, (_, r) in zip(y, low.iterrows()):
        if r.key in dark.index:
            d = dark.loc[r.key]
            a.barh(yi, 2 * d.rsd, left=d["median"] - d.rsd, color="0.85", height=0.6)
    a.errorbar(low["median"], y, xerr=[np.nan_to_num(low["median"] - low.ci_lo),
                                       np.nan_to_num(low.ci_hi - low["median"])],
               fmt="o", color="tab:orange", mec="k", ms=7, capsize=3,
               label="low-sun nights: median (95 % CI)")
    a.plot([], [], "s", color="0.85", ms=10, label="dark nights (leave-one-out): median +-1 robust sd")
    for yi, (_, r) in zip(y, low.iterrows()):
        a.annotate(f"n={int(r.n)}", (r["median"], yi), xytext=(6, 5), textcoords="offset points",
                   fontsize=7)
    a.set_yticks(y)
    a.set_yticklabels([("ALL units" if r.key == "ALL" else f"{_site(r.site)} {r.lat:.1f} N")
                       for _, r in low.iterrows()], fontsize=8)
    a.axvline(0, color="0.4", lw=1)
    a.set_xlim(-40, 40)
    a.set_xlabel("ln C - median of the dark nights within +-45 d [%]\n(network seasonal course removed)")
    a.set_title("Continuity at the edges of the no-night season")
    a.legend(fontsize=7.5, loc="lower right")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "p3_continuity.png", dpi=130)
    plt.close(fig)


def fig_qc(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dark = df[(df.run == "op") & df.ok]
    low = df[(df.run == LOWSUN_RUN) & df.ok & (df.flag == LOW_SUN)]
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.4))
    for a, est, lab in ((ax[0], "resid_pct", "night-mean median (8-12 km)"),
                        (ax[1], "resid_tmed", "median over profiles (8-12 km)")):
        bins = np.linspace(-30, 60, 91)
        for g, col, name in ((dark, "0.4", "dark nights (op)"), (low, "tab:orange", "low-sun nights")):
            x = g[est].values
            x = x[np.isfinite(x)]
            a.hist(np.clip(x, bins[0], bins[-1]), bins=bins, density=True, histtype="step", lw=1.8,
                   color=col, label=f"{name}: n={x.size}, median {np.median(x):+.1f} %, "
                                    f"|x|>15 %: {100 * np.mean(np.abs(x) > 15):.0f} %")
        a.set_xlabel(f"far-range residual as C_L bias at the window centre [%] ({lab})")
        a.set_ylabel("density")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    fig.suptitle("Solar-background QC statistic: dark versus low-sun nights")
    fig.tight_layout()
    fig.savefig(OUT / "p3_qc.png", dpi=130)
    plt.close(fig)


def main():
    df = load_runs()
    cen = census()
    kap = kappas()
    say(f"{len(df)} rows, {df.key.nunique()} units, {df.groupby('key').date.nunique().sum()} unit-nights")
    validation(df, cen)
    av, gain = availability(df, cen)
    bp, bs = bias_vs_sza(df, kap, cen)
    threshold_bias(df, kap)
    window_variants(df, kap)
    trigger_study(df, kap)
    cont, series = continuity(df, cen, kap)
    season = network_harmonic()
    lc, ls = local_continuity(df, cen, season, kap)
    residual_qc(df)
    fig_availability(av)
    fig_bias(bp, bs)
    fig_continuity(series, ls, cen)
    fig_qc(df)
    (OUT / "p3_digest.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
