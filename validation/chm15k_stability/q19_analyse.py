# -*- coding: utf-8 -*-
"""q19 - final analysis: offsets, scatter/amplitude, twins, availability, injection, figures.

Inputs (all under C:/DATA/Projects/202610_CHM15k_stability/q19/):
  record/  (q19_gradient_gate.py record: op = operational v2.2, fix = measured noise in the tier-2 gates)
  gate/    (q19_gradient_gate.py gate: exact constants with each gate variant of q19_gate_variants.json)
  inject/  (q19_inject.py), q19_blocks_windows.pkl (q19_offline.py), q19_unit_gradient.csv
Series compared per unit (tier 1 is bit-identical in all of them):
  t1 = tier-1 nights only; v22 = operational eprof_v2.2 (published); fix = measured noise, no new gate;
  <variant> = tier 1 + the tier-2 nights the gate variant keeps (with best-first re-selection).
Offsets of tier-2 constants: ln C - expected tier-1 level at that date (local line through the unit's
tier-1 nights within +-45 d; the plain +-45 d median of q16e is also given).
Writes q19_*.csv, q19_tables.md and landscape figures (altitude on the Y axis).
"""
from __future__ import annotations

import glob
import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q = L.Q
Q19 = Q.Q19
warnings.filterwarnings("ignore")
MD = []
SEASONS = (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)), ("JJA", (6, 7, 8)), ("SON", (9, 10, 11)))
LABEL = {"t1": "tier 1 only", "v22": "eprof_v2.2 (published)", "fix": "v2.2, measured noise",
         "fix_vc": "measured noise + vertical gate (k=2)", "fix_vc_k3": "measured noise + vertical gate (k=3)",
         "op_vc": "operational noise + vertical gate", "fix_blk": "measured noise + cleanest-block only",
         "op_gext": "linear gradient over [zb-1 km, z_hi] (proposal as briefed)",
         "fix_vc2": "FINAL: measured noise + vertical gate (k=2)",
         "fix_vc2_k3": "measured noise + vertical gate (k=3)",
         "op_vc2": "operational noise + vertical gate (k=2)"}
FIRST_DESIGN = ("fix_vc", "fix_vc_k3", "op_vc")   # 'no block above' rejected (superseded, see summary)


def md_table(df, fmt=None):
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                f = (fmt or {}).get(c, "{:.1f}")
                cells.append("–" if not np.isfinite(v) else f.format(v))
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def load_gate():
    fr = []
    for f in sorted(Q.GATE_DIR.glob("*/*.csv")):
        try:
            fr.append(pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str}))
        except Exception:  # noqa: BLE001
            continue
    g = pd.concat(fr, ignore_index=True)
    g["key"] = g.wmo + "_" + g.ident
    g["t"] = pd.to_datetime(g.date)
    g["month"] = g.t.dt.month
    g["ok"] = g.flag.isin(Q.SUCCESS) & (g.c > 0)
    g["lnc"] = np.where(g.ok, np.log(g.c.where(g.c > 0)), np.nan)
    return g


def build_series(rec, gate, variants):
    """Successful nights of each series: key, date, t, month, lnc, tier src, type."""
    op = rec[rec.variant == "op"]
    fx = rec[rec.variant == "fix"]
    t1 = op[(op.tier == 1) & op.ok].assign(src="t1")
    S = {"t1": t1,
         "v22": op[op.ok].assign(src=lambda x: np.where(x.tier == 1, "t1", "t2")),
         "fix": fx[fx.ok].assign(src=lambda x: np.where(x.tier == 1, "t1", "t2"))}
    for name in variants:
        g = gate[(gate.variant == name) & gate.ok & (gate.tier == 2)].assign(src="t2")
        S[name] = pd.concat([t1, g], ignore_index=True)
    keep = ["key", "type", "date", "t", "month", "lnc", "c", "zb", "zt", "src"]
    return {k: v[[c for c in keep if c in v.columns]].copy() for k, v in S.items()}


def tier2_offsets(S, rec, keys, ver=None):
    """Per-night offsets of the tier-2 nights of every series (units in keys): plain +-45 d median
    (q16e), local line through all tier-1 nights (main), and local line through VERIFIED tier-1 nights
    only (ver = {(key, date): bool}; tier-1 nights whose own window passes the vertical test)."""
    op = rec[rec.variant == "op"]
    vmask = None
    if ver is not None:
        vmask = pd.Series([ver.get((k, d), False) for k, d in zip(op.key, op.date)], index=op.index)
    out = []
    for name, s in S.items():
        if name == "t1":
            continue
        t2 = s[(s.src == "t2") & s.key.isin(keys)]
        for key, g in t2.groupby("key"):
            gop = op[op.key == key]
            med, loc, nref = L.unit_offsets(gop, g.t.values, g.lnc.values)
            vloc = np.full(len(g), np.nan)
            if vmask is not None:
                _, vloc, _ = L.unit_offsets(gop, g.t.values, g.lnc.values, ref_mask=vmask)
            out.append(pd.DataFrame(dict(series=name, key=key, date=g.date.values, month=g.month.values,
                                         zb=g.zb.values, off_med=med, off_loc=loc, off_ver=vloc, nref=nref)))
    return pd.concat(out, ignore_index=True)


def tier1_loo(rec, keys):
    """Leave-one-out offsets of tier-1 nights (the night-to-night floor of the offset metric)."""
    op = rec[(rec.variant == "op") & rec.key.isin(keys)]
    out = []
    for key, g in op.groupby("key"):
        t1 = g[(g.tier == 1) & g.ok].sort_values("t")
        for i in range(len(t1)):
            m = np.ones(len(t1), bool)
            m[i] = False
            ref, _ = L.local_ref(t1.t.values[m], t1.lnc.values[m], [t1.t.values[i]])
            out.append(dict(key=key, date=t1.date.values[i], month=t1.month.values[i],
                            off_loc=t1.lnc.values[i] - ref[0]))
    return pd.DataFrame(out)


def stats_block(o):
    o = np.asarray(o, float)
    o = o[np.isfinite(o)]
    if not o.size:
        return dict(n=0, median=np.nan, mean=np.nan, sd=np.nan, p90=np.nan)
    return dict(n=int(o.size), median=100 * np.median(o), mean=100 * np.mean(o),
                sd=100 * L.robust_sd(o), p90=100 * np.percentile(np.abs(o), 90))


def availability(rec, S, keys):
    op = rec[(rec.variant == "op") & rec.key.isin(keys)]
    clear = op[~op.flag.isin([0.0, -1.0, -98.0, -99.0])]
    n_clear = clear.groupby("key").size()
    rows = []
    for name, s in S.items():
        n = s[s.key.isin(keys)].groupby("key").size().reindex(n_clear.index).fillna(0)
        a = 100 * n / n_clear
        rows.append(dict(series=name, nights=int(n.sum()), avail_pooled=100 * n.sum() / n_clear.sum(),
                         avail_unit_median=float(np.median(a)), avail_unit_min=float(np.min(a))))
    return pd.DataFrame(rows), n_clear


def verified_tier1(rec):
    """(key, date) -> True when the night's TIER-1 window passes the main vertical-consistency rule."""
    W = pickle.load(open(Q19 / "q19_blocks_windows.pkl", "rb"))
    G = Q.load_gate_variants()[0][1]["gate"]
    out = {}
    for w in W:
        if w["variant"] == "op" and w["which"] == "sel":
            fn = L.vc_verdict2 if not G.get("require_above", True) else L.vc_verdict
            out[(w["key"], w["date"])] = fn(w, k=G["k"], eps_tol=G["eps_tol"],
                                            tol_lo=G.get("tol_lo", np.inf)) == "accept"
    return out


def twins(S, rec, ver=None):
    """Same-night twin offsets of each series' tier-2 nights (Schiphol A-D, Budapest A/B, Magurele A/B).
    ver: optional {(key, date): bool}; when given only VERIFIED tier-1 nights of the twin serve as the
    reference (on a night that needs tier 2 the twin's own low tier-1 window is often in aerosol too)."""
    op = rec[rec.variant == "op"]
    t1ok = {}
    for k, g in op.groupby("key"):
        g = g[(g.tier == 1) & g.ok]
        if ver is not None:
            g = g[np.array([ver.get((k, d), False) for d in g.date], dtype=bool)]
        t1ok[k] = g[["date", "lnc"]]
    rows = []
    for name, s in S.items():
        if name == "t1":
            continue
        t2 = s[(s.src == "t2") & s.key.str.split("_").str[0].isin(Q.TWINS)]
        ev = t2.assign(lnc_eval=t2.lnc, label=name)[["key", "date", "lnc_eval", "label"]]
        rows.append(L.twin_offsets(None, t1ok, ev))
    T = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    # tier-1 vs tier-1 same-night floor (leave-one-out over the pair reference)
    fl = []
    for wmo, ids in Q.TWINS.items():
        for a in ids:
            for b in ids:
                if a >= b:
                    continue
                ka, kb = f"{wmo}_{a}", f"{wmo}_{b}"
                if ka not in t1ok or kb not in t1ok:
                    continue
                m = t1ok[ka].merge(t1ok[kb], on="date", suffixes=("_x", "_y"))
                if len(m) < 10:
                    continue
                d = (m.lnc_x - m.lnc_y).values
                tt = pd.to_datetime(m.date).values.astype("datetime64[D]").astype(float)
                for i in range(len(m)):
                    w = (np.abs(tt - tt[i]) <= 45)
                    w[i] = False
                    if w.sum() >= 3:
                        fl.append(d[i] - np.median(d[w]))
    return T, np.array(fl)


def unit_series(S, keys):
    rows = []
    for name, s in S.items():
        for key, g in s[s.key.isin(keys)].groupby("key"):
            rows.append(dict(series=name, key=key, **L.series_stats(g.t.values, g.lnc.values)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------
# injection
# ------------------------------------------------------------------------------------------------
def load_inject():
    """Pass 1 (inject/: ungated op/fix records + first-design gates), pass 2 (inject2/: revised gates) and
    pass 3 (inject3/: aerosol and uniform cases re-injected ADDITIVELY). The aerosol/uniform rows of passes
    1-2 are dropped: their multiplicative injection also scaled the night's boundary-layer aerosol without
    the matching extinction (the Klett correction then inflated C) -- see q19_inject.aerosol_additive."""
    fr = [pd.read_csv(f, dtype={"date": str}) for f in glob.glob(str(Q19 / "inject" / "*.csv"))]
    fr += [pd.read_csv(f, dtype={"date": str}) for f in glob.glob(str(Q19 / "inject2" / "*.csv"))]
    fr = [f for f in fr if len(f)]
    if not fr:
        return pd.DataFrame()
    d = pd.concat(fr, ignore_index=True)
    d = d[~d.kind.isin(["aerosol", "uniform"])]
    f3 = [pd.read_csv(f, dtype={"date": str}) for f in glob.glob(str(Q19 / "inject3" / "*.csv"))]
    f3 = [f for f in f3 if len(f)]
    if f3:
        d = pd.concat([d] + f3, ignore_index=True)
    d["key"] = d.wmo + "_" + d.ident.astype(str)
    d["ok"] = d.flag.isin(Q.SUCCESS) & (d.c > 0)
    d["series"] = np.where(d.gate.fillna("") == "", d.variant, d.gate)
    return d


def inject_tables(I):
    out = {}
    # first-choice verdict of the FINAL rule on the operational selection (record rows carry the block
    # statistics of the selected window computed with k=2, eps_tol=0.10, tol_lo=0.10)
    I = I.copy()
    I["verdict_final"] = np.where(I.lcb > 0.10, "aerosol", np.where(I.ucb_lo < -0.10, "tail", "accept"))
    I.loc[~np.isfinite(I.lcb), "verdict_final"] = "unverifiable"
    # (0) forced tier 2 on the clean night (native noise)
    c = I[I.case == "clean"]
    out["clean_forced"] = c.groupby("series").agg(nights=("date", "size"), kept=("ok", "mean"),
                                                  bias_med=("bias", lambda x: 100 * np.nanmedian(x)),
                                                  bias_p90=("bias", lambda x: 100 * np.nanpercentile(np.abs(x), 90)))
    # (i) noise ladder: nights rescued by tier 2 under operational v2.2 at each level
    n = I[I.kind == "noise"].copy()
    trig = n[(n.series == "op") & (n.tier == 2) & n.ok][["key", "date", "amp"]]
    first = trig.groupby(["key", "date"]).amp.min().rename("a_trig").reset_index()
    nn = n.merge(trig.assign(t2=True), on=["key", "date", "amp"], how="inner")
    rows = []
    for (a, sname), g in nn.groupby(["amp", "series"]):
        rows.append(dict(a=a, noise_factor=np.sqrt(1 + a * a), series=sname, n=len(g), kept=100 * g.ok.mean(),
                         bias_med=100 * np.nanmedian(g.bias), bias_sd=100 * L.robust_sd(g.bias.values)))
    out["noise"] = pd.DataFrame(rows)
    allt = []
    for sname, g in nn.groupby("series"):
        allt.append(dict(series=sname, n=len(g), kept=100 * g.ok.mean(),
                         false_rejection=100 * (1 - g.ok.mean()),
                         bias_med=100 * np.nanmedian(g.bias), bias_mean=100 * np.nanmean(g.bias),
                         bias_sd=100 * L.robust_sd(g.bias.values)))
    out["noise_all"] = pd.DataFrame(allt)
    out["a_trig"] = first
    # (ii)-(iv) per case and condition
    rows = []
    for (kind, amp, H, cond, sname), g in I[I.kind.isin(["aerosol", "baseline", "uniform"])].groupby(
            ["kind", "amp", "H", "cond", "series"], dropna=False):
        ref = I[(I.kind == kind) & (I.amp == amp) & ((I.H == H) | (I.H.isna() & pd.isna(H))) & (I.cond == cond)
                & (I.series == "op")]
        rr = ref.set_index(["key", "date"])
        gg = g.set_index(["key", "date"])
        common = gg.index.intersection(rr.index)
        # nights where the OPERATIONAL v2.2 produced a (biased) tier-2 constant
        t2 = rr.loc[common]
        sel = t2[(t2.tier == 2) & t2.ok].index
        gs = gg.loc[sel] if len(sel) else gg.iloc[:0]
        det = np.nan
        if sname == "op" and len(sel):
            det = 100 * np.mean(gs.verdict_final.astype(str) != "accept")
        rows.append(dict(kind=kind, amp=amp, H=H, cond=cond, series=sname, n_t2=len(sel),
                         kept=100 * gs.ok.mean() if len(gs) else np.nan,
                         bias_med=100 * np.nanmedian(gs.bias) if len(gs) else np.nan,
                         bias_p90=100 * np.nanpercentile(np.abs(gs.bias.dropna()), 90) if gs.bias.notna().any() else np.nan,
                         s_win_med=100 * np.nanmedian(t2.loc[sel].s_win) if len(sel) else np.nan,
                         detect_first=det))
    out["cases"] = pd.DataFrame(rows)
    return out


# ------------------------------------------------------------------------------------------------
# figures
# ------------------------------------------------------------------------------------------------
def fig_offsets(P, loo, order, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    a = ax[0]
    data = [100 * loo.off_loc.dropna().values] + [100 * P[P.series == s].off_loc.dropna().values for s in order]
    labs = ["tier 1\n(leave-one-out)"] + [s for s in order]
    bp = a.boxplot(data, tick_labels=labs, whis=(10, 90), showfliers=False, patch_artist=True)
    for b, col in zip(bp["boxes"], ["0.85"] + ["tab:red"] + ["tab:orange"] + ["tab:green"] * (len(order) - 2)):
        b.set_facecolor(col)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("offset vs the unit's tier-1 level (local line, +-45 d) [%]")
    a.set_title("tier-2 constants vs tier-1 (boxes 25-75 %, whiskers 10-90 %)")
    a.tick_params(axis="x", labelsize=7, rotation=20)
    a.grid(alpha=0.3, axis="y")
    a = ax[1]
    for s, col in zip(order, ("tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:brown", "0.4")):
        g = P[P.series == s].groupby("month").off_loc
        m = g.median()
        n = g.size()
        m = m[n >= 10]
        a.plot(m.index, 100 * m.values, "o-", color=col, label=s, ms=4)
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel("month")
    a.set_ylabel("median offset [%] (>= 10 nights)")
    a.set_title("tier-2 offset by month")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    a = ax[2]
    for s, col in zip(order, ("tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:brown", "0.4")):
        n = P[P.series == s].groupby("month").size()
        a.plot(n.index, n.values, "o-", color=col, label=s, ms=4)
    a.set_xlabel("month")
    a.set_ylabel("tier-2 nights (26 CHM15k units, 2025-01 to 2026-10)")
    a.set_title("tier-2 nights kept, by month")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def fig_injection(T, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(23, 5.6))
    cols = {"op": "tab:red", "fix": "tab:orange", "fix_vc2": "tab:green", "fix_vc2_k3": "tab:olive",
            "op_vc2": "tab:blue", "fix_blk": "tab:purple", "op_gext": "tab:brown", "fix_vc": "0.65"}
    lab = {"op": "v2.2 (published)", "fix": "v2.2, measured noise", "fix_vc2": "final gate (k=2)",
           "fix_vc2_k3": "final gate (k=3)", "op_vc2": "gate, operational noise", "op_gext": "linear gradient (as briefed)",
           "fix_vc": "first design (no block above = reject)"}
    a = ax[0]
    N = T["noise"]
    for s_ in ("op", "fix", "fix_vc2", "fix_vc2_k3", "op_vc2", "op_gext", "fix_vc"):
        g = N[N.series == s_]
        if len(g):
            a.plot(g.noise_factor, g.kept, "o-", color=cols[s_], label=lab[s_], ms=5)
    a.set_xlabel("noise factor (total / native)")
    a.set_ylabel("% of the v2.2-rescued clean nights kept")
    a.set_title("(i) white noise on verified-clean nights")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    C = T["cases"]
    a = ax[1]
    for cond, mk in (("forced", "o"), ("natural", "s")):
        for H, ls in ((1000.0, "-"), (2000.0, "--")):
            g = C[(C.kind == "aerosol") & (C.H == H) & (C.cond == cond) & (C.series == "op")].sort_values("amp")
            a.plot(100 * g.amp, g.detect_first, ls, marker=mk, color="k" if cond == "forced" else "tab:gray",
                   label=f"H = {H / 1e3:g} km, {'native noise, tier 1 off' if cond == 'forced' else 'tier-2 noise level'}")
    a.set_ylim(0, 102)
    a.set_xlabel("injected aerosol R-1 at 4 km [%]")
    a.set_ylabel("% of v2.2 tier-2 windows rejected (first candidate)")
    a.set_title("(ii) detection power, decaying aerosol layer")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    a = ax[2]
    for s_ in ("op", "fix_vc2"):
        for cond, mk in (("forced", "o"), ("natural", "s")):
            for H, ls in ((1000.0, "-"), (2000.0, "--")):
                g = C[(C.kind == "aerosol") & (C.H == H) & (C.cond == cond) & (C.series == s_)].sort_values("amp")
                a.plot(100 * g.amp, g.bias_med, ls, marker=mk, color=cols[s_], ms=4,
                       label=f"{lab[s_]}, H={H / 1e3:g} km, {cond}")
    a.plot([0, 30], [0, 30], color="0.6", lw=0.8)
    a.set_xlabel("injected aerosol R-1 at 4 km [%]")
    a.set_ylabel("median bias of the constants kept [%]")
    a.set_title("(ii) bias after the gate (re-selection keeps most nights)")
    a.legend(fontsize=6)
    a.grid(alpha=0.3)
    a = ax[3]
    for kind, mk in (("baseline", "s"), ("uniform", "^")):
        for s_ in ("op", "fix_vc2"):
            g = C[(C.kind == kind) & (C.cond == "forced") & (C.series == s_)].sort_values("amp")
            a.plot(100 * g.amp, g.bias_med, mk + ("-" if s_ == "op" else "--"), color=cols[s_], ms=6,
                   label=f"{kind}: {lab[s_]}")
    a.axhline(0, color="0.6", lw=0.8)
    a.set_xlabel("baseline B / signal(4 km) [%]  |  uniform R-1 [%]")
    a.set_ylabel("median bias of the constants kept [%]")
    a.set_title("(iii) additive baseline, (iv) uniform aerosol: the gate is blind")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=115)
    plt.close(fig)


def fig_profiles(examples, path):
    """Landscape panel of example R(z) profiles with the window, the blocks and the verdict."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import q19_offline as O
    fig, axs = plt.subplots(1, len(examples), figsize=(4.0 * len(examples), 6.0), squeeze=False)
    for a, (title, key, date, s_idx, e_idx) in zip(axs[0], examples):
        P = L.profile_of(key, date)
        if P is None:
            continue
        z, sig, p, sg, _ = P
        pp = Q.prep_profile(sig, p, z, sg)
        b = O.block_arrays(pp, s_idx, e_idx)
        R = sig / p / b["Rw"]
        sR = sg / p / b["Rw"]
        nb = 10
        m = (z > 300) & (z < 10000) & np.isfinite(R)
        zz, RR, ss = z[m] / 1e3, R[m], sR[m]
        k = len(zz) // nb
        zb_ = zz[:k * nb].reshape(k, nb).mean(1)
        Rb = np.nanmean(RR[:k * nb].reshape(k, nb), 1)
        sb = np.sqrt(np.nansum(ss[:k * nb].reshape(k, nb) ** 2, 1)) / nb
        a.errorbar(Rb, zb_, xerr=sb, fmt="o-", ms=2.5, color="k", lw=0.7, elinewidth=0.6,
                   label="R(z)/R(window), 300 m")
        a.axhspan(b["zs"] / 1e3, b["ze"] / 1e3, color="tab:orange", alpha=0.35, label="tier-2 window")
        for exc, sgm, zc, up in zip(b["exc"], b["sig"], b["zc"], b["is_up"]):
            col = "tab:red" if exc - 2 * sgm > 0.10 else ("tab:purple" if (not up and exc + 2 * sgm < -0.10)
                                                          else "tab:green")
            a.errorbar(1 / (1 + exc), zc / 1e3, xerr=sgm / (1 + exc) ** 2, fmt="s", color=col, ms=6,
                       capsize=2)
        a.axvline(1, color="tab:red", lw=0.7)
        a.axvspan(1 / 1.1, 1.0 / 0.9, color="tab:green", alpha=0.07)
        a.set_xlim(0.2, 2.5)
        a.set_ylim(0.3, 10)
        v = L.vc_verdict(b, k=2.0, eps_tol=0.10, tol_lo=0.10)
        a.set_title(f"{title}\n{key} {date}\nverdict: {v}", fontsize=8)
        a.set_xlabel("R(z) / R(window)")
        a.grid(alpha=0.3)
    axs[0][0].set_ylabel("range above the instrument [km]")
    axs[0][0].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    rec = L.load()
    variants = [v[0] for v in Q.load_gate_variants()]
    gate = load_gate()
    S = build_series(rec, gate, variants)
    chm = sorted(rec[(rec.type == "CHM15k") & (rec.key != L.PAYERNE_KEY)].key.unique())
    order = ["v22", "fix"] + variants
    fig_order = ["v22", "fix"] + [v for v in variants if v not in FIRST_DESIGN]
    # ---------------- offsets ----------------
    ver = verified_tier1(rec)
    P = tier2_offsets(S, rec, chm, ver)
    P.to_csv(Q19 / "q19_offsets_nights.csv", index=False)
    loo = tier1_loo(rec, chm)
    loo.to_csv(Q19 / "q19_tier1_loo.csv", index=False)
    av, n_clear = availability(rec, S, chm)
    rows = []
    for name in ["t1"] + order:
        r = dict(series=name, label=LABEL.get(name, name))
        if name == "t1":
            r.update({f"off_{k}": v for k, v in stats_block(loo.off_loc).items()})
        else:
            p = P[P.series == name]
            r.update({f"off_{k}": v for k, v in stats_block(p.off_loc).items()})
            r["off_med_q16e"] = stats_block(p.off_med)["median"]
            r["off_vs_verified_median"] = stats_block(p.off_ver)["median"]
            r["off_vs_verified_mean"] = stats_block(p.off_ver)["mean"]
            r["JJA_median"] = stats_block(p[p.month.isin([6, 7, 8])].off_loc)["median"]
        a = av[av.series == name].iloc[0]
        r.update(nights=a.nights, tier2=int(((S[name].src == "t2") & S[name].key.isin(chm)).sum()),
                 avail_pooled=a.avail_pooled, avail_unit_median=a.avail_unit_median)
        rows.append(r)
    T = pd.DataFrame(rows)
    n_t2_v22 = T.loc[T.series == "v22", "tier2"].values[0]
    T["tier2_kept_pct"] = 100 * T.tier2 / n_t2_v22
    T.to_csv(Q19 / "q19_offsets_summary.csv", index=False)
    MD.append("## Offsets of tier-2 constants (26 CHM15k units, Payerne A excluded)")
    MD.append(md_table(T[["series", "label", "tier2", "tier2_kept_pct", "off_n", "off_median", "off_mean",
                          "off_sd", "off_p90", "JJA_median", "off_med_q16e", "off_vs_verified_median",
                          "off_vs_verified_mean", "avail_pooled", "avail_unit_median"]]))
    # first verdict breakdown of the operational tier-2 nights for the main gate
    main_name = variants[0]
    G = gate[(gate.variant == main_name)]
    op = rec[rec.variant == "op"]
    t2op = op[(op.tier == 2) & op.ok & op.key.isin(chm)]
    gg = G.set_index(["key", "date"])
    br = []
    for (key, d), r in t2op.set_index(["key", "date"]).iterrows():
        if (key, d) not in gg.index:
            continue
        x = gg.loc[(key, d)]
        if isinstance(x, pd.DataFrame):
            x = x.iloc[0]
        br.append(dict(key=key, date=d, first=x.first_verdict, final=("kept" if x.ok else "rejected"),
                       moved=x.moved, lnc_v22=np.log(r.c), lnc_gate=x.lnc))
    BR = pd.DataFrame(br)
    if len(BR):
        BR["t"] = pd.to_datetime(BR.date)
        offs = []
        for key, g in BR.groupby("key"):
            _, loc, _ = L.unit_offsets(op[op.key == key], g.t.values, g.lnc_v22.values)
            offs.append(pd.Series(loc, index=g.index))
        BR["off_v22"] = pd.concat(offs)
        tab = BR.groupby(["first"]).agg(n=("off_v22", "size"), off_v22_median=("off_v22", lambda x: 100 * np.nanmedian(x)),
                                        off_v22_mean=("off_v22", lambda x: 100 * np.nanmean(x)),
                                        kept_after_reselection=("final", lambda x: 100 * np.mean(x == "kept")))
        tab = tab.reset_index()
        tab.to_csv(Q19 / "q19_first_verdicts.csv", index=False)
        MD.append(f"## Operational tier-2 nights by the {main_name} first-candidate verdict (offset of the PUBLISHED constant)")
        MD.append(md_table(tab))
    # published tier-2 constants split by the FINAL outcome under the main gate
    G1 = gate[gate.variant == main_name].drop_duplicates(["key", "date"]).set_index(["key", "date"])
    outc = []
    for key, g in t2op.groupby("key"):
        _, loc, _ = L.unit_offsets(op[op.key == key], g.t.values, g.lnc.values)
        for d, lnc, o in zip(g.date, g.lnc, loc):
            x = G1.loc[(key, d)] if (key, d) in G1.index else None
            if x is None or not x.ok:
                outc.append(dict(outcome="rejected", off_pub=o, off_new=np.nan))
            elif abs(x.lnc - lnc) < 1e-9:
                outc.append(dict(outcome="kept, same window", off_pub=o, off_new=o))
            else:
                outc.append(dict(outcome="kept, window re-selected", off_pub=o, off_new=o + (x.lnc - lnc)))
    OC = pd.DataFrame(outc)
    rows = []
    for name, g in OC.groupby("outcome"):
        a, b = stats_block(g.off_pub), stats_block(g.off_new)
        rows.append(dict(outcome=name, nights=len(g), published_median=a["median"], published_mean=a["mean"],
                         final_median=b["median"], final_mean=b["mean"]))
    OCT = pd.DataFrame(rows)
    OCT.to_csv(Q19 / "q19_outcomes.csv", index=False)
    MD.append(f"## Published tier-2 constants by final {main_name} outcome (offset vs tier-1 local line)")
    MD.append(md_table(OCT, fmt={"nights": "{:.0f}"}))
    # ---------------- per unit ----------------
    U = unit_series(S, chm)
    U.to_csv(Q19 / "q19_unit_series_stats.csv", index=False)
    piv = U.pivot_table(index="key", columns="series", values=["amp", "n2n", "n"])
    net = U.groupby("series")[["amp", "n2n", "n"]].median().reindex(["t1"] + order).reset_index()
    # paired differences vs v22
    dd = []
    for s in ["t1"] + order:
        if s == "v22":
            continue
        d_amp = (piv["amp"][s] - piv["amp"]["v22"]).dropna()
        d_n2n = (piv["n2n"][s] - piv["n2n"]["v22"]).dropna()
        dd.append(dict(series=s, d_amp_median=np.median(d_amp), d_n2n_median=np.median(d_n2n),
                       units=len(d_amp)))
    net = net.merge(pd.DataFrame(dd), on="series", how="left")
    MD.append("## Night-to-night scatter and annual amplitude (network medians over units; paired differences vs v2.2)")
    MD.append(md_table(net, fmt={"n": "{:.0f}", "units": "{:.0f}"}))
    # per-unit offsets for the main gate + unit gradient
    ug = pd.read_csv(Q19 / "q19_unit_gradient.csv") if (Q19 / "q19_unit_gradient.csv").exists() else None
    pu = P.groupby(["series", "key"]).off_loc.agg(["median", "size"]).reset_index()
    pu = pu.pivot_table(index="key", columns="series", values="median")
    pu = (100 * pu).round(1)
    if ug is not None:
        pu = pu.join(ug.set_index("key")[["g_med"]].round(1).rename(columns={"g_med": "unit_gradient_%/km"}))
    pu.to_csv(Q19 / "q19_unit_offsets.csv")
    # ---------------- twins ----------------
    for tag, vv in (("all tier-1 twin nights", None), ("verified tier-1 twin nights", ver)):
        TW, floor = twins(S, rec, vv)
        if not len(TW):
            continue
        # Budapest A instrument event (2026-04-05..04-19: C -40 % on tier-1 nights too, then recovers)
        ev = (TW.key == "0-20000-0-12843_A") & (TW.date >= "20260403") & (TW.date <= "20260422")
        if ev.any():
            tw2 = TW[~ev].groupby("label").offset.apply(lambda x: pd.Series(stats_block(x))).unstack().reset_index()
            MD.append(f"### twins ({tag}), Budapest A event 2026-04-03..22 excluded")
            MD.append(md_table(tw2, fmt={"n": "{:.0f}"}))
        TW.to_csv(Q19 / f"q19_twins_{'verified' if vv else 'all'}.csv", index=False)
        tw = TW.groupby("label").offset.apply(lambda x: pd.Series(stats_block(x))).unstack().reset_index()
        tw = pd.concat([pd.DataFrame([dict(label="tier1/tier1 floor", **stats_block(floor))]), tw])
        MD.append(f"## Co-located twins ({tag}): same-night offset = ln C_X - ln C_Y minus the pair's median "
                  "ratio on nights where both are tier 1 (+-45 d)")
        MD.append(md_table(tw, fmt={"n": "{:.0f}"}))
    # ---------------- availability per instrument type / Payerne ----------------
    other = []
    for typ, keys in (("Payerne A (CHM15k, abnormal)", [L.PAYERNE_KEY]),
                      ("CL61", sorted(rec[rec.type == "CL61"].key.unique())),
                      ("CL31", sorted(rec[rec.type == "CL31"].key.unique()))):
        Po = tier2_offsets(S, rec, keys)
        a2, _ = availability(rec, S, keys)
        for name in ["t1"] + order:
            r = dict(group=typ, series=name)
            if name != "t1":
                r.update({f"off_{k}": v for k, v in stats_block(Po[Po.series == name].off_loc).items()})
            r.update(nights=int(a2[a2.series == name].nights.values[0]),
                     avail=float(a2[a2.series == name].avail_pooled.values[0]))
            other.append(r)
    OT = pd.DataFrame(other)
    OT.to_csv(Q19 / "q19_other_instruments.csv", index=False)
    MD.append("## Payerne A, CL61, CL31")
    MD.append(md_table(OT, fmt={"off_n": "{:.0f}", "nights": "{:.0f}"}))
    # ---------------- injection ----------------
    I = load_inject()
    setup = Q19 / "inject2" / "inject_setup.json"
    n_exp = 2 * len(json.load(open(setup))["nights"]) if setup.exists() else -1
    n_exp = 3 * n_exp // 2 if n_exp > 0 else n_exp
    n_got = (len(glob.glob(str(Q19 / "inject" / "*.csv"))) + len(glob.glob(str(Q19 / "inject2" / "*.csv")))
             + len(glob.glob(str(Q19 / "inject3" / "*.csv"))))
    if n_got < n_exp:
        print(f"injection incomplete ({n_got}/{n_exp} nights): skipped")
        I = I.iloc[:0]
    if len(I):
        IT = inject_tables(I)
        for k, v in IT.items():
            if isinstance(v, pd.DataFrame):
                v.to_csv(Q19 / f"q19_inject_{k}.csv", index=k == "clean_forced")
        MD.append("## Injection (i): white noise on clean nights, nights rescued by tier 2 under v2.2")
        MD.append(md_table(IT["noise_all"], fmt={"n": "{:.0f}"}))
        MD.append(md_table(IT["noise"], fmt={"n": "{:.0f}", "a": "{:.0f}", "noise_factor": "{:.2f}"}))
        MD.append("## Injection (0): forced tier 2 on clean nights, native noise")
        MD.append(md_table(IT["clean_forced"].reset_index(), fmt={"nights": "{:.0f}", "kept": "{:.2f}"}))
        MD.append("## Injection (ii)-(iv): aerosol / baseline / uniform")
        MD.append(md_table(IT["cases"], fmt={"n_t2": "{:.0f}", "amp": "{:.2f}", "H": "{:.0f}"}))
        fig_injection(IT, Q19 / "q19_injection.png")
    fig_offsets(P, loo, fig_order, Q19 / "q19_offsets.png")
    (Q19 / "q19_tables.md").write_text("\n\n".join(MD), encoding="utf-8")
    print("\n\n".join(MD))


if __name__ == "__main__":
    main()
