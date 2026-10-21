# -*- coding: utf-8 -*-
"""ref23 step 6: write the report ref23_agreement.md from the analysis outputs.

Every number in the report is read from the CSV outputs of ref23_analysis.py / ref23_cloud.py
(no number is typed by hand), so a re-run of the chain regenerates a consistent report.

Usage: python validation/v23_reference/ref23_report.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402
import ref23_stats as S  # noqa: E402
from ref23_analysis import pair_label, SITE_LABEL, TYPE  # noqa: E402

OUT = RC.OUT
FIG = {k: (OUT / f"ref23_{k}.png").as_posix() for k in
       ("timeseries_cl61", "timeseries_twins", "pulls", "reliability", "profiles", "lambda",
        "cloud_uncertainty", "budget")}


def f1(x, sgn=False, nd=1):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{nd}f}" if sgn else f"{x:.{nd}f}"


def load():
    d = {}
    for name in ("summary", "summary_tier", "summary_season", "summary_pooled", "summary_cl31_era",
                 "wv_transmission", "reproduction", "reproduction_nights",
                 "cloud_uncertainty", "cloud_vs_published", "budget_units", "lambda_alpha",
                 "reliability", "fate_v22_to_v23", "pairs"):
        p = OUT / f"{name}.csv"
        d[name] = pd.read_csv(p, dtype={"date": str}, low_memory=False) if p.exists() else pd.DataFrame()
    # every Rayleigh night (both versions): night counts per site and budget components per tier
    rp = OUT / "tables" / "rays.pkl"
    cols = ["site", "ident", "date", "variant", "c", "tier", "zb", "zt", "g", "sig_ens", "sig_vert",
            "sig_noise", "sig_aer", "sig_mol"]
    d["rays"] = pd.read_pickle(rp)[cols] if rp.exists() else pd.DataFrame(columns=cols)
    cp = OUT / "tables" / "clouds.pkl"
    d["clouds"] = (pd.read_pickle(cp)[["site", "ident", "date", "c"]] if cp.exists()
                   else pd.DataFrame(columns=["site", "ident", "date", "c"]))
    return d


def cl31_step(clouds, days=30):
    """Ratio of the Payerne CL31 cloud constant after / before the 2026-07-07 block swap (medians of
    the successful days within `days` days on either side, the swap day excluded)."""
    g = clouds[(clouds.site == "payerne") & (clouds.ident == "B") & np.isfinite(clouds.c.astype(float))]
    t = pd.to_datetime(g.date, format="%Y%m%d")
    swap = pd.Timestamp("2026-07-07")
    before = g.c[(t < swap) & (t >= swap - pd.Timedelta(days=days))]
    after = g.c[(t > swap) & (t <= swap + pd.Timedelta(days=days))]
    if len(before) < 3 or len(after) < 3:
        return np.nan, 0, 0
    return float(np.median(after) / np.median(before)), len(before), len(after)


def tier_budget(rays, site, ident):
    """Median window base [km], gradient [%/km] and budget components [%] of the v2.3 nights of one
    unit, per tier ({1: {...}, 2: {...}})."""
    g = rays[(rays.site == site) & (rays.ident == ident) & (rays.variant == "v23") & rays.c.notna()]
    out = {}
    for t, gg in g.groupby("tier"):
        out[int(t)] = dict(n=len(gg), zb=np.nanmedian(gg.zb) / 1e3, g=100 * np.nanmedian(gg.g),
                           **{k: 100 * np.nanmedian(gg[k]) for k in
                              ("sig_ens", "sig_vert", "sig_noise", "sig_aer", "sig_mol")})
    return out


def get(st, site, pair, version, variant="main"):
    m = st[(st.site == site) & (st.pair == pair) & (st.version == version) & (st.variant == variant)]
    if not len(m) or m.iloc[0].n < 8:
        return None
    return m.iloc[0].to_dict()


def pooled(pp, group, version, variant="main"):
    m = pp[(pp.group == group) & (pp.version == version) & (pp.variant == variant)]
    return m.iloc[0].to_dict() if len(m) and m.iloc[0].n >= 8 else None


def summary_table(st):
    head = ("| site | pair X / R | n nights (v2.3) | bias: median ln ρ [%] (95 % CI) | "
            "observed σ [%] | predicted σ [%] v2.2 → v2.3 | centred pull v2.2 → v2.3 "
            "(pre-review budget) | uncentred pull v2.2 → v2.3 | % within 2 σ v2.2 → v2.3 | "
            "Spearman v2.3 | verdict (v2.3) |")
    lines = [head, "|" + "---|" * 11]
    for (site, pair), g in st[st.variant == "main"].groupby(["site", "pair"], sort=False):
        kind = g.kind.iloc[0]
        if kind == "CC":
            s = get(st, site, pair, "cloud")
            s22 = s
            spre = None
        else:
            s = get(st, site, pair, "v23")
            s22 = get(st, site, pair, "v22pub")
            spre = get(st, site, pair, "v23pre")
        if s is None:
            continue
        pre = f" ({f1(spre['pull_c'], nd=2)})" if spre else ""
        lines.append(
            f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {s['n']:.0f} | "
            f"{f1(s['bias'], True)} ({f1(s['ci_lo'], True)}, {f1(s['ci_hi'], True)}) | "
            f"{f1(s['scatter'])} | {f1(s22['sig_pred']) if s22 else '–'} → {f1(s['sig_pred'])} | "
            f"{f1(s22['pull_c'], nd=2) if s22 else '–'} → {f1(s['pull_c'], nd=2)}{pre} | "
            f"{f1(s22['pull_u'], nd=2) if s22 else '–'} → {f1(s['pull_u'], nd=2)} | "
            f"{f1(s22['in2_u'], nd=0) if s22 else '–'} → {f1(s['in2_u'], nd=0)} | "
            f"{f1(s['spearman'], True, 2)} | {S.verdict(s)} |")
    return "\n".join(lines)


def sensitivity_table(st, pairs):
    head = ("| site | pair | bias 500–3000 m | bias 1500–3000 m | centred pull: main | with mol | "
            "no conv | 1500–3000 m | 60-min SNR | cloud sigma empirical | local cloud reference |")
    lines = [head, "|" + "---|" * 11]
    for site, pair in pairs:
        kind_v = "cloud" if pair.startswith("B:cloud") else "v23"
        rows = {v: get(st, site, pair, kind_v, v) for v in
                ("main", "with_mol", "no_conv", "hi_band", "w60", "cloud_empirical", "cloud_local")}
        m = rows["main"]
        if m is None:
            continue
        hb = rows["hi_band"]

        def pc(k):
            r = rows[k]
            return f1(r["pull_c"], nd=2) if r else "–"
        lines.append(f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {f1(m['bias'], True)} | "
                     f"{f1(hb['bias'], True) if hb else '–'} | {pc('main')} | {pc('with_mol')} | "
                     f"{pc('no_conv')} | {pc('hi_band')} | {pc('w60')} | {pc('cloud_empirical')} | "
                     f"{pc('cloud_local')} |")
    return "\n".join(lines)


def tier_table(tt):
    lines = ["| site | pair | tier of X | n | bias [%] | observed sigma [%] | predicted sigma [%] | "
             "centred pull | % within 2 sigma |", "|" + "---|" * 9]
    for (site, pair), g in tt[tt.version == "v23"].groupby(["site", "pair"], sort=False):
        for r in g.itertuples():
            if r.split not in ("strict", "fallback"):
                continue
            lines.append(f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {r.split} | {r.n:.0f} | "
                         f"{f1(r.bias, True)} | {f1(r.scatter)} | {f1(r.sig_pred)} | "
                         f"{f1(r.pull_c, nd=2)} | {f1(r.in2_u, nd=0)} |")
    return "\n".join(lines)


def season_table(ss, pairs):
    """v2.3 bias / centred pull per season for the given pairs (n >= 6)."""
    lines = ["| site | pair | DJF | MAM | JJA | SON |", "|---|---|---|---|---|---|"]
    s3 = ss[ss.version == "v23"]
    for site, pair in pairs:
        g = s3[(s3.site == site) & (s3.pair == pair)].set_index("split")
        if not len(g):
            continue
        cells = []
        for sea in ("DJF", "MAM", "JJA", "SON"):
            if sea in g.index and g.loc[sea].n >= 6:
                r = g.loc[sea]
                cells.append(f"n {r.n:.0f}: {f1(r.bias, True)} %, pull {f1(r.pull_c, nd=2)}")
            else:
                cells.append("–")
        lines.append(f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def pooled_table(pp, variant="main"):
    lines = ["| group | version | n | observed sigma (pair-centred) [%] | predicted sigma [%] | "
             "centred pull | % within 2 sigma (centred) | uncentred pull | % within 2 sigma | "
             "Spearman (p) |", "|" + "---|" * 10]
    for grp, g in pp[pp.variant == variant].groupby("group", sort=False):
        for r in g.itertuples():
            if not np.isfinite(r.n) or r.n < 8:
                continue
            lines.append(f"| {grp} | {S.VERSION_LABEL.get(r.version, r.version)} | {r.n:.0f} | "
                         f"{f1(r.scatter)} | {f1(r.sig_pred)} | {f1(r.pull_c, nd=2)} | "
                         f"{f1(r.in2_c, nd=0)} | {f1(r.pull_u, nd=2)} | {f1(r.in2_u, nd=0)} | "
                         f"{f1(r.spearman, True, 2)} ({r.spearman_p:.2g}) |")
    return "\n".join(lines)


def fate_md(fate):
    lines = ["| site | pair | v2.2 night → v2.3 | n | median ln rho with the v2.2 constant [%] | "
             "robust sigma [%] | median ln rho with the v2.3 constant [%] |", "|" + "---|" * 7]
    order = ["strict, kept", "strict, rejected by u/C > 60 %", "fallback, same window",
             "fallback, window moved", "fallback, rejected"]
    for (site, pair), g in fate.groupby(["site", "pair"], sort=False):
        g = g.set_index("fate")
        for c in order:
            if c not in g.index:
                continue
            r = g.loc[c]
            lines.append(f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {c} | {r.n:.0f} | "
                         f"{f1(r.med_v22, True)} | {f1(r.scat_v22)} | {f1(r.med_v23, True)} |")
    return "\n".join(lines)


def repro_md(rep):
    lines = ["| unit | nights rerun | calibrated (v2.2 rerun) | published | common | constant "
             "bit-identical | uncertainty identical |", "|" + "---|" * 7]
    for r in rep.itertuples():
        lines.append(f"| {SITE_LABEL[r.site]} {r.type} {r.ident} | {r.nights_rerun} | {r.ok_rerun} | "
                     f"{r.published} | {r.both} | {f1(r.identical_pct, nd=1)} % | "
                     f"{f1(r.u_identical_pct, nd=1)} % |")
    return "\n".join(lines)


def cloud_md(ct, cmp):
    lines = ["| channel | days | profiles / day (median) | scenes / day | reported 1 sigma "
             "(std/median) [%] | standard error of the median [%] | day-to-day sigma_w (fit) [%] | "
             "drift over 7 d [%] | consecutive-day pull with the reported sigma |",
             "|" + "---|" * 9]
    for r in ct.itertuples():
        lines.append(f"| {r.label} | {r.days} | {r.n_prof_med:.0f} | {r.n_scenes_med:.0f} | "
                     f"{f1(r.sig_rep_med)} | {f1(r.sig_sem_med)} | {f1(r.sig_w_fit)} | "
                     f"{f1(r.drift_7d)} | {f1(r.pull_rep_consec, nd=2)} |")
    return "\n".join(lines)


def main():
    d = load()
    st, pp = d["summary"], d["summary_pooled"]
    branch = subprocess.run(["git", "-C", str(RC.REPO), "rev-parse", "--abbrev-ref", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    text = TEMPLATE(d, st, pp, branch)
    (OUT / "ref23_agreement.md").write_text(text, encoding="utf-8")
    print(f"wrote {OUT / 'ref23_agreement.md'} ({len(text)} chars)")


def counts(d):
    """Night counts per unit and period, from the pair table and the reproduction table."""
    rep = d["reproduction"]
    out = {}
    for r in rep.itertuples():
        out[(r.site, r.ident)] = dict(rerun=r.nights_rerun, v22=r.ok_rerun, pub=r.published,
                                      both=r.both, ident=r.identical_pct)
    return out


def bridging(pt):
    """Share of same-night vs bridged cloud references per cloud-referenced pair (v2.3 nights)."""
    rows = []
    sel = pt[pt.kind.isin(["RC", "self"]) & np.isfinite(pt.get("lnr_v23", np.nan))]
    for (site, pair), g in sel.groupby(["site", "pair"], sort=False):
        k = g.cr_kind.value_counts()
        rows.append(dict(site=site, pair=pair, n=len(g), same=int(k.get("same", 0)),
                         bridge=int(k.get("bridge", 0)),
                         lag_med=float(np.nanmedian(g.loc[g.cr_kind == "bridge", "cr_lag"]))
                         if (g.cr_kind == "bridge").any() else np.nan,
                         sig_med=100 * float(np.nanmedian(g.cr_sig))))
    return pd.DataFrame(rows)


def _s(st, site, pair, version="v23", variant="main"):
    return get(st, site, pair, version, variant) or {}


def _b(s, key="bias"):
    """'+x.x % (CI a, b)' from a stats dict."""
    if not s:
        return "–"
    return f"{f1(s[key], True)} % ({f1(s['ci_lo'], True)}, {f1(s['ci_hi'], True)})"


def bias_table(st, pairs):
    """Bias of each pair in the two bands, in clean air and its Angstrom sensitivity."""
    lines = ["| site | pair | n | 500–3000 m | 1500–3000 m | clean air 500–3000 m | "
             "clean air 1500–3000 m | median sigma_conv: 500–3000 m / clean air |",
             "|" + "---|" * 8]
    for site, pair in pairs:
        v = "cloud" if pair.startswith("B:cloud") else "v23"
        m = get(st, site, pair, v, "main")
        if m is None:
            continue
        hb, cl, ch = (get(st, site, pair, v, k) for k in ("hi_band", "clean", "clean_hi"))
        lines.append(f"| {SITE_LABEL[site]} | {pair_label(site, pair)} | {m['n']:.0f} | {_b(m)} | "
                     f"{_b(hb)} | {_b(cl)} | {_b(ch)} | {f1(m.get('sig_conv'))} / "
                     f"{f1(cl.get('sig_conv')) if cl else '–'} % |")
    return "\n".join(lines)


MAIN_PAIRS = [("payerne", "A:ray/C:cloud"), ("payerne", "A:ray/C:ray"),
              ("payerne", "C:ray/C:cloud"), ("payerne", "B:cloud/C:cloud"),
              ("lindenberg", "0:ray/C:cloud"), ("lindenberg", "0:ray/C:ray"),
              ("lindenberg", "C:ray/C:cloud"), ("aosta", "0:ray/B:cloud"), ("aosta", "0:ray/B:ray"),
              ("aosta", "B:ray/B:cloud")]
# CHM15k (1064 nm) vs CL61 (910 nm): the pairs that need the wavelength conversion
MIXED_PAIRS = [("payerne", "A:ray/C:cloud"), ("payerne", "A:ray/C:ray"),
               ("lindenberg", "0:ray/C:cloud"), ("lindenberg", "0:ray/C:ray"),
               ("aosta", "0:ray/B:cloud"), ("aosta", "0:ray/B:ray")]


def pct(x, nd=0):
    return "–" if x is None or not np.isfinite(x) else f"{x:.{nd}f} %"


def conv_summary(d, st):
    """Size of the 910 -> 1064 nm conversion term over the CHM15k / CL61 pairs: the change of the
    median ln rho per unit of Angstrom exponent (alpha 0.5 -> 1.5), the per-night sigma_conv in the
    500-3000 m band and in clean air, and the total predicted sigma for scale."""
    out = dict(slope=(np.nan, np.nan), conv=(np.nan, np.nan), conv_clean=(np.nan, np.nan),
               pred=(np.nan, np.nan))
    la = d["lambda_alpha"]
    if len(la):
        piv = la.pivot_table(index=["site", "pair"], columns="alpha", values="med")
        # CHM15k / CL61 only (the Payerne CHM15k / CL31 pair is left out)
        piv = piv[[not (s_ == "payerne" and p.endswith("B:cloud")) for s_, p in piv.index]]
        if 1.5 in piv and 0.5 in piv:
            dd = (piv[1.5] - piv[0.5]).abs()
            out["slope"] = (dd.min(), dd.max())

    def rng(vals):
        v = [x for x in vals if x is not None and np.isfinite(x)]
        return (min(v), max(v)) if v else (np.nan, np.nan)
    out["conv"] = rng([_s(st, *k).get("sig_conv") for k in MIXED_PAIRS])
    out["conv_clean"] = rng([_s(st, *k, variant="clean").get("sig_conv") for k in MIXED_PAIRS])
    out["pred"] = rng([_s(st, *k).get("sig_pred") for k in MIXED_PAIRS])
    return out


def rng_txt(r, nd=1):
    return f"{f1(r[0], nd=nd)}–{f1(r[1], nd=nd)}"


def fate_overview(fate, n_min=3):
    """v2.2 fallback nights that v2.3 no longer calibrates vs the strict nights, per pair: offset of
    their median ln rho from the strict median and their robust spreads (v2.2 constants)."""
    rows = []
    for (site, pair), g in fate.groupby(["site", "pair"], sort=False):
        g = g.set_index("fate")
        if "fallback, rejected" in g.index and "strict, kept" in g.index:
            rj, sk = g.loc["fallback, rejected"], g.loc["strict, kept"]
            if rj.n >= n_min:
                rows.append(dict(site=site, pair=pair, n=rj.n, off=abs(rj.med_v22 - sk.med_v22),
                                 s_rej=rj.scat_v22, s_str=sk.scat_v22))
    return pd.DataFrame(rows)


def build_text(d, st, pp, branch):
    """The report.  Numbers come from the analysis CSVs; the wording around them was written
    after reading the results of the run of 2026-10-21 and must be re-read if the data change."""
    pt = d["pairs"]
    rep = d["reproduction"]
    ct, cmp = d["cloud_uncertainty"], d["cloud_vs_published"]
    fate = d["fate_v22_to_v23"]
    tt = d["summary_tier"]
    br = bridging(pt)
    P = {k: _s(st, *k) for k in MAIN_PAIRS}
    P22 = {k: _s(st, *k, version="v22pub") for k in MAIN_PAIRS}
    pc = {g: pooled(pp, g, "v23") for g in pp.group.unique()}
    pc22 = {g: pooled(pp, g, "v22pub") for g in pp.group.unique()}
    pcb = {g: pooled(pp, g, "v22b") for g in pp.group.unique()}
    twin_pairs = [(s_, p_) for s_, p_ in st[st.kind == "RR_twin"][["site", "pair"]]
                  .drop_duplicates().itertuples(index=False)]
    T = {k: _s(st, *k) for k in twin_pairs}
    T22 = {k: _s(st, *k, version="v22pub") for k in twin_pairs}
    L = []
    w = L.append
    # ------------------------------------------------------------------------------- header
    w("# Co-located instruments at the E-PROFILE reference sites: agreement and the v2.3 "
      "calibration uncertainty\n")
    w(f"*Study ref23, 2026-10-21. Branch `{branch}` (on top of `v2.3-dev`, with the reviewed v2.3 "
      "uncertainty budget merged). Scripts "
      "`validation/v23_reference/ref23_*.py`; every number below is written by "
      "`ref23_report.py` from the CSV outputs in `C:/DATA/Projects/202610_v23/reference_agreement/`.*\n")
    w(SUMMARY(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, br, rep))
    w(METHODS(d, st, rep, ct, cmp, br))
    w(RESULTS(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, tt))
    w(CAVEATS(d, st, P, rep, br))
    w(REPRO())
    return "\n".join(L)


def PULL_SENTENCE(g1, g2, g3, g4):
    """One sentence on the v2.3 pooled centred pulls (written from the values, not by hand)."""
    def cls(x):
        if not np.isfinite(x):
            return "not evaluable"
        if x < S.PULL_LO:
            return "over-covered"
        if x > S.PULL_HI:
            return "under-covered"
        return "calibrated"
    parts = [f"CHM15k vs CL61 cloud {cls(g1.get('pull_c', np.nan))}",
             f"CHM15k vs CL61 Rayleigh {cls(g2.get('pull_c', np.nan))}",
             f"CL61 Rayleigh vs cloud {cls(g3.get('pull_c', np.nan))}",
             f"CHM15k twins {cls(g4.get('pull_c', np.nan))}"]
    return ("The v2.2 ensemble under-stated every comparison; with v2.3 (by the [0.75, 1.33] rule): "
            + ", ".join(parts) + ". The cloud calibration is the cleanest test: an independent method "
            "that shares no error term with the Rayleigh budget. The over-coverage of the Rayleigh / "
            "Rayleigh pairs (twins most) is expected: the ratio of two calibrated profiles cancels the "
            "atmospheric part of their calibration errors (aerosol and molecular model of similar "
            "windows the same night), which a single-unit budget must cover (§4.4).")


def SPEARMAN_SENTENCE(gs):
    """How well the per-night sigma ranks the per-night disagreement (Spearman, four groups)."""
    v = [g.get("spearman", np.nan) for g in gs]
    vals = ", ".join(f1(x, True, 2) for x in v)
    fin = [x for x in v if np.isfinite(x)]
    if not fin:
        return ""
    word = "weakly" if max(fin) < 0.3 else "weakly to moderately"
    nil = "; nil for the CHM15k against the cloud reference" if abs(v[0]) < 0.1 else ""
    return (f" Night by night, a larger stated σ flags a larger disagreement only {word} (Spearman "
            f"of |ln ρ − median| with σ_pred {vals}, same order{nil}).")


# The prose sections are filled in below (SUMMARY, METHODS, RESULTS, CAVEATS, REPRO).
def SUMMARY(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, br, rep):
    pt = d["pairs"]
    L = ["## Summary\n"]
    w = L.append
    n_sn = len(pd.unique(pt.site + pt.date))
    ry = d["rays"]
    n_read = len(ry[["site", "date"]].drop_duplicates()) if len(ry) else 0
    last = ry.groupby("site").date.max() if len(ry) else pd.Series(dtype=str)
    ends = []
    for day in sorted(last.unique()):
        names = ", ".join(sorted(SITE_LABEL[s_] for s_ in last[last == day].index))
        ends.append(f"{day[:4]}-{day[4:6]}-{day[6:]} at {names}")
    nrep = int(rep.both.sum())
    nid = int(np.round(np.sum(rep.both * rep.identical_pct / 100.0)))
    w(f"* **What was done.** Every night from 2025-01-01 to the last L1 ({'; '.join(ends)}) of six "
      f"E-PROFILE reference sites — {n_read} site-nights, every unit's L1 read once — was "
      "recalibrated in memory: v2.2 (the published constant — reproduced bit-identically on "
      f"{nid} of {nrep} common nights), and the release candidate v2.3 (`eprof_v2.3` + the v2.3 "
      "budget as revised by the independent review, every component stored); the cloud "
      "calibrations of the CL61 and CL31 reproduce the published ones. The observable is ln ρ, the "
      "night median of the ratio of the two units' calibrated attenuated backscatter over their "
      "common clear, detected pixels of the calibration night (v3 dashboard conventions), each "
      f"unit divided by its own constant of that night; {n_sn} site-nights have at least one pair "
      "value.\n")
    s, s22 = P[("payerne", "A:ray/C:cloud")], P22[("payerne", "A:ray/C:cloud")]
    sr = P[("payerne", "A:ray/C:ray")]
    ss = P[("payerne", "C:ray/C:cloud")]
    tt = d["summary_tier"]

    def tier(site, pair, split, version="v23"):
        m = tt[(tt.site == site) & (tt.pair == pair) & (tt.split == split) & (tt.version == version)]
        return m.iloc[0].to_dict() if len(m) else {}
    ss1, ss2 = tier("payerne", "C:ray/C:cloud", "strict"), tier("payerne", "C:ray/C:cloud", "fallback")
    sa1 = tier("payerne", "A:ray/C:cloud", "strict")
    w(f"* **Payerne, CL61 reference.** CHM15k A / CL61 cloud: {s.get('n', 0):.0f} v2.3 nights, "
      f"bias {_b(s)}, observed σ {f1(s.get('scatter'))} % vs predicted {f1(s.get('sig_pred'))} %, "
      f"{f1(s.get('in2_u'), nd=0)} % of the nights within 2 σ (v2.2: {f1(s22.get('in2_u'), nd=0)} % "
      f"on its {s22.get('n', 0):.0f} nights, bias {f1(s22.get('bias'), True)} %), centred pull "
      f"{_pl(s, s22)}: **{_within(s)} the v2.3 uncertainty**"
      f"{'' if s.get('pull_c', 9) <= S.PULL_HI else ' — the scatter is under-predicted'}; this "
      f"CHM15k is abnormal (noisy receiver; strict-tier nights alone: {sa1.get('n', 0):.0f} nights, "
      f"bias {f1(sa1.get('bias'), True)} %, centred pull {f1(sa1.get('pull_c'), nd=2)}). CHM15k A / "
      f"CL61 Rayleigh: {sr.get('n', 0):.0f} nights, bias {_b(sr)}, centred pull "
      f"{f1(sr.get('pull_c'), nd=2)}, {f1(sr.get('in2_u'), nd=0)} % within 2 σ. The CL61's own two "
      f"calibrations (Rayleigh vs cloud, same photons): bias {_b(ss)}; on strict-tier nights "
      f"{f1(ss1.get('bias'), True)} % with centred pull {f1(ss1.get('pull_c'), nd=2)} "
      f"({f1(ss1.get('in2_u'), nd=0)} % within 2 σ), on its fallback-tier nights — which v2.3 keeps, "
      f"the vertical-consistency test not being applied to the CL61 — {f1(ss2.get('bias'), True)} %.\n")
    # biases (v2.3): Rayleigh-Rayleigh CHM15k / CL61 in the two bands and in clean air
    def bb(site, pair, var):
        r = _s(st, site, pair, "v23", var)
        return f"{f1(r.get('bias'), True)} %" if r else "–"
    w("* **Bias.** In the v3 band (500–3000 m) the CHM15k reads lower than the CL61 Rayleigh at "
      f"Lindenberg ({bb('lindenberg', '0:ray/C:ray', 'main')}) and Aosta "
      f"({bb('aosta', '0:ray/B:ray', 'main')}), but most of it is near-range shape (overlap) and "
      "wavelength conversion: above 1500 m in clean air the two Rayleigh calibrations agree to "
      f"{bb('lindenberg', '0:ray/C:ray', 'clean_hi')} and {bb('aosta', '0:ray/B:ray', 'clean_hi')}. "
      "Against the CL61 cloud calibration: Lindenberg "
      f"{bb('lindenberg', '0:ray/C:cloud', 'main')} (clean air above 1500 m "
      f"{bb('lindenberg', '0:ray/C:cloud', 'clean_hi')}), Aosta {bb('aosta', '0:ray/B:cloud', 'main')} "
      f"({bb('aosta', '0:ray/B:cloud', 'clean_hi')}); the CL61 Rayleigh / cloud ratio is "
      f"{bb('lindenberg', 'C:ray/C:cloud', 'main')} (Lindenberg) and "
      f"{bb('aosta', 'B:ray/B:cloud', 'main')} (Aosta), the size of the cloud method's own "
      "S = 18.8 ± 0.8 sr floor (±4.3 %). Twins carry static, instrumental biases (Budapest B / A "
      f"{bb('budapest', 'B:ray/A:ray', 'main')}, 1500–3000 m {bb('budapest', 'B:ray/A:ray', 'hi_band')}; "
      f"Schiphol B, overlap deficit, {bb('schiphol', 'B:ray/A:ray', 'main')} → "
      f"{bb('schiphol', 'B:ray/A:ray', 'hi_band')} above 1500 m) that no per-night uncertainty "
      "covers. The Payerne CL31 / CL61 rows are dominated by the CL31's electronic dark before its "
      "July 2026 block swap (§4.2).\n")
    keys = ("CHM15k Rayleigh vs CL61 cloud", "CHM15k Rayleigh vs CL61 Rayleigh",
            "CL61 Rayleigh vs CL61 cloud", "CHM15k twins")
    g1, g2, g3, g4 = (pc.get(k) or {} for k in keys)
    h1, h2, h3, h4 = (pc22.get(k) or {} for k in keys)
    q1, q2, q3, q4 = (pooled(pp, k, "v23pre") or {} for k in keys)
    w("* **Is the uncertainty now calibrated?** Centred pull (robust std, 1 = calibrated) v2.2 → "
      f"v2.3 (pre-review budget in brackets): CHM15k vs CL61 cloud {_pl(g1, h1)} "
      f"({f1(q1.get('pull_c'), nd=2)}), CHM15k vs CL61 Rayleigh {_pl(g2, h2)} "
      f"({f1(q2.get('pull_c'), nd=2)}), CL61 Rayleigh vs cloud {_pl(g3, h3)} "
      f"({f1(q3.get('pull_c'), nd=2)}), CHM15k twins {_pl(g4, h4)} ({f1(q4.get('pull_c'), nd=2)}); "
      f"nights within 2 σ (uncentred) {f1(h1.get('in2_u'), nd=0)} → {f1(g1.get('in2_u'), nd=0)} %, "
      f"{f1(h2.get('in2_u'), nd=0)} → {f1(g2.get('in2_u'), nd=0)} %, {f1(h3.get('in2_u'), nd=0)} → "
      f"{f1(g3.get('in2_u'), nd=0)} %, {f1(h4.get('in2_u'), nd=0)} → {f1(g4.get('in2_u'), nd=0)} %. "
      + PULL_SENTENCE(g1, g2, g3, g4) + SPEARMAN_SENTENCE([g1, g2, g3, g4]) + "\n")
    fz = fate.set_index(["site", "pair", "fate"])

    def fv(site, pair, cls_, col="med_v22"):
        k = (site, pair, cls_)
        return (f"{f1(fz.loc[k, col], True)} % (n = {fz.loc[k, 'n']:.0f})" if k in fz.index
                else "–")
    fo = fate_overview(fate)
    good = fo[fo.n >= 5].sort_values("off").iloc[0] if (fo.n >= 5).any() else None
    w("* **The repaired fallback tier, judged by an independent instrument.** The v2.2 fallback-tier "
      "nights that v2.3 no longer calibrates were mostly bad nights: over the "
      f"{len(fo)} pairs with ≥ 3 of them, their median ln ρ (v2.2 constant) sits a median "
      f"{f1(fo.off.median())} % (range {f1(fo.off.min())}–{f1(fo.off.max())} %) from that of the "
      f"strict nights and they scatter more (robust σ, median over the pairs: "
      f"{f1(fo.s_rej.median())} % vs {f1(fo.s_str.median())} %). Rejected vs strict: Payerne "
      f"CHM15k A / CL61 cloud {fv('payerne', 'A:ray/C:cloud', 'fallback, rejected')} vs "
      f"{fv('payerne', 'A:ray/C:cloud', 'strict, kept')}; Lindenberg CHM15k / CL61 cloud "
      f"{fv('lindenberg', '0:ray/C:cloud', 'fallback, rejected')} vs "
      f"{fv('lindenberg', '0:ray/C:cloud', 'strict, kept')}; Budapest B / A "
      f"{fv('budapest', 'B:ray/A:ray', 'fallback, rejected')} vs "
      f"{fv('budapest', 'B:ray/A:ray', 'strict, kept')}"
      + (f". Not every pair shows a bias ({SITE_LABEL[good.site]} "
         f"{pair_label(good.site, good.pair)}: {fv(good.site, good.pair, 'fallback, rejected')} vs "
         f"{fv(good.site, good.pair, 'strict, kept')}): v2.3 also drops some unbiased nights"
         if good is not None else "")
      + " (§4.6). The CL61's own fallback nights, which v2.3 does not test, keep a low Rayleigh "
      "constant at Payerne (§4.2).\n")
    cv = conv_summary(d, st)
    w("* **Caveats** (§5): abnormal Payerne CHM15k; short Payerne CL61 window; the CL31 dark and "
      f"swap; the 910 / 1064 nm conversion (σ_conv {rng_txt(cv['conv'])} % per night in the "
      "500–3000 m band) and the near-range shape; the cloud reference is mostly bridged and its "
      "reported uncertainty is a profile spread; the ratio of two Rayleigh calibrations cancels "
      "their common atmospheric errors.\n")
    return "\n".join(L)


def METHODS(d, st, rep, ct, cmp, br):
    pt = d["pairs"]
    L = []
    w = L.append
    nn = pt.groupby("site").date.nunique().to_dict()
    w("## 1. Question and design\n")
    w("Co-located instruments see the same atmosphere. If each is calibrated correctly, the ratio of "
      "their calibrated attenuated backscatter is 1, and its night-to-night spread should be what "
      "their stated calibration uncertainties predict. This is an independent truth for the new "
      "v2.3 per-night uncertainty (`uncertainty_model = \"v2.3\"`, doc/v2.3/uncertainty.md), which "
      "was validated on the time series of single units and on CHM15k twins only. Three questions:\n")
    w("1. Are the calibrated instruments consistent with each other within their stated "
      "uncertainties — first of all with the CL61 at Payerne?\n2. What is the bias between them?\n"
      "3. Does the v2.3 uncertainty match the observed disagreement, where the v2.2 ensemble "
      "uncertainty did not?\n")
    w("## 2. Data and calibrations\n")
    w("### 2.1 Sites, units, periods\n")
    w("| site | WIGOS id | units | nights read (first → last) | nights with a pair value "
      "(first → last) |\n|---|---|---|---|---|")
    ry = d["rays"]
    read = {}
    for sk, g in ry.groupby("site"):
        read[sk] = f"{g.date.nunique()} ({g.date.min()} → {g.date.max()})"
    units = {"payerne": "A CHM15k (Rayleigh; module TUB140016, abnormal noisy receiver), B CL31 "
                        "(cloud; optical block replaced 2026-07-07), C CL61 (cloud + Rayleigh, L1 "
                        "from 2026-02-24)",
             "lindenberg": "0 CHM15k (Rayleigh), C CL61 (cloud + Rayleigh, L1 to 2026-06-30)",
             "aosta": "0 CHM15k (Rayleigh), B CL61 (cloud + Rayleigh, L1 from 2025-10-31)",
             "schiphol": "A, B, C, D CHM15k (Rayleigh), a few km apart on the airport",
             "budapest": "A, B CHM15k (Rayleigh)", "magurele": "A, B CHM15k (Rayleigh)"}
    first_last = {}
    for sk in RC.SITES:
        ds = sorted(pd.unique(pt[pt.site == sk].date)) if len(pt) else []
        first_last[sk] = (ds[0], ds[-1]) if ds else ("–", "–")
    for sk, s in RC.SITES.items():
        w(f"| {s['name']} | {s['wmo']} | {units[sk]} | {read.get(sk, '–')} | {nn.get(sk, 0)} "
          f"({first_last[sk][0]} → {first_last[sk][1]}) |")
    w("\nThe CL61 of Payerne has daily L1 files from 2026-02-24 on (the 5-min granule layout noted "
      "earlier is no longer the only source), so its paired window is 2026-02-24 → 2026-10-08.\n")
    w("### 2.2 One read per night, everything in memory\n")
    w("`ref23_run.py` reads, for every night D and every unit of a site, the L1 files of D−1 and D "
      "ONCE (`load_instrument_day`, the runner's read-once path) and computes in memory: the two "
      "Rayleigh calibrations, the cloud calibration of day D and the comparison ingredients "
      "(12 workers, disk-bound: 46 min for the whole set). Nothing operational is written "
      "(the NetCDF write is patched out, as in the q16a/u23 harnesses).\n")
    w("**Rayleigh** (CHM15k, CL61), `calibrate_rayleigh` unmodified, `options.json`, CAMS "
      "`D:/CAMS_Monthly_04;D:/CAMS_daily` (0.4° only), CL61 with the operational water-vapour line "
      "`ALC_WV_SPECTRUM={\"CL61\": [910.55, 0.188]}`, dark-free as operational:\n")
    w("* **v2.2** — `eprof_v2.2` with the v2.3 budget computed on its window and the relative QC "
      "limit switched off (`max_rel_uncertainty = 1e9`); the v2.2 `u > C` gate is re-applied "
      "offline. This gives the published constant and window, its published ensemble uncertainty "
      "(= the budget's ensemble term) AND the v2.3 budget on the same window "
      "(\"v2.2 constant + v2.3 budget\", which isolates the effect of the uncertainty model).\n"
      "* **v2.3** — the release candidate: `eprof_v2.3` (repaired tier 2) + "
      "`uncertainty_model = \"v2.3\"` with the defaults **as revised by the independent review** "
      "(merged 2026-10-21, doc/v2.3/uncertainty.md §8): photon noise measured as the noise of the "
      "window mean on the binned stack (`noise_estimator = \"window_mean\"`), the vertical term "
      "saturating at 15 % (`vertical_cap_rel`), the molecular climatology evaluated at the window "
      "height above sea level; relative QC limit 60 % at 2 σ (it now almost never fires). This is the "
      "main v2.3 column. The **pre-review** budget (per-gate noise, no cap, window height AGL — what "
      "was first shipped) comes from the first pass of the harness and is the secondary column "
      "(\"v2.3 pre-review\"); its constants are the same, its σ and its QC-rejected nights differ.\n")
    w("Every budget component is stored as a 1-σ relative value (ensemble, vertical consistency — "
      "capped and uncapped —, molecular model, residual aerosol, photon noise), with the sub-window "
      "constants, the height gradient, the window and the tier, so the totals can be recomputed if "
      "the budget is revised again (`tables/rays.pkl` reviewed, `tables/rays_prereview.pkl` "
      "pre-review). The reviewed budget needed a second pass (`ref23_run.py run --rayleigh-only`): "
      "its noise term reads the binned profile stack and cannot be derived from the stored "
      "per-gate term. The tier is not readable from the v2.2 selection message (a fallback pick "
      "and a strict pick both say \"max composite quality\"); it is obtained by re-running the "
      "selector with the noise-aware tier disabled: strict-only success = tier 1. On every night "
      "where both versions calibrate, tier 1 of v2.3 equals tier 1 of v2.2 (same constant, same "
      "window) and the tiers agree.\n")
    w("**Reproduction of the published constants** (dashboard NetCDFs fetched 2026-10-09, "
      "`calibration_version` 220):\n")
    w(repro_md(rep))
    nonid = []
    for r in rep.itertuples():
        k = int(round(r.both * (1 - r.identical_pct / 100.0)))
        if k:
            nonid.append(f"{SITE_LABEL[r.site]} {r.type} {r.ident}: {k}")
    rn = d.get("reproduction_nights", pd.DataFrame())
    shared_txt = ""
    if len(rn):
        per_day = rn.groupby("date").site.nunique()
        shared = sorted(per_day[per_day >= 2].index)
        n_sh = int(rn.date.isin(shared).sum())
        if shared:
            days = ", ".join(f"{x[:4]}-{x[4:6]}-{x[6:]}" for x in shared)
            shared_txt = (f" {n_sh} of the {len(rn)} fall on {len(shared)} dates shared by two or more "
                          f"sites ({days}), which points to an input of those days that the "
                          "operational run saw differently (for instance the CAMS file available at "
                          "the time) rather than to the calibration code.")
    w("\nNon-identical common nights: " + ("; ".join(nonid) if nonid else "none") + "."
      + shared_txt + " For the three Payerne CL61 nights of June 2026, switching to daily-only CAMS "
      "does not reproduce the published value either. The operational run also applies the "
      "Cloudnet-classification screen, which this harness — like the q16a/u23 harnesses — does not. "
      "On the identical nights the published uncertainty equals the rerun's ensemble term.\n")
    w("### 2.3 Cloud references — what the published cloud uncertainty is\n")
    w("The liquid-cloud (O'Connor) calibration was rerun for every UTC day with the operational code "
      "path (`scripts/run_network_calibration._do_cloud`: the shared read sliced to day D, WV "
      "correction with the unit's line, PVC η tables, aerosol transmission correction). It "
      "reproduces the published cloud constants (`calibration_version` 1000) on every common day:\n")
    lines = ["| channel | days rerun | days published | common | identical (1e-6) | "
             "published u / rerun std·C/median |", "|---|---|---|---|---|---|"]
    from ref23_cloud import LABEL
    for r in cmp.itertuples():
        lines.append(f"| {LABEL[(r.site, r.ident)]} | {r.n_rerun} | {r.n_pub} | {r.n_both} | "
                     f"{f1(r.within_1e6, nd=1)} % | {f1(r.u_pub_vs_rerun, nd=3)} |")
    w("\n".join(lines) + "\n")
    w("**What the published cloud uncertainty is.** `lidar_constant_uncertainty` of a cloud day = "
      "C_L × std(C_i)/median(C_i), the C_i being the per-profile O'Connor coefficients of the "
      "day's accepted liquid-cloud profiles (`calibration/cloud/calibration.py`: `cal_std = "
      "np.std(valid_C, ddof=1)`; `run_network_calibration._do_cloud`: `cl_unc = cl × std/coef`). "
      "It is therefore (i) a **1-σ** number, whereas the Rayleigh product is 2 σ; (ii) the "
      "**spread of single profiles**, not the uncertainty of the day's median that is used; "
      "(iii) free of any systematic term (S = 18.8 sr, multiple-scattering η, water vapour, "
      "aerosol transmission). It does not decrease with the number of profiles or scenes "
      "(figure, right).\n")
    w("**Is it meaningful?** The test is the day-to-day reproducibility of the cloud constant: the "
      "structure function D(L) = robust variance of ln C(d+L) − ln C(d) over pairs of successful "
      "days L apart is fitted as 2σ_w² + q·L over L ≤ 7 d (σ_w = random error of one day value, an "
      "upper bound since it also holds real day-to-day change; q·L = real drift over L days).\n")
    w(cloud_md(ct, cmp))
    w(f"\n![Cloud-calibration uncertainty: day-to-day reproducibility against the reported spread "
      f"(left), reported values (middle), reported spread against the number of cloud scenes "
      f"(right)]({FIG['cloud_uncertainty']})\n")
    pr = ct.pull_rep_consec if len(ct) else pd.Series(dtype=float)
    sem_low = ct[ct.sig_sem_med < ct.sig_w_fit].label.tolist() if len(ct) else []
    sem_ok = ct[ct.sig_sem_med >= ct.sig_w_fit].label.tolist() if len(ct) else []
    w("The reported spread is of the size of the observed day-to-day reproducibility or larger "
      f"(consecutive-day pulls of the day difference with the reported σ {f1(pr.min(), nd=2)}–"
      f"{f1(pr.max(), nd=2)}, i.e. ≤ 1): as a per-day 1 σ it is **usable and conservative**, although "
      "for the wrong reason — the day median is dominated by scene-to-scene variability (cloud "
      "type, droplet size through η, aerosol below cloud), which the profile spread happens to "
      "sample, while the standard error of the median over independent scenes (counted as runs of "
      "accepted profiles separated by more than 10 min) is smaller than the fitted day-to-day σ_w "
      f"for {', '.join(sem_low) if sem_low else 'no channel'}"
      + (f" (not for {', '.join(sem_ok)})" if sem_ok else "") + ". The reported spread is used "
      "per day as the cloud 1 σ; the empirical constant σ_w and the standard error are "
      "sensitivities.\n")
    w("**Bridging.** A clear night — what a Rayleigh calibration needs — rarely has a liquid-cloud "
      "calibration on the same days. The cloud value of night D (the night D−1 → D) is the cloud "
      "calibration of day D and/or D−1 (log-mean when both exist, σ/√k); otherwise the nearest "
      "successful day(s) within ±7 days of the night (same era for the CL31), with "
      "σ = √(σ_day² + q·L), q the channel's drift rate from the structure function, L the distance "
      "in days. Sensitivity: the median of every successful cloud day within ±7 days "
      "(σ = 1.2533 σ_day/√n ⊕ drift). The Payerne CL31 is never bridged across its 2026-07-07 swap "
      "(the swap day is excluded).\n")
    if len(br):
        lines = ["| pair | v2.3 nights | same-night cloud value | bridged | median bridge lag [d] | "
                 "median σ of the cloud reference [%] |", "|---|---|---|---|---|---|"]
        for r in br.itertuples():
            lines.append(f"| {SITE_LABEL[r.site]} {pair_label(r.site, r.pair)} | {r.n} | {r.same} | "
                         f"{r.bridge} | {f1(r.lag_med, nd=0)} | {f1(r.sig_med)} |")
        w("\n".join(lines) + "\n")
    w("## 3. Observable and statistics\n")
    w("### 3.1 ln ρ\n")
    w("For night n and a pair X / R (R = the reference: the CL61 at the mixed sites):\n")
    w("  ρ(n) = median over the admitted pixels of β_X(h, z) / β_R(h, z),   "
      "β = rcs_0 × 10⁶ / C_L(n)\n")
    w("* **each unit divided by ITS OWN constant of night n** (the cloud reference: its night value, "
      "bridged if missing — §2.3);\n"
      "* β of the 910 nm units divided by their two-way water-vapour transmission (CAMS, the "
      "unit's own line: CL61 910.55 nm / FWHM 0.188 nm, CL31 909.7 / 6.0 nm); a CL61 or CL31 "
      "compared with a CHM15k is converted to 1064 nm with the v3 dashboard's component-separated "
      "rule β₁₀₆₄ = β_mol,1064 + (β₉₁₀ − β_mol,910)·(910/1064)^α (CAMS T/p molecular attenuated "
      "backscatter, α = 1); two 910 nm units are compared at 910 nm;\n"
      "* pixels = hourly medians of rcs_0 on the v3 display grid (30 m up to 3 km), hours of the "
      "**calibration night** (SZA > 100°, the night the constants come from), with ≥ 30 min of "
      "samples; **clear for every co-located unit** (union of the cloud bases and of the Vaisala "
      "vertical visibility of all units, ±15 min, the v3 screen); **both units detected at "
      "SNR ≥ 3** with the v3 `common` rule (5-min windows, an hour admitted when ≥ 50 % of its "
      "data are detected by both);\n"
      "* bands: **500–3000 m** (the v3 dashboard band) and **1500–3000 m** (above full overlap of "
      "every unit); **clean air**: only the gates whose night-median scattering ratio at 1064 nm "
      "is below 1.5, where the wavelength conversion hardly depends on α;\n"
      "* ln ρ(z): the same ratio per gate, median over the night's admitted hours (figure 4);\n"
      "* dark-free, as operational. At Payerne the measured covered-telescope dark is subtracted "
      "from the profiles of the cloud-calibrated units as a sensitivity (their constants are "
      "dark-immune, so this is the coherent state of a cloud channel).\n")
    w("### 3.2 Predicted uncertainty and common terms\n")
    w("  σ_pred(n) = √(σ_X² + σ_R² + σ_conv²), all 1 σ relative, from the per-night components.\n")
    wv = d.get("wv_transmission", pd.DataFrame())
    c61 = wv[wv.type == "CL61"] if len(wv) else wv
    wv_txt = (f"The two-way water-vapour transmission of the narrow CL61 line is ≥ "
              f"{f1(c61.t2_30_p05.min(), nd=2)} at 3 km (5th percentile of the nights); a ±20 % "
              f"humidity error moves it by ≤ {f1(c61.dlnr_20pct_30_p05.max())} % there (§4.8), which "
              "is neglected." if len(c61) else "")
    w("* **Rayleigh vs Rayleigh** (CHM15k twins, CHM15k vs CL61 Rayleigh): the molecular-model "
      "term is dropped for both units — both use US-1976 instead of the real T/p over similar "
      "windows the same night, so its error is common and cancels in the ratio (doc/v2.3 drops it "
      "for twins); \"with mol\" is the sensitivity. The other terms are treated as independent, "
      "although the vertical and residual-aerosol terms are partly atmospheric and therefore "
      "partly common too (see §4.4).\n"
      "* **Rayleigh vs cloud**: the full Rayleigh budget (the cloud method has no molecular model) "
      "+ the cloud reference σ (§2.3).\n"
      "* **σ_conv** — the wavelength conversion: |ln ρ(α = 1.5) − ln ρ(α = 0.5)| / 2 per night, "
      "i.e. α = 1 ± 0.5 (1 σ) for the 910 → 1064 nm Ångström exponent of the boundary-layer "
      f"aerosol. {wv_txt}\n"
      "* v2.2: σ_X = published `lidar_constant_uncertainty` / (2 C_L) (the ensemble), same σ_R "
      "(cloud) and σ_conv.\n")
    w("### 3.3 Statistics\n")
    w("Per site × pair × version: **n** nights; **bias** = median ln ρ with a 95 % CI from a "
      "block bootstrap over calendar weeks (consecutive nights share aerosol and gain state); "
      "**observed σ** = 1.4826 × MAD of ln ρ; **predicted σ** = median σ_pred; **pulls** "
      "z = ln ρ / σ_pred (uncentred: bias and scatter together) and z_c = (ln ρ − median) / σ_pred "
      "(centred: scatter only), each summarised by its robust std about zero (1.4826 × median|z|; "
      "1 for a calibrated Gaussian uncertainty) and the fraction |z| < 2 (95.4 %); **Spearman** "
      "rank correlation of |ln ρ − median| with σ_pred (does a larger stated uncertainty flag a "
      "larger disagreement?).\n")
    w("Verdict rules: **consistent** = no significant bias (CI contains 0), ≥ 90 % of the nights "
      "within 2 σ_pred (uncentred) and a centred pull robust std in [0.75, 1.33]; otherwise "
      "**significant bias** (CI excludes 0; \"within / beyond the stated u\" according to the "
      "≥ 90 % coverage), **u under-stated** (centred pull > 1.33) or **u over-stated** (< 0.75).\n")
    return "\n".join(L)


def _pl(s, s22=None, key="pull_c"):
    """'a → b' for a pull from the v2.2 and v2.3 stats."""
    return f"{f1(s22.get(key), nd=2) if s22 else '–'} → {f1(s.get(key), nd=2) if s else '–'}"


def RESULTS(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, tt):
    pt = d["pairs"]
    L = ["## 4. Results\n"]
    w = L.append
    w("### 4.1 Summary\n")
    w("One row per site × pair; v2.3 statistics, with the published v2.2 values where an arrow "
      "is shown (each version on its own nights). Band 500–3000 m, α = 1, 5-min SNR windows.\n")
    w(summary_table(st))
    w("\n![ln ρ per night with ±2 σ_pred: CL61 sites (v2.3 red, filled = strict tier, open = "
      f"fallback tier; v2.2 published grey)]({FIG['timeseries_cl61']})\n")
    w(f"![ln ρ per night with ±2 σ_pred: CHM15k twins]({FIG['timeseries_twins']})\n")
    w("Bias by band, in clean air and its wavelength-conversion sensitivity (v2.3):\n")
    w(bias_table(st, MAIN_PAIRS))
    w("\nSensitivity of the centred pull (v2.3) to the analysis choices:\n")
    w(sensitivity_table(st, MAIN_PAIRS + sorted(T)))
    w("\n" + RESULTS_TEXT(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, tt))
    return "\n".join(L)


def _within(s):
    """'within' / 'not within' the stated uncertainty, by the uncentred coverage rule."""
    if not s:
        return "not evaluable"
    return "within" if s.get("in2_u", 0) >= S.COVER_OK else "not within"


def RESULTS_TEXT(d, st, pp, P, P22, pc, pc22, pcb, T, T22, fate, tt):
    pt = d["pairs"]
    L = []
    w = L.append
    k_ac, k_acr, k_self = ("payerne", "A:ray/C:cloud"), ("payerne", "A:ray/C:ray"), \
        ("payerne", "C:ray/C:cloud")
    s_ac, s_acr, s_self = P[k_ac], P[k_acr], P[k_self]
    s22_ac, s22_acr, s22_self = P22[k_ac], P22[k_acr], P22[k_self]
    dk = _s(st, *k_ac, variant="dark")
    dkhi = _s(st, *k_ac, variant="dark_hi")
    hi = _s(st, *k_ac, variant="hi_band")
    tier_ac = tt[(tt.site == "payerne") & (tt.pair == "A:ray/C:cloud") & (tt.version == "v23")]
    t1 = tier_ac[tier_ac.split == "strict"]
    t1 = t1.iloc[0].to_dict() if len(t1) else {}
    w("### 4.2 Payerne — agreement with the CL61\n")
    w(f"**CHM15k A (Rayleigh) / CL61 C (cloud).** v2.3: {s_ac.get('n', 0):.0f} nights, bias "
      f"{_b(s_ac)}, observed σ {f1(s_ac.get('scatter'))} % against a predicted "
      f"{f1(s_ac.get('sig_pred'))} %, centred pull {_pl(s_ac, s22_ac)} (v2.2 → v2.3), "
      f"{f1(s22_ac.get('in2_u'), nd=0)} → {f1(s_ac.get('in2_u'), nd=0)} % of the nights within 2 σ. "
      f"The agreement is **{_within(s_ac)} the v2.3 uncertainty** by the coverage rule; the scatter "
      f"{'is' if s_ac.get('pull_c', 9) <= S.PULL_HI else 'is still not'} predicted "
      f"(centred pull {f1(s_ac.get('pull_c'), nd=2)}). The v2.2 series ({s22_ac.get('n', 0):.0f} "
      f"nights, bias {_b(s22_ac)}) contained many fallback-tier nights that v2.3 removes (§4.6). "
      f"On the strict-tier nights alone: {t1.get('n', 0):.0f} nights, bias {f1(t1.get('bias'), True)} %, "
      f"centred pull {f1(t1.get('pull_c'), nd=2)}. Above full overlap (1500–3000 m) the bias is "
      f"{_b(hi)}; with the CL61's measured dark subtracted from its profile "
      f"{f1(dk.get('bias'), True)} % (500–3000 m) and {f1(dkhi.get('bias'), True)} % (1500–3000 m).\n")
    w(f"**CHM15k A / CL61 C (Rayleigh).** {s_acr.get('n', 0):.0f} nights, bias {_b(s_acr)}, "
      f"observed σ {f1(s_acr.get('scatter'))} %, predicted {f1(s_acr.get('sig_pred'))} %, centred "
      f"pull {_pl(s_acr, s22_acr)}, within 2 σ {f1(s22_acr.get('in2_u'), nd=0)} → "
      f"{f1(s_acr.get('in2_u'), nd=0)} %.\n")
    ts = tt[(tt.site == "payerne") & (tt.pair == "C:ray/C:cloud") & (tt.version == "v23")]
    ts = ts.set_index("split")
    tb = tier_budget(d["rays"], "payerne", "C")
    if {"strict", "fallback"} <= set(ts.index) and {1, 2} <= set(tb):
        s1_, s2_, b1_, b2_ = ts.loc["strict"], ts.loc["fallback"], tb[1], tb[2]
        tier_txt = (
            f"The CL61's own fallback-tier nights, which v2.3 keeps (the vertical-consistency test is "
            f"not applied to the CL61), sit at {f1(s2_.bias, True)} % ({s2_.n:.0f} nights; strict-tier "
            f"nights {f1(s1_.bias, True)} %, {s1_.n:.0f} nights). They use higher windows (median base "
            f"{f1(b2_['zb'])} km vs {f1(b1_['zb'])} km) with a steep negative gradient of the "
            f"sub-window constants (median {f1(b2_['g'], True)} %/km vs {f1(b1_['g'], True)} %/km): "
            "the signature of the range-growing CL61 dark baseline measured at Payerne (about −5 % of "
            "the molecular signal at 3 km, −25 % at 7 km; doc/reports/cl61_cloud_vs_rayleigh_origin.md "
            "§2.4), which lowers the dark-free Rayleigh constant of a high window. The v2.3 budget "
            f"gives these nights a larger σ (median vertical {f1(b2_['sig_vert'])} %, noise "
            f"{f1(b2_['sig_noise'])} %, ensemble {f1(b2_['sig_ens'])} %, against "
            f"{f1(b1_['sig_vert'])} / {f1(b1_['sig_noise'])} / {f1(b1_['sig_ens'])} % on the strict "
            "nights), not enough to cover a bias of that size.")
        tl = tt[(tt.site == "lindenberg") & (tt.pair == "C:ray/C:cloud") & (tt.version == "v23")]
        tl = tl.set_index("split")
        if {"strict", "fallback"} <= set(tl.index):
            tier_txt += (f" At Lindenberg the CL61's fallback nights sit at "
                         f"{f1(tl.loc['fallback'].bias, True)} % ({tl.loc['fallback'].n:.0f} nights) "
                         f"against {f1(tl.loc['strict'].bias, True)} % on its strict nights: the "
                         "bias belongs to the Payerne unit, not to the CL61 fallback tier as such.")
    else:
        tier_txt = ""
    w(f"**CL61 C Rayleigh / CL61 C cloud** (same photons: ln ρ = ln C_cloud − ln C_Rayleigh). "
      f"{s_self.get('n', 0):.0f} nights, bias {_b(s_self)}, centred pull {_pl(s_self, s22_self)}, "
      f"within 2 σ {f1(s22_self.get('in2_u'), nd=0)} → {f1(s_self.get('in2_u'), nd=0)} %. "
      + tier_txt + "\n")
    ce = d["summary_cl31_era"]

    def era(pair, e, var):
        m = ce[(ce.pair == pair) & (ce.era == e) & (ce.variant == var)]
        return m.iloc[0].to_dict() if len(m) else {}
    b0, b0d = era("B:cloud/C:cloud", "before 2026-07-07", "main"), \
        era("B:cloud/C:cloud", "before 2026-07-07", "dark")
    b1 = era("B:cloud/C:cloud", "after 2026-07-07", "main")
    a0, a0d = era("A:ray/B:cloud", "before 2026-07-07", "main"), \
        era("A:ray/B:cloud", "before 2026-07-07", "dark")
    w("**CL31 B (cloud) / CL61 C (cloud)** — two 910 nm units compared at 910 nm, both "
      "cloud-calibrated. Before the CL31's optical-block swap (2026-07-07): "
      f"{b0.get('n', 0):.0f} nights, {f1(b0.get('bias'), True)} % (observed σ "
      f"{f1(b0.get('scatter'))} %) dark-free, {f1(b0d.get('bias'), True)} % ({f1(b0d.get('scatter'))} %, "
      f"centred pull {f1(b0d.get('pull_c'), nd=2)}) with the measured covered-telescope dark of "
      f"both units subtracted. After the swap: {b1.get('n', 0):.0f} nights, "
      f"{f1(b1.get('bias'), True)} % (observed σ {f1(b1.get('scatter'))} %, centred pull "
      f"{f1(b1.get('pull_c'), nd=2)}, {f1(b1.get('in2_u'), nd=0)} % within 2 σ) without any dark "
      "correction — the new block's dark has not been measured, and the old profile is not applied "
      "to it. The CL31's cloud calibration is thus consistent with the CL61's once its profile is "
      "clean; its old block carried a dark large enough to make the dark-free comparison "
      f"meaningless. CHM15k A / CL31 B before the swap: {f1(a0.get('bias'), True)} % dark-free, "
      f"{f1(a0d.get('bias'), True)} % with the CL31 dark.\n")
    # --- Lindenberg / Aosta ---------------------------------------------------------------------
    w("### 4.3 Lindenberg and Aosta\n")
    for site in ("lindenberg", "aosta"):
        c61 = "C" if site == "lindenberg" else "B"
        kc, kr, ks = (site, f"0:ray/{c61}:cloud"), (site, f"0:ray/{c61}:ray"), \
            (site, f"{c61}:ray/{c61}:cloud")
        sc, sr, ss = P.get(kc, {}), P.get(kr, {}), P.get(ks, {})
        sc22, sr22, ss22 = P22.get(kc, {}), P22.get(kr, {}), P22.get(ks, {})
        rc, rch = _s(st, *kr, variant="clean"), _s(st, *kr, variant="clean_hi")
        cc, cch = _s(st, *kc, variant="clean"), _s(st, *kc, variant="clean_hi")
        w(f"**{SITE_LABEL[site]}.** CHM15k / CL61 Rayleigh: {sr.get('n', 0):.0f} nights, bias "
          f"{_b(sr)} in 500–3000 m but {_b(rch)} in clean air above 1500 m; observed σ "
          f"{f1(sr.get('scatter'))} % vs predicted {f1(sr.get('sig_pred'))} %; centred pull "
          f"{_pl(sr, sr22)}; within 2 σ {f1(sr22.get('in2_u'), nd=0)} → {f1(sr.get('in2_u'), nd=0)} %. "
          f"CHM15k / CL61 cloud: {sc.get('n', 0):.0f} nights, bias {_b(sc)} (clean air above "
          f"1500 m {_b(cch)}); centred pull {_pl(sc, sc22)}; within 2 σ "
          f"{f1(sc22.get('in2_u'), nd=0)} → {f1(sc.get('in2_u'), nd=0)} %. CL61 Rayleigh / CL61 "
          f"cloud: {ss.get('n', 0):.0f} nights, bias {_b(ss)}, centred pull {_pl(ss, ss22)}, within "
          f"2 σ {f1(ss22.get('in2_u'), nd=0)} → {f1(ss.get('in2_u'), nd=0)} %.\n")
    # --- twins ------------------------------------------------------------------------------------
    w("### 4.4 CHM15k twins (Schiphol, Budapest, Măgurele)\n")
    lines = ["| pair | n | bias 500–3000 m [%] | bias 1500–3000 m [%] | observed σ [%] | "
             "predicted σ v2.2 → v2.3 [%] | centred pull v2.2 → v2.3 | within 2 σ v2.2 → v2.3 |",
             "|" + "---|" * 8]
    for k in sorted(T):
        s, s22 = T[k], T22[k]
        if not s:
            continue
        h = _s(st, *k, variant="hi_band")
        lines.append(f"| {SITE_LABEL[k[0]]} {pair_label(*k)} | {s['n']:.0f} | {_b(s)} | {_b(h)} | "
                     f"{f1(s['scatter'])} | {f1(s22.get('sig_pred'))} → {f1(s['sig_pred'])} | "
                     f"{_pl(s, s22)} | {f1(s22.get('in2_u'), nd=0)} → {f1(s['in2_u'], nd=0)} % |")
    w("\n".join(lines) + "\n")
    tw, tw22 = pc.get("CHM15k twins") or {}, pc22.get("CHM15k twins") or {}
    w(f"Pooled over the {len([k for k in T if T[k]])} twin pairs the centred pull goes from "
      f"{f1(tw22.get('pull_c'), nd=2)} (v2.2) to {f1(tw.get('pull_c'), nd=2)} (v2.3): the v2.2 "
      "ensemble under-states the twin scatter, the v2.3 budget over-covers it "
      f"({f1(tw.get('in2_c'), nd=0)} % of the pair-centred nights within 2 σ). This is expected and "
      "is not a defect for the single-unit use of the uncertainty. The ratio of two calibrated "
      "profiles measures only the part of the two calibration errors that differs between the "
      "units: the atmospheric part (aerosol in the window, molecular model — the same air and "
      "similar windows the same night) cancels, and so does any real gain change of either unit, "
      "which its constant and its signal follow alike. The single-unit budget must cover the "
      "atmospheric part, so it over-covers a twin ratio. The review of the budget puts the twin "
      "pull of a perfectly calibrated budget anywhere in 0.40–0.77 on its twin-constant truth "
      f"(doc/v2.3/uncertainty.md §4); {f1(tw.get('pull_c'), nd=2)} is "
      f"{'inside' if 0.40 <= tw.get('pull_c', np.nan) <= 0.77 else 'outside'} that range. "
      "The static twin biases are "
      "instrumental: Schiphol B, whose applied overlap is ~25 % too low at 500–1000 m, reads high "
      "in the 500–3000 m band and much less above 1500 m.\n")
    # --- calibration of the uncertainty -------------------------------------------------------
    w("### 4.5 Is the stated uncertainty calibrated?\n")
    w("Pooled groups (each pair centred on its own median for the centred statistics; Payerne "
      "CHM15k A separately):\n")
    w(pooled_table(pp))
    parts = []
    for grp in ("CHM15k Rayleigh vs CL61 cloud", "CHM15k Rayleigh vs CL61 Rayleigh",
                "CL61 Rayleigh vs CL61 cloud", "CHM15k twins", "Payerne CHM15k A (abnormal) vs CL61"):
        a_, b_ = pooled(pp, grp, "v23pre"), pooled(pp, grp, "v23")
        if a_ and b_:
            parts.append(f"{grp} {f1(a_['pull_c'], nd=2)} → {f1(b_['pull_c'], nd=2)} "
                         f"(median σ_pred {f1(a_['sig_pred'])} → {f1(b_['sig_pred'])} %, n "
                         f"{a_['n']:.0f} → {b_['n']:.0f})")
    w("\nPre-review → reviewed v2.3 budget (centred pull; the n changes because the relative QC "
      "limit fired on more nights with the pre-review budget): " + "; ".join(parts) + ".\n")
    w("\nLike for like — only the nights both the published v2.2 and the v2.3 series calibrate "
      "(on strict-tier nights the two constants are identical, so this isolates the uncertainty "
      "model):\n")
    w(pooled_table(pp, "common"))
    w("\nSensitivity of the pooled centred pull (v2.3) to the common-term choice and to the "
      "conversion term:\n")
    lines = ["| group | main | with the molecular term | without σ_conv | 1500–3000 m |",
             "|---|---|---|---|---|"]
    for grp in pp.group.unique():
        cells = []
        for var in ("main", "with_mol", "no_conv", "hi_band"):
            r = pooled(pp, grp, "v23", var)
            cells.append(f1(r["pull_c"], nd=2) if r else "–")
        lines.append(f"| {grp} | " + " | ".join(cells) + " |")
    w("\n".join(lines))
    w(f"\n![Pull distributions, v2.2 published (grey) and v2.3 (red), uncentred (top) and centred "
      f"on each pair's median (bottom)]({FIG['pulls']})\n")
    w(f"![Reliability: observed robust RMS of ln ρ (pair-centred) against the binned predicted σ; "
      f"1:1 = calibrated]({FIG['reliability']})\n")
    w("Split by the tier of the Rayleigh unit X (v2.3):\n")
    w(tier_table(tt))
    w("\nBy season (v2.3; bias of ln ρ and centred pull; only seasons with ≥ 6 nights):\n")
    w(season_table(d["summary_season"], MAIN_PAIRS))
    w("\n")
    w(f"![Median v2.3 budget components of the Rayleigh units of this study]({FIG['budget']})\n")
    # --- fate -----------------------------------------------------------------------------------
    w("### 4.6 What v2.3 did to the v2.2 nights — judged by the independent reference\n")
    w("Each v2.2 night of the Rayleigh unit X is classed by its fate in v2.3; the median ln ρ is "
      "given with the v2.2 constant and, for the kept nights, with the v2.3 constant.\n")
    w(fate_md(fate))
    w("\n")
    # --- height ---------------------------------------------------------------------------------
    w("### 4.7 Height dependence\n")
    w(f"![ln ρ(z): median and inter-quartile range over the nights, per pair (height on the Y "
      f"axis; dotted lines: 0.5, 1.5 and 3 km)]({FIG['profiles']})\n")
    w("The ratio profiles are not flat below ~1.5 km: the CHM15k overlap correction, the CL61 "
      "near-range response and, at Schiphol B, the known overlap deficit shape the lowest kilometre, "
      "and the wavelength conversion weighs most where the aerosol is (the boundary layer). The "
      "500–3000 m band therefore mixes a near-range shape difference into the bias; the 1500–3000 m "
      "band is the one that tests the calibration constants (table of §4.1).\n")
    # --- lambda -----------------------------------------------------------------------------------
    w("### 4.8 The wavelength conversion\n")
    cv = conv_summary(d, st)
    if np.isfinite(cv["slope"][0]):
        w(f"Between α = 0.5 and α = 1.5 the median ln ρ of the CHM15k / CL61 pairs moves by "
          f"{rng_txt(cv['slope'])} % (500–3000 m): the absolute bias of a 910 / 1064 nm comparison "
          "in the boundary layer is only defined to that order unless the Ångström exponent is "
          f"known. Per night, σ_conv (α = 1 ± 0.5) has a median of {rng_txt(cv['conv'])} % in the "
          f"500–3000 m band, against a total predicted σ of {rng_txt(cv['pred'])} %: it **must "
          "enter** the predicted uncertainty of a 910 / 1064 nm comparison (without it the centred "
          "pulls rise, see the sensitivity tables). In clean air (scattering ratio < 1.5) it falls "
          f"to {rng_txt(cv['conv_clean'])} % and is negligible there.\n")
    wv = d.get("wv_transmission", pd.DataFrame())
    if len(wv):
        c61 = wv[wv.type == "CL61"]
        w(f"The two-way water-vapour transmission of the narrow CL61 line (910.55 nm, FWHM 0.188 nm) "
          f"is {f1(c61.t2_15_med.min(), nd=3)}–{f1(c61.t2_15_med.max(), nd=3)} at 1.5 km and "
          f"{f1(c61.t2_30_med.min(), nd=3)}–{f1(c61.t2_30_med.max(), nd=3)} at 3 km (median nights; "
          f"5th percentile at 3 km ≥ {f1(c61.t2_30_p05.min(), nd=3)}): a ±20 % humidity error moves the "
          f"CL61 backscatter by at most {f1(c61.dlnr_20pct_30_p05.max())} % at 3 km, less in the "
          "band — neglected. (The CL31's wide line absorbs more: "
          + "; ".join(f"T² {f1(r.t2_30_med, nd=3)} at 3 km" for r in wv[wv.type == 'CL31'].itertuples())
          + ".)\n")
    w(f"![Wavelength conversion: median ln ρ against the Ångström exponent (left), the per-night "
      f"conversion term (middle) and against the CHM15k calibration uncertainty (right)]"
      f"({FIG['lambda']})\n")
    return "\n".join(L)


def CAVEATS(d, st, P, rep, br):
    pt = d["pairs"]
    pay_c_nights = int(pt[(pt.site == "payerne") & (pt.pair == "A:ray/C:cloud")]
                       ["lnr_v23"].notna().sum()) if "lnr_v23" in pt else 0
    ce = d["summary_cl31_era"]

    def era(e, var):
        m = ce[(ce.pair == "B:cloud/C:cloud") & (ce.era == e) & (ce.variant == var)]
        return m.iloc[0].to_dict() if len(m) else {}
    bc, bcd = era("before 2026-07-07", "main"), era("before 2026-07-07", "dark")
    L = ["## 5. Caveats\n"]
    w = L.append
    bu = d["budget_units"]
    fa = bu[(bu.site == "payerne") & (bu.ident == "A")]
    fa = float(fa.frac_tier2.iloc[0]) if len(fa) else np.nan
    w("* **Payerne CHM15k A is abnormal.** Its receiver (module TUB140016) is noisy (night "
      "background 5–25 × the network) and loses its temperature regulation on hot afternoons; "
      f"{pct(fa)} of its v2.3 nights (and more of its v2.2 nights) come from the fallback tier, "
      "whose constants sit far from the CL61 (§4.6). It is reported separately and kept out of the "
      "pooled groups; the vertical-consistency study recommends keeping it on tier 1 only "
      "(doc/v2.3/vertical_consistency.md).\n")
    w(f"* **Short Payerne CL61 window.** The CL61 exists from 2026-02-24; {pay_c_nights} nights carry "
      "a v2.3 CHM15k / CL61-cloud value. Its diode temperature regulation was lost on 2026-06-11 "
      "(doc/reports/cl61_cloud_vs_rayleigh_origin.md), so a summer spectral drift of the 910 nm "
      "line cannot be excluded; the summer nights are kept.\n")
    stp, nb_, na_ = cl31_step(d["clouds"])
    ba = era("after 2026-07-07", "main")
    w("* **Payerne CL31: two eras and a large dark.** The optical block was replaced on 2026-07-07 "
      f"(cloud constant × {f1(stp, nd=2)}, medians of the {nb_} / {na_} successful days within 30 "
      "days before / after): the swap day is excluded and no cloud value is carried across it. Its "
      "cloud uncertainty is usable (§2.3), but the old block's profile carried a large electronic "
      "dark (−14 … −60 % of the molecular signal in the covered-telescope campaign): before the "
      f"swap the CL31 / CL61 ratio is {f1(bc.get('bias'), True)} % (observed σ "
      f"{f1(bc.get('scatter'))} %) without the measured dark and {f1(bcd.get('bias'), True)} % "
      f"({f1(bcd.get('scatter'))} %) with it; after the swap it is {f1(ba.get('bias'), True)} % "
      f"({f1(ba.get('scatter'))} %) without any dark correction (the new block's dark has not been "
      "measured). The CL31 summary rows mix both eras and the dark-free state; §4.2 gives them per "
      "era.\n")
    cv = conv_summary(d, st)
    w("* **Wavelength conversion.** In the 500–3000 m band the CHM15k / CL61 ratio moves by "
      f"{rng_txt(cv['slope'])} % per unit of Ångström exponent (boundary-layer aerosol): with "
      "α = 1 ± 0.5 the absolute bias of the mixed-wavelength pairs is only defined to about half "
      f"that there. It enters σ_pred per night (median σ_conv {rng_txt(cv['conv'])} %); in clean "
      f"air (scattering ratio < 1.5) it falls to {rng_txt(cv['conv_clean'])} % (§4.8).\n")
    w("* **Overlap and near range.** ln ρ(z) is not flat below ~1.5 km (CHM15k overlap, the known "
      "25 % overlap deficit of Schiphol B at 500–1000 m, CL61 near-range response): the 500–3000 m "
      "band mixes a near-range shape difference into the bias. The 1500–3000 m band is the "
      "calibration-relevant one (§4.7).\n")
    w("* **The cloud reference is mostly bridged.** Clear nights rarely have liquid clouds on the "
      "same days; the bridged value carries the measured drift of the cloud constant. The cloud "
      "method's own systematic floor (S = 18.8 ± 0.8 sr → ±4.3 %, the η table) is not in its "
      "reported uncertainty, so a bias of that size between a Rayleigh and a cloud calibration is "
      "not a disagreement.\n")
    w("* **Co-located Rayleigh calibrations are not independent.** Co-located units see the same "
      "aerosol and pick similar windows the same night, so the atmospheric part of the budget "
      "(vertical term, residual aerosol, molecular model) is largely common and cancels in their "
      "ratio; the budget treats it as independent and must over-cover twins (a perfectly "
      "calibrated budget gives a twin pull of 0.40–0.77, doc/v2.3/uncertainty.md §4) and, to a "
      "lesser degree, the CHM15k / CL61 Rayleigh pairs. Only the comparison with the cloud "
      "calibration is free of this. The Schiphol units are a few km apart, which adds real "
      "aerosol differences.\n")
    w("* **Classification screen and published nights.** The harness does not apply the operational "
      "Cloudnet-classification screen (a handful of published nights are not reproduced, §2.2), "
      "and the published v2.2 series contains nights that v2.3 rejects; both versions are therefore "
      "also compared on their common nights.\n")
    m0, mh = _s(st, "payerne", "A:ray/C:cloud"), _s(st, "payerne", "A:ray/C:cloud", variant="hi_band")
    k0, kh = (_s(st, "payerne", "A:ray/C:cloud", variant="dark"),
              _s(st, "payerne", "A:ray/C:cloud", variant="dark_hi"))

    def dlt(a, b):
        return (f1(b["bias"] - a["bias"], True) if a and b and np.isfinite(a.get("bias", np.nan))
                and np.isfinite(b.get("bias", np.nan)) else "–")
    w("* **Dark.** Everything is dark-free, as operational. The covered-telescope campaign measured, "
      "in the 3–7 km Rayleigh windows, a baseline of −16 % (CL61, growing with range) and −17 % "
      "(CHM15k, flat pedestal) of the night molecular signal at Payerne "
      "(doc/reports/cl61_cloud_vs_rayleigh_origin.md §2.4). It biases the dark-free Rayleigh "
      "constants of these units (in v2.2 and v2.3 alike) and was not removed here; no other "
      "site has a dark measurement. Subtracting the CL61 dark from its cloud-calibrated profile "
      "(whose constant is dark-immune) moves the CHM15k A / CL61-cloud ratio by "
      f"{dlt(m0, k0)} points in 500–3000 m and {dlt(mh, kh)} points in 1500–3000 m: in the "
      "aerosol-rich band the signal is large against the dark; in clean air aloft it is not.\n")
    w("* **SNR conditioning and the site-wide screen.** Only pixels both units detect (SNR ≥ 3, "
      "5 min) enter, and an hour is used only if every co-located unit is cloud-free: the sample is "
      "common to the pair, but conditioned on signal; the 60-min SNR window is a sensitivity.\n")
    return "\n".join(L)


def REPRO():
    return """## 6. Reproduction

```
# one read per night of every unit; Rayleigh v2.2 + v2.3 (budget components), cloud, ingredients
# (12 workers; ~45 min on the home PC, disk-bound):  runs/<site>/<YYYYMM>.pkl  (pre-review budget)
python validation/v23_reference/ref23_run.py run --workers 12
# the Rayleigh calibrations again with the REVIEWED budget (merged 2026-10-21; ~15 min):
#   runs_rayleigh_reviewed/
python validation/v23_reference/ref23_run.py run --workers 12 --rayleigh-only
python validation/v23_reference/ref23_check.py        # every night complete vs the L1 archive now
python validation/v23_reference/ref23_tables.py       # tables/rays.pkl, clouds.pkl, published.pkl, ...
python validation/v23_reference/ref23_analysis.py     # pairs.csv, summary*.csv, fate, figures
python validation/v23_reference/ref23_report.py       # this report
```

Data integrity: `ref23_check.py` confirms that every stream with an L1 file now was read and processed
on every night. The ten Schiphol months that were running while the Amsterdam L1 folder was being moved
to another disk (2026-10-21, 20:45) were re-run (`runs_superseded/` keeps the first version): identical
on all 297 nights.

| file | content |
|---|---|
| `pairs.csv` | one row per site × pair × night: ln ρ per calibration state, band, α, SNR window, clean air, dark; every σ term; tiers; the cloud reference (kind, lag, σ) |
| `summary.csv` | statistics per site × pair × version × variant (main, with mol, no conv, 1500–3000 m, 60-min SNR, clean air, dark, cloud-σ models, local cloud reference, common nights) |
| `summary_tier.csv`, `summary_season.csv`, `summary_pooled.csv`, `reliability.csv` | splits, pooled groups, reliability bins |
| `fate_v22_to_v23.csv` | what v2.3 did to each v2.2 night, judged by the reference |
| `reproduction.csv`, `reproduction_nights.csv`, `cloud_vs_published.csv` | reruns against the published constants (the non-identical nights listed) |
| `wv_transmission.csv`, `summary_cl31_era.csv` | CL61 / CL31 water-vapour transmission; Payerne CL31 pairs per optical-block era |
| `cloud_uncertainty.csv`, `cloud_structure.csv` | the cloud uncertainty check |
| `budget_units.csv`, `lambda_alpha.csv`, `profiles_median.csv` | figure data |
| `tables/rays.pkl` | every Rayleigh night, both versions, every budget component (recompute the totals from here if the budget is revised) |
"""


def TEMPLATE(d, st, pp, branch):
    """The report text (filled from the analysis outputs)."""
    return build_text(d, st, pp, branch)


if __name__ == "__main__":
    main()
