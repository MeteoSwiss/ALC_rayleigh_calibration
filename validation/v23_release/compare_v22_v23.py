# -*- coding: utf-8 -*-
"""Calibration release v2.3 against the v2.2 release: what changed, stream by stream and night by night.

Inputs (local copies of the CSCS release trees, CSVs + NetCDFs only):
  --v22  C:/DATA/Projects/202611_v23_release/calout_v22_rel   (alc_v22_release.sbatch, 2025-01-01..2026-08-13)
  --v23  C:/DATA/Projects/202611_v23_release/calout_v23_rel   (alc_v23_release_array.sbatch, 2025-01-01..archive end)
Per stream: <key>/<key>_cal.csv (date, method, flag, cal_value, uncertainty, n_profiles, bottom_height,
top_height, [version], message), <key>/<key>_kalman.csv (method, date, kalman, kalman_std), the yearly NetCDFs
<key>/<year>/ALC_calibration_<key><year>.nc (calibration_version, calibration_method, uncertainty components).

What is compared (common window = the v2.2 window unless --end is given):
  * completeness: streams and days in each tree
  * Rayleigh nights: identical (tier 1, bit-identical by design), changed (tier-2 window replaced), dropped,
    added (new tier-2 rescues, low-sun nights), per type and per month; the ratio v2.3/v2.2 on changed nights;
    window tops (CL61 cap); low-sun nights per station and their level against the dark nights
  * uncertainty: relative 2-sigma uncertainty on identical nights, v2.2 ensemble vs v2.3 budget, and the QC
  * Kalman: v2.3 / v2.2 on common dates per stream (weighting + tier-2 + eta)
  * cloud: day ratio v2.3 / v2.2 (eta tables v2 vs 2026-07) per type, day counts
  * NetCDF: version codes and the uncertainty components

Outputs: --out (default C:/DATA/Projects/202611_v23_release/compare/): tables (CSV), summary.txt, and the
figures in doc/reports/figs_v23_release/ (landscape).

Run:  python validation/v23_release/compare_v22_v23.py [--v22 DIR --v23 DIR --out DIR --end YYYYMMDD]
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
FIGS = REPO / "doc" / "reports" / "figs_v23_release"
SUCCESS = {1.0, 0.5, 0.25}
TYPE_COL = {"CHM15k": "#d62728", "CL61": "#2ca02c", "Mini-MPL": "#8c564b", "CL31": "#ff7f0e", "CL51": "#9467bd"}
LOW_SUN = 0.25


def census():
    c = json.loads((REPO / "validation" / "scope_l1_2026_census.json").read_text(encoding="utf-8"))
    rows = c if isinstance(c, list) else c.get("streams", [])
    return {f"{r['wmo']}_{r['ident']}": r for r in rows}


def read_cal(tree):
    """{key: DataFrame(date:int, method, flag:float, C, u, bot, top, version, msg)}"""
    out = {}
    for f in sorted(glob.glob(str(Path(tree) / "*" / "*_cal.csv"))):
        key = Path(f).parent.name
        try:
            d = pd.read_csv(f, dtype=str, keep_default_na=False)
        except Exception:
            continue
        if d.empty or "method" not in d:
            continue
        d["date"] = pd.to_numeric(d["date"], errors="coerce")
        d["flag"] = pd.to_numeric(d["flag"], errors="coerce")
        d["C"] = pd.to_numeric(d["cal_value"], errors="coerce")
        d["u"] = pd.to_numeric(d["uncertainty"], errors="coerce")
        d["bot"] = pd.to_numeric(d.get("bottom_height", ""), errors="coerce")
        d["top"] = pd.to_numeric(d.get("top_height", ""), errors="coerce")
        d["version"] = pd.to_numeric(d.get("version", ""), errors="coerce")
        d["msg"] = d.get("message", "")
        d = d.dropna(subset=["date"])
        d["date"] = d["date"].astype(int)
        out[key] = d[["date", "method", "flag", "C", "u", "bot", "top", "version", "msg"]]
    return out


def read_kalman(tree):
    out = {}
    for f in sorted(glob.glob(str(Path(tree) / "*" / "*_kalman.csv"))):
        key = Path(f).parent.name
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        if d.empty or "kalman" not in d:
            continue
        d["date"] = pd.to_numeric(d["date"], errors="coerce")
        out[key] = d.dropna(subset=["date"])
    return out


def nc_summary(tree):
    import netCDF4
    vers, comps, n = Counter(), Counter(), 0
    for f in glob.glob(str(Path(tree) / "*" / "*" / "*.nc")):
        try:
            ds = netCDF4.Dataset(f)
        except Exception:
            continue
        v = np.asarray(ds["calibration_version"][:]).astype(int) if "calibration_version" in ds.variables else np.array([])
        for x in v:
            vers[int(x)] += 1
        comps[tuple(sorted(x for x in ds.variables if x.startswith("lidar_constant_uncertainty_")))] += 1
        n += 1
        ds.close()
    return n, vers, comps


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser()
    ap.add_argument("--v22", default="C:/DATA/Projects/202611_v23_release/calout_v22_rel")
    ap.add_argument("--v23", default="C:/DATA/Projects/202611_v23_release/calout_v23_rel")
    ap.add_argument("--out", default="C:/DATA/Projects/202611_v23_release/compare")
    ap.add_argument("--end", default="20260813", help="end of the common window (v2.2 archive end)")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    cen = census()
    itype = lambda k: cen.get(k, {}).get("type", "?")
    site = lambda k: cen.get(k, {}).get("site", k)
    end = int(a.end)
    L = []
    say = L.append

    c22, c23 = read_cal(a.v22), read_cal(a.v23)
    k22, k23 = read_kalman(a.v22), read_kalman(a.v23)
    keys = sorted(set(c22) | set(c23))
    say(f"streams: v2.2 {len(c22)}, v2.3 {len(c23)}, union {len(keys)}; common window to {end}")
    missing23 = [k for k in c22 if k not in c23]
    missing22 = [k for k in c23 if k not in c22]
    say(f"  only in v2.2: {len(missing23)} {missing23[:8]}   only in v2.3: {len(missing22)} {missing22[:8]}")
    last23 = {k: int(d.date.max()) for k, d in c23.items()}
    beyond = sorted((k, v) for k, v in last23.items() if v > end)
    say(f"  v2.3 streams with days beyond {end}: {len(beyond)} -> " + ", ".join(f"{k}:{v}" for k, v in beyond[:12]))

    # ------------------------------------------------------------------ Rayleigh nights
    rows = []
    lowsun = []
    for k in keys:
        t = itype(k)
        d22 = c22.get(k)
        d23 = c23.get(k)
        r22 = d22[(d22.method == "rayleigh") & (d22.date <= end)].set_index("date") if d22 is not None else None
        r23 = d23[(d23.method == "rayleigh") & (d23.date <= end)].set_index("date") if d23 is not None else None
        if r22 is None or r23 is None:
            continue
        s22 = r22[r22.flag.isin(SUCCESS) & (r22.C > 0)]
        s23 = r23[r23.flag.isin(SUCCESS) & (r23.C > 0)]
        common = s22.index.intersection(s23.index)
        for dte in common:
            a_, b_ = s22.loc[dte], s23.loc[dte]
            same = abs(b_.C / a_.C - 1) < 1e-9
            rows.append(dict(key=k, itype=t, date=int(dte), kind="identical" if same else "changed",
                             C22=a_.C, C23=b_.C, u22=a_.u / a_.C, u23=b_.u / b_.C, top22=a_.top, top23=b_.top,
                             flag23=b_.flag))
        for dte in s22.index.difference(s23.index):
            a_ = s22.loc[dte]
            rows.append(dict(key=k, itype=t, date=int(dte), kind="dropped", C22=a_.C, C23=np.nan, u22=a_.u / a_.C,
                             u23=np.nan, top22=a_.top, top23=np.nan, flag23=r23.flag.get(dte, np.nan)))
        for dte in s23.index.difference(s22.index):
            b_ = s23.loc[dte]
            rows.append(dict(key=k, itype=t, date=int(dte), kind="lowsun" if b_.flag == LOW_SUN else "added",
                             C22=np.nan, C23=b_.C, u22=np.nan, u23=b_.u / b_.C, top22=np.nan, top23=b_.top,
                             flag23=b_.flag))
        # low-sun nights over the WHOLE v2.3 record (not only the common window), vs the dark-night level
        r23all = d23[(d23.method == "rayleigh") & d23.flag.isin(SUCCESS) & (d23.C > 0)]
        ls = r23all[r23all.flag == LOW_SUN]
        if len(ls):
            dark = r23all[r23all.flag != LOW_SUN]
            lowsun.append(dict(key=k, site=site(k), itype=t, lat=cen.get(k, {}).get("lat", np.nan), n_lowsun=len(ls),
                               n_dark=len(dark), months=",".join(sorted({str(x)[:6] for x in ls.date})),
                               level_ratio=float(np.median(ls.C) / np.median(dark.C)) if len(dark) else np.nan))
    R = pd.DataFrame(rows)
    R.to_csv(out / "rayleigh_nights.csv", index=False, float_format="%.6g")
    say("\n== Rayleigh nights in the common window, per type (counts) ==")
    piv = R.pivot_table(index="itype", columns="kind", values="date", aggfunc="count", fill_value=0)
    say(piv.to_string())
    for t in sorted(R.itype.unique()):
        s = R[(R.itype == t)]
        ch = s[s.kind == "changed"]
        n22 = int((s.kind.isin(["identical", "changed", "dropped"])).sum())
        n23 = int((s.kind.isin(["identical", "changed", "added", "lowsun"])).sum())
        msg = f"  {t}: successes v2.2 {n22} -> v2.3 {n23} ({100 * (n23 / n22 - 1) if n22 else np.nan:+.1f} %); "
        if len(ch):
            r = ch.C23 / ch.C22
            msg += f"changed nights ratio v2.3/v2.2 median {np.median(r):.3f} (P10-P90 {np.percentile(r, 10):.3f}-{np.percentile(r, 90):.3f}, n {len(ch)}); "
        if t == "CL61":
            msg += f"window top v2.2 median {np.nanmedian(s[s.kind != 'added'].top22):.0f} m, v2.3 {np.nanmedian(s[s.kind != 'dropped'].top23):.0f} m; "
        say(msg)
    # per-stream changes (for the appendix table)
    per = R.groupby(["key", "itype", "kind"]).size().unstack(fill_value=0)
    per["site"] = [site(k) for k in per.index.get_level_values(0)]
    per.to_csv(out / "rayleigh_nights_per_stream.csv")
    LS = pd.DataFrame(lowsun).sort_values("lat", ascending=False) if lowsun else pd.DataFrame()
    LS.to_csv(out / "lowsun_nights.csv", index=False, float_format="%.4g")
    say("\n== polar-summer low-sun nights (flag 0.25, whole v2.3 record) ==")
    say(LS.to_string(index=False) if len(LS) else "  none")

    # Payerne A
    pa = R[R.key == "0-20000-0-06610_A"]
    if len(pa):
        say(f"\nPayerne CHM15k A (tier 1 only in v2.3): v2.2 successes {int(pa.kind.isin(['identical', 'changed', 'dropped']).sum())}, "
            f"v2.3 {int(pa.kind.isin(['identical', 'changed', 'added']).sum())}")

    # ------------------------------------------------------------------ uncertainty on identical nights
    idn = R[R.kind == "identical"]
    say("\n== relative 2-sigma uncertainty on identical (tier-1) nights: v2.2 ensemble vs v2.3 budget ==")
    for t in sorted(idn.itype.unique()):
        s = idn[idn.itype == t]
        say(f"  {t}: n {len(s)}  v2.2 median {100 * s.u22.median():.1f} % (P10-P90 {100 * s.u22.quantile(.1):.1f}-{100 * s.u22.quantile(.9):.1f})  "
            f"v2.3 median {100 * s.u23.median():.1f} % (P10-P90 {100 * s.u23.quantile(.1):.1f}-{100 * s.u23.quantile(.9):.1f})")
    # QC rejections per flag in v2.3 vs v2.2 (Rayleigh rows, common window)
    say("\n== Rayleigh flags, common window (rows per flag) ==")
    fl = []
    for ver, cc in (("v2.2", c22), ("v2.3", c23)):
        cnt = Counter()
        for k, d in cc.items():
            r = d[(d.method == "rayleigh") & (d.date <= end)]
            cnt.update(r.flag.round(2).tolist())
        fl.append((ver, cnt))
    allflags = sorted(set(fl[0][1]) | set(fl[1][1]), key=lambda x: (np.isnan(x), x))
    say("  flag   " + "  ".join(f"{f:>7}" for f in allflags))
    for ver, cnt in fl:
        say(f"  {ver}  " + "  ".join(f"{cnt.get(f, 0):>7d}" for f in allflags))

    # ------------------------------------------------------------------ Kalman
    say("\n== Kalman best estimate v2.3 / v2.2 on common dates (per stream median), per type and method ==")
    krows = []
    for k in keys:
        if k not in k22 or k not in k23:
            continue
        for m in ("rayleigh", "cloud"):
            a_ = k22[k][k22[k].method == m].set_index("date").kalman
            b_ = k23[k][k23[k].method == m].set_index("date").kalman
            cm = a_.index.intersection(b_.index)
            cm = cm[cm <= end]
            if len(cm) < 10:
                continue
            r = (b_.loc[cm] / a_.loc[cm]).replace([np.inf, -np.inf], np.nan).dropna()
            if len(r):
                krows.append(dict(key=k, itype=itype(k), method=m, n=len(r), ratio_med=float(r.median()),
                                  ratio_p05=float(r.quantile(.05)), ratio_p95=float(r.quantile(.95))))
    K = pd.DataFrame(krows)
    K.to_csv(out / "kalman_ratio_per_stream.csv", index=False, float_format="%.5g")
    for (t, m), s in K.groupby(["itype", "method"]):
        say(f"  {t:8s} {m:8s}: streams {len(s):3d}  median of per-stream median ratio {s.ratio_med.median():.4f}  "
            f"(streams P10-P90 {s.ratio_med.quantile(.1):.4f}-{s.ratio_med.quantile(.9):.4f})")

    # ------------------------------------------------------------------ cloud days
    crows = []
    for k in keys:
        if k not in c22 or k not in c23:
            continue
        a_ = c22[k]
        b_ = c23[k]
        a_ = a_[(a_.method == "cloud") & (a_.date <= end)].set_index("date")
        b_ = b_[(b_.method == "cloud") & (b_.date <= end)].set_index("date")
        sa = a_[a_.flag.isin(SUCCESS | {0.25}) & (a_.C > 0)]
        sb = b_[b_.flag.isin(SUCCESS | {0.25}) & (b_.C > 0)]
        cm = sa.index.intersection(sb.index)
        for dte in cm:
            crows.append(dict(key=k, itype=itype(k), date=int(dte), ratio=sb.loc[dte].C / sa.loc[dte].C))
        crows.append(dict(key=k, itype=itype(k), date=-1, ratio=np.nan, n22=len(sa), n23=len(sb), ncommon=len(cm)))
    Cd = pd.DataFrame(crows)
    Cd.to_csv(out / "cloud_days.csv", index=False, float_format="%.6g")
    say("\n== cloud days in the common window: v2.3 / v2.2 constant ratio (eta tables v2 vs 2026-07) ==")
    for t in ("CL31", "CL51", "CL61"):
        s = Cd[(Cd.itype == t) & (Cd.date > 0)]
        cnt = Cd[(Cd.itype == t) & (Cd.date < 0)]
        if len(s):
            say(f"  {t}: common days {len(s)} (v2.2 {int(cnt.n22.sum())}, v2.3 {int(cnt.n23.sum())}); ratio median {s.ratio.median():.4f} "
                f"P5-P95 {s.ratio.quantile(.05):.4f}-{s.ratio.quantile(.95):.4f}")

    # ------------------------------------------------------------------ NetCDF
    say("\n== NetCDF ==")
    for ver, tree in (("v2.2", a.v22), ("v2.3", a.v23)):
        n, vers, comps = nc_summary(tree)
        say(f"  {ver}: {n} files; version codes {dict(sorted(vers.items()))}; uncertainty components {[(len(k), v) for k, v in comps.items()]}")

    # ------------------------------------------------------------------ figures
    # F1 availability per month per type (Rayleigh)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), sharey=False)
    for ax, t in zip(axes, ("CHM15k", "CL61", "Mini-MPL")):
        s = R[R.itype == t].copy()
        s["month"] = (s.date // 100).astype(int)
        m22 = s[s.kind.isin(["identical", "changed", "dropped"])].groupby("month").size()
        m23 = s[s.kind.isin(["identical", "changed", "added", "lowsun"])].groupby("month").size()
        months = sorted(set(m22.index) | set(m23.index))
        x = np.arange(len(months))
        ax.bar(x - 0.2, [m22.get(m, 0) for m in months], 0.4, color="#9e9e9e", label="v2.2")
        ax.bar(x + 0.2, [m23.get(m, 0) for m in months], 0.4, color=TYPE_COL[t], label="v2.3")
        ax.set_xticks(x[::2])
        ax.set_xticklabels([str(m) for m in months][::2], rotation=60, fontsize=8)
        ax.set_title(f"{t}: Rayleigh calibrations per month, common window", fontsize=10)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("calibrated nights (all streams)")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_availability_month.png", dpi=130)
    plt.close(fig)
    # F2 changed-night ratio + window tops (CL61) + uncertainty
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    ax = axes[0]
    bins = np.linspace(0.5, 1.5, 51)
    for t in ("CHM15k", "CL61", "Mini-MPL"):
        ch = R[(R.itype == t) & (R.kind == "changed")]
        if len(ch):
            ax.hist((ch.C23 / ch.C22).clip(0.5, 1.5), bins, histtype="step", lw=2, color=TYPE_COL[t],
                    label=f"{t}: n {len(ch)}, median {np.median(ch.C23 / ch.C22):.3f}")
    ax.axvline(1, color="grey", lw=1)
    ax.set_xlabel("C_L v2.3 / v2.2 on changed (tier-2) nights")
    ax.set_title("tier-2 nights the repaired fallback moved", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    s = R[(R.itype == "CL61")]
    ax.hist(s.top22.dropna() / 1e3, np.arange(2, 9, 0.25), histtype="step", lw=2, color="#9e9e9e", label=f"v2.2 (n {s.top22.notna().sum()})")
    ax.hist(s.top23.dropna() / 1e3, np.arange(2, 9, 0.25), histtype="step", lw=2, color=TYPE_COL["CL61"], label=f"v2.3 (n {s.top23.notna().sum()})")
    ax.set_xlabel("CL61 molecular-window top [km ASL]")
    ax.set_title("CL61 windows: tier-2 cap at 5 km AGL", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    for t in ("CHM15k", "CL61", "Mini-MPL"):
        s = idn[idn.itype == t]
        if len(s):
            ax.hist(100 * s.u22.clip(0, 80), np.arange(0, 82, 2), histtype="step", lw=1.2, ls="--", color=TYPE_COL[t])
            ax.hist(100 * s.u23.clip(0, 80), np.arange(0, 82, 2), histtype="step", lw=2, color=TYPE_COL[t],
                    label=f"{t}: v2.2 {100 * s.u22.median():.0f} % -> v2.3 {100 * s.u23.median():.0f} % (n {len(s)})")
    ax.set_xlabel("relative uncertainty of the nightly constant, 2 sigma [%] (dashed v2.2 ensemble, solid v2.3 budget)")
    ax.set_title("identical tier-1 nights", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_rayleigh_changes.png", dpi=130)
    plt.close(fig)
    # F3 cloud ratio + Kalman ratios + low-sun
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    ax = axes[0]
    for t in ("CL31", "CL51", "CL61"):
        s = Cd[(Cd.itype == t) & (Cd.date > 0)]
        if len(s):
            ax.hist(s.ratio.clip(0.85, 1.05), np.linspace(0.85, 1.05, 81), histtype="step", lw=2, color=TYPE_COL[t],
                    label=f"{t}: median {s.ratio.median():.4f} (n {len(s)})")
    ax.axvline(1, color="grey", lw=1)
    ax.set_xlabel("cloud C_L v2.3 / v2.2 per day")
    ax.set_title("eta tables v2 (a_G 3.3 um, alpha 20 /km) vs 2026-07", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    for (t, m), s in K.groupby(["itype", "method"]):
        ax.hist(s.ratio_med.clip(0.7, 1.3), np.linspace(0.7, 1.3, 61), histtype="step", lw=2 if m == "rayleigh" else 1.4,
                ls="-" if m == "rayleigh" else "--", color=TYPE_COL.get(t, "k"), label=f"{t} {m}: median {s.ratio_med.median():.3f} ({len(s)} streams)")
    ax.axvline(1, color="grey", lw=1)
    ax.set_xlabel("per-stream median of Kalman v2.3 / v2.2")
    ax.set_title("Kalman best estimate (weighting + tier 2 + eta)", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax = axes[2]
    if len(LS):
        y = np.arange(len(LS))
        ax.barh(y, LS.n_lowsun, color="#1f77b4")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{s[:18]} ({t})" for s, t in zip(LS.site, LS.itype)], fontsize=7)
        for yy, (n, r) in enumerate(zip(LS.n_lowsun, LS.level_ratio)):
            ax.text(n + 0.3, yy, f"{n} nights, level {r:.2f} x dark" if np.isfinite(r) else f"{n}", va="center", fontsize=7)
        ax.invert_yaxis()
    ax.set_xlabel("polar-summer low-sun nights (flag 0.25) in the v2.3 record")
    ax.set_title("nights that had no calibration before", fontsize=10)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_cloud_kalman_lowsun.png", dpi=130)
    plt.close(fig)
    txt = "\n".join(L)
    print(txt)
    (out / "summary.txt").write_text(txt, encoding="utf-8")
    print("figures ->", FIGS)


if __name__ == "__main__":
    main()
