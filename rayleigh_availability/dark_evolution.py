# -*- coding: utf-8 -*-
"""Payerne dark-measurement campaign (May-September 2026): per-session dark profiles of the CHM15k (A),
CL31 (B) and CL61 (C) and their evolution between sessions -- paper Appendix 5, Figure A5-2.

Per session (covered intervals verified in the L1 data, dark_windows_detected.json):
  * dark profile = time median of rcs_0 over the covered interval, smoothed in range as for the
    calibration (dark_profiles._smooth_range), converted to attenuated-backscatter units with the
    OPERATIONAL lidar constant of that day (dashboard Kalman: CHM15k Rayleigh, CL31/CL61 cloud), i.e.
    the backscatter the dark adds to the distributed product;
  * layer values in 0.5-1.5 km (boundary layer) and 2.5-6.5 km (Rayleigh fit range) expressed in % of
    the molecular attenuated backscatter at the instrument wavelength (US standard atmosphere, two-way
    molecular transmission, no water vapour), with the standard error from the per-profile scatter;
  * evolution: weighted mean, between-session scatter vs within-session error (reduced chi^2) and the
    linear trend per month with its standard error.
The CL31 optical module was replaced on 2026-07-07 (gain change), so its sessions before and after are
two different instruments and are analysed separately.

Run:  python rayleigh_availability/dark_evolution.py
"""
from __future__ import annotations
import csv
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dark_profiles as DP                                     # noqa: E402
from validation.paper import intercompare as IC                # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUTDIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v08")
CALIB = OUTDIR / "calib_dashboard"
PAY = "0-20000-0-06610"
STATION_ALT = 491.0
NAMES = {"A": "CHM15k", "B": "CL31", "C": "CL61"}
WL = {"A": 1064.47, "B": 909.7, "C": 910.55}
CALFILE = {"A": f"{PAY}_A_rayleigh.csv", "B": f"{PAY}_B_cloud.csv", "C": f"{PAY}_C_cloud.csv"}
LAYERS = {"BL": (500.0, 1500.0), "fit": (2500.0, 6500.0)}
SWAP = datetime(2026, 7, 7, 13, 0)


def kalman(ident):
    rows = [r for r in csv.DictReader(open(CALIB / CALFILE[ident], encoding="utf-8"))
            if r["C_kalman"] not in ("", "nan")]
    return (np.array([np.datetime64(r["time"]) for r in rows]),
            np.array([float(r["C_kalman"]) for r in rows]))


def c_at(dd, cc, t):
    return float(np.interp(np.datetime64(t.date()).astype("f8"), dd.astype("datetime64[D]").astype("f8"), cc))


def session(ident, w0, w1):
    t0 = datetime.strptime(w0, "%Y-%m-%d %H:%M:%S")
    t1 = datetime.strptime(w1, "%Y-%m-%d %H:%M:%S")
    blocks, rng = [], None
    day = t0.date()
    while day <= t1.date():
        got = DP._read_day(ident, datetime.combine(day, datetime.min.time()))
        if got is not None:
            times, r, x = got
            rng = r if rng is None else rng
            if r.size == rng.size:
                m = (times >= t0) & (times <= t1)
                if m.any():
                    blocks.append(x[m])
        day += timedelta(days=1)
    X = np.vstack(blocks)
    prof = DP._smooth_range(rng, np.nanmedian(X, axis=0))
    out = dict(start=t0, rng=rng, prof=prof, n=X.shape[0])
    for name, (lo, hi) in LAYERS.items():
        lay = (rng >= lo) & (rng <= hi)
        per = np.nanmean(X[:, lay], axis=1)                    # per-profile layer mean
        out[name] = float(np.nanmedian(per))
        out[name + "_se"] = float(1.2533 * np.nanstd(per) / np.sqrt(np.isfinite(per).sum()))
    return out


def main():
    win = json.loads((DP.DATA / "dark_windows_detected.json").read_text())
    res = {}
    for ident in ("A", "B", "C"):
        dd, cc = kalman(ident)
        rows = []
        for w0, w1 in win[ident]:
            s = session(ident, w0, w1)
            cl = c_at(dd, cc, s["start"])
            z = s["rng"]
            mol = IC._molecular_beta(z, STATION_ALT, WL[ident])          # Mm-1 sr-1, incl. T^2_mol
            s["beta"] = s["prof"] / cl * 1e6
            for name, (lo, hi) in LAYERS.items():
                lay = (z >= lo) & (z <= hi)
                m = float(np.nanmean(mol[lay]))
                s[name + "_beta"] = s[name] / cl * 1e6
                s[name + "_pct"] = 100.0 * s[name + "_beta"] / m
                s[name + "_pct_se"] = 100.0 * s[name + "_se"] / cl * 1e6 / m
            s["C_L"] = cl
            s["era"] = "post" if (ident == "B" and s["start"] > SWAP) else "pre" if ident == "B" else "all"
            rows.append(s)
        res[ident] = rows
    stats = summarize(res)
    (OUTDIR / "dark_evolution_stats.json").write_text(json.dumps(stats, indent=1, default=str))
    print(json.dumps(stats, indent=1, default=str))
    figure(res)


def summarize(res):
    out = {}
    for ident, rows in res.items():
        groups = {"all": rows} if ident != "B" else {
            "pre-swap": [r for r in rows if r["era"] == "pre"], "post-swap": [r for r in rows if r["era"] == "post"]}
        for g, rr in groups.items():
            key = NAMES[ident] + ("" if g == "all" else f" {g}")
            out[key] = {"n_sessions": len(rr), "dates": [r["start"].strftime("%Y-%m-%d") for r in rr]}
            for lay in LAYERS:
                v = np.array([r[lay + "_pct"] for r in rr])
                e = np.array([r[lay + "_pct_se"] for r in rr])
                w = 1 / np.maximum(e, 1e-12) ** 2
                mean = float(np.sum(w * v) / np.sum(w))
                chi2 = float(np.sum(w * (v - mean) ** 2) / max(len(v) - 1, 1))
                tt = np.array([(r["start"] - rr[0]["start"]).total_seconds() / 86400 / 30.44 for r in rr])
                if len(v) >= 3:
                    A = np.vstack([tt, np.ones_like(tt)]).T
                    coef, res_, *_ = np.linalg.lstsq(A, v, rcond=None)
                    resid = v - A @ coef
                    s2 = np.sum(resid ** 2) / max(len(v) - 2, 1)
                    se_slope = float(np.sqrt(s2 / np.sum((tt - tt.mean()) ** 2)))
                    slope = float(coef[0])
                else:
                    slope = se_slope = float("nan")
                out[key][lay] = dict(values_pct=[round(x, 2) for x in v], se_pct=[round(x, 2) for x in e],
                                     median_pct=float(np.median(v)), min_pct=float(v.min()),
                                     max_pct=float(v.max()), std_between_pct=float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                                     median_se_pct=float(np.median(e)), weighted_mean_pct=mean,
                                     reduced_chi2=chi2, trend_pct_per_month=slope,
                                     trend_se_pct_per_month=se_slope)
    return out


def figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.gridspec import GridSpec
    plt.rcParams.update({"font.size": 10})
    fig = plt.figure(figsize=(18, 10.5))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1.25, 1], hspace=0.38, wspace=0.75)
    d0 = min(r["start"] for rr in res.values() for r in rr)
    d1 = max(r["start"] for rr in res.values() for r in rr)
    cmap = plt.get_cmap("viridis")
    norm = lambda t: (t - d0).total_seconds() / max((d1 - d0).total_seconds(), 1)
    xlim = {"A": (-0.03, 0.01), "B": (-0.6, 0.6), "C": (-0.03, 0.01)}
    for k, ident in enumerate(("A", "B", "C")):
        ax = fig.add_subplot(gs[0, 2 * k:2 * k + 2])
        for r in res[ident]:
            m = (r["rng"] > 150) & (r["rng"] <= 8000)
            ls = "--" if r["era"] == "pre" else "-"
            ax.plot(r["beta"][m], r["rng"][m], ls, color=cmap(norm(r["start"])), lw=1.2,
                    label=r["start"].strftime("%d %b") + (" (pre-swap)" if r["era"] == "pre" else ""))
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlim(*xlim[ident])
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_ylim(0, 8000)
        ax.set_xlabel(r"dark / C$_L$  [Mm$^{-1}$ sr$^{-1}$]")
        ax.set_ylabel("Range [m]")
        ax.grid(alpha=0.3)
        lab = "abc"[k]
        title = f"({lab}) {NAMES[ident]} ({WL[ident]:.1f} nm): dark profile per session"
        if ident == "B":
            title += "\n(dashed: before the optical-module replacement of 7 July)"
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=2, loc="lower left")
    cols = {"A": "#d62728", "B": "#ff7f0e", "C": "#1f77b4"}
    # (d) boundary layer, all instruments; (e) fit range, CHM15k + CL61; (f) fit range, CL31 alone
    # (its dark is comparable to the molecular signal and would hide the other two on a common axis)
    panels = [("BL", ("A", "B", "C"), "d", None), ("fit", ("A", "C"), "e", (-45, 15)), ("fit", ("B",), "f", None)]
    for k, (lay, idents, letter, ylim) in enumerate(panels):
        lo, hi = LAYERS[lay]
        ax = fig.add_subplot(gs[1, 2 * k:2 * k + 2])
        for ident in idents:
            for era in ("pre", "post", "all"):
                rr = [r for r in res[ident] if r["era"] == era]
                if not rr:
                    continue
                t = [r["start"] for r in rr]
                v = [r[lay + "_pct"] for r in rr]
                e = [r[lay + "_pct_se"] for r in rr]
                lab = NAMES[ident] + {"pre": " (before 7 Jul)", "post": " (after 7 Jul)", "all": ""}[era]
                ax.errorbar(t, v, yerr=e, fmt="o" if era != "pre" else "s", ls="-" if era != "pre" else "--",
                            color=cols[ident], mfc=cols[ident] if era != "pre" else "white", ms=5, capsize=2,
                            lw=1, label=lab)
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(SWAP, color="0.5", ls=":", lw=1)
        ax.set_ylabel(r"dark / molecular $\beta_{att}$  [%]")
        which = "boundary layer" if lay == "BL" else "Rayleigh fit range"
        ax.set_title(f"({letter}) {' + '.join(NAMES[i] for i in idents)}: dark in {lo / 1000:.1f}–"
                     f"{hi / 1000:.1f} km ({which})", fontsize=10)
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        if ylim:
            ax.set_ylim(*ylim)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Payerne dark-measurement campaign, 12 May \u2013 29 September 2026 (covered telescope)",
                 fontsize=12.5, fontweight="bold", y=0.995)
    p = OUTDIR / "fig_v08_dark_evolution.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    print(f"-> {p}")


if __name__ == "__main__":
    main()
