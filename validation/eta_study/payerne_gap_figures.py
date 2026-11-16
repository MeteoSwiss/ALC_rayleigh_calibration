# -*- coding: utf-8 -*-
"""Figures of the Payerne CL61-vs-CHM15k gap dissection (report 04 v2, section 5.3-5.4).

  fig_payerne_gap_monthly.png   medrel of the CL61 channels vs the CHM15k (v2.2+dark, v2.2) per month,
                                June-October 2026, and the uncalibrated CHM15k/CL61 and CL31/CL61
                                signal ratios (ext/monthly_gap.csv, ext/raw_ratio_monthly.csv from
                                payerne_extension.py)
  fig_payerne_cl61_nights.png   per-night CL61 Rayleigh constants (hood dark v3 / no dark) against the
                                molecular-window centre, and the monthly medians of the CL61 and CHM15k
                                series (the calibration CSVs read by payerne_extension.series)
  fig_cl61_cloud_rayleigh_units.png   cloud / Rayleigh ratio per CL61 unit, operational vs proposed eta
                                (step4/ratio_cloud_rayleigh.csv), with the opaque-low-cloud dark estimate
                                of each unit (cl61_cloud_dark.log, remote_dark 2026-10-09)

Run:  python validation/eta_study/payerne_gap_figures.py [--out doc/reports/figs_04v2]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from validation.eta_study.payerne_extension import series, EXT, MONTHS   # noqa: E402

STUDY = Path("C:/DATA/Projects/202611_eta_study")
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
COL = {"raydark": "#1f77b4", "raynodark": "#7f7f7f", "cloud": "#404040", "cloud_prop": "#2ca02c",
       "chm": "#d62728", "cl31": "#ff7f0e"}


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def fig_monthly(out, plt):
    gp = EXT / "monthly_gap.csv"
    rp = EXT / "raw_ratio_monthly.csv"
    if not (gp.exists() and rp.exists()):
        print("monthly_gap.csv / raw_ratio_monthly.csv missing - run payerne_extension.py first")
        return
    gap = _read_csv(gp)
    rr = _read_csv(rp)
    curs = [("C Rayleigh+dark (v3)", "CL61 Rayleigh + hood dark", COL["raydark"]),
            ("C cloud (oper eta)", "CL61 cloud, operational eta", COL["cloud"]),
            ("C cloud (proposed eta, x0.90)", "CL61 cloud, proposed eta (x 0.90)", COL["cloud_prop"]),
            ("C Rayleigh no dark", "CL61 Rayleigh, no dark", COL["raynodark"])]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))
    x = np.arange(len(MONTHS))
    for ax, ref, title in ((axes[0], "A v2.2+dark (v3)", "vs CHM15k, v2.2 + hood dark (v3)"),
                           (axes[1], "A v2.2", "vs CHM15k, v2.2 (no dark)")):
        for i, (ck, lab, c) in enumerate(curs):
            y = [next((float(r["medrel"]) for r in gap if r["cur"] == ck and r["ref"] == ref and r["month"] == m
                       and r["subset"] == "all" and r["band"] == "0.5-3 km"), np.nan) for m in MONTHS]
            yn = [next((float(r["medrel"]) for r in gap if r["cur"] == ck and r["ref"] == ref and r["month"] == m
                        and r["subset"] == "night" and r["band"] == "0.5-3 km"), np.nan) for m in MONTHS]
            ax.bar(x + (i - 1.5) * 0.2, y, 0.18, color=c, label=lab)
            ax.plot(x + (i - 1.5) * 0.2, yn, "k_", ms=9, mew=1.4, label="night only" if i == 0 else None)
        ax.axhline(0, color="grey", lw=1)
        hrs = [next((r["hours"] for r in gap if r["ref"] == ref and r["month"] == m and r["subset"] == "all"
                     and r["band"] == "0.5-3 km"), "0") for m in MONTHS]
        ax.set_xticks(x)
        ax.set_xticklabels([f"{m[5:]}\n{h} h" for m, h in zip(MONTHS, hrs)])
        ax.set_ylabel("median relative difference over 500-3000 m [%]")
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3, axis="y")
        ax.set_ylim(-46 if ref == "A v2.2" else -32, 28)
        ax.legend(fontsize=8, loc="lower left")
    ax = axes[2]
    ac = np.array([float(r["AC_12_night"]) for r in rr])
    ac_d = np.array([float(r["AC_12_day"]) for r in rr])
    bc = np.array([float(r["BC_12_night"]) for r in rr])
    ref_i = MONTHS.index("2026-07")
    ax.plot(x, 100 * (ac / ac[ref_i] - 1), "o-", color=COL["chm"], lw=2, ms=8, label="CHM15k / CL61, night, 1-2 km")
    ax.plot(x, 100 * (ac_d / ac_d[ref_i] - 1), "o--", color=COL["chm"], lw=1.4, ms=6, alpha=0.6, label="CHM15k / CL61, day, 1-2 km")
    bc_rel = 100 * (bc / bc[ref_i] - 1)
    bc_rel[np.array(MONTHS) < "2026-07"] = np.nan        # CL31 optical block replaced 2026-07-07
    ax.plot(x, bc_rel, "s-", color=COL["cl31"], lw=2, ms=8, label="CL31 / CL61, night, 1-2 km (new CL31 optics from 2026-07-07)")
    ax.axhline(0, color="grey", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([m[5:] for m in MONTHS])
    ax.set_ylabel("uncalibrated signal ratio, change vs July [%]")
    ax.set_title("uncalibrated rcs_0 ratios, 1-2 km\n(hood dark subtracted, no constants)", fontsize=11)
    ax.set_ylim(-44, 8)
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("Payerne 2026: CL61 attenuated backscatter (910.55 nm line, molecular-aware conversion to 1064 nm)"
                 " against the co-located CHM15k, per month", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "fig_payerne_gap_monthly.png", dpi=140)
    plt.close(fig)
    print("->", out / "fig_payerne_gap_monthly.png")


def fig_nights(out, plt):
    from scipy.stats import theilslopes
    ser = series()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))
    ax = axes[0]
    for key, lab, c, mk in (("C_raydark", "Rayleigh + hood dark (v3)", COL["raydark"], "o"),
                            ("C_raynodark", "Rayleigh, no dark", COL["raynodark"], "^")):
        d = ser[key]
        zc = np.array([0.5 * (v[2] + v[3]) for v in d.values()]) / 1000
        C = np.array([v[0] for v in d.values()])
        ok = np.isfinite(zc) & (C > 0)
        sl = theilslopes(np.log(C[ok]), zc[ok])
        ax.plot(zc[ok], C[ok], mk, color=c, ms=5, alpha=0.7,
                label=f"{lab}: {100 * sl[0]:+.1f} %/km [{100 * sl[2]:+.1f}, {100 * sl[3]:+.1f}], n={ok.sum()}")
        zz = np.linspace(3, 7, 10)
        ax.plot(zz, np.exp(sl[1] + sl[0] * zz), "-", color=c, lw=2)
    ax.set_xlabel("molecular-window centre [km AGL]")
    ax.set_ylabel("CL61 Rayleigh constant C_L [Vaisala units]")
    ax.set_title("Payerne CL61, Feb-Oct 2026: the constant depends on where the window sits\n"
                 "(Theil-Sen slope of ln C_L; the hood dark removes most of the dependence)", fontsize=10.5)
    ax.set_ylim(0.3, 1.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    months = ["2026-02", "2026-03", "2026-04", "2026-05"] + MONTHS
    x = np.arange(len(months))

    def monthly(d, m):
        v = [val[0] for day, val in d.items() if day.startswith(m.replace("-", ""))]
        return (np.median(v), len(v)) if v else (np.nan, 0)

    ax = axes[1]
    for key, lab, c in (("C_raydark", "Rayleigh + hood dark (v3)", COL["raydark"]),
                        ("C_raynodark", "Rayleigh, no dark", COL["raynodark"]),
                        ("C_cloud", "cloud, operational eta", COL["cloud"])):
        y = [monthly(ser[key], m) for m in months]
        ax.plot(x, [v[0] for v in y], "o-", color=c, lw=2, ms=7, label=lab)
        for xx, (v, n) in zip(x, y):
            if n:
                ax.text(xx, v + 0.02, str(n), ha="center", fontsize=7, color=c)
    y = [monthly(ser["C_cloud"], m)[0] * 0.90 for m in months]
    ax.plot(x, y, "o--", color=COL["cloud_prop"], lw=1.6, ms=6, label="cloud, proposed eta (x 0.90)")
    ax.set_xticks(x)
    ax.set_xticklabels([m[5:] for m in months])
    ax.set_ylabel("monthly median of the per-night C_L [Vaisala units]")
    ax.set_title("CL61 monthly medians (number of nights/days above the points)", fontsize=10.5)
    ax.set_ylim(0.6, 1.3)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    ax = axes[2]
    for key, lab, c in (("A_v22dark", "CHM15k v2.2 + hood dark (v3)", COL["chm"]),
                        ("A_v22", "CHM15k v2.2, no dark", "#e8a0a0")):
        y = [monthly(ser[key], m) for m in months]
        ax.plot(x, [v[0] * 1e-11 for v in y], "o-", color=c, lw=2, ms=7, label=lab)
        for xx, (v, n) in zip(x, y):
            if n:
                ax.text(xx, v * 1e-11 + 0.1, str(n), ha="center", fontsize=7, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels([m[5:] for m in months])
    ax.set_ylabel("monthly median C_L [1e11 rcs units]")
    ax.set_title("CHM15k TUB140016 monthly medians\n(the distributed L2 used 5.5e11 in Jun-Jul and 5.9e11 in Aug)",
                 fontsize=10.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "fig_payerne_cl61_nights.png", dpi=140)
    plt.close(fig)
    print("->", out / "fig_payerne_cl61_nights.png")


def fig_units(out, plt):
    rows = _read_csv(STUDY / "step4" / "ratio_cloud_rayleigh.csv")
    darks = {}
    rx = re.compile(r"^(\S+)\s+(\S+).*dark/mol 3-6 km\s+([+-]?[\d.]+|\+nan)\s+\+/-\s+(\S+)\s+% \(raw\s+([+-]?[\d.]+|\+nan)")
    for line in (DATA / "cl61_cloud_dark.log").read_text(encoding="utf-8").splitlines():
        m = rx.match(line)
        if m:
            try:
                darks[m.group(1)] = float(m.group(5))
            except ValueError:
                pass
    labs = []
    for r in rows:
        if r["variant"] != "a":
            continue
        labs.append((r["key"], r["site"].title().replace("_", " "), r["reference"]))
    labs.sort(key=lambda t: float(next(r["ratio"] for r in rows if r["key"] == t[0] and r["reference"] == t[2] and r["variant"] == "a")))
    fig, ax = plt.subplots(figsize=(14, 6.2))
    y = np.arange(len(labs))
    for i, (key, site, ref) in enumerate(labs):
        for var, c, mk, lab in (("a", COL["cloud"], "o", "operational eta (PVC 5.5 um, 10 /km, S 18.8 sr)"),
                                ("e", COL["cloud_prop"], "D", "proposed (PVC 3.3 um, alpha_k, Mie S)")):
            r = next(r for r in rows if r["key"] == key and r["reference"] == ref and r["variant"] == var)
            v, lo, hi = float(r["ratio"]), float(r["ci_lo"]), float(r["ci_hi"])
            ax.errorbar(v, i + (0.12 if var == "a" else -0.12), xerr=[[v - lo], [hi - v]], fmt=mk, color=c, ms=7,
                        capsize=3, lw=1.2, label=lab if i == 0 else None)
        d = darks.get(key)
        tag = "hood dark corrected" if ref == "rayleigh+dark" else ("no dark" + (f", low-cloud dark {d:+.0f} %" if d is not None else ""))
        ax.text(1.42, i, tag, va="center", fontsize=8, color="#555555")
    ax.axvline(1, color="grey", lw=1)
    ax.axvspan(0.95, 1.05, color="grey", alpha=0.12, lw=0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{site}" for _, site, _ in labs])
    ax.set_xlabel("cloud C_L / Rayleigh C_L (ratio of medians, 95 % block-bootstrap CI)")
    ax.set_xlim(0.78, 1.72)
    ax.set_title("Liquid-cloud vs Rayleigh lidar constant of the CL61 units: the proposed eta shifts every unit by -9 %\n"
                 "(the Rayleigh reference is dark-free except Payerne; low-cloud dark = remote estimate at 3-6 km,"
                 " % of the molecular signal)", fontsize=10.5)
    ax.grid(alpha=0.3, axis="x")
    ax.legend(fontsize=9, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_cl61_cloud_rayleigh_units.png", dpi=140)
    plt.close(fig)
    print("->", out / "fig_cl61_cloud_rayleigh_units.png")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "doc" / "reports" / "figs_04v2"))
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    fig_nights(out, plt)
    fig_units(out, plt)
    fig_monthly(out, plt)


if __name__ == "__main__":
    main()
