# -*- coding: utf-8 -*-
"""Amsterdam / Lindenberg : la soustraction du dark ESTIME ameliore-t-elle la comparaison ?

Compare deux runs Rayleigh locaux strictement identiques sauf ALC_DARK_PROFILE :
  * Amsterdam (0-20000-0-06240, 4 CHM15k co-localises) : dispersion inter-unites de C_L par
    nuit commune (CV robuste sur >= 3 unites), avec vs sans soustraction.
  * Lindenberg (0-20000-0-10393, CHM15k "0" vs CL61 "C") : stabilite du ratio C_0/C_C par
    nuit commune, avec vs sans.

Run : python rayleigh_availability/compare_dark_runs.py --nodark <dir> --dark <dir> --station amst|lind
Sortie : stats imprimees + figure doc/reports/figs_dark_clearsky/compare_<station>.png
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "doc" / "reports" / "figs_dark_clearsky"

STATIONS = {
    "amst": dict(wmo="0-20000-0-06240", idents=["A", "B", "C", "D"], mode="cv"),
    "lind": dict(wmo="0-20000-0-10393", idents=["0", "C"], mode="ratio"),
}


def read_success(outdir: Path, wmo: str, ident: str) -> dict[str, float]:
    """date -> cal_value des nuits Rayleigh reussies (flag 1/0.5)."""
    key = f"{wmo}_{ident}"
    csv = outdir / key / f"{key}_cal.csv"
    out = {}
    if not csv.exists():
        return out
    with open(csv, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip().split(",")
        idx = {c: header.index(c) for c in ("date", "method", "flag", "cal_value")}
        for line in fh:
            p = line.rstrip("\n").split(",")
            if len(p) <= idx["cal_value"] or p[idx["method"]] != "rayleigh":
                continue
            try:
                if float(p[idx["flag"]]) in (1.0, 0.5) and float(p[idx["cal_value"]]) > 0:
                    out[p[idx["date"]]] = float(p[idx["cal_value"]])
            except ValueError:
                continue
    return out


def robust_cv(vals: np.ndarray) -> float:
    med = np.median(vals)
    return float(1.4826 * np.median(np.abs(vals - med)) / med) if med > 0 else np.nan


def stats_cv(series: dict[str, dict[str, float]], idents) -> tuple[list[str], list[float]]:
    """Par nuit avec >= 3 unites : CV robuste inter-unites des series NORMALISEES.

    Chaque unite est d'abord divisee par sa propre mediane : les C_L absolus different
    legitimement entre unites (materiel), et un simple realignement de niveaux par la
    soustraction ne doit PAS compter comme une amelioration. Ce CV mesure la divergence
    des co-fluctuations nuit-a-nuit -- ce que des instruments co-localises DOIVENT partager."""
    med = {i: np.median(list(series[i].values())) for i in idents if series[i]}
    norm = {i: {d: v / med[i] for d, v in series[i].items()}
            for i in idents if i in med and med[i] > 0}
    dates = sorted(set().union(*[set(s) for s in norm.values()]))
    out_d, out_v = [], []
    for d in dates:
        vals = np.array([norm[i][d] for i in norm if d in norm[i]])
        if vals.size >= 3:
            out_d.append(d)
            out_v.append(robust_cv(vals))
    return out_d, out_v


def stats_ratio(series, idents):
    """Par nuit commune aux 2 flux : ln(C_a / C_b)."""
    a, b = idents
    dates = sorted(set(series[a]) & set(series[b]))
    return dates, [float(np.log(series[a][d] / series[b][d])) for d in dates]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodark", required=True)
    ap.add_argument("--dark", required=True)
    ap.add_argument("--station", choices=list(STATIONS), required=True)
    args = ap.parse_args()
    spec = STATIONS[args.station]
    wmo, idents = spec["wmo"], spec["idents"]

    res = {}
    for lab, outdir in (("sans dark", Path(args.nodark)), ("avec dark", Path(args.dark))):
        series = {i: read_success(outdir, wmo, i) for i in idents}
        n_per = {i: len(series[i]) for i in idents}
        if spec["mode"] == "cv":
            dates, vals = stats_cv(series, idents)
            metric = "CV robuste inter-unites"
        else:
            dates, vals = stats_ratio(series, idents)
            metric = f"ln(C_{idents[0]}/C_{idents[1]})"
        res[lab] = (dates, np.array(vals))
        v = np.array(vals)
        print(f"[{lab}] nuits/flux: {n_per} ; {len(dates)} nuits comparables")
        if v.size:
            if spec["mode"] == "cv":
                print(f"   {metric}: mediane = {np.median(v) * 100:.2f} %  "
                      f"p75 = {np.percentile(v, 75) * 100:.2f} %")
            else:
                print(f"   {metric}: mediane = {np.median(v):+.4f} "
                      f"(ratio {np.exp(np.median(v)):.3f})  "
                      f"disp. robuste = {1.4826 * np.median(np.abs(v - np.median(v))):.4f}")

    # test apparie sur les nuits communes aux deux runs
    d0 = dict(zip(*res["sans dark"])) if res["sans dark"][0] else {}
    d1 = dict(zip(*res["avec dark"])) if res["avec dark"][0] else {}
    common = sorted(set(d0) & set(d1))
    if common:
        a = np.array([d0[d] for d in common])
        b = np.array([d1[d] for d in common])
        if spec["mode"] == "cv":
            delta = (np.median(b) - np.median(a)) * 100
            better = int(np.sum(b < a))
            print(f"\nAPPARIE ({len(common)} nuits) : delta CV median = {delta:+.2f} points de %"
                  f" ; le dark reduit le CV sur {better}/{len(common)} nuits")
        else:
            print(f"\nAPPARIE ({len(common)} nuits) : |ln ratio| median "
                  f"{np.median(np.abs(a)):.4f} -> {np.median(np.abs(b)):.4f} ; "
                  f"dispersion {1.4826 * np.median(np.abs(a - np.median(a))):.4f} -> "
                  f"{1.4826 * np.median(np.abs(b - np.median(b))):.4f}")

    # figure (paysage) : serie temporelle de la metrique, les deux runs superposes
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(13.5, 5.2))
    for lab, col in (("sans dark", "#1f77b4"), ("avec dark", "#d62728")):
        dates, vals = res[lab]
        if len(dates):
            t = [datetime.strptime(d, "%Y%m%d") for d in dates]
            y = np.array(vals) * (100 if spec["mode"] == "cv" else 1)
            ax.plot(t, y, ".", ms=4, alpha=0.6, color=col, label=lab)
            ax.axhline(np.median(y), color=col, ls="--", lw=1)
    ax.set_ylabel("CV inter-unites [%]" if spec["mode"] == "cv"
                  else f"ln(C_{idents[0]}/C_{idents[1]})")
    if spec["mode"] == "ratio":
        ax.axhline(0, color="#444", ls=":", lw=1)
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(f"{wmo} — effet de la soustraction du dark estime (nuits Rayleigh reussies)",
                 fontweight="bold")
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / f"compare_{args.station}.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
