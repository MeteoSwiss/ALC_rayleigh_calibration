# -*- coding: utf-8 -*-
"""Synthese du scan dark nuits-claires : theta par flux vs le detecteur FE de reference.

Pour chaque flux scanne (les *_diag.json de dark_clearsky/), assemble :
  * theta, sem_theta, theta_corr (calibration affine du selftest), split-half r ;
  * le gradient across-night FE ANNEE-MOIS du meme flux, recalcule du CSV reseau exactement
    comme en phase 4 par 4.2 (ln cal_value vs mi-hauteur de fenetre, demoyennage par annee-mois)
    -- c'est le detecteur de reference valide (plancher ~4 %/km) que theta doit au moins egaler.

La question a trancher : les 5 flux "dark fort" connus (Messina, Payerne, Montsec, Bern,
Twenthe) se separent-ils des controles sains sur l'axe theta ? Si oui la methode nuits-claires
est un ecran reseau utilisable ; sinon le FE reste le seul detecteur sans capot.

Run : python rayleigh_availability/dark_clearsky_summary.py
Sorties : <DATA>/dark_clearsky/summary.csv + doc/reports/figs_dark_clearsky/scan_discrimination.png
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_clearsky")
CSV_ROOT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04")
FIG = REPO / "doc" / "reports" / "figs_dark_clearsky"

# Les 5 flux "dark fort" + cas limites de la phase 4 (verite de terrain du test de separation)
STRONG = {"0-20000-0-00203_A", "0-20000-0-06610_A", "0-20008-0-MSA_A",
          "0-20008-0-BRN_A", "0-20000-0-06290_A"}
BORDER = {"0-528-0-06233_A", "0-20000-0-06260_A"}


def fe_gradient(key: str, start: str = "00000000") -> tuple[float, float, int]:
    """Gradient across-night ln(C) vs mi-hauteur de fenetre en effets fixes annee-mois,
    replique de la phase 4 par 4.2. Retourne (pente %/km, SE, n)."""
    csv = CSV_ROOT / key / f"{key}_cal.csv"
    if not csv.exists():
        return np.nan, np.nan, 0
    rows = []
    with open(csv, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip().split(",")
        idx = {c: header.index(c) for c in
               ("date", "method", "flag", "cal_value", "bottom_height", "top_height")}
        for line in fh:
            p = line.rstrip("\n").split(",")
            if len(p) <= idx["top_height"] or p[idx["method"]] != "rayleigh":
                continue
            try:
                if float(p[idx["flag"]]) not in (1.0, 0.5):
                    continue
                c = float(p[idx["cal_value"]])
                z = 0.5 * (float(p[idx["bottom_height"]]) + float(p[idx["top_height"]])) / 1000.0
                if c > 0 and np.isfinite(z) and p[idx["date"]] >= start:
                    rows.append((p[idx["date"]][:6], z, np.log(c)))
            except ValueError:
                continue
    if len(rows) < 15:
        return np.nan, np.nan, len(rows)
    ym = np.array([r[0] for r in rows])
    z = np.array([r[1] for r in rows])
    y = np.array([r[2] for r in rows])
    zd, yd = z.copy(), y.copy()
    kept = np.zeros(z.size, bool)
    for g in np.unique(ym):
        m = ym == g
        if m.sum() >= 3:
            zd[m] -= z[m].mean()
            yd[m] -= y[m].mean()
            kept |= m
    if kept.sum() < 15:
        return np.nan, np.nan, int(kept.sum())
    zd, yd = zd[kept], yd[kept]
    vz = np.sum(zd ** 2)
    if vz <= 0:
        return np.nan, np.nan, int(kept.sum())
    slope = np.sum(zd * yd) / vz
    resid = yd - slope * zd
    n = zd.size
    dof = max(n - 1 - len(np.unique(ym[kept])), 1)
    se = float(np.sqrt(np.sum(resid ** 2) / dof / vz))
    return 100.0 * float(slope), 100.0 * se, n


def main():
    rows = []
    for dj in sorted(DATA.glob("*_diag.json")):
        if "_noaer" in dj.name:            # tests de sensibilite, pas des flux
            continue
        d = json.load(open(dj, encoding="utf-8"))
        key = f"{d['wmo']}_{d['ident']}"
        fe, fe_se, fe_n = fe_gradient(key, start=d.get("start", "00000000"))
        st = d.get("selftest") or {}
        theta0 = (st.get("zero") or {}).get("theta_hat", np.nan)
        split_r = d.get("split_r", np.nan)
        n_used = d.get("n_used", 0)
        # Garde-fou : un plancher de faux positif |theta0| > 0.5, une reproductibilite
        # split-half < 0.3 ou < 30 nuits = detecteur NON fiable sur ce flux (a griser).
        reliable = (np.isfinite(theta0) and abs(theta0) <= 0.5
                    and np.isfinite(split_r) and split_r >= 0.3 and n_used >= 30)
        rows.append(dict(
            key=key, itype=d.get("itype", "?"),
            theta=d.get("theta", np.nan), sem_theta=d.get("sem_theta", np.nan),
            theta_corr=d.get("theta_corr", np.nan),
            theta0=theta0,
            lam=(st.get("capot") or {}).get("theta_hat", np.nan),
            split_r=split_r, n_used=n_used,
            reliable=bool(reliable),
            cv_A=d.get("cv_A", np.nan), bias_tpl=d.get("bias_tpl", np.nan),
            fe=fe, fe_se=fe_se, fe_n=fe_n,
            group=("dark_fort" if key in STRONG else
                   "limite" if key in BORDER else "controle"),
        ))
    if not rows:
        print("aucun diagnostic trouve dans", DATA)
        return

    out_csv = DATA / "summary.csv"
    cols = list(rows[0].keys())
    with open(out_csv, "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"-> {out_csv}  ({len(rows)} flux)")

    hdr = (f"{'flux':26s} {'type':7s} {'grp':9s} {'theta':>7s} {'+-':>5s} {'corr':>6s} "
           f"{'FE %/km':>8s} {'+-':>5s} {'r_sh':>5s} {'n':>4s}")
    print(hdr)
    for r in sorted(rows, key=lambda r: (r["itype"], r["group"], -abs(r["theta"] or 0))):
        print(f"{r['key']:26s} {r['itype']:7s} {r['group']:9s} "
              f"{r['theta']:+7.2f} {r['sem_theta']:5.2f} {r['theta_corr']:+6.2f} "
              f"{r['fe']:+8.2f} {r['fe_se']:5.2f} {r['split_r']:5.2f} {r['n_used']:4d}")

    # Figure de discrimination : theta vs FE, colore par groupe (paysage)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = {"dark_fort": "#d62728", "limite": "#e8871a", "controle": "#1f77b4"}
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))
    for pane, typ in ((0, "CHM15k"), (1, "CL61")):
        ax = axes[pane]
        sub = [r for r in rows if r["itype"] == typ or (typ == "CL61" and r["itype"] == "CL31")]
        for r in sub:
            mk = "s" if r["itype"] == "CL31" else "o"
            c = col[r["group"]] if r["reliable"] else "#bbbbbb"
            ax.errorbar(r["fe"], r["theta"], xerr=r["fe_se"], yerr=r["sem_theta"],
                        fmt=mk, ms=7, color=c, capsize=2, lw=1)
            ax.annotate(r["key"].split("_")[0].split("-")[-1], (r["fe"], r["theta"]),
                        fontsize=7, xytext=(4, 4), textcoords="offset points",
                        color="#333" if r["reliable"] else "#999")
        ax.axhline(0, color="#444", ls=":", lw=1)
        ax.axvline(0, color="#444", ls=":", lw=1)
        ax.set_xlabel("gradient FE annee-mois [%/km] (detecteur de reference)")
        ax.set_ylabel(r"$\theta$ nuits-claires [unites dark Payerne]")
        ax.set_title(f"{typ} — separation dark fort / controles")
        ax.grid(alpha=0.3)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=l)
               for l, c in col.items()]
    axes[0].legend(handles=handles, fontsize=9)
    fig.suptitle("Detecteur dark nuits-claires vs gradient FE — scan reseau",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / "scan_discrimination.png"
    fig.savefig(p, dpi=140)
    print(f"figure -> {p}")


if __name__ == "__main__":
    main()
