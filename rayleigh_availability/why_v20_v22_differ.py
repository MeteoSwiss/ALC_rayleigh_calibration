# -*- coding: utf-8 -*-
"""Where the v2.0 -> v2.2 difference comes from, in one figure per site.

The question this answers: v2.2 is a strict SUPERSET of v2.0 -- on every night v2.0 calibrates,
v2.2 returns the bit-identical constant -- so the whole difference between the two versions is the
ADDED nights, and nothing else. This draws that decomposition directly:

  row 1  C_L(t): v2.0 nights (grey) vs v2.2-recovered nights (blue), same axis
  row 2  the fit window each night used (bottom..top), same time axis

Reading it: recovered nights cluster where v2.0 has NO nights at all (summer at Payerne), and they
sit systematically ~1-2 km higher, because on those nights no low window passes the gates -- the
air below is loaded. Whether that lifts or lowers C_L is the site's own C(z) signature.

Run:  python rayleigh_availability/why_v20_v22_differ.py
Out:  rayleigh_availability/figs/why_v20_v22_differ.png
"""
from __future__ import annotations
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

REPO = Path(__file__).resolve().parents[1]
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
RAW = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/inter-comparison_dashboard/calib_raw")
FIG = REPO / "rayleigh_availability" / "figs"
MAN = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
VALID = (1.0, 0.5)

# (title, v2.0 record, v2.2 record, station altitude m)
def corpus(lab):
    inst = next(i for i in MAN if i["label"] == lab)
    return (json.loads((DATA / "baselines" / f"base_eprof_v2_{lab}.json").read_text()),
            json.loads((DATA / "candidates" / f"cand_N2.5_{lab}.json").read_text()),
            float(inst["alt"]))


def payerne_cl61():
    return (json.loads((RAW / "payerne_C_v2.0.json").read_text()),
            json.loads((RAW / "payerne_C_v2.2.json").read_text()), 490.0)


SITES = [
    ("Payerne CHM15k", lambda: corpus("PAYERNE_CHM15k_A")),
    ("Payerne CL61", payerne_cl61),
    ("Lindenberg CHM15k", lambda: corpus("LINDENBERG_CHM15k_0")),
    ("Aosta CHM15k", lambda: corpus("SAINT-CHRISTOPHE_A_CHM15k_0")),
    ("Aosta CL61", lambda: corpus("SAINT-CHRISTOPHE_A_CL61_B")),
    ("SIRTA/Palaiseau CHM15k", lambda: corpus("PALAISEAU_CHM15k_B")),
]


def split(base, cand):
    """(kept, recovered) as {date: record}; kept = valid under BOTH (identical constants)."""
    bv = {d for d, v in base.items() if v[0] in VALID and v[1]}
    kept = {d: v for d, v in cand.items() if v[0] in VALID and v[1] and d in bv}
    rec = {d: v for d, v in cand.items() if v[0] in VALID and v[1] and d not in bv}
    return kept, rec


def dt(d):
    return datetime.strptime(d, "%Y%m%d")


def main():
    n = len(SITES)
    fig, axes = plt.subplots(2, n, figsize=(3.5 * n, 6.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.35, 1]))
    for k, (title, loader) in enumerate(SITES):
        try:
            base, cand, alt = loader()
        except (FileNotFoundError, StopIteration):
            axes[0, k].set_title(f"{title}\n(no data)", fontsize=9)
            continue
        kept, rec = split(base, cand)
        ax, axw = axes[0, k], axes[1, k]

        for grp, col, lab in ((kept, "#777777", "v2.0 (= v2.2, identique)"),
                              (rec, "#1f77b4", "ajoutees par v2.2")):
            if not grp:
                continue
            t = [dt(d) for d in sorted(grp)]
            c = [grp[d][1] for d in sorted(grp)]
            ax.plot(t, c, "o", ms=3.2, color=col, label=f"{lab} ({len(grp)})")
            # window bottom..top as a vertical bar per night, AGL km
            for d in sorted(grp):
                v = grp[d]
                if v[3] and v[4]:
                    axw.plot([dt(d), dt(d)], [(v[3] - alt) / 1000, (v[4] - alt) / 1000],
                             "-", color=col, lw=0.8, alpha=0.55)
        if kept:
            m = np.median([v[1] for v in kept.values()])
            ax.axhline(m, color="#d62728", lw=1.1, ls=":", zorder=1)
            if rec:
                mr = np.median([v[1] for v in rec.values()])
                ax.set_title(f"{title}\najoutees: {mr / m - 1:+.0%} vs v2.0", fontsize=9.5)
            else:
                ax.set_title(title, fontsize=9.5)
        else:
            ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.6, loc="upper left")
        axw.grid(alpha=0.3)
        axw.set_ylim(0, 8)
        axw.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
        axw.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%y"))
        # keep the C axis readable: scale from the v2.0 population, not the outliers
        vals = [v[1] for v in kept.values()] or [v[1] for v in rec.values()]
        if vals:
            lo, hi = np.percentile(vals, [2, 98])
            pad = 0.9 * (hi - lo) if hi > lo else 0.4 * hi
            ax.set_ylim(max(0, lo - pad), hi + pad)
    axes[0, 0].set_ylabel(r"$C_L$")
    axes[1, 0].set_ylabel("fenetre de fit\n[km AGL]")
    fig.suptitle("D'ou vient l'ecart v2.0 -> v2.2 : uniquement des nuits AJOUTEES "
                 "(les nuits communes sont bit-identiques)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "why_v20_v22_differ.png"
    fig.savefig(out, dpi=140)
    print(f"-> {out}")

    # the same decomposition as printed numbers
    print(f"\n{'site':26s} {'v2.0':>6s} {'ajout':>6s} {'C v2.0':>10s} {'C ajout':>10s} "
          f"{'ecart':>7s} {'z v2.0':>7s} {'z ajout':>8s}")
    for title, loader in SITES:
        try:
            base, cand, alt = loader()
        except (FileNotFoundError, StopIteration):
            continue
        kept, rec = split(base, cand)
        if not kept:
            continue
        ck = np.median([v[1] for v in kept.values()])
        zk = np.median([(v[3] + v[4]) / 2 - alt for v in kept.values() if v[3]]) / 1000
        if rec:
            cr = np.median([v[1] for v in rec.values()])
            zr = np.median([(v[3] + v[4]) / 2 - alt for v in rec.values() if v[3]]) / 1000
            print(f"{title:26s} {len(kept):6d} {len(rec):6d} {ck:10.3g} {cr:10.3g} "
                  f"{cr / ck - 1:+6.1%} {zk:6.1f}km {zr:7.1f}km")
        else:
            print(f"{title:26s} {len(kept):6d} {0:6d} {ck:10.3g} {'-':>10s} {'-':>7s} {zk:6.1f}km")


if __name__ == "__main__":
    main()
