# -*- coding: utf-8 -*-
"""CL61 Rayleigh-vs-cloud agreement, v2.0 against v2.2, across every CL61 in the corpus.

The CL61 carries two independent retrievals of the SAME constant: the Rayleigh fit (clear nights)
and the O'Connor liquid-cloud integration (cloudy days). They should agree; how far the Rayleigh
value sits from the same-station cloud Kalman is therefore an absolute per-night error meter that
window games cannot fool -- the cloud method never sees the molecular window.

Reference = the fullcal cloud Kalman (legacy MS-eta tables, pre-2026-07). That vintage biases the
ABSOLUTE ratio (current PVC tables lower the cloud constant ~19% at Payerne, lifting Ray/cloud
accordingly), but it divides out of every v2.0-vs-v2.2 and kept-vs-recovered comparison, which is
what this script is for.

Out:  <DATA>/cl61_ray_vs_cloud.json  +  figs/cl61_ray_vs_cloud.png  + printed table
"""
from __future__ import annotations
import csv
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rayleigh_availability"))

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
FULLCAL = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/fullcal_l1_2026")
MAN = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
FIG = REPO / "rayleigh_availability" / "figs"
VALID = (1.0, 0.5)


def cloud_kalman(key):
    """{date: cloud kalman value} from the fullcal archive."""
    f = FULLCAL / key / f"{key}_kalman.csv"
    if not f.exists():
        return {}
    out = {}
    with open(f, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["method"] == "cloud":
                try:
                    out[row["date"]] = float(row["kalman"])
                except ValueError:
                    pass
    return out


def match(dstr, ref, tol_days=5):
    """Cloud reference on this date, else nearest within tol."""
    if dstr in ref:
        return ref[dstr]
    d0 = datetime.strptime(dstr, "%Y%m%d")
    best, bv = None, None
    for k in range(1, tol_days + 1):
        for dd in (d0 - timedelta(days=k), d0 + timedelta(days=k)):
            v = ref.get(dd.strftime("%Y%m%d"))
            if v is not None:
                return v
    return None


def med(v):
    v = [x for x in v if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else np.nan


def main():
    rows = []
    for inst in MAN:
        if inst["group"] != "CL61":
            continue
        lab, key, alt = inst["label"], f"{inst['wmo']}_{inst['ident']}", float(inst["alt"])
        base = json.loads((DATA / "baselines" / f"base_eprof_v2_{lab}.json").read_text())
        cand = json.loads((DATA / "candidates" / f"cand_N2.5_{lab}.json").read_text())
        cloud = cloud_kalman(key)
        if not cloud:
            print(f"{lab}: no cloud series"); continue

        def valid(rec):
            return {d: v for d, v in rec.items() if v[0] in VALID and v[1]}

        bv, cv = valid(base), valid(cand)
        pop = {}
        for tag, days in (("v2.0", bv), ("kept", {d: cv[d] for d in cv if d in bv}),
                          ("recovered", {d: cv[d] for d in cv if d not in bv})):
            m = []
            for d, v in days.items():
                cw = match(d, cloud)
                if cw and cw > 0:
                    mid = (0.5 * (v[3] + v[4]) - alt) if (v[3] and v[4]) else np.nan
                    m.append((d, v[1] / cw, mid))
            pop[tag] = m
        if len(pop["v2.0"]) < 3 and len(pop["recovered"]) < 3:
            print(f"{lab}: too few matched nights"); continue

        # in-population height slope of ln(ratio-to-cloud): the origin discriminator
        allm = pop["kept"] + pop["recovered"]
        slope = np.nan
        if len(allm) >= 6:
            z = np.array([m[2] for m in allm], float) / 1000.0
            r = np.array([m[1] for m in allm], float)
            ok = np.isfinite(z) & np.isfinite(r) & (r > 0)
            if ok.sum() >= 6 and np.ptp(z[ok]) > 0.5:
                slope = float(np.polyfit(z[ok], np.log(r[ok]), 1)[0] * 100)

        rows.append(dict(
            label=lab.replace("_CL61", "").replace("_A", "").replace("_B", "").replace("_C", ""),
            n20=len(pop["v2.0"]), nk=len(pop["kept"]), nr=len(pop["recovered"]),
            r20=med([m[1] for m in pop["v2.0"]]),
            rk=med([m[1] for m in pop["kept"]]),
            rr=med([m[1] for m in pop["recovered"]]),
            r22=med([m[1] for m in pop["kept"] + pop["recovered"]]),
            zk=med([m[2] for m in pop["kept"]]), zr=med([m[2] for m in pop["recovered"]]),
            slope=slope,
        ))

    print(f"\n{'station':22s} {'n20/nk/nr':>10s} {'v2.0/cld':>9s} {'v2.2/cld':>9s} "
          f"{'rec/cld':>8s} {'rec-v20':>8s} {'z_k':>5s} {'z_r':>5s} {'%/km':>6s}")
    for r in rows:
        d = (r["rr"] / r["r20"] - 1) * 100 if np.isfinite(r["rr"]) and np.isfinite(r["r20"]) else np.nan
        print(f"{r['label']:22s} {r['n20']:3d}/{r['nk']:3d}/{r['nr']:3d} "
              f"{r['r20']:9.3f} {r['r22']:9.3f} {r['rr']:8.3f} {d:+7.1f}% "
              f"{r['zk']/1000 if np.isfinite(r['zk']) else np.nan:5.1f} "
              f"{r['zr']/1000 if np.isfinite(r['zr']) else np.nan:5.1f} {r['slope']:+6.1f}")

    (DATA / "cl61_ray_vs_cloud.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")

    # pooled: recovered-vs-v2.0 shift and the height slope
    dd = [(r["rr"] / r["r20"] - 1) * 100 for r in rows
          if np.isfinite(r.get("rr", np.nan)) and np.isfinite(r.get("r20", np.nan))]
    sl = [r["slope"] for r in rows if np.isfinite(r.get("slope", np.nan))]
    print(f"\npooled: recovered-vs-v2.0 shift median {np.median(dd):+.1f}% over {len(dd)} stations; "
          f"height slope median {np.median(sl):+.1f} %/km over {len(sl)}")

    _figure(rows)


def _figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [r for r in rows if np.isfinite(r.get("r20", np.nan)) or np.isfinite(r.get("r22", np.nan))]
    if not rows:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 4.8),
                                   gridspec_kw=dict(width_ratios=[1.7, 1]))
    y = np.arange(len(rows))
    for r, yy in zip(rows, y):
        ax1.plot([r["r20"]], [yy], "o", color="#1f77b4", ms=7)
        ax1.plot([r["r22"]], [yy], "s", color="#d62728", ms=6)
        if np.isfinite(r.get("rr", np.nan)):
            ax1.plot([r["rr"]], [yy], "^", color="#ff7f0e", ms=6)
        lo = np.nanmin([r["r20"], r["r22"], r.get("rr", np.nan)])
        hi = np.nanmax([r["r20"], r["r22"], r.get("rr", np.nan)])
        ax1.plot([lo, hi], [yy, yy], "-", color="#bbb", lw=1, zorder=0)
    ax1.axvline(1.0, color="#444", lw=1, ls=":")
    ax1.set_yticks(y, [r["label"] for r in rows], fontsize=9)
    ax1.set_xlabel("median  $C_{Rayleigh}$ / $C_{cloud}^{Kalman}$   (legacy-η cloud reference)")
    ax1.plot([], [], "o", color="#1f77b4", label="v2.0 nights")
    ax1.plot([], [], "s", color="#d62728", label="v2.2 all")
    ax1.plot([], [], "^", color="#ff7f0e", label="v2.2 recovered only")
    ax1.legend(fontsize=8, loc="lower right")
    ax1.grid(alpha=0.3, axis="x")
    ax1.set_title("CL61 Rayleigh vs same-station cloud calibration")

    ok = [r for r in rows if np.isfinite(r.get("slope", np.nan))
          and np.isfinite(r.get("zr", np.nan)) and np.isfinite(r.get("zk", np.nan))]
    for r in ok:
        dz = (r["zr"] - r["zk"]) / 1000.0
        d = (r["rr"] / r["r20"] - 1) * 100 if np.isfinite(r.get("rr", np.nan)) else np.nan
        pred = r["slope"] * dz
        ax2.plot(pred, d, "o", ms=7, color="#1f77b4")
        ax2.annotate(r["label"][:9], (pred, d), fontsize=7, xytext=(3, 3),
                     textcoords="offset points")
    lim = ax2.get_xlim() + ax2.get_ylim()
    lo, hi = min(lim), max(lim)
    ax2.plot([lo, hi], [lo, hi], ":", color="#d62728", lw=1, label="1:1")
    ax2.set_xlabel("height slope × Δz(rec−kept)  [%]")
    ax2.set_ylabel("recovered vs v2.0, rel. to cloud  [%]")
    ax2.set_title("Is the shift the height effect?")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.suptitle("v2.2 recovered nights against the height-independent cloud reference — all corpus CL61s",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "cl61_ray_vs_cloud.png"
    fig.savefig(out, dpi=140)
    print(f"figure -> {out}")


if __name__ == "__main__":
    main()
