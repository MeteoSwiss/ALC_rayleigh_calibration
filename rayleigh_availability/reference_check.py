# -*- coding: utf-8 -*-
"""Are the RECOVERED nights trustworthy? — inter-unit test on the Amsterdam quad.

The availability guards (canary / continuity / sigma_SD / Kalman outliers) can all be satisfied by
a change that adds *self-consistently wrong* nights. The only way to catch that is an independent
reference, and the cleanest one in the network needs no external instrument at all:

  0-20000-0-06240 A/B/C/D — FOUR co-located CHM15k at Amsterdam Schiphol.

Same site, same night, same atmosphere, same wavelength: no water-vapour or wavelength conversion
is involved, so any disagreement between the units is calibration error. Each unit has its own true
constant, so the statistic is the RATIO between a pair, normalised by that pair's own median over
the record. A perfect calibration gives ratio 1 every night; the scatter of the ratio is the joint
calibration error of the pair.

The test then asks the decisive question:

    is the pair-ratio scatter on nights that only v2.2 recovers WORSE than on the nights v2
    already had?

If the recovered nights are sound, the two scatters are comparable. If the recovery is admitting
junk (or altitude-shifted constants that no longer describe the same instrument), the recovered
nights are visibly noisier or offset.

Run:  python rayleigh_availability/reference_check.py [config]
"""
from __future__ import annotations
import itertools
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                       # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
QUAD_WMO = "0-20000-0-06240"
REF = "eprof_v2"
FIG = REPO / "rayleigh_availability" / "figs"


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def series(rec):
    """{date: C_L} for valid nights."""
    return {d: v[1] for d, v in rec.items()
            if IND.is_valid(v[0]) and v[1] and np.isfinite(v[1]) and v[1] > 0}


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "N1.5"
    quad = [i for i in MANIFEST if i["wmo"] == QUAD_WMO and i["group"] == "CHM15k"]
    if len(quad) < 2:
        print("Amsterdam quad not in the corpus"); return

    ref_s, new_s = {}, {}
    for i in quad:
        r = load(BASE / f"base_{REF}_{i['label']}.json")
        n = load(CAND / f"cand_{cfg}_{i['label']}.json")
        if r and n:
            ref_s[i["ident"]] = series(r)
            new_s[i["ident"]] = series(n)
    if len(ref_s) < 2:
        print("need at least two units with both runs"); return
    print(f"Amsterdam quad, config {cfg}: units {sorted(ref_s)}")
    print(f"  valid nights per unit  v2: "
          f"{ {k: len(v) for k, v in sorted(ref_s.items())} }")
    print(f"                       v2.2: {  {k: len(v) for k, v in sorted(new_s.items())} }\n")

    print(f"  {'pair':6s} {'n both-old':>10s} {'n any-new':>9s} | "
          f"{'scatter old':>12s} {'scatter new-nights':>18s} {'offset new':>11s}")
    rows = []
    for a, b in itertools.combinations(sorted(ref_s), 2):
        # nights where BOTH units were already valid under v2 -> the reference scatter
        old_nights = sorted(set(ref_s[a]) & set(ref_s[b]))
        # nights where both units are valid under v2.2 and AT LEAST ONE was recovered
        new_nights = sorted((set(new_s[a]) & set(new_s[b]))
                            - (set(ref_s[a]) & set(ref_s[b])))
        if len(old_nights) < 8 or len(new_nights) < 5:
            continue
        r_old = np.array([ref_s[a][d] / ref_s[b][d] for d in old_nights])
        r_new = np.array([new_s[a][d] / new_s[b][d] for d in new_nights])
        base = np.median(r_old)                       # the pair's own calibration offset
        # robust relative scatter of the ratio, normalised by the pair's median
        s_old = 1.4826 * np.median(np.abs(r_old - base)) / base * 100
        s_new = 1.4826 * np.median(np.abs(r_new - base)) / base * 100
        off = (np.median(r_new) / base - 1) * 100
        rows.append((f"{a}/{b}", len(old_nights), len(new_nights), s_old, s_new, off))
        print(f"  {a}/{b:4s} {len(old_nights):10d} {len(new_nights):9d} | "
              f"{s_old:11.1f}% {s_new:17.1f}% {off:+10.1f}%")

    if rows:
        so = np.median([r[3] for r in rows]); sn = np.median([r[4] for r in rows])
        of = np.median([abs(r[5]) for r in rows])
        print(f"\n  median over {len(rows)} pairs: scatter {so:.1f}% (v2 nights) -> "
              f"{sn:.1f}% (nights v2.2 adds), |offset| {of:.1f}%")
        print(f"  ratio new/old scatter = {sn/so:.2f}x   "
              f"[gate: <= 1.2x means the recovered nights are as trustworthy as the existing ones]")
        verdict = "PASS" if (sn <= 1.2 * so and of <= 5.0) else "CONCERN"
        print(f"  -> {verdict}")
        _figure(rows, cfg)


def _figure(rows, cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    lbl = [r[0] for r in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.bar(x - 0.2, [r[3] for r in rows], 0.4, label="nights v2 already had", color="#888")
    ax.bar(x + 0.2, [r[4] for r in rows], 0.4, label="nights v2.2 adds", color="#1f77b4")
    ax.set_xticks(x); ax.set_xticklabels(lbl)
    ax.set_xlabel("Amsterdam co-located CHM15k pair")
    ax.set_ylabel("pair-ratio scatter [%]")
    ax.set_title(f"Are recovered nights as consistent between co-located units?  ({cfg})")
    ax.legend(); ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    out = FIG / "phase2_amsterdam_interunit.png"
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
