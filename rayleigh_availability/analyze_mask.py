# -*- coding: utf-8 -*-
"""Phase 3 verdict — what the target-classification pre-fit mask changes.

Phase 2 left one open mechanism: the nights v2.2 recovers are fitted 1.4-1.9 km HIGHER than the
ones v2 already had, and C_L is not independent of fit height, so the recovered nights inherit
(gradient x height shift). The obvious candidate cause is residual aerosol in the LOW windows --
steady elevated layers that the temporal MAD screen cannot see, because a layer present all night
defines each altitude's median instead of standing out from it.

If that is the cause, masking those cells with an independent classification should:
  (a) make low windows eligible again on recovered nights  -> the height shift shrinks
  (b) flatten dC_L/dz                                      -> the offset mechanism disappears
  (c) move the recovered nights TOWARDS the co-located reference
If the gradient survives the mask, the cause is not aerosol and the Phase-2 caveat stands as
written (it would then be instrumental -- overlap/afterpulse -- and a two-pass correction is the
only route).

Compares, on the classified streams only and on nights where a classification exists:
    v2  vs  Mv2   the mask alone, on the CURRENT operational gates
    N2.5 vs M2.5  the mask added to the noise-aware gates

Run:  python rayleigh_availability/analyze_mask.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                        # noqa: E402
from canaries import CANARIES                                                   # noqa: E402
from run_classification import SUBSET                                           # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
PAIRS = [("eprof_v2", "Mv2", "operational gates"), ("N2.5", "M2.5", "noise-aware gates")]


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def series(cfg, label):
    return load(BASE / f"base_{cfg}_{label}.json") if cfg.startswith("eprof") else \
        load(CAND / f"cand_{cfg}_{label}.json")


def common(a, b):
    """Restrict both series to the dates they share (the mask run skips unclassified nights)."""
    keys = set(a) & set(b)
    return {d: a[d] for d in keys}, {d: b[d] for d in keys}


def valid(rec):
    return {d: v for d, v in rec.items() if IND.is_valid(v[0]) and v[1]}


def mid_height(v):
    """Window mid-height AGL. Records store ASL heights; the offset cancels in a difference."""
    if v[3] is None or v[4] is None:
        return np.nan
    return 0.5 * (v[3] + v[4])


def summarise_pair(ref_cfg, new_cfg, insts):
    rows = []
    for inst in insts:
        a0, b0 = series(ref_cfg, inst["label"]), series(new_cfg, inst["label"])
        v2 = series("eprof_v2", inst["label"])
        if not a0 or not b0:
            continue
        a, b = common(a0, b0)
        va, vb = valid(a), valid(b)
        v2v = valid(v2 or {})
        both = sorted(set(va) & set(vb))
        # The mask is NOT a superset construction -- it changes the profile the fit sees, so it can
        # move a night the operational series already has. Split the change by that, because moving
        # existing nights is a far bigger deployment decision than adding new ones.
        keptn = [d for d in both if d in v2v]
        recn = [d for d in both if d not in v2v]

        def _dc(ds):
            return np.array([vb[d][1] / va[d][1] - 1.0 for d in ds]) * 100 if ds else np.array([])

        dc_k, dc_r = _dc(keptn), _dc(recn)
        dh = np.array([mid_height(vb[d]) - mid_height(va[d]) for d in both]) if both else np.array([])
        gained = sorted(set(vb) - set(va))
        # what the mask converted: the flag those nights carried under the reference config
        from_flag = {}
        for d in gained:
            from_flag[a[d][0]] = from_flag.get(a[d][0], 0) + 1
        rows.append(dict(
            label=inst["label"], role=SUBSET.get(inst["label"], ""), n=len(a),
            n_a=len(va), n_b=len(vb),
            lost=sorted(set(va) - set(vb)), gained=gained, from_flag=from_flag,
            med_dc_kept=float(np.median(dc_k)) if dc_k.size else np.nan,
            moved_kept=int(np.sum(np.abs(dc_k) > 1.0)), n_kept=len(keptn),
            med_dc_rec=float(np.median(dc_r)) if dc_r.size else np.nan, n_rec=len(recn),
            med_dh=float(np.median(dh)) if dh.size else np.nan,
            grad_a=IND.altitude_gradient(a)["slope_pct_per_km"],
            grad_b=IND.altitude_gradient(b)["slope_pct_per_km"],
            m11_a=sum(1 for v in a.values() if v[0] == -11),
            m11_b=sum(1 for v in b.values() if v[0] == -11)))
    return rows


def main():
    insts = [i for i in MANIFEST if i["label"] in SUBSET]
    have = [i for i in insts if (DATA / "classification" / i["label"]).is_dir()]
    print(f"classified streams: {len(have)}/{len(insts)}\n")

    for ref_cfg, new_cfg, what in PAIRS:
        rows = summarise_pair(ref_cfg, new_cfg, have)
        if not rows:
            print(f"== {ref_cfg} -> {new_cfg} ({what}): no output yet\n")
            continue
        print(f"== {ref_cfg} -> {new_cfg}: the classification mask on the {what} ==")
        print(f"{'stream':28s} {'nights':>6s} {'valid':>12s} {'lost':>5s} {'gain':>5s} "
              f"{'dC kept%':>9s} {'moved':>7s} {'dC rec%':>8s} {'dh m':>6s} "
              f"{'grad %/km':>18s} {'-11':>9s}")
        for r in rows:
            print(f"{r['label'][:28]:28s} {r['n']:6d} {r['n_a']:5d}->{r['n_b']:5d}  "
                  f"{len(r['lost']):5d} {len(r['gained']):5d} "
                  f"{r['med_dc_kept']:9.2f} {r['moved_kept']:3d}/{r['n_kept']:<3d} "
                  f"{r['med_dc_rec']:8.2f} {r['med_dh']:6.0f} "
                  f"{r['grad_a']:8.1f}->{r['grad_b']:8.1f} "
                  f"{r['m11_a']:4d}->{r['m11_b']:3d}")
        tot_a = sum(r["n_a"] for r in rows)
        tot_b = sum(r["n_b"] for r in rows)
        lost = sum(len(r["lost"]) for r in rows)
        gained = sum(len(r["gained"]) for r in rows)
        moved = sum(r["moved_kept"] for r in rows)
        n_kept = sum(r["n_kept"] for r in rows)
        dcs = [r["med_dc_kept"] for r in rows if np.isfinite(r["med_dc_kept"])]
        ga = np.array([r["grad_a"] for r in rows], float)
        gb = np.array([r["grad_b"] for r in rows], float)
        ok = np.isfinite(ga) & np.isfinite(gb)
        print(f"{'TOTAL':28s} {sum(r['n'] for r in rows):6d} {tot_a:5d}->{tot_b:5d}  "
              f"{lost:5d} {gained:5d} {np.median(dcs) if dcs else np.nan:9.2f} "
              f"{moved:3d}/{n_kept:<3d}")
        conv = {}
        for r in rows:
            for f, n in r["from_flag"].items():
                conv[f] = conv.get(f, 0) + n
        if conv:
            print("  nights gained, by the flag they carried without the mask: "
                  + ", ".join(f"flag {int(f)}={n}"
                              for f, n in sorted(conv.items(), key=lambda kv: -kv[1])))
        if ok.any():
            print(f"  median |dC_L/dz|: {np.median(np.abs(ga[ok])):.1f} -> "
                  f"{np.median(np.abs(gb[ok])):.1f} %/km   "
                  f"(flattened on {int(np.sum(np.abs(gb[ok]) < np.abs(ga[ok])))}/{int(ok.sum())} streams)")
        n11 = sum(r["m11_b"] for r in rows)
        print(f"  flag -11 fired on {n11} nights "
              f"({100.0 * n11 / max(sum(r['n'] for r in rows), 1):.1f}% of classified nights)")
        for wmo, ident, date, why in CANARIES:
            lab = next((i["label"] for i in have
                        if i["wmo"] == wmo and i["ident"] == ident), None)
            if lab:
                rec = series(new_cfg, lab) or {}
                v = rec.get(date)
                state = "MISSING" if v is None else ("REJECTED" if not IND.is_valid(v[0])
                                                     else f"ADMITTED C_L={v[1]:.3g}")
                print(f"  canary {lab} {date}: {state}  [{why}]")
        print()


if __name__ == "__main__":
    main()
