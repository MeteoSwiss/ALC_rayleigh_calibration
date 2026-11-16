# -*- coding: utf-8 -*-
"""q19 - exact constants for the REVISED rule variants (require_above=False) with minimal L1 reads.

The first-design variants (fix_vc, fix_vc_k3, op_vc) were run exactly on every tier-2 night. The revised
rule differs only for windows with no measurable block above them. The offline best-first simulation over
the stored relaxed cells reproduces the exact gate-run selections 3689/3689, so the selection changes
under the revised rule are known in advance (q19_rerun_nights.json, 122 nights). Only those nights are
re-run (each L1 night read once, all revised variants in memory); for every other tier-2 night the
revised variant selects the same window as its first-design counterpart, so its exact result is reused.
Writes gate/<unit>/revised.csv (variants fix_vc2, fix_vc2_k3, op_vc2).
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_gradient_gate as Q  # noqa: E402

PAIRS = {"fix_vc2": "fix_vc", "fix_vc2_k3": "fix_vc_k3", "op_vc2": "op_vc"}


def main():
    warnings.filterwarnings("ignore")
    V = dict(Q.load_gate_variants())
    variants = [(n, V[n]) for n in PAIRS]
    need = [tuple(x) for x in json.load(open(Q.Q19 / "q19_rerun_nights.json"))]
    rr = Q.Q19 / "gate_rerun"
    tasks = []
    by_unit = {}
    for wmo, ident, d in need:
        by_unit.setdefault((wmo, ident), []).append(d)
    for (wmo, ident), ds in by_unit.items():
        tasks.append((wmo, ident, sorted(ds), variants, rr / f"{wmo}_{ident}.csv"))
    print(f"{len(tasks)} units, {len(need)} nights to re-run", flush=True)
    Q._pool(Q.gate_chunk, tasks, 8)
    # merge: re-run rows where available, otherwise the first-design rows renamed
    fr = [pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str}) for f in rr.glob("*.csv")]
    rerun = pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()
    done = set(zip(rerun.wmo, rerun.ident, rerun.date)) if len(rerun) else set()
    old = []
    for f in sorted(Q.GATE_DIR.glob("*/chunk*.csv")):
        old.append(pd.read_csv(f, dtype={"date": str, "ident": str, "wmo": str}))
    old = pd.concat(old, ignore_index=True)
    out = []
    for new, src in PAIRS.items():
        o = old[old.variant == src].copy()
        o = o[[(w, i, d) not in done for w, i, d in zip(o.wmo, o.ident, o.date)]]
        o["variant"] = new
        out.append(o)
        if len(rerun):
            out.append(rerun[rerun.variant == new])
    M = pd.concat(out, ignore_index=True)
    for (wmo, ident), g in M.groupby(["wmo", "ident"]):
        g.to_csv(Q.GATE_DIR / f"{wmo}_{ident}" / "revised.csv", index=False)
    print("revised rows:", M.groupby("variant").size().to_dict())


if __name__ == "__main__":
    main()
