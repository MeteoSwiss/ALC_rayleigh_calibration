# -*- coding: utf-8 -*-
"""q19 - the vertical-consistency statistic computed 'both ways': with the measured night-mean noise
(default) and with eprof_v2.2's propagated sigma_signal (sqrt(300 s/dt_native) too large).

Offline (stored profiles, q19_offline.py and q19_offline.py op): verdict of the main rule on the
operational first-choice tier-2 window, and the best-first re-selection over the relaxed cells
(C of a re-selected window estimated as C x ratio_med(new)/ratio_med(old)). Offsets vs the unit's
tier-1 level (local line). Writes q19_sigma_ways.csv.
"""
from __future__ import annotations

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q19 = L.Q.Q19
warnings.filterwarnings("ignore")


def evaluate(tag, base, G):
    W = pickle.load(open(Q19 / f"q19_blocks_windows{tag}.pkl", "rb"))
    C = pickle.load(open(Q19 / f"q19_blocks_cells{tag}.pkl", "rb"))
    Wd = {(w["key"], w["date"], w["variant"], w["which"]): w for w in W}
    t2 = pd.read_csv(Q19 / f"q19_design_tier2_{base}.csv", dtype={"date": str})
    first, kept, off = [], [], []
    for r in t2.itertuples():
        v = L.vc_verdict2(Wd.get((r.key, r.date, base, "sel")), k=G["k"], eps_tol=G["eps_tol"], tol_lo=G["tol_lo"])
        first.append(v)
        ok, dl = False, np.nan
        cc = C.get((r.key, r.date, base))
        if cc is not None:
            cells, bl = cc
            for i in np.argsort(-np.nan_to_num(cells["quality"].values, nan=-np.inf), kind="stable"):
                if L.vc_verdict2(bl[i], k=G["k"], eps_tol=G["eps_tol"], tol_lo=G["tol_lo"]) == "accept":
                    a, b = cells["ratio_med"].values[i], r.sel_ratio_med
                    dl = np.log(a / b) if (a > 0 and b > 0) else 0.0
                    ok = True
                    break
        kept.append(ok)
        off.append(r.off_loc + dl if ok else np.nan)
    first, kept, off = np.array(first), np.array(kept), np.array(off)
    fin = np.isfinite(t2.off_loc.values)
    o = off[fin & kept]
    return dict(sigma="measured" if tag == "" else "propagated (v2.2)", base=base, n=int(fin.sum()),
                first_accept=100 * np.mean(first[fin] == "accept"),
                first_aerosol=100 * np.mean(first[fin] == "aerosol"),
                first_tail=100 * np.mean(first[fin] == "tail"),
                kept_resel=100 * np.mean(kept[fin]), off_median=100 * np.median(o), off_mean=100 * np.mean(o),
                off_sd=100 * L.robust_sd(o), off_p90=100 * np.percentile(np.abs(o), 90))


def main():
    G = L.Q.load_gate_variants()[0][1]["gate"]
    rows = [evaluate(tag, base, G) for tag in ("", "_opsigma") for base in ("op", "fix")]
    R = pd.DataFrame(rows)
    R.to_csv(Q19 / "q19_sigma_ways.csv", index=False)
    pd.set_option("display.width", 200)
    print(R.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
