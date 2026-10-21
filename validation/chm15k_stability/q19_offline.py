# -*- coding: utf-8 -*-
"""q19 - offline block statistics for every recorded window (no recalibration).

For each night of the record run (both variants) the stored profile (time-cleaned night-mean signal,
p_mol, measured night-mean noise on the fit grid) is re-read and, for
  * the selected window (tier 1 or tier 2),
  * the window the noise-aware tier WOULD select on a tier-1 night (t2sel),
  * every relaxed (tier-2 eligible) cell of the tier-2 nights,
the per-block window excesses R_window / R_block - 1 and their noise sigma are computed (750 m blocks
from the window top to z_hi, and down to max(1.5 km, zb - 1 km) below the window). Any k / tolerance /
power rule can then be applied offline (q19_lib.block_verdict).
Writes q19_blocks_windows.pkl and q19_blocks_cells.pkl in the q19 folder.
"""
from __future__ import annotations

import pickle
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q19_lib as L  # noqa: E402

Q = L.Q
warnings.filterwarnings("ignore")


def block_arrays(P, s, e, block_m=750.0, min_gates=10):
    """Window level and per-block excesses (above and below the window) for one window."""
    st = Q.window_stats(P, s, e)
    if P is None or not np.isfinite(st["Rw"]):
        return None
    z, ok, R, wR = P["z"], P["ok"], P["R"], P["wR"]
    Rw, sRw, ze, zs, z_hi, z_lo = st["Rw"], st["sRw"], st["ze"], st["zs"], st["z_hi"], st["z_lo"]
    nb = max(1, int(round(block_m / P["dz"])))
    up = np.where(ok & (z > ze) & (z <= z_hi))[0]
    gu = [up[i:i + nb] for i in range(0, up.size, nb)]
    if len(gu) > 1 and gu[-1].size < min_gates:
        gu[-2] = np.concatenate([gu[-2], gu[-1]])
        gu = gu[:-1]
    lo = np.where(ok & (z >= z_lo) & (z < zs))[0][::-1]
    gl = [lo[i:i + nb] for i in range(0, lo.size, nb)]
    rows = []
    for side, gs in (("up", gu), ("lo", gl)):
        for g in gs:
            if g.size < min_gates:
                continue
            w = wR[g]
            sw = float(np.sum(w))
            if sw <= 0:
                continue
            Rb = float(np.sum(w * R[g]) / sw)
            if not Rb > 0:
                continue
            sRb = float(1.0 / np.sqrt(sw))
            r = Rw / Rb
            rows.append((side == "up", r - 1.0, r * np.hypot(sRw / Rw, sRb / Rb), float(np.mean(z[g]))))
    a = np.array(rows, float) if rows else np.zeros((0, 4))
    return dict(zs=zs, ze=ze, z_hi=z_hi, z_lo=z_lo, Rw=Rw, sRw=sRw, g_up=st["g_up"], sg_up=st["sg_up"],
                is_up=a[:, 0].astype(bool), exc=a[:, 1], sig=a[:, 2], zc=a[:, 3])


SIGMA_ROW = 3          # 3 = measured night-mean noise (default); 4 = eprof_v2.2's propagated sigma_signal


def unit_job(key, sigma_row=3):
    df = L.load_unit(key)
    outw, outc = [], {}
    for month, g in df.groupby(df.date.str[:6]):
        f = Q.REC_DIR / key / f"{month}.npz"
        if not f.exists():
            continue
        with np.load(f) as zf:
            files = set(zf.files)
            prof = {k[2:]: zf[k] for k in files if k.startswith("p_")}
            cells = {k[2:]: zf[k] for k in files if k.startswith("c_")}
        for r in g.itertuples():
            pa = prof.get(r.date)
            if pa is None:
                continue
            z, sig, p, sg = pa[0], pa[1], pa[2], pa[sigma_row]
            P = Q.prep_profile(sig, p, z, sg)
            for which, s_, e_ in (("sel", r.sel_s, r.sel_e), ("t2", r.t2_s, r.t2_e)):
                if not (np.isfinite(s_) and np.isfinite(e_)):
                    continue
                b = block_arrays(P, int(s_), int(e_))
                if b is not None:
                    outw.append(dict(key=key, date=r.date, variant=r.variant, which=which, **b))
            ck = f"{r.variant}_{r.date}"
            if ck in cells:
                c = pd.DataFrame(cells[ck].astype(float), columns=L.CELL_COLS)
                lst = []
                for cr in c.itertuples():
                    b = block_arrays(P, int(cr.s), int(cr.e))
                    lst.append(b)
                outc[(key, r.date, r.variant)] = (c, lst)
    return outw, outc


def main():
    """python q19_offline.py [op]  ('op' = use eprof_v2.2's propagated sigma instead of the measured one)."""
    use_op = len(sys.argv) > 1 and sys.argv[1] == "op"
    row = 4 if use_op else 3
    tag = "_opsigma" if use_op else ""
    keys = sorted(L.load().key.unique())
    t0 = time.time()
    W, C = [], {}
    with ProcessPoolExecutor(max_workers=12) as ex:
        for k, (w, c) in zip(keys, ex.map(unit_job, keys, [row] * len(keys))):
            W += w
            C.update(c)
            print(f"{k}: {len(w)} windows, {len(c)} cell sets ({time.time() - t0:.0f} s)", flush=True)
    pickle.dump(W, open(Q.Q19 / f"q19_blocks_windows{tag}.pkl", "wb"))
    pickle.dump(C, open(Q.Q19 / f"q19_blocks_cells{tag}.pkl", "wb"))


if __name__ == "__main__":
    main()
