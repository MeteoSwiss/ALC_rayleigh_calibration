# -*- coding: utf-8 -*-
"""q19 - the Payerne CHM15k 2026-09-07 case (junk tier-2 rescue) through the vertical-consistency gate.

Operationally (daily runner with --classify) the Rayleigh night D = 2026-09-07 of 0-20000-0-06610_A is
calibrated with the Cloudnet target-classification pre-fit mask + the -11 veto: eprof_v2.2 returns flag
0.5, C_L = 1.0219e12 (2x the unit median 5.4e11), window 6064-6544 m ASL, u/C 42 %. Without the screen
the same night is flag -2.

This script regenerates the two classification curtains (D-1, D) exactly as the runner does
(scripts/run_network_calibration._do_classification, CAMS daily temperature), into the q19 folder
only, and then runs the night through q19_gradient_gate in every variant:
  op / fix (record)  : operational v2.2 / v2.2 with the measured night-mean noise
  gate_op / gate_fix : the same + the vertical-consistency gate on tier 2
each with and without the classification screen. It also writes a landscape figure of R(z).
Outputs only under C:/DATA/Projects/202610_CHM15k_stability/q19/.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "scripts"))

from common import OUT  # noqa: E402

Q19 = OUT / "q19"
CLS_ROOT = Q19 / "cls"
WMO, IDENT = "0-20000-0-06610", "A"
NIGHT = "20260907"
NEIGHBOURS = ["20260904", "20260905", "20260906", "20260910", "20260911", "20260912"]


def make_classification(days):
    """Write the per-day classification NetCDFs with the operational runner code (q19 folder only)."""
    os.environ["ALC_FULLCAL_DIR"] = str(CLS_ROOT)
    os.environ["ALC_CAMS_DIR"] = "D:/CAMS_daily"
    import run_network_calibration as RNC          # reads the env at import
    assert Path(RNC.OUT) == CLS_ROOT, RNC.OUT
    s = next(x for x in json.load(open(REPO / "validation/scope_l1_2026_census.json", encoding="utf-8"))
             if x["wmo"] == WMO and x["ident"] == IDENT)
    for d in days:
        dt = datetime.strptime(d, "%Y%m%d")
        n = RNC._do_classification(s, dt, dt, shared=None, force=False)
        print("classification", d, n, flush=True)


def main():
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    from datetime import timedelta
    days = sorted({(datetime.strptime(d, "%Y%m%d") + timedelta(days=k)).strftime("%Y%m%d")
                   for d in [NIGHT] + NEIGHBOURS for k in (-1, 0)})
    make_classification(days)

    import q19_gradient_gate as Q
    from calibration.io.classification import classification_files_for_night, read_classification_curtain
    variants = Q.load_gate_variants()
    gv = variants[0][1]["gate"]
    s = Q.census_entry(WMO, IDENT)
    key = f"{WMO}_{IDENT}"
    rows = []
    profiles = {}
    for d in [NIGHT] + NEIGHBOURS:
        native = Q.load_native(s, d)
        cls = read_classification_curtain(classification_files_for_night(CLS_ROOT, key, WMO, d))
        for screen in (True, False):
            runs = [("op", "record", False, None), ("fix", "record", True, None)]
            runs += [(name, "gate", spec.get("base", "op") == "fix", spec["gate"]) for name, spec in variants]
            for variant, mode, fix, gpar in runs:
                res, cap = Q.run_one(s, d, native, mode, gate=gpar, classification=cls if screen else None,
                                     cells=False, sigma_fix=fix)
                sel = cap.get("gsel") if mode == "gate" and cap.get("gsel") else cap.get("sel") or {}
                gd = cap.get("gate") or {}
                verdict = ""
                st = {}
                if mode == "record" and cap.get("sel") and cap.get("profile") and cap["profile"][3] is not None:
                    sig, p, z, sg = cap["profile"]
                    P = Q.prep_profile(sig, p, z, sg)
                    st = Q.candidate_stats(P, cap["sel"]["s"], cap["sel"]["e"], gv)
                    verdict = Q.gate_decision(st, gv)[1]
                rows.append(dict(date=d, screen=screen, variant=variant, flag=res["flag"], c=res["c"],
                                 u_rel=res["u"] / res["c"] if res["c"] == res["c"] else np.nan,
                                 zb_agl=res["zb"], zt_agl=res["zt"], tier=cap.get("tier"),
                                 verdict_selected=verdict, gate_verdict=gd.get("verdict", ""),
                                 g_up=sel.get("g_up"), sg_up=sel.get("sg_up"),
                                 lcb=st.get("lcb"), best_s=st.get("best_s"), n_blk=st.get("n_blk"),
                                 z_hi=sel.get("z_hi"), cls_ok=cls is not None,
                                 sigma_ratio=cap.get("sigma_ratio"), message=res["message"]))
                if variant == "op" and "profile" in cap:
                    profiles[(d, screen)] = (cap["profile"], cap.get("sel"))
                print({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in rows[-1].items()},
                      flush=True)
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(Q19 / "q19_payerne_20260907.csv", index=False)
    plot(profiles, df)


def plot(profiles, df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(16, 6.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    for a, (d, screen) in zip(axs, ((NIGHT, True), (NIGHT, False), ("20260905", True))):
        if (d, screen) not in profiles:
            continue
        (sig, p, z, sg), sel = profiles[(d, screen)]
        R = sig / p
        sR = sg / p
        ref = np.nanmedian(R[(z > 2000) & (z < 4000)]) if sel is None else sel["Rw"]
        nb = 10
        m = (z > 300) & (z < 10000) & np.isfinite(R)
        zz, RR, ss = z[m] / 1e3, R[m] / ref, sR[m] / ref
        k = len(zz) // nb
        zb_ = zz[:k * nb].reshape(k, nb).mean(1)
        Rb = RR[:k * nb].reshape(k, nb).mean(1)
        sb = np.sqrt((ss[:k * nb].reshape(k, nb) ** 2).sum(1)) / nb
        a.plot(RR, zz, ",", color="0.75")
        a.errorbar(Rb, zb_, xerr=sb, fmt="o-", ms=3, color="k", lw=0.8, elinewidth=0.7,
                   label="R(z)/R(window), 300 m blocks, measured noise")
        if sel:
            a.axhspan(sel["zs"] / 1e3, sel["ze"] / 1e3, color="tab:orange", alpha=0.35,
                      label="v2.2 window")
            a.axhspan(sel["z_lo"] / 1e3, sel["z_hi"] / 1e3, color="tab:blue", alpha=0.08,
                      label="extended band")
        a.axvline(1, color="tab:red", lw=0.8)
        a.set_xlim(0, 3)
        a.set_ylim(0.3, 10)
        r = df[(df.date == d) & (df.screen == screen) & (df.variant == "op")]
        tt = f"Payerne CHM15k {d} ({'with' if screen else 'without'} classification screen)"
        if len(r):
            r = r.iloc[0]
            tt += f"\nv2.2 flag {r.flag:g}, C={r.c:.3g}, tier {r.tier}"
            if r.verdict_selected:
                tt += f"\ng_up={r.g_up:+.1f}+/-{r.sg_up:.1f} %/km, gate verdict: {r.verdict_selected}"
        a.set_title(tt, fontsize=9)
        a.set_xlabel("R(z) = signal / p_mol, normalised to the window")
        a.grid(alpha=0.3)
    axs[0].set_ylabel("range above the instrument [km]")
    axs[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(Q19 / "q19_payerne_20260907.png", dpi=120)


if __name__ == "__main__":
    main()
