# -*- coding: utf-8 -*-
"""Review attack: in-sample tuning. Do the design choices hold OUT OF SAMPLE?

The sub-window range (1 km below, 4 km floor, 1 km above), D = 1 km, the 1 % weight floor, the 2 %
aerosol allowance and the molecular climatology were chosen on the 36 units the doc reports. Here the
same candidates are evaluated on
  * a time hold-out of the design units (2025 nights vs 2026 nights), and
  * the FRESH units (27 CHM15k, 8 CL61, 2 Mini-MPL never used in the design),
with the exact shipped budget recomputed offline from the recorded wide-range sub-windows (identical
machinery to u23_analysis; the exact code value is also stored by the review harness).

Usage: python u23r_holdout.py [eprof_v2.2|eprof_v2.3]
"""
from __future__ import annotations

import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402
import u23_analysis as A  # noqa: E402
from u23r_truth import prepare  # noqa: E402

warnings.filterwarnings("ignore")

RANGES = {"window only": dict(below_m=0.0, above_m=0.0, min_alt_m=1500.0),
          "window + 1 km above": dict(below_m=0.0, above_m=1000.0, min_alt_m=1500.0),
          "+1 km below, floor 4 km (shipped)": dict(below_m=1000.0, above_m=1000.0, min_alt_m=4000.0),
          "+1 km below, floor 3 km": dict(below_m=1000.0, above_m=1000.0, min_alt_m=3000.0),
          "+1 km below, floor 1.5 km": dict(below_m=1000.0, above_m=1000.0, min_alt_m=1500.0)}
FLOORS = (0.005, 0.01, 0.02)


def _fit_worker(args):
    name, floor, recs = args
    cfg = dict(A.BASE_CFG, snr_min=3.0, floor=floor, **RANGES[name])
    out = np.full((len(recs), 4), np.nan)
    for i, (z, r, s, zb, zt) in enumerate(recs):
        out[i] = A.line_fit(z, r, s, zb, zt, cfg)
    return name, floor, out


def fits(df):
    recs = [(z, r, s, float(zb), float(zt)) for z, r, s, zb, zt in
            zip(df._z, df._r, df._s, df.zb_agl, df.zt_agl)]
    jobs = [(n, f, recs) for n in RANGES for f in FLOORS]
    res = {}
    with ProcessPoolExecutor(max_workers=24) as ex:
        for name, floor, arr in ex.map(_fit_worker, jobs):
            res[(name, floor)] = arr
    return res


def sig_total(df, arr, D, aer=None, mol="shipped"):
    g, gs, res = arr[:, 0], arr[:, 1], arr[:, 2]
    sv = A.vertical_from_fit(g, gs, res, D, "abs", True)
    sv = np.where(np.isfinite(sv), sv, 2.0 * df.sig_ens.values / 0.3 * D)   # module fallback
    s_aer = df.def_sig_aer.values if aer is None else aer * (
        df.type.map({"CHM15k": 1064.0, "CL61": 910.0, "Mini-MPL": 532.0}).values / 1064.0) ** 3
    if mol == "shipped":
        s_mol = df.def_sig_mol.values
    elif mol == "residual only":
        s_mol = np.full(len(df), 0.0135)
    else:
        s_mol = np.zeros(len(df))
    return np.sqrt(df.sig_ens.values ** 2 + sv ** 2 + np.nan_to_num(df.def_sig_noise.values) ** 2
                   + s_mol ** 2 + s_aer ** 2), np.sqrt(df.sig_ens.values ** 2 + sv ** 2
                                                       + np.nan_to_num(df.def_sig_noise.values) ** 2)


def main():
    import pandas as pd
    method = sys.argv[1] if len(sys.argv) > 1 else "eprof_v2.2"
    po = prepare("orig", method)
    pf = prepare("fresh", method)
    sets = {"orig": po["df"], "fresh": pf["df"]}
    tw = po["tw"]
    rows = []
    for set_name, df in sets.items():
        d = df[df.ok & np.isfinite(df.sig23)].copy().reset_index(drop=True)
        d["year"] = [int(x[:4]) for x in d.date]
        F = fits(d)
        chm = (d.type == "CHM15k") & (d.key != C.PAYERNE) & np.isfinite(d.anom)
        splits = {"all": chm, "2025": chm & (d.year == 2025), "2026": chm & (d.year == 2026)}
        if set_name == "fresh":
            splits["CL61"] = (d.type == "CL61") & np.isfinite(d.anom)
            splits["Mini-MPL"] = (d.type == "Mini-MPL") & np.isfinite(d.anom)
        keys = list(zip(d.key, d.date))
        for (name, floor), arr in F.items():
            for D in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
                sig, sig_nomol = sig_total(d, arr, D)
                for sp, m in splits.items():
                    st = C.pull_stats(d.anom[m] / sig[m])
                    r = dict(set=set_name, split=sp, design=name, floor=floor, D=D, n=st["n"],
                             ts_pull=st["rstd"], ts_out2=st["out2"],
                             med_sig=100 * np.median(sig[m]) if m.any() else np.nan)
                    if set_name == "orig" and sp in ("all", "2025", "2026"):
                        t = tw if sp == "all" else tw[[x[:4] == sp for x in tw.date]]
                        smap = dict(zip(keys, sig_nomol))
                        sa = np.array([smap.get((a, dd), np.nan) for a, dd in zip(t.ka, t.date)])
                        sb = np.array([smap.get((b, dd), np.nan) for b, dd in zip(t.kb, t.date)])
                        tp = C.pull_stats(t.anom.values / np.hypot(sa, sb))
                        r.update(twin_n=tp["n"], twin_pull=tp["rstd"], twin_out2=tp["out2"])
                    rows.append(r)
        # aerosol allowance and molecular term, shipped design (fresh and orig)
        arr = F[("+1 km below, floor 4 km (shipped)", 0.01)]
        for aer in (0.0, 0.01, 0.02, 0.03):
            for mol in ("shipped", "residual only", "none"):
                sig, _ = sig_total(d, arr, 1.0, aer=aer, mol=mol)
                for sp, m in splits.items():
                    st = C.pull_stats(d.anom[m] / sig[m])
                    rows.append(dict(set=set_name, split=sp, design="shipped range, D 1 km",
                                     floor=0.01, D=1.0, aerosol=aer, molecular=mol, n=st["n"],
                                     ts_pull=st["rstd"], ts_out2=st["out2"],
                                     med_sig=100 * np.median(sig[m]) if m.any() else np.nan))
        # exact code values of the shipped budget (no emulation)
        for sp, m in splits.items():
            st = C.pull_stats(d.anom[m] / d.sig23[m])
            rows.append(dict(set=set_name, split=sp, design="shipped (exact code value)", floor=0.01,
                             D=1.0, n=st["n"], ts_pull=st["rstd"], ts_out2=st["out2"],
                             v22_pull=C.pull_stats(d.anom[m] / d.sig_ens[m])["rstd"],
                             v22_out2=C.pull_stats(d.anom[m] / d.sig_ens[m])["out2"],
                             obs_sigma=100 * C.rstd(d.anom[m]),
                             med_sig=100 * np.median(d.sig23[m]) if m.any() else np.nan))
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / f"u23r_holdout_{method}.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)
    sel = out[(out.floor == 0.01) & (out.D == 1.0) & out.get("aerosol", np.nan).isna()]
    print(sel.round(3).to_string(index=False))
    print("\n--- D sensitivity, shipped range, weight floor 1 % ---")
    dd = out[(out.design == "+1 km below, floor 4 km (shipped)") & (out.floor == 0.01)
             & out.get("aerosol", np.nan).isna()]
    print(dd.pivot_table(index=["set", "split"], columns="D", values="ts_pull").round(3).to_string())
    print(dd[dd.set == "orig"].pivot_table(index="split", columns="D", values="twin_pull").round(3).to_string())
    print("\n--- weight floor, shipped range, D 1 km ---")
    ff = out[(out.design == "+1 km below, floor 4 km (shipped)") & (out.D == 1.0)
             & out.get("aerosol", np.nan).isna()]
    print(ff.pivot_table(index=["set", "split"], columns="floor", values="ts_pull").round(3).to_string())
    print("\n--- aerosol / molecular, shipped design ---")
    aa = out[out.design == "shipped range, D 1 km"]
    print(aa.pivot_table(index=["set", "split", "molecular"], columns="aerosol",
                         values="ts_pull").round(3).to_string())


if __name__ == "__main__":
    main()
