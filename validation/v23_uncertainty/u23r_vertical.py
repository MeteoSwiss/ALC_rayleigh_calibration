# -*- coding: utf-8 -*-
"""Review: which vertical-term design is CONDITIONALLY calibrated (not only on average)?

For each candidate (sub-window range x D) the vertical term is recomputed offline from the recorded
wide-range sub-windows (u23_analysis.line_fit, the module's own fit) and combined with the other
components as stored by the review harness ("shipped": per-gate noise; "revised": the window-mean
noise and the molecular climatology at the window centre ASL, exact code values rev_*). Reported per
candidate, CHM15k (Payerne excluded): pooled time-series pull, the pull in the bottom / top decile
of the reported sigma and in the top octile of |g| (the time-series truth is an UPPER bound of the
calibration error, so a pull well below 1 there is a demonstrated over-estimate), the spread of
log(pull) over the deciles (0 = conditionally calibrated), and the twin pull.

Usage: python u23r_vertical.py <set> <method> [noise: shipped|revised]
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
RANGES = {"shipped: 1 km below (floor 4 km) + 1 km above": dict(below_m=1000.0, above_m=1000.0,
                                                                 min_alt_m=4000.0),
          "window + 1 km above": dict(below_m=0.0, above_m=1000.0, min_alt_m=4000.0),
          "window only": dict(below_m=0.0, above_m=0.0, min_alt_m=4000.0)}
DS = (0.5, 0.75, 1.0)


def _fit(args):
    name, recs, sub_sig = args
    cfg = dict(A.BASE_CFG, snr_min=3.0, floor=0.01, **RANGES[name])
    out = np.full((len(recs), 4), np.nan)
    for i, (z, r, s, zb, zt) in enumerate(recs):
        ss = s if sub_sig is None or sub_sig[i].size != z.size else np.where(
            np.isfinite(sub_sig[i]), sub_sig[i], s)
        out[i] = A.line_fit(z, r, ss, zb, zt, cfg)
    return name, out


def evaluate(d, tw, noise="shipped"):
    import pandas as pd
    recs = [(z, r, s, float(zb), float(zt)) for z, r, s, zb, zt in
            zip(d._z, d._r, d._s, d.zb_agl, d.zt_agl)]
    sub_sig = None
    if noise == "revised" and "_sub_sig_wm" in d:
        sub_sig = list(d._sub_sig_wm)
    with ProcessPoolExecutor(max_workers=8) as ex:
        fits = dict(ex.map(_fit, [(n, recs, sub_sig) for n in RANGES]))
    n_col = "def_sig_noise" if noise == "shipped" else "rev_sig_noise"
    m_col = "def_sig_mol" if noise == "shipped" else "rev_sig_mol"
    chm = ((d.type == "CHM15k") & (d.key != C.PAYERNE) & np.isfinite(d.anom)).values
    keys = list(zip(d.key, d.date))
    rows = []
    for name, arr in fits.items():
        for D in DS:
            sv = A.vertical_from_fit(arr[:, 0], arr[:, 1], arr[:, 2], D, "abs", True)
            sv = np.where(np.isfinite(sv), sv, 2.0 * d.sig_ens.values / 0.3 * D)
            base = (d.sig_ens.values ** 2 + np.nan_to_num(d[n_col].values) ** 2)
            sig = np.sqrt(base + sv ** 2 + np.nan_to_num(d[m_col].values) ** 2
                          + np.nan_to_num(d.def_sig_aer.values) ** 2)
            a, s = d.anom.values[chm], sig[chm]
            st = C.pull_stats(a / s)
            q = pd.qcut(s, 10, labels=False)
            dec = np.array([C.pull_stats(a[q == k] / s[q == k])["rstd"] for k in range(10)])
            gq = pd.qcut(np.abs(arr[chm, 0]), 8, labels=False, duplicates="drop")
            top_g = C.pull_stats(a[gq == gq.max()] / s[gq == gq.max()])["rstd"]
            r = dict(range=name, D=D, noise=noise, n=st["n"], ts_pull=st["rstd"], ts_out2=st["out2"],
                     pull_bottom_decile=dec[0], pull_top_decile=dec[-1], top_g_octile=top_g,
                     decile_log_sd=float(np.std(np.log(dec))), med_sig=100 * np.median(s),
                     qc06_rejected=float(np.mean(2 * sig[np.isfinite(d.anom.values)] > 0.6)))
            if tw is not None and len(tw):
                nomol = np.sqrt(base + sv ** 2)
                smap = dict(zip(keys, nomol))
                sa = np.array([smap.get((x, y), np.nan) for x, y in zip(tw.ka, tw.date)])
                sb = np.array([smap.get((x, y), np.nan) for x, y in zip(tw.kb, tw.date)])
                tp = C.pull_stats(tw.anom.values / np.hypot(sa, sb))
                r.update(twin_pull=tp["rstd"], twin_out2=tp["out2"])
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    import pandas as pd
    set_name = sys.argv[1] if len(sys.argv) > 1 else "orig"
    method = sys.argv[2] if len(sys.argv) > 2 else "eprof_v2.2"
    noise = sys.argv[3] if len(sys.argv) > 3 else "shipped"
    p = prepare(set_name, method)
    d = p["df"]
    d = d[d.ok & np.isfinite(d.sig23)].reset_index(drop=True)
    out = evaluate(d, p.get("tw"), noise)
    out.insert(0, "method", method)
    out.insert(0, "set", set_name)
    out.to_csv(C.OUT / f"u23r_vertical_{set_name}_{method}_{noise}.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
