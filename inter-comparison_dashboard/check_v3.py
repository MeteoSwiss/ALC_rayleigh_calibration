# -*- coding: utf-8 -*-
"""Non-regression + approximation check for the v3 payload.

Decodes v3_<site>.json exactly the way the browser does (base64 -> log-quantised uint16 -> the
affine transform chain) and recomputes the statistics table, then compares it to the v2 page's
precomputed numbers in data_<site>.json.  Any deviation here is a deviation the operator would see.

Also quantifies the approximation the v3 SPEC allowed but v3 does NOT use: dividing a window-median
profile by the window-median C_L instead of dividing per hour.

Run:  python inter-comparison_dashboard/check_v3.py [site] [combo]
"""
from __future__ import annotations
import base64
import json
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from validation.paper import intercompare as IC
import sites

DATA = sites.DATA_ROOT
TARGET = 1064.0


def dequant(enc):
    """base64 log-quantised uint16 -> float64 (n_h x n_z).  The JS `decodeStream` twin."""
    q = np.frombuffer(base64.b64decode(enc["b"]), "<u2").reshape(enc["nh"], enc["nz"]).astype("i4")
    lo = np.asarray(enc["lo"], "f8")[None, :]
    hi = np.asarray(enc["hi"], "f8")[None, :]
    sign = np.where(q >= 32768, -1.0, 1.0)
    qq = q & 32767
    with np.errstate(all="ignore"):
        mag = lo * np.power(hi / lo, (qq - 1) / 32766.0)
    out = sign * mag
    out[qq == 0] = 0.0
    out[q == 0] = np.nan
    return out


def defloat(blk):
    if blk is None:
        return None
    return np.frombuffer(base64.b64decode(blk["b"]), "<f4").reshape(blk["n"], blk["m"]).astype("f8")


def transform(raw, c, corr, wv, wl, hday, dark=None, calibrated=True):
    v = raw.copy()
    if dark is not None:
        v = v - np.asarray(dark, "f8")[None, :]
    if calibrated:
        v = v * 1e6 / np.asarray(c, "f8")[:, None]
    if wv and corr.get("wv") is not None:
        v = v / defloat(corr["wv"])[hday]
    need = corr.get("lam") is not None and abs(corr["lam"] - TARGET) >= 1.0
    if need and wl == "molecular" and corr.get("bmt") is not None:
        v = defloat(corr["bmt"])[hday] + (v - defloat(corr["bml"])[hday]) * corr["f"]
    elif need and wl == "angstrom":
        v = v * corr["f"]
    return v


def main():
    site = sys.argv[1] if len(sys.argv) > 1 else "payerne"
    P = json.loads((DATA / f"v3_{site}.json").read_text(encoding="utf-8"))
    V2 = json.loads((DATA / f"data_{site}.json").read_text(encoding="utf-8"))
    hday = np.asarray(P["hour_day"])
    b0, b1 = P["band"]
    band = slice(b0, b1 + 1)
    nz = len(P["z"])
    print(f"{P['name']}: {len(P['hours'])} paired hours, {nz} display gates, "
          f"band gates {b0}..{b1} ({b1-b0+1})")

    # v3 state == the v2 default combo: v2.2 constants, WV comparison ON, molecular wavelength.
    v2combo = V2["combos"]["v2.2_wv1_molecular"]["r0_" + str(len(V2["months"]) - 1)]
    # v3 channel -> (ident, method, variant) reproducing that combo
    sel = {"payerne": [("A", "rayleigh", "v2.2"), ("B", "cloud", "cloudWV"),
                       ("C", "cloud", "cloudWV"), ("C", "rayleigh", "v2.2")],
           "amsterdam": [(i, "rayleigh", "v2.2") for i in "ABCD"],
           "lindenberg": [("0", "rayleigh", "v2.2"), ("C", "cloud", "cloudWV"),
                          ("C", "rayleigh", "v2.2")]}[site]

    vals = {}
    for src in P["sources"]:
        for k, (ident, method, variant) in enumerate(sel):
            raw = dequant(P["streams"][f"{src}|{ident}"])
            corr = P["corr"][P["corr_of"][f"{src}|{ident}"]]
            rec = P["calib"][f"{ident}|{method}|{variant}"]
            vals[(src, k)] = transform(raw, rec["c"], corr, True, "molecular", hday,
                                       calibrated=(src == "L1"))

    # paired mask (v2 `have`): every channel finite somewhere in the band, in every source
    have = np.logical_and.reduce([np.any(np.isfinite(vals[(s, k)][:, band]), axis=1)
                                  for s in P["sources"] for k in range(len(sel))])
    print(f"paired rows recomputed: {int(have.sum())}  (v2 n_hours = {v2combo['n_hours']})")

    iref = 0
    zmask = np.zeros(nz, bool)
    zmask[band] = True
    worst = 0.0
    print(f"\n{'src':4s} {'channel':22s} {'medrelbias %':>26s} {'r_log':>16s} "
          f"{'relbias %':>22s} {'N':>14s}")
    for src in P["sources"]:
        for k, (ident, method, variant) in enumerate(sel):
            if k == iref:
                continue
            new = IC._stats(vals[(src, k)][have], vals[(src, iref)][have], zmask)
            old = v2combo[src]["stats"][k]
            row = []
            for key, fmt in (("medrelbias_pct", "%+.4f"), ("r_log", "%.4f"),
                             ("relbias_pct", "%+.4f"), ("n", "%d")):
                o, n = old[key], new[key]
                row.append(f"{fmt % n} / {fmt % o}")
                if key != "n" and o:
                    worst = max(worst, abs(n - o) / max(abs(o), 1e-9) * 100)
            lbl = f"{ident} {method} {variant}"
            print(f"{src:4s} {lbl:22s} {row[0]:>26s} {row[1]:>16s} {row[2]:>22s} {row[3]:>14s}")
    print(f"\nworst relative deviation new-vs-v2 on the tabulated statistics: {worst:.4f} %")

    # ---- the approximation the spec allowed for the profile panel, measured but NOT used --------
    # Well-posed form: run the SAME transform chain but with C_L frozen at its window median,
    # i.e. "divide the window's median profile by the median C_L".  (Dividing an *uncalibrated*
    # transformed profile is not even defined once the molecular wavelength conversion adds its
    # additive term, which lives in calibrated units — that is one reason v3 does it per hour.)
    print("\nprofile panel: exact per-hour C_L division vs a frozen window-median C_L")
    mx = 0.0
    for k, (ident, method, variant) in enumerate(sel):
        raw = dequant(P["streams"][f"L1|{ident}"])
        corr = P["corr"][P["corr_of"][f"L1|{ident}"]]
        c = np.asarray(P["calib"][f"{ident}|{method}|{variant}"]["c"], "f8")
        cfix = np.full_like(c, float(np.median(c[have])))
        exact = np.nanmedian(transform(raw, c, corr, True, "molecular", hday)[have], axis=0)
        approx = np.nanmedian(transform(raw, cfix, corr, True, "molecular", hday)[have], axis=0)
        m = np.isfinite(exact) & np.isfinite(approx) & (np.abs(exact) > 1e-3)
        d = np.abs(approx[m] - exact[m]) / np.abs(exact[m]) * 100
        if d.size:
            print(f"  {ident} {method:8s} {variant:11s} max {d.max():7.3f} %   "
                  f"median {np.median(d):6.3f} %   (C_L spread "
                  f"{100*(c[have].max()/c[have].min()-1):.1f} %)")
            mx = max(mx, d.max())
    print(f"  -> worst-case profile approximation error {mx:.2f} % — v3 does NOT use it")


if __name__ == "__main__":
    main()
