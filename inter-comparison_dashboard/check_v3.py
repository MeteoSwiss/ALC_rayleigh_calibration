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
import re
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from validation.paper import intercompare as IC
import sites
import variants_v3 as V3

DATA = sites.DATA_ROOT
TARGET = 1064.0

# Mechanical gate (2026-08-16 review: the check used to have no threshold, no exit code, and a
# headline metric normalised by a near-zero denominator).  Healthy worst on the gated sweep is a
# fraction of a point (quantisation jitter); a real perturbation of the wiring moves rows by
# several points, so these absolute tolerances separate the two cleanly.
TOL_PP = 0.5      # percentage points, on medrelbias_pct / relbias_pct
# r_log is fragile on the CL31 with both corrections off (near-zero in-band signal: positive-pair
# selection flips on tiny numerical differences; healthy worst 0.12 on that one corner while
# medrel matches to 0.001 pt) — the tolerance sits above that but far below a wiring swap.
TOL_R = 0.15

# L1 rows of the dark-run variants are INFORMATIONAL, not gated: the v2 snapshot predates the
# 2026-08-16 dark refit + window extension of diag_v22_dark, so their deviation measures data
# vintage, not wiring (their L2 rows use no constants and stay gated — they match exactly).
STALE_CALS = {"v2.2dark", "v2.2dark_est"}

# v2 swept ONE variant dimension across all channels; v3 deliberately dropped the redundant
# cloud x pipeline-version corner (the cloud method does not depend on the molecular method), so
# only these cloud mappings are expressible from the v3 payload — the rest is skipped with a note.
CAL_MAP_CLOUD = {"v2.2": "cloudWV", "v2.2sansWV": "cloudNoWV"}
# v2's estimated-dark variant kept its own id; v3 collapsed the naming (the page now discloses
# the estimated provenance via dark_kind instead)
CAL_MAP_RAYLEIGH = {"v2.2dark_est": "v2.2dark"}
# for a 1064 nm channel the WV correction is an exact no-op, so v2's 'sansWV' variant is the
# ordinary v2.2 run — v3 dropped that physically meaningless entry (operator request), map it back
CAL_MAP_1064 = {"v2.2sansWV": "v2.2"}


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


def _combo_channels(site, cal, P):
    """v3 (ident, method, variant, dark_vec|None) list reproducing the v2 combo's channels, in the
    v2 channel order.  variant=None marks a channel the v3 payload cannot express for this cal
    (stats row skipped; the channel still feeds the paired mask via its canonical variant,
    which is exact because dividing by a positive constant never changes finiteness)."""
    base = {"payerne": [("A", "rayleigh"), ("B", "cloud"), ("C", "cloud"), ("C", "rayleigh")],
            "amsterdam": [(i, "rayleigh") for i in "ABCD"],
            "lindenberg": [("0", "rayleigh"), ("C", "cloud"), ("C", "rayleigh")]}[site]
    out = []
    for ident, method in base:
        if method == "rayleigh":
            variant = CAL_MAP_RAYLEIGH.get(cal, cal)
            corr = P["corr"][P["corr_of"][f"L1|{ident}"]]
            if corr.get("wv") is None:              # 1064 nm channel: WV variants collapse
                variant = CAL_MAP_1064.get(cal, variant)
        else:
            variant = CAL_MAP_CLOUD.get(cal)
        if variant is not None:
            rec = P["calib"].get(f"{ident}|{method}|{variant}")
            if not (rec and rec.get("ok")):
                variant = None
        # v2 subtracted the measured dark from every L1 channel under a dark variant
        db = P.get("dark", {}).get(ident) if cal in V3.DARK_RUNS else None
        out.append((ident, method, variant, db))
    return out


def main():
    site = sys.argv[1] if len(sys.argv) > 1 else "payerne"
    only = sys.argv[2] if len(sys.argv) > 2 else None
    P = json.loads((DATA / f"v3_{site}.json").read_text(encoding="utf-8"))
    V2 = json.loads((DATA / f"data_{site}.json").read_text(encoding="utf-8"))
    hday = np.asarray(P["hour_day"])
    b0, b1 = P["band"]
    band = slice(b0, b1 + 1)
    nz = len(P["z"])
    rng = "r0_" + str(len(V2["months"]) - 1)
    print(f"{P['name']}: {len(P['hours'])} paired hours, {nz} display gates, "
          f"band gates {b0}..{b1} ({b1-b0+1})")

    iref = 0
    zmask = np.zeros(nz, bool)
    zmask[band] = True
    worst_pp, worst_r, skipped = 0.0, 0.0, []
    combos = [c for c in V2["combos"] if only is None or c == only]
    print(f"sweeping {len(combos)} v2 combos against the v3 wiring "
          f"(gate: |d| <= {TOL_PP} pt / {TOL_R} on r_log)")
    for ckey in combos:
        m = re.match(r"(.+)_wv([01])_(molecular|none|angstrom)$", ckey)
        if not m:
            skipped.append((ckey, "unparseable combo key"))
            continue
        cal, wv, wl = m.group(1), m.group(2) == "1", m.group(3)
        v2combo = V2["combos"][ckey].get(rng)
        if not v2combo:
            skipped.append((ckey, "no full-range entry"))
            continue
        chans = _combo_channels(site, cal, P)
        if chans[iref][2] is None:
            skipped.append((ckey, "reference channel not expressible"))
            continue
        vals = {}
        for src in P["sources"]:
            for k, (ident, method, variant, db) in enumerate(chans):
                raw = dequant(P["streams"][f"{src}|{ident}"])
                corr = P["corr"][P["corr_of"][f"{src}|{ident}"]]
                # canonical variant stands in for an inexpressible one, ONLY for the mask
                vmask = variant or ("v2.2" if method == "rayleigh" else "cloudWV")
                rec = P["calib"].get(f"{ident}|{method}|{vmask}")
                if not (rec and rec.get("ok")):
                    rec = None
                vals[(src, k)] = transform(raw, rec["c"] if rec else np.ones(len(P["hours"])),
                                           corr, wv, wl, hday,
                                           dark=(db if src == "L1" else None),
                                           calibrated=(src == "L1" and rec is not None))
        have = np.logical_and.reduce([np.any(np.isfinite(vals[(s, k)][:, band]), axis=1)
                                      for s in P["sources"] for k in range(len(chans))])
        print(f"\n== {ckey}: rows {int(have.sum())} (v2 n_hours = {v2combo['n_hours']})")
        print(f"{'src':4s} {'channel':26s} {'medrelbias %':>26s} {'r_log':>16s} "
              f"{'relbias %':>22s} {'N':>14s}")
        for src in P["sources"]:
            for k, (ident, method, variant, db) in enumerate(chans):
                if k == iref:
                    continue
                if variant is None:
                    skipped.append((ckey, f"{src} {ident} {method}: no v3 record for '{cal}'"))
                    continue
                gated = not (src == "L1" and cal in STALE_CALS)
                new = IC._stats(vals[(src, k)][have], vals[(src, iref)][have], zmask)
                old = v2combo[src]["stats"][k]
                row = []
                for key, fmt in (("medrelbias_pct", "%+.4f"), ("r_log", "%.4f"),
                                 ("relbias_pct", "%+.4f"), ("n", "%d")):
                    o, n = old[key], new[key]
                    row.append(f"{fmt % n} / {fmt % o}")
                    if o is None or n is None or not gated:
                        continue
                    if key in ("medrelbias_pct", "relbias_pct"):
                        worst_pp = max(worst_pp, abs(n - o))          # absolute, in points
                    elif key == "r_log":
                        worst_r = max(worst_r, abs(n - o))
                lbl = f"{ident} {method} {variant}" + ("" if gated else "  [info]")
                print(f"{src:4s} {lbl:26s} {row[0]:>26s} {row[1]:>16s} {row[2]:>22s} "
                      f"{row[3]:>14s}")
    if skipped:
        print(f"\nnot expressible from the v3 payload ({len(skipped)}; by design, "
              f"see CAL_MAP_CLOUD):")
        for ckey, why in skipped:
            print(f"  - {ckey}: {why}")
    # ---- the approximation the spec allowed for the profile panel, measured but NOT used --------
    # Well-posed form: run the SAME transform chain but with C_L frozen at its window median,
    # i.e. "divide the window's median profile by the median C_L".  (Dividing an *uncalibrated*
    # transformed profile is not even defined once the molecular wavelength conversion adds its
    # additive term, which lives in calibrated units — that is one reason v3 does it per hour.)
    print("\nprofile panel: exact per-hour C_L division vs a frozen window-median C_L")
    mx = 0.0
    chans = _combo_channels(site, "v2.2", P)
    vals = {}
    for k, (ident, method, variant, db) in enumerate(chans):
        vals[k] = transform(dequant(P["streams"][f"L1|{ident}"]),
                            P["calib"][f"{ident}|{method}|{variant}"]["c"],
                            P["corr"][P["corr_of"][f"L1|{ident}"]], True, "molecular", hday)
    have = np.logical_and.reduce([np.any(np.isfinite(vals[k][:, band]), axis=1)
                                  for k in range(len(chans))])
    for k, (ident, method, variant, db) in enumerate(chans):
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

    # ---- the mechanical gate, LAST so every section above always prints -------------------------
    print(f"\nworst absolute deviation new-vs-v2: {worst_pp:.4f} pt "
          f"(medrelbias/relbias), {worst_r:.5f} (r_log)")
    if worst_pp > TOL_PP or worst_r > TOL_R:
        print(f"FAIL: deviation exceeds tolerance ({TOL_PP} pt / {TOL_R})")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
