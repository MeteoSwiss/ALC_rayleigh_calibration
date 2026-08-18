# -*- coding: utf-8 -*-
"""Phase 0 — build the tuning/validation corpus for the Rayleigh availability work.

Why a new corpus: the C8 gate tuning (2026-06-20) used Feb-May 2026 only, so it never saw the
summer failure mode that the forensics exposed — flag -2 rejects 30-47 % of CLEAR CHM15k nights in
Jun-Aug versus 8-13 % in winter. This corpus therefore spans winter + shoulder + summer of BOTH
2025 and 2026, and deliberately over-samples the deficit tail (stations whose -2 rate is far above
the network median) while keeping healthy stations and the other instrument types as controls.

Stream selection
  * the 24 streams of validation/scope_l1_2026.json  (C8's corpus: 10 CHM15k, 10 CL61, 4 Mini-MPL)
    -> continuity with the previous tuning, and the non-CHM15k types guard against regressions
  * the worst CHM15k streams by -2 rate per CLEAR night, measured on the local v2 archive
  * the Amsterdam quad (0-20000-0-06240 A/B/C/D) -- four co-located CHM15k, the inter-unit
    reference for "did the newly admitted nights stay consistent?"

Split (declared here, enforced by the analysis):
  * tune     = streams marked "tune", 2025 nights ONLY
  * holdout  = every other stream, plus EVERY stream's 2026 nights
    i.e. both a station holdout and a year holdout, which is what the C8 review demanded before
    thresholds tuned on one period are treated as permanent.

Out: rayleigh_availability/scope_availability.json
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CENSUS = json.loads((REPO / "validation" / "scope_l1_2026_census.json").read_text())
BASE24 = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
# v2 per-night archive used ONLY to rank the deficit tail (its flags are a mixed code vintage at
# 20260601, so it is never used as a baseline -- Phase 0 re-runs those).
V2_ARCHIVE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/fullcal_l1_2026")
L1_ROOT = Path("D:/E-PROFILE_L1_2026")
OUT = REPO / "rayleigh_availability" / "scope_availability.json"

# Contiguous ranges covering winter + shoulder + summer. The real archive extent is ~all of 2025
# and 2026 to mid-July (the census first/last is only its scan window, not the data extent).
# 2025 runs to 31 December: an earlier cut at 30 September left Nov/Dec out of BOTH years, so the
# "winter" season was represented by Jan/Feb alone and every time series carried a 3-month hole
# that looked like a calibration failure rather than a corpus boundary.
PERIODS_2025 = [("20250101", "20251231")]
PERIODS_2026 = [("20260101", "20260713")]

N_TAIL = 8              # deficit-tail CHM15k streams to add
MIN_CLEAR = 50          # clear nights required before a stream can be ranked
AMSTERDAM = "0-20000-0-06240"


def minus2_rate(key):
    """(-2 count, clear-night count) for one stream from the v2 archive; clear = flag not in {0,-1}."""
    f = V2_ARCHIVE / key / f"{key}_cal.csv"
    if not f.exists():
        return 0, 0
    n2 = nclear = 0
    with open(f, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("method") != "rayleigh":
                continue
            try:
                fl = float(r["flag"])
            except (TypeError, ValueError):
                continue
            if fl in (0.0, -1.0):
                continue
            nclear += 1
            if fl == -2.0:
                n2 += 1
    return n2, nclear


def l1_years(wmo, ident):
    """Years for which this stream actually has L1 files on disk."""
    out = []
    for year in (2025, 2026):
        n = 0
        for mm in range(1, 13):
            d = L1_ROOT / wmo / str(year) / f"{mm:02d}"
            if d.is_dir():
                n += sum(1 for p in d.glob(f"L1_{wmo}_{ident}{year}*.nc"))
            if n:
                break
        if n:
            out.append(year)
    return out


def main():
    cen = {(c["wmo"], c["ident"]): c for c in CENSUS}
    chosen, why = {}, {}

    for b in BASE24:                                    # C8's corpus
        chosen[(b["wmo"], b["ident"])] = b
        why[(b["wmo"], b["ident"])] = "base24"

    ranked = []                                         # deficit tail
    for d in sorted(V2_ARCHIVE.glob("*_*")):
        if not d.is_dir():
            continue
        wmo, _, ident = d.name.rpartition("_")
        c = cen.get((wmo, ident))
        if not c or c["type"] != "CHM15k":
            continue
        n2, nclear = minus2_rate(d.name)
        if nclear >= MIN_CLEAR:
            ranked.append((100.0 * n2 / nclear, wmo, ident, c))
    ranked.sort(reverse=True)
    added = 0
    for rate, wmo, ident, c in ranked:
        if added >= N_TAIL:
            break
        if (wmo, ident) in chosen:
            continue
        chosen[(wmo, ident)] = c
        why[(wmo, ident)] = f"tail({rate:.0f}%)"
        added += 1

    for c in CENSUS:                                    # Amsterdam quad (inter-unit reference)
        if c["wmo"] == AMSTERDAM and c["type"] == "CHM15k":
            chosen.setdefault((c["wmo"], c["ident"]), c)
            why.setdefault((c["wmo"], c["ident"]), "amsterdam_quad")

    # Stratified tune/holdout: within each group, every other stream (sorted) goes to tune, so both
    # splits carry deficit-tail AND healthy streams of each instrument type.
    out = []
    by_group = {}
    for k, c in sorted(chosen.items()):
        by_group.setdefault(c["type"], []).append((k, c))
    for group, items in by_group.items():
        for i, ((wmo, ident), c) in enumerate(items):
            years = l1_years(wmo, ident)
            if not years:
                continue
            periods = ([p for p in PERIODS_2025 if 2025 in years]
                       + [p for p in PERIODS_2026 if 2026 in years])
            n2, nclear = minus2_rate(f"{wmo}_{ident}")
            out.append(dict(
                wmo=wmo, ident=ident, type=c["type"], site=c["site"],
                lat=c["lat"], lon=c["lon"], alt=c["alt"],
                label=f"{c['site'][:18]}_{c['type']}_{ident}".replace(" ", "_"),
                group=c["type"], split=("tune" if i % 2 == 0 else "holdout"),
                reason=why[(wmo, ident)], years=years, periods=periods,
                v2_minus2_pct_clear=(round(100.0 * n2 / nclear, 1) if nclear else None),
                v2_clear_nights=nclear))
    out.sort(key=lambda r: (r["group"], r["wmo"], r["ident"]))
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")

    print(f"-> {OUT}  ({len(out)} streams)")
    for g in sorted({r['group'] for r in out}):
        rows = [r for r in out if r["group"] == g]
        nt = sum(r["split"] == "tune" for r in rows)
        print(f"  {g:9s} {len(rows):3d} streams ({nt} tune / {len(rows)-nt} holdout)")
    print(f"\n  {'key':32s} {'site':20s} {'grp':8s} {'split':8s} {'reason':14s} {'-2%clr':>7s} yrs")
    for r in out:
        print(f"  {r['wmo']+'_'+r['ident']:32s} {r['site'][:20]:20s} {r['group']:8s} "
              f"{r['split']:8s} {r['reason']:14s} "
              f"{('%.0f' % r['v2_minus2_pct_clear']) if r['v2_minus2_pct_clear'] is not None else '-':>7s} "
              f"{','.join(str(y) for y in r['years'])}")


if __name__ == "__main__":
    main()
