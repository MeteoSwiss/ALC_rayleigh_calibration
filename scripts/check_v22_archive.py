#!/usr/bin/env python3
"""Acceptance gates for the v2.2 release archive, run BEFORE the operational cutover.

The flip (uncommenting the block in ops/config.sh) is allowed only when this exits 0. Each gate
answers one question that, left unchecked, has a way of destroying something irreversible:

  G1 vintage      Is the new tree homogeneously v2.2, one row per (stream, method, date), and is
                  the OLD tree still free of v2.2 rows?  A mixed archive cannot be rolled back.
  G2 continuity   Do the constants match the already-validated v2.2 archive everywhere the release
                  was not supposed to change them?  The ONLY intended change is the CL61 water-
                  vapour spectrum, so any moved constant on another type is a regression.
                  Includes the canaries: two Payerne nights that must stay rejected.
  G3 completeness Does every stream carry every artefact the operational server needs -- including
                  both .npz caches, whose absence makes the daily runner skip OmB/sensitivity
                  FOREVER rather than rebuild them, and the per-year NetCDF in its operational
                  place rather than under parts_nc/.
  G4 cloudcover   Are the new status columns present and physical (0..8 octas)?

Gates that cannot be evaluated (a tree not supplied) are reported SKIPPED, never passed.

Run (on balfrin, where the trees live):
  python scripts/check_v22_archive.py \
      --new /scratch/mch/mhrvo/E_PROFILE_calout_v22_rel \
      --ref /scratch/mch/mhrvo/E_PROFILE_calout_v22_04 \
      --census validation/scope_l1_2026_census.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

V22_VERSION_CODE = 220
#: Molecular vintages this release supersedes (calibration/io/output.py VERSION_CODES). A row
#: carrying one of these in the NEW tree means an old result survived the recompute. The cloud code
#: (1000) is deliberately absent: it is a different retrieval, not an older Rayleigh, and it shares
#: the yearly NetCDF with the Rayleigh rows.
_SUPERSEDED_MOLECULAR = {25, 90, 95, 100, 110, 120, 200, 205, 300, 310}

#: Nights the calibration must keep REJECTING. Both are genuine Payerne outliers that v1 accepted
#: and v2 correctly throws out (6x and 1.9x the surrounding level); if a release starts accepting
#: them, its gates have been loosened by accident.
CANARIES = [("0-20000-0-06610_A", "rayleigh", "20260531"),
            ("0-20000-0-06610_A", "rayleigh", "20260608")]

#: The one instrument type whose constants this release is MEANT to move (constructor WV spectrum).
EXPECTED_TO_CHANGE = "CL61"


def _read_cal(path: Path) -> dict:
    """{(method, date): (flag, cal_value)} from a <key>_cal.csv."""
    out = {}
    try:
        with path.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                try:
                    flag = float(r["flag"])
                except (TypeError, ValueError):
                    continue
                try:
                    val = float(r["cal_value"])
                except (TypeError, ValueError):
                    val = float("nan")
                out[(r["method"], r["date"])] = (flag, val)
    except OSError:
        return {}
    return out


def _streams(tree: Path) -> list:
    return sorted(p.name for p in tree.iterdir() if p.is_dir() and not p.name.startswith("_"))


def _types_from_census(path) -> dict:
    if not path:
        return {}
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
    except OSError:
        return {}
    rows = raw if isinstance(raw, list) else raw.get("streams", [])
    return {f"{r['wmo']}_{r['ident']}": str(r.get("type", "")) for r in rows}


def gate_vintage(new: Path, old: Path | None) -> tuple:
    """G1 -- homogeneous v2.2, no duplicate rows, and the old tree untouched by v2.2."""
    msgs, bad = [], 0
    try:
        import netCDF4
    except ImportError:
        netCDF4 = None

    dup_streams = 0
    for k in _streams(new):
        cal = new / k / f"{k}_cal.csv"
        if not cal.exists():
            continue
        seen = Counter()
        with cal.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                seen[(r["method"], r["date"])] += 1
        d = sum(1 for v in seen.values() if v > 1)
        if d:
            dup_streams += 1
            if dup_streams <= 3:
                msgs.append(f"  duplicate (method,date) rows: {k} x{d}")
    if dup_streams:
        msgs.append(f"  {dup_streams} stream(s) carry duplicate rows")
        bad += 1
    else:
        msgs.append(f"  no duplicate (method,date) rows in {len(_streams(new))} streams")

    if netCDF4 is None:
        msgs.append("  netCDF4 unavailable -> version check SKIPPED")
    else:
        checked = wrong = seen_files = noversion = 0
        for k in _streams(new)[:40]:          # a sample is enough: the code stamps every row alike
            for nc in sorted((new / k).glob("*/ALC_calibration_*.nc")):
                seen_files += 1
                try:
                    with netCDF4.Dataset(nc) as ds:
                        if "calibration_version" not in ds.variables:
                            noversion += 1
                            continue
                        vals = set(int(v) for v in ds.variables["calibration_version"][:].ravel())
                except OSError:
                    continue
                checked += 1
                # The yearly NetCDF carries BOTH methods -- the run is --methods rayleigh,cloud and
                # they append to the same file -- so 1000 (cloud_oconnor) beside 220 is correct, not
                # a mixed vintage. Demanding {220} alone failed every CL31/CL51 stream, where
                # Rayleigh never succeeds at 910 nm and only cloud rows exist. What actually must
                # hold is that no row carries a SUPERSEDED molecular vintage.
                stale = vals & _SUPERSEDED_MOLECULAR
                if stale:
                    wrong += 1
                    msgs.append(f"  {nc.name}: superseded molecular version(s) {sorted(stale)}")
        msgs.append(f"  no superseded molecular vintage in {checked - wrong}/{checked} sampled "
                    f"NetCDFs (220 = v2.2 Rayleigh, 1000 = cloud; both are expected)")
        if wrong:
            bad += 1
        # A vintage gate that verified nothing must not report success: either there are no
        # NetCDFs at all (the archive is unusable) or they carry no version variable (the vintage
        # is unverifiable) -- both are failures, not passes.
        if seen_files == 0:
            msgs.append("  NO calibration NetCDF found in the sampled streams -> cannot verify "
                        "the vintage, and the archive is missing its exchange product")
            bad += 1
        elif checked == 0:
            msgs.append(f"  {noversion} NetCDF(s) carry no calibration_version variable -> vintage "
                        f"unverifiable")
            bad += 1
        # leftovers of the chunked layout mean the promotion step did not run
        left = [k for k in _streams(new) if (new / k / "parts_nc").exists()]
        if left:
            msgs.append(f"  {len(left)} stream(s) still have parts_nc/ -> NetCDFs not promoted "
                        f"into <key>/<year>/ (e.g. {left[0]})")
            bad += 1
        else:
            msgs.append("  no parts_nc/ leftovers: NetCDFs are in the operational layout")

    if old is None:
        msgs.append("  old tree not supplied -> 'old tree free of v2.2' SKIPPED")
    elif netCDF4 is not None:
        hits = 0
        for k in _streams(old)[:40]:
            for nc in sorted((old / k).glob("*/ALC_calibration_*.nc")):
                try:
                    with netCDF4.Dataset(nc) as ds:
                        if "calibration_version" in ds.variables:
                            if V22_VERSION_CODE in set(
                                    int(v) for v in ds.variables["calibration_version"][:].ravel()):
                                hits += 1
                except OSError:
                    continue
        msgs.append(f"  old tree: {hits} sampled NetCDF(s) contain v2.2 rows (must be 0)")
        if hits:
            bad += 1
    return bad == 0, msgs


def gate_continuity(new: Path, ref: Path | None, types: dict) -> tuple:
    """G2 -- only the CL61 constants may have moved, and the canaries stay rejected."""
    msgs, bad = [], 0
    if ref is None:
        msgs.append("  reference tree not supplied -> SKIPPED")
        return None, msgs

    common = [k for k in _streams(new) if (ref / k).is_dir()]
    moved = defaultdict(int)
    compared = defaultdict(int)
    worst = {}
    for k in common:
        t = types.get(k, "?")
        a = _read_cal(new / k / f"{k}_cal.csv")
        b = _read_cal(ref / k / f"{k}_cal.csv")
        for key in set(a) & set(b):
            fa, va = a[key]
            fb, vb = b[key]
            compared[t] += 1
            if fa != fb or not _same(va, vb):
                moved[t] += 1
                if va == va and vb == vb and vb:
                    rel = abs(va - vb) / abs(vb) * 100.0
                    if rel > worst.get(t, (0.0, ""))[0]:
                        worst[t] = (rel, f"{k} {key[0]} {key[1]}")
    for t in sorted(compared):
        n, m = compared[t], moved[t]
        tag = "expected" if t == EXPECTED_TO_CHANGE else "MUST BE 0"
        extra = f", worst {worst[t][0]:.2f}% at {worst[t][1]}" if t in worst else ""
        msgs.append(f"  {t:8s}: {m:6d}/{n:6d} rows moved ({tag}){extra}")
        if m and t != EXPECTED_TO_CHANGE:
            bad += 1
    if not compared:
        msgs.append("  no rows in common with the reference -> cannot judge continuity")
        bad += 1

    for key, method, date in CANARIES:
        cal = _read_cal(new / key / f"{key}_cal.csv")
        hit = cal.get((method, date))
        if hit is None:
            msgs.append(f"  canary {key} {date}: ABSENT (cannot confirm rejection)")
            bad += 1
        elif hit[0] > 0:
            msgs.append(f"  canary {key} {date}: ACCEPTED with flag {hit[0]} -- must be rejected")
            bad += 1
        else:
            msgs.append(f"  canary {key} {date}: still rejected (flag {hit[0]:g})")
    return bad == 0, msgs


def _same(a, b) -> bool:
    if a != a and b != b:          # both NaN
        return True
    if a != a or b != b:
        return False
    return abs(a - b) <= 1e-12 * max(1.0, abs(a), abs(b))


def gate_completeness(new: Path, ref: Path | None, census: dict) -> tuple:
    """G3 -- every stream carries every artefact the operational server needs."""
    msgs, bad = [], 0
    expected = set(_streams(ref)) if ref is not None else set()
    expected |= set(census)
    if not expected:
        expected = set(_streams(new))
        msgs.append("  no reference/census -> completeness judged against the new tree itself")
    have = set(_streams(new))
    missing = sorted(expected - have)
    if missing:
        msgs.append(f"  {len(missing)} stream(s) MISSING from the new tree "
                    f"(e.g. {', '.join(missing[:4])})")
        bad += 1
    else:
        msgs.append(f"  all {len(expected)} expected streams present")

    required = ["{k}_cal.csv", "{k}_kalman.csv", "{k}_hk.csv", "{k}_status.csv",
                "{k}_sens.csv", "{k}_omb.csv", "_sens_cache.npz", "_omb_cache.npz"]
    counts = Counter()
    for k in have:
        for pat in required:
            if (new / k / pat.format(k=k)).exists():
                counts[pat] += 1
        if list((new / k).glob("*/ALC_calibration_*.nc")):
            counts["<year>/*.nc"] += 1
    n = len(have)
    for pat in required + ["<year>/*.nc"]:
        got = counts[pat]
        flag = "" if got == n else "   <-- INCOMPLETE"
        msgs.append(f"  {pat:20s} {got:4d}/{n}{flag}")
        if got != n:
            bad += 1
    # the caches must match their CSVs one-for-one, else the regression guard latches
    if counts["_omb_cache.npz"] != counts["{k}_omb.csv"]:
        msgs.append("  _omb_cache.npz count != _omb.csv count -> the daily run would SKIP OmB")
        bad += 1
    if counts["_sens_cache.npz"] != counts["{k}_sens.csv"]:
        msgs.append("  _sens_cache.npz count != _sens.csv count -> the daily run would SKIP sens")
        bad += 1
    return bad == 0, msgs


def gate_cloudcover(new: Path) -> tuple:
    """G4 -- the new status columns exist and are physical."""
    msgs, bad = [], 0
    have = _streams(new)
    withcol = withval = 0
    srcs = Counter()
    out_of_range = 0
    for k in have:
        p = new / k / f"{k}_status.csv"
        if not p.exists():
            continue
        with p.open(newline="", encoding="utf-8") as fh:
            rd = csv.DictReader(fh)
            if "mean_cloud_cover" not in (rd.fieldnames or []):
                continue
            withcol += 1
            any_val = False
            for r in rd:
                v = r.get("mean_cloud_cover", "")
                if v == "":
                    continue
                any_val = True
                srcs[r.get("cloud_src", "")] += 1
                try:
                    f = float(v)
                except ValueError:
                    out_of_range += 1
                    continue
                if not (0.0 <= f <= 8.0):
                    out_of_range += 1
            withval += 1 if any_val else 0
    msgs.append(f"  mean_cloud_cover column present in {withcol}/{len(have)} status files")
    msgs.append(f"  streams with at least one value: {withval}")
    msgs.append(f"  sources: {dict(srcs)}")
    if withcol == 0:
        msgs.append("  no stream carries the column -> the producer did not run")
        bad += 1
    if out_of_range:
        msgs.append(f"  {out_of_range} value(s) outside 0..8 octas -> screening is wrong")
        bad += 1
    else:
        msgs.append("  every value inside 0..8 octas")
    return bad == 0, msgs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", type=Path, required=True, help="the release archive under test")
    ap.add_argument("--ref", type=Path, default=None,
                    help="the already-validated v2.2 archive (continuity reference)")
    ap.add_argument("--old", type=Path, default=None,
                    help="the live operational tree (checked for v2.2 contamination)")
    ap.add_argument("--census", default=None, help="census JSON (expected stream list)")
    ap.add_argument("--json", type=Path, default=None, help="write a machine-readable report here")
    args = ap.parse_args()

    types = _types_from_census(args.census)
    gates = [
        ("G1 vintage", gate_vintage(args.new, args.old)),
        ("G2 continuity", gate_continuity(args.new, args.ref, types)),
        ("G3 completeness", gate_completeness(args.new, args.ref, types)),
        ("G4 cloudcover", gate_cloudcover(args.new)),
    ]
    report, failed, skipped = {}, 0, 0
    for name, (ok, msgs) in gates:
        state = "SKIPPED" if ok is None else ("PASS" if ok else "FAIL")
        print(f"\n=== {name}: {state} ===")
        for m in msgs:
            print(m)
        report[name] = {"state": state, "detail": msgs}
        failed += 1 if ok is False else 0
        skipped += 1 if ok is None else 0

    print("\n" + "=" * 60)
    if failed:
        print(f"VERDICT: {failed} gate(s) FAILED -> do NOT flip ops/config.sh")
    elif skipped:
        print(f"VERDICT: all evaluated gates passed, {skipped} SKIPPED -> supply the missing tree "
              f"before flipping")
    else:
        print("VERDICT: all gates passed -> the archive may be flipped (see "
              "doc/OPERATIONS_v22_cutover.md)")
    if args.json:
        args.json.write_text(json.dumps(report, indent=1), encoding="utf-8")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
