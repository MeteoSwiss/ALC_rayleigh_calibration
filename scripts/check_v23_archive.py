#!/usr/bin/env python3
"""Acceptance gates for the v2.3 operational archive, run BEFORE the cutover (exit 0 = go).

  G1 vintage      one row per (stream, method, date); every calibration NetCDF carries only the v2.3
                  codes 230 (Rayleigh), 231 (Rayleigh, polar-summer low-sun night) and 1100 (cloud,
                  eta tables v2); no parts_nc/ leftovers.
  G2 continuity   against the release tree (--ref): inside the release window every calibration row
                  and every Kalman row is bit-identical, except where v2.3.0 is MEANT to differ:
                  the cloud rows of the four CL31 -> CL51 swap streams (per-day type) and the Kalman of
                  the streams with a kalman_restart event or a station_kalman_params entry.
  G3 content      the five uncertainty components in every NetCDF with Rayleigh rows; flag 0.25 only
                  north of 57 N; no stream with flag -99 on more than 1 % of its rows (the first-pass
                  failure mode: every 910 nm night -99); isolated -99 days (a corrupt L1 file) are
                  listed; cloud CSV rows carry 1100.
  G4 completeness every stream has cal / kalman / status CSVs (hk only where the L1 has housekeeping
                  fields: listed, not failed); OmB and sensitivity outputs and both caches except the
                  known exceptions (out-of-domain stations, dead streams).
  G5 canaries     two Payerne nights that every v2 release must keep rejecting.

Gates that cannot be evaluated are reported SKIPPED, never passed.

    python scripts/check_v23_archive.py --new <v2.3 tree> --ref <release tree> \\
        --census validation/scope_l1_2026_census.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V23_CODES = {230, 231, 1100}
CLOUD_CODE = 1100
RELEASE_END = "20261101"
UNCERTAINTY_COMPONENTS = ["lidar_constant_uncertainty_ensemble", "lidar_constant_uncertainty_vertical",
                          "lidar_constant_uncertainty_molecular", "lidar_constant_uncertainty_aerosol",
                          "lidar_constant_uncertainty_noise"]
CANARIES = [("0-20000-0-06610_A", "rayleigh", "20260531"),
            ("0-20000-0-06610_A", "rayleigh", "20260608")]
#: Release report v2.3 item 12: the OmB reader has only the Europe CAMS monthlies.
NO_OMB_KNOWN = {"0-20008-0-EDT_A", "0-20008-0-EDT_B", "0-20008-0-LAU_A", "0-554-10-UAAONIWACL31_A"}
POLAR_MIN_LAT = 57.0


def _streams(tree: Path) -> list:
    return sorted(p.name for p in tree.iterdir() if p.is_dir() and not p.name.startswith("_"))


def _rows(path: Path) -> list:
    try:
        with path.open(newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def _census(path) -> dict:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = raw if isinstance(raw, list) else raw.get("streams", [])
    return {f"{r['wmo']}_{r['ident']}": r for r in rows}


def _expected_to_move() -> tuple[set, set]:
    """(streams whose cloud rows may move: the CL31/CL51 swaps, streams whose Kalman may move)."""
    ev = json.loads((REPO / "monitoring" / "instrument_events.json").read_text(encoding="utf-8"))
    restarts = {k for k, v in ev.items() if isinstance(v, list) and any(e.get("kalman_restart") for e in v)}
    swaps = {k for k, v in ev.items() if isinstance(v, list)
             and any(e.get("kalman_restart") and "CL51" in e.get("label", "") for e in v)}
    opts = json.loads((REPO / "options.json").read_text(encoding="utf-8"))
    kalman = restarts | set(opts.get("station_kalman_params", {}))
    return swaps, kalman | swaps


def gate_vintage(new: Path) -> tuple:
    import netCDF4
    msgs, ok = [], True
    dup = [k for k in _streams(new)
           if any(n > 1 for n in Counter((r["method"], r["date"]) for r in _rows(new / k / f"{k}_cal.csv")).values())]
    msgs.append(f"  duplicate (method,date) rows: {len(dup)} streams {dup[:3]}")
    ok &= not dup
    codes, bad, files = Counter(), [], 0
    for nc in sorted(new.glob("*/*/ALC_calibration_*.nc")):
        files += 1
        with netCDF4.Dataset(nc) as ds:
            vals = {int(v) for v in ds["calibration_version"][:].compressed()} if "calibration_version" in ds.variables else set()
        codes.update(vals)
        if not vals or vals - V23_CODES:
            bad.append(f"{nc.name}: {sorted(vals)}")
    msgs.append(f"  {files} NetCDF, version codes {dict(sorted(codes.items()))}; outside {sorted(V23_CODES)}: {len(bad)} {bad[:3]}")
    ok &= files > 0 and not bad
    left = [k for k in _streams(new) if (new / k / "parts_nc").exists()]
    msgs.append(f"  parts_nc/ leftovers: {len(left)}")
    ok &= not left
    return ok, msgs


def gate_continuity(new: Path, ref: Path | None) -> tuple:
    if ref is None:
        return None, ["  --ref not supplied -> SKIPPED"]
    swaps, kalman_moves = _expected_to_move()
    msgs, ok = [], True
    n_cmp = n_moved = n_expected = 0
    k_cmp = k_moved = 0
    unexpected = []
    for k in _streams(ref):
        if not (new / k).is_dir():
            unexpected.append(f"{k}: missing in new tree")
            continue
        a = {(r["method"], r["date"]): (r["flag"], r["cal_value"]) for r in _rows(new / k / f"{k}_cal.csv")}
        b = {(r["method"], r["date"]): (r["flag"], r["cal_value"]) for r in _rows(ref / k / f"{k}_cal.csv")}
        for key in b:
            if key[1] > RELEASE_END:
                continue
            n_cmp += 1
            if a.get(key) != b[key]:
                if k in swaps and key[0] == "cloud":
                    n_expected += 1
                else:
                    n_moved += 1
                    if len(unexpected) < 8:
                        unexpected.append(f"{k} {key}: {b[key]} -> {a.get(key)}")
        ka = {(r["method"], r["date"]): r["kalman"] for r in _rows(new / k / f"{k}_kalman.csv")}
        kb = {(r["method"], r["date"]): r["kalman"] for r in _rows(ref / k / f"{k}_kalman.csv")}
        if k in kalman_moves:
            continue
        for key in kb:
            if key[1] > RELEASE_END:
                continue
            k_cmp += 1
            if ka.get(key) != kb[key]:
                k_moved += 1
                if len(unexpected) < 8:
                    unexpected.append(f"{k} kalman {key}: {kb[key]} -> {ka.get(key)}")
    msgs.append(f"  calibration rows in the release window: {n_cmp} compared, {n_expected} moved as intended "
                f"(cloud rows of {sorted(swaps)}), {n_moved} moved unexpectedly")
    msgs.append(f"  Kalman rows outside {len(kalman_moves)} intended streams: {k_cmp} compared, {k_moved} moved")
    msgs += [f"    {u}" for u in unexpected]
    ok &= n_moved == 0 and k_moved == 0 and not any("missing" in u for u in unexpected)
    return ok, msgs


def gate_content(new: Path, census: dict) -> tuple:
    import netCDF4
    msgs, ok = [], True
    missing_u, checked = [], 0
    for nc in sorted(new.glob("*/*/ALC_calibration_*.nc")):
        with netCDF4.Dataset(nc) as ds:
            if "calibration_version" not in ds.variables:
                continue
            vals = {int(v) for v in ds["calibration_version"][:].compressed()}
            if not vals & {230, 231}:
                continue
            checked += 1
            miss = [c for c in UNCERTAINTY_COMPONENTS if c not in ds.variables]
            if miss:
                missing_u.append(f"{nc.name}: {miss}")
    msgs.append(f"  uncertainty components: {checked} Rayleigh NetCDF, {len(missing_u)} incomplete {missing_u[:2]}")
    ok &= not missing_u
    low_sun_south, minus99, cloud_noversion, nrows = [], Counter(), Counter(), Counter()
    for k in _streams(new):
        lat = float(census.get(k, {}).get("lat", 90.0))
        for r in _rows(new / k / f"{k}_cal.csv"):
            nrows[k] += 1
            try:
                f = float(r["flag"])
            except (TypeError, ValueError):
                continue
            if f == 0.25 and lat < POLAR_MIN_LAT:
                low_sun_south.append(f"{k} {r['date']}")
            if f == -99:
                minus99[k] += 1
            if r["method"] == "cloud" and str(r.get("version", "")).strip() != str(CLOUD_CODE):
                cloud_noversion[k] += 1
    msgs.append(f"  flag 0.25 south of {POLAR_MIN_LAT} N: {len(low_sun_south)} {low_sun_south[:3]}")
    systematic = {k: n for k, n in minus99.items() if n > 0.01 * nrows[k]}
    msgs.append(f"  flag -99: {sum(minus99.values())} isolated rows in {len(minus99)} streams {sorted(minus99.items())}; "
                f"systematic (> 1 % of a stream): {len(systematic)} {sorted(systematic.items())[:3]}")
    msgs.append(f"  cloud rows without version {CLOUD_CODE}: {sum(cloud_noversion.values())} in {len(cloud_noversion)} streams")
    ok &= not low_sun_south and not systematic and not cloud_noversion
    return ok, msgs


def gate_completeness(new: Path, census: dict) -> tuple:
    msgs, ok = [], True
    core = ["_cal.csv", "_kalman.csv", "_status.csv"]
    addons = ["_omb.csv", "_sens.csv"]
    caches = ["_omb_cache.npz", "_sens_cache.npz"]
    no_core, no_add, dead, no_hk = [], [], [], []
    for k in _streams(new):
        d = new / k
        if [s for s in core if not (d / f"{k}{s}").exists()]:
            no_core.append(k)
        if not (d / f"{k}_hk.csv").exists():
            no_hk.append(k)
        kal = _rows(d / f"{k}_kalman.csv")
        if not kal:
            dead.append(k)                       # no Kalman -> nothing for the add-ons to use
            continue
        miss = [s for s in addons if not (d / f"{k}{s}").exists()]
        miss += [c for c in caches if not (d / c).exists()]
        if k in NO_OMB_KNOWN:
            miss = [m for m in miss if "omb" not in m]
        if miss:
            no_add.append(f"{k} {miss}")
    msgs.append(f"  {len(_streams(new))} streams; missing core CSVs: {len(no_core)} {no_core[:3]}")
    msgs.append(f"  no _hk.csv (L1 without housekeeping fields, as in v2.2): {len(no_hk)} {no_hk[:5]}")
    msgs.append(f"  no Kalman (dead or too few nights, add-ons not expected): {len(dead)} {dead[:5]}")
    msgs.append(f"  missing add-on outputs (known out-of-domain OmB excused): {len(no_add)} {no_add[:3]}")
    not_in_census = [k for k in _streams(new) if k not in census]
    msgs.append(f"  streams not in the census: {len(not_in_census)} {not_in_census[:3]}")
    ok &= not no_core and not no_add
    return ok, msgs


def gate_canaries(new: Path) -> tuple:
    msgs, ok = [], True
    for k, m, d in CANARIES:
        row = next((r for r in _rows(new / k / f"{k}_cal.csv") if r["method"] == m and r["date"] == d), None)
        flag = float(row["flag"]) if row else None
        good = row is not None and flag not in (1.0, 0.5, 0.25)
        msgs.append(f"  {k} {m} {d}: flag {flag} -> {'rejected' if good else 'ACCEPTED or missing'}")
        ok &= good
    return ok, msgs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True, type=Path)
    ap.add_argument("--ref", type=Path, default=None, help="the v2.3 release tree (continuity gate)")
    ap.add_argument("--census", type=Path, default=REPO / "validation" / "scope_l1_2026_census.json")
    a = ap.parse_args()
    census = _census(a.census)
    gates = [("G1 vintage", gate_vintage(a.new)), ("G2 continuity", gate_continuity(a.new, a.ref)),
             ("G3 content", gate_content(a.new, census)), ("G4 completeness", gate_completeness(a.new, census)),
             ("G5 canaries", gate_canaries(a.new))]
    failed = 0
    for name, (ok, msgs) in gates:
        state = "SKIPPED" if ok is None else ("PASS" if ok else "FAIL")
        failed += ok is False
        print(f"{name}: {state}")
        print("\n".join(msgs))
    print("RESULT:", "GO" if failed == 0 else f"NO-GO ({failed} gate(s) failed)")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
