"""
scope_l1_2026.py — enumerate the E-PROFILE L1 2026 archive and classify every instrument by
reading the NetCDF metadata (the `instrument_type` global attribute + station coordinates).

instruments.json is NOT used for typing: it predates the CL61 fleet (0 CL61 entries), so it
cannot identify the very instruments we need. The L1 files carry a reliable `instrument_type`
attribute ('CHM15k' | 'CL61' | 'CL51' | 'CL31' | 'Mini-MPL' | ...) and station_lat/lon/alt.

L1 layout: <root>/<WMO>/YYYY/MM/L1_<WMO>_<id><YYYYMMDD>.nc  (one daily file per instrument).
One <WMO> dir can host several instruments distinguished by the 1-char identifier.

Writes scope_l1_2026.json: the run manifest (selected CHM15k + all Mini-MPL + all CL61 with
enough coverage), each {label, group, wmo, ident, type, site, lat, lon, alt, n_days, ...}.
"""
from __future__ import annotations
import glob
import json
import os
from collections import defaultdict
from pathlib import Path

from netCDF4 import Dataset

ROOT = Path("D:/E-PROFILE_L1_2026")
REPO = Path(__file__).resolve().parents[1]
OUT_JSON = REPO / "validation" / "scope_l1_2026.json"
CACHE = REPO / "validation" / "scope_l1_2026_census.json"   # full classified census (cache)

MIN_DAYS = 20      # need enough nights for the drift-insensitive variability metrics

# The "selected" CHM15k from the long-run method study (longrun_methods.py SEL): exact
# (WMO, identifier) so we reproduce the same 10 instruments.
SELECTED_CHM15K = {
    ("0-20000-0-06610", "A"): "Payerne",   ("0-20000-0-10393", "0"): "Lindenberg",
    ("0-380-5-1", "0"): "Aosta",           ("0-250-1001-07151", "B"): "Palaiseau",
    ("0-20008-0-UGR", "A"): "Granada",     ("0-20008-0-INO", "A"): "Magurele",
    ("0-20000-0-01311", "A"): "Bergen",    ("0-20000-0-01492", "A"): "Oslo",
    ("0-20000-0-10140", "0"): "Hamburg",   ("0-20000-0-10962", "0"): "Hohenpeiss",
}


def parse_l1_name(fname: str):
    """'L1_0-20000-0-06610_A20260201.nc' -> ('0-20000-0-06610', 'A', '20260201')."""
    core = fname[3:-3]
    date = core[-8:]
    wmo_ident = core[:-8]
    if "_" not in wmo_ident or not date.isdigit():
        return None
    wmo, ident = wmo_ident.rsplit("_", 1)
    return wmo, ident, date


def scan_archive():
    """{(wmo, ident): (sorted dates, one representative file path)}."""
    present = defaultdict(set)
    one_file = {}
    for wmo_dir in sorted(ROOT.iterdir()):
        if not wmo_dir.is_dir():
            continue
        for f in glob.glob(str(wmo_dir / "2026" / "*" / "L1_*.nc")):
            parsed = parse_l1_name(os.path.basename(f))
            if parsed is None:
                continue
            wmo, ident, date = parsed
            present[(wmo, ident)].add(date)
            one_file.setdefault((wmo, ident), f)
    return {k: (sorted(v), one_file[k]) for k, v in present.items()}


def read_meta(path):
    """(instrument_type, site, lat, lon, alt) from one L1 file."""
    with Dataset(path, "r") as d:
        itype = getattr(d, "instrument_type", "?")
        title = getattr(d, "title", "")
        site = title.split(" ")[0] if title else "?"
        def v(name, default):
            if name in d.variables:
                try:
                    return float(d.variables[name][:].flatten()[0])
                except Exception:
                    return default
            return default
        lat = v("station_latitude", 0.0)
        lon = v("station_longitude", 0.0)
        alt = v("station_altitude", 0.0)
    return itype, site, lat, lon, alt


def build_census(streams):
    rows = []
    for (wmo, ident), (dates, f) in streams.items():
        try:
            itype, site, lat, lon, alt = read_meta(f)
        except Exception as e:
            print(f"  WARN read {wmo}_{ident}: {e}")
            continue
        rows.append(dict(wmo=wmo, ident=ident, type=itype, site=site,
                         lat=lat, lon=lon, alt=alt, n_days=len(dates),
                         first=dates[0], last=dates[-1]))
    return rows


def main():
    streams = scan_archive()
    print(f"L1 archive: {len(streams)} (WMO, identifier) instrument-streams present")
    rows = build_census(streams)
    CACHE.write_text(json.dumps(rows, indent=1), encoding="utf-8")

    census = defaultdict(int)
    for r in rows:
        census[r["type"]] += 1
    print("\nType census (by L1 `instrument_type` attribute):")
    for t, n in sorted(census.items(), key=lambda kv: -kv[1]):
        print(f"  {t:12s} {n}")

    def desc(r):
        return f"{r['site']} [{r['wmo']}_{r['ident']}] {r['type']}  {r['n_days']}d {r['first']}-{r['last']}"

    by_key = {(r["wmo"], r["ident"]): r for r in rows}

    # selected CHM15k (exact WMO,ident; fall back to any CHM15k at that WMO)
    sel_chm = []
    for (wmo, ident), name in SELECTED_CHM15K.items():
        r = by_key.get((wmo, ident))
        if r is None or r["type"] != "CHM15k":
            alt = [x for x in rows if x["wmo"] == wmo and x["type"] == "CHM15k"]
            r = alt[0] if alt else None
        if r is None:
            print(f"   MISSING selected CHM15k: {name} ({wmo}_{ident})")
            continue
        sel_chm.append({**r, "name": name})

    mpl = [r for r in rows if r["type"] in ("Mini-MPL", "MPL")]
    cl61 = [r for r in rows if r["type"] == "CL61"]

    print(f"\nSelected CHM15k ({len(sel_chm)}):")
    for r in sel_chm:
        print("   ", desc(r))
    print(f"\nMini-MPL present ({len(mpl)}):")
    for r in sorted(mpl, key=lambda r: -r["n_days"]):
        print("   ", desc(r), "" if r["n_days"] >= MIN_DAYS else "  (<20d skip)")
    print(f"\nAll CL61 present ({len(cl61)}; {sum(r['n_days']>=MIN_DAYS for r in cl61)} with >= {MIN_DAYS}d):")
    for r in sorted(cl61, key=lambda r: -r["n_days"]):
        print("   ", desc(r), "" if r["n_days"] >= MIN_DAYS else "  (<20d skip)")

    # ---- run manifest ----------------------------------------------------
    def slabel(site, suffix):
        base = site.split(",")[0].strip().replace(" ", "").replace("_", "-")
        return f"{base}_{suffix}"

    manifest = []
    for r in sel_chm:
        manifest.append({**{k: r[k] for k in ("wmo", "ident", "type", "site", "lat", "lon",
                                              "alt", "n_days", "first", "last")},
                         "label": slabel(r.get("name", r["site"]), "CHM15k"), "group": "CHM15k"})
    for r in sorted(mpl, key=lambda r: r["site"]):
        if r["n_days"] >= MIN_DAYS:
            manifest.append({**{k: r[k] for k in ("wmo", "ident", "type", "site", "lat", "lon",
                                                  "alt", "n_days", "first", "last")},
                             "label": slabel(r["site"], "MPL"), "group": "Mini-MPL"})
    for r in sorted(cl61, key=lambda r: r["site"]):
        if r["n_days"] >= MIN_DAYS:
            manifest.append({**{k: r[k] for k in ("wmo", "ident", "type", "site", "lat", "lon",
                                                  "alt", "n_days", "first", "last")},
                             "label": slabel(r["site"], "CL61"), "group": "CL61"})

    # de-duplicate labels (two streams from the same town)
    seen = defaultdict(int)
    for m in manifest:
        seen[m["label"]] += 1
        if seen[m["label"]] > 1:
            m["label"] = f"{m['label']}-{m['ident']}"

    OUT_JSON.write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    ng = defaultdict(int)
    for m in manifest:
        ng[m["group"]] += 1
    print(f"\nWrote {OUT_JSON.name}: {len(manifest)} instruments "
          f"(CHM15k={ng['CHM15k']}, Mini-MPL={ng['Mini-MPL']}, CL61={ng['CL61']})")


if __name__ == "__main__":
    main()
