#!/usr/bin/env python3
"""Build the {station-key -> CEDA L2 URL} map used by the dashboard's per-page "CEDA L2 data" link.

The CEDA E-PROFILE archive lays daily L2 files out as::

    https://data.ceda.ac.uk/badc/eprofile/data/daily_files/<country>/<station>/<institution>-<maker>-<model>_<channel>

e.g. Payerne -> ``meteoswiss-lufft-chm15k_A``, ``meteoswiss-vaisala-cl31_B``,
``bern-university-vaisala-cl61_C``. The institution / station / country slugs are CEDA-curated and
NOT reliably derivable from our metadata, so this script **crawls the live tree** and matches each of
our station keys to a real folder by (country, maker/model, channel, fuzzy station name). Because every
emitted URL corresponds to a folder actually found in the listing, the links are validated by
construction ("double-check all the links"); unmatched keys are logged, never guessed.

Run occasionally (needs internet); the dashboard reads the JSON via ``--ceda-links`` / ``$ALC_CEDA_LINKS``.

  python scripts/build_ceda_links.py --manifest validation/scope_l1_2026_census.json \
      --l2dir A:/E-PROFILE_L2_monthly --out validation/ceda_links.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from monitoring import config, index  # noqa: E402

CEDA_ROOT = "https://data.ceda.ac.uk/badc/eprofile/data/daily_files"

# Instrument type -> the "<maker>-<model>" part of the CEDA folder name. Mini-MPL/MPL are not part of
# the E-PROFILE CEDA distribution, so they simply go unmatched (no link).
MAKER_MODEL = {
    "CHM15k": "lufft-chm15k", "CHM8k": "lufft-chm8k",
    "CL31": "vaisala-cl31", "CL51": "vaisala-cl51", "CL61": "vaisala-cl61",
}


def slug(s: str) -> str:
    """Lowercase, non-alphanumerics -> single hyphen (matches CEDA's folder-naming convention)."""
    return re.sub(r"[^a-z0-9]+", "-", str(s or "").strip().lower()).strip("-")


def _list_dir(url: str, timeout: int = 30) -> list[str]:
    """Child folder names listed at a CEDA browse URL (parsed from the directory-listing anchors)."""
    try:
        req = urllib.request.Request(url + "/", headers={"User-Agent": "alc-ceda-linker/1.0"})
        html = urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8", "ignore")
    except Exception as e:  # noqa: BLE001
        print(f"  ! fetch failed: {url} ({e})", flush=True)
        return []
    # Anchors point at child paths '.../daily_files/<...>/<name>'; take the last path segment.
    base = url.split("data.ceda.ac.uk", 1)[-1].rstrip("/")
    names = set()
    for href in re.findall(r'href="([^"]+)"', html):
        h = href.split("data.ceda.ac.uk", 1)[-1].rstrip("/")
        if h.startswith(base + "/"):
            tail = h[len(base) + 1:]
            if tail and "/" not in tail:
                names.add(tail)
    return sorted(names)


def _match_station(station_slug: str, candidates: list[str]) -> str | None:
    """Best CEDA station-folder match for our station-name slug: exact, then prefix/substring."""
    if station_slug in candidates:
        return station_slug
    for c in candidates:
        if c == station_slug:
            return c
    for c in candidates:                       # our name is a prefix of theirs or vice-versa
        if c.startswith(station_slug) or station_slug.startswith(c):
            return c
    for c in candidates:                       # loose substring (e.g. "london" in "western-london")
        if station_slug and (station_slug in c or c in station_slug):
            return c
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=config.DEFAULT_MANIFEST)
    ap.add_argument("--l2dir", type=Path, default=config.DEFAULT_L2_DIR,
                    help="L2 archive for station name/country (same as the dashboard build)")
    ap.add_argument("--out", type=Path, default=Path("validation/ceda_links.json"))
    ap.add_argument("--only", default=None,
                    help="comma-separated country substrings to restrict the crawl (debug/testing)")
    args = ap.parse_args()

    # Reuse the dashboard's metadata (key, itype, station name, country from the L2 attrs).
    st = index._load_manifest(args.manifest)
    st = index._enrich_metadata(st, args.l2dir)
    only = [o.strip().lower() for o in args.only.split(",")] if args.only else None

    # Crawl only the countries our stations sit in (map our country string -> a CEDA country slug).
    ceda_countries = _list_dir(CEDA_ROOT)
    print(f"CEDA: {len(ceda_countries)} country folders", flush=True)
    want_countries = {}
    for _, r in st.iterrows():
        cslug = slug(r.get("country"))
        if not cslug:
            continue
        if only and not any(o in cslug for o in only):
            continue
        match = cslug if cslug in ceda_countries else _match_station(cslug, ceda_countries)
        if match:
            want_countries[cslug] = match

    # station-folder listing per needed country (cached), then instrument folders per station.
    station_cache: dict = {}
    links, unmatched = {}, []
    for _, r in st.iterrows():
        key, itype = str(r["key"]), str(r.get("itype", ""))
        mm = MAKER_MODEL.get(itype)
        channel = key.rsplit("_", 1)[-1] if "_" in key else ""
        cslug = slug(r.get("country"))
        if not mm or not channel or not cslug or cslug not in want_countries:
            if only is None or (cslug and any(o in cslug for o in only)):
                unmatched.append((key, "no maker/model or country" if not mm else "country not on CEDA"))
            continue
        country = want_countries[cslug]
        if country not in station_cache:
            station_cache[country] = _list_dir(f"{CEDA_ROOT}/{country}")
        stations = station_cache[country]
        st_match = _match_station(slug(r.get("name")), stations)
        if not st_match:
            unmatched.append((key, f"no station match in {country} for '{slug(r.get('name'))}'"))
            continue
        insts = _list_dir(f"{CEDA_ROOT}/{country}/{st_match}")
        want_suffix = f"-{mm}_{channel}"
        folder = next((f for f in insts if f.endswith(want_suffix)), None)
        if not folder:
            unmatched.append((key, f"no {want_suffix} in {country}/{st_match} ({insts})"))
            continue
        links[key] = f"{CEDA_ROOT}/{country}/{st_match}/{folder}"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(links, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"\nMatched {len(links)} / {len(st)} keys -> {args.out}", flush=True)
    if unmatched:
        print(f"Unmatched ({len(unmatched)}):", flush=True)
        for k, why in unmatched[:40]:
            print(f"  {k}: {why}", flush=True)
        if len(unmatched) > 40:
            print(f"  ... and {len(unmatched) - 40} more", flush=True)


if __name__ == "__main__":
    main()
