"""
download_cloudnet_cl61.py — download CL61 raw NetCDF from the Cloudnet open API into the
RAW per-day-folder layout (<root>/<site>/YYYYMMDD/<filename>.nc) used by data_level=RAW.

Usage:  python download_cloudnet_cl61.py <site> <dateFrom> <dateTo>
        site in {lindenberg, hyytiala};  dates YYYY-MM-DD.

Cloudnet raw-files API: GET /api/raw-files?site=&instrumentPid=&dateFrom=&dateTo=
Each result has 'filename', 'measurementDate', 'downloadUrl' (open, no auth).
Handles both daily and 5-min files (all files for a day land in the same day-folder; the
RAW reader concatenates them).
"""
from __future__ import annotations
import sys, json, urllib.request, urllib.parse
from pathlib import Path

API = "https://cloudnet.fmi.fi/api/raw-files"
SITES = {
    "lindenberg": dict(site="lindenberg", pid="https://hdl.handle.net/21.12132/3.695573e5981845d9"),
    "hyytiala":   dict(site="hyytiala",   pid="https://hdl.handle.net/21.12132/3.241bda142975460b"),
}
OUT = Path("R:/CL61/RAW_cloudnet_dl")


def query(site, pid, dfrom, dto):
    q = {"site": site, "instrumentPid": pid, "dateFrom": dfrom, "dateTo": dto}
    url = API + "?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def main():
    if len(sys.argv) < 4:
        print("usage: download_cloudnet_cl61.py <site> <dateFrom YYYY-MM-DD> <dateTo>")
        return
    key, dfrom, dto = sys.argv[1], sys.argv[2], sys.argv[3]
    s = SITES[key]
    files = query(s["site"], s["pid"], dfrom, dto)
    # Each day usually has BOTH a daily concatenated file (e.g. 20240929_lindenberg_cl61_*.nc,
    # nt=8640) AND the 288 5-min 'live_*' files (same data). Prefer the daily file; fall back
    # to the 5-min files only on days with no daily (the RAW reader handles either, but never
    # both, to avoid double-counting).
    from collections import defaultdict
    byday = defaultdict(list)
    for f in files:
        byday[str(f["measurementDate"]).replace("-", "")[:8]].append(f)
    wanted = []
    for day, fs in byday.items():
        daily = [f for f in fs if not f["filename"].lower().startswith("live")]
        wanted.extend(daily if daily else fs)
    print(f"{key}: {len(files)} files -> {len(wanted)} to fetch ({len(byday)} days; daily-preferred) {dfrom}..{dto}")
    n_dl, n_skip, n_err, nbytes = 0, 0, 0, 0
    for f in wanted:
        fn = f["filename"]
        ymd = str(f["measurementDate"]).replace("-", "")[:8]
        dest = OUT / key / ymd / fn
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and dest.stat().st_size > 0:
            n_skip += 1
            continue
        try:
            urllib.request.urlretrieve(f["downloadUrl"], dest)
            n_dl += 1
            nbytes += dest.stat().st_size
        except Exception as e:
            n_err += 1
            print(f"  ERR {fn}: {e}")
    print(f"  downloaded={n_dl} skipped={n_skip} errors={n_err} "
          f"({nbytes/1e6:.0f} MB);  days={len(set((OUT/key).glob('*')))}")
    print("DOWNLOAD_DONE")


if __name__ == "__main__":
    main()
