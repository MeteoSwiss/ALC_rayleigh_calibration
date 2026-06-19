"""
download_cloudnet_longrun.py — sampled CL61 raw download from Cloudnet for Lindenberg +
Hyytiala, into the RAW per-day-folder layout (<root>/<folder>/YYYYMMDD/*.nc).

Samples 5 nights/month (+ the previous day, so each night window is complete) over the
period. Per day it prefers the single daily concatenated file; if a site has ONLY 1-min
"live" files (Hyytiala has no daily concat), it subsamples them to ~5 min (SUBSAMPLE_LIVE).

Downloads MISSING DAYS ONLY: a day-folder that already holds >=1 .nc is left untouched
(no redownload, no daily+live duplication). Resumable. Note Hyytiala's disk folder is "hyy"
(the API site name is "hyytiala") — kept consistent with the data already on disk.
"""
from __future__ import annotations
import json, calendar, urllib.request, urllib.parse
from collections import defaultdict
from pathlib import Path

API = "https://cloudnet.fmi.fi/api/raw-files"
SITES = [
    dict(site="lindenberg", folder="lindenberg", pid="https://hdl.handle.net/21.12132/3.695573e5981845d9"),
    dict(site="hyytiala",   folder="hyy",        pid="https://hdl.handle.net/21.12132/3.241bda142975460b"),
]
OUT = Path("R:/CL61/RAW_cloudnet_dl")
SAMPLE = [3, 9, 15, 21, 27]
WANT = sorted(set(SAMPLE + [d - 1 for d in SAMPLE]))   # include previous day for the night window
SUBSAMPLE_LIVE = 5     # sites with no daily concat (Hyytiala): keep every Nth 1-min file (~5 min)
START = (2024, 9)     # earliest gap for either site (hyy ≈ Sep 2024, Lindenberg ≈ Oct 2024)
END = (2026, 6)


def months(s, e):
    y, m = s
    while (y, m) <= e:
        yield y, m
        m += 1
        if m > 12:
            m = 1
            y += 1


def query(site, pid, dfrom, dto):
    url = API + "?" + urllib.parse.urlencode({"site": site, "instrumentPid": pid,
                                              "dateFrom": dfrom, "dateTo": dto})
    with urllib.request.urlopen(url, timeout=180) as r:
        return json.load(r)


def main():
    tot_dl, tot_mb = 0, 0.0
    for s in SITES:
        site, folder, pid = s["site"], s["folder"], s["pid"]
        for y, m in months(START, END):
            dfrom = f"{y}-{m:02d}-01"
            dto = f"{y}-{m:02d}-{calendar.monthrange(y, m)[1]:02d}"
            try:
                files = query(site, pid, dfrom, dto)
            except Exception as e:
                print(f"  query err {site} {dfrom}: {e}")
                continue
            byday = defaultdict(list)
            for f in files:
                byday[str(f["measurementDate"]).replace("-", "")[:8]].append(f)
            new_days = skipped = 0
            for d in WANT:
                ymd = f"{y}{m:02d}{d:02d}"
                if ymd not in byday:
                    continue
                dest_daily = OUT / folder / f"{ymd}.nc"   # canonical single daily file
                dest_dir   = OUT / folder / ymd           # fallback folder (live files)
                # Skip if already have a valid daily file OR a non-empty day folder
                if dest_daily.exists() and dest_daily.stat().st_size > 0:
                    skipped += 1
                    continue
                if any(dest_dir.glob("*.nc")):
                    skipped += 1
                    continue
                fs = byday[ymd]
                daily = [f for f in fs if "live" not in f["filename"].lower()]
                if daily:
                    # Site has a pre-concatenated daily file — download it directly to YYYYMMDD.nc
                    f = daily[0]
                    tmp = str(dest_daily) + ".tmp"
                    try:
                        urllib.request.urlretrieve(f["downloadUrl"], tmp)
                        Path(tmp).rename(dest_daily)
                        tot_dl += 1
                        tot_mb += dest_daily.stat().st_size / 1e6
                        new_days += 1
                    except Exception as e:
                        print(f"  dl err {f['filename']}: {e}")
                        Path(tmp).unlink(missing_ok=True)
                else:
                    # Live files only (Hyytiala) — download subsampled into a day folder;
                    # homogenize_cl61_daily.py will concatenate them into a daily file.
                    live = sorted(fs, key=lambda f: f["filename"])
                    to_get = live[::SUBSAMPLE_LIVE]
                    dest_dir.mkdir(parents=True, exist_ok=True)
                    got = 0
                    for f in to_get:
                        dest = dest_dir / f["filename"]
                        if dest.exists() and dest.stat().st_size > 0:
                            continue
                        try:
                            urllib.request.urlretrieve(f["downloadUrl"], dest)
                            tot_dl += 1
                            got += 1
                            tot_mb += dest.stat().st_size / 1e6
                        except Exception as e:
                            print(f"  dl err {f['filename']}: {e}")
                    if got:
                        new_days += 1
            print(f"  {site}->{folder} {y}-{m:02d}: {new_days} new days, {skipped} already present  "
                  f"(running total {tot_dl} files, {tot_mb/1e3:.1f} GB)")
    print(f"downloaded {tot_dl} files, {tot_mb/1e3:.1f} GB total")
    print("DOWNLOAD_LONGRUN_DONE")


if __name__ == "__main__":
    main()
