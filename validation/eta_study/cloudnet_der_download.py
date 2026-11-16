"""Inventory and download of the Cloudnet droplet effective radius (`der`) product, all sites.

Companion of cloudnet_sites_droplet_size.py (eta study, extension of step 1d to more Cloudnet sites).

1. Query the Cloudnet portal API (https://cloudnet.fmi.fi/api/files?product=der&dateFrom=..&dateTo=..)
   once for ALL sites (latest version of every file only, which is the API default).
2. Write the availability per site (files, MB, first/last day, coordinates) and the full file list to
   C:/DATA/Projects/202611_eta_study/cloudnet_sites/.
3. Download what is not on disk yet into D:/Cloudnet_Droplet_radius/<site>/<file>, verifying the
   SHA-256 checksum the portal publishes (written to <file>.part first, renamed when verified).

Files of the July collection 952a7e80 (Payerne, Palaiseau, Lindenberg to 2026-06-30, Chilbolton
Jan-Mar 2025) sit flat in D:/Cloudnet_Droplet_radius/ and are never moved or overwritten. A flat file
whose checksum equals the portal's counts as on disk; one the portal has reprocessed since (newer
version, different checksum) is downloaded again into the per-site folder, so that every site is
analysed with the latest portal processing while the July files stay as step 1d read them.

Data licence: ACTRIS Cloudnet, CC-BY 4.0 (cite the portal and the site PIs; see the README of the
collection in D:/Cloudnet_Droplet_radius/).

Run:  python validation/eta_study/cloudnet_der_download.py [--dry-run] [--threads 8]
                                                          [--sites a,b,c] [--from 2025-01-01 --to 2026-10-31]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time as time_mod
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://cloudnet.fmi.fi/api/files"
ROOT = Path("D:/Cloudnet_Droplet_radius")
OUT = Path("C:/DATA/Projects/202611_eta_study/cloudnet_sites")


def query_files(date_from: str, date_to: str) -> list[dict]:
    """All `der` files in the date range (every site), latest version of each."""
    url = f"{API}?product=der&dateFrom={date_from}&dateTo={date_to}"
    with urllib.request.urlopen(url, timeout=600) as r:
        return json.loads(r.read().decode("utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def local_status(f: dict) -> str:
    """'site' (per-site copy, verified at download), 'flat_same' (July collection, identical to the
    portal), 'flat_outdated' (July collection, portal has a newer version) or 'missing'."""
    if (ROOT / f["site"]["id"] / f["filename"]).exists():
        return "site"
    flat = ROOT / f["filename"]
    if flat.exists():
        return "flat_same" if sha256(flat) == f["checksum"] else "flat_outdated"
    return "missing"


def download_one(f: dict, retries: int = 3) -> tuple[str, str]:
    """Download one file to ROOT/<site>/<filename>, verified against the portal SHA-256."""
    dest = ROOT / f["site"]["id"] / f["filename"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    err = ""
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(f["downloadUrl"], timeout=300) as r, open(part, "wb") as out:
                while True:
                    block = r.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
            if sha256(part) != f["checksum"]:
                err = "checksum mismatch"
                continue
            part.replace(dest)
            return f["filename"], "ok"
        except Exception as exc:  # noqa: BLE001 - network hiccups are retried
            err = f"{type(exc).__name__}: {exc}"
            time_mod.sleep(2 * attempt)
    part.unlink(missing_ok=True)
    return f["filename"], f"FAILED ({err})"


def write_inventory(files: list[dict], path_sites: Path, path_files: Path) -> None:
    agg = defaultdict(lambda: dict(n=0, bytes=0, first="9999", last="0000"))
    for f in files:
        a = agg[f["site"]["id"]]
        s = f["site"]
        a.update(name=s["humanReadableName"], country=s["countryCode"], lat=s["latitude"],
                 lon=s["longitude"], alt=s["altitude"], site_type="+".join(s["type"]))
        a["n"] += 1
        a["bytes"] += int(f["size"])
        a["first"] = min(a["first"], f["measurementDate"])
        a["last"] = max(a["last"], f["measurementDate"])
    path_sites.parent.mkdir(parents=True, exist_ok=True)
    with open(path_sites, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["site", "name", "country", "site_type", "lat", "lon", "alt_m", "n_files", "MB",
                    "first", "last"])
        for site, a in sorted(agg.items(), key=lambda kv: -kv[1]["n"]):
            w.writerow([site, a["name"], a["country"], a["site_type"], a["lat"], a["lon"], a["alt"],
                        a["n"], f"{a['bytes'] / 1e6:.1f}", a["first"], a["last"]])
    with open(path_files, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["site", "date", "filename", "bytes", "sha256", "uuid", "pid", "updatedAt",
                    "volatile", "downloadUrl"])
        for f in sorted(files, key=lambda f: (f["site"]["id"], f["measurementDate"])):
            w.writerow([f["site"]["id"], f["measurementDate"], f["filename"], f["size"], f["checksum"],
                        f["uuid"], f["pid"], f["updatedAt"], f["volatile"], f["downloadUrl"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="date_from", default="2025-01-01")
    ap.add_argument("--to", dest="date_to", default="2026-10-31")
    ap.add_argument("--sites", default="", help="comma-separated Cloudnet site ids (default: all)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true", help="inventory + volume only, no download")
    a = ap.parse_args()

    files = query_files(a.date_from, a.date_to)
    print(f"portal: {len(files)} der files, {sum(int(f['size']) for f in files) / 1e9:.2f} GB, "
          f"{len({f['site']['id'] for f in files})} sites ({a.date_from}..{a.date_to})", flush=True)
    write_inventory(files, OUT / "der_availability.csv", OUT / "der_filelist.csv")
    if a.sites:
        keep = {s.strip() for s in a.sites.split(",") if s.strip()}
        files = [f for f in files if f["site"]["id"] in keep]

    todo, n_site, flat_same, flat_newer = [], 0, 0, []
    for f in files:
        st = local_status(f)
        if st == "site":
            n_site += 1
        elif st == "flat_same":
            flat_same += 1
        else:
            todo.append(f)
            if st == "flat_outdated":          # the July file stays; the newer version goes per-site
                flat_newer.append(f["filename"])
    per_site = defaultdict(lambda: [0, 0])
    for f in todo:
        per_site[f["site"]["id"]][0] += 1
        per_site[f["site"]["id"]][1] += int(f["size"])
    vol = sum(v[1] for v in per_site.values())
    print(f"on disk: {n_site} per-site files, {flat_same} July-collection files identical to the "
          f"portal; {len(flat_newer)} July-collection files have a newer portal version (re-downloaded "
          f"per-site, the July file is kept)", flush=True)
    print(f"to download: {len(todo)} files, {vol / 1e9:.2f} GB into {ROOT}/<site>/", flush=True)
    for site, (n, b) in sorted(per_site.items(), key=lambda kv: -kv[1][1]):
        print(f"  {site:22s} {n:4d} files {b / 1e6:7.1f} MB", flush=True)
    if flat_newer:
        (OUT / "july_collection_outdated.txt").write_text("\n".join(sorted(flat_newer)) + "\n",
                                                          encoding="utf-8")
    if a.dry_run or not todo:
        return

    t0 = time_mod.time()
    done, failed, got = 0, [], 0
    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        futs = {ex.submit(download_one, f): f for f in todo}
        for fut in as_completed(futs):
            name, status = fut.result()
            done += 1
            if status == "ok":
                got += int(futs[fut]["size"])
            else:
                failed.append(f"{name} {status}")
            if done % 500 == 0 or done == len(todo):
                dt = time_mod.time() - t0
                print(f"  {done}/{len(todo)} files, {got / 1e9:.2f} GB, {dt:.0f} s "
                      f"({got / 1e6 / max(dt, 1e-3):.1f} MB/s), {len(failed)} failed", flush=True)
    manifest = ROOT / f"download_manifest_{time_mod.strftime('%Y%m%d')}.csv"
    with open(manifest, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["site", "date", "filename", "bytes", "sha256", "pid", "updatedAt"])
        for f in sorted(todo, key=lambda f: (f["site"]["id"], f["measurementDate"])):
            if (ROOT / f["site"]["id"] / f["filename"]).exists():
                w.writerow([f["site"]["id"], f["measurementDate"], f["filename"], f["size"],
                            f["checksum"], f["pid"], f["updatedAt"]])
    print(f"manifest -> {manifest}", flush=True)
    if failed:
        (OUT / "download_failed.txt").write_text("\n".join(failed) + "\n", encoding="utf-8")
        print(f"{len(failed)} FAILED -> {OUT / 'download_failed.txt'} (re-run to retry)", flush=True)


if __name__ == "__main__":
    main()
