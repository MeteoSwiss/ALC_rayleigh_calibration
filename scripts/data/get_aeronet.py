#!/usr/bin/env python3
"""
Fetch AERONET sun-photometer AOD from the AERONET v3 web API.

Python port of get_aeronet_from_api.m (M. Hervo). Uses the print_web_data_v3 endpoint:
  https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site=<SITE>&year=..&AOD15=1&AVG=10

Returns a list of dict rows with a parsed UTC datetime and the AOD columns. Helper
aod_at_532() interpolates to 532 nm (the miniMPL wavelength) via the Angstrom exponent
from the 500/440 nm pair when 532 is absent.

Usage (CLI):  python get_aeronet.py --site Toulouse --start 20250110 --end 20250116
"""
from __future__ import annotations

import argparse
import csv
import io
import urllib.request
from datetime import datetime

BASE = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"


def fetch_aeronet(site: str, start: str, end: str, level: int = 15) -> list[dict]:
    """site: AERONET site name; start/end: YYYYMMDD; level: 15 (1.5) or 20 (2.0)."""
    s, e = datetime.strptime(start, "%Y%m%d"), datetime.strptime(end, "%Y%m%d")
    q = (f"{BASE}?site={site}"
         f"&year={s.year}&month={s.month}&day={s.day}"
         f"&year2={e.year}&month2={e.month}&day2={e.day}"
         f"&AOD{level}=1&AVG=10")
    with urllib.request.urlopen(q, timeout=120) as r:
        text = r.read().decode("utf-8", errors="ignore")

    # The CSV header is the line beginning with "AERONET_Site,Date(dd:mm:yyyy),..."
    lines = text.replace("<br>", "\n").splitlines()
    hdr_i = next((i for i, ln in enumerate(lines)
                  if ln.startswith("AERONET_Site,") and "Date" in ln), None)
    if hdr_i is None:
        return []
    reader = csv.DictReader(io.StringIO("\n".join(lines[hdr_i:])))
    out = []
    for row in reader:
        date_k = next((k for k in row if k.startswith("Date(")), None)
        time_k = next((k for k in row if k.startswith("Time(")), None)
        if not date_k or not row.get(date_k) or ":" not in str(row[date_k]):
            continue
        try:
            dt = datetime.strptime(f"{row[date_k]} {row[time_k]}", "%d:%m:%Y %H:%M:%S")
        except (ValueError, TypeError):
            continue

        def num(key):
            for k in row:
                if k.replace("-", "_").startswith(key):
                    try:
                        v = float(row[k])
                        return None if v <= -999 else v
                    except (ValueError, TypeError):
                        return None
            return None

        out.append({
            "dt": dt,
            "AOD_500": num("AOD_500"), "AOD_440": num("AOD_440"),
            "AOD_532": num("AOD_532"), "AOD_675": num("AOD_675"),
            "AOD_1020": num("AOD_1020"),
            "angstrom_440_870": num("440-870_Angstrom") or num("Angstrom_Exponent"),
        })
    return out


def aod_at_532(row: dict) -> float | None:
    """Best estimate of AOD at 532 nm for one AERONET record."""
    if row.get("AOD_532") is not None:
        return row["AOD_532"]
    a500, a440 = row.get("AOD_500"), row.get("AOD_440")
    ae = row.get("angstrom_440_870")
    if a500 is not None and ae is not None:
        return a500 * (532.0 / 500.0) ** (-ae)          # Angstrom interpolation
    if a500 is not None and a440 is not None and a440 > 0 and a500 > 0:
        ae2 = -__import__("math").log(a500 / a440) / __import__("math").log(500.0 / 440.0)
        return a500 * (532.0 / 500.0) ** (-ae2)
    return a500


def nightly_aod_532(site: str, date_yyyymmdd: str) -> float | None:
    """Median 532 nm AOD over the day-of and day-after-dawn closest to a night.

    AERONET only measures in daylight, so for the night labelled D (spanning D-1 evening
    to D morning) we take the photometer points from the daytime of D-1 and D and use the
    median, mirroring the report's per-night comparison.
    """
    from datetime import timedelta
    d = datetime.strptime(date_yyyymmdd, "%Y%m%d")
    start = (d - timedelta(days=1)).strftime("%Y%m%d")
    rows = fetch_aeronet(site, start, date_yyyymmdd, level=15)
    vals = [aod_at_532(r) for r in rows]
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    import statistics
    return statistics.median(vals)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default="Toulouse")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    args = ap.parse_args()
    rows = fetch_aeronet(args.site, args.start, args.end)
    print(f"{len(rows)} AERONET records for {args.site} {args.start}..{args.end}")
    for r in rows[:8]:
        print(f"  {r['dt']}  AOD500={r['AOD_500']}  AOD532~{aod_at_532(r)}  AE={r['angstrom_440_870']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
