"""Eta study: in-situ accumulation-mode particle number (N50, N100) near the Cloudnet sites, from EBAS.

An independent, in-situ test of whether the Cloudnet N_scaled (and through it a_G and alpha) follows the aerosol.
N100 (dry particles larger than 100 nm) is the usual proxy of the cloud condensation nuclei active at the
0.2-0.3 % supersaturations of stratocumulus bases; N50 is given for lower activation diameters.

Data: the public EBAS THREDDS server (https://thredds.nilu.no/thredds/catalog/ebas/catalog.html, no login),
particle number size distributions (DMPS / SMPS / TDMPS / TSMPS / MPSS, level 2, matrix aerosol / pm10 / pm1;
APS, optical counters and volatility (non-volatile) matrices excluded). Most level-2 records on that server end
between 2006 and 2024, so this is a CLIMATOLOGY of the station, not the 2025-2026 period of the Cloudnet data:
for each station the most recent file of at least 26 weeks is used, at most its last 3 years.

Matching: the EBAS station nearest to each Cloudnet site of the E-PROFILE box within 100 km and within 1000 m of
altitude (mountain observatories above valley sites excluded). N50 / N100 are integrated per hourly spectrum
(dN/dlogD over the bins above the cut, a spectrum counting only if every bin above the cut is valid), then the
median over the hours. A station whose median hourly dN/dlogD maximum is below 100 cm-3 is flagged
(not a dN/dlogD record) and not used.

Run:  python validation/eta_study/cloudnet_sites_ebas_n100.py
Outputs (C:/DATA/Projects/202611_eta_study/cloudnet_sites/): external/ebas/ebas_catalog.xml,
external/ebas/ebas_pnsd_stations.csv, external/ebas/n100_<station>.csv (hourly N50/N100), ebas_n100_sites.csv.
No personal metadata (creator names, e-mails) is stored.
"""
from __future__ import annotations

import re
import sys
import time as time_mod
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cloudnet_sites_droplet_size as dsz  # noqa: E402

OUT = dsz.OUT
EBAS = OUT / "external" / "ebas"
CATALOG_URL = "https://thredds.nilu.no/thredds/catalog/ebas/catalog.xml"
DODS = "https://thredds.nilu.no/thredds/dodsC/ebas/"
INSTRUMENTS = ("dmps", "smps", "tdmps", "tsmps", "mpss")
MATRICES = ("aerosol", "pm10", "pm1", "pm25", "")
MIN_WEEKS = 26
MAX_YEARS = 3
MAX_DIST_KM = 100.0
MAX_DALT_M = 1000.0
MIN_PEAK = 100.0           # cm-3: median hourly dN/dlogD maximum; the cleanest station here (Arctic) has 214


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2.0 * 6371.0 * np.arcsin(np.sqrt(a))


def duration_days(code: str) -> float:
    m = re.fullmatch(r"(\d+)(y|mo|w|d|h)", code)
    if not m:
        return 0.0
    n, u = int(m.group(1)), m.group(2)
    return n * {"y": 365.25, "mo": 30.44, "w": 7.0, "d": 1.0, "h": 1 / 24}[u]


def catalog_files():
    """Level-2 particle number size distribution files of the EBAS THREDDS catalogue."""
    import pandas as pd
    import requests
    fp = EBAS / "ebas_catalog.xml"
    if not fp.exists():
        EBAS.mkdir(parents=True, exist_ok=True)
        fp.write_bytes(requests.get(CATALOG_URL, timeout=600).content)
    names = re.findall(r'<dataset name="([^"]+\.nc)"', fp.read_text(encoding="utf-8"))
    rows = []
    for n in names:
        f = n.split(".")
        if len(f) < 12 or f[4] != "particle_number_size_distribution" or f[3] not in INSTRUMENTS:
            continue
        if f[5] not in MATRICES or "lev2" not in f[-2]:
            continue
        start = pd.Timestamp(f[1][:8])
        dur = duration_days(f[6])
        rows.append(dict(station=f[0], name=n, instrument=f[3], matrix=f[5], start=start,
                         end=start + pd.Timedelta(days=dur), weeks=dur / 7.0))
    return pd.DataFrame(rows)


def station_coordinates(files):
    """Latitude, longitude, altitude of every station (geospatial attributes of one file, via OPeNDAP)."""
    import pandas as pd
    from netCDF4 import Dataset
    fp = EBAS / "ebas_pnsd_stations.csv"
    known = pd.read_csv(fp) if fp.exists() else pd.DataFrame(columns=["station", "lat", "lon", "alt_m"])
    rows = known.to_dict("records")
    for st, g in files.groupby("station"):
        if st in set(known["station"]):
            continue
        try:
            with Dataset(DODS + g["name"].iloc[0]) as nc:
                rows.append(dict(station=st, lat=float(nc.geospatial_lat_min), lon=float(nc.geospatial_lon_min),
                                 alt_m=float(getattr(nc, "geospatial_vertical_min", np.nan))))
        except Exception as exc:  # noqa: BLE001 - one unreachable file must not stop the inventory
            print(f"  {st}: {type(exc).__name__}: {exc}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(fp, index=False)
    return out


def pnsd_variable(nc) -> str:
    """The mean size distribution of the file: '..._amean' where statistics are given, else the plain one
    (never the percentiles, precision, flags or metadata companions)."""
    cands = [k for k, v in nc.variables.items() if k.startswith("particle_number_size_distribution")
             and v.dimensions[:1] == ("D",) and v.dimensions[-1:] == ("time",) and v.ndim == 2]
    for pref in ("particle_number_size_distribution_amean", "particle_number_size_distribution"):
        if pref in cands:
            return pref
    plain = [k for k in cands if not re.search(r"(perc|prec|_qc|metadata)", k)]
    if not plain:
        raise KeyError(f"no size-distribution variable among {list(nc.variables)}")
    return plain[0]


def hourly_n(name: str, station: str):
    """Hourly N50 / N100 (cm-3) of one EBAS file, last MAX_YEARS years only; cached as CSV. None when the file
    does not resolve the accumulation mode (fewer than 10 size bins, or no bin above 300 nm)."""
    import pandas as pd
    from netCDF4 import Dataset, num2date
    fp = EBAS / f"n100_{station}.csv"
    if fp.exists():
        return pd.read_csv(fp, parse_dates=["time"])
    with Dataset(DODS + name) as nc:
        if nc.variables["D"].size < 10 or float(np.max(nc.variables["D"][:])) < 300.0:
            return None
        tv = nc.variables["time"]
        t = np.asarray(tv[:], float)
        i0 = int(np.searchsorted(t, t[-1] - MAX_YEARS * 365.25))
        D = np.asarray(nc.variables["D"][:], float)                        # nm
        x = np.ma.filled(np.ma.asarray(nc.variables[pnsd_variable(nc)][:, i0:], dtype=float),
                         np.nan)                                           # dN/dlogD, cm-3, (D, time)
        when = num2date(t[i0:], tv.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    lD = np.log10(D)
    edges = np.r_[lD[0] - (lD[1] - lD[0]) / 2, (lD[1:] + lD[:-1]) / 2, lD[-1] + (lD[-1] - lD[-2]) / 2]
    dlog = np.diff(edges)
    out = {"time": pd.to_datetime(when)}
    for cut in (50, 100):
        m = D >= cut
        n = np.sum(x[m] * dlog[m][:, None], axis=0)                       # NaN if any bin above the cut is NaN
        out[f"N{cut}"] = n
    with np.errstate(all="ignore"):
        out["peak"] = np.nanmax(np.where(np.isfinite(x), x, -np.inf), axis=0)   # dN/dlogD maximum of the hour
    out["peak"] = np.where(np.isfinite(out["peak"]), out["peak"], np.nan)
    df = pd.DataFrame(out)
    df.attrs["D_max_nm"] = float(D.max())
    df.to_csv(fp, index=False, float_format="%.5g")
    return df


def main():
    import pandas as pd
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time_mod.time()
    files = catalog_files()
    coords = station_coordinates(files)
    tab = pd.read_csv(OUT / "site_table_extinction.csv")
    tab = tab[tab["eprofile_box"]]
    rows = []
    for _, r in tab.iterrows():
        d = haversine_km(r["lat"], r["lon"], coords["lat"].to_numpy(float), coords["lon"].to_numpy(float))
        dalt = np.abs(coords["alt_m"].to_numpy(float) - r["alt_m"])
        ok = (d <= MAX_DIST_KM) & (dalt <= MAX_DALT_M)
        row = dict(site=r["site"])
        if ok.any():
            i = int(np.argmin(np.where(ok, d, np.inf)))
            st = coords.iloc[i]
            g = files[(files["station"] == st["station"]) & (files["weeks"] >= MIN_WEEKS)].sort_values("end")
            f, h = None, None
            for _, cand in g.iloc[::-1].iterrows():          # the most recent usable file
                h = hourly_n(cand["name"], st["station"])
                if h is not None:
                    f = cand
                    break
            if f is not None:
                v = h["N100"].to_numpy(float)
                ok_h = np.isfinite(v)
                row.update(ebas_station=st["station"], ebas_lat=st["lat"], ebas_lon=st["lon"], ebas_alt_m=st["alt_m"],
                           ebas_dist_km=float(d[i]), ebas_instrument=f["instrument"], ebas_matrix=f["matrix"],
                           ebas_period=f"{h['time'][ok_h].min():%Y-%m}..{h['time'][ok_h].max():%Y-%m}" if ok_h.any() else "",
                           ebas_n_hours=int(ok_h.sum()), N100_med=float(np.nanmedian(v)) if ok_h.any() else np.nan,
                           N100_mean=float(np.nanmean(v)) if ok_h.any() else np.nan,
                           N50_med=float(np.nanmedian(h["N50"])) if np.isfinite(h["N50"]).any() else np.nan,
                           ebas_peak_dNdlogD_med=float(np.nanmedian(h["peak"])),
                           ebas_file=f["name"])
                # plausibility: a dN/dlogD spectrum peaks at >= 100 cm-3 even in the Arctic; a far lower peak
                # means another convention (e.g. dN per bin) -> the station is flagged and not used
                row["ebas_qc"] = "ok" if row["ebas_peak_dNdlogD_med"] >= MIN_PEAK else "implausible spectrum level"
        rows.append(row)
        print(f"  {r['site']:22s} {row.get('ebas_station', '-'):8s} {row.get('ebas_dist_km', np.nan):6.1f} km "
              f"{row.get('ebas_period', ''):16s} N100 {row.get('N100_med', np.nan):7.0f}  N50 "
              f"{row.get('N50_med', np.nan):7.0f}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "ebas_n100_sites.csv", index=False, float_format="%.6g")
    print(f"total runtime {time_mod.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
