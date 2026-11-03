"""Eta study: what explains the differences of droplet size (a_G) and extinction (alpha) between Cloudnet sites?

Inputs (cloudnet_sites_droplet_size.py, cloudnet_sites_extinction.py, cloudnet_sites_ebas_n100.py;
C:/DATA/Projects/202611_eta_study/cloudnet_sites/): site_table_extinction.csv (one row per site: a_G, alpha, N_scaled
site medians, liquid-layer depth, base height, retrieval-status fractions), site_table.csv (radar, band, LWP source,
a_G of the 50 / 150 m layers), site_radar_table.csv, season tables, extinction_scenes/<site>.npz, ebas_n100_sites.csv.

Candidate explanations, one column each of site_explain_table.csv:

  aerosol, downloaded without login (raw files cached under cloudnet_sites/external/, re-used when present):
    * AERONET version 3 direct-sun daily averages, level 1.5 (cloud-screened, quality-controlled; level 2.0 is
      counted and compared but not used: many 2025-2026 records are not yet promoted), 2025-01-01 .. 2026-10-31, at
      the AERONET site nearest to the Cloudnet site within 100 km (several within 10 km: the longest record; a
      photometer within 500 m of the site altitude preferred): AOD 500 nm, Angstrom exponent 440-870 nm, aerosol
      index AI = AOD x AE (Nakajima et al. 2001), and the SDA fine-mode AOD at 500 nm (O'Neill et al. 2003).
      https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3
    * surface PM2.5 and PM10 from the EEA air-quality download service (https://eeadmz1-downloads-api-appservice.
      azurewebsites.net/), verified E1a data for 2025 completed with the up-to-date E2a data for 2026: mean of the
      daily means 2025-01-01 .. 2026-10-31 at the nearest BACKGROUND sampling point within 100 km (within 500 m of
      the site altitude first, then rural or suburban before urban), station metadata from
      https://discomap.eea.europa.eu/map/fme/metadata/PanEuropean_metadata.csv. A station more than 500 m below or
      above the site is flagged and not used.
    * in-situ N50 / N100 climatology from EBAS (cloudnet_sites_ebas_n100.py), merged when present.
  thermodynamics (local CAMS 0.4 deg model levels, the calibration's water-vapour input):
    * cloud-base temperature of every scene (CAMS temperature at the liquid-base altitude, nearest 3-hourly
      step, D:/CAMS_Monthly_04 for 2025-01 .. 2026-06 and D:/CAMS_daily for 2026-07 .. 2026-10), and the adiabatic
      liquid-water lapse rate Gamma_ad(T_base, p_base) derived from it
  cloud and retrieval (in the der files):
    * liquid-layer depth and base height above the site; the share of window scenes whose lowest-100-m pixels carry
      retrieval status 3 (rain flag of the profile) or 2 (droplet pixel that also carries the categorize 'falling'
      bit: drizzle, ice, or radar echo continuing above the lidar-detected liquid top; CloudnetPy der.py /
      falling.py); the near-base gradient a_G(50 m) / a_G(150 m); the liquid water content that the retrieval puts
      in the lowest 100 m and its ratio to the adiabatic value (retrieval_consistency); radar band and LWP source
      (multi-channel radiometer vs the 89-GHz channel of the RPG radar); site altitude, |latitude|, record length
      and season mix.
  instrument tests:
    * sites whose radar changed (radar_changes: a_G, N, alpha per radar, against the reflectivity-offset signature
      a_G ~ k^(1/3), N ~ 1/k, alpha ~ k^(-1/3));
    * neighbouring sites 4-7 km apart on their common days (neighbour_pairs): same cloud field, different instruments.

Statistics: Spearman rank correlations of the site medians of a_G, alpha and N_scaled with every candidate, on
the 30 Cloudnet sites inside the E-PROFILE box (all have >= 20 days), and three sensitivity sets.

Run:  python validation/eta_study/cloudnet_sites_explain.py [--workers 12] [--offline] [--plot-only]
Outputs (C:/DATA/Projects/202611_eta_study/cloudnet_sites/): external/aeronet/*, external/eea/*, cams_tbase/<site>.npz,
aeronet_sites.csv, eea_pm_sites.csv, cams_tbase_sites.csv, retrieval_consistency_sites.csv,
season_ratio_extinction.csv, radar_change_extinction.csv, neighbour_pairs_extinction.csv, site_explain_table.csv,
correlations_explain.csv, explain_tables.md; doc/reports/figs_04v2/fig_cloudnet_sites_vs_aerosol.png.
Read-only on every other input.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import io
import re
import sys
import time as time_mod
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cloudnet_sites_droplet_size as dsz  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = dsz.OUT
EXT = OUT / "external"
SCENE_DIR = OUT / "extinction_scenes"
TBASE_DIR = OUT / "cams_tbase"
FIG_DIR = REPO / "doc" / "reports" / "figs_04v2"

PERIOD = ("2025-01-01", "2026-10-31")
N_PERIOD_DAYS = 669
MAX_DIST_KM = 100.0
COLOCATED_KM = 10.0                 # several AERONET sites this close: take the longest record
MIN_AERONET_DAYS = 30
MAX_DELEV_M = 500.0                 # AERONET / PM station within 500 m of the site altitude preferred
MIN_PM_FRAC = 0.4                   # PM: at least 40 % of the period's days with a valid daily mean
MIN_HOURS_PER_DAY = 18              # EEA rule for a daily mean from hourly values (75 %)

AERONET_URL = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
AERONET_SITES_URL = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
EEA_API = "https://eeadmz1-downloads-api-appservice.azurewebsites.net"
EEA_META_URL = "https://discomap.eea.europa.eu/map/fme/metadata/PanEuropean_metadata.csv"

CAMS_MONTHLY = Path("D:/CAMS_Monthly_04")
CAMS_DAILY = Path("D:/CAMS_daily")
CAMS_EPOCH_OFFSET_H = 613608.0      # hours from 1900-01-01 to 1970-01-01
CAMS_MAX_DT_H = 1.6                 # 3-hourly steps: nearest within 1.6 h
LAPSE_BELOW = 6.5e-3                # K/m, below the lowest model level (base under the model orography)


# =============================================================================
#  Small helpers
# =============================================================================
def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2.0 * 6371.0 * np.arcsin(np.sqrt(a))


def http_get(url, params=None, timeout=180, tries=3):
    import requests
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            last = f"HTTP {r.status_code}"
        except Exception as exc:  # noqa: BLE001 - network hiccup: retry
            last = f"{type(exc).__name__}: {exc}"
        time_mod.sleep(3 * (k + 1))
    raise RuntimeError(f"GET {url} failed: {last}")


def _q(x, p):
    return dsz._q(x, p)


# =============================================================================
#  AERONET
# =============================================================================
def aeronet_locations(offline: bool):
    import pandas as pd
    fp = EXT / "aeronet" / "aeronet_locations_v3.txt"
    if not fp.exists():
        if offline:
            raise FileNotFoundError(fp)
        fp.parent.mkdir(parents=True, exist_ok=True)
        fp.write_bytes(http_get(AERONET_SITES_URL).content)
    loc = pd.read_csv(fp, skiprows=1)
    loc.columns = ["aeronet_site", "lon", "lat", "elev_m"]
    return loc.drop_duplicates("aeronet_site")


def aeronet_file(site: str, product: str, offline: bool) -> Path | None:
    """Daily averages (AVG=20) of `product` (AOD15, AOD20, SDA15) for the whole period; cached. The 'Contact'
    header line (PI name and e-mail) is not stored: no personal data in the archive."""
    fp = EXT / "aeronet" / f"{site}_{product}_daily.csv"
    if fp.exists():
        return fp
    if offline:
        return None
    y1, m1, d1 = (int(x) for x in PERIOD[0].split("-"))
    y2, m2, d2 = (int(x) for x in PERIOD[1].split("-"))
    params = {"site": site, "year": y1, "month": m1, "day": d1, "year2": y2, "month2": m2, "day2": d2,
              product: 1, "AVG": 20, "if_no_html": 1}
    text = http_get(AERONET_URL, params=params).text
    lines = [ln for ln in text.splitlines() if not ln.startswith("Contact")]
    fp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    time_mod.sleep(0.5)                  # be gentle with the AERONET server
    return fp


def read_aeronet(fp: Path | None):
    import pandas as pd
    if fp is None or not fp.exists():
        return None
    lines = fp.read_text(encoding="utf-8", errors="replace").splitlines()
    hdr = [k for k, ln in enumerate(lines) if ln.startswith("AERONET_Site")]
    if not hdr or len(lines) <= hdr[0] + 1:
        return None
    df = pd.read_csv(io.StringIO("\n".join(lines[hdr[0]:])), low_memory=False)
    df = df.replace([-999.0, -999], np.nan)
    dcol = [c for c in df.columns if c.startswith("Date")][0]
    df["date"] = pd.to_datetime(df[dcol], format="%d:%m:%Y", errors="coerce")
    return df


def aod500(df):
    """AOD at 500 nm; where the instrument has no 500-nm channel, from 440 nm and the 440-675 Angstrom exponent."""
    a = df.get("AOD_500nm", pd_nan(df))
    fill = df.get("AOD_440nm", pd_nan(df)) * (500.0 / 440.0) ** (-df.get("440-675_Angstrom_Exponent", pd_nan(df)))
    return a.where(np.isfinite(a), fill)


def pd_nan(df):
    import pandas as pd
    return pd.Series(np.nan, index=df.index)


def aeronet_for_sites(tab, offline: bool):
    """Per Cloudnet site: the AERONET site used and its statistics over the period (daily level-1.5 values)."""
    import pandas as pd
    loc = aeronet_locations(offline)
    rows = []
    for _, r in tab.iterrows():
        d = haversine_km(r["lat"], r["lon"], loc["lat"].to_numpy(), loc["lon"].to_numpy())
        cand = loc.assign(dist_km=d)[d <= MAX_DIST_KM].sort_values("dist_km")
        row = dict(site=r["site"], aeronet_site="", aeronet_dist_km=np.nan, aeronet_elev_m=np.nan,
                   aeronet_n_candidates=len(cand))
        chosen = None
        # co-located candidates (<= 10 km): the longest record; otherwise the nearest with enough days. A
        # candidate within MAX_DELEV_M of the site altitude is preferred (a valley photometer below a mountain
        # site adds the boundary layer under the station to the column)
        cand = cand.assign(far_elev=(np.abs(cand["elev_m"] - r["alt_m"]) > MAX_DELEV_M).astype(int))
        cand = cand.sort_values(["far_elev", "dist_km"])
        near = cand[cand["dist_km"] <= COLOCATED_KM]
        if (near["far_elev"] == 0).any():
            near = near[near["far_elev"] == 0]
        best_n = -1
        for _, c in near.iterrows():
            df = read_aeronet(aeronet_file(c["aeronet_site"], "AOD15", offline))
            n = 0 if df is None else int(np.isfinite(aod500(df)).sum())
            if n >= MIN_AERONET_DAYS and n > best_n:
                chosen, best_n = (c, df), n
        if chosen is None:
            for _, c in cand[cand["dist_km"] > COLOCATED_KM].iterrows():
                df = read_aeronet(aeronet_file(c["aeronet_site"], "AOD15", offline))
                n = 0 if df is None else int(np.isfinite(aod500(df)).sum())
                if n >= MIN_AERONET_DAYS:
                    chosen = (c, df)
                    break
        if chosen is not None:
            c, df = chosen
            a5 = aod500(df)
            ae = df.get("440-870_Angstrom_Exponent", pd_nan(df))
            ok = np.isfinite(a5)
            row.update(aeronet_site=c["aeronet_site"], aeronet_dist_km=float(c["dist_km"]),
                       aeronet_elev_m=float(c["elev_m"]), aod_n_days=int(ok.sum()),
                       aod500_med=_q(a5, 50), aod500_mean=float(np.nanmean(a5)), aod500_p75=_q(a5, 75),
                       ae_med=_q(ae, 50), ai_med=_q(a5 * ae, 50),
                       aod_frac_days_dust=float(np.mean((a5[ok] > 0.2) & (ae[ok] < 0.6))))
            l2 = read_aeronet(aeronet_file(c["aeronet_site"], "AOD20", offline))
            row["aod_l2_n_days"] = 0 if l2 is None else int(np.isfinite(aod500(l2)).sum())
            if l2 is not None and np.isfinite(aod500(l2)).sum() >= MIN_AERONET_DAYS:
                row["aod500_l2_med"] = _q(aod500(l2), 50)
            sda = read_aeronet(aeronet_file(c["aeronet_site"], "SDA15", offline))
            if sda is not None:
                fine = sda.get("Fine_Mode_AOD_500nm[tau_f]", pd_nan(sda))
                row.update(sda_n_days=int(np.isfinite(fine).sum()), fine_aod500_med=_q(fine, 50),
                           fmf500_med=_q(sda.get("FineModeFraction_500nm[eta]", pd_nan(sda)), 50),
                           coarse_aod500_med=_q(sda.get("Coarse_Mode_AOD_500nm[tau_c]", pd_nan(sda)), 50))
        rows.append(row)
        print(f"  AERONET {r['site']:22s} -> {row['aeronet_site'] or '-':22s} "
              f"{row['aeronet_dist_km']:6.1f} km  AOD500 {row.get('aod500_med', np.nan):.3f}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "aeronet_sites.csv", index=False, float_format="%.6g")
    return out


# =============================================================================
#  EEA surface PM
# =============================================================================
RURAL_SUBURBAN = ("rural", "rural-regional", "rural-nearcity", "rural-remote", "suburban")
EEA_POLS = {"PM2.5": "6001", "PM10": "5"}


def eea_metadata(offline: bool):
    import pandas as pd
    fp = EXT / "eea" / "PanEuropean_metadata.csv"
    if not fp.exists():
        if offline:
            raise FileNotFoundError(fp)
        fp.parent.mkdir(parents=True, exist_ok=True)
        fp.write_bytes(http_get(EEA_META_URL, timeout=600).content)
    md = pd.read_csv(fp, sep="\t", low_memory=False)
    md["pol"] = md["AirPollutantCode"].astype(str).str.rsplit("/", n=1).str[-1]
    md = md[md["pol"].isin(EEA_POLS.values()) & (md["AirQualityStationType"] == "background")]
    end = md["ObservationDateEnd"].astype(str)
    md = md[md["ObservationDateEnd"].isna() | (end >= "2025-01-01")]
    return md.drop_duplicates("SamplingPoint")


def eea_url_list(country: str, pollutant: str, dataset: int, offline: bool) -> list[str]:
    import requests
    fp = EXT / "eea" / "urls" / f"{country}_{pollutant}_ds{dataset}.txt"
    if not fp.exists():
        if offline:
            return []
        fp.parent.mkdir(parents=True, exist_ok=True)
        body = {"countries": [country], "cities": [], "pollutants": [pollutant], "dataset": dataset,
                "dateTimeStart": f"{PERIOD[0]}T00:00:00Z", "dateTimeEnd": f"{PERIOD[1]}T23:59:59Z"}
        r = requests.post(EEA_API + "/ParquetFile/urls", json=body, timeout=300)
        r.raise_for_status()
        urls = [u.strip() for u in r.text.splitlines() if u.strip().startswith("http")]
        fp.write_text("\n".join(urls) + "\n", encoding="utf-8")
    return [u for u in fp.read_text(encoding="utf-8").splitlines() if u.strip()]


def eea_daily(spo: str, country: str, pollutant: str, offline: bool):
    """Daily means of one sampling point over the period: verified E1a (dataset 2) first, the up-to-date E2a
    (dataset 1) for the days E1a does not cover yet."""
    import pandas as pd
    series = []
    for ds in (2, 1):
        fname = re.sub(r"[^A-Za-z0-9._-]", "_", spo)        # the service writes ':' of the id as '_'
        url = next((u for u in eea_url_list(country, pollutant, ds, offline) if u.endswith(f"/{fname}.parquet")), None)
        if url is None:
            continue
        fp = EXT / "eea" / "parquet" / f"ds{ds}" / country / f"{fname}.parquet"
        if not fp.exists():
            if offline:
                continue
            fp.parent.mkdir(parents=True, exist_ok=True)
            fp.write_bytes(http_get(url, timeout=300).content)
        df = pd.read_parquet(fp)
        df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
        df = df[df["Validity"].isin([1, 2, 3]) & (df["Value"] >= 0)]
        df = df[(df["Start"] >= PERIOD[0]) & (df["Start"] < pd.Timestamp(PERIOD[1]) + pd.Timedelta(days=1))]
        if df.empty:
            continue
        day = df["Start"].dt.floor("D")
        hourly = df["AggType"].astype(str).str.lower().eq("hour")
        parts = []
        if hourly.any():
            g = df[hourly].groupby(day[hourly])["Value"]
            m = g.mean()[g.count() >= MIN_HOURS_PER_DAY]
            parts.append(m)
        daily = df["AggType"].astype(str).str.lower().eq("day")
        if daily.any():
            parts.append(df[daily].groupby(day[daily])["Value"].mean())
        if parts:
            s = pd.concat(parts)
            series.append(s[~s.index.duplicated()])
    if not series:
        return None
    merged = series[0]
    for s in series[1:]:
        merged = pd.concat([merged, s[~s.index.isin(merged.index)]])
    return merged.sort_index()


def eea_for_sites(tab, offline: bool):
    import pandas as pd
    md = eea_metadata(offline)
    rows = []
    for _, r in tab[tab["eprofile_box"]].iterrows():
        row = dict(site=r["site"])
        for pol, code in EEA_POLS.items():
            m = md[md["pol"] == code]
            d = haversine_km(r["lat"], r["lon"], m["Latitude"].to_numpy(float), m["Longitude"].to_numpy(float))
            c = m.assign(dist_km=d)[d <= MAX_DIST_KM].copy()
            c["area_rank"] = np.where(c["AirQualityStationArea"].isin(RURAL_SUBURBAN), 0, 1)
            c["far_alt"] = (np.abs(c["Altitude"].astype(float) - r["alt_m"]) > MAX_DELEV_M).astype(int)
            c = c.sort_values(["far_alt", "area_rank", "dist_km"])
            key = pol.replace(".", "")
            row[f"{key}_n_candidates"] = len(c)
            for _, s in c.iterrows():
                spo = str(s["SamplingPoint"])
                try:
                    daily = eea_daily(spo, str(s["Countrycode"]), pol, offline)
                except Exception as exc:  # noqa: BLE001 - one bad file must not stop the site loop
                    print(f"    EEA {spo}: {type(exc).__name__}: {exc}", flush=True)
                    daily = None
                if daily is None or len(daily) < MIN_PM_FRAC * N_PERIOD_DAYS:
                    continue
                row.update({f"{key}_station": str(s["AirQualityStationEoICode"]), f"{key}_spo": spo,
                            f"{key}_area": str(s["AirQualityStationArea"]), f"{key}_dist_km": float(s["dist_km"]),
                            f"{key}_alt_m": float(s["Altitude"]), f"{key}_n_days": int(len(daily)),
                            f"{key}_mean": float(daily.mean()), f"{key}_median": float(daily.median())})
                break
        rows.append(row)
        print(f"  EEA {r['site']:22s} PM2.5 {row.get('PM25_mean', np.nan):5.1f} ({row.get('PM25_station', '-')}, "
              f"{row.get('PM25_dist_km', np.nan):.0f} km)  PM10 {row.get('PM10_mean', np.nan):5.1f} "
              f"({row.get('PM10_station', '-')}, {row.get('PM10_dist_km', np.nan):.0f} km)", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "eea_pm_sites.csv", index=False, float_format="%.6g")
    return out


# =============================================================================
#  CAMS cloud-base temperature
# =============================================================================
def cams_files():
    """(path, keep-date or None): monthly 0.4 deg files 2025-01 .. 2026-06, then the daily files (which span
    D-1 03 UTC .. D+1 00 UTC) for 2026-07 .. 2026-10, keeping only their own date."""
    import pandas as pd
    jobs = [(str(p), None) for p in sorted(CAMS_MONTHLY.glob("CAMS_Beta_2*.nc"))]
    last_month = max(Path(p).stem[-6:] for p, _ in jobs)
    for p in sorted(CAMS_DAILY.glob("CAMS_Beta_2*.nc")):
        ds = p.stem[-8:]
        if ds[:6] > last_month and ds <= PERIOD[1].replace("-", ""):
            jobs.append((str(p), ds))
    return jobs


def cams_columns(job):
    """Temperature columns at the grid points of the sites from one CAMS file. The level altitudes are taken as
    the file's time-mean (they move by a few tens of metres at most, < 0.3 K at the cloud base)."""
    from netCDF4 import Dataset
    path, keep, points = job               # points: {site: (iy, ix)}
    out = {}
    with Dataset(path) as nc:
        t_h = np.asarray(nc.variables["time"][:], float)
        t_s = (t_h - CAMS_EPOCH_OFFSET_H) * 3600.0
        if keep is not None:
            day0 = np.datetime64(f"{keep[:4]}-{keep[4:6]}-{keep[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
            tsel = np.flatnonzero((t_s >= day0) & (t_s < day0 + 86400))
        else:
            tsel = np.arange(t_s.size)
        if tsel.size == 0:
            return out
        tv, zv = nc.variables["t"], nc.variables["altitude"]
        cy, cx = tv.chunking()[2], tv.chunking()[3]
        groups = {}
        for s, (iy, ix) in points.items():
            groups.setdefault((iy // cy, ix // cx), []).append((s, iy, ix))
        t0, t1 = int(tsel[0]), int(tsel[-1]) + 1
        tm = t_s[t0:t1]
        for members in groups.values():
            ys = [m[1] for m in members]
            xs = [m[2] for m in members]
            y0, y1, x0, x1 = min(ys), max(ys) + 1, min(xs), max(xs) + 1
            tt = np.ma.filled(np.ma.asarray(tv[t0:t1, :, y0:y1, x0:x1], dtype=float), np.nan)
            zz = np.ma.filled(np.ma.asarray(zv[t0:t1:4, :, y0:y1, x0:x1], dtype=float), np.nan)
            for s, iy, ix in members:
                out[s] = (tm, tt[:, :, iy - y0, ix - x0].astype(np.float32),
                          np.nanmean(zz[:, :, iy - y0, ix - x0], axis=0).astype(np.float32))
    return out


def gamma_ad(T, p):
    """Adiabatic liquid-water-content lapse rate (g m-3 per km) at temperature T (K) and pressure p (Pa)."""
    g, rd, rv, cp, lv, eps = 9.81, 287.05, 461.5, 1004.0, 2.501e6, 0.622
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))          # Bolton 1980, over water
    qs = eps * es / (p - (1 - eps) * es)
    gm = g * (1 + lv * qs / (rd * T)) / (cp + lv ** 2 * qs * eps / (rd * T ** 2))   # moist-adiabatic lapse rate
    rho = p / (rd * T * (1 + 0.61 * qs))
    return rho * qs * (lv * gm / (rv * T ** 2) - g / (rd * T)) * 1e3 * 1e3          # kg m-4 -> g m-3 km-1


def cams_tbase(tab, workers: int):
    """Per scene: CAMS temperature at the liquid-base altitude (box sites). Cached per site."""
    import pandas as pd
    from netCDF4 import Dataset
    box = tab[tab["eprofile_box"]]
    TBASE_DIR.mkdir(parents=True, exist_ok=True)
    todo = [s for s in box["site"] if not (TBASE_DIR / f"{s}.npz").exists()]
    if todo:
        jobs0 = cams_files()
        with Dataset(jobs0[0][0]) as nc:
            lat = np.asarray(nc.variables["latitude"][:], float)
            lon = np.asarray(nc.variables["longitude"][:], float)
        points = {}
        for _, r in box[box["site"].isin(todo)].iterrows():
            points[r["site"]] = (int(np.argmin(np.abs(lat - r["lat"]))), int(np.argmin(np.abs(lon - r["lon"]))))
        jobs = [(p, k, points) for p, k in jobs0]
        print(f"CAMS: {len(jobs)} files, {len(points)} sites -> {workers} processes", flush=True)
        cols = {s: [] for s in points}
        t0 = time_mod.time()
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for k, res in enumerate(ex.map(cams_columns, jobs), 1):
                for s, v in res.items():
                    cols[s].append(v)
                if k % 20 == 0:
                    print(f"  {k}/{len(jobs)} CAMS files ({time_mod.time() - t0:.0f} s)", flush=True)
        for s in points:
            tm = np.concatenate([c[0] for c in cols[s]])
            T = np.concatenate([c[1] for c in cols[s]])
            Z = np.concatenate([np.repeat(c[2][None, :], c[0].size, axis=0) for c in cols[s]])
            order = np.argsort(tm, kind="stable")
            tm, T, Z = tm[order], T[order], Z[order]
            keep = np.r_[True, np.diff(tm) > 0]
            tm, T, Z = tm[keep], T[keep], Z[keep]
            z = np.load(SCENE_DIR / f"{s}.npz")
            alt = float(box.loc[box["site"] == s, "alt_m"].iloc[0])
            zb = z["base_agl"].astype(float) + alt
            ts = z["t_epoch"].astype(float)
            idx = np.clip(np.searchsorted(tm, ts), 1, tm.size - 1)
            idx = np.where(np.abs(tm[idx - 1] - ts) <= np.abs(tm[idx] - ts), idx - 1, idx)
            ok_t = np.abs(tm[idx] - ts) <= CAMS_MAX_DT_H * 3600.0
            tb = np.full(ts.size, np.nan)
            for i in np.unique(idx[ok_t]):
                sel = ok_t & (idx == i)
                zc, tc = Z[i][::-1].astype(float), T[i][::-1].astype(float)     # ascending altitude
                good = np.isfinite(zc) & np.isfinite(tc)
                zc, tc = zc[good], tc[good]
                v = np.interp(zb[sel], zc, tc)
                below = zb[sel] < zc[0]
                v[below] = tc[0] + LAPSE_BELOW * (zc[0] - zb[sel][below])
                tb[sel] = v
            np.savez_compressed(TBASE_DIR / f"{s}.npz", t_epoch=z["t_epoch"], T_base=tb.astype(np.float32),
                                z_base_msl=zb.astype(np.float32))
    rows = []
    for _, r in box.iterrows():
        fp = TBASE_DIR / f"{r['site']}.npz"
        if not fp.exists():
            continue
        tb = np.load(fp)["T_base"].astype(float)
        zb = np.load(fp)["z_base_msl"].astype(float)
        sc = np.load(SCENE_DIR / f"{r['site']}.npz")
        ok = np.isfinite(sc["alpha_100"]) & np.isfinite(tb)
        p = 101325.0 * (1 - 2.25577e-5 * zb) ** 5.25588                     # standard atmosphere
        gad = gamma_ad(tb, p)
        rows.append(dict(site=r["site"], tbase_n=int(ok.sum()), tbase_med_C=_q(tb[ok] - 273.15, 50),
                         tbase_p25_C=_q(tb[ok] - 273.15, 25), tbase_p75_C=_q(tb[ok] - 273.15, 75),
                         frac_tbase_below0=float(np.mean(tb[ok] < 273.15)),
                         frac_tbase_below_m10=float(np.mean(tb[ok] < 263.15)),
                         gamma_ad_med=_q(gad[ok], 50), zbase_msl_med=_q(zb[ok], 50)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "cams_tbase_sites.csv", index=False, float_format="%.6g")
    return out


# =============================================================================
#  Assembly, correlations, figure
# =============================================================================
def season_mix(s: str):
    z = np.load(SCENE_DIR / f"{s}.npz")
    ok = np.isfinite(z["alpha_100"])
    month = (z["date"][ok] // 100) % 100
    return float(np.mean(np.isin(month, (12, 1, 2)))), float(np.mean(np.isin(month, (6, 7, 8))))


RHO_W = 1000.0          # kg m-3
H_EFF_M = 50.0          # m: median height of the lowest-100-m layer above the base


def retrieval_consistency(tab):
    """Liquid water content that the retrieval places in the lowest 100 m, LWC = 4/3 pi rho_w N r_v^3 with the
    volume-mean radius r_v = r_e exp(-sigma^2) of the lognormal, against the adiabatic value at the median height
    of the layer, LWC_ad = Gamma_ad(T_base, p_base) x min(50 m, depth / 2): f_ad = LWC / LWC_ad per scene.
    N r^3 does not depend on the radar calibration (r ~ k^(1/3), N ~ 1/k): f_ad is set by the LWP and by how the
    Frisch scaling distributes it over the lidar-detected liquid layer."""
    import pandas as pd
    rows = []
    for _, r in tab[tab["eprofile_box"]].iterrows():
        s = r["site"]
        if not (TBASE_DIR / f"{s}.npz").exists():
            continue
        z = np.load(SCENE_DIR / f"{s}.npz")
        tb = np.load(TBASE_DIR / f"{s}.npz")
        rv = z["re_s_100"].astype(float) * np.exp(-dsz.SIGMA_LN ** 2)
        lwc = 4.0 / 3.0 * np.pi * RHO_W * z["N_100"].astype(float) * rv ** 3 * 1e3            # g m-3
        zb = tb["z_base_msl"].astype(float)
        p = 101325.0 * (1 - 2.25577e-5 * zb) ** 5.25588
        h = np.minimum(H_EFF_M, z["depth"].astype(float) / 2.0)
        lwc_ad = gamma_ad(tb["T_base"].astype(float), p) * h * 1e-3                             # g m-3
        ok = np.isfinite(z["alpha_100"]) & np.isfinite(lwc) & np.isfinite(lwc_ad) & (lwc_ad > 0)
        f = lwc[ok] / lwc_ad[ok]
        rows.append(dict(site=s, lwc100_med=_q(lwc[ok], 50), lwc_ad_med=_q(lwc_ad[ok], 50), f_ad_med=_q(f, 50),
                         f_ad_p25=_q(f, 25), f_ad_p75=_q(f, 75), frac_superadiabatic=float(np.mean(f > 1.0))))
    out = pd.DataFrame(rows)
    # how far below the Cloudnet liquid base the condensation level would have to lie for the median scene to be
    # adiabatic (the layer's median height grows from 50 m to 50 m x f_ad): 50 m x (f_ad - 1), where f_ad > 1
    out["base_offset_for_adiabatic_m"] = np.where(out["f_ad_med"] > 1, H_EFF_M * (out["f_ad_med"] - 1.0), 0.0)
    out.to_csv(OUT / "retrieval_consistency_sites.csv", index=False, float_format="%.6g")
    return out


def ebas_sites():
    """In-situ N50 / N100 climatology from cloudnet_sites_ebas_n100.py (flagged stations removed)."""
    import pandas as pd
    fp = OUT / "ebas_n100_sites.csv"
    if not fp.exists():
        return None
    e = pd.read_csv(fp)
    bad = e.get("ebas_qc", pd.Series("ok", index=e.index)).astype(str) != "ok"
    e.loc[bad, ["N100_med", "N50_med", "N100_mean"]] = np.nan
    return e


def season_ratios():
    """Winter / summer (DJF / JJA) day-weighted site medians of alpha and a_G (>= 10 days in each season)."""
    import pandas as pd
    se = pd.read_csv(OUT / "season_table_extinction.csv")
    sa = pd.read_csv(OUT / "season_table.csv")
    rows = []
    for s, g in se.groupby("site"):
        ga = sa[sa["site"] == s].set_index("season")
        g = g.set_index("season")
        if not {"DJF", "JJA"} <= set(g.index):
            continue
        ok = (g.loc["DJF", "n_days"] >= dsz.MIN_DAYS_SEASON) and (g.loc["JJA", "n_days"] >= dsz.MIN_DAYS_SEASON)
        rows.append(dict(site=s, n_days_DJF=int(g.loc["DJF", "n_days"]), n_days_JJA=int(g.loc["JJA", "n_days"]),
                         alpha_DJF=g.loc["DJF", "alpha_daymed"], alpha_JJA=g.loc["JJA", "alpha_daymed"],
                         aG_DJF=ga.loc["DJF", "aG_daymed"], aG_JJA=ga.loc["JJA", "aG_daymed"],
                         N_DJF=g.loc["DJF", "N_med"], N_JJA=g.loc["JJA", "N_med"], both_seasons=bool(ok)))
    out = pd.DataFrame(rows)
    out["alpha_ratio"] = out["alpha_DJF"] / out["alpha_JJA"]
    out["aG_ratio"] = out["aG_DJF"] / out["aG_JJA"]
    out["N_ratio"] = out["N_DJF"] / out["N_JJA"]
    out.to_csv(OUT / "season_ratio_extinction.csv", index=False, float_format="%.6g")
    return out


def radar_changes():
    """Sites whose radar changed during the period: a_G, N_scaled and alpha per radar (same scenes, the radar of
    each day from cloudnet_sites_droplet_size's per-day table). A reflectivity offset k multiplies a_G by k^(1/3),
    N_scaled by 1/k and alpha by k^(-1/3) (Frisch scaled retrieval): the 'expected' columns give, from the a_G
    ratio to the site's main radar, the N and alpha ratios such an offset would produce."""
    import pandas as pd
    rt = pd.read_csv(OUT / "site_radar_table.csv")
    rows = []
    for s, g in rt.groupby("site"):
        if len(g) < 2:
            continue
        days = pd.read_csv(OUT / "scenes_all" / f"{s}_days.csv", dtype={"date": str}, keep_default_na=False)
        d2r = dict(zip(days["date"].astype(int), days["radar"]))
        z = np.load(SCENE_DIR / f"{s}.npz")
        rad = np.array([d2r.get(int(d), "") for d in z["date"]])
        a = z["alpha_100"].astype(float) * 1e3
        ag = z["re_s_100"].astype(float) * 1e6 * dsz.AG_FACTOR
        n = z["N_100"].astype(float) * 1e-6
        month = (z["date"] // 100) % 100
        for r in g["radar"]:
            m = (rad == r) & np.isfinite(a)
            if m.sum() < 1000:
                continue
            rows.append(dict(site=s, radar=r, band=g.loc[g["radar"] == r, "band"].iloc[0], n_scenes=int(m.sum()),
                             n_days=int(np.unique(z["date"][m]).size),
                             months=",".join(str(x) for x in sorted(set(month[m].tolist()))),
                             aG_med=_q(ag[m], 50), N_med=_q(n[m], 50), alpha_med=_q(a[m], 50)))
    out = pd.DataFrame(rows)
    if len(out):
        main_r = out.loc[out.groupby("site")["n_scenes"].idxmax(), ["site", "aG_med", "N_med", "alpha_med"]]
        out = out.merge(main_r, on="site", suffixes=("", "_main"))
        out["aG_ratio"] = out["aG_med"] / out["aG_med_main"]
        out["N_ratio"] = out["N_med"] / out["N_med_main"]
        out["alpha_ratio"] = out["alpha_med"] / out["alpha_med_main"]
        out["N_ratio_if_offset"] = out["aG_ratio"] ** -3.0
        out["alpha_ratio_if_offset"] = out["aG_ratio"] ** -1.0
        out["offset_dB_from_aG"] = 30.0 * np.log10(out["aG_ratio"])
    out.to_csv(OUT / "radar_change_extinction.csv", index=False, float_format="%.6g")
    return out


NEIGHBOURS = [("lindenberg", "falkenberg"), ("leipzig", "leipzig-lim"), ("kenttarova", "pallasjarvi"),
              ("kenttarova", "sammaltunturinlaatas"), ("pallasjarvi", "sammaltunturinlaatas")]
MIN_SCENES_DAY = 20


def neighbour_pairs():
    """Cloudnet sites 4-7 km apart observe the same cloud field: on the days both have scenes, any systematic
    difference of their a_G, N_scaled or alpha comes from the instruments and the retrieval, not from the aerosol
    or the cloud climate. Medians over the common days of the daily medians (days with >= 20 scenes at both)."""
    import pandas as pd
    t = pd.read_csv(OUT / "site_table_extinction.csv").set_index("site")
    rows = []
    for a, b in NEIGHBOURS:
        daily = {}
        for s in (a, b):
            z = np.load(SCENE_DIR / f"{s}.npz")
            ok = np.isfinite(z["alpha_100"])
            df = pd.DataFrame(dict(date=z["date"][ok], alpha=z["alpha_100"][ok] * 1e3,
                                   aG=z["re_s_100"][ok] * 1e6 * dsz.AG_FACTOR, N=z["N_100"][ok] * 1e-6))
            g = df.groupby("date")
            med = g.median()
            daily[s] = med[g.size() >= MIN_SCENES_DAY]
        common = daily[a].index.intersection(daily[b].index)
        row = dict(site_a=a, site_b=b, distance_km=float(haversine_km(t.loc[a, "lat"], t.loc[a, "lon"],
                                                                        t.loc[b, "lat"], t.loc[b, "lon"])),
                   radar_a=t.loc[a, "radar"], radar_b=t.loc[b, "radar"], lwp_a=t.loc[a, "lwp_source"],
                   lwp_b=t.loc[b, "lwp_source"], n_common_days=int(common.size),
                   first_common=str(int(common.min())) if common.size else "",
                   last_common=str(int(common.max())) if common.size else "")
        if common.size:
            for v in ("aG", "N", "alpha"):
                va, vb = daily[a].loc[common, v], daily[b].loc[common, v]
                row[f"{v}_a"] = float(np.median(va))
                row[f"{v}_b"] = float(np.median(vb))
                row[f"{v}_ratio_b_over_a"] = float(np.median(vb / va))       # median of the same-day ratios
                row[f"{v}_frac_days_b_larger"] = float(np.mean(vb > va))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "neighbour_pairs_extinction.csv", index=False, float_format="%.6g")
    return out


def assemble(tab, aer, pm, tb, rc=None, eb=None):
    import pandas as pd
    st = pd.read_csv(OUT / "site_table.csv")[["site", "aG50_med", "aG150_med", "re_der_med"]]
    t = tab.merge(st, on="site", how="left").merge(aer, on="site", how="left")
    for extra in (pm, tb, rc, eb):
        if extra is not None and len(extra):
            t = t.merge(extra, on="site", how="left")
    if "N100_med" in t:
        t["N_over_N100"] = t["N_med"] / t["N100_med"]
        t["N_over_N50"] = t["N_med"] / t["N50_med"]
    for key in ("PM25", "PM10"):
        if f"{key}_alt_m" in t:
            far = (t[f"{key}_alt_m"] - t["alt_m"]).abs() > MAX_DELEV_M
            t[f"{key}_alt_flag"] = np.where(far, "station > 500 m from site altitude: not used", "")
            t[f"{key}_mean_used"] = t[f"{key}_mean"].where(~far)
    t["abs_lat"] = t["lat"].abs()
    t["aG50_over_aG150"] = t["aG50_med"] / t["aG150_med"]
    t["band_W"] = t["band"].astype(str).str.startswith("W").astype(float)
    t["lwp_radiometer"] = (~t["lwp_source"].astype(str).str.contains("radar 89")).astype(float)
    mix = [season_mix(s) for s in t["site"]]
    t["frac_scenes_DJF"] = [m[0] for m in mix]
    t["frac_scenes_JJA"] = [m[1] for m in mix]
    t.to_csv(OUT / "site_explain_table.csv", index=False, float_format="%.6g")
    return t


CANDIDATES = [  # (column, label, group)
    ("N_med", "N_scaled (cm-3)", "retrieval / cloud"),
    ("aod500_med", "AERONET AOD 500 nm", "aerosol"),
    ("fine_aod500_med", "AERONET fine-mode AOD 500 nm (SDA)", "aerosol"),
    ("ae_med", "AERONET Angstrom exponent 440-870", "aerosol"),
    ("ai_med", "AERONET aerosol index AOD x AE", "aerosol"),
    ("PM25_mean_used", "surface PM2.5 (ug m-3)", "aerosol"),
    ("PM10_mean_used", "surface PM10 (ug m-3)", "aerosol"),
    ("N100_med", "in-situ N100, EBAS climatology (cm-3)", "aerosol"),
    ("N50_med", "in-situ N50, EBAS climatology (cm-3)", "aerosol"),
    ("tbase_med_C", "cloud-base temperature (C, CAMS)", "thermodynamics"),
    ("gamma_ad_med", "adiabatic LWC lapse rate (g m-3 km-1)", "thermodynamics"),
    ("alt_m", "site altitude (m)", "geography"),
    ("abs_lat", "absolute latitude", "geography"),
    ("depth_med", "liquid-layer depth (m)", "retrieval / cloud"),
    ("base_agl_med", "liquid-base height above site (m)", "retrieval / cloud"),
    ("frac_st3_100", "share of scenes with status 3 (rain in profile)", "retrieval / cloud"),
    ("frac_st2_100", "share of scenes with status 2 (droplets + falling bit)", "retrieval / cloud"),
    ("frac_st1_100", "share of scenes with reliable pixels", "retrieval / cloud"),
    ("aG50_over_aG150", "a_G(50 m) / a_G(150 m)", "retrieval / cloud"),
    ("lwc100_med", "retrieved LWC in the lowest 100 m (g m-3)", "retrieval / cloud"),
    ("f_ad_med", "retrieved / adiabatic LWC in the lowest 100 m", "retrieval / cloud"),
    ("band_W", "W-band radar (1) vs Ka (0)", "instrument"),
    ("lwp_radiometer", "LWP from radiometer (1) vs radar 89-GHz channel (0)", "instrument"),
    ("n_days", "days of data", "sampling"),
    ("frac_scenes_DJF", "share of scenes in DJF", "sampling"),
]
TARGETS = [("aG_med", "a_G"), ("alpha_med", "alpha"), ("N_med", "N_scaled")]


def correlation_sets(t):
    box = t[t["eprofile_box"] & (t["n_days"] >= dsz.MIN_DAYS_SITE)]
    return {"E-PROFILE box (30 sites)": box,
            "E-PROFILE box without Schneefernerhaus": box[box["site"] != "schneefernerhaus"],
            "E-PROFILE box without L'Aquila and Cabauw": box[~box["site"].isin(["laquila", "cabauw"])],
            "all sites >= 20 days except ATTO-Campina": t[(t["n_days"] >= dsz.MIN_DAYS_SITE) & (t["site"] != "campina")]}


def correlations(t):
    import pandas as pd
    from scipy.stats import spearmanr
    rows = []
    for set_name, g in correlation_sets(t).items():
        for tcol, tlab in TARGETS:
            for col, lab, grp in CANDIDATES:
                if col == tcol or col not in g:
                    continue
                x, y = g[col].to_numpy(float), g[tcol].to_numpy(float)
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 6 or np.nanstd(x[ok]) == 0:
                    continue
                rho, p = spearmanr(x[ok], y[ok])
                rows.append(dict(set=set_name, target=tlab, candidate=col, label=lab, group=grp, n=int(ok.sum()),
                                 spearman=float(rho), p=float(p)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "correlations_explain.csv", index=False, float_format="%.6g")
    return out


def _f(x, fmt="{:.2f}"):
    try:
        return "" if not np.isfinite(float(x)) else fmt.format(float(x))
    except (TypeError, ValueError):
        return str(x)


def write_tables(t, cor):
    """The report's comparison and correlation tables, generated (never typed by hand)."""
    box = t[t["eprofile_box"] & (t["n_days"] >= dsz.MIN_DAYS_SITE)].sort_values("aG_med").reset_index(drop=True)
    L = ["| site | alt (m) | a_G (um) | alpha (/km) | N (cm-3) | T_base (C) | depth (m) | status 2 / 3 | LWC/LWC_ad "
         "| AERONET site (km) | AOD 500 | AE | fine AOD | PM2.5 (station, km) | PM10 | N100 in situ (station, years) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in box.iterrows():
        pm_st = f"{_f(r.get('PM25_mean'), '{:.1f}')} ({r.get('PM25_station', '')}, {_f(r.get('PM25_dist_km'), '{:.0f}')})" \
            if np.isfinite(float(r.get("PM25_mean", np.nan))) else ""
        if pm_st and isinstance(r.get("PM25_alt_flag"), str) and r.get("PM25_alt_flag"):
            pm_st += f", station {_f(r.get('PM25_alt_m') - r['alt_m'], '{:+.0f}')} m: not used"
        n100 = f"{_f(r.get('N100_med'), '{:.0f}')} ({r.get('ebas_station', '')}, {r.get('ebas_period', '')})" \
            if np.isfinite(float(r.get("N100_med", np.nan))) else ""
        aer_site = f"{r.get('aeronet_site', '')} ({_f(r.get('aeronet_dist_km'), '{:.0f}')})" \
            if isinstance(r.get("aeronet_site"), str) and r.get("aeronet_site") else ""
        L.append(f"| {r['name']} | {r['alt_m']:.0f} | {_f(r['aG_med'])} | {_f(r['alpha_med'], '{:.1f}')} | "
                 f"{_f(r['N_med'], '{:.0f}')} | {_f(r.get('tbase_med_C'), '{:.1f}')} | {_f(r['depth_med'], '{:.0f}')} | "
                 f"{_f(100 * r['frac_st2_100'], '{:.0f}')} / {_f(100 * r['frac_st3_100'], '{:.0f}')} % | "
                 f"{_f(r.get('f_ad_med'))} | {aer_site} | {_f(r.get('aod500_med'), '{:.3f}')} | "
                 f"{_f(r.get('ae_med'))} | {_f(r.get('fine_aod500_med'), '{:.3f}')} | {pm_st} | "
                 f"{_f(r.get('PM10_mean'), '{:.1f}')} | {n100} |")
    L += ["", "| candidate | a_G: rho (p), n | alpha: rho (p), n | N_scaled: rho (p), n |", "|---|---|---|---|"]
    main = cor[cor["set"] == "E-PROFILE box (30 sites)"]
    for col, lab, _grp in CANDIDATES:
        cells = []
        for _tc, tl in TARGETS:
            m = main[(main["candidate"] == col) & (main["target"] == tl)]
            cells.append("" if m.empty else f"{m['spearman'].iloc[0]:+.2f} ({m['p'].iloc[0]:.2g}), {m['n'].iloc[0]}")
        L.append(f"| {lab} | " + " | ".join(cells) + " |")
    for set_name in [s for s in cor["set"].unique() if s != "E-PROFILE box (30 sites)"]:
        L += ["", f"Sensitivity set: {set_name}", "", "| candidate | a_G | alpha | N_scaled |", "|---|---|---|---|"]
        g = cor[cor["set"] == set_name]
        for col, lab, _grp in CANDIDATES:
            cells = []
            for _tc, tl in TARGETS:
                m = g[(g["candidate"] == col) & (g["target"] == tl)]
                cells.append("" if m.empty else f"{m['spearman'].iloc[0]:+.2f} ({m['p'].iloc[0]:.2g}), {m['n'].iloc[0]}")
            if any(cells):
                L.append(f"| {lab} | " + " | ".join(cells) + " |")
    (OUT / "explain_tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def plot_vs(t, cor, panels, path):
    """Rows: a_G and alpha; columns: the candidate variables of `panels` [(column, axis label, log x?, short name)].
    Column 'N_med' carries the constant-water-content lines (a_G ~ N^-1/3, alpha ~ N^+1/3) through the median
    site; column 'f_ad_med' marks the adiabatic limit (f_ad = 1: larger values are physically impossible)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, LogLocator
    box = t[t["eprofile_box"] & (t["n_days"] >= dsz.MIN_DAYS_SITE)]
    main = cor[cor["set"] == "E-PROFILE box (30 sites)"]
    fig, axes = plt.subplots(2, len(panels), figsize=(6.6 * len(panels), 10.2), dpi=150)
    band_col = {True: "#1f77b4", False: "#ff7f0e"}
    plain = FuncFormatter(lambda v, _pos: f"{v:g}")
    letters = iter("abcdefghij")
    rows_def = [("aG_med", "a_G", "site median a$_G$ [µm]", False, "a$_G$"),
                ("alpha_med", "alpha", "site median $\\alpha$ [km$^{-1}$]", True, "$\\alpha$")]
    for i, (tcol, tlab, ylab, logy, tshort) in enumerate(rows_def):
        for j, (col, xlab, logx, xshort) in enumerate(panels):
            ax = axes[i, j]
            g = box[np.isfinite(box[col].to_numpy(float))]
            for _, r in g.iterrows():
                w = str(r["band"]).startswith("W")
                mk = "o" if "radar 89" not in str(r["lwp_source"]) else "s"
                ax.plot(r[col], r[tcol], mk, ms=7.5, color=band_col[w], mec="k", mew=0.6, zorder=3)
                ax.annotate(r["name"][:13], (r[col], r[tcol]), xytext=(4, 3), textcoords="offset points",
                            fontsize=6.8, zorder=4)
            m = main[(main["candidate"] == col) & (main["target"] == tlab)]
            stat = "" if m.empty else (f"Spearman {m['spearman'].iloc[0]:+.2f} (p = {m['p'].iloc[0]:.1g}, "
                                       f"n = {m['n'].iloc[0]})")
            ax.set_title(f"({next(letters)}) {tshort} vs {xshort}: {stat}", fontsize=10.5)
            if col == "N_med":           # same liquid water content: a_G ~ N^(-1/3), alpha ~ N^(+1/3)
                nn = np.logspace(np.log10(0.9 * g[col].min()), np.log10(1.1 * g[col].max()), 50)
                n0, y0 = float(np.median(g[col])), float(np.median(g[tcol]))
                expo = -1.0 / 3.0 if tcol == "aG_med" else 1.0 / 3.0
                ax.plot(nn, y0 * (nn / n0) ** expo, color="0.35", ls="-.", lw=1.1, zorder=2)
            if col == "f_ad_med":
                ax.axvspan(1.0, 3.0, color="tab:red", alpha=0.07, zorder=1)
                ax.axvline(1.0, color="tab:red", lw=1.0, zorder=2)
            if logx:
                ax.set_xscale("log")
                ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
                ax.xaxis.set_major_formatter(plain)
                ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _pos: ""))
            if logy:
                ax.set_yscale("log")
                ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
                ax.yaxis.set_major_formatter(plain)
                ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _pos: ""))
            ax.axhline(20.0 if tcol == "alpha_med" else 3.3, color="k", ls="--", lw=1.0, zorder=2)
            ax.set_xlabel(xlab, fontsize=9.5)
            ax.set_ylabel(ylab)
            ax.grid(alpha=0.3, which="both")
    handles = [Line2D([], [], marker="o", ls="", color="#1f77b4", mec="k", label="W-band (94 GHz) radar"),
               Line2D([], [], marker="o", ls="", color="#ff7f0e", mec="k", label="Ka-band (35 GHz) radar"),
               Line2D([], [], marker="s", ls="", color="0.6", mec="k", label="LWP from the 89-GHz channel of the radar"),
               Line2D([], [], ls="--", color="k", label="table values a$_G$ = 3.3 µm, $\\alpha$ = 20 km$^{-1}$"),
               Line2D([], [], ls="-.", color="0.35", label="same water content: a$_G$ $\\propto$ N$^{-1/3}$, "
                                                            "$\\alpha$ $\\propto$ N$^{+1/3}$"),
               Patch(color="tab:red", alpha=0.15, label="retrieved LWC above the adiabatic value 50 m above the base")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=9, frameon=False)
    fig.suptitle("Site medians of the Cloudnet droplet radius a$_G$ and extinction $\\alpha$ (lowest 100 m above the "
                 "liquid base, 2025-01 .. 2026-10) against candidate explanations: the 30 sites of the E-PROFILE box",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.035, 1, 0.965))
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--offline", action="store_true", help="use only the cached external files")
    ap.add_argument("--plot-only", action="store_true", help="tables + figure from site_explain_table.csv")
    a = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    import pandas as pd
    t0 = time_mod.time()
    if a.plot_only:
        t = pd.read_csv(OUT / "site_explain_table.csv")
    else:
        tab = pd.read_csv(OUT / "site_table_extinction.csv")
        aer = aeronet_for_sites(tab, a.offline)
        pm = eea_for_sites(tab, a.offline)
        tb = cams_tbase(tab, a.workers)
        rc = retrieval_consistency(tab)
        eb = ebas_sites()
        season_ratios()
        radar_changes()
        neighbour_pairs()
        t = assemble(tab, aer, pm, tb, rc, eb)
    cor = correlations(t)
    write_tables(t, cor)
    panels = [("N_med", "site median N$_{scaled}$ [cm$^{-3}$] (Cloudnet retrieval, same scenes)", True,
               "N$_{scaled}$"),
              ("fine_aod500_med", "AERONET fine-mode AOD at 500 nm (SDA level 1.5, median of the days)", True,
               "fine-mode AOD"),
              ("f_ad_med", "retrieved / adiabatic LWC in the lowest 100 m (median of the scenes)", True,
               "LWC / LWC$_{ad}$")]
    plot_vs(t, cor, panels, FIG_DIR / "fig_cloudnet_sites_vs_aerosol.png")
    with pd.option_context("display.width", 250, "display.max_rows", 400, "display.max_columns", 20):
        print(cor[cor["set"] == "E-PROFILE box (30 sites)"].round(3).to_string(index=False))
    print(f"total runtime {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
