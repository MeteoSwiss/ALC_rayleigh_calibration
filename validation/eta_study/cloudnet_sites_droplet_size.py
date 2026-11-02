"""Eta study, extension of step 1d: Cloudnet droplet size at the calibration cloud base, all `der` sites.

Step 1d (cloudnet_microphysics_match.py) fixed a_G = 3.3 um from three Cloudnet sites. This script asks
whether that value holds at every Cloudnet site that produces the droplet effective radius (`der`)
product, 2025-01-01..2026-10-31 (files from cloudnet_der_download.py). Per site it computes:

  A. all liquid-cloud scenes: every 30-s Cloudnet profile whose LOWEST liquid layer -- the layer a
     ceilometer reports as its first cloud base -- has its base inside the calibration cloud-base
     window 0.25-2.4 km above the site (the window applies to the Cloudnet liquid base, which sits
     60-95 m below the ceilometer CBH at the step-1d sites);
  B. where an E-PROFILE Vaisala ceilometer stands at (<= 1 km) or near (<= 11 km) the site: the same
     quantities at that unit's accepted cloud-calibration profiles. The nearest Cloudnet profile within
     15 min is taken and the liquid base is searched next to the ceilometer CBH, exactly as in step 1d.

Conventions shared with step 1d:
  * liquid pixel = der_retrieval_status 1..3 (status 4 "surrounding ice" never carries a radius);
  * liquid base = lowest gate of a contiguous liquid run (B: the lowest liquid gate within +-150 m of
    the ceilometer CBH, followed down through contiguous liquid); base edge = gate centre minus half
    the local gate spacing (step 1d used the median spacing: identical on regular grids, but several
    RPG-FMCW-94 grids here change spacing with the chirp, e.g. 9.6 / 34 m at Hyytiala);
  * a layer of depth D = the run's gates whose centre lies in [base_edge, base_edge + D), D = 50 /
    100 / 150 m; medians over the status == 1 (reliable) pixels of the layer only;
  * r_e = der_scaled (Frisch method 2, scaled to the LWP), a_G = r_e exp(-1.5 sigma^2) = 0.832 r_e
    (lognormal sigma = 0.35, the CloudnetPy der default with N = 200 cm-3); der (method 1, fixed N)
    and N_scaled are carried along. Values outside 0 < r < 1 mm and N <= 0 are discarded.

Scene sources for B (time, CBH of each accepted calibration profile):
  msfov : C:/DATA/Projects/202606_E-PROFILE_calibration/old/ms_fov/results.jsonl (Payerne B/C,
          Palaiseau A, Lindenberg C >= 2025-05; the list behind the step-1d 3.3 um);
  dump  : C:/DATA/Projects/202611_eta_study/cloudnet_sites/calib_dump/<key>_profiles.npz
          (cloudnet_sites_calib_dump.py: Chilbolton, L'Aquila and Norunda CL51, Middle Wallop and
          Akrotiri CL31).

--files latest (default) reads the latest portal version of every file (per-site folders of
D:/Cloudnet_Droplet_radius, falling back to the flat July collection where the portal file is
identical); --files july reads only the flat July collection (step-1d inputs) to size the effect of
the portal reprocessing on the original sites.

Run:  python validation/eta_study/cloudnet_sites_droplet_size.py [--workers 8] [--files latest|july]
                                                                [--sites a,b] [--plot-only]
Outputs (C:/DATA/Projects/202611_eta_study/cloudnet_sites/[july/]): scenes_all/<site>.npz,
calib_scenes.csv, site_table.csv/.md (all tables), site_radar_table.csv, season_table.csv,
calib_table.csv, across_sites.csv, version_effect.md (latest run, once the july run exists),
fig_aG_sites.png (one panel per site), fig_aG_summary.png (sites, latitude, a_G vs N_scaled).
Read-only on every input.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import json
import sys
import time as time_mod
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CLOUDNET = Path("D:/Cloudnet_Droplet_radius")
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
MSFOV_JSONL = DATA / "old" / "ms_fov" / "results.jsonl"
OUT = Path("C:/DATA/Projects/202611_eta_study/cloudnet_sites")
DUMP_DIR = OUT / "calib_dump"
AVAIL_CSV = OUT / "der_availability.csv"
CENSUS = REPO / "validation" / "scope_l1_2026_census.json"
REPORT_FIG_DIR = REPO / "doc" / "reports" / "figs_eta_cloudnet_sites"

SIGMA_LN = 0.35
AG_FACTOR = float(np.exp(-1.5 * SIGMA_LN ** 2))       # 0.832
LAYERS_M = (50, 100, 150)
CBH_WINDOW_M = (250.0, 2400.0)                         # calibration cloud-base window, AGL
BASE_SEARCH_M = 150.0                                  # B: liquid gate searched within +-150 m of CBH
MAX_DT_H = 0.25                                        # B: nearest Cloudnet profile within 15 min
AG_STEP1D = 3.3                                        # um, the value under test (never tuned here)
# E-PROFILE domain = the operational CAMS Europe+Arctic box (ALC_CAMS_AREA N80/W-30/S27/E45)
EPROFILE_BOX = dict(north=80.0, west=-30.0, south=27.0, east=45.0)
MIN_DAYS_SITE = 20                                     # a site enters the across-site statistics
MIN_DAYS_SEASON = 10                                   # a season counts in the season-balanced value

# unit key -> (label, Cloudnet site, type, distance km, scene source, earliest day)
UNITS = {
    "0-20000-0-06610_C": ("Payerne CL61", "payerne", "CL61", 0.1, "msfov", None),
    "0-20000-0-06610_B": ("Payerne CL31 (distorted unit)", "payerne", "CL31", 0.1, "msfov", None),
    "0-250-1001-07151_A": ("Palaiseau CL31", "palaiseau", "CL31", 0.2, "msfov", None),
    "0-20000-0-10393_C": ("Lindenberg CL61 (>=2025-05)", "lindenberg", "CL61", 0.3, "msfov", "20250501"),
    "0-826-300-3_A": ("Chilbolton CL51", "chilbolton", "CL51", 0.1, "dump", None),
    "0-20008-0-CEO_A": ("L'Aquila CL51", "laquila", "CL51", 0.1, "dump", None),
    "0-20008-0-NOR_A": ("Norunda CL51", "norunda", "CL51", 0.1, "dump", None),
    "0-20000-0-03749_A": ("Middle Wallop CL31 (9 km from Chilbolton)", "chilbolton", "CL31", 9.2, "dump", None),
    "0-20000-0-17601_A": ("Akrotiri CL31 (11 km from Limassol)", "limassol", "CL31", 10.9, "dump", None),
}
UNIT_COLOURS = {"0-20000-0-06610_C": "C3", "0-20000-0-06610_B": "C7", "0-250-1001-07151_A": "C0",
                "0-20000-0-10393_C": "C1", "0-826-300-3_A": "C2", "0-20008-0-CEO_A": "C4",
                "0-20008-0-NOR_A": "C5", "0-20000-0-03749_A": "C6", "0-20000-0-17601_A": "C8"}


# =============================================================================
#  Inputs: files, sites, instruments
# =============================================================================
def der_files(site: str, mode: str) -> dict[str, Path]:
    """date (YYYYMMDD) -> der file. latest: per-site folder first (portal version), then the flat
    July collection (identical to the portal wherever no per-site copy was needed). july: flat only."""
    files = {p.name[:8]: p for p in CLOUDNET.glob(f"*_{site}_der.nc")}
    if mode == "latest":
        files.update({p.name[:8]: p for p in (CLOUDNET / site).glob(f"*_{site}_der.nc")})
    return files


def load_sites() -> dict[str, dict]:
    with open(AVAIL_CSV, newline="", encoding="utf-8") as f:
        rows = {r["site"]: r for r in csv.DictReader(f)}
    for r in rows.values():
        r["lat"], r["lon"] = float(r["lat"]), float(r["lon"])
        b = EPROFILE_BOX
        r["eprofile_box"] = (b["south"] <= r["lat"] <= b["north"]) and (b["west"] <= r["lon"] <= b["east"])
    return rows


RADAR_KEYS = ("FMCW", "MIRA", "BASTA", "KAZR", "Copernicus", "Galileo", "WACR", "radar")
LIDAR_KEYS = ("CHM15k", "CL31", "CL51", "CL61", "CT25k", "DA10", "PollyXT", "CS135", "Halo", "minimpl",
              "MPL")
MWR_KEYS = ("HATPRO", "MWR-3C", "WVR-1100", "Radiometrics", "MP3000", "HAMSTRAD", "LHATPRO")


def instruments(source: str) -> tuple[str, str, str]:
    """(radar, lidar, LWP source) from the `source` global attribute (one instrument per line)."""
    lines = [s.strip() for s in source.split("\n") if s.strip()]

    def pick(keys):
        hits = [s for s in lines if any(k.lower() in s.lower() for k in keys)]
        return hits[0] if hits else ""
    radar, lidar, mwr = pick(RADAR_KEYS), pick(LIDAR_KEYS), pick(MWR_KEYS)
    for prefix in ("RPG-Radiometer Physics ", "METEK ", "ProSensing ", "Lufft ", "Vaisala ",
                   "Radiometrics ", "TROPOS ", "Campbell Scientific "):
        radar, lidar, mwr = (x.replace(prefix, "") for x in (radar, lidar, mwr))
    if not mwr:   # no radiometer listed: the LWP comes from the 89-GHz channel of an RPG radar
        mwr = "radar 89 GHz channel" if "FMCW" in radar else "none"
    return radar, lidar, mwr


def radar_band(radar: str) -> str:
    r = radar.upper()
    if any(k in r for k in ("FMCW-94", "BASTA", "WACR", "-94")):
        return "W (94 GHz)"
    if any(k in r for k in ("MIRA", "KAZR", "COPERNICUS", "FMCW-35", "-35")):
        return "Ka (35 GHz)"
    return "?"


# =============================================================================
#  Core: liquid runs and layer medians (shared by A and B)
# =============================================================================
def gate_edges(h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lower / upper edge of every gate from the LOCAL spacing (grids change with the RPG chirp)."""
    d = np.diff(h)
    lower = np.empty_like(h)
    upper = np.empty_like(h)
    lower[1:] = h[1:] - 0.5 * d
    lower[0] = h[0] - 0.5 * d[0]
    upper[:-1] = h[:-1] + 0.5 * d
    upper[-1] = h[-1] + 0.5 * d[-1]
    return lower, upper


def lowest_liquid_run(liq: np.ndarray):
    """Per profile: (any liquid, base gate, top gate) of the LOWEST contiguous liquid run."""
    nz = liq.shape[1]
    anyl = liq.any(axis=1)
    jb = np.argmax(liq, axis=1)
    after = ~liq & (np.arange(nz)[None, :] > jb[:, None])
    jt = np.where(after.any(axis=1), np.argmax(after, axis=1) - 1, nz - 1)
    return anyl, jb, jt


def run_near(liq_row: np.ndarray, h: np.ndarray, cloud_msl: float):
    """B: lowest liquid gate within +-BASE_SEARCH_M of the ceilometer cloud base, followed DOWN to the
    base of its contiguous liquid run and UP to its top (step-1d definition). None if no liquid."""
    cand = np.flatnonzero(liq_row & (np.abs(h - cloud_msl) <= BASE_SEARCH_M))
    if cand.size == 0:
        return None
    jb = int(cand[0])
    while jb > 0 and liq_row[jb - 1]:
        jb -= 1
    jt = jb
    while jt + 1 < liq_row.size and liq_row[jt + 1]:
        jt += 1
    return jb, jt


def layer_medians(st, ders, der, nsc, h, lower, upper, rows, jb, jt, kmax):
    """Medians over the status == 1 pixels of the lowest 50/100/150 m of the liquid runs (rows, jb, jt).
    Returns a dict of (n,) arrays: n_liq_D, n_st1_D, re_s_D (der_scaled), re_d_100 (der), N_100, and
    eff_100 = the depth the "100 m" layer really spans (whole gates: the gate-centre rule makes it
    ~90-110 m depending on the grid, and a sub-metre shift of the base edge can add or drop a gate)."""
    nz = h.size
    g = jb[:, None] + np.arange(kmax)[None, :]
    inrun = g <= jt[:, None]
    gc = np.minimum(g, nz - 1)
    rel = h[gc] - lower[jb][:, None]
    r2 = rows[:, None]
    s, vs, vd, vn = st[r2, gc], ders[r2, gc], der[r2, gc], nsc[r2, gc]
    ok_s = (vs > 0) & (vs < 1e-3)
    ok_d = (vd > 0) & (vd < 1e-3)
    ok_n = vn > 0
    st1 = s == 1
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN rows -> NaN, on purpose
        for D in LAYERS_M:
            m = inrun & (rel < D)
            m1 = m & st1
            out[f"n_liq_{D}"] = m.sum(axis=1).astype(np.int16)
            out[f"n_st1_{D}"] = m1.sum(axis=1).astype(np.int16)
            out[f"re_s_{D}"] = np.nanmedian(np.where(m1 & ok_s, vs, np.nan), axis=1)
            if D == 100:
                out["re_d_100"] = np.nanmedian(np.where(m1 & ok_d, vd, np.nan), axis=1)
                out["N_100"] = np.nanmedian(np.where(m1 & ok_n, vn, np.nan), axis=1)
                last = np.minimum(jb + np.maximum(m.sum(axis=1), 1) - 1, nz - 1)   # layer gates are contiguous
                out["eff_100"] = upper[last] - lower[jb]
    return out


def _filled(var, fill=np.nan, dtype=float):
    """Read a variable with its _FillValue turned into `fill` (the masked-array trap of step 1d)."""
    return np.ma.filled(np.ma.asarray(var[:], dtype=dtype), fill)


# =============================================================================
#  Worker: one (site, day) file -> A arrays + B records
# =============================================================================
A_FIELDS = ("t_epoch", "base_agl", "depth", "re_s_50", "re_s_100", "re_s_150", "re_d_100", "N_100",
            "eff_100", "n_st1_100", "n_liq_100")


def process_day(job):
    """job = (site, ds, path, scenes) with scenes = [(row_id, unit_key, t_epoch_s, cbh_m, alt_m), ...]."""
    from netCDF4 import Dataset
    site, ds, path, scenes = job
    res = dict(site=site, ds=ds, ok=False)
    try:
        with Dataset(path) as nc:
            t_h = np.asarray(nc.variables["time"][:], float)
            h = np.asarray(nc.variables["height"][:], float)
            salt = float(np.nanmedian(_filled(nc.variables["altitude"])))
            st = _filled(nc.variables["der_retrieval_status"], fill=0, dtype=int)
            der = _filled(nc.variables["der"])
            ders = _filled(nc.variables["der_scaled"])
            nsc = _filled(nc.variables["N_scaled"])
            source = getattr(nc, "source", "")
    except Exception as exc:  # noqa: BLE001 - a corrupt file must not kill the run
        res["error"] = f"{type(exc).__name__}: {exc}"
        return res
    if h.size < 2 or t_h.size == 0:
        res["error"] = "empty"
        return res
    radar, lidar, mwr = instruments(source)
    lower, upper = gate_edges(h)
    kmax = int(np.ceil(max(LAYERS_M) / max(float(np.min(np.diff(h))), 1.0))) + 2
    k1 = int(np.argmin(np.abs(h - salt - 1000.0)))
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    liq = (st >= 1) & (st <= 3)
    res.update(ok=True, radar=radar, lidar=lidar, mwr=mwr, salt=salt, n_prof=int(t_h.size),
               dz_1km=float(upper[k1] - lower[k1]), has_der_scaled=bool(np.isfinite(ders).any()))

    # --- A: all profiles, lowest liquid layer in the calibration cloud-base window ---
    anyl, jb, jt = lowest_liquid_run(liq)
    base_agl = lower[jb] - salt
    win = anyl & (base_agl >= CBH_WINDOW_M[0]) & (base_agl <= CBH_WINDOW_M[1])
    rows = np.flatnonzero(win)
    res["n_liq_any"] = int(anyl.sum())
    res["n_win"] = int(rows.size)
    a = {}
    if rows.size:
        lm = layer_medians(st, ders, der, nsc, h, lower, upper, rows, jb[rows], jt[rows], kmax)
        res["n_st1_100_any"] = int((lm["n_st1_100"] > 0).sum())
        keep = np.isfinite(lm["re_s_150"])          # any reliable radius within 150 m -> stored
        r = rows[keep]
        a = dict(t_epoch=(day0 + np.round(t_h[r] * 3600.0)).astype(np.int64),
                 base_agl=base_agl[r].astype(np.float32),
                 depth=(upper[jt[r]] - lower[jb[r]]).astype(np.float32))
        for k in ("re_s_50", "re_s_100", "re_s_150", "re_d_100", "N_100", "eff_100"):
            a[k] = lm[k][keep].astype(np.float32)
        a["n_st1_100"] = lm["n_st1_100"][keep]
        a["n_liq_100"] = lm["n_liq_100"][keep]
    else:
        res["n_st1_100_any"] = 0
    res["A"] = a

    # --- B: accepted calibration profiles of the units at this site ---
    recs = []
    if scenes:
        hrs = np.array([(t - day0) / 3600.0 for _, _, t, _, _ in scenes])
        tis = np.clip(np.searchsorted(t_h, hrs), 1, t_h.size - 1)
        tis = np.where(np.abs(t_h[tis - 1] - hrs) <= np.abs(t_h[tis] - hrs), tis - 1, tis)
        if t_h.size == 1:
            tis = np.zeros_like(tis)
        b_rows, b_jb, b_jt, b_idx = [], [], [], []
        for (rid, key, _t, cbh, alt), hh, ti in zip(scenes, hrs, tis):
            rec = dict(rid=rid, dt_s=float((t_h[ti] - hh) * 3600.0))
            if abs(t_h[ti] - hh) > MAX_DT_H:
                rec["status"] = "notime"
                recs.append(rec)
                continue
            cloud_msl = cbh + (salt if alt is None else alt)
            rec.update(status="matched", cloud_msl=cloud_msl)
            rn = run_near(liq[ti], h, cloud_msl)
            if rn is None:
                # where is the profile's lowest liquid base instead (NaN: no liquid at all)?
                rec.update(liquid=False,
                           lowest_offset=float(lower[jb[ti]] - cloud_msl) if anyl[ti] else np.nan)
            else:
                rec.update(liquid=True, base_offset=float(lower[rn[0]] - cloud_msl),
                           liq_depth=float(upper[rn[1]] - lower[rn[0]]), base_status=int(st[ti, rn[0]]))
                b_rows.append(ti)
                b_jb.append(rn[0])
                b_jt.append(rn[1])
                b_idx.append(len(recs))
            recs.append(rec)
        if b_rows:
            lm = layer_medians(st, ders, der, nsc, h, lower, upper, np.array(b_rows), np.array(b_jb),
                               np.array(b_jt), kmax)
            for j, i in enumerate(b_idx):
                for k, v in lm.items():
                    recs[i][k] = float(v[j])
    res["B"] = recs
    return res


# =============================================================================
#  Calibration scenes (B)
# =============================================================================
def scenes_msfov(path: Path):
    """(unit key, YYYYMMDD, t_epoch_s, cbh_m) from the July ms_fov per-profile results."""
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            key = r.get("key")
            if key not in UNITS or UNITS[key][4] != "msfov" or "time" not in r:
                continue
            for t, z in zip(r["time"], r["cbh"]):
                out.append((key, r["date"], int(t), float(z)))
    return out


def scenes_dump(folder: Path):
    out = []
    for key, u in UNITS.items():
        fp = folder / f"{key}_profiles.npz"
        if u[4] != "dump" or not fp.exists():
            continue
        z = np.load(fp, allow_pickle=True)
        t_s = (z["time_ns"] // 1_000_000_000).astype(np.int64)
        for d, t, cbh in zip(z["date"], t_s, z["cbh_m"]):
            out.append((key, str(int(d)), int(t), float(cbh)))
    return out


def unit_altitudes() -> dict[str, float | None]:
    """Altitude used for the cloud base MSL: the Cloudnet site's (None) for co-located units, as in
    step 1d; the ceilometer's own (census) for the units that stand a few km away."""
    census = {f"{c['wmo']}_{c['ident']}": c for c in json.load(open(CENSUS, encoding="utf-8"))}
    return {k: (None if u[3] <= 1.0 else float(census[k]["alt"])) for k, u in UNITS.items()}


# =============================================================================
#  Driver
# =============================================================================
def run(args, out_dir: Path):
    import pandas as pd
    sites = load_sites()
    use = sorted(sites) if args.sites == "all" else [s.strip() for s in args.sites.split(",")]
    if args.files == "july":
        use = [s for s in use if any(CLOUDNET.glob(f"*_{s}_der.nc"))]
    alts = unit_altitudes()

    # calibration scenes grouped per (Cloudnet site, day)
    sc = pd.DataFrame(scenes_msfov(MSFOV_JSONL) + scenes_dump(DUMP_DIR),
                      columns=["key", "date", "t_epoch", "cbh_m"])
    earliest = sc["key"].map(lambda k: UNITS[k][5] or "00000000")
    sc = sc[sc["date"] >= earliest].reset_index(drop=True)
    sc["cn_site"] = sc["key"].map(lambda k: UNITS[k][1])
    by_day = defaultdict(list)
    for rid, k, d, t, z, s in zip(sc.index, sc["key"], sc["date"], sc["t_epoch"], sc["cbh_m"], sc["cn_site"]):
        by_day[(s, d)].append((int(rid), k, int(t), float(z), alts[k]))

    jobs, n_missing_cal = [], 0
    for s in use:
        files = der_files(s, args.files)
        for ds, p in sorted(files.items()):
            jobs.append((s, ds, str(p), by_day.pop((s, ds), [])))
    n_missing_cal = sum(len(v) for (s, _), v in by_day.items() if s in use)
    print(f"[{args.files}] {len(jobs)} site-days on {len(use)} sites, {len(sc)} calibration profiles "
          f"({n_missing_cal} on days without a der file) -> {args.workers} processes", flush=True)

    acc = defaultdict(lambda: defaultdict(list))
    day_meta = defaultdict(list)
    b_recs = {}
    errors = []
    t0 = time_mod.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for k, res in enumerate(ex.map(process_day, jobs, chunksize=4), 1):
            s = res["site"]
            if not res["ok"]:
                errors.append(f"{s} {res['ds']} {res.get('error')}")
            else:
                day_meta[s].append(dict(date=res["ds"], radar=res["radar"], lidar=res["lidar"],
                                        mwr=res["mwr"], salt=res["salt"], n_prof=res["n_prof"],
                                        dz_1km=res["dz_1km"], n_liq_any=res["n_liq_any"],
                                        n_win=res["n_win"], n_st1_100_any=res["n_st1_100_any"],
                                        has_der_scaled=res["has_der_scaled"]))
                for f, v in res["A"].items():
                    acc[s][f].append(v)
                n = len(res["A"].get("t_epoch", ()))
                acc[s]["date"].append(np.full(n, int(res["ds"]), np.int32))
                acc[s]["radar_day"].append(np.full(n, len(day_meta[s]) - 1, np.int32))
                for rec in res["B"]:
                    b_recs[rec.pop("rid")] = rec
            if k % 1000 == 0:
                print(f"  {k}/{len(jobs)} site-days ({time_mod.time() - t0:.0f} s)", flush=True)
    print(f"[{args.files}] read in {time_mod.time() - t0:.0f} s, {len(errors)} unreadable files", flush=True)
    (out_dir / "scenes_all").mkdir(parents=True, exist_ok=True)
    if errors:
        (out_dir / "unreadable_files.txt").write_text("\n".join(errors) + "\n", encoding="utf-8")
    for s in use:
        meta = pd.DataFrame(day_meta[s])
        if meta.empty:
            continue
        meta.to_csv(out_dir / "scenes_all" / f"{s}_days.csv", index=False)
        arrays = {f: (np.concatenate(acc[s][f]) if acc[s][f] else np.array([]))
                  for f in A_FIELDS + ("date", "radar_day")}
        np.savez_compressed(out_dir / "scenes_all" / f"{s}.npz", **arrays)

    # B: per calibration profile
    extra = pd.DataFrame.from_dict(b_recs, orient="index")
    b = sc.join(extra)
    b["status"] = b["status"].fillna("nofile")
    b["label"] = b["key"].map(lambda k: UNITS[k][0])
    b["itype"] = b["key"].map(lambda k: UNITS[k][2])
    b.to_csv(out_dir / "calib_scenes.csv", index=False, float_format="%.5g")
    return use


# =============================================================================
#  Summaries
# =============================================================================
def _q(x, p):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if x.size else np.nan


def _season(month: np.ndarray) -> np.ndarray:
    return np.array(["DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF"])[month - 1]


def load_site(out_dir: Path, s: str):
    import pandas as pd
    fp = out_dir / "scenes_all" / f"{s}.npz"
    if not fp.exists():
        return None, None
    z = dict(np.load(fp))
    # keep_default_na=False: an instrument missing from `source` stays "" instead of becoming NaN
    meta = pd.read_csv(out_dir / "scenes_all" / f"{s}_days.csv", dtype={"date": str}, keep_default_na=False)
    return z, meta


def daily_median(date, x):
    """Median over days of the per-day medians (day weighting: one long stratocumulus day = one vote)."""
    ok = np.isfinite(x)
    if not ok.any():
        return np.nan, 0
    d, x = date[ok], x[ok]
    order = np.argsort(d, kind="stable")
    d, x = d[order], x[order]
    cut = np.flatnonzero(np.diff(d)) + 1
    meds = [np.median(v) for v in np.split(x, cut)]
    return float(np.median(meds)), len(meds)


def summarise_sites(out_dir: Path, use):
    import pandas as pd
    sites = load_sites()
    rows, season_rows, radar_rows = [], [], []
    for s in use:
        z, meta = load_site(out_dir, s)
        if z is None:
            continue
        info = sites[s]
        radars = Counter(meta["radar"])
        main_radar = radars.most_common(1)[0][0]
        ag = z["re_s_100"] * 1e6 * AG_FACTOR
        ok = np.isfinite(ag)
        date = z["date"]
        ag_daymed, n_days_sc = daily_median(date, ag)
        row = dict(site=s, name=info["name"], country=info["country"], lat=info["lat"], lon=info["lon"],
                   alt_m=float(np.nanmedian(meta["salt"])), eprofile_box=info["eprofile_box"],
                   radar=" / ".join(f"{r} ({n} d)" for r, n in radars.most_common()),
                   band=radar_band(main_radar),
                   lidar=Counter(meta["lidar"]).most_common(1)[0][0],
                   lwp_source=Counter(meta["mwr"]).most_common(1)[0][0],
                   dz_1km_m=float(np.nanmedian(meta["dz_1km"])),
                   n_files=len(meta), first=str(meta["date"].min()), last=str(meta["date"].max()),
                   n_profiles=int(meta["n_prof"].sum()), n_liq_window=int(meta["n_win"].sum()),
                   frac_st1_100=float(meta["n_st1_100_any"].sum() / max(meta["n_win"].sum(), 1)),
                   n_scenes=int(ok.sum()), n_days=n_days_sc,
                   aG_med=_q(ag, 50), aG_p25=_q(ag, 25), aG_p75=_q(ag, 75), aG_p10=_q(ag, 10),
                   aG_p90=_q(ag, 90), aG_daymed=ag_daymed,
                   aG50_med=_q(z["re_s_50"] * 1e6 * AG_FACTOR, 50),
                   aG150_med=_q(z["re_s_150"] * 1e6 * AG_FACTOR, 50),
                   re_med=_q(z["re_s_100"][ok] * 1e6, 50), re_der_med=_q(z["re_d_100"][ok] * 1e6, 50),
                   N_med=_q(z["N_100"][ok] * 1e-6, 50), N_p25=_q(z["N_100"][ok] * 1e-6, 25),
                   N_p75=_q(z["N_100"][ok] * 1e-6, 75),
                   eff_100_med=_q(z["eff_100"][ok], 50),
                   base_agl_med=_q(z["base_agl"][ok], 50))
        # day-to-day spread of the site (IQR of the daily medians): the yardstick for site differences
        if ok.any():
            dd = pd.Series(ag[ok]).groupby(date[ok]).median()
            row["aG_day_p25"], row["aG_day_p75"] = _q(dd, 25), _q(dd, 75)
        month = (date // 100) % 100
        seas = _season(month) if date.size else np.array([])
        seasonal = []
        for sea in ("DJF", "MAM", "JJA", "SON"):
            m = seas == sea
            dm, nd = daily_median(date[m], ag[m])
            season_rows.append(dict(site=s, season=sea, aG_daymed=dm, n_days=nd, n_scenes=int(np.isfinite(ag[m]).sum()),
                                    aG_med=_q(ag[m], 50)))
            if nd >= MIN_DAYS_SEASON:
                seasonal.append(dm)
        # season-balanced value: records cover different seasons (Chilbolton = one winter), so sites are
        # also compared on the median of their four seasonal day-medians (all four seasons required)
        row["n_seasons"] = len(seasonal)
        row["aG_season_bal"] = float(np.median(seasonal)) if len(seasonal) == 4 else np.nan
        rows.append(row)
        rid = z["radar_day"]
        for rname in radars:
            idx = np.flatnonzero(meta["radar"].to_numpy() == rname)
            m = np.isin(rid, idx)
            dm, nd = daily_median(date[m], ag[m])
            radar_rows.append(dict(site=s, radar=rname, band=radar_band(rname), n_files=int(idx.size),
                                   n_days=nd, n_scenes=int(np.isfinite(ag[m]).sum()), aG_med=_q(ag[m], 50),
                                   aG_daymed=dm, N_med=_q(z["N_100"][m] * 1e-6, 50)))
    tab = pd.DataFrame(rows).sort_values("lat", ascending=False).reset_index(drop=True)
    tab.to_csv(out_dir / "site_table.csv", index=False, float_format="%.4g")
    pd.DataFrame(season_rows).to_csv(out_dir / "season_table.csv", index=False, float_format="%.4g")
    pd.DataFrame(radar_rows).to_csv(out_dir / "site_radar_table.csv", index=False, float_format="%.4g")
    return tab, pd.DataFrame(season_rows), pd.DataFrame(radar_rows)


def summarise_calib(out_dir: Path):
    """B per unit, next to the all-scene value of its site (all days, and only the unit's own days)."""
    import pandas as pd
    b = pd.read_csv(out_dir / "calib_scenes.csv", dtype={"date": str}, low_memory=False)
    rows = []
    for key, (label, site, itype, dist, src, _e) in UNITS.items():
        g = b[b["key"] == key]
        if g.empty:
            continue
        m = g[g["status"] == "matched"]
        liq = m[m["liquid"] == True]  # noqa: E712 - object column
        unm = m[m["liquid"] == False]  # noqa: E712
        off = unm["lowest_offset"].to_numpy(float) if "lowest_offset" in unm else np.full(len(unm), np.nan)
        ag = liq["re_s_100"].to_numpy(float) * 1e6 * AG_FACTOR
        ok = np.isfinite(ag)
        dmed, nd = daily_median(liq["date"].astype(int).to_numpy(), ag)
        row = dict(key=key, unit=label, type=itype, site=site, distance_km=dist, source=src,
                   n_profiles=len(g), n_days=g["date"].nunique(), n_matched=len(m),
                   frac_matched=len(m) / max(len(g), 1), dt_abs_med_s=_q(np.abs(m["dt_s"]), 50),
                   frac_liquid=len(liq) / max(len(m), 1), base_offset_med_m=_q(liq["base_offset"], 50),
                   # profiles WITHOUT liquid near the CBH: share with no liquid at all in the Cloudnet profile,
                   # and where the lowest liquid base lies otherwise (P10/P50/P90 of base - CBH, m)
                   unmatched_frac_no_liquid=float(np.mean(~np.isfinite(off))) if off.size else np.nan,
                   unmatched_lowest_offset_p10=_q(off, 10), unmatched_lowest_offset_p50=_q(off, 50),
                   unmatched_lowest_offset_p90=_q(off, 90),
                   n_st1=int(ok.sum()), n_days_st1=nd, frac_st1_of_matched=int(ok.sum()) / max(len(m), 1),
                   aG_med=_q(ag, 50), aG_p25=_q(ag, 25), aG_p75=_q(ag, 75), aG_p10=_q(ag, 10),
                   aG_p90=_q(ag, 90), aG_daymed=dmed,
                   aG50_med=_q(liq["re_s_50"].to_numpy(float) * 1e6 * AG_FACTOR, 50),
                   aG150_med=_q(liq["re_s_150"].to_numpy(float) * 1e6 * AG_FACTOR, 50),
                   re_der_med=_q(liq["re_d_100"].to_numpy(float)[ok] * 1e6, 50),
                   N_med=_q(liq["N_100"].to_numpy(float)[ok] * 1e-6, 50),
                   eff_100_med=_q(liq["eff_100"].to_numpy(float)[ok], 50))
        z, _meta = load_site(out_dir, site)
        if z is not None:
            ag_all = z["re_s_100"] * 1e6 * AG_FACTOR
            row["site_all_aG_med"] = _q(ag_all, 50)
            on_days = np.isin(z["date"], liq["date"].astype(int).unique())
            row["site_samedays_aG_med"] = _q(ag_all[on_days], 50)
            row["site_samedays_n"] = int(np.isfinite(ag_all[on_days]).sum())
            row["ratio_cal_to_all"] = row["aG_med"] / row["site_all_aG_med"]
            row["ratio_cal_to_samedays"] = row["aG_med"] / row["site_samedays_aG_med"]
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(out_dir / "calib_table.csv", index=False, float_format="%.4g")

    # pooled over units: every accepted profile counts (as the step-1d "all Vaisala" row)
    step1d = ["0-20000-0-06610_C", "0-20000-0-06610_B", "0-250-1001-07151_A", "0-20000-0-10393_C"]
    coloc = [k for k, u in UNITS.items() if u[3] <= 1.0]
    pools = [("step-1d units (Payerne B/C, Palaiseau A, Lindenberg C)", step1d),
             ("step-1d units excl. Payerne CL31 (distorted)", [k for k in step1d if k != "0-20000-0-06610_B"]),
             ("all co-located units (<= 1 km)", coloc),
             ("all co-located units excl. Payerne CL31", [k for k in coloc if k != "0-20000-0-06610_B"])]
    m = b[(b["status"] == "matched") & (b["liquid"] == True)]  # noqa: E712
    ag_all = m["re_s_100"].to_numpy(float) * 1e6 * AG_FACTOR
    prow = []
    for name, keys in pools:
        sel = m["key"].isin(keys).to_numpy() & np.isfinite(ag_all)
        unit_meds = [float(np.median(ag_all[sel & (m["key"] == k).to_numpy()])) for k in keys
                     if (sel & (m["key"] == k).to_numpy()).any()]
        prow.append(dict(pool=name, n_units=len(unit_meds), n_profiles=int(sel.sum()),
                         n_days=int(m.loc[sel, "date"].nunique()), aG_med=_q(ag_all[sel], 50),
                         aG_p25=_q(ag_all[sel], 25), aG_p75=_q(ag_all[sel], 75),
                         median_of_unit_medians=float(np.median(unit_meds))))
    pooled = pd.DataFrame(prow)
    pooled.to_csv(out_dir / "calib_pooled.csv", index=False, float_format="%.4g")
    return tab, pooled


def across_sites(tab, season, radar):
    """Spread between sites and its dependence on latitude, altitude and radar band (statistics only)."""
    from scipy.stats import spearmanr
    L = []
    for name, sel in (("E-PROFILE domain (CAMS box N27-80, W30-E45)", tab["eprofile_box"]),
                      ("all sites", tab["eprofile_box"] | ~tab["eprofile_box"])):
        g0 = tab[sel & (tab["n_days"] >= MIN_DAYS_SITE)]
        if len(g0) < 3:
            continue
        for col, lab in (("aG_med", "scene-weighted"), ("aG_daymed", "day-weighted"),
                         ("aG_season_bal", "season-balanced (4 seasons)")):
            g = g0[np.isfinite(g0[col].to_numpy(float))]
            if len(g) < 3:
                continue
            v = g[col].to_numpy(float)
            rho, p = spearmanr(g["lat"].abs(), v)
            rho_a, p_a = spearmanr(g["alt_m"], v)
            L.append(dict(group=name, weighting=lab, n_sites=len(g), median_of_sites=float(np.median(v)),
                          p25=_q(v, 25), p75=_q(v, 75), min=float(v.min()), max=float(v.max()),
                          sd=float(np.std(v, ddof=1)), min_site=g.loc[g[col].idxmin(), "site"],
                          max_site=g.loc[g[col].idxmax(), "site"],
                          n_within_2p5_5=int(((v >= 2.5) & (v <= 5.0)).sum()),
                          frac_within_10pct_of_3p3=float(np.mean(np.abs(v / AG_STEP1D - 1) <= 0.10)),
                          frac_within_20pct_of_3p3=float(np.mean(np.abs(v / AG_STEP1D - 1) <= 0.20)),
                          n_below_3p3=int((v < AG_STEP1D).sum()),
                          spearman_abslat=float(rho), p_abslat=float(p), spearman_alt=float(rho_a),
                          p_alt=float(p_a),
                          median_day_iqr=float(np.median(g["aG_day_p75"] - g["aG_day_p25"]))))
        for band in ("W (94 GHz)", "Ka (35 GHz)"):
            gb = g0[(g0["band"] == band) & np.isfinite(g0["aG_daymed"].to_numpy(float))]
            if len(gb):
                L.append(dict(group=name, weighting=f"day-weighted, {band} only", n_sites=len(gb),
                              median_of_sites=float(np.median(gb["aG_daymed"])),
                              min=float(gb["aG_daymed"].min()), max=float(gb["aG_daymed"].max())))
        # seasonal cycle inside the sites: winter (DJF) vs summer (JJA) day-weighted medians
        sp = season[season["site"].isin(g0["site"])].pivot(index="site", columns="season")
        both = (sp[("n_days", "DJF")] >= MIN_DAYS_SEASON) & (sp[("n_days", "JJA")] >= MIN_DAYS_SEASON)
        ratio = (sp[("aG_daymed", "DJF")] / sp[("aG_daymed", "JJA")])[both]
        L.append(dict(group=name, weighting="DJF / JJA day-weighted, per site", n_sites=int(both.sum()),
                      median_of_sites=float(ratio.median()), min=float(ratio.min()), max=float(ratio.max()),
                      n_djf_larger=int((ratio > 1).sum())))
    return L


# =============================================================================
#  Markdown tables
# =============================================================================
def _f(x, fmt="{:.2f}"):
    try:
        return "" if not np.isfinite(float(x)) else fmt.format(float(x))
    except (TypeError, ValueError):
        return str(x)


def write_markdown(out_dir, tab, season, radar, cal, pooled, stats):
    L = ["# Cloudnet droplet size at the calibration cloud base, all `der` sites", "",
         f"a_G = {AG_FACTOR:.3f} r_e(der_scaled), status-1 pixels of the lowest 100 m above the lowest liquid "
         f"base, base {CBH_WINDOW_M[0] / 1e3:.2f}-{CBH_WINDOW_M[1] / 1e3:.1f} km AGL. Sites north to south.", "",
         "| site | lat | alt (m) | E-PROFILE box | radar (days) | LWP | dz@1km (m) | '100 m' spans (m) | days | "
         "scenes | frac status-1 | a_G med [P25-P75] | a_G day-med | a_G season-bal. | a_G 50 / 150 m | "
         "r_e der (fixed N) | N_scaled cm-3 med [P25-P75] |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tab.iterrows():
        L.append(f"| {r['name']} | {r['lat']:.1f} | {r['alt_m']:.0f} | {'yes' if r['eprofile_box'] else 'no'} | "
                 f"{r['radar']} | {r['lwp_source']} | {r['dz_1km_m']:.0f} | {_f(r['eff_100_med'], '{:.0f}')} | "
                 f"{r['n_days']} | {r['n_scenes']} | "
                 f"{100 * r['frac_st1_100']:.0f} % | {_f(r['aG_med'])} [{_f(r['aG_p25'], '{:.1f}')}-"
                 f"{_f(r['aG_p75'], '{:.1f}')}] | {_f(r['aG_daymed'])} | {_f(r['aG_season_bal'])} | "
                 f"{_f(r['aG50_med'])} / {_f(r['aG150_med'])} | {_f(r['re_der_med'])} | {_f(r['N_med'], '{:.0f}')} "
                 f"[{_f(r['N_p25'], '{:.0f}')}-{_f(r['N_p75'], '{:.0f}')}] |")
    L += ["", "## Across sites (sites with >= %d days)" % MIN_DAYS_SITE, "",
          "| group | weighting | sites | median of site medians [P25-P75] | min (site) - max (site) | SD | "
          "within 2.5-5 um | below 3.3 | within +-10 / +-20 % of 3.3 | Spearman abs(lat) (p) | Spearman alt (p) | "
          "median within-site day IQR |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in stats:
        if "spearman_abslat" in r:
            L.append(f"| {r['group']} | {r['weighting']} | {r['n_sites']} | {r['median_of_sites']:.2f} "
                     f"[{r['p25']:.2f}-{r['p75']:.2f}] | {r['min']:.2f} ({r['min_site']}) - {r['max']:.2f} "
                     f"({r['max_site']}) | {r['sd']:.2f} | {r['n_within_2p5_5']}/{r['n_sites']} | "
                     f"{r['n_below_3p3']}/{r['n_sites']} | {100 * r['frac_within_10pct_of_3p3']:.0f} / "
                     f"{100 * r['frac_within_20pct_of_3p3']:.0f} % | {r['spearman_abslat']:+.2f} "
                     f"({r['p_abslat']:.2g}) | {r['spearman_alt']:+.2f} ({r['p_alt']:.2g}) | "
                     f"{r['median_day_iqr']:.2f} |")
        elif "n_djf_larger" in r:
            L.append(f"| {r['group']} | {r['weighting']} | {r['n_sites']} | ratio {r['median_of_sites']:.2f} | "
                     f"{r['min']:.2f} - {r['max']:.2f} | | | DJF larger at {r['n_djf_larger']}/{r['n_sites']} | | | | |")
        else:
            L.append(f"| {r['group']} | {r['weighting']} | {r['n_sites']} | {r['median_of_sites']:.2f} | "
                     f"{r['min']:.2f} - {r['max']:.2f} | | | | | | | |")
    L += ["", "## Calibration profiles of the E-PROFILE Vaisala units (B)", "",
          "| unit | km | profiles (days) | matched | liquid near CBH | '100 m' spans (m) | status-1 in 100 m (days) | "
          "a_G med [P25-P75] | a_G day-med | a_G 50 / 150 m | N cm-3 | site all scenes | site, same days | "
          "cal / all |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in cal.iterrows():
        L.append(f"| {r['unit']} | {r['distance_km']:.1f} | {r['n_profiles']} ({r['n_days']}) | "
                 f"{100 * r['frac_matched']:.0f} % | {100 * r['frac_liquid']:.0f} % | "
                 f"{_f(r['eff_100_med'], '{:.0f}')} | {r['n_st1']} "
                 f"({r['n_days_st1']}) | {_f(r['aG_med'])} [{_f(r['aG_p25'], '{:.1f}')}-"
                 f"{_f(r['aG_p75'], '{:.1f}')}] | {_f(r['aG_daymed'])} | {_f(r['aG50_med'])} / "
                 f"{_f(r['aG150_med'])} | {_f(r['N_med'], '{:.0f}')} | {_f(r.get('site_all_aG_med'))} | "
                 f"{_f(r.get('site_samedays_aG_med'))} | {_f(r.get('ratio_cal_to_all'))} |")
    L += ["", "Profiles without Cloudnet liquid within +-150 m of the CBH:", "",
          "| unit | share of matched | of which no liquid at all | lowest liquid base - CBH, P10 / P50 / P90 (m) |",
          "|---|---|---|---|"]
    for _, r in cal.iterrows():
        L.append(f"| {r['unit']} | {100 * (1 - r['frac_liquid']):.0f} % | {100 * r['unmatched_frac_no_liquid']:.0f} % | "
                 f"{_f(r['unmatched_lowest_offset_p10'], '{:+.0f}')} / {_f(r['unmatched_lowest_offset_p50'], '{:+.0f}')} / "
                 f"{_f(r['unmatched_lowest_offset_p90'], '{:+.0f}')} |")
    L += ["", "Pooled over units (every accepted profile counts):", "",
          "| pool | units | profiles (days) | a_G med [P25-P75] | median of the unit medians |",
          "|---|---|---|---|---|"]
    for _, r in pooled.iterrows():
        L.append(f"| {r['pool']} | {r['n_units']} | {r['n_profiles']} ({r['n_days']}) | {_f(r['aG_med'])} "
                 f"[{_f(r['aG_p25'], '{:.1f}')}-{_f(r['aG_p75'], '{:.1f}')}] | {_f(r['median_of_unit_medians'])} |")
    L += ["", "## Seasonal day-weighted medians of a_G (days in brackets)", "",
          "| site | DJF | MAM | JJA | SON |", "|---|---|---|---|---|"]
    for s in tab["site"]:
        g = season[season["site"] == s].set_index("season")
        cells = [f"{_f(g.loc[x, 'aG_daymed'])} ({int(g.loc[x, 'n_days'])})" if x in g.index and g.loc[x, "n_days"] else ""
                 for x in ("DJF", "MAM", "JJA", "SON")]
        L.append(f"| {tab.loc[tab['site'] == s, 'name'].iloc[0]} | " + " | ".join(cells) + " |")
    multi = radar.groupby("site").filter(lambda g: len(g) > 1)
    if len(multi):
        L += ["", "## Sites whose radar changed (day-weighted a_G per radar)", "",
              "| site | radar | band | days | a_G day-med | N cm-3 |", "|---|---|---|---|---|---|"]
        for _, r in multi.iterrows():
            L.append(f"| {r['site']} | {r['radar']} | {r['band']} | {r['n_days']} | {_f(r['aG_daymed'])} | "
                     f"{_f(r['N_med'], '{:.0f}')} |")
    (out_dir / "site_table.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# =============================================================================
#  Figures (landscape)
# =============================================================================
def plot_panels(out_dir, tab, cal_scenes, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(tab)
    ncol = 7
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(26, 2.55 * nrow + 1.2), sharex=True)
    axes = np.atleast_2d(axes)
    bins = np.arange(0.0, 10.01, 0.25)
    for ax, (_, r) in zip(axes.ravel(), tab.iterrows()):
        z, _meta = load_site(out_dir, r["site"])
        ag = z["re_s_100"] * 1e6 * AG_FACTOR
        ag = ag[np.isfinite(ag)]
        if ag.size:
            ax.hist(np.clip(ag, bins[0], bins[-1] - 1e-6), bins=bins, weights=np.full(ag.size, 1.0 / ag.size),
                    histtype="stepfilled", color="0.78", ec="0.45", lw=0.8, label="all scenes")
            ax.axvline(np.median(ag), color="k", lw=1.6)
        for key, (label, site, *_rest) in UNITS.items():
            if site != r["site"]:
                continue
            g = cal_scenes[(cal_scenes["key"] == key) & (cal_scenes["liquid"] == True)]  # noqa: E712
            x = g["re_s_100"].to_numpy(float) * 1e6 * AG_FACTOR
            x = x[np.isfinite(x)]
            if x.size:
                ax.hist(np.clip(x, bins[0], bins[-1] - 1e-6), bins=bins, weights=np.full(x.size, 1.0 / x.size),
                        histtype="step", lw=1.6, color=UNIT_COLOURS[key],
                        label=f"{label.split(' (')[0]} cal. ({np.median(x):.2f})")
        ax.axvline(AG_STEP1D, color="C3", ls=":", lw=1.4)
        thin = " - THIN" if r["n_days"] < MIN_DAYS_SITE else ""
        ax.set_title(f"{r['name']} ({r['country']}) {r['lat']:.1f}°, {r['alt_m']:.0f} m\n"
                     f"{r['band'].split(' ')[0]}-band, {r['lidar']}{thin}", fontsize=8.5)
        bal = f"\n4-season {r['aG_season_bal']:.2f}" if np.isfinite(r["aG_season_bal"]) else ""
        ax.text(0.98, 0.95, f"med {r['aG_med']:.2f} [{r['aG_p25']:.1f}-{r['aG_p75']:.1f}]\n"
                f"day-med {r['aG_daymed']:.2f}{bal}\n{r['n_scenes']} sc. / {r['n_days']} d",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
        handles = [h for h in ax.get_legend_handles_labels()[1] if "cal." in h]
        if handles:      # the small-radius corner is nearly empty at every site
            ax.legend(fontsize=6.3, loc="upper left", frameon=False)
        ax.tick_params(labelsize=7.5)
        ax.grid(alpha=0.25)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("a$_G$ = 0.832 r$_e$(der_scaled) [$\\mu$m]", fontsize=8.5)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of scenes / bin", fontsize=8.5)
    fig.suptitle("Cloudnet droplet size in the lowest 100 m above the lowest liquid base (base 0.25-2.4 km AGL, "
                 "status-1 pixels), 2025-01..2026-10, sites north to south.  Grey: all scenes (black line = median); "
                 "coloured: accepted cloud-calibration profiles of E-PROFILE Vaisala units; red dotted: step-1d "
                 "a$_G$ = 3.3 $\\mu$m.  End bins hold the clipped tails.", fontsize=10.5, y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(path, dpi=110)
    plt.close(fig)


BAND_COLOURS = {"W (94 GHz)": "C0", "Ka (35 GHz)": "C1", "?": "0.5"}


def plot_summary(tab, cal, stats, path):
    """(a) E-PROFILE-domain sites north to south: median, P25-P75, day-weighted, season-balanced and the
    calibration profiles of the Vaisala units; (b) all sites vs latitude; (c) site a_G vs N_scaled."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=(24, 9.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.28)
    ax_f, ax_l, ax_n = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) forest plot of the sites inside the E-PROFILE CAMS box, north at the top
    box = tab[tab["eprofile_box"]].sort_values("lat").reset_index(drop=True)
    for i, r in box.iterrows():
        thin = r["n_days"] < MIN_DAYS_SITE
        c = BAND_COLOURS.get(r["band"], "0.5")
        a = 0.35 if thin else 1.0
        ax_f.errorbar(r["aG_med"], i, xerr=[[r["aG_med"] - r["aG_p25"]], [r["aG_p75"] - r["aG_med"]]],
                      fmt="o", ms=7, color=c, alpha=a, elinewidth=1.2, capsize=2.5)
        ax_f.plot(r["aG_daymed"], i, "x", color=c, ms=7, mew=1.6, alpha=a)
        if np.isfinite(r["aG_season_bal"]):
            ax_f.plot(r["aG_season_bal"], i, "|", color="k", ms=12, mew=2.0, alpha=a)
        for j, (_, u) in enumerate(cal[cal["site"] == r["site"]].iterrows()):
            if np.isfinite(u["aG_med"]):
                ax_f.plot(u["aG_med"], i + 0.28 + 0.12 * j, "D", ms=6, color=UNIT_COLOURS[u["key"]], mec="k", mew=0.6)
    ax_f.set_yticks(range(len(box)))
    ax_f.set_yticklabels([f"{r['name']} ({r['lat']:.0f}°N, {r['alt_m']:.0f} m, {r['band'].split(' ')[0]})"
                          for _, r in box.iterrows()], fontsize=8)
    ax_f.set_ylim(-0.7, len(box) - 0.3)
    st = [s for s in stats if s.get("group", "").startswith("E-PROFILE") and s.get("weighting") == "scene-weighted"]
    if st:
        ax_f.axvline(st[0]["median_of_sites"], color="0.3", ls="--", lw=1.2)
    ax_f.set_xlabel("a$_G$ in the lowest 100 m [$\\mu$m]")
    ax_f.set_title("(a) sites inside the E-PROFILE CAMS box, north at the top", fontsize=10)

    # (b) all sites vs latitude; labels only where they do not pile up (outside Europe, extremes)
    med_box = float(np.median(box.loc[box["n_days"] >= MIN_DAYS_SITE, "aG_med"]))
    for _, r in tab.iterrows():
        thin = r["n_days"] < MIN_DAYS_SITE
        c = BAND_COLOURS.get(r["band"], "0.5")
        ax_l.errorbar(r["lat"], r["aG_med"], yerr=[[r["aG_med"] - r["aG_p25"]], [r["aG_p75"] - r["aG_med"]]],
                      fmt="o", ms=7, color=c, mfc=c if r["eprofile_box"] else "white", mew=1.5,
                      alpha=0.35 if thin else 1.0, elinewidth=1.0, capsize=2)
        if (not r["eprofile_box"] or abs(r["aG_med"] - med_box) > 1.2) and not thin:
            ax_l.annotate(r["name"][:14], (r["lat"], r["aG_med"]), xytext=(4, 4), textcoords="offset points",
                          fontsize=7.5)
    ax_l.set_xlabel("site latitude [°N]")
    ax_l.set_ylabel("a$_G$ in the lowest 100 m [$\\mu$m]  (median, bar = P25-P75)")
    ax_l.set_title("(b) all 41 sites vs latitude", fontsize=10)

    # (c) a_G vs N_scaled per site: the retrieval trades r_e against N through the LWP scaling
    for _, r in tab.iterrows():
        thin = r["n_days"] < MIN_DAYS_SITE
        c = BAND_COLOURS.get(r["band"], "0.5")
        ax_n.plot(r["N_med"], r["aG_med"], "o", ms=8, color=c, mfc=c if r["eprofile_box"] else "white", mew=1.5,
                  alpha=0.35 if thin else 1.0)
        if not thin and (r["N_med"] < 60 or r["N_med"] > 900 or r["aG_med"] > 5 or r["aG_med"] < 2.5):
            ax_n.annotate(r["name"][:14], (r["N_med"], r["aG_med"]), xytext=(5, 3), textcoords="offset points",
                          fontsize=7.5)
    nn = np.logspace(1, 6, 50)
    ref = tab[tab["eprofile_box"] & (tab["n_days"] >= MIN_DAYS_SITE)]
    n0, a0 = float(np.median(ref["N_med"])), float(np.median(ref["aG_med"]))
    ax_n.plot(nn, a0 * (nn / n0) ** (-1.0 / 3.0), color="0.4", lw=1, ls="-.")
    ax_n.set_xscale("log")
    ax_n.set_ylim(0.0, 7.5)
    ax_n.set_xlabel("N$_{scaled}$ in the lowest 100 m [cm$^{-3}$]  (site median)")
    ax_n.set_ylabel("a$_G$ in the lowest 100 m [$\\mu$m]  (site median)")
    ax_n.set_title("(c) site a$_G$ vs N$_{scaled}$ (dash-dot: a$_G$ $\\propto$ N$^{-1/3}$, fixed water content)",
                   fontsize=10)

    for ax in (ax_f, ax_l, ax_n):
        ax.grid(alpha=0.3)
        if ax is ax_n:
            ax.axhline(AG_STEP1D, color="C3", ls=":", lw=1.5)
        elif ax is ax_f:
            ax.axvline(AG_STEP1D, color="C3", ls=":", lw=1.5)
            ax.axvspan(2.5, 5.0, color="C3", alpha=0.06)
        else:
            ax.axhline(AG_STEP1D, color="C3", ls=":", lw=1.5)
            ax.axhspan(2.5, 5.0, color="C3", alpha=0.06)
    handles = [Line2D([], [], marker="o", ls="", color="C0", label="W-band (94 GHz) radar"),
               Line2D([], [], marker="o", ls="", color="C1", label="Ka-band (35 GHz) radar"),
               Line2D([], [], marker="o", ls="", color="k", mfc="white", label="outside the E-PROFILE CAMS box"),
               Line2D([], [], marker="o", ls="", color="0.5", alpha=0.35, label=f"faint: < {MIN_DAYS_SITE} days"),
               Line2D([], [], marker="x", ls="", color="k", label="(a) day-weighted median"),
               Line2D([], [], marker="|", ls="", color="k", ms=10, mew=2, label="(a) season-balanced median"),
               Line2D([], [], marker="D", ls="", color="C2", mec="k", label="(a) calibration profiles of a Vaisala unit"),
               Line2D([], [], ls="--", color="0.3", label="(a) median of the site medians (>= 20 days)"),
               Line2D([], [], ls=":", color="C3", label="step-1d a$_G$ = 3.3 $\\mu$m; shaded 2.5-5 $\\mu$m")]
    fig.legend(handles=handles, fontsize=8.5, loc="lower center", ncol=5, frameon=False)
    fig.suptitle("Cloudnet droplet size at the calibration cloud base across sites (lowest liquid base 0.25-2.4 km "
                 "AGL, status-1 pixels, scene-weighted medians unless stated, 2025-01..2026-10)", fontweight="bold")
    fig.subplots_adjust(left=0.13, right=0.99, top=0.92, bottom=0.14)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def write_compact_tables(out_dir, tab, cal):
    """The two tables of doc/reports/eta_cloudnet_sites_2026-11.md, generated (never typed by hand)."""
    L = ["| site | lat | alt (m) | radar (main), band | LWP from | days | scenes | a_G median [P25-P75] (um) | "
         "day-weighted | season-balanced | N_scaled median (cm-3) |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tab.iterrows():
        main = r["radar"].split(" (")[0]
        box = "" if r["eprofile_box"] else " *"
        thin = " (thin)" if r["n_days"] < MIN_DAYS_SITE else ""
        L.append(f"| {r['name']} ({r['country']}){box}{thin} | {r['lat']:.1f} | {r['alt_m']:.0f} | {main}, "
                 f"{r['band'].split(' ')[0]} | {r['lwp_source']} | {r['n_days']} | {r['n_scenes']} | "
                 f"**{_f(r['aG_med'])}** [{_f(r['aG_p25'], '{:.1f}')}-{_f(r['aG_p75'], '{:.1f}')}] | "
                 f"{_f(r['aG_daymed'])} | {_f(r['aG_season_bal'])} | {_f(r['N_med'], '{:.0f}')} |")
    L += ["", "| unit | distance (km) | accepted profiles (days) | with status-1 data (days) | "
          "a_G median [P25-P75] (um) | site, all scenes | site, same days | calibration / same days |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in cal.iterrows():
        L.append(f"| {r['unit']} | {r['distance_km']:.1f} | {r['n_profiles']} ({r['n_days']}) | {r['n_st1']} "
                 f"({r['n_days_st1']}) | **{_f(r['aG_med'])}** [{_f(r['aG_p25'], '{:.1f}')}-{_f(r['aG_p75'], '{:.1f}')}] | "
                 f"{_f(r.get('site_all_aG_med'))} | {_f(r.get('site_samedays_aG_med'))} | "
                 f"{_f(r.get('ratio_cal_to_samedays'))} |")
    (out_dir / "note_tables.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def version_effect(out_dir: Path):
    """Latest portal processing vs the July collection (step-1d files) on the sites present in both:
    all-scene site medians and calibration-profile medians, same code, same days."""
    import pandas as pd
    jd = out_dir / "july"
    if out_dir != OUT or not (jd / "site_table.csv").exists():
        return None
    new_s, old_s = pd.read_csv(out_dir / "site_table.csv"), pd.read_csv(jd / "site_table.csv")
    new_c, old_c = pd.read_csv(out_dir / "calib_table.csv"), pd.read_csv(jd / "calib_table.csv")
    L = ["# Effect of the portal reprocessing (latest files vs the July collection of step 1d)", "",
         "Same code. Days differ only where the latest record is longer (Jul-Oct 2026 added).", "",
         "| what | July files: a_G med (n) | latest files: a_G med (n) | latest, July days only | change |",
         "|---|---|---|---|---|"]
    for _, o in old_s.iterrows():
        n = new_s[new_s["site"] == o["site"]]
        if n.empty:
            continue
        n = n.iloc[0]
        z, _m = load_site(out_dir, o["site"])
        zj, _mj = load_site(jd, o["site"])
        ag = z["re_s_100"] * 1e6 * AG_FACTOR
        same = np.isin(z["date"], np.unique(zj["date"]))
        m_same = _q(ag[same], 50)
        L.append(f"| {o['name']}, all scenes | {o['aG_med']:.3f} ({o['n_scenes']}) | {n['aG_med']:.3f} "
                 f"({n['n_scenes']}) | {m_same:.3f} ({int(np.isfinite(ag[same]).sum())}) | "
                 f"{100 * (m_same / o['aG_med'] - 1):+.1f} % |")
    for _, o in old_c.iterrows():
        n = new_c[new_c["key"] == o["key"]]
        if n.empty or not np.isfinite(o["aG_med"]):
            continue
        n = n.iloc[0]
        L.append(f"| {o['unit']}, calibration profiles | {o['aG_med']:.3f} ({o['n_st1']}) | {n['aG_med']:.3f} "
                 f"({n['n_st1']}) | (same scenes) | {100 * (n['aG_med'] / o['aG_med'] - 1):+.1f} % |")
    (out_dir / "version_effect.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return L


def summarise_and_plot(out_dir: Path, use):
    import pandas as pd
    tab, season, radar = summarise_sites(out_dir, use)
    cal, pooled = summarise_calib(out_dir)
    stats = across_sites(tab, season, radar)
    pd.DataFrame(stats).to_csv(out_dir / "across_sites.csv", index=False, float_format="%.4g")
    write_markdown(out_dir, tab, season, radar, cal, pooled, stats)
    write_compact_tables(out_dir, tab, cal)
    version_effect(out_dir)
    cal_scenes = pd.read_csv(out_dir / "calib_scenes.csv", dtype={"date": str}, low_memory=False)
    plot_panels(out_dir, tab, cal_scenes, out_dir / "fig_aG_sites.png")
    plot_summary(tab, cal, stats, out_dir / "fig_aG_summary.png")
    print((out_dir / "site_table.md").read_text(encoding="utf-8"), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--files", choices=("latest", "july"), default="latest")
    ap.add_argument("--sites", default="all")
    ap.add_argument("--plot-only", action="store_true", help="tables + figures from the saved scenes")
    ap.add_argument("--copy-figs", action="store_true", help=f"also copy the figures to {REPORT_FIG_DIR}")
    a = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")       # site names such as Galați on a cp1252 Windows console
    out_dir = OUT if a.files == "latest" else OUT / "july"
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time_mod.time()
    if a.plot_only:
        use = sorted(p.stem for p in (out_dir / "scenes_all").glob("*.npz"))
        if a.sites != "all":
            use = [s for s in use if s in a.sites.split(",")]
    else:
        use = run(a, out_dir)
    summarise_and_plot(out_dir, use)
    if a.copy_figs:
        import shutil
        REPORT_FIG_DIR.mkdir(parents=True, exist_ok=True)
        for f in ("fig_aG_sites.png", "fig_aG_summary.png"):
            shutil.copy2(out_dir / f, REPORT_FIG_DIR / f)
    print(f"total runtime {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
