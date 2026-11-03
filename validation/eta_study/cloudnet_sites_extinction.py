"""Eta study: Cloudnet cloud extinction at the calibration cloud base, all `der` sites.

Companion of cloudnet_sites_droplet_size.py, which measured the droplet equivalent-area radius a_G at the
41 Cloudnet sites that produce the droplet effective radius product (`der`). This script measures the
second input of the multiple-scattering (eta) tables, the in-cloud extinction alpha, on EXACTLY the same
scenes, so that the two parameters are comparable site by site:

  * scene selection, liquid runs, layer rule, reliability flags and validity masks are imported from
    cloudnet_sites_droplet_size (`dsz`): every 30-s profile whose LOWEST liquid layer has its base
    0.25-2.4 km above the site (A, all scenes), and the accepted cloud-calibration profiles of the
    E-PROFILE Vaisala units at or near a Cloudnet site (B);
  * per pixel (der_retrieval_status == 1, 0 < der_scaled < 1 mm, N_scaled > 0) the extinction of a
    lognormal size distribution with Q_ext = 2,
        alpha = 2 pi N_scaled r_e^2 exp(-3 sigma^2) = 2 pi N_scaled a_G^2,   sigma = 0.35,
    then the median over the pixels of the lowest 50 / 100 / 150 m above the liquid base. This is the
    step-1d definition (cloudnet_microphysics_match._layer_stats, columns alpha_km_50/100/150 of
    C:/DATA/Projects/202611_eta_study/step1d/scenes_msfov.csv). The median of the per-pixel product is not
    the product of the layer medians of N and a_G; the latter is carried along (alpha_prodmed) only to
    size the difference.

The per-scene a_G and N_scaled of the re-read are checked against the droplet script's cache
(cloudnet_sites/scenes_all/<site>.npz): they must be identical, which proves that the scene selection
is the same. The re-read also counts, for every scene of the cloud-base window, the retrieval status of
the lowest-100-m pixels (1 reliable; 2 droplet pixel that also carries the categorize 'falling' bit, i.e.
drizzle, ice or radar echo above the lidar-detected liquid top; 3 rain flag of the profile), used in the
explanation of the site differences (cloudnet_sites_explain.py).

Run:  python validation/eta_study/cloudnet_sites_extinction.py [--workers 16] [--sites a,b] [--plot-only]
Outputs (C:/DATA/Projects/202611_eta_study/cloudnet_sites/):
  extinction_scenes/<site>.npz, <site>_days.csv   per-scene alpha / a_G / N, per-day status counts
  calib_scenes_extinction.csv                      B: per calibration profile
  site_table_extinction.csv                        one row per site
  season_table_extinction.csv, calib_table_extinction.csv, across_sites_extinction.csv,
  note_tables_extinction.md                        the report tables (generated)
  doc/reports/figs_04v2/fig_cloudnet_sites_extinction.png  (a) calibration units, (b) site medians
  doc/reports/figs_04v2/fig_cloudnet_sites_maps.png        maps of the site medians of a_G and alpha
Read-only on every input.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import sys
import time as time_mod
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cloudnet_sites_droplet_size as dsz  # noqa: E402  (same scene selection, same constants)

REPO = Path(__file__).resolve().parents[2]
OUT = dsz.OUT                                           # C:/DATA/Projects/202611_eta_study/cloudnet_sites
SCENE_DIR = OUT / "extinction_scenes"
STEP1D_CSV = Path("C:/DATA/Projects/202611_eta_study/step1d/scenes_msfov.csv")
FIG_DIR = REPO / "doc" / "reports" / "figs_04v2"

ALPHA_FACTOR = float(np.exp(-3.0 * dsz.SIGMA_LN ** 2))   # 0.692: a_G^2 = r_e^2 exp(-3 sigma^2)
ALPHA_ADOPTED = 20.0                                       # /km, generation-2 eta tables (decision note s.9)
# extinction profile of the calibration clouds measured in decision note section 9 (km-1)
ALPHA_PROFILE_S9 = {"onset": 5.0, "signal peak": 13.0, "100 m higher": 30.0}
STEP1D_UNITS = [("0-20000-0-06610_C", "Payerne CL61", "#7f3c8d"),
                ("0-20000-0-06610_B", "Payerne CL31", "#1f77b4"),
                ("0-250-1001-07151_A", "Palaiseau CL31", "#ff7f0e"),
                ("0-20000-0-10393_C", "Lindenberg CL61", "#2ca02c")]

A_FIELDS = ("t_epoch", "date", "base_agl", "depth", "re_s_100", "N_100", "alpha_50", "alpha_100",
            "alpha_150", "n_st1_100", "n_liq_100", "n_st3_100")


# =============================================================================
#  Core: layer statistics with the per-pixel extinction
# =============================================================================
def layer_stats(st, ders, nsc, h, lower, rows, jb, jt, kmax):
    """Medians over the status == 1 pixels of the lowest 50/100/150 m of the liquid runs (rows, jb, jt).

    Same gates, same masks as dsz.layer_medians (gate centre in [base_edge, base_edge + D), truncated at
    the top of the contiguous liquid run; 0 < r < 1 mm; N > 0). The extinction is formed PER PIXEL and its
    median taken over the layer, as in step 1d. Also counts the status-2 / status-3 pixels of the layer.
    """
    nz = h.size
    g = jb[:, None] + np.arange(kmax)[None, :]
    inrun = g <= jt[:, None]
    gc = np.minimum(g, nz - 1)
    rel = h[gc] - lower[jb][:, None]
    r2 = rows[:, None]
    s, vs, vn = st[r2, gc], ders[r2, gc], nsc[r2, gc]
    ok_s = (vs > 0) & (vs < 1e-3)
    ok_n = vn > 0
    st1 = s == 1
    with np.errstate(all="ignore"):
        alpha_px = 2.0 * np.pi * vn * vs ** 2 * ALPHA_FACTOR          # m-1
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)              # all-NaN rows -> NaN, on purpose
        for D in dsz.LAYERS_M:
            m = inrun & (rel < D)
            m1 = m & st1
            out[f"n_liq_{D}"] = m.sum(axis=1).astype(np.int16)
            out[f"n_st1_{D}"] = m1.sum(axis=1).astype(np.int16)
            out[f"re_s_{D}"] = np.nanmedian(np.where(m1 & ok_s, vs, np.nan), axis=1)
            out[f"alpha_{D}"] = np.nanmedian(np.where(m1 & ok_s & ok_n, alpha_px, np.nan), axis=1)
            if D == 100:
                out["N_100"] = np.nanmedian(np.where(m1 & ok_n, vn, np.nan), axis=1)
                out["n_st2_100"] = (m & (s == 2)).sum(axis=1).astype(np.int16)
                out["n_st3_100"] = (m & (s == 3)).sum(axis=1).astype(np.int16)
    return out


# =============================================================================
#  Worker: one (site, day) file -> A arrays, per-day status counts, B records
# =============================================================================
def process_day(job):
    """job = (site, ds, path, scenes), scenes = [(row_id, unit_key, t_epoch_s, cbh_m, alt_m), ...]
    (exactly the jobs of dsz.run)."""
    from netCDF4 import Dataset
    site, ds, path, scenes = job
    res = dict(site=site, ds=ds, ok=False)
    try:
        with Dataset(path) as nc:
            t_h = np.asarray(nc.variables["time"][:], float)
            h = np.asarray(nc.variables["height"][:], float)
            salt = float(np.nanmedian(dsz._filled(nc.variables["altitude"])))
            st = dsz._filled(nc.variables["der_retrieval_status"], fill=0, dtype=int)
            ders = dsz._filled(nc.variables["der_scaled"])
            nsc = dsz._filled(nc.variables["N_scaled"])
    except Exception as exc:  # noqa: BLE001 - a corrupt file must not kill the run
        res["error"] = f"{type(exc).__name__}: {exc}"
        return res
    if h.size < 2 or t_h.size == 0:
        res["error"] = "empty"
        return res
    lower, upper = dsz.gate_edges(h)
    kmax = int(np.ceil(max(dsz.LAYERS_M) / max(float(np.min(np.diff(h))), 1.0))) + 2
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    liq = (st >= 1) & (st <= 3)
    res.update(ok=True, salt=salt)

    # --- A: all profiles whose lowest liquid layer has its base in the calibration window ---
    anyl, jb, jt = dsz.lowest_liquid_run(liq)
    base_agl = lower[jb] - salt
    win = anyl & (base_agl >= dsz.CBH_WINDOW_M[0]) & (base_agl <= dsz.CBH_WINDOW_M[1])
    rows = np.flatnonzero(win)
    counts = dict(n_win=int(rows.size), n_win_st1=0, n_win_st2=0, n_win_st3=0, n_win_all3=0)
    a = {}
    if rows.size:
        ls = layer_stats(st, ders, nsc, h, lower, rows, jb[rows], jt[rows], kmax)
        # status of the lowest-100-m pixels of EVERY window scene (also those without a reliable pixel)
        counts.update(n_win_st1=int((ls["n_st1_100"] > 0).sum()), n_win_st2=int((ls["n_st2_100"] > 0).sum()),
                      n_win_st3=int((ls["n_st3_100"] > 0).sum()),
                      n_win_all3=int(((ls["n_st3_100"] > 0) & (ls["n_st3_100"] == ls["n_liq_100"])).sum()))
        keep = np.isfinite(ls["re_s_150"])       # dsz rule: any reliable radius within 150 m -> stored
        r = rows[keep]
        a = dict(t_epoch=(day0 + np.round(t_h[r] * 3600.0)).astype(np.int64),
                 date=np.full(r.size, int(ds), np.int32),
                 base_agl=base_agl[r].astype(np.float32),
                 depth=(upper[jt[r]] - lower[jb[r]]).astype(np.float32))
        for k in ("re_s_100", "N_100", "alpha_50", "alpha_100", "alpha_150"):
            a[k] = ls[k][keep].astype(np.float32)
        for k in ("n_st1_100", "n_liq_100", "n_st3_100"):
            a[k] = ls[k][keep]
    res["A"] = a
    res["counts"] = counts

    # --- B: accepted calibration profiles of the units at this site (dsz.process_day, B part) ---
    recs = []
    if scenes:
        hrs = np.array([(t - day0) / 3600.0 for _, _, t, _, _ in scenes])
        tis = np.clip(np.searchsorted(t_h, hrs), 1, t_h.size - 1)
        tis = np.where(np.abs(t_h[tis - 1] - hrs) <= np.abs(t_h[tis] - hrs), tis - 1, tis)
        if t_h.size == 1:
            tis = np.zeros_like(tis)
        b_rows, b_jb, b_jt, b_idx = [], [], [], []
        for (rid, _key, _t, cbh, alt), hh, ti in zip(scenes, hrs, tis):
            rec = dict(rid=rid)
            if abs(t_h[ti] - hh) > dsz.MAX_DT_H:
                rec["status"] = "notime"
                recs.append(rec)
                continue
            cloud_msl = cbh + (salt if alt is None else alt)
            rec["status"] = "matched"
            rn = dsz.run_near(liq[ti], h, cloud_msl)
            if rn is None:
                rec["liquid"] = False
            else:
                rec["liquid"] = True
                b_rows.append(ti)
                b_jb.append(rn[0])
                b_jt.append(rn[1])
                b_idx.append(len(recs))
            recs.append(rec)
        if b_rows:
            ls = layer_stats(st, ders, nsc, h, lower, np.array(b_rows), np.array(b_jb), np.array(b_jt), kmax)
            for j, i in enumerate(b_idx):
                for k in ("re_s_100", "N_100", "alpha_50", "alpha_100", "alpha_150", "n_st1_100"):
                    recs[i][k] = float(ls[k][j])
    res["B"] = recs
    return res


# =============================================================================
#  Driver
# =============================================================================
def build_jobs(sites_arg: str):
    """The jobs of dsz.run (same files, same calibration scenes grouped per site-day)."""
    import pandas as pd
    sites = dsz.load_sites()
    use = sorted(sites) if sites_arg == "all" else [s.strip() for s in sites_arg.split(",")]
    alts = dsz.unit_altitudes()
    sc = pd.DataFrame(dsz.scenes_msfov(dsz.MSFOV_JSONL) + dsz.scenes_dump(dsz.DUMP_DIR),
                      columns=["key", "date", "t_epoch", "cbh_m"])
    earliest = sc["key"].map(lambda k: dsz.UNITS[k][5] or "00000000")
    sc = sc[sc["date"] >= earliest].reset_index(drop=True)
    sc["cn_site"] = sc["key"].map(lambda k: dsz.UNITS[k][1])
    by_day = defaultdict(list)
    for rid, k, d, t, z, s in zip(sc.index, sc["key"], sc["date"], sc["t_epoch"], sc["cbh_m"], sc["cn_site"]):
        by_day[(s, d)].append((int(rid), k, int(t), float(z), alts[k]))
    jobs = []
    for s in use:
        for ds, p in sorted(dsz.der_files(s, "latest").items()):
            jobs.append((s, ds, str(p), by_day.pop((s, ds), [])))
    return use, sc, jobs


def run(args):
    import pandas as pd
    use, sc, jobs = build_jobs(args.sites)
    print(f"{len(jobs)} site-days on {len(use)} sites, {len(sc)} calibration profiles -> {args.workers} processes",
          flush=True)
    acc = defaultdict(lambda: defaultdict(list))
    day_rows = defaultdict(list)
    b_recs, errors = {}, []
    t0 = time_mod.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for k, res in enumerate(ex.map(process_day, jobs, chunksize=4), 1):
            s = res["site"]
            if not res["ok"]:
                errors.append(f"{s} {res['ds']} {res.get('error')}")
            else:
                day_rows[s].append(dict(date=res["ds"], salt=res["salt"], **res["counts"]))
                for f, v in res["A"].items():
                    acc[s][f].append(v)
                for rec in res["B"]:
                    b_recs[rec.pop("rid")] = rec
            if k % 2000 == 0:
                print(f"  {k}/{len(jobs)} site-days ({time_mod.time() - t0:.0f} s)", flush=True)
    print(f"read in {time_mod.time() - t0:.0f} s, {len(errors)} unreadable files", flush=True)
    SCENE_DIR.mkdir(parents=True, exist_ok=True)
    if errors:
        (SCENE_DIR / "unreadable_files.txt").write_text("\n".join(errors) + "\n", encoding="utf-8")
    for s in use:
        if not day_rows[s]:
            continue
        pd.DataFrame(day_rows[s]).to_csv(SCENE_DIR / f"{s}_days.csv", index=False)
        arrays = {f: (np.concatenate(acc[s][f]) if acc[s][f] else np.array([])) for f in A_FIELDS}
        np.savez_compressed(SCENE_DIR / f"{s}.npz", **arrays)
    b = sc.join(pd.DataFrame.from_dict(b_recs, orient="index"))
    b["status"] = b["status"].fillna("nofile")
    b["label"] = b["key"].map(lambda k: dsz.UNITS[k][0])
    b.to_csv(OUT / "calib_scenes_extinction.csv", index=False, float_format="%.5g")
    return use


# =============================================================================
#  Checks and summaries
# =============================================================================
def check_against_droplet_cache(use) -> list[str]:
    """The re-read must reproduce the droplet script's per-scene a_G and N exactly (same scenes)."""
    lines = []
    for s in use:
        fz, fd = SCENE_DIR / f"{s}.npz", OUT / "scenes_all" / f"{s}.npz"
        if not (fz.exists() and fd.exists()):
            continue
        z, d = np.load(fz), np.load(fd)
        same_n = z["t_epoch"].size == d["t_epoch"].size
        ok = same_n and np.array_equal(z["t_epoch"], d["t_epoch"]) \
            and np.allclose(z["re_s_100"], d["re_s_100"], equal_nan=True, rtol=0, atol=0) \
            and np.allclose(z["N_100"], d["N_100"], equal_nan=True, rtol=0, atol=0)
        lines.append(f"{s}: {z['t_epoch'].size} scenes, droplet cache {d['t_epoch'].size} -> "
                     f"{'IDENTICAL' if ok else 'DIFFERENT'}")
    return lines


def _q(x, p):
    return dsz._q(x, p)


def load_scenes(s: str):
    fp = SCENE_DIR / f"{s}.npz"
    return dict(np.load(fp)) if fp.exists() else None


def summarise_sites(use):
    """One row per site: extinction (km-1) next to a_G and N_scaled of the same scenes."""
    import pandas as pd
    meta = pd.read_csv(OUT / "site_table.csv", dtype={"first": str, "last": str})
    meta = meta.set_index("site")
    rows, season_rows = [], []
    for s in use:
        z = load_scenes(s)
        if z is None or s not in meta.index:
            continue
        m = meta.loc[s]
        days = pd.read_csv(SCENE_DIR / f"{s}_days.csv")
        alpha = z["alpha_100"] * 1e3                         # km-1
        ag = z["re_s_100"] * 1e6 * dsz.AG_FACTOR             # um
        n_cm3 = z["N_100"] * 1e-6
        ok = np.isfinite(alpha)
        date = z["date"]
        a_day, n_days = dsz.daily_median(date, alpha)
        n_win = max(int(days["n_win"].sum()), 1)
        row = dict(site=s, name=m["name"], country=m["country"], lat=m["lat"], lon=m["lon"], alt_m=m["alt_m"],
                   eprofile_box=bool(m["eprofile_box"]), radar=m["radar"], band=m["band"], lidar=m["lidar"],
                   lwp_source=m["lwp_source"], first=m["first"], last=m["last"],
                   n_days=n_days, n_scenes=int(ok.sum()),
                   alpha_med=_q(alpha, 50), alpha_p25=_q(alpha, 25), alpha_p75=_q(alpha, 75),
                   alpha_p10=_q(alpha, 10), alpha_p90=_q(alpha, 90), alpha_daymed=a_day,
                   alpha50_med=_q(z["alpha_50"] * 1e3, 50), alpha150_med=_q(z["alpha_150"] * 1e3, 50),
                   aG_med=_q(ag[ok], 50), aG_p25=_q(ag[ok], 25), aG_p75=_q(ag[ok], 75),
                   N_med=_q(n_cm3[ok], 50), N_p25=_q(n_cm3[ok], 25), N_p75=_q(n_cm3[ok], 75),
                   aG_med_site_table=m["aG_med"], N_med_site_table=m["N_med"],
                   depth_med=_q(z["depth"][ok], 50), base_agl_med=_q(z["base_agl"][ok], 50),
                   n_win=int(days["n_win"].sum()),
                   frac_st1_100=days["n_win_st1"].sum() / n_win, frac_st2_100=days["n_win_st2"].sum() / n_win,
                   frac_st3_100=days["n_win_st3"].sum() / n_win, frac_all3_100=days["n_win_all3"].sum() / n_win)
        # product of the layer medians: 2 pi N a_G^2 from the SITE medians (not the step-1d definition)
        row["alpha_prodmed"] = 2.0 * np.pi * row["N_med"] * 1e6 * (row["aG_med"] * 1e-6) ** 2 * 1e3
        month = (date // 100) % 100
        seas = dsz._season(month) if date.size else np.array([])
        seasonal = []
        for sea in ("DJF", "MAM", "JJA", "SON"):
            sel = seas == sea
            dm, nd = dsz.daily_median(date[sel], alpha[sel])
            season_rows.append(dict(site=s, season=sea, alpha_daymed=dm, n_days=nd,
                                    n_scenes=int(np.isfinite(alpha[sel]).sum()), alpha_med=_q(alpha[sel], 50),
                                    aG_med=_q(ag[sel & ok], 50), N_med=_q(n_cm3[sel & ok], 50)))
            if nd >= dsz.MIN_DAYS_SEASON:
                seasonal.append(dm)
        row["n_seasons"] = len(seasonal)
        row["alpha_season_bal"] = float(np.median(seasonal)) if len(seasonal) == 4 else np.nan
        rows.append(row)
    tab = pd.DataFrame(rows).sort_values("lat", ascending=False).reset_index(drop=True)
    tab.to_csv(OUT / "site_table_extinction.csv", index=False, float_format="%.6g")
    season = pd.DataFrame(season_rows)
    season.to_csv(OUT / "season_table_extinction.csv", index=False, float_format="%.6g")
    return tab, season


def summarise_calib(tab):
    """B per unit (41-site framework), next to the all-scene value of its site on the same days, and the
    step-1d values of the four original units (the numbers of decision note section 9)."""
    import pandas as pd
    b = pd.read_csv(OUT / "calib_scenes_extinction.csv", dtype={"date": str}, low_memory=False)
    s1 = pd.read_csv(STEP1D_CSV, usecols=["key", "date", "n_st1_100", "aG_um_100", "alpha_km_50", "alpha_km_100",
                                          "alpha_km_150"], low_memory=False)
    s1 = s1[(s1["n_st1_100"] >= 1) & np.isfinite(s1["aG_um_100"])]
    rows = []
    for key, (label, site, itype, dist, src, _e) in dsz.UNITS.items():
        g = b[(b["key"] == key) & (b["status"] == "matched") & (b["liquid"] == True)]  # noqa: E712
        alpha = g["alpha_100"].to_numpy(float) * 1e3
        ok = np.isfinite(alpha)
        ag = g["re_s_100"].to_numpy(float) * 1e6 * dsz.AG_FACTOR
        dmed, nd = dsz.daily_median(g["date"].astype(int).to_numpy(), alpha)
        row = dict(key=key, unit=label, type=itype, site=site, distance_km=dist, n_profiles=int(ok.sum()),
                   n_days=nd, alpha_med=_q(alpha, 50), alpha_p25=_q(alpha, 25), alpha_p75=_q(alpha, 75),
                   alpha_p10=_q(alpha, 10), alpha_p90=_q(alpha, 90), alpha_daymed=dmed,
                   alpha50_med=_q(g["alpha_50"].to_numpy(float) * 1e3, 50),
                   alpha150_med=_q(g["alpha_150"].to_numpy(float) * 1e3, 50),
                   aG_med=_q(ag[ok], 50), N_med=_q(g["N_100"].to_numpy(float)[ok] * 1e-6, 50))
        z = load_scenes(site)
        if z is not None:
            a_all = z["alpha_100"] * 1e3
            on_days = np.isin(z["date"], g.loc[ok, "date"].astype(int).unique())
            row["site_all_alpha_med"] = _q(a_all, 50)
            row["site_samedays_alpha_med"] = _q(a_all[on_days], 50)
            row["ratio_cal_to_samedays"] = row["alpha_med"] / row["site_samedays_alpha_med"]
        g1 = s1[s1["key"] == key]
        if len(g1):
            row.update(step1d_n=len(g1), step1d_alpha50_med=float(np.nanmedian(g1["alpha_km_50"])),
                       step1d_alpha100_med=float(np.nanmedian(g1["alpha_km_100"])),
                       step1d_alpha150_med=float(np.nanmedian(g1["alpha_km_150"])),
                       step1d_alpha100_p10=_q(g1["alpha_km_100"], 10), step1d_alpha100_p90=_q(g1["alpha_km_100"], 90),
                       step1d_aG100_med=float(np.nanmedian(g1["aG_um_100"])))
        rows.append(row)
    rows.append(dict(key="step1d_pooled", unit="step 1d, four units pooled", n_profiles=len(s1),
                     step1d_n=len(s1), step1d_alpha50_med=float(np.nanmedian(s1["alpha_km_50"])),
                     step1d_alpha100_med=float(np.nanmedian(s1["alpha_km_100"])),
                     step1d_alpha150_med=float(np.nanmedian(s1["alpha_km_150"])),
                     step1d_alpha100_p10=_q(s1["alpha_km_100"], 10), step1d_alpha100_p90=_q(s1["alpha_km_100"], 90),
                     step1d_aG100_med=float(np.nanmedian(s1["aG_um_100"]))))
    cal = pd.DataFrame(rows)
    cal.to_csv(OUT / "calib_table_extinction.csv", index=False, float_format="%.6g")
    return cal


def across_sites(tab):
    """Spread of the site medians of alpha (E-PROFILE domain, >= 20 days; and all sites)."""
    import pandas as pd
    from scipy.stats import spearmanr
    L = []
    for name, sel in (("E-PROFILE domain (CAMS box N27-80, W30-E45)", tab["eprofile_box"]),
                      ("all sites", tab["eprofile_box"] | ~tab["eprofile_box"])):
        g0 = tab[sel & (tab["n_days"] >= dsz.MIN_DAYS_SITE)]
        for col, lab in (("alpha_med", "scene-weighted"), ("alpha_daymed", "day-weighted"),
                         ("alpha_season_bal", "season-balanced (4 seasons)"), ("alpha50_med", "lowest 50 m"),
                         ("alpha150_med", "lowest 150 m"), ("alpha_prodmed", "2 pi N_med a_G,med^2")):
            g = g0[np.isfinite(g0[col].to_numpy(float))]
            if len(g) < 3:
                continue
            v = g[col].to_numpy(float)
            rho_l, p_l = spearmanr(g["lat"].abs(), v)
            rho_a, p_a = spearmanr(g["alt_m"], v)
            rho_g, p_g = spearmanr(g["aG_med"], v)
            rho_n, p_n = spearmanr(g["N_med"], v)
            L.append(dict(group=name, weighting=lab, n_sites=len(g), median_of_sites=float(np.median(v)),
                          p25=_q(v, 25), p75=_q(v, 75), min=float(v.min()), max=float(v.max()),
                          min_site=g.loc[g[col].idxmin(), "site"], max_site=g.loc[g[col].idxmax(), "site"],
                          geo_sd_factor=float(np.exp(np.std(np.log(v), ddof=1))),
                          n_below_10=int((v < 10).sum()), n_10_30=int(((v >= 10) & (v <= 30)).sum()),
                          n_above_30=int((v > 30).sum()), n_below_20=int((v < ALPHA_ADOPTED).sum()),
                          frac_within_25pct_of_20=float(np.mean(np.abs(v / ALPHA_ADOPTED - 1) <= 0.25)),
                          spearman_abslat=float(rho_l), p_abslat=float(p_l), spearman_alt=float(rho_a),
                          p_alt=float(p_a), spearman_aG=float(rho_g), p_aG=float(p_g), spearman_N=float(rho_n),
                          p_N=float(p_n)))
    out = pd.DataFrame(L)
    out.to_csv(OUT / "across_sites_extinction.csv", index=False, float_format="%.6g")
    return out


# =============================================================================
#  Figures (landscape)
# =============================================================================
def plot_extinction(tab, path):
    """(a) alpha at the calibration profiles of the four step-1d Vaisala units (step-1d match, as panel (a)
    of the paper's droplet figure); (b) site medians and interquartile ranges of the E-PROFILE-domain
    sites with >= 20 days."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    sc = pd.read_csv(STEP1D_CSV, usecols=["key", "n_st1_100", "aG_um_100", "alpha_km_100"], low_memory=False)
    sc = sc[(sc["n_st1_100"] >= 1) & np.isfinite(sc["aG_um_100"])]
    pooled = sc["alpha_km_100"].to_numpy(float)
    pooled = pooled[np.isfinite(pooled)]
    st = tab[tab["eprofile_box"] & (tab["n_days"] >= dsz.MIN_DAYS_SITE)].sort_values("alpha_med").reset_index(drop=True)
    med_sites = float(np.median(st["alpha_med"]))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(17, 5.6), dpi=200, gridspec_kw=dict(width_ratios=[1, 1.5]))
    bins = np.logspace(np.log10(0.5), np.log10(2000.0), 37)       # 0.1 decade per bin
    for key, lab, col in STEP1D_UNITS:
        v = sc.loc[sc["key"] == key, "alpha_km_100"].to_numpy(float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        # fraction per bin, not a density: log-spaced bins would otherwise distort the tails
        ax.hist(np.clip(v, bins[0], bins[-1] * 0.999), bins, histtype="step", lw=1.8, color=col,
                weights=np.full(v.size, 1.0 / v.size),
                label=f"{lab} (n = {v.size}, median {np.median(v):.0f} km$^{{-1}}$)")
    ax.axvline(ALPHA_ADOPTED, color="k", ls="--", lw=1.4,
               label=f"adopted $\\alpha$ = {ALPHA_ADOPTED:.0f} km$^{{-1}}$ (four units pooled: median "
                     f"{np.median(pooled):.0f})")
    for i, (where, val) in enumerate(ALPHA_PROFILE_S9.items()):
        ax.axvline(val, color="0.45", ls=":", lw=1.1,
                   label="ceilometer profile (s. 9): 5 / 13 / 30 km$^{-1}$" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 0.22)
    ax.set_xlabel("Cloudnet extinction $\\alpha$ = 2$\\pi$ N$_{scaled}$ $a_G^2$ in the lowest 100 m above the "
                  "liquid base [km$^{-1}$]")
    ax.set_ylabel("fraction of calibration profiles per 0.1-decade bin")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.2, loc="upper left", handlelength=1.6)
    ax.set_title("(a) Calibration profiles at the three Cloudnet sites with a co-located ceilometer", fontsize=11)

    y = np.arange(len(st))
    left = 2.0
    bx.barh(y, st["alpha_med"] - left, left=left, color="#9ecae1", edgecolor="0.3", lw=0.5)
    bx.errorbar(st["alpha_med"], y, xerr=[st["alpha_med"] - st["alpha_p25"], st["alpha_p75"] - st["alpha_med"]],
                fmt="none", ecolor="0.25", elinewidth=0.8, capsize=2)
    bx.axvline(ALPHA_ADOPTED, color="k", ls="--", lw=1.4, label=f"adopted $\\alpha$ = {ALPHA_ADOPTED:.0f} km$^{{-1}}$")
    bx.axvline(med_sites, color="tab:red", ls=":", lw=1.6,
               label=f"median of the site medians ({med_sites:.1f} km$^{{-1}}$)")
    bx.set_xscale("log")
    bx.set_xlim(left, 300.0)
    bx.set_yticks(y)
    bx.set_yticklabels([f"{n} ({int(d)} d)" for n, d in zip(st["name"], st["n_days"])], fontsize=7)
    bx.set_ylim(-0.7, len(st) - 0.3)
    bx.set_xlabel("site median of $\\alpha$ (all liquid scenes, lowest 100 m above the base) [km$^{-1}$], "
                  "bars: interquartile range")
    bx.grid(alpha=0.3, axis="x", which="both")
    bx.legend(fontsize=8, loc="lower right")
    bx.set_title(f"(b) {len(st)} Cloudnet sites of the E-PROFILE domain, 2025–2026")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return dict(med_sites=med_sites, pooled_med=float(np.median(pooled)), n_pooled=int(pooled.size))


# Maps: each panel = the whole E-PROFILE box + a zoom on the crowded central-European sites.
ZOOM = dict(west=-2.5, east=17.5, south=45.5, north=54.0)
# sites less than 10 km apart would hide each other's colour: displaced (degrees lon, lat) for visibility
DISPLACE_MAIN = {"pallasjarvi": (-2.4, 0.0), "sammaltunturinlaatas": (2.4, 0.0)}
DISPLACE_ZOOM = {"leipzig": (-0.28, 0.15), "leipzig-lim": (0.0, -0.28), "lindenberg": (0.0, 0.19),
                 "falkenberg": (0.0, -0.22)}
LABEL_MAIN = {"ny-alesund": "bottom center", "snaefellsnes": "top center", "pallasjarvi": "middle left",
              "kenttarova": "top center", "sammaltunturinlaatas": "middle right", "hyytiala": "top right",
              "norunda": "middle left", "cluj": "top center", "galati": "middle right", "bucharest": "bottom center",
              "potenza": "middle right", "laquila": "top left", "lampedusa": "bottom center",
              "granada": "bottom center", "limassol": "bottom center", "arm-graciosa": "bottom center"}
LABEL_ZOOM = {"leipzig": "middle left", "leipzig-lim": "bottom center", "melpitz": "middle right",
              "lindenberg": "middle right", "falkenberg": "middle right", "juelich": "middle left",
              "cabauw": "top center", "kappelrodeck": "middle left", "munich": "middle right",
              "schneefernerhaus": "bottom center", "payerne": "middle left", "palaiseau": "middle left",
              "chilbolton": "top center", "rzecin": "top center"}


def _in_zoom(lon, lat):
    return (ZOOM["west"] <= lon <= ZOOM["east"]) and (ZOOM["south"] <= lat <= ZOOM["north"])


def plot_maps(tab, path):
    """Two panels side by side: site median of a_G (a) and of alpha (b) on the E-PROFILE box, each with a
    zoom on central Europe. plotly scattergeo + kaleido (cartopy is not installed)."""
    import plotly.graph_objects as go
    box = dsz.EPROFILE_BOX
    sel = tab[tab["eprofile_box"]].copy().reset_index(drop=True)
    alpha_ticks = np.array([5, 7, 10, 15, 20, 30, 50, 70])
    # (geo id, zoom?, values, colour scale, range, ticks, colour-bar x, colour-bar title)
    a_val = sel["aG_med"].to_numpy(float)
    l_val = np.log10(sel["alpha_med"].to_numpy(float))
    layout_geo = {}
    panels = [("geo", False, a_val, "Viridis", (1.8, 5.6), None, None, None),
              ("geo2", True, a_val, "Viridis", (1.8, 5.6), None, 0.472, "a<sub>G</sub> (µm)"),
              ("geo3", False, l_val, "Plasma", (np.log10(5.0), np.log10(70.0)), alpha_ticks, None, None),
              ("geo4", True, l_val, "Plasma", (np.log10(5.0), np.log10(70.0)), alpha_ticks, 0.99,
               "α (km<sup>-1</sup>)")]
    domains = {"geo": (0.0, 0.27), "geo2": (0.275, 0.465), "geo3": (0.515, 0.785), "geo4": (0.79, 0.98)}
    fig = go.Figure()
    for gid, zoom, val, cmap, (vmin, vmax), ticks, xbar, ctitle in panels:
        lon = sel["lon"].to_numpy(float).copy()
        lat = sel["lat"].to_numpy(float).copy()
        disp = DISPLACE_ZOOM if zoom else DISPLACE_MAIN
        for i, s in enumerate(sel["site"]):
            if s in disp:
                lon[i] += disp[s][0]
                lat[i] += disp[s][1]
        inside = np.array([_in_zoom(x, y) for x, y in zip(sel["lon"], sel["lat"])])
        show = inside if zoom else np.ones(len(sel), bool)
        labels = [n if (zoom or not z) else "" for n, z in zip(sel["name"], inside)]
        pos_map = LABEL_ZOOM if zoom else LABEL_MAIN
        marker = dict(size=13 if zoom else 11, color=val[show], colorscale=cmap, cmin=vmin, cmax=vmax,
                      line=dict(width=1.0, color="black"), showscale=xbar is not None)
        if xbar is not None:
            cbar = dict(title=dict(text=ctitle, side="right"), x=xbar, len=0.8, thickness=14, y=0.48)
            if ticks is not None:
                cbar.update(tickvals=np.log10(ticks).tolist(), ticktext=[str(t) for t in ticks])
            marker["colorbar"] = cbar
        fig.add_trace(go.Scattergeo(lon=lon[show], lat=lat[show], mode="markers+text",
                                    text=[t for t, k in zip(labels, show) if k],
                                    textposition=[pos_map.get(s, "top center") for s, k in zip(sel["site"], show) if k],
                                    textfont=dict(size=10 if zoom else 9, color="#1a1a1a"), marker=marker,
                                    geo=gid, showlegend=False))
        if not zoom:     # outline of the zoom area on the main map
            zl = np.r_[np.linspace(ZOOM["west"], ZOOM["east"], 30), np.full(10, ZOOM["east"]),
                       np.linspace(ZOOM["east"], ZOOM["west"], 30), np.full(10, ZOOM["west"])]
            za = np.r_[np.full(30, ZOOM["south"]), np.linspace(ZOOM["south"], ZOOM["north"], 10),
                       np.full(30, ZOOM["north"]), np.linspace(ZOOM["north"], ZOOM["south"], 10)]
            fig.add_trace(go.Scattergeo(lon=zl, lat=za, mode="lines", line=dict(width=1.2, color="#333333", dash="dot"),
                                        geo=gid, showlegend=False, hoverinfo="skip"))
        rng = ZOOM if zoom else box
        layout_geo[gid] = dict(
            domain=dict(x=list(domains[gid]), y=[0.02, 0.92]),
            projection=dict(type="conic conformal", parallels=[35, 65], rotation=dict(lon=7.5)),
            lonaxis=dict(range=[rng["west"], rng["east"]], showgrid=True, dtick=5 if zoom else 10, gridcolor="#d0d0d0"),
            lataxis=dict(range=[rng["south"], rng["north"]], showgrid=True, dtick=5 if zoom else 10, gridcolor="#d0d0d0"),
            resolution=50, showcoastlines=True, coastlinecolor="#555555", coastlinewidth=0.7, showcountries=True,
            countrycolor="#999999", countrywidth=0.5, showland=True, landcolor="#f2f0eb", showocean=True,
            oceancolor="#e6eef5", showlakes=False, showframe=True, framecolor="#666666")
    fig.update_layout(**layout_geo, width=2800, height=770, margin=dict(l=10, r=10, t=55, b=30),
                      paper_bgcolor="white", font=dict(family="Arial", size=12),
                      title=dict(text="Cloudnet droplet radius and extinction in the lowest 100 m above the liquid "
                                      "cloud base (base 0.25-2.4 km above the site, reliable pixels, 2025-01 .. "
                                      "2026-10): the 30 sites inside the E-PROFILE box", x=0.5, y=0.975, yref="container",
                                 font=dict(size=17)))
    for x, txt in ((0.232, "(a) site median of the droplet radius a<sub>G</sub> (µm)"),
                   (0.747, "(b) site median of the extinction α (km<sup>-1</sup>, log colour scale)")):
        fig.add_annotation(text=txt, x=x, y=0.96, xref="paper", yref="paper", showarrow=False, font=dict(size=15))
    fig.add_annotation(text="Dotted frame on each map = the zoom to its right. Markers of sites less than 10 km "
                            "apart (Kenttärova / Pallasjärvi / Sammaltunturinlaatas, Leipzig / Leipzig LIM, "
                            "Lindenberg / Falkenberg) are displaced slightly so that each colour is visible. "
                            "The 11 Cloudnet sites outside the box are in the table only.",
                       x=0.5, y=-0.03, xref="paper", yref="paper", showarrow=False, font=dict(size=12))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_image(str(path), width=2800, height=770, scale=1)


def _fmt(x, fmt="{:.1f}"):
    try:
        return "" if not np.isfinite(float(x)) else fmt.format(float(x))
    except (TypeError, ValueError):
        return str(x)


def write_note_tables(tab, cal):
    """The tables of doc/reports/eta_cloudnet_sites_extinction_2026-11.md, generated (never typed by hand)."""
    box = tab[tab["eprofile_box"] & (tab["n_days"] >= dsz.MIN_DAYS_SITE)].sort_values("lat", ascending=False)
    L = ["| site | lat | alt (m) | days | scenes | alpha median [P25-P75] (/km) | day-weighted | season-balanced | "
         "lowest 50 / 150 m | 2 pi N_med a_G,med^2 | a_G (um) | N_scaled (cm-3) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in box.iterrows():
        L.append(f"| {r['name']} ({r['country']}) | {r['lat']:.1f} | {r['alt_m']:.0f} | {r['n_days']} | {r['n_scenes']} | "
                 f"**{_fmt(r['alpha_med'])}** [{_fmt(r['alpha_p25'])}-{_fmt(r['alpha_p75'])}] | {_fmt(r['alpha_daymed'])} | "
                 f"{_fmt(r['alpha_season_bal'])} | {_fmt(r['alpha50_med'])} / {_fmt(r['alpha150_med'])} | "
                 f"{_fmt(r['alpha_prodmed'])} | {_fmt(r['aG_med'], '{:.2f}')} | {_fmt(r['N_med'], '{:.0f}')} |")
    out = tab[~tab["eprofile_box"]].sort_values("lat", ascending=False)
    L += ["", "| site outside the E-PROFILE box | lat | alt (m) | days | alpha median [P25-P75] (/km) | a_G (um) | N (cm-3) |",
          "|---|---|---|---|---|---|---|"]
    for _, r in out.iterrows():
        L.append(f"| {r['name']} ({r['country']}) | {r['lat']:.1f} | {r['alt_m']:.0f} | {r['n_days']} | "
                 f"{_fmt(r['alpha_med'])} [{_fmt(r['alpha_p25'])}-{_fmt(r['alpha_p75'])}] | {_fmt(r['aG_med'], '{:.2f}')} | "
                 f"{_fmt(r['N_med'], '{:.0f}')} |")
    L += ["", "| unit | profiles with status-1 data (days) | alpha median [P25-P75] (/km) | lowest 50 / 150 m | "
          "site, all scenes | site, same days | calibration / same days | step 1d: 50 / 100 / 150 m (profiles) |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in cal.iterrows():
        if r["key"] == "step1d_pooled":
            L.append(f"| {r['unit']} | | | | | | | {_fmt(r['step1d_alpha50_med'])} / {_fmt(r['step1d_alpha100_med'])} / "
                     f"{_fmt(r['step1d_alpha150_med'])} ({int(r['step1d_n'])}) |")
            continue
        s1 = "" if not np.isfinite(float(r.get("step1d_n", np.nan))) else \
            (f"{_fmt(r['step1d_alpha50_med'])} / {_fmt(r['step1d_alpha100_med'])} / {_fmt(r['step1d_alpha150_med'])} "
             f"({int(r['step1d_n'])})")
        L.append(f"| {r['unit']} | {r['n_profiles']} ({_fmt(r['n_days'], '{:.0f}')}) | **{_fmt(r['alpha_med'])}** "
                 f"[{_fmt(r['alpha_p25'])}-{_fmt(r['alpha_p75'])}] | {_fmt(r['alpha50_med'])} / {_fmt(r['alpha150_med'])} | "
                 f"{_fmt(r.get('site_all_alpha_med'))} | {_fmt(r.get('site_samedays_alpha_med'))} | "
                 f"{_fmt(r.get('ratio_cal_to_samedays'), '{:.2f}')} | {s1} |")
    (OUT / "note_tables_extinction.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def summarise_and_plot(use):
    tab, _season = summarise_sites(use)
    cal = summarise_calib(tab)
    acr = across_sites(tab)
    write_note_tables(tab, cal)
    fig_stats = plot_extinction(tab, FIG_DIR / "fig_cloudnet_sites_extinction.png")
    plot_maps(tab, FIG_DIR / "fig_cloudnet_sites_maps.png")
    cols = ["name", "lat", "alt_m", "n_days", "n_scenes", "alpha_med", "alpha_p25", "alpha_p75", "alpha_daymed",
            "alpha_season_bal", "alpha50_med", "alpha150_med", "alpha_prodmed", "aG_med", "N_med", "frac_st3_100"]
    import pandas as pd
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_rows", 100):
        print(tab[cols].round(2).to_string(index=False))
        print(cal.round(2).to_string(index=False))
        print(acr.round(3).to_string(index=False))
    print(fig_stats)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--sites", default="all")
    ap.add_argument("--plot-only", action="store_true", help="tables + figures from the saved scenes")
    a = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")       # site names such as Galați on a cp1252 Windows console
    t0 = time_mod.time()
    if a.plot_only:
        use = sorted(p.stem for p in SCENE_DIR.glob("*.npz"))
        if a.sites != "all":
            use = [s for s in use if s in a.sites.split(",")]
    else:
        use = run(a)
        for line in check_against_droplet_cache(use):
            print(line, flush=True)
    summarise_and_plot(use)
    print(f"total runtime {time_mod.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
