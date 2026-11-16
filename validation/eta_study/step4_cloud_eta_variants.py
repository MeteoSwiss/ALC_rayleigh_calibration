"""Step 4 of the eta study: which multiple-scattering correction makes the CL61 LIQUID-CLOUD
lidar constant agree with the CL61's own RAYLEIGH constant?

For every CL61 unit (Payerne, Lindenberg from 2025-05, and the 9 others as a secondary check) the
operational cloud calibration is re-run EXACTLY (same loader, WV 910.55/0.188 on 0.4 deg CAMS, aerosol
transmission correction, day aggregation through calibration.flags.cloud_flag and the runner's success
set) once per eta-table variant, on the unit's known cloud days:
  a      operational PVC table (a_G 5.5 um, alpha 10 /km)                    = _filters._ETA_CL61
  a_re   the same PVC recomputed here (a_G 5.5, alpha 10)                    (sanity: must match a)
  b      legacy Hopkin ladder (git 83f451f^; the CL61 borrowed the CL51 table)
  c      PVC a_G 3.3 um, alpha = CL61 type median alpha_k (step 2 tail slope)
  d      PVC a_G 3.3 um, alpha = this unit's median alpha_k
  e      = c with the Mie lidar ratio S(r_e = 3.3/0.832 um, sigma 0.35, 910.55 nm) instead of 18.8 sr
         (post hoc: C_L,e = C_L,c * S_mie / 18.8)
  f_lo / f_hi   PVC a_G 3.3, alpha = CL61 alpha_k P10 / P90
  f_ag25 / f_ag50  PVC a_G 2.5 / 5.0 um, alpha = CL61 alpha_k median
Reference: the CL61 Rayleigh v2.2 at 910.55 nm. Payerne: WITH the measured dark (diag_v22_l55s008dark,
primary) and without (diag_v22_l55s008, secondary). Lindenberg: without dark only (no hood measurement
there). Other CL61: without dark (step4/rayleigh_v22_l55_cl61, rerun of scripts/run_streams_parallel.py).
Statistic: ratio of medians of the daily cloud C_L and the daily Rayleigh C_L over the common period,
with a day-block bootstrap (7-day blocks, each series resampled independently) 95 % CI.

Run:  python validation/eta_study/step4_cloud_eta_variants.py [--workers 16]
Outputs: C:/DATA/Projects/202611_eta_study/step4/ (daily_cloud_variants.csv, ratio_cloud_rayleigh.csv,
         mie_lidar_ratio.csv, eta_tables.csv, fig_step4_timeseries.png, summary.txt)
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
os.environ["ALC_WV_SPECTRUM"] = '{"CL61": [910.55, 0.188]}'
for _v in ("ALC_ETA_TABLES", "ALC_WV_DISABLE", "ALC_DARK_PROFILE"):
    os.environ.pop(_v, None)

import argparse
import copy
import math
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
import extract_incloud_profiles as E  # noqa: E402  (paths, unit/day selection, CAMS guard)

OUT = Path("C:/DATA/Projects/202611_eta_study/step4")
STEP2 = Path("C:/DATA/Projects/202611_eta_study/step2/profiles")
DIAG = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
RAY_DARK = {"0-20000-0-06610_C": DIAG / "diag_v22_l55s008dark"}
RAY_NODARK = {"0-20000-0-06610_C": DIAG / "diag_v22_l55s008", "0-20000-0-10393_C": DIAG / "diag_v22_l55s008"}
RAY_OTHERS = OUT / "rayleigh_v22_l55_cl61"
CLOUD_OPER = DIAG / "diag_v22_l55s008"           # operational-config cloud rows, for the check of variant a
S_THEO = 18.8
SUCCESS = (1, 1.0, 0.5, 0.25)                    # scripts/run_network_calibration._SUCCESS
SIGMA_LN = 0.35
BLOCK_DAYS = 7
N_BOOT = 2000


# =============================================================================
#  eta tables and Mie lidar ratio
# =============================================================================
def pvc_table(a_um, alpha_km):
    from validation.multiple_scattering_eta import CBH_GRID_KM, eta_table
    import calibration.cloud._filters as CF
    # "operational" of this study = the 2026-07 table; since 2026-11 _ETA_CL61 holds the v2 table
    _ETA_CL61 = getattr(CF, "_ETA_PVC55_CL61", CF._ETA_CL61)
    n = _ETA_CL61.shape[0]                       # same 10 nodes (0.25-2.375 km) as the operational table
    return np.column_stack([CBH_GRID_KM[:n], eta_table("CL61", a_um, alpha=alpha_km * 1e-3)[:n]])


def mie_lidar_ratio(r_e_um, wavelength_um=0.91055, m=1.3275 - 4.86e-7j, sigma=SIGMA_LN):
    """Extinction-to-backscatter ratio [sr] of a lognormal water-droplet population (number
    distribution, ln-std sigma) with effective radius r_e. miepython's qback is the radar-style
    backscatter efficiency (4 pi dsigma/dOmega / (pi r^2)) -> S = 4 pi <Q_ext r^2> / <Q_back r^2>.
    Check: r_e 10 um gives 18.9 sr (O'Connor 2004: 18.8 at 905 nm). Water n, k at ~0.91 um:
    Hale & Querry (1973)."""
    import miepython
    r = np.logspace(np.log10(0.05), np.log10(80.0), 6000)
    qe, _qs, qb, _g = miepython.efficiencies_mx(m, 2 * np.pi * r / wavelength_um)
    r_m = r_e_um * np.exp(-2.5 * sigma ** 2)
    n = np.exp(-(np.log(r / r_m)) ** 2 / (2 * sigma ** 2)) / r
    w = n * r ** 2
    return float(4 * np.pi * np.trapezoid(qe * w, r) / np.trapezoid(qb * w, r))


# =============================================================================
#  Worker: one unit-day, every variant
# =============================================================================
def run_day(job):
    """job = (unit dict, YYYYMMDD, {variant: table}) -> list of row dicts."""
    from calibration.cloud import CloudCalConfig
    from calibration.cloud.calibration import (
        liquid_cloud_calibration_from_data, set_defaults, build_cloud_input_from_day)
    import calibration.cloud._filters as CF
    from calibration.flags import cloud_flag
    from calibration.io.instrument_day import load_instrument_day

    u, ds, tables = job
    wmo, ident = u["key"].rsplit("_", 1)
    fp = E.L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"
    if not fp.exists():
        return []
    cfg = set_defaults(CloudCalConfig(
        nc_file=str(fp), instrument="CL61", apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=E.CAMS_FOLDERS, abs_cs_lookup_table=str(E.WV_LUT), cams_folder_fallback="",
        wv_source="cams", era5_cache="", station_latitude=u["lat"], station_longitude=u["lon"],
        average_time_s=0.0, average_range_m=0.0))
    try:
        idd = load_instrument_day([str(fp)], "CL61", E.CAMS_FOLDERS, cams_folder_fallback="",
                                  read_cams=False, build_working=False)
        if idd is None:
            return []
        data0 = build_cloud_input_from_day(idd.slice_to_date(datetime.strptime(ds, "%Y%m%d").date()), cfg)
    except Exception:  # noqa: BLE001
        return []
    oper = CF._ETA_CL61.copy()                   # restored after the loop (the module global the correction reads)
    rows = []
    for var, tab in tables.items():
        CF._ETA_CL61 = np.asarray(tab)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = liquid_cloud_calibration_from_data(copy.deepcopy(data0), cfg)
            n = int(getattr(res, "n_profiles", 0))
            coef, std = float(res.cal_median), float(res.cal_std)
            flag = cloud_flag(n, coef, std)
            cl = float(res.lidar_constant)
            ok = (flag in SUCCESS) and math.isfinite(cl) and cl > 0
            rows.append(dict(key=u["key"], date=int(ds), variant=var, flag=flag, n_profiles=n,
                             cal_median=coef, cal_std=std, C_L=(cl if ok else np.nan)))
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(key=u["key"], date=int(ds), variant=var, flag=-99, n_profiles=0,
                             cal_median=np.nan, cal_std=np.nan, C_L=np.nan, err=str(exc)[:80]))
        finally:
            CF._ETA_CL61 = oper
    return rows


# =============================================================================
#  Reference + statistics
# =============================================================================
def rayleigh_series(key, folder):
    import pandas as pd
    f = folder / key / f"{key}_cal.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d = d[(d.method == "rayleigh") & (d.flag > 0) & (d.cal_value > 0)]
    return d.set_index("date").cal_value.astype(float)


def block_bootstrap_ratio(cloud, ray, rng):
    """Ratio of medians with a 7-day block bootstrap (each series resampled independently)."""
    def blocks(s):
        days = np.array([datetime.strptime(str(int(d)), "%Y%m%d").toordinal() for d in s.index])
        b = (days - days.min()) // BLOCK_DAYS
        return [s.values[b == k] for k in np.unique(b)]
    bc, br = blocks(cloud), blocks(ray)
    out = np.empty(N_BOOT)
    for i in range(N_BOOT):
        c = np.concatenate([bc[j] for j in rng.integers(0, len(bc), len(bc))])
        r = np.concatenate([br[j] for j in rng.integers(0, len(br), len(br))])
        out[i] = np.median(c) / np.median(r)
    return float(np.median(cloud.values) / np.median(ray.values)), np.percentile(out, [2.5, 97.5])


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--plot-only", action="store_true", help="redraw the figure from daily_cloud_variants.csv")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.plot_only:
        units = {k: u for k, u in E.select_units(network=True).items() if u["itype"] == "CL61"}
        figure(pd.read_csv(OUT / "daily_cloud_variants.csv"), units, None)
        return
    log = open(OUT / "summary.txt", "w", encoding="utf-8")

    def say(s=""):
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    from calibration.cloud._filters import _ETA_CL61, _ETA_LEGACY_CL51
    # ---- alpha_k from step 2 (tail slope, CL61 working gates)
    units = {k: u for k, u in E.select_units(network=True).items() if u["itype"] == "CL61"}
    ak = {}
    for k in units:
        p = STEP2 / f"{k}_working.parquet"
        if p.exists():
            d = pd.read_parquet(p, columns=["valid", "alpha_k"])
            ak[k] = d.loc[d.valid, "alpha_k"].to_numpy(float)
    allk = np.concatenate(list(ak.values()))
    a10, a50, a90 = (float(np.nanpercentile(allk, q)) for q in (10, 50, 90))
    say(f"CL61 alpha_k (tail slope, {len(ak)} units, {np.isfinite(allk).sum()} profiles): "
        f"P10 {a10:.1f}  median {a50:.1f}  P90 {a90:.1f} /km")

    # ---- Mie lidar ratio
    mie = pd.DataFrame([dict(r_e_um=r, S_sr=mie_lidar_ratio(r)) for r in (3.5, 4.0, 4.5, 5.5, 8.0, 10.0)])
    re_c = 3.3 / np.exp(-1.5 * SIGMA_LN ** 2)
    S_e = mie_lidar_ratio(re_c)
    mie.loc[len(mie)] = dict(r_e_um=re_c, S_sr=S_e)
    mie.to_csv(OUT / "mie_lidar_ratio.csv", index=False, float_format="%.4g")
    say("Mie lidar ratio (lognormal sigma 0.35, 910.55 nm, water 1.3275-4.86e-7i):\n" + mie.round(3).to_string(index=False))
    say(f"variant e uses S = {S_e:.2f} sr at r_e = {re_c:.2f} um (a_G 3.3 um) -> C_L x {S_e / S_THEO:.4f}")

    # ---- eta tables
    base = {"a": _ETA_CL61, "a_re": pvc_table(5.5, 10.0), "b": _ETA_LEGACY_CL51,
            "c": pvc_table(3.3, a50), "f_lo": pvc_table(3.3, a10), "f_hi": pvc_table(3.3, a90),
            "f_ag25": pvc_table(2.5, a50), "f_ag50": pvc_table(5.0, a50)}
    trow = [dict(variant=v, cbh_km=c, eta=e) for v, t in base.items() for c, e in t]
    unit_tab = {}
    for k in units:
        if k in ak and np.isfinite(ak[k]).any():
            unit_tab[k] = pvc_table(3.3, float(np.nanmedian(ak[k])))
            trow += [dict(variant=f"d:{k}", cbh_km=c, eta=e) for c, e in unit_tab[k]]
    pd.DataFrame(trow).to_csv(OUT / "eta_tables.csv", index=False, float_format="%.5f")
    say("eta at 0.5 / 1 / 2 km: " + "; ".join(
        f"{v} {np.interp(0.5, t[:, 0], t[:, 1]):.3f}/{np.interp(1.0, t[:, 0], t[:, 1]):.3f}/"
        f"{np.interp(2.0, t[:, 0], t[:, 1]):.3f}" for v, t in base.items()))
    say(f"max |a_re - a| on the table nodes: {np.max(np.abs(base['a_re'][:, 1] - _ETA_CL61[:, 1])):.2e}")

    # ---- re-run the calibration per unit-day and variant
    jobs = []
    for k, u in units.items():
        tabs = dict(base)
        if k in unit_tab:
            tabs["d"] = unit_tab[k]
        jobs += [(u, d, tabs) for d in u["days"]]
    jobs.sort(key=lambda j: j[0]["key"] != "0-20000-0-10393_C")      # largest unit first
    say(f"{len(units)} CL61 units, {len(jobs)} unit-days x {len(base) + 1} variants -> {a.workers} processes")
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(run_day, jobs, chunksize=2), 1):
            rows += r
            if i % 100 == 0:
                print(f"  {i}/{len(jobs)} unit-days", flush=True)
    df = pd.DataFrame(rows)
    ce = df[df.variant == "c"].copy()
    ce["variant"] = "e"
    ce["C_L"] = ce["C_L"] * S_e / S_THEO
    df = pd.concat([df, ce], ignore_index=True)
    df.to_csv(OUT / "daily_cloud_variants.csv", index=False, float_format="%.6g")

    # ---- check: variant a vs the operational-config cloud rows (diag_v22_l55s008)
    for k in ("0-20000-0-06610_C", "0-20000-0-10393_C"):
        f = CLOUD_OPER / k / f"{k}_cal.csv"
        if not f.exists():
            continue
        op = pd.read_csv(f)
        op = op[(op.method == "cloud") & (op.cal_value > 0)].set_index("date").cal_value
        mine = df[(df.key == k) & (df.variant == "a")].set_index("date").C_L.dropna()
        j = pd.concat([op.rename("oper"), mine.rename("mine")], axis=1, join="inner")
        say(f"check variant a vs {CLOUD_OPER.name} cloud ({k}): {len(j)} common days, "
            f"median ratio {np.median(j.mine / j.oper):.6f}, max |rel diff| {np.max(np.abs(j.mine / j.oper - 1)):.2e}")

    # ---- ratios cloud / Rayleigh
    rng = np.random.default_rng(0)
    rrow = []
    variants = ["a", "a_re", "b", "c", "d", "e", "f_lo", "f_hi", "f_ag25", "f_ag50"]
    for k, u in units.items():
        refs = []
        if k in RAY_DARK:
            refs.append(("rayleigh+dark", rayleigh_series(k, RAY_DARK[k])))
        if k in RAY_NODARK:
            refs.append(("rayleigh", rayleigh_series(k, RAY_NODARK[k])))
        else:
            refs.append(("rayleigh", rayleigh_series(k, RAY_OTHERS)))
        for ref_name, ray in refs:
            if ray is None or ray.empty:
                continue
            for var in variants:
                cl = df[(df.key == k) & (df.variant == var)].set_index("date").C_L.dropna()
                if cl.empty:
                    continue
                lo, hi = max(cl.index.min(), ray.index.min()), min(cl.index.max(), ray.index.max())
                c2, r2 = cl[(cl.index >= lo) & (cl.index <= hi)], ray[(ray.index >= lo) & (ray.index <= hi)]
                if len(c2) < 5 or len(r2) < 5:
                    continue
                ratio, ci = block_bootstrap_ratio(c2, r2, rng)
                rrow.append(dict(key=k, site=u["site"], reference=ref_name, variant=var, period=f"{lo}-{hi}",
                                 n_cloud_days=len(c2), n_ray_days=len(r2), ratio=ratio, ci_lo=ci[0], ci_hi=ci[1],
                                 cloud_med=float(np.median(c2)), ray_med=float(np.median(r2))))
    rdf = pd.DataFrame(rrow)
    rdf.to_csv(OUT / "ratio_cloud_rayleigh.csv", index=False, float_format="%.4f")
    for key in ("0-20000-0-06610_C", "0-20000-0-10393_C"):
        sub = rdf[rdf.key == key]
        if len(sub):
            say(f"\n== {key} ({sub.site.iloc[0]}): cloud C_L / Rayleigh C_L, ratio of medians [95 % block-bootstrap CI] ==")
            say(sub.pivot_table(index="variant", columns="reference", values="ratio").reindex(variants).round(4).to_string())
            say(sub[["reference", "variant", "period", "n_cloud_days", "n_ray_days", "ratio", "ci_lo", "ci_hi"]]
                .round(4).to_string(index=False))
    oth = rdf[~rdf.key.isin(["0-20000-0-06610_C", "0-20000-0-10393_C"])]
    if len(oth):
        say("\n== other CL61 (secondary, Rayleigh without dark): ratio per variant ==")
        say(oth.pivot_table(index=["key", "site"], columns="variant", values="ratio")[
            [v for v in variants if v in oth.variant.unique()]].round(3).to_string())
        say("median over units: " + oth.groupby("variant").ratio.median().reindex(variants).round(4).to_string())
    figure(df, units, rdf)
    log.close()


def figure(df, units, rdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    keys = ["0-20000-0-06610_C", "0-20000-0-10393_C"]
    fig, axes = plt.subplots(1, 2, figsize=(19, 6))
    for ax, k in zip(axes, keys):
        if k not in units:
            continue
        def ser(var):
            s = df[(df.key == k) & (df.variant == var)].set_index("date").C_L.dropna()
            s.index = pd.to_datetime(s.index.astype(str))
            return s
        band_lo = pd.concat([ser(v) for v in ("f_lo", "f_hi", "f_ag25", "f_ag50")], axis=1).min(1)
        band_hi = pd.concat([ser(v) for v in ("f_lo", "f_hi", "f_ag25", "f_ag50")], axis=1).max(1)
        ax.fill_between(band_lo.index, band_lo, band_hi, color="C2", alpha=0.15, step="mid",
                        label="PVC a_G 3.3 band (alpha P10-P90, a_G 2.5-5)")
        for var, c, lab in (("a", "k", "a operational PVC 5.5 um / 10 /km"), ("b", "C1", "b legacy Hopkin"),
                            ("c", "C2", "c PVC 3.3 um / alpha_k type"), ("e", "C4", "e = c + Mie S")):
            s = ser(var)
            ax.plot(s.index, s.values, ".", ms=4, color=c, alpha=0.6, label=lab)
            ax.plot(s.rolling("30D").median().index, s.rolling("30D").median(), "-", color=c, lw=1.6)
        a_ser = ser("a")
        t0, t1 = a_ser.index.min() - pd.Timedelta(days=7), a_ser.index.max() + pd.Timedelta(days=7)
        for ref, folder, mk in (("Rayleigh + dark", RAY_DARK.get(k), "s"), ("Rayleigh (no dark)", RAY_NODARK.get(k), "^")):
            if folder is None:
                continue
            r = rayleigh_series(k, folder)
            if r is not None and not r.empty:
                r.index = pd.to_datetime(r.index.astype(str))
                r = r[(r.index >= t0) & (r.index <= t1)]           # the compared (cloud) period only
                ax.plot(r.index, r.values, mk, ms=4, mfc="none", color="C3" if "dark" in ref and "no" not in ref else "C0",
                        label=f"{ref} (median {r.median():.3f})")
        lo, hi = np.nanpercentile(a_ser.values, [1, 99])
        ax.set_ylim(0.6 * lo, 1.45 * hi)                          # robust: a few Rayleigh outliers fall outside
        ax.set_title(f"CL61 {units[k]['site']} ({k})")
        ax.set_ylabel("lidar constant C_L (CL61 L1 scale)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Step 4: daily cloud C_L per eta variant vs the CL61 Rayleigh C_L (dots: days; lines: 30-day "
                 "rolling median)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_step4_timeseries.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
