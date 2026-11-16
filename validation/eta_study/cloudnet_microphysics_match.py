"""Step 1d of the eta study: Cloudnet droplet microphysics at the cloud-calibration scenes.

For every accepted liquid-cloud calibration profile ("scene") of the four Vaisala units at the
three Cloudnet sites, take the nearest 30-s Cloudnet `der` profile (within 6 min, as in
validation/_ms_cloudnet_tests.match_one), locate the Cloudnet LIQUID BASE next to the
ceilometer cloud base (cloud base MSL = ceilometer CBH + Cloudnet site altitude), and take the
median microphysics in the lowest 50 / 100 / 150 m of that liquid layer, using ONLY pixels with
der_retrieval_status == 1 (reliable retrieval).

Per layer the scene gets:
  r_e (der_scaled, LWP-scaled, Frisch method 2)   -> the r_e used for a_G and alpha
  r_e (der, fixed-N, Frisch method 1)             -> for comparison only
  N_scaled
  a_G   = r_e * exp(-1.5 sigma^2)                  (lognormal, sigma = 0.35 -> 0.832 r_e)
  alpha = 2 pi N_scaled r_e^2 exp(-3 sigma^2)      (extinction of a lognormal, Q_ext = 2)
The OLD script (_ms_cloudnet_tests.match_one) computed alpha = 2 pi N_scaled der^2 * 1.4, i.e.
a factor 1.4 where the lognormal moment is exp(-3 sigma^2) = 0.69 (x2.02 too large). It is
reported in two ways so the bug can be sized: (i) the old formula on the SAME status==1 layer
pixels, (ii) an exact replica of match_one (window CBH-30..CBH+300 m, every liquid status, and
N_scaled read WITHOUT masking its _FillValue, as the old code did).

Liquid-base definition: liquid pixels = der_retrieval_status 1..3 (droplets diagnosed; status 4
"surrounding ice" marks the ICE pixels next to droplets and never carries a der value).
The lowest liquid gate within +-150 m of the ceilometer cloud base is taken, then followed DOWN
through contiguous liquid gates to the base of that liquid run (and UP to its top). A layer of
depth D holds the gates whose centre lies in [base_edge, base_edge + D), base_edge = centre of
the base gate - dz/2, truncated at the top of the contiguous liquid run.

Scene sources (time, CBH of each accepted profile):
  msfov : old/ms_fov/results.jsonl -- the per-profile list behind _ms_cloudnet_tests (July 2026
          run; CL31 rows are WV-spectrum independent and are what the _l55 re-runs seed verbatim)
  dump  : cloud_profile_dump_l55s008/<key>_profiles.npz -- operational configuration
          (910.55/0.188 nm, PVC eta, 0.4 deg CAMS); Payerne B/C and Lindenberg C only.

Run:  python validation/eta_study/cloudnet_microphysics_match.py [--workers 12] [--source msfov,dump]
Outputs (C:/DATA/Projects/202611_eta_study/step1d/): scenes_<source>.csv, summary_<source>.csv,
summary_<source>.md, fig_microphysics_<source>.png.  Read-only on every input.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import time as time_mod
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

CLOUDNET = Path("D:/Cloudnet_Droplet_radius")
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
MSFOV_JSONL = DATA / "old" / "ms_fov" / "results.jsonl"
DUMP_DIR = DATA / "cloud_profile_dump_l55s008"
OUT = Path("C:/DATA/Projects/202611_eta_study/step1d")

SIGMA_LN = 0.35
AG_FACTOR = float(np.exp(-1.5 * SIGMA_LN ** 2))       # 0.832
ALPHA_FACTOR = float(np.exp(-3.0 * SIGMA_LN ** 2))    # 0.692
OLD_FACTOR = 1.4                                       # the bug in _ms_cloudnet_tests.match_one
RHO_W = 1.0e6                                          # g m-3 (old-script replica only)
MAX_DT_H = 0.1                                         # 6 min
BASE_SEARCH_M = 150.0                                  # liquid gate searched within +-150 m of CBH
LAYERS_M = (50, 100, 150)
CN_LAST_DAY = "20260630"                               # last day of the Cloudnet collection

# unit key -> (site, Cloudnet site name, instrument type, distorted flag, earliest date)
UNITS = {
    "0-20000-0-06610_C": ("Payerne", "payerne", "CL61", False, None),
    "0-20000-0-06610_B": ("Payerne", "payerne", "CL31", True, None),    # distorted unit
    "0-250-1001-07151_A": ("Palaiseau", "palaiseau", "CL31", False, None),
    "0-20000-0-10393_C": ("Lindenberg", "lindenberg", "CL61", False, "20250501"),  # 2x hardware step
}


# =============================================================================
#  Scene sources
# =============================================================================
def scenes_from_msfov(path: Path):
    """(key, date, t_epoch_s, cbh_m, s_app) from the ms_fov per-profile results."""
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            if r.get("key") not in UNITS or "time" not in r:
                continue
            for t, z, s in zip(r["time"], r["cbh"], r["S_noeta"]):
                out.append((r["key"], r["date"], int(t), float(z), float(s)))
    return out


def scenes_from_dump(folder: Path):
    """(key, date, t_epoch_s, cbh_m, s_app) from the operational per-profile dump."""
    out = []
    for key in UNITS:
        fp = folder / f"{key}_profiles.npz"
        if not fp.exists():
            continue
        z = np.load(fp, allow_pickle=True)
        t_s = (z["time_ns"] // 1_000_000_000).astype(np.int64)
        for d, t, cbh, s in zip(z["date"], t_s, z["cbh_m"], z["s_apparent"]):
            out.append((key, str(int(d)), int(t), float(cbh), float(s)))
    return out


def keep_scene(key: str, ds: str) -> bool:
    earliest = UNITS[key][4]
    return (earliest is None or ds >= earliest) and ds <= CN_LAST_DAY


# =============================================================================
#  Worker: one (Cloudnet site, date) -> per-scene microphysics
# =============================================================================
def _filled(var, sl, fill=np.nan, dtype=float):
    """Read a slice with its _FillValue turned into `fill` (the masked-array trap)."""
    x = var[sl]
    return np.ma.filled(np.ma.asarray(x, dtype=dtype), fill)


def _layer_stats(idx, st, der, ders, nsc):
    """Medians over the status==1 pixels of one layer (gate indices `idx`)."""
    res = dict(n_liq=int(idx.size), n_st1=0, re_scaled=np.nan, re_der=np.nan, N_scaled=np.nan,
               alpha_new=np.nan, alpha_der_corr=np.nan, alpha_old_same=np.nan)
    if idx.size == 0:
        return res
    ok1 = st[idx] == 1
    res["n_st1"] = int(ok1.sum())
    if not ok1.any():
        return res
    rs = ders[idx][ok1]
    rd = der[idx][ok1]
    nn = nsc[idx][ok1]
    rs = np.where((rs > 0) & (rs < 1e-3), rs, np.nan)
    rd = np.where((rd > 0) & (rd < 1e-3), rd, np.nan)
    nn = np.where(nn > 0, nn, np.nan)
    with np.errstate(all="ignore"):
        if np.isfinite(rs).any():
            res["re_scaled"] = float(np.nanmedian(rs))
        if np.isfinite(rd).any():
            res["re_der"] = float(np.nanmedian(rd))
        if np.isfinite(nn).any():
            res["N_scaled"] = float(np.nanmedian(nn))
        a_new = 2.0 * np.pi * nn * rs ** 2 * ALPHA_FACTOR
        a_dc = 2.0 * np.pi * nn * rd ** 2 * ALPHA_FACTOR
        a_old = 2.0 * np.pi * nn * rd ** 2 * OLD_FACTOR
        for name, arr in (("alpha_new", a_new), ("alpha_der_corr", a_dc), ("alpha_old_same", a_old)):
            if np.isfinite(arr).any():
                res[name] = float(np.nanmedian(arr))
    return res


def match_day(job):
    """job = (cn_site, ds, [(row_id, t_epoch_s, cbh_m), ...]) -> list of (row_id, dict)."""
    from netCDF4 import Dataset
    cn_site, ds, scenes = job
    fp = CLOUDNET / f"{ds}_{cn_site}_der.nc"
    if not fp.exists():
        return [(rid, dict(status="nofile")) for rid, _, _ in scenes]
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    hrs = np.array([(t - day0) / 3600.0 for _, t, _ in scenes])
    try:
        with Dataset(fp) as nc:
            cnt = np.asarray(nc.variables["time"][:], float)
            h = np.asarray(nc.variables["height"][:], float)
            salt = float(np.asarray(nc.variables["altitude"][:]).ravel()[0])
            tis = np.clip(np.searchsorted(cnt, hrs), 1, cnt.size - 1)
            tis = np.where(np.abs(cnt[tis - 1] - hrs) <= np.abs(cnt[tis] - hrs), tis - 1, tis)
            ti0, ti1 = int(tis.min()), int(tis.max())
            sl = (slice(ti0, ti1 + 1), slice(None))
            st = _filled(nc.variables["der_retrieval_status"], sl, fill=0, dtype=int)
            der = _filled(nc.variables["der"], sl)
            ders = _filled(nc.variables["der_scaled"], sl)
            nsc = _filled(nc.variables["N_scaled"], sl)
            # old-script replica: np.asarray keeps the 9.97e36 _FillValue of N_scaled
            der_old = np.asarray(nc.variables["der"][sl], float)
            n_old = np.asarray(nc.variables["N_scaled"][sl], float)
    except Exception as exc:  # noqa: BLE001 - a corrupt file must not kill the run
        return [(rid, dict(status=f"readerr {type(exc).__name__}")) for rid, _, _ in scenes]
    der_old[~np.isfinite(der_old) | (der_old > 1e-3)] = np.nan
    dz = float(np.median(np.diff(h)))

    out = []
    for (rid, _t, cbh), hh, ti in zip(scenes, hrs, tis):
        dt_s = (cnt[ti] - hh) * 3600.0
        if abs(cnt[ti] - hh) > MAX_DT_H:
            out.append((rid, dict(status="notime", dt_s=dt_s)))
            continue
        r = ti - ti0
        cloud_msl = cbh + salt
        rec = dict(status="matched", dt_s=dt_s, site_alt=salt, dz=dz, cloud_msl=cloud_msl)

        # --- exact replica of the old match_one value ---
        hm = (h >= cloud_msl - 30) & (h <= cloud_msl + 300)
        re_col = der_old[r, hm]
        if np.isfinite(re_col).any():
            re_o = float(np.nanmedian(re_col))
            rec["old_re_um"] = re_o * 1e6
            ncol = n_old[r, hm]
            if np.isfinite(ncol).any():
                nm = float(np.nanmedian(ncol))
                rec["old_N_fill_contaminated"] = bool(nm > 1e30)
                if nm > 0:
                    lwc = (4.0 / 3.0) * np.pi * RHO_W * nm * (re_o ** 3) * OLD_FACTOR
                    rec["old_alpha_script"] = 3.0 * lwc / (2.0 * RHO_W * re_o)

        # --- Cloudnet liquid base next to the ceilometer cloud base ---
        # status 4 ("surrounding ice") flags the ICE pixels next to droplets and never carries a
        # der value -> liquid = status 1..3 (checked: der finite on every 1-3 pixel, none on 4)
        liq = (st[r] >= 1) & (st[r] <= 3)
        cand = np.where(liq & (np.abs(h - cloud_msl) <= BASE_SEARCH_M))[0]
        if cand.size == 0:
            rec["liquid"] = False
            out.append((rid, rec))
            continue
        jb = int(cand.min())
        while jb > 0 and liq[jb - 1]:
            jb -= 1
        jt = jb
        while jt + 1 < liq.size and liq[jt + 1]:
            jt += 1
        base_edge = h[jb] - dz / 2.0
        rec.update(liquid=True, base_msl=float(base_edge), base_offset=float(base_edge - cloud_msl),
                   liq_depth=float((jt - jb + 1) * dz), base_status=int(st[r, jb]))
        for D in LAYERS_M:
            idx = np.arange(jb, jt + 1)
            idx = idx[(h[idx] - base_edge) < D]
            ls = _layer_stats(idx, st[r], der[r], ders[r], nsc[r])
            for k, v in ls.items():
                rec[f"{k}_{D}"] = v
        out.append((rid, rec))
    return out


# =============================================================================
#  Driver
# =============================================================================
def run_source(name: str, scenes, workers: int):
    import pandas as pd
    df = pd.DataFrame(scenes, columns=["key", "date", "t_epoch", "cbh_m", "s_app"])
    df = df[[keep_scene(k, d) for k, d in zip(df["key"], df["date"])]].reset_index(drop=True)
    df["site"] = [UNITS[k][0] for k in df["key"]]
    df["itype"] = [UNITS[k][2] for k in df["key"]]
    df["distorted_unit"] = [UNITS[k][3] for k in df["key"]]
    df["time_utc"] = pd.to_datetime(df["t_epoch"], unit="s")
    grp = defaultdict(list)
    for rid, k, d, t, z in zip(df.index, df["key"], df["date"], df["t_epoch"], df["cbh_m"]):
        grp[(UNITS[k][1], d)].append((rid, int(t), float(z)))
    jobs = [(s, d, v) for (s, d), v in sorted(grp.items())]
    print(f"[{name}] {len(df)} scenes on {len(jobs)} site-days -> {workers} processes", flush=True)
    t0 = time_mod.time()
    recs = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for k, res in enumerate(ex.map(match_day, jobs, chunksize=4), 1):
            for rid, rec in res:
                recs[rid] = rec
            if k % 200 == 0:
                print(f"  {k}/{len(jobs)} site-days ({time_mod.time() - t0:.0f} s)", flush=True)
    print(f"[{name}] matched in {time_mod.time() - t0:.1f} s", flush=True)
    extra = pd.DataFrame.from_dict(recs, orient="index").sort_index()
    df = df.join(extra)
    # physical units for the outputs: um, cm-3, km-1
    df["old_alpha_script_km"] = df["old_alpha_script"] * 1e3 if "old_alpha_script" in df else np.nan
    for D in LAYERS_M:
        df[f"re_scaled_um_{D}"] = df[f"re_scaled_{D}"] * 1e6
        df[f"re_der_um_{D}"] = df[f"re_der_{D}"] * 1e6
        df[f"aG_um_{D}"] = df[f"re_scaled_um_{D}"] * AG_FACTOR
        df[f"N_scaled_cm3_{D}"] = df[f"N_scaled_{D}"] * 1e-6
        df[f"alpha_km_{D}"] = df[f"alpha_new_{D}"] * 1e3
        df[f"alpha_der_corr_km_{D}"] = df[f"alpha_der_corr_{D}"] * 1e3
        df[f"alpha_old_same_km_{D}"] = df[f"alpha_old_same_{D}"] * 1e3
        df = df.drop(columns=[f"re_scaled_{D}", f"re_der_{D}", f"N_scaled_{D}", f"alpha_new_{D}",
                              f"alpha_der_corr_{D}", f"alpha_old_same_{D}"])
    df = df.drop(columns=["old_alpha_script"], errors="ignore")
    return df


GROUPS = [
    ("Payerne CL61", lambda d: (d.site == "Payerne") & (d.itype == "CL61")),
    ("Payerne CL31 (distorted unit)", lambda d: (d.site == "Payerne") & (d.itype == "CL31")),
    ("Palaiseau CL31", lambda d: (d.site == "Palaiseau") & (d.itype == "CL31")),
    ("Lindenberg CL61 (>=2025-05)", lambda d: (d.site == "Lindenberg") & (d.itype == "CL61")),
    ("CL31 all", lambda d: d.itype == "CL31"),
    ("CL61 all", lambda d: d.itype == "CL61"),
    ("Vaisala all", lambda d: d.itype.isin(["CL31", "CL61"])),
    ("Vaisala excl. Payerne CL31", lambda d: ~d.distorted_unit.astype(bool)),
]


def _q(x, p):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if x.size else np.nan


def summarise(df):
    import pandas as pd
    rows = []
    for gname, sel in GROUPS:
        g = df[sel(df)]
        if g.empty:
            continue
        m = g[g["status"] == "matched"]
        liq = m[m["liquid"] == True]  # noqa: E712 - pandas object column
        row = {"group": gname, "n_scenes": len(g), "n_days": g["date"].nunique(),
               "n_matched": len(m), "n_days_matched": m["date"].nunique(),
               "n_liquid": len(liq), "frac_liquid": len(liq) / max(len(m), 1),
               "base_offset_med_m": _q(liq["base_offset"], 50),
               "liq_depth_med_m": _q(liq["liq_depth"], 50),
               "frac_base_status1": float((liq["base_status"] == 1).mean()) if len(liq) else np.nan}
        old = m["old_alpha_script_km"]
        contam = m.get("old_N_fill_contaminated")
        row["old_script_n"] = int(np.isfinite(old).sum())
        row["old_script_frac_fill"] = float(contam.fillna(False).astype(bool).sum() / max(np.isfinite(old).sum(), 1)) if contam is not None else np.nan
        clean = old[(old > 0) & (old < 1e3)]
        row["old_script_alpha_med_km"] = _q(clean, 50)
        row["old_script_re_med_um"] = _q(m["old_re_um"], 50)
        for D in LAYERS_M:
            s1 = liq[liq[f"n_st1_{D}"] > 0]
            row[f"frac_st1_{D}_of_matched"] = len(s1) / max(len(m), 1)
            row[f"frac_st1_{D}_of_liquid"] = len(s1) / max(len(liq), 1)
            row[f"n_st1_{D}"] = len(s1)
            row[f"n_days_st1_{D}"] = s1["date"].nunique()
            for var in ("aG_um", "alpha_km"):
                col = s1[f"{var}_{D}"]
                row[f"{var}_{D}_med"] = _q(col, 50)
                row[f"{var}_{D}_p10"] = _q(col, 10)
                row[f"{var}_{D}_p90"] = _q(col, 90)
                daymed = s1.groupby("date")[f"{var}_{D}"].median()
                row[f"{var}_{D}_daymed"] = _q(daymed, 50)
            row[f"re_scaled_um_{D}_med"] = _q(s1[f"re_scaled_um_{D}"], 50)
            row[f"re_der_um_{D}_med"] = _q(s1[f"re_der_um_{D}"], 50)
            row[f"N_scaled_cm3_{D}_med"] = _q(s1[f"N_scaled_cm3_{D}"], 50)
            row[f"alpha_der_corr_km_{D}_med"] = _q(s1[f"alpha_der_corr_km_{D}"], 50)
            row[f"alpha_old_same_km_{D}_med"] = _q(s1[f"alpha_old_same_km_{D}"], 50)
        rows.append(row)
    return pd.DataFrame(rows)


def write_markdown(summ, name, path):
    def f(x, fmt="{:.1f}"):
        return "" if not np.isfinite(x) else fmt.format(x)
    L = [f"# Step 1d - Cloudnet microphysics at calibration scenes (source: {name})", ""]
    L.append("| group | scenes | days | matched | liquid near CBH | base-CBH (m) | liquid depth (m) |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in summ.iterrows():
        L.append(f"| {r.group} | {r.n_scenes} | {r.n_days} | {r.n_matched} | "
                 f"{r.n_liquid} ({100 * r.frac_liquid:.0f} %) | {f(r.base_offset_med_m, '{:+.0f}')} | "
                 f"{f(r.liq_depth_med_m, '{:.0f}')} |")
    for D in LAYERS_M:
        L += ["", f"## Lowest {D} m above the Cloudnet liquid base (status == 1 pixels only)", "",
              "| group | n (days) | frac status1 | a_G um med [P10-P90] | alpha /km med [P10-P90] | "
              "a_G day-med | alpha day-med | r_e scaled | r_e der | N cm-3 | alpha old-formula same px |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for _, r in summ.iterrows():
            L.append(
                f"| {r.group} | {r[f'n_st1_{D}']} ({r[f'n_days_st1_{D}']}) | "
                f"{100 * r[f'frac_st1_{D}_of_matched']:.0f} % | "
                f"{f(r[f'aG_um_{D}_med'], '{:.2f}')} [{f(r[f'aG_um_{D}_p10'])}-{f(r[f'aG_um_{D}_p90'])}] | "
                f"{f(r[f'alpha_km_{D}_med'])} [{f(r[f'alpha_km_{D}_p10'])}-{f(r[f'alpha_km_{D}_p90'])}] | "
                f"{f(r[f'aG_um_{D}_daymed'], '{:.2f}')} | {f(r[f'alpha_km_{D}_daymed'])} | "
                f"{f(r[f're_scaled_um_{D}_med'], '{:.2f}')} | {f(r[f're_der_um_{D}_med'], '{:.2f}')} | "
                f"{f(r[f'N_scaled_cm3_{D}_med'], '{:.0f}')} | {f(r[f'alpha_old_same_km_{D}_med'])} |")
    L += ["", "## Old-script replica (window CBH-30..+300 m, all statuses, der + unmasked N_scaled, x1.4)", "",
          "| group | n | frac N fill-contaminated | r_e der med um | alpha med /km (fill excluded) |",
          "|---|---|---|---|---|"]
    for _, r in summ.iterrows():
        L.append(f"| {r.group} | {r.old_script_n} | {100 * r.old_script_frac_fill:.0f} % | "
                 f"{f(r.old_script_re_med_um, '{:.2f}')} | {f(r.old_script_alpha_med_km)} |")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def plot(df, name, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    units = [("Payerne CL61", "C3"), ("Lindenberg CL61 (>=2025-05)", "C1"),
             ("Palaiseau CL31", "C0"), ("Payerne CL31 (distorted unit)", "C7")]
    sel = dict(GROUPS)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    for lab, c in units:
        g = df[sel[lab](df) & (df["status"] == "matched")]
        g = g[g["liquid"] == True]  # noqa: E712
        g = g[g["n_st1_100"] > 0]
        if g.empty:
            continue
        ls = "--" if "distorted" in lab else "-"
        for ax, col, bins in ((axes[0], "aG_um_100", np.arange(0, 16, 0.4)),
                              (axes[1], "alpha_km_100", np.logspace(-2, 2.5, 46)),
                              (axes[2], "re_der_um_100", np.arange(0, 16, 0.4))):
            x = g[col].to_numpy(float)
            x = x[np.isfinite(x)]
            if x.size:
                # fraction per bin (not a density): log-spaced alpha bins would otherwise
                # exaggerate the low tail
                ax.hist(np.clip(x, bins[0], bins[-1]), bins=bins, weights=np.full(x.size, 1.0 / x.size),
                        histtype="step", lw=2, ls=ls, color=c,
                        label=f"{lab} (n={x.size}, med {np.median(x):.2f})")
    axes[0].axvline(5.5, color="k", ls=":", lw=1.5, label="PVC operational a$_G$ = 5.5 $\\mu$m")
    axes[1].axvline(10.0, color="k", ls=":", lw=1.5, label="PVC operational $\\alpha$ = 10 km$^{-1}$")
    axes[1].set_xscale("log")
    axes[0].set_xlabel("a$_G$ = 0.832 r$_e$(der_scaled) [$\\mu$m]")
    axes[1].set_xlabel("$\\alpha$ = 2$\\pi$ N$_{scaled}$ r$_e^2$ e$^{-3\\sigma^2}$ [km$^{-1}$]")
    axes[2].set_xlabel("r$_e$ from der (fixed-N, comparison only) [$\\mu$m]")
    for ax in axes:
        ax.set_ylabel("fraction of scenes per bin (end bins hold the clipped tails)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5)
    fig.suptitle(f"Cloudnet microphysics in the lowest 100 m above the liquid base at the "
                 f"calibration scenes (status==1 pixels; scenes: {name})", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--source", default="msfov,dump")
    ap.add_argument("--msfov", default=str(MSFOV_JSONL))
    ap.add_argument("--dump", default=str(DUMP_DIR))
    ap.add_argument("--plot-only", action="store_true", help="redraw the figures from scenes_<source>.csv")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t_all = time_mod.time()
    if a.plot_only:
        import pandas as pd
        for name in a.source.split(","):
            plot(pd.read_csv(OUT / f"scenes_{name}.csv", low_memory=False), name,
                 OUT / f"fig_microphysics_{name}.png")
        return
    for name in a.source.split(","):
        scenes = scenes_from_msfov(Path(a.msfov)) if name == "msfov" else scenes_from_dump(Path(a.dump))
        df = run_source(name, scenes, a.workers)
        df.to_csv(OUT / f"scenes_{name}.csv", index=False, float_format="%.5g")
        summ = summarise(df)
        summ.to_csv(OUT / f"summary_{name}.csv", index=False, float_format="%.5g")
        write_markdown(summ, name, OUT / f"summary_{name}.md")
        plot(df, name, OUT / f"fig_microphysics_{name}.png")
        print((OUT / f"summary_{name}.md").read_text(encoding="utf-8"), flush=True)
    print(f"total runtime {time_mod.time() - t_all:.1f} s", flush=True)


if __name__ == "__main__":
    main()
