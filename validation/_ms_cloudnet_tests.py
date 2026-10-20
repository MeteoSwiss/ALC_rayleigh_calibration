"""Tests 2 & 3 (Cloudnet measured-microphysics) with PARALLEL per-day file reads.

The Cloudnet matching (one netCDF per site-day, ~1500 files) is embarrassingly parallel:
each (site, date) is matched independently in a worker, then aggregated. Produces:
  fig_test2_cloudnet.png    first-principles etaS with measured droplet radius+extinction,
                            residual height-trend legacy vs PVC (lower = better).
  fig_test3_dropletsize.png observed apparent etaS vs measured droplet effective radius.
Run after _ms_fov_experiment.py. Independent of the Test-1/C_L figures.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validation.multiple_scattering_eta import eta_for_cloud_base  # noqa: E402

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/ms_fov")   # moved from .../ms_fov (2026-10)
REPO_FIGS = Path(__file__).resolve().parents[1] / "doc" / "reports" / "figs_multiple_scattering"
CLOUDNET = Path("D:/Cloudnet_Droplet_radius")   # moved from A: (2026-10)
S_THEO = 18.8
RHO_W = 1.0e6  # g m-3
WORKERS = 30   # home gaming PC has 32 cores

OPTICS = {"CHM15k": (0.23e-3, 0.15e-3, 1064e-9),
          "CL61":   (0.56e-3, 0.28e-3, 910.55e-9),
          "CL31":   (0.83e-3, 0.57e-3, 910e-9)}
SITE_CN = {"Payerne": "payerne", "Lindenberg": "lindenberg", "Palaiseau": "palaiseau"}
COLOR = {"CHM15k": "C2", "CL61": "C3", "CL31": "C0"}


def _legacy_eta(itype, z_km):
    tabs = np.load(OUT / "tables.npz")
    if itype == "CHM15k":
        return np.ones_like(z_km)
    tab = tabs["leg31"] if itype == "CL31" else tabs["leg51"]
    return np.interp(z_km, tab[:, 0], tab[:, 1])


def match_one(job):
    """Match one (site, date) to its Cloudnet der file. job = (site, date, profiles) where
    profiles = {itype: [(time_epoch_s, cbh_m, etaS), ...]}. Returns {itype: [(cbh_km, etaS,
    re_um, alpha), ...]}.

    I/O-lean: the der/N arrays are (2880 x 554); we read ONLY the contiguous time-slice that
    the day's cloud profiles actually span (usually a small fraction of the day) instead of the
    whole array -- the A: drive is the bottleneck, not CPU."""
    from netCDF4 import Dataset
    site, ds, profiles = job
    fp = CLOUDNET / f"{ds}_{SITE_CN[site]}_der.nc"
    if not fp.exists():
        return None
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    try:
        with Dataset(fp) as nc:
            cnt = np.asarray(nc.variables["time"][:], float)     # hours since midnight (small)
            cnh = np.asarray(nc.variables["height"][:], float)   # MSL (small)
            salt = float(np.asarray(nc.variables["altitude"][:]).ravel()[0])
            # which time rows do we actually need?
            hrs = np.array([(t - day0) / 3600.0 for pl in profiles.values() for (t, _, _) in pl])
            if hrs.size == 0:
                return None
            tis = np.searchsorted(cnt, hrs)
            tis = np.clip(tis, 0, cnt.size - 1)
            ti0, ti1 = int(tis.min()), int(tis.max())
            der = np.asarray(nc.variables["der"][ti0:ti1 + 1, :], float)      # only the needed span
            hasN = "N_scaled" in nc.variables
            N = np.asarray(nc.variables["N_scaled"][ti0:ti1 + 1, :], float) if hasN else None
    except Exception:
        return None
    der[~np.isfinite(der) | (der > 1e-3)] = np.nan
    out = defaultdict(list)
    for it, plist in profiles.items():
        for tsec, cbh_m, s in plist:
            hh = (tsec - day0) / 3600.0
            ti = int(np.argmin(np.abs(cnt - hh)))
            if abs(cnt[ti] - hh) > 0.1:
                continue
            r = ti - ti0                                    # row within the read slice
            if r < 0 or r >= der.shape[0]:
                continue
            cloud_msl = cbh_m + salt
            hm = (cnh >= cloud_msl - 30) & (cnh <= cloud_msl + 300)
            re_col = der[r, hm]
            if not np.isfinite(re_col).any():
                continue
            re = float(np.nanmedian(re_col))
            alpha = np.nan
            if N is not None:
                Ncol = N[r, hm]
                if np.isfinite(Ncol).any():
                    Nm = float(np.nanmedian(Ncol))
                    if Nm > 0:
                        lwc = (4.0 / 3.0) * np.pi * RHO_W * Nm * (re ** 3) * 1.4
                        alpha = 3.0 * lwc / (2.0 * RHO_W * re)
            out[it].append((cbh_m / 1000.0, s, re * 1e6, alpha))
    return {it: v for it, v in out.items()}


def build_jobs(rows):
    grp = defaultdict(lambda: defaultdict(list))   # (site,date) -> itype -> list
    for r in rows:
        if "S_noeta" not in r:
            continue
        key = (r["site"], r["date"])
        for t, z, s in zip(r["time"], r["cbh"], r["S_noeta"]):
            grp[key][r["itype"]].append((int(t), float(z), float(s)))
    return [(site, ds, {it: v for it, v in d.items()}) for (site, ds), d in grp.items()]


def main():
    REPO_FIGS.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = ["Payerne", "Lindenberg", "Palaiseau"]
    cache = OUT / "match_cache.npz"
    if cache.exists() and "--rematch" not in sys.argv:
        # instant: figures regenerate from the cached match without touching A: again
        z = np.load(cache, allow_pickle=False)
        match = {s: defaultdict(dict) for s in sites}
        for key in z.files:
            s, it, fld = key.split("__")
            match[s][it][fld] = z[key]
        print(f"loaded match cache ({cache.name}); pass --rematch to rebuild", flush=True)
    else:
        rows = [json.loads(l) for l in open(OUT / "results.jsonl", encoding="utf-8")]
        rows = [r for r in rows if "err" not in r and "empty" not in r]
        jobs = build_jobs(rows)
        print(f"{len(jobs)} site-days to match (I/O-bound -> {WORKERS} THREADS)", flush=True)
        # ThreadPool, not ProcessPool: the work is reading ~1500 small netCDF files. netCDF/HDF5
        # release the GIL during I/O, so threads read concurrently with ZERO spawn cost and no
        # pickling -- the Windows ProcessPool spawn (30x re-import of the calibration stack) left
        # both CPU and disk idle. Threads saturate the disk queue instead.
        import time as _t
        agg = {s: defaultdict(lambda: defaultdict(list)) for s in SITE_CN}
        done = 0
        t0 = _t.time()
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            for job, res in zip(jobs, ex.map(match_one, jobs)):   # map preserves order
                done += 1
                if done % 300 == 0:
                    print(f"  matched {done}/{len(jobs)}  ({done/(_t.time()-t0):.0f}/s)", flush=True)
                if not res:
                    continue
                site = job[0]
                for it, tuples in res.items():
                    for cbh, s, re, al in tuples:
                        agg[site][it]["cbh"].append(cbh); agg[site][it]["etas"].append(s)
                        agg[site][it]["re"].append(re); agg[site][it]["alpha"].append(al)
        match = {s: {it: {k: np.asarray(v, float) for k, v in d.items()} for it, d in sd.items()}
                 for s, sd in agg.items()}
        flat = {f"{s}__{it}__{fld}": arr for s, sd in match.items() for it, d in sd.items()
                for fld, arr in d.items()}
        np.savez(cache, **flat)
        print(f"matched in {_t.time()-t0:.1f}s -> cached {cache.name}", flush=True)

    def _binned(x, y, edges):
        xc, ym = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            m = (x >= a) & (x < b) & np.isfinite(y)
            if m.sum() >= 15:
                xc.append(0.5 * (a + b)); ym.append(np.median(y[m]))
        return np.array(xc), np.array(ym)

    # ---- Test 3: apparent etaS vs measured droplet effective radius ----
    reb = np.array([2, 4, 6, 8, 10, 12, 16])
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
    print("\nTest 3 — apparent etaS vs measured droplet radius (slope over 4-12um, "
          "normalised; PVC predicts DECREASE, legacy predicts flat):")
    for ax, site in zip(axes, sites):
        for it, d in match[site].items():
            if "re" not in d or d["re"].size < 50:
                continue
            band = (d["cbh"] > 0.6) & (d["cbh"] < 1.6)
            re, s = d["re"][band], d["etas"][band]
            xc, ym = _binned(re, s, reb)
            if xc.size < 3:
                continue
            ymn = ym / np.median(ym)
            sl = np.polyfit(xc, ymn, 1)[0]
            print(f"  {site:11s} {it:7s} n={re.size:5d}  slope={sl:+.4f} /um")
            ax.plot(xc, ymn, "o-", color=COLOR.get(it, "C4"), lw=1.8,
                    label=f"{it} (n={re.size})")
        ax.set_title(f"{site} (cloud base 0.6-1.6 km)")
        ax.set_xlabel("Cloudnet droplet eff. radius [$\\mu$m]")
        ax.set_ylabel("normalised apparent $\\eta$S"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Test 3 — apparent $\\eta$S vs Cloudnet-measured droplet effective radius "
                 "(fixed cloud-height band)\nPVC/physics predicts a decrease with droplet size; "
                 "the legacy table has NO droplet dependence", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for dd in (OUT, REPO_FIGS):
        fig.savefig(dd / "fig_test3_dropletsize.png", dpi=150)
    plt.close(fig)
    print("  fig_test3_dropletsize.png")

    # ---- Test 2: first-principles closure with measured microphysics ----
    print("\nTest 2 — first-principles closure with measured microphysics "
          "(|corr(theoretical/observed, height)|; lower=better):")
    print(f"{'site/instr':22s} {'n':>6s} {'legacy':>8s} {'PVC(meas)':>10s}  winner")
    fig, ax = plt.subplots(figsize=(10, 5.6))
    labels, legc, pvcc = [], [], []
    for site in sites:
        for it, d in match[site].items():
            if "re" not in d:
                continue
            m = np.isfinite(d["re"]) & np.isfinite(d["etas"]) & np.isfinite(d["cbh"])
            m &= (d["cbh"] > 0.3) & (d["cbh"] < 2.6)
            if m.sum() < 80:
                continue
            cbh, es, re = d["cbh"][m], d["etas"][m], np.maximum(d["re"][m], 2.0)
            al = d["alpha"][m]
            al = np.where(np.isfinite(al) & (al > 1e-3) & (al < 0.2), al, 10e-3)
            rt, rl, wl = OPTICS[it]
            eta_meas = np.array([eta_for_cloud_base(z * 1000, rt, rl, wl / (np.pi * a * 1e-6),
                                                    alpha=ai, dr=10.0)
                                 for z, a, ai in zip(cbh, re, al)])
            eta_leg = _legacy_eta(it, cbh)
            rP = (eta_meas * S_THEO) / es
            rL = (eta_leg * S_THEO) / es
            cP = abs(np.corrcoef(cbh, rP / np.median(rP))[0, 1])
            cL = abs(np.corrcoef(cbh, rL / np.median(rL))[0, 1])
            win = "PVC" if cP < cL else "legacy"
            print(f"{site+'/'+it:22s} {m.sum():6d} {cL:8.3f} {cP:10.3f}  {win}")
            labels.append(f"{site[:3]}/{it}"); legc.append(cL); pvcc.append(cP)
    x = np.arange(len(labels)); w = 0.38
    ax.bar(x - w / 2, legc, w, color="orange", label="legacy table")
    ax.bar(x + w / 2, pvcc, w, color="red", label="PVC (Cloudnet-measured microphysics)")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("|corr(theoretical/observed $\\eta$S, height)|  (lower = better)")
    ax.set_title("Test 2 — first-principles closure with Cloudnet-measured droplet size + "
                 "extinction (no fitting)\nresidual height-trend after each MS model; lower bar wins",
                 fontweight="bold")
    ax.grid(alpha=0.3, axis="y"); ax.legend()
    fig.tight_layout()
    for dd in (OUT, REPO_FIGS):
        fig.savefig(dd / "fig_test2_cloudnet.png", dpi=150)
    plt.close(fig)
    print("  fig_test2_cloudnet.png")


if __name__ == "__main__":
    main()
